import argparse
import array
import collections
import gc
import hashlib
import mmap
import os
import pickle as pkl
import signal
import struct
import sys
import tempfile
import time

from barcode_store import BarcodeStore


_NEED_SWAP = sys.byteorder == "little"
_EMPTY_ARRAY = array.array("Q")
_MIB = 1024 * 1024

DEFAULT_CANDIDATE_CACHE_MIB = 256
DEFAULT_HEAD_CACHE_MIB = 64
DEFAULT_PROGRESS_EVERY = 100000
CHECKPOINT_SCHEMA_VERSION = 1
DEFAULT_BASE_DIR = "path/to/your/project_directory"  # Replace with your MTPLite project directory

# head_index.idx record: 11 x uint32 kmer | uint64 data-offset | uint32 count
_KMER_LEN = 11
_KMER_BYTES = _KMER_LEN * 4
_IDX_RECORD = struct.Struct(f">{_KMER_LEN}IQI")
_IDX_ENTRY_SIZE = _IDX_RECORD.size
_KMER_FMT = struct.Struct(f">{_KMER_LEN}I")
_TAIL_FMT = struct.Struct(">QI")

# Conservative allowance for an OrderedDict node and its object references.
_CACHE_ENTRY_OVERHEAD = 224
_MAX_PARENT_INDEX = (1 << 32) - 1


def log(message):
    """Emit a timestamped message immediately, including under nohup."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def load_pickle(in_file):
    log(f"Loading {in_file}...")
    with open(in_file, "rb") as stream:
        return pkl.load(stream)


def atomic_write_pickle(data, out_file):
    """Atomically replace *out_file* with a protocol-4 pickle."""
    destination = os.path.abspath(os.path.expanduser(os.fspath(out_file)))
    directory = os.path.dirname(destination)
    os.makedirs(directory, exist_ok=True)

    prefix = f".{os.path.basename(destination)}."
    fd, temporary = tempfile.mkstemp(
        dir=directory,
        prefix=prefix,
        suffix=".tmp",
    )
    try:
        with os.fdopen(fd, "wb") as stream:
            fd = None
            pkl.dump(data, stream, protocol=4)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    except BaseException:
        if fd is not None:
            os.close(fd)
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_pickle(data, out_file):
    log(f"Writing results to {out_file}...")
    atomic_write_pickle(data, out_file)


class ByteLRUCache:
    """An LRU whose caller accounts for each entry in bytes."""

    __slots__ = (
        "max_bytes",
        "current_bytes",
        "_entries",
        "hits",
        "misses",
        "evictions",
        "puts",
        "rejections",
    )

    def __init__(self, max_bytes):
        max_bytes = int(max_bytes)
        if max_bytes < 0:
            raise ValueError("max_bytes must be non-negative")
        self.max_bytes = max_bytes
        self.current_bytes = 0
        self._entries = collections.OrderedDict()
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self.puts = 0
        self.rejections = 0

    def get(self, key):
        try:
            value, _size = self._entries[key]
        except KeyError:
            self.misses += 1
            return None
        self._entries.move_to_end(key)
        self.hits += 1
        return value

    def put(self, key, value, size):
        size = int(size)
        if size < 0:
            raise ValueError("cache entry size must be non-negative")
        if self.max_bytes == 0 or size > self.max_bytes:
            self.rejections += 1
            return False

        previous = self._entries.pop(key, None)
        if previous is not None:
            self.current_bytes -= previous[1]

        while self._entries and self.current_bytes + size > self.max_bytes:
            _, (_, evicted_size) = self._entries.popitem(last=False)
            self.current_bytes -= evicted_size
            self.evictions += 1

        self._entries[key] = (value, size)
        self.current_bytes += size
        self.puts += 1
        return True

    def clear(self):
        self._entries.clear()
        self.current_bytes = 0

    def stats(self):
        return {
            "hits": self.hits,
            "misses": self.misses,
            "evictions": self.evictions,
            "puts": self.puts,
            "rejections": self.rejections,
            "entries": len(self._entries),
            "current_bytes": self.current_bytes,
            "max_bytes": self.max_bytes,
        }


def _advise_random(mm):
    if mm is None:
        return
    madvise = getattr(mm, "madvise", None)
    advice = getattr(mmap, "MADV_RANDOM", None)
    if madvise is None or advice is None:
        return
    try:
        madvise(advice)
    except (OSError, ValueError):
        pass


def _advise_barcode_store_random(store):
    _advise_random(getattr(store, "_idx_mm", None))
    _advise_random(getattr(store, "_data_mm", None))


class HeadIndexQuery:
    """Binary-search lookup over the sorted head-index flat files.

    Only compact (offset, count) metadata is cached. Posting arrays are decoded
    for the caller and immediately become reclaimable.
    """

    __slots__ = (
        "_idx_mm",
        "_data_mm",
        "_idx_fd",
        "_data_fd",
        "_n_records",
        "_data_size",
        "_cache",
    )

    def __init__(
        self,
        idx_path,
        data_path,
        cache_bytes=DEFAULT_HEAD_CACHE_MIB * _MIB,
    ):
        self._idx_fd = None
        self._data_fd = None
        self._idx_mm = None
        self._data_mm = None
        self._data_size = 0
        self._cache = ByteLRUCache(cache_bytes)

        try:
            self._idx_fd = os.open(os.fspath(idx_path), os.O_RDONLY)
            idx_size = os.fstat(self._idx_fd).st_size
            if idx_size == 0 or idx_size % _IDX_ENTRY_SIZE:
                raise ValueError(
                    "head index size must be a non-zero multiple of "
                    f"{_IDX_ENTRY_SIZE} bytes"
                )
            self._idx_mm = mmap.mmap(
                self._idx_fd,
                idx_size,
                access=mmap.ACCESS_READ,
            )
            self._n_records = idx_size // _IDX_ENTRY_SIZE

            self._data_fd = os.open(os.fspath(data_path), os.O_RDONLY)
            self._data_size = os.fstat(self._data_fd).st_size
            if self._data_size:
                self._data_mm = mmap.mmap(
                    self._data_fd,
                    self._data_size,
                    access=mmap.ACCESS_READ,
                )
            _advise_random(self._idx_mm)
            _advise_random(self._data_mm)
        except BaseException:
            self.close()
            raise

    @staticmethod
    def _metadata_cache_size(packed_key, metadata):
        return (
            sys.getsizeof(packed_key)
            + sys.getsizeof(metadata)
            + sum(sys.getsizeof(value) for value in metadata)
            + _CACHE_ENTRY_OVERHEAD
        )

    def _find_metadata(self, packed_key):
        idx_mm = self._idx_mm
        lo, hi = 0, self._n_records
        while lo < hi:
            mid = (lo + hi) >> 1
            offset = mid * _IDX_ENTRY_SIZE
            mid_key = idx_mm[offset : offset + _KMER_BYTES]
            if mid_key < packed_key:
                lo = mid + 1
            elif mid_key > packed_key:
                hi = mid
            else:
                return _TAIL_FMT.unpack_from(idx_mm, offset + _KMER_BYTES)
        return (-1, 0)

    def get_packed(self, packed_key):
        if not isinstance(packed_key, bytes):
            packed_key = bytes(packed_key)
        if len(packed_key) != _KMER_BYTES:
            raise ValueError(
                f"packed head-index key must contain {_KMER_BYTES} bytes"
            )

        metadata = self._cache.get(packed_key)
        if metadata is None:
            metadata = self._find_metadata(packed_key)
            self._cache.put(
                packed_key,
                metadata,
                self._metadata_cache_size(packed_key, metadata),
            )

        data_offset, count = metadata
        if data_offset < 0 or count == 0:
            return _EMPTY_ARRAY

        data_end = data_offset + count * 8
        if self._data_mm is None or data_end > self._data_size:
            raise ValueError(
                "head-index posting range lies outside head_index.data"
            )

        hits = array.array("Q")
        view = memoryview(self._data_mm)[data_offset:data_end]
        try:
            hits.frombytes(view)
        finally:
            view.release()
        if _NEED_SWAP:
            hits.byteswap()
        return hits

    def get(self, kmer_tuple):
        return self.get_packed(_KMER_FMT.pack(*kmer_tuple))

    def cache_stats(self):
        return self._cache.stats()

    def close(self):
        cache = getattr(self, "_cache", None)
        if cache is not None:
            cache.clear()

        for mm_attr in ("_idx_mm", "_data_mm"):
            mm = getattr(self, mm_attr, None)
            if mm is not None:
                mm.close()
                setattr(self, mm_attr, None)

        for fd_attr in ("_idx_fd", "_data_fd"):
            fd = getattr(self, fd_attr, None)
            if fd is not None:
                os.close(fd)
                setattr(self, fd_attr, None)


def read_memory_status(status_path="/proc/self/status"):
    """Return selected Linux process-memory counters in KiB."""
    wanted = {"RssAnon", "RssFile", "VmRSS", "VmHWM"}
    counters = {}
    try:
        with open(status_path, "r") as stream:
            for line in stream:
                name, separator, remainder = line.partition(":")
                if not separator or name not in wanted:
                    continue
                fields = remainder.split()
                if fields:
                    counters[name] = int(fields[0])
    except (OSError, ValueError):
        return {}
    return counters


def _format_mib(kib):
    return f"{kib / 1024.0:.1f}MiB"


def memory_summary():
    counters = read_memory_status()
    if not counters:
        return "memory=unavailable"
    return " ".join(
        f"{name}={_format_mib(counters[name])}"
        for name in ("RssAnon", "RssFile", "VmRSS", "VmHWM")
        if name in counters
    )


def _cache_summary(name, stats):
    return (
        f"{name}=entries:{stats['entries']:,},"
        f"bytes:{stats['current_bytes'] / _MIB:.1f}/{stats['max_bytes'] / _MIB:.1f}MiB,"
        f"hits:{stats['hits']:,},misses:{stats['misses']:,},"
        f"evictions:{stats['evictions']:,},rejections:{stats['rejections']:,}"
    )


class ProjectAssembler:
    def __init__(
        self,
        barcode_query,
        head_query,
        int_to_rid=None,
        k_mer_size=15,
        end_fraction=0.15,
        end_cap=2500,
        query_stride=2,
        candidate_cache_mib=DEFAULT_CANDIDATE_CACHE_MIB,
        progress_every=DEFAULT_PROGRESS_EVERY,
    ):
        if not 0 <= int(candidate_cache_mib):
            raise ValueError("candidate_cache_mib must be non-negative")
        if int(progress_every) < 0:
            raise ValueError("progress_every must be non-negative")
        self.barcode_query = barcode_query
        self.head_query = head_query
        self.int_to_rid = int_to_rid or {}
        self.k_mer_size = int(k_mer_size)
        self.end_fraction = float(end_fraction)
        self.end_cap = int(end_cap)
        self.query_stride = max(1, int(query_stride))
        self.progress_every = int(progress_every)
        self._candidate_cache = ByteLRUCache(
            int(candidate_cache_mib) * _MIB
        )

        if array.array("Q").itemsize != 8:
            raise RuntimeError("this assembler requires 64-bit array('Q')")
        if array.array("I").itemsize != 4:
            raise RuntimeError("this assembler requires 32-bit array('I')")

    def _end_window(self, length):
        return min(int(length * self.end_fraction), self.end_cap)

    def _get_candidate_ids_from_suffix(self, barcodes):
        kmer_size = self.k_mer_size
        length = len(barcodes)
        window = self._end_window(length)
        if window < kmer_size:
            return array.array("Q")

        candidates = set()
        candidates_update = candidates.update
        start_scan = length - window
        get_packed = getattr(self.head_query, "get_packed", None)
        if get_packed is not None:
            # Convert the suffix window to big-endian bytes once. Each lookup
            # then copies only its 44-byte key instead of allocating an array
            # slice plus 11 boxed Python integers.
            packed_barcodes = barcodes[start_scan:]
            if not (
                isinstance(packed_barcodes, array.array)
                and packed_barcodes.typecode == "I"
            ):
                packed_barcodes = array.array("I", packed_barcodes)
            if _NEED_SWAP:
                packed_barcodes.byteswap()
            packed_view = memoryview(packed_barcodes).cast("B")
            try:
                last_relative_start = window - kmer_size
                for relative_start in range(
                    0,
                    last_relative_start + 1,
                    self.query_stride,
                ):
                    byte_start = relative_start * 4
                    packed_key = bytes(
                        packed_view[
                            byte_start : byte_start + kmer_size * 4
                        ]
                    )
                    hits = get_packed(packed_key)
                    if hits:
                        candidates_update(hits)
            finally:
                packed_view.release()
        else:
            head_get = self.head_query.get
            for start in range(
                start_scan,
                length - kmer_size + 1,
                self.query_stride,
            ):
                key = tuple(barcodes[start : start + kmer_size])
                hits = head_get(key)
                if hits:
                    candidates_update(hits)

        # Snapshot the set's iteration order. Sorting this array would change
        # which paths a capped legacy search sees.
        return array.array("Q", candidates)

    @staticmethod
    def _candidate_cache_size(read_int, candidates):
        return (
            sys.getsizeof(read_int)
            + sys.getsizeof(candidates)
            + _CACHE_ENTRY_OVERHEAD
        )

    def _candidate_ids_for_read(self, read_int):
        cached = self._candidate_cache.get(read_int)
        if cached is not None:
            return cached

        barcode_result = self.barcode_query.get(read_int)
        # BarcodeStore returns (read_length, barcodes), while the lightweight
        # query doubles used by callers historically return barcodes directly.
        if isinstance(barcode_result, tuple) and len(barcode_result) == 2:
            _, barcodes = barcode_result
        else:
            barcodes = barcode_result
        if not barcodes:
            candidates = array.array("Q")
        else:
            candidates = self._get_candidate_ids_from_suffix(barcodes)

        self._candidate_cache.put(
            read_int,
            candidates,
            self._candidate_cache_size(read_int, candidates),
        )
        return candidates

    def candidate_cache_stats(self):
        return self._candidate_cache.stats()

    def head_cache_stats(self):
        cache_stats = getattr(self.head_query, "cache_stats", None)
        if cache_stats is not None:
            return cache_stats()
        return {
            "hits": 0,
            "misses": 0,
            "evictions": 0,
            "puts": 0,
            "rejections": 0,
            "entries": 0,
            "current_bytes": 0,
            "max_bytes": 0,
        }

    @staticmethod
    def _ancestry_contains(node_ids, parent_indexes, node_index, read_int):
        while True:
            if node_ids[node_index] == read_int:
                return True
            if node_index == 0:
                return False
            node_index = parent_indexes[node_index]

    @staticmethod
    def _reconstruct_child_path(
        node_ids,
        parent_indexes,
        node_index,
        child_id,
    ):
        reversed_path = [child_id]
        while True:
            reversed_path.append(node_ids[node_index])
            if node_index == 0:
                break
            node_index = parent_indexes[node_index]
        reversed_path.reverse()
        return tuple(reversed_path)

    @staticmethod
    def _retain_top_path(retained, path, limit):
        if limit <= 0:
            return
        retained.append(path)
        retained.sort(key=lambda item: (-len(item), item))
        if len(retained) > limit:
            retained.pop()

    def _search_source(
        self,
        source_id,
        target_ids,
        max_depth,
        top_k_per_anchor,
        max_expansions_per_source,
        source_label,
        ordinal,
        total_anchors,
    ):
        if max_depth > 255:
            raise ValueError("max_depth must be at most 255")

        node_ids = array.array("Q", [source_id])
        parent_indexes = array.array("I", [0])
        depths = array.array("B", [1])
        cursor = 0
        expanded = 0
        reached_cap = False
        retained = []
        frontier_high_water = 1
        started = time.monotonic()

        log(
            f"ANCHOR START {ordinal}/{total_anchors} {source_label} "
            f"{memory_summary()}"
        )

        while cursor < len(node_ids) and not reached_cap:
            node_index = cursor
            cursor += 1
            depth = depths[node_index]
            if depth > max_depth:
                continue

            right_id = node_ids[node_index]
            candidates = self._candidate_ids_for_read(right_id)
            for candidate in candidates:
                next_id = int(candidate)
                if self._ancestry_contains(
                    node_ids,
                    parent_indexes,
                    node_index,
                    next_id,
                ):
                    continue

                next_depth = depth + 1
                if next_id in target_ids and next_id != source_id:
                    path = self._reconstruct_child_path(
                        node_ids,
                        parent_indexes,
                        node_index,
                        next_id,
                    )
                    self._retain_top_path(
                        retained,
                        path,
                        top_k_per_anchor,
                    )

                if next_depth < max_depth:
                    if len(node_ids) > _MAX_PARENT_INDEX:
                        raise OverflowError(
                            "compact BFS parent index exceeded uint32 capacity"
                        )
                    node_ids.append(next_id)
                    parent_indexes.append(node_index)
                    depths.append(next_depth)
                    frontier = len(node_ids) - cursor
                    if frontier > frontier_high_water:
                        frontier_high_water = frontier

                expanded += 1
                if (
                    self.progress_every
                    and expanded % self.progress_every == 0
                ):
                    log(
                        f"ANCHOR PROGRESS {ordinal}/{total_anchors} "
                        f"{source_label} expansions={expanded:,} "
                        f"frontier={len(node_ids) - cursor:,} "
                        f"frontier_high_water={frontier_high_water:,} "
                        f"retained={len(retained):,} {memory_summary()} "
                        f"{_cache_summary('candidate_cache', self.candidate_cache_stats())} "
                        f"{_cache_summary('head_cache', self.head_cache_stats())}"
                    )

                if (
                    max_expansions_per_source is not None
                    and expanded >= max_expansions_per_source
                ):
                    reached_cap = True
                    log(
                        f"ANCHOR CAP {ordinal}/{total_anchors} {source_label} "
                        f"expansions={expanded:,} "
                        f"frontier={len(node_ids) - cursor:,}"
                    )
                    break

        selected = [list(path) for path in retained]
        elapsed = time.monotonic() - started
        metrics = {
            "expanded": expanded,
            "frontier": len(node_ids) - cursor,
            "frontier_high_water": frontier_high_water,
            "reached_cap": reached_cap,
            "elapsed": elapsed,
        }
        return selected, metrics

    def find_top_paths_per_anchor(
        self,
        anchor_ids,
        max_depth=10,
        top_k_per_anchor=5,
        max_expansions_per_source=None,
        target_ids=None,
        on_anchor_complete=None,
        anchor_offset=0,
        total_anchors=None,
    ):
        anchor_ids = list(anchor_ids)
        if target_ids is None:
            target_ids = anchor_ids
        target_ids = set(target_ids)
        if total_anchors is None:
            total_anchors = anchor_offset + len(anchor_ids)

        results = {}
        for local_index, source_id in enumerate(anchor_ids):
            ordinal = anchor_offset + local_index + 1
            source_label = self.int_to_rid.get(source_id, source_id)
            selected, metrics = self._search_source(
                source_id=source_id,
                target_ids=target_ids,
                max_depth=max_depth,
                top_k_per_anchor=top_k_per_anchor,
                max_expansions_per_source=max_expansions_per_source,
                source_label=source_label,
                ordinal=ordinal,
                total_anchors=total_anchors,
            )
            results[source_id] = selected

            # A main-pipeline callback atomically checkpoints this anchor before
            # DONE is emitted. A killed in-progress anchor is therefore retried.
            if on_anchor_complete is not None:
                on_anchor_complete(source_id, selected, ordinal, total_anchors)

            if len(selected) == top_k_per_anchor:
                log(
                    f"ANCHOR SATURATED {source_label}: reached "
                    f"{top_k_per_anchor} paths"
                )
            log(
                f"ANCHOR DONE {ordinal}/{total_anchors} {source_label} "
                f"expansions={metrics['expanded']:,} "
                f"frontier_high_water={metrics['frontier_high_water']:,} "
                f"retained={len(selected):,} "
                f"capped={metrics['reached_cap']} "
                f"elapsed={metrics['elapsed']:.1f}s {memory_summary()} "
                f"{_cache_summary('candidate_cache', self.candidate_cache_stats())} "
                f"{_cache_summary('head_cache', self.head_cache_stats())}"
            )

        return results


class CheckpointError(RuntimeError):
    pass


def _anchor_fingerprint(anchor_ids):
    digest = hashlib.sha256()
    digest.update(struct.pack(">Q", len(anchor_ids)))
    for anchor_id in anchor_ids:
        digest.update(struct.pack(">Q", int(anchor_id)))
    return digest.hexdigest()


def _normalize_input_paths(input_paths):
    if hasattr(input_paths, "items"):
        items = input_paths.items()
    else:
        items = ((os.fspath(path), path) for path in input_paths)
    normalized = {}
    for role, path in items:
        normalized[str(role)] = os.path.abspath(
            os.path.expanduser(os.fspath(path))
        )
    return normalized


def _input_signatures(input_paths):
    signatures = {}
    for role, path in _normalize_input_paths(input_paths).items():
        stat_result = os.stat(path)
        mtime_ns = getattr(
            stat_result,
            "st_mtime_ns",
            int(stat_result.st_mtime * 1000000000),
        )
        signatures[role] = {
            "path": path,
            "size": stat_result.st_size,
            "mtime_ns": mtime_ns,
        }
    return signatures


def build_checkpoint_metadata(anchor_ids, input_paths, search_config):
    ordered_anchor_ids = [int(anchor_id) for anchor_id in anchor_ids]
    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "python_version": tuple(sys.version_info[:2]),
        "anchor_count": len(ordered_anchor_ids),
        "anchor_fingerprint": _anchor_fingerprint(ordered_anchor_ids),
        "anchor_ids": ordered_anchor_ids,
        "input_signatures": _input_signatures(input_paths),
        "search_config": dict(search_config),
    }


def _checkpoint_failure(message):
    return CheckpointError(f"{message}; use --restart to start a fresh run")


def _validate_checkpoint_state(checkpoint, expected_metadata):
    if not isinstance(checkpoint, dict):
        raise _checkpoint_failure("checkpoint is not a dictionary")

    metadata_fields = (
        "schema_version",
        "python_version",
        "anchor_count",
        "anchor_fingerprint",
        "anchor_ids",
        "input_signatures",
        "search_config",
    )
    for field in metadata_fields:
        if checkpoint.get(field) != expected_metadata.get(field):
            raise _checkpoint_failure(f"checkpoint {field} does not match")

    completed = checkpoint.get("completed_anchor_ids")
    paths_by_anchor = checkpoint.get("paths_by_anchor")
    if not isinstance(completed, list):
        raise _checkpoint_failure("checkpoint completed anchors are invalid")
    if not isinstance(paths_by_anchor, dict):
        raise _checkpoint_failure("checkpoint paths are invalid")

    anchors = expected_metadata["anchor_ids"]
    if completed != anchors[: len(completed)]:
        raise _checkpoint_failure(
            "checkpoint completed anchors are not an exact ordered prefix"
        )
    if len(completed) != len(set(completed)):
        raise _checkpoint_failure("checkpoint completed anchors are duplicated")
    if set(paths_by_anchor) != set(completed):
        raise _checkpoint_failure(
            "checkpoint result keys do not match completed anchors"
        )

    for source_id in completed:
        paths = paths_by_anchor[source_id]
        if not isinstance(paths, list):
            raise _checkpoint_failure("checkpoint source paths are invalid")
        for path in paths:
            if not isinstance(path, list) or not path:
                raise _checkpoint_failure("checkpoint contains an invalid path")
            if path[0] != source_id or not all(
                isinstance(read_id, int) for read_id in path
            ):
                raise _checkpoint_failure("checkpoint contains an invalid path")


def write_checkpoint(
    checkpoint_path,
    metadata,
    completed_anchor_ids,
    paths_by_anchor,
):
    checkpoint = dict(metadata)
    checkpoint["completed_anchor_ids"] = list(completed_anchor_ids)
    checkpoint["paths_by_anchor"] = dict(paths_by_anchor)
    _validate_checkpoint_state(checkpoint, metadata)
    atomic_write_pickle(checkpoint, checkpoint_path)


def load_checkpoint(checkpoint_path, expected_metadata):
    try:
        with open(checkpoint_path, "rb") as stream:
            checkpoint = pkl.load(stream)
    except Exception as error:
        raise _checkpoint_failure(
            f"checkpoint could not be read ({error})"
        ) from error
    _validate_checkpoint_state(checkpoint, expected_metadata)
    return checkpoint


def _assert_inputs_unchanged(expected_metadata, input_paths):
    current = _input_signatures(input_paths)
    expected = expected_metadata["input_signatures"]
    if current != expected:
        changed = sorted(
            set(current) | set(expected),
            key=str,
        )
        changed = [
            role
            for role in changed
            if current.get(role) != expected.get(role)
        ]
        raise RuntimeError(
            "assembler inputs changed during the run: "
            + ", ".join(changed)
        )


def _positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def _nonnegative_int(value):
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def create_argument_parser():
    parser = argparse.ArgumentParser(
        description="Find barcode paths between anchor reads."
    )
    parser.add_argument(
        "--run-dir",
        "--base-dir",
        dest="base_dir",
        default=DEFAULT_BASE_DIR,
        help="Project artifact root (default: edit DEFAULT_BASE_DIR in this file)",
    )
    parser.add_argument(
        "--max-expansions-per-source",
        type=_positive_int,
        default=10000000,
    )
    parser.add_argument(
        "--candidate-cache-mib",
        type=_nonnegative_int,
        default=DEFAULT_CANDIDATE_CACHE_MIB,
    )
    parser.add_argument(
        "--head-cache-mib",
        type=_nonnegative_int,
        default=DEFAULT_HEAD_CACHE_MIB,
    )
    parser.add_argument(
        "--progress-every",
        type=_nonnegative_int,
        default=DEFAULT_PROGRESS_EVERY,
        help="Expansion interval for diagnostics; 0 disables periodic logs.",
    )
    parser.add_argument(
        "--checkpoint-path",
        default=None,
        help="Checkpoint pickle path (default: bridges_v1.1/barcode_assembler.checkpoint.pkl)",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        help="Replace any existing checkpoint and start at the first anchor.",
    )
    return parser


class TerminationRequested(BaseException):
    def __init__(self, signum):
        super().__init__(signum)
        self.signum = signum


def _signal_name(signum):
    try:
        return signal.Signals(signum).name
    except (ValueError, AttributeError):
        return str(signum)


def _termination_handler(signum, _frame):
    raise TerminationRequested(signum)


def _install_signal_handlers():
    previous = {}
    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous[signum] = signal.getsignal(signum)
            signal.signal(signum, _termination_handler)
    except ValueError:
        # signal.signal is limited to the main thread.
        for signum, handler in previous.items():
            try:
                signal.signal(signum, handler)
            except ValueError:
                pass
        return {}
    return previous


def _restore_signal_handlers(previous):
    for signum, handler in previous.items():
        signal.signal(signum, handler)


def _checkpoint_path(args):
    if args.checkpoint_path:
        return os.path.abspath(
            os.path.expanduser(os.fspath(args.checkpoint_path))
        )
    return os.path.join(
        args.base_dir,
        "bridges_v1.1",
        "barcode_assembler.checkpoint.pkl",
    )


def _convert_paths_to_read_ids(paths_by_anchor, ordered_anchor_ids, int_to_rid):
    converted = {}
    for source_id in ordered_anchor_ids:
        if source_id not in paths_by_anchor:
            raise RuntimeError(
                f"assembler has no completed result for anchor integer {source_id}"
            )
        try:
            source_read_id = int_to_rid[source_id]
        except KeyError as error:
            raise RuntimeError(
                f"int_to_rid.pkl has no label for anchor integer {source_id}"
            ) from error
        if not isinstance(source_read_id, str):
            raise RuntimeError(
                f"read label for integer {source_id} is not a string"
            )

        converted_paths = []
        for path in paths_by_anchor[source_id]:
            converted_path = []
            for read_int in path:
                try:
                    read_id = int_to_rid[read_int]
                except KeyError as error:
                    raise RuntimeError(
                        f"int_to_rid.pkl has no label for read integer {read_int}"
                    ) from error
                if not isinstance(read_id, str):
                    raise RuntimeError(
                        f"read label for integer {read_int} is not a string"
                    )
                converted_path.append(read_id)
            converted_paths.append(converted_path)
        converted[source_read_id] = converted_paths
    return converted


def run_pipeline(args):
    args.base_dir = os.path.abspath(os.path.expanduser(os.fspath(args.base_dir)))
    checkpoint_path = _checkpoint_path(args)

    k = 21
    k_mer_size = 11
    end_fraction = 0.25
    end_cap = 4000
    query_stride = 1
    max_depth = 15
    top_k_per_anchor = 10

    rid_map_pkl = os.path.join(args.base_dir, "rid_maps", "rid_to_int.pkl")
    int_rid_map_pkl = os.path.join(
        args.base_dir,
        "rid_maps",
        "int_to_rid.pkl",
    )
    anchor_pkl = os.path.join(args.base_dir, "output", "anchors.pkl")
    bridges_pkl = os.path.join(
        args.base_dir,
        "bridges_v1.1",
        "bridges.pkl",
    )
    store_dir = os.path.join(args.base_dir, f"barcode_{k}mers")
    idx_path = os.path.join(store_dir, "read_unikmer_map.idx")
    data_path = os.path.join(store_dir, "read_unikmer_map.data")
    head_idx_path = os.path.join(store_dir, "head_index", "head_index.idx")
    head_data_path = os.path.join(
        store_dir,
        "head_index",
        "head_index.data",
    )

    input_paths = {
        "rid_to_int": rid_map_pkl,
        "int_to_rid": int_rid_map_pkl,
        "anchors": anchor_pkl,
        "barcode_idx": idx_path,
        "barcode_data": data_path,
        "head_idx": head_idx_path,
        "head_data": head_data_path,
    }
    checkpoint_target = os.path.realpath(checkpoint_path)
    protected_paths = dict(input_paths)
    protected_paths["final bridges output"] = bridges_pkl
    for role, protected_path in protected_paths.items():
        if checkpoint_target == os.path.realpath(protected_path):
            raise ValueError(
                f"checkpoint path must not overwrite {role}: {protected_path}"
            )
    search_config = {
        "k": k,
        "k_mer_size": k_mer_size,
        "end_fraction": end_fraction,
        "end_cap": end_cap,
        "query_stride": query_stride,
        "max_depth": max_depth,
        "top_k_per_anchor": top_k_per_anchor,
        "max_expansions_per_source": args.max_expansions_per_source,
    }

    log(
        "START barcode assembler "
        f"python={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro} "
        f"base_dir={args.base_dir} max_expansions={args.max_expansions_per_source:,} "
        f"candidate_cache={args.candidate_cache_mib}MiB "
        f"head_cache={args.head_cache_mib}MiB "
        f"progress_every={args.progress_every:,} checkpoint={checkpoint_path}"
    )
    if tuple(sys.version_info[:2]) != (3, 8):
        log(
            "WARNING documented runtime is Python 3.8; checkpoints are "
            "restricted to this runtime's major/minor version"
        )
    log(
        "SIGKILL cannot be logged; use scheduler accounting or "
        f"/usr/bin/time -v. {memory_summary()}"
    )

    rid_to_int = load_pickle(rid_map_pkl)
    anchors = load_pickle(anchor_pkl)
    ordered_anchor_ids = []
    anchor_labels = {}
    for read_id in anchors.keys():
        read_int = rid_to_int.get(read_id)
        if read_int is None:
            continue
        read_int = int(read_int)
        ordered_anchor_ids.append(read_int)
        anchor_labels[read_int] = read_id
    if not ordered_anchor_ids:
        raise ValueError("No valid anchor IDs found in anchors.pkl.")

    metadata = build_checkpoint_metadata(
        ordered_anchor_ids,
        input_paths,
        search_config,
    )

    # These two objects are large and are not needed during graph traversal.
    del anchors
    del rid_to_int
    gc.collect()
    log(
        f"Prepared {len(ordered_anchor_ids):,} ordered anchors and released "
        f"startup pickles. {memory_summary()}"
    )

    if args.restart:
        if os.path.exists(bridges_pkl):
            log(
                f"RESTART requested; existing {bridges_pkl} remains stale "
                "until this run completes"
            )
        write_checkpoint(checkpoint_path, metadata, [], {})
        completed_anchor_ids = []
        paths_by_anchor = {}
    elif os.path.exists(checkpoint_path):
        checkpoint = load_checkpoint(checkpoint_path, metadata)
        completed_anchor_ids = list(checkpoint["completed_anchor_ids"])
        paths_by_anchor = dict(checkpoint["paths_by_anchor"])
        log(
            f"Loaded checkpoint with {len(completed_anchor_ids):,}/"
            f"{len(ordered_anchor_ids):,} completed anchors"
        )
    else:
        completed_anchor_ids = []
        paths_by_anchor = {}
        write_checkpoint(checkpoint_path, metadata, [], {})
        log(f"Initialized checkpoint {checkpoint_path}")

    for ordinal, source_id in enumerate(completed_anchor_ids, start=1):
        log(
            f"ANCHOR RESUMED {ordinal}/{len(ordered_anchor_ids)} "
            f"{anchor_labels.get(source_id, source_id)}"
        )

    remaining = ordered_anchor_ids[len(completed_anchor_ids) :]
    if remaining:
        barcode_query = None
        head_query = None
        assembler = None
        try:
            barcode_query = BarcodeStore(idx_path, data_path)
            _advise_barcode_store_random(barcode_query)
            head_query = HeadIndexQuery(
                head_idx_path,
                head_data_path,
                cache_bytes=args.head_cache_mib * _MIB,
            )
            assembler = ProjectAssembler(
                barcode_query=barcode_query,
                head_query=head_query,
                int_to_rid=anchor_labels,
                k_mer_size=k_mer_size,
                end_fraction=end_fraction,
                end_cap=end_cap,
                query_stride=query_stride,
                candidate_cache_mib=args.candidate_cache_mib,
                progress_every=args.progress_every,
            )

            def checkpoint_anchor(source_id, selected, _ordinal, _total):
                paths_by_anchor[source_id] = selected
                completed_anchor_ids.append(source_id)
                _assert_inputs_unchanged(metadata, input_paths)
                write_checkpoint(
                    checkpoint_path,
                    metadata,
                    completed_anchor_ids,
                    paths_by_anchor,
                )

            assembler.find_top_paths_per_anchor(
                anchor_ids=remaining,
                max_depth=max_depth,
                top_k_per_anchor=top_k_per_anchor,
                max_expansions_per_source=args.max_expansions_per_source,
                target_ids=ordered_anchor_ids,
                on_anchor_complete=checkpoint_anchor,
                anchor_offset=len(completed_anchor_ids),
                total_anchors=len(ordered_anchor_ids),
            )
            log(
                "SEARCH COMPLETE "
                f"{_cache_summary('candidate_cache', assembler.candidate_cache_stats())} "
                f"{_cache_summary('head_cache', head_query.cache_stats())} "
                f"{memory_summary()}"
            )
        except MemoryError:
            details = ["MEMORY ERROR during graph search", memory_summary()]
            if assembler is not None:
                details.append(
                    _cache_summary(
                        "candidate_cache",
                        assembler.candidate_cache_stats(),
                    )
                )
            if head_query is not None:
                details.append(
                    _cache_summary("head_cache", head_query.cache_stats())
                )
            log(" ".join(details))
            raise
        finally:
            try:
                if head_query is not None:
                    head_query.close()
            finally:
                try:
                    if barcode_query is not None:
                        barcode_query.close()
                finally:
                    assembler = None
                    gc.collect()

    _assert_inputs_unchanged(metadata, input_paths)
    if completed_anchor_ids != ordered_anchor_ids:
        raise RuntimeError("search ended without completing every anchor")

    # Load this large map only after search arrays, caches, and mmaps are gone.
    int_to_rid = load_pickle(int_rid_map_pkl)
    read_paths_by_anchor = _convert_paths_to_read_ids(
        paths_by_anchor,
        ordered_anchor_ids,
        int_to_rid,
    )
    del int_to_rid
    _assert_inputs_unchanged(metadata, input_paths)

    total_paths = sum(len(paths) for paths in read_paths_by_anchor.values())
    write_pickle(read_paths_by_anchor, bridges_pkl)
    log(
        f"Saved {total_paths:,} total paths for "
        f"{len(read_paths_by_anchor):,} anchors to {bridges_pkl}"
    )
    if total_paths == 0:
        log("No path found.")
    return 0


def main(argv=None):
    parser = create_argument_parser()
    args = parser.parse_args(argv)
    previous_handlers = _install_signal_handlers()
    try:
        return run_pipeline(args)
    except TerminationRequested as error:
        log(
            f"TERMINATED by {_signal_name(error.signum)}; the last completed "
            f"checkpoint, if initialized, is at {_checkpoint_path(args)}"
        )
        return 128 + error.signum
    except KeyboardInterrupt:
        log(
            f"INTERRUPTED; the last completed anchor remains in "
            f"{_checkpoint_path(args)}"
        )
        return 130
    except MemoryError:
        log(
            "ABORTED after MemoryError; the last completed checkpoint, if "
            f"initialized, is at {_checkpoint_path(args)}. "
            f"{memory_summary()}"
        )
        return 1
    except CheckpointError as error:
        log(f"CHECKPOINT ERROR: {error}")
        return 2
    finally:
        _restore_signal_handlers(previous_handlers)


if __name__ == "__main__":
    raise SystemExit(main())
