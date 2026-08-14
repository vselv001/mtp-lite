import array
import collections
import functools
import mmap
import os
import pickle as pkl
import struct
import sys

from barcode_store import BarcodeStore

_NEED_SWAP = sys.byteorder == "little"

_EMPTY_ARRAY = array.array("Q")
BARCODE_CACHE_SIZE = 20000
HEAD_INDEX_CACHE_SIZE = 50000

# head_index.idx record: 11 × uint32 kmer | uint64 data-offset | uint32 count
_KMER_LEN = 11
_KMER_BYTES = _KMER_LEN * 4                        # 44
_IDX_RECORD = struct.Struct(f">{_KMER_LEN}IQI")    # 56 bytes
_IDX_ENTRY_SIZE = _IDX_RECORD.size                  # 56
_KMER_FMT = struct.Struct(f">{_KMER_LEN}I")
_TAIL_FMT = struct.Struct(">QI")                    # offset + count


def load_pickle(in_file):
    print(f"Loading {in_file}...")
    with open(in_file, "rb") as f:
        return pkl.load(f)


def write_pickle(data, out_file):
    print(f"Writing results to {out_file}...")
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, 'wb') as f:
        pkl.dump(data, f)


class CachedBarcodeStore:
    """BarcodeStore wrapper with LRU cache, compatible with assembler's .get(int) interface."""
    __slots__ = ("_store",)

    def __init__(self, idx_path, data_path):
        self._store = BarcodeStore(idx_path, data_path)

    @functools.lru_cache(maxsize=BARCODE_CACHE_SIZE)
    def get(self, read_int):
        _, barcodes = self._store.get(read_int)
        return barcodes

    def close(self):
        self.get.cache_clear()
        self._store.close()


class HeadIndexQuery:
    """Binary-search lookup on the sorted head_index .idx/.data flat files."""
    __slots__ = ("_idx_mm", "_data_mm", "_idx_fd", "_data_fd", "_n_records")

    def __init__(self, idx_path, data_path):
        self._idx_fd = os.open(idx_path, os.O_RDONLY)
        idx_size = os.fstat(self._idx_fd).st_size
        self._idx_mm = mmap.mmap(self._idx_fd, idx_size, access=mmap.ACCESS_READ)
        self._n_records = idx_size // _IDX_ENTRY_SIZE

        self._data_fd = os.open(data_path, os.O_RDONLY)
        data_size = os.fstat(self._data_fd).st_size
        self._data_mm = mmap.mmap(self._data_fd, data_size, access=mmap.ACCESS_READ) if data_size else None

    @functools.lru_cache(maxsize=HEAD_INDEX_CACHE_SIZE)
    def get(self, kmer_tuple):
        idx_mm = self._idx_mm
        lo, hi = 0, self._n_records
        while lo < hi:
            mid = (lo + hi) >> 1
            off = mid * _IDX_ENTRY_SIZE
            mid_kmer = _KMER_FMT.unpack_from(idx_mm, off)
            if mid_kmer < kmer_tuple:
                lo = mid + 1
            elif mid_kmer > kmer_tuple:
                hi = mid
            else:
                data_offset, count = _TAIL_FMT.unpack_from(idx_mm, off + _KMER_BYTES)
                if count == 0:
                    return _EMPTY_ARRAY
                raw = self._data_mm[data_offset:data_offset + count * 8]
                hits = array.array("Q")
                hits.frombytes(raw)
                if _NEED_SWAP:
                    hits.byteswap()
                return hits
        return _EMPTY_ARRAY

    def close(self):
        self.get.cache_clear()
        if self._idx_mm:
            self._idx_mm.close()
            self._idx_mm = None
        if self._data_mm:
            self._data_mm.close()
            self._data_mm = None
        for fd_attr in ("_idx_fd", "_data_fd"):
            fd = getattr(self, fd_attr, None)
            if fd is not None:
                os.close(fd)
                setattr(self, fd_attr, None)


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
    ):
        self.barcode_query = barcode_query
        self.head_query = head_query
        self.int_to_rid = int_to_rid or {}
        self.k_mer_size = k_mer_size
        self.end_fraction = end_fraction
        self.end_cap = end_cap
        self.query_stride = max(1, int(query_stride))

    def _end_window(self, n):
        return min(int(n * self.end_fraction), self.end_cap)

    def _get_candidate_ids_from_suffix(self, barcodes):
        k = self.k_mer_size
        n = len(barcodes)
        window = self._end_window(n)
        if window < k:
            return set()
        start_scan = n - window
        stride = self.query_stride
        head_get = self.head_query.get
        candidates = set()
        candidates_update = candidates.update
        for i in range(start_scan, n - k + 1, stride):
            kmer = tuple(barcodes[i : i + k])
            hits = head_get(kmer)
            if hits:
                candidates_update(hits)
        return candidates

    def find_top_paths_per_anchor(
        self,
        anchor_ids,
        max_depth=10,
        top_k_per_anchor=5,
        max_expansions_per_source=None,
    ):
        anchor_ids = list(anchor_ids)
        target_ids = set(anchor_ids)
        results = {}
        barcode_get = self.barcode_query.get
        get_candidates = self._get_candidate_ids_from_suffix

        for source_id in anchor_ids:
            source_label = self.int_to_rid.get(source_id, source_id)

            queue = collections.deque([(source_id, [source_id], frozenset([source_id]))])
            seen_paths = {(source_id,)}
            source_paths = []
            expanded = 0
            reached_expansion_cap = False

            while queue and not reached_expansion_cap:
                _, path, path_set = queue.popleft()
                if len(path) > max_depth:
                    continue

                right_id = path[-1]
                right_barcodes = barcode_get(right_id)
                if not right_barcodes:
                    continue

                candidates = get_candidates(right_barcodes)
                for next_id in candidates:
                    if next_id in path_set:
                        continue

                    next_path = path + [next_id]
                    next_path_key = tuple(next_path)
                    if next_path_key in seen_paths:
                        continue
                    seen_paths.add(next_path_key)

                    if next_id in target_ids and next_id != source_id:
                        source_paths.append(next_path)

                    if len(next_path) < max_depth:
                        queue.append((next_id, next_path, path_set | {next_id}))

                    expanded += 1
                    if (
                        max_expansions_per_source is not None
                        and expanded >= max_expansions_per_source
                    ):
                        reached_expansion_cap = True
                        break

            source_paths.sort(key=lambda p: (-len(p), tuple(p)))
            selected = source_paths[:top_k_per_anchor]
            results[source_id] = selected
            if len(selected) == top_k_per_anchor:
                print(
                    f"ANCHOR SATURATED {source_label}: reached {top_k_per_anchor} paths"
                )

        return results


def main():
    base_dir = "path/to/your/project_directory"  # Replace with your MTPLite project directory
    k = 21  # Replace with the k-mer size used for unikmer extraction
    k_mer_size = 11
    end_fraction = 0.25
    end_cap = 4000
    query_stride = 1
    max_depth = 15
    top_k_per_anchor = 10
    max_expansions_per_source = 10000000

    rid_map_pkl = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    anchor_pkl = os.path.join(base_dir, "output", "anchors.pkl")
    bridges_pkl = os.path.join(base_dir, "bridges_v1.1", "bridges.pkl")
    store_dir = os.path.join(base_dir, f"barcode_{k}mers")
    idx_path = os.path.join(store_dir, "read_unikmer_map.idx")
    data_path = os.path.join(store_dir, "read_unikmer_map.data")
    head_idx_path = os.path.join(base_dir, f"barcode_{k}mers", "head_index", "head_index.idx")
    head_data_path = os.path.join(base_dir, f"barcode_{k}mers", "head_index", "head_index.data")

    rid_to_int = load_pickle(rid_map_pkl)
    anchors = load_pickle(anchor_pkl)
    int_to_rid = {int(v): k for k, v in rid_to_int.items()}
    start_ids = [int(rid_to_int[rid]) for rid in anchors.keys() if rid in rid_to_int]
    if not start_ids:
        raise ValueError("No valid anchor IDs found in anchors.pkl.")

    barcode_query = CachedBarcodeStore(idx_path, data_path)
    head_query = HeadIndexQuery(head_idx_path, head_data_path)

    try:
        assembler = ProjectAssembler(
            barcode_query=barcode_query,
            head_query=head_query,
            int_to_rid=int_to_rid,
            k_mer_size=k_mer_size,
            end_fraction=end_fraction,
            end_cap=end_cap,
            query_stride=query_stride,
        )
        paths_by_anchor = assembler.find_top_paths_per_anchor(
            anchor_ids=start_ids,
            max_depth=max_depth,
            top_k_per_anchor=top_k_per_anchor,
            max_expansions_per_source=max_expansions_per_source,
        )
    finally:
        barcode_query.close()
        head_query.close()

    read_paths_by_anchor = {}
    total_paths = 0
    for source_id, source_paths in paths_by_anchor.items():
        source_read_id = int_to_rid.get(source_id, source_id)
        converted = [[int_to_rid.get(rid, rid) for rid in path] for path in source_paths]
        read_paths_by_anchor[source_read_id] = converted
        total_paths += len(converted)

    write_pickle(read_paths_by_anchor, bridges_pkl)
    print(
        f"Saved {total_paths:,} total paths "
        f"for {len(read_paths_by_anchor):,} anchors to {bridges_pkl}"
    )

    if total_paths == 0:
        print("No path found.")


if __name__ == "__main__":
    main()
