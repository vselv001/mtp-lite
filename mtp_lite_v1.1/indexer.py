import array
import os
import pickle
import shutil
import struct
import mmap

from barcode_store import BarcodeStore

DEFAULT_PROGRESS_EVERY = 100000

# ── flat-file head index reader ────────────────────────────────────


class HeadIndexStore:
    """Read-only binary-search accessor for head_index .idx / .data files.

    .idx layout (sorted by kmer key):
        [kmer_bytes  (k*4 B)] [data_offset (8 B)] [num_read_ids (4 B)]

    .data layout:
        Packed uint64 read-IDs, referenced by offset from .idx entries.
    """

    __slots__ = (
        "k_mer_size", "_key_size", "_entry_size", "_tail_st",
        "_idx_mm", "_data_mm", "_idx_fd", "_data_fd", "_num_entries",
    )

    def __init__(self, idx_path, data_path, k_mer_size):
        self.k_mer_size = k_mer_size
        self._key_size = k_mer_size * 4
        self._entry_size = self._key_size + 12  # +8 offset +4 count
        self._tail_st = struct.Struct(">QI")

        self._idx_fd = os.open(idx_path, os.O_RDONLY)
        idx_size = os.fstat(self._idx_fd).st_size
        self._idx_mm = mmap.mmap(self._idx_fd, idx_size, access=mmap.ACCESS_READ)
        self._num_entries = idx_size // self._entry_size

        self._data_fd = os.open(data_path, os.O_RDONLY)
        data_size = os.fstat(self._data_fd).st_size
        self._data_mm = (
            mmap.mmap(self._data_fd, data_size, access=mmap.ACCESS_READ)
            if data_size else None
        )

    def query(self, kmer_key_bytes):
        """Binary search for a kmer key (raw bytes). Return list of read IDs."""
        lo, hi = 0, self._num_entries
        es = self._entry_size
        ks = self._key_size
        mm = self._idx_mm

        while lo < hi:
            mid = (lo + hi) >> 1
            pos = mid * es
            mid_key = mm[pos : pos + ks]
            if mid_key < kmer_key_bytes:
                lo = mid + 1
            elif mid_key > kmer_key_bytes:
                hi = mid
            else:
                offset, count = self._tail_st.unpack_from(mm, pos + ks)
                if count == 0:
                    return []
                raw = self._data_mm[offset : offset + count * 8]
                arr = array.array("Q")
                arr.frombytes(raw)
                return list(arr)
        return []

    def query_tuple(self, kmer_tuple):
        """Convenience: accept a tuple of barcode ints."""
        key = struct.pack(f">{self.k_mer_size}I", *kmer_tuple)
        return self.query(key)

    def close(self):
        self._idx_mm.close()
        os.close(self._idx_fd)
        if self._data_mm is not None:
            self._data_mm.close()
        os.close(self._data_fd)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# ── indexer (write path) ───────────────────────────────────────────


class BarcodeIndexer:
    def __init__(
        self,
        db_path,
        mode,
        k_mer_size,
        index_stride=3,
        end_fraction=0.25,
        end_cap=4000,
        **_ignored,
    ):
        self.db_path = db_path
        self.k_mer_size = k_mer_size
        self.index_stride = index_stride
        self.end_fraction = end_fraction
        self.end_cap = end_cap
        self._kmer_key_struct = struct.Struct(f">{k_mer_size}I")
        self._pack_q = struct.Struct(">Q").pack

        self.idx_path = os.path.join(db_path, "head_index.idx")
        self.data_path = os.path.join(db_path, "head_index.data")

        if mode == "w":
            if os.path.exists(db_path):
                shutil.rmtree(db_path)
            os.makedirs(db_path, exist_ok=True)
            self._kmer_acc = {}  # kmer_bytes -> bytearray of packed uint64 read-IDs
        else:
            self._store = HeadIndexStore(self.idx_path, self.data_path, k_mer_size)

    def close(self):
        store = getattr(self, "_store", None)
        if store is not None:
            store.close()

    def _kmer_key(self, barcode_arr, start):
        k = self.k_mer_size
        return self._kmer_key_struct.pack(*barcode_arr[start : start + k])

    def _get_starts(self, n, use_full_depth):
        k = self.k_mer_size
        if n < k:
            return None
        if use_full_depth:
            return range(0, n - k + 1)
        stride = self.index_stride
        window = min(int(n * self.end_fraction), self.end_cap)
        if window < k:
            return None
        head_starts = range(0, window - k + 1, stride)
        tail_start = n - window
        tail_starts = range(tail_start, n - k + 1, stride)
        seen = set()
        starts = []
        for i in head_starts:
            if i not in seen:
                seen.add(i)
                starts.append(i)
        for i in tail_starts:
            if i not in seen:
                seen.add(i)
                starts.append(i)
        return starts

    # ── write: accumulate in memory ────────────────────────────────

    def add_reads_from_barcode_store(
        self,
        store,
        rid_to_int,
        read_ids,
        anchor_ids=None,
        progress_every=DEFAULT_PROGRESS_EVERY,
        **_ignored,
    ):
        total = len(read_ids)
        anchor_ids = anchor_ids or set()
        print(
            f"Indexing {total:,} reads from BarcodeStore "
            f"(anchor reads at full depth: {len(anchor_ids):,})..."
        )

        processed = 0
        missing = 0
        full_depth_reads = 0
        store_get = store.get
        kmer_key = self._kmer_key
        pack_q = self._pack_q
        acc = self._kmer_acc

        for read_id in read_ids:
            read_int = rid_to_int.get(read_id)
            if read_int is None:
                missing += 1
                continue

            _, barcode_arr = store_get(read_int)
            if not barcode_arr:
                missing += 1
                continue

            use_full_depth = read_id in anchor_ids
            starts = self._get_starts(len(barcode_arr), use_full_depth)
            if starts is None:
                missing += 1
                continue

            rid_packed = pack_q(int(read_int))
            for i in starts:
                key = kmer_key(barcode_arr, i)
                buf = acc.get(key)
                if buf is None:
                    acc[key] = bytearray(rid_packed)
                else:
                    buf += rid_packed

            if use_full_depth:
                full_depth_reads += 1

            processed += 1
            if processed % progress_every == 0:
                print(f"  indexed {processed:,}/{total:,} reads  (unique kmers: {len(acc):,})")

        print(
            f"Done accumulating. Indexed: {processed:,}, missing/skipped: {missing:,}, "
            f"full-depth anchor reads: {full_depth_reads:,}, "
            f"unique k-mers: {len(acc):,}"
        )

        self._write_flat_files()

    def add_reads_bulk(self, reads_dict):
        """Index reads from an in-memory dict (legacy interface)."""
        print(f"Indexing {len(reads_dict)} reads...")
        kmer_key = self._kmer_key
        pack_q = self._pack_q
        acc = self._kmer_acc

        for idx, (r_id, barcode_arr) in enumerate(reads_dict.items(), start=1):
            starts = self._get_starts(len(barcode_arr), use_full_depth=False)
            if starts is None:
                continue
            rid_packed = pack_q(r_id)
            for i in starts:
                key = kmer_key(barcode_arr, i)
                buf = acc.get(key)
                if buf is None:
                    acc[key] = bytearray(rid_packed)
                else:
                    buf += rid_packed
            if idx % DEFAULT_PROGRESS_EVERY == 0:
                print(f"  indexed {idx:,} reads")

        self._write_flat_files()

    def _write_flat_files(self):
        """Sort accumulated k-mers and write .idx + .data flat files."""
        acc = self._kmer_acc
        print(f"Sorting {len(acc):,} unique k-mers...")
        sorted_keys = sorted(acc)

        tail_st = struct.Struct(">QI")  # data_offset (8B) + count (4B)
        key_size = self.k_mer_size * 4
        entry_size = key_size + tail_st.size

        print(f"Writing {self.data_path} and {self.idx_path} ...")
        data_offset = 0
        with open(self.data_path, "wb", buffering=16 * 1024 * 1024) as data_f, \
             open(self.idx_path, "wb", buffering=16 * 1024 * 1024) as idx_f:
            for key in sorted_keys:
                rid_bytes = acc[key]
                count = len(rid_bytes) // 8
                # write index entry: kmer_key + offset + count
                idx_f.write(key)
                idx_f.write(tail_st.pack(data_offset, count))
                # write data
                data_f.write(rid_bytes)
                data_offset += len(rid_bytes)

        self._kmer_acc = None  # free memory
        print(f"Wrote {len(sorted_keys):,} index entries, {data_offset:,} bytes of read-ID data.")

    # ── read: query via HeadIndexStore ─────────────────────────────

    def query_head_index(self, kmer_tuple):
        return self._store.query_tuple(kmer_tuple)


def load_pickle(path):
    print(f"Loading {path}...")
    with open(path, "rb") as f:
        return pickle.load(f)


def main():
    base_dir = "path/to/your/project_directory"  # Replace with your MTPLite project directory
    k = 21  # Replace with the k-mer size used for unikmer extraction

    k_mer_size = 11
    progress_every = DEFAULT_PROGRESS_EVERY

    store_dir = os.path.join(base_dir, f"barcode_{k}mers")
    idx_path = os.path.join(store_dir, "read_unikmer_map.idx")
    data_path = os.path.join(store_dir, "read_unikmer_map.data")
    rid_map_pkl = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    read_stats_pkl = os.path.join(base_dir, "stats", "read_stats.pkl")
    anchors_pkl = os.path.join(base_dir, "output", "anchors.pkl")
    out_dir = os.path.join(base_dir, f"barcode_{k}mers", "head_index")

    rid_to_int = load_pickle(rid_map_pkl)
    read_stats = load_pickle(read_stats_pkl)
    read_ids = read_stats.keys()

    anchor_ids = set()
    if os.path.exists(anchors_pkl):
        anchors = load_pickle(anchors_pkl)
        anchor_ids = set(anchors.keys())
        print(f"Loaded {len(anchor_ids):,} anchor IDs for full-depth indexing.")
    else:
        print(f"Anchor file not found at {anchors_pkl}. Using shallow indexing for all reads.")

    indexer = BarcodeIndexer(
        db_path=out_dir,
        mode="w",
        k_mer_size=k_mer_size,
    )
    try:
        print(f"Opening barcode store: {store_dir}")
        with BarcodeStore(idx_path, data_path) as store:
            indexer.add_reads_from_barcode_store(
                store=store,
                rid_to_int=rid_to_int,
                read_ids=read_ids,
                anchor_ids=anchor_ids,
                progress_every=progress_every,
            )
    finally:
        indexer.close()


if __name__ == "__main__":
    main()
