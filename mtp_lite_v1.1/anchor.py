import pickle as pkl
import os
import sys
import heapq
import array
import mmap
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from barcode_store import IDX_ENTRY

_EMPTY = np.empty(0, dtype=np.uint32)


class FlatBarcodeQuery:
    """Numpy-cached mmap reader for the .idx/.data flat files."""

    __slots__ = ("_idx_mm", "_data_mm", "_idx_fd", "_data_fd",
                 "_entry_size", "_unpack", "_cache")

    def __init__(self, idx_path, data_path):
        self._idx_fd = os.open(idx_path, os.O_RDONLY)
        idx_size = os.fstat(self._idx_fd).st_size
        self._idx_mm = mmap.mmap(self._idx_fd, idx_size, access=mmap.ACCESS_READ)

        self._data_fd = os.open(data_path, os.O_RDONLY)
        data_size = os.fstat(self._data_fd).st_size
        self._data_mm = mmap.mmap(self._data_fd, data_size, access=mmap.ACCESS_READ) if data_size else None

        self._entry_size = IDX_ENTRY.size
        self._unpack = IDX_ENTRY.unpack_from
        self._cache = {}

    def get(self, read_int):
        cache = self._cache
        result = cache.get(read_int)
        if result is not None:
            return result

        off = read_int * self._entry_size
        data_offset, num_barcodes, read_len = self._unpack(self._idx_mm, off)
        if num_barcodes == 0:
            cache[read_int] = _EMPTY
            return _EMPTY

        raw = self._data_mm[data_offset:data_offset + num_barcodes * 4]
        barcodes = np.frombuffer(raw, dtype=np.uint32).copy()
        cache[read_int] = barcodes
        return barcodes

    def get_stats(self, read_int):
        """Return (num_barcodes, read_len) from the .idx entry."""
        off = read_int * self._entry_size
        _, num_barcodes, read_len = self._unpack(self._idx_mm, off)
        return num_barcodes, read_len

    def close(self):
        self._cache.clear()
        self._idx_mm.close()
        os.close(self._idx_fd)
        if self._data_mm is not None:
            self._data_mm.close()
        os.close(self._data_fd)


def loadPickle(inFile):
    print(f"Loading {inFile}...")
    with open(inFile, 'rb') as f:
        return pkl.load(f)


def writePickle(data, outFile):
    print(f"Writing results to {outFile}...")
    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    with open(outFile, 'wb') as f:
        pkl.dump(data, f)


def select_anchors(universe, rid_to_int, ctx, update_every=100000):
    get_stats = ctx.get_stats
    heap = []
    for rid, read_int in rid_to_int.items():
        num_barcodes, _ = get_stats(read_int)
        if num_barcodes > 0:
            heap.append((-num_barcodes, read_int, rid))
    heapq.heapify(heap)

    print(f"Heap initialized with {len(heap)} valid reads.")

    selected = {}
    processed = 0
    get_barcodes = ctx.get
    u = universe
    _count_nonzero = np.count_nonzero

    while heap:
        processed += 1
        if processed % update_every == 0:
            print(
                f"[{processed}] iterations | "
                f"heap size: {len(heap)} | "
                f"selected: {len(selected)}"
            )

        _, read_int, best_id = heapq.heappop(heap)

        barcodes = get_barcodes(read_int)
        if len(barcodes) == 0:
            continue

        mask = u[barcodes] != 0
        actual_score = int(_count_nonzero(mask))
        if actual_score == 0:
            continue

        next_best = -heap[0][0] if heap else -1
        if actual_score >= next_best:
            kept = barcodes[mask]
            kept_arr = array.array('I')
            kept_arr.frombytes(kept.tobytes())
            _, read_len = get_stats(read_int)
            selected[best_id] = {
                "barcodes": kept_arr,
                "length": read_len,
            }
            u[kept] = 0
        else:
            heapq.heappush(heap, (-actual_score, read_int, best_id))

    print(f"Finished after {processed} iterations.")
    return selected


def anchor_executor():
    base_dir = "path/to/your/project_directory"  # Replace with your MTPLite project directory
    k = 21  # Replace with the k-mer size used for unikmer extraction

    universe_pkl = os.path.join(base_dir, "universe", "universe.pkl")
    rid_int_map_pkl = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    out_pkl = os.path.join(base_dir, "output", "anchors.pkl")
    idx_path = os.path.join(base_dir, f"barcode_{k}mers", "read_unikmer_map.idx")
    data_path = os.path.join(base_dir, f"barcode_{k}mers", "read_unikmer_map.data")

    universe_list = loadPickle(universe_pkl)

    universe_arr = np.fromiter(universe_list, dtype=np.uint32, count=len(universe_list))
    max_barcode = int(universe_arr.max())
    universe = np.zeros(max_barcode + 1, dtype=np.uint8)
    universe[universe_arr] = 1
    universe_size = len(universe_arr)

    rid_to_int = loadPickle(rid_int_map_pkl)

    print(f"Finished loading. Universe size: {universe_size}, Reads: {len(rid_to_int)}")

    ctx = FlatBarcodeQuery(idx_path, data_path)
    try:
        selected_anchors = select_anchors(universe, rid_to_int, ctx)

        remaining = int(np.count_nonzero(universe))
        covered = universe_size - remaining
        coverage_pct = 100.0 * (covered / universe_size) if universe_size else 0.0

        print(
            f"Unikmers covered: {covered}/{universe_size} "
            f"({coverage_pct:.2f}%)"
        )

        writePickle(selected_anchors, out_pkl)
        print(f"Saved {len(selected_anchors)} anchors to {out_pkl}")
    finally:
        ctx.close()

if __name__ == "__main__":
    anchor_executor()
