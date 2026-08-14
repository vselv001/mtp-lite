import os
import pickle as pkl

from barcode_store import BarcodeStore


class RollingHash:
    __slots__ = ('MASK_64', 'MOD_2', 'BASE', 'min_overlap',
                 'id_map', 'main_hashes_1', 'main_hashes_2', 'main_strings',
                 'powers_1', 'powers_2')

    def __init__(self, anchors, min_overlap):
        self.MASK_64 = (1 << 64) - 1
        self.MOD_2 = 10**9 + 7
        self.BASE = 4_294_967_291
        self.min_overlap = min_overlap

        self.id_map = {}
        self.main_hashes_1 = {}
        self.main_hashes_2 = {}
        self.main_strings = {}

        self.powers_1 = [1]
        self.powers_2 = [1]

        self._preprocess(anchors)

    def _preprocess(self, anchors):
        id_map = self.id_map
        MASK_64 = self.MASK_64
        MOD_2 = self.MOD_2
        BASE = self.BASE

        for main_id, info in anchors.items():
            barcodes = info["barcodes"] if isinstance(info, dict) else info
            self.main_strings[main_id] = barcodes

            h1 = h2 = 0
            row_h1 = []
            row_h2 = []
            append_h1 = row_h1.append
            append_h2 = row_h2.append

            for char_idx, char_id in enumerate(barcodes):
                id_map.setdefault(char_id, []).append((main_id, char_idx))
                h1 = (h1 * BASE + char_id) & MASK_64
                append_h1(h1)
                h2 = (h2 * BASE + char_id) % MOD_2
                append_h2(h2)

            self.main_hashes_1[main_id] = row_h1
            self.main_hashes_2[main_id] = row_h2

    def _ensure_powers(self, n):
        current_len = len(self.powers_1)
        if n < current_len:
            return
        p1, p2 = self.powers_1[-1], self.powers_2[-1]
        BASE = self.BASE
        MASK_64 = self.MASK_64
        MOD_2 = self.MOD_2
        pw1_append = self.powers_1.append
        pw2_append = self.powers_2.append
        for _ in range(n - current_len + 1):
            p1 = (p1 * BASE) & MASK_64
            p2 = (p2 * BASE) % MOD_2
            pw1_append(p1)
            pw2_append(p2)

    def process_candidate(self, cand):
        cand_len = len(cand)
        if cand_len < self.min_overlap:
            return False

        id_map = self.id_map
        first_entries = id_map.get(cand[0])
        if first_entries is None:
            return False
        last_entries = id_map.get(cand[-1])
        if last_entries is None:
            return False

        self._ensure_powers(cand_len)

        # Precompute candidate prefix hashes once — O(n) total,
        # then any substring hash is O(1) instead of O(n) per anchor.
        MASK_64 = self.MASK_64
        MOD_2 = self.MOD_2
        BASE = self.BASE
        powers_1 = self.powers_1
        powers_2 = self.powers_2

        cand_h1 = [0] * cand_len
        cand_h2 = [0] * cand_len
        h1 = h2 = 0
        for i, c in enumerate(cand):
            h1 = (h1 * BASE + c) & MASK_64
            h2 = (h2 * BASE + c) % MOD_2
            cand_h1[i] = h1
            cand_h2[i] = h2

        main_strings = self.main_strings
        main_hashes_1 = self.main_hashes_1
        main_hashes_2 = self.main_hashes_2
        min_overlap = self.min_overlap

        # cand prefix == anchor suffix
        prefix_match = False
        for main_id, m_pos in first_entries:
            main_len = len(main_strings[main_id])
            req_len = main_len - m_pos
            if req_len < min_overlap or cand_len < req_len:
                continue

            # anchor hash for [m_pos .. main_len-1]
            mh1 = main_hashes_1[main_id]
            mh2 = main_hashes_2[main_id]
            end = main_len - 1
            if m_pos == 0:
                target_h1, target_h2 = mh1[end], mh2[end]
            else:
                target_h1 = (mh1[end] - mh1[m_pos - 1] * powers_1[req_len]) & MASK_64
                target_h2 = (mh2[end] - mh2[m_pos - 1] * powers_2[req_len]) % MOD_2

            # candidate prefix hash [0 .. req_len-1] — O(1) lookup
            if cand_h1[req_len - 1] == target_h1 and cand_h2[req_len - 1] == target_h2:
                prefix_match = True
                break

        if not prefix_match:
            return False

        # cand suffix == anchor prefix
        last_idx = cand_len - 1
        for main_id, m_pos in last_entries:
            req_len = m_pos + 1
            if req_len < min_overlap or cand_len < req_len:
                continue

            target_h1 = main_hashes_1[main_id][m_pos]
            target_h2 = main_hashes_2[main_id][m_pos]

            # candidate suffix hash [cand_len-req_len .. cand_len-1] — O(1) lookup
            start_idx = cand_len - req_len
            if start_idx == 0:
                curr_h1 = cand_h1[last_idx]
                curr_h2 = cand_h2[last_idx]
            else:
                curr_h1 = (cand_h1[last_idx] - cand_h1[start_idx - 1] * powers_1[req_len]) & MASK_64
                curr_h2 = (cand_h2[last_idx] - cand_h2[start_idx - 1] * powers_2[req_len]) % MOD_2

            if curr_h1 == target_h1 and curr_h2 == target_h2:
                return True

        return False


def load_pickle(in_file):
    print(f"Loading {in_file}...")
    with open(in_file, "rb") as f:
        return pkl.load(f)


def write_pickle(data, out_file):
    print(f"Writing results to {out_file}...")
    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_file, "wb") as f:
        pkl.dump(data, f)


def run_rolling_hash():
    base_dir = "path/to/your/project_directory"  # Replace with your MTPLite project directory
    k = 21  # Replace with the k-mer size used for unikmer extraction

    anchor_pkl = os.path.join(base_dir, "output", "anchors.pkl")
    rid_map_pkl = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    store_dir = os.path.join(base_dir, f"barcode_{k}mers")
    idx_path = os.path.join(store_dir, "read_unikmer_map.idx")
    data_path = os.path.join(store_dir, "read_unikmer_map.data")

    out_pkl = os.path.join(base_dir, "bridges_v1.1", "direct_bridges.pkl")
    progress_every = 50000
    min_overlap = 100

    anchors = load_pickle(anchor_pkl)
    print(f"Finished loading {len(anchors)} anchor IDs")

    solver = RollingHash(anchors, min_overlap)
    print("Finished hashing anchors.")

    anchor_ids = set(anchors.keys())

    rid_int_map = load_pickle(rid_map_pkl)
    # Filter out anchors up front; build list of (read_id, read_int) to iterate
    candidate_reads = [
        (read_id, read_int)
        for read_id, read_int in rid_int_map.items()
        if read_id not in anchor_ids
    ]
    print(f"Reads after removing anchors: {len(candidate_reads):,}")

    bridge_reads = {}
    processed = 0

    print(f"Opening barcode store: {store_dir}")
    with BarcodeStore(idx_path, data_path) as store:
        for read_id, read_int in candidate_reads:
            processed += 1
            if processed % progress_every == 0:
                print(f"Processed {processed:,} reads | bridges found: {len(bridge_reads)}")

            read_len, barcodes = store.get(read_int)
            if not barcodes:
                continue

            if solver.process_candidate(barcodes):
                bridge_reads[read_id] = {
                    "barcodes": barcodes,
                    "length": read_len,
                }

    print("Finished selecting bridge reads.")
    write_pickle(bridge_reads, out_pkl)
    print(f"Saved {len(bridge_reads)} bridge reads.")

if __name__ == "__main__":
    run_rolling_hash()
