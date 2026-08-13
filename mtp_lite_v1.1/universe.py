import pickle as pkl
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from barcode_store import BarcodeStore

def build_universe(candidates, rid_to_int, ctx):
    universe = set()

    print("Processing candidates...")
    for i, read_id in enumerate(candidates, 1):
        read_int = rid_to_int.get(read_id)
        if read_int is None:
            continue
        _, barcodes = ctx.get(read_int)
        if not barcodes:
            continue
        universe.update(barcodes)

        if i % 100000 == 0:
            print(f"  Processed {i} reads...")

    return universe

def loadPickle(inFile):
    print(f"Loading {inFile}...")
    with open(inFile, 'rb') as f:
        return pkl.load(f)

def writePickle(data, outFile):
    print(f"Writing results to {outFile}...")
    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    with open(outFile, 'wb') as f:
        pkl.dump(data, f)

def main():
    base_dir = os.path.join(os.path.expanduser("~"), "MTPLite")
    k = 21

    idx_path = os.path.join(base_dir, f"barcode_{k}mers", "read_unikmer_map.idx")
    data_path = os.path.join(base_dir, f"barcode_{k}mers", "read_unikmer_map.data")
    candidates_pkl = os.path.join(base_dir, "stats", "read_stats.pkl")
    rid_int_map_pkl = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    out_pkl = os.path.join(base_dir, "universe", "universe.pkl")

    # 1. Load data
    candidates = loadPickle(candidates_pkl)
    rid_to_int = loadPickle(rid_int_map_pkl)

    # 2. Initialize barcode query context
    ctx = BarcodeStore(idx_path, data_path)

    # 3. Execution
    try:
        print("Building universe...")
        universe = build_universe(candidates, rid_to_int, ctx)
        writePickle(universe, out_pkl)
    finally:
        # 4. Manual Cleanup
        # Using finally ensures the connection closes even if the build fails
        ctx.close()

    print("Done.")

if __name__ == "__main__":
    main()
