import pickle as pkl
import os
import csv
import numpy as np
from collections import defaultdict
from operator import itemgetter


def loadPickle(inFile):
    print(f"Loading {inFile}...")
    with open(inFile, 'rb') as f:
        return pkl.load(f)

def writePickle(data, outFile):
    print(f"Writing results to {outFile}...")
    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    with open(outFile, 'wb') as f:
        pkl.dump(data, f)

def writeCSV(rows, outFile):
    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    with open(outFile, mode='w', newline='') as f:
        csv.writer(f).writerows(rows)


# ── distribution helpers ────────────────────────────────────────────

def _aggregate_distribution(read_stats):
    """Return {unikmer_count: [num_reads, base_sum]} from read_stats."""
    agg = defaultdict(lambda: [0, 0])
    for unikmer_count, read_len in read_stats.values():
        entry = agg[unikmer_count]
        entry[0] += 1
        entry[1] += read_len
    return agg


def _compute_bounds(agg, start_index=50):
    """Compute mean +/- 3 sigma bounds from aggregated distribution."""
    X, F = [], []
    for count in sorted(agg):
        if count >= start_index:
            X.append(count)
            F.append(agg[count][0])

    X = np.array(X, dtype=np.float64)
    F = np.array(F, dtype=np.float64)
    total = F.sum()
    mean = np.dot(X, F) / total
    stdev = np.sqrt(np.dot(F, (X - mean) ** 2) / total)

    lower = max(int(np.floor(mean - 3 * stdev)), start_index)
    upper = int(np.ceil(mean + 3 * stdev))
    print(f"mean: {mean:.1f}  std: {stdev:.1f}  bounds: [{lower}, {upper}]")
    return lower, upper


# ── main pipeline ───────────────────────────────────────────────────

def bin_reads_pipeline(read_stats, coverage, genome_size, dist_csv, bin_pkl, bin_csv):

    # 1. Aggregate distribution (two counters per bucket, no read_id storage)
    print("1. aggregating unikmer distribution ...")
    agg = _aggregate_distribution(read_stats)
    counts_sorted = sorted(agg)

    # Write distribution CSV
    header = ["# of unikmers", "# of reads", "sum of bases (bp)"]
    rows = [[c, agg[c][0], agg[c][1]] for c in counts_sorted]
    rows.insert(0, header)
    writeCSV(rows, dist_csv)

    # 2. Compute bounds (in-memory, no CSV re-parse)
    print("2. computing distribution bounds ...")
    lower, upper = _compute_bounds(agg)
    adjusted_upper = int(np.round(upper / 500)) * 500
    print(f"   adjusted upper: {adjusted_upper}")

    # 3. Coverage threshold per bin
    print("3. determining coverage per bin ...")
    if genome_size < 500_000:
        print("   genome too small for binning")
        return
    elif genome_size < 5_000_000:
        div = 100
    elif genome_size < 50_000_000:
        div = 1_000
    elif genome_size < 500_000_000:
        div = 10_000
    else:
        div = 100_000

    min_bases = int(genome_size * coverage / div)

    # 4. Build bin boundaries
    print("4. creating bins ...")
    unikmer_to_bin = {}
    bin_ids = []

    # Fixed bins 0-5
    idx = 0
    for i in range(6):
        if i in agg:
            unikmer_to_bin[i] = i
            bin_ids.append(i)
            idx += 1

    # Dynamic bins up to adjusted_upper
    window_start = counts_sorted[idx]
    window_bases = 0

    while idx < len(counts_sorted) and counts_sorted[idx] <= adjusted_upper:
        c = counts_sorted[idx]
        unikmer_to_bin[c] = window_start
        window_bases += agg[c][1]
        idx += 1

        if window_bases >= min_bases:
            bin_ids.append(window_start)
            window_bases = 0
            if idx < len(counts_sorted):
                window_start = counts_sorted[idx]

    if window_start not in bin_ids:
        bin_ids.append(window_start)

    # Overflow bin
    overflow = adjusted_upper + 1
    bin_ids.append(overflow)

    # 5. Assign reads (collect into lists for a single fast sort)
    print("5. assigning reads to bins ...")
    bin_lists = {b: [] for b in bin_ids}

    for read_id, (unikmer_count, read_len) in read_stats.items():
        if unikmer_count <= 5:
            target = unikmer_count
        elif unikmer_count > adjusted_upper:
            target = overflow
        else:
            target = unikmer_to_bin[unikmer_count]
        bin_lists[target].append((read_id, read_len))

    # 6. Sort by length desc, convert to ordered dicts
    print("6. sorting reads in bins by length ...")
    getter = itemgetter(1)
    binned = {}
    for b in bin_ids:
        bin_lists[b].sort(key=getter, reverse=True)
        binned[b] = dict(bin_lists[b])

    # 7. Write outputs
    print("7. writing outputs ...")
    writePickle(binned, bin_pkl)

    summary = [["unikmer bin", "# of reads", "sum of bases (bp)"]]
    for b, reads in binned.items():
        summary.append([b, len(reads), sum(reads.values())])
    writeCSV(summary, bin_csv)

    print("done")


if __name__ == "__main__":
    base_dir = "path/to/your/project_directory"  # Replace with your MTPLite project directory

    read_stats_pkl = os.path.join(base_dir, "stats", "read_stats.pkl")
    dist_csv = os.path.join(base_dir, "stats", "unikmer_distribution.csv")
    bin_pkl = os.path.join(base_dir, "bins", "binned_reads.pkl")
    bin_csv = os.path.join(base_dir, "bins", "binned_reads.csv")

    coverage = 1000  # Replace with the estimated sequencing coverage
    genome_size = 14_550_000  # Replace with the estimated genome size in base pairs

    read_stats = loadPickle(read_stats_pkl)
    bin_reads_pipeline(read_stats, coverage, genome_size, dist_csv, bin_pkl, bin_csv)
