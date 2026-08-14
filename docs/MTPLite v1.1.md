### Overview

MTPLite v1.1 is a **targeted metagenomic read selection pipeline** that identifies the most informative reads from a high-coverage HiFi dataset for assembling a specific genomic region. It reduces **1.31M reads (1000x coverage) down to ~14K reads**, which are then assembled with hifiasm. The target here is _C. elegans_ chromosome 1 (~15 Mbp).

Unikmers are **single-copy k-mers** — k-mers that appear exactly once in the genome. They are identified operationally from the sequencing reads themselves by taking k-mers whose frequency falls within the interval [μ − 3σ, μ + 3σ] of the overall k-mer count distribution. K-mers appearing fewer than 5 times are excluded first as likely sequencing errors. With k=21, the k-mers within that window are single-copy with high probability.

The central design idea is the **barcode representation**: every read is compressed into an ordered list of unikmer IDs (its "barcode"). This collapses a raw 5–30 kb HiFi sequence into a sparse integer array, where two reads overlapping on the genome share a run of identical IDs in order. All downstream overlap, coverage, and pathfinding logic operates on these integer arrays rather than raw nucleotide sequences — making set operations, hashing, and graph traversal cheap.

---

### Tools and Software Used

This inventory covers the software invoked directly by the v1.1 scripts and `run_order.txt`, the third-party Python packages imported by the pipeline, and the aligner invoked internally by QUAST. Versions are taken from the original `awink` Conda environment and the saved QUAST log; where no version was recorded, this is stated explicitly.

|Category|Tool|Recorded version|Use in MTPLite v1.1|
|---|---|---:|---|
|Runtime and environment|**Conda**|Not recorded|Managed the `awink` software environment and installed the pipeline dependencies|
|Runtime and environment|**Python**|3.8.20|Ran the barcode construction, binning, anchor selection, indexing, bridge finding, path finding, and final read-selection programs|
|Runtime and environment|**Bash**|Not recorded|Ran `extractUnikmers.sh` and `assembly.sh` and orchestrated command-line stages|
|K-mer processing|**Jellyfish**|2.2.10|`count` counted canonical 21-mers, `histo` generated the frequency histogram, and `dump` extracted k-mers in the selected frequency range|
|Assembly|**hifiasm**|0.25.0|Assembled the selected HiFi reads and produced GFA assembly graphs|
|Assembly evaluation|**QUAST**|5.3.0|Compared the primary assembly with the chromosome 1 reference and generated assembly statistics, plots, and reports|
|Alignment used by QUAST|**minimap2**|Bundled QUAST version not separately recorded|Aligned the assembly to the reference during QUAST evaluation; it was invoked internally by QUAST rather than directly by the MTPLite scripts|
|Python package|**NumPy**|1.24.3|Performed weighted distribution calculations and vectorized unikmer-coverage operations|
|Python package|**Biopython**|1.78|Read FASTQ/FASTA records and wrote the final selected-read FASTA via `Bio.SeqIO` and `Bio.SeqRecord`|
|Python package|**XlsxWriter**|3.1.1|Wrote the Jellyfish k-mer frequency distribution to `prefix.xlsx` for inspection|
|Shell utilities|**GNU Awk (`awk`)**|Not recorded|Extracted unikmer sequences from `jellyfish dump` output and primary-contig sequences from hifiasm GFA output|
|Shell utilities|**`sed`, `cut`, and `xargs`**|Not recorded|Parsed the lower and upper unikmer-frequency bounds printed by `freq_distribution_kmers.py`|
|Execution and file utilities|**`nohup` and `mkdir`**|Not recorded|Ran long Python stages in the background and created output directories|

The pipeline also uses Python standard-library modules rather than additional third-party packages for its core data structures and file handling: `array`, `collections`, `csv`, `functools`, `gzip`, `heapq`, `mmap`, `operator`, `os`, `pickle`, `shutil`, `struct`, and `sys`.

The main intermediate and interchange formats are FASTQ/FASTA, GFA, CSV, XLSX, Python pickle files, and custom memory-mapped `.idx`/`.data` binary stores.

---

### Pipeline Steps

#### 1. `extractUnikmers.sh` — K-mer Frequency Filtering

- Uses **Jellyfish** to count all 21-mers in the input reads (32 threads)
- Computes mean/stddev of k-mer frequencies via `freq_distribution_kmers.py`, then extracts k-mers in the **[mean − 3σ, mean + 3σ]** range
- These "unikmers" are k-mers at expected single-copy frequency — not repetitive, not erroneous
- **Output:** `prefix.unikmers` (13,468,121 unique k-mers)

#### 2. `read_unikmer_map.py` — Barcode Each Read

- Builds a canonical 2-bit integer encoding for each unikmer (min of forward and reverse-complement hash)
- Scans every read using a **rolling hash** with simultaneous forward + reverse complement tracking
- Assigns each read a **barcode** = ordered list of unikmer IDs found in it
- Stores barcodes in a **memory-mapped flat file** (`.idx` + `.data`) for O(1) random access
  - `.idx`: 16 bytes per read (data offset 8B + num_barcodes 4B + read_len 4B)
  - `.data`: packed uint32 barcode arrays
- Also builds `rid_to_int.pkl` and `int_to_rid.pkl` mappings
- **Output:** `barcode_21mers/read_unikmer_map.{idx,data}`, `rid_maps/rid_to_int.pkl`, `rid_maps/int_to_rid.pkl`

#### 3. `bin_reads.py` — Adaptive Read Binning

Reads with very few unikmers are almost always low-information (short reads, repeat-dense reads, or error-heavy reads), while reads with an extreme number of unikmers are typically either exceptionally long HiFi reads or reads that happen to straddle unusually single-copy regions. Binning separates these populations so that downstream steps can treat them differently and so that the distribution's tails don't dominate coverage estimates.

- **Step 1 — Aggregate distribution.** Iterates over `read_stats` (`read_id → (unikmer_count, read_len)`) and builds a dict keyed by unikmer count, storing `[num_reads, base_sum]` per bucket. This is written to `unikmer_distribution.csv`.
- **Step 2 — Compute bounds.** From buckets with count ≥ 50 (to ignore the degenerate low tail), computes the frequency-weighted mean and standard deviation of the unikmer-count distribution, then derives `[μ − 3σ, μ + 3σ]`. The upper bound is rounded up to the nearest 500 to produce `adjusted_upper`. Observed: mean ≈ 9,860, σ ≈ 3,330, adjusted upper = 20,000.
- **Step 3 — Per-bin coverage budget.** Chooses a divisor based on genome size (100 / 1,000 / 10,000 / 100,000 for progressively larger genomes) and computes `min_bases = genome_size * coverage / div`. Each dynamic bin must accumulate at least this many bases before it closes — this keeps every bin roughly coverage-equivalent rather than equal-width in unikmer count.
- **Step 4 — Build bin boundaries.** Fixed singleton bins for unikmer counts 0–5 (these are small enough that a single count is its own bin). From count 6 up to `adjusted_upper`, walks sorted counts left-to-right; accumulates `agg[c][1]` (base sum) into a rolling window anchored at `window_start`; closes the bin and opens a new one whenever the window's base sum reaches `min_bases`. Everything above `adjusted_upper` falls into a single **overflow bin** (`adjusted_upper + 1`).
- **Step 5 — Assign reads.** Each read is routed to its bin: counts ≤ 5 map to the singleton bins, counts > `adjusted_upper` go to overflow, everything else goes via the `unikmer_to_bin` lookup table built in step 4.
- **Step 6 — Sort within bin.** Each bin's read list is sorted by read length descending using `operator.itemgetter(1)` and then materialised as a `dict` (preserving insertion order). Downstream steps that draw reads from "low-complexity bins" therefore always see the longest reads first.
- **Step 7 — Write outputs.** `binned_reads.pkl` holds `{bin_id: {read_id: read_len}}`; `binned_reads.csv` is a human-readable summary.

- **Output:** `bins/binned_reads.pkl`, `bins/binned_reads.csv`, `stats/unikmer_distribution.csv`

#### 4. `universe.py` — Universal K-mer Set

- Iterates all candidate reads via `BarcodeStore` flat-file accessor
- Collects every unikmer ID that appears in any read → the "universe" to be covered
- **Output:** `universe/universe.pkl` (13,468,121 unikmer IDs)

#### 5. `anchor.py` — Greedy Set-Cover Anchor Selection

Anchor selection is framed as **weighted set cover**: the "universe" is the set of all unikmer IDs, each read is a "set" (its barcode), and the goal is to pick the smallest collection of reads whose union covers every unikmer. Exact set cover is NP-hard; the classical greedy approximation is used, accelerated by a lazy-evaluated heap.

**Data structures.**
- `universe` is loaded from `universe.pkl` (a list of unikmer IDs) and reshaped into a **dense `uint8` bit-vector** of length `max_barcode + 1`. A value of `1` means "this unikmer is still uncovered", `0` means "already covered". This allows `universe[barcodes] != 0` to produce a boolean mask in a single vectorised numpy operation, and `universe[kept] = 0` to mark a whole batch of unikmers as covered in one shot.
- `FlatBarcodeQuery` wraps the barcode `.idx` + `.data` flat files with mmap. `get(read_int)` returns a `numpy.uint32` array of the read's barcode (copied out of the mmap so the caller can index safely); repeat calls are served from a per-instance dict cache (`self._cache`). `get_stats(read_int)` returns `(num_barcodes, read_len)` read directly from the `.idx` entry without touching the data file — used by the heap-seeding loop to avoid loading full barcode arrays.

**Algorithm — lazy greedy with a max-heap.**
1. **Seed.** For every mapped read, call `get_stats` to obtain `num_barcodes`. Push `(-num_barcodes, read_int, rid)` onto a list and `heapq.heapify` it. The negative sign turns Python's min-heap into a max-heap on unikmer count.
2. **Pop.** Extract the read with the largest *initial* unikmer count (i.e. its best-case score).
3. **Recompute actual score.** Fetch the read's barcode via `get_barcodes`. Compute `mask = u[barcodes] != 0` and `actual_score = count_nonzero(mask)` — the number of *currently uncovered* unikmers. Initial counts are optimistic upper bounds; the actual score monotonically decreases as the universe shrinks.
4. **Lazy check.** Peek at the next-best initial score in the heap (`-heap[0][0]`). If `actual_score ≥ next_best`, no other read can possibly beat this one, so accept it immediately: mark `u[kept] = 0`, stash the kept barcodes + read length in the `selected` dict. If `actual_score < next_best`, **re-insert** the read with its updated score and continue — do *not* recompute every other read's score.
5. **Terminate.** When the heap empties or every remaining read has `actual_score == 0`, stop.

This is the "lazy greedy" variant of Nemhauser et al.'s submodular maximisation — correctness follows from monotonicity (actual scores only ever decrease), and empirically only a small fraction of reads need their scores recomputed.

**Bookkeeping.** Each selected anchor is stored as `{"barcodes": array.array('I'), "length": read_len}` where `barcodes` is the subset of the anchor's barcode that was *newly* covered (i.e. the intersection with the uncovered universe at the moment of selection). Downstream stages use this to know which unikmers each anchor uniquely contributes.

**Result on this dataset:** 1,540 anchors cover all 13.47M unikmers (100%).

- **Output:** `output/anchors.pkl`

#### 6. `direct_bridge.py` — Rolling-Hash Overlap Detection

- For each non-anchor read, checks if its prefix overlaps an anchor's suffix AND its suffix overlaps an anchor's prefix (minimum overlap: **100 barcodes**)
- Precomputes **prefix hash arrays** for all anchors and each candidate, enabling O(1) substring hash lookups
- Uses dual rolling hashes (64-bit masked + modular 10^9+7) for collision avoidance
- **Result:** 24 direct bridge reads
- **Output:** `bridges_v1.1/direct_bridges.pkl`

#### 7. `indexer.py` — Build Head Index

The head index is an **inverted index over barcode k-mers**: given an 11-long tuple of unikmer IDs, return every read whose barcode contains that tuple contiguously. The assembler uses it to find candidate overlap partners without scanning all 1.31M reads.

**Key design choices.**
- **k = 11 over barcodes.** Note this k is *not* the nucleotide k of 21 — it is a k-mer *over the unikmer-ID sequence*. Eleven consecutive unikmer IDs is an extremely specific signature (roughly 13M^11 possible keys), so hash collisions are negligible and a match is strong evidence of a genomic overlap.
- **Asymmetric indexing depth.** Anchors are indexed at **full depth** (every valid start position, stride 1), since anchors are the source/destination of every BFS path and must be reachable from anywhere. Non-anchor reads are indexed only in their **head and tail 25%** (capped at 4,000 positions per end) with **stride 3** — sufficient to support overlap chaining while keeping the index size tractable. This asymmetry is driven by `anchor_ids` passed in from `anchors.pkl`.
- **Flat-file format, sorted by key.** Keys are variable-length tuples serialised to a fixed 44-byte big-endian struct (`>11I`), then followed by an 8-byte data offset and a 4-byte count per entry (56 bytes per entry total). All entries are sorted by the packed key bytes so that lookups become a binary search over the `.idx` mmap.

**Write path (`BarcodeIndexer` in `mode="w"`).**
1. **Accumulate.** `_kmer_acc` is a dict mapping each 44-byte key to a `bytearray` of packed uint64 read-IDs (`>Q`). `add_reads_from_barcode_store` iterates `read_ids`, fetches each barcode from `BarcodeStore`, chooses `starts` via `_get_starts(n, use_full_depth)`, and for each start position extends the bytearray under that key with the read's ID. Using `bytearray` extension avoids Python-list overhead; 8-byte aligned packing lets the data file be read back directly into an `array.array("Q")`.
2. **`_get_starts` logic.**
   - If the barcode is shorter than k, skip.
   - Full-depth (anchors): `range(0, n - k + 1)`.
   - Shallow (non-anchors): `window = min(n * 0.25, 4000)`; emit `range(0, window - k + 1, 3)` for the head and `range(n - window, n - k + 1, 3)` for the tail; dedupe with a `seen` set so short reads whose head and tail overlap don't double-count positions.
3. **Sort and write.** `_write_flat_files` sorts keys lexicographically (the big-endian packing makes byte-order sort match tuple-order sort), then streams each key → `(data_offset, count)` into `head_index.idx` and the corresponding packed read-IDs into `head_index.data`. Buffered writes at 16 MiB keep I/O efficient.

**Read path (`HeadIndexStore`).**
- Mmaps both `.idx` and `.data`. `query(kmer_key_bytes)` runs a **binary search** over the `.idx` mmap by computing midpoint byte offsets (`mid * entry_size`) and slicing out the 44-byte key for comparison. On a hit, it unpacks the trailing `>QI` to get `(offset, count)`, slices `count * 8` bytes from `.data`, and loads them into an `array.array("Q")`. Zero Python-level iteration per probe.

**Stats on this dataset:** 71.2M unique 11-mer keys, ~17.2 GB of packed read-ID data.

#### 8. `barcode_assembler.py` — BFS Path Finding Between Anchors

The assembler treats the anchor set as nodes in an implicit graph where an edge `A → B` exists if some read's **tail barcode k-mers** match the head index at a position present in `B`. Rather than materialise the graph, it performs **breadth-first search from each anchor** independently, collecting up to `top_k_per_anchor` completed paths that end at some other anchor.

**Core abstractions.**
- `CachedBarcodeStore` — a thin `BarcodeStore` wrapper with `functools.lru_cache(maxsize=20,000)` on `get(read_int)`. Because BFS revisits the same reads across many anchor searches, this cache materially reduces mmap reads.
- `HeadIndexQuery` — reimplementation of `HeadIndexStore` tailored for the assembler. It mmaps `head_index.idx`/`.data`, performs binary search over 56-byte records, and wraps the result in another `lru_cache(maxsize=50,000)` keyed on the 11-tuple. Handles endian conversion via `array.byteswap()` on little-endian machines since the on-disk packed uint64s are big-endian (`>Q`).

**Candidate generation — `_get_candidate_ids_from_suffix`.** For a given read (barcode array of length `n`):
- `window = min(n * end_fraction, end_cap) = min(n * 0.25, 4000)` — the tail window size.
- Slides a length-k (k=11) window over the last `window` positions with `stride = query_stride = 1`, building a tuple of unikmer IDs at each step.
- Queries the head index for each tuple and unions the resulting read IDs into a `candidates` set.
- Returns the full union. These are all the reads whose barcodes share any 11-consecutive-unikmer tail signature with the source read's tail — i.e. plausible downstream overlap partners.

**BFS driver — `find_top_paths_per_anchor`.**
- `target_ids` is the set of all anchors; the search succeeds whenever it reaches any anchor that isn't the source.
- For each anchor `source_id`:
  - Initialise `queue = deque([(source_id, [source_id], frozenset([source_id]))])`. The frozenset is the path's node set, used for O(1) cycle checking.
  - `seen_paths` — a set of tuple-keyed paths to avoid re-enqueuing identical extensions via different orders.
  - Pop left → level-order BFS. If `len(path) > max_depth (=15)`, skip.
  - Fetch the right endpoint's barcode via `barcode_get`; compute candidates via `get_candidates` on its tail window.
  - For each candidate not already in `path_set` and not already seen: append to `next_path`; if `next_id` is another anchor, record the path; if the path is still shorter than `max_depth`, enqueue it for further extension.
  - Count every extension toward `expanded`. Stop when `expanded ≥ max_expansions_per_source (=10,000,000)` — a safety rail against combinatorial blowup on hub-like anchors.
- Per-source results are sorted by `(-len(path), tuple(path))` (longer paths first, then lexicographic) and truncated to `top_k_per_anchor = 10`. Anchors that hit that cap are logged as `ANCHOR SATURATED`.

**Why tail-anchored search.** HiFi reads naturally overlap in a head-meets-tail fashion, so a read `r` is a plausible successor of the current endpoint `e` exactly when `r` contains some k-mer from `e`'s tail — which is what the head index (containing every read's head k-mers, plus full-depth anchor k-mers) was designed to answer. Querying `e`'s tail against an index of everyone's head region is therefore the right directionality for forward chaining.

**Output conversion.** Internal paths use integer read IDs for speed; before writing, they are remapped to the original string read IDs via `int_to_rid` and saved as `{source_read_id: [[read_id, ...], ...]}`.

**Result on this dataset:** 15,355 paths across 1,540 anchors (≈10 paths/anchor average).

- **Output:** `bridges_v1.1/bridges.pkl`

#### 9. `final_read_selection.py` — Merge & Extract

- Unions: anchors (1,540) + direct bridges (24) + reads from BFS paths + low-complexity bin reads (bins 0–200)
- Extracts sequences from original FASTQ → output FASTA
- **Result:** **14,176 selected reads** from 1.31M input (~1.08%)
- **Output:** `output/mtpv1.1.fasta`

#### 10. `assembly.sh` — HiFi Assembly

- Runs **hifiasm** on the selected reads (40 threads)
- Extracts primary contigs from GFA → final FASTA
- **Output:** `output/mtpv1.1/mtpv1.1.asm.fasta`

#### 11. QUAST — Assembly Quality Assessment

- Runs QUAST against the reference genome (`chr1.fasta`)
- **Output:** `output/quast_mtpv1.1/`

---

### Data Flow

```
hifi_1000x_chr1.fastq
       │
  [extractUnikmers] → 13.47M unikmers
       │
  [read_unikmer_map] → flat barcode store + rid maps (1.31M reads)
       │
       ├── [bin_reads] → binned reads by complexity
       ├── [universe] → 13.47M unikmer IDs
       │        │
       │   [anchor] → 1,540 anchors (100% coverage)
       │        │
       │        ├── [direct_bridge] → 24 direct bridges
       │        │
       │        └── [indexer] → head index (71.2M k-mers)
       │                │
       │           [barcode_assembler] → 15,355 BFS paths
       │
  [final_read_selection] → 14,176 reads (mtpv1.1.fasta)
       │
  [assembly.sh / hifiasm] → mtpv1.1.asm.fasta
       │
  [QUAST] → quality report
```

---

### QUAST Results (vs. reference NC_003279.8, chr1)

|Metric|Value|
|---|---|
|**Contigs**|**3**|
|**Largest contig**|**12,142,537 bp**|
|**Total length**|**15,079,659 bp**|
|**Reference length**|15,072,434 bp|
|**N50**|**12,142,537 bp**|
|**Genome fraction**|**99.993%**|
|**Misassemblies**|**0**|
|**Local misassemblies**|**0**|
|Mismatches per 100 kbp|0.14|
|Indels per 100 kbp|1.36|
|GC content|35.74% (ref: 35.75%)|
|N's per 100 kbp|0.00|
|Duplication ratio|1.001|
|Unaligned contigs|0|

---

### Assessment

**Excellent assembly quality.** From 1.31M reads, the pipeline selected 14,176 (1.08%) and produced an assembly with:

- **Zero misassemblies**
- **99.993% genome fraction** covered
- **N50 of 12.1 Mbp** (over 80% of the chromosome in one contig)
- Only 3 contigs total, 2 of which are ≥50 kb
- Extremely low error rates (0.14 mismatches/100 kbp, 1.36 indels/100 kbp)

The key algorithmic insight is the **barcode-based read representation**: by converting raw sequences into ordered unikmer-ID lists, the pipeline can reason about read overlaps and coverage in a compressed, noise-filtered space — enabling efficient greedy set cover and BFS-based path finding to select a minimal, high-quality read set for assembly.

Key v1.1 changes that improved assembly quality:
- **Larger BFS search space:** depth 15 (was 10), 10 paths per anchor (was 5), 10M expansions (was 5M) → 15,355 paths (was 7,631)
- **Shorter index k-mer size:** k=11 (was 15) with wider indexing window (25% vs 15%) → better connectivity in the overlap graph
- **Stricter direct bridging:** minimum overlap 100 barcodes (was 25) → higher-confidence bridges
- **More selected reads:** 14,176 (was 9,998) → better coverage, fewer contigs (3 vs 5), larger N50 (12.1 Mbp vs 8.0 Mbp)
