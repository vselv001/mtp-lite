# How to run MTPLite v1.1

Follow this guide after completing the [setup guide](SETUP.md). The commands
below reproduce the order implemented by the pipeline; run each stage only
after its prerequisites finish successfully.

## Before you start

Confirm all of the following:

- The Conda environment in [`environment.yml`](../environment.yml) is active.
- `input/hifi_1000x_chr1.fastq` and `reference/chr1.fasta` exist at the
  project root expected by your configured scripts.
- The project-root paths inside the Python and Bash entry points have been
  updated consistently if you are not using the original environment.
- You have sufficient local disk and RAM for the barcode and head-index files.
- You have read the [re-run safety notes](SETUP.md#6-re-run-safety).

## Execute the pipeline

Run the following from the repository root. The `mkdir` command creates a
local log directory that is ignored by Git. Each redirection preserves the
program's exit code, unlike a backgrounded pipeline whose failure may be easy
to miss.

```bash
cd mtp_lite_v1.1
mkdir -p logs

# 1. Count canonical 21-mers and write the frequency-filtered unikmer list.
bash extractUnikmers.sh prefix 21 > logs/extract_unikmers.log 2>&1

# 2. Build read-ID maps, a memory-mapped barcode store, and read statistics.
python -u read_unikmer_map.py > logs/read_unikmer_map.log 2>&1

# 3. Bin reads by their unikmer-barcode counts.
python -u bin_reads.py > logs/bin_reads.log 2>&1

# 4. Build the observed unikmer universe.
python -u universe.py > logs/universe.log 2>&1

# 5. Choose greedy set-cover anchors.
python -u anchor.py > logs/anchor.log 2>&1

# 6. Detect reads directly bridging anchor ends.
python -u direct_bridge.py > logs/direct_bridge.log 2>&1

# 7. Build the barcode head index.
python -u indexer.py > logs/indexer.log 2>&1

# 8. Search for bounded paths between anchors.
python -u barcode_assembler.py > logs/barcode_assembler.log 2>&1

# 9. Merge selected read IDs and write the assembly FASTA input.
python -u final_read_selection.py > logs/final_read_selection.log 2>&1

# 10. Assemble selected reads and extract primary contigs.
bash assembly.sh > logs/assembly.log 2>&1
```

Use `tail -f logs/<stage>.log` in another terminal to monitor a stage. If a
command fails, inspect its log and resolve the underlying issue before moving
to the next numbered command.

## Evaluate the assembly

Return to the repository root and run QUAST against the target reference:

```bash
cd ..
quast.py output/mtpv1.1/mtpv1.1.asm.fasta \
  -r reference/chr1.fasta \
  -o output/quast_mtpv1.1
```

QUAST writes its HTML report, text reports, and plots to
`output/quast_mtpv1.1/`. The v1.1 evidence report committed to this repository
is [report.pdf](../report.pdf).

## Outputs and checkpoints

| Step | Confirm this output exists before continuing |
| --- | --- |
| 1 | `jellyfish_data/prefix.unikmers` |
| 2 | `barcode_21mers/read_unikmer_map.idx`, `barcode_21mers/read_unikmer_map.data`, `stats/read_stats.pkl` |
| 3 | `bins/binned_reads.pkl` and `stats/unikmer_distribution.csv` |
| 4 | `universe/universe.pkl` |
| 5 | `output/anchors.pkl` |
| 6 | `bridges_v1.1/direct_bridges.pkl` |
| 7 | `barcode_21mers/head_index/head_index.idx` and `head_index.data` |
| 8 | `bridges_v1.1/bridges.pkl` |
| 9 | `output/mtpv1.1.fasta` |
| 10 | `output/mtpv1.1/mtpv1.1.asm.fasta` |
| QUAST | `output/quast_mtpv1.1/` |

## Interpreting the reference experiment

For the recorded *C. elegans* chromosome I run, MTPLite selected 14,176 of
1.31 million input reads and produced a three-contig assembly. QUAST reported
99.993% genome fraction, an N50 of 12,142,537 bp, and zero misassemblies. See
[report.pdf](../report.pdf) for the complete result tables and plots.

These values validate the supplied dataset and parameterization only. Record
the reference version, software versions, input checksums, command logs, and
any configuration changes before comparing a new run with this result.

## Common failures

| Symptom | Likely cause and response |
| --- | --- |
| `FileNotFoundError` for an input or intermediate | Check that every `base_dir`, `DIR`, `READ_FILE`, and `OUT_DIR` value points to the same project root; then rerun the prerequisite stage. |
| `jellyfish`, `hifiasm`, or `quast.py` is not found | Activate the MTPLite Conda environment and rerun the preflight version checks. |
| The indexer is killed or runs out of memory | Run on a machine with more RAM and local disk; the head-index accumulator is intentionally memory-intensive. |
| No anchors or paths are found | Inspect the unikmer-frequency window, barcode outputs, and k-mer/coverage settings before changing downstream thresholds. |
| Unexpected results after a rerun | Check whether the head-index or assembly-output directories were cleared, then rebuild downstream stages in dependency order. |
