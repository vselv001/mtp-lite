# Run MTPLite v1.1

Follow this guide only after completing [Setup](SETUP.md). It describes the
dependency order of the research workflow, not a single command that can run
without configuration. Replace every placeholder in the scripts and save a
run record before beginning.

## Before running

Confirm that:

- the `mtplite-v1.1` Conda environment is active;
- all paths, prefixes, k-mer sizes, and resource settings are configured
  consistently across the entry points;
- input reads and, if evaluation is planned, a target reference are readable;
- local storage and RAM are adequate for barcode and head-index construction;
- the intended data, software, parameter, and hardware record has been
  captured as described in [Reproducibility](REPRODUCIBILITY.md).

Run the placeholder scan from the setup guide before proceeding. All commands
below are launched from the script directory. The shell scripts use their
configured variables; they do not take the historical positional arguments
shown in [run_order.md](run_order.md).

Before stage 1, verify the two required handoffs: `final_read_selection.py`
uses `<project-root>/output` for `output_dir`, and `assembly.sh` uses the
resulting `<project-root>/output/<prefix>.fasta` as `READ_FILE`. Its `OUT_DIR`
is the separate directory for hifiasm outputs.

```bash
cd mtp_lite_v1.1
mkdir -p logs

# 1. Count canonical nucleotide k-mers and extract the unikmer list.
bash extractUnikmers.sh > logs/extract_unikmers.log 2>&1

# 2. Build read-ID maps, a memory-mapped barcode store, and read statistics.
python -u read_unikmer_map.py > logs/read_unikmer_map.log 2>&1

# 3. Bin reads by barcode complexity.
python -u bin_reads.py > logs/bin_reads.log 2>&1

# 4. Construct the observed unikmer universe.
python -u universe.py > logs/universe.log 2>&1

# 5. Select greedy set-cover anchors.
python -u anchor.py > logs/anchor.log 2>&1

# 6. Identify high-confidence direct bridge reads.
python -u direct_bridge.py > logs/direct_bridge.log 2>&1

# 7. Build the on-disk barcode head index.
python -u indexer.py > logs/indexer.log 2>&1

# 8. Search bounded paths between anchors.
python -u barcode_assembler.py > logs/barcode_assembler.log 2>&1

# 9. Merge selected read IDs and emit the assembly FASTA.
python -u final_read_selection.py > logs/final_read_selection.log 2>&1

# 10. Assemble the selected reads and extract primary contigs.
bash assembly.sh > logs/assembly.log 2>&1
```

Execute one stage at a time. Inspect the corresponding log and checkpoint
before advancing. If an upstream stage is rerun or its configuration changes,
rebuild every dependent output.

## Stage checkpoints

The directory and filename prefixes below are examples from the recorded v1.1
configuration. Substitute the values you configured.

| Step | Purpose | Expected artifact |
| --- | --- | --- |
| 1 | Count k-mers and filter unikmers | `jellyfish_data/<prefix>.unikmers` |
| 2 | Build barcode store and mappings | `barcode_<k>mers/read_unikmer_map.{idx,data}`, `stats/read_stats.pkl` |
| 3 | Bin reads | `bins/binned_reads.pkl`, `stats/unikmer_distribution.csv` |
| 4 | Build unikmer universe | `universe/universe.pkl` |
| 5 | Select anchors | `<project-root>/output/anchors.pkl` |
| 6 | Detect direct bridges | `bridges_v1.1/direct_bridges.pkl` |
| 7 | Build head index | `barcode_<k>mers/head_index/head_index.{idx,data}` |
| 8 | Find bridge paths | `bridges_v1.1/bridges.pkl` |
| 9 | Write selected reads | `<project-root>/output/<prefix>.fasta` |
| 10 | Assemble reads | `<assembly_dir>/<prefix>.asm.fasta` |

Use `tail -f logs/<stage>.log` from another terminal to monitor long-running
stages. A missing intermediate normally indicates either an incomplete
prerequisite or inconsistent configuration—not a downstream algorithm error.

## Evaluate an assembly

After assembly, return to the repository root and adapt the paths to your
configured assembly and reference:

```bash
cd ..
quast.py <assembly_dir>/<prefix>.asm.fasta \
  -r reference/target.fasta \
  -o <output_dir>/quast_<prefix>
```

QUAST writes HTML, text reports, and plots under the output directory. Record
the command, QUAST version, reference identifier/checksum, and full report
alongside the run. The included [report.pdf](../report.pdf) is the result for
the reference v1.1 experiment only.

## Troubleshooting

| Observation | Check first |
| --- | --- |
| A path, input, or intermediate cannot be found | Compare every configured root, input path, prefix, and output directory; then rebuild the missing prerequisite. |
| Jellyfish, hifiasm, or QUAST is unavailable | Activate the Conda environment and repeat the version checks in [Setup](SETUP.md). |
| The indexer is killed or exhausts memory | Stop before running downstream stages; use a machine with more RAM and fast local disk, then record the changed hardware context. |
| No anchors or bridge paths are produced | Inspect the unikmer histogram/window and barcode artifacts before relaxing downstream bridge or path thresholds. |
| A rerun produces inconsistent outputs | Confirm which output directories were cleared and rebuild all downstream stages with a single, documented configuration. |

## Interpreting the v1.1 run

The reference experiment selected 14,176 of 1.31 million input reads and
assembled three contigs. Against the recorded chromosome I reference, QUAST
reported 99.993% genome fraction and no reported misassemblies. These values
are evidence for the captured configuration, not acceptance criteria for a new
dataset. Use the [technical notes](<MTPLite v1.1.md#reference-experiment-results>)
and [report.pdf](../report.pdf) when comparing results.
