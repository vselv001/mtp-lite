# Historical run-order record

This page preserves the concise sequence associated with the original v1.1
experiment. It is included as provenance, not as an executable protocol. The
historical shell invocation used positional arguments and an absolute project
path; the current mirrored shell scripts use configured variables instead.

Use [How to run](HOW_TO_RUN.md) for current execution guidance and
[Setup](SETUP.md) to configure the entry points first.

## Recorded stage order

1. Extract frequency-filtered unikmers with `extractUnikmers.sh`.
2. Build read-ID maps, read barcodes, and per-read statistics with
   `read_unikmer_map.py`.
3. Bin reads with `bin_reads.py`.
4. Build the observed unikmer universe with `universe.py`.
5. Select anchors with `anchor.py`.
6. Detect direct bridges with `direct_bridge.py`.
7. Build the head index with `indexer.py`.
8. Search anchor-to-anchor paths with `barcode_assembler.py`.
9. Select and export reads with `final_read_selection.py`.
10. Assemble selected reads with `assembly.sh`.
11. Evaluate the primary-contig assembly with QUAST against the target
    reference.

The original run launched several Python stages with `nohup`, but the stages
are dependency-ordered and should not be treated as parallel tasks. Confirm
each output and its log before continuing to the next stage.

## Historical context

The original notes describe the goal as reducing contig count in the selected
read assembly. The documented final configuration retained 14,176 reads and
produced three primary contigs under the recorded reference-based evaluation.
See [methods and reference results](<MTPLite v1.1.md>) and
[report.pdf](../report.pdf) for the evidence and limitations.
