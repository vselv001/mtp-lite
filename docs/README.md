# MTPLite documentation

This directory contains the documentation needed to reproduce or adapt the
MTPLite v1.1 research workflow. Read the pages in the following order for a
new installation.

| Document | When to use it |
| --- | --- |
| [Setup guide](SETUP.md) | Install the recorded software environment, prepare inputs, and configure the current path-bound scripts. |
| [How to run](HOW_TO_RUN.md) | Execute the pipeline in dependency order and evaluate an assembly with QUAST. |
| [Technical design notes](<MTPLite v1.1.md>) | Understand the barcode representation, selection algorithms, and v1.1 parameters. |
| [Original run order](run_order.md) | Consult the concise command sequence preserved from the initial experiment. |
| [QUAST result report](../report.pdf) | Review the evidence for the reported v1.1 assembly results. |

## Documentation scope

The checked-in code is a mirror of a research workflow. It is intentionally
documented as implemented: entry-point scripts use dataset-specific values and
hard-coded project paths rather than command-line options. The setup guide
identifies every path that must be adapted before running on another machine.

Raw reads, reference sequences, intermediate binary stores, and assembly
outputs are excluded from version control by the root [`.gitignore`](../.gitignore).
Keep their provenance, checksums, commands, and software versions alongside
each experiment outside the repository or in a tracked experiment manifest.
