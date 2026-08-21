# MTPLite documentation

This documentation describes the MTPLite v1.1 research workflow as it is
mirrored in this repository. It distinguishes the recorded *C. elegans*
chromosome I experiment from guidance for adapting the workflow to a new
dataset.

## Reading paths

### Reproduce or adapt a run

1. [Setup](SETUP.md) — create the recorded environment and replace every
   configuration placeholder.
2. [Reproducibility guide](REPRODUCIBILITY.md) — capture the data, parameter,
   software, hardware, and command record for the run.
3. [How to run](HOW_TO_RUN.md) — execute the stages in dependency order and
   check the expected artifacts.
4. [Technical notes](<MTPLite v1.1.md>) — review the algorithm and parameter
   rationale before changing a setting.

### Understand the reference study

- [Technical notes and results](<MTPLite v1.1.md>) explain the barcode-based
  method and place the reported result in context.
- [report.pdf](../report.pdf) is the retained QUAST output for the v1.1
  assembly evaluated against *C. elegans* chromosome I reference NC_003279.8.
- [Original run order](run_order.md) preserves a concise historical command
  record. It is not current configuration guidance.

## Documentation principles

- **Configured, not turnkey.** The scripts contain intentionally visible
  placeholders for paths, filenames, prefixes, and resource settings. Complete
  the setup checklist before running any stage. In particular, use one project
  artifact root and pass the selected-read FASTA from final selection to
  `assembly.sh` exactly as written.
- **Evidence-scoped.** The reported metrics apply only to the recorded input,
  reference, toolchain, and parameters.
- **Reproducible by record.** Large data and generated artifacts are excluded
  from version control. Preserve their provenance, checksums, commands, logs,
  tool versions, and output checksums with each experiment.
- **Mirror-aware.** This repository is not the source of truth for the
  implementation. See [CONTRIBUTING.md](../CONTRIBUTING.md) for maintenance
  expectations.
