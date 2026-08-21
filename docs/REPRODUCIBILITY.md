# Reproducibility record

MTPLite is a computational research workflow. A meaningful rerun needs more
than source files: it needs a durable record of the exact data, reference,
configuration, environment, commands, resources, and results. This page lists
the minimum material to retain for every MTPLite experiment.

Raw reads, reference sequences, large binary indices, and assembly outputs are
intentionally ignored by Git. Keep their provenance in a secure experiment
location or a tracked manifest that does not expose restricted data.

## Before execution

Create a run identifier and capture the following before the first pipeline
stage.

| Category | Record |
| --- | --- |
| Repository | Remote URL, commit SHA, branch or tag, and whether this mirror differs from the canonical repository |
| Input reads | Source/accession, acquisition date, sequencing platform and chemistry if known, file format, byte size, and SHA-256 checksum |
| Target reference | Source/accession, assembly/version, sequence identifier, byte size, and SHA-256 checksum |
| Study intent | Target region, scientific question, and rationale for selection/assembly evaluation |
| Environment | `environment.yml`, resolved package list, operating system, and executable version output |
| Compute | Host or scheduler context, CPU model/count, memory, scratch/storage type, and available disk |
| Configuration | Every edited placeholder plus all numerical method settings and output prefixes |

Capture configuration before running rather than trying to reconstruct it from
changed scripts. The commands below are useful starting points:

```bash
git rev-parse HEAD
git remote -v
conda list --explicit
sha256sum input/reads.fastq reference/target.fasta
python --version
jellyfish --version
hifiasm --version
quast.py --version
```

Adapt filenames to the configured layout. If access controls prohibit sharing
input data, retain stable accessions and checksums rather than copying data into
the repository.

## During execution

Run the stages in [dependency order](HOW_TO_RUN.md). Preserve:

- the exact command for every stage, including working directory and output
  redirection;
- unedited stdout/stderr logs, job-scheduler records, start/end times, and exit
  statuses;
- the derived unikmer histogram bounds and selected k-mer count;
- stage checkpoints, selected read count, and any warnings or failures;
- resource observations for the barcode store and head-index build.

If a stage is rerun, record why, which outputs were discarded, and which
downstream stages were rebuilt. Do not overwrite a completed run's evidence
without first archiving it under a distinct run identifier.

## After execution

For the final selected-read FASTA, assembly, and QUAST result, retain paths,
byte sizes, SHA-256 checksums, and a copy of the generated QUAST report. Record
the reference-evaluation command verbatim. Summarize at least:

| Result group | Suggested measures |
| --- | --- |
| Selection | Input reads, retained unikmers, anchors, direct bridges, paths, final selected reads, and selected-read fraction |
| Assembly | Contig count, total length, largest contig, N50, and output checksum |
| Reference evaluation | Reference identifier, genome fraction, misassemblies, local misassemblies, mismatch/indel rates, duplication ratio, and QUAST version |

Interpret evaluation results in light of the target and reference. For a new
dataset, do not state that the v1.1 metrics were reproduced unless the input,
reference, configuration, and evaluation procedure are comparable and the
evidence is retained.

## Minimal manifest template

Copy this template into an experiment record and complete it before publishing
or comparing a result:

```text
run_id:
date_started:
repository_remote:
git_commit:
canonical-source relationship:

input_reads:
  source_or_accession:
  platform_and_chemistry:
  path_or_access_method:
  sha256:

reference:
  source_or_accession:
  assembly_and_sequence_id:
  sha256:

configuration:
  project_root:
  read_path:
  output_dir_and_prefix:
  nucleotide_kmer_size:
  unikmer_frequency_bounds:
  estimated_coverage:
  genome_size_bp:
  bridge_and_index_settings:
  path_search_settings:
  assembly_threads:

environment:
  conda_environment:
  resolved_packages_file:
  executable_versions:
  operating_system:

compute:
  cpu:
  memory:
  storage:
  scheduler_or_host:

artifacts:
  logs:
  selected_reads_sha256:
  assembly_sha256:
  quast_output:

results_summary:
  selected_reads:
  contigs:
  n50_bp:
  genome_fraction:
  misassemblies:
  notes:
```

## What this mirror can and cannot reproduce

This repository preserves code, environment specifications, the recorded
parameter rationale, and a saved QUAST report. It does not distribute raw
reads, reference FASTA, generated barcode/index stores, run logs, or a single
central configuration file. It supports a documented rerun after the user
supplies and configures those inputs, but cannot by itself establish bitwise
reproduction of the historical run.
