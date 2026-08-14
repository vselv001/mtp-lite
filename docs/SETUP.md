# Setup guide

This guide prepares a system to reproduce the MTPLite v1.1 workflow. It
covers the recorded environment and the configuration required by the current
scripts; it does not change pipeline parameters or claim portability that the
implementation does not yet provide.

## 1. System requirements

The recorded v1.1 run used Linux, Bash, Conda, fast local storage, and the
following software.

| Component | Recorded version | Purpose |
| --- | ---: | --- |
| Python | 3.8.20 | Runs the pipeline programs. |
| NumPy | 1.24.3 | Vectorized unikmer-coverage operations. |
| Biopython | 1.78 | FASTA/FASTQ parsing and FASTA writing. |
| XlsxWriter | 3.1.1 | Jellyfish histogram workbook output. |
| Jellyfish | 2.2.10 | Canonical k-mer counting and unikmer extraction. |
| hifiasm | 0.25.0 | HiFi read assembly. |
| QUAST | 5.3.0 | Reference-based assembly evaluation. |

The workflow also calls `awk`, `sed`, `cut`, and `xargs`. QUAST invokes
minimap2 internally; its separately installed version was not recorded.

The head-index construction is the main resource bottleneck. In the supplied
experiment it produced roughly 71.2 million unique 11-ID keys and 17.2 GB of
packed read-ID data, in addition to the in-memory accumulator used during
index construction. Provision substantial RAM and fast local scratch/storage
before starting a deep-coverage run.

## 2. Create the recorded software environment

From the repository root, create the Conda environment described in the
version-controlled [`environment.yml`](../environment.yml):

```bash
conda env create -f environment.yml
conda activate mtplite-v1.1
```

For an existing environment, install the same packages and versions from that
file. The versions reflect the successful v1.1 experiment; validate any
version changes before treating their results as comparable.

Verify the executable and Python dependencies before investing in a long run:

```bash
python --version
jellyfish --version
hifiasm --version
quast.py --version
python -c "import Bio, numpy, xlsxwriter; print('Python dependencies available')"
```

## 3. Prepare the project and data layout

The current configuration expects this project-root layout:

```text
MTPLite/
├── input/
│   └── hifi_1000x_chr1.fastq
├── reference/
│   └── chr1.fasta
├── mtp_lite_v1.1/
└── ... pipeline-generated directories ...
```

`hifi_1000x_chr1.fastq` is the PacBio HiFi read input. `chr1.fasta` is used
only by the final QUAST evaluation. The read-mapping stage accepts
uncompressed FASTQ or FASTA according to the file extension; the final
selection stage also supports gzip-compressed FASTQ/FASTA. The configured
v1.1 run uses an uncompressed FASTQ.

The root `.gitignore` excludes `input/`, `reference/`, and every large
generated output directory. Do not commit raw sequencing data or binary
indices. For each experiment, retain a provenance record with the input source,
download date, checksums, reference identifier, and any parameter changes.

## 4. Configure the current scripts

MTPLite v1.1 is not yet parameterized by a config file or command-line
arguments. The Python entry points derive their root directory as
`~/MTPLite`, while the shell scripts retain the original absolute path
`/home/vselv001/MTPLite`. Configure all of them to one consistent project root
before running on a different account or checkout path.

| Files | What to update |
| --- | --- |
| `read_unikmer_map.py`, `bin_reads.py`, `universe.py`, `anchor.py`, `direct_bridge.py`, `indexer.py`, `barcode_assembler.py`, `final_read_selection.py` | The `base_dir` assignment in each entry point. |
| `extractUnikmers.sh` | `DIR`, which determines input and Jellyfish-output paths. |
| `assembly.sh` | `READ_FILE` and `OUT_DIR`. |

Keep the configured filename prefixes consistent across stages. In particular,
the first stage writes `jellyfish_data/prefix.unikmers`, and
`read_unikmer_map.py` expects that exact filename.

The documented v1.1 parameters target a 14.55 Mbp region at 1,000× coverage.
Before using another target, review the k-mer length, expected genome size,
coverage, binning settings, bridge threshold, and index/BFS limits in the
[technical design notes](<MTPLite v1.1.md>). Changing those values makes a new
experiment and should be recorded with its output.

## 5. Preflight check

Before a full run, confirm that the configured input is readable and that the
repo will not accidentally stage large artifacts:

```bash
git status --short
git check-ignore -v input/hifi_1000x_chr1.fastq reference/chr1.fasta
```

Then move to the script directory. The Python programs import local modules,
so run them from `mtp_lite_v1.1/` as shown in the [run guide](HOW_TO_RUN.md).

## 6. Re-run safety

Most stages overwrite their own files. Two operations remove existing
directories before recreating them:

- `indexer.py` removes `barcode_21mers/head_index/`.
- `final_read_selection.py` empties `output/mtpv1.1/`.

Archive results you want to preserve before rerunning either stage. Do not run
dependent stages concurrently: the numbered sequence in the run guide is a
dependency order.
