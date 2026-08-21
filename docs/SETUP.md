# Setup

This guide prepares a system to run or adapt the MTPLite v1.1 research
workflow. It documents the environment that supported the recorded
*C. elegans* chromosome I study. It does not turn the mirrored scripts into a
general-purpose command-line application.

## 1. Assess suitability

MTPLite is intended for high-coverage PacBio HiFi data and targeted assembly.
The v1.1 workflow operates over large, on-disk barcode and head-index files.
The recorded index contained roughly 71.2 million unique 11-ID keys and
17.2 GB of packed read-ID data. Index construction also accumulates keys in
memory. Use fast local storage and provision RAM conservatively before starting
a deep-coverage run.

The source code is a mirror and contains user-editable placeholders. Treat a
modified parameter set or another organism as a separate experiment, not as a
direct reproduction of the v1.1 result.

## 2. Create the recorded environment

From the repository root, create the version-pinned Conda environment:

```bash
conda env create -f environment.yml
conda activate mtplite-v1.1
```

The recorded toolchain is Python 3.8.20, NumPy 1.24.3, Biopython 1.78,
XlsxWriter 3.1.1, Jellyfish 2.2.10, hifiasm 0.25.0, and QUAST 5.3.0. Bash and
standard Unix utilities (`awk`, `sed`, `cut`, and `xargs`) are also required.
QUAST invokes minimap2 internally; a standalone minimap2 version was not
recorded.

Verify the available environment before a long run:

```bash
python --version
jellyfish --version
hifiasm --version
quast.py --version
python -c "import Bio, numpy, xlsxwriter; print('Python dependencies available')"
```

If the solved environment differs from the recorded versions, save the
resolved `conda list --explicit` output with the experiment record.

## 3. Prepare data outside version control

Create or link data directories at the project root. The names below are a
convention only; they must agree with the values you place in the scripts.

```text
MTPLite/
├── input/
│   └── reads.fastq              # PacBio HiFi FASTQ or FASTA
├── reference/
│   └── target.fasta             # used only for reference-based evaluation
├── mtp_lite_v1.1/
└── ... generated workflow directories
```

The root `.gitignore` excludes raw reads, reference sequences, intermediate
binary stores, and assembly outputs. Before running, capture the input source,
accession or download date, file checksum, and reference identifier in the
[reproducibility record](REPRODUCIBILITY.md).

`read_unikmer_map.py` accepts uncompressed reads: a path ending in `.fastq` is
parsed as FASTQ and any other filename is parsed as FASTA. In contrast,
`final_read_selection.py` also supports `.fastq.gz`, `.fasta.gz`, and `.fa.gz`.
Use an uncompressed `.fastq` or `.fasta` path consistently throughout the
workflow.

## 4. Configure every entry point

The repository does not provide a central configuration file. Replace the
placeholder strings and study-specific values in the following files before
execution. Use absolute paths to avoid ambiguity. Set `DIR` in
`extractUnikmers.sh`, every Python `base_dir`, and
`barcode_assembler.py`'s `DEFAULT_BASE_DIR` to the same project artifact root.
The `--base-dir` option remains available as an explicit override for
`barcode_assembler.py`.

| File(s) | Required configuration |
| --- | --- |
| `extractUnikmers.sh` | `DIR`, `READ_FILE`, `PREFIX`, `KMER_SIZE`, and `THREADS` |
| `read_unikmer_map.py` | `base_dir`, `k`, `unikmer_file`, and `read_file` |
| `bin_reads.py` | `base_dir`, estimated `coverage`, and `genome_size` |
| `universe.py`, `anchor.py`, `direct_bridge.py`, `indexer.py` | `base_dir` and the applicable k-mer / path-search settings |
| `barcode_assembler.py` | `DEFAULT_BASE_DIR` and the applicable path-search settings; optionally override the root with `--base-dir` |
| `final_read_selection.py` | `base_dir`, `output_dir`, `output_prefix`, and `read_file` |
| `assembly.sh` | `READ_FILE`, `OUT_DIR`, `PREFIX`, and `THREADS` |

For the recorded v1.1 experiment, the key values were a nucleotide k-mer size
of 21, expected coverage of 1,000×, target genome size of 14,550,000 bp,
11-ID barcode index keys, a 100-ID direct-bridge minimum overlap, and 40
hifiasm threads. The [technical notes](<MTPLite v1.1.md#recorded-v11-configuration>)
list the complete recorded settings and their roles.

The `PREFIX` in `extractUnikmers.sh` determines the unikmer filename. Set
`read_unikmer_map.py`'s `unikmer_file` to that exact generated path. Similarly,
the final-selection and assembly prefixes must agree with each other and with
the later QUAST command. Set `final_read_selection.py`'s `output_dir` to
`<project-root>/output`: this is where `anchor.py` creates `anchors.pkl`. The
selected reads will therefore be `<project-root>/output/<prefix>.fasta`; copy
that exact path into `assembly.sh` as `READ_FILE`. `OUT_DIR` may be a separate
assembly-results directory.

## 5. Preflight review

Before starting a full run, verify that no placeholder remains:

```bash
rg -n '^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*[[:space:]]*=[[:space:]]*"(path/to/your|prefix_for_output_files|number_of_)' mtp_lite_v1.1
git status --short
git check-ignore -v input/reads.fastq reference/target.fasta
```

The first command should return no active configuration placeholders. The last
command should show that research data are ignored rather than staged. Run the
Python programs from `mtp_lite_v1.1/` so their local module imports resolve.

## 6. Re-run safety

The workflow is staged but not transaction-safe. Keep an experiment-specific
directory or archive results before re-running. In particular:

- `indexer.py` removes and rebuilds `barcode_21mers/head_index/`.
- `final_read_selection.py` empties its configured assembly-output directory
  before writing the selected-read FASTA.

Do not run dependent stages concurrently. Follow the
[dependency-ordered run guide](HOW_TO_RUN.md), and rebuild downstream outputs
after modifying an upstream configuration or input.
