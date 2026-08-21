# MTPLite v1.1: methods and reference experiment

## Scope

MTPLite is a targeted metagenomic read-selection workflow for assembling a
specified genomic region from deep PacBio HiFi sequencing data. Rather than
assembling every input read, it identifies a small, information-rich subset and
uses hifiasm to assemble that selection.

This document describes the v1.1 implementation preserved in this mirror and
the recorded experiment on *Caenorhabditis elegans* chromosome I. It is a
methods record, not a claim that the same thresholds or performance will hold
for other sequencing technologies, genomes, coverage levels, or targets.

## Terminology

| Term | Meaning in MTPLite |
| --- | --- |
| **Unikmer** | A canonical nucleotide k-mer retained because its observed read-count frequency falls within a data-derived window expected to enrich for single-copy sequence. It is an operational proxy, not a proof of genomic uniqueness. |
| **Barcode** | The ordered sequence of unikmer IDs found in a read. |
| **Universe** | The set of unikmer IDs observed across candidate reads and targeted by anchor selection. |
| **Anchor** | A read selected by greedy set cover because it contributes previously uncovered unikmers. |
| **Direct bridge** | A non-anchor read meeting the configured barcode-overlap criterion at anchor ends. |
| **Head index** | A disk-backed inverted index from barcode-ID tuples to read IDs, used to generate overlap candidates during path search. |

## Workflow

### 1. K-mer counting and unikmer selection

`extractUnikmers.sh` uses Jellyfish to count canonical nucleotide k-mers in the
input reads and write a histogram. `freq_distribution_kmers.py` saves that
histogram to an XLSX workbook and estimates a frequency window from the
observed distribution. The lower bound is clamped to 5; selected k-mers are
then dumped from the Jellyfish count database.

In the v1.1 run, `k = 21`. Frequency filtering suppresses very-low-frequency
k-mers that are more likely to reflect sequencing error and high-frequency
k-mers that are more likely to occur in repetitive sequence. It does not
independently establish copy number, so inspect the histogram and record the
derived bounds for every new run.

### 2. Read barcoding and storage

`read_unikmer_map.py` canonicalizes selected k-mers with two-bit encoding and
scans each read with forward and reverse-complement rolling encodings. Each
read is represented by the ordered IDs of its retained k-mers.

To avoid keeping all barcodes as Python objects, the workflow writes a
memory-mapped barcode store:

```text
read_unikmer_map.idx: uint64 data offset | uint32 barcode count | uint32 read length
read_unikmer_map.data: packed uint32 barcode IDs
```

The index record is 16 bytes per read. Read-ID mapping files and per-read
barcode/read-length statistics are saved alongside the store.

### 3. Read complexity and universe construction

`bin_reads.py` groups reads by the number of retained unikmers. It calculates
distribution bounds from sufficiently populated frequency buckets and creates
coverage-aware bins; reads within a bin are ordered by length. The final
selection later includes configured low-complexity bins, so this stage retains
read populations that anchor selection alone may underrepresent.

`universe.py` scans the barcode store and collects the observed unikmer IDs.
This universe is the coverage target for anchor selection.

### 4. Greedy anchor selection

`anchor.py` frames anchor selection as set cover: each read barcode is a set of
unikmers, and the objective is to cover the observed universe with a compact
read subset. Exact minimum set cover is intractable at this scale, so the code
uses a lazy greedy approximation.

The implementation seeds a max-heap with each read's barcode count. When a
read is popped, it recomputes the number of currently uncovered unikmers using
a dense NumPy coverage vector. If its updated score is still competitive, the
read is accepted and its newly covered IDs are marked. Otherwise, it is
reinserted with its updated score. This exploits the fact that scores only
decrease as coverage accumulates.

### 5. Connectivity recovery

MTPLite adds two classes of reads to connect selected anchors.

- **Direct bridges.** `direct_bridge.py` evaluates non-anchor reads against
  anchor barcode ends. It uses dual rolling hashes and a minimum barcode
  overlap threshold to identify reads that bridge compatible anchor ends.
- **Path-derived bridges.** `indexer.py` builds a sorted, disk-backed inverted
  index over contiguous barcode-ID tuples. Anchors are indexed at full depth;
  non-anchors are indexed in bounded head and tail windows. `barcode_assembler.py`
  queries an endpoint's tail against this index and performs bounded
  breadth-first searches until another anchor is reached.

These steps generate overlap candidates from barcode agreement. They do not
replace alignment-based overlap validation or a formal assembly-graph model.

### 6. Final selection, assembly, and evaluation

`final_read_selection.py` unions anchor IDs, direct bridges, reads from
anchor-to-anchor paths, and reads from the configured low-complexity bins. It
extracts those sequences from the original read file into a FASTA file.

`assembly.sh` runs hifiasm on the selected reads and extracts primary-contig
segments from its GFA output. QUAST can then compare the output FASTA with a
target reference. Reference-based metrics should always be reported with the
exact reference version and evaluation command.

## Recorded v1.1 configuration

The values below describe the reference experiment. Values marked
"data-derived" are computed during the run and must be retained in its logs or
manifest.

| Component | Recorded setting | Role |
| --- | ---: | --- |
| Input technology | PacBio HiFi | Long-read input for targeted selection |
| Target | *C. elegans* chromosome I | Reference-evaluated region |
| Estimated coverage | 1,000× | Read-binning coverage budget |
| Estimated target size | 14,550,000 bp | Read-binning coverage budget |
| Nucleotide unikmer size | 21 | Jellyfish counting and read barcoding |
| Unikmer frequency window | Data-derived, mean ± 3 SD; lower bound ≥ 5 | Candidate unikmer filter |
| Jellyfish threads | 32 | K-mer counting |
| Direct-bridge minimum overlap | 100 barcode IDs | Direct bridge criterion |
| Barcode index key | 11 barcode IDs | Candidate overlap signature |
| Non-anchor index stride | 3 | Index-size control |
| Indexed/query end window | 25% of barcode, capped at 4,000 IDs | Candidate overlap window |
| Path-search query stride | 1 | Tail-window queries |
| Maximum path depth | 15 reads | Bounded search |
| Paths retained per anchor | 10 | Bounded search output |
| Expansion cap per anchor | 10,000,000 | Bounded search safety limit |
| Final bin inclusion | Bins 0–200 | Supplemental read selection |
| hifiasm threads | 40 | Assembly |

The recorded Conda environment is specified in
[`environment.yml`](../environment.yml). The implementation now exposes paths
and most study-specific settings as placeholders, so a new run must explicitly
configure and document all of them.

## Reference experiment results

The following result is from the included QUAST report for `mtpv1.1.asm`,
evaluated against *C. elegans* chromosome I reference NC_003279.8
(15,072,434 bp). It is a single reference-based evaluation, not an independent
benchmark or a comparison with alternative methods.

| Measure | Recorded value |
| --- | ---: |
| Input reads | 1.31 million |
| Unikmers retained | 13,468,121 |
| Anchor reads | 1,540; 100% of observed unikmer universe covered |
| Direct bridge reads | 24 |
| Anchor-to-anchor paths | 15,355 |
| Final selected reads | 14,176 (≈1.08% of input) |
| Assembly contigs | 3 |
| Largest contig / N50 | 12,142,537 bp |
| Assembly length | 15,079,659 bp |
| Genome fraction | 99.993% |
| Misassemblies / local misassemblies | 0 / 0 |
| Mismatches / indels per 100 kbp | 0.14 / 1.36 |
| Duplication ratio | 1.001 |
| N bases per 100 kbp | 0.00 |

The saved [QUAST report](../report.pdf) contains the source tables and plots
for these values, including 21 mismatches and 205 indels across the aligned
assembly. Its reported absence of misassemblies is specific to that reference,
alignment, and QUAST version.

## Interpretation and limitations

The reference experiment supports the following limited observation: for the
documented chromosome I input and configuration, MTPLite retained a small
fraction of reads while producing an assembly with high reference coverage.
It does not establish general superiority, optimality of the selected subset,
or performance on non-reference sequence.

Important limitations for new studies include:

- The unikmer criterion is frequency-based and may behave differently with
  repeat content, ploidy, heterozygosity, read quality, coverage, and target
  size.
- Parameter changes alter the computational experiment. Validate the selection
  and assembly with appropriate orthogonal evidence for the biological question.
- The head index is both storage- and memory-intensive; resource limits may
  affect feasible parameters or completion.
- Raw inputs, references, intermediate stores, execution logs, and complete
  provenance for the reference run are not distributed in this mirror.

Use [Reproducibility](REPRODUCIBILITY.md) to preserve the evidence needed for a
new or adapted run.
