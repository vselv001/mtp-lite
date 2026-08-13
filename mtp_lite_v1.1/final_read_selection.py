import os
import pickle as pkl
import gzip
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_pickle(path):
    print(f"Loading {path}...")
    with open(path, "rb") as handle:
        return pkl.load(handle)

def resolve_existing_path(*paths):
    for path in paths:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"None of these files exist: {paths}")

def write_pickle(data, path):
    print(f"Writing results to {path}...")
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "wb") as handle:
        pkl.dump(data, handle)

def _open_reads_file(reads_path):
    if reads_path.endswith(".gz"):
        return gzip.open(reads_path, "rt")
    return open(reads_path, "r")


def _detect_seqio_format(reads_path):
    lower = reads_path.lower()
    if lower.endswith(".fastq") or lower.endswith(".fastq.gz"):
        return "fastq"
    if lower.endswith(".fasta") or lower.endswith(".fa") or lower.endswith(".fasta.gz") or lower.endswith(".fa.gz"):
        return "fasta"
    raise ValueError(f"Unsupported reads file format: {reads_path}")


def ensure_empty_dir(path):
    os.makedirs(path, exist_ok=True)
    for name in os.listdir(path):
        target = os.path.join(path, name)
        if os.path.isdir(target):
            shutil.rmtree(target)
        else:
            os.remove(target)


def write_selected_reads_fasta(reads_path, selected_reads, out_fasta):
    print(f"Writing selected reads FASTA to {out_fasta}...")
    out_dir = os.path.dirname(out_fasta)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    selected_set = set(selected_reads)
    written = 0
    seqio_fmt = _detect_seqio_format(reads_path)

    with _open_reads_file(reads_path) as handle, open(out_fasta, "w") as out_handle:
        for record in SeqIO.parse(handle, seqio_fmt):
            if record.id in selected_set:
                SeqIO.write(SeqRecord(record.seq, id=record.id, description=""), out_handle, "fasta")
                written += 1

    print(f"Wrote {written:,} selected reads to FASTA.")


def collect_from_anchors(anchor_obj):
    if not isinstance(anchor_obj, dict):
        return set()
    return set(anchor_obj.keys())


def collect_from_direct_bridges(direct_obj):
    if not isinstance(direct_obj, dict):
        return set()
    return set(direct_obj.keys())


def collect_from_bridges(bridges_obj):
    selected = set()
    if not isinstance(bridges_obj, dict):
        return selected

    for anchor_id, paths in bridges_obj.items():
        if isinstance(anchor_id, str):
            selected.add(anchor_id)

        if not isinstance(paths, list):
            continue

        for path in paths:
            if not isinstance(path, (list, tuple)):
                continue
            for read_id in path:
                if isinstance(read_id, str):
                    selected.add(read_id)

    return selected


def collect_from_bins(binned_obj, bins=range(6)):
    selected = set()
    if not isinstance(binned_obj, dict):
        return selected

    for bin_id in bins:
        reads_in_bin = binned_obj.get(bin_id, {})
        if not isinstance(reads_in_bin, dict):
            continue
        for read_id in reads_in_bin.keys():
            if isinstance(read_id, str):
                selected.add(read_id)

    return selected


def main():
    base_dir = os.path.join(os.path.expanduser("~"), "MTPLite")
    output_dir = os.path.join(base_dir, "output")
    mtp_output_dir = os.path.join(output_dir, "mtpv1.1")
    helper_dir = os.path.join(base_dir, "helper")

    read_file = os.path.join(base_dir, "input", "hifi_1000x_chr1.fastq")
    anchors_pkl = os.path.join(output_dir, "anchors.pkl")
    bridges_dir = os.path.join(base_dir, "bridges_v1.1")
    bridges_pkl = os.path.join(bridges_dir, "bridges.pkl")
    direct_bridges_pkl = os.path.join(bridges_dir, "direct_bridges.pkl")
    binned_reads_pkl = os.path.join(base_dir, "bins", "binned_reads.pkl")
    final_reads_fasta = os.path.join(output_dir, "mtpv1.1.fasta")

    ensure_empty_dir(mtp_output_dir)

    anchors = load_pickle(anchors_pkl)
    bridges = load_pickle(bridges_pkl)
    direct_bridges = load_pickle(direct_bridges_pkl)
    binned_reads = load_pickle(binned_reads_pkl)

    final_reads = set()
    final_reads.update(collect_from_anchors(anchors))
    final_reads.update(collect_from_bridges(bridges))
    final_reads.update(collect_from_direct_bridges(direct_bridges))
    final_reads.update(collect_from_bins(binned_reads, bins=range(201)))

    final_reads = sorted(final_reads)
    write_selected_reads_fasta(read_file, final_reads, final_reads_fasta)

    print(f"Selected total reads: {len(final_reads):,}")
    print(f"FASTA saved at: {final_reads_fasta}")


if __name__ == "__main__":
    main()
