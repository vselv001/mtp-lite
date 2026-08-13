from Bio import SeqIO
import pickle as pkl
import os
import array
import mmap
from barcode_store import IDX_ENTRY

# ── 2-bit kmer encoding ────────────────────────────────────────────
# A=0  C=1  G=2  T=3  → complement is XOR 3 (A↔T, C↔G)
_ENC = bytearray(128)
_ENC[ord('A')] = 0; _ENC[ord('a')] = 0
_ENC[ord('C')] = 1; _ENC[ord('c')] = 1
_ENC[ord('G')] = 2; _ENC[ord('g')] = 2
_ENC[ord('T')] = 3; _ENC[ord('t')] = 3

_IS_ACGT = bytearray(128)
for _c in b'ACGTacgt':
    _IS_ACGT[_c] = 1

# ── helpers ─────────────────────────────────────────────────────────

def loadPickle(inFile):
    print(f"Loading {inFile}...")
    with open(inFile, 'rb') as f:
        return pkl.load(f)

def writePickle(data, outFile):
    print(f"Writing results to {outFile}...")
    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    with open(outFile, 'wb') as f:
        pkl.dump(data, f)

def create_directory(path):
    os.makedirs(path, exist_ok=True)

def build_read_id_maps(read_file, out_str_to_int=None, out_int_to_str=None):
    str_to_int = {}
    int_to_str = {}
    next_id = 0
    fmt = "fastq" if read_file.endswith(".fastq") else "fasta"
    with open(read_file) as handle:
        for record in SeqIO.parse(handle, fmt):
            rid = record.id
            if rid not in str_to_int:
                str_to_int[rid] = next_id
                int_to_str[next_id] = rid
                next_id += 1
            if next_id % 100000 == 0:
                print(f"Mapped {next_id:,} read IDs")
    print(f"Total reads mapped: {next_id:,}")
    if out_str_to_int is not None:
        writePickle(str_to_int, out_str_to_int)
    if out_int_to_str is not None:
        writePickle(int_to_str, out_int_to_str)
    return str_to_int, int_to_str

def get_unikmer_map(unikmer_file, k):
    """Build unikmer dict keyed by canonical 2-bit integer.

    Each kmer is stored under min(forward_int, rc_int) so the inner loop
    only needs one lookup per position.
    """
    mask = (1 << (2 * k)) - 1
    enc = _ENC
    unikmer_map = {}
    with open(unikmer_file) as f:
        for idx, line in enumerate(f):
            kmer_str = line.strip()
            # encode forward
            fwd = 0
            for ch in kmer_str:
                fwd = (fwd << 2) | enc[ord(ch)]
            # encode RC via complement + reverse of 2-bit groups
            comp = fwd ^ mask
            rev = 0
            tmp = comp
            for _ in range(k):
                rev = (rev << 2) | (tmp & 3)
                tmp >>= 2
            unikmer_map[min(fwd, rev)] = idx
    return unikmer_map

# ── writer ──────────────────────────────────────────────────────────

_DATA_BUF_LIMIT = 16 * 1024 * 1024   # flush every 16 MB

def barcode_reads_flat(unikmer_map, rid_to_int, read_file, k, idx_path, data_path):
    """Barcode every read and write results into .idx + .data flat files."""

    num_reads = len(rid_to_int)
    read_stats = {}
    entry_size = IDX_ENTRY.size            # 16
    pack_entry = IDX_ENTRY.pack

    # Pre-allocate index (zeroed → unfilled entries read as offset=0, count=0, len=0)
    idx_size = num_reads * entry_size
    with open(idx_path, 'wb') as f:
        f.truncate(idx_size)

    idx_fd = os.open(idx_path, os.O_RDWR)
    idx_mm = mmap.mmap(idx_fd, idx_size)

    data_f = open(data_path, 'wb', buffering=_DATA_BUF_LIMIT)
    data_offset = 0

    # Rolling-hash constants
    mask = (1 << (2 * k)) - 1
    rc_shift = 2 * (k - 1)
    enc = _ENC
    is_acgt = _IS_ACGT

    uniq_get = unikmer_map.get
    rid_lookup = rid_to_int.__getitem__
    fmt = "fastq" if read_file.endswith(".fastq") else "fasta"

    with open(read_file, buffering=4 * 1024 * 1024) as handle:
        counter = 0
        for record in SeqIO.parse(handle, fmt):
            counter += 1
            read_id = record.id
            read_int = rid_lookup(read_id)

            seq_bytes = record.seq.encode('ascii') if isinstance(record.seq, str) else bytes(str(record.seq), 'ascii')
            read_len = len(seq_bytes)

            # ── rolling hash with rolling RC ────────────────────────
            barcodes = []
            bc_append = barcodes.append
            fwd = 0
            rev = 0
            valid_run = 0

            for i in range(read_len):
                b = seq_bytes[i]
                if is_acgt[b]:
                    base = enc[b]
                    fwd = ((fwd << 2) | base) & mask
                    rev = (rev >> 2) | ((3 ^ base) << rc_shift)
                    valid_run += 1
                else:
                    fwd = 0
                    rev = 0
                    valid_run = 0

                if valid_run >= k:
                    hit = uniq_get(fwd if fwd <= rev else rev)
                    if hit is not None:
                        bc_append(hit)

            num_barcodes = len(barcodes)

            # write barcode data
            if num_barcodes:
                barcode_arr = array.array('I', barcodes)
                data_f.write(barcode_arr.tobytes())

            # write index entry (direct mmap store)
            pos = read_int * entry_size
            idx_mm[pos:pos + entry_size] = pack_entry(data_offset, num_barcodes, read_len)
            data_offset += num_barcodes * 4

            read_stats[read_id] = (num_barcodes, read_len)

            if counter % 100000 == 0:
                print(f"processed reads: {counter}")

    data_f.close()
    idx_mm.flush()
    idx_mm.close()
    os.close(idx_fd)

    return read_stats

# ── entry point ─────────────────────────────────────────────────────

def mapper_executor():
    k = 21
    base_dir = os.path.join(os.path.expanduser("~"), "MTPLite")

    unikmer_file = os.path.join(base_dir, "jellyfish_data", "prefix.unikmers")
    read_file = os.path.join(base_dir, "input", "hifi_1000x_chr1.fastq")
    rid_map_file = os.path.join(base_dir, "rid_maps", "rid_to_int.pkl")
    int_rid_map_file = os.path.join(base_dir, "rid_maps", "int_to_rid.pkl")
    stats_pkl = os.path.join(base_dir, "stats", "read_stats.pkl")

    out_dir = os.path.join(base_dir, f"barcode_{k}mers")
    create_directory(out_dir)
    idx_path = os.path.join(out_dir, "read_unikmer_map.idx")
    data_path = os.path.join(out_dir, "read_unikmer_map.data")

    print("loading unikmers")
    unikmer_map = get_unikmer_map(unikmer_file, k)

    if os.path.exists(rid_map_file) and os.path.exists(int_rid_map_file):
        print("loading rid_to_int")
        rid_to_int = loadPickle(rid_map_file)
    else:
        print("building read id maps")
        rid_to_int, _ = build_read_id_maps(read_file, rid_map_file, int_rid_map_file)

    print("processing reads")
    read_stats = barcode_reads_flat(
        unikmer_map,
        rid_to_int,
        read_file,
        k,
        idx_path,
        data_path,
    )

    print("saving read stats pickle")
    writePickle(read_stats, stats_pkl)
    print("DONE")

if __name__ == "__main__":
    mapper_executor()
