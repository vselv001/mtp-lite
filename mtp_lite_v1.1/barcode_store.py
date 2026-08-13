import os
import array
import struct
import mmap

# Index entry: offset into .data (8 B) | num_barcodes (4 B) | read_len (4 B)
IDX_ENTRY = struct.Struct(">QII")          # 16 bytes per read


class BarcodeStore:
    """Read-only O(1) accessor for the .idx / .data flat files."""

    __slots__ = ("_idx_mm", "_data_mm", "_idx_fd", "_data_fd")

    def __init__(self, idx_path, data_path):
        self._idx_fd = os.open(idx_path, os.O_RDONLY)
        idx_size = os.fstat(self._idx_fd).st_size
        self._idx_mm = mmap.mmap(self._idx_fd, idx_size, access=mmap.ACCESS_READ)

        self._data_fd = os.open(data_path, os.O_RDONLY)
        data_size = os.fstat(self._data_fd).st_size
        self._data_mm = mmap.mmap(self._data_fd, data_size, access=mmap.ACCESS_READ) if data_size else None

    def get(self, read_int):
        """Return (read_len, array('I', barcodes)) for a given integer read id."""
        off = read_int * IDX_ENTRY.size
        data_offset, num_barcodes, read_len = IDX_ENTRY.unpack_from(self._idx_mm, off)
        if num_barcodes == 0:
            return read_len, array.array('I')
        raw = self._data_mm[data_offset:data_offset + num_barcodes * 4]
        barcodes = array.array('I')
        barcodes.frombytes(raw)
        return read_len, barcodes

    def close(self):
        self._idx_mm.close()
        os.close(self._idx_fd)
        if self._data_mm is not None:
            self._data_mm.close()
        os.close(self._data_fd)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
