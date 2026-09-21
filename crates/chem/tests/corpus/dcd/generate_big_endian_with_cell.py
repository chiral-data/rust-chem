#!/usr/bin/env python3
"""Generates big_endian_with_cell.dcd (#341).

No tool anywhere in this repo's oracle image can write a big-endian DCD:
MDAnalysis's `DCDWriter` (`tools/oracle/oracles/mdanalysis.py`) has no
endianness parameter in its own C extension and always writes native/
little-endian. CHARMM/NAMD-produced big-endian files exist in the wild
(the whole reason `io::dcd`'s reader auto-detects endianness at all — see
its module doc), but nothing here can *author* one, so this hand-assembles
the bytes directly, mirroring the exact layout `io/dcd.rs`'s own test-only
`RawDcdBuilder` already builds in Rust (`test_a_hand_built_big_endian_
file_with_a_unit_cell_reads_correctly`) — just committed as a real file
and a real script instead of thrown away after one test run.

Layout (CHARMM DCD, Fortran unformatted records: each record is its own
byte length as an i32, the body, then that same length again):
  header record (84 bytes): b"CORD", NSET, ISTART, NSAVC, 4 zeros, NAMNF,
    DELTA (f32) + has-cell flag (i32), 8 zeros, VERSION=24
  title record: NTITLE=0
  atoms record: NATOM
  per frame: [cell record: 6 f64, A/gamma/B/beta/alpha/C order, if any],
    then one record per axis (x, then y, then z), each NATOM f32 values

Run: python3 generate_big_endian_with_cell.py
"""

import struct
from pathlib import Path

BIG_ENDIAN = ">"
NATOMS = 4
# [A, gamma, B, beta, alpha, C], angles as cosines (all 90 degrees -> 0.0),
# the modern convention -- the same six numbers the Rust test already
# hand-verified.
CELL = [30.0, 0.0, 25.0, 0.0, 0.0, 20.0]
POSITIONS = [(float(i), 0.0, 0.0) for i in range(NATOMS)]


def record(body: bytes) -> bytes:
    length = struct.pack(BIG_ENDIAN + "i", len(body))
    return length + body + length


def header_bytes() -> bytes:
    hdr = b"CORD"
    hdr += struct.pack(BIG_ENDIAN + "i", 1)  # NSET
    hdr += struct.pack(BIG_ENDIAN + "i", 0)  # ISTART
    hdr += struct.pack(BIG_ENDIAN + "i", 1)  # NSAVC
    hdr += struct.pack(BIG_ENDIAN + "5i", 0, 0, 0, 0, 0)
    hdr += struct.pack(BIG_ENDIAN + "i", 0)  # NAMNF
    hdr += struct.pack(BIG_ENDIAN + "f", 1.0)  # DELTA (CHARMM: f32)
    hdr += struct.pack(BIG_ENDIAN + "i", 1)  # has a unit cell
    hdr += struct.pack(BIG_ENDIAN + "8i", *([0] * 8))
    hdr += struct.pack(BIG_ENDIAN + "i", 24)  # CHARMM version
    assert len(hdr) == 84, len(hdr)
    return hdr


def build() -> bytes:
    out = record(header_bytes())
    out += record(struct.pack(BIG_ENDIAN + "i", 0))  # NTITLE = 0
    out += record(struct.pack(BIG_ENDIAN + "i", NATOMS))

    out += record(struct.pack(BIG_ENDIAN + "6d", *CELL))
    for axis in range(3):
        values = [p[axis] for p in POSITIONS]
        out += record(struct.pack(BIG_ENDIAN + f"{NATOMS}f", *values))
    return out


if __name__ == "__main__":
    Path(__file__).with_name("big_endian_with_cell.dcd").write_bytes(build())
