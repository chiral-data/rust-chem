#!/usr/bin/env python3
"""Generates binary_big_endian.ply (#341).

`io/ply.rs`'s `ascii` and `binary_little_endian` cases are already
exercised cheaply by round-tripping through this crate's own writer, which
never needs a generator: it's already the reproducible source. Only
`binary_big_endian` was a real gap -- this crate never writes it (module
doc: "This crate never writes big-endian PLY"), and the one existing test
for it (`test_a_real_plyfile_produced_big_endian_fixture_reads_correctly`)
embedded a 218-byte literal directly in the test source, produced once by
hand and never committed as a file or a script.

Confirmed empirically (not assumed): `trimesh` (already pinned for #340)
cannot write `binary_big_endian` PLY at all -- its exporter's only
encodings are `ascii`/`binary`/`binary_little_endian` (read directly from
its own source). `plyfile` is the tool that actually produced the original
byte literal, and is the one used here -- a second, oracle-only pinned
dependency (GPL-3), the same "external tool in tools/oracle/'s own image,
never vendored into chem" posture OpenBabel (GPL-2) already established.

Geometry matches the byte literal it replaces exactly: 3 vertices at
(0,0,0)/(1,0,0)/(0,1,0), a shared (0,0,1) normal, red/green/blue vertex
colors, one triangle face.

Run: python3 generate_binary_big_endian.py
"""

from pathlib import Path

import numpy as np
from plyfile import PlyData, PlyElement

VERTEX_DTYPE = [
    ("x", "f4"), ("y", "f4"), ("z", "f4"),
    ("nx", "f4"), ("ny", "f4"), ("nz", "f4"),
    ("red", "u1"), ("green", "u1"), ("blue", "u1"),
]
VERTICES = np.array(
    [
        (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 255, 0, 0),
        (1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0, 255, 0),
        (0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0, 0, 255),
    ],
    dtype=VERTEX_DTYPE,
)
FACES = np.array([([0, 1, 2],)], dtype=[("vertex_indices", "i4", (3,))])


def build() -> PlyData:
    return PlyData(
        [PlyElement.describe(VERTICES, "vertex"), PlyElement.describe(FACES, "face")],
        text=False,
        byte_order=">",
    )


if __name__ == "__main__":
    build().write(str(Path(__file__).with_name("binary_big_endian.ply")))
