#!/usr/bin/env python3
"""Generates permuted_axes.ccp4 (#341).

Confirmed empirically (not assumed): gemmi's Python API can write a
genuinely permuted `MAPC`/`MAPR`/`MAPS` order, via `Ccp4Map.set_header_i32`
after `update_ccp4_header()` -- no hand-assembled bytes needed, unlike DCD.
`update_ccp4_header()` first sets sane defaults (MAPC=1, MAPR=2, MAPS=3,
i.e. columns/rows/sections mapped straight onto crystallographic a/b/c);
patching those three words afterwards only changes what the file's
already-fixed storage order *means*, not the bytes underneath -- exactly
`io/ccp4.rs`'s own documented contract ("MAPC/MAPR/MAPS only permute which
canonical axis each file dimension is").

This reproduces the exact case `io/ccp4.rs`'s own hand-computed test
(`test_a_permuted_axis_non_cubic_cell_matches_independently_computed_
values`) already pins: file dims NX=2 (columns), NY=1 (rows), NZ=3
(sections); MAPC=3 (columns represent c), MAPR=1 (rows -> a), MAPS=2
(sections -> b); cell (10, 20, 30, 90, 90, 90); values 1..6 stored
columns-fastest.

Run: python3 generate_permuted_axes.py
"""

from pathlib import Path

import gemmi

grid = gemmi.FloatGrid(2, 1, 3)  # (columns, rows, sections) = (NX, NY, NZ)
grid.set_unit_cell(gemmi.UnitCell(10.0, 20.0, 30.0, 90.0, 90.0, 90.0))
grid.spacegroup = gemmi.SpaceGroup("P1")
# Columns-fastest: (col, sec) = (0,0),(1,0),(0,1),(1,1),(0,2),(1,2) -> 1..6.
grid.set_value(0, 0, 0, 1.0)
grid.set_value(1, 0, 0, 2.0)
grid.set_value(0, 0, 1, 3.0)
grid.set_value(1, 0, 1, 4.0)
grid.set_value(0, 0, 2, 5.0)
grid.set_value(1, 0, 2, 6.0)

m = gemmi.Ccp4Map()
m.grid = grid
m.update_ccp4_header()
m.set_header_i32(17, 3)  # MAPC = c
m.set_header_i32(18, 1)  # MAPR = a
m.set_header_i32(19, 2)  # MAPS = b

if __name__ == "__main__":
    m.write_ccp4_map(str(Path(__file__).with_name("permuted_axes.ccp4")))
