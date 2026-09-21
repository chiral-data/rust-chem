#!/usr/bin/env python3
"""Generates triclinic_scaled.lammpstrj (#341).

No tool needed to author this one -- LAMMPS dump is plain ASCII text, so
the only question is getting the numbers right, not the bytes. Two
existing `io/lammpstrj.rs` tests each hand-verify one half of this trap in
isolation (a triclinic box alone; scaled coordinates alone); this combines
both into one fixture, since a scaled coordinate's conversion to Cartesian
genuinely depends on the box being triclinic, not just orthogonal --
`real = origin + xs*a + ys*b + zs*c` where `a/b/c` are the triclinic edge
vectors, not `real = origin + frac * (hi - lo)`.

The box is the exact one `test_a_triclinic_box_matches_independently_
verified_values` already hand-verified: true edges lx=ly=lz=10, tilts
xy=2, xz=1, yz=0.5 (bounding-box-shifted per LAMMPS's own convention into
the `0 13 2 / 0 10.5 1 / 0 10 0.5` stated below). The fractional position
(0.5, 0.25, 0.1) is the same one `test_scaled_coordinates_convert_using_
the_box` already used for an orthogonal box; here it resolves to
`origin + 0.5*(10,0,0) + 0.25*(2,10,0) + 0.1*(1,0.5,10)` =
`(5 + 0.5 + 0.1, 2.5 + 0.05, 1.0)` = `(5.6, 2.55, 1.0)` -- hand-derived
here and cross-checked against MDAnalysis's `DumpReader` during this
story's research, the same discipline both source tests already used.

Run: python3 generate_triclinic_scaled.py
"""

from pathlib import Path

TEXT = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS xy xz yz pp pp pp
0 13 2
0 10.5 1
0 10 0.5
ITEM: ATOMS id type xs ys zs
1 1 0.5 0.25 0.1
"""

if __name__ == "__main__":
    Path(__file__).with_name("triclinic_scaled.lammpstrj").write_text(TEXT)
