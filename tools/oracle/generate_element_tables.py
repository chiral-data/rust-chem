#!/usr/bin/env python3
"""Prints `core::elements`'s COVALENT_RADII and VDW_RADII tables (#388).

Generated from the pinned OpenBabel rather than typed in, so the Rust table is
reproducible and its provenance is one command:

    docker run --rm --entrypoint python3 -v "$PWD":/src chem-oracle \\
        /src/tools/oracle/generate_element_tables.py

OpenBabel fills a gap with a placeholder rather than saying "unknown": a
covalent radius past curium (Cordero 2008 stops at Z = 96) and a vdW radius of
exactly 2.0. Both become `None`, since an invented number is worse than none.
"""

from openbabel import openbabel as ob

LAST_Z = 118
LAST_CORDERO_Z = 96
VDW_UNKNOWN = 2.0


def covalent(z: int):
    return ob.GetCovalentRad(z) if z <= LAST_CORDERO_Z else None


def vdw(z: int):
    r = ob.GetVdwRad(z)
    return None if r == VDW_UNKNOWN else r


def rust_array(name: str, radius) -> str:
    # Index 0 is the placeholder every `core::elements` table carries for
    # `Element::UNKNOWN`.
    cells = ["None"] + [
        "None" if (r := radius(z)) is None else f"Some({r!r})" for z in range(1, LAST_Z + 1)
    ]
    body = "\n".join(
        "    " + ", ".join(cells[i : i + 8]) + "," for i in range(0, len(cells), 8)
    )
    return f"pub(crate) const {name}: [Option<f64>; {LAST_Z + 1}] = [\n{body}\n];"


if __name__ == "__main__":
    print(rust_array("COVALENT_RADII", covalent))
    print()
    print(rust_array("VDW_RADII", vdw))
