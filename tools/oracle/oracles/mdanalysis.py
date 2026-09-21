"""MDAnalysis as the trajectory oracle for XTC, TRR, DCD and NCTRAJ (#340).

Like `gemmi`/`meeko`, this does not implement the `Oracle` dataclass and is
not part of `load()`'s list -- it answers "do these positions, frame by
frame, agree within this format's own real precision", not a SMILES-shaped
question. It is also the harness's first *tolerance* comparison rather than
an identity one: `crates/chem/src/io/format.rs` solved this exact problem
for its own internal fidelity matrix (`frame_tolerance`/`positions_match`,
#339), and `run.py`'s own copies of those two functions mirror the same two
numbers (`0.01` Angstrom for XTC, `1e-3` elsewhere) so the Rust and Python
sides of this question can't quietly diverge.

MDAnalysis can both read and write XTC, TRR, DCD and the Amber NetCDF
(NCTRAJ) format -- confirmed by survey, and empirically here via
`load_mdanalysis`'s own round trip. LAMMPS trajectory is the fifth
`Kind::Frames` format this milestone adds, and is handled entirely in
`run.py`: MDAnalysis's `DumpReader` is read-only (no `DumpWriter` exists),
so that one fixture is generated as plain text there instead of authored
through this module.
"""

import tempfile
from pathlib import Path
from typing import NamedTuple, Optional

import MDAnalysis as mda
import numpy as np

#: One reader class per format this module can read. Low-level
#: (`mda.coordinates.<FMT>.<FMT>Reader`) rather than a full `Universe` --
#: confirmed by survey that these need no topology at all to yield
#: positions and a box, which is all this check ever asks of them.
#:
#: `lammpstrj` is read-only here on purpose: MDAnalysis's `DumpReader`
#: exists but has no `DumpWriter` counterpart (confirmed by survey), so it
#: belongs in this dict (`read_frames` needs it) but not in
#: `WRITABLE_FORMATS` below (`write_fixture` cannot support it).
_READER = {
    "xtc": mda.coordinates.XTC.XTCReader,
    "trr": mda.coordinates.TRR.TRRReader,
    "dcd": mda.coordinates.DCD.DCDReader,
    "nctraj": mda.coordinates.TRJ.NCDFReader,
    "lammpstrj": mda.coordinates.LAMMPS.DumpReader,
}

#: The four formats MDAnalysis can write. LAMMPS trajectory is deliberately
#: absent -- see the module doc.
WRITABLE_FORMATS = ("xtc", "trr", "dcd", "nctraj")

#: `mda.Writer`'s generic factory picks a writer class from the *filename's
#: extension*, not from a format string -- and does not recognise `chem`'s
#: own `.nctraj` extension (confirmed by survey: `TypeError: No trajectory
#: or frame writer for format 'NCTRAJ'`). `write_fixture` passes this
#: explicitly instead, so the fixture's on-disk suffix can stay `.nctraj`
#: (matching what a caller later hands to `chem convert --from`/`--to`,
#: which cares about neither suffix) rather than needing a second,
#: MDAnalysis-only naming convention.
_WRITE_FORMAT = {"xtc": "XTC", "trr": "TRR", "dcd": "DCD", "nctraj": "NCDF"}


class Frames(NamedTuple):
    """Per-frame positions plus the (single, constant) box this check's own
    fixtures always use. A real trajectory's box can vary frame to frame;
    nothing here needs that generality.
    """

    positions: tuple[tuple[tuple[float, float, float], ...], ...]
    box: Optional[tuple[float, float, float, float, float, float]]


def read_frames(path: Path, fmt: str) -> Optional[Frames]:
    """Reads any of the four writable formats back, or `None` if MDAnalysis
    could not read this file at all.
    """
    reader_cls = _READER[fmt]
    try:
        positions = []
        box = None
        with reader_cls(str(path)) as reader:
            for timestep in reader:
                positions.append(
                    tuple(tuple(float(c) for c in p) for p in timestep.positions)
                )
                box = tuple(float(c) for c in timestep.dimensions)
    except Exception:
        return None
    if not positions:
        return None
    return Frames(tuple(positions), box)


def write_fixture(
    path: Path,
    fmt: str,
    positions_per_frame: tuple[tuple[tuple[float, float, float], ...], ...],
    box: tuple[float, float, float, float, float, float],
) -> None:
    """Authors a small trajectory in `fmt` via MDAnalysis's own writer.

    `fmt` must be one of `WRITABLE_FORMATS`. No topology is built beyond a
    bare atom count -- every one of this check's comparisons is positions
    and box only, the same floor `crates/chem/src/io/format.rs`'s own
    `Kind::Frames` fixtures use.
    """
    n_atoms = len(positions_per_frame[0])
    universe = mda.Universe.empty(n_atoms, trajectory=True)
    with mda.Writer(str(path), n_atoms=n_atoms, format=_WRITE_FORMAT[fmt]) as writer:
        for positions in positions_per_frame:
            universe.atoms.positions = np.array(positions, dtype=np.float32)
            universe.atoms.dimensions = np.array(box, dtype=np.float32)
            writer.write(universe.atoms)


#: A trivial two-atom, one-frame fixture -- enough to prove MDAnalysis can
#: write its own XTC and read it back, the same "broken rather than strict"
#: gate `load_gemmi`/`load_meeko` already apply to their own toolkits.
_SANITY_POSITIONS = (((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),)
_SANITY_BOX = (10.0, 10.0, 10.0, 90.0, 90.0, 90.0)


def load_mdanalysis() -> None:
    """Confirms MDAnalysis can write and read back its own XTC before its
    answers are trusted. Unlike `load_gemmi`/`load_meeko`, there is no
    single callable to hand back -- callers use `write_fixture`/
    `read_frames` directly -- so this is a pure gate, raising rather than
    returning on success.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "sanity.xtc"
        write_fixture(path, "xtc", _SANITY_POSITIONS, _SANITY_BOX)
        back = read_frames(path, "xtc")
    if back is None or back.positions != _SANITY_POSITIONS:
        raise RuntimeError(
            "MDAnalysis cannot write and read back its own two-atom XTC "
            "file, so it is broken rather than strict — refusing to report "
            "its answers as findings"
        )
