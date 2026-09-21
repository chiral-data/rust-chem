"""cpptraj as the NCTRAJ reference oracle (#340).

Amber's own analysis tool, not a second opinion: the issue's own framing is
"the definition of whether a file is a valid Amber trajectory," the same
status RDKit had for commonchem in #229. Used here purely as a validity
check on `chem`'s NCTRAJ *output* -- does Amber's own tool accept it and
report back the same frame count and coordinates -- not as a fixture
author; `oracles/mdanalysis.py` already authors this format's fixtures
with a real writer, and cpptraj's own interface is a command script, not a
Python library, better spent judging output than generating it.

No standalone `cpptraj` package exists on any conda channel (confirmed by
survey) -- only the full `ambertools` conda-forge bundle provides it,
installed into its own environment by the Dockerfile's dedicated stage. A
real, disclosed ~2.3GB image-size cost, paid once at build time, isolated
into its own layer so an unrelated pip pin bump does not repeat it.

cpptraj also refuses to treat a coordinate-only NCTRAJ as its own topology
(confirmed by survey: "Could not determine format of topology" against the
trajectory file itself) -- every real Amber workflow pairs a trajectory
with a separate topology, so this hands it a minimal PDB naming only atom
count, which is all `read_frames` needs and all `chem`'s own NCTRAJ writer
ever states in the first place.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import NamedTuple, Optional

import gemmi


def binary() -> Path:
    """Where the AmberTools conda environment's `cpptraj` lives.

    The Dockerfile sets `CPPTRAJ_BIN`; falls back to the environment's
    default path for a developer who built that stage locally under the
    same name.
    """
    return Path(os.environ.get("CPPTRAJ_BIN", "/opt/conda/envs/ambertools/bin/cpptraj"))


def _topology_pdb(n_atoms: int) -> str:
    """A minimal PDB cpptraj accepts as a topology -- one residue, `n_atoms`
    generic carbons. States nothing this check needs beyond atom count,
    the same floor a coordinate-only NCTRAJ trajectory itself states.
    """
    lines = [
        f"ATOM  {i + 1:>5}  C{i + 1:<3}RES A   1       0.000   0.000   0.000  1.00  0.00           C"
        for i in range(n_atoms)
    ]
    return "\n".join(lines) + "\nEND\n"


class Frames(NamedTuple):
    positions: tuple[tuple[tuple[float, float, float], ...], ...]
    box: Optional[tuple[float, float, float, float, float, float]]


def read_frames(nctraj_path: Path, n_atoms: int) -> Optional[Frames]:
    """Asks cpptraj to read `nctraj_path` and report every frame's atoms
    back as a multi-`MODEL` PDB, or `None` if cpptraj declined the file
    outright (not valid Amber NetCDF by its own definition) or produced
    nothing gemmi can read.

    PDB, not cpptraj's own `mdcrd`, because gemmi already has a robust PDB
    reader this module can reuse (`oracles/gemmi.py` already trusts it for
    `check_pdb`) rather than a second, hand-rolled parser for `mdcrd`'s
    flat, wrapped-at-ten-columns text stream.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        topology_path = tmp / "topology.pdb"
        topology_path.write_text(_topology_pdb(n_atoms))
        out_path = tmp / "out.pdb"
        script = (
            f"parm {topology_path}\n"
            f"trajin {nctraj_path}\n"
            f"trajout {out_path} pdb\n"
            "run\n"
        )
        result = subprocess.run(
            [str(binary())],
            input=script,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not out_path.exists():
            return None
        pdb_text = out_path.read_text()

    try:
        structure = gemmi.read_structure_string(pdb_text, format=gemmi.CoorFormat.Pdb)
    except Exception:
        return None
    if len(structure) == 0:
        return None

    positions = tuple(
        tuple((atom.pos.x, atom.pos.y, atom.pos.z) for chain in model for residue in chain for atom in residue)
        for model in structure
    )
    cell = structure.cell
    box = (
        (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
        if cell.a > 0
        else None
    )
    return Frames(positions, box)


def load_cpptraj() -> None:
    """Confirms the `cpptraj` binary actually runs inside the image before
    its answers are trusted -- the same discipline as every other oracle's
    `load_*` gate, adapted for a CLI tool rather than a Python import.
    """
    try:
        result = subprocess.run(
            [str(binary()), "--version"], capture_output=True, text=True, timeout=30
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        raise RuntimeError(
            f"cpptraj did not run ({e}) — refusing to report its answers as findings"
        ) from e
    if result.returncode != 0 or "CPPTRAJ" not in result.stdout:
        raise RuntimeError(
            "cpptraj --version did not identify itself, so it is broken "
            "rather than strict — refusing to report its answers as findings"
        )
