"""gemmi as a structural oracle for mmCIF/PDB (#224).

Unlike RDKit/OpenBabel (see `oracles/__init__.py`), gemmi does not answer a
SMILES-shaped question — there is no bare-SMILES `parses`/`identity` here, so
it does not implement the `Oracle` dataclass and is not part of `load()`'s
list. It answers a structural question instead: given mmCIF or PDB text,
what atoms, unit cell, chains and residues does a correct, independent
reader see. That is also why it is loaded and sanity-checked separately
(`load_gemmi` below) rather than folded into `load()`'s SMILES-sanity-gated
list — `Oracle.parses("CCO")` means nothing to a toolkit that only reads
structures.

RDKit has no mmCIF/PDB support at all, and OpenBabel's mmCIF writer emits a
file this very module reads **zero atoms** from — neither can judge this
format, which is the whole reason this module exists.

That second claim used to read "a different, incompatible small-molecule
crystallography dialect (fractional coordinates, no chains)". Measured against
the pinned openbabel-wheel 3.1.1.21 (#258), with and without a `CRYST1` cell,
it writes `Cartn_x`/`Cartn_y`/`Cartn_z` both times — the fractional half was
never true. What it actually does is drop the chain, occupancy and
`B_iso_or_equiv` columns, and the result is not merely awkward to read:

    source pdb        3 atoms, chain A, b=[42.5, 37.25, 55.0]
    obabel -> mmcif   0 atoms, no chains, no b
    chem   -> mmcif   3 atoms, chain A, b=[42.5, 37.25, 55.0]

Its **PDB** writer is a separate defect, and the one `check_pdb` pins: it
zeroes `B_iso_or_equiv`'s PDB equivalent on every write while preserving the
occupancy column beside it, which destroys the per-atom confidence a predicted
structure stores there.
"""

from typing import Callable, NamedTuple, Optional

import gemmi as _gemmi

#: A minimal but complete `_atom_site` loop. Verified against the installed
#: gemmi (0.7.5): a loop missing `auth_asym_id`/`auth_seq_id` reads with
#: zero models rather than an error — silently, so this fixture is exactly
#: as complete as a real file, not the smallest string that merely avoids
#: an exception.
SANITY = """data_sanity
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 C C . LIG A 1 A LIG . 0.000 0.000 0.000 1.00 20.00 1
"""

_FORMAT = {".cif": _gemmi.CoorFormat.Mmcif, ".pdb": _gemmi.CoorFormat.Pdb}


class Summary(NamedTuple):
    """The structural facts this oracle compares against `chem`'s own read.

    Residue identity is `(chain, seqid, name)` rather than anything
    positional — atom *order* within a residue is not a claim either reader
    makes, only which residues exist and what they are named/numbered.
    """

    atom_count: int
    cell: Optional[tuple[float, float, float, float, float, float]]
    chain_ids: tuple[str, ...]
    residues: tuple[tuple[str, int, str], ...]


def summarize(text: str, suffix: str = ".cif") -> Optional[Summary]:
    """Reads `text` (mmCIF by default; pass `suffix=".pdb"` for PDB) and
    summarises what gemmi saw, or `None` if gemmi could not read it at all
    (either it raised, or it read zero models — gemmi does both depending
    on exactly what's missing, so both count as "could not read this").
    """
    try:
        structure = _gemmi.read_structure_string(text, format=_FORMAT[suffix])
    except Exception:
        return None
    if len(structure) == 0:
        return None

    model = structure[0]
    atom_count = sum(len(residue) for chain in model for residue in chain)
    cell = (
        (
            structure.cell.a,
            structure.cell.b,
            structure.cell.c,
            structure.cell.alpha,
            structure.cell.beta,
            structure.cell.gamma,
        )
        if structure.cell.a > 0
        else None
    )
    chain_ids = tuple(chain.name for chain in model)
    residues = tuple(
        (chain.name, residue.seqid.num, residue.name)
        for chain in model
        for residue in chain
    )
    return Summary(atom_count, cell, chain_ids, residues)


class Sites(NamedTuple):
    """The per-atom values [`Summary`] deliberately leaves out.

    Separate from `Summary` rather than more fields on it: `check_mmcif`
    compares whole tuples with `!=`, so widening `Summary` would silently
    change a check this does not belong to.

    Rounded to the two decimals both formats store, so the comparison is
    about the value and not about float formatting.
    """

    occupancies: tuple[float, ...]
    b_factors: tuple[float, ...]


def sites(text: str, suffix: str = ".pdb") -> Optional[Sites]:
    """Per-atom occupancy and B-factor, in file order.

    The check that matters for #258: `Carries` is presence-based, so a writer
    emitting a constant for every atom satisfies every mask assertion in the
    crate. Only comparing values catches it -- which is exactly OpenBabel's
    PDB defect.
    """
    try:
        structure = _gemmi.read_structure_string(text, format=_FORMAT[suffix])
    except Exception:
        return None
    if len(structure) == 0:
        return None

    atoms = [atom for chain in structure[0] for residue in chain for atom in residue]
    return Sites(
        tuple(round(atom.occ, 2) for atom in atoms),
        tuple(round(atom.b_iso, 2) for atom in atoms),
    )


#: The PDB counterpart of `SANITY`, complete in the same way -- fixed columns
#: through the element symbol, so a column-offset bug in the fixture cannot be
#: mistaken for one in what is being tested.
SANITY_PDB = """\
ATOM      1  N   ALA A   1      11.104  13.207   2.428  0.80 42.50           N
END
"""


def load_gemmi() -> Callable[..., Optional[Summary]]:
    """Confirms gemmi can read something, the same way `oracles.load()`
    confirms RDKit/OpenBabel can before trusting either. Returns
    `summarize` itself (so a caller need not re-import this module), or
    raises if gemmi is broken rather than merely unable to parse a name it
    doesn't recognise.
    """
    if summarize(SANITY) is None:
        raise RuntimeError(
            "gemmi cannot read a trivial mmCIF fixture, so it is broken "
            "rather than strict — refusing to report its answers as findings"
        )
    # The PDB path was in `_FORMAT` from the start but never exercised until
    # #258. Gated too, so a `check_pdb` failure is unambiguously about chem
    # rather than about this reader.
    if summarize(SANITY_PDB, ".pdb") is None:
        raise RuntimeError(
            "gemmi cannot read a trivial PDB fixture, so it is broken "
            "rather than strict — refusing to report its answers as findings"
        )
    measured = sites(SANITY_PDB)
    if measured != Sites((0.80,), (42.50,)):
        raise RuntimeError(
            f"gemmi reports {measured} for a fixture stating 0.80/42.50, so its "
            "per-atom values cannot be trusted as a reference"
        )
    return summarize
