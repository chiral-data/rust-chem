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
different, incompatible small-molecule crystallography dialect (fractional
coordinates, no chains) — neither can judge this format, which is the whole
reason this module exists.
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
    return summarize
