"""gemmi as a structural oracle for mmCIF, PDB, CIF-core (#320) and, as of
#340, CCP4/MRC (#224 for the pattern).

BinaryCIF was surveyed for #340 and stays out of scope: the pinned (and
latest installable) `gemmi==0.7.5` raises on both `gemmi.cif.read` and
`gemmi.read_structure` for the already-committed `bcif/*.bcif` fixtures
(#319) -- confirmed directly, not assumed from the version number alone.
The gap `crates/chem/tests/corpus/README.md` already discloses for
BinaryCIF stands; nothing here reads it.

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

import numpy as _np

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


class SmallMoleculeSummary(NamedTuple):
    """The structural facts [`check_cif_core`] compares against `chem`'s own
    read of the small-molecule dictionary (#320) -- fractional sites rather
    than [`Summary`]'s Cartesian ones, and no chain/residue concept at all in
    this dictionary.
    """

    cell: Optional[tuple[float, float, float, float, float, float]]
    spacegroup_hm: Optional[str]
    sites: tuple[tuple[str, str, float, float, float, float, Optional[float]], ...]


def summarize_small_molecule(text: str) -> Optional[SmallMoleculeSummary]:
    """The CIF-core counterpart of [`summarize`] (#320) -- a different gemmi
    entry point entirely (`make_small_structure_from_block`, not
    `read_structure_string`), since this is a different dictionary, not a
    parameter variant of mmCIF/PDB reading.
    """
    try:
        doc = _gemmi.cif.Document()
        doc.parse_string(text)
        block = doc.sole_block()
        st = _gemmi.make_small_structure_from_block(block)
    except Exception:
        return None
    if len(st.sites) == 0:
        return None

    cell = (
        (st.cell.a, st.cell.b, st.cell.c, st.cell.alpha, st.cell.beta, st.cell.gamma)
        if st.cell.a > 0
        else None
    )
    sites = tuple(
        (
            site.label,
            site.type_symbol,
            round(site.fract.x, 4),
            round(site.fract.y, 4),
            round(site.fract.z, 4),
            round(site.occ, 2),
            round(site.u_iso, 4) if site.u_iso else None,
        )
        for site in st.sites
    )
    return SmallMoleculeSummary(cell, st.spacegroup_hm or None, sites)


#: A minimal but complete small-molecule CIF -- one atom, a cell, a space
#: group, esd on a cell length, so the sanity check exercises the same three
#: things `check_cif_core` actually relies on.
SANITY_CIF_CORE = """data_sanity
_cell_length_a 4.9134(2)
_cell_length_b 4.9134
_cell_length_c 5.4052
_cell_angle_alpha 90.00
_cell_angle_beta 90.00
_cell_angle_gamma 120.00
_symmetry_space_group_name_H-M 'P 32 2 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0.4697 0.0000 0.3333 1.0
"""


class VolumeSummary(NamedTuple):
    """The structural facts [`check_ccp4`] compares against `chem`'s own
    read of a CCP4/MRC map (#340) -- grid dimensions and cell exactly,
    density values within a small `float32` tolerance, since that is the
    map's own storage precision and not a claim either reader makes about
    exactness.
    """

    dims: tuple[int, int, int]
    cell: tuple[float, float, float, float, float, float]
    values: tuple[float, ...]


def summarize_ccp4(path) -> Optional[VolumeSummary]:
    """Reads a CCP4/MRC map and summarises what gemmi saw, or `None` if
    gemmi could not read it. A different entry point again
    (`gemmi.read_ccp4_map`, not `read_structure_string`/
    `make_small_structure_from_block`) -- CCP4 is a grid, not atoms, so
    there is no structure to build here at all.
    """
    try:
        m = _gemmi.read_ccp4_map(str(path))
        m.setup(0.0)
    except Exception:
        return None
    grid = m.grid
    if grid.nu * grid.nv * grid.nw == 0:
        return None
    cell = grid.unit_cell
    values = tuple(round(float(v), 4) for v in _np.asarray(grid.array).flatten())
    return VolumeSummary(
        (grid.nu, grid.nv, grid.nw),
        (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma),
        values,
    )


#: A tiny (2x2x2) but complete map -- non-cubic cell, non-constant values,
#: so the sanity check cannot pass by accident the way an all-zero grid
#: could.
def write_ccp4_fixture(path) -> None:
    grid = _gemmi.FloatGrid(2, 2, 2)
    grid.set_unit_cell(_gemmi.UnitCell(4, 5, 6, 90, 90, 90))
    grid.spacegroup = _gemmi.SpaceGroup("P1")
    grid.array[:] = _np.arange(8, dtype=_np.float32).reshape(2, 2, 2)
    m = _gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header()
    m.write_ccp4_map(str(path))


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
    # The small-molecule path (#320) is a different gemmi entry point
    # entirely (`make_small_structure_from_block`) -- gated the same way the
    # PDB path is, so a `check_cif_core` failure is unambiguously about chem.
    if summarize_small_molecule(SANITY_CIF_CORE) is None:
        raise RuntimeError(
            "gemmi cannot read a trivial small-molecule CIF fixture, so it "
            "is broken rather than strict — refusing to report its answers "
            "as findings"
        )
    # The CCP4 path (#340) is a third gemmi entry point (`read_ccp4_map`) --
    # gated the same way, over a real temp file since `read_ccp4_map` (unlike
    # `read_structure_string`) has no in-memory-text overload.
    import tempfile as _tempfile
    from pathlib import Path as _Path

    with _tempfile.TemporaryDirectory() as _tmp:
        _path = _Path(_tmp) / "sanity.ccp4"
        write_ccp4_fixture(_path)
        if summarize_ccp4(_path) is None:
            raise RuntimeError(
                "gemmi cannot read a trivial CCP4 map it just wrote itself, "
                "so it is broken rather than strict — refusing to report "
                "its answers as findings"
            )
    return summarize
