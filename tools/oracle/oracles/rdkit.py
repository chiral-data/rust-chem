"""RDKit as an oracle.

The primary one: it is the only toolkit here that produces InChI, which is what
every identity comparison rests on. See `oracles/__init__.py`.
"""

from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

from . import Oracle

# RDKit logs parse failures to stderr at load. The rejection corpus makes that
# deliberate and voluminous, and a harness whose own expected cases look like
# errors is unreadable.
RDLogger.DisableLog("rdApp.*")


def _mol(smiles: str):
    return Chem.MolFromSmiles(smiles)


def parses(smiles: str) -> bool:
    return _mol(smiles) is not None


def identity(smiles: str) -> Optional[str]:
    mol = _mol(smiles)
    return Chem.MolToInchi(mol) if mol is not None else None


def identity_of_sdf(text: str) -> Optional[str]:
    supplier = Chem.SDMolSupplier()
    supplier.SetData(text)
    for mol in supplier:
        if mol is not None:
            return Chem.MolToInchi(mol)
    return None


def fingerprint(smiles: str, radius: int, nbits: int) -> Optional[list[int]]:
    """Morgan bits.

    `radius` and `nbits` are honoured here. The harness this replaces took both
    and then hardcoded `radius=2, fpSize=2048` in the generator, so the
    parameters were a latent bug the moment anyone passed anything else.
    """
    mol = _mol(smiles)
    if mol is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return sorted(generator.GetFingerprint(mol).GetOnBits())


def identity_of_commonchem(text: str) -> Optional[str]:
    """InChI of the first molecule in a commonchem/rdkitjson document.

    RDKit defines this format, so unlike every other check here the comparison
    is against the reference implementation rather than against a toolkit that
    merely happens to support the format.
    """
    try:
        mols = Chem.JSONToMols(text)
    except Exception:
        return None
    for mol in mols:
        if mol is not None:
            return Chem.MolToInchi(mol)
    return None


def bond_stereo_of_commonchem(text: str) -> Optional[list]:
    """`(begin, end, stereo, stereo_atoms)` for every double bond, in order.

    Read off the bonds rather than through `MolToSmiles`, which ignores the
    stereo atoms — see the field's comment in `oracles/__init__.py`.
    """
    try:
        mols = Chem.JSONToMols(text)
    except Exception:
        return None
    if not mols or mols[0] is None:
        return None
    out = []
    for bond in mols[0].GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            out.append(
                (
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    str(bond.GetStereo()),
                    sorted(bond.GetStereoAtoms()),
                )
            )
    return out


def commonchem_of_smiles(smiles: str) -> Optional[str]:
    """RDKit's own commonchem document for a SMILES — the read direction's input."""
    mol = _mol(smiles)
    return Chem.MolToJSON(mol) if mol is not None else None


def formula(smiles: str) -> Optional[str]:
    mol = _mol(smiles)
    return rdMolDescriptors.CalcMolFormula(mol) if mol is not None else None


def formula_of_commonchem(text: str) -> Optional[str]:
    """Molecular formula of the first molecule in a commonchem document.

    The hydrogen-count check InChI cannot perform — see the field comment in
    `oracles/__init__.py`.
    """
    try:
        mols = Chem.JSONToMols(text)
    except Exception:
        return None
    for mol in mols:
        if mol is not None:
            return rdMolDescriptors.CalcMolFormula(mol)
    return None


def formula_of_pdb(text: str) -> Optional[str]:
    """Molecular formula of a PDB, as RDKit reads it.

    The hydrogen question, asked of a format that states no count of its own. A
    reader has to imply one from the bonds, and PDB's own columns cannot show
    whether it did -- an implicit hydrogen creates no atom and fills no field,
    so `check_pdb`'s structural and per-atom comparisons are blind to it. That
    is how #285 stayed invisible here for a whole milestone.

    Only meaningful where `CONECT` covers every atom. RDKit infers the missing
    bonds from geometry and this crate deliberately does not, so on a partly
    connected file the formulae differ by bond perception rather than by
    hydrogen counting -- `check_pdb` notes those files rather than comparing
    them.

    `sanitize=False` because a partial CONECT block leaves valences RDKit would
    reject outright, and a formula needs none of that analysis.
    """
    try:
        mol = Chem.MolFromPDBBlock(text, removeHs=False, sanitize=False)
    except Exception:
        return None
    if mol is None:
        return None
    try:
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return None


def oracle() -> Oracle:
    return Oracle(
        name="rdkit",
        parses=parses,
        identity=identity,
        identity_of_sdf=identity_of_sdf,
        fingerprint=fingerprint,
        identity_of_commonchem=identity_of_commonchem,
        bond_stereo_of_commonchem=bond_stereo_of_commonchem,
        commonchem_of_smiles=commonchem_of_smiles,
        formula=formula,
        formula_of_commonchem=formula_of_commonchem,
        formula_of_pdb=formula_of_pdb,
    )
