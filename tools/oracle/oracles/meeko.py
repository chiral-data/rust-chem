"""Meeko as the third PDBQT dialect (#258).

`io/pdbqt.rs` targets AutoDock's own documented spec rather than any
toolkit's quirks, precisely because obabel and Meeko disagree with each
other — #173 records that downstream code parses both. That stance makes
"read both" a promise, and until this module nothing held it to anything:
Meeko was not installed, so only OpenBabel's dialect was ever exercised.

Like `gemmi`, this does not implement the `Oracle` dataclass and is not in
`load()`'s list — it answers "what does a PDBQT written by the reference
ligand-preparation tool look like", not a SMILES-shaped question.

The two dialects differ in more than whitespace, which is the point:

    meeko   ATOM      1  C   UNL     1  ...  1.00  0.00    +0.034 C
    obabel  ATOM      1  C   UNL     1  ...  0.00  0.00    +0.000 C
    chem    ATOM      1 C    LIG     1  ...  1.00  0.00     0.000 C

Meeko computes real Gasteiger charges and keeps polar hydrogens as `HD`;
obabel writes zeros for both charge and occupancy; chem writes an unsigned
charge, `LIG` rather than `UNL`, and starts the atom name one column
earlier. All three are read back correctly — the PDBQT type column at 78-79
is authoritative, so the name column's offset changes no meaning.
"""

from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from meeko import MoleculePreparation, PDBQTWriterLegacy

# Importing meeko turns RDKit's logger back on, undoing what `oracles/rdkit.py`
# did for the reason stated there: the rejection corpus makes parse-failure
# logging deliberate and voluminous. Silence it again, after the import rather
# than before, or meeko simply re-enables it.
RDLogger.DisableLog("rdApp.*")

#: Enough to exercise a torsion tree: ethanol's C-O bond is rotatable, so
#: Meeko emits ROOT/BRANCH/ENDBRANCH rather than a flat single-fragment file.
SANITY = "CCO"


def pdbqt_of_smiles(smiles: str) -> Optional[str]:
    """Meeko's own PDBQT for a molecule, or `None` if it declined to write one.

    A fixed embedding seed, so a failing run is reproducible: the coordinates
    are arbitrary but they must not differ between two runs of the same check.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xF00D) != 0:
        return None

    # Meeko refuses a disconnected molecule outright ("Must have 1" fragment),
    # which is a ligand-preparation tool declining something that is not a
    # ligand rather than a defect. A caller treats that like any other "this
    # writer had no opinion" and moves on.
    try:
        setups = MoleculePreparation()(mol)
    except ValueError:
        return None
    if not setups:
        return None
    # `write_string` returns (text, ok, error_message) -- the bool is not
    # redundant with an empty string, so it is what gets checked.
    text, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    return text if ok else None


def load_meeko():
    """Confirms Meeko can write a PDBQT before its output is trusted.

    The same discipline as `oracles.load()` and `load_gemmi`: a toolkit that
    is broken rather than strict must fail loudly here instead of quietly
    reporting that our reader disagrees with nothing.
    """
    text = pdbqt_of_smiles(SANITY)
    if text is None or "ROOT" not in text:
        raise RuntimeError(
            "meeko cannot write a PDBQT for ethanol, so it is broken rather "
            "than strict — refusing to report its answers as findings"
        )
    return pdbqt_of_smiles
