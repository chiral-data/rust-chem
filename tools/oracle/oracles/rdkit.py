"""RDKit as an oracle.

The primary one: it is the only toolkit here that produces InChI, which is what
every identity comparison rests on. See `oracles/__init__.py`.
"""

from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

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


def oracle() -> Oracle:
    return Oracle(
        name="rdkit",
        parses=parses,
        identity=identity,
        identity_of_sdf=identity_of_sdf,
        fingerprint=fingerprint,
    )
