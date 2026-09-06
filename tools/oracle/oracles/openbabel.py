"""OpenBabel as an oracle.

Runs as a separate process in a dev image and is never linked, vendored, or
named in Cargo.toml — which is what keeps a GPL-2 toolkit judging an MIT crate
a matter of mere aggregation.

It has no InChI here, so identity comes from its own canonical SMILES. That is
a weaker key than RDKit's InChI: it is canonical within OpenBabel but says
nothing about layers, so a mismatch it reports localises less precisely. It
earns its place by disagreeing with RDKit — aromaticity perception especially —
and those disagreements are findings rather than failures.
"""

from typing import Optional

from openbabel import openbabel as ob

from . import Oracle

# OpenBabel narrates every parse failure to stderr. The rejection corpus makes
# that deliberate and noisy.
ob.obErrorLog.SetOutputLevel(0)


def _read(text: str, fmt: str):
    conversion = ob.OBConversion()
    conversion.SetInFormat(fmt)
    mol = ob.OBMol()
    if not conversion.ReadString(mol, text):
        return None
    return mol if mol.NumAtoms() > 0 else None


def parses(smiles: str) -> bool:
    return _read(smiles, "smi") is not None


def _canonical(mol) -> str:
    conversion = ob.OBConversion()
    conversion.SetOutFormat("can")
    return conversion.WriteString(mol).split()[0]


def identity(smiles: str) -> Optional[str]:
    mol = _read(smiles, "smi")
    return _canonical(mol) if mol is not None else None


def identity_of_sdf(text: str) -> Optional[str]:
    mol = _read(text, "sdf")
    return _canonical(mol) if mol is not None else None


def round_trip_pdb(text: str) -> Optional[str]:
    """PDB in, OpenBabel's own PDB out.

    Exists so `check_pdb` can record what this toolkit does to the B-factor
    column rather than assume it -- #173 names the behaviour, and #258 measured
    it: the column is zeroed while the occupancy beside it survives.
    """
    mol = _read(text, "pdb")
    if mol is None:
        return None
    conversion = ob.OBConversion()
    conversion.SetInAndOutFormats("pdb", "pdb")
    return conversion.WriteString(mol)


def pdbqt_of_smiles(smiles: str) -> Optional[str]:
    """OpenBabel's own PDBQT for a molecule.

    The dialect chem is already known to read; paired with Meeko's in
    `check_pdbqt` so "read both" means both rather than whichever one was
    installed.
    """
    mol = _read(smiles, "smi")
    if mol is None:
        return None
    conversion = ob.OBConversion()
    conversion.SetInAndOutFormats("smi", "pdbqt")
    return conversion.WriteString(mol)


def oracle() -> Oracle:
    return Oracle(
        name="openbabel",
        parses=parses,
        identity=identity,
        identity_of_sdf=identity_of_sdf,
        fingerprint=None,
        round_trip_pdb=round_trip_pdb,
        pdbqt_of_smiles=pdbqt_of_smiles,
    )
