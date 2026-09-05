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


def oracle() -> Oracle:
    return Oracle(
        name="openbabel",
        parses=parses,
        identity=identity,
        identity_of_sdf=identity_of_sdf,
        fingerprint=None,
    )
