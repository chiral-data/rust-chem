"""The oracle registry.

An oracle answers three questions about a molecule and nothing else: can you
read this, what is its identity, and what bits does its fingerprint set. Keeping
the surface that narrow is what lets a third one be added without touching any
check.

gemmi (#224) answers a structurally different question — atoms/cell/chains/
residues as deposited, for mmCIF and PDB, which `chem` reads starting with
#223/#224 — and does not fit `Oracle`'s SMILES-shaped contract at all: there
is no bare-SMILES `parses`/`identity` for a toolkit that only reads
structures, so it is not part of `load()`'s list below. See `gemmi.py` and
its own `load_gemmi()`, sanity-gated the same way, against a trivial mmCIF
fixture instead of `SANITY` (`"CCO"` means nothing to gemmi).
"""

from dataclasses import dataclass
from typing import Callable, Optional


@dataclass
class Oracle:
    name: str
    #: Whether this toolkit accepts the SMILES at all.
    parses: Callable[[str], bool]
    #: A canonical identity, or None if the input could not be read.
    #:
    #: InChI where the toolkit can produce it — canonical independently of
    #: anyone's SMILES writer, and layered, so a mismatch says *what* moved
    #: rather than merely that something did.
    identity: Callable[[str], Optional[str]]
    #: Identity computed from a block of SDF text rather than SMILES.
    identity_of_sdf: Callable[[str], Optional[str]]
    #: Set bits of a Morgan fingerprint, or None if unsupported.
    fingerprint: Optional[Callable[[str, int, int], Optional[list[int]]]] = None


#: Something every cheminformatics toolkit can read. If one cannot, it is
#: broken rather than opinionated.
SANITY = "CCO"


def load() -> list[Oracle]:
    """Every oracle available in this image, in a stable order.

    An import failure is fatal rather than skipped: a harness that quietly runs
    with one fewer oracle reports green for a reason nobody asked for.

    So is an oracle that imports and then cannot read ethanol. OpenBabel does
    exactly that when its plugin loader hits a missing system library — every
    format silently fails to register, `SetInFormat` returns false, and the
    toolkit reports that no molecule on earth parses. It answers every question
    confidently and wrongly, which produced 31 false mismatches the first time
    this ran. Checking costs one parse and turns a silent lie into a loud stop.
    """
    from . import openbabel, rdkit

    oracles = [rdkit.oracle(), openbabel.oracle()]
    for oracle in oracles:
        if not oracle.parses(SANITY) or oracle.identity(SANITY) is None:
            raise RuntimeError(
                f"oracle {oracle.name!r} cannot read {SANITY!r}, so it is broken "
                f"rather than strict — refusing to report its answers as findings"
            )
    return oracles
