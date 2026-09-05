"""Drives the `chem` binary as a subprocess.

The interface an oracle harness should use, and the one thing worth keeping
from the `compare.py` this replaces: no bindings, no linkage, nothing that
could accidentally make a toolkit a dependency of the crate. What is tested is
the shipped binary, exactly as a user runs it.
"""

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def binary() -> Path:
    """Where the release binary lives.

    The image sets `CHEM_BIN` and builds it at image-build time. `compare.py`
    used to shell out to `cargo build` on every run, which meant a container
    that had to reach the network to do its job.
    """
    return Path(os.environ.get("CHEM_BIN", REPO_ROOT / "target/release/chem"))


@dataclass
class Run:
    stdout: str
    stderr: str
    code: int


def run(args: list[str], stdin: str = "") -> Run:
    result = subprocess.run(
        [str(binary()), *args],
        input=stdin,
        capture_output=True,
        text=True,
    )
    return Run(result.stdout, result.stderr, result.returncode)


@dataclass
class Record:
    name: str
    smiles: str


def read_corpus(path: Path) -> list[Record]:
    """Parses a corpus file the way `chem` does: `SMILES<space>name`.

    Deliberately not via the binary. The harness has to be able to say "chem
    disagrees about this line", which means knowing what the line said
    independently of chem's opinion of it.
    """
    records = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"{path.stem}:{smiles}"
        records.append(Record(name=name, smiles=smiles))
    return records


def is_known_gap(name: str) -> bool:
    """Whether this entry is a pinned gap.

    Marked by the entry's *name*, not a nearby comment. A comment-scanning
    version armed itself on the phrase appearing in a file's own header and
    swept up every line below it; a name carries the marker into every listing
    and cannot mis-arm.

    Gaps run both ways: input that should be rejected and is not, and valid
    input the parser cannot read yet. Both are pinned rather than deleted so
    the corpus describes the parser as it behaves.
    """
    return name.startswith("known-gap-")


def parse(smiles: str) -> bool:
    """Whether `chem` accepts this SMILES.

    `chem info` reports a skipped record on stderr and reads zero molecules,
    so acceptance is "did a row come back", not the exit code — a file of
    mixed good and bad lines exits 0 while skipping some.
    """
    result = run(["info", "-"], stdin=f"{smiles} probe\n")
    rows = [
        line
        for line in result.stdout.splitlines()
        if line and not line.startswith("name\t")
    ]
    return len(rows) > 0


def write_smiles(smiles: str) -> str | None:
    """Round-trips through chem's SMILES writer, returning what it emitted."""
    result = run(["aromatic", "-", "--out-format", "smiles"], stdin=f"{smiles} probe\n")
    for line in result.stdout.splitlines():
        if line.strip():
            return line.split(None, 1)[0]
    return None


def write_sdf(smiles: str) -> str | None:
    """SMILES in, SDF out, via the command that computes coordinates."""
    result = run(["coords", "-", "--out-format", "sdf"], stdin=f"{smiles} probe\n")
    return result.stdout if result.stdout.strip() else None


def fingerprint(smiles: str, radius: int, nbits: int) -> list[int] | None:
    """The set bits of a Morgan fingerprint.

    `chem fp` writes `name<TAB>hex`, least-significant-bit-first within each
    nibble, so bit *n* is character *n // 4* at position *n % 4*.

    Piped through `chem aromatic` first, matching the documented
    `chem aromatic x.smi | chem fp` pipeline: `chem fp` deliberately does not
    perceive aromaticity on its own (#192), and comparing against RDKit's
    fingerprint — which perceives aromaticity as part of RDKit's own
    sanitization — only makes sense against the same pipeline a real user
    would run for that comparison.
    """
    perceived = run(["aromatic", "-", "--format", "smiles"], stdin=f"{smiles} probe\n")
    if perceived.code != 0:
        return None
    rows = [
        line
        for line in perceived.stdout.splitlines()
        if line and not line.startswith("#")
    ]
    if not rows:
        return None
    perceived_smiles = rows[0].split()[0]

    result = run(
        ["fp", "-", "--radius", str(radius), "--size", str(nbits), "-b", "cpu"],
        stdin=f"{perceived_smiles} probe\n",
    )
    if result.code != 0:
        return None
    rows = [
        line
        for line in result.stdout.splitlines()
        if line and not line.startswith("#") and not line.startswith("name\t")
    ]
    if not rows:
        return None

    digits = rows[0].split("\t")[1].strip()
    bits = []
    for index, char in enumerate(digits):
        nibble = int(char, 16)
        for offset in range(4):
            if nibble & (1 << offset):
                bits.append(index * 4 + offset)
    return sorted(bit for bit in bits if bit < nbits)
