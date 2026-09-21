#!/usr/bin/env python3
"""Differential harness: RDKit and OpenBabel as oracles for `chem`.

    docker build -t chem-oracle -f tools/oracle/Dockerfile .
    docker run --rm chem-oracle
    docker run --rm chem-oracle python3 tools/oracle/run.py --check write --verbose

Nothing here is a workspace member, a dev-dependency, or named in Cargo.toml.
It drives the release binary as a subprocess, so a developer with neither
Python nor a cheminformatics toolkit still runs the whole `cargo` gate.

# What a finding is, and is not

A **mismatch** is `chem` disagreeing with an oracle: a failure, and the point.

A **disagreement** is the oracles differing from *each other* on the same input
— aromaticity perception, most often. That is a fact about the toolkits, not
about us, and the milestone wants it recorded in the fidelity table rather than
resolved by picking a favourite. Reported, never fatal.

An **expected loss** is a mismatch the crate already knows about, because a
format's `Carries` mask says the attribute does not survive. Reported as
context, never fatal: a drop the mask predicted is a pass. That is the only
parity rule that works for lossy formats — demanding byte-identical round trips
produces a permanently red suite that tells you nothing.
"""

import argparse
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import chem
from oracles import Oracle, load

CORPUS = chem.REPO_ROOT / "crates/chem/tests/corpus"


@dataclass
class Report:
    passed: int = 0
    mismatches: list[str] = field(default_factory=list)
    disagreements: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def ok(self):
        self.passed += 1

    def mismatch(self, line: str):
        self.mismatches.append(line)

    def disagree(self, line: str):
        self.disagreements.append(line)

    def note(self, line: str):
        self.notes.append(line)


def inchi_layers(a: str, b: str) -> str:
    """Which InChI layers differ, by name.

    The reason identity is InChI rather than a canonical SMILES string. A
    harness that reports "differs" and stops costs a debugging session per
    finding; one that says the stereo layer moved has already done the work.
    """
    names = {
        "c": "connectivity",
        "h": "hydrogens",
        "q": "charge",
        "p": "protonation",
        "b": "double-bond stereo",
        "t": "tetrahedral stereo",
        "m": "stereo parity",
        "s": "stereo type",
        "i": "isotope",
    }

    def layers(text: str) -> dict[str, str]:
        out = {}
        for part in text.split("/")[1:]:
            if part and part[0] in names:
                out[part[0]] = part[1:]
            else:
                out.setdefault("formula", part)
        return out

    left, right = layers(a), layers(b)
    moved = []
    for key in sorted(set(left) | set(right)):
        if left.get(key) != right.get(key):
            label = names.get(key, key)
            moved.append(f"{label} ({left.get(key, '-')} -> {right.get(key, '-')})")
    return "; ".join(moved) if moved else "identical layers, different strings"


def difference(oracle: Oracle, before: str, after: str) -> str:
    """How two identities differ, as precisely as this oracle allows."""
    if before.startswith("InChI=") and after.startswith("InChI="):
        return inchi_layers(before, after)
    # OpenBabel's key is its own canonical SMILES: canonical within OpenBabel,
    # but unlayered, so this is the best it can say.
    return f"{before} -> {after}"


def check_parse(oracles: list[Oracle], verbose: bool) -> Report:
    """Does `chem` accept what the oracles accept, and refuse what they refuse?

    The rejection corpus is where this earns its keep. A suite of valid
    molecules proves what a parser accepts and never what it refuses.
    """
    report = Report()
    for path in sorted(CORPUS.glob("*.smi")):
        rejecting = path.stem == "invalid"

        for record in chem.read_corpus(path):
            ours = chem.parse(record.smiles)
            theirs = {o.name: o.parses(record.smiles) for o in oracles}

            # A pinned gap is a known deviation, reported as context and never
            # fatal — in either direction. In `invalid.smi` it is input chem
            # wrongly accepts (#190); elsewhere it is valid input chem cannot
            # read yet (#191).
            if chem.is_known_gap(record.name):
                report.note(
                    f"{path.name}: {record.name} — chem "
                    f"{'accepted' if ours else 'rejected'}, pinned gap; "
                    f"oracles {theirs}"
                )
                continue

            if rejecting:
                if ours:
                    report.mismatch(
                        f"{path.name}: chem ACCEPTED {record.smiles!r} "
                        f"({record.name}), which must be rejected; oracles {theirs}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    reject ok  {record.name:<34} {record.smiles}")
                continue

            for name, accepted in theirs.items():
                if accepted != ours:
                    report.mismatch(
                        f"{path.name}: chem {'accepted' if ours else 'rejected'} "
                        f"{record.smiles!r} ({record.name}), {name} did the opposite"
                    )
                    break
            else:
                report.ok()
                if verbose:
                    print(f"    ok         {record.name:<34} {record.smiles}")

            if len(set(theirs.values())) > 1:
                report.disagree(
                    f"{path.name}: {record.name} — oracles split on whether "
                    f"{record.smiles!r} parses: {theirs}"
                )
    return report


def check_write(oracles: list[Oracle], verbose: bool) -> Report:
    """Does a round trip through chem's SMILES writer preserve the molecule?"""
    report = Report()
    for path in sorted(CORPUS.glob("*.smi")):
        if path.stem == "invalid":
            continue
        for record in chem.read_corpus(path):
            if chem.is_known_gap(record.name):
                continue  # chem cannot read it yet; the parse check reports it
            written = chem.write_smiles(record.smiles)
            if written is None:
                report.mismatch(f"{path.name}: chem wrote nothing for {record.name}")
                continue

            identities = {}
            for oracle in oracles:
                before = oracle.identity(record.smiles)
                after = oracle.identity(written)
                if before is None:
                    continue
                identities[oracle.name] = (before, after)

                if after is None:
                    report.mismatch(
                        f"{path.name}: {record.name} — chem wrote {written!r}, "
                        f"which {oracle.name} cannot read back"
                    )
                elif before != after:
                    report.mismatch(
                        f"{path.name}: {record.name} — {oracle.name} says "
                        f"{difference(oracle, before, after)}"
                        f"  [{record.smiles} -> {written}]"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok         {record.name:<34} {record.smiles}")

            verdicts = {
                name: before == after for name, (before, after) in identities.items()
            }
            if len(set(verdicts.values())) > 1:
                report.disagree(
                    f"{path.name}: {record.name} — oracles split on the same "
                    f"round trip: {verdicts}"
                )
    return report


def check_sdf(oracles: list[Oracle], verbose: bool) -> Report:
    """SMILES in, SDF out, oracle reads it back — is it the same molecule?

    The path that caught three nitro compounds losing their formal charges in
    #183. Losses a `Carries` mask predicted are expected and reported as such.
    """
    report = Report()
    for path in sorted(CORPUS.glob("*.smi")):
        if path.stem == "invalid":
            continue
        for record in chem.read_corpus(path):
            if chem.is_known_gap(record.name):
                continue  # chem cannot read it yet; the parse check reports it
            sdf = chem.write_sdf(record.smiles)
            if sdf is None:
                report.mismatch(f"{path.name}: chem wrote no SDF for {record.name}")
                continue

            for oracle in oracles:
                before = oracle.identity(record.smiles)
                after = oracle.identity_of_sdf(sdf)
                if before is None:
                    continue
                if after is None:
                    report.mismatch(
                        f"{path.name}: {record.name} — {oracle.name} cannot read "
                        f"the SDF chem wrote"
                    )
                elif before != after:
                    report.mismatch(
                        f"{path.name}: {record.name} — SDF round trip, "
                        f"{oracle.name} says {difference(oracle, before, after)}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok         {record.name:<34} {record.smiles}")
    return report


def check_json(oracles: list[Oracle], verbose: bool) -> Report:
    """commonchem JSON, both directions, against the toolkit that defines it (#229).

    The only format in scope where the oracle is the reference implementation
    rather than a second opinion, so both directions are worth checking:

    - **write**: chem emits a document, RDKit reads it, and the InChI must
      match RDKit's InChI for the original SMILES. This is what proves the
      `rdkitRepresentation` extension is emitted correctly — aromaticity lives
      only there, so a document missing it reads back kekulised and non-
      aromatic, with a different InChI.
    - **read**: RDKit emits a document (in its own `rdkitjson` v12 dialect,
      which it always uses and we never write), chem reads it and writes its
      canonical SMILES, and RDKit's InChI of *that* must match. This is what
      proves the dialect is accepted and the extension is applied on the way
      in.

    Only RDKit takes part. OpenBabel has no reader for this format, and gemmi
    is structural — an oracle without `identity_of_commonchem` is skipped, the
    same way `check_fp` skips one without `fingerprint`.
    """
    report = Report()
    for path in sorted(CORPUS.glob("*.smi")):
        if path.stem == "invalid":
            continue
        for record in chem.read_corpus(path):
            if chem.is_known_gap(record.name):
                continue  # chem cannot read it yet; the parse check reports it

            for oracle in oracles:
                if oracle.identity_of_commonchem is None:
                    continue
                before = oracle.identity(record.smiles)
                if before is None:
                    continue

                # --- write direction -------------------------------------
                written = chem.write_commonchem(record.smiles)
                if written is None:
                    report.mismatch(
                        f"{path.name}: chem wrote no commonchem for {record.name}"
                    )
                    continue
                after = oracle.identity_of_commonchem(written)
                if after is None:
                    report.mismatch(
                        f"{path.name}: {record.name} — {oracle.name} cannot read "
                        f"the commonchem chem wrote"
                    )
                elif before != after:
                    report.mismatch(
                        f"{path.name}: {record.name} — commonchem write, "
                        f"{oracle.name} says {difference(oracle, before, after)}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok  write  {record.name:<34} {record.smiles}")

                # Hydrogen counts, which the InChI comparison above cannot
                # see: its `/p` layer normalises mobile protons away, so a
                # carboxylate carrying two impossible hydrogens has the same
                # InChI as a correct one. commonchem is the first format here
                # that states a per-atom hydrogen count, so this is the only
                # check that looks at what it actually wrote.
                if oracle.formula is not None and oracle.formula_of_commonchem is not None:
                    want = oracle.formula(record.smiles)
                    got = oracle.formula_of_commonchem(written)
                    if want is not None and got is not None:
                        if want != got:
                            report.mismatch(
                                f"{path.name}: {record.name} — commonchem write, "
                                f"formula {want} -> {got}"
                            )
                        else:
                            report.ok()

                # --- read direction --------------------------------------
                if oracle.commonchem_of_smiles is None:
                    continue
                theirs = oracle.commonchem_of_smiles(record.smiles)
                if theirs is None:
                    continue
                ours = chem.read_commonchem(theirs)
                if ours is None:
                    report.mismatch(
                        f"{path.name}: {record.name} — chem cannot read the "
                        f"commonchem {oracle.name} wrote"
                    )
                    continue
                back = oracle.identity(ours)
                if back is None:
                    report.mismatch(
                        f"{path.name}: {record.name} — {oracle.name} cannot read "
                        f"the SMILES chem wrote from their commonchem"
                    )
                elif before != back:
                    report.mismatch(
                        f"{path.name}: {record.name} — commonchem read, "
                        f"{oracle.name} says {difference(oracle, before, back)}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok  read   {record.name:<34} {record.smiles}")

                # --- bond stereo, read off the bonds ---------------------
                # Not through any SMILES writer: RDKit's ignores the
                # `stereoAtoms` its JSON reader preserves, so two chemically
                # opposite documents render identically. See the field comment
                # in `oracles/__init__.py`.
                if oracle.bond_stereo_of_commonchem is not None:
                    mine = oracle.bond_stereo_of_commonchem(written)
                    theirs_stereo = oracle.bond_stereo_of_commonchem(theirs)
                    if mine is not None and theirs_stereo is not None:
                        # Only bonds that actually assert a configuration.
                        # Counting `STEREONONE` would compare how many double
                        # bonds each side *has*, which differs legitimately
                        # whenever aromaticity perception does -- `chem` does
                        # not perceive outside `chem aromatic` (#192), so it
                        # writes `C1=CC=CC=C1` with three plain double bonds
                        # where RDKit writes an aromatic ring with none. That
                        # is a different question, and the `write`/`read`
                        # comparisons above already answer it.
                        mine_kinds = sorted(
                            s for _, _, s, _ in mine if s != "STEREONONE"
                        )
                        their_kinds = sorted(
                            s for _, _, s, _ in theirs_stereo if s != "STEREONONE"
                        )
                        if mine_kinds != their_kinds:
                            report.mismatch(
                                f"{path.name}: {record.name} — bond stereo, "
                                f"chem {mine_kinds} vs {oracle.name} {their_kinds}"
                            )
                        else:
                            report.ok()
    return report


def tanimoto(a: set[int], b: set[int]) -> float:
    union = len(a | b)
    return len(a & b) / union if union else 1.0


#: Fold width for counting atom environments, which is deliberately not the
#: width being tested. Collisions vanish across the whole corpus by 8192; this
#: leaves an order of magnitude of headroom as it grows, and costs 0.2s against
#: 8192's 0.1s. A megabit costs 1.7s and buys no further agreement, because the
#: hex payload is 256KB per molecule.
ENVIRONMENT_WIDTH = 65536


def check_fp(oracles: list[Oracle], verbose: bool, radius: int, nbits: int) -> Report:
    """Do our fingerprints rank molecules the way RDKit's do?

    Not a bit-for-bit comparison, which cannot pass and would be noise. Morgan
    bit *positions* are the output of a particular hash of a particular
    invariant packing: for methane we set one bit and so does RDKit, ours at
    168 and theirs at 1264. Matching those would mean cloning RDKit's hash,
    which is a decision nobody has taken (#192).

    What is comparable across implementations — and what a fingerprint is
    actually *for* — is similarity ordering. So: with every corpus molecule as
    a query, does the nearest neighbour agree? That property holds for any
    chemically equivalent fingerprint regardless of hash, and breaks the moment
    ours stops describing the same environments.

    It holds only where the ranking is decided by more than one bit, which is
    why a third of this corpus is skipped with a note. Two molecules sharing a
    single bit out of 2048 have told you nothing: the bit is as likely to be a
    collision as a shared environment, so the "nearest" neighbour is whichever
    collision each implementation happened to get. Bit *count* is the other
    thing comparable across hashes, and it is checked here too — at
    `ENVIRONMENT_WIDTH` rather than the width under test, since only a fold
    wide enough to avoid collisions turns a bit count into an environment
    count. That is what would catch our enumeration drifting, which ranking
    structurally cannot (#253).
    """
    report = Report()
    molecules: list[tuple[str, str]] = []
    for path in sorted(CORPUS.glob("*.smi")):
        if path.stem == "invalid":
            continue
        for record in chem.read_corpus(path):
            if not chem.is_known_gap(record.name):
                molecules.append((record.name, record.smiles))

    for oracle in oracles:
        if oracle.fingerprint is None:
            continue

        ours, theirs = {}, {}
        for name, smiles in molecules:
            mine = chem.fingerprint(smiles, radius, nbits)
            yours = oracle.fingerprint(smiles, radius, nbits)
            if mine is None or yours is None:
                continue
            ours[name] = set(mine)
            theirs[name] = set(yours)

            # How many atom environments did each side enumerate?
            #
            # Counted at a width where a count means that, and not at the width
            # being tested. Folded into 2048 bits a smaller count means either
            # a collision or a missing environment, and nothing distinguishes
            # them: this comparison used to run at the requested width and
            # report "different environment counts, not just a different hash",
            # which was exactly backwards for the only molecule it ever fired
            # on. RDKit puts alanine's carboxyl carbon and its hydroxyl oxygen
            # -- a C and an O -- on bit 807, so it shows 12 bits over 13
            # environments while chem shows 13 over 13 (#253).
            #
            # Worth asserting rather than noting, because a bit count is the
            # one quantity comparable between two Morgan implementations
            # without cloning the hash, and the ranking comparison below
            # cannot see it: an implementation enumerating a different set of
            # environments could still rank neighbours identically.
            mine_envs = chem.fingerprint(smiles, radius, ENVIRONMENT_WIDTH)
            their_envs = oracle.fingerprint(smiles, radius, ENVIRONMENT_WIDTH)
            if mine_envs is not None and their_envs is not None:
                if len(mine_envs) != len(their_envs):
                    report.mismatch(
                        f"{name} — {len(mine_envs)} atom environments vs "
                        f"{oracle.name}'s {len(their_envs)}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok  envs   {name:<34} {len(mine_envs)}")

        names = sorted(ours)
        for name in names:
            others = [other for other in names if other != name]
            if not others:
                continue
            mine = max(others, key=lambda o: tanimoto(ours[name], ours[o]))
            yours = max(others, key=lambda o: tanimoto(theirs[name], theirs[o]))

            # A nearest neighbour reached by a single shared bit is not a
            # neighbour, it is a coincidence. In a 2048-bit space one bit is as
            # likely to be a hash collision as a shared environment, and
            # tracing every deciding bit in this corpus found four of five were
            # collisions: RDKit called hydroxide — one bit, an anionic oxygen —
            # the nearest thing to cyclopentadiene, while `chem` answered
            # cyclohexane on a bit carried by every saturated carbocycle here.
            # The check reported that as *us* diverging.
            #
            # Zero shared bits is the same failure at its extreme: every
            # candidate ties and `max` returns whichever it saw first, so the
            # comparison is of iteration orders. That case was skipped before
            # this rule generalised it, and it produced findings that moved
            # whenever the corpus grew.
            #
            # Skipped with a note rather than silently, so the gap in coverage
            # stays visible.
            mine_shared = len(ours[name] & ours[mine])
            their_shared = len(theirs[name] & theirs[yours])
            if mine_shared < 2 or their_shared < 2:
                report.note(
                    f"{name}: nearest neighbour rests on "
                    f"{min(mine_shared, their_shared)} shared bit(s), too few to "
                    f"tell chemistry from a hash collision — not compared"
                )
                continue

            if mine == yours:
                report.ok()
                if verbose:
                    print(f"    ok         {name:<34} nearest {mine}")
            else:
                report.mismatch(
                    f"{name} — nearest neighbour disagrees with {oracle.name}: "
                    f"ours {mine} (t={tanimoto(ours[name], ours[mine]):.3f}), "
                    f"theirs {yours} (t={tanimoto(theirs[name], theirs[yours]):.3f})"
                )
    return report


MMCIF_CORPUS = CORPUS / "mmcif"


def check_mmcif(oracles: list[Oracle], verbose: bool) -> Report:
    """Does `chem`'s mmCIF round trip agree with gemmi's independent read?

    Unlike every other check, this ignores `oracles` entirely (kept as a
    parameter only so it fits `main()`'s generic dispatch) and drives
    gemmi directly — neither RDKit nor OpenBabel can judge this format (see
    `oracles/gemmi.py`'s module doc). Comparison is structural — atom
    count, cell, chain ids, residue identities — computed by gemmi on both
    the original fixture and on what `chem convert --from mmcif --to
    mmcif` wrote back, rather than a text diff: this format has no
    canonical spelling to hold either side to.
    """
    from oracles import gemmi as gemmi_oracle

    summarize = gemmi_oracle.load_gemmi()
    report = Report()
    for path in sorted(MMCIF_CORPUS.glob("*.cif")):
        original = path.read_text()
        reference = summarize(original)
        if reference is None:
            report.mismatch(f"{path.name}: gemmi itself could not read this fixture")
            continue

        written = chem.convert_mmcif(original, "mmcif")
        if written is None:
            report.mismatch(f"{path.name}: chem could not round-trip this file")
            continue

        ours = summarize(written)
        if ours is None:
            report.mismatch(f"{path.name}: gemmi cannot read what chem wrote back")
        elif ours != reference:
            report.mismatch(f"{path.name}: chem's round trip disagrees with gemmi — {reference} vs {ours}")
        else:
            report.ok()
            if verbose:
                print(f"    ok         {path.name:<34} {reference.atom_count} atoms")
    return report


CIF_CORE_CORPUS = CORPUS / "cif_core"


def check_cif_core(oracles: list[Oracle], verbose: bool) -> Report:
    """Does `chem`'s CIF-core round trip agree with gemmi's independent
    read (#320)?

    Unlike BinaryCIF (#319), gemmi genuinely reads this dictionary
    (`gemmi.make_small_structure_from_block`, a different entry point than
    mmCIF/PDB's `read_structure_string`) -- confirmed by `load_gemmi`'s own
    sanity gate. Otherwise this mirrors `check_mmcif` exactly: structural
    comparison (cell, space group, fractional sites) between gemmi's read
    of the original fixture and gemmi's read of what `chem convert --from
    cif-core --to cif-core` wrote back.
    """
    from oracles import gemmi as gemmi_oracle

    gemmi_oracle.load_gemmi()
    report = Report()
    for path in sorted(CIF_CORE_CORPUS.glob("*.cif")):
        original = path.read_text()
        reference = gemmi_oracle.summarize_small_molecule(original)
        if reference is None:
            report.mismatch(f"{path.name}: gemmi itself could not read this fixture")
            continue

        written = chem.convert_cif_core(original, "cif-core")
        if written is None:
            report.mismatch(f"{path.name}: chem could not round-trip this file")
            continue

        ours = gemmi_oracle.summarize_small_molecule(written)
        if ours is None:
            report.mismatch(f"{path.name}: gemmi cannot read what chem wrote back")
        elif ours != reference:
            report.mismatch(
                f"{path.name}: chem's round trip disagrees with gemmi — {reference} vs {ours}"
            )
        else:
            report.ok()
            if verbose:
                print(f"    ok         {path.name:<34} {len(reference.sites)} sites")
    return report


PDB_CORPUS = CORPUS / "pdb"


def _pdb_states_every_bond(text: str) -> bool:
    """Whether every atom in this PDB appears in a `CONECT` record.

    The question decides whether a formula comparison means anything: RDKit
    infers the bonds a `CONECT` block leaves out, from geometry, and this crate
    deliberately does not (`io/pdb.rs`'s module doc). So on a partly connected
    file the two disagree about *bonds*, and comparing formulae there would
    measure bond perception while claiming to measure hydrogen counts.

    Serial numbers only -- deliberately not a PDB parser. A harness that parsed
    the format properly could be wrong in the same way the crate is, and a
    fixture's connectivity is a property of the text rather than of anyone's
    reading of it.
    """
    atoms, connected = set(), set()
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms.add(line[6:11].strip())
        elif line.startswith("CONECT"):
            connected.add(line[6:11].strip())
    return bool(atoms) and atoms <= connected


def check_pdb(oracles: list[Oracle], verbose: bool) -> Report:
    """Does `chem`'s PDB round trip keep what gemmi reads -- values included?

    Three questions, and the second is why this exists separately from
    `check_mmcif` rather than as more fixtures for it:

    1. The structural summary, as `check_mmcif` does it.
    2. **Per-atom occupancy and B-factor values.** `Carries` is presence-based,
       so a writer emitting a constant for every atom passes every mask
       assertion in the crate (#257 hands this over explicitly). Only a value
       comparison can see it.
    3. OpenBabel through the same path -- recorded as a `note`, never a
       mismatch. Its PDB writer zeroes the B-factor column while preserving the
       occupancy beside it, and being *different* from that is the correct
       behaviour (#173), so agreement here would be the bug.
    4. **The hydrogen count**, which none of the above can see. PDB states no
       count, so a reader implies one from the bonds -- and an implicit hydrogen
       creates no atom and fills no column, so questions 1 and 2 are identical
       whether the count was implied or left blank. #285 was invisible here for
       a whole milestone because of it (#293).

    Uses `oracles` unlike `check_mmcif`, which ignores it: OpenBabel is the
    subject of question 3 rather than a judge of questions 1 and 2.
    """
    from oracles import gemmi as gemmi_oracle

    summarize = gemmi_oracle.load_gemmi()
    report = Report()
    for path in sorted(PDB_CORPUS.glob("*.pdb")):
        original = path.read_text()
        reference = summarize(original, ".pdb")
        reference_sites = gemmi_oracle.sites(original)
        if reference is None or reference_sites is None:
            report.mismatch(f"{path.name}: gemmi itself could not read this fixture")
            continue

        written = chem.convert_pdb(original, "pdb")
        if written is None:
            report.mismatch(f"{path.name}: chem could not round-trip this file")
            continue

        ours = summarize(written, ".pdb")
        if ours is None:
            report.mismatch(f"{path.name}: gemmi cannot read what chem wrote back")
            continue
        if ours != reference:
            report.mismatch(
                f"{path.name}: chem's round trip disagrees with gemmi — {reference} vs {ours}"
            )
            continue

        our_sites = gemmi_oracle.sites(written)
        if our_sites != reference_sites:
            report.mismatch(
                f"{path.name}: per-atom values moved — "
                f"occupancy {reference_sites.occupancies} -> {our_sites.occupancies}, "
                f"b-factor {reference_sites.b_factors} -> {our_sites.b_factors}"
            )
            continue

        report.ok()
        if verbose:
            print(
                f"    ok         {path.name:<34} "
                f"{reference.atom_count} atoms, {len(set(reference_sites.b_factors))} distinct b"
            )

        # Question 4. The hydrogen count, via our own `pdb -> smi` output.
        #
        # Both formulae come from RDKit -- one read from the file, one computed
        # from the SMILES we wrote -- so what is compared is hydrogen counting
        # rather than canonical-string convention. OpenBabel judges nothing
        # here: it perceives aromaticity, so it would give a third answer for
        # the same file and disagreeing with it would mean nothing.
        rdkit_pdb_formula = next(
            (
                f
                for f in (getattr(o, "formula_of_pdb", None) for o in oracles)
                if f is not None
            ),
            None,
        )
        rdkit_smiles_formula = next(
            (f for f in (getattr(o, "formula", None) for o in oracles) if f is not None),
            None,
        )
        if rdkit_pdb_formula is not None and rdkit_smiles_formula is not None:
            if not _pdb_states_every_bond(original):
                # Not a mismatch and not silence: a reader of this output should
                # be able to tell a fixture that was skipped from one that was
                # never looked at.
                report.note(
                    f"{path.name}: hydrogen counts not compared -- CONECT does not "
                    "cover every atom, so the oracle infers bonds from geometry "
                    "where chem does not"
                )
            else:
                theirs = rdkit_pdb_formula(original)
                as_smiles = chem.convert_pdb(original, "smi")
                ours = (
                    rdkit_smiles_formula(as_smiles.split()[0])
                    if as_smiles and as_smiles.split()
                    else None
                )
                if theirs is None or ours is None:
                    report.mismatch(
                        f"{path.name}: could not compare hydrogen counts -- "
                        f"rdkit {theirs!r}, ours {ours!r}"
                    )
                elif theirs != ours:
                    report.mismatch(
                        f"{path.name}: hydrogen counts disagree with rdkit -- "
                        f"{theirs} vs our {ours}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok         {path.name:<34} formula {ours}")

        # Question 3. Not a mismatch: this is the oracle being wrong, recorded
        # so the divergence is visible rather than assumed.
        for oracle in oracles:
            theirs = getattr(oracle, "round_trip_pdb", None)
            if theirs is None:
                continue
            their_text = theirs(original)
            their_sites = gemmi_oracle.sites(their_text) if their_text else None
            if their_sites is None:
                report.note(f"{path.name}: {oracle.name} wrote a PDB gemmi cannot read")
            elif their_sites.b_factors != reference_sites.b_factors:
                kept = (
                    "occupancy preserved"
                    if their_sites.occupancies == reference_sites.occupancies
                    else f"occupancy also moved to {their_sites.occupancies}"
                )
                report.note(
                    f"{path.name}: {oracle.name} rewrote b-factors "
                    f"{reference_sites.b_factors} -> {their_sites.b_factors} ({kept}); "
                    "chem keeps them, which is the point"
                )
    return report


#: AutoDock types that name an element other than themselves. Spelled out
#: here rather than asked of `chem`, for the reason `read_corpus` exists: a
#: check has to know what a line said independently of chem's reading of it.
AUTODOCK_ELEMENT = {"A": "C", "OA": "O", "NA": "N", "SA": "S", "HD": "H", "HS": "H"}


def _pdbqt_atoms(text: str) -> list[tuple[str, float]]:
    """The element and partial charge of each atom of a PDBQT, in file order.

    Fixed columns, because that is what the format is: Meeko writes `+0.034`
    where chem writes ` 0.034`, so whitespace splitting would shift the fields
    apart on one dialect and not the other.
    """
    atoms = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        charge = line[70:76].strip()
        atom_type = line[77:79].strip()
        element = AUTODOCK_ELEMENT.get(atom_type, atom_type)
        atoms.append((element, round(float(charge), 4) if charge else 0.0))
    return atoms


def _mol2_atoms(text: str) -> list[tuple[str, float]]:
    """The same, read out of a Mol2 atom block."""
    atoms = []
    in_block = False
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            in_block = line.strip() == "@<TRIPOS>ATOM"
            continue
        if not in_block or not line.strip():
            continue
        fields = line.split()
        element = fields[5].split(".")[0]
        atoms.append((element, round(float(fields[8]), 4)))
    return atoms


def check_pdbqt(oracles: list[Oracle], verbose: bool) -> Report:
    """Does `chem` read every PDBQT dialect, not just its own?

    `io/pdbqt.rs` targets AutoDock's documented spec rather than obabel's or
    Meeko's quirks -- the two disagree with each other, and #173 records that
    downstream code parses both. That makes "read both" a promise, and until
    Meeko joined the image (#258) only one dialect was ever exercised.

    **Writes Mol2, not PDBQT.** A PDBQT round trip cannot measure the reader
    while #259 is open: these files carry no bonds, so chem reads N one-atom
    fragments and its writer keeps only the largest -- the first version of
    this check reported twelve findings that were all that one defect, with
    every atom read perfectly. Mol2 carries element and partial charge and
    drops no components, so what is compared is the read.

    Not a text comparison either: the dialects differ deliberately (`UNL` vs
    `LIG`, an explicit `+`, the atom-name column) and none of it changes what
    the file means.
    """
    from oracles import meeko as meeko_oracle

    write_pdbqt = meeko_oracle.load_meeko()
    report = Report()

    writers: list[tuple[str, Callable[[str], Optional[str]]]] = [("meeko", write_pdbqt)]
    for oracle in oracles:
        theirs = getattr(oracle, "pdbqt_of_smiles", None)
        if theirs is not None:
            writers.append((oracle.name, theirs))

    # `hard.smi` only: this is about dialects, not breadth, and every writer
    # here runs a 3D embedding per molecule.
    for record in chem.read_corpus(CORPUS / "hard.smi"):
        name = record.name
        if chem.is_known_gap(name):
            continue
        for writer_name, write in writers:
            theirs = write(record.smiles)
            if theirs is None:
                continue
            expected = _pdbqt_atoms(theirs)
            if not expected:
                report.note(f"{name}: {writer_name} wrote a PDBQT with no atoms")
                continue

            ours = chem.convert_pdbqt(theirs, "mol2")
            if ours is None:
                report.mismatch(f"{name}: chem could not read {writer_name}'s PDBQT")
                continue

            measured = _mol2_atoms(ours)
            if measured != expected:
                report.mismatch(
                    f"{name}: chem's read of {writer_name}'s PDBQT disagrees — "
                    f"{expected} vs {measured}"
                )
            else:
                report.ok()
                if verbose:
                    print(
                        f"    ok         {name:<26} {writer_name:<9} "
                        f"{len(expected)} atoms"
                    )
    return report


def frame_tolerance(fmt: str) -> float:
    """Mirrors `crates/chem/src/io/format.rs`'s own `frame_tolerance`
    (#339) exactly, so the Rust and Python sides of "is this trajectory
    close enough" never quietly diverge: `0.01` Angstrom for XTC's real
    quantization step, `1e-3` elsewhere for `f32` rounding noise.
    """
    return 0.01 if fmt == "xtc" else 1e-3


def positions_match(
    a: tuple[tuple[float, float, float], ...],
    b: tuple[tuple[float, float, float], ...],
    tol: float,
) -> bool:
    """The Python counterpart of `format.rs`'s `positions_match` (#339):
    per-axis absolute difference, not a vector norm, for the same reason
    the Rust version chose it -- simpler to reason about at a boundary.
    """
    if len(a) != len(b):
        return False
    return all(
        abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol and abs(x[2] - y[2]) <= tol
        for x, y in zip(a, b)
    )


def boxes_match(
    a: Optional[tuple[float, ...]], b: Optional[tuple[float, ...]], tol: float
) -> bool:
    """Cell lengths at `tol`, angles at a fixed, generous 0.1 degree --
    angles are not subject to a format's own coordinate quantization the
    way lengths are, so they get one tolerance regardless of `fmt`.
    """
    if a is None or b is None:
        return a is b
    lengths_ok = all(abs(a[i] - b[i]) <= tol for i in range(3))
    angles_ok = all(abs(a[i] - b[i]) <= 0.1 for i in range(3, 6))
    return lengths_ok and angles_ok


#: One canonical (positions, box) pair every trajectory-format fixture in
#: this section is authored from or graded against -- ten atoms, not two:
#: `io::xtc`'s own coordinate block stores `size <= 9` atoms as raw,
#: lossless floats (see its module doc) and only compresses above that
#: threshold, so a two-atom fixture would never exercise XTC's real lossy
#: path at all (the exact gap #339 found and fixed for the Rust-side
#: fixture, reused here for the same reason). Positions are not
#: grid-aligned (`i * 1.0033`, not `i * 1.0`) for the same reason: a round
#: number sits exactly on XTC's 0.001nm quantization step and would never
#: exercise real rounding error either.
FIXTURE_POSITIONS = (
    tuple((i * 1.0033, 0.0, 0.0) for i in range(10)),
    tuple((i * 1.0033 + 0.5017, 0.0, 0.0) for i in range(10)),
)
FIXTURE_BOX = (20.0, 20.0, 20.0, 90.0, 90.0, 90.0)


def _render_lammpstrj(
    positions_per_frame: tuple[tuple[tuple[float, float, float], ...], ...],
    box: tuple[float, float, float, float, float, float],
) -> str:
    """Builds a minimal LAMMPS dump text fixture from the same canonical
    (positions, box) every other trajectory format's fixture is authored
    from (#340). Generated, not hand-typed, so twenty numeric lines can
    never silently drift from the tuple every comparison here is graded
    against -- MDAnalysis has no LAMMPS-dump *writer* (confirmed by
    survey: `DumpReader` exists, `DumpWriter` does not), so this is the
    one trajectory format here `chem`, not an external tool, originates.
    """
    lines = []
    for step, positions in enumerate(positions_per_frame):
        lines += [
            "ITEM: TIMESTEP",
            str(step),
            "ITEM: NUMBER OF ATOMS",
            str(len(positions)),
            "ITEM: BOX BOUNDS pp pp pp",
            f"0.0 {box[0]}",
            f"0.0 {box[1]}",
            f"0.0 {box[2]}",
            "ITEM: ATOMS id type x y z",
        ]
        lines += [
            f"{i + 1} 1 {x:.6f} {y:.6f} {z:.6f}" for i, (x, y, z) in enumerate(positions)
        ]
    return "\n".join(lines) + "\n"


LAMMPSTRJ_FIXTURE = _render_lammpstrj(FIXTURE_POSITIONS, FIXTURE_BOX)

TRAJECTORY_FORMATS = ("xtc", "trr", "dcd", "nctraj", "lammpstrj")


def check_trajectory(oracles: list[Oracle], verbose: bool) -> Report:
    """Does chem's `Kind::Frames` round trip agree with MDAnalysis, within
    each format's own real precision (#340)?

    Ignores `oracles` entirely (kept only so it fits `main()`'s generic
    dispatch), the same posture `check_mmcif`/`check_cif_core` already
    take for a structurally different oracle. Comparison is a numeric
    positions/box check with tolerance, not identity -- the first check
    in this harness that isn't, mirroring `format.rs`'s own
    `frame_tolerance`/`positions_match` (#339).

    Shape, for every format: an original fixture with known positions and
    box (MDAnalysis-authored for xtc/trr/dcd/nctraj -- it can write all
    four; a generated text literal for lammpstrj, since MDAnalysis's
    `DumpReader` has no writer counterpart) is round-tripped through
    `chem convert <fmt> --to <fmt>`, then MDAnalysis reads `chem`'s output
    back and the result is compared to the original known values -- the
    same "oracle authors, chem round-trips, oracle re-reads and compares"
    shape `check_ccp4`/`check_mesh` use.
    """
    from oracles import mdanalysis as mda_oracle

    mda_oracle.load_mdanalysis()
    report = Report()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for fmt in TRAJECTORY_FORMATS:
            source = tmp / f"source.{fmt}"
            if fmt == "lammpstrj":
                source.write_text(LAMMPSTRJ_FIXTURE)
            else:
                mda_oracle.write_fixture(source, fmt, FIXTURE_POSITIONS, FIXTURE_BOX)

            written = tmp / f"written.{fmt}"
            if not chem.convert_file(source, fmt, fmt, written):
                report.mismatch(f"{fmt}: chem could not round-trip its own format")
                continue

            back = mda_oracle.read_frames(written, fmt)
            if back is None:
                report.mismatch(f"{fmt}: MDAnalysis could not read what chem wrote back")
                continue

            tol = frame_tolerance(fmt)
            if len(back.positions) != len(FIXTURE_POSITIONS):
                report.mismatch(
                    f"{fmt}: {len(FIXTURE_POSITIONS)} frames in, {len(back.positions)} out"
                )
            elif not all(
                positions_match(expected, got, tol)
                for expected, got in zip(FIXTURE_POSITIONS, back.positions)
            ):
                report.mismatch(f"{fmt}: positions moved by more than {tol} Angstrom")
            elif not boxes_match(FIXTURE_BOX, back.box, tol):
                report.mismatch(f"{fmt}: box disagrees — {FIXTURE_BOX} vs {back.box}")
            else:
                report.ok()
                if verbose:
                    print(f"    ok         {fmt:<10} {len(back.positions)} frames")

        # The named traps a fixture has to be committed, not generated, to
        # pin (#341): a big-endian DCD with a unit cell, and a LAMMPS dump
        # combining a triclinic box with scaled coordinates. MDAnalysis can
        # *read* both (its DCD reader auto-detects endianness; its
        # DumpReader is read-only regardless of format) even though it
        # cannot *author* either -- so unlike the generated fixtures above,
        # these compare MDAnalysis's own read of the original committed
        # file against MDAnalysis's read of chem's round-trip, the same
        # "oracle reads both sides" shape `check_mmcif` already uses,
        # rather than a Python-side hardcoded expectation duplicating what
        # the fixture's own generator script already states.
        for fmt, corpus_dir in (("dcd", CORPUS / "dcd"), ("lammpstrj", CORPUS / "lammpstrj")):
            for path in sorted(corpus_dir.glob(f"*.{fmt}")):
                reference = mda_oracle.read_frames(path, fmt)
                if reference is None:
                    report.mismatch(f"{path.name}: MDAnalysis itself could not read this fixture")
                    continue

                written = tmp / f"corpus_written_{path.stem}.{fmt}"
                if not chem.convert_file(path, fmt, fmt, written):
                    report.mismatch(f"{path.name}: chem could not round-trip this fixture")
                    continue

                ours = mda_oracle.read_frames(written, fmt)
                if ours is None:
                    report.mismatch(f"{path.name}: MDAnalysis cannot read what chem wrote back")
                    continue

                tol = frame_tolerance(fmt)
                if len(ours.positions) != len(reference.positions):
                    report.mismatch(
                        f"{path.name}: {len(reference.positions)} frames in, "
                        f"{len(ours.positions)} out of chem's round trip"
                    )
                elif not all(
                    positions_match(expected, got, tol)
                    for expected, got in zip(reference.positions, ours.positions)
                ):
                    report.mismatch(f"{path.name}: chem's round trip moved positions by more than {tol} Angstrom")
                elif not boxes_match(reference.box, ours.box, tol):
                    report.mismatch(
                        f"{path.name}: box disagrees after chem's round trip — "
                        f"{reference.box} vs {ours.box}"
                    )
                else:
                    report.ok()
                    if verbose:
                        print(f"    ok         {path.name:<30} {len(ours.positions)} frames")
    return report


def check_nctraj_reference(oracles: list[Oracle], verbose: bool) -> Report:
    """Does Amber's own tool accept what `chem` writes as NCTRAJ, and read
    back the same coordinates (#340)?

    A different question than `check_trajectory` asks of NCTRAJ, which is
    why it is a separate check rather than a sixth branch there: this one
    is validity, not agreement with a second implementation. cpptraj is
    "the definition of whether a file is a valid Amber trajectory," the
    same status RDKit had for commonchem in #229 -- so the interesting
    failure here is cpptraj refusing the file outright, not a numeric
    disagreement.
    """
    from oracles import cpptraj as cpptraj_oracle

    cpptraj_oracle.load_cpptraj()
    report = Report()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        source = tmp / "source.lammpstrj"
        source.write_text(LAMMPSTRJ_FIXTURE)
        written = tmp / "written.nctraj"
        if not chem.convert_file(source, "lammpstrj", "nctraj", written):
            report.mismatch("nctraj: chem could not write NCTRAJ from the fixture")
            return report

        n_atoms = len(FIXTURE_POSITIONS[0])
        back = cpptraj_oracle.read_frames(written, n_atoms)
        if back is None:
            report.mismatch(
                "nctraj: cpptraj refused to read what chem wrote — not a valid "
                "Amber trajectory by Amber's own tool"
            )
            return report

        tol = frame_tolerance("nctraj")
        if len(back.positions) != len(FIXTURE_POSITIONS):
            report.mismatch(
                f"nctraj: {len(FIXTURE_POSITIONS)} frames in, {len(back.positions)} "
                "out of cpptraj"
            )
        elif not all(
            positions_match(expected, got, tol)
            for expected, got in zip(FIXTURE_POSITIONS, back.positions)
        ):
            report.mismatch(f"nctraj: cpptraj's read moved positions by more than {tol} Angstrom")
        else:
            report.ok()
            if verbose:
                print(f"    ok         nctraj     {len(back.positions)} frames, cpptraj-validated")
    return report


CCP4_TOLERANCE = 1e-3


def check_ccp4(oracles: list[Oracle], verbose: bool) -> Report:
    """Does chem's CCP4/MRC round trip agree with gemmi's independent read
    (#340)?

    The same shape `check_mmcif` established for gemmi (#224): ignores
    `oracles` entirely, drives gemmi directly on both sides of a round
    trip gemmi itself authors. Grid dimensions and cell compared exactly;
    density values within `CCP4_TOLERANCE`, since the map's own on-disk
    storage is `float32` and neither reader makes an exactness claim
    beyond that.
    """
    from oracles import gemmi as gemmi_oracle

    gemmi_oracle.load_gemmi()
    report = Report()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        source = tmp / "source.ccp4"
        gemmi_oracle.write_ccp4_fixture(source)
        reference = gemmi_oracle.summarize_ccp4(source)
        if reference is None:
            report.mismatch("ccp4: gemmi itself could not read its own fixture")
            return report

        written = tmp / "written.ccp4"
        if not chem.convert_file(source, "ccp4", "ccp4", written):
            report.mismatch("ccp4: chem could not round-trip this file")
            return report

        ours = gemmi_oracle.summarize_ccp4(written)
        problem = _ccp4_disagreement(reference, ours)
        if problem is not None:
            report.mismatch(f"ccp4: {problem}")
        else:
            report.ok()
            if verbose:
                print(f"    ok         ccp4       dims={ours.dims}")

        # The named trap a fixture has to be committed, not generated, to
        # pin (#341): a permuted MAPC/MAPR/MAPS with a non-cubic cell.
        # Same "gemmi reads both sides of chem's round trip" shape as
        # above, just against a permanent file instead of one this run
        # authored and will discard.
        for path in sorted((CORPUS / "ccp4").glob("*.ccp4")):
            reference = gemmi_oracle.summarize_ccp4(path)
            if reference is None:
                report.mismatch(f"{path.name}: gemmi itself could not read this fixture")
                continue

            written = tmp / f"corpus_written_{path.stem}.ccp4"
            if not chem.convert_file(path, "ccp4", "ccp4", written):
                report.mismatch(f"{path.name}: chem could not round-trip this fixture")
                continue

            ours = gemmi_oracle.summarize_ccp4(written)
            problem = _ccp4_disagreement(reference, ours)
            if problem is not None:
                report.mismatch(f"{path.name}: {problem}")
            else:
                report.ok()
                if verbose:
                    print(f"    ok         {path.name:<30} dims={ours.dims}")
    return report


def _ccp4_disagreement(reference, ours) -> Optional[str]:
    """The comparison `check_ccp4` runs on both its generated and its
    committed fixture, factored out so the two don't drift apart.
    """
    if ours is None:
        return "gemmi cannot read what chem wrote back"
    if ours.dims != reference.dims:
        return f"dims disagree — {reference.dims} vs {ours.dims}"
    if not boxes_match(reference.cell, ours.cell, CCP4_TOLERANCE):
        return f"cell disagrees — {reference.cell} vs {ours.cell}"
    if not all(abs(a - b) <= CCP4_TOLERANCE for a, b in zip(reference.values, ours.values)):
        return "density values disagree beyond float32 tolerance"
    return None


MESH_TOLERANCE = 1e-3


def check_mesh(oracles: list[Oracle], verbose: bool) -> Report:
    """Does chem's OBJ/PLY round trip agree with trimesh's independent
    read (#340)?

    The lowest bar the issue names for these two formats: "geometry
    either matches or does not." trimesh authors a closed, watertight
    fixture (a cube) once; chem round-trips it through both OBJ and PLY;
    trimesh re-reads each and the result is compared to trimesh's own
    read of the original, the same "oracle authors, chem round-trips,
    oracle re-reads" shape `check_ccp4` uses.
    """
    from oracles import mesh as mesh_oracle

    mesh_oracle.load_trimesh()
    report = Report()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        source = tmp / "source.obj"
        mesh_oracle.write_fixture(source)
        reference = mesh_oracle.summarize(source)
        if reference is None:
            report.mismatch("mesh: trimesh itself could not read its own fixture")
            return report

        for fmt in ("obj", "ply"):
            written = tmp / f"written.{fmt}"
            if not chem.convert_file(source, "obj", fmt, written):
                report.mismatch(f"{fmt}: chem could not write this format")
                continue

            ours = mesh_oracle.summarize(written)
            problem = _mesh_disagreement(reference, ours)
            if problem is not None:
                report.mismatch(f"{fmt}: {problem}")
            else:
                report.ok()
                if verbose:
                    print(f"    ok         {fmt:<10} {ours.vertex_count}v {ours.face_count}f")

        # The named trap a fixture has to be committed, not generated, to
        # pin (#341): PLY's `binary_big_endian` encoding, the one trimesh
        # cannot itself write (confirmed by survey) -- read only, both
        # sides, the same "oracle reads both sides of chem's round trip"
        # shape the trajectory/CCP4 corpus checks above use.
        for path in sorted((CORPUS / "ply").glob("*.ply")):
            reference = mesh_oracle.summarize(path)
            if reference is None:
                report.mismatch(f"{path.name}: trimesh itself could not read this fixture")
                continue

            written = tmp / f"corpus_written_{path.stem}.ply"
            if not chem.convert_file(path, "ply", "ply", written):
                report.mismatch(f"{path.name}: chem could not round-trip this fixture")
                continue

            ours = mesh_oracle.summarize(written)
            problem = _mesh_disagreement(reference, ours)
            if problem is not None:
                report.mismatch(f"{path.name}: {problem}")
            else:
                report.ok()
                if verbose:
                    print(f"    ok         {path.name:<30} {ours.vertex_count}v {ours.face_count}f")
    return report


def _mesh_disagreement(reference, ours) -> Optional[str]:
    """The comparison `check_mesh` runs on both its generated and its
    committed fixture, factored out so the two don't drift apart.
    """
    if ours is None:
        return "trimesh cannot read what chem wrote"
    if ours.vertex_count != reference.vertex_count:
        return f"{reference.vertex_count} vertices in, {ours.vertex_count} out"
    if ours.face_count != reference.face_count:
        return f"{reference.face_count} faces in, {ours.face_count} out"
    if abs(ours.volume - reference.volume) > MESH_TOLERANCE:
        return f"volume disagrees — {reference.volume} vs {ours.volume}"
    return None


CHECKS = {
    "parse": check_parse,
    "write": check_write,
    "sdf": check_sdf,
    "fp": check_fp,
    "mmcif": check_mmcif,
    "cif_core": check_cif_core,
    "pdb": check_pdb,
    "pdbqt": check_pdbqt,
    "json": check_json,
    "trajectory": check_trajectory,
    "nctraj_reference": check_nctraj_reference,
    "ccp4": check_ccp4,
    "mesh": check_mesh,
}


REGRESSIONS = CORPUS / "regressions"


def baseline(check: str) -> set[str]:
    """Findings already recorded for this check.

    Without this the first run makes CI red forever: the harness found 39 real
    divergences on day one, every one of them a bug in `chem` rather than in
    the harness, and none of them fixable inside the story that built the
    machine. A baseline turns "we are not yet at parity" into "we got no worse",
    which is the only question a CI job can usefully ask before the fidelity
    work is done.

    New findings still fail. Recorded ones are reported as context, so they
    stay visible rather than becoming invisible debt.
    """
    path = REGRESSIONS / f"{check}.md"
    if not path.exists():
        return set()
    return {
        line[2:].strip()
        for line in path.read_text().splitlines()
        if line.startswith("- ")
    }


def promote(check: str, report: Report) -> int:
    """Records this check's findings as the new baseline.

    The loop the harness exists for: a divergence an oracle finds becomes a
    committed file, reviewed in the diff like anything else, and `cargo test`
    keeps its fixtures loadable long after the container is forgotten.
    """
    REGRESSIONS.mkdir(exist_ok=True)
    path = REGRESSIONS / f"{check}.md"
    body = [
        f"# `{check}` — recorded divergences",
        "",
        "Written by `python3 tools/oracle/run.py --promote`. Each line is a way",
        "`chem` differs from an oracle today. They are recorded rather than",
        "fixed: this is the baseline a run compares against, so a *new*",
        "divergence fails CI while these do not.",
        "",
        "Shrinking this file is the fidelity work. Never edit it by hand —",
        "re-run with `--promote` so it matches what the oracles actually say.",
        "",
    ]
    body += [f"- {line}" for line in report.mismatches]
    path.write_text("\n".join(body) + "\n")
    print(f"  recorded {len(report.mismatches)} finding(s) -> {path}")
    return len(report.mismatches)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", choices=sorted(CHECKS), action="append")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--promote", action="store_true")
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--nbits", type=int, default=2048)
    args = parser.parse_args()

    binary = chem.binary()
    if not binary.exists():
        print(f"no chem binary at {binary}", file=sys.stderr)
        print("build it: cargo build --release -p chem --features cli", file=sys.stderr)
        return 2

    oracles = load()
    selected = args.check or sorted(CHECKS)
    print(f"chem     {binary}")
    print(f"oracles  {', '.join(o.name for o in oracles)}")
    print(f"corpus   {CORPUS}")
    print()

    new_findings = 0
    recorded = 0
    for name in selected:
        print(f"[{name}]")
        if name == "fp":
            report = CHECKS[name](oracles, args.verbose, args.radius, args.nbits)
        else:
            report = CHECKS[name](oracles, args.verbose)

        if args.promote:
            for line in report.mismatches:
                print(f"  finding     {line}")
            promote(name, report)
            print()
            continue

        known = baseline(name)
        fresh = [line for line in report.mismatches if line not in known]
        already = [line for line in report.mismatches if line in known]
        fixed = known - set(report.mismatches)

        for line in report.notes:
            print(f"  note        {line}")
        for line in report.disagreements:
            print(f"  DISAGREE    {line}")
        for line in already:
            print(f"  recorded    {line}")
        for line in sorted(fixed):
            print(f"  FIXED       {line}")
        for line in fresh:
            print(f"  NEW         {line}")
        print(
            f"  {report.passed} ok, {len(fresh)} new, {len(already)} recorded, "
            f"{len(fixed)} fixed, {len(report.disagreements)} oracle "
            f"disagreement(s), {len(report.notes)} note(s)"
        )
        if fixed:
            print(
                "  -> a recorded divergence stopped happening; re-run with "
                "--promote to shrink the baseline"
            )
        print()

        new_findings += len(fresh)
        recorded += len(already)

    if args.promote:
        return 0

    # Only a *new* divergence fails. Oracle-vs-oracle disagreements, pinned
    # known gaps and recorded findings are all context.
    print(f"total: {new_findings} new, {recorded} recorded")
    return 1 if new_findings else 0


if __name__ == "__main__":
    sys.exit(main())
