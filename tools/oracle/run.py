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
from dataclasses import dataclass, field
from pathlib import Path

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
    thing comparable across hashes — 32 of 35 corpus molecules match RDKit's
    exactly — and it is reported as a note, since ranking cannot see it (#253).
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

            if len(mine) != len(yours):
                report.note(
                    f"{name}: {len(mine)} bits vs {oracle.name}'s {len(yours)} "
                    f"— different environment counts, not just a different hash"
                )

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


CHECKS = {
    "parse": check_parse,
    "write": check_write,
    "sdf": check_sdf,
    "fp": check_fp,
    "mmcif": check_mmcif,
    "json": check_json,
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
