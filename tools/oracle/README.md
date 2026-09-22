# The differential harness

RDKit and OpenBabel as oracles for `chem`, so the milestone's two promises can
both hold: behavioural parity with those toolkits, in pure Rust, with neither
of them anywhere in the dependency graph. gemmi (#224) joins them for mmCIF
and PDB specifically, where neither RDKit (no support at all) nor OpenBabel
(a different, incompatible mmCIF dialect) can judge the result — see
"Adding an oracle" below for why it isn't a third `Oracle` implementation.

v0.9.0 (#340) added four more: MDAnalysis for the five trajectory formats,
gemmi again for CCP4/MRC, trimesh for OBJ/PLY, and cpptraj — Amber's own
tool, not a second opinion — as NCTRAJ's reference oracle specifically. None
of RDKit/OpenBabel/gemmi/Meeko can read any of these; without new oracles
they would ship on self-consistency alone, exactly the standard #173
rejected. See "Tolerance, not identity" below for what's different about
judging a trajectory.

```sh
docker build -t chem-oracle -f tools/oracle/Dockerfile .
docker run --rm chem-oracle                                           # every check
docker run --rm chem-oracle python3 tools/oracle/run.py --check write --verbose
```

Nothing here is a workspace member, a dev-dependency, or named in `Cargo.toml`.
It drives `target/release/chem` as a subprocess, so a developer with neither
Python nor a toolkit installed still runs the whole `cargo` gate.

## The checks

| Check | Question |
|---|---|
| `parse` | Does `chem` accept what the oracles accept, and refuse what they refuse? |
| `write` | Does a round trip through our SMILES writer preserve the molecule? |
| `sdf` | Does a molecule survive `chem coords` to SDF and back? |
| `fp` | Do our fingerprints rank molecules the way RDKit's do? |
| `mmcif` | Does `chem`'s mmCIF round trip agree with gemmi's independent read? |
| `cif_core` | Does `chem`'s CIF-core round trip agree with gemmi's independent read? |
| `pdb` | Does `chem`'s PDB round trip agree with gemmi, including per-atom occupancy/B-factor? |
| `pdbqt` | Does `chem` read what OpenBabel and Meeko each write as PDBQT? |
| `json` | Does commonchem JSON agree with RDKit, which defines the format? |
| `trajectory` | Does chem's XTC/TRR/DCD/NCTRAJ/LAMMPS-trajectory round trip agree with MDAnalysis, within each format's own real precision? |
| `nctraj_reference` | Does Amber's own tool accept what `chem` writes as NCTRAJ, and read back the same coordinates? |
| `ccp4` | Does chem's CCP4/MRC round trip agree with gemmi's independent read? |
| `mesh` | Does chem's OBJ/PLY round trip agree with trimesh's independent read? |

`json` (#229) is the only check whose oracle is the format's *reference
implementation* rather than a second opinion, so it runs both directions: what
`chem` writes must read back as the same molecule, and what RDKit writes (in
its own `rdkitjson` dialect, which `chem` accepts and never emits) must survive
being read. `nctraj_reference` (#340) is the same shape for a different
reason: it is not a second opinion at all, but a validity question — is this
a file Amber's own tool accepts — so the interesting failure there is cpptraj
refusing the file outright, not a numeric disagreement.

BinaryCIF was surveyed for #340 and has no check: the pinned (and latest
installable) `gemmi==0.7.5` cannot read it at all, confirmed directly against
the already-committed `bcif/*.bcif` fixtures (#319). The gap
`crates/chem/tests/corpus/README.md` already discloses stands.

## Identity, not strings

`chem` has no canonical SMILES writer yet, so comparing our output text against
an oracle's would differ cosmetically on nearly every molecule and mean nothing.
Atom order and ring-closure digits are free variables.

The comparison key is **InChI**, canonical independently of anyone's writer and
*layered* — so a mismatch says `protonation (+1 -> -)` or
`isotope (1+1 -> -)` rather than "differs". Localising a finding is most of the
work of acting on it.

OpenBabel has no InChI here, so its key is its own canonical SMILES: canonical
within OpenBabel, unlayered, and therefore less precise. It earns its place by
disagreeing with RDKit rather than by being exact.

`json` needs a second key alongside InChI, for a reason worth stating: InChI's
`/p` layer normalises mobile protons away, so `CC(=O)[O-]` and the impossible
`CC(=O)[OH2-]` have the *same* InChI. commonchem is the first format `chem`
writes that states a per-atom hydrogen count, so an InChI-only comparison would
be blind to the one thing it newly carries. That check therefore also compares
molecular formula, which is order-independent and counts every hydrogen — and
which is what caught #240.

`fp` is the exception, deliberately. Morgan bit *positions* are the output of a
particular hash, and ours differ from RDKit's (#192) — that check compares
nearest-neighbour agreement instead, which holds for any chemically equivalent
fingerprint.

## Tolerance, not identity

Every check above compares by identity: two things either produce the same
InChI, the same structural summary, the same set of atoms, or they don't.
`trajectory`, `nctraj_reference` and `ccp4` (#340) are the first checks that
can't be — a trajectory format's own stored precision is real, not a defect,
so a byte- or even float-identical round trip is the wrong question to ask.

The tolerance is not invented per check. `crates/chem/src/io/format.rs`
already solved this exact problem for its own internal fidelity matrix
(`frame_tolerance`/`positions_match`, #339): `0.01` Angstrom for XTC's real
quantization step, `1e-3` elsewhere for `f32` rounding noise. `run.py`'s
`frame_tolerance`/`positions_match` are the Python side of the same two
numbers, so "is this trajectory close enough" can't quietly answer
differently depending on which language is asking. `ccp4` uses a plain
`1e-3` float32 tolerance for density values, for the same reason: the map's
own on-disk storage is `float32`, and neither reader claims exactness beyond
that.

## Fixtures authored by the oracle, not committed

None of `trajectory`/`nctraj_reference`/`ccp4`/`mesh`'s eight format families
have a real fixture in `crates/chem/tests/corpus/` today — #341 ("Binary
fixtures we can still explain") is the later story about committing real,
explainable binaries there. Rather than wait for it or duplicate it, every
one of these checks authors its own small fixture at run time, using the
oracle tool itself as the author wherever it can write one (MDAnalysis for
XTC/TRR/DCD/NCTRAJ, gemmi for CCP4, trimesh for OBJ/PLY) — which has the
added benefit of making the fixture's origin genuinely independent of
`chem`, so `chem` reading it back is a real test of the reader rather than a
circular one. The one exception is LAMMPS trajectory: MDAnalysis's
`DumpReader` is read-only (no `DumpWriter` exists), so that one fixture is
generated as a plain text literal instead, from the same canonical
positions every other trajectory format's fixture uses.

## What fails a run, and what does not

Only a **new** divergence fails.

- **new** — `chem` disagrees with an oracle in a way not already recorded. Exit 1.
- **recorded** — a divergence already in `crates/chem/tests/corpus/regressions/`.
  Reported, never fatal.
- **oracle disagreement** — the oracles differing from *each other*. A fact about
  the toolkits, belonging in the fidelity table rather than resolved by picking a
  favourite.
- **pinned gap** — a corpus entry named `known-gap-*`, recording something `chem`
  wrongly accepts (#190) or cannot yet read (#191).
- **fixed** — a recorded divergence that stopped happening. Re-run with
  `--promote` to shrink the baseline.

The baseline exists because the harness found 23 real divergences the day it was
built, none of them fixable inside the story that built the machine. Without it
CI would be red forever, which asks no useful question. With it, CI asks "did we
get worse?" — and shrinking `regressions/` is the fidelity work itself.

## Recording a new baseline

The container writes into the repo, so the corpus has to be mounted:

```sh
docker run --rm -v "$PWD/crates/chem/tests/corpus:/usr/src/rust-chem/crates/chem/tests/corpus" \
    chem-oracle python3 tools/oracle/run.py --promote
```

Then rebuild the image so it carries the new baseline, and commit the result —
the files are reviewed in the diff like any other change. Never edit them by
hand.

## Adding an oracle

Implement `parses`, `identity`, `identity_of_sdf` and optionally `fingerprint`,
then add it to `oracles/load()`. That narrow surface is what lets a third slot in
without touching a check — as long as the new oracle answers a SMILES-shaped
question at all.

gemmi (#224) doesn't: it has no bare-SMILES `parses`/`identity`, since it only
reads structures (mmCIF, PDB), not small molecules from a string. Forcing it
into `Oracle`'s shape would be the same mistake `load()`'s sanity check exists
to catch — instead it's `oracles/gemmi.py`, a separate structural interface
(`summarize(text) -> atom count, cell, chains, residues`), loaded and
sanity-gated on its own (`load_gemmi()`, against a trivial mmCIF fixture
instead of `"CCO"`), and used by its own check (`mmcif`) rather than folded
into `parse`/`write`/`sdf`/`fp`. A future oracle that answers the *same*
SMILES-shaped question `Oracle` already asks is the common case and belongs
in `load()`'s list directly; one that answers a structurally different
question, the way gemmi does, earns its own path instead.

`load()` refuses to return an oracle that cannot read ethanol. OpenBabel does
exactly that when a system library is missing — its plugin loader aborts, every
format silently fails to register, and it reports that no molecule on earth
parses. It answers every question confidently and wrongly, which produced 31
false mismatches the first time this ran.

`oracles/mdanalysis.py` and `oracles/mesh.py` (#340) follow gemmi's path, not
`load()`'s: neither answers a SMILES-shaped question, so each is its own
module with its own `load_*` gate (`load_mdanalysis`, `load_trimesh`),
sanity-checked by writing and reading back its own trivial fixture before
its answers are trusted.

`oracles/cpptraj.py` (#340) is a third shape again: cpptraj is a CLI tool,
not a Python library, so `load_cpptraj` runs `cpptraj --version` as a
subprocess rather than importing anything — the same discipline as every
other `load_*` gate, adapted to how this one oracle is actually reached.
