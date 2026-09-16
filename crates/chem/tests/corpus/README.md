# The differential corpus

Fixtures for the oracle harness in `tools/oracle/`, and for the hermetic
`tests/corpus.rs` that consumes the same files without needing a toolkit.

## Why these are ours

Every molecule here was written for this repository. RDKit's suite is BSD-3 and
could have been vendored with attribution; OpenBabel's is GPL-2 and could not,
this crate being MIT. Authoring our own skips that review entirely and means
every fixture is one we can explain — at the cost of inheriting none of the
breadth those suites have. Widening the corpus is cheap; the promotion loop
below does it automatically whenever an oracle finds something.

## Layout

One file per theme, in the `SMILES<space>name` format `chem` already reads, so
they double as ordinary input. `#` comments carry the reasoning.

| File | What it probes |
|---|---|
| `hard.smi` | The cases named in #4: buckyball, tetrahedral stereo, cis/trans, disconnected fragments, aromatic N-H, isotope |
| `rings.smi` | Ring-closure bookkeeping — digit reuse, fused and bridged systems, the two-digit `%NN` form |
| `charges.smi` | Formal charge, including a zwitterion where dropping charges yields an impossible molecule |
| `aromatics.smi` | Aromaticity perception, the axis where the two oracles most visibly disagree |
| `stereo.smi` | Stereo layers our SMILES writer does not emit yet, so these are expected losses |
| `invalid.smi` | Input that must be **rejected** |
| `regressions/` | Promoted automatically — see below |

The structure formats need whole files rather than a line of SMILES, so they sit
in their own directories. One molecule per file, deliberately: RDKit reads only
the first `MODEL` of a framed PDB, so a multi-record fixture would compare one
molecule against several and pass for the wrong reason.

| File | What it probes |
|---|---|
| `pdb/dipeptide-with-ligand.pdb` | Cells, chains, residues, and per-atom occupancy and B-factor values. `CONECT` is deliberately **partial** — two of eight atoms — which is what real PDB files look like |
| `pdb/water-no-cell.pdb` | The contrast to the file beside it: no `CRYST1`, no `CONECT`, explicit hydrogen *atoms* |
| `pdb/ligand-fully-connected.pdb` | The **hydrogen count**. Every atom in a `CONECT` record, which is what lets an oracle compare a count PDB never states (#293) |
| `mmcif/` | The mmCIF analogues of the first two |
| `sdf/ethanol.mol` | A real molfile as it actually looks: no `$$$$` at all, which is an SDF multi-record separator a single-molecule `.mol` file never has (#318) |
| `bcif/` | The BinaryCIF encoding of the two `mmcif/` fixtures beside them, generated with `chem convert *.cif --to bcif` -- this crate's own writer, not vendored bytes, so what's tested is that its encoder and decoder agree on real structure-sized input (#319). No oracle here: gemmi has no BinaryCIF support at all (confirmed from its source -- `CoorFormat` only knows `Pdb`/`Mmcif`/`Mmjson`, and anything else is silently misread as PDB text rather than rejected), so unlike `pdb/`/`mmcif/` this format is checked at the unit level and against a real, independently-produced `.bcif` file during development, not vendored here |
| `cif_core/quartz.cif` | The small-molecule crystallography dictionary, not mmCIF's: fractional coordinates, an estimated standard deviation on nearly every number (`4.9134(2)`), both a space-group symbol and its International Tables number (#320). Unlike BinaryCIF, gemmi genuinely reads this dictionary (`gemmi.read_small_structure`), so it is oracle-checked |
| `cif_core/symmetry-operators.cif` | A `_symmetry_equiv_pos_as_xyz` loop present and deliberately unused: this crate reads and writes only the asymmetric unit a file states, never expanding symmetry, so this fixture is what proves the operator loop is parsed-and-discarded rather than choking the reader |
| `psf/water.psf` | A small, complete topology: bonds, one angle, donors, and acceptors, with zero exclusions (`!NNB` present but empty) (#321) |
| `psf/exclusions.psf` | The `!NNB` two-array reconstruction specifically: more than one atom has exclusion partners, which an off-by-one in the `IBLO14` cumulative-pointer walk would get wrong while still parsing successfully |

## Pinned gaps

An entry **named** `known-gap-*` records a place `chem` does not behave as the
corpus otherwise promises. They run in both directions: input that should be
rejected and is not, and valid input the parser cannot read yet.

Pinned rather than deleted, so the corpus describes the parser as it behaves.
Closing a gap fails the test that pins it, which is the prompt to rename the
entry and drop its comment — that is how #191's nine entries left this
file, and #190's two.

None are currently pinned.

The marker is the entry's *name*, not a nearby comment. The first version
scanned for a `KNOWN GAP` comment, armed itself on the phrase appearing in a
file's own header, and swept up every line below it.

## `regressions/`

One file per check, recording every way `chem` differs from an oracle today.
Written by `python3 tools/oracle/run.py --promote`, reviewed in the diff like
any other change, never edited by hand.

They are a **baseline, not a bug list**: a run compares against them, so a *new*
divergence fails while these do not. The harness found 23 the day it was built,
none of them fixable inside the story that built it, and without a baseline CI
would be red forever — which asks no useful question. With one it asks "did we
get worse?"

Shrinking these files is the fidelity work. `sdf.md` is mostly formal charge,
exactly as the `Carries` mask predicts; `fp.md` is aromaticity perception (#192).
