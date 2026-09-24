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
| `prmtop/water.prmtop` | A small, complete topology: bonds, one angle, an `ATOMIC_NUMBER` section, and a charge stated in Amber units to prove the 18.2223 scaling round-trips (#322) |
| `prmtop/exclusions.prmtop` | `NUMBER_EXCLUDED_ATOMS`/`EXCLUDED_ATOMS_LIST`'s sequential-consume shape specifically: more than one atom has real exclusions, and more than one uses the `1`/`0` "no exclusions" sentinel, which a miscounted consume would get wrong while still parsing successfully |
| `prmtop/no_atomic_number.prmtop` | No `ATOMIC_NUMBER` section at all -- the mass-based element-inference fallback every real PRMTOP reader (this one included) needs for pre-AmberTools-12 files |
| `top/urea_like.top` | A small, complete single-moleculetype topology: bonds, one angle, one proper (funct 9) and one improper (funct 4) dihedral, an `[ exclusions ]` entry, and an unresolved `#include` that this reader records but never follows (#323) |
| `top/multi_moleculetype.top` | Two `[ moleculetype ]` blocks in one file -- normal GROMACS practice, not a round-trip artifact -- proving one `Molecule` per block, plus a `[ molecules ]` count (`Water 1000`) this crate never expands |
| `top/conditional.top` | `#ifdef`/`#else` driven by an in-file `#define`: both branches state a same-atoms `[ angles ]` line in a different order, so reading the wrong branch (or both) is directly visible in the result |
| `lammps/atomic.data` | A small, topology-free `atomic`-style system with a nonzero box origin (`xlo`/`ylo`/`zlo` all `10.0`), proving both that a bondless file reads fine and that the origin-loss decision is real (#324) |
| `lammps/full.data` | A `full`-style system (bonds, one angle, a proper and an improper dihedral) with an explicit `Atoms # full` comment, exercising the comment-based atom-style resolution path |
| `lammps/ambiguous_style.data` | A 6-column `Atoms` section with no style comment -- ambiguous between `charge` and `molecular`/`bond`/`angle` -- proving `AmbiguousAtomStyle` fires by default and `LammpsReadOptions::atom_style` is the escape hatch |
| `lammps/coarse_grained.data` | A `Masses` entry stating `1.0`, a reduced-unit bead mass matching no real element -- proving this reads as `Element::UNKNOWN`, not a failure or a guess |
| `lammps/triclinic.data` | A genuinely non-orthogonal box (`xy xz yz` all nonzero) -- pins the closed-form box-to-`UnitCell` conversion's shape, which a sign or swapped cosine would still "succeed" while getting wrong |
| `edr/em-first-3-frames.edr` | Real GROMACS 2022 output: the names block and first three frames of an energy minimisation, cut at a frame boundary so it is still a whole file (#399). Its expected values are what `gmx energy` prints for the same bytes |
| `dcd/big_endian_with_cell.dcd` | Big-endian and a unit cell at once (#341), the byte-swap path and the cell heuristic exercised together. Generated -- see "Fixtures we generate" below |
| `ccp4/permuted_axes.ccp4` | A permuted `MAPC`/`MAPR`/`MAPS` order with a non-cubic cell (#341) -- proves the axis relabelling is honoured independently of file storage order. Generated |
| `ply/binary_big_endian.ply` | The one PLY encoding this crate never writes (#341); `ascii`/`binary_little_endian` need no fixture at all, since round-tripping through this crate's own writer already covers them. Generated |
| `lammpstrj/triclinic_scaled.lammpstrj` | Scaled coordinates *in* a triclinic box (#341) -- the conversion genuinely depends on the box's own edge vectors, not just its bounding lengths, unlike an orthogonal box |

Two traps named alongside these four turned out to already be pinned:
`prmtop/water.prmtop`'s stated charge already exercises the 18.2223 Amber
scaling factor (#322, above), and `sdf/ethanol.mol` already has zero `$$$$`
occurrences (#318, above). Neither needed new work for #341.

## Fixtures we generate

A fixture is easy to hand-author for a `.smi` line and hard for a
big-endian DCD or a CCP4 map with permuted axes — nobody hand-types those
bytes. The resolution that keeps this corpus's own "ours by construction"
rule (above) without asking anyone to: a **script committed beside the
fixture it produces** (`generate_*.py` next to the file it writes). The
bytes are still ours, and reading the script is how you explain them,
rather than either vendoring someone else's file or hand-transcribing a
byte array into a test.

Generators may use the same toolkits `tools/oracle/`'s Docker image
already pins (MDAnalysis, gemmi, trimesh; `plyfile` is one exception, a
second, generator-only pin — see `ply/generate_binary_big_endian.py` for
why) — none of that makes a toolkit a dependency of `chem` itself, the
same posture the oracle harness already takes. Where a needed capability
doesn't exist in any pinned tool (no tool here can *write* a big-endian
DCD, confirmed by survey), the generator hand-assembles the bytes
directly instead, the same way a Rust test's own in-process builder would
have — just committed, not thrown away after one run.

There is no automated check that re-runs a generator and diffs its output
against the committed file on every push; that would be real, but is a
different, heavier story than this one. The discipline instead mirrors
`--promote`'s own: if you touch a generator, re-run it, diff the result
against what's committed, and commit both together, reviewed like any
other change.

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
