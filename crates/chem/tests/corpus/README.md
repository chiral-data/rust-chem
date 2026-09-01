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

## Pinned gaps

An entry **named** `known-gap-*` records a place `chem` does not behave as the
corpus otherwise promises. They run in both directions: input that should be
rejected and is not, and valid input the parser cannot read yet.

Pinned rather than deleted, so the corpus describes the parser as it behaves.
Closing a gap fails the test that pins it, which is the prompt to rename the
entry and drop its comment — that is how #191's nine entries left this file.

The current ones are in `invalid.smi`: #190, where a run of bond symbols is
silently collapsed so `C##C` parses as ethyne.

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
