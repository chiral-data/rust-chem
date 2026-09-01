# The differential harness

RDKit and OpenBabel as oracles for `chem`, so the milestone's two promises can
both hold: behavioural parity with those toolkits, in pure Rust, with neither
of them anywhere in the dependency graph.

```sh
docker build -t chem-oracle -f tools/oracle/Dockerfile .
docker run --rm chem-oracle                                           # every check
docker run --rm chem-oracle python3 tools/oracle/run.py --check write --verbose
```

Nothing here is a workspace member, a dev-dependency, or named in `Cargo.toml`.
It drives `target/release/chem` as a subprocess, so a developer with neither
Python nor a toolkit installed still runs the whole `cargo` gate.

## The four checks

| Check | Question |
|---|---|
| `parse` | Does `chem` accept what the oracles accept, and refuse what they refuse? |
| `write` | Does a round trip through our SMILES writer preserve the molecule? |
| `sdf` | Does a molecule survive `chem coords` to SDF and back? |
| `fp` | Do our fingerprints rank molecules the way RDKit's do? |

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

`fp` is the exception, deliberately. Morgan bit *positions* are the output of a
particular hash, and ours differ from RDKit's (#192) — that check compares
nearest-neighbour agreement instead, which holds for any chemically equivalent
fingerprint.

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
without touching a check.

gemmi belongs here and is deliberately absent: its value is mmCIF and PDB, and
`chem` reads neither yet, so it would be a pinned dependency exercising nothing.

`load()` refuses to return an oracle that cannot read ethanol. OpenBabel does
exactly that when a system library is missing — its plugin loader aborts, every
format silently fails to register, and it reports that no molecule on earth
parses. It answers every question confidently and wrongly, which produced 31
false mismatches the first time this ran.
