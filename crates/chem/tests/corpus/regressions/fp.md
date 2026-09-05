# `fp` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- kekule-cyclopentadiene — nearest neighbour disagrees with rdkit: ours cyclohexane (t=0.100), theirs hydroxide (t=0.125)
- sodium-chloride — nearest neighbour disagrees with rdkit: ours alanine (t=0.071), theirs pyridine (t=0.100)
