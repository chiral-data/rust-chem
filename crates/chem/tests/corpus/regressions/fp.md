# `fp` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- benzene — nearest neighbour disagrees with rdkit: ours pyridine (t=0.333), theirs benzene-kekule (t=1.000)
- benzene-kekule — nearest neighbour disagrees with rdkit: ours kekule-cyclopentadiene (t=0.500), theirs benzene (t=1.000)
- glycine-zwitterion — nearest neighbour disagrees with rdkit: ours trans-cyclohexanediol (t=0.308), theirs alanine (t=0.150)
- kekule-cyclopentadiene — nearest neighbour disagrees with rdkit: ours benzene-kekule (t=0.500), theirs hydroxide (t=0.125)
- pyridine — nearest neighbour disagrees with rdkit: ours pyrrole (t=0.889), theirs benzene (t=0.333)
- pyrrole — nearest neighbour disagrees with rdkit: ours pyridine (t=0.889), theirs imidazole (t=0.235)
- trans-cyclohexanediol — nearest neighbour disagrees with rdkit: ours glycine-zwitterion (t=0.308), theirs norbornane (t=0.231)
