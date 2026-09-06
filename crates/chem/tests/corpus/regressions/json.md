# `json` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- hard.smi: dimethylcarbene — commonchem write, rdkit says formula (C3H6 -> C3H8); hydrogens (1-2H3 -> 3H2,1-2H3)
- hard.smi: carbon-13-atom — commonchem write, rdkit says formula (C -> CH4); hydrogens (- -> 1H4)
