# `sdf` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- hard.smi: cis-difluoroethene — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- hard.smi: cis-difluoroethene — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
- stereo.smi: cis — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- stereo.smi: cis — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
- stereo.smi: trans-reversed-spelling — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- stereo.smi: trans-reversed-spelling — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
