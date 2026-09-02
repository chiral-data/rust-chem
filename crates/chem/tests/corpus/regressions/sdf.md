# `sdf` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- aromatics.smi: imidazole — rdkit cannot read the SDF chem wrote
- charges.smi: ammonium — SDF round trip, rdkit says protonation (+1 -> -)
- charges.smi: ammonium — SDF round trip, openbabel says [NH4+] -> N
- charges.smi: hydroxide — SDF round trip, rdkit says protonation (-1 -> -)
- charges.smi: hydroxide — SDF round trip, openbabel says [OH-] -> O
- charges.smi: nitrobenzene — rdkit cannot read the SDF chem wrote
- charges.smi: nitrobenzene — SDF round trip, openbabel says [O-][N+](=O)c1ccccc1 -> ON(=O)c1ccccc1
- charges.smi: magnesium-ion — SDF round trip, rdkit says formula (Mg -> Mg.2H); charge (+2 -> -)
- charges.smi: magnesium-ion — SDF round trip, openbabel says [Mg+2] -> [MgH2]
- charges.smi: tetramethylammonium — rdkit cannot read the SDF chem wrote
- charges.smi: tetramethylammonium — SDF round trip, openbabel says C[N+](C)(C)C -> CN(C)(C)C
- charges.smi: glycine-zwitterion — SDF round trip, openbabel says [O-]C(=O)C[NH3+] -> NCC(=O)O
- hard.smi: alanine — SDF round trip, rdkit says stereo parity (0 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- hard.smi: alanine — SDF round trip, openbabel says C[C@@H](C(=O)O)N -> CC(C(=O)O)N
- hard.smi: cis-difluoroethene — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- hard.smi: cis-difluoroethene — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
- hard.smi: sodium-chloride — SDF round trip, rdkit says formula (ClH.Na -> ClH.Na.H); hydrogens (1H; -> 1H;;); protonation (-1 -> -); charge (;+1 -> -)
- hard.smi: sodium-chloride — SDF round trip, openbabel says [Na+].[Cl-] -> [NaH].Cl
- hard.smi: pyrrole — rdkit cannot read the SDF chem wrote
- hard.smi: pyrrole — SDF round trip, openbabel says c1ccc[nH]1 -> C1C=CC=N1
- hard.smi: carbon-13-methane — SDF round trip, rdkit says isotope (1+1 -> -)
- hard.smi: carbon-13-methane — SDF round trip, openbabel says [13CH4] -> C
- stereo.smi: l-alanine — SDF round trip, rdkit says stereo parity (0 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- stereo.smi: l-alanine — SDF round trip, openbabel says C[C@@H](C(=O)O)N -> CC(C(=O)O)N
- stereo.smi: d-alanine — SDF round trip, rdkit says stereo parity (1 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- stereo.smi: d-alanine — SDF round trip, openbabel says C[C@H](C(=O)O)N -> CC(C(=O)O)N
- stereo.smi: cis — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- stereo.smi: cis — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
- stereo.smi: trans-reversed-spelling — SDF round trip, rdkit says double-bond stereo (2-1- -> 2-1+)
- stereo.smi: trans-reversed-spelling — SDF round trip, openbabel says F/C=C\F -> F/C=C/F
- stereo.smi: trans-cyclohexanediol — SDF round trip, rdkit says tetrahedral stereo (5-,6+ -> -)
- stereo.smi: trans-cyclohexanediol — SDF round trip, openbabel says O[C@@H]1CC[C@@H](CC1)O -> OC1CCC(CC1)O
