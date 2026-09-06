# `json` — recorded divergences

Written by `python3 tools/oracle/run.py --promote`. Each line is a way
`chem` differs from an oracle today. They are recorded rather than
fixed: this is the baseline a run compares against, so a *new*
divergence fails CI while these do not.

Shrinking this file is the fidelity work. Never edit it by hand —
re-run with `--promote` so it matches what the oracles actually say.

- aromatics.smi: imidazole — rdkit cannot read the SMILES chem wrote from their commonchem
- charges.smi: ammonium — commonchem read, rdkit says formula (H3N -> N); hydrogens (1H3 -> -); protonation (+1 -> -); charge (- -> +1)
- charges.smi: hydroxide — commonchem read, rdkit says formula (H2O -> O); hydrogens (1H2 -> -); protonation (-1 -> -); charge (- -> -1)
- charges.smi: nitrobenzene — commonchem write, formula C6H5NO2 -> C6H7NO2
- charges.smi: glycine-zwitterion — commonchem write, formula C2H5NO2 -> C2H7NO2
- charges.smi: glycine-zwitterion — commonchem read, rdkit says formula (C2H5NO2 -> C2H3NO2); hydrogens (1,3H2,(H,4,5) -> 1H2,(H,4,5)); protonation (- -> -1); charge (- -> +1)
- hard.smi: alanine — commonchem read, rdkit says formula (C3H7NO2 -> C3H6NO2); hydrogens (2H,4H2,1H3,(H,5,6) -> 4H2,1H3,(H,5,6)); stereo parity (0 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- hard.smi: sodium-chloride — commonchem write, formula ClNa -> H2ClNa
- hard.smi: pyrrole — rdkit cannot read the SMILES chem wrote from their commonchem
- hard.smi: carbon-13-methane — commonchem read, rdkit says formula (CH4 -> C); hydrogens (1H4 -> -)
- stereo.smi: l-alanine — commonchem read, rdkit says formula (C3H7NO2 -> C3H6NO2); hydrogens (2H,4H2,1H3,(H,5,6) -> 4H2,1H3,(H,5,6)); stereo parity (0 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- stereo.smi: d-alanine — commonchem read, rdkit says formula (C3H7NO2 -> C3H6NO2); hydrogens (2H,4H2,1H3,(H,5,6) -> 4H2,1H3,(H,5,6)); stereo parity (1 -> -); stereo type (1 -> -); tetrahedral stereo (2- -> -)
- stereo.smi: trans-cyclohexanediol — commonchem read, rdkit says formula (C6H12O2 -> C6H10O2); hydrogens (5-8H,1-4H2 -> 7-8H,1-4H2); tetrahedral stereo (5-,6+ -> -)
