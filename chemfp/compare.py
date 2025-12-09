"""Validate Rust Morgan fingerprints against RDKit with side-by-side comparison."""

import subprocess
import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFingerprintGenerator
except ImportError:
    print("ERROR: RDKit not installed")
    sys.exit(1)


class Validator:
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.binary_path = self.project_root / "target" / "release" / "cli"

    def build(self):
        print("Building Rust binary...")
        result = subprocess.run(
            ["cargo", "build", "--release", "--bin", "cli"],
            cwd=self.project_root,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print("Build failed:", result.stderr)
            sys.exit(1)
        print(f"Built: {self.binary_path}\n")

    def get_rdkit_data(self, smiles, radius=2, nbits=2048):
        """Get complete RDKit data."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Initial invariants (bit-packed)
        invariants = []
        for atom in mol.GetAtoms():
            Z = atom.GetAtomicNum()
            deg = atom.GetDegree()
            H = atom.GetTotalNumHs()
            val = deg + H
            chg = atom.GetFormalCharge()
            arom = atom.GetIsAromatic()

            code = Z % 128
            code |= min(val, 15) << 8
            code |= min(H, 8) << 12
            code |= min(deg, 15) << 16
            code |= ((chg + 5) & 0xF) << 20
            if arom:
                code |= 1 << 24

            invariants.append(code)

        # All environments
        info = {}
        AllChem.GetMorganFingerprint(mol, radius=radius, bitInfo=info)

        envs = []
        for env_hash, features in info.items():
            for atom_idx, rad in features:
                envs.append((atom_idx, rad, env_hash))
        envs.sort(key=lambda x: (x[1], x[0]))

        # Final fingerprint
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = gen.GetFingerprint(mol)
        # fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nbits)
        bits = sorted(fp.GetOnBits())

        return {
            "invariants": invariants,
            "environments": envs,
            "bits": bits,
            "mol": mol,
        }

    def get_rust_fp(self, smiles, radius=2, nbits=2048):
        """Get Rust fingerprint bits."""
        result = subprocess.run(
            [str(self.binary_path), smiles, str(radius), str(nbits)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        output = result.stdout.strip()
        if not output:
            return []
        return sorted(map(int, output.split(",")))

    def compare(self, smiles, name, radius=2, nbits=2048):
        """Compare side-by-side."""
        print(f"\n{'=' * 100}")
        print(f"{name} ({smiles})")
        print(f"{'=' * 100}")

        rdkit = self.get_rdkit_data(smiles, radius, nbits)
        if not rdkit:
            print("ERROR: Failed to parse SMILES")
            return False

        rust_bits = self.get_rust_fp(smiles, radius, nbits)
        if rust_bits is None:
            print("ERROR: Rust failed")
            return False

        # 1. Initial Invariants
        print(f"\n{'STEP 1: INITIAL ATOM INVARIANTS (Bit-Packed)':^100}")
        print(f"{'-' * 100}")
        print(f"{'Atom':<6} {'Symbol':<8} {'Z':<4} {'Deg':<5} {'H':<4} {'Val':<5}")
        print(f"{'Chg':<5} {'Arom':<6} {'RDKit Invariant':<20}")
        print(f"{'-' * 100}")

        for i, inv in enumerate(rdkit["invariants"]):
            atom = rdkit["mol"].GetAtoms()[i]
            print(
                f"{i:<6} {atom.GetSymbol():<8} {atom.GetAtomicNum():<4} "
                f"{atom.GetDegree():<5} {atom.GetTotalNumHs():<4} "
                f"{atom.GetDegree() + atom.GetTotalNumHs():<5} "
                f"{atom.GetFormalCharge():<5} {str(atom.GetIsAromatic()):<6} "
                f"{inv:<10} (0x{inv:08x})"
            )

        # 2. Bonds
        if rdkit["mol"].GetNumBonds() > 0:
            print(f"\n{'STEP 2: BONDS':^100}")
            print(f"{'-' * 100}")
            print(f"{'Bond':<8} {'Atoms':<12} {'Type':<10}")
            print(f"{'-' * 100}")
            for bond in rdkit["mol"].GetBonds():
                print(
                    f"{bond.GetIdx():<8} "
                    f"{bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx():<8} "
                    f"{int(bond.GetBondType()):<10}"
                )

        # 3. Environments - Side by Side
        print(f"\n{'STEP 3: ENVIRONMENTS (All Radius Layers)':^100}")
        print(f"{'-' * 100}")
        print(f"{'RDKit':<50} | {'Rust (Expected)':<49}")
        print(f"{'-' * 100}")
        print(
            f"{'Atom':<6} {'Radius':<8} {'Hash':<18} {'→ Bit':<10}",
            "|",
            f"{'Should Match':<49}",
        )
        print(f"{'-' * 100}")

        for atom_idx, rad, env_hash in rdkit["environments"]:
            bit_id = env_hash % nbits
            print(
                f"{atom_idx:<6} {rad:<8} {env_hash:<10} (0x{env_hash:08x})",
                f" → {bit_id:<10} | ",
                f"{'✓ If your hash matches' if bit_id in rust_bits else '✗ Missing in Rust'}",
            )

        # 4. Final Fingerprint
        print(f"\n{'STEP 4: FINAL FINGERPRINT':^100}")
        print(f"{'-' * 100}")
        print(f"{'Source':<15} {'Bits Set':<12} {'Bit Positions':<73}")
        print(f"{'-' * 100}")
        print(f"{'RDKit':<15} {len(rdkit['bits']):<12} {str(rdkit['bits'][:20]):<73}")
        print(f"{'Rust':<15} {len(rust_bits):<12} {str(rust_bits[:20]):<73}")

        # 5. Comparison
        rdkit_set = set(rdkit["bits"])
        rust_set = set(rust_bits)
        match = rdkit_set == rust_set

        print(f"\n{'RESULT':^100}")
        print(f"{'-' * 100}")
        print(f"Match: {'✓ PASS' if match else '✗ FAIL'}")

        if not match:
            only_rdkit = sorted(rdkit_set - rust_set)
            only_rust = sorted(rust_set - rdkit_set)
            both = sorted(rdkit_set & rust_set)

            print(f"Both have:     {len(both)} bits → {both[:10]}")
            print(f"Only RDKit:    {len(only_rdkit)} bits → {only_rdkit[:10]}")
            print(f"Only Rust:     {len(only_rust)} bits → {only_rust[:10]}")

        return match

    def run(self):
        """Run validation."""
        tests = [
            ("C", "Methane"),
            ("CC", "Ethane"),
            # ("CCC", "Propane"),
        ]

        results = [self.compare(smiles, name) for smiles, name in tests]

        print(f"\n{'=' * 100}")
        print(f"{'SUMMARY':^100}")
        print(f"{'=' * 100}")
        print(f"Passed: {sum(results)}/{len(results)}")

        return 0 if all(results) else 1


if __name__ == "__main__":
    validator = Validator()
    validator.build()
    sys.exit(validator.run())
