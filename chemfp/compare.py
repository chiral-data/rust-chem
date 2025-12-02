"""
Validate Rust Morgan fingerprints against RDKit.
Builds Rust binary, generates fingerprints, and compares results.
"""

import subprocess
import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
except ImportError:
    print("ERROR: RDKit not installed")
    print("Install with: pip install -r requirements.txt")
    sys.exit(1)


class RustFingerprintValidator:
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.binary_name = "cli"
        self.binary_path = None

    def build_rust_binary(self):
        """Build the Rust CLI tool."""
        print("Building Rust fingerprint CLI...")

        result = subprocess.run(
            ["cargo", "build", "--bin", self.binary_name],
            cwd=self.project_root,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print("ERROR: Failed to build Rust binary")
            print(result.stderr)
            sys.exit(1)

        self.binary_path = self.project_root / "target" / "debug" / self.binary_name

        if not self.binary_path.exists():
            print(f"ERROR: Binary not found at {self.binary_path}")
            sys.exit(1)

        print(f"Built binary: {self.binary_path}")

    def get_rust_fingerprint(self, smiles, radius=2, nbits=2048):
        """Get fingerprint from Rust implementation."""
        if not self.binary_path:
            self.build_rust_binary()

        result = subprocess.run(
            [str(self.binary_path), smiles, str(radius), str(nbits)],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"ERROR: Rust fingerprint failed for {smiles}")
            print(result.stderr)
            return None

        # Parse output: comma-separated list of bit positions
        try:
            output = result.stdout.strip()
            if not output:
                return set()
            bit_positions = set(map(int, output.split(",")))
            return bit_positions
        except Exception as e:
            print(f"ERROR parsing Rust output: {e}")
            return None

    def get_rdkit_fingerprint(self, smiles, radius=2, nbits=2048):
        """Get fingerprint from RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # ✅ Use GetMorganFingerprint with bitInfo to get feature hashes
        from rdkit.Chem import AllChem

        bit_info = {}
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=radius, nBits=nbits, bitInfo=bit_info
        )

        # bit_info format: {bit_position: [(atom_idx, radius), ...]}
        # We want to count unique features across all radii
        unique_features = set()
        for bit_pos, atom_list in bit_info.items():
            for atom_idx, feature_radius in atom_list:
                unique_features.add((atom_idx, feature_radius))

        print(f"  RDKit unique features: {len(unique_features)}")
        print(f"  RDKit bit positions: {len(bit_info)}")

        # Get set bits
        bit_positions = set(fp.GetOnBits())
        return bit_positions

    # def get_rdkit_fingerprint(self, smiles, radius=2, nbits=2048):
    #     """Get fingerprint from RDKit."""
    #     mol = Chem.MolFromSmiles(smiles)
    #     if mol is None:
    #         return None
    #
    #     gen = rdFingerprintGenerator.GetMorganGenerator(
    #         radius=radius, fpSize=nbits)
    #
    #     info = rdFingerprintGenerator.AdditionalOutput()
    #     fp = gen.GetFingerprint(mol, additionalOutput=info)
    #
    #     print(f"  RDKit atom counts: {sorted(info.GetAtomCounts().keys())}")
    #     print(f"  RDKit num atoms: {len(info.GetAtomCounts())}")
    #
    #     # Get set bits
    #     bit_positions = set(fp.GetOnBits())
    #     return bit_positions

    def compare_fingerprints(self, rust_bits, rdkit_bits):
        """Compare two fingerprints and return statistics."""
        if rust_bits is None or rdkit_bits is None:
            return None

        intersection = rust_bits & rdkit_bits
        union = rust_bits | rdkit_bits

        only_rust = rust_bits - rdkit_bits
        only_rdkit = rdkit_bits - rust_bits

        if len(union) == 0:
            similarity = 1.0
        else:
            similarity = len(intersection) / len(union)

        return {
            "rust_bits": len(rust_bits),
            "rdkit_bits": len(rdkit_bits),
            "intersection": len(intersection),
            "union": len(union),
            "similarity": similarity,
            "only_rust": len(only_rust),
            "only_rdkit": len(only_rdkit),
            "match": similarity == 1.0,
        }

    def validate_test_set(self):
        """Run validation on a test set of molecules."""
        test_molecules = [
            ("C", "Methane"),
            ("CC", "Ethane"),
            ("CCC", "Propane"),
            ("C=C", "Ethene"),
            ("C#C", "Ethyne"),
            ("CCO", "Ethanol"),
            ("CC(C)C", "Isobutane"),
            ("c1ccccc1", "Benzene"),
            ("CC(=O)O", "Acetic acid"),
            ("CC(C)(C)C", "Neopentane"),
            ("C1CCCCC1", "Cyclohexane"),
            ("c1ccc(O)cc1", "Phenol"),
            ("CCN", "Ethylamine"),
            ("O", "Water"),
            ("N", "Ammonia"),
        ]

        print("\n" + "=" * 80)
        print("VALIDATION RESULTS")
        print("=" * 80)

        results = []

        for smiles, name in test_molecules:
            print(f"\nTesting: {name} ({smiles})")

            rust_fp = self.get_rust_fingerprint(smiles)
            rdkit_fp = self.get_rdkit_fingerprint(smiles)

            if rust_fp is None or rdkit_fp is None:
                print("  SKIP: Failed to generate fingerprint")
                continue

            stats = self.compare_fingerprints(rust_fp, rdkit_fp)
            results.append((name, smiles, stats))

            print(f"  Rust bits:   {stats['rust_bits']}")
            print(f"  RDKit bits:  {stats['rdkit_bits']}")
            print(f"  Similarity:  {stats['similarity']:.4f}")
            print(f"  Match:       {'PASS' if stats['match'] else 'FAIL'}")

            if not stats["match"]:
                print(f"  Only Rust:   {stats['only_rust']} bits")
                print(f"  Only RDKit:  {stats['only_rdkit']} bits")

        # Summary
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)

        total = len(results)
        matches = sum(1 for _, _, s in results if s["match"])
        avg_similarity = (
            sum(s["similarity"]
                for _, _, s in results) / total if total > 0 else 0
        )

        print(f"Total tested:     {total}")
        print(f"Exact matches:    {matches}")
        print(f"Match rate:       {matches / total * 100:.1f}%")
        print(f"Avg similarity:   {avg_similarity:.4f}")

        if matches == total:
            print("\nRESULT: ALL TESTS PASSED")
            return 0
        else:
            print("\nRESULT: SOME TESTS FAILED")
            return 1


def main():
    print("RDKit Morgan Fingerprint Validator")
    print("=" * 80)

    validator = RustFingerprintValidator()
    exit_code = validator.validate_test_set()

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
