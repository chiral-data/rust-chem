set -e

# Validates this crate's Morgan fingerprints against RDKit's, which is the only
# third-party check in the repository. Needs Docker; the image installs RDKit.

# build the Docker image:
docker build -t chem-rdkit-validator -f crates/chem/Dockerfile .

# run the validation script inside the container:
docker run --rm chem-rdkit-validator
