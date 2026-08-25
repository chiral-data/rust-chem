set -e

# build the Docker image:
docker build -t chemfp-validator -f crates/chem-fp/Dockerfile .

# run the validation script inside the container:
docker run --rm chemfp-validator
