# Use an official Rust image as a parent image.
# Using a specific version for reproducibility.
FROM rust:1.89-slim-bullseye

# Install Python and pip for the validation script
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory to the project root
WORKDIR /usr/src/rust-chem

# Install Python dependencies
RUN pip3 install --no-cache-dir "rdkit>=2023.9.1" "numpy>=1.24.0"

# Copy the entire project context into the container
COPY . .


# The compare.py script runs 'cargo build' itself.
# To make the image self-contained, we can pre-build the binary.
# The script is smart enough not to rebuild if it's already there.
# This ensures all dependencies are fetched when the image is built.
RUN cargo build --release -p chem --features cli

# Set the default command to run the Python validation script
CMD ["python3", "crates/chem/compare.py"]
