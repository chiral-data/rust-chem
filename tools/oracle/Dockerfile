# The differential harness: RDKit and OpenBabel as oracles for `chem`.
#
# Its own image rather than anything in the cargo workspace, so that a
# developer with neither Python nor a cheminformatics toolkit installed still
# runs the whole `cargo` gate. Nothing here is ever named in Cargo.toml.
#
#   docker build -t chem-oracle -f tools/oracle/Dockerfile .
#   docker run --rm chem-oracle
#
# Build context is the repository root, not this directory — the image needs
# the crate to build the binary it drives.

# --- deps ---------------------------------------------------------------------
# The toolkits install before any source is copied, so editing Rust or the
# harness re-runs neither pip nor apt. That ordering is the difference between
# a five-second rebuild and a three-minute one.
FROM rust:1.89-slim-bookworm AS deps

# libxrender1 and friends are not optional. Without them OpenBabel's plugin
# loader aborts partway through registering formats, `SetInFormat("smi")`
# returns false, and every molecule looks unparseable — a blind oracle that
# still answers, which is worse than one that fails to import.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 python3-pip \
        libxrender1 libxext6 libsm6 libice6 && \
    rm -rf /var/lib/apt/lists/*

COPY tools/oracle/requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /tmp/requirements.txt

# --- harness ------------------------------------------------------------------
FROM deps AS harness

WORKDIR /usr/src/rust-chem
COPY . .

# Built once here rather than by the harness at run time. `compare.py` used to
# shell out to `cargo build` on every run, which meant a container that had to
# reach the network to do its job.
RUN cargo build --release -p chem --features cli

ENV CHEM_BIN=/usr/src/rust-chem/target/release/chem
CMD ["python3", "tools/oracle/run.py"]
