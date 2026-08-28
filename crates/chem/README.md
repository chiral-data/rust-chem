# chem

Cheminformatics in Rust: molecules, SMILES and SDF, Morgan fingerprints,
Tanimoto similarity search, and 2D depiction — with optional GPU kernels and a
command-line tool.

```toml
[dependencies]
chem = "0.6"
```

```rust
use chem::io::smiles::parse_smiles;
use chem::fp::morgan::MorganFingerprint;

let mol = parse_smiles("c1ccccc1O")?;
let fp = MorganFingerprint::get_fingerprint_as_bitvec(
    &mol, 2, 2048, None, None, false, true, false,
)?;
```

## Modules

| module | what it holds |
| --- | --- |
| [`chem::core`](https://docs.rs/chem/latest/chem/core/) | atoms, bonds, molecules, the molecular graph, ring perception, 2D layout |
| [`chem::io`](https://docs.rs/chem/latest/chem/io/) | SMILES and SDF reading and writing, aromaticity perception, whole-file readers |
| [`chem::fp`](https://docs.rs/chem/latest/chem/fp/) | Morgan fingerprints and Tanimoto similarity, on the CPU |
| [`chem::draw`](https://docs.rs/chem/latest/chem/draw/) | 2D depiction as drawable primitives, and an SVG writer |
| [`chem::search`](https://docs.rs/chem/latest/chem/search/) | fingerprint generation and similarity search over a chosen backend |
| [`chem::gpu`](https://docs.rs/chem/latest/chem/gpu/) | wgpu compute kernels — **requires the `gpu` feature** |

Longer notes: [core data structures](docs/core.md), [the Morgan
algorithm](docs/morgan.md).

## Examples

Runnable, one idea each. `cargo run --example <name>`, and docs.rs renders them
on the crate page.

| example | shows |
| --- | --- |
| `build_a_molecule` | atoms and bonds by hand, implicit hydrogens, computing 2D coordinates |
| `read_and_convert` | reading a `.smi` or `.sdf`, what happens to records that fail, perceiving aromaticity, writing SDF back |
| `rank_by_similarity` | fingerprinting a small library and ranking it by Tanimoto |
| `smiles_to_svg` | parse, lay out, draw — takes a SMILES and an output path |
| `hello_gpu` | the GPU context and a compute pipeline (needs `--features gpu`) |
| `search_cpu_vs_gpu` | where the GPU starts paying for itself (needs `--features gpu`) |

## Features

`gpu` is **off by default**. `wgpu` is a large dependency tree and irrelevant to
anyone parsing SMILES or fingerprinting on a CPU, which is what every GPU-less
machine and every WebAssembly build already do.

```toml
chem = { version = "0.6", features = ["gpu"] }
```

Without it, `chem::search` still offers its whole API — it simply never finds a
device, and says so through `FingerprintSearch::gpu_init_error()`.

The GPU is worth reaching for when a dataset stays resident and is queried
repeatedly: warm similarity search over 200,000 targets runs ~56× faster than the
CPU and is nearly flat in dataset size. It does **not** pay for a single
one-shot query, where the upload dominates. `cargo run --release --example
search_cpu_vs_gpu --features gpu` measures this on your own hardware.

## The command-line tool

```sh
cargo install chem --features cli
```

```sh
chem info mols.smi                       # what is in this file
chem fp mols.smi -o fps.tsv              # Morgan fingerprints
chem search fps.tsv --query 'c1ccccc1O'  # rank by similarity
chem aromatic mols.smi                   # perceive aromatic rings
chem coords mols.smi -o mols.sdf         # 2D layout, written as SDF
chem draw mols.smi --outdir svg/         # one SVG per molecule
```

Reads a named file or standard input, writes to `-o` or standard output, and puts
progress and warnings on standard error — so `chem aromatic x.smi | chem fp`
composes. Exit codes distinguish "the file was unusable" from "nothing matched".

## Known limitations

- Aromaticity perception is a Hückel heuristic over the SSSR, not a full model.
  It handles the common heteroaromatics; exotic systems and charged aromatics are
  untested.
- SDF output writes atoms, bonds and 2D coordinates. Charges, isotopes, chirality
  and bond stereo are not yet written.
- Aromatic bonds are written to SDF as type 4, which round-trips here but is a
  query bond type in strict molfile.

## Licence

MIT.
