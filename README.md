# rust-chem

Cheminformatics in Rust. Two crates in this workspace:

- **[`chem`](crates/chem/)** — the library and the `chem` command-line tool.
  Published to [crates.io](https://crates.io/crates/chem); documentation at
  [docs.rs/chem](https://docs.rs/chem). See [its README](crates/chem/README.md)
  for the API and the feature flags.
- **`chem-app`** — a desktop and browser workbench built on it. Not a published
  crate, but the browser build runs at
  [chem.chiral.one](https://chem.chiral.one); see
  [the user guide](crates/chem-app/USER_GUIDE.md),
  [testing notes](crates/chem-app/docs/E2E-TESTING.md) and
  [how it is published](crates/chem-app/docs/DEPLOY.md).

## Working on it

```sh
cargo test --all-features                    # everything
cargo test -p chem --no-default-features     # the lean library, no wgpu or clap
cargo run --release -p chem-app              # the desktop workbench
crates/chem-app/e2e.sh                       # the browser build, served locally
```

`chem`'s default build deliberately excludes the GPU, so `--all-features` covers
one of four configurations. CI runs the lean one separately; run it locally
before changing anything under `crates/chem/src/gpu` or `crates/chem/src/search`.

## Roadmap

| Wave | Release | Delivers | New |
| --- | --- | --- | --- |
| IR and the registry | v0.7.0 ✅ | The data model, the format registry, the oracle harness | 2 |
| Conversion surface | v0.8.0 ✅ | `chem convert`, option bags, `Supplier`/`Writer`, `-L`/`-H` | 9 |
| Records that are not molecules | v0.9.0 ✅ | Binary encodings, and record kinds beyond the molecule: trajectories, volume grids, force-field topologies, meshes, tabular | 18 |
| 3 — crystals and comp-chem input | v0.10.0 | Crystallography and the comp-chem input writers — templated text | ~35 |
| 4 — comp-chem output parsers | v0.11.0 | The log scrapers, and the auto-detecting "Generic Output" dispatcher | ~25 |
| 5 — reactions, biopolymers, 2D drawing | v0.12.0 | Reactions, query molecules, the binary ChemDraw parsers | ~40 |
| 6 — parity certification | v1.0.0 | Multithreaded parse, the published fidelity matrix, fingerprint formats | ~10 |

## Inspiration

[sdfrust](https://github.com/pdbabin/sdfrust) ·
[smilesDrawer](https://github.com/reymond-group/smilesDrawer)  ·
[RDKit](https://github.com/rdkit/rdkit) ·
[OpenBabel](https://github.com/openbabel/openbabel) ·
[gemmi](https://github.com/project-gemmi/gemmi) ·
[Meeko](https://github.com/forlilab/Meeko) ·
[MDAnalysis](https://github.com/MDAnalysis/mdanalysis) ·
[cpptraj](https://github.com/Amber-MD/cpptraj) ·
[trimesh](https://github.com/mikedh/trimesh) ·
[Jmol](https://jmol.sourceforge.net/)

## Licence

MIT.
