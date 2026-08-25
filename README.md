# rust-chem

Cheminformatics in Rust. Two crates in this workspace:

- **[`chem`](crates/chem/)** — the library and the `chem` command-line tool.
  Published to [crates.io](https://crates.io/crates/chem); documentation at
  [docs.rs/chem](https://docs.rs/chem). See [its README](crates/chem/README.md)
  for the API and the feature flags.
- **`chem-app`** — a desktop and browser workbench built on it. Not published;
  see [the user guide](crates/chem-app/USER_GUIDE.md) and
  [testing notes](crates/chem-app/docs/E2E-TESTING.md).

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

## Inspiration

[sdfrust](https://github.com/pdbabin/sdfrust) ·
[smilesDrawer](https://github.com/reymond-group/smilesDrawer)

## Licence

MIT.
