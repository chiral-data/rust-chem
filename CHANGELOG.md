# Changelog

All notable changes to this project are documented here.

## [0.3.0] - 2026-08-17

Generalized cheminformatics workbench UI — `chem_app` moves from a single-purpose fingerprint-search demo toward a general workbench.

### Features

- Multi-dataset Data window: loading a file or the example set adds a named entry to a "Loaded Files" list rather than overwriting the active dataset; click any entry to swap back to it (#62)
- SDF file support alongside SMILES, with multi-record splitting since the parser only handles one `$$$$`-terminated record at a time (#63)
- Per-molecule computed-property table (Name, SMILES, Formula, MW, Fingerprint, Aromatic) replacing the old plain molecule list; clicking a row opens a floating, viewport-capped molecule detail window (#64)
- Aromaticity detection wired into the UI for the first time, as a direct "Detect Aromaticity" action alongside fingerprint computation (#66)
- SMILES writer (`chemio::smiles_writer`) — the reverse of SMILES parsing, adapted from a previously-orphaned, unbuilt algorithm into a real, tested part of `chemio` (#70)

### Bug Fixes

- `chemio::aromaticity` no longer flags fully saturated rings (e.g. cyclohexane) as aromatic — wiring aromaticity detection into the UI surfaced a real correctness bug in the underlying algorithm (#72)

### Chores

- Removed orphaned workspace-root `src/lib.rs`/`smiles_writer.rs`, safe now that the real SMILES-writer content lives in `chemio` (#71)

## [0.2.1] - 2026-08-17

### CI

- Add CI job to build the web (WASM) bundle (72b5448)


### Documentation

- Add user guide (#45) (d4fd5aa)


### Features

- Make dataset file loading platform-agnostic (#43) (cc4e43e)

- Add WASM/browser build target (#44) (f7444ad)

- Async GPU dispatch for full web GPU parity (#47) (631677d)

- Reduce Morgan batch shader to 6 storage-buffer bindings (#51) (81d5f0e)


### Testing

- Add deterministic CPU-path search tests (#46) (2f5d72d)



