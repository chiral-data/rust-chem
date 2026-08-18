# Changelog

All notable changes to this project are documented here.

## [0.4.0] - 2026-08-18

2D structure depiction — `chem_app` draws molecules rather than only describing them in text. Inspired by [smilesDrawer](https://github.com/reymond-group/smilesDrawer), whose option surface informs the renderer's configuration.

### Features

- 2D coordinate storage on `chemcore::Molecule`, with a `Point2`/`BoundingBox` geometry module (#87)
- SDF atom coordinates are preserved instead of being parsed and discarded (#88)
- Structure renderer: bonds by order, atom labels, with bond endpoints inset to the edge of the label (#89)
- Sizing options — fixed `scale`, a target `bond_length`, or fit-to-rect (#90)
- Ring perception (SSSR) in `chemcore`, replacing a private finder that was capped at 7-membered rings; also populates `Bond::in_ring`, which nothing previously set (#91)
- 2D coordinate generation for molecules without any, so SMILES-derived structures can be drawn: ring systems as regular polygons, fused rings across their shared edge, chains at ~120° (#92)
- Per-element colours and light/dark theming (#96)
- Atom display options — `show_carbons`, `explicit_hydrogens`, `atom_visualization`, plus charge and isotope annotations (#97)
- Structures in the query panel and, optionally, per row of the dataset table (#98)

### Bug Fixes

- SMILES ring-closure bonds between aromatic atoms parsed as `Single`, so benzene had five aromatic bonds and one single one. Because Morgan invariants weight bond order, the same molecule fingerprinted differently depending on how its SMILES was written — two spellings of phenol scored 0.27 Tanimoto against each other instead of 1.0 (#95)
- The second line of a ring bond was drawn on an arbitrary side, landing outside the ring as often as inside (#93)
- A lone carbon with no bonds rendered as a blank panel, having neither a label nor a bond to imply a vertex (#97)

### Refactor

- The deferred-work plumbing shared by the GPU-capable operations is extracted into `Task<T>`, removing 10 of `app.rs`'s 23 `cfg(target_arch)` splits. An `Operation` trait was evaluated and deliberately not adopted: with all four operations in view, their signatures share nothing and a trait would add a layer without removing one (#99)

### Known limitations

- No force-directed layout refinement (smilesDrawer's Kamada-Kawai pass), so heavily bridged, spiro and cage systems can still overlap
- No stereochemistry rendering or condensed pseudo-element groups
- Reaction rendering is not implemented (tracked separately)

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



