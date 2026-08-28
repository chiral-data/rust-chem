# Changelog

All notable changes to this project are documented here.

## [0.6.0] - 2026-08-28

**`chem` is on crates.io.** Everything the workbench could do is now available as a library and a command-line tool: `chem = "0.6"` for the library, `cargo install chem --features cli` for the tool.

The libraries were seven crates. They are one. crates.io has a flat namespace, so seven crates meant claiming seven names permanently for one project — and nothing had been published yet, so the choice was still free. One crate is the reversible direction: splitting later is easy, un-publishing seven abandoned names is not. `chem::core`, `chem::io`, `chem::fp`, `chem::gpu`, `chem::draw`, `chem::search`.

### Features

- **A `chem` command-line tool** with six subcommands over the same operations the workbench runs: `info`, `fp`, `search`, `aromatic`, `coords`, `draw`. It reads a named file or standard input, writes to `-o` or standard output, and puts progress and warnings on standard error — so `chem aromatic x.smi | chem fp` composes. Exit codes distinguish an unusable file from an empty result (#139, #140, #141, #142)
- **Depiction became a library.** `describe_structure`, the shape types and the SVG writer left the GUI crate, so anything can draw a molecule without linking a window and a GPU surface (#138)
- **File-level reading became a library.** Splitting a `.smi` by line and an `.sdf` on `$$$$` lived in the app; now both front ends agree, by construction, on what a file contains. Records that fail come back to the caller rather than into a log, which is what lets a batch tool count them and exit on them (#139)
- **An SDF writer.** There was none, so 2D coordinates had nowhere to go — SMILES cannot carry them, and `chem coords` would have discarded the thing it just computed (#161)
- **Backend selection became a library**, so the CLI reaches the GPU kernels with no window server (#140)
- The GPU is behind an off-by-default `gpu` feature: **32 dependencies without it, 177 with**. `wgpu` is irrelevant to anyone parsing SMILES or fingerprinting on a CPU, which is what every GPU-less machine and every WebAssembly build already do (#144)
- Six runnable examples, rendered by docs.rs (#145)

### Bug Fixes

Four in the SMILES parser, all found by pointing a batch tool at whole files rather than one molecule at a time.

- **An unspecified bond to an aromatic atom parsed as aromatic even when the other atom was not.** The default consulted only the atom being added, so the answer depended on which end the author wrote first: `Oc1ccccc1` got an aromatic C–O bond and `c1ccccc1O` did not, and the two spellings of phenol scored 0.467 against each other. **32 of the 122 molecules in `test.smi` were affected** — every toluene, cresol, xylene and aniline in the file had a wrong fingerprint (#167)
- **Aromaticity detection missed pyridine-type nitrogen entirely**, returning zero aromatic atoms for pyridine, pyrimidine, pyrazine and imidazole, and 6 of 10 for quinoline. The π-electron count was inverted and keyed off atom degree, which cannot distinguish the two nitrogen types. A single pass also under-detected fused rings order-dependently — naphthalene needed two passes while anthracene needed one, so detection depended on how the SMILES was written (#160)
- **Back-to-back ring-closure digits were read as one number.** A bare ring label is exactly one digit, so `c12` closes rings 1 and 2 — but `digit1` consumed both, leaving them open. Anthracene and quinoline failed to parse at all while their isomers succeeded, purely because of where the digits landed (#69)
- **A string of bond characters parsed as a molecule.** Every bond character is legal in isolation, so `$$$$` — the SDF record terminator — built an atomless molecule reported as success. Reading an SDF as SMILES exited 0 claiming to have read molecules (#151)

### Corrections to v0.5.0 and earlier

- `chemfp` shipped three binaries from a library crate, a dead 79-line `BitVec` nothing imported, a commented-out copy of its own module list, and a `rayon` dependency used only in a doc comment
- `gpu::init_logging()` called `env_logger::builder().init()`, which panics if a logger already exists — a landmine in a library, since the caller cannot know whether something got there first
- `crates/chem-app/README` and the root README described the repository around v0.1.0 and cited paths that had moved twice

### Testing

278 tests at the start of this milestone, **353** at the end. Three findings worth recording:

- Two tests were written to match the behaviour rather than the specification. `test_parse_ring_number` had a comment reading "only takes first digit" above an assertion expecting `123`; the polycyclic aromaticity test used lowercase SMILES and admitted in its own comment that the parser had already set every flag before detection ran — it would have passed with the function deleted
- `build_molecule` produced four separate bugs this milestone, so the aromatic-bond test asserts the *invariant* — no aromatic bond may have a non-aromatic endpoint — rather than a fifth list of cases
- Asserting on fingerprints from inside the SMILES parser's own tests only became possible once the crates merged

### Known limitations

- Aromaticity perception is a Hückel heuristic over the SSSR. Exotic systems and charged aromatics are untested
- SDF output writes atoms, bonds and 2D coordinates; charges, isotopes, chirality and bond stereo are not yet written, and aromatic bonds are written as type 4, a query bond type in strict molfile
- The GPU path panics on a batch large enough to outgrow `max_storage_buffer_binding_size` (#158), and `--backend auto` prefers the GPU even for a one-shot search where it is slower (#159)

## [0.5.0] - 2026-08-22

Floating-window workbench UI — `chem_app` moves from four fixed panels to a menu bar over a free canvas, with **Datasets**, **Operations**, **Inspector** and **Settings** as windows whose arrangement survives a restart. Plus the visualization work deferred from v0.4.0: SVG export, layout refinement for bridged systems, and structures in search results.

The app is now called **Chem Workbench**. "ChemFP Demo — Molecular Fingerprint Search" described what v0.2.0 was.

### Features

- Menu bar over an empty workspace, with a registry holding which windows exist and whether they are open; the View menu is generated from it (#115)
- **Datasets** window: a virtualised molecule table, so a ten-thousand-molecule dataset scrolls rather than being silently cut to its first twenty rows, with resizable columns and each loaded file labelled with its format and molecule count (#124)
- Loaded datasets can be removed. Removing one you are not looking at leaves its neighbour's fingerprints and results alone; removing the active one discards what belonged to it (#131)
- **Operations** window: the four operations as collapsing sections, each reporting what it produced, how long it took and which backend ran it. 2D coordinate generation becomes an action for the first time, having previously only happened as a side effect of opening a structure (#116)
- **Inspector** window: the parsed query drawn with its fingerprint, and the ranked results (#119)
- Search results draw the molecule they found, and expanding one shows the score's arithmetic — shared bits, bits in either, and their division — over a single grid coloured by which fingerprint holds each bit (#130)
- **Settings** window, reachable from the menu bar, holding a light/dark/system theme picker the app never had, and the structure display options that govern every structure it draws (#126)
- Several molecule detail windows open at once, so two molecules can be compared rather than remembered. Eight at a time; the ninth closes the oldest (#125)
- Window geometry, which windows were open, the theme, the display options and the operation parameters survive a restart, and a reload on web. **View → Reset layout** recovers a window dragged somewhere unreachable, which is newly necessary now that the position persists (#128)
- Structures export as SVG, from any detail window or the query (#132)
- Layout refinement for bridged and cage systems, by stress majorization on the Kamada-Kawai energy (#129)

### Bug Fixes

- The dataset table widened its own window until it covered its neighbour instead of scrolling. `ScrollArea::vertical` leaves the horizontal axis disabled, and egui sizes a disabled axis to its content, so a seven-column table stretched its container without limit (#115)
- The two fingerprints in a search result could not be compared: grid width was a per-call-site argument, so the same 2048 bits formed a landscape grid in one view and a portrait one in another, with no bit in the same position in both — and two side by side overflowed the window (#130)
- Unset fingerprint cells were a hardcoded near-white: invisible on a light background and glaring on a dark one. Harmless until this release made the theme a choice (#130)
- `StructureOptions::atom_visualization` had been renamed to `atom_inspector` by a text replace made for an unrelated window, leaving a field name related to neither its type nor its meaning (#132)
- The dataset remove control used `✕` (U+2715), which is not in egui's bundled font subset and drew as a missing-glyph box. It is painted from two line segments, as egui does for its own window close button (#131)

### Refactor

- `describe_structure` separates what to draw from painting it, which is what made SVG export possible — the renderer previously computed geometry and emitted egui shapes in one pass, with nothing in between to serialise. The one part a backend supplies is measuring a label, and that is the point rather than the inconvenience: bonds inset to the edge of their labels, so insets measured with one font and drawn with another leave bonds striking through text. It also gave the renderer its first tests of what it assembles rather than of its helpers (#132)
- The app's single 30-field struct is split into shared state and per-view state, and each panel moved to its own file. Every dataset-changing path previously ended in the same six-line reset, duplicated three times, half of it clearing state the shared half can no longer reach; an epoch counter replaces it, and a view holding stale indices notices for itself (#113)
- Each operation reports its own outcome. One shared status string served fourteen writers and one reader, so running aromaticity wiped the fingerprint result off the screen (#116)
- A uniform `dyn Window` trait for the registry was evaluated and deliberately not adopted: it would unify four call sites, the molecule detail windows cannot join it, and the menu bar needs named access to one window. A layer added and none removed (#127)

### Documentation

- `chem_app/e2e.sh` builds and serves the web build for hands-on testing, refusing to serve a bundle older than the run and stamping the commit into the page — CI cannot see what the app does, and a stale bundle once made a change look ineffective across eleven minutes of commits. `docs/E2E-TESTING.md` lists what to check per story (#117)

### Testing

- `cargo test --workspace` no longer deadlocks on a machine with a GPU. Several tests each requested their own adapter and device, and concurrent creation against one physical GPU hangs in the driver — test binaries stayed alive for half an hour holding GPU handles. Two crates were affected: `chemgpu`'s Morgan tests and `chem_app`'s search tests now take a lazily-initialised shared instance, the pattern [#19](https://github.com/chiral-data/rust-chem/issues/19) established for the rest of `chemgpu`. The suite goes from hanging indefinitely to 3.6 seconds (#135)

### CI

- The pipeline runs on pull requests into milestone branches, not only `main`. It previously ran once per milestone, so a break could sit unattributed behind a dozen merged PRs. Linux only off the path to `main`, where the three-OS sweep belongs (#114)

### Corrections to v0.4.0

- v0.4.0's known limitations say bridged, spiro and cage systems "can still overlap". **They do not overlap.** Measurement found no clashing atoms in any of them, and spiro and steroid systems were already laid out exactly. The defect was stretched bonds, in bridged and cage systems only — adamantane drew a bond 4.4× its target length, rendering as a long spurious line across the molecule. That is what #129 fixes.

### Known limitations

- A cage cannot have uniform bond lengths *and* no crossings in two dimensions, so refinement is a trade rather than a fix. It runs only where bond lengths are bad, and is kept only where it improves them without introducing an overlap — so a cage comes out compact and tight rather than spacious.
- An atom closing two rings with consecutive bare digits (`C34`) fails to parse ([#69](https://github.com/chiral-data/rust-chem/issues/69)). That is how any polycyclic where three or more rings meet at a shared atom is written, so **steroids cannot be loaded from SMILES at all** — only from SDF, where connectivity is explicit.
- An exported SVG is internally consistent but not pixel-identical to the screen: it estimates metrics for the font it declares, so a viewer substituting a metrically different face shows slight drift. Converting glyphs to outlines would remove it, at the cost of a font-parsing dependency and labels that are no longer text.
- Persistence covers layout and preferences, not data. Loaded datasets, computed fingerprints, search results and open molecule windows are not restored — a fingerprint cache restored under settings no longer visible is worse than none.
- Reaction rendering is still not implemented ([#85](https://github.com/chiral-data/rust-chem/issues/85)), and stereochemistry and condensed pseudo-element groups remain unaddressed.

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



