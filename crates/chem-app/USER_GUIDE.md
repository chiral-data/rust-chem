# Chem Workbench — User Guide

`chem-app` is a cheminformatics workbench built on `chem::core`, `chem::io`, `chem::fp`, and `chem::gpu`: it loads SMILES and SDF datasets, draws structures, computes Morgan fingerprints, detects aromaticity, and searches by Tanimoto similarity. It runs both as a native desktop app and as a browser (WASM) build.

## Running it

### Native

From the repo root:

```bash
cargo run --release -p chem-app
```

### Web (browser)

Requires the `wasm32-unknown-unknown` target and [`trunk`](https://trunkrs.dev/):

```bash
rustup target add wasm32-unknown-unknown
cargo install trunk
```

Then, from `crates/chem-app/`:

```bash
trunk serve --address 0.0.0.0 --port 8080
```

This builds, serves, and live-reloads on file changes at `http://<this-machine>:8080`. For a one-off static build instead (e.g. to check what a real deployment would look like):

```bash
trunk build --release
cd dist && python3 -m http.server 8080 --bind 0.0.0.0
```

The web build starts up on the CPU and upgrades to GPU acceleration via WebGPU a moment later if the browser supports it (see the menu bar's **🚀 GPU** / **💻 CPU** indicator). Browsers without WebGPU, or whose adapter can't fit the fingerprinting shader's buffer requirements, stay on CPU automatically — no user action needed either way.

## Using it

A menu bar sits above a workspace that three floating windows sit on:

- **Datasets** — the datasets you have loaded, and the active one's molecules
- **Operations** — everything you run against the dataset
- **Inspector** — the query and ranked results
- **Settings** — preferences: theme, and how structures are drawn (closed until you open it)

Each is movable and resizable, and each can be closed from its own **✕** or toggled from the **View** menu. Close the last one and the empty workspace tells you where to get them back. Window content is still being reorganised across v0.5.0, so a control described under one window may move to a neighbouring one in a later release.

### The menu bar

- **File** — Load File, Load Examples (the same actions as the buttons in Datasets).
- **View** — one checkbox per window, **Reset layout**, and **Close all molecule windows**.
- **Settings** — opens and closes the Settings window. It stays lit while the window is open.
- **Right-hand side** — current FPS, and the **GPU** / **CPU** chips described below.

### 1. Load a dataset

In the **Datasets** window, or from the **File** menu:

- **📋 Load Examples** — loads 15 built-in molecules (methane, benzene, phenol, aniline, aspirin-adjacent structures, etc.) instantly. Good default if you just want to try things out.
- **📂 Load File** — pick any format the library reads, from disk. As of v0.8.0 that is eleven: SMILES (`.smi`/`.smiles`/`.txt`), SDF (`.sdf`), CXSMILES (`.cxsmiles`), XYZ (`.xyz`), PDB (`.pdb`/`.ent`), mmCIF (`.cif`/`.mmcif`), Mol2 (`.mol2`), PDBQT (`.pdbqt`), GRO (`.gro`), CML (`.cml`) and commonchem JSON (`.json`). The dialog's list is generated from the format registry, so it stays in step with what the library supports.
  - The file's extension decides how it is read; an unrecognised one is treated as SMILES.
  - Formats that carry structures rather than SMILES strings — PDB, mmCIF, Mol2 and the rest — show their own name in the SMILES column, e.g. `(PDB)`.
  - A structure file's 3D coordinates are not a 2D drawing: the depiction is laid out from the connectivity rather than flattened from the conformer, so it shows what is bonded to what and not the real geometry.
  - SMILES format is one molecule per line: `SMILES [optional name]`. Lines starting with `#` are comments and blank lines are skipped. If no name is given, molecules are auto-named `Molecule_<line number>`.
  - SDF files can hold multiple `$$$$`-terminated molecule records; each is parsed independently, using the record's own name field if present.
  - Each load adds a new entry to the **Files** list rather than replacing what's already there — click any entry to switch back to it. Each is shown with its format and molecule count, so two SMILES files are told apart without switching between them. Loading a file with the same name as an existing entry (e.g. reloading the same path) updates that entry in place instead of adding a duplicate.
  - The **✕** beside an entry removes one you've finished with. The last one can't be removed, since the app always has a dataset. Removing the one you're looking at discards its fingerprints and search results and switches to a neighbour; removing any other leaves your work untouched.

The table below lists the active dataset's molecules and scrolls through all of them, however many there are — only the rows on screen are drawn. Click a name to open that molecule's detail window; click it again to close it. A row stays highlighted while its window is open.

Each detail window has an **Export SVG** button, which saves that molecule's structure as a vector file — a file on the desktop build, a download in the browser. Exports always use the light palette, whatever theme the app is showing, since a structure bound for a printed page shouldn't carry a dark background's colours.

Several detail windows can be open at once, so two molecules can be compared side by side rather than remembered. Eight is the limit: opening a ninth closes the oldest, since rows can be clicked much faster than windows can be closed. **View → Close all molecule windows** clears them.

**Show structures in table** adds a structure column. Every molecule is drawn: a file that brought its own layout is shown as the file drew it, and anything else is laid out on demand. That layout is shared with the detail windows and the search results, so one molecule looks the same everywhere it appears.

### 2. Run an operation

The **Operations** window holds everything that computes, one collapsing section each. Each section's header carries what happened last time it ran — what it produced, how long it took, and whether GPU or CPU did it — so a collapsed section still reports itself. Switching datasets clears those, since they describe data that is no longer active.

**Backend**, at the top, picks **🚀 GPU** or **💻 CPU** for the operations that can use either (fingerprints and search). If GPU initialisation failed, the reason is shown here in full, with a **Retry**.

**Fingerprints**

- **Radius** — Morgan fingerprint radius (0–5). Higher values capture larger structural neighborhoods around each atom.
- **Size** — fingerprint length in bits (512–4096, logarithmic slider).
- **⚡ Compute Fingerprints** — generates a fingerprint for every molecule in the loaded dataset.

You need to do this at least once before you can search.

**Aromaticity**

- **🔬 Detect Aromaticity** — runs ring perception across the dataset and flags aromatic atoms, which the dataset table's *Aromatic* column then reflects.

**2D Coordinates**

- **📐 Generate Coordinates** — lays out every molecule that doesn't already have coordinates, so it can be drawn. Molecules whose coordinates came from an SDF file keep them, and the section says how many it generated against how many it kept. Structures are also laid out on demand when you open one, so this is for doing the whole dataset at once.

**Convert**

- **Write as** — picks the output format, from everything this build can write.
- **What it will cost, before you run it.** Formats hold different things: XYZ has no bond block, SMILES has no coordinates, PDB has no isotopes. The section lists what this dataset would lose to the format you picked, and how many molecules lose each thing, as soon as you pick it. A conversion that keeps everything says so.
- **⟳ Convert** — writes the dataset in the chosen format, reads it back, and adds *that* as a new dataset, then switches to it. The one you converted from stays in the Files list, so you can click between the two and see what changed. Converting benzene from a CML file to SMILES, for instance, shows `C1CCCCC1` where the original had `c1ccccc1` — cyclohexane, drawn without the aromatic ring. The report predicts that; the new dataset is it.
- Converting again to the same format replaces the earlier result rather than adding another entry.
- The report accounts for losses that come from the *pair* of formats rather than the target alone — a few conversions lose something both formats otherwise carry, and the reason is named when so. `chem convert` on the command line reports the same thing, from the same code.

**Export**

- **💾 Export…** — writes the active dataset in its own format: a file on the desktop build, a download in the browser. Convert first if you want a different format.

### 3. Look at the results

The **Inspector** window has two sections.

**Query** — the parsed query molecule drawn, its details, and its fingerprint as a bit grid. **Export SVG** saves the structure as a vector file you can drop into a document or a slide. It is labelled with the SMILES it was parsed from, which is not necessarily what is currently in the box in Operations.

**Results** — each hit drawn, with its rank, name, SMILES and similarity score. Seeing the molecule is the point: two structures can score 0.9 for reasons obvious in a drawing and invisible in a SMILES string.

A structure always appears, laid out on demand where the file carried no coordinates of its own.

Click **▼ Why?** on a result to see how the score was arrived at. One grid holds both fingerprints, each bit coloured by which of them has it set: **blue** for bits in both, **amber** for bits only in this molecule, **violet** for bits only in your query, and background for bits in neither.

That is the score, drawn. Tanimoto similarity is the count of shared bits divided by the count set in either, so the blue cells are the numerator and blue plus amber plus violet is the denominator — both printed above the grid. For a molecule's atoms and bonds, click its name in the Datasets table to open a detail window.

### What is remembered

The workspace comes back the way you left it. Window positions and sizes, which windows were open, the theme, the structure display options, and the fingerprint radius, size and top-k all survive a restart — on the desktop build, and on web across a reload of the same browser.

What is *not* remembered is the data: loaded files, computed fingerprints, search results, the query text, and which molecule windows were open. Those are re-made each session rather than restored, so a fingerprint is never shown that was computed under settings you can no longer see.

**View → Reset layout** puts the windows back where a first launch has them, which is the way out if one has been dragged somewhere unreachable.

### Settings

Open it with **Settings** in the menu bar, or from the **View** menu. It stays open while you work, so a change can be watched taking effect rather than applied blind.

- **Theme** — Light, Dark, or follow the system. Structure colours follow it: the same molecule is drawn with a light or dark palette to match.
- **Structures** — which carbons are labelled, how atoms are annotated, whether hydrogens are explicit. These apply to every structure the app draws, which is why they are here rather than beside any one of them.
- **Show structures in the dataset table** — the thumbnail column.

### Backend chips

Top-right of the menu bar shows the current FPS and a **CPU** / **GPU** selector — click either one to switch, from anywhere, without opening a window. The one currently in use is highlighted. The Operations window's **Backend** section is the same setting with more detail.

- **💻 CPU** — always available; click it any time to force CPU-only fingerprinting/search.
- **🚀 GPU** (green) — GPU is available. Click it to switch to (or back to) GPU acceleration — switching back after having used it once is instant, no re-init needed.
- **⚠ GPU** (red) — a GPU init attempt actually failed. Hover it to see why (e.g. "GPU unavailable: Morgan: No suitable GPU adapter found" or a storage-buffer-limit message). Click it to retry — useful if you've since enabled WebGPU in your browser, for instance.

## Known limitations (web build)

- GPU acceleration on web depends on the browser and its adapter — Chrome/Edge with WebGPU enabled are the best bet today; browsers without WebGPU, or whose adapter reports too few storage-buffer bindings for the fingerprinting shader, run on CPU instead. Either way the app falls back automatically rather than failing.
- Search re-uploads the target dataset to the GPU on every query in the web build, rather than reusing a cached upload the way native does — imperceptible at demo-scale dataset sizes, but a difference worth knowing about.
- Very large datasets may be slower to fingerprint on CPU than they would be on a native GPU-enabled build.
