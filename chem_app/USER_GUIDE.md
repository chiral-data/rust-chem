# Chem Workbench — User Guide

`chem_app` is a cheminformatics workbench built on `chemcore`, `chemio`, `chemfp`, and `chemgpu`: it loads SMILES and SDF datasets, draws structures, computes Morgan fingerprints, detects aromaticity, and searches by Tanimoto similarity. It runs both as a native desktop app and as a browser (WASM) build.

## Running it

### Native

From the repo root:

```bash
cargo run --release -p chem_app
```

### Web (browser)

Requires the `wasm32-unknown-unknown` target and [`trunk`](https://trunkrs.dev/):

```bash
rustup target add wasm32-unknown-unknown
cargo install trunk
```

Then, from `chem_app/`:

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
- **View** — one checkbox per window, and **Close all molecule windows**.
- **⚙** — opens and closes **Settings**, at the right-hand end.
- **Right-hand side** — current FPS, and the **GPU** / **CPU** chips described below.

### 1. Load a dataset

In the **Datasets** window, or from the **File** menu:

- **📋 Load Examples** — loads 15 built-in molecules (methane, benzene, phenol, aniline, aspirin-adjacent structures, etc.) instantly. Good default if you just want to try things out.
- **📂 Load File** — pick a `.smi` / `.smiles` / `.txt` file, or a `.sdf` file, from disk.
  - SMILES format is one molecule per line: `SMILES [optional name]`. Lines starting with `#` are comments and blank lines are skipped. If no name is given, molecules are auto-named `Molecule_<line number>`.
  - SDF files can hold multiple `$$$$`-terminated molecule records; each is parsed independently, using the record's own name field if present.
  - Each load adds a new entry to the **Files** list rather than replacing what's already there — click any entry to switch back to it. Each is shown with its format and molecule count, so two SMILES files are told apart without switching between them. Loading a file with the same name as an existing entry (e.g. reloading the same path) updates that entry in place instead of adding a duplicate.

The table below lists the active dataset's molecules and scrolls through all of them, however many there are — only the rows on screen are drawn. Click a name to open that molecule's detail window; click it again to close it. A row stays highlighted while its window is open.

Several detail windows can be open at once, so two molecules can be compared side by side rather than remembered. Eight is the limit: opening a ninth closes the oldest, since rows can be clicked much faster than windows can be closed. **View → Close all molecule windows** clears them.

**Show structures in table** adds a structure column. A structure is drawn where the molecule has coordinates: SDF files bring their own, and for anything parsed from SMILES, run **2D Coordinates** in Operations to generate them. Rows without coordinates show a dash rather than a blank cell.

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

**Similarity Search**

- Type a SMILES string in the text box, e.g. `c1ccccc1O` (phenol) or `CC(=O)Oc1ccccc1C(=O)O` (aspirin-like).
- Click **Parse**, or just stop typing — it auto-parses after a short idle delay (debounced so it doesn't re-parse on every keystroke).
- On success the parsed molecule is drawn in the **Inspector** window's *Query* section, with its details and its fingerprint. Invalid SMILES shows an error here, beside the box you typed it in.
- **Top K** — how many ranked results to return.
- **🔍 Search** — ranked by Tanimoto similarity. It needs both a parsed query and computed dataset fingerprints; if either is missing, the section says which.

### 3. Look at the results

The **Inspector** window has two sections.

**Query** — the parsed query molecule drawn, its details, and its fingerprint as a bit grid. It is labelled with the SMILES it was parsed from, which is not necessarily what is currently in the box in Operations.

**Results** — each hit with its rank, summary and similarity score. Click **▼ Show** on any result to expand it in place: atom list, bond list, and a side-by-side fingerprint comparison against your query. Click **▲ Hide** to collapse it again.

### Settings

Open it with the **⚙** at the right of the menu bar, or from the **View** menu. It stays open while you work, so a change can be watched taking effect rather than applied blind.

- **Theme** — Light, Dark, or follow the system. Structure colours follow it: the same molecule is drawn with a light or dark palette to match.
- **Structures** — which carbons are labelled, how atoms are annotated, whether hydrogens are explicit. These apply to every structure the app draws, which is why they are here rather than beside any one of them.
- **Show structures in the dataset table** — the thumbnail column. A structure is drawn where the molecule has coordinates; run **2D Coordinates** in Operations for anything parsed from SMILES.

### Backend chips

Top-right of the menu bar shows the current FPS and a **CPU** / **GPU** selector — click either one to switch, from anywhere, without opening a window. The one currently in use is highlighted. The Operations window's **Backend** section is the same setting with more detail.

- **💻 CPU** — always available; click it any time to force CPU-only fingerprinting/search.
- **🚀 GPU** (green) — GPU is available. Click it to switch to (or back to) GPU acceleration — switching back after having used it once is instant, no re-init needed.
- **⚠ GPU** (red) — a GPU init attempt actually failed. Hover it to see why (e.g. "GPU unavailable: Morgan: No suitable GPU adapter found" or a storage-buffer-limit message). Click it to retry — useful if you've since enabled WebGPU in your browser, for instance.

## Known limitations (web build)

- GPU acceleration on web depends on the browser and its adapter — Chrome/Edge with WebGPU enabled are the best bet today; browsers without WebGPU, or whose adapter reports too few storage-buffer bindings for the fingerprinting shader, run on CPU instead. Either way the app falls back automatically rather than failing.
- Search re-uploads the target dataset to the GPU on every query in the web build, rather than reusing a cached upload the way native does — imperceptible at demo-scale dataset sizes, but a difference worth knowing about.
- Very large datasets may be slower to fingerprint on CPU than they would be on a native GPU-enabled build.
