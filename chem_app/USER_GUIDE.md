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

- **Data Sources** — loaded files and the dataset table
- **Operations** — everything you run against the dataset
- **Data Visualization** — ranked search results

Each is movable and resizable, and each can be closed from its own **✕** or toggled from the **View** menu. Close the last one and the empty workspace tells you where to get them back. Window content is still being reorganised across v0.5.0, so a control described under one window may move to a neighbouring one in a later release.

### The menu bar

- **File** — Load File, Load Examples (the same actions as the buttons in Data Sources).
- **View** — one checkbox per window.
- **Right-hand side** — current FPS, and the **GPU** / **CPU** chips described below.

### 1. Load a dataset

In the **Data Sources** window, or from the **File** menu:

- **📋 Load Examples** — loads 15 built-in molecules (methane, benzene, phenol, aniline, aspirin-adjacent structures, etc.) instantly. Good default if you just want to try things out.
- **📂 Load File** — pick a `.smi` / `.smiles` / `.txt` file, or a `.sdf` file, from disk.
  - SMILES format is one molecule per line: `SMILES [optional name]`. Lines starting with `#` are comments and blank lines are skipped. If no name is given, molecules are auto-named `Molecule_<line number>`.
  - SDF files can hold multiple `$$$$`-terminated molecule records; each is parsed independently, using the record's own name field if present.
  - Each load adds a new entry to the **Loaded Files** list rather than replacing what's already there — click any entry to switch back to it. Loading a file with the same name as an existing entry (e.g. reloading the same path) updates that entry in place instead of adding a duplicate.

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
- On success you'll see the parsed molecule drawn, its info, and its fingerprint as a bit grid. Invalid SMILES shows an error instead.
- **Top K** — how many ranked results to return.
- **🔍 Search** — ranked by Tanimoto similarity. It needs both a parsed query and computed dataset fingerprints; if either is missing, the section says which.

### 3. Explore results

In the **Data Visualization** window, each hit shows its rank, structure summary, and similarity score. Click **▼ Show** on any result to expand it in place: atom list, bond list, and a side-by-side fingerprint comparison against your query molecule. Click **▲ Hide** to collapse it again.

### Backend chips

Top-right of the menu bar shows the current FPS and a **CPU** / **GPU** selector — click either one to switch, from anywhere, without opening a window. The one currently in use is highlighted. The Operations window's **Backend** section is the same setting with more detail.

- **💻 CPU** — always available; click it any time to force CPU-only fingerprinting/search.
- **🚀 GPU** (green) — GPU is available. Click it to switch to (or back to) GPU acceleration — switching back after having used it once is instant, no re-init needed.
- **⚠ GPU** (red) — a GPU init attempt actually failed. Hover it to see why (e.g. "GPU unavailable: Morgan: No suitable GPU adapter found" or a storage-buffer-limit message). Click it to retry — useful if you've since enabled WebGPU in your browser, for instance.

## Known limitations (web build)

- GPU acceleration on web depends on the browser and its adapter — Chrome/Edge with WebGPU enabled are the best bet today; browsers without WebGPU, or whose adapter reports too few storage-buffer bindings for the fingerprinting shader, run on CPU instead. Either way the app falls back automatically rather than failing.
- Search re-uploads the target dataset to the GPU on every query in the web build, rather than reusing a cached upload the way native does — imperceptible at demo-scale dataset sizes, but a difference worth knowing about.
- Very large datasets may be slower to fingerprint on CPU than they would be on a native GPU-enabled build.
