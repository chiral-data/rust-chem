# ChemFP Demo — User Guide

`chem_app` is a molecular fingerprint search demo built on `chemcore`, `chemio`, `chemfp`, and `chemgpu`. It runs both as a native desktop app and as a browser (WASM) build.

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

The web build starts up on the CPU and upgrades to GPU acceleration via WebGPU a moment later if the browser supports it (see the status bar's **🚀 GPU** / **💻 CPU** indicator). Browsers without WebGPU, or whose adapter can't fit the fingerprinting shader's buffer requirements, stay on CPU automatically — no user action needed either way.

## Using the demo

The window is split into four panels: **Dataset** (left), **Query Molecule** (right), **Search Results** (center), and a status bar (top).

### 1. Load a dataset

In the **Dataset** panel:

- **📋 Load Examples** — loads 15 built-in molecules (methane, benzene, phenol, aniline, aspirin-adjacent structures, etc.) instantly. Good default if you just want to try things out.
- **📂 Load File** — pick a `.smi` / `.smiles` / `.txt` file from disk. Format is one molecule per line: `SMILES [optional name]`. Lines starting with `#` are comments and blank lines are skipped. If no name is given, molecules are auto-named `Molecule_<line number>`.

### 2. Compute fingerprints

Still in the **Dataset** panel:

- **Radius** — Morgan fingerprint radius (0–5). Higher values capture larger structural neighborhoods around each atom.
- **Size** — fingerprint length in bits (512–4096, logarithmic slider).
- **⚡ Compute Fingerprints** — generates a fingerprint for every molecule in the loaded dataset. The status line reports how long it took and whether it ran on GPU or CPU.

You need to do this at least once before you can search.

### 3. Enter a query molecule

In the **Query Molecule** panel:

- Type a SMILES string in the text box, e.g. `c1ccccc1O` (phenol) or `CC(=O)Oc1ccccc1C(=O)O` (aspirin-like).
- Click **Parse**, or just stop typing — it auto-parses after a short idle delay (debounced so it doesn't re-parse on every keystroke).
- On success, you'll see the parsed molecule's info and its fingerprint rendered as a bit grid.
- Invalid SMILES shows an error message instead.

### 4. Search

- Set **Top K** — how many ranked results to return.
- Click **🔍 Search** (enabled once you have both a parsed query and computed dataset fingerprints). Results are ranked by Tanimoto similarity.

### 5. Explore results

In the **Search Results** panel, each hit shows its rank, structure summary, and similarity score. Click **▼ Show** on any result to expand it in place: atom list, bond list, and a side-by-side fingerprint comparison against your query molecule. Click **▲ Hide** to collapse it again.

### Status bar

Top-right of the window shows the current FPS and whether fingerprint search is running on **🚀 GPU** or **💻 CPU**, plus a button next to it:

- **🚀 GPU** (green) — GPU is active. Click **Use CPU** to switch to the CPU path without losing the GPU context (switching back is instant, no re-init needed).
- **💻 CPU** (yellow) — on CPU by choice (either GPU hasn't been tried yet, or you switched away from it manually). The button reads **Use GPU** if a GPU context is already available (instant switch back) or **Try GPU** if one needs to be initialized first.
- **⚠ CPU** (red) — a GPU init attempt actually failed. Hover the indicator to see why (e.g. "GPU unavailable: Morgan: No suitable GPU adapter found" or a storage-buffer-limit message). Click **Retry GPU** to try again — useful if you've since enabled WebGPU in your browser, for instance.

## Known limitations (web build)

- GPU acceleration on web depends on the browser and its adapter — Chrome/Edge with WebGPU enabled are the best bet today; browsers without WebGPU, or whose adapter reports too few storage-buffer bindings for the fingerprinting shader, run on CPU instead. Either way the app falls back automatically rather than failing.
- Search re-uploads the target dataset to the GPU on every query in the web build, rather than reusing a cached upload the way native does — imperceptible at demo-scale dataset sizes, but a difference worth knowing about.
- Very large datasets may be slower to fingerprint on CPU than they would be on a native GPU-enabled build.
