# Testing chem_app by hand

CI proves the code compiles on both targets, that the tests pass, that clippy is
quiet at `-D warnings`, that rustdoc has no broken links, and that the web bundle
builds. It proves nothing whatsoever about what the app *does*, because there is
no display on a runner and no one looking at it.

So everything below has to be checked by a person. The list is not
box-ticking: each entry exists because it is a behaviour a test could not have
caught, and several of them are behaviours that were broken at some point during
v0.5.0 and only found by looking.

## Running it

```bash
chem_app/e2e.sh                        # builds, then serves on 127.0.0.1:8080
chem_app/e2e.sh --port 8081
chem_app/e2e.sh --address 0.0.0.0      # reachable from another machine
```

Native, where there is a display:

```bash
cargo run --release -p chem_app
```

Use `--release` for anything involving timings or the FPS readout. A debug wasm
build makes fingerprint generation and layout slow enough that the numbers the
app reports are meaningless.

## First: confirm what you are testing

The app prints its build id at the top right, next to the FPS readout, and
`e2e.sh` prints the same string when it starts. **If they disagree, stop** — the
browser is showing you a cached page. Hard-refresh (Ctrl+Shift+R).

This is not hypothetical. During v0.5.0 an interrupted build left an older
bundle in place, and the preview looked unchanged across eleven minutes of
commits; the conclusion drawn from it ("the change had no effect") was wrong.
`e2e.sh` now rebuilds from scratch, refuses to serve artefacts older than the
run, and refuses to start if another server holds the port — but the browser's
cache is outside its reach, which is what the id in the corner is for.

A build id ending `-dirty` means the tree had uncommitted changes, so the bundle
matches no commit and "it works on `abc1234`" is not reproducible. Fine while
iterating, not fine as evidence.

## What to check, by story

### Windows and the shell (#103)

- **View menu** toggles each of the three windows; unchecking closes, rechecking
  reopens where it was.
- **Close all three.** The canvas should say so and point at the View menu. A
  blank canvas with no way back is the failure this guards against.
- **Default layout**: three windows, none overlapping, visible canvas between
  and around them. Check on a *narrow* window too — the layout is derived from
  the viewport, so resize the browser and reload.
- **Drag a window mostly off-screen.** It should be held back. There is no
  "reset layout" yet (#108), so a window that escapes is unrecoverable.
- **Narrow Datasets until the table won't fit.** The table must scroll
  sideways *inside* the window. If the window widens to fit the table instead,
  that is the `ScrollArea` regression from #115 returning: `ScrollArea::vertical`
  sizes a disabled axis to its content, and only `both()` clamps it. **This one
  cannot be tested headlessly** — it is egui runtime layout behaviour.
- Title bar and browser tab both read **Chem Workbench**.

### Operations (#105)

- **Each section's header reports its last run** — what it produced, how long,
  which backend. Collapse a section and the header still says it.
- **Run Fingerprints, then Aromaticity.** The fingerprint result must still be
  there. It used to be wiped: one shared status string, last writer wins.
- **Run 2D Coordinates twice.** First run generates; the second should report
  everything as *kept from file*, which is the same branch that preserves
  coordinates an SDF supplied (#88).
- **Press Search with nothing set up.** It should name the missing prerequisite,
  not just grey out.
- **Backend radio and the menu bar chips must agree**, in both directions.
- **The query debounce survives a closed window.** Type a SMILES, close the
  Operations window inside 300ms, wait, reopen it. The query should have parsed.
  The debounce ticks from the frame loop precisely so that closing its window
  doesn't stop it, and nothing else verifies that.

### Datasets (#104)

- **Load the examples, then scroll the table.** All 15 rows reachable. The table
  used to stop at 20 whatever the dataset held, so on a large file the check is
  that scrolling keeps going — only the visible rows are drawn, and a bug there
  shows as blank or repeated rows while scrolling fast.
- **Each file in the list shows its format and molecule count.** Load a `.smi`
  and a `.sdf` and confirm both are labelled correctly.
- **Drag a column edge.** Columns are resizable now; SMILES is the one worth
  widening.
- **Turn on "Show structures in table" with a SMILES dataset.** Cells show a dash,
  not a structure — coordinates don't exist yet. Run **2D Coordinates** in
  Operations and they appear. With an SDF file they should be there immediately,
  since the file supplies them.
- **Collapse the Files section.** The table takes the whole window. This is the
  escape hatch when the window is too short for both.

### Molecule detail windows (#107)

- **Click three different molecule names.** Three windows, cascaded so each is
  visible rather than exactly beneath the last, and all three rows highlighted.
  Opening a second used to close the first.
- **Click an open row again.** Its window closes and the row un-highlights. The
  others are untouched.
- **Open nine.** The oldest closes rather than the click being refused, and its
  row un-highlights — which is what tells you which one went.
- **View → Close all molecule windows** shows the count and clears them.
- **Switch dataset with several open.** All of them close: they are keyed by row
  index, which means nothing against a different dataset.
- **Drag two windows apart, then close and reopen one.** It returns where you
  left it, not to the cascade position — each window keeps its own geometry
  because each has its own id.

### Settings (#121)

- **Settings sits beside File and View**, plain text like they are, and lights
  while the window is open.
  It is not a menu: a menu would overlay the canvas and close the moment you
  clicked a control, which is the opposite of what these settings need. It was
  a bare gear at the far right first, which nobody would have found.
- **Settings also appears in the View menu.** A window unreachable from there
  would be the worse inconsistency; the gear is a shortcut, not the only door.
- **Change the theme.** The whole app switches, and structure colours follow —
  check a drawn molecule in the Inspector, not just the chrome.
- **Change a structure option with a structure visible.** That is the entire
  reason this is a window rather than a modal: the query structure, the table
  thumbnails and any open detail window should all follow, while you watch.
- **Close every content window but leave Settings open.** The canvas must still
  say where windows come from — Settings holding preferences rather than data is
  why it doesn't count as something to look at.
- **The thumbnail toggle is here now**, back alongside the structure options
  rather than beside the table.

### Inspector (#106)

- **Change a Display option** — turn on explicit hydrogens, say. Every structure
  in the app should follow: the query, the table thumbnails, and any open detail
  window. That is the whole reason the options moved here from beside the table.
- **The thumbnail toggle stayed behind**, in Datasets, next to the table it
  governs. Confirm it still works from there.
- **The Query section labels itself with what was parsed**, not with what is in
  the box. Parse a SMILES, then type something different without pressing Parse
  and before the idle delay elapses — the label should still name the parsed one.
- **Collapse each section.** The header keeps its summary — the hit count, the
  parsed SMILES — so a collapsed section still says what is in it.

### Dataset invalidation (#102)

Do this in one sequence, because it is one behaviour:

1. Load Examples, compute fingerprints, run a search, expand a result row, and
   click a molecule name to open its detail window.
2. Now switch to a different entry in **Loaded Files**.

Everything derived from the old dataset must go: detail window closed, results
gone, expanded row collapsed, Fingerprint column back to "No", and the operation
sections back to "—". Three separate views notice this independently, via an
epoch on the shared state, so a partial reset is the interesting failure.

### Structure rendering, whenever it changes

- Methane draws (a lone carbon has no bond to imply a vertex — it drew blank
  once).
- Benzene's inner lines sit *inside* the ring, not outside.
- Two SMILES spellings of the same molecule score 1.0 against each other. A
  ring-closure aromaticity bug (#95) made phenol score 0.27 against itself, and
  the text views hid it completely.

### GPU / CPU

- Chrome or Edge for WebGPU; Firefox needs `dom.webgpu.enabled`.
- No GPU is a *pass*, not a failure: a red ⚠ chip with a reason on hover, and
  everything still works on CPU. Check the reason reads sensibly.
- Switch backend and re-run fingerprints; the section should report the backend
  that actually ran it.

## Adding a story to this list

Ask what the story does that a unit test cannot see, and write that down. Good
candidates: anything about layout, scrolling, focus, or hover; anything that only
happens on one platform; anything whose failure mode is "looks fine, wrong
result". If a story has nothing in that category, say so in its PR rather than
inventing a checklist item.

## Traps

- **`index.html` is the one file with a stable name.** The js and wasm are
  content-hashed, so those can't be served stale — but a cached `index.html`
  points at a hash that no longer exists, which shows as a stuck loading
  overlay. Hard-refresh.
- **A long-lived `trunk serve` on another port** rebuilds into `chem_app/dist/`
  whenever any file changes, in debug unless it was started with `--release`.
  `e2e.sh` builds into `dist-e2e/` to stay out of its way. Check for one before
  concluding anything about which build you are seeing.
- **`cargo test --workspace` twice at once will hang.** `chem_app/src/search.rs`
  has two tests that call `FingerprintSearch::new()`, which initialises a real
  GPU; two concurrent runs deadlock inside the driver, holding GPU handles until
  killed. CI never sees it because runners have no GPU. Run one at a time.
