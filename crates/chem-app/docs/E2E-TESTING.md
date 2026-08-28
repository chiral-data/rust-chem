# Testing chem-app by hand

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
crates/chem-app/e2e.sh                 # builds, then serves on 127.0.0.1:8080
crates/chem-app/e2e.sh --port 8081
crates/chem-app/e2e.sh --address 0.0.0.0  # reachable from another machine
```

Native, where there is a display:

```bash
cargo run --release -p chem-app
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
- **Load three files, then remove the *first* while looking at the third.** The
  third must still be the active one, with its fingerprints and results intact.
  This is the case the feature can get quietly wrong: every entry after the
  removed one shifts down a slot, so an unadjusted index would leave you looking
  at a different dataset with nothing invalidated.
- **Now remove the one you're looking at.** A neighbour takes over, and its
  fingerprints, results and open molecule windows are gone — they belonged to
  the dataset that just went.
- **Remove down to one.** The remove control greys out, and hovering says why.
- **Look at the remove control itself.** It is two painted lines, not a
  character. `✕` (U+2715) was tried and rendered as a missing-glyph box — egui
  paints its own window close button the same way rather than using a character,
  which is the hint. **Any new glyph needs looking at on screen**: the bundled
  font carries a subset, so a character that renders in your terminal and in this
  file may still be a box in the app.
- **Drag a column edge.** Columns are resizable now; SMILES is the one worth
  widening.
- **Turn on "Show structures in table" with a SMILES dataset.** Cells show a dash,
  not a structure — coordinates don't exist yet. Run **2D Coordinates** in
  Operations and they appear. With an SDF file they should be there immediately,
  since the file supplies them.
- **Collapse the Files section.** The table takes the whole window. This is the
  escape hatch when the window is too short for both.

### SVG export (#109)

- **Export a molecule from its detail window**, then open the file. Bonds should
  stop at the edge of atom labels, not strike through them — that inset is
  computed from an estimate for the font the file declares, not from egui's
  metrics, and getting it wrong is the specific failure this design avoids.
- **Export the query from the Inspector** too; both surfaces have a button,
  since with several detail windows open there is no single "current" structure.
- **Export something aromatic** (benzene, phenol). The inner ring lines must be
  dashed in the file, not solid — the dash survives as an attribute rather than
  being cut into segments.
- **Export while the app is in dark theme.** The file must still come out with
  dark strokes on nothing, ready for a light document. It does not inherit the
  live palette.
- **Open the file in a browser and in an editor** if you can. It is valid XML —
  verified with a real parser — but how a viewer substitutes fonts is exactly
  what cannot be checked here.
- **Try a molecule with an awkward name** — anything with slashes or spaces. The
  suggested filename should be sanitised, and never start with a dot.

### Structures in result rows (#111)

- **Run a search on a dataset with coordinates.** Every hit draws its molecule at
  96x72, larger than the table's 64x48 because a result is being compared
  against the query rather than scanned.
- **Search a freshly parsed SMILES dataset** without running 2D Coordinates
  first: rows show a dash, not a blank or a broken box. Then run it and they
  fill in — the same convention as the dataset table.
- **The expansion is fingerprints only now.** Click **▼ Why?** and you get the
  shared-bit count, the either-bit count, their division, and both grids stacked.
  The atom and bond lists moved out; they are in a detail window, opened from
  the Datasets table.
- **Check the division matches the row's score.** The panel shows
  `intersection/union = x.xxx`, and it should equal the similarity printed in the
  row above — it comes from the same function that ranked the results, so a
  disagreement means a real bug.
- **One grid, three colours.** Blue for bits in both fingerprints, amber for this
  molecule only, violet for the query only. Count the blue against the "shared"
  figure if you want to be sure the legend describes the grid.
- **Search a molecule against itself** — the top hit at 1.000 should be entirely
  blue, with no amber or violet at all. Any other colour on a perfect match means
  the grid and the score disagree.
- **Switch the theme with the panel open.** Unset cells follow it. They were a
  hardcoded near-white, which only became wrong once the theme was a choice.
- **Check phenol against itself** — both spellings are in the examples. Scoring
  1.000 with two visibly identical drawings side by side is the check that a
  fingerprinting bug would now be obvious rather than hidden behind text.

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

### Persistence (#108)

This is the one section that cannot be checked in a single sitting — every item
needs a restart, and on web a reload of the same browser.

- **Move and resize the windows, close one, change the theme and the fingerprint
  radius. Restart.** All of it should come back: geometry, which windows were
  open, the theme, the display options, radius, size and top-k.
- **Load a file and compute fingerprints, then restart.** The dataset is *not*
  restored, and neither are the fingerprints or results. That is deliberate — a
  fingerprint cache restored under settings you can no longer see is worse than
  none.
- **Drag a window mostly off-screen, then View → Reset layout.** It comes back.
  Without the reset this state would persist across restarts, which is what
  makes the menu item necessary rather than a convenience.
- **On web, check the same browser** — persistence is `localStorage`, so a
  different browser or a private window starts fresh, correctly.
- **Then check it survives a version change.** Not something to do by hand every
  time, but worth knowing the shape: a saved file from an older build must be
  ignored rather than crash, and the app should open on defaults.

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
- **A long-lived `trunk serve` on another port** rebuilds into `crates/chem-app/dist/`
  whenever any file changes, in debug unless it was started with `--release`.
  `e2e.sh` builds into `dist-e2e/` to stay out of its way. Check for one before
  concluding anything about which build you are seeing.
- **GPU tests share one device now.** `cargo test --workspace` used to hang on a
  machine with a GPU — several tests each requested their own adapter and device,
  and concurrent creation against one GPU deadlocks in the driver. Fixed in
  #134 by sharing one instance per crate, the pattern #19 established. If a run
  ever hangs again, look for a test calling a `*::new()` that builds its own
  context instead of taking the shared one, and check `nvidia-smi -pm 1` is on.
