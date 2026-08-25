#!/usr/bin/env bash
#
# Builds chem-app's web bundle and serves it for hands-on testing.
#
# The point of this script rather than `trunk build` plus a web server by hand:
# it is very easy to end up serving a bundle that isn't the code you think you
# are testing. That happened during v0.5.0 — a build was interrupted, an older
# bundle was left in the directory, and the preview looked unchanged for
# eleven minutes' worth of commits. So this script
#
#   * refuses to reuse a directory: the dist is removed and rebuilt every run;
#   * fails if the build did not actually produce fresher artefacts;
#   * stamps the commit into the page itself, so the app tells you which build
#     you are looking at (top right, next to the FPS readout);
#   * refuses to start if something else already holds the port, naming it.
#
# Usage:
#   crates/chem-app/e2e.sh                   # 127.0.0.1:8080
#   crates/chem-app/e2e.sh --port 8081
#   crates/chem-app/e2e.sh --address 0.0.0.0 # reachable from other machines
#
# For the native build there is nothing to serve: `cargo run --release -p chem-app`.

set -euo pipefail

ADDRESS="127.0.0.1"
PORT="8080"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --address) ADDRESS="$2"; shift 2 ;;
        --port)    PORT="$2";    shift 2 ;;
        -h|--help) sed -n '3,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# Derived from the script's own location rather than named, so moving the crate
# again cannot leave these disagreeing.
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$APP_DIR/../.." && pwd)"
DIST="$APP_DIR/dist-e2e"

command -v trunk >/dev/null || { echo "trunk not installed: cargo install trunk" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 needed to serve the bundle" >&2; exit 1; }

# --- what exactly are we about to test ---------------------------------------
BRANCH="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD)"
COMMIT="$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
BUILD_ID="$COMMIT"
if ! git -C "$REPO_ROOT" diff-index --quiet HEAD -- 2>/dev/null; then
    # Worth shouting about: a dirty tree means the bundle matches no commit, so
    # "it works on abc1234" is not a claim anyone can reproduce.
    BUILD_ID="$COMMIT-dirty"
    echo "WARNING: uncommitted changes — serving $BUILD_ID, which matches no commit"
fi

# --- refuse to fight over the port -------------------------------------------
if ss -ltnH 2>/dev/null | grep -q ":$PORT "; then
    echo "port $PORT is already in use:" >&2
    ss -ltnp 2>/dev/null | grep ":$PORT " >&2 || true
    echo "stop it, or pass --port" >&2
    exit 1
fi

# --- build, and prove it built ------------------------------------------------
STARTED_AT="$(date +%s)"
rm -rf "$DIST"
echo "building $BRANCH @ $BUILD_ID (release, wasm32) …"
( cd "$APP_DIR" && CHEM_BUILD_ID="$BUILD_ID" trunk build --release --dist "$DIST" )

WASM="$(find "$DIST" -name '*.wasm' -print -quit)"
[[ -n "$WASM" ]] || { echo "build produced no wasm in $DIST" >&2; exit 1; }
BUILT_AT="$(stat -c %Y "$WASM")"
if (( BUILT_AT < STARTED_AT )); then
    echo "artefact in $DIST predates this run — refusing to serve a stale bundle" >&2
    exit 1
fi

cat <<EOF

  serving  $BRANCH @ $BUILD_ID
  bundle   $(basename "$WASM")
  built    $(date -d "@$BUILT_AT" '+%H:%M:%S')
  open     http://$ADDRESS:$PORT/

  The build id also appears in the app, top right. If it does not match the
  line above, you are looking at a cached page — hard-refresh (Ctrl+Shift+R).
  Ctrl-C to stop.

EOF

exec python3 -m http.server "$PORT" --bind "$ADDRESS" --directory "$DIST"
