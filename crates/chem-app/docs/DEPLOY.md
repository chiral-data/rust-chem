# Publishing the workbench

The browser build is served at **https://chem.chiral.one**.

GitHub Actions builds it; Vercel only hosts. Vercel never sees a Rust toolchain
— `.github/workflows/deploy.yml` runs `trunk build --release`, arranges the
three files into Vercel's prebuilt output layout, and uploads that.

## What triggers a deploy

**Push to `v080`.** That is the only live trigger today.

The workflow also declares `workflow_dispatch`, and it does **not work yet**:
GitHub dispatches only workflows that exist on the default branch, and this file
lives on `v080`. So there is no "Run workflow" button and `gh workflow run`
returns 404 until `main` carries it. To redeploy by hand before then, push a
commit or re-run the last job from the Actions tab.

`main` is not a trigger yet. Its HEAD is the v0.7.0 release, 34 commits behind
`v080`, so deploying from it would publish an app with none of the v0.8.0 work.
When v0.8.0 releases, **replace** `v080` with `main` rather than listing both:
production is a single alias, so two branches deploying to it means last push
wins, and a push to a stale milestone branch could overwrite the release.

The deploy is its own only gate. `ci.yml` does not run on pushes to a milestone
branch, so there is no CI run to depend on, and `v080` is unprotected — a direct
push can reach production if it compiles. The wasm clippy step and the bundle
assertions are what stand in the way; the build id and `vercel rollback` are what
you have afterwards.

The deploy does not wait for tests. Everything reaching `v080` came through a PR
that CI checked, and gating would re-run exactly that; a direct push that was
never PR'd deploys unchecked, which is the trade.

## The three secrets

| secret | where it comes from |
|---|---|
| `VERCEL_TOKEN` | vercel.com/account/tokens, **Scope set to the `chiral` team**, not a personal account. The CLI cannot create one |
| `VERCEL_ORG_ID` | the `chiral` team id, `vercel project inspect` or `.vercel/project.json` |
| `VERCEL_PROJECT_ID` | the same two places |

The workflow's first step fails with the names of any that are missing, rather
than letting the Vercel CLI fail somewhere less obvious. The token goes through
the environment, never into a command line: a `${{ }}` interpolation lands in
the script the runner writes to disk, and this repository is public.

To rotate: create the new token, `gh secret set VERCEL_TOKEN`, then delete the
old one in the dashboard. Nothing else changes.

**Write the token's expiry date here when you create it.** A secret that expires
with nobody expecting it is the likeliest cause of a future mystery 403.

### Reading a credential failure

The three messages mean different things, and only the third is about scope:

| message | cause |
|---|---|
| `The token provided via VERCEL_TOKEN ... is not valid` | wrong, truncated or revoked token |
| `Project not found ({"VERCEL_PROJECT_ID":…,"VERCEL_ORG_ID":…})` | the two ids do not belong together |
| `Could not retrieve Project Settings` | the token is valid and **cannot see this project** — nearly always scoped to a personal account instead of the team |

The Scope dropdown on the token form defaults to your personal account, which
is how the third one happens. The workflow now prints what the token can reach
before it depends on it, and fails naming the scope if the team is missing.

## Checking a deploy actually landed

```sh
curl -sI https://chem.chiral.one | grep -iE 'HTTP|cache-control'
curl -sI https://chem.chiral.one/chem-app-<hash>_bg.wasm \
  | grep -iE 'HTTP|content-type|content-encoding|cache-control'
```

- `index.html` → `cache-control: no-cache`
- the hashed js and wasm → `public, max-age=31536000, immutable`
- the wasm → `content-type: application/wasm`, and compressed in transit

The in-app build id is compiled into the wasm and drawn on a canvas, so no HTTP
client can read it. `build-info.json` is its machine-readable twin, written by
the same job from the same variable:

```sh
curl -s https://chem.chiral.one/build-info.json
```

It names the commit, the ref, and the Actions run that built it. Open the page
and the id in the **top right** must read the same seven characters — that pair
is how you tell what is live from what your browser kept. The workflow asserts
the match on every deploy and prints the id in its last step.

## Undoing one

These commands resolve the project from a `.vercel/` link, which is gitignored
and absent in a fresh clone — so give them the ids, or they will not find it:

```sh
export VERCEL_ORG_ID=team_SV5XJeJeURf8s8XEhKhAsH27
export VERCEL_PROJECT_ID=prj_rT8mRl2SE04api6dn8A58OhpQGiY

vercel rollback --scope chiral-b84fb3de                    # to the previous production deploy
vercel rollback status --scope chiral-b84fb3de             # watch it land
vercel promote <deployment-url> --scope chiral-b84fb3de    # or forward to a specific one
```

Both re-point the alias and take effect immediately; neither rebuilds anything.
The ids are identifiers rather than credentials, and they survive a project
rename — only the token needs care.

## Traps

- **A stale `index.html`.** It is the only file with a stable name, so a cached
  copy points at a content hash that no longer exists — which shows as a stuck
  loading overlay, not an error. That is what its `no-cache` header is for;
  hard-refresh if you meet it anyway.
- **No build id in the corner** means `CHEM_BUILD_ID` was not set for
  `trunk build`. It is compiled in through `option_env!`, so it cannot be
  injected into the page afterwards, and when unset the app displays nothing
  rather than a placeholder. A bundle in that state is unattributable.
- **`integrity` hashes.** `index.html` ships SRI hashes and `<base href="/" />`
  with absolute asset paths, so the app must be served from a domain root and
  the host must not rewrite the js or wasm bytes. Transport compression is fine,
  being decoded before the hash is checked. Serving under a path prefix needs
  `trunk build --public-url`.
- **Deployment protection.** A Vercel project can be set to answer visitors with
  a login instead of the app — a site that looks deployed and is unreachable.
  `curl -sI` from a shell with no Vercel session is what settles it.
- **`.vercel/`** is written by the CLI and by the output layout, and is ignored
  in `.gitignore`. Keep it that way: #271 was a 4.4 MB wasm committed by
  accident, and the blob is still in the object store.

- **Do not connect the GitHub repository to the Vercel project.** One click in
  the dashboard makes every push to every branch trigger a Vercel-side build —
  on a builder with no Rust toolchain, so a wall of failures, and failing checks
  on other people's PRs since this repo is public. The project is deliberately
  unconnected and its build settings blank; Actions pushes a finished directory.
- **Do not add an SPA fallback rewrite.** There is no client-side router (one
  canvas, state in `localStorage`), and a catch-all to `index.html` would answer
  a stale client's request for a vanished `chem-app-<oldhash>_bg.wasm` with HTML
  and a 200 — turning a clean 404 into a MIME error and breaking the one
  diagnosis the trap above teaches. Unknown paths should 404.
- **Do not serve pre-compressed `.br` sidecars as identity.** That rewrites the
  bytes SRI covers and the browser refuses the module. Transport compression is
  fine, being decoded before the hash is checked.

## Headers, and why they are not in a `vercel.json`

They are in `crates/chem-app/vercel-output.json`, which the workflow copies to
`.vercel/output/config.json`. The upload is that output directory rather than the
repository, so a `vercel.json` in the tree would never be read.

## Running it locally instead

`crates/chem-app/e2e.sh` builds the same bundle and serves it, stamping the same
build id from your working tree — with `-dirty` appended when the tree has
uncommitted changes. See [E2E-TESTING.md](E2E-TESTING.md).
