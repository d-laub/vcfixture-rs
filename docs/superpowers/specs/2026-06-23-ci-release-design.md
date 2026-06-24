# CI & release pipeline design

Add GitHub Actions to `vcfixture` (repo `d-laub/vcfixture-rs`): a CI workflow that
lints, checks commit format, and tests pull requests across an OS/toolchain
matrix; and a manually-triggered release workflow that uses commitizen to bump
the version and regenerate the changelog, cuts a GitHub release from that
changelog, and publishes to crates.io via Trusted Publishing (OIDC).

## Goals

- Every PR (and push to `main`) is linted (`cargo fmt`, `cargo clippy`), has its
  commits validated as conventional-commits, and is tested on Linux + macOS
  across stable and the project's MSRV.
- A single, deliberate human action (`workflow_dispatch`) produces a release:
  version bump, changelog, git tag, GitHub release, and crates.io publish — in
  that one run, gated behind a passing test + package check.
- No long-lived crates.io credential stored in the repo: publishing authenticates
  with crates.io Trusted Publishing (OIDC).
- Local tooling (prek hooks, commitizen, pixi tasks) and CI stay conceptually
  aligned, but CI installs toolchains natively for speed (it does not run pixi).

## Non-goals

- **No fully-automatic release on merge.** Releases are intentional and
  human-triggered. Every qualifying merge does not ship.
- **No independently re-runnable publish step.** Bump, release, and publish live
  in one workflow run. If `cargo publish` fails after the tag/release already
  exist, recovery is a manual re-run (documented in the runbook below). This is
  an accepted trade-off for the simplicity of one workflow.
- **No pixi in CI.** pixi.lock remains the local-dev source of truth; CI uses
  marketplace actions for Rust and `uv` for commitizen.

## File layout

```
.github/workflows/ci.yml        # PR + push-to-main: lint, commit-check, test matrix
.github/workflows/release.yml   # workflow_dispatch: bump -> tag -> release -> publish
.cz.toml                        # commitizen config (new — none exists today)
Cargo.toml                      # add rust-version = "<detected MSRV>"
```

## Decisions (locked during brainstorming)

| Question | Decision |
|----------|----------|
| Release trigger / flow | One manual-dispatch workflow does bump → tag → test → GH release → publish |
| crates.io auth | Trusted Publishing (OIDC), no stored token |
| CI matrix | `{ubuntu-latest, macos-latest} × {stable, MSRV}` |
| CI tooling install | Native actions (`dtolnay/rust-toolchain`, `Swatinem/rust-cache`) + `uv` for commitizen |
| MSRV value | Auto-detected with `cargo-msrv`, pinned as `rust-version` in `Cargo.toml` |
| PR commit linting | `cz check` over every commit in the PR range |
| Bump mechanism | `commitizen-tools/commitizen-action` in the release workflow |

## Component 1 — `.cz.toml` (commitizen config)

No commitizen config exists yet (the `cz bump` pixi task and prek `commit-msg`
hook currently run against defaults). Add `.cz.toml`:

```toml
[tool.commitizen]
name = "cz_conventional_commits"
version_provider = "cargo"      # read/write [package] version in Cargo.toml
tag_format = "v$version"
update_changelog_on_bump = true
major_version_zero = true       # crate is 0.x: feat/fix bump the minor/patch, not to 1.0
```

`version_provider = "cargo"` means commitizen reads and writes the version directly
in `Cargo.toml`; there is no duplicated `version =` to keep in sync inside
`.cz.toml`.

## Component 2 — MSRV in `Cargo.toml`

Run `cargo-msrv` locally to find the lowest Rust toolchain that compiles the crate
with its current dependency tree (`noodles-*`, `ndarray`, `indexmap`, `rand`,
`proptest`, …) and pin it:

```toml
[package]
# ...
rust-version = "<detected>"     # e.g. "1.74"
```

The same value feeds the CI matrix's MSRV row. If `cargo-msrv` and the CI matrix
ever disagree, CI is the source of truth and the pin is corrected.

## Component 3 — `ci.yml` (PR + push-to-main)

Triggers: `pull_request` and `push` to `main`. A `concurrency` group keyed on the
ref cancels superseded in-progress runs.

Three jobs:

- **`lint`** — `ubuntu-latest`, stable. `cargo fmt --all -- --check` then
  `cargo clippy --all-features -- -D warnings`. Mirrors the prek `cargo-fmt` /
  `cargo-clippy` hooks.
- **`commit-check`** — `ubuntu-latest`, `pull_request` only. Checkout with
  `fetch-depth: 0`, `astral-sh/setup-uv`, `uv tool install commitizen`, then
  `cz check --rev-range origin/${{ github.base_ref }}..HEAD`. Validates every
  commit in the PR is conventional-commits so a later `cz bump` computes the
  correct version. (The commitizen-action only *bumps*; it has no check mode, so
  this job uses commitizen directly via `uv`.)
- **`test`** — matrix `os = [ubuntu-latest, macos-latest] × toolchain =
  [stable, <MSRV>]` (4 jobs). Each: `dtolnay/rust-toolchain@<toolchain>`,
  `Swatinem/rust-cache`, `cargo test --all-features --locked`. Clippy is **not**
  run on the MSRV rows — lint output drifts between toolchain versions, so MSRV
  rows only build and test.

## Component 4 — `release.yml` (manual release)

Trigger: `workflow_dispatch` only (no recursion guard needed, since it is never
fired by a push). Permissions: `contents: write` (push the bump commit + tag,
create the release) and `id-token: write` (crates.io OIDC). Single job, ordered so
nothing is committed, tagged, or published unless the crate is green:

1. **Checkout** with `fetch-depth: 0` (full history for changelog) and the default
   token persisted (used for the push in step 4).
2. **Setup** Rust stable (`dtolnay/rust-toolchain@stable`) + `Swatinem/rust-cache`.
3. **Gate** — `cargo test --all-features --locked` and `cargo publish --dry-run`.
   The dry-run validates packaging before any irreversible step; packaging
   structure does not depend on the version number, so running it pre-bump is
   sufficient.
4. **Bump** — `commitizen-tools/commitizen-action` with:
   - `push: false` (we own the push so the Cargo.lock sync stays in one commit),
   - `changelog_increment_filename: .release-notes.md` (writes just the new
     version's changelog section, for the GitHub release body),
   - `github_token: ${{ secrets.GITHUB_TOKEN }}`.

   The action bumps `Cargo.toml`, regenerates `CHANGELOG.md`, creates the bump
   commit, tags `v$X.Y.Z`, and exposes the new version as the `REVISION` env var
   (and a `version` step output).
5. **Lock sync + push** — the action bumps `Cargo.toml` but not `Cargo.lock`, so
   the committed lock would otherwise be stale (and break the next `--locked`
   build). Run `cargo update --workspace` to sync the package version in
   `Cargo.lock`, fold it into the bump commit (`git add Cargo.lock`,
   `git commit --amend --no-edit`, `git tag -f "v$REVISION"`), then push the
   branch and tag. Because the bump commit was never pushed, the branch push is a
   fast-forward and the tag is new to the remote — no force-push of shared history.
6. **GitHub release** — `gh release create "v$REVISION" --notes-file
   .release-notes.md --title "v$REVISION"`.
7. **Publish** — `rust-lang/crates-io-auth-action` (exchanges the OIDC token for a
   short-lived crates.io token) then `cargo publish`.

## Manual / one-time prerequisites

These are outside the workflow files and are done once during rollout:

1. **Reserve the crate name + first publish.** crates.io Trusted Publishing is
   configured per-crate, which requires the crate to already exist. The user is
   authenticated with `cargo` on this machine and has authorized running the first
   `cargo publish` (current version `0.1.0`) by hand during rollout to reserve the
   name. All subsequent releases go through `release.yml` + OIDC.
2. **Configure Trusted Publishing** on crates.io for `vcfixture`: repository
   `d-laub/vcfixture-rs`, workflow `release.yml`. No token is stored in GitHub.
3. **Branch protection** on `main`, if enabled, must allow the release job's push
   of the bump commit + tag (e.g. allow the `github-actions` actor to bypass, or
   leave `main` unprotected for direct pushes by Actions).

## Release runbook (operating the pipeline, post-rollout)

1. Merge conventional-commit PRs into `main` as usual.
2. When ready to release, run the **Release** workflow from the Actions tab
   (`workflow_dispatch`).
3. The run bumps the version from the accumulated commits, updates the changelog,
   tags, creates the GitHub release, and publishes to crates.io.
4. **If `cargo publish` (step 7) fails** after the tag and GitHub release already
   exist: the version is bumped and tagged but unpublished. Fix the cause, then
   re-publish manually from the tagged commit (`git checkout v$X.Y.Z &&
   cargo publish`) or re-run just the publish step. Do not re-run the whole
   workflow, which would attempt a second bump.

## Testing / validation

CI workflows are validated by their own first runs:

- Open the implementation PR; confirm `lint`, `commit-check`, and all four `test`
  matrix jobs pass. This exercises `ci.yml` end-to-end.
- `release.yml` is validated by the rollout: after the manual first publish and
  Trusted Publishing setup, trigger it once and confirm it bumps, tags, releases,
  and publishes the next version.
- Locally: `cz bump --dry-run` (already wired as the `bump-dry` pixi task) and
  `cargo publish --dry-run` confirm the bump and packaging before any CI run.
