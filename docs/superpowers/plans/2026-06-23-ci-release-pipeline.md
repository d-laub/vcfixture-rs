# CI & Release Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add GitHub Actions CI (lint + commit-check + test matrix on PRs) and a manually-triggered commitizen-driven release workflow that bumps, changelogs, tags, GitHub-releases, and publishes `vcfixture` to crates.io via Trusted Publishing.

**Architecture:** Two workflows under `.github/workflows/` — `ci.yml` (runs on PRs and pushes to `main`) and `release.yml` (`workflow_dispatch` only). Releases are driven by `commitizen-tools/commitizen-action` with `push: false` so we own the push and can fold the `Cargo.lock` sync into the single bump commit. crates.io publishing uses OIDC (no stored token). commitizen reads/writes the version directly in `Cargo.toml` (`version_provider = "cargo"`).

**Tech Stack:** GitHub Actions; `dtolnay/rust-toolchain`, `Swatinem/rust-cache`, `astral-sh/setup-uv`, `commitizen-tools/commitizen-action`, `rust-lang/crates-io-auth-action`; commitizen (conventional commits); cargo-msrv; actionlint; pixi (local dev only).

## Global Constraints

- Repository: `d-laub/vcfixture-rs`. Crate name: `vcfixture`. Current version: `0.1.0`.
- crate is **0.x** — `major_version_zero = true`; feat/fix bump minor/patch, never to 1.0 automatically.
- crates.io auth is **Trusted Publishing (OIDC) only** — never add a `CARGO_REGISTRY_TOKEN` secret.
- CI does **not** run pixi; it installs toolchains via marketplace actions + `uv`.
- Test matrix: `{ubuntu-latest, macos-latest} × {stable, MSRV}`. Tests run `cargo test --all-features --locked`. Clippy runs **only** on stable/Linux, never on MSRV rows.
- tag format is `v$version` (e.g. `v0.1.0`).
- All commits made while implementing this plan MUST be conventional-commits (`feat:`, `fix:`, `ci:`, `docs:`, `chore:`, `build:`), because the `commit-check` CI job validates them.
- Spec: `docs/superpowers/specs/2026-06-23-ci-release-design.md`.

---

### Task 1: commitizen config (`.cz.toml`)

No commitizen config exists yet; `cz` currently runs on defaults. This adds a config so `version_provider = "cargo"` reads/writes the version in `Cargo.toml` and bumps generate the changelog with the `v$version` tag format.

**Files:**
- Create: `.cz.toml`

**Interfaces:**
- Produces: a valid commitizen config consumed by Task 4 (`commitizen-action`) and Task 5 (`cz changelog`). Version source of truth is `Cargo.toml [package] version`.

- [ ] **Step 1: Create `.cz.toml`**

```toml
[tool.commitizen]
name = "cz_conventional_commits"
version_provider = "cargo"
tag_format = "v$version"
update_changelog_on_bump = true
major_version_zero = true
```

- [ ] **Step 2: Verify the config parses and finds the version**

Run: `pixi run bump-dry`
(this is the existing `bump-dry = "cz bump --dry-run"` pixi task)
Expected: command exits without a config/parse error and either prints a planned bump (e.g. `bump: version 0.1.0 → 0.x.0`) or reports there are no eligible commits. Either outcome confirms `.cz.toml` is valid and the cargo version provider resolved `0.1.0`. A traceback or "config file not found / invalid" is a FAIL.

- [ ] **Step 3: Commit**

```bash
git add .cz.toml
git commit -m "build: add commitizen config with cargo version provider"
```

---

### Task 2: `Cargo.toml` crates.io metadata + MSRV pin

Enrich the package metadata so the first crates.io publish is clean, and pin a detected MSRV that the CI matrix will use. Both edits touch `Cargo.toml`, so they are one task.

**Files:**
- Modify: `Cargo.toml` (`[package]` table)

**Interfaces:**
- Produces: `MSRV` — the detected minimum supported Rust version string (e.g. `1.74`). Task 3's test matrix consumes this exact value. Record it in the commit message and in the matrix.

- [ ] **Step 1: Add crates.io metadata to `[package]`**

Add these keys to the existing `[package]` table in `Cargo.toml` (keep the existing `name`, `version`, `edition`, `license`, `description`):

```toml
repository = "https://github.com/d-laub/vcfixture-rs"
readme = "README.md"
keywords = ["vcf", "testing", "fixtures", "bioinformatics", "proptest"]
categories = ["development-tools::testing", "science"]
```

- [ ] **Step 2: Verify packaging is valid with the new metadata**

Run: `cargo publish --dry-run`
Expected: PASS — `Packaging vcfixture v0.1.0`, `Verifying vcfixture v0.1.0`, finishes with no errors. (`keywords` ≤ 5 and ≤ 20 chars each; `categories` must be valid crates.io slugs — the dry-run rejects invalid ones, so this is the real check.)

- [ ] **Step 3: Detect the MSRV**

Install and run cargo-msrv (note: `find` bisects by building repeatedly with the full dependency tree, so this can take several minutes):

```bash
cargo install cargo-msrv --locked
cargo msrv find -- cargo check --all-features
```
Expected: prints `Minimum Supported Rust Version: <X.Y.Z>`. Use the `X.Y` (major.minor) form as the MSRV value below. If `cargo msrv` is unavailable or fails, fall back: try `rustup toolchain install 1.74 && cargo +1.74 check --all-features`, adjusting the number upward until it compiles; the lowest that compiles is the MSRV.

- [ ] **Step 4: Pin `rust-version` in `[package]`**

Add to the `[package]` table (substitute the value detected in Step 3):

```toml
rust-version = "<MSRV from Step 3, e.g. 1.74>"
```

- [ ] **Step 5: Verify the pinned toolchain compiles the crate**

```bash
rustup toolchain install <MSRV>
cargo +<MSRV> test --all-features --locked
```
Expected: PASS — builds and all tests pass on the pinned toolchain. If it fails, raise `rust-version` to the lowest version that passes and re-run.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml
git commit -m "build: add crates.io metadata and pin MSRV <MSRV>"
```

---

### Task 3: CI workflow (`ci.yml`) + actionlint tooling

Add the PR/push CI workflow and a local `actionlint` task (via pixi) to validate workflow YAML before pushing.

**Files:**
- Create: `.github/workflows/ci.yml`
- Modify: `pixi.toml` (add `actionlint` dependency + `lint-actions` task)

**Interfaces:**
- Consumes: `MSRV` from Task 2 — substitute the literal value into the matrix `toolchain` list.
- Produces: a CI workflow whose `commit-check` job assumes the conventional-commit constraint from Global Constraints.

- [ ] **Step 1: Add actionlint to pixi**

In `pixi.toml`, add to `[dependencies]`:

```toml
actionlint = "*"
```

And add to `[tasks]`:

```toml
lint-actions = "actionlint"
```

- [ ] **Step 2: Install the new tool**

Run: `pixi install`
Expected: resolves and installs `actionlint` with no error. (If `actionlint` is not on the configured conda channels, instead install the binary per https://github.com/rhysd/actionlint and run it directly in Step 5; remove the pixi changes.)

- [ ] **Step 3: Create `.github/workflows/ci.yml`**

Substitute `<MSRV>` (two places: the matrix entry) with the value from Task 2.

```yaml
name: CI

on:
  pull_request:
  push:
    branches: [main]

concurrency:
  group: ci-${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

permissions:
  contents: read

jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt, clippy
      - uses: Swatinem/rust-cache@v2
      - name: rustfmt
        run: cargo fmt --all -- --check
      - name: clippy
        run: cargo clippy --all-features -- -D warnings

  commit-check:
    if: github.event_name == 'pull_request'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0
      - uses: astral-sh/setup-uv@v5
      - name: install commitizen
        run: uv tool install commitizen
      - name: check commit messages
        run: cz check --rev-range origin/${{ github.base_ref }}..HEAD

  test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
        toolchain: [stable, "<MSRV>"]
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@master
        with:
          toolchain: ${{ matrix.toolchain }}
      - uses: Swatinem/rust-cache@v2
      - name: test
        run: cargo test --all-features --locked
```

- [ ] **Step 4: Add the pixi changes and workflow, then lint the workflow**

Run: `pixi run lint-actions`
Expected: no output / exit 0 (actionlint reports no problems). Any reported error (unknown action input, bad expression, YAML issue) is a FAIL — fix inline and re-run.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/ci.yml pixi.toml pixi.lock
git commit -m "ci: add PR lint, commit-check, and test matrix workflow"
```

---

### Task 4: Release workflow (`release.yml`)

Add the manual-dispatch release workflow: gate on tests, bump+changelog+tag via commitizen-action (no push), sync the lockfile into the bump commit, push, create the GitHub release, and publish to crates.io via OIDC.

**Files:**
- Create: `.github/workflows/release.yml`

**Interfaces:**
- Consumes: `.cz.toml` from Task 1 (commitizen-action reads it); crates.io Trusted Publishing config from Task 5 (must exist before the publish step succeeds, but the workflow file can be committed first).
- Produces: a tag `v$version`, a GitHub release, and a crates.io publish per run.

- [ ] **Step 1: Create `.github/workflows/release.yml`**

```yaml
name: Release

on:
  workflow_dispatch:

permissions:
  contents: write
  id-token: write

jobs:
  release:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0

      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2

      # Gate: nothing is bumped, tagged, or published unless these pass.
      - name: test
        run: cargo test --all-features --locked
      - name: package check
        run: cargo publish --dry-run

      # Bump version + changelog + tag. push:false so we own the push and can
      # fold the Cargo.lock sync into the single bump commit.
      - name: Bump version
        id: cz
        uses: commitizen-tools/commitizen-action@0.23.1
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          push: false
          changelog_increment_filename: .release-notes.md

      # commitizen bumps Cargo.toml but not Cargo.lock; sync it and amend it into
      # the bump commit so the committed lockfile is never stale.
      - name: Sync lockfile and push
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          cargo update --workspace
          git add Cargo.lock
          git commit --amend --no-edit
          git tag -f "v${{ steps.cz.outputs.version }}"
          git push origin "HEAD:${{ github.ref_name }}"
          git push origin "v${{ steps.cz.outputs.version }}"

      - name: Create GitHub release
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          gh release create "v${{ steps.cz.outputs.version }}" \
            --title "v${{ steps.cz.outputs.version }}" \
            --notes-file .release-notes.md

      - name: Authenticate to crates.io
        id: auth
        uses: rust-lang/crates-io-auth-action@v1

      - name: Publish
        run: cargo publish
        env:
          CARGO_REGISTRY_TOKEN: ${{ steps.auth.outputs.token }}
```

- [ ] **Step 2: Lint the workflow**

Run: `pixi run lint-actions`
Expected: exit 0, no problems reported. Fix any reported issue inline and re-run.
(Note: actionlint validates syntax/expressions/known action schemas. The `commitizen-action` version `0.23.1` and `crates-io-auth-action@v1` tags should be confirmed current at implementation time; if actionlint or the marketplace flags an unknown ref, bump to the latest released tag.)

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/release.yml
git commit -m "ci: add manual release workflow with commitizen bump and OIDC publish"
```

- [ ] **Step 4: Push the branch and open a PR to validate CI end-to-end**

```bash
git push -u origin HEAD
gh pr create --fill
```
Expected: on the PR, the `lint`, `commit-check`, and all four `test` matrix jobs run and pass. This is the real validation of `ci.yml` (Task 3) and confirms the workflow YAML is accepted by GitHub. `release.yml` does not run on PRs (it is `workflow_dispatch` only) and is validated during rollout (Task 5).

- [ ] **Step 5: Merge the PR**

After CI is green and review (if any) is complete, merge to `main`.

---

### Task 5: Rollout — first publish, baseline tag, and Trusted Publishing

One-time steps to make the pipeline operational. The first `cargo publish` is irreversible (a published version cannot be deleted, only yanked) — confirm the crate is in final shape before running it. The user is authenticated with `cargo` on this machine and has authorized the first publish.

**Files:**
- Create: `CHANGELOG.md` (generated)

**Interfaces:**
- Consumes: everything from Tasks 1–4 merged to `main`.
- Produces: the reserved crate name on crates.io, baseline tag `v0.1.0`, initial `CHANGELOG.md`, and a configured Trusted Publisher — after which `release.yml` is fully operational.

- [ ] **Step 1: Confirm clean state on `main`**

```bash
git checkout main && git pull
git status
```
Expected: on `main`, up to date, working tree clean, Tasks 1–4 merged.

- [ ] **Step 2: First manual publish (IRREVERSIBLE — confirm before running)**

```bash
cargo publish --dry-run   # final sanity check
cargo publish             # reserves the name; publishes vcfixture 0.1.0
```
Expected: `Uploading vcfixture v0.1.0` then success. After this, `vcfixture` exists on crates.io.

- [ ] **Step 3: Create the baseline tag so commitizen has a starting point**

commitizen computes the next version from commits since the last `v*` tag. Without a baseline tag, the first `cz bump` would consider the entire history. Tag the published commit:

```bash
git tag v0.1.0
git push origin v0.1.0
```

- [ ] **Step 4: Generate the initial changelog**

```bash
cz changelog
```
Expected: writes `CHANGELOG.md` containing a `0.1.0` section derived from history. Review it, then commit:

```bash
git add CHANGELOG.md
git commit -m "docs: add initial changelog"
git push origin main
```
(If `main` is protected against direct pushes, open a tiny PR instead.)

- [ ] **Step 5: Configure crates.io Trusted Publishing**

In the crates.io web UI: `vcfixture` → Settings → Trusted Publishing → Add. Set:
- Repository owner: `d-laub`
- Repository name: `vcfixture-rs`
- Workflow filename: `release.yml`
- (Environment: leave blank unless you add a GitHub Actions environment to `release.yml`.)

Expected: a trusted publisher entry now exists; no `CARGO_REGISTRY_TOKEN` secret is created in GitHub.

- [ ] **Step 6: Confirm branch protection allows the release push**

If `main` has branch protection, ensure the release job can push the bump commit + tag (allow the `github-actions` actor / bot to bypass required reviews and push, or leave the relevant restrictions off for Actions). If `main` is unprotected, no action needed.

- [ ] **Step 7: Validate the release pipeline end-to-end**

When there is at least one releasable conventional commit on `main` after `v0.1.0` (a `fix:` or `feat:`), trigger the **Release** workflow from the GitHub Actions tab (`Run workflow`). Confirm in order:
- the run bumps the version (Cargo.toml + CHANGELOG.md updated),
- a new `v$version` tag and matching GitHub release appear with notes from the changelog,
- `cargo publish` succeeds and the new version shows on crates.io.

Note: this triggers a real publish of the next version. If you have no changes to release yet, defer this step until you do — the pipeline is otherwise ready.

**Recovery (if `cargo publish` fails after tag/release exist):** the version is bumped, tagged, and GitHub-released but not on crates.io. Fix the cause, then publish manually from the tagged commit — `git checkout v$version && cargo publish` — rather than re-running the whole workflow (which would attempt a second bump).

---

## Self-Review notes

- **Spec coverage:** `.cz.toml` (Task 1) ↔ spec §Component 1; MSRV pin (Task 2) ↔ §Component 2; `ci.yml` lint/commit-check/test-matrix (Task 3) ↔ §Component 3; `release.yml` gate→bump→lock-sync→release→OIDC-publish (Task 4) ↔ §Component 4; first publish + Trusted Publishing + branch protection + baseline tag + changelog (Task 5) ↔ §Manual prerequisites and §Release runbook. crates.io metadata enrichment (Task 2) is an addition beyond the spec, justified to make the first irreversible publish clean.
- **Version coupling:** the MSRV value is detected at implementation time (Task 2 Step 3) and substituted into `Cargo.toml` (Task 2) and the CI matrix (Task 3) — a genuine produce/consume dependency, not a placeholder.
- **Action version tags** (`commitizen-action@0.23.1`, `crates-io-auth-action@v1`, `setup-uv@v5`) should be confirmed current when implementing; Task 4 Step 2 calls this out.
