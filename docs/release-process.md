# Releasing dottir

This is the human-facing checklist for cutting a release. The
`.github/workflows/release.yml` workflow does the heavy lifting once
a `v*.*.*` tag is pushed.

## What a release produces

For every tag `vX.Y.Z`, the workflow publishes:

1. A **GitHub Release** at
   `https://github.com/kavonrtep/dottir/releases/tag/vX.Y.Z`, with:
   - `dottir-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz` (CLI + GUI,
     dynamic ELF, README/LICENSE/CHANGELOG included)
   - `dottir-vX.Y.Z-x86_64-pc-windows-msvc.zip` (same, MSVC build)
   - A `.sha256` file for each archive
   - Auto-generated release notes from commits since the previous tag
2. Two **conda packages** on the `petrnovak` Anaconda channel:
   - `linux-64::dottir=X.Y.Z`
   - `win-64::dottir=X.Y.Z`

   Install with:

       conda install -c petrnovak dottir

## Cutting a release

1. **Update the changelog.** Add a section for the new version to
   `docs/CHANGELOG.md`.

2. **Bump the workspace version.**

       cargo set-version 0.2.0    # one-shot, requires cargo-edit
       # OR hand-edit [workspace.package].version in the root Cargo.toml

3. **Commit and push to `main`.**

       git commit -am "release: v0.2.0"
       git push

4. **Tag and push the tag.**

       git tag v0.2.0
       git push origin v0.2.0

5. **Watch the workflow** at
   <https://github.com/kavonrtep/dottir/actions/workflows/release.yml>.
   Expect ~10 min total — the conda Windows build is the slow leg.

## What can go wrong

### Tag and Cargo.toml disagree

The `check-tag` job fails immediately. To recover:

    git tag -d v0.2.0
    git push --delete origin v0.2.0
    # fix the version, re-commit, re-tag, re-push

### A `build` job fails

Builds are independent — Linux can succeed while Windows fails (or
vice-versa). The release page is not created until **all** build
jobs succeed (the `publish` job has `needs: build` without a matrix
filter), so a partial failure is safe: no half-published release.

Re-trigger by deleting and re-pushing the same tag.

### Conda upload fails after the GitHub release has been published

This is the messy case. The GitHub release exists; the conda
package does not. Options:

1. Build and upload locally:

       conda install -n base conda-build anaconda-client
       DOTTIR_VERSION=0.2.0 conda build recipe/ \
           --output-folder conda-out -c conda-forge
       anaconda --token "$ANACONDA_API_TOKEN" upload \
           --user petrnovak --force \
           conda-out/linux-64/dottir-0.2.0-*.conda

2. Re-run only the failed `conda` job from the GitHub Actions UI.
   `--force` on `anaconda upload` makes this idempotent: re-running
   a successful upload is a no-op.

## Secrets required

| Secret | Used by | Purpose |
|--------|---------|---------|
| `GITHUB_TOKEN` | `publish` job (auto-provided) | `gh release create` |
| `ANACONDA_API_TOKEN` | `conda` job | `anaconda upload --user petrnovak` |

`GITHUB_TOKEN` is set up automatically by GitHub Actions.
`ANACONDA_API_TOKEN` lives in repo Settings → Secrets and variables
→ Actions.

## Local sanity-checks before tagging

Before pushing the tag, it's worth running locally:

    cargo fmt --all --check
    cargo clippy --workspace --all-targets --locked -- -D warnings
    cargo test --workspace --locked

The same checks run in `ci.yml` on every push, but catching them
before tagging saves a wasted release attempt.
