# Making a PAMTRA release

This describes how to cut a new PAMTRA version and get it published on
conda-forge. It's maintainer-facing (parallel to [AI.md](AI.md), which covers
day-to-day build/test); PyPI publishing is intentionally out of scope (see
"Why not PyPI" below).

## 1. Bump the version number

Three places currently need to agree on the version — bump all of them
together:

- [pyproject.toml](pyproject.toml) — `[project] version = "..."`
- [meson.build](meson.build) — `project(..., version: '...', ...)`
- [conda-recipe/recipe.yaml](conda-recipe/recipe.yaml) — `context.version`

Use [semantic versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`).

## 2. Tag and cut a GitHub release

```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```

Then create a GitHub Release from that tag (`gh release create vX.Y.Z` or via
the GitHub UI) — GitHub's auto-generated source tarball for the tag is what
the conda-forge recipe points at, so the release doesn't need any attached
build artifacts, just the tag itself.

Grab the tarball's checksum for the next step:

```bash
curl -sL https://github.com/igmk/pamtra/archive/refs/tags/vX.Y.Z.tar.gz | shasum -a 256
```

## 3. Update the conda-forge recipe

[conda-recipe/recipe.yaml](conda-recipe/recipe.yaml) uses the [v1 recipe
format](https://conda-forge.org/blog/2025/02/27/conda-forge-v1-recipe-support/)
(rattler-build) — this is what `conda-forge/staged-recipes` requires for new
submissions; the old conda-build `meta.yaml` format is deprecated there.

Update `source:` from the local-path form used for dev testing to the tagged
release:

```yaml
source:
  url: https://github.com/igmk/pamtra/archive/refs/tags/vX.Y.Z.tar.gz
  sha256: <checksum from step 2>
```

Validate locally before pushing anything upstream — this is the fast
feedback loop, catches recipe mistakes before conda-forge's CI does:

```bash
# one-time: pixi global install rattler-build
rattler-build build --recipe conda-recipe/recipe.yaml -c conda-forge
```

This builds the package and runs its embedded test (the `pyPamtra` import
check, plus the full `pytest tests/` suite) inside a fresh conda environment
— a good end-to-end check independent of the pip-based CI.

## 4a. First-ever conda-forge submission (one-time)

Only needed once, to get PAMTRA onto conda-forge in the first place:

1. Fork `conda-forge/staged-recipes`, branch from `main`.
2. Add the recipe at `recipes/pamtra/recipe.yaml` (their required layout).
3. Open a PR there, following their PR checklist, and request review from
   `@conda-forge/help-python` (or whichever team their checklist points to).
4. Once merged, conda-forge's automation creates `conda-forge/pamtra-feedstock`
   and does the first build/publish. Maintainers listed in the recipe's
   `extra.recipe-maintainers` get commit rights on that feedstock repo.

## 4b. Every subsequent release

Since PAMTRA isn't on PyPI, conda-forge's auto-tick bot (which mostly tracks
PyPI/CRAN releases) won't reliably notice new GitHub-only tags. So each
release needs a manual PR against `conda-forge/pamtra-feedstock`:

1. Clone/update your fork of `pamtra-feedstock`.
2. Edit `recipe/recipe.yaml`: bump `version`, update `sha256` for the new
   tarball (same command as step 2 above), reset `build.number` to `0`.
3. Open a PR against the feedstock; its own CI (linux/osx/win matrix) builds
   and tests it.
4. Merge once green — conda-forge publishes the new version automatically.

## Verification checklist

- [ ] `rattler-build build --recipe conda-recipe/recipe.yaml -c conda-forge`
      succeeds locally against the real tagged tarball (not the local `path:`
      source used for day-to-day recipe edits)
- [ ] The feedstock PR's own CI (linux-64/osx-64/osx-arm64/win-64) is green
      before merging
- [ ] After merge, `conda install -c conda-forge pamtra` works in a clean
      environment

## Why not PyPI

PAMTRA links against netCDF (C + Fortran bindings, which pull in HDF5/zlib/
curl), FFTW, and OpenBLAS. A PyPI wheel can't depend on system packages the
way a conda package can — those libraries would have to be bundled *inside*
each wheel (via `cibuildwheel` + `auditwheel`/`delocate`/`delvewheel`), which
is a second, non-trivial CI pipeline. conda-forge gets this for free because
conda already manages those dependencies as packages. If PyPI becomes
worthwhile later, an sdist-only release (no compiled wheel, `pip install`
compiles from source using the user's local toolchain — same as `pip
install .` today) would be the low-effort first step, and would also let
conda-forge's auto-tick bot pick up new versions automatically.
