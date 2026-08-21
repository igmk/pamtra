# Making a PAMTRA release

This describes how to cut a new PAMTRA version and get it published on both
conda-forge and PyPI. It's maintainer-facing (parallel to [AI.md](AI.md),
which covers day-to-day build/test).

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

Pushing the tag alone (before creating any GitHub Release) already triggers
[wheels.yml](.github/workflows/wheels.yml)'s `publish-pypi` job: it builds
wheels for every supported platform/Python version and publishes them to
PyPI via trusted publishing (OIDC, no stored token) — no separate manual
step. See "PyPI" below for the one-time setup this depends on, and check
[the Actions tab](https://github.com/igmk/pamtra/actions/workflows/wheels.yml)
to confirm it went green before moving on.

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
- [ ] `wheels.yml`'s `publish-pypi` job (triggered by the tag push in step 2)
      is green, and `pip install pamtra` works in a clean environment/venv
      with no system libraries preinstalled

## PyPI

`pip install pamtra` works via [wheels.yml](.github/workflows/wheels.yml)'s
`build` + `publish-pypi` jobs (`cibuildwheel`, triggered on `vX.Y.Z` tags).
No conda/system libraries needed at install time — PAMTRA's own C/Fortran
dependencies are bundled into the wheel itself:

- **OpenBLAS**: `tools/build_openblas_static.sh` builds it single-threaded
  (`USE_THREAD=0` — PAMTRA only uses it for small per-particle T-matrix
  solves, not large GEMMs, so this costs nothing in practice) and statically,
  and `meson.build`'s `pyPamtraLib` target link-time-restricts its exported
  symbol table to just `PyInit_pyPamtraLib`. Both matter because numpy/scipy
  wheels already bundle their own dynamically-linked OpenBLAS, and two
  copies loaded into the same process can collide (duplicate global symbols,
  shared thread-pool state) regardless of whether the versions match —
  matching versions alone doesn't rename or namespace anything.
- **FFTW**: `tools/build_fftw_static.sh`, static for the same
  one-less-shared-object reason as OpenBLAS, though it has no collision risk
  of its own (nothing else commonly bundles it).
- **netCDF-C/-Fortran/HDF5**: only needed by the standalone `pamtra` CLI
  executable, not by `pyPamtraLib` (`import pyPamtra` gets its NetCDF I/O
  through the pure-Python `netCDF4` package instead) — so wheel builds pass
  `-Dbuild_cli=false` (`meson_options.txt`) and skip this dependency chain
  entirely rather than bundle it. `tools/build_netcdf_stack.sh` still exists
  for `pip install .`/conda builds, which build the CLI by default.

**One-time setup this depends on**: a
[trusted publisher](https://pypi.org/manage/account/publishing/) registered
on pypi.org for project `pamtra`, GitHub repo `igmk/pamtra`, workflow
`wheels.yml`, environment `pypi` — no API token stored anywhere. Without
this, `publish-pypi` fails at the trusted-publishing handshake even though
the build itself succeeds.
