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

## 4a. First-ever conda-forge submission (one-time, done)

Already completed — `conda-forge/pamtra-feedstock` exists. Left here for
reference in case the feedstock ever needs to be recreated from scratch:

1. Fork `conda-forge/staged-recipes`, branch from `main`.
2. Add the recipe at `recipes/pamtra/recipe.yaml` (their required layout).
3. Open a PR there, following their PR checklist, and request review from
   `@conda-forge/help-python` (or whichever team their checklist points to).
4. Once merged, conda-forge's automation creates `conda-forge/pamtra-feedstock`
   and does the first build/publish. Maintainers listed in the recipe's
   `extra.recipe-maintainers` get commit rights on that feedstock repo.

## 4b. Every subsequent release

Now that PyPI releases are live (see "PyPI" below), conda-forge's autotick
bot (`regro-cf-autotick-bot`) notices the new PyPI release on its own and
opens a version-bump PR against `conda-forge/pamtra-feedstock` — no manual
PR needed to get one started. But the bot only bumps `recipe/recipe.yaml`'s
`version`/`sha256`/`build.number`; it does **not** diff against whatever
changed in [conda-recipe/recipe.yaml](conda-recipe/recipe.yaml) since the
last feedstock sync. Concretely, this bit us going from 1.0.3 to 1.1.0: the
bot's PR built fine but failed at test time with
`ModuleNotFoundError: No module named 'meteo_si'`, because `meteo_si`
had been added as a run dependency in the interim and the bot had no way to
know that.

1. **Diff `conda-recipe/recipe.yaml` against the bot PR's
   `recipe/recipe.yaml`.** Apply any dependency added/removed/changed in the
   former to the latter by hand — this is the step that's easy to skip, and
   the failure mode (an import error, or an unsolvable environment) doesn't
   point back at "the recipe is out of sync" on its own.
2. Push the fix straight to the bot's PR branch — maintainers listed in
   `extra.recipe-maintainers` have push access to
   `regro-cf-autotick-bot/pamtra-feedstock`, same as any fork you have
   commit rights on:
   ```bash
   git clone --branch <bot-branch-name> https://github.com/regro-cf-autotick-bot/pamtra-feedstock.git
   cd pamtra-feedstock
   # edit recipe/recipe.yaml to match conda-recipe/recipe.yaml
   git commit -am "..."
   git push origin <bot-branch-name>
   ```
3. Once the feedstock's own CI (linux-64/osx-64/win-64) is green, merge —
   conda-forge publishes the new version automatically.

If no bot PR shows up within a few days of the PyPI release, fall back to
opening one manually: clone/update your fork of `pamtra-feedstock`, apply
the same `recipe/recipe.yaml` edits, and open the PR yourself.

Separately, conda-forge sometimes opens **migration PRs** (e.g. a new
Python version) on their own bot branch, independent of version bumps —
these only touch `.ci_support/*.yaml`/`.azure-pipelines/*.yml`, not
`recipe.yaml`. If one is open when a version-bump PR merges, it'll still be
testing the old version; merge `main` into the migration PR's branch (same
clone-push pattern as above, plus `git remote add upstream
https://github.com/conda-forge/pamtra-feedstock.git && git fetch upstream
main && git merge upstream/main`) to pick up the bump before merging it —
the two touch disjoint files, so this is conflict-free in practice.

## Verification checklist

- [ ] `rattler-build build --recipe conda-recipe/recipe.yaml -c conda-forge`
      succeeds locally against the real tagged tarball (not the local `path:`
      source used for day-to-day recipe edits)
- [ ] `conda-recipe/recipe.yaml` and the feedstock PR's `recipe/recipe.yaml`
      actually match (see 4b — the bot doesn't sync dependency changes)
- [ ] The feedstock PR's own CI (linux-64/osx-64/win-64) is green before
      merging
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
