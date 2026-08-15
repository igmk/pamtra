"""Run the standalone pamtra CLI binary end-to-end and sanity-check it
against the pyPamtra Python API, on the same fixed, no-external-data
scenario used by the golden regression tests (see _scenario.py).

The CLI and Python API are two separate entry points into (mostly) the
same Fortran core, driven by two separate "fill settings" paths (see
AI.md), so their outputs are not expected to match exactly -- notably,
the CLI's netCDF output keeps `noutlevels` worth of extra padding in its
layer dimension that the Python API's in-memory result does not, and
this scenario's second hydrometeor layer does not come back with a
Ze value from the CLI (radar_mode 'simple' picks up only the lowest
hydrometeor-bearing layer) even though the Python path can see both.
That gap is a pre-existing quirk of the CLI path, not something this
test tries to fix -- it instead checks the CLI runs, produces a valid
netCDF file, and that where the two do report a value, it agrees.
"""

import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from _scenario import build_pamtra  # noqa: E402

FREQ = 35.5
REPO_ROOT = Path(__file__).parent.parent


def find_pamtra_binary():
    found = shutil.which("pamtra")
    if found is not None:
        return found
    candidate = REPO_ROOT / "build" / "pamtra"
    if candidate.is_file():
        return str(candidate)
    return None


@pytest.fixture(scope="module")
def pamtra_binary():
    binary = find_pamtra_binary()
    if binary is None:
        pytest.skip(
            "pamtra binary not found on PATH or in build/ -- build it with "
            "`pip install .` or `pixi run build-cli` first")
    return binary


@pytest.fixture(scope="module")
def python_result():
    pam = build_pamtra()
    pam.runPamtra([FREQ])
    return pam


@pytest.fixture(scope="module")
def cli_result(pamtra_binary, tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("cli_binary")

    pam = build_pamtra()
    descriptor = tmp_path / "descriptor.txt"
    profile = tmp_path / "profile.lev"
    namelist = tmp_path / "settings.nml"
    pam.df.writeFile(str(descriptor))
    pam.writePamtraProfile(str(profile))
    # Written by hand rather than via pam.writeNmlFile(): that dumps every
    # nmlSet key, including some the Fortran namelist reader rejects.
    # randomseed is the only setting _scenario.build_pamtra() changes from
    # the Fortran defaults.
    namelist.write_text("&settings\nrandomseed=10\n/\n")

    subprocess.run(
        [pamtra_binary, "-f", str(FREQ), "-d", str(descriptor),
         "-p", str(profile), "-n", str(namelist), "-o", str(tmp_path)],
        check=True, capture_output=True, text=True, cwd=tmp_path)

    ncfiles = glob.glob(str(tmp_path / "profile*.nc"))
    assert len(ncfiles) == 1, f"expected one output file, found {ncfiles}"
    return netCDF4.Dataset(ncfiles[0])


def test_cli_runs_and_writes_output(cli_result):
    assert "tb" in cli_result.variables
    assert "Ze" in cli_result.variables


def test_brightness_temperature_matches_python(cli_result, python_result):
    cli_tb = cli_result.variables["tb"][:]
    py_tb = python_result.r["tb"]
    assert cli_tb.shape == py_tb.shape
    assert np.isfinite(cli_tb).all()
    # Loose tolerance: the two entry points fill settings via separate code
    # paths (see module docstring), so small systematic differences are
    # expected -- this is a sanity check, not a bit-exact regression test.
    np.testing.assert_allclose(cli_tb, py_tb, rtol=0.01)


def test_radar_reflectivity_peak_matches_python(cli_result, python_result):
    cli_ze = np.ma.masked_equal(cli_result.variables["Ze"][:], -9999)
    py_ze = np.ma.masked_equal(python_result.r["Ze"], -9999)
    assert cli_ze.count() >= 1
    assert py_ze.count() >= 1
    np.testing.assert_allclose(cli_ze.max(), py_ze.max(), rtol=0.01)
