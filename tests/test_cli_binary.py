"""Run the standalone pamtra CLI binary end-to-end and check it against
the pyPamtra Python API, on the same fixed, no-external-data scenario
used by the golden regression tests (see _scenario.py).

Writing this test is what turned up the writePamtraProfile bugs fixed
alongside it (see git history) -- notably a header field (the
per-gridpoint layer/level count) that undercounted .lev files by one,
silently dropping the top atmospheric layer from the Fortran read. With
that fixed, the CLI and Python outputs agree closely (limited mainly by
the ASCII profile format's precision, e.g. temperature/pressure are only
written to 2/1 decimal places), so this compares them at a fairly tight
tolerance rather than just smoke-testing that the binary runs. The CLI's
netCDF output does still carry `noutlevels` worth of extra padding in its
layer dimension that the Python API's in-memory result does not; the
comparisons below account for that (matching layer counts, or comparing
only the non-fill-value entries) rather than expecting identical shapes
throughout.
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
    # Limited mainly by the ASCII profile format's precision (temperature
    # is only written to 2 decimal places, pressure to 1), not by any
    # remaining difference between the two entry points.
    np.testing.assert_allclose(cli_tb, py_tb, rtol=1e-3)


def test_radar_reflectivity_matches_python(cli_result, python_result):
    cli_ze = np.ma.masked_equal(cli_result.variables["Ze"][:], -9999)
    py_ze = np.ma.masked_equal(python_result.r["Ze"], -9999)
    # Both hydrometeor layers in this scenario should come back with a Ze
    # value from both paths, in the same order -- the CLI's array is just
    # the Python one plus a few extra always-masked noutlevels slots.
    assert py_ze.count() == 2
    assert cli_ze.count() == py_ze.count()
    np.testing.assert_allclose(cli_ze.compressed(), py_ze.compressed(), rtol=1e-3)
