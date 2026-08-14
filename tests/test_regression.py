"""Golden-output regression tests.

The Fortran core is treated as a black box here: run a small, fixed
scenario end-to-end and compare selected outputs against values stored in
tests/golden/. This is the primary safety net against silent numerical
regressions in the RT/scattering code, since the Fortran core is built
around shared global module state rather than pure functions and isn't a
realistic target for per-subroutine unit tests (see AI.md).

To intentionally update the reference values after a real physics change,
run: python tests/generate_golden_data.py
"""

from pathlib import Path

import numpy as np
import pytest

from _scenario import run_scenario

GOLDEN_DIR = Path(__file__).parent / "golden"


@pytest.fixture(scope="module")
def result():
    return run_scenario()


@pytest.fixture(scope="module")
def reference():
    return np.load(GOLDEN_DIR / "rain_snow_reference.npz")


def test_radar_reflectivity(result, reference):
    np.testing.assert_allclose(result.r["Ze"], reference["Ze"], rtol=1e-5)


def test_hydrometeor_attenuation(result, reference):
    np.testing.assert_allclose(result.r["Att_hydro"], reference["Att_hydro"], rtol=1e-5)


def test_brightness_temperature(result, reference):
    np.testing.assert_allclose(result.r["tb"], reference["tb"], rtol=1e-5)


def test_no_missing_values(result):
    # layer index 1 (rain) and 2 (snow) both carry hydrometeor mass and
    # should produce a real (non -9999 fill value) reflectivity.
    assert (result.r["Ze"][0, 0, 1:3] > -9999).all()
