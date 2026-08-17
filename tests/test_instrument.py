"""Tests for pyPamtra.PamtraInstrument / pyPamtra.addInstrument().

PamtraInstrument is a thin convenience wrapper around the existing
runPamtra()/to_xarray(): it temporarily overrides pam.nmlSet, calls
runPamtra(), and snapshots pam.r into its own .results. The tests here
focus on exactly the things that wrapper has to get right: nmlSet is
restored (even on failure), each instrument's results are independent
of the others sharing the same pyPamtra object, and the error paths are
clear.
"""

import numpy as np
import pytest

pytest.importorskip("xarray")

import pyPamtra
from _scenario import build_pamtra


def test_addInstrument_runs_and_populates_results():
    pam = build_pamtra()
    instrument = pyPamtra.PamtraInstrument("simple", [10.65, 35.5], radar_mode="simple")

    returned = pam.addInstrument(instrument)

    assert returned is instrument
    assert pam.instruments["simple"] is instrument
    assert instrument.results is not None
    assert "Ze" in instrument.results
    np.testing.assert_allclose(instrument.results["Ze"].values, pam.r["Ze"])


def test_addInstrument_run_false_defers_execution():
    pam = build_pamtra()
    instrument = pyPamtra.PamtraInstrument("simple", [35.5], radar_mode="simple")

    pam.addInstrument(instrument, run=False)
    assert instrument.results is None

    instrument.run()
    assert instrument.results is not None


def test_two_instruments_share_profile_but_not_results():
    pam = build_pamtra()
    simple = pam.addInstrument(pyPamtra.PamtraInstrument("simple", [35.5], radar_mode="simple"))
    spectrum = pam.addInstrument(pyPamtra.PamtraInstrument("spectrum", [35.5], radar_mode="spectrum"))

    assert "radar_spectra" not in simple.results or simple.results["radar_spectra"].size <= 1
    assert "radar_spectra" in spectrum.results
    assert spectrum.results["radar_spectra"].size > 1
    assert simple.results is not spectrum.results


def test_nmlSet_restored_after_run():
    pam = build_pamtra()
    original = dict(pam.nmlSet)

    pam.addInstrument(pyPamtra.PamtraInstrument("simple", [35.5], radar_mode="simple"))

    assert pam.nmlSet == original


def test_nmlSet_restored_after_failed_run():
    pam = build_pamtra()
    original = dict(pam.nmlSet)

    with pytest.raises(TypeError, match="unknown nmlSet override"):
        pam.addInstrument(pyPamtra.PamtraInstrument("bad", [35.5], not_a_real_setting=True))

    assert pam.nmlSet == original
    assert "bad" not in pam.instruments


def test_run_without_parent_raises():
    instrument = pyPamtra.PamtraInstrument("orphan", [35.5])
    with pytest.raises(RuntimeError, match="no parent"):
        instrument.run()


def test_to_netcdf_before_run_raises():
    instrument = pyPamtra.PamtraInstrument("norun", [35.5])
    with pytest.raises(RuntimeError, match="no results yet"):
        instrument.to_netcdf("unused.nc")


def test_to_netcdf_writes_file(tmp_path):
    pam = build_pamtra()
    instrument = pam.addInstrument(pyPamtra.PamtraInstrument("simple", [35.5], radar_mode="simple"))

    out_file = tmp_path / "simple.nc"
    instrument.to_netcdf(str(out_file))

    assert out_file.exists()
    assert out_file.stat().st_size > 0


def test_scalar_frequency_is_wrapped_in_a_list():
    instrument = pyPamtra.PamtraInstrument("simple", 35.5)
    assert instrument.frequencies == [35.5]
