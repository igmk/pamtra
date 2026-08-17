"""Tests for pyPamtra.pyPamtraXr.

pyPamtraXr is a pure delegation layer over pyPamtra/to_xarray()/
from_xarray()/addInstrument() -- these tests focus on the delegation
being correct (same result as calling the wrapped pam directly), the
escape hatch (.pam) being the real, same object, and that there is no
separate "live" state to get out of sync.
"""

import numpy as np
import pytest

pytest.importorskip("xarray")

import pyPamtra
from _scenario import FREQUENCIES, build_pamtra

HYDROMETEOR_KWARGS = dict(
    hydro_name="rwc_q", liq_ice=1, moment_in=3, nbin=30,
    dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
    scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
)


def build_pamxr():
    pamxr = pyPamtra.pyPamtraXr()
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)

    reference = build_pamtra()  # legacy pam, only used here to source a valid profile
    kwargs = {
        k: reference.p[k]
        for k in ["hgt_lev", "temp_lev", "press_lev", "relhum_lev"]
    }
    pamxr.set_profile(**kwargs)
    pamxr.pam.p["hydro_q"][0, 0, 1, 0] = 1e-3
    return pamxr


def test_construction_creates_fresh_pam():
    pamxr = pyPamtra.pyPamtraXr()
    assert isinstance(pamxr.pam, pyPamtra.pyPamtra)


def test_construction_wraps_existing_pam():
    existing = pyPamtra.pyPamtra()
    pamxr = pyPamtra.pyPamtraXr(pam=existing)
    assert pamxr.pam is existing


def test_add_and_remove_hydrometeor():
    pamxr = pyPamtra.pyPamtraXr()
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    assert pamxr.pam.df.nhydro == 1
    assert pamxr.pam.df.data["hydro_name"][0] == "rwc_q"

    pamxr.remove_hydrometeor("rwc_q")
    assert pamxr.pam.df.nhydro == 0


def test_set_profile_matches_direct_createProfile():
    pamxr = build_pamxr()
    reference = build_pamtra()
    np.testing.assert_allclose(pamxr.pam.p["hgt_lev"], reference.p["hgt_lev"])
    np.testing.assert_allclose(pamxr.pam.p["temp_lev"], reference.p["temp_lev"])


def test_set_profile_from_xarray_round_trip():
    pamxr = build_pamxr()
    ds = pamxr.to_xarray(source="p")

    pamxr2 = pyPamtra.pyPamtraXr()
    pamxr2.add_hydrometeor(**HYDROMETEOR_KWARGS)
    pamxr2.set_profile_from_xarray(ds)

    np.testing.assert_allclose(pamxr2.pam.p["hgt_lev"], pamxr.pam.p["hgt_lev"])
    np.testing.assert_allclose(pamxr2.pam.p["hydro_q"], pamxr.pam.p["hydro_q"])


def test_run_returns_dataset_and_populates_pam_r():
    pamxr = build_pamxr()
    results = pamxr.run(FREQUENCIES)

    assert "tb" in results
    assert "Ze" in results
    assert "tb" in pamxr.pam.r
    np.testing.assert_allclose(results["tb"].values, pamxr.pam.r["tb"])


def test_add_instrument_and_instruments_property():
    pamxr = build_pamxr()
    instrument = pamxr.add_instrument(
        pyPamtra.PamtraInstrument("simple", FREQUENCIES[0], radar_mode="simple")
    )

    assert pamxr.instruments["simple"] is instrument
    assert instrument.results is not None
    assert "Ze" in instrument.results


def test_to_netcdf_writes_results_by_default(tmp_path):
    pamxr = build_pamxr()
    pamxr.run(FREQUENCIES)

    out_file = tmp_path / "results.nc"
    pamxr.to_netcdf(str(out_file))

    assert out_file.exists()
    assert out_file.stat().st_size > 0


def test_to_netcdf_can_write_profile_instead(tmp_path):
    pamxr = build_pamxr()

    out_file = tmp_path / "profile.nc"
    pamxr.to_netcdf(str(out_file), source="p")

    assert out_file.exists()
    assert out_file.stat().st_size > 0


def test_escape_hatch_mutation_is_visible_through_pamxr():
    pamxr = build_pamxr()
    pamxr.pam.p["hydro_q"][0, 0, 0, 0] = 7e-3
    assert pamxr.pam.p["hydro_q"][0, 0, 0, 0] == 7e-3
