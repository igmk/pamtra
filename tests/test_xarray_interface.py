"""Tests for pyPamtra.to_xarray()/from_xarray().

These are purely additive on top of self.p/self.r/self.df.data -- they
never read a live view of them and never write to them, so the tests
here focus on: (1) round-trip fidelity through from_xarray(), (2) that
mutating the returned Dataset can never affect the pyPamtra object it
came from, and (3) the optional-dependency contract.
"""

import sys

import numpy as np
import pytest

xr = pytest.importorskip("xarray")

import pyPamtra
from _scenario import FREQUENCIES, build_pamtra, run_scenario

P_KEYS_TO_CHECK = [
    "hgt_lev", "temp_lev", "press_lev", "relhum_lev",
    "hydro_q", "lat", "lon", "unixtime", "sfc_type", "obs_height",
]


def test_to_xarray_missing_xarray_gives_clear_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "xarray", None)
    pam = build_pamtra()
    with pytest.raises(ImportError, match="optional 'xarray' package"):
        pam.to_xarray()


def test_to_xarray_invalid_source_raises():
    pam = build_pamtra()
    with pytest.raises(ValueError, match="source must be"):
        pam.to_xarray(source="bogus")


def test_to_xarray_p_round_trip_via_from_xarray():
    pam = build_pamtra()
    ds = pam.to_xarray(source="p")

    pam2 = pyPamtra.pyPamtra()
    pam2.df = pam.df  # same hydrometeors/order as the source profile
    pam2.from_xarray(ds)

    for key in P_KEYS_TO_CHECK:
        np.testing.assert_allclose(pam2.p[key], pam.p[key], err_msg=key)


def test_to_xarray_p_units_and_dims():
    pam = build_pamtra()
    ds = pam.to_xarray(source="p")

    assert ds["hgt_lev"].attrs["units"] == "m"
    assert ds["hgt_lev"].dims == ("grid_x", "grid_y", "heightbins_plus1")
    assert ds["hydro_q"].dims == ("grid_x", "grid_y", "heightbins", "hydrometeor")
    assert ds.sizes["hydrometeor"] == pam.df.nhydro


def test_to_xarray_mutation_does_not_affect_pam():
    pam = build_pamtra()
    original = pam.p["temp_lev"].copy()

    ds = pam.to_xarray(source="p")
    ds["temp_lev"].values[:] = -1

    np.testing.assert_allclose(pam.p["temp_lev"], original)


def test_to_xarray_results_after_run():
    pam = run_scenario()
    ds = pam.to_xarray(source="r")

    assert "tb" in ds
    assert "Ze" in ds
    assert ds.attrs["pamtraVersion"] == pam.r["pamtraVersion"]
    np.testing.assert_allclose(ds["tb"].values, pam.r["tb"])


def test_to_xarray_results_side_keys_get_real_dims_not_synthesized():
    # self.dimensions only covered Ze/Att_hydro/Att_atmo/tb on the results
    # side; everything else (angles_deg, emissivity, radar_snr/quality/
    # moments/slopes/edges/vel, psd_*) fell back to synthesized
    # key_dim0/key_dim1/... names
    pam = run_scenario()
    ds = pam.to_xarray(source="r")

    # emissivity's angle axis is half the length of tb's/angles_deg's --
    # a genuinely different set, so it must not share their "angles" dim
    assert ds["emissivity"].dims == ("grid_x", "grid_y", "passive_polarisation", "frequency", "emis_angles")
    assert ds.sizes["emis_angles"] != ds.sizes["angles"]


def test_to_xarray_promotes_angles_deg_and_freqs_to_coordinates():
    # angles_deg/self.set["freqs"] are the real values the "angles"/
    # "frequency" dimensions index into -- they should be genuine xarray
    # coordinates (named to match their dimension), not separately-named
    # data variables, so .sel(angles=..., frequency=...) works
    pam = run_scenario()
    ds = pam.to_xarray(source="r")

    assert "angles_deg" not in ds.data_vars
    assert "angles" in ds.coords
    np.testing.assert_allclose(ds.coords["angles"].values, pam.r["angles_deg"])

    assert "frequency" in ds.coords
    np.testing.assert_allclose(ds.coords["frequency"].values, pam.set["freqs"])

    np.testing.assert_allclose(
        ds["tb"].sel(angles=180.0, frequency=pam.set["freqs"][0]).values,
        pam.r["tb"][:, :, :, 0, 0, :],
    )


def test_to_xarray_radar_spectra_has_all_six_dims():
    # self.dimensions["radar_spectra"] declared only 5 of its 6 actual
    # dims when radar_mode == "spectrum", so it silently fell back to
    # synthesized names even when correctly shaped, real data was present
    pam = build_pamtra()
    pam.nmlSet["radar_mode"] = "spectrum"
    pam.runPamtra(FREQUENCIES)
    ds = pam.to_xarray(source="r")

    assert ds["radar_spectra"].dims == (
        "grid_x", "grid_y", "heightbins", "frequency", "radar_polarisation", "radar_nfft",
    )
    assert ds["radar_vel"].dims == ("frequency", "radar_nfft")


def test_to_xarray_save_ssp_and_save_psd_do_not_collide_dim_names():
    # scatter_matrix/extinct_matrix/emis_vector each reuse the same
    # "stokes"/"angle" axis size twice (in/out) -- xarray raises if two
    # axes of one variable share a dim name, so these need genuinely
    # distinct names, not just correct lengths
    pam = build_pamtra()
    pam.nmlSet["save_ssp"] = True
    pam.nmlSet["save_psd"] = True
    pam.runPamtra(FREQUENCIES)
    ds = pam.to_xarray(source="r")  # raises ValueError on a dim-name collision

    assert len(set(ds["scatter_matrix"].dims)) == len(ds["scatter_matrix"].dims)
    assert len(set(ds["extinct_matrix"].dims)) == len(ds["extinct_matrix"].dims)
    assert ds["psd_d"].dims == ("grid_x", "grid_y", "heightbins", "hydrometeor", "sizebin")
    assert ds["psd_d"].attrs["units"] == "m"


def test_to_xarray_previously_ungapped_keys_get_real_dims_and_units():
    # groundtemp/airturb/hgt/press_lev previously fell back to synthesized
    # key_dim0/key_dim1/... names because self.dimensions had no entry for
    # them (or, for press_lev, had one under the wrong key: "p_lev", which
    # self.p never actually uses)
    pam = build_pamtra()
    ds = pam.to_xarray(source="p")

    assert ds["groundtemp"].dims == ("grid_x", "grid_y")
    assert ds["groundtemp"].attrs["units"] == "K"
    assert ds["airturb"].dims == ("grid_x", "grid_y", "heightbins")
    assert ds["hgt"].dims == ("grid_x", "grid_y", "heightbins")
    assert ds["press_lev"].dims == ("grid_x", "grid_y", "heightbins_plus1")
    assert ds["press_lev"].attrs["units"] == "Pa"


def test_to_xarray_sanitizes_non_identifier_dim_names():
    # self.dimensions["radar_prop"] declares a literal "2" as its last
    # dim name, which is not a valid dimension identifier
    pam = build_pamtra()
    ds = pam.to_xarray(source="p")

    assert "2" not in ds["radar_prop"].dims
    assert all(d.isidentifier() for d in ds["radar_prop"].dims)


def test_to_xarray_outer_dims_round_trip():
    pam = build_pamtra()
    ds = pam.to_xarray(source="p", outer_dims={"grid_x": "lat", "grid_y": "lon_idx"})
    assert "lat" in ds["hgt_lev"].dims
    assert "grid_x" not in ds.dims

    pam2 = pyPamtra.pyPamtra()
    pam2.df = pam.df
    pam2.from_xarray(ds, outer_dims={"lat": "grid_x", "lon_idx": "grid_y"})

    np.testing.assert_allclose(pam2.p["hgt_lev"], pam.p["hgt_lev"])
