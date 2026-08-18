"""Tests for pyPamtra.pyPamtraXr.

.p/.r/.df/.df_4d/.df_full_spec are real, persistent attributes here (not
hidden/export-on-demand): .p/.r are xarray.Dataset, .df is a
pandas.DataFrame indexed by hydro_name, .df_4d/.df_full_spec are
xarray.Dataset. self.pam is the one place the RT engine actually runs;
these tests check that translation is correct in both directions
(pam -> pamxr via from_pam()/refresh(), pamxr -> pam via run()), and that
self.df stays eagerly in sync with self.pam.df on every hydrometeor call
(reusing pamDescriptorFile's own tested methods, not a reimplementation).
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


def build_pamxr(outer_dims=None):
    pamxr = pyPamtra.pyPamtraXr(outer_dims=outer_dims)
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)

    reference = build_pamtra()  # legacy pam, only used here to source a valid profile
    kwargs = {
        k: reference.p[k]
        for k in ["hgt_lev", "temp_lev", "press_lev", "relhum_lev"]
    }
    pamxr.set_profile(**kwargs)
    pamxr.pam.p["hydro_q"][0, 0, 1, 0] = 1e-3
    pamxr.refresh()
    return pamxr


def test_construction_creates_fresh_empty_state():
    pamxr = pyPamtra.pyPamtraXr()
    assert isinstance(pamxr.pam, pyPamtra.pyPamtra)
    assert len(pamxr.p.data_vars) == 0
    assert len(pamxr.r.data_vars) == 0
    assert len(pamxr.df) == 0
    assert len(pamxr.df_4d.data_vars) == 0
    assert len(pamxr.df_full_spec.data_vars) == 0


def test_construction_wraps_existing_pam():
    existing = pyPamtra.pyPamtra()
    pamxr = pyPamtra.pyPamtraXr(pam=existing)
    assert pamxr.pam is existing


def test_from_pam_bridges_importer_output():
    reference = build_pamtra()  # built the classic way, like any pyPamtra.importer function would
    pamxr = pyPamtra.pyPamtraXr.from_pam(reference)

    assert pamxr.pam is reference
    assert len(pamxr.df) == 2
    np.testing.assert_allclose(pamxr.p["hgt_lev"].values, reference.p["hgt_lev"])


def test_add_hydrometeor_updates_df_and_pam_df():
    pamxr = pyPamtra.pyPamtraXr()
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)

    assert list(pamxr.df.index) == ["rwc_q"]
    assert pamxr.df.loc["rwc_q", "dist_name"] == "exp"
    assert pamxr.pam.df.nhydro == 1
    assert pamxr.pam.df.data["hydro_name"][0] == "rwc_q"


def test_remove_hydrometeor():
    pamxr = pyPamtra.pyPamtraXr()
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    pamxr.remove_hydrometeor("rwc_q")

    assert len(pamxr.df) == 0
    assert pamxr.pam.df.nhydro == 0


def test_add_4d_and_remove_4d():
    pamxr = build_pamxr()

    arr = np.zeros(pamxr.pam._shape4D)
    arr[:] = 0.5
    pamxr.add_4d("a_ms", arr)

    assert "a_ms" in pamxr.df_4d
    assert "a_ms" not in pamxr.df.columns
    assert pamxr.df_4d["a_ms"].dims == ("grid_x", "grid_y", "heightbins", "hydrometeor")

    pamxr.remove_4d("a_ms", np.array([-99.0]))
    assert "a_ms" not in pamxr.df_4d
    assert "a_ms" in pamxr.df.columns


def test_add_full_spectra():
    pamxr = build_pamxr()
    pamxr.add_full_spectra()

    assert "d_ds" in pamxr.df_full_spec
    assert pamxr.df_full_spec["d_ds"].dims == ("grid_x", "grid_y", "heightbins", "hydrometeor", "sizebin")
    # d_bound_ds has one more bin than the rest: distinct dim, not just distinct size
    assert pamxr.df_full_spec["d_bound_ds"].dims[-1] == "sizebin_plus1"


def test_set_profile_matches_direct_createProfile():
    pamxr = build_pamxr()
    reference = build_pamtra()
    np.testing.assert_allclose(pamxr.p["hgt_lev"].values, reference.p["hgt_lev"])
    np.testing.assert_allclose(pamxr.p["temp_lev"].values, reference.p["temp_lev"])


def test_set_profile_from_xarray_round_trip():
    pamxr = build_pamxr()
    ds = pamxr.p

    pamxr2 = pyPamtra.pyPamtraXr()
    pamxr2.add_hydrometeor(**HYDROMETEOR_KWARGS)
    pamxr2.set_profile_from_xarray(ds)

    np.testing.assert_allclose(pamxr2.p["hgt_lev"].values, pamxr.p["hgt_lev"].values)
    np.testing.assert_allclose(pamxr2.p["hydro_q"].values, pamxr.p["hydro_q"].values)


def test_run_returns_dataset_and_populates_r():
    pamxr = build_pamxr()
    results = pamxr.run(FREQUENCIES)

    assert "tb" in results
    assert "Ze" in results
    assert results is pamxr.r
    np.testing.assert_allclose(results["tb"].values, pamxr.pam.r["tb"])


def test_run_translates_p_into_pam_before_running():
    pamxr = build_pamxr()
    # mutate .p directly (not through set_profile), then run(): pam.p must
    # reflect this, not whatever pam.p happened to hold before
    pamxr.p["hydro_q"].values[0, 0, 1, 0] = 9e-3

    pamxr.run(FREQUENCIES)

    np.testing.assert_allclose(pamxr.pam.p["hydro_q"][0, 0, 1, 0], 9e-3)


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


def test_to_netcdf_invalid_source_raises():
    pamxr = build_pamxr()
    with pytest.raises(ValueError, match="source must be"):
        pamxr.to_netcdf("unused.nc", source="bogus")


def test_refresh_picks_up_escape_hatch_mutation():
    pamxr = build_pamxr()

    pamxr.pam.df.addHydrometeor(
        hydro_name="swc_q", liq_ice=-1, moment_in=3, nbin=30,
        dist_name="exp_field_t", d_1=0.51e-10, d_2=2.0e-2,
        scat_name="mie-sphere", vel_size_mod="heymsfield10_particles",
    )
    assert len(pamxr.df) == 1  # not yet reflected

    pamxr.refresh()
    assert len(pamxr.df) == 2
    assert "swc_q" in pamxr.df.index


def test_escape_hatch_is_same_object():
    pamxr = build_pamxr()
    pamxr.pam.p["hydro_q"][0, 0, 0, 0] = 7e-3
    assert pamxr.pam.p["hydro_q"][0, 0, 0, 0] == 7e-3


def test_repr_reflects_state():
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["lat", "lon"])
    assert "no hydrometeors" in repr(pamxr)
    assert "no profile set" in repr(pamxr)
    assert "no results yet" in repr(pamxr)
    assert "outer_dims=['lat', 'lon']" in repr(pamxr)

    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    assert "1 hydrometeor (rwc_q)" in repr(pamxr)

    reference = build_pamtra()
    pamxr.set_profile(**{k: reference.p[k] for k in ["hgt_lev", "temp_lev", "press_lev", "relhum_lev"]})
    assert "no profile set" not in repr(pamxr)

    pamxr.run(FREQUENCIES)
    assert "results populated" in repr(pamxr)


def test_outer_dims_rejects_dict():
    # list(a_dict) silently returns its keys -- reject explicitly rather
    # than accepting the pre-N-D {"grid_x": "lat", ...} rename dict as a
    # (wrong) list of names
    with pytest.raises(TypeError, match="not a dict"):
        pyPamtra.pyPamtraXr(outer_dims={"grid_x": "lat", "grid_y": "lon"})


def test_outer_dims_two_names_relabels_p_and_r():
    pamxr = build_pamxr(outer_dims=["lat", "lon"])
    assert pamxr.p["hgt_lev"].dims == ("lat", "lon", "heightbins_plus1")

    results = pamxr.run(FREQUENCIES)
    assert "lat" in results["tb"].dims
    assert "grid_x" not in results.dims


def test_outer_dims_run_matches_plain_run_numerically():
    plain = build_pamxr()
    renamed = build_pamxr(outer_dims=["lat", "lon"])

    results_plain = plain.run(FREQUENCIES)
    results_renamed = renamed.run(FREQUENCIES)

    np.testing.assert_allclose(results_renamed["tb"].values, results_plain["tb"].values)
    np.testing.assert_allclose(renamed.pam.p["hydro_q"], plain.pam.p["hydro_q"])


def test_outer_dims_applies_to_instrument_results():
    pamxr = build_pamxr(outer_dims=["lat", "lon"])
    instrument = pamxr.add_instrument(
        pyPamtra.PamtraInstrument("simple", FREQUENCIES[0], radar_mode="simple")
    )
    assert "lat" in instrument.results["Ze"].dims


def test_changing_outer_dims_and_refresh_relabels_p():
    pamxr = build_pamxr()
    assert pamxr.p["hgt_lev"].dims == ("grid_x", "grid_y", "heightbins_plus1")

    pamxr.outer_dims = ["time", "grid_y"]
    pamxr.refresh()

    assert pamxr.p["hgt_lev"].dims == ("time", "grid_y", "heightbins_plus1")


def test_outer_dims_single_name():
    # 1 outer dim: a stateless rename + squeeze, no reshape needed since
    # pyPamtra already treats a bare leading axis as ngridx with ngridy=1
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["scan"])
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    reference = build_pamtra()
    pamxr.set_profile(
        hgt_lev=reference.p["hgt_lev"][:, 0],
        temp_lev=reference.p["temp_lev"][:, 0],
        press_lev=reference.p["press_lev"][:, 0],
        relhum_lev=reference.p["relhum_lev"][:, 0],
    )
    assert pamxr.p["hgt_lev"].dims == ("scan", "heightbins_plus1")

    results = pamxr.run(FREQUENCIES)
    assert results["tb"].dims[0] == "scan"
    np.testing.assert_allclose(results["tb"].values, pamxr.pam.r["tb"][:, 0])


def test_outer_dims_three_names_round_trip():
    # 3+ outer dims: a real reshape (pyPamtra's own grid is always 2D),
    # exercising the actual "arbitrary number of dimensions" ask
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["time", "lat", "lon"])
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)

    reference = build_pamtra()  # shape (1, 2, nlyr[+1])
    T, La, Lo = 2, 1, 2
    def tile(key):
        base = reference.p[key][0, 0]  # single profile, shape (nlyr[+1],)
        return np.broadcast_to(base, (T, La, Lo) + base.shape).copy()

    pamxr.set_profile(
        hgt_lev=tile("hgt_lev"), temp_lev=tile("temp_lev"),
        press_lev=tile("press_lev"), relhum_lev=tile("relhum_lev"),
    )
    assert pamxr.p["hgt_lev"].dims == ("time", "lat", "lon", "heightbins_plus1")
    assert pamxr.p["hgt_lev"].shape == (T, La, Lo, reference.p["hgt_lev"].shape[-1])

    # give one (time, lat, lon) gridpoint some liquid water, run, and check
    # that exact gridpoint's result -- confirms the reshape didn't shuffle data
    pamxr.p["hydro_q"].values[1, 0, 0, 0, 0] = 1e-3
    results = pamxr.run(FREQUENCIES)
    assert results["tb"].dims[:3] == ("time", "lat", "lon")
    assert results["tb"].shape[:3] == (T, La, Lo)

    other_pamxr = pyPamtra.pyPamtraXr()
    other_pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    other_pamxr.set_profile(
        hgt_lev=reference.p["hgt_lev"][0, 0], temp_lev=reference.p["temp_lev"][0, 0],
        press_lev=reference.p["press_lev"][0, 0], relhum_lev=reference.p["relhum_lev"][0, 0],
    )
    other_pamxr.p["hydro_q"].values[0, 0, 0, 0] = 1e-3
    single_results = other_pamxr.run(FREQUENCIES)

    np.testing.assert_allclose(
        results["tb"].values[1, 0, 0], single_results["tb"].values[0, 0]
    )


def test_outer_dims_three_names_refresh_without_prior_shape_raises():
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["time", "lat", "lon"])
    pamxr.add_hydrometeor(**HYDROMETEOR_KWARGS)
    # bypass pyPamtraXr entirely: build a profile directly on the escape
    # hatch, so pyPamtraXr never learns the (time, lat, lon) shape
    reference = build_pamtra()
    pamxr.pam.createProfile(
        hgt_lev=reference.p["hgt_lev"], temp_lev=reference.p["temp_lev"],
        press_lev=reference.p["press_lev"], relhum_lev=reference.p["relhum_lev"],
    )
    with pytest.raises(RuntimeError, match="no profile has been set"):
        pamxr.refresh()
