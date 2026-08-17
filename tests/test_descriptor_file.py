import tempfile
from pathlib import Path

import numpy as np
import pytest

import pyPamtra
from pyPamtra.descriptorFile import (
    riming_dependent_area_size,
    riming_dependent_mass_size,
    riming_dependent_ssrga_name,
    ssrga_parameter,
)

REPO_ROOT = Path(__file__).parent.parent
DESCRIPTOR_FILE = REPO_ROOT / "descriptorfiles" / "descriptor_file_2m_icon.txt"

HYDROMETEOR = (
    "rwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
    3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
    "mie-sphere", "khvorostyanov01_drops", -99.0,
)


def make_pam_with_hydrometeor():
    # add4D/remove4D/removeHydrometeor all key off self.parent.p's grid
    # shape (set via pam.p directly here -- no need to build a full
    # atmospheric profile just to exercise the descriptor-file side).
    pam = pyPamtra.pyPamtra()
    pam.p["ngridx"] = 1
    pam.p["ngridy"] = 1
    pam.p["max_nlyrs"] = 2
    pam.df.addHydrometeor(HYDROMETEOR)
    return pam


def test_read_file():
    pam = pyPamtra.pyPamtra()
    pam.df.readFile(str(DESCRIPTOR_FILE))
    assert pam.df.nhydro == 6
    assert list(pam.df.data["hydro_name"]) == [
        "cwc_q", "iwc_q", "rwc_q", "swc_q", "gwc_q", "hwc_q",
    ]


def test_write_read_roundtrip(tmp_path):
    pam = pyPamtra.pyPamtra()
    pam.df.readFile(str(DESCRIPTOR_FILE))

    out_file = tmp_path / "roundtrip.txt"
    pam.df.writeFile(str(out_file))

    pam2 = pyPamtra.pyPamtra()
    pam2.df.readFile(str(out_file))

    assert (pam.df.data == pam2.df.data).all()


def test_add_hydrometeor():
    pam = pyPamtra.pyPamtra()
    pam.df.addHydrometeor((
        "rwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
        3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
        "mie-sphere", "khvorostyanov01_drops", -99.0,
    ))
    assert pam.df.nhydro == 1
    assert pam.df.data["hydro_name"][0] == "rwc_q"


def test_add_hydrometeor_kwargs_matches_tuple():
    # kwargs form should produce byte-identical data to the equivalent
    # positional tuple, including the defaulted (-99.) fields
    pam_tuple = pyPamtra.pyPamtra()
    pam_tuple.df.addHydrometeor(HYDROMETEOR)

    pam_kwargs = pyPamtra.pyPamtra()
    pam_kwargs.df.addHydrometeor(
        hydro_name="rwc_q", liq_ice=1, moment_in=3, nbin=30,
        dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
        scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
    )

    assert (pam_tuple.df.data == pam_kwargs.df.data).all()


def test_add_hydrometeor_kwargs_unknown_field_raises():
    pam = pyPamtra.pyPamtra()
    with pytest.raises(TypeError, match="unexpected keyword"):
        pam.df.addHydrometeor(
            hydro_name="rwc_q", liq_ice=1, moment_in=3, nbin=30,
            dist_name="exp", d_1=1.0e-5, scat_name="mie-sphere",
            vel_size_mod="khvorostyanov01_drops",
            not_a_real_field=1.0,
        )


def test_add_hydrometeor_kwargs_missing_required_raises():
    pam = pyPamtra.pyPamtra()
    with pytest.raises(TypeError, match="missing required"):
        pam.df.addHydrometeor(hydro_name="rwc_q", liq_ice=1)


@pytest.mark.parametrize("field,value", [("liq_ice", 0), ("moment_in", 99)])
def test_add_hydrometeor_kwargs_invalid_enum_raises(field, value):
    pam = pyPamtra.pyPamtra()
    kwargs = dict(
        hydro_name="rwc_q", liq_ice=1, moment_in=3, nbin=30,
        dist_name="exp", d_1=1.0e-5, scat_name="mie-sphere",
        vel_size_mod="khvorostyanov01_drops",
    )
    kwargs[field] = value
    with pytest.raises(ValueError):
        pam.df.addHydrometeor(**kwargs)


def test_add_hydrometeor_tuple_and_kwargs_together_raises():
    pam = pyPamtra.pyPamtra()
    with pytest.raises(TypeError, match="either hydroTuple or keyword"):
        pam.df.addHydrometeor(HYDROMETEOR, hydro_name="something_else")


def test_remove_hydrometeor():
    pam = make_pam_with_hydrometeor()
    pam.df.removeHydrometeor("rwc_q")
    assert pam.df.nhydro == 0
    assert len(pam.df.data) == 0


def test_remove_hydrometeor_missing_name_raises():
    pam = make_pam_with_hydrometeor()
    with pytest.raises(ValueError):
        pam.df.removeHydrometeor("does_not_exist")


def test_remove_hydrometeor_all():
    pam = make_pam_with_hydrometeor()
    pam.df.removeHydrometeor("all")
    assert pam.df.nhydro == 0
    assert len(pam.df.data) == 0


def test_add4d_remove4d_roundtrip():
    pam = make_pam_with_hydrometeor()
    assert "a_ms" in pam.df.data.dtype.names

    arr = np.arange(np.prod(pam._shape4D), dtype=float).reshape(pam._shape4D)
    pam.df.add4D("a_ms", arr)
    assert "a_ms" not in pam.df.data.dtype.names
    np.testing.assert_array_equal(pam.df.data4D["a_ms"], arr)

    pam.df.remove4D("a_ms", np.array([-99.0]))
    assert "a_ms" in pam.df.data.dtype.names
    np.testing.assert_allclose(pam.df.data["a_ms"], [-99.0])
    assert "a_ms" not in pam.df.data4D


def test_add_full_spectra():
    pam = make_pam_with_hydrometeor()  # nbin=30 in HYDROMETEOR
    pam.df.addFullSpectra()

    assert pam.df.fs_nbin == 30
    for key in ["rho_ds", "d_ds", "n_ds", "mass_ds", "area_ds",
                "as_ratio", "canting", "fallvelocity"]:
        assert pam.df.dataFullSpec[key].shape == pam._shape4D + (30,)
    # bin boundaries: one more than bin centers
    assert pam.df.dataFullSpec["d_bound_ds"].shape == pam._shape4D + (31,)


def test_ssrga_parameter():
    kappa, beta, gamma, zeta1, alpha_eff = ssrga_parameter(0.0)
    # M=0 should reduce to the p3 (constant) coefficients
    np.testing.assert_allclose([alpha_eff, kappa, beta, gamma, zeta1],
                                [0.575, 0.194, 5.42, 2.76, 0.067])


def test_riming_dependent_ssrga_name():
    name = riming_dependent_ssrga_name(0.5)
    assert name.startswith("ss-rayleigh-gans_")
    assert len(name.split("_")) == 5


def test_riming_dependent_ssrga_name_invalid_input_raises():
    with pytest.raises(ValueError):
        riming_dependent_ssrga_name("not a number")


@pytest.mark.parametrize("monomer", ["column", "dendrite", "needle", "plate", "rosette", "mean"])
def test_riming_dependent_mass_size(monomer):
    a_m, b_m = riming_dependent_mass_size(0.05, monomer)
    assert np.isfinite(a_m)
    assert np.isfinite(b_m)


def test_riming_dependent_mass_size_clamps_above_range():
    # M beyond the table's last entry (0.8155) is clamped to it, not extrapolated
    at_max, _ = riming_dependent_mass_size(0.8155, "column")
    beyond_max, _ = riming_dependent_mass_size(10.0, "column")
    np.testing.assert_allclose(beyond_max, at_max)


def test_riming_dependent_mass_size_invalid_monomer_raises():
    with pytest.raises(ValueError):
        riming_dependent_mass_size(0.05, "not_a_monomer")


@pytest.mark.parametrize("monomer", ["column", "dendrite", "needle", "plate", "rosette", "mean"])
def test_riming_dependent_area_size(monomer):
    a_A, b_A = riming_dependent_area_size(0.05, monomer)
    assert np.isfinite(a_A)
    assert np.isfinite(b_A)


def test_riming_dependent_area_size_invalid_monomer_raises():
    with pytest.raises(ValueError):
        riming_dependent_area_size(0.05, "not_a_monomer")
