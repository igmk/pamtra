import numpy as np

import pyPamtra


def test_e_sat_gg_water_at_freezing():
    # Saturation vapor pressure of water at 0 C is a well known reference
    # point (~611 Pa, e.g. Goff-Gratch or Clausius-Clapeyron).
    e_sat = pyPamtra.meteoSI.e_sat_gg_water(273.15)
    np.testing.assert_allclose(e_sat, 610.34, rtol=1e-3)


def test_moist_rho_q_dry_air_ideal_gas():
    # At q=0 (dry air), moist_rho_q must reduce to the ideal gas law:
    # rho = p / (R_specific_dry_air * T).
    rho = pyPamtra.meteoSI.moist_rho_q(101325.0, 273.15, 0.0)
    np.testing.assert_allclose(rho, 1.2923, rtol=1e-3)


def test_rh2q_q2rh_roundtrip():
    rh = 0.5
    q = pyPamtra.meteoSI.rh2q(rh, 280.0, 90000.0)
    rh_back = pyPamtra.meteoSI.q2rh(q, 280.0, 90000.0)
    np.testing.assert_allclose(rh_back, rh, rtol=1e-6)


def test_rh2q_rejects_percent_input():
    # rh must be given as a fraction (0-1), not percent; the function
    # raises rather than silently misinterpreting e.g. rh=50 as 5000%.
    try:
        pyPamtra.meteoSI.rh2q(50.0, 280.0, 90000.0)
    except TypeError:
        pass
    else:
        raise AssertionError("expected TypeError for rh > 5")
