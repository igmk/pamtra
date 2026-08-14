"""Shared golden-output regression scenario.

Used by both test_regression.py and generate_golden_data.py, so the two
can never drift apart. Deliberately uses only Mie-sphere scattering and
pyPamtra.importer.createUsStandardProfile, so it needs no external
scattering database or emissivity atlas and works with PAMTRA_DATADIR="".
"""

import pyPamtra

FREQUENCIES = [10.65, 35.5]


def build_pamtra():
    pam = pyPamtra.pyPamtra()

    # Rain-like liquid hydrometeor: exponential PSD, mass concentration only.
    pam.df.addHydrometeor((
        "rwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
        3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
        "mie-sphere", "khvorostyanov01_drops", -99.0,
    ))
    # Snow-like ice hydrometeor: exp_field_t PSD, mass concentration only.
    pam.df.addHydrometeor((
        "swc_q", -99.0, -1, -99.0, 0.038, 2.0, 0.3971, 1.88,
        3, 30, "exp_field_t", -99.0, -99.0, -99.0, -99.0, 0.51e-10, 2.0e-2,
        "mie-sphere", "heymsfield10_particles", -99.0,
    ))

    pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0.0, 1000.0, 2000.0, 3000.0])
    pam.set["verbose"] = 0
    pam.set["pyVerbose"] = 0
    pam.nmlSet["randomseed"] = 10

    pam.p["hydro_q"][0, 0, 1, 0] = 1e-3  # rwc_q, layer index 1
    pam.p["hydro_q"][0, 0, 2, 1] = 5e-4  # swc_q, layer index 2

    return pam


def run_scenario():
    pam = build_pamtra()
    pam.runPamtra(FREQUENCIES)
    return pam
