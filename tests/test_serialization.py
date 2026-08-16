import numpy as np

from _scenario import run_scenario


def test_write_load_results_numpy_roundtrip(tmp_path):
    # writeResultsToNumpy/loadResultsFromNumpy pickle the session; both
    # used to open their file in text mode, which pickle rejects in
    # Python 3 (TypeError: write() argument must be str, not bytes).
    pam = run_scenario()

    pkl_file = tmp_path / "session.pkl"
    pam.writeResultsToNumpy(str(pkl_file))

    pam2 = type(pam)()
    pam2.loadResultsFromNumpy(str(pkl_file))

    np.testing.assert_array_equal(pam2.r["tb"], pam.r["tb"])
    np.testing.assert_array_equal(pam2.r["Ze"], pam.r["Ze"])
