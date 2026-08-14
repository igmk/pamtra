"""End-to-end test of runPicklePamtra: it writes .job pickles to disk and
waits for a separate consumer process to write back .result pickles (see
tools/picklePam/ for the reference consumer, which is unrelated legacy
Python 2 code and out of scope here). This test acts as that consumer in
a background thread, so the actual runPicklePamtra code path -- job
writing, result waiting/reading, and _joinResults -- is exercised for
real rather than mocked.
"""

import glob
import os
import pickle
import threading
import time

import numpy as np

from _scenario import FREQUENCIES, build_pamtra
from pyPamtra.libWrapper import parallelPamtraFortranWrapper


def _consumer(pickle_path, stop_event, poll_interval=0.05):
    while not stop_event.is_set():
        for job_file in glob.glob(os.path.join(pickle_path, "*.job")):
            claimed = job_file + ".taken"
            try:
                os.rename(job_file, claimed)
            except OSError:
                continue  # another worker got there first

            with open(claimed, "rb") as f:
                indices, settings, nmlSet, dfPart, dfPart4D, dfPartFS, profilePart = pickle.load(f)

            result = parallelPamtraFortranWrapper(
                indices, settings, nmlSet, dfPart, dfPart4D, dfPartFS, profilePart,
                returnModule=False,
            )

            result_file = job_file[: -len(".job")] + ".result"
            with open(result_file + ".tmp", "wb") as f:
                pickle.dump(result, f)
            os.rename(result_file + ".tmp", result_file)
            os.remove(claimed)
        time.sleep(poll_interval)


def test_run_pickle_pamtra_matches_run_pamtra(tmp_path):
    pam = build_pamtra()

    stop_event = threading.Event()
    consumer = threading.Thread(target=_consumer, args=(str(tmp_path), stop_event), daemon=True)
    consumer.start()
    try:
        pam.runPicklePamtra(FREQUENCIES, picklePath=str(tmp_path), maxWait=30)
    finally:
        stop_event.set()
        consumer.join(timeout=5)

    reference = build_pamtra()
    reference.runPamtra(FREQUENCIES)

    np.testing.assert_allclose(pam.r["Ze"], reference.r["Ze"], rtol=1e-5)
    np.testing.assert_allclose(pam.r["tb"], reference.r["tb"], rtol=1e-5)
