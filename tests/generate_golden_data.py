"""Regenerate tests/golden/rain_snow_reference.npz.

Only run this deliberately, after confirming a change to the RT/scattering
code is an intentional physics change and not a regression:

    python tests/generate_golden_data.py
"""

import os
from pathlib import Path

import numpy as np

os.environ.setdefault("PAMTRA_DATADIR", "")

from _scenario import run_scenario  # noqa: E402

GOLDEN_DIR = Path(__file__).parent / "golden"


def main():
    pam = run_scenario()
    GOLDEN_DIR.mkdir(exist_ok=True)
    np.savez(
        GOLDEN_DIR / "rain_snow_reference.npz",
        Ze=pam.r["Ze"],
        Att_hydro=pam.r["Att_hydro"],
        tb=pam.r["tb"],
    )
    print(f"Wrote {GOLDEN_DIR / 'rain_snow_reference.npz'}")


if __name__ == "__main__":
    main()
