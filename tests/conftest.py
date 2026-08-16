import os

# pyPamtra raises RuntimeError at import time if PAMTRA_DATADIR is unset (see
# python/pyPamtra/core.py). An empty string is an explicitly supported value
# for the parts of PAMTRA exercised here (Mie-sphere scattering, no surface
# emissivity atlas / external scattering database lookups).
os.environ.setdefault("PAMTRA_DATADIR", "")
