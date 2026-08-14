"""Download and cache PAMTRA's external data (scattering databases and
surface emissivity atlases) and point PAMTRA_DATADIR at it.

Deliberately kept outside the ``pyPamtra`` package: ``import pyPamtra``
raises at import time if ``PAMTRA_DATADIR`` isn't set (see
python/pyPamtra/core.py), so a helper meant to be run *before* that env
var exists can't live inside the package it's setting up for.

This data is optional -- pyPamtra works without it (Mie-sphere
scattering and the built-in surface emissivity defaults need no external
files at all, see doc/source/installation.rst). It's needed for the
Liu/Hong DDA scattering databases and the FASTEM/TESSEM2/TELSEM2 surface
emissivity atlases.

Usage from Python, before importing pyPamtra::

    import pamtra_data
    pamtra_data.fetch_data()
    import pyPamtra

Usage from the shell, via the console script this module registers::

    export PAMTRA_DATADIR=$(pamtra-fetch-data)
"""

import argparse
import os
import sys

import pooch

DATA_URL = "https://github.com/igmk/pamtra/releases/download/data-v1/pamtra_data.tar.bz2"
DATA_HASH = "sha256:ae647647e634a8d26c8ac6906950c64a4dc6d203493fef4a759a15039dd6d018"
ARCHIVE_NAME = "pamtra_data.tar.bz2"
EXTRACT_DIR = "pamtra_data"


def _pooch():
    return pooch.create(
        path=pooch.os_cache("pamtra"),
        base_url="",
        registry={ARCHIVE_NAME: DATA_HASH},
        urls={ARCHIVE_NAME: DATA_URL},
    )


def fetch_data(progressbar=True):
    """
    Download (once) and cache PAMTRA's scattering-database and surface
    emissivity data, and set the ``PAMTRA_DATADIR`` environment variable
    to point at it for the current process.

    On repeat calls this returns the cached path immediately without
    re-downloading, unless the cached copy is missing or its checksum no
    longer matches.

    Parameters
    ----------
    progressbar : bool, optional
        Show a download progress bar (default True).

    Returns
    -------
    str
        Path to the extracted data directory.
    """
    p = _pooch()
    p.fetch(
        ARCHIVE_NAME,
        processor=pooch.Untar(extract_dir=EXTRACT_DIR),
        progressbar=progressbar,
    )
    # pooch.Untar(extract_dir=...) always unpacks into <pooch cache>/<extract_dir>,
    # regardless of whether the archive itself has a single wrapping directory
    # (this one doesn't -- it unpacks straight into emissivity/, hongdb/, scatdb/).
    datadir = os.path.join(p.path, EXTRACT_DIR)
    os.environ["PAMTRA_DATADIR"] = datadir
    return datadir


def _cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Download and cache PAMTRA's scattering/emissivity data, "
                     "print the resulting PAMTRA_DATADIR path to stdout."
    )
    parser.add_argument(
        "--quiet", action="store_true", help="suppress the download progress bar"
    )
    args = parser.parse_args(argv)

    datadir = fetch_data(progressbar=not args.quiet)
    print(datadir)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
