..  _quickstart:


Quickstart
==========

This is the shortest path from an installed :ref:`pyPamtra` to a first
brightness temperature and radar reflectivity. It's meant as a starting
point for both new users and AI assistants working in this repository --
see :ref:`profiles`, :ref:`running`, and :ref:`results` for the details
behind each step.

A PAMTRA run always needs three things: a **descriptor file** (what
hydrometeors are present and how their scattering is computed), a
**profile** (the atmospheric state), and a call to one of the ``run*``
methods.

::

    import pyPamtra

    pam = pyPamtra.pyPamtra()

    # One liquid hydrometeor ("cwc_q"), exponential PSD, Mie-sphere
    # scattering. See :ref:`descriptorFile` for what each field means.
    pam.df.addHydrometeor((
        "cwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
        3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
        "mie-sphere", "khvorostyanov01_drops", -99.0,
    ))

    # A built-in US Standard Atmosphere profile needs no external data.
    pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0.0, 1000.0, 2000.0])

    # Put some liquid water content into the lowest layer.
    pam.p["hydro_q"][0, 0, 0, 0] = 1e-3

    pam.runPamtra([35.5])  # frequency in GHz

    print(pam.r["tb"][0, 0, 0, 0])  # brightness temperature, [angle, V/H]
    print(pam.r["Ze"][0, 0, 0, 0])  # radar reflectivity, [layer]

This runs entirely without ``PAMTRA_DATADIR`` pointing at real data, because
Mie-sphere scattering and :any:`pyPamtra.importer.createUsStandardProfile`
need no external database -- set ``export PAMTRA_DATADIR=""`` first if you
want to skip PAMTRA's one-time automatic data download for this example
(see :ref:`installation`). Real work usually starts from an existing descriptor
file in ``descriptorfiles/`` (:ref:`descriptorFile`) and a profile built
from your own model output via one of the importers in
:any:`pyPamtra.importer` (:ref:`profiles`), rather than
``createUsStandardProfile`` and a hand-built hydrometeor.

See ``examples/`` in the repository for complete, runnable scripts and
notebooks, e.g. ``examples/simple_passive.ipynb`` and
``examples/radar_example.py``.
