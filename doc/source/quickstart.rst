..  _quickstart:


Quickstart
==========

This is the shortest path from an installed PAMTRA to a first brightness
temperature and radar reflectivity. It's meant as a starting point for
both new users and AI assistants working in this repository -- see
:ref:`profiles`, :ref:`running`, and :ref:`results` for the details
behind each step.

A PAMTRA run always needs three things: a **descriptor file** (what
hydrometeors are present and how their scattering is computed), a
**profile** (the atmospheric state), and a call to one of the ``run*``
methods.

There are two ways to drive this: :ref:`pyPamtraXr`, a small curated
facade recommended for new code, and :ref:`pyPamtra` itself, the
lower-level dict-of-arrays interface every part of PAMTRA is ultimately
built on (needed for anything :ref:`pyPamtraXr` doesn't cover, and for
extending an existing script that already uses it).

Recommended: pyPamtraXr
****************************

::

    import pyPamtra

    pamxr = pyPamtra.pyPamtraXr()

    # One liquid hydrometeor ("cwc_q"), exponential PSD, Mie-sphere
    # scattering. See :ref:`descriptorFile` for what each field means.
    pamxr.add_hydrometeor(
        hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
        dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
        scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
    )

    # A built-in US Standard Atmosphere profile needs no external data.
    # pamxr.p is now a real, labeled xarray.Dataset (see :ref:`pyPamtraXr`).
    pamxr.set_profile(hgt_lev=[0.0, 1000.0, 2000.0], temp_lev=[288.15, 281.5, 275.0],
                       press_lev=[101300.0, 89000.0, 79000.0], relhum_lev=[80.0, 70.0, 60.0])

    # Put some liquid water content into the lowest layer.
    pamxr.p["hydro_q"].values[0, 0, 0, 0] = 1e-3

    results = pamxr.run(35.5)  # frequency in GHz

    print(results["tb"].values)  # brightness temperature, [gridx, gridy, outlevel, angle, frequency, V/H]
    print(results["Ze"].values)  # radar reflectivity, [gridx, gridy, layer, frequency, radar_npol, radar_npeaks]

``results`` is a labeled ``xarray.Dataset`` (units in ``.attrs``, see
:ref:`results`'s xarray section) -- ``pamxr.run()`` returns it directly,
and it's the same object as ``pamxr.r`` afterwards. Write it out with
``pamxr.to_netcdf("output.nc")``.

Escape hatch
    ``pamxr.pam`` is the underlying :ref:`pyPamtra` object -- the same
    one, always, not a copy -- for anything not covered by
    :ref:`pyPamtraXr`'s curated attributes/methods (``nmlSet``/``set``
    tweaks, or any other :ref:`pyPamtra` method). It's also how existing
    :any:`pyPamtra.importer` functions plug in, via
    ``pyPamtra.pyPamtraXr.from_pam(pam)`` -- see :ref:`pyPamtraXr`.

Multiple instrument configurations
    :ref:`pyPamtraXr`'s ``add_instrument()`` runs several
    frequency/``nmlSet`` configurations against the same profile, each
    keeping its own results -- see :ref:`running`'s instrument section.

Lower-level: pyPamtra
****************************

The same run, built directly on :ref:`pyPamtra`'s dict-of-arrays
interface -- what ``pyPamtraXr`` is a facade over, and what's needed for
anything outside its curated surface anyway::

    import pyPamtra

    pam = pyPamtra.pyPamtra()

    # (addHydrometeor also accepts a positional tuple in field order --
    # see :ref:`descriptorFile` -- but the keyword form below is safer:
    # unrecognized or missing fields raise instead of silently landing
    # in the wrong column.)
    pam.df.addHydrometeor(
        hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
        dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
        scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
    )

    pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0.0, 1000.0, 2000.0])
    pam.p["hydro_q"][0, 0, 0, 0] = 1e-3

    pam.runPamtra([35.5])  # frequency in GHz

    print(pam.r["tb"][0, 0, 0, 0])  # brightness temperature, [angle, V/H]
    print(pam.r["Ze"][0, 0, 0, 0])  # radar reflectivity, [layer]

Both examples run entirely without ``PAMTRA_DATADIR`` pointing at real
data, because Mie-sphere scattering and
:any:`pyPamtra.importer.createUsStandardProfile` need no external
database -- set ``export PAMTRA_DATADIR=""`` first if you want to skip
PAMTRA's one-time automatic data download for this example (see
:ref:`installation`). Real work usually starts from an existing
descriptor file in ``descriptorfiles/`` (:ref:`descriptorFile`) and a
profile built from your own model output via one of the importers in
:any:`pyPamtra.importer` (:ref:`profiles`), rather than
``createUsStandardProfile`` and a hand-built hydrometeor.

See ``examples/`` in the repository for complete, runnable scripts and
notebooks, e.g. ``examples/pyPamtraXr_quickstart.py``,
``examples/simple_passive.ipynb``, and ``examples/radar_example.py``.
