..  _pyPamtraXr:


pyPamtraXr
==========

:any:`pyPamtra.pamtra_xr.pyPamtraXr` is a small, curated facade over
:ref:`pyPamtra` for users who'd rather not touch
``pam.p["key"][x, y, z]``-style dict/array indexing or ``pam.df``'s
positional-tuple recarray rows directly. It is recommended for new code
-- see :ref:`quickstart` for a full first example.

Unlike a plain :ref:`pyPamtra` object, ``.p``/``.r``/``.df``/``.df_4d``/
``.df_full_spec`` are real, persistent, labeled attributes here:

.. list-table::
   :header-rows: 1

   * - Attribute
     - Type
     - What it replaces
   * - ``pamxr.p``
     - ``xarray.Dataset``
     - ``pam.p`` (the atmospheric profile)
   * - ``pamxr.r``
     - ``xarray.Dataset``
     - ``pam.r`` (results, empty before the first ``run()``)
   * - ``pamxr.df``
     - ``pandas.DataFrame``, indexed by ``hydro_name``
     - ``pam.df.data`` (per-hydrometeor scalar properties)
   * - ``pamxr.df_4d``
     - ``xarray.Dataset``
     - ``pam.df.data4D`` (per-gridpoint property overrides)
   * - ``pamxr.df_full_spec``
     - ``xarray.Dataset``
     - ``pam.df.dataFullSpec`` (full particle-size-spectrum input)

``.df`` is a ``pandas.DataFrame``, not an ``xarray.Dataset`` --
it's a small, mixed-dtype (strings, ints, floats), few-row table (one
row per hydrometeor), which is what pandas is for; xarray is built for
labeled N-D arrays like ``.p``/``.r``/``.df_4d``/``.df_full_spec``, and
is a comparatively awkward fit for a handful of heterogeneous scalar
columns. This costs nothing extra: pandas is already a hard dependency
of xarray.

There is still exactly one place the RT engine actually runs:
``pamxr.pam``, a real :ref:`pyPamtra` object, always accessible as the
escape hatch for anything this class doesn't cover. ``run()`` (and
``add_instrument()``) translate ``.p`` into ``pamxr.pam.p`` right before
calling ``pamxr.pam.runPamtra()``/``addInstrument()``, via the existing
:any:`pyPamtra.core.pyPamtra.from_xarray` -- not a new, independently
maintained reimplementation of ``createProfile()``'s shape/defaulting
logic. ``.df``/``.df_4d``/``.df_full_spec`` instead stay eagerly synced
from ``pamxr.pam.df`` on every ``add_hydrometeor()``/
``remove_hydrometeor()``/``add_4d()``/``remove_4d()``/
``add_full_spectra()`` call, by calling the existing, already-tested
:any:`pyPamtra.descriptorFile.pamDescriptorFile` methods and re-deriving
the pandas/xarray mirror from the result -- so there is exactly one
implementation of what a valid hydrometeor row looks like, not two.

Escape hatch
****************

``pamxr.pam`` is the underlying :ref:`pyPamtra` object -- always the
same one, not a copy. Anything not covered by ``pyPamtraXr``'s curated
methods (``nmlSet``/``set`` tweaks, or any other :ref:`pyPamtra` method)
is reachable through it::

    pamxr = pyPamtra.pyPamtraXr()
    pamxr.pam.nmlSet["hydro_threshold"] = 1e-8   # not covered by pyPamtraXr's own methods

If you mutate ``pamxr.pam`` directly, ``.p``/``.r``/``.df``/``.df_4d``/
``.df_full_spec`` go stale until the next ``run()`` (which only
refreshes ``.p``/``.r``) or an explicit ``pamxr.refresh()`` call --
e.g. after set_profile() has already populated ``pamxr.pam.p["hydro_q"]``::

    pamxr.pam.p["hydro_q"][0, 0, 0, 0] = 1e-3   # bypasses pamxr.p
    pamxr.refresh()                             # pamxr.p now reflects it

You can also wrap an already-built :ref:`pyPamtra` object instead of
creating a fresh one -- most usefully, to reuse any of the existing
:any:`pyPamtra.importer` functions (which all build and return a plain
:ref:`pyPamtra` object) without modifying them::

    pam = pyPamtra.importer.createUsStandardProfile(pyPamtra.pyPamtra(), hgt_lev=[0.0, 1000.0])
    pamxr = pyPamtra.pyPamtraXr.from_pam(pam)   # pam becomes pamxr.pam, not a copy

Naming (and counting) the grid dimensions
***********************************************

The leading dimensions of every array are a generic 2D grid index
(``grid_x``/``grid_y``) by default, since that's all a plain
:ref:`pyPamtra` profile natively is. ``outer_dims`` is an ordered list
naming what your data actually varies along instead -- a latitude and
longitude, a time series, a flight track -- and, unlike a plain
:ref:`pyPamtra` object, it doesn't have to be exactly 2 names::

    pamxr = pyPamtra.pyPamtraXr(outer_dims=["lat", "lon"])   # 2: a direct rename
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["scan"])          # 1
    pamxr = pyPamtra.pyPamtraXr(outer_dims=["time", "lat", "lon"])  # 3+

:ref:`pyPamtra`'s own core is genuinely, always a 2D grid (its Fortran
layer is set up with exactly ``atmo_ngridx``/``atmo_ngridy``, nothing
else) -- there's no way around that, and ``pyPamtraXr`` doesn't try to.
2 names is therefore always a direct, lossless rename. Any other count
is handled by folding the named axes into ``pyPamtra``'s single grid
axis (with the second one left at size 1) right before handing off to
``self.pam``, and splitting it back apart when reading ``.p``/``.r``/
instrument results back out -- entirely within ``pyPamtraXr``, never
touching :any:`pyPamtra.core.pyPamtra.to_xarray`/``from_xarray`` (which
stay genuinely, natively 2D) or the Fortran boundary.

This applies to every ``xarray.Dataset`` the object builds afterwards --
``.p``, ``.r``, and any :any:`pyPamtra.instrument.PamtraInstrument`
results -- not just one call's output. It's stored as ``pamxr.outer_dims``;
change it directly and call ``refresh()`` to re-label the current state::

    pamxr.outer_dims = ["time", "grid_y"]
    pamxr.refresh()

Every method that takes its own ``outer_dims`` argument (``run()``,
``run_parallel()``, ``add_instrument()``, ``set_profile_from_xarray()``)
defaults to ``pamxr.outer_dims`` and accepts a per-call override
instead.

With 3 or more names, ``refresh()`` (and anything built on it) needs to
know the original sizes to split the flattened grid back apart -- that
information only exists once a profile has actually been set through
``set_profile()``/``set_profile_from_xarray()``. Calling ``refresh()``
before that (e.g. after building a profile directly on ``pamxr.pam``,
bypassing ``pyPamtraXr`` entirely) raises ``RuntimeError`` rather than
guessing.

Class reference
********************

.. automodule:: pyPamtra.pamtra_xr
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
