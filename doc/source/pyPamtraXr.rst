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

Class reference
********************

.. automodule:: pyPamtra.pamtra_xr
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
