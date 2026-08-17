..  _pyPamtraXr:


pyPamtraXr
==========

:any:`pyPamtra.pamtra_xr.pyPamtraXr` is a small, curated facade over :ref:`pyPamtra`
for users who'd rather not touch ``pam.p["key"][x, y, z]``-style
dict/array indexing directly. It is recommended for new code -- see
:ref:`quickstart` for a full first example.

There is no separate, "live" xarray-native copy of the profile or
results behind ``pyPamtraXr``: every array still lives in
``pamxr.pam.p``/``pamxr.pam.r`` exactly as in plain :ref:`pyPamtra`,
unchanged. ``pyPamtraXr``'s methods are verbs that delegate to it, and
the only way to get a labeled ``xarray.Dataset`` out is to explicitly
ask for one, via ``to_xarray()`` or the return value of ``run()`` --
never a stored, directly-settable ``.profile``/``.results`` attribute.
This keeps the actual state in exactly one place and matches the
copy/snapshot semantics :any:`pyPamtra.core.pyPamtra.to_xarray` /
:any:`pyPamtra.core.pyPamtra.from_xarray` already use.

Escape hatch
****************

``pamxr.pam`` is the underlying :ref:`pyPamtra` object -- always the
same one, not a copy. Anything not covered by ``pyPamtraXr``'s curated
methods below (``nmlSet``/``set`` tweaks, per-grid-point array edits, or
any other :ref:`pyPamtra` method) is reachable through it, and any
change made through it is immediately visible to ``pyPamtraXr`` too,
since it's the same object::

    pamxr = pyPamtra.pyPamtraXr()
    pamxr.pam.nmlSet["hydro_threshold"] = 1e-8   # not covered by pyPamtraXr's own methods
    pamxr.pam.p["hydro_q"][0, 0, 0, 0] = 1e-3    # per-grid-point edit

You can also wrap an already-built :ref:`pyPamtra` object instead of
creating a fresh one, e.g. to migrate an existing script incrementally::

    pam = pyPamtra.pyPamtra()
    # ... existing script builds pam as usual ...
    pamxr = pyPamtra.pyPamtraXr(pam=pam)

Class reference
********************

.. automodule:: pyPamtra.pamtra_xr
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
