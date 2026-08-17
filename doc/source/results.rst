..  _results:


Reading Results
=================

After a ``run*`` call (:ref:`running`), results live in ``pam.r``, a dict
of numpy arrays keyed by output variable. The leading dimensions are
always the profile grid (``ngridx``, ``ngridy``); the remaining
dimensions depend on the variable -- e.g. output level and angle for
passive output, layer for active output, frequency and polarization
where applicable. Check ``pam.r[var].shape`` and the corresponding
:ref:`settings` (``nfreqs``, ``noutlevels``, number of hydrometeors in
:ref:`descriptorFile`, ``radar_npol``) to make sense of a given array
rather than assuming a fixed axis order across variables.

The most commonly used variables:

.. list-table::
   :header-rows: 1

   * - Variable
     - Contents
   * - ``tb``
     - Brightness temperature [K], passive output
   * - ``Ze``
     - Radar reflectivity [dBZ], active output
   * - ``Att_hydro``
     - Two-way hydrometeor attenuation [dB]
   * - ``Att_atmo``
     - Two-way atmospheric gas attenuation [dB]
   * - ``radar_spectra``
     - Full Doppler spectrum (only if ``nmlSet["radar_mode"]`` is
       ``"spectrum"`` or ``"moments"``)
   * - ``radar_vel``
     - Doppler velocity bins for ``radar_spectra``
   * - ``radar_moments``
     - Radar moments (mean Doppler velocity, skewness, kurtosis, ...)
       derived from the spectrum
   * - ``radar_slopes``
     - Left/right spectral slopes
   * - ``emissivity``
     - Surface emissivity actually used

``pam.r["nmlSettings"]`` carries a copy of the ``nmlSet`` the run was made
with, so a saved result is self-describing.


Writing output
****************

:any:`pyPamtra.core.pyPamtra.writeResultsToNetCDF` writes ``pam.r`` (plus
selected profile variables) to a NetCDF file::

    pam.writeResultsToNetCDF("output.nc", ncForm="NETCDF4")

:any:`pyPamtra.core.pyPamtra.writeResultsToNumpy` /
:any:`pyPamtra.core.pyPamtra.loadResultsFromNumpy` round-trip the *entire*
session (``pam.r``, ``pam.p``, ``pam.nmlSet``, ``pam.set``, and the
descriptor file) through a single pickle file instead -- lighter-weight
than NetCDF and useful for resuming/inspecting a run later, but Python
pickles rather than NetCDF, so they aren't meant for sharing outside
Python or across pyPamtra versions. Pass ``seperateFiles=True`` to write
one plain ``.npy`` file per variable into a directory instead, useful for
very large result sets.


Quick-look plotting
**********************

:any:`pyPamtra.plot` has two plotting helpers built on ``pam.r``:
:any:`pyPamtra.plot.plotTB` (a 2D brightness-temperature map over the
profile grid) and :any:`pyPamtra.plot.plotTBLine` (brightness temperature
along one grid axis, e.g. for a flight track or transect). Both operate
directly on ``pam.r`` after a run.


xarray interface
**********************

:any:`pyPamtra.core.pyPamtra.to_xarray` builds an ``xarray.Dataset``
snapshot of ``pam.p`` or ``pam.r``, with dimension names and units taken
from PAMTRA's own metadata instead of the bare, unlabeled arrays in the
dicts::

    ds = pam.to_xarray(source="r")   # or source="p" for the profile
    ds["tb"]                         # labeled dims, units in .attrs
    ds.to_netcdf("output.nc")        # xarray's own writer

This is an additive, read-only snapshot: it never holds a live view into
``pam.p``/``pam.r``, so mutating the returned ``Dataset`` never affects
the ``pyPamtra`` object, and nothing about ``pam.p``/``pam.r`` themselves
changes. Requires the optional ``xarray`` package
(``pip install pamtra[xarray]``).

Pass ``outer_dims`` to relabel the leading grid dimensions to something
more meaningful than the generic grid index, e.g.
``pam.to_xarray(outer_dims={"grid_x": "lat"})``.

:any:`pyPamtra.core.pyPamtra.from_xarray` is the inverse for profiles:
it builds a profile from an ``xarray.Dataset`` (e.g. one assembled
independently, not necessarily round-tripped through ``to_xarray``) by
feeding its variables into :any:`pyPamtra.core.pyPamtra.createProfile`,
so all of that method's usual defaulting and validation still applies::

    pam2 = pyPamtra.pyPamtra()
    pam2.df.addHydrometeor(...)   # descriptor file is not part of the Dataset
    pam2.from_xarray(ds)
