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
