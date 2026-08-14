..  _profiles:


Building a Profile
===================

The profile is PAMTRA's description of the atmospheric column(s) to
simulate: temperature, pressure, humidity, hydrometeor content, and
surface properties on a (x, y, z) grid. It's stored in ``pam.p`` (a dict
of numpy arrays) and is separate from the :ref:`descriptorFile`, which
describes *how* each hydrometeor scatters rather than *how much* of it is
present at each grid point.

There are two ways to build one: :any:`pyPamtra.core.pyPamtra.createProfile`
for data you already have in memory, or one of the importers in
:any:`pyPamtra.importer` for reading a specific model's native output
format directly.


createProfile
*************

:any:`pyPamtra.core.pyPamtra.createProfile` takes keyword arguments named
after the profile variables. Variables ending in ``_lev`` are given at
layer *boundaries* (one more entry than the number of layers); variables
without the suffix are layer-average values. Everything is in SI units,
except ``relhum``, which is in percent.

Mandatory: ``hgt_lev``, (``temp_lev`` or ``temp``), (``press_lev`` or
``press``), and (``relhum_lev`` or ``relhum``).

Optional, guessed with a warning if not provided: ``timestamp``, ``lat``,
``lon``, ``wind10u``, ``wind10v``, ``hydro_q``, ``hydro_n``,
``hydro_reff``, ``obs_height``, ``sfc_type``, ``sfc_model``, ``sfc_refl``,
``sfc_salinity``.

If hydrometeors are present, ``hydro_q`` (mass mixing ratio, kg/kg) and/or
``hydro_n`` / ``hydro_reff`` need a last dimension matching the number of
hydrometeors already added to ``pam.df`` (see :ref:`descriptorFile`) --
add hydrometeors to the descriptor file *before* calling ``createProfile``
with hydrometeor data. Which of ``hydro_q``/``hydro_n``/``hydro_reff`` is
actually required for a given hydrometeor depends on that hydrometeor's
``moment_in`` field in the descriptor file.

::

    import numpy as np
    import pyPamtra

    pam = pyPamtra.pyPamtra()
    pam.df.readFile("descriptorfiles/descriptor_file_2m_icon.txt")

    pam.createProfile(
        hgt_lev=np.array([[[0., 500., 1000.]]]),
        temp_lev=np.array([[[290., 285., 280.]]]),
        press_lev=np.array([[[101300., 95000., 89000.]]]),
        relhum_lev=np.array([[[80., 90., 95.]]]),
        hydro_q=np.zeros((1, 1, 2, 6)),  # (x, y, layers, n_hydrometeors)
    )

:any:`pyPamtra.core.pyPamtra.createFullProfile` is a variant that additionally
lets you pass all the ``sfc_*``/wind/etc. arguments positionally instead
of relying on the create-if-missing defaults.


Importing model output
***********************

:any:`pyPamtra.importer` has 17 functions, each reading one model's
native output format and calling ``createProfile`` for you. They all take
a ``descriptorFile`` argument (a path or an already-loaded
:any:`pyPamtra.descriptorFile.pamDescriptorFile`) so hydrometeor moments
in the model output can be matched to the right descriptor rows.

.. list-table::
   :header-rows: 1

   * - Function
     - Model / format
   * - ``readWrfDataset``
     - WRF
   * - ``readCosmoDe1MomDataset``
     - COSMO-DE, 1-moment microphysics
   * - ``readCosmoDe2MomDataset``
     - COSMO-DE, 2-moment microphysics
   * - ``readCosmoDe2MomDatasetOnFlightTrack``
     - COSMO-DE 2-moment, sampled along a flight track
   * - ``readCosmoReAn2km`` / ``readCosmoReAn6km``
     - COSMO reanalysis, 2 km / 6 km grid
   * - ``readIconNwp1MomDataset``
     - ICON-NWP, 1-moment microphysics
   * - ``readIconNwp2MomDataset``
     - ICON-NWP, 2-moment microphysics
   * - ``readIconNwp1MomDataset_cells``
     - ICON-NWP 1-moment, unstructured cell grid
   * - ``readIcon1momMeteogram`` / ``readIcon2momMeteogram``
     - ICON meteogram output, 1-mom / 2-mom
   * - ``readICON2mom`` / ``readICON1momNWP``
     - ICON 2-moment / 1-moment NWP output
   * - ``readIcon2momOnFlightTrack``
     - ICON 2-moment, sampled along a flight track
   * - ``readHIRHAM``
     - HIRHAM regional climate model
   * - ``readECMWF``
     - ECMWF/IFS
   * - ``readMesoNH``
     - Meso-NH
   * - ``createUsStandardProfile``
     - Synthetic US Standard Atmosphere (no external data needed)
   * - ``ncToDict``
     - Generic NetCDF-to-dict reader, not model-specific

Each function's own docstring (see the :ref:`pyPamtra` API reference, or
read ``python/pyPamtra/importer.py`` directly) documents its specific
arguments -- they differ per model (grid subsetting, forecast index,
tarball handling, etc.).


Parameterized PSD vs. hydro_fullspec
**************************************

By default, each hydrometeor's particle size distribution is computed
from the analytic form named in its descriptor file row (``dist_name``,
``p_1``..``p_4``, see :ref:`descriptorFile`) plus whatever moment(s) are
in the profile (``hydro_q``/``hydro_n``/``hydro_reff``).

Setting ``pam.nmlSet["hydro_fullspec"] = True`` switches to passing
explicit, discrete size-bin arrays (diameter, mass, number concentration,
density, aspect ratio per bin) directly from Python instead, bypassing
the analytic PSD entirely. Build these with
:any:`pyPamtra.descriptorFile.pamDescriptorFile.addFullSpectra` after
adding the corresponding hydrometeors to ``pam.df`` -- useful for
bin-resolved microphysics schemes or in situ particle spectra that don't
fit a closed-form PSD. See ``examples/pyPamTestFullSpec.py``.
