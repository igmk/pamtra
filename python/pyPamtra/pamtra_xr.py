# -*- coding: utf-8 -*-
'''
pyPamtraXr: a curated, xarray-native facade over pyPamtra, for users who
would rather never touch pam.p["key"][x,y,z]-style dict/array indexing.

Deliberately does not expose the profile or results as a public,
directly-settable attribute -- there is no "live" xarray Dataset backing
this class. Every raw array still lives in self.pam.p/self.pam.r exactly
as in plain pyPamtra (see pyPamtra.core), unchanged; pyPamtraXr's methods
are verbs (add_hydrometeor, set_profile, run, ...) that delegate to it,
and the only way to get a Dataset out is to explicitly ask for one via
to_xarray() or the return value of run(). This keeps the actual state in
one place (self.pam), avoids a second xarray-native copy of the data
getting out of sync with it, and matches the copy/snapshot semantics
pyPamtra.to_xarray()/from_xarray() already use.

self.pam is always the same underlying pyPamtra object and is fully
accessible for anything not covered by this curated surface (e.g.
nmlSet/set tweaks, or legacy methods like writeResultsToNumpy()).
'''

from .core import pyPamtra


class pyPamtraXr(object):
  '''
  Curated, xarray-native facade over a pyPamtra object.

  Parameters
  ----------
  pam : pyPamtra, optional
      Wrap an existing pyPamtra object instead of creating a fresh one --
      useful to migrate an existing script incrementally.

  Attributes
  ----------
  pam : pyPamtra
      The wrapped pyPamtra object. Always fully accessible: anything not
      covered by this class's curated methods (nmlSet/set tweaks, or any
      other pyPamtra method) is reachable via pamxr.pam.*.

  Examples
  --------
  >>> pamxr = pyPamtra.pyPamtraXr()
  >>> pamxr.add_hydrometeor(
  ...     hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
  ...     dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
  ...     scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
  ... )
  >>> pamxr.set_profile(hgt_lev=[0., 1000., 2000.], temp_lev=[288., 281., 275.],
  ...                    press_lev=[101300., 90000., 80000.], relhum_lev=[80., 70., 60.])
  >>> results = pamxr.run(35.5)
  >>> results["tb"]
  '''

  def __init__(self,pam=None):
    self.pam = pam if pam is not None else pyPamtra()

  # -- hydrometeors --------------------------------------------------

  def add_hydrometeor(self,**kwargs):
    '''
    Add a hydrometeor by keyword arguments. See
    pyPamtra.descriptorFile.pamDescriptorFile.addHydrometeor (the
    keyword form) for the field reference.
    '''
    return self.pam.df.addHydrometeor(**kwargs)

  def remove_hydrometeor(self,name):
    '''Remove a hydrometeor previously added via add_hydrometeor().'''
    return self.pam.df.removeHydrometeor(name)

  # -- profile ---------------------------------------------------------

  def set_profile(self,**kwargs):
    '''
    Define the atmospheric profile. See
    pyPamtra.core.pyPamtra.createProfile for the accepted keyword
    arguments (hgt_lev/temp_lev/press_lev/relhum_lev are mandatory,
    everything else is optional and defaulted with a warning).
    '''
    return self.pam.createProfile(**kwargs)

  def set_profile_from_xarray(self,ds,outer_dims=None):
    '''
    Define the atmospheric profile from an xarray.Dataset (e.g. one
    assembled from a model's own xarray-native output). See
    pyPamtra.core.pyPamtra.from_xarray. Hydrometeors must already be
    added via add_hydrometeor() first if the Dataset includes hydro_q/
    hydro_n/hydro_reff.
    '''
    return self.pam.from_xarray(ds,outer_dims=outer_dims)

  # -- running -----------------------------------------------------------

  def run(self,freqs,outer_dims=None,**kwargs):
    '''
    Run the RT engine (pyPamtra.core.pyPamtra.runPamtra) and return the
    results as an xarray.Dataset (also reachable afterwards via
    to_xarray(source="r")).

    Parameters
    ----------
    freqs : float or list of float
        Frequencies in GHz.
    outer_dims : dict, optional
        Passed through to to_xarray().
    **kwargs
        Passed through to runPamtra() (e.g. checkData=False).

    Returns
    -------
    xarray.Dataset
    '''
    self.pam.runPamtra(freqs,**kwargs)
    return self.to_xarray(source="r",outer_dims=outer_dims)

  def run_parallel(self,freqs,outer_dims=None,**kwargs):
    '''
    Run the RT engine across local CPU cores
    (pyPamtra.core.pyPamtra.runParallelPamtra) and return the results as
    an xarray.Dataset. See run() for the parameters shared with the
    serial variant; **kwargs also accepts runParallelPamtra's
    pp_local_workers/pp_deltaF/pp_deltaX/pp_deltaY/timeout.

    Returns
    -------
    xarray.Dataset
    '''
    self.pam.runParallelPamtra(freqs,**kwargs)
    return self.to_xarray(source="r",outer_dims=outer_dims)

  # -- instruments -----------------------------------------------------

  def add_instrument(self,instrument,run=True,outer_dims=None):
    '''
    Register and (by default) run a
    pyPamtra.instrument.PamtraInstrument -- see
    pyPamtra.core.pyPamtra.addInstrument. Access it afterwards via
    self.instruments[instrument.name].
    '''
    return self.pam.addInstrument(instrument,run=run,outer_dims=outer_dims)

  @property
  def instruments(self):
    '''dict of name -> PamtraInstrument, see add_instrument().'''
    return self.pam.instruments

  # -- export ------------------------------------------------------------

  def to_xarray(self,source="p",outer_dims=None):
    '''
    Explicit snapshot of the profile (source="p") or results
    (source="r") as an xarray.Dataset. See pyPamtra.core.pyPamtra.to_xarray.
    '''
    return self.pam.to_xarray(source=source,outer_dims=outer_dims)

  def to_netcdf(self,fname,source="r",outer_dims=None,**kwargs):
    '''
    Write to_xarray(source, outer_dims) to fname via
    xarray.Dataset.to_netcdf(). Defaults to the results (source="r");
    pass source="p" to write the profile instead.
    '''
    return self.to_xarray(source=source,outer_dims=outer_dims).to_netcdf(fname,**kwargs)
