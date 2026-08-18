# -*- coding: utf-8 -*-
'''
pyPamtraXr: a curated, xarray/pandas-native facade over pyPamtra, for
users who would rather not touch pam.p["key"][x,y,z]-style dict/array
indexing or pam.df's positional-tuple recarray rows directly.

Unlike an earlier version of this class, .p/.r/.df/.df_4d/.df_full_spec
are real, persistent attributes here -- not hidden, not export-on-demand.
There is still exactly one place the RT engine itself runs: self.pam, a
real pyPamtra object, always accessible as the escape hatch for anything
this class doesn't cover. run() (and add_instrument()) translate .p into
self.pam.p right before calling self.pam.runPamtra()/addInstrument(), via
the existing, already-tested pyPamtra.from_xarray() -- not a new,
independently-maintained reimplementation of createProfile's shape/
defaulting logic. .df/.df_4d/.df_full_spec instead stay eagerly synced
from self.pam.df on every add_hydrometeor()/remove_hydrometeor()/
add_4d()/remove_4d()/add_full_spectra() call, by calling the existing,
already-tested pamDescriptorFile methods and re-deriving the pandas/xr
mirror from the result -- so there is exactly one implementation of
"what does a valid hydrometeor row look like", not two.

If self.pam is mutated directly through the escape hatch, .p/.r/.df/
.df_4d/.df_full_spec go stale until refresh() (or the next run()) is
called -- the same category of caveat as any cached view of mutable
state.
'''

import numpy as np

from .core import pyPamtra

# canonical dims, matching pyPamtra.core._CANONICAL_DIM_NAMES
_DF_4D_DIMS = ["grid_x", "grid_y", "heightbins", "hydrometeor"]
_DF_FULL_SPEC_DIMS = ["grid_x", "grid_y", "heightbins", "hydrometeor", "sizebin"]


def _require_xarray_pandas():
  try:
    import xarray  # noqa: F401
    import pandas  # noqa: F401
  except ImportError as e:
    raise ImportError("pyPamtraXr requires the optional 'xarray' and 'pandas' packages: pip install pamtra[xarray]") from e


def _descriptor_data_to_dataframe(data):
  import pandas as pd
  return pd.DataFrame(data).set_index("hydro_name")


def _dict_to_dataset(d, dims):
  import xarray as xr
  data_vars = {}
  for key, arr in d.items():
    arr = np.asarray(arr)
    these_dims = list(dims) if len(dims) == arr.ndim else ["%s_dim%i" % (key,i) for i in range(arr.ndim)]
    data_vars[key] = (these_dims,arr.copy())
  return xr.Dataset(data_vars)


def _full_spec_to_dataset(d):
  import xarray as xr
  data_vars = {}
  for key, arr in d.items():
    arr = np.asarray(arr)
    dims = list(_DF_FULL_SPEC_DIMS)
    if key == "d_bound_ds":
      # one bin longer than the rest: bin boundaries, not bin centers
      dims[-1] = "sizebin_plus1"
    if len(dims) != arr.ndim:
      dims = ["%s_dim%i" % (key,i) for i in range(arr.ndim)]
    data_vars[key] = (dims,arr.copy())
  return xr.Dataset(data_vars)


class pyPamtraXr(object):
  '''
  Curated, xarray/pandas-native facade over a pyPamtra object.

  Requires the optional 'xarray' package (which pulls in pandas as a
  dependency of its own).

  Parameters
  ----------
  pam : pyPamtra, optional
      Wrap an existing pyPamtra object instead of creating a fresh one
      (see also from_pam(), the named constructor for this).
  outer_dims : dict, optional
      Rename map applied to the leading grid dimensions of every
      Dataset this object builds (.p, .r, and instrument results), e.g.
      {"grid_x": "lat", "grid_y": "lon"} to label them something more
      meaningful than the generic grid index -- see
      pyPamtra.core.pyPamtra.to_xarray. Stored as self.outer_dims;
      change it directly and call refresh() to re-label existing state.

  Attributes
  ----------
  pam : pyPamtra
      The wrapped pyPamtra object -- always the same one, never a copy.
      The escape hatch for anything not covered by this class's curated
      methods (nmlSet/set tweaks, or any other pyPamtra method), and the
      only place the RT engine actually runs.
  outer_dims : dict
      See the outer_dims parameter above.
  p : xarray.Dataset
      The atmospheric profile. See set_profile()/set_profile_from_xarray()
      to populate it with validation, or assign/mutate it directly for
      full control (no validation gate on direct assignment).
  r : xarray.Dataset
      Results, populated by run()/run_parallel() (empty before the first
      run).
  df : pandas.DataFrame
      Per-hydrometeor scalar properties, indexed by hydro_name -- see
      doc/source/descriptorFile.rst for what each column means. Kept in
      sync with self.pam.df by add_hydrometeor()/remove_hydrometeor();
      read-mutate this directly only if you also keep self.pam.df in
      sync yourself (or call refresh() afterwards).
  df_4d : xarray.Dataset
      Per-gridpoint hydrometeor property overrides, see
      pyPamtra.descriptorFile.pamDescriptorFile.add4D. Kept in sync by
      add_4d()/remove_4d().
  df_full_spec : xarray.Dataset
      Full particle-size-spectrum arrays, see
      pyPamtra.descriptorFile.pamDescriptorFile.addFullSpectra. Kept in
      sync by add_full_spectra().

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

  def __init__(self,pam=None,outer_dims=None):
    _require_xarray_pandas()
    self.pam = pam if pam is not None else pyPamtra()
    self.outer_dims = outer_dims or {}
    self.refresh()

  @classmethod
  def from_pam(cls,pam,outer_dims=None):
    '''
    Build a pyPamtraXr from an existing pyPamtra object -- e.g. one
    built by one of the pyPamtra.importer functions, so the existing
    importers keep working unmodified. pam becomes the new object's
    escape hatch (self.pam), not a copy.
    '''
    return cls(pam=pam,outer_dims=outer_dims)

  def refresh(self):
    '''
    Resync .p/.r/.df/.df_4d/.df_full_spec from the current state of
    .pam (using the current self.outer_dims). Only needed after
    mutating .pam directly through the escape hatch, or after changing
    self.outer_dims -- every method on this class keeps .p/.r in sync
    with both on its own otherwise.

    Returns
    -------
    self
    '''
    self.p = self.pam.to_xarray(source="p",outer_dims=self.outer_dims)
    self.r = self.pam.to_xarray(source="r",outer_dims=self.outer_dims)
    self._sync_df()
    return self

  def _sync_df(self):
    self.df = _descriptor_data_to_dataframe(self.pam.df.data)
    self.df_4d = _dict_to_dataset(self.pam.df.data4D,_DF_4D_DIMS)
    self.df_full_spec = _full_spec_to_dataset(self.pam.df.dataFullSpec)

  # -- hydrometeors --------------------------------------------------

  def add_hydrometeor(self,**kwargs):
    '''
    Add a hydrometeor by keyword arguments. See
    pyPamtra.descriptorFile.pamDescriptorFile.addHydrometeor for the
    field reference. Updates self.df (and self.pam.df, which is what
    actually gets used at run() time).
    '''
    self.pam.df.addHydrometeor(**kwargs)
    self._sync_df()
    return self.df

  def remove_hydrometeor(self,name):
    '''Remove a hydrometeor previously added via add_hydrometeor().'''
    self.pam.df.removeHydrometeor(name)
    self._sync_df()
    return self.df

  def add_4d(self,key,arr):
    '''
    Replace hydrometeor property `key` with per-gridpoint values. See
    pyPamtra.descriptorFile.pamDescriptorFile.add4D. Updates
    self.df_4d.
    '''
    self.pam.df.add4D(key,arr)
    self._sync_df()
    return self.df_4d

  def remove_4d(self,key,val):
    '''Undo add_4d(): restore `key` to a single scalar value per hydrometeor.'''
    self.pam.df.remove4D(key,val)
    self._sync_df()
    return self.df_4d

  def add_full_spectra(self):
    '''
    Switch to full particle-size-spectrum input. See
    pyPamtra.descriptorFile.pamDescriptorFile.addFullSpectra. Updates
    self.df_full_spec (fill it in before calling run()).
    '''
    self.pam.df.addFullSpectra()
    self._sync_df()
    return self.df_full_spec

  # -- profile ---------------------------------------------------------

  def set_profile(self,**kwargs):
    '''
    Define the atmospheric profile with validation/defaulting. See
    pyPamtra.core.pyPamtra.createProfile for the accepted keyword
    arguments (hgt_lev/temp_lev/press_lev/relhum_lev are mandatory,
    everything else is optional and defaulted with a warning). Updates
    self.p (and self.pam.p).
    '''
    self.pam.createProfile(**kwargs)
    self.p = self.pam.to_xarray(source="p",outer_dims=self.outer_dims)
    return self.p

  def set_profile_from_xarray(self,ds,outer_dims=None):
    '''
    Define the atmospheric profile from an xarray.Dataset (e.g. one
    assembled from a model's own xarray-native output), with the same
    validation/defaulting as set_profile(). Hydrometeors must already
    be added via add_hydrometeor() first if ds includes hydro_q/
    hydro_n/hydro_reff. Updates self.p (and self.pam.p).

    To skip validation/defaulting entirely (e.g. ds is already known to
    have every field set), assign self.p = ds directly instead.

    Parameters
    ----------
    ds : xarray.Dataset
    outer_dims : dict, optional
        How ds's own dimensions relate to the canonical names (undone
        before handing off to pyPamtra.from_xarray()) -- not the same
        thing as self.outer_dims, which is how the *resulting* self.p
        gets labeled once ds has been read in.
    '''
    self.pam.from_xarray(ds,outer_dims=outer_dims)
    self.p = self.pam.to_xarray(source="p",outer_dims=self.outer_dims)
    return self.p

  # -- running -----------------------------------------------------------

  def run(self,freqs,outer_dims=None,**kwargs):
    '''
    Run the RT engine (pyPamtra.core.pyPamtra.runPamtra). self.p is
    translated into self.pam.p first (via pyPamtra.from_xarray(), so
    the same defaulting/validation as set_profile_from_xarray() applies
    here too); self.pam.df is assumed already in sync (add_hydrometeor()
    etc. keep it that way). Populates and returns self.r.

    Parameters
    ----------
    freqs : float or list of float
        Frequencies in GHz.
    outer_dims : dict, optional
        Overrides self.outer_dims for this call, for both reading self.p
        back in and labeling the returned results.
    **kwargs
        Passed through to runPamtra() (e.g. checkData=False).

    Returns
    -------
    xarray.Dataset
    '''
    outer_dims = self.outer_dims if outer_dims is None else outer_dims
    self.pam.from_xarray(self.p,outer_dims=outer_dims)
    self.pam.runPamtra(freqs,**kwargs)
    self.r = self.pam.to_xarray(source="r",outer_dims=outer_dims)
    return self.r

  def run_parallel(self,freqs,outer_dims=None,**kwargs):
    '''
    Run the RT engine across local CPU cores
    (pyPamtra.core.pyPamtra.runParallelPamtra). See run() for how self.p
    is translated in and what outer_dims does; the keyword arguments
    also accept runParallelPamtra's
    pp_local_workers/pp_deltaF/pp_deltaX/pp_deltaY/timeout. Populates
    and returns self.r.

    Returns
    -------
    xarray.Dataset
    '''
    outer_dims = self.outer_dims if outer_dims is None else outer_dims
    self.pam.from_xarray(self.p,outer_dims=outer_dims)
    self.pam.runParallelPamtra(freqs,**kwargs)
    self.r = self.pam.to_xarray(source="r",outer_dims=outer_dims)
    return self.r

  # -- instruments -----------------------------------------------------

  def add_instrument(self,instrument,run=True,outer_dims=None):
    '''
    Register and (by default) run a
    pyPamtra.instrument.PamtraInstrument against the current self.p --
    see pyPamtra.core.pyPamtra.addInstrument. Access it afterwards via
    self.instruments[instrument.name]. self.p is translated into
    self.pam.p first, same as run(); outer_dims defaults to
    self.outer_dims, same meaning as in run().
    '''
    outer_dims = self.outer_dims if outer_dims is None else outer_dims
    self.pam.from_xarray(self.p,outer_dims=outer_dims)
    return self.pam.addInstrument(instrument,run=run,outer_dims=outer_dims)

  @property
  def instruments(self):
    '''dict of name -> PamtraInstrument, see add_instrument().'''
    return self.pam.instruments

  # -- export ------------------------------------------------------------

  def to_netcdf(self,fname,source="r",**kwargs):
    '''
    Write self.p or self.r to fname via xarray.Dataset.to_netcdf().
    Defaults to the results (source="r"); pass source="p" to write the
    profile instead.
    '''
    if source not in ("p","r"):
      raise ValueError("source must be 'p' or 'r', got %r" % (source,))
    return getattr(self,source).to_netcdf(fname,**kwargs)
