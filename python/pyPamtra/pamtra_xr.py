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

self.outer_dims (an ordered list of names, e.g. ["time","lat","lon"])
lets the leading "grid" dimensions be named -- and counted -- arbitrarily,
even though pyPamtra's own core is genuinely always a 2D (ngridx, ngridy)
grid. Two names are a direct, lossless rename (pyPamtra's grid already
is 2D, nothing needs reshaping). Any other count is handled by folding
the named axes into pyPamtra's single ngridx axis (ngridy=1) before
handing off to self.pam, and splitting it back apart when reading .p/.r/
instrument results back out -- entirely within this class, never
touching pyPamtra.core.pyPamtra.to_xarray()/from_xarray() or the Fortran
boundary, both of which stay genuinely, natively 2D.
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


def _normalize_outer_dims(outer_dims):
  '''
  outer_dims is an ordered list of names, e.g. ["time","lat","lon"] --
  not the {"grid_x": "lat", "grid_y": "lon"} rename dict the pre-N-D
  version of this class used. A dict is rejected explicitly rather than
  silently accepted: list(a_dict) returns its keys, so passing one by
  mistake would otherwise quietly turn into a no-op self-rename instead
  of raising.
  '''
  if isinstance(outer_dims,dict):
    raise TypeError(
      "outer_dims must be an ordered list of dimension names (e.g. "
      "['lat','lon']), not a dict -- got %r" % (outer_dims,)
    )
  return list(outer_dims)


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


# -- outer-dims flatten/unflatten -------------------------------------
#
# pyPamtra.core.pyPamtra is genuinely, always a 2D (ngridx, ngridy) grid
# (see PamtraFortranWrapper setting vars_atmosphere.atmo_ngridx/atmo_ngridy
# directly on the Fortran module) -- there is no way to give it more or
# fewer than 2 outer dims. len(outer_dims) == 2 is therefore always a
# direct, lossless rename. Anything else (1, or 3+) is a real reshape:
# the outer_dims axes get folded into a single leading axis (pyPamtra's
# ngridx, with ngridy left at 1) on the way in, and split back apart
# (using outer_shape, the sizes remembered from whichever flatten call
# established them) on the way out.


def _flatten_dataset(ds,outer_dims):
  '''
  Fold ds's outer_dims (an ordered list of dimension names) into the
  literal grid_x[,grid_y] dims pyPamtra.core.pyPamtra actually
  understands.

  Returns (flattened_ds, outer_shape), where outer_shape is the tuple
  of original sizes along outer_dims (needed to reverse this via
  _unflatten_dataset() later), or None when len(outer_dims) <= 2 (a
  pure rename needs no shape to reverse).
  '''
  import xarray as xr
  outer_dims = list(outer_dims)
  if len(outer_dims) == 2:
    rename = {d: n for d,n in zip(outer_dims,["grid_x","grid_y"]) if d in ds.dims}
    return ds.rename(rename), None
  if len(outer_dims) == 0:
    return ds, None
  if len(outer_dims) == 1:
    rename = {outer_dims[0]: "grid_x"} if outer_dims[0] in ds.dims else {}
    return ds.rename(rename), None

  outer_shape = tuple(ds.sizes[d] for d in outer_dims)
  n = int(np.prod(outer_shape,dtype=int))
  data_vars = {}
  for key in ds.data_vars:
    da = ds[key]
    if not all(d in da.dims for d in outer_dims):
      data_vars[key] = da
      continue
    otherDims = [d for d in da.dims if d not in outer_dims]
    daT = da.transpose(*outer_dims,*otherDims)
    flatValues = daT.values.reshape((n,)+daT.shape[len(outer_dims):])
    data_vars[key] = (["grid_x"]+otherDims,flatValues,da.attrs)
  return xr.Dataset(data_vars), outer_shape


def _flatten_kwargs(kwargs,outer_dims):
  '''
  Same idea as _flatten_dataset(), but for plain-array kwargs (e.g. as
  passed to createProfile()/set_profile()), which have no dimension
  names to match against -- only the leading len(outer_dims) axes of
  each array (by position) are folded together.
  '''
  n = len(outer_dims)
  if n <= 2:
    return kwargs, None
  outer_shape = None
  flat_kwargs = {}
  for key,val in kwargs.items():
    arr = np.asarray(val)
    if arr.ndim < n:
      flat_kwargs[key] = val
      continue
    thisShape = arr.shape[:n]
    if outer_shape is None:
      outer_shape = thisShape
    elif thisShape != outer_shape:
      raise ValueError(
        "all profile arrays must share the same leading %i-dim shape for "
        "outer_dims %r; got %r for '%s', expected %r" % (n,outer_dims,thisShape,key,outer_shape)
      )
    flat_kwargs[key] = arr.reshape((int(np.prod(thisShape,dtype=int)),)+arr.shape[n:])
  return flat_kwargs, outer_shape


def _unflatten_dataset(ds,outer_dims,outer_shape):
  '''Inverse of _flatten_dataset()/_flatten_kwargs().'''
  import xarray as xr
  outer_dims = list(outer_dims)
  if len(outer_dims) == 2:
    rename = {d: n for d,n in zip(["grid_x","grid_y"],outer_dims) if d in ds.dims}
    return ds.rename(rename)
  if len(outer_dims) <= 1:
    rename = {"grid_x": outer_dims[0]} if outer_dims and "grid_x" in ds.dims else {}
    out = ds.rename(rename) if rename else ds
    return out.squeeze("grid_y",drop=True) if "grid_y" in out.dims else out

  hasGridX = any("grid_x" in ds[key].dims for key in ds.data_vars)
  if not hasGridX:
    return ds
  if outer_shape is None:
    raise RuntimeError(
      "cannot restore %i-dimensional outer_dims %r: no profile has been "
      "set with this many outer dims yet -- call set_profile()/"
      "set_profile_from_xarray() first" % (len(outer_dims),outer_dims)
    )
  outer_shape = tuple(int(s) for s in outer_shape)

  data_vars = {}
  for key in ds.data_vars:
    da = ds[key]
    if "grid_x" not in da.dims:
      data_vars[key] = da
      continue
    restDims = [d for d in da.dims if d not in ("grid_x","grid_y")]
    restShape = tuple(s for d,s in zip(da.dims,da.shape) if d not in ("grid_x","grid_y"))
    values = da.values.reshape(outer_shape+restShape)
    data_vars[key] = (list(outer_dims)+restDims,values,da.attrs)
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
  outer_dims : list of str, optional
      Names for the leading "grid" dimensions of every Dataset this
      object builds (.p, .r, and instrument results) -- e.g.
      ["time", "lat", "lon"] instead of the generic default
      ["grid_x", "grid_y"]. Any number of names is allowed: 2 is always
      a direct rename (pyPamtra's own grid is genuinely 2D); any other
      count folds/splits the extra dimensionality around that 2D core
      (see the module docstring). Stored as self.outer_dims; change it
      directly and call refresh() to re-label existing state.

  Attributes
  ----------
  pam : pyPamtra
      The wrapped pyPamtra object -- always the same one, never a copy.
      The escape hatch for anything not covered by this class's curated
      methods (nmlSet/set tweaks, or any other pyPamtra method), and the
      only place the RT engine actually runs.
  outer_dims : list of str
      See the outer_dims parameter above.

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
    self.outer_dims = _normalize_outer_dims(outer_dims) if outer_dims is not None else ["grid_x","grid_y"]
    self._outer_shape = None
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

    Raises RuntimeError if self.outer_dims has more than 2 names and no
    profile has ever been set with that many outer dims (there is then
    no way to know how to split .pam's plain 2D grid back apart).

    Returns
    -------
    self
    '''
    self.p = _unflatten_dataset(self.pam.to_xarray(source="p"),self.outer_dims,self._outer_shape)
    self.r = _unflatten_dataset(self.pam.to_xarray(source="r"),self.outer_dims,self._outer_shape)
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
    self.df_4d. arr must already be shaped for pyPamtra's native 2D
    grid (self.pam._shape4D) -- unlike .p, self.outer_dims folding is
    not applied here.
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

    Arrays may have len(self.outer_dims) leading axes instead of the 1
    or 2 pyPamtra.core.pyPamtra.createProfile natively understands --
    see the module docstring for how that's reconciled.
    '''
    kwargs,outer_shape = _flatten_kwargs(kwargs,self.outer_dims)
    if outer_shape is not None:
      self._outer_shape = outer_shape
    self.pam.createProfile(**kwargs)
    self.p = _unflatten_dataset(self.pam.to_xarray(source="p"),self.outer_dims,self._outer_shape)
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
    outer_dims : list of str, optional
        ds's own outer dimension names, if different from
        self.outer_dims (defaults to self.outer_dims).
    '''
    outer_dims = self.outer_dims if outer_dims is None else _normalize_outer_dims(outer_dims)
    flatDs,outerShape = _flatten_dataset(ds,outer_dims)
    if outerShape is not None:
      self._outer_shape = outerShape
    self.pam.from_xarray(flatDs)
    self.p = _unflatten_dataset(self.pam.to_xarray(source="p"),self.outer_dims,self._outer_shape)
    return self.p

  # -- running -----------------------------------------------------------

  def run(self,freqs,outer_dims=None,**kwargs):
    '''
    Run the RT engine (pyPamtra.core.pyPamtra.runPamtra). self.p is
    translated into self.pam.p first (folding self.outer_dims down to
    pyPamtra's native 2D grid, then via pyPamtra.from_xarray(), so the
    same defaulting/validation as set_profile_from_xarray() applies
    here too); self.pam.df is assumed already in sync (add_hydrometeor()
    etc. keep it that way). Populates and returns self.r.

    Parameters
    ----------
    freqs : float or list of float
        Frequencies in GHz.
    outer_dims : list of str, optional
        Overrides self.outer_dims for this call, for both reading self.p
        back in and labeling the returned results.
    **kwargs
        Passed through to runPamtra() (e.g. checkData=False).

    Returns
    -------
    xarray.Dataset
    '''
    outer_dims = self.outer_dims if outer_dims is None else _normalize_outer_dims(outer_dims)
    flatP,outerShape = _flatten_dataset(self.p,outer_dims)
    if outerShape is not None:
      self._outer_shape = outerShape
    self.pam.from_xarray(flatP)
    self.pam.runPamtra(freqs,**kwargs)
    self.r = _unflatten_dataset(self.pam.to_xarray(source="r"),outer_dims,self._outer_shape)
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
    outer_dims = self.outer_dims if outer_dims is None else _normalize_outer_dims(outer_dims)
    flatP,outerShape = _flatten_dataset(self.p,outer_dims)
    if outerShape is not None:
      self._outer_shape = outerShape
    self.pam.from_xarray(flatP)
    self.pam.runParallelPamtra(freqs,**kwargs)
    self.r = _unflatten_dataset(self.pam.to_xarray(source="r"),outer_dims,self._outer_shape)
    return self.r

  # -- instruments -----------------------------------------------------

  def add_instrument(self,instrument,run=True,outer_dims=None):
    '''
    Register and (by default) run a
    pyPamtra.instrument.PamtraInstrument against the current self.p --
    see pyPamtra.core.pyPamtra.addInstrument. Access it afterwards via
    self.instruments[instrument.name]. self.p is translated into
    self.pam.p first, same as run(); outer_dims defaults to
    self.outer_dims, same meaning as in run(). If run=True,
    instrument.results is also unfolded back to self.outer_dims before
    returning (if run=False, the deferred instrument.run() call will
    leave instrument.results in pyPamtra's native grid_x/grid_y form --
    call refresh()-style unflatten yourself, or just call run() now).
    '''
    outer_dims = self.outer_dims if outer_dims is None else _normalize_outer_dims(outer_dims)
    flatP,outerShape = _flatten_dataset(self.p,outer_dims)
    if outerShape is not None:
      self._outer_shape = outerShape
    self.pam.from_xarray(flatP)
    self.pam.addInstrument(instrument,run=run,outer_dims=None)
    if run and instrument.results is not None:
      instrument.results = _unflatten_dataset(instrument.results,outer_dims,self._outer_shape)
    return instrument

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
