# -*- coding: utf-8 -*-
'''
Thin composition layer for running several differently-configured
"instruments" (a frequency set + nmlSet overrides) against the same
pyPamtra profile and descriptor file, each getting its own labeled
xarray.Dataset of results instead of sharing one flat pam.r dict that
gets overwritten by every run.

This is purely a convenience wrapper: PamtraInstrument.run() just calls
the existing pyPamtra.runPamtra() with a temporarily-overridden nmlSet
and snapshots pam.r via pyPamtra.to_xarray() -- no Fortran boundary or
storage changes.
'''


class PamtraInstrument(object):
  '''
  A named frequency set + nmlSet overrides, run against a pyPamtra
  profile via pam.addInstrument(instrument) (or instrument.run(pam)
  directly). Results are stored per-instrument in self.results as an
  xarray.Dataset (see pyPamtra.to_xarray), instead of overwriting the
  single shared pam.r dict.

  Parameters
  ----------
  name : str
      Name this instrument is registered under in pam.instruments.
  frequencies : float or list of float
      Frequencies in GHz to run this instrument at.
  **nmlOverrides
      nmlSet overrides specific to this instrument (e.g.
      radar_mode="moments", active=True, passive=False). Validated
      against the target pyPamtra object's own namelist defaults when
      run() is called, so a typo raises immediately instead of being
      silently ignored.

  Attributes
  ----------
  parent : pyPamtra or None
      The pyPamtra object this instrument is attached to. Set by
      pyPamtra.addInstrument(); can also be set directly.
  results : xarray.Dataset or None
      Populated by run(). Requires the optional 'xarray' package (see
      pyPamtra.to_xarray).

  Examples
  --------
  >>> pam = pyPamtra.pyPamtra()
  >>> pam.df.addHydrometeor(...)
  >>> pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0., 1000.])
  >>> pam.addInstrument(pyPamtra.PamtraInstrument("simpleRadar", 35.5, radar_mode="simple"))
  >>> pam.addInstrument(pyPamtra.PamtraInstrument("dopplerRadar", 35.5, radar_mode="spectrum"))
  >>> pam.instruments["simpleRadar"].results["Ze"]
  >>> pam.instruments["dopplerRadar"].results["radar_spectra"]
  '''

  def __init__(self,name,frequencies,**nmlOverrides):
    if not hasattr(frequencies,"__iter__"):
      frequencies = [frequencies]
    self.name = name
    self.frequencies = list(frequencies)
    self.nmlOverrides = nmlOverrides
    self.parent = None
    self.results = None

  def run(self,parent=None,outer_dims=None):
    '''
    Run this instrument's frequencies against parent (or self.parent,
    e.g. as set by pyPamtra.addInstrument()), with this instrument's
    nmlSet overrides applied only for the duration of this call --
    parent.nmlSet is restored afterwards even if the run fails, so
    instruments sharing one pyPamtra object never leak settings into
    each other. Populates and returns self.results.

    Parameters
    ----------
    parent : pyPamtra, optional
        Overrides self.parent for this call.
    outer_dims : dict, optional
        Passed through to pyPamtra.to_xarray().

    Returns
    -------
    xarray.Dataset
    '''
    if parent is not None:
      self.parent = parent
    if self.parent is None:
      raise RuntimeError(
        "instrument '%s' has no parent pyPamtra object -- call "
        "pam.addInstrument(instrument) or instrument.run(pam)" % self.name
      )

    unknown = set(self.nmlOverrides) - set(self.parent._nmlDefaultSet.keys())
    if unknown:
      raise TypeError(
        "instrument '%s' got unknown nmlSet override(s): %s"
        % (self.name,", ".join(sorted(unknown)))
      )

    savedNmlSet = dict(self.parent.nmlSet)
    try:
      self.parent.nmlSet.update(self.nmlOverrides)
      self.parent.runPamtra(self.frequencies)
      self.results = self.parent.to_xarray(source="r",outer_dims=outer_dims)
    finally:
      self.parent.nmlSet.clear()
      self.parent.nmlSet.update(savedNmlSet)

    return self.results

  def to_netcdf(self,fname,**kwargs):
    '''
    Write self.results to fname via xarray.Dataset.to_netcdf(). Requires
    run() to have been called first (either directly or via
    pam.addInstrument(instrument, run=True), the default).
    '''
    if self.results is None:
      raise RuntimeError("instrument '%s' has no results yet -- call run() first" % self.name)
    return self.results.to_netcdf(fname,**kwargs)
