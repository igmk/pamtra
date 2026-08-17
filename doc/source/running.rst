..  _running:


Running a Simulation
======================

Once a profile (:ref:`profiles`) and descriptor file (:ref:`descriptorFile`)
are set up, and :ref:`settings` (``pam.nmlSet``/``pam.set``) are configured,
one of four methods on :any:`pyPamtra.core.pyPamtra` actually runs the RT engine
and populates the result dict ``pam.r`` (:ref:`results`). All four take a
``freqs`` argument -- a single frequency or a sorted list of frequencies in
GHz.

``runPamtra``
****************

The default, serial entry point::

    pam.runPamtra([10.65, 35.5, 94.0])

Use this unless the run is large enough that parallelism actually matters
-- it's the simplest to debug.

``runParallelPamtra``
************************

Splits the run across local CPU cores with :py:mod:`multiprocessing`, by
frequency and/or grid point::

    pam.runParallelPamtra(
        [10.65, 35.5, 94.0],
        pp_local_workers="auto",  # or an explicit int
        pp_deltaF=1,   # frequencies per worker chunk (0 = all in one)
        pp_deltaX=0,   # x-grid points per chunk (0 = all in one)
        pp_deltaY=0,   # y-grid points per chunk (0 = all in one)
        timeout=None,  # seconds per chunk, or None
    )

Use this for a multicore workstation and a profile with enough grid
points/frequencies to be worth splitting. The ``pp_delta*`` arguments
control the granularity of the split (small chunks -> more, smaller jobs).

``runPicklePamtra`` / ``runPicklePamtraSFTP``
*************************************************

Cluster-oriented variants of ``runParallelPamtra``: instead of running
in-process, each chunk is pickled to ``picklePath`` and expected to be
picked up and executed by separate jobs (e.g. submitted to a batch
scheduler), with ``runPicklePamtraSFTP`` transferring the pickles to a
remote host over SFTP first via :any:`pyPamtra.tools.sftp2Cluster`. Use
these only if a single machine's cores (``runParallelPamtra``) aren't
enough -- they add real operational complexity (job submission, pickle
lifecycle) that plain parallel execution doesn't need.

Checking input first
***********************

All four methods run ``pam._checkData()`` first unless
called with ``checkData=False``, which validates array shapes and value
ranges in the profile before handing off to the compiled extension --
worth leaving on, since a shape mismatch caught here is much easier to
debug than the same problem surfacing as a Fortran-side crash or garbage
output.

Running several instrument configurations
*********************************************

Each ``run*`` call above overwrites the single, shared ``pam.r`` --
running a second configuration (different frequencies, or a different
``radar_mode``/``active``/``passive`` combination) against the same
profile discards the first result unless you save it yourself first.
:any:`pyPamtra.PamtraInstrument` wraps that bookkeeping: a name, a
frequency list, and any ``nmlSet`` overrides specific to that
configuration, run via :any:`pyPamtra.core.pyPamtra.addInstrument`::

    pam.addInstrument(pyPamtra.PamtraInstrument("simple", 35.5, radar_mode="simple"))
    pam.addInstrument(pyPamtra.PamtraInstrument("spectrum", 35.5, radar_mode="spectrum"))

    pam.instruments["simple"].results["Ze"]
    pam.instruments["spectrum"].results["radar_spectra"]

Each instrument's ``nmlSet`` overrides apply only for the duration of
its own run (``pam.nmlSet`` is restored afterwards, even if the run
fails), and each instrument keeps its own :any:`xarray.Dataset` of
results (see :ref:`results`'s xarray section) rather than sharing
``pam.r`` -- so the two calls above don't interfere with each other.
Pass ``run=False`` to register an instrument without running it yet, and
call ``instrument.run()`` explicitly later.
