..  _fortran_core:


Fortran Core
=============

The RT engine itself -- gas absorption, hydrometeor scattering (Mie,
T-matrix, self-similar Rayleigh-Gans), the RT4 passive solver, and the
radar Doppler spectrum simulator -- lives in ``src/`` as ~107 Fortran 90
files, compiled either into the standalone ``pamtra`` CLI binary or into
the ``pyPamtraLib`` extension that :ref:`pyPamtra` calls into (see
:ref:`installation` for how the two builds relate). It predates the
Python wrapper and is considered legacy: new functionality is added in
Python where possible, and the Fortran core is treated mostly as a fixed
computational backend rather than something to extend routine by
routine.

For the physics itself, the reference is Mech et al. (2020), *PAMTRA
1.0: the Passive and Active Microwave radiative TRAnsfer tool*, Geosci.
Model Dev., 13, 4229-4251, https://doi.org/10.5194/gmd-13-4229-2020 --
Sect. 2 in particular describes the model framework each src/ subsystem
implements.

For where in ``src/`` a given subsystem actually lives, and the (fairly
important) caveat that most of it is built around shared module-level
state rather than pure functions -- which shapes both how it's called
from Python (``pyPamtra.libWrapper``) and how it can realistically be
tested -- see `AI.md <https://github.com/igmk/pamtra/blob/master/AI.md>`_
in the repository root, which documents the Fortran source layout in
more detail than is useful to duplicate here.
