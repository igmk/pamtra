"""Direct tests of a couple of Fortran routines that f2py already exposes
on the compiled pyPamtraLib extension (see f2py_sources in meson.build).

Most of src/ is built around shared vars_*/settings module state and isn't
a realistic unit-test target (see AI.md) -- but for the handful of
routines meson.build already wires up for f2py, testing them straight
through the compiled extension needs no new build infrastructure: no
separate driver programs, no extra meson test() targets, just calling the
same pyPamtraLib Python module the rest of pyPamtra already depends on.
"""

import numpy as np

from pyPamtra import pyPamtraLib


def test_viscosity_air():
    # Sutherland's law, no shared module state involved at all.
    eta = pyPamtraLib.viscosity_air(273.15)
    np.testing.assert_allclose(eta, 1.716743876889605e-05, rtol=1e-9)


def test_eps_water_tkc():
    # get_eps_water reads settings.liq_mod, so set it the same way
    # libWrapper.py does before calling into the Fortran core.
    pyPamtraLib.settings.liq_mod = "TKC"
    pyPamtraLib.report_module.verbose = 0

    errorstatus, eps = pyPamtraLib.eps_water.get_eps_water(0.0, 20.0, 10.0)

    assert errorstatus == 0
    np.testing.assert_allclose(eps.real, 60.6084, rtol=1e-3)
    np.testing.assert_allclose(eps.imag, 32.9184, rtol=1e-3)
