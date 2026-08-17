"""
pyPamtraXr quickstart: the recommended entry point for new code, a small
curated facade over pyPamtra. See doc/source/pyPamtraXr.rst and
doc/source/quickstart.rst for the reference docs this mirrors, and
examples/xarray_and_instruments_howto.py for the lower-level pieces
(addHydrometeor(**kwargs), to_xarray()/from_xarray(), PamtraInstrument)
this is built on.

Needs PAMTRA_DATADIR set (can be empty, e.g. `export PAMTRA_DATADIR=""`)
and the optional 'xarray' package (`pip install pamtra[xarray]`).
"""

import pyPamtra

pamxr = pyPamtra.pyPamtraXr()

# One liquid hydrometeor ("cwc_q"), exponential PSD, Mie-sphere
# scattering. See doc/source/descriptorFile.rst for what each field means.
pamxr.add_hydrometeor(
    hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
    dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
    scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
)

# A hand-built profile needs no external data (real work usually starts
# from a model output importer instead, see doc/source/profiles.rst).
pamxr.set_profile(
    hgt_lev=[0.0, 1000.0, 2000.0],
    temp_lev=[288.15, 281.5, 275.0],
    press_lev=[101300.0, 89000.0, 79000.0],
    relhum_lev=[80.0, 70.0, 60.0],
)

# Put some liquid water content into the lowest layer. Not yet covered by
# pyPamtraXr's curated surface, so this reaches into the wrapped pyPamtra
# object directly via the "escape hatch" -- always the same object, never
# a copy, so the edit is immediately visible to pamxr too.
pamxr.pam.p["hydro_q"][0, 0, 0, 0] = 1e-3

results = pamxr.run(35.5)  # frequency in GHz

print("brightness temperature [K]:", results["tb"].values.ravel())
print("radar reflectivity [dBz]:", results["Ze"].values.ravel())

# --------------------------------------------------------------------
# Multiple instrument configurations against the same profile, each
# keeping its own results (see doc/source/running.rst).
# --------------------------------------------------------------------

simple = pamxr.add_instrument(
    pyPamtra.PamtraInstrument("simpleRadar", 35.5, radar_mode="simple")
)
spectrum = pamxr.add_instrument(
    pyPamtra.PamtraInstrument("dopplerRadar", 35.5, radar_mode="spectrum")
)

print()
print("simpleRadar Ze [dBz]:", pamxr.instruments["simpleRadar"].results["Ze"].values.ravel())
print("dopplerRadar spectrum shape:", pamxr.instruments["dopplerRadar"].results["radar_spectra"].shape)

pamxr.to_netcdf("pyPamtraXr_quickstart_results.nc")
print()
print("wrote pyPamtraXr_quickstart_results.nc")
