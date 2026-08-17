"""
pyPamtraXr quickstart: the recommended entry point for new code, a small
curated facade over pyPamtra. .p/.r are real xarray.Dataset, .df is a
pandas.DataFrame, .df_4d/.df_full_spec are xarray.Dataset -- see
doc/source/pyPamtraXr.rst and doc/source/quickstart.rst for the
reference docs this mirrors, and examples/xarray_and_instruments_howto.py
for the lower-level pieces (addHydrometeor(**kwargs), to_xarray()/
from_xarray(), PamtraInstrument) this is built on.

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
print("pamxr.df (per-hydrometeor scalar properties):")
print(pamxr.df)

# A hand-built profile needs no external data (real work usually starts
# from a model output importer instead, see "from_pam" below).
pamxr.set_profile(
    hgt_lev=[0.0, 1000.0, 2000.0],
    temp_lev=[288.15, 281.5, 275.0],
    press_lev=[101300.0, 89000.0, 79000.0],
    relhum_lev=[80.0, 70.0, 60.0],
)

# pamxr.p is a real xarray.Dataset -- mutate it in place like any other.
pamxr.p["hydro_q"].values[0, 0, 0, 0] = 1e-3

results = pamxr.run(35.5)  # frequency in GHz

print()
print("brightness temperature [K]:", results["tb"].values.ravel())
print("radar reflectivity [dBz]:", results["Ze"].values.ravel())

# --------------------------------------------------------------------
# from_pam(): bridge an existing pyPamtra.importer function's output in
# unmodified -- every importer builds and returns a plain pyPamtra
# object, and from_pam() wraps it as pamxr.pam (not a copy).
# --------------------------------------------------------------------

pam = pyPamtra.importer.createUsStandardProfile(pyPamtra.pyPamtra(), hgt_lev=[0.0, 500.0, 1500.0])
pam.df.addHydrometeor(
    hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
    dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
    scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
)
pam.p["hydro_q"][0, 0, 0, 0] = 1e-3

pamxr2 = pyPamtra.pyPamtraXr.from_pam(pam)
print()
print("from_pam: pamxr2.pam is pam:", pamxr2.pam is pam)
print("from_pam: pamxr2.p has the imported profile:", "hgt_lev" in pamxr2.p)

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
