"""
How-to for three additive pyPamtra features that don't change any
existing behavior: keyword-argument hydrometeor definitions, the
xarray interface, and running several instrument configurations
against one profile. See doc/source/descriptorFile.rst,
doc/source/results.rst, and doc/source/running.rst for the reference
docs this mirrors.

Needs PAMTRA_DATADIR set (can be empty, e.g. `export PAMTRA_DATADIR=""`)
and the optional 'xarray' package (`pip install pamtra[xarray]`).
"""

import pyPamtra

# --------------------------------------------------------------------
# 1. addHydrometeor(**kwargs): same descriptor file row as the old
#    positional-tuple form, but fields are named, so a typo or a
#    transposed value raises instead of silently landing in the wrong
#    column. See doc/source/descriptorFile.rst for what each field
#    means.
# --------------------------------------------------------------------

pam = pyPamtra.pyPamtra()

# old style, still fully supported (equivalent to the addHydrometeor() call below):
# pam.df.addHydrometeor((
#     "cwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
#     3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
#     "mie-sphere", "khvorostyanov01_drops", -99.0,
# ))

pam.df.addHydrometeor(
    hydro_name="cwc_q", liq_ice=1, moment_in=3, nbin=30,
    dist_name="exp", p_3=8.0e6, d_1=1.0e-5, d_2=6.0e-3,
    scat_name="mie-sphere", vel_size_mod="khvorostyanov01_drops",
)

pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0.0, 1000.0, 2000.0])
pam.p["hydro_q"][0, 0, 0, 0] = 1e-3

# --------------------------------------------------------------------
# 2. to_xarray()/from_xarray(): a labeled, unit-aware xarray.Dataset
#    snapshot of pam.p or pam.r, built from metadata pyPamtra already
#    carries (self.dimensions/self.units). Purely a snapshot -- it
#    never holds a live view into pam.p/pam.r, so nothing you do to
#    the returned Dataset can affect the pyPamtra object.
# --------------------------------------------------------------------

profile_ds = pam.to_xarray(source="p")
print("profile as xarray.Dataset:")
print(profile_ds)
print()
print("hgt_lev units:", profile_ds["hgt_lev"].attrs["units"])

# and back, e.g. after building/editing a profile the xarray way:
pam_from_ds = pyPamtra.pyPamtra()
pam_from_ds.df = pam.df  # descriptor file is not part of the Dataset
pam_from_ds.from_xarray(profile_ds)

# --------------------------------------------------------------------
# 3. PamtraInstrument: run several frequency/nmlSet configurations
#    against the same profile, each keeping its own results instead of
#    overwriting the single shared pam.r. nmlSet overrides apply only
#    for the duration of that instrument's run and are restored
#    afterwards, even if the run fails.
# --------------------------------------------------------------------

simple = pam.addInstrument(
    pyPamtra.PamtraInstrument("simpleRadar", 35.5, radar_mode="simple")
)
spectrum = pam.addInstrument(
    pyPamtra.PamtraInstrument("dopplerRadar", 35.5, radar_mode="spectrum")
)

print()
print("simpleRadar Ze [dBz]:", simple.results["Ze"].values.ravel())
print("dopplerRadar spectrum shape:", spectrum.results["radar_spectra"].shape)

simple.to_netcdf("instrument_howto_simpleRadar.nc")
print("wrote instrument_howto_simpleRadar.nc")
