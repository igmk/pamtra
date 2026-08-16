"""Round-trip tests for writePamtraProfile / readPamtraProfile (the .lev/.lay
ASCII profile format).

Pure Python, no compiled pamtra binary or PAMTRA_DATADIR needed -- unlike
test_cli_binary.py, which only exercises this pair for the .lev format
with moment_in==3, these run unconditionally and cover every moment_in
value (0, 1, 2, 3, 12, 13, 23) in both .lev and .lay, since writing the
CLI test exposed real bugs here (see git history) that had gone
undetected for exactly the cases this file doesn't cover.
"""

import numpy as np
import pytest

import pyPamtra

MOMENT_CASES = [0, 1, 2, 3, 12, 13, 23]

N_VALUE = 1.23e6
REFF_VALUE = 4.56e-4
Q_VALUE = 7.89e-4


def make_pam(moment_in):
  pam = pyPamtra.pyPamtra()
  pam.df.addHydrometeor((
      "test_hyd", -99.0, -1, -99.0, -99.0, -99.0, -99.0, -99.0,
      moment_in, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
      "mie-sphere", "khvorostyanov01_drops", -99.0,
  ))
  pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=[0.0, 1000.0, 2000.0])
  pam.set["pyVerbose"] = 0
  # createUsStandardProfile only supplies level-centered fields; the .lay
  # writer needs the layer-centered ones too (same midpoint average
  # createProfile itself uses for "hgt" -- see core.py).
  pam.p["hgt"] = (pam.p["hgt_lev"][..., 1:] + pam.p["hgt_lev"][..., :-1]) / 2.
  pam.p["press"] = (pam.p["press_lev"][..., 1:] + pam.p["press_lev"][..., :-1]) / 2.
  pam.p["temp"] = (pam.p["temp_lev"][..., 1:] + pam.p["temp_lev"][..., :-1]) / 2.
  pam.p["relhum"] = (pam.p["relhum_lev"][..., 1:] + pam.p["relhum_lev"][..., :-1]) / 2.
  # createProfile defaults this regardless of radar settings, but
  # writePamtraProfile never writes a column for it -- drop it so
  # readPamtraProfile doesn't go looking for a column that isn't there.
  del pam.p["airturb"]
  return pam


@pytest.mark.parametrize("levLay", ["lev", "lay"])
@pytest.mark.parametrize("moment_in", MOMENT_CASES)
def test_profile_roundtrip(tmp_path, moment_in, levLay):
  pam = make_pam(moment_in)
  # layer index 0 for both .lev and .lay (see the module docstring on
  # writePamtraProfile's .lev branch for why level-file layer indexing is
  # zz-1, not zz)
  if moment_in in (1, 12, 13):
    pam.p["hydro_n"][0, 0, 0, 0] = N_VALUE
  if moment_in in (2, 12, 23):
    pam.p["hydro_reff"][0, 0, 0, 0] = REFF_VALUE
  if moment_in in (3, 13, 23):
    pam.p["hydro_q"][0, 0, 0, 0] = Q_VALUE

  profileFile = str(tmp_path / ("profile." + levLay))
  pam.writePamtraProfile(profileFile)

  readBack = make_pam(moment_in)
  readBack.readPamtraProfile(profileFile)

  if moment_in in (1, 12, 13):
    np.testing.assert_allclose(readBack.p["hydro_n"][0, 0, 0, 0], N_VALUE, rtol=1e-6)
  if moment_in in (2, 12, 23):
    np.testing.assert_allclose(readBack.p["hydro_reff"][0, 0, 0, 0], REFF_VALUE, rtol=1e-6)
  if moment_in in (3, 13, 23):
    np.testing.assert_allclose(readBack.p["hydro_q"][0, 0, 0, 0], Q_VALUE, rtol=1e-6)


@pytest.mark.parametrize("levLay", ["lev", "lay"])
def test_profile_roundtrip_last_layer(tmp_path, levLay):
  # Regression test for the off-by-one that indexed hydro_* (0-based, size
  # nlyrs) by the .lev branch's 1-based level loop variable directly:
  # reading/writing the *last* layer is exactly where that ran off the end
  # of the array.
  pam = make_pam(3)
  pam.p["hydro_q"][0, 0, 1, 0] = Q_VALUE  # last of 2 layers (hgt_lev has 3 entries)

  profileFile = str(tmp_path / ("profile." + levLay))
  pam.writePamtraProfile(profileFile)

  readBack = make_pam(3)
  readBack.readPamtraProfile(profileFile)
  np.testing.assert_allclose(readBack.p["hydro_q"][0, 0, 1, 0], Q_VALUE, rtol=1e-6)
