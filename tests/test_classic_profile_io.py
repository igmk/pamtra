"""Tests for readClassicPamtraProfile (the old '.cla'-style classic
profile format, distinct from the newer .lev/.lay format covered by
test_profile_io.py).

There is no example/fixture file for this format anywhere in the repo,
no writer to round-trip against, and (see the docstring on
readClassicPamtraProfile) its header doesn't match what the standalone
CLI's own '.cla' reader expects -- so unlike test_profile_io.py, this
can't cross-check against an independent implementation. These tests
instead pin down readClassicPamtraProfile's own documented column
layout with a hand-written fixture, so a change to it doesn't silently
break without at least one test noticing.
"""

import numpy as np
import pytest

import pyPamtra

# 1 gridpoint, 2 layers. Field-by-field:
#   header:        year month day time ngridx ngridy nlyrs deltax deltay
#   per gridpoint: model_i model_j
#                  lat lon lfrac wind10u wind10v
#                  iwv cwp iwp rwp swp gwp [hwp for n_moments==2]
#                  hgt_lev[0] press_lev[0] temp_lev[0] relhum_lev[0]
#   per layer:     hgt_lev press_lev temp_lev relhum_lev cwc_q iwc_q rwc_q
#                  swc_q gwc_q [hwc_q cwc_n iwc_n rwc_n swc_n gwc_n hwc_n
#                  for n_moments==2]
PROFILE_1MOM = """\
2020 01 01 1200 1 1 2 1000 1000
1 1
50.0 10.0 1 0.0 0.0
10.0 0.0 0.0 0.0 0.0 0.0
0.0 101325.0 288.15 50.0
1000.0 89000.0 281.0 60.0 0.0 0.0 0.001 0.0 0.0
2000.0 79000.0 275.0 70.0 0.0 0.0 0.0 0.0005 0.0
"""

PROFILE_2MOM = """\
2020 01 01 1200 1 1 2 1000 1000
1 1
50.0 10.0 1 0.0 0.0
10.0 0.0 0.0 0.0 0.0 0.0 0.0
0.0 101325.0 288.15 50.0
1000.0 89000.0 281.0 60.0 0.0 0.0 0.001 0.0 0.0 0.0 0.0 0.0 1.5e6 0.0 0.0 0.0
2000.0 79000.0 275.0 70.0 0.0 0.0 0.0 0.0005 0.0 0.0 0.0 0.0 0.0 2.5e5 0.0 0.0
"""

# hydro_q/hydro_n column order is cwc, iwc, rwc, swc, gwc[, hwc]
RWC_INDEX = 2
SWC_INDEX = 3


def test_read_classic_profile_1mom(tmp_path):
  profileFile = tmp_path / "profile.cla"
  profileFile.write_text(PROFILE_1MOM)

  pam = pyPamtra.pyPamtra()
  pam.readClassicPamtraProfile(str(profileFile))

  assert pam.p["ngridx"] == 1
  assert pam.p["ngridy"] == 1
  assert pam.p["nlyrs"][0, 0] == 2
  assert pam.p["sfc_type"][0, 0] == 1
  np.testing.assert_allclose(pam.p["lat"][0, 0], 50.0)
  np.testing.assert_allclose(pam.p["lon"][0, 0], 10.0)
  np.testing.assert_allclose(pam.p["iwv"][0, 0], 10.0)
  np.testing.assert_allclose(pam.p["hgt_lev"][0, 0], [0.0, 1000.0, 2000.0])
  np.testing.assert_allclose(pam.p["press_lev"][0, 0], [101325.0, 89000.0, 79000.0])
  np.testing.assert_allclose(pam.p["temp_lev"][0, 0], [288.15, 281.0, 275.0])
  # relhum_lev is converted from % to a fraction
  np.testing.assert_allclose(pam.p["relhum_lev"][0, 0], [0.5, 0.6, 0.7])
  np.testing.assert_allclose(pam.p["hydro_q"][0, 0, 0, RWC_INDEX], 0.001)
  np.testing.assert_allclose(pam.p["hydro_q"][0, 0, 1, SWC_INDEX], 0.0005)


def test_read_classic_profile_2mom(tmp_path):
  profileFile = tmp_path / "profile.cla"
  profileFile.write_text(PROFILE_2MOM)

  pam = pyPamtra.pyPamtra()
  pam.readClassicPamtraProfile(str(profileFile), n_moments=2)

  assert pam.p["sfc_type"][0, 0] == 1
  np.testing.assert_allclose(pam.p["hydro_q"][0, 0, 0, RWC_INDEX], 0.001)
  np.testing.assert_allclose(pam.p["hydro_n"][0, 0, 0, RWC_INDEX], 1.5e6)
  np.testing.assert_allclose(pam.p["hydro_q"][0, 0, 1, SWC_INDEX], 0.0005)
  np.testing.assert_allclose(pam.p["hydro_n"][0, 0, 1, SWC_INDEX], 2.5e5)


def test_read_classic_profile_invalid_n_moments(tmp_path):
  profileFile = tmp_path / "profile.cla"
  profileFile.write_text(PROFILE_1MOM)

  pam = pyPamtra.pyPamtra()
  with pytest.raises(IOError):
    pam.readClassicPamtraProfile(str(profileFile), n_moments=3)
