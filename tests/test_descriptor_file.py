import tempfile
from pathlib import Path

import pyPamtra

REPO_ROOT = Path(__file__).parent.parent
DESCRIPTOR_FILE = REPO_ROOT / "descriptorfiles" / "descriptor_file_2m_icon.txt"


def test_read_file():
    pam = pyPamtra.pyPamtra()
    pam.df.readFile(str(DESCRIPTOR_FILE))
    assert pam.df.nhydro == 6
    assert list(pam.df.data["hydro_name"]) == [
        "cwc_q", "iwc_q", "rwc_q", "swc_q", "gwc_q", "hwc_q",
    ]


def test_write_read_roundtrip(tmp_path):
    pam = pyPamtra.pyPamtra()
    pam.df.readFile(str(DESCRIPTOR_FILE))

    out_file = tmp_path / "roundtrip.txt"
    pam.df.writeFile(str(out_file))

    pam2 = pyPamtra.pyPamtra()
    pam2.df.readFile(str(out_file))

    assert (pam.df.data == pam2.df.data).all()


def test_add_hydrometeor():
    pam = pyPamtra.pyPamtra()
    pam.df.addHydrometeor((
        "rwc_q", -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0,
        3, 30, "exp", -99.0, -99.0, 8.0e6, -99.0, 1.0e-5, 6.0e-3,
        "mie-sphere", "khvorostyanov01_drops", -99.0,
    ))
    assert pam.df.nhydro == 1
    assert pam.df.data["hydro_name"][0] == "rwc_q"
