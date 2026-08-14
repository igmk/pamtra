import pyPamtra


def test_write_read_nml_roundtrip(tmp_path):
    pam = pyPamtra.pyPamtra()
    pam.nmlSet["emissivity"] = 0.85
    pam.nmlSet["active"] = False
    pam.nmlSet["gas_mod"] = "L93"

    nml_file = tmp_path / "test.nml"
    pam.writeNmlFile(str(nml_file))

    pam2 = pyPamtra.pyPamtra()
    pam2.readNmlFile(str(nml_file))

    assert pam2.nmlSet["emissivity"] == 0.85
    assert pam2.nmlSet["active"] is False
    assert pam2.nmlSet["gas_mod"] == "L93"
    # a key not touched above should still round-trip with its default
    assert pam2.nmlSet["passive"] is True
