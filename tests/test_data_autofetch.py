"""pyPamtra/core.py's PAMTRA_DATADIR handling, exercised in subprocesses
so each case gets a genuinely fresh interpreter (import side effects only
run once per process, so this can't be tested by re-importing pyPamtra
within the main test session).
"""

import os
import subprocess
import sys

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")

_SETUP = """
import functools, hashlib, http.server, os, sys, tarfile, tempfile, threading

tmp = tempfile.mkdtemp()
src = os.path.join(tmp, "src")
for d in ("emissivity", "hongdb", "scatdb"):
    os.makedirs(os.path.join(src, d))
with open(os.path.join(src, "emissivity", "x.dat"), "w") as f:
    f.write("x")

serve_dir = os.path.join(tmp, "serve")
os.makedirs(serve_dir)
archive = os.path.join(serve_dir, "pamtra_data.tar.bz2")
with tarfile.open(archive, "w:bz2") as tf:
    for d in ("emissivity", "hongdb", "scatdb"):
        tf.add(os.path.join(src, d), arcname=d)
sha256 = hashlib.sha256(open(archive, "rb").read()).hexdigest()

handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=serve_dir)
server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
threading.Thread(target=server.serve_forever, daemon=True).start()
port = server.server_address[1]

sys.path.insert(0, "python")
import pamtra_data
pamtra_data.DATA_URL = f"http://127.0.0.1:{port}/pamtra_data.tar.bz2"
pamtra_data.DATA_HASH = f"sha256:{sha256}"

import pooch
pooch.os_cache = lambda name: os.path.join(tmp, "cache")

os.environ.pop("PAMTRA_DATADIR", None)
"""


def _run(script, timeout=90):
    return subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=REPO_ROOT,
    )


def test_unset_datadir_autofetches_and_imports():
    result = _run(_SETUP + """
import pyPamtra
assert os.environ["PAMTRA_DATADIR"], "PAMTRA_DATADIR should be set after import"
assert os.path.isdir(os.environ["PAMTRA_DATADIR"])
print("AUTOFETCH_OK", os.environ["PAMTRA_DATADIR"])
""")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "AUTOFETCH_OK" in result.stdout


def test_explicit_empty_datadir_skips_autofetch():
    result = _run("""
import os
os.environ["PAMTRA_DATADIR"] = ""
import pyPamtra
print("SKIPPED_OK")
""")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "SKIPPED_OK" in result.stdout
    # must not have touched pooch/network at all
    assert "Downloading" not in result.stdout


def test_unreachable_url_raises_clear_runtime_error():
    result = _run("""
import os, sys
sys.path.insert(0, "python")
import pamtra_data
pamtra_data.DATA_URL = "http://127.0.0.1:1/nonexistent.tar.bz2"  # nothing listening
os.environ.pop("PAMTRA_DATADIR", None)
import pyPamtra
""")
    assert result.returncode != 0
    assert "RuntimeError" in result.stderr
    assert "PAMTRA_DATADIR not set" in result.stderr
