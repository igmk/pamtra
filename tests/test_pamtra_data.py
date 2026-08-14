"""Network-free test of pamtra_data.fetch_data(): serves a small fixture
archive (same shape as the real one -- no single wrapping directory) over
a local HTTP server instead of hitting the real ~250 MB GitHub release download.
"""

import functools
import hashlib
import http.server
import os
import sys
import tarfile
import threading

import pooch
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import pamtra_data  # noqa: E402


@pytest.fixture
def fixture_archive(tmp_path):
    src = tmp_path / "src"
    (src / "dirA").mkdir(parents=True)
    (src / "dirB").mkdir(parents=True)
    (src / "dirA" / "a.txt").write_text("hello a")
    (src / "dirB" / "b.txt").write_text("hello b")

    serve_dir = tmp_path / "serve"
    serve_dir.mkdir()
    archive = serve_dir / "fixture.tar.bz2"
    with tarfile.open(archive, "w:bz2") as tf:
        tf.add(src / "dirA", arcname="dirA")
        tf.add(src / "dirB", arcname="dirB")

    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    return serve_dir, sha256


@pytest.fixture
def http_server(fixture_archive):
    serve_dir, sha256 = fixture_archive
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(serve_dir))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        yield f"http://127.0.0.1:{port}/fixture.tar.bz2", sha256
    finally:
        server.shutdown()
        thread.join(timeout=5)


def test_fetch_data(monkeypatch, tmp_path, http_server):
    url, sha256 = http_server

    monkeypatch.setattr(pamtra_data, "DATA_URL", url)
    monkeypatch.setattr(pamtra_data, "DATA_HASH", f"sha256:{sha256}")
    monkeypatch.setattr(pooch, "os_cache", lambda name: str(tmp_path / "cache"))
    monkeypatch.delenv("PAMTRA_DATADIR", raising=False)

    datadir = pamtra_data.fetch_data(progressbar=False)

    assert datadir == str(tmp_path / "cache" / pamtra_data.EXTRACT_DIR)
    assert os.environ["PAMTRA_DATADIR"] == datadir
    assert sorted(os.listdir(datadir)) == ["dirA", "dirB"]
    assert open(os.path.join(datadir, "dirA", "a.txt")).read() == "hello a"


def test_fetch_data_rejects_corrupted_cache(monkeypatch, tmp_path, http_server):
    # A hash mismatch (corrupted/tampered download) must raise, not
    # silently hand back a broken PAMTRA_DATADIR.
    url, _real_sha256 = http_server

    monkeypatch.setattr(pamtra_data, "DATA_URL", url)
    monkeypatch.setattr(pamtra_data, "DATA_HASH", "sha256:" + "0" * 64)
    monkeypatch.setattr(pooch, "os_cache", lambda name: str(tmp_path / "cache"))
    monkeypatch.delenv("PAMTRA_DATADIR", raising=False)

    with pytest.raises(Exception):
        pamtra_data.fetch_data(progressbar=False)
