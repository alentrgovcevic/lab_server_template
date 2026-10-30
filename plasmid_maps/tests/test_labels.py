"""Label images for a lab database or label printer: /api/qr/{text}.png and /api/code128/{text}.png, and the
bearer token that lets a lab database fetch them without a sign-in."""

import struct
import urllib.error
import urllib.parse
import urllib.request
import zlib

import numpy as np
import pytest
import zxingcpp  # pinned in requirements.txt: an independent decoder, so a missing one fails loudly

from plasmidbank import labels
from tests.test_security import _login, local_app  # noqa: F401  (local_app is a fixture)


TOKEN = "t" * 43


def _pixels(png: bytes) -> np.ndarray:
    """Our 1-bit greyscale PNG back to a uint8 array (0 black, 255 white), with stdlib + numpy."""
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    w, h, depth, colour = struct.unpack(">IIBB", png[16:26])
    assert (depth, colour) == (1, 0)
    i, idat = 8, b""
    while i < len(png):
        n = struct.unpack(">I", png[i : i + 4])[0]
        if png[i + 4 : i + 8] == b"IDAT":
            idat += png[i + 8 : i + 8 + n]
        i += 12 + n
    raw = zlib.decompress(idat)
    stride = (w + 7) // 8 + 1
    rows = [
        np.unpackbits(
            np.frombuffer(raw[r * stride + 1 : (r + 1) * stride], dtype=np.uint8)
        )[:w]
        for r in range(h)
    ]
    return (np.array(rows) * 255).astype(np.uint8)


def _decode(png: bytes):
    found = zxingcpp.read_barcodes(_pixels(png))
    assert len(found) == 1
    return found[0]


CASES = [
    ("qr", "010079198"),
    ("qr", "PMC-01234"),
    ("qr", "a.b_c-D"),
    ("qr", "Z" * 40),
    ("code128", "010079198"),
    ("code128", "020067850"),
    ("code128", "PMS-00567"),
    ("code128", "x.y_Z"),
]
MAKE = {"qr": labels.qr_png, "code128": labels.code128_png}
FORMAT = {
    "qr": zxingcpp.BarcodeFormat.QRCode,
    "code128": zxingcpp.BarcodeFormat.Code128,
}


@pytest.mark.parametrize("kind,text", CASES)
def test_label_decodes_back_to_the_text(kind, text):
    got = _decode(MAKE[kind](text))
    assert (got.text, got.format) == (text, FORMAT[kind])
    if kind == "qr":
        assert got.ec_level == "M"


@pytest.mark.parametrize(
    "kind,px,size",
    [
        ("qr", None, (600, 600)),
        ("qr", 50, (200, 200)),
        ("qr", 333, (333, 333)),
        ("qr", 5000, (1200, 1200)),
        ("code128", None, (900, 200)),
        ("code128", 100, (400, 89)),
        ("code128", 9999, (1800, 400)),
    ],
)
def test_size_is_clamped_and_still_decodes(kind, px, size):
    png = MAKE[kind]("010079198", px)
    a = _pixels(png)
    assert (a.shape[1], a.shape[0]) == size
    assert _decode(png).text == "010079198"


# the standards' minimum quiet zones, in modules; and how many modules the first dark run spans
# (a QR finder pattern's top edge is 7, Code 128's start code begins with a 2-module bar)
@pytest.mark.parametrize("kind,quiet,first_run", [("qr", 4, 7), ("code128", 10, 2)])
def test_quiet_zone_is_white(kind, quiet, first_run):
    a = _pixels(MAKE[kind]("PMC-01234"))
    dark_cols = np.where(a.min(axis=0) == 0)[0]
    dark_rows = np.where(a.min(axis=1) == 0)[0]
    row = a[dark_rows[0]]
    run = np.argmax(row[dark_cols[0]:] == 255)
    module = run / first_run
    assert module >= 1 and module == int(module)  # whole pixels per module, never resampled
    assert dark_cols[0] >= quiet * module
    assert a.shape[1] - 1 - dark_cols[-1] >= quiet * module
    if kind == "qr":
        assert dark_rows[0] >= quiet * module and a.shape[0] - 1 - dark_rows[-1] >= quiet * module


# ── through the server: token, session, validation ─────────────────────────────


def _get(port, path, headers=None):
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", headers=headers or {}
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.headers.get("content-type"), r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers.get("content-type"), e.read()


BEARER = {"Authorization": f"Bearer {TOKEN}"}


@pytest.mark.parametrize("kind,text", [("qr", "010079198"), ("code128", "PMC-01234")])
def test_token_fetches_a_png(local_app, monkeypatch, kind, text):
    monkeypatch.setenv("PLASMIDBANK_LABEL_TOKEN", TOKEN)
    port = local_app[3]
    status, ctype, body = _get(port, f"/api/{kind}/{text}.png", BEARER)
    assert (status, ctype) == (200, "image/png")
    assert _decode(body).text == text


@pytest.mark.parametrize("text", ["a b", "a*b", "x" * 41, "æble", "a;b"])
def test_bad_text_is_400(local_app, monkeypatch, text):
    monkeypatch.setenv("PLASMIDBANK_LABEL_TOKEN", TOKEN)
    port = local_app[3]
    for kind in ("qr", "code128"):
        status, _, _ = _get(port, f"/api/{kind}/{urllib.parse.quote(text)}.png", BEARER)
        assert status == 400, (kind, text)


@pytest.mark.parametrize(
    "env,headers,status",
    [
        (TOKEN, {}, 401),  # no token and no session: the normal sign-in gate
        (TOKEN, {"Authorization": "Bearer wrong"}, 401),
        (TOKEN, {"Authorization": f"Bearer {TOKEN}x"}, 401),
        (TOKEN, {"Authorization": f"Basic {TOKEN}"}, 401),
        (None, BEARER, 503),  # fail closed: no token configured on the server
        ("", BEARER, 503),
        ("", {"Authorization": "Bearer "}, 503),
    ],
)
def test_refused_without_the_right_token(local_app, monkeypatch, env, headers, status):
    if env is None:
        monkeypatch.delenv("PLASMIDBANK_LABEL_TOKEN", raising=False)
    else:
        monkeypatch.setenv("PLASMIDBANK_LABEL_TOKEN", env)
    port = local_app[3]
    for kind in ("qr", "code128"):
        got, _, body = _get(port, f"/api/{kind}/010079198.png", headers)
        assert got == status, (kind, body)
        assert not body.startswith(b"\x89PNG")


@pytest.mark.parametrize(
    "path",
    ["/api/maps", "/api/me", "/api/qr/010079198.svg", "/api/qr/010079198", "/api/stats"],
)
def test_token_opens_nothing_else(local_app, monkeypatch, path):
    monkeypatch.setenv("PLASMIDBANK_LABEL_TOKEN", TOKEN)
    status, _, _ = _get(local_app[3], path, BEARER)
    assert status == 401, path


def test_signed_in_browser_session_needs_no_token(local_app, monkeypatch):
    monkeypatch.delenv("PLASMIDBANK_LABEL_TOKEN", raising=False)
    s, a, c, port = local_app
    a.set_password("maryna.label", "a long enough password")
    _, cookie = _login(c, "maryna.label", "a long enough password")
    status, ctype, body = _get(port, "/api/qr/020067850.png", {"Cookie": cookie})
    assert (status, ctype) == (200, "image/png")
    assert _decode(body).text == "020067850"
