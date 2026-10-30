"""Label images for a lab database or label printer: a QR code or a Code 128 barcode of a plant, seed-bag or
cryo number, as a 1-bit black-on-white PNG. A lab database fetches them (e.g. FileMaker's Insert from URL) into a
container field and prints the number itself, so the barcode carries no human-readable text.

Bars and modules are scaled by a whole number of pixels and never resampled; whatever is left over
of the requested width goes into the quiet zones, so the image is exactly the size asked for.
segno encodes the QR, python-barcode the Code 128 pattern; the PNG writer is stdlib (no Pillow).
"""

from __future__ import annotations

import re
import struct
import zlib

import barcode
import segno

TEXT = re.compile(r"^[A-Za-z0-9._-]{1,40}$")
QR_PX, QR_MIN, QR_MAX = 600, 200, 1200
BAR_PX, BAR_MIN, BAR_MAX = 900, 400, 1800  # ~74 mm at 300 dpi is 874 px
BAR_HEIGHT_RATIO = 200 / 900
QR_QUIET = 4  # modules, the QR minimum
BAR_QUIET = 10  # modules, the Code 128 minimum


def clamp(px: int | None, default: int, lo: int, hi: int) -> int:
    return default if px is None else min(max(px, lo), hi)


def _png(rows: list[list[bool]]) -> bytes:
    """1-bit greyscale PNG; True is black."""
    h, w = len(rows), len(rows[0])
    raw = bytearray()
    for row in rows:
        raw.append(0)  # filter: none
        byte, n = 0, 0
        for black in row:
            byte = (byte << 1) | (0 if black else 1)  # in greyscale 1 is white
            n += 1
            if n == 8:
                raw.append(byte)
                byte, n = 0, 0
        if n:
            raw.append((byte << (8 - n)) | ((1 << (8 - n)) - 1))  # pad with white

    def chunk(kind: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", zlib.crc32(kind + data))
        )

    ihdr = struct.pack(">IIBBBBB", w, h, 1, 0, 0, 0, 0)
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", ihdr)
        + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
        + chunk(b"IEND", b"")
    )


def _place(modules: list[bool], scale: int, width: int) -> list[bool]:
    """One pixel row: modules scaled by `scale`, centred in `width` white pixels."""
    left = (width - len(modules) * scale) // 2
    row = [False] * width
    for i, black in enumerate(modules):
        if black:
            row[left + i * scale : left + (i + 1) * scale] = [True] * scale
    return row


def qr_png(text: str, px: int | None = None) -> bytes:
    """Square QR code of exactly `text`, error correction M, at least a 4-module quiet zone."""
    side = clamp(px, QR_PX, QR_MIN, QR_MAX)
    # make_qr, not make: never a Micro QR for short numbers, and never boosted above level M
    matrix = [
        [bool(v) for v in r]
        for r in segno.make_qr(text, error="m", boost_error=False).matrix
    ]
    n = len(matrix)
    scale = side // (n + 2 * QR_QUIET)
    top = (side - n * scale) // 2
    blank = [False] * side
    rows = [blank] * top
    for r in matrix:
        line = _place(r, scale, side)
        rows.extend([line] * scale)
    rows.extend([blank] * (side - len(rows)))
    return _png(rows)


def code128_png(text: str, px: int | None = None) -> bytes:
    """Code 128 of `text`, bars only (no printed text), at least a 10-module quiet zone each side."""
    width = clamp(px, BAR_PX, BAR_MIN, BAR_MAX)
    pattern = barcode.get_barcode_class("code128")(text).build()[0]
    modules = [c == "1" for c in pattern]
    scale = width // (len(modules) + 2 * BAR_QUIET)
    if (
        scale < 1
    ):  # 40 mixed characters need ~495 modules: widen rather than drop the quiet zones
        scale, width = 1, len(modules) + 2 * BAR_QUIET
    line = _place(modules, scale, width)
    return _png([line] * round(width * BAR_HEIGHT_RATIO))
