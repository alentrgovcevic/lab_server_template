"""Local, content-addressed copies of the readable maps: mirror/<kind>/<id>__<sha8><ext>.

- Only text formats (GenBank, FASTA, bare sequence). `.clc` binaries are never copied.
- The name carries the first 8 hex digits of the raw file's sha256, so an edited file lands in a
  NEW mirror file and every older version stays. Nothing here ever overwrites or deletes.
- The extension follows the sniffed format, not the source name (a Vector NTI `.seq` that is
  really GenBank is mirrored as `.gb`).
- Writes go to a temp name and are renamed into place, so a killed run never leaves a half file.
"""

from __future__ import annotations

import os
from pathlib import Path

from .index import DEFAULT_HOME

DEFAULT_MIRROR = DEFAULT_HOME / "mirror"
MIRRORED_FORMATS = {"genbank": ".gb", "fasta": ".fasta", "raw": ".seq"}


def mirror_name(kind: str, plasmid_id: str | None, sha256: str, fmt: str) -> Path:
    kind_dir = kind if kind in ("PMC", "EC") else "UNIDENTIFIED"
    return (
        Path(kind_dir)
        / f"{plasmid_id or 'UNIDENTIFIED'}__{sha256[:8]}{MIRRORED_FORMATS[fmt]}"
    )


def mirror_bytes(
    data: bytes,
    *,
    kind: str,
    plasmid_id: str | None,
    sha256: str,
    fmt: str,
    mirror_dir: str | Path,
) -> tuple[Path | None, bool]:
    """Write `data` into the mirror unless that version is already there.

    Returns (mirror path, newly written); (None, False) for formats that are not mirrored.
    """
    if fmt not in MIRRORED_FORMATS:
        return None, False
    dest = Path(mirror_dir).expanduser() / mirror_name(kind, plasmid_id, sha256, fmt)
    if dest.exists() and dest.stat().st_size == len(data):
        return dest, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.{os.getpid()}.tmp")
    tmp.write_bytes(data)
    os.replace(tmp, dest)
    return dest, True
