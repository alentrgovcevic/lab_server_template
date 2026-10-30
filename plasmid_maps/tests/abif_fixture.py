"""Write minimal ABIF (.ab1) files for tests: real Sanger files with a peak table are not on hand."""

import struct


def make_ab1(seq: str, qual: list[int], peaks: list[int], channels: dict[str, list[int]], order="GATC") -> bytes:
    """A minimal ABIF file: base calls, qualities, peak table, four processed traces."""
    tags = [
        (b"PBAS", 2, 2, seq.encode()),
        (b"PCON", 2, 2, bytes(qual)),
        (b"PLOC", 2, 4, struct.pack(f">{len(peaks)}h", *peaks)),
        (b"FWO_", 1, 2, order.encode()),
        (b"SMPL", 1, 18, bytes([len(b"read")]) + b"read"),
    ] + [
        (b"DATA", 9 + i, 4, struct.pack(f">{len(channels[b])}h", *channels[b])) for i, b in enumerate(order)
    ]
    data, entries, off = b"", [], 128 + 28 * len(tags)
    for name, num, etype, payload in tags:
        esize = {2: 1, 4: 2, 18: 1}[etype]
        n = len(payload) // esize
        if len(payload) <= 4:
            field = payload.ljust(4, b"\0")
        else:
            field = struct.pack(">i", off + len(data))
            data += payload
        entries.append(name + struct.pack(">ihhii", num, etype, esize, n, len(payload)) + field + b"\0\0\0\0")
    root = b"tdir" + struct.pack(">ihhiii", 1, 1023, 28, len(tags), 28 * len(tags), 128) + b"\0\0\0\0"
    header = (b"ABIF" + struct.pack(">h", 101) + root).ljust(128, b"\0")
    return header + b"".join(entries) + data
