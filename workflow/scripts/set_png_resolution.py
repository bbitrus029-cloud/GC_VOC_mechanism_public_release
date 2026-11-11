#!/usr/bin/env python3
"""Add PNG physical-resolution metadata without decoding or changing pixels."""
from pathlib import Path
import hashlib
import struct
import sys
import zlib


def set_resolution(path, dpi):
    path = Path(path)
    raw = path.read_bytes()
    assert raw[:8] == b'\x89PNG\r\n\x1a\n'
    chunks = []
    pos = 8
    while pos < len(raw):
        size = struct.unpack('>I', raw[pos:pos + 4])[0]
        chunks.append((raw[pos + 4:pos + 8], raw[pos:pos + size + 12]))
        pos += size + 12
    pixels = b''.join(c for kind, c in chunks if kind == b'IDAT')
    ppm = round(float(dpi) / .0254)
    data = b'pHYs' + struct.pack('>IIB', ppm, ppm, 1)
    physical = struct.pack('>I', 9) + data + struct.pack('>I', zlib.crc32(data))
    output = raw[:8]
    for kind, chunk in chunks:
        if kind == b'pHYs':
            continue
        output += chunk
        if kind == b'IHDR':
            output += physical
    check = []
    pos = 8
    while pos < len(output):
        size = struct.unpack('>I', output[pos:pos + 4])[0]
        if output[pos + 4:pos + 8] == b'IDAT':
            check.append(output[pos:pos + size + 12])
        pos += size + 12
    assert hashlib.sha256(pixels).digest() == hashlib.sha256(b''.join(check)).digest()
    path.write_bytes(output)


if __name__ == '__main__':
    set_resolution(sys.argv[1], float(sys.argv[2]))
