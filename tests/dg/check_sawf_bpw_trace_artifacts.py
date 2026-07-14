#!/usr/bin/env python3
import argparse
import math
import struct
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("run_directory", type=Path)
args = parser.parse_args()

fragments = args.run_directory / "data_dcdft" / "fragments"
if not fragments.is_dir():
    fragments = args.run_directory / "fragments"
trace_files = sorted(fragments.glob("*/buffer_periodic_wannier_trace.bin"))
basis_files = sorted(fragments.glob("*/buffer_periodic_wannier_basis.bin"))
if len(trace_files) != 8 or len(basis_files) != 8:
    raise SystemExit(f"expected 8 BPW and trace files, got {len(basis_files)} and {len(trace_files)}")

expected_magic = -22022218
for path in trace_files:
    payload = memoryview(path.read_bytes())
    offset = 0

    def unpack(fmt):
        nonlocal_offset = offset
        size = struct.calcsize(fmt)
        if nonlocal_offset + size > len(payload):
            raise SystemExit(f"truncated trace payload: {path}")
        values = struct.unpack_from(fmt, payload, nonlocal_offset)
        return values, nonlocal_offset + size

    (magic, version), offset = unpack("<2i")
    header, offset = unpack("<10i")
    fragment = header[0]
    domain, buffer, box = header[1:4], header[4:7], header[7:10]
    metrics, offset = unpack("<4d")
    (nkeep,), offset = unpack("<i")
    if magic != expected_magic or version != 1 or nkeep != 32:
        raise SystemExit(f"invalid trace header: {path}")
    if box != tuple(domain[i] + 2 * buffer[i] for i in range(3)):
        raise SystemExit(f"inconsistent buffer box: {path}")
    if not all(math.isfinite(x) and x > 0.0 for x in metrics):
        raise SystemExit(f"invalid trace metric: {path}")

    seen = set()
    max_abs = 0.0
    for _ in range(6):
        (axis, side, npoints), offset = unpack("<3i")
        (area_weight, alpha), offset = unpack("<2d")
        if axis not in (1, 2, 3) or side not in (-1, 1) or (axis, side) in seen:
            raise SystemExit(f"invalid or duplicate face: {path}")
        seen.add((axis, side))
        expected_points = math.prod(domain[i] for i in range(3) if i != axis - 1)
        if npoints != expected_points or area_weight <= 0.0 or alpha <= 0.0:
            raise SystemExit(f"invalid face metadata: {path}")
        count = npoints * nkeep
        values, offset = unpack(f"<{2 * count}d")
        if not all(math.isfinite(x) for x in values):
            raise SystemExit(f"non-finite face data: {path}")
        max_abs = max(max_abs, max(abs(x) for x in values))
    if offset != len(payload) or max_abs == 0.0 or fragment < 1:
        raise SystemExit(f"invalid trace length or zero trace: {path}")

print("PASS 8 fragment BPW artifacts contain finite six-face symmetry-closed traces")
