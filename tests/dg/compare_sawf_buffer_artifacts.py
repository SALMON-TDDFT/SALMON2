#!/usr/bin/env python3
import argparse
import struct
from pathlib import Path

import numpy as np


def take(payload, offset, dtype, count):
    values = np.frombuffer(payload, dtype=dtype, count=count, offset=offset)
    return values, offset + values.nbytes


def read_bpw(path):
    payload = path.read_bytes()
    offset = 0
    ints, offset = take(payload, offset, "<i4", 2)
    if tuple(ints) != (-22022215, 3):
        raise ValueError(f"invalid BPW header: {path}")
    header, offset = take(payload, offset, "<i4", 10)
    fragment = int(header[0])
    dims, offset = take(payload, offset, "<i4", 4)
    _nspin, nbasis, nproj, nkeep = map(int, dims)
    _, offset = take(payload, offset, "<i4", nproj)
    _, offset = take(payload, offset, "<i4", nproj)
    lambdas, offset = take(payload, offset, "<f8", nkeep)
    _, offset = take(payload, offset, "<i4", nkeep)
    spreads, offset = take(payload, offset, "<f8", nkeep)
    _, offset = take(payload, offset, "<f8", nkeep)
    _, offset = take(payload, offset, "<f8", nbasis * nkeep)
    r, offset = take(payload, offset, "<f8", 3 * nkeep * nkeep)
    r = r.reshape((3, nkeep, nkeep), order="F")
    centers, offset = take(payload, offset, "<f8", 3 * nkeep)
    centers = centers.reshape((3, nkeep), order="F")
    h, offset = take(payload, offset, "<f8", nkeep * nkeep)
    h = h.reshape((nkeep, nkeep), order="F")
    v, offset = take(payload, offset, "<f8", 3 * nkeep * nkeep)
    v = v.reshape((3, nkeep, nkeep), order="F")
    aa, offset = take(payload, offset, "<f8", 3 * nkeep * nkeep)
    aa = aa.reshape((3, nkeep, nkeep), order="F")
    _, offset = take(payload, offset, "<i4", 1)
    if offset != len(payload):
        raise ValueError(f"unexpected BPW payload length: {path}")
    return {"fragment": fragment, "nkeep": nkeep, "lambda": lambdas,
            "spread": spreads, "r": r, "center": centers, "h": h, "v": v, "aa": aa}


def read_trace(path):
    payload = path.read_bytes()
    offset = 0
    header, offset = take(payload, offset, "<i4", 12)
    if tuple(header[:2]) != (-22022218, 1):
        raise ValueError(f"invalid trace header: {path}")
    fragment = int(header[2])
    domain = tuple(map(int, header[3:6]))
    metrics, offset = take(payload, offset, "<f8", 4)
    nkeep_values, offset = take(payload, offset, "<i4", 1)
    nkeep = int(nkeep_values[0])
    faces = []
    features = []
    for _ in range(6):
        face_header, offset = take(payload, offset, "<i4", 3)
        axis, side, npoints = map(int, face_header)
        weights, offset = take(payload, offset, "<f8", 2)
        area, alpha = map(float, weights)
        u, offset = take(payload, offset, "<f8", npoints * nkeep)
        dn, offset = take(payload, offset, "<f8", npoints * nkeep)
        u = u.reshape((npoints, nkeep), order="F")
        dn = dn.reshape((npoints, nkeep), order="F")
        faces.append((axis, side, area, alpha, u, dn))
        features.extend((np.sqrt(area) * u, np.sqrt(area) * dn / alpha))
    if offset != len(payload):
        raise ValueError(f"unexpected trace payload length: {path}")
    return {"fragment": fragment, "domain": domain, "hgs": metrics[:3],
            "nkeep": nkeep, "faces": faces, "features": np.vstack(features)}


def rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), np.finfo(float).tiny)


parser = argparse.ArgumentParser()
parser.add_argument("smaller", type=Path)
parser.add_argument("larger", type=Path)
args = parser.parse_args()

maxima = {name: 0.0 for name in ("gauge", "overlap", "ww", "position", "face")}
for fragment in range(1, 9):
    name = f"{fragment:06d}"
    a = read_bpw(args.smaller / "fragments" / name / "buffer_periodic_wannier_basis.bin")
    b = read_bpw(args.larger / "fragments" / name / "buffer_periodic_wannier_basis.bin")
    ta = read_trace(args.smaller / "fragments" / name / "buffer_periodic_wannier_trace.bin")
    tb = read_trace(args.larger / "fragments" / name / "buffer_periodic_wannier_trace.bin")
    if a["nkeep"] != b["nkeep"] or ta["domain"] != tb["domain"]:
        raise SystemExit(f"incompatible buffer artifacts for fragment {fragment}")
    cross = ta["features"].T @ tb["features"]
    left, _, right_h = np.linalg.svd(cross, full_matrices=False)
    gauge = left @ right_h
    maxima["gauge"] = max(maxima["gauge"], rel(ta["features"] @ gauge, tb["features"]))
    maxima["overlap"] = max(maxima["overlap"], rel(gauge.T @ gauge, np.eye(a["nkeep"])))
    maxima["ww"] = max(maxima["ww"], rel(gauge.T @ a["h"] @ gauge, b["h"]))
    for axis in range(3):
        maxima["position"] = max(maxima["position"], rel(gauge.T @ a["aa"][axis] @ gauge, b["aa"][axis]))
    for fa, fb in zip(ta["faces"], tb["faces"]):
        if fa[:2] != fb[:2]:
            raise SystemExit(f"face ordering mismatch for fragment {fragment}")
        _, _, area_a, alpha_a, ua, dna = fa
        _, _, area_b, alpha_b, ub, dnb = fb
        block_a = area_a * (-0.5 * (dna.T @ ua + ua.T @ dna) + alpha_a * ua.T @ ua)
        block_b = area_b * (-0.5 * (dnb.T @ ub + ub.T @ dnb) + alpha_b * ub.T @ ub)
        maxima["face"] = max(maxima["face"], rel(gauge.T @ block_a @ gauge, block_b))

print(" ".join(f"{key}={value:.8e}" for key, value in maxima.items()))
