"""Read-only numerical kernel for matched, frozen-potential density snapshots.

Caller must establish identical physical grid/order, DC seed, frozen potential,
pseudopotentials, electron count, temperature, and MPI/rank-fragment mapping.
This kernel cannot authenticate provenance or certify physical accuracy.
"""
from math import fsum, isfinite, sqrt
from pathlib import Path
import hashlib
import json
import argparse


def compare_densities(dc, conventional, dg, weights, *, relaxed=None):
    """Report weighted L2 differences with one common DC-density denominator.

    conventional and dg must both come from H[rho_DC], before density feedback.
    relaxed is optional and is never substituted for the frozen DG snapshot.
    Differences from conventional LCFO are not errors against an exact solution.
    """
    arrays = [list(values) for values in (dc, conventional, dg, weights)]
    if relaxed is not None:
        arrays.append(list(relaxed))
    n = len(arrays[0])
    if not n or any(len(values) != n for values in arrays):
        raise ValueError("density snapshots must share a nonempty physical grid")
    if any(not isfinite(x) for values in arrays for x in values):
        raise ValueError("nonfinite density or integration weight")
    dc, conventional, dg, weights = arrays[:4]
    if any(w <= 0 for w in weights):
        raise ValueError("integration weights must be positive")
    denominator = fsum(w*x*x for w, x in zip(weights, dc))
    if not isfinite(denominator) or denominator <= 0:
        raise ValueError("DC density must have a finite positive weighted norm")
    baseline = [b-a for a, b in zip(dc, conventional)]
    increment = [b-a for a, b in zip(conventional, dg)]

    def norm(delta):
        value = fsum(w*x*x for w, x in zip(weights, delta)) / denominator
        if not isfinite(value):
            raise ValueError("density difference norm overflow")
        return sqrt(value)

    cross = 2*fsum(w*a*b for w, a, b in zip(weights, baseline, increment))/denominator
    if not isfinite(cross):
        raise ValueError("density cross term overflow")
    return {
        "dc_lcfo_correction": norm(baseline),
        "dg_increment": norm(increment),
        "frozen_total_displacement": norm([b-a for a, b in zip(dc, dg)]),
        "relaxation_displacement": None if relaxed is None else
            norm([b-a for a, b in zip(dg, arrays[4])]),
        "cross_term": cross,
        "absolute_dc_error_available": False,
    }


def read_snapshot_set(prefix):
    """Validate a complete one-rank-per-fragment export, without mutating it.

    SHA256 values identify the analyzed bytes, not an authenticated restart.
    No pass/fail physical accuracy threshold is inferred here.
    """
    paths = sorted(Path(prefix).parent.glob(Path(prefix).name + ".rank-*"))
    if not paths:
        raise ValueError("density snapshots unavailable")
    rows, hashes, reference = {}, {}, None
    for rank, path in enumerate(paths):
        raw = path.read_bytes()
        lines = raw.decode("ascii").splitlines()
        if len(lines) < 6 or lines[0] != "SALMON_DG_DENSITY_DIAGNOSTIC_V1":
            raise ValueError("unknown density snapshot format")
        nproc, stored_rank, fragment, count, grid_count = map(int, lines[1].split())
        provenance = tuple(map(int, lines[2].split()))
        electron_count, temperature = map(float, lines[3].split())
        if (nproc != len(paths) or stored_rank != rank or fragment != rank+1 or
                path.name != Path(prefix).name + f".rank-{rank:08d}" or
                count < 1 or grid_count < count or len(provenance) != 3 or
                any(x == 0 for x in provenance) or not isfinite(electron_count) or
                electron_count <= 0 or not isfinite(temperature) or temperature < 0):
            raise ValueError("invalid snapshot identity or controls")
        identity = (nproc, grid_count, provenance, electron_count, temperature)
        if reference is not None and identity != reference:
            raise ValueError("rank-disagreeing density snapshot provenance")
        reference = identity
        if len(lines) != count+5 or lines[-1] != "END_DENSITY_DIAGNOSTIC":
            raise ValueError("truncated or extended density snapshot")
        for line in lines[4:-1]:
            tokens = line.split()
            if len(tokens) != 7:
                raise ValueError("invalid density snapshot row")
            point = int(tokens[0])
            values = tuple(map(float, tokens[1:]))
            if (point < 1 or point > grid_count or point in rows or
                    any(not isfinite(x) for x in values) or values[0] <= 0 or
                    any(x < 0 for x in values[1:5])):
                raise ValueError("invalid or duplicate physical grid row")
            rows[point] = values
        hashes[path.name] = hashlib.sha256(raw).hexdigest()
    if len(rows) != reference[1]:
        raise ValueError("incomplete physical grid coverage")
    weights, dc, conventional, dg, relaxed, potential = zip(*(rows[i] for i in sorted(rows)))
    return {
        "metrics": compare_densities(dc, conventional, dg, weights, relaxed=relaxed),
        "electron_counts": {name: fsum(w*r for w, r in zip(weights, density))
                            for name, density in (("dc", dc), ("conventional", conventional),
                                                  ("dg_frozen", dg), ("dg_final", relaxed))},
        "target_electrons": reference[3], "temperature_au": reference[4],
        "mpi_ranks": reference[0], "source_fingerprints": reference[2],
        "sha256": hashes,
        "accuracy_certified": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path)
    args = parser.parse_args()
    print(json.dumps(read_snapshot_set(args.prefix), indent=2, allow_nan=False))
