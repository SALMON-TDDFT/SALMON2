#!/usr/bin/env python3
"""Structural and generator-reference checks; NOT a scattering/ghost certification."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def validate(path, log):
    lines = path.read_text().splitlines()
    n, h, lmax, z = map(float, lines[0].split())
    n, lmax = int(n), int(lmax)
    assert lines[1] == 'KY_SOC_NLCC_V1'
    assert lines[2] == 'ANGSTROM EV U_R_INV_SQRT_ANGSTROM RHO_INV_ANGSTROM3'
    ipc, lloc, nc = map(int, lines[3].split())
    assert (n, lmax, z, nc) == (2000, 2, 16, 5)
    assert lines[-1] == 'END_KY_SOC_NLCC_V1'
    channels = [list(map(float, x.split())) for x in lines[4:4+nc]]
    assert [(int(x[0]), int(x[1])) for x in channels] == [(0,1),(1,1),(1,3),(2,3),(2,5)]
    assert abs(sum(x[3] for x in channels)-z) < 1e-12
    rows = [list(map(float,x.split())) for x in lines[4+nc:-1]]
    assert len(rows) == n+1
    for i, row in enumerate(rows):
        assert len(row) == 5+2*nc
        assert all(math.isfinite(x) for x in row)
        assert abs(row[0]-i*h) < 1e-10
        assert row[2] >= 0
    assert rows[0][3] == 0
    residual = 0.
    for i in range(2,n-1):
        rho = lambda j: rows[j][2]
        d = (rho(i-2)-8*rho(i-1)+8*rho(i+1)-rho(i+2))/(12*h)
        dd = (-rho(i-2)+16*rho(i-1)-30*rho(i)+16*rho(i+1)-rho(i+2))/(12*h*h)
        residual = max(residual,abs(d-rows[i][3]),abs(dd-rows[i][4]))
    assert residual < 1e-9
    text = log.read_text()
    roots = [list(map(float,line.split()[1:])) for line in text.splitlines() if line.strip().startswith('C2_ROOT')]
    assert len(roots) == nc and 'C2_RANGE_REJECT' not in text
    assert max(abs(x[2]) for x in roots) < 1e-10
    assert 'CRPS_SOC_NLCC_WRITTEN' in text
    # These are the AE reference eigenvalues carried in the generator header.
    splits = {f'l{l}': channels[b][4]-channels[a][4] for l,a,b in [(1,1,2),(2,3,4)]}
    assert all(abs(x)>1e-5 for x in splits.values())
    assert max(abs(row[7]-row[9]) for row in rows) > 1e-5
    assert max(abs(row[11]-row[13]) for row in rows) > 1e-5
    return dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), valence=z,
                core_radius_angstrom=ipc*h, local_channel=lloc,
                reference_j_splitting_eV=splits, nlcc_derivative_residual=residual,
                maximum_norm_constraint_residual=max(abs(x[2]) for x in roots),
                scattering_validated=False, separable_projector_ghost_validated=False,
                solid_or_time_evolution_validated=False)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pseudo',type=Path)
    parser.add_argument('log',type=Path)
    args=parser.parse_args()
    print(json.dumps(validate(args.pseudo,args.log),indent=2))
