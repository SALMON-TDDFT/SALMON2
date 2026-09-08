#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np


def ints(f, n):
    return np.fromfile(f, dtype="<i4", count=n)


def reals(f, n):
    return np.fromfile(f, dtype="<f8", count=n)


def read_wavefunctions(path):
    with open(path, "rb") as f:
        nfrag, nspin, nstate_frag, nstate_tot = ints(f, 4)
        nmat = ints(f, nspin)
        nbasis = ints(f, nfrag * nspin).reshape((nfrag, nspin), order="F")
        index_basis = ints(f, nstate_frag * nfrag * nspin).reshape(
            (nstate_frag, nfrag, nspin), order="F"
        )
        coef = reals(f, nstate_frag * nstate_tot * nspin).reshape(
            (nstate_frag, nstate_tot, nspin), order="F"
        )
    return nfrag, nspin, nstate_frag, nstate_tot, nmat, nbasis, index_basis, coef


def read_basis_functions(path):
    with open(path, "rb") as f:
        nxyz = ints(f, 3)
        nspin, nstate_frag = ints(f, 2)
        nbasis = ints(f, nspin)
        phi = reals(f, int(np.prod(nxyz)) * nspin * nstate_frag).reshape(
            (nxyz[0], nxyz[1], nxyz[2], nspin, nstate_frag), order="F"
        )
    return nxyz, nspin, nstate_frag, nbasis, phi


def read_rgrid_index(path):
    with open(path, "rb") as f:
        local_num = ints(f, 3)
        total_num = ints(f, 3)
        idx = [ints(f, int(local_num[a])) for a in range(3)]
    return local_num, total_num, idx


def read_eigen(path):
    vals = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))
    return np.array(vals, dtype=float)


def read_atoms(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.replace("'", " ").split()
        if len(parts) >= 4:
            atoms.append((parts[0], np.array([float(parts[1]), float(parts[2]), float(parts[3])])))
    return atoms


def pdelta(delta, length):
    return delta - np.rint(delta / length) * length


def c_sp3_value(x, y, z, center, ih, sigma, al):
    signs = [(1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0)]
    sx, sy, sz = signs[ih]
    dx = pdelta(x - center[0], al[0])
    dy = pdelta(y - center[1], al[1])
    dz = pdelta(z - center[2], al[2])
    r2 = dx * dx + dy * dy + dz * dz
    if r2 > (8.0 * sigma) ** 2:
        return 0.0
    pi = math.pi
    gaussian = math.exp(-0.5 * r2 / (sigma * sigma))
    ns = 1.0 / (pi**0.75 * sigma**1.5)
    npart = math.sqrt(2.0) / (pi**0.75 * sigma**2.5)
    return 0.5 * (
        ns * gaussian
        + sx * npart * dx * gaussian
        + sy * npart * dy * gaussian
        + sz * npart * dz * gaussian
    )


def reciprocal_from_index(g, al):
    return 2.0 * math.pi * np.array([g[0] / al[0], g[1] / al[1], g[2] / al[2]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--atom", required=True)
    ap.add_argument("--al", nargs=3, type=float, required=True)
    ap.add_argument("--num-rgrid", nargs=3, type=int, required=True)
    ap.add_argument("--num-wann", type=int, default=512)
    ap.add_argument("--num-bands", type=int, default=512)
    ap.add_argument("--projection-width", type=float, default=1.0)
    args = ap.parse_args()

    data = Path(args.data)
    total = data / "total"
    fragments = data / "fragments"
    seed = args.seed
    al = np.array(args.al, dtype=float)
    hvol = float(np.prod(al / np.array(args.num_rgrid, dtype=float)))
    atoms = read_atoms(args.atom)
    c_atoms = [pos for sym, pos in atoms if sym == "C"]
    if len(c_atoms) * 4 > args.num_wann:
        raise SystemExit("C:sp3 projection count exceeds num_wann")

    eig = read_eigen(total / f"{seed}_eigen.data")
    nband = args.num_bands
    nwann = args.num_wann

    first_wf = fragments / "000001" / "wavefunctions.bin"
    nfrag, nspin, nstate_frag, nstate_tot, nmat, nbasis_all, index_basis, _ = read_wavefunctions(first_wf)
    if nspin != 1:
        raise SystemExit("spin-polarized mode is not supported")

    # Full-state matrices are small enough for this diagnostic path.
    amn = np.zeros((nband, nwann), dtype=float)
    mmn = np.zeros((6, nband, nband), dtype=np.complex128)
    gvecs = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0), (-1, 0, 0)]

    nproj_csp3 = len(c_atoms) * 4
    for ifrag in range(1, nfrag + 1):
        fdir = fragments / f"{ifrag:06d}"
        _, _, _, _, _, _, _, coef_all = read_wavefunctions(fdir / "wavefunctions.bin")
        nxyz, _, _, nbasis_frag, phi = read_basis_functions(fdir / "basis_functions.bin")
        _, _, idx = read_rgrid_index(fdir / "rgrid_index.bin")
        nbasis = int(nbasis_frag[0])
        coef = coef_all[:nbasis, :nband, 0]
        phi_core = phi[:, :, :, 0, :nbasis]

        # AMN: C:sp3 projections plus cyclic shells if num_wann > 4*nC.
        chunk_proj = np.zeros((nbasis, nwann), dtype=float)
        for ip in range(nwann):
            ip_base = ip % nproj_csp3
            shell = ip // nproj_csp3
            atom = c_atoms[ip_base // 4]
            ih = ip_base % 4
            sigma = args.projection_width / math.sqrt(1.0 + 0.75 * shell)
            norm = 0.0
            for ix, gx_i in enumerate(idx[0]):
                x = (int(gx_i) - 1) * al[0] / args.num_rgrid[0]
                for iy, gy_i in enumerate(idx[1]):
                    y = (int(gy_i) - 1) * al[1] / args.num_rgrid[1]
                    for iz, gz_i in enumerate(idx[2]):
                        z = (int(gz_i) - 1) * al[2] / args.num_rgrid[2]
                        gv = c_sp3_value(x, y, z, atom, ih, sigma, al)
                        if gv == 0.0:
                            continue
                        norm += gv * gv * hvol
                        chunk_proj[:, ip] += phi_core[ix, iy, iz, :] * gv * hvol
            if norm > 1.0e-300:
                chunk_proj[:, ip] /= math.sqrt(norm)
        amn += coef.T @ chunk_proj

        # MMN: <psi_m|exp(-i G.r)|psi_n> for Gamma neighbor shells.
        for ig, g in enumerate(gvecs):
            kvec = reciprocal_from_index(g, al)
            sre = np.zeros((nbasis, nbasis), dtype=float)
            sim = np.zeros((nbasis, nbasis), dtype=float)
            for ix, gx_i in enumerate(idx[0]):
                x = (int(gx_i) - 1) * al[0] / args.num_rgrid[0]
                for iy, gy_i in enumerate(idx[1]):
                    y = (int(gy_i) - 1) * al[1] / args.num_rgrid[1]
                    for iz, gz_i in enumerate(idx[2]):
                        z = (int(gz_i) - 1) * al[2] / args.num_rgrid[2]
                        phase = float(kvec @ np.array([x, y, z]))
                        c = math.cos(phase)
                        s = -math.sin(phase)
                        v = phi_core[ix, iy, iz, :]
                        outer = np.outer(v, v) * hvol
                        sre += outer * c
                        sim += outer * s
            tmp = coef.T @ (sre + 1j * sim) @ coef
            mmn[ig, :, :] += tmp

    win = total / f"{seed}.win"
    with open(win, "w") as f:
        f.write(f"num_bands = {nband}\n")
        f.write(f"num_wann = {nwann}\n")
        f.write("num_iter = 100\nmp_grid = 1 1 1\ngamma_only = false\nwrite_hr = true\nwrite_rmn = true\n\n")
        f.write("begin unit_cell_cart\nbohr\n")
        f.write(f"{al[0]:23.15e} 0.0 0.0\n0.0 {al[1]:23.15e} 0.0\n0.0 0.0 {al[2]:23.15e}\n")
        f.write("end unit_cell_cart\n\nbegin atoms_cart\nbohr\n")
        for sym, pos in atoms:
            f.write(f"{sym} {pos[0]:23.15e} {pos[1]:23.15e} {pos[2]:23.15e}\n")
        f.write("end atoms_cart\n\nbegin kpoints\n0.0 0.0 0.0\nend kpoints\n\n")
        f.write("begin projections\n")
        if nproj_csp3 < nwann:
            f.write("random\n")
        f.write("C:sp3\nend projections\n")

    with open(total / f"{seed}.eig", "w") as f:
        for i in range(nband):
            f.write(f"{i+1:8d}{1:8d} {eig[i] * 27.211386245988:23.15e}\n")

    with open(total / f"{seed}.amn", "w") as f:
        f.write("SALMON diagnostic LCFO generated projections\n")
        f.write(f"{nband:10d}{1:10d}{nwann:10d}\n")
        for ip in range(nwann):
            for ib in range(nband):
                f.write(f"{ib+1:8d}{ip+1:8d}{1:8d} {amn[ib, ip]:23.15e} {0.0:23.15e}\n")

    with open(total / f"{seed}.mmn", "w") as f:
        f.write("SALMON diagnostic LCFO generated overlaps\n")
        f.write(f"{nband:10d}{1:10d}{6:10d}\n")
        for ig, g in enumerate(gvecs):
            f.write(f"{1:8d}{1:8d}{g[0]:8d}{g[1]:8d}{g[2]:8d}\n")
            for j in range(nband):
                for i in range(nband):
                    z = mmn[ig, i, j]
                    f.write(f" {z.real:23.15e} {z.imag:23.15e}\n")

    print(f"wrote Wannier90 seed files for {seed} in {total}")


if __name__ == "__main__":
    main()
