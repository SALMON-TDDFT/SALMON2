#!/usr/bin/env python3
"""Export Gamma-only Wannier90 seed files from a conventional SALMON restart WF."""

import argparse
from pathlib import Path

import numpy as np

HARTREE_TO_EV = 27.211386245988


def read_atoms(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        fields = line.replace("'", " ").split()
        if len(fields) >= 4:
            atoms.append((fields[0], np.asarray(fields[1:4], dtype=float)))
    return atoms


def read_energies(path, count):
    values = []
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit():
            values.append(float(fields[1]))
    if len(values) < count:
        raise ValueError(f"eigenvalue file has {len(values)} states; expected {count}")
    return np.asarray(values[:count])


def grid_coordinates(mesh, cell):
    axes = [np.arange(n, dtype=float) * length / n for n, length in zip(mesh, cell)]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    return np.column_stack((x.ravel(order="F"), y.ravel(order="F"), z.ravel(order="F")))


def real_spd_shell(delta, sigma):
    q = delta / sigma
    x, y, z = q.T
    gaussian = np.exp(-0.5 * np.einsum("ij,ij->i", q, q))
    return np.column_stack(
        (
            gaussian,
            z * gaussian,
            x * gaussian,
            y * gaussian,
            (2.0 * z * z - x * x - y * y) / np.sqrt(6.0) * gaussian,
            np.sqrt(2.0) * z * x * gaussian,
            np.sqrt(2.0) * y * z * gaussian,
            (x * x - y * y) / np.sqrt(2.0) * gaussian,
            np.sqrt(2.0) * x * y * gaussian,
        )
    )


def build_amn(psi, coordinates, atoms, cell, sigma, hvol, chunk_atoms, shell):
    nband = psi.shape[1]
    nchan = 4 if shell == "sp" else 9
    nproj = nchan * len(atoms)
    amn = np.empty((nband, nproj), dtype=float)
    for first in range(0, len(atoms), chunk_atoms):
        last = min(len(atoms), first + chunk_atoms)
        projectors = np.empty((len(coordinates), nchan * (last - first)), dtype=float)
        for local, (_, center) in enumerate(atoms[first:last]):
            delta = coordinates - center
            delta -= np.rint(delta / cell) * cell
            projectors[:, nchan * local : nchan * (local + 1)] = real_spd_shell(delta, sigma)[:, :nchan]
        norms = np.sqrt(hvol * np.einsum("ij,ij->j", projectors, projectors))
        if np.any(norms <= 1.0e-14):
            raise ValueError("zero-norm pseudo-channel projector")
        projectors /= norms
        amn[:, nchan * first : nchan * last] = hvol * (psi.T @ projectors)
    return amn


def build_mmn(psi, coordinates, cell, hvol):
    result = []
    for axis in range(3):
        phase = np.exp(-2j * np.pi * coordinates[:, axis] / cell[axis])
        result.append(hvol * (psi.T @ (phase[:, None] * psi)))
    return result


def read_symmetry_operations(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            rows.append([float(value) for value in line.split()])
    if len(rows) % 3:
        raise ValueError("symmetry file must contain three rows per operation")
    operations = []
    for first in range(0, len(rows), 3):
        block = np.asarray(rows[first : first + 3])
        rotation = np.rint(block[:, :3]).astype(int)
        if np.max(np.abs(block[:, :3] - rotation)) > 1.0e-10:
            raise ValueError("symmetry rotation is not integer")
        operations.append((rotation, np.mod(block[:, 3], 1.0)))
    return operations


def atom_permutation(atoms, cell, rotation, translation, tolerance=1.0e-8):
    fractional = np.asarray([position / cell for _, position in atoms])
    symbols = [symbol for symbol, _ in atoms]
    permutation = np.full(len(atoms), -1, dtype=int)
    used = set()
    for source, position in enumerate(fractional):
        transformed = rotation @ position + translation
        best = None
        for target, candidate in enumerate(fractional):
            if target in used or symbols[target] != symbols[source]:
                continue
            delta = transformed - candidate
            delta -= np.rint(delta)
            distance = np.linalg.norm(delta * cell)
            if best is None or distance < best[0]:
                best = (distance, target)
        if best is None or best[0] > tolerance:
            raise ValueError(f"symmetry operation does not map atom {source + 1}")
        permutation[source] = best[1]
        used.add(best[1])
    return permutation


def spd_rotation_representation(rotation):
    p_basis = np.column_stack((np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])))
    p_rep = p_basis.T @ rotation @ p_basis
    # Evaluate the five real harmonics on deterministic sample points and solve
    # for their representation. This follows the same z2,zx,yz,x2-y2,xy order.
    samples = np.asarray(((1, 2, 3), (2, -1, 1), (-2, 1, 3), (1, -3, 2), (3, 1, -2)), dtype=float)
    def d_values(points):
        x, y, z = points.T
        return np.column_stack(((2*z*z-x*x-y*y)/np.sqrt(6), np.sqrt(2)*z*x, np.sqrt(2)*y*z,
                                (x*x-y*y)/np.sqrt(2), np.sqrt(2)*x*y))
    d_rep = np.linalg.solve(d_values(samples), d_values(samples @ rotation))
    block = np.zeros((9, 9))
    block[0, 0] = 1.0
    block[1:4, 1:4] = p_rep
    block[4:9, 4:9] = d_rep
    return block


def build_dmn(psi, atoms, mesh, cell, operations, amn, energies, hvol, shell):
    nband, nwann = amn.shape
    nchan = 4 if shell == "sp" else 9
    indices = np.indices(tuple(mesh), dtype=int).reshape(3, -1, order="F")
    d_wann_all, d_band_all, diagnostics = [], [], []
    for rotation, translation in operations:
        shift = translation * mesh
        if np.max(np.abs(shift - np.rint(shift))) > 1.0e-10:
            raise ValueError("Full-WF DMN currently requires grid-commensurate translations")
        target = rotation @ indices + np.rint(shift).astype(int)[:, None]
        target %= mesh[:, None]
        target_flat = np.ravel_multi_index(tuple(target), tuple(mesh), order="F")
        d_band = hvol * (psi[target_flat, :].T @ psi)
        permutation = atom_permutation(atoms, cell, rotation, translation)
        orbital_rep = spd_rotation_representation(rotation.astype(float))[:nchan, :nchan]
        d_wann = np.zeros((nwann, nwann))
        for source, target_atom in enumerate(permutation):
            d_wann[nchan * target_atom : nchan * (target_atom + 1),
                   nchan * source : nchan * (source + 1)] = orbital_rep
        band_unitarity = np.max(np.abs(d_band.T @ d_band - np.eye(nband)))
        wann_unitarity = np.max(np.abs(d_wann.T @ d_wann - np.eye(nwann)))
        energy_scale = max(np.ptp(energies), 1.0e-30)
        h_cov = np.max(np.abs(d_band.T @ (energies[:, None] * d_band) - np.diag(energies))) / energy_scale
        amn_cov = np.max(np.abs(d_band @ amn - amn @ d_wann)) / max(1.0, np.max(np.abs(amn)))
        diagnostics.append((band_unitarity, wann_unitarity, h_cov, amn_cov))
        d_wann_all.append(d_wann.astype(complex))
        d_band_all.append(d_band.astype(complex))
    return d_wann_all, d_band_all, diagnostics


def write_dmn(path, d_wann, d_band):
    nband, nwann = d_band[0].shape[0], d_wann[0].shape[0]
    with path.open("w") as out:
        out.write("SALMON conventional Full-WF SAWF symmetry data\n")
        out.write(f"{nband:9d}{len(d_band):9d}{1:9d}{1:9d}\n\n{1:9d}\n\n{1:9d}\n\n")
        out.write("".join(f"{1:9d}" for _ in d_band) + "\n")
        for matrix in d_wann:
            out.write("\n")
            for value in matrix.ravel(order="F"):
                out.write(f" ({value.real:18.10E},{value.imag:18.10E})")
            out.write("\n")
        for matrix in d_band:
            out.write("\n")
            for value in matrix.ravel(order="F"):
                out.write(f" ({value.real:18.10E},{value.imag:18.10E})")
            out.write("\n")


def write_win(path, nband, nwann, cell, atoms, num_iter, shell):
    with path.open("w") as out:
        out.write(f"num_bands = {nband}\nnum_wann = {nwann}\nnum_iter = {num_iter}\n")
        out.write("mp_grid = 1 1 1\ngamma_only = true\nwrite_hr = true\nwrite_rmn = true\n\n")
        out.write("begin unit_cell_cart\nbohr\n")
        out.write(f"{cell[0]:23.15e} 0.0 0.0\n0.0 {cell[1]:23.15e} 0.0\n0.0 0.0 {cell[2]:23.15e}\n")
        out.write("end unit_cell_cart\n\nbegin atoms_cart\nbohr\n")
        for symbol, position in atoms:
            out.write(f"{symbol} {position[0]:23.15e} {position[1]:23.15e} {position[2]:23.15e}\n")
        out.write("end atoms_cart\n\nbegin kpoints\n0.0 0.0 0.0\nend kpoints\n\n")
        out.write("begin projections\n")
        for symbol, position in atoms:
            channels = "s;p" if shell == "sp" else "s;p;d"
            out.write(f"f={position[0]/cell[0]:.15f},{position[1]/cell[1]:.15f},{position[2]/cell[2]:.15f}:{channels}\n")
        out.write("end projections\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wfn", required=True, type=Path)
    parser.add_argument("--eigen", required=True, type=Path)
    parser.add_argument("--atom", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", required=True)
    parser.add_argument("--mesh", nargs=3, required=True, type=int)
    parser.add_argument("--cell", nargs=3, required=True, type=float)
    parser.add_argument("--num-bands", required=True, type=int)
    parser.add_argument("--wfn-states", type=int, help="number of states stored in wfn.bin")
    parser.add_argument("--num-wann", required=True, type=int)
    parser.add_argument("--sigma", default=1.0, type=float)
    parser.add_argument("--num-iter", default=100, type=int)
    parser.add_argument("--chunk-atoms", default=4, type=int)
    parser.add_argument("--symmetry", type=Path)
    parser.add_argument("--shell", choices=("sp", "spd"), default="spd")
    args = parser.parse_args()

    mesh = np.asarray(args.mesh)
    cell = np.asarray(args.cell)
    npoint = int(np.prod(mesh))
    wfn_states = args.wfn_states or args.num_bands
    if args.num_bands > wfn_states:
        raise ValueError("num-bands cannot exceed wfn-states")
    expected_bytes = npoint * wfn_states * 8
    if args.wfn.stat().st_size != expected_bytes:
        raise ValueError(f"wfn size is {args.wfn.stat().st_size}; expected {expected_bytes}")
    atoms = read_atoms(args.atom)
    nchan = 4 if args.shell == "sp" else 9
    if nchan * len(atoms) != args.num_wann:
        raise ValueError(f"complete {args.shell} shells require num_wann = {nchan} * number of atoms")

    psi_all = np.memmap(args.wfn, dtype=np.float64, mode="r", shape=(npoint, wfn_states), order="F")
    psi = psi_all[:, : args.num_bands]
    coordinates = grid_coordinates(mesh, cell)
    hvol = float(np.prod(cell / mesh))
    overlap = hvol * (psi.T @ psi)
    orth_residual = float(np.max(np.abs(overlap - np.eye(args.num_bands))))
    amn = build_amn(psi, coordinates, atoms, cell, args.sigma, hvol, args.chunk_atoms, args.shell)
    singular = np.linalg.svd(amn, compute_uv=False)
    mmn = build_mmn(psi, coordinates, cell, hvol)
    energies = read_energies(args.eigen, args.num_bands)

    args.output.mkdir(parents=True, exist_ok=True)
    base = args.output / args.seed
    write_win(base.with_suffix(".win"), args.num_bands, args.num_wann, cell, atoms, args.num_iter, args.shell)
    with base.with_suffix(".eig").open("w") as out:
        for band, energy in enumerate(energies, 1):
            out.write(f"{band:8d}{1:8d} {energy * HARTREE_TO_EV:23.15e}\n")
    with base.with_suffix(".amn").open("w") as out:
        out.write("SALMON conventional Full-WF pseudo-channel projections\n")
        out.write(f"{args.num_bands:10d}{1:10d}{args.num_wann:10d}\n")
        for projector in range(args.num_wann):
            for band in range(args.num_bands):
                out.write(f"{band+1:8d}{projector+1:8d}{1:8d} {amn[band,projector]:23.15e} {0.0:23.15e}\n")
    with base.with_suffix(".mmn").open("w") as out:
        out.write("SALMON conventional Full-WF Gamma overlaps\n")
        out.write(f"{args.num_bands:10d}{1:10d}{3:10d}\n")
        for matrix, gvec in zip(mmn, ((1, 0, 0), (0, 1, 0), (0, 0, 1))):
            out.write(f"{1:8d}{1:8d}{gvec[0]:8d}{gvec[1]:8d}{gvec[2]:8d}\n")
            for column in range(args.num_bands):
                for row in range(args.num_bands):
                    value = matrix[row, column]
                    out.write(f" {value.real:23.15e} {value.imag:23.15e}\n")
    np.savetxt(args.output / "projection_singular_values.dat", singular)
    dmn_text = ""
    if args.symmetry:
        operations = read_symmetry_operations(args.symmetry)
        d_wann, d_band, dmn_diagnostics = build_dmn(
            psi, atoms, mesh, cell, operations, amn, energies, hvol, args.shell
        )
        write_dmn(base.with_suffix(".dmn"), d_wann, d_band)
        with base.with_suffix(".win").open("a") as out:
            out.write("site_symmetry = true\nsymmetrize_eps = 1.0d-6\n")
        dmn_text = "".join(
            f"dmn_op={index} band_unitarity={values[0]:.15e} wann_unitarity={values[1]:.15e} "
            f"h_covariance={values[2]:.15e} amn_covariance={values[3]:.15e}\n"
            for index, values in enumerate(dmn_diagnostics, 1)
        )
    (args.output / "export_diagnostics.txt").write_text(
        f"orthogonality_max_abs={orth_residual:.15e}\n"
        f"amn_smax={singular[0]:.15e}\n"
        f"amn_smin={singular[-1]:.15e}\n"
        f"amn_condition={singular[0]/singular[-1]:.15e}\n" + dmn_text
    )
    print((args.output / "export_diagnostics.txt").read_text(), end="")


if __name__ == "__main__":
    main()
