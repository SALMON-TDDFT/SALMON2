#!/usr/bin/env python3
"""Project polished conventional WFs and a Wannier90 gauge onto DC fragment bases."""

import argparse
import re
from pathlib import Path

import numpy as np

AU_LENGTH_AA = 0.529177210903
GLOBAL_MAGIC = -22022216
FRAGMENT_MAGIC = -22022218


def write_fortran(stream, array, dtype):
    np.asarray(array, dtype=dtype).ravel(order="F").tofile(stream)


def read_i4(stream, count):
    return np.fromfile(stream, dtype="<i4", count=count)


def read_f8(stream, count):
    return np.fromfile(stream, dtype="<f8", count=count)


def read_basis(path):
    with path.open("rb") as stream:
        shape = read_i4(stream, 3)
        nspin, nstate = read_i4(stream, 2)
        nbasis = read_i4(stream, nspin)
        phi = read_f8(stream, int(np.prod(shape)) * nspin * nstate).reshape(
            (*shape, nspin, nstate), order="F"
        )
    return shape, nbasis, phi


def read_grid(path):
    with path.open("rb") as stream:
        local = read_i4(stream, 3)
        total = read_i4(stream, 3)
        indices = [read_i4(stream, int(n)) - 1 for n in local]
    return total, indices


def read_wavefunction_header(path):
    with path.open("rb") as stream:
        nfrag, nspin, nstate_frag, _ = read_i4(stream, 4)
        nmat = read_i4(stream, nspin)
        nbasis = read_i4(stream, nfrag * nspin)
        index_basis = read_i4(stream, nstate_frag * nfrag * nspin)
    return nfrag, nspin, nstate_frag, nmat, nbasis, index_basis


def read_energies(path, count):
    values = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit():
            values.append(float(fields[1]))
    if len(values) < count:
        raise ValueError("eigenvalue file has too few states")
    return np.asarray(values[:count], dtype="<f8")


def read_u_matrix(path):
    lines = path.read_text().splitlines()
    dimensions = [int(value) for value in lines[1].split()]
    if dimensions[0] != 1:
        raise ValueError("only Gamma-point u matrices are supported")
    ncolumn, nrow = dimensions[1], dimensions[2]
    values = np.loadtxt(lines[4:])
    if values.shape != (nrow * ncolumn, 2):
        raise ValueError(f"unexpected u.mat payload shape {values.shape}")
    return (values[:, 0] + 1j * values[:, 1]).reshape((nrow, ncolumn), order="F")


def read_u(path, nband, nwann, u_dis_path=None):
    gauge = read_u_matrix(path)
    if gauge.shape == (nband, nwann):
        return gauge
    if gauge.shape != (nwann, nwann) or u_dis_path is None:
        raise ValueError(f"unexpected u.mat shape {gauge.shape}")
    disentanglement = read_u_matrix(u_dis_path)
    if disentanglement.shape != (nband, nwann):
        raise ValueError(f"unexpected u_dis.mat shape {disentanglement.shape}")
    return disentanglement @ gauge


def read_rdat(path, nwann):
    with path.open() as stream:
        stream.readline()
        if int(stream.readline()) != nwann:
            raise ValueError("r.dat num_wann mismatch")
        nrpts = int(stream.readline())
        position = np.zeros((3, nwann, nwann), dtype=np.complex128)
        for _ in range(nrpts * nwann * nwann):
            fields = stream.readline().split()
            rvec = tuple(int(value) for value in fields[:3])
            if rvec != (0, 0, 0):
                continue
            row, column = int(fields[3]) - 1, int(fields[4]) - 1
            values = [float(value) for value in fields[5:11]]
            position[:, row, column] = [values[0] + 1j * values[1],
                                               values[2] + 1j * values[3],
                                               values[4] + 1j * values[5]]
    return position / AU_LENGTH_AA


def read_spreads(path, nwann):
    pattern = re.compile(r"WF centre and spread\s+\d+\s+\(.*?\)\s+([0-9.Ee+-]+)")
    values = [float(match.group(1)) for match in pattern.finditer(path.read_text())]
    if len(values) < nwann:
        raise ValueError("could not find final Wannier spreads in wout")
    return np.asarray(values[-nwann:]) / (AU_LENGTH_AA**2)


def fragment_owner(centers, cell, fragment_grid):
    wrapped = np.mod(centers, cell[:, None])
    coordinate = np.floor(wrapped / (cell / fragment_grid)[:, None]).astype(int)
    coordinate = np.minimum(coordinate, fragment_grid[:, None] - 1)
    return ((coordinate[0] * fragment_grid[1] + coordinate[1]) * fragment_grid[2] + coordinate[2] + 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wfn", type=Path, required=True)
    parser.add_argument("--wfn-states", type=int, required=True)
    parser.add_argument("--u-mat", type=Path, required=True)
    parser.add_argument("--u-dis-mat", type=Path)
    parser.add_argument("--r-dat", type=Path, required=True)
    parser.add_argument("--wout", type=Path, required=True)
    parser.add_argument("--eigen", type=Path, required=True)
    parser.add_argument("--data-dcdft", type=Path, required=True)
    parser.add_argument("--mesh", nargs=3, type=int, required=True)
    parser.add_argument("--cell", nargs=3, type=float, required=True)
    parser.add_argument("--fragment-grid", nargs=3, type=int, required=True)
    parser.add_argument("--num-bands", type=int, required=True)
    parser.add_argument("--num-wann", type=int, required=True)
    args = parser.parse_args()

    mesh = np.asarray(args.mesh)
    cell = np.asarray(args.cell)
    fragment_grid = np.asarray(args.fragment_grid)
    npoint = int(np.prod(mesh))
    psi = np.memmap(args.wfn, dtype="<f8", mode="r", shape=(npoint, args.wfn_states), order="F")[:, : args.num_bands]
    gauge = read_u(args.u_mat, args.num_bands, args.num_wann, args.u_dis_mat)
    energies = read_energies(args.eigen, args.num_bands)
    position = read_rdat(args.r_dat, args.num_wann)
    centers = np.real(np.column_stack([np.diag(position[axis]) for axis in range(3)]).T)
    spreads = read_spreads(args.wout, args.num_wann)
    owners = fragment_owner(centers, cell, fragment_grid).astype("<i4")
    hvol = float(np.prod(cell / mesh))

    fragments = sorted((args.data_dcdft / "fragments").glob("[0-9][0-9][0-9][0-9][0-9][0-9]"))
    if len(fragments) != int(np.prod(fragment_grid)):
        raise ValueError("fragment count does not match fragment-grid")
    projected_overlap = np.zeros((args.num_wann, args.num_wann), dtype=np.complex128)
    fragment_coefficients = []
    state_projections = []
    projected_state_overlap = np.zeros((args.num_bands, args.num_bands), dtype=np.float64)
    for fragment in fragments:
        shape, nbasis_by_spin, phi = read_basis(fragment / "basis_functions.bin")
        total, indices = read_grid(fragment / "rgrid_index.bin")
        if np.any(total != mesh):
            raise ValueError(f"global grid mismatch in {fragment}")
        core_indices = []
        for axis, core_size in zip(indices, shape):
            excess = len(axis) - int(core_size)
            if excess < 0 or excess % 2:
                raise ValueError(f"buffer/core grid mismatch in {fragment}")
            edge = excess // 2
            core_indices.append(axis[edge : edge + int(core_size)])
        nbasis = int(nbasis_by_spin[0])
        global_indices = np.ravel_multi_index(
            np.meshgrid(*core_indices, indexing="ij"), tuple(mesh), order="F"
        ).ravel(order="F")
        projection = hvol * (phi[:, :, :, 0, :nbasis].reshape((-1, nbasis), order="F").T @ psi[global_indices, :])
        coefficient = np.asfortranarray(projection @ gauge, dtype="<c16")
        projected_overlap += coefficient.conj().T @ coefficient
        fragment_coefficients.append((fragment, nbasis, coefficient))
        source_seed = fragment / "wavefunctions.bin"
        nfrag, nspin, nstate_frag, nmat, nbasis_all, index_basis = read_wavefunction_header(source_seed)
        if nspin != 1 or nstate_frag < nbasis:
            raise ValueError(f"unsupported wavefunction header in {fragment}")
        projected_state_overlap += projection.T @ projection
        state_projections.append((fragment, nfrag, nspin, nstate_frag, nbasis, nmat, nbasis_all, index_basis, projection))

    state_eigenvalues, state_eigenvectors = np.linalg.eigh(projected_state_overlap)
    if state_eigenvalues[0] <= 1.0e-10:
        raise ValueError("projected Full-state metric is rank deficient")
    state_lowdin = (state_eigenvectors / np.sqrt(state_eigenvalues)) @ state_eigenvectors.T
    for fragment, nfrag, nspin, nstate_frag, nbasis, nmat, nbasis_all, index_basis, projection in state_projections:
        state_coef = np.zeros((nstate_frag, args.num_bands), dtype="<f8", order="F")
        state_coef[:nbasis, :] = projection @ state_lowdin
        with (fragment / "wavefunctions_full_polish_seed.bin").open("wb") as stream:
            np.asarray([nfrag, nspin, nstate_frag, args.num_bands], dtype="<i4").tofile(stream)
            nmat.astype("<i4").tofile(stream)
            nbasis_all.astype("<i4").tofile(stream)
            index_basis.astype("<i4").tofile(stream)
            write_fortran(stream, state_coef, "<f8")
            write_fortran(stream, energies[:, None], "<f8")

    overlap_eigenvalues, overlap_eigenvectors = np.linalg.eigh(projected_overlap)
    if overlap_eigenvalues[0] <= 1.0e-10:
        raise ValueError("projected Wannier metric is rank deficient")
    lowdin = (overlap_eigenvectors / np.sqrt(overlap_eigenvalues)) @ overlap_eigenvectors.conj().T
    for fragment, nbasis, coefficient in fragment_coefficients:
        coefficient = np.asfortranarray(coefficient @ lowdin, dtype="<c16")
        with (fragment / "full_polish_wannier_coef.bin").open("wb") as stream:
            np.asarray([FRAGMENT_MAGIC, 1, nbasis, args.num_wann, 1], dtype="<i4").tofile(stream)
            write_fortran(stream, coefficient, "<c16")
    gauge = gauge @ lowdin
    position = np.einsum("ij,ajk,kl->ail", lowdin.conj().T, position, lowdin)
    centers = np.real(np.column_stack([np.diag(position[axis]) for axis in range(3)]).T)
    owners = fragment_owner(centers, cell, fragment_grid).astype("<i4")

    total_dir = args.data_dcdft / "total"
    total_dir.mkdir(exist_ok=True)
    with (total_dir / "wannier90_global_basis.bin").open("wb") as stream:
        np.asarray([GLOBAL_MAGIC, 2, args.num_bands, args.num_wann, len(fragments)], dtype="<i4").tofile(stream)
        owners.tofile(stream)
        write_fortran(stream, centers, "<f8")
        np.asarray(spreads, dtype="<f8").tofile(stream)
        write_fortran(stream, gauge, "<c16")
        np.asarray([1], dtype="<i4").tofile(stream)
        write_fortran(stream, position, "<c16")
    print(f"wrote {len(fragments)} fragment coefficient files and global metadata")
    overlap_residual = np.linalg.norm(projected_overlap - np.eye(args.num_wann), ord=2)
    print(f"projected_wannier_overlap_residual_2norm={overlap_residual:.8e}")
    print(f"projected_wannier_overlap_eigenvalue_min={overlap_eigenvalues[0]:.8e}")
    print(f"projected_wannier_overlap_eigenvalue_max={overlap_eigenvalues[-1]:.8e}")
    print(f"projected_state_overlap_eigenvalue_min={state_eigenvalues[0]:.8e}")
    print(f"projected_state_overlap_eigenvalue_max={state_eigenvalues[-1]:.8e}")


if __name__ == "__main__":
    main()
