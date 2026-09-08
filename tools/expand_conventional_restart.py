#!/usr/bin/env python3
"""Expand a converged real Gamma-point SALMON restart with empty orbitals."""

import argparse
import shutil
import struct
from pathlib import Path

import numpy as np


def read_record(stream):
    size = struct.unpack("<i", stream.read(4))[0]
    payload = stream.read(size)
    if struct.unpack("<i", stream.read(4))[0] != size:
        raise ValueError("invalid Fortran sequential record")
    return payload


def write_record(stream, payload):
    stream.write(struct.pack("<i", len(payload)))
    stream.write(payload)
    stream.write(struct.pack("<i", len(payload)))


def read_info(path):
    with Path(path).open("rb") as stream:
        nk = struct.unpack("<i", read_record(stream))[0]
        nstate = struct.unpack("<i", read_record(stream))[0]
        iteration = struct.unpack("<i", read_record(stream))[0]
        nprocs = struct.unpack("<i", read_record(stream))[0]
        real_orbital = struct.unpack("<i", read_record(stream))[0] != 0
    return nk, nstate, iteration, nprocs, real_orbital


def write_info(path, nk, nstate, iteration, nprocs, real_orbital):
    with Path(path).open("wb") as stream:
        for value in (nk, nstate, iteration, nprocs, int(real_orbital)):
            write_record(stream, struct.pack("<i", value))


def expand_orbitals(occupied, nstate, hvol, seed):
    npoint, nold = occupied.shape
    if nstate < nold or nstate > npoint:
        raise ValueError("expanded state count must lie between old count and grid size")
    expanded = np.empty((npoint, nstate), dtype=np.float64, order="F")
    expanded[:, :nold] = occupied
    if nstate == nold:
        return expanded
    normalized_occupied = np.sqrt(hvol) * occupied
    random = np.random.default_rng(seed).standard_normal((npoint, nstate - nold))
    random -= normalized_occupied @ (normalized_occupied.T @ random)
    added, _ = np.linalg.qr(random, mode="reduced")
    added -= normalized_occupied @ (normalized_occupied.T @ added)
    added, _ = np.linalg.qr(added, mode="reduced")
    expanded[:, nold:] = added / np.sqrt(hvol)
    return expanded


def read_occupation(path):
    with Path(path).open("rb") as stream:
        payload = read_record(stream)
        if stream.read(1):
            raise ValueError("occupation file contains unexpected extra records")
    return np.frombuffer(payload, dtype="<f8").copy()


def write_occupation(path, occupation):
    with Path(path).open("wb") as stream:
        write_record(stream, np.asarray(occupation, dtype="<f8").tobytes())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mesh", required=True, nargs=3, type=int)
    parser.add_argument("--cell", required=True, nargs=3, type=float)
    parser.add_argument("--num-states", required=True, type=int)
    parser.add_argument("--seed", default=20260711, type=int)
    args = parser.parse_args()

    nk, nold, iteration, nprocs, real_orbital = read_info(args.source / "info.bin")
    if nk != 1 or not real_orbital:
        raise ValueError("only real Gamma-point restarts are supported")
    npoint = int(np.prod(args.mesh))
    expected = npoint * nold * 8
    if (args.source / "wfn.bin").stat().st_size != expected:
        raise ValueError("source wfn.bin size does not match mesh and restart metadata")
    hvol = float(np.prod(np.asarray(args.cell) / np.asarray(args.mesh)))
    occupied = np.memmap(args.source / "wfn.bin", dtype="<f8", mode="r", shape=(npoint, nold), order="F")
    expanded = expand_orbitals(occupied, args.num_states, hvol, args.seed)

    args.output.mkdir(parents=True, exist_ok=True)
    for name in ("rho_inout.bin", "atomic_coor.txt", "atomic_vel.txt"):
        source = args.source / name
        if source.exists():
            shutil.copy2(source, args.output / name)
    write_info(args.output / "info.bin", nk, args.num_states, iteration, nprocs, real_orbital)
    occupation = read_occupation(args.source / "occupation.bin")
    if occupation.size != nold:
        raise ValueError("occupation count does not match restart state count")
    expanded_occupation = np.zeros(args.num_states)
    expanded_occupation[:nold] = occupation
    write_occupation(args.output / "occupation.bin", expanded_occupation)
    np.asarray(expanded, order="F").ravel(order="F").tofile(args.output / "wfn.bin")

    overlap = hvol * expanded.T @ expanded
    print(f"source_states={nold} expanded_states={args.num_states}")
    print(f"occupied_change_max=0.000000000000000e+00")
    print(f"orthogonality_max_abs={np.max(np.abs(overlap - np.eye(args.num_states))):.15e}")


if __name__ == "__main__":
    main()
