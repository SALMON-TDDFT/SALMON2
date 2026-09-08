#!/usr/bin/env python3
from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]


def lapack_libraries():
    configured = os.environ.get("SALMON_LAPACK_LIBS")
    if configured:
        return shlex.split(configured)
    if shutil.which("pkg-config"):
        for package in ("openblas", "lapack"):
            probe = subprocess.run(
                ["pkg-config", "--libs", package], capture_output=True, text=True
            )
            if probe.returncode == 0:
                return shlex.split(probe.stdout)
    if shutil.which("brew"):
        probe = subprocess.run(
            ["brew", "--prefix", "openblas"], capture_output=True, text=True
        )
        if probe.returncode == 0:
            return [f"-L{probe.stdout.strip()}/lib", "-lopenblas"]
    return ["-llapack", "-lblas"]

operator_source = (root / "src/gs/dc/dg_hybrid_divided_operator.f90").read_text()
payload_source = (root / "src/gs/dc/dg_hybrid_variational_payload.f90").read_text()
map_source = (root / "src/common/dg_hybrid_wannier_complement.f90").read_text()
shape_gate = operator_source.find("terminal transform shape differs between ranks")
transform_broadcast = operator_source.find("MPI_Bcast(reference_transform")
assert shape_gate >= 0 and transform_broadcast >= 0 and shape_gate < transform_broadcast, \
    "terminal transform shape agreement must precede every count-dependent broadcast"
compose_start = payload_source.index("subroutine compose_dg_hybrid_variational_hamiltonian")
compose_end = payload_source.index("end subroutine compose_dg_hybrid_variational_hamiltonian", compose_start)
compose_source = payload_source[compose_start:compose_end]
assert "stat=allocation_status" in compose_source and "move_alloc" in compose_source
assert "allocate(iterate%" not in compose_source, \
    "variational composition must publish allocations only after collective success"
assert "compute_dg_hybrid_union_to_complete_binding" in map_source
assert "recomputed_binding_fingerprint==complete_transform_binding_fingerprint" in operator_source, \
    "terminal composition must authenticate the raw Task 5 transform binding"

with tempfile.TemporaryDirectory(prefix="hybrid-divided-operator-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "hybrid_divided_operator"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
        "-ffree-line-length-none", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/common/dg_hybrid_wannier_complement.f90"),
        str(root / "src/gs/dc/dg_hybrid_variational_payload.f90"),
        str(root / "src/gs/dc/dg_hybrid_fragment_basis.f90"),
        str(root / "src/gs/dc/dg_hybrid_divided_operator.f90"),
        str(root / "tests/dg/test_dg_hybrid_divided_operator_mpi.f90"),
        *lapack_libraries(),
        "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    identity_fingerprints = []
    rectangular_fingerprints = []
    for nrank in (1, 2, 4, 8):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(nrank), str(exe)],
            capture_output=True, text=True, env=env, timeout=30,
        )
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)
        assert f"PASS hybrid divided operator on {nrank} ranks" in run.stdout
        match = re.search(
            r"HYBRID_DIVIDED_OPERATOR ranks=\d+ identity=(-?\d+) rectangular=(-?\d+)",
            run.stdout,
        )
        assert match, run.stdout
        identity_fingerprints.append(int(match.group(1)))
        rectangular_fingerprints.append(int(match.group(2)))
    assert len(set(identity_fingerprints)) == 1, identity_fingerprints
    assert len(set(rectangular_fingerprints)) == 1, rectangular_fingerprints
    assert identity_fingerprints[0] != rectangular_fingerprints[0]

print("PASS hybrid divided operator on 1, 2, 4, and 8 ranks")
