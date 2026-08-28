#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
scf_source = (root / "src/gs/dc/dg_hybrid_continuation_scf.f90").read_text().lower()
assert "s_dg_hybrid_continuation_callbacks" not in scf_source, (
    "the rejected public continuation callback bundle is still present"
)
assert "subroutine run_dg_hybrid_continuation_scf_fixture" in scf_source, (
    "the synthetic fixture entry point is missing"
)
with tempfile.TemporaryDirectory(prefix="hybrid-continuation-scf-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "hybrid_continuation_scf"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
        "-ffree-line-length-none", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/common/dg_hybrid_continuation_residuals.f90"),
        str(root / "src/common/dg_hybrid_continuation_acceptance.f90"),
        str(root / "src/gs/dc/dg_hybrid_continuation_state.f90"),
        str(root / "src/gs/dc/dg_hybrid_continuation_controller.f90"),
        str(root / "src/gs/dc/dg_hybrid_continuation_scf.f90"),
        str(root / "tests/dg/test_dg_hybrid_continuation_scf_mpi.f90"),
        "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for nrank in (1, 2, 4, 8):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(nrank), str(exe)],
                             capture_output=True, text=True, env=env)
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)
        assert f"PASS hybrid continuation SCF on {nrank} ranks" in run.stdout
print("PASS hybrid continuation SCF on 1, 2, 4, and 8 ranks")
