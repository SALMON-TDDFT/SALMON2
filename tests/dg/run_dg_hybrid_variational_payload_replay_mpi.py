#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-variational-replay-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "hybrid_variational_replay"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008", "-ffree-line-length-none",
        "-I", str(build), "-J", str(build), "-fcheck=all", "-fbacktrace",
        str(root / "src/gs/dc/dg_hybrid_variational_payload.f90"),
        str(root / "tests/dg/test_dg_hybrid_variational_payload_replay_mpi.f90"), "-o", str(exe),
    ], check=True)
    env = os.environ.copy(); env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for nrank in (1, 2, 4):
        prefix = build / f"payload-{nrank}"
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(nrank), str(exe), str(prefix)],
                             capture_output=True, text=True, env=env)
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)
        assert f"PASS variational payload replay on {nrank} ranks" in run.stdout
print("PASS variational payload replay on 1, 2, and 4 ranks")
