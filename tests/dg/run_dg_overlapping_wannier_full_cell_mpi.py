#!/usr/bin/env python3
from pathlib import Path
import os, re, shutil, subprocess, tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-full-cell-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "full_cell"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/gs/dc/dg_overlapping_wannier_full_cell.f90"),
        str(root / "tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90"),
        "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    signatures = []
    for nrank in (1, 2, 4, 8):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(nrank), str(exe)],
                             capture_output=True, text=True, env=env)
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)
        assert f"PASS full-cell tiled projection on {nrank} ranks" in run.stdout
        match = re.search(r"FULL_CELL ranks=\d+ values=([^\n]+)", run.stdout)
        assert match, run.stdout
        signatures.append([float(value) for value in match.group(1).split()])
    reference = signatures[0]
    assert all(max(abs(a-b) for a,b in zip(reference,row)) < 1e-12
               for row in signatures[1:])
print("PASS full-cell tiled projection on 1, 2, 4, and 8 ranks")
