#!/usr/bin/env python3
"""Exercise the real SALMON PP/Hartree/XC hybrid callback at t=0 and one step."""

from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
salmon = root / "build-hybrid-commit" / "salmon"
assert salmon.exists(), "build-hybrid-commit/salmon is required"

if os.environ.get("SALMON_LAPACK_LIBS"):
    libs = shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(
    ["brew", "--prefix", "openblas"], capture_output=True
).returncode == 0:
    prefix = subprocess.check_output(["brew", "--prefix", "openblas"], text=True).strip()
    libs = [f"-L{prefix}/lib", "-lopenblas"]
else:
    libs = ["-llapack", "-lblas"]

env = os.environ.copy()
env["OMP_NUM_THREADS"] = "1"
env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
mpiexec = shutil.which("mpiexec")
mpifort = shutil.which("mpifort")

with tempfile.TemporaryDirectory(prefix="hybrid-production-smoke-") as name:
    work = Path(name)
    (work / "config.h").write_text("")
    writer = work / "write_hybrid_checkpoint"
    sources = [
        "src/common/dg_hybrid_sparse_metric.f90",
        "src/common/dg_hybrid_sparse_operators.f90",
        "src/rt/dg/rt_dg_hybrid_checkpoint.f90",
        "src/rt/dg/rt_dg_hybrid_initialization.f90",
        "src/rt/dg/rt_dg_hybrid_density_update.f90",
        "tests/dg/test_rt_dg_hybrid_initialization_mpi.f90",
    ]
    subprocess.run(
        [
            mpifort,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(work),
            "-J",
            str(work),
            *[str(root / source) for source in sources],
            *libs,
            "-o",
            str(writer),
        ],
        check=True,
    )
    shutil.copy2(root / "testsuites/pseudo/Si.cpi", work / "Si.cpi")
    (work / "atom.dat").write_text("  'Si' 0.0 0.0 0.0 1\n")

    for nrank in (1, 2, 4):
        checkpoint = work / "hybrid_dg_ground_state.chk"
        written = subprocess.run(
            [mpiexec, "-n", str(nrank), str(writer), str(checkpoint), "write_production", "2560"],
            cwd=work,
            env=env,
            capture_output=True,
            text=True,
        )
        assert written.returncode == 0, (nrank, written.stdout, written.stderr)
        input_text = f"""
&calculation
 theory='tddft_response'
/
&control
 sysname='hybrid_production_smoke'
/
&units
 unit_system='a.u.'
/
&parallel
 nproc_k=1
 nproc_ob=1
 nproc_rgrid={nrank},1,1
 yn_eigenexa='n'
/
&system
 yn_periodic='y'
 al=40d0,8d0,8d0
 nelem=1
 nstate=2
 nelec=4
 natom=1
 file_atom_coor='atom.dat'
/
&pseudo
 izatom(1)=14
 file_pseudo(1)='Si.cpi'
 lloc_ps(1)=2
/
&functional
 xc='PZ'
/
&rgrid
 num_rgrid=40,8,8
/
&kgrid
 num_kgrid=1,1,1
/
&tgrid
 dt=0.02d0
 nt=1
/
&propagation
 yn_rt_dg_hybrid_continuation='y'
 yn_dg_length_gauge='y'
/
&emfield
 ae_shape1='none'
/
"""
        run = subprocess.run(
            [mpiexec, "-n", str(nrank), str(salmon)],
            input=input_text,
            cwd=work,
            env=env,
            capture_output=True,
            text=True,
            timeout=180,
        )
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)

print("PASS production hybrid GS-to-RT smoke on 1, 2, and 4 ranks")
