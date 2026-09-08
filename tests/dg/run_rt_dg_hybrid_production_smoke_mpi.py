#!/usr/bin/env python3
"""Exercise the real SALMON PP/Hartree/XC hybrid callback at t=0 and one step."""

from pathlib import Path
import os
import math
import re
import shlex
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
main_source = (root / "src/rt/main_tddft.f90").read_text()


def require_hybrid_route_contract(source: str) -> None:
    continuation = source.split("subroutine run_dg_hybrid_continuation_rt()", 1)[1].split(
        "end subroutine run_dg_hybrid_continuation_rt", 1)[0]
    projection = source.split("subroutine project_salmon_local_rows", 1)[1].split(
        "end subroutine project_salmon_local_rows", 1)[0]
    assert "[HYBRID-RT-ROUTE] propagator=EXP potential=PP+HARTREE+XC" in continuation
    assert "[HYBRID-RT-STEP] step=" in continuation and "path=PP+HARTREE+XC" in continuation
    assert "call propagate_rt_dg_hybrid_length_gauge" in continuation
    for token in ("call hartree", "call exchange_correlation_density", "call update_vlocal"):
        assert token in projection.lower(), token


require_hybrid_route_contract(main_source)
for old, replacement in (
    ("[HYBRID-RT-ROUTE]", "[REMOVED-HYBRID-RT-ROUTE]"),
    ("[HYBRID-RT-STEP]", "[REMOVED-HYBRID-RT-STEP]"),
    ("potential=PP+HARTREE+XC", "potential=CONVENTIONAL"),
):
    mutated = main_source.replace(old, replacement)
    try:
        require_hybrid_route_contract(mutated)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"route contract mutation survived: {old}")

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
    configured_build = os.environ.get("SALMON_DG_PRODUCTION_BUILD")
    if configured_build:
        salmon_build = Path(configured_build).resolve()
        cache = (salmon_build / "CMakeCache.txt").read_text()
        assert "USE_MPI:BOOL=ON" in cache and "USE_SCALAPACK:BOOL=ON" in cache
    else:
        salmon_build = work / "salmon-build"
        subprocess.run(
            [
                "cmake", "-S", str(root), "-B", str(salmon_build),
                "-DUSE_MPI=ON", "-DUSE_SCALAPACK=ON", "-DUSE_EIGENEXA=OFF", "-DUSE_WANNIER90=OFF",
                f"-DCMAKE_Fortran_COMPILER={mpifort}", "-DCMAKE_BUILD_TYPE=Debug",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            ["cmake", "--build", str(salmon_build), "-j2"],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    salmon = salmon_build / "salmon"
    assert salmon.exists(), "current-source SALMON build did not produce the executable"
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
 yn_scalapack='y'
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
        assert run.stdout.count("[HYBRID-RT-ROUTE] propagator=EXP potential=PP+HARTREE+XC") == 1, run.stdout
        assert run.stdout.count("[HYBRID-RT-HANDOFF]") == 1, run.stdout
        assert run.stdout.count("[HYBRID-RT-POTENTIAL] update=0 path=PP+HARTREE+XC") == 1, run.stdout
        steps = re.findall(
            r"\[HYBRID-RT-STEP\] step=(\d+) metric_norm=\s*([^ ]+) orbital_energy=\s*([^ ]+) "
            r"polarization_norm=\s*([^ ]+) density_norm=\s*([^ ]+) update_count=(\d+) path=PP\+HARTREE\+XC",
            run.stdout,
        )
        assert len(steps) == 1 and steps[0][0] == "1" and steps[0][5] == "2", run.stdout
        observables = [float(value) for value in steps[0][1:5]]
        assert all(math.isfinite(value) for value in observables), observables
        assert observables[0] > 0.0 and observables[3] >= 0.0, observables
        assert "hybrid DG RT Hartree/XC update failed" not in run.stdout + run.stderr
        assert "production DG requires a build with MPI and ScaLAPACK support" not in run.stdout + run.stderr
        assert "[DG-OW-RT]" not in run.stdout, "wrong coefficient-only RT route"
        for line in run.stdout.splitlines():
            if line.startswith("[HYBRID-RT-"):
                print(f"ranks={nrank} {line}")

print("PASS production hybrid GS-to-RT smoke on 1, 2, and 4 ranks")
