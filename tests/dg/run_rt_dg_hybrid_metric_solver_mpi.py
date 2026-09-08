#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
if os.environ.get("SALMON_LAPACK_LIBS"):
  lapack_libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("pkg-config") and subprocess.run(["pkg-config","--exists","openblas"],check=False).returncode==0:
  lapack_libs=shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
elif shutil.which("brew"):
  probe=subprocess.run(["brew","--prefix","openblas"],capture_output=True,text=True)
  lapack_libs=[f"-L{probe.stdout.strip()}/lib","-lopenblas"] if probe.returncode==0 else ["-llapack","-lblas"]
else:
  lapack_libs=["-llapack","-lblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-metric-solver-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_metric_solver"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_metric_solver.f90"),
    str(root/"tests/dg/test_rt_dg_hybrid_metric_solver_mpi.f90"),*lapack_libs,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid metric solver on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_METRIC_SOLVER ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid metric solver on 1, 2, 4, and 8 ranks")
