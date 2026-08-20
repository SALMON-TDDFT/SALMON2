#!/usr/bin/env python3
from pathlib import Path
import os,shlex,shutil,subprocess,tempfile,re
root=Path(__file__).resolve().parents[2]
if os.environ.get("SALMON_LAPACK_LIBS"): libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
else: libs=["-llapack","-lblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-length-gauge-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_length_gauge"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero","-fbacktrace",
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),str(root/"src/common/dg_hybrid_sparse_operators.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_metric_solver.f90"),str(root/"src/rt/dg/rt_dg_hybrid_length_gauge.f90"),
    str(root/"tests/dg/test_rt_dg_hybrid_length_gauge_mpi.f90"),*libs,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");fps=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid length gauge on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_LENGTH_GAUGE ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout;fps.append(int(match.group(1)))
  assert len(set(fps))==1,fps
print("PASS hybrid length gauge on 1, 2, 4, and 8 ranks")
