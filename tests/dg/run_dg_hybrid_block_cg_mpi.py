#!/usr/bin/env python3
from pathlib import Path
import os,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
if shutil.which("pkg-config"):
  lapack=shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
else:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();lapack=[f"-L{prefix}/lib","-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-block-cg-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_block_cg"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dg_hybrid_block_cg.f90"),str(root/"tests/dg/test_dg_hybrid_block_cg_mpi.f90"),*lapack,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid block CG on {nrank} ranks" in run.stdout
print("PASS hybrid block CG on 1, 2, 4, and 8 ranks")
