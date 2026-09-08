#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-divided-scf-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_divided_scf"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dc_scf_convergence.f90"),
    str(root/"src/gs/dc/dg_hybrid_divided_scf.f90"),
    str(root/"tests/dg/test_dg_hybrid_divided_scf_mpi.f90"),"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  iterations=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env,timeout=60)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid divided SCF on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_DIVIDED_SCF ranks=\d+ iterations=(\d+)",run.stdout);assert match,run.stdout
    iterations.append(int(match.group(1)))
  assert len(set(iterations))==1,iterations
print("PASS hybrid divided SCF on 1, 2, 4, and 8 ranks")
