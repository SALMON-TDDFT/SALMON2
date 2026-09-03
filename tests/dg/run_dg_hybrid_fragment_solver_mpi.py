#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
if shutil.which("pkg-config"):
  libs=shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
else:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip()
  libs=[f"-L{prefix}/lib","-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-fragment-solver-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_fragment_solver"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dg_hybrid_fragment_basis.f90"),str(root/"src/gs/dc/dg_hybrid_fragment_solver.f90"),
    str(root/"src/gs/occupation_kernel.f90"),str(root/"src/gs/dc/dc_fragment_occupation.f90"),
    str(root/"tests/dg/test_dg_hybrid_fragment_solver_mpi.f90"),*libs,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  fingerprints=[];workspaces=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],
      capture_output=True,text=True,env=env,timeout=30)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid fragment solver on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_FRAGMENT_SOLVER ranks=\d+ fingerprint=(-?\d+) workspace=(\d+)",run.stdout);assert match
    fingerprints.append(int(match.group(1)))
    workspaces.append(int(match.group(2)))
  assert len(set(fingerprints))==1,fingerprints
  assert len(set(workspaces))==1,workspaces
print("PASS hybrid fragment solver on 1, 2, 4, and 8 ranks")
