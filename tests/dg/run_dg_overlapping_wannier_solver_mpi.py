#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-solver-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"solver"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_solver.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_density.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_solver_mpi.f90"),
      "-llapack","-lblas","-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    signatures=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env,timeout=30)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier solver and density on {n} ranks" in p.stdout
      match=re.search(r"SOLVER ranks=\d+ signature=(-?\d+)",p.stdout)
      assert match,p.stdout
      signatures.append(int(match.group(1)))
    assert len(set(signatures))==1,signatures
print("PASS overlapping-Wannier solver and density fixture on 1, 2, 4, and 8 ranks")
