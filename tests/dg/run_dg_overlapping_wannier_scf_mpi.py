#!/usr/bin/env python3
from pathlib import Path
import os,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-scf-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"scf"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_solver.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_density.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_scf.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_scf_mpi.f90"),
      "-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env,timeout=30)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier SCF on {n} ranks" in p.stdout
print("PASS overlapping-Wannier SCF fixture on 1, 2, 4, and 8 ranks")
