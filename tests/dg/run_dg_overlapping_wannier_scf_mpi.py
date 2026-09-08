#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-scf-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"scf"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_solver.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_density.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_scf.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_scf_mpi.f90"),
      "-llapack","-lblas",
      "-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    fingerprints=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env,timeout=30)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier SCF on {n} ranks" in p.stdout
      match=re.search(r"fingerprint=(-?\d+)",p.stdout)
      assert match,p.stdout
      fingerprints.append(int(match.group(1)))
    assert len(set(fingerprints))==1,fingerprints
print("PASS overlapping-Wannier SCF fixture on 1, 2, 4, and 8 ranks")
