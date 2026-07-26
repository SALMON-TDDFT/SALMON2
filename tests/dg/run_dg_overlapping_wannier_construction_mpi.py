#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-construction-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"construction"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_metric.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_construction.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_construction_mpi.f90"),
      "-llapack","-lblas","-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    signatures=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier construction on {n} ranks" in p.stdout
      match=re.search(r"CONSTRUCTION ranks=\d+ fingerprint=(-?\d+) centers=([^\n]+)",p.stdout)
      assert match,p.stdout
      signatures.append((int(match.group(1)),tuple(int(x) for x in match.group(2).split())))
    assert len(set(signatures))==1,signatures
print("PASS overlapping-Wannier construction fixture on 1, 2, 4, and 8 ranks")
