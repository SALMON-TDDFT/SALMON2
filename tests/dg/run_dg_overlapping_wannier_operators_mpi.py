#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-operators-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"operators"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_operators.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_operators_mpi.f90"),
      "-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    signatures=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier weak operators on {n} ranks" in p.stdout
      match=re.search(r"OPERATORS ranks=\d+ values=([^\n]+)",p.stdout)
      assert match,p.stdout
      signatures.append([float(value) for value in match.group(1).split()])
    reference=signatures[0]
    assert all(max(abs(a-b) for a,b in zip(reference,row))<1e-13 for row in signatures[1:])
print("PASS overlapping-Wannier weak operators fixture on 1, 2, 4, and 8 ranks")
