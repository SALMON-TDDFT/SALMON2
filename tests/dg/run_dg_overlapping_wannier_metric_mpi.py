#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-metric-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"metric"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow",
      str(root/"src/gs/dc/dg_overlapping_wannier_metric.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_metric_mpi.f90"),
      "-llapack","-lblas","-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    matrices=[];ownership=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier metric on {n} ranks" in p.stdout
      match=re.search(r"METRIC ranks=\d+ ownership=(\d+)\n([^\n]+)",p.stdout)
      assert match,p.stdout
      ownership.append(int(match.group(1)))
      matrices.append([float(value) for value in match.group(2).split()])
    assert ownership==[4,4,4,4],ownership
    reference=matrices[0]
    assert all(len(row)==len(reference) for row in matrices)
    assert all(max(abs(a-b) for a,b in zip(reference,row))<1e-13 for row in matrices[1:])
print("PASS overlapping-Wannier metric fixture on 1, 2, 4, and 8 ranks")
