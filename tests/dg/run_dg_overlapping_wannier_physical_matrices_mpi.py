#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-physical-") as name:
    build=Path(name);(build/"config.h").write_text("")
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    fixtures=[
      ("nonlocal","dg_overlapping_wannier_nonlocal.f90","test_dg_overlapping_wannier_nonlocal_mpi.f90"),
      ("observables","dg_overlapping_wannier_observables.f90","test_dg_overlapping_wannier_observables_mpi.f90"),
    ]
    for label,module,test in fixtures:
      exe=build/label
      subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
        "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
        str(root/"src/gs/dc"/module),str(root/"tests/dg"/test),"-o",str(exe)],check=True)
      signatures=[]
      for n in (1,2,4,8):
        p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
        assert p.returncode==0,(label,n,p.stdout,p.stderr)
        assert f"PASS overlapping-Wannier {label} on {n} ranks" in p.stdout
        marker="NONLOCAL" if label=="nonlocal" else "OBSERVABLES"
        match=re.search(rf"{marker} ranks=\d+ values=([^\n]+)",p.stdout)
        assert match,p.stdout
        signatures.append([float(value) for value in match.group(1).split()])
      reference=signatures[0]
      assert all(max(abs(a-b) for a,b in zip(reference,row))<1e-13 for row in signatures[1:])
print("PASS overlapping-Wannier physical matrices on 1, 2, 4, and 8 ranks")
