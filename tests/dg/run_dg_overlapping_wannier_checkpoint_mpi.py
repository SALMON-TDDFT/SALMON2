#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-checkpoint-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"checkpoint"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_checkpoint.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90"),"-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    env["OW_CHECKPOINT_PREFIX"]=str(build/"route-checkpoint")
    fingerprints=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env,timeout=30)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier checkpoint on {n} ranks" in p.stdout
      fingerprints.append(re.search(r"fingerprint=([0-9A-F]+)",p.stdout).group(1))
    assert len(set(fingerprints))==1,fingerprints
print("PASS overlapping-Wannier checkpoint fixture on 1, 2, 4, and 8 ranks")
