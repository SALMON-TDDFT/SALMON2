#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-metric-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_metric"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),
    str(root/"tests/dg/test_dg_hybrid_sparse_metric_mpi.f90"),"-llapack","-lblas","-o",str(exe)],check=True)
  env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid sparse metric on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_METRIC ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid sparse metric on 1, 2, 4, and 8 ranks")
