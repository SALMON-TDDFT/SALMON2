#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-fragment-basis-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_fragment_basis"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dg_hybrid_fragment_basis.f90"),
    str(root/"tests/dg/test_dg_hybrid_fragment_basis_mpi.f90"),"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  fingerprints=[]
  for nrank in (1,2,4):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid fragment basis on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_FRAGMENT_BASIS ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid fragment basis on 1, 2, and 4 ranks")
