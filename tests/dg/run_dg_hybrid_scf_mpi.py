#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
if shutil.which("pkg-config"):
  libraries=shlex.split(subprocess.check_output(["pkg-config","--libs","scalapack"],text=True))+\
    shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
else:
  scalapack=subprocess.check_output(["brew","--prefix","scalapack"],text=True).strip()
  openblas=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip()
  libraries=[f"-L{scalapack}/lib","-lscalapack",f"-L{openblas}/lib","-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-scf-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_scf"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-DUSE_SCALAPACK","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dg_hybrid_generalized_eigensystem.f90"),str(root/"src/gs/dc/dg_hybrid_block_cg.f90"),
    str(root/"src/gs/dc/dg_hybrid_scf.f90"),str(root/"tests/dg/test_dg_hybrid_scf_mpi.f90"),*libraries,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid SCF on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_SCF ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid SCF on 1, 2, 4, and 8 ranks")
