#!/usr/bin/env python3
from pathlib import Path
import os,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-fragment-stream-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_fragment_stream"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-std=f2008","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/dc/dg_hybrid_fragment_basis.f90"),
    str(root/"src/gs/dc/dg_hybrid_fragment_basis_stream.f90"),
    str(root/"tests/dg/test_dg_hybrid_fragment_basis_stream_mpi.f90"),"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (2,4):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid fragment basis stream on {nrank} ranks" in run.stdout
print("PASS hybrid fragment basis stream on 2 and 4 ranks")
