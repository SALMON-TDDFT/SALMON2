#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
source=(root/"src/gs/dc/dg_hybrid_generalized_eigensystem.f90").read_text().lower()
backend=source[source.index("subroutine solve_dg_hybrid_generalized_scalapack"):]
backend=backend[:backend.index("end subroutine solve_dg_hybrid_generalized_scalapack")]
assert "call mpi_allreduce(local_dot,global_dot,1" not in backend,(
  "complete eigensystem diagnostics retain an n-state squared scalar reduction loop")
assert "call mpi_bcast(remote_row,nstate" not in backend,(
  "complete eigensystem fingerprint retains an n-state squared vector broadcast loop")
complete=source[source.index("subroutine solve_dg_hybrid_generalized_complete_once"):]
complete=complete[:complete.index("end subroutine solve_dg_hybrid_generalized_complete_once")]
normalized=re.sub(r"\s+","",complete)
assert complete.count("call solver(")==1,"complete eigensystem wrapper must invoke its backend exactly once"
assert "global_count,global_count" in normalized,"complete eigensystem wrapper did not request the full basis"
if shutil.which("pkg-config"):
  scalapack=shlex.split(subprocess.check_output(["pkg-config","--libs","scalapack"],text=True))
  openblas=shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
else:
  scalapack_prefix=subprocess.check_output(["brew","--prefix","scalapack"],text=True).strip()
  openblas_prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip()
  scalapack=[f"-L{scalapack_prefix}/lib","-lscalapack"]
  openblas=[f"-L{openblas_prefix}/lib","-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-generalized-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_generalized"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-DUSE_SCALAPACK","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/occupation_kernel.f90"),
    str(root/"src/gs/dc/dg_hybrid_occupation_policy.f90"),
    str(root/"src/gs/dc/dg_hybrid_ground_state_types.f90"),
    str(root/"src/gs/dc/dg_hybrid_generalized_eigensystem.f90"),
    str(root/"tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90"),*scalapack,*openblas,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid generalized eigensystem on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_GENERALIZED_EIGENSYSTEM ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid generalized eigensystem on 1, 2, 4, and 8 ranks")
