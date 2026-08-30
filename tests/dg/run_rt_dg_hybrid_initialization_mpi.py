#!/usr/bin/env python3
from pathlib import Path
import os,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
initialization_source=(root/"src/rt/dg/rt_dg_hybrid_initialization.f90").read_text().lower()
initializer=initialization_source.split("subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[1].split(
  "end subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[0]
for forbidden in ("allocate(metric(n,n)","position_rows(3,n,n)","metric_graph(n,n)","operator_graph(n,n)"):
  assert forbidden not in initializer, f"RT initialization still replicates a dense global object: {forbidden}"
assert "redistribute_ground_state_rows" in initializer
if os.environ.get("SALMON_LAPACK_LIBS"): libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
else: libs=["-llapack","-lblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-rt-init-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_rt_init"
  sources=["src/common/dg_hybrid_sparse_metric.f90","src/common/dg_hybrid_sparse_operators.f90",
    "src/rt/dg/rt_dg_hybrid_checkpoint.f90","src/rt/dg/rt_dg_hybrid_initialization.f90",
    "src/rt/dg/rt_dg_hybrid_density_update.f90",
    "tests/dg/test_rt_dg_hybrid_initialization_mpi.f90"]
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",*[str(root/s) for s in sources],*libs,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (1,2,4):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),str(build/f"state-{nrank}.chk")],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid RT initialization on {nrank} ranks" in run.stdout
  cross=build/"cross-rank.chk"
  write=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),str(cross),"write_only"],capture_output=True,text=True,env=env)
  assert write.returncode==0,(write.stdout,write.stderr)
  read=subprocess.run([shutil.which("mpiexec"),"-n","4",str(exe),str(cross),"read_only"],capture_output=True,text=True,env=env)
  assert read.returncode==0,(read.stdout,read.stderr)
  write_fp=[line for line in write.stdout.splitlines() if line.startswith("HYBRID_RT_INIT_FINGERPRINT=")][0]
  read_fp=[line for line in read.stdout.splitlines() if line.startswith("HYBRID_RT_INIT_FINGERPRINT=")][0]
  assert write_fp==read_fp,(write_fp,read_fp)
  coalesced=build/"coalesced-rank.chk"
  write4=subprocess.run([shutil.which("mpiexec"),"-n","4",str(exe),str(coalesced),"write_only"],capture_output=True,text=True,env=env)
  assert write4.returncode==0,(write4.stdout,write4.stderr)
  write4_fp=[line for line in write4.stdout.splitlines() if line.startswith("HYBRID_RT_INIT_FINGERPRINT=")][0]
  for nrank in (2,1):
    read_small=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),str(coalesced),"read_only"],capture_output=True,text=True,env=env)
    assert read_small.returncode==0,(nrank,read_small.stdout,read_small.stderr)
    read_small_fp=[line for line in read_small.stdout.splitlines() if line.startswith("HYBRID_RT_INIT_FINGERPRINT=")][0]
    assert write4_fp==read_small_fp,(write4_fp,read_small_fp)
print("PASS hybrid RT initialization on 1, 2, and 4 ranks")
