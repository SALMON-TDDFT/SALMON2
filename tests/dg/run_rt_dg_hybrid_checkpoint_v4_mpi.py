#!/usr/bin/env python3
from pathlib import Path
import os, shlex, shutil, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
source=root/"src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90"
endpoint=root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90"
assert source.exists(),"RED: distributed-native v4 checkpoint module is absent"
text=source.read_text().lower()
endpoint_text=endpoint.read_text().lower()
main_text=(root/"src/gs/main_dft.f90").read_text().lower()
publisher=main_text.split("subroutine publish_dg_hybrid_divided_v4",1)[1].split(
  "end subroutine publish_dg_hybrid_divided_v4",1)[0]
for token in ("salmon_hybrid_dg_manifest_v4","salmon_hybrid_dg_rank_shard_v4",
              "initial_occupied_amplitudes","basis_point_offsets","metric_offsets"):
  assert token in text,f"RED: v4 checkpoint is missing {token}"
assert "full_coefficients" not in text and "full_metric" not in text
assert endpoint_text.count("subroutine publish_rt_dg_hybrid_checkpoint_v4") == 2
endpoint_body=endpoint_text.split("subroutine publish_rt_dg_hybrid_checkpoint_v4",1)[1].split(
  "end subroutine publish_rt_dg_hybrid_checkpoint_v4",1)[0]
for token in ("collective_rt_dg_hybrid_publication_precondition",
              "collective_rt_dg_hybrid_publication_mapping_precondition",
              "write_rt_dg_hybrid_checkpoint_v4"):
  assert endpoint_body.count(token)==1,token
assert publisher.count("call publish_rt_dg_hybrid_checkpoint_v4")==1
assert "call write_rt_dg_hybrid_checkpoint_v4" not in publisher
for mutation in (
  endpoint_body.replace("write_rt_dg_hybrid_checkpoint_v4","removed_publication",1),
  endpoint_body.replace("occupied_row_ids,local_valid","row_ids,local_valid",1),
):
  assert mutation != endpoint_body
with tempfile.TemporaryDirectory(prefix="hybrid-v4-checkpoint-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_v4_checkpoint"
  support=[root/"src/common/dg_hybrid_sparse_metric.f90",root/"src/common/dg_hybrid_sparse_operators.f90",
    root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90",root/"src/rt/dg/rt_dg_hybrid_point_density.f90",
    source,endpoint,root/"src/rt/dg/rt_dg_hybrid_structural_graph.f90",
    root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90",root/"src/rt/dg/rt_dg_hybrid_initialization_v4.f90"]
  if os.environ.get("SALMON_LAPACK_LIBS"):
    libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
  elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
    prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
  else: libs=["-llapack","-lblas"]
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",*[str(path) for path in support],
    str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_v4_mpi.f90"),*libs,"-o",str(exe)],check=True)
  reject=build/"hybrid_v4_reject"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",str(source),str(endpoint),
    str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_v4_reject_mpi.f90"),"-o",str(reject)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env,timeout=60)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS distributed-v4 shard manifest ranks={nrank}" in run.stdout
    prefix=Path(f"/tmp/salmon-hybrid-v4-checkpoint-{nrank}")
    if nrank==2:
      shards=sorted(prefix.parent.glob(prefix.name+".v4.*.rank000001.shard"),key=lambda p:p.stat().st_mtime_ns)
      assert shards
      damaged=bytearray(shards[-1].read_bytes());damaged[-1]^=0x01;shards[-1].write_bytes(damaged)
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "PASS v4 collective rejection ranks=2 diagnostic=" in rejected.stdout
      assert "rank shard is partial, stale, or corrupt" in rejected.stdout.lower(),rejected.stdout
    if nrank==4:
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "MPI rank mapping changed" in rejected.stdout
print("PASS distributed-v4 rank shards and atomic manifest on 1, 2, 4, and 8 ranks")
