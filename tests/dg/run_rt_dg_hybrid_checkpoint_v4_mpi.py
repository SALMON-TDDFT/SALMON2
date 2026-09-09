#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
source=root/"src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90"
assert source.exists(),"RED: distributed-native v4 checkpoint module is absent"
text=source.read_text().lower()
for token in ("salmon_hybrid_dg_manifest_v4","salmon_hybrid_dg_rank_shard_v4",
              "initial_occupied_amplitudes","basis_point_offsets","metric_offsets"):
  assert token in text,f"RED: v4 checkpoint is missing {token}"
assert "full_coefficients" not in text and "full_metric" not in text
with tempfile.TemporaryDirectory(prefix="hybrid-v4-checkpoint-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_v4_checkpoint"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",str(source),
    str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_v4_mpi.f90"),"-o",str(exe)],check=True)
  reject=build/"hybrid_v4_reject"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",str(source),
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
