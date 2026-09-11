#!/usr/bin/env python3
from pathlib import Path
import os, shlex, shutil, struct, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
source=root/"src/rt/dg/rt_dg_hybrid_checkpoint_v5.f90"
endpoint=source
assert source.exists(),"RED: distributed-native v5 checkpoint module is absent"
text=source.read_text().lower()
endpoint_text=endpoint.read_text().lower()
main_text=(root/"src/gs/main_dft.f90").read_text().lower()
publisher=main_text.split("subroutine publish_dg_hybrid_divided_v5",1)[1].split(
  "end subroutine publish_dg_hybrid_divided_v5",1)[0]
for token in ("salmon_hybrid_dg_manifest_v5","salmon_hybrid_dg_rank_shard_v5",
              "initial_occupied_amplitudes","basis_point_offsets","metric_offsets",
              "system_fingerprint","pseudopotential_fingerprint","pseudopotential_digest"):
  assert token in text,f"RED: v5 checkpoint is missing {token}"
assert "full_coefficients" not in text and "full_metric" not in text
reader=text.split("subroutine read_rt_dg_hybrid_checkpoint_v5",1)[1].split(
  "end subroutine read_rt_dg_hybrid_checkpoint_v5",1)[0]
assert "valid_read_dimensions" in reader,"RED: shard dimensions are allocated before checked extent validation"
assert "unit=-1" in reader,"RED: checkpoint reader may close an uninitialized unit"
writer=text.split("subroutine write_rt_dg_hybrid_checkpoint_v5",1)[1].split(
  "end subroutine write_rt_dg_hybrid_checkpoint_v5",1)[0]
assert writer.count("unit=-1") >= 2,"RED: writer does not reset its unit before both OPEN operations"
assert "close_if_open" in writer,"RED: writer closes a failed/uninitialized unit"
collective_diagnostics=("transaction broadcast failed","shard-size gather failed","shard-digest gather failed",
                        "fragment-map gather failed","manifest metadata")
collective_diagnostic_counts={diagnostic:text.count(diagnostic) for diagnostic in collective_diagnostics}
def require_collective_safety(candidate: str) -> None:
  assert ";call mpi_" not in candidate, \
    "RED: v5 transaction starts another MPI phase before checking the preceding collective"
  assert "ierr/=mpi_success.or." not in candidate, \
    "RED: v5 I/O reads an undefined collective output when MPI returns an error"
  for diagnostic in collective_diagnostics:
    assert candidate.count(diagnostic)>=collective_diagnostic_counts[diagnostic], \
      f"RED: missing collective-safe v5 diagnostic: {diagnostic}"
require_collective_safety(text)
for diagnostic in collective_diagnostics:
  mutation=text.replace(diagnostic,"removed collective diagnostic",1)
  try: require_collective_safety(mutation)
  except AssertionError: pass
  else: raise AssertionError(f"collective guard mutation survived: {diagnostic}")
assert "expected_shard_extent" in text and "valid_read_dimensions" in reader, \
  "RED: reader lacks exact overflow-safe serialized extent validation"
assert "manifest_size/=" in reader and "actual_size<" in reader,"RED: fixed headers are read without exact/minimum-size checks"
assert "certified_rank<nocc" in text,"RED: v5 checkpoint permits certification below occupied rank"
for field in ("global_grid_count","certified_rank","operator_structure_fingerprint","scope_fingerprint","payload_fingerprint",
              "system_fingerprint","pseudopotential_fingerprint"):
  assert field in reader,f"RED: shard {field} is not compared with the manifest"
assert "dg_sha256_final(hash,digest)" in text,"RED: checkpoint does not retain a full streaming SHA-256 digest"
assert "shard_digests(4,nproc)" in text,"RED: manifest truncates shard SHA-256"
assert endpoint_text.count("subroutine publish_rt_dg_hybrid_checkpoint_v5") == 2
endpoint_body=endpoint_text.split("subroutine publish_rt_dg_hybrid_checkpoint_v5",1)[1].split(
  "end subroutine publish_rt_dg_hybrid_checkpoint_v5",1)[0]
assert endpoint_body.count("collective_rt_dg_hybrid_publication_precondition")==2
for token in ("collective_rt_dg_hybrid_publication_mapping_precondition",
              "write_rt_dg_hybrid_checkpoint_v5"):
  assert endpoint_body.count(token)==1,token
assert "authorization%valid" in endpoint_body and "endpoint authorization failed" in endpoint_body
for field in ("checkpoint_version","published_rank","basis_fingerprint","operator_fingerprint"):
  assert f"authorization%{field}" in endpoint_body
assert publisher.count("call publish_rt_dg_hybrid_checkpoint_v5")==1
assert "call write_rt_dg_hybrid_checkpoint_v5" not in publisher
for mutation in (
  endpoint_body.replace("write_rt_dg_hybrid_checkpoint_v5","removed_publication",1),
  endpoint_body.replace("occupied_row_ids,local_valid","row_ids,local_valid",1),
):
  assert mutation != endpoint_body
with tempfile.TemporaryDirectory(prefix="hybrid-v5-checkpoint-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_v5_checkpoint"
  support=[root/"src/common/dg_portable_sha256.f90",root/"src/common/dg_hybrid_sparse_metric.f90",root/"src/common/dg_hybrid_sparse_operators.f90",
    root/"src/gs/dc/dg_hybrid_publication_policy.f90",
    root/"tests/dg/legacy_support/dg_hybrid_continuation_controller.f90",
    root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90",root/"src/rt/dg/rt_dg_hybrid_point_density.f90",
    source,root/"src/rt/dg/rt_dg_hybrid_structural_graph.f90",
    root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90",root/"src/rt/dg/rt_dg_hybrid_initialization_v5.f90"]
  if os.environ.get("SALMON_LAPACK_LIBS"):
    libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
  elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
    prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
  else: libs=["-llapack","-lblas"]
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",*[str(path) for path in support],
    str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_v5_mpi.f90"),*libs,"-o",str(exe)],check=True)
  reject=build/"hybrid_v5_reject"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",str(root/"src/common/dg_portable_sha256.f90"),str(source),
    str(root/"tests/dg/test_rt_dg_hybrid_v5_reader_legacy_v4_reject_mpi.f90"),"-o",str(reject)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (1,2,4,8):
    shutil.rmtree(Path(f"/tmp/salmon-hybrid-v5-open-failure-{nrank}"),ignore_errors=True)
    manifest_open_failure=Path(f"/tmp/salmon-hybrid-v5-manifest-open-failure-{nrank}.manifest.temporary")
    if manifest_open_failure.is_dir(): shutil.rmtree(manifest_open_failure)
    elif manifest_open_failure.exists(): manifest_open_failure.unlink()
    manifest_open_failure.mkdir();(manifest_open_failure/"open-must-fail").write_text("sentinel")
    open_failure_env=env.copy();open_failure_env["SALMON_TEST_V5_MANIFEST_OPEN_FAILURE"]="1"
    failed_open=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],
      capture_output=True,text=True,env=open_failure_env,timeout=60)
    assert failed_open.returncode==0,(nrank,failed_open.stdout,failed_open.stderr)
    assert f"PASS v5 collective manifest OPEN failure ranks={nrank}" in failed_open.stdout
    shutil.rmtree(manifest_open_failure)
    manifest_failure=Path(f"/tmp/salmon-hybrid-v5-manifest-failure-{nrank}.manifest")
    if manifest_failure.is_dir(): shutil.rmtree(manifest_failure)
    elif manifest_failure.exists(): manifest_failure.unlink()
    manifest_failure.mkdir();(manifest_failure/"rename-must-fail").write_text("sentinel")
    failure_env=env.copy();failure_env["SALMON_TEST_V5_MANIFEST_FAILURE"]="1"
    failed_manifest=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],
      capture_output=True,text=True,env=failure_env,timeout=60)
    assert failed_manifest.returncode==0,(nrank,failed_manifest.stdout,failed_manifest.stderr)
    assert f"PASS v5 collective manifest failure ranks={nrank}" in failed_manifest.stdout
    shutil.rmtree(manifest_failure)
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env,timeout=60)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS distributed-v5 shard manifest ranks={nrank}" in run.stdout
    old_prefix=Path(f"/tmp/salmon-hybrid-old-v4-{nrank}")
    old_manifest=bytearray(120+20*nrank)
    old_manifest[:32]=b"SALMON_HYBRID_DG_MANIFEST_V4".ljust(32,b" ")
    struct.pack_into("=i",old_manifest,32,4)
    Path(str(old_prefix)+".manifest").write_bytes(old_manifest)
    old_reject=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(reject),str(old_prefix)],
      capture_output=True,text=True,env=env,timeout=30)
    assert old_reject.returncode==0,(nrank,old_reject.stdout,old_reject.stderr)
    assert "unsupported distributed checkpoint schema v4" in old_reject.stdout.lower(),old_reject.stdout
    Path(str(old_prefix)+".manifest").unlink()
    corrupt_prefix=Path(f"/tmp/salmon-hybrid-corrupt-v5-{nrank}")
    corrupt_manifest=bytearray(176+44*nrank)
    corrupt_manifest[:32]=b"SALMON_HYBRID_DG_MANIFEST_BAD".ljust(32,b" ")
    struct.pack_into("=i",corrupt_manifest,32,5)
    Path(str(corrupt_prefix)+".manifest").write_bytes(corrupt_manifest)
    corrupt_reject=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(reject),str(corrupt_prefix)],
      capture_output=True,text=True,env=env,timeout=30)
    assert corrupt_reject.returncode==0,(nrank,corrupt_reject.stdout,corrupt_reject.stderr)
    assert "distributed-v5 manifest magic/schema/extent is corrupt" in corrupt_reject.stdout.lower(),corrupt_reject.stdout
    Path(str(corrupt_prefix)+".manifest").unlink()
    prefix=Path(f"/tmp/salmon-hybrid-v5-checkpoint-{nrank}")
    if nrank==1:
      shard=sorted(prefix.parent.glob(prefix.name+".v5.*.rank000000.shard"),key=lambda p:p.stat().st_mtime_ns)[-1]
      original=shard.read_bytes()
      for label,damaged,diagnostic in (
        ("truncated",original[:100],"truncated or has an invalid fixed header"),
        ("negative",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("huge",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("operator-count",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("point-count",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("coefficient-product",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("large-consistent-product",bytearray(original),"negative, overflowing, or invalid dimensions"),
        ("certified-below-occupied",bytearray(original),"disagrees with manifest common metadata"),
      ):
        if label=="negative": struct.pack_into("=i",damaged,216,-1)
        if label=="huge": struct.pack_into("=i",damaged,216,2**31-1)
        if label=="operator-count": struct.pack_into("=i",damaged,232,struct.unpack_from("=i",damaged,232)[0]+1)
        if label=="point-count": struct.pack_into("=i",damaged,236,struct.unpack_from("=i",damaged,236)[0]+1)
        if label=="coefficient-product": struct.pack_into("=i",damaged,248,2**31-1)
        if label=="large-consistent-product":
          for offset,value in ((216,10**9),(220,10**9+1),(228,10**9+1),(244,10**9)):
            struct.pack_into("=i",damaged,offset,value)
        if label=="certified-below-occupied":
          nocc=struct.unpack_from("=i",damaged,56)[0];struct.pack_into("=i",damaged,60,nocc-1)
        shard.write_bytes(damaged)
        rejected=subprocess.run([shutil.which("mpiexec"),"-n","1",str(reject),str(prefix)],
          capture_output=True,text=True,env=env,timeout=30)
        assert rejected.returncode==0,(label,rejected.stdout,rejected.stderr)
        assert diagnostic in rejected.stdout.lower(),(label,rejected.stdout)
        shard.write_bytes(original)
    if nrank==2:
      manifest=Path(str(prefix)+".manifest")
      original_manifest=manifest.read_bytes();damaged_manifest=bytearray(original_manifest)
      manifest_grid=struct.unpack_from("=i",damaged_manifest,44)[0]
      struct.pack_into("=i",damaged_manifest,44,manifest_grid+1)
      manifest.write_bytes(damaged_manifest)
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "disagrees with manifest common metadata" in rejected.stdout.lower(),rejected.stdout
      manifest.write_bytes(original_manifest)
      damaged_manifest=bytearray(original_manifest)
      manifest_nocc=struct.unpack_from("=i",damaged_manifest,48)[0]
      struct.pack_into("=i",damaged_manifest,52,manifest_nocc-1)
      manifest.write_bytes(damaged_manifest)
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "disagrees with manifest common metadata" in rejected.stdout.lower(),rejected.stdout
      manifest.write_bytes(original_manifest)
      manifest.write_bytes(original_manifest+b"X")
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "manifest magic/schema/extent is corrupt" in rejected.stdout.lower(),rejected.stdout
      manifest.write_bytes(original_manifest)
      shards=sorted(prefix.parent.glob(prefix.name+".v5.*.rank000001.shard"),key=lambda p:p.stat().st_mtime_ns)
      assert shards
      original_shard=shards[-1].read_bytes();damaged=bytearray(original_shard);damaged[-1]^=0x01;shards[-1].write_bytes(damaged)
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "PASS v5 collective rejection ranks=2 diagnostic=" in rejected.stdout
      assert "rank shard is partial, stale, or corrupt" in rejected.stdout.lower(),rejected.stdout
      shards[-1].write_bytes(original_shard)
    if nrank==4:
      rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(reject),str(prefix)],
        capture_output=True,text=True,env=env,timeout=30)
      assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
      assert "MPI rank mapping changed" in rejected.stdout
print("PASS distributed-v5 rank shards and atomic manifest on 1, 2, 4, and 8 ranks")
