#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,struct,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
checkpoint_source=(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text().lower()
writer=checkpoint_source.split("subroutine write_rt_dg_hybrid_checkpoint",1)[1].split(
  "end subroutine write_rt_dg_hybrid_checkpoint",1)[0]
assert "checkpoint_version=2" in checkpoint_source
assert "operators%metric_values" not in writer
assert "operator_metric" not in writer
if os.environ.get("SALMON_LAPACK_LIBS"):
  lapack_libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("pkg-config") and subprocess.run(["pkg-config","--exists","openblas"],check=False).returncode==0:
  lapack_libs=shlex.split(subprocess.check_output(["pkg-config","--libs","openblas"],text=True))
elif shutil.which("brew"):
  probe=subprocess.run(["brew","--prefix","openblas"],capture_output=True,text=True)
  lapack_libs=[f"-L{probe.stdout.strip()}/lib","-lopenblas"] if probe.returncode==0 else ["-llapack","-lblas"]
else:
  lapack_libs=["-llapack","-lblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-checkpoint-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_checkpoint";occupied_exe=build/"hybrid_occupied_checkpoint";checkpoint=build/"state.chk"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),str(root/"src/common/dg_hybrid_sparse_operators.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90"),str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90"),
    *lapack_libs,"-o",str(exe)],check=True)
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),str(root/"src/common/dg_hybrid_sparse_operators.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90"),str(root/"tests/dg/test_rt_dg_hybrid_occupied_checkpoint_mpi.f90"),
    *lapack_libs,"-o",str(occupied_exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  legacy_checkpoint=build/"legacy-v1.chk"
  write_legacy=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"write_legacy",str(legacy_checkpoint)],capture_output=True,text=True,env=env)
  assert write_legacy.returncode==0,(write_legacy.stdout,write_legacy.stderr);legacy_fingerprints=[]
  for nrank in (1,2,4,8):
    legacy=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read",str(legacy_checkpoint)],capture_output=True,text=True,env=env)
    assert legacy.returncode==0,(nrank,legacy.stdout,legacy.stderr);assert f"PASS hybrid checkpoint on {nrank} ranks" in legacy.stdout
    match=re.search(r"HYBRID_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",legacy.stdout);assert match,legacy.stdout;legacy_fingerprints.append(int(match.group(1)))
  assert len(set(legacy_fingerprints))==1,legacy_fingerprints
  write=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"write",str(checkpoint)],capture_output=True,text=True,env=env)
  assert write.returncode==0,(write.stdout,write.stderr);fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read",str(checkpoint)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr);assert f"PASS hybrid checkpoint on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout;fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
  for mode in ("read_stale","read_stale_state","read_stale_selection","read_stale_window","read_stale_packet",
               "read_stale_complement","read_stale_position","read_stale_operator"):
    stale=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),mode,str(checkpoint)],capture_output=True,text=True,env=env)
    assert stale.returncode==0,(mode,stale.stdout,stale.stderr)
  rank_stale=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"read_rank_stale",str(checkpoint)],capture_output=True,text=True,env=env)
  assert rank_stale.returncode==0,(rank_stale.stdout,rank_stale.stderr)
  incomplete=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"write_incomplete",str(build/"incomplete.chk")],capture_output=True,text=True,env=env)
  assert incomplete.returncode==0,(incomplete.stdout,incomplete.stderr)
  corrupt=build/"corrupt.chk";payload=bytearray(checkpoint.read_bytes());payload[-1]^=0x5A;corrupt.write_bytes(payload)
  bad=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"read_corrupt",str(corrupt)],capture_output=True,text=True,env=env)
  assert bad.returncode==0,(bad.stdout,bad.stderr)
  for label,data in (("truncated_header",checkpoint.read_bytes()[:24]),
                     ("truncated_payload",checkpoint.read_bytes()[:len(checkpoint.read_bytes())//2])):
    broken=build/f"{label}.chk";broken.write_bytes(data)
    rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),"read_corrupt",str(broken)],capture_output=True,text=True,env=env)
    assert rejected.returncode==0,(label,rejected.stdout,rejected.stderr)
  occupied=build/"occupied.chk"
  write_occupied=subprocess.run([shutil.which("mpiexec"),"-n","2",str(occupied_exe),"write",str(occupied)],capture_output=True,text=True,env=env)
  assert write_occupied.returncode==0,(write_occupied.stdout,write_occupied.stderr);occupied_fingerprints=[]
  for nrank in (1,2,4,8):
    read_occupied=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(occupied_exe),"read",str(occupied)],capture_output=True,text=True,env=env)
    assert read_occupied.returncode==0,(nrank,read_occupied.stdout,read_occupied.stderr)
    match=re.search(r"HYBRID_OCCUPIED_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",read_occupied.stdout);assert match,read_occupied.stdout
    occupied_fingerprints.append(int(match.group(1)))
  assert len(set(occupied_fingerprints))==1,occupied_fingerprints
  for mode,target in (("stale",occupied),("stale_provenance",occupied),("changed_occupation",occupied),("old",checkpoint)):
    rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(occupied_exe),mode,str(target)],capture_output=True,text=True,env=env)
    assert rejected.returncode==0,(mode,rejected.stdout,rejected.stderr)
  for mode in ("write_incomplete","write_bad_occupation","write_stale_scf"):
    rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(occupied_exe),mode,str(build/f"{mode}.chk")],capture_output=True,text=True,env=env)
    assert rejected.returncode==0,(mode,rejected.stdout,rejected.stderr)
  hostile=build/"occupied_hostile_extent.chk";hostile_payload=bytearray(occupied.read_bytes())
  struct.pack_into("=i",hostile_payload,24,1_000_000_000);hostile.write_bytes(hostile_payload[:128])
  rejected=subprocess.run([shutil.which("mpiexec"),"-n","2",str(occupied_exe),"old",str(hostile)],capture_output=True,text=True,env=env)
  assert rejected.returncode==0,(rejected.stdout,rejected.stderr)
print("PASS hybrid checkpoint on 1, 2, 4, and 8 ranks")
