#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,struct,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
checkpoint_source=(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text().lower()
main_source=(root/"src/gs/main_dft.f90").read_text().lower()
writer=checkpoint_source.split("subroutine write_rt_dg_hybrid_checkpoint",1)[1].split(
  "end subroutine write_rt_dg_hybrid_checkpoint",1)[0]
assert "checkpoint_version=2" in checkpoint_source
assert "operators%metric_values" not in writer
assert "operator_metric" not in writer
assert "write_rt_dg_hybrid_ground_state_checkpoint" in checkpoint_source
assert "read_rt_dg_hybrid_ground_state_checkpoint" in checkpoint_source
continuation=main_source.split("subroutine run_dg_hybrid_concrete_continuation",1)[1].split(
  "end subroutine run_dg_hybrid_concrete_continuation",1)[0]
assert "write_rt_dg_hybrid_ground_state_checkpoint" in continuation
assert "hybrid_dg_ground_state.chk" in continuation
for component in ("energy%e_tot","energy%e_kin","energy%e_h","energy%e_xc","energy%e_ion_ion",
                  "energy%e_ion_loc","energy%e_ion_nloc"):
  assert component in continuation, f"complete checkpoint omits energy provenance: {component}"
for component in ("pp%zion","pp%lmax","pp%nrmax","ppg%nlma"):
  assert component in continuation, f"complete checkpoint omits pseudopotential provenance: {component}"
assert continuation.index("if(.not.final_refresh_performed)") < continuation.index(
  "write_rt_dg_hybrid_ground_state_checkpoint")
assert "write_rt_dg_hybrid_occupied_checkpoint" not in continuation
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
  for nrank in (1,2,4):
    complete=build/f"complete-{nrank}.chk"
    complete_write=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"write_complete",str(complete)],capture_output=True,text=True,env=env)
    assert complete_write.returncode==0,(nrank,complete_write.stdout,complete_write.stderr)
    accepted_bytes=complete.read_bytes()
    complete_read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read_complete",str(complete)],capture_output=True,text=True,env=env)
    assert complete_read.returncode==0,(nrank,complete_read.stdout,complete_read.stderr)
    assert f"PASS complete hybrid checkpoint on {nrank} ranks" in complete_read.stdout
    invalid=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"write_bad_complete",str(complete)],capture_output=True,text=True,env=env)
    assert invalid.returncode==0,(nrank,invalid.stdout,invalid.stderr)
    assert complete.read_bytes()==accepted_bytes,"rejected complete checkpoint replaced the accepted file"
    incomplete=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"write_incomplete_complete",str(complete)],capture_output=True,text=True,env=env)
    assert incomplete.returncode==0,(nrank,incomplete.stdout,incomplete.stderr)
    assert complete.read_bytes()==accepted_bytes,"incomplete complete checkpoint replaced the accepted file"
    for mode in ("write_bad_grid_complete","write_out_of_range_grid_complete","write_missing_position_convention_complete"):
      invalid_catalog=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),mode,str(complete)],capture_output=True,text=True,env=env)
      assert invalid_catalog.returncode==0,(nrank,mode,invalid_catalog.stdout,invalid_catalog.stderr)
      assert complete.read_bytes()==accepted_bytes,f"{mode} replaced the accepted file"
    interrupted=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"write_interrupted_complete",str(complete)],capture_output=True,text=True,env=env)
    assert interrupted.returncode==0,(nrank,interrupted.stdout,interrupted.stderr)
    assert complete.read_bytes()==accepted_bytes,"interrupted complete checkpoint replaced the accepted file"
    for label,fraction in (("metadata",0.08),("basis",0.45),("matrix",0.67),("state",0.90)):
      corrupt_complete=build/f"complete-{nrank}-{label}-corrupt.chk";corrupt_payload=bytearray(accepted_bytes)
      corrupt_payload[int(len(corrupt_payload)*fraction)]^=0x31;corrupt_complete.write_bytes(corrupt_payload)
      corrupt_read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read_complete_corrupt",str(corrupt_complete)],capture_output=True,text=True,env=env)
      assert corrupt_read.returncode==0,(nrank,label,corrupt_read.stdout,corrupt_read.stderr)
    for label,offset,value,kind in (("operation-count",48,3,"i"),("identity-flag",36,1,"i"),
                                    ("kinetic-fingerprint",96,999999,"q"),
                                    ("payload-fingerprint",200,999999,"q")):
      corrupt_complete=build/f"complete-{nrank}-{label}.chk";corrupt_payload=bytearray(accepted_bytes)
      struct.pack_into(f"={kind}",corrupt_payload,offset,value);corrupt_complete.write_bytes(corrupt_payload)
      corrupt_read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read_complete_corrupt",str(corrupt_complete)],capture_output=True,text=True,env=env)
      assert corrupt_read.returncode==0,(nrank,label,corrupt_read.stdout,corrupt_read.stderr)
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
