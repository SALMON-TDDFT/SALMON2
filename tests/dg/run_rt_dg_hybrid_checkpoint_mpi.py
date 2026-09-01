#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,struct,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
checkpoint_source=(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text().lower()
main_source=(root/"src/gs/main_dft.f90").read_text().lower()
coalesced_redistribution=checkpoint_source.split("subroutine redistribute_coalesced_ground_state",1)[1].split(
  "end subroutine redistribute_coalesced_ground_state",1)[0]
ownership_validation=checkpoint_source.split("subroutine validate_ground_state_global_ownership",1)[1].split(
  "end subroutine validate_ground_state_global_ownership",1)[0]
ground_state_writer=checkpoint_source.split("subroutine write_rt_dg_hybrid_ground_state_checkpoint",1)[1].split(
  "end subroutine write_rt_dg_hybrid_ground_state_checkpoint",1)[0]
stream_readers="\n".join(re.findall(
  r"subroutine stream_read_[\s\S]*?end subroutine(?:\s+stream_read_\w+)?",checkpoint_source))
stream_writers="\n".join(re.findall(
  r"subroutine stream_write_[\s\S]*?end subroutine(?:\s+stream_write_\w+)?",checkpoint_source))
assert not re.search(r"call\s+mpi_bcast\s*\(\s*buffer\s*,",stream_readers),(
  "checkpoint reader broadcasts each dense shard to every reader rank")
assert not re.search(r"call\s+mpi_bcast\s*\(\s*buffer\s*,",stream_writers),(
  "checkpoint writer broadcasts each dense shard to every writer rank")
for redistributed_buffer in (
  "metric_columns_buffer","operator_columns_buffer","row_components","coefficient_buffer",
  "position_buffer","c_buffer","b_buffer","grid_integer_buffer","grid_real_buffer",
  "grid_basis_buffer","rt_grid_basis_buffer","u_buffer","rt_components","rt_scalar_buffer",
  "rt_vector_buffer","rt_tensor_buffer","point_ids_buffer","basis_ids_buffer","metadata_buffer",
  "normal_buffer","face_weight_buffer","face_weights_buffer","face_values_buffer","interface_buffer",
  "record_owner_buffer","nonlocal_values_buffer",
):
  assert not re.search(rf"call\s+mpi_bcast\s*\(\s*{redistributed_buffer}\s*,",coalesced_redistribution),(
    f"checkpoint redistribution broadcasts {redistributed_buffer} to every reader rank")
assert "call mpi_send" in coalesced_redistribution and "call mpi_recv" in coalesced_redistribution,(
  "checkpoint redistribution does not transfer records directly to their canonical owner")
assert "mpi_allgatherv" not in coalesced_redistribution,(
  "checkpoint redistribution replicates a distributed record catalog on every reader rank")
assert "subroutine sort_face_records" not in coalesced_redistribution and "subroutine sort_id_owners" not in (
  coalesced_redistribution),"checkpoint redistribution retains quadratic insertion sorting"
assert coalesced_redistribution.count("call mpi_alltoallv")>=10,(
  "checkpoint redistribution does not route variable face/nonlocal payloads directly by owner")
assert "subroutine scale_counts_checked" in coalesced_redistribution and "send_data_counts64" in (
  coalesced_redistribution),"checkpoint redistribution does not guard classic-MPI count overflow"
assert "mpi_allgatherv(local_ids" not in ownership_validation,(
  "checkpoint ownership validation replicates every face/nonlocal ID on every rank")
assert "call mpi_alltoallv" in ownership_validation,(
  "checkpoint ownership validation does not distribute ID uniqueness buckets")
assert "mpi_allgatherv" not in ground_state_writer,(
  "ground-state writer still replicates a distributed ownership catalog")
writer=checkpoint_source.split("subroutine write_rt_dg_hybrid_checkpoint",1)[1].split(
  "end subroutine write_rt_dg_hybrid_checkpoint",1)[0]
assert re.search(r"checkpoint_version\s*=\s*2",checkpoint_source)
assert re.search(r"occupied_version\s*=\s*2",checkpoint_source)
assert re.search(r"ground_state_version\s*=\s*3",checkpoint_source)
assert "operators%metric_values" not in writer
assert "operator_metric" not in writer
assert "write_rt_dg_hybrid_ground_state_checkpoint" in checkpoint_source
assert "read_rt_dg_hybrid_ground_state_checkpoint" in checkpoint_source
for named_contract in (
  "type,public::s_rt_dg_hybrid_construction_catalog",
  "type,public::s_rt_dg_hybrid_certified_basis_payload",
  "type,public::s_rt_dg_hybrid_electron_count_receipt",
  "type,public::s_rt_dg_hybrid_rt_space_payload",
  "type,public::s_rt_dg_hybrid_energy_window_receipt",
  "type,public::s_rt_dg_hybrid_symmetry_receipt",
  "type,public::s_rt_dg_hybrid_handoff_receipts",
  "authenticate_rt_dg_hybrid_ground_state_payload",
):
  assert named_contract in checkpoint_source, f"missing v3 checkpoint contract: {named_contract}"
continuation=main_source.split("subroutine run_dg_hybrid_concrete_continuation",1)[1].split(
  "end subroutine run_dg_hybrid_concrete_continuation",1)[0]
assert "write_rt_dg_hybrid_ground_state_checkpoint" in continuation
assert "hybrid_dg_ground_state.chk" in continuation
for component in ("calc_total_energy_periodic(dc%mg_tot,ewald,dc%system_tot","fixed_payload%kinetic_rows",
                  "fixed_payload%nonlocal_rows","eexc_tmp(energy_ix,energy_iy,energy_iz)",
                  "checkpoint_energy%e_ion_ion"):
  assert component in continuation, f"complete checkpoint omits energy provenance: {component}"
assert "checkpoint_payload%energy_receipt=[energy%e_tot" not in continuation
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
  nonmpi_probe=build/"nonmpi_checkpoint_probe.f90";nonmpi_probe.write_text("""program nonmpi_checkpoint_probe
  use rt_dg_hybrid_checkpoint, only: rt_dg_hybrid_ground_state_checkpoint_version
  if (rt_dg_hybrid_ground_state_checkpoint_version /= 3) error stop
end program nonmpi_checkpoint_probe
""")
  subprocess.run([shutil.which("gfortran"),"-cpp","-I",str(build),"-J",str(build),
    str(root/"src/common/dg_hybrid_sparse_metric.f90"),str(root/"src/common/dg_hybrid_sparse_operators.f90"),
    str(root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90"),str(nonmpi_probe),*lapack_libs,
    "-o",str(build/"nonmpi_checkpoint_probe")],check=True)
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
  complete_fingerprints=[]
  for nrank in (1,2,4,8):
    complete=build/f"complete-{nrank}.chk"
    complete_write=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"write_complete",str(complete)],capture_output=True,text=True,env=env)
    assert complete_write.returncode==0,(nrank,complete_write.stdout,complete_write.stderr)
    accepted_bytes=complete.read_bytes()
    assert struct.unpack_from("=i",accepted_bytes,16)[0]==3,"complete GS checkpoint is not version 3"
    match=re.search(r"HYBRID_GS_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",complete_write.stdout)
    assert match,complete_write.stdout;complete_fingerprints.append(int(match.group(1)))
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
    for label,fraction in ((("metadata",0.08),("basis",0.45),("state",0.90)) if nrank==1 else ()):
      corrupt_complete=build/f"complete-{nrank}-{label}-corrupt.chk";corrupt_payload=bytearray(accepted_bytes)
      corrupt_payload[int(len(corrupt_payload)*fraction)]^=0x31;corrupt_complete.write_bytes(corrupt_payload)
      corrupt_read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read_complete_corrupt",str(corrupt_complete)],capture_output=True,text=True,env=env)
      assert corrupt_read.returncode==0,(nrank,label,corrupt_read.stdout,corrupt_read.stderr)
    for label,offset,value,kind in ((("operation-count",96,3,"i"),("identity-flag",36,1,"i"),
                                    ("kinetic-fingerprint",224,999999,"q"),
                                    ("payload-fingerprint",336,999999,"q")) if nrank==1 else ()):
      corrupt_complete=build/f"complete-{nrank}-{label}.chk";corrupt_payload=bytearray(accepted_bytes)
      struct.pack_into(f"={kind}",corrupt_payload,offset,value);corrupt_complete.write_bytes(corrupt_payload)
      corrupt_read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),"read_complete_corrupt",str(corrupt_complete)],capture_output=True,text=True,env=env)
      assert corrupt_read.returncode==0,(nrank,label,corrupt_read.stdout,corrupt_read.stderr)
  assert len(set(complete_fingerprints))==1,complete_fingerprints
  complete1=build/"complete-1.chk"
  between_levels=build/"complete-between-levels.chk"
  between_levels_write=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "write_complete_between_levels",str(between_levels)],capture_output=True,text=True,env=env)
  assert between_levels_write.returncode==0,(between_levels_write.stdout,between_levels_write.stderr)
  signed_spread=build/"complete-signed-spread.chk"
  signed_spread_write=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "write_complete_signed_spread",str(signed_spread)],capture_output=True,text=True,env=env)
  assert signed_spread_write.returncode==0,(signed_spread_write.stdout,signed_spread_write.stderr)
  distinct_selection=build/"complete-distinct-selection.chk"
  distinct_selection_write=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "write_complete_distinct_selection",str(distinct_selection)],capture_output=True,text=True,env=env)
  assert distinct_selection_write.returncode==0,(distinct_selection_write.stdout,distinct_selection_write.stderr)
  bad_proof=build/"complete-bad-proof-below-requested.chk"
  bad_proof_write=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "write_bad_proof_below_requested_complete",str(bad_proof)],capture_output=True,text=True,env=env)
  assert bad_proof_write.returncode==0,(bad_proof_write.stdout,bad_proof_write.stderr)
  legacy_dynamic=build/"complete-legacy-dynamic.chk"
  legacy_dynamic_write=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "write_complete_legacy_dynamic",str(legacy_dynamic)],capture_output=True,text=True,env=env)
  assert legacy_dynamic_write.returncode==0,(legacy_dynamic_write.stdout,legacy_dynamic_write.stderr)
  write_match=re.search(r"HYBRID_GS_LEGACY_DYNAMIC phase=write fingerprint=(-?\d+)",legacy_dynamic_write.stdout)
  assert write_match,legacy_dynamic_write.stdout
  legacy_dynamic_read=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "read_complete_legacy_dynamic",str(legacy_dynamic)],capture_output=True,text=True,env=env)
  assert legacy_dynamic_read.returncode==0,(legacy_dynamic_read.stdout,legacy_dynamic_read.stderr)
  read_match=re.search(r"HYBRID_GS_LEGACY_DYNAMIC phase=read fingerprint=(-?\d+)",legacy_dynamic_read.stdout)
  assert read_match and read_match.group(1)==write_match.group(1),legacy_dynamic_read.stdout
  hostile_shape=build/"complete-1-hostile-shape.chk";hostile_payload=bytearray(complete1.read_bytes())
  first_shard_arrays=24+15*4+25*4+56*8+39*8
  leading_arrays=((4,8),(4,8),(2,8),(10,8),(2,8),(4,4),(5,4),(16,4),(5,4),(4,4))
  first_two_dimensional_array=first_shard_arrays+sum(4+length*element_bytes for length,element_bytes in leading_arrays)
  struct.pack_into("=ii",hostile_payload,first_two_dimensional_array,-2,-2);hostile_shape.write_bytes(hostile_payload)
  rejected_hostile=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),"read_complete_corrupt",str(hostile_shape)],capture_output=True,text=True,env=env)
  assert rejected_hostile.returncode==0,(rejected_hostile.stdout,rejected_hostile.stderr)
  overflow_shape=build/"complete-1-overflow-shape.chk";overflow_payload=bytearray(complete1.read_bytes())
  struct.pack_into("=ii",overflow_payload,first_two_dimensional_array,50000,50000);overflow_shape.write_bytes(overflow_payload)
  rejected_overflow=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),"read_complete_corrupt",str(overflow_shape)],capture_output=True,text=True,env=env)
  assert rejected_overflow.returncode==0,(rejected_overflow.stdout,rejected_overflow.stderr)
  nonfinite_header=build/"complete-1-nonfinite-header.chk";nonfinite_payload=bytearray(complete1.read_bytes())
  first_shard_reals=24+15*4+25*4+56*8
  struct.pack_into("=d",nonfinite_payload,first_shard_reals+(17-1)*8,float("nan"));nonfinite_header.write_bytes(nonfinite_payload)
  rejected_nonfinite=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),"read_complete_corrupt",str(nonfinite_header)],capture_output=True,text=True,env=env)
  assert rejected_nonfinite.returncode==0,(rejected_nonfinite.stdout,rejected_nonfinite.stderr)
  complete8=build/"complete-8.chk"
  first_shard_fingerprints=24+15*4+25*4
  coalesced_common_corrupt=build/"complete-8-common-corrupt.chk"
  coalesced_common_payload=bytearray(complete8.read_bytes())
  struct.pack_into("=q",coalesced_common_payload,first_shard_fingerprints+(22-1)*8,999999)
  coalesced_common_corrupt.write_bytes(coalesced_common_payload)
  rejected_coalesced_common=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "read_complete_coalesced_corrupt",str(coalesced_common_corrupt)],capture_output=True,text=True,env=env)
  assert rejected_coalesced_common.returncode==0,(rejected_coalesced_common.stdout,rejected_coalesced_common.stderr)
  fingerprint_bytes=struct.pack("=q",complete_fingerprints[-1]);checkpoint_bytes=complete8.read_bytes();positions=[];start=0
  while (position:=checkpoint_bytes.find(fingerprint_bytes,start))>=0:
    positions.append(position);start=position+len(fingerprint_bytes)
  assert len(positions)==8,positions
  coalesced_nonfirst_fingerprint=build/"complete-8-nonfirst-payload-fingerprint-corrupt.chk"
  coalesced_nonfirst_payload=bytearray(checkpoint_bytes)
  struct.pack_into("=q",coalesced_nonfirst_payload,positions[1],999999)
  coalesced_nonfirst_fingerprint.write_bytes(coalesced_nonfirst_payload)
  rejected_nonfirst_fingerprint=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),
    "read_complete_coalesced_corrupt",str(coalesced_nonfirst_fingerprint)],capture_output=True,text=True,env=env)
  assert rejected_nonfirst_fingerprint.returncode==0,(rejected_nonfirst_fingerprint.stdout,rejected_nonfirst_fingerprint.stderr)
  coalesced_4_to_2=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),
    "read_complete_coalesced",str(build/"complete-4.chk")],capture_output=True,text=True,env=env)
  assert coalesced_4_to_2.returncode==0,(coalesced_4_to_2.stdout,coalesced_4_to_2.stderr)
  coalesced_4_to_3=subprocess.run([shutil.which("mpiexec"),"-n","3",str(exe),
    "read_complete_coalesced",str(build/"complete-4.chk")],capture_output=True,text=True,env=env)
  assert coalesced_4_to_3.returncode==0,(coalesced_4_to_3.stdout,coalesced_4_to_3.stderr)
  match=re.search(r"HYBRID_GS_CHECKPOINT ranks=3 fingerprint=(-?\d+)",coalesced_4_to_3.stdout)
  assert match and int(match.group(1))==complete_fingerprints[2],coalesced_4_to_3.stdout
  for reader_ranks in (1,2,4):
    coalesced=subprocess.run([shutil.which("mpiexec"),"-n",str(reader_ranks),str(exe),
      "read_complete_coalesced",str(complete8)],capture_output=True,text=True,env=env)
    assert coalesced.returncode==0,(reader_ranks,coalesced.stdout,coalesced.stderr)
    match=re.search(r"HYBRID_GS_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",coalesced.stdout)
    assert match and int(match.group(1))==complete_fingerprints[-1],coalesced.stdout
  for reader_ranks in (1,2,4,8):
    expanded=subprocess.run([shutil.which("mpiexec"),"-n",str(reader_ranks),str(exe),
      "read_complete_coalesced",str(complete1)],capture_output=True,text=True,env=env)
    assert expanded.returncode==0,(reader_ranks,expanded.stdout,expanded.stderr)
    match=re.search(r"HYBRID_GS_CHECKPOINT ranks=\d+ fingerprint=(-?\d+)",expanded.stdout)
    assert match and int(match.group(1))==complete_fingerprints[0],expanded.stdout
  authenticated1=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),"read_complete_auth",str(complete1)],capture_output=True,text=True,env=env)
  assert authenticated1.returncode==0,(authenticated1.stdout,authenticated1.stderr)
  for mode in ("write_bad_explicit_compat_complete","write_bad_explicit_proof_complete",
               "write_bad_extension_complete","write_bad_boundary_complete","write_bad_window_mode_complete",
               "write_bad_face_id_complete","write_bad_nonlocal_id_complete",
               "write_bad_electron_defect_complete","write_bad_omitted_tail_complete",
               "write_bad_face_point_offset_complete","write_bad_face_weight_offset_complete",
               "write_bad_face_basis_offset_complete","write_bad_face_value_offset_complete",
               "write_bad_face_observable_offset_complete","write_bad_face_point_tail_complete",
               "write_bad_face_weight_tail_complete","write_bad_face_basis_tail_complete",
               "write_bad_face_value_tail_complete","write_bad_face_observable_tail_complete"):
    invalid_semantic=subprocess.run([shutil.which("mpiexec"),"-n","1",str(exe),mode,str(complete1)],capture_output=True,text=True,env=env)
    assert invalid_semantic.returncode==0,(mode,invalid_semantic.stdout,invalid_semantic.stderr)
  for mode in ("write_duplicate_face_id_complete","write_duplicate_nonlocal_id_complete"):
    duplicate_id=subprocess.run([shutil.which("mpiexec"),"-n","2",str(exe),mode,
      str(build/f"{mode}.chk")],capture_output=True,text=True,env=env)
    assert duplicate_id.returncode==0,(mode,duplicate_id.stdout,duplicate_id.stderr)
  authenticated=subprocess.run([shutil.which("mpiexec"),"-n","8",str(exe),"read_complete_auth",str(complete8)],capture_output=True,text=True,env=env)
  assert authenticated.returncode==0,(authenticated.stdout,authenticated.stderr)
  wrong_complete_version=build/"complete-v2-rejected.chk";wrong_payload=bytearray(complete8.read_bytes())
  struct.pack_into("=i",wrong_payload,16,2);wrong_complete_version.write_bytes(wrong_payload)
  rejected_version=subprocess.run([shutil.which("mpiexec"),"-n","8",str(exe),"read_complete_corrupt",str(wrong_complete_version)],capture_output=True,text=True,env=env)
  assert rejected_version.returncode==0,(rejected_version.stdout,rejected_version.stderr)
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
  assert struct.unpack_from("=i",checkpoint.read_bytes(),16)[0]==2,"generic checkpoint version changed"
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
  assert struct.unpack_from("=i",occupied.read_bytes(),16)[0]==2,"occupied checkpoint version changed"
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
