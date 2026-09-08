#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,struct,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
initialization_source=(root/"src/rt/dg/rt_dg_hybrid_initialization.f90").read_text().lower()
density_update_source=(root/"src/rt/dg/rt_dg_hybrid_density_update.f90").read_text().lower()
structural_source=(root/"src/rt/dg/rt_dg_hybrid_structural_graph.f90").read_text().lower()
projection_source=(root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90").read_text().lower()
rt_environment_source=(root/"src/rt/initialization_rt.f90").read_text().lower()
main_rt_source=(root/"src/rt/main_tddft.f90").read_text().lower()
main_rt_code="\n".join(line.split("!",1)[0] for line in main_rt_source.splitlines())
initializer=initialization_source.split("subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[1].split(
  "end subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[0]
compact=re.sub(r"\s+|&", "", initialization_source)
initializer_compact=re.sub(r"\s+|&", "", initializer)
density_compact=re.sub(r"\s+|&", "", density_update_source)
main_compact=re.sub(r"\s+|&", "", main_rt_code)

assert "type,public::s_rt_dg_hybrid_v3_startup_receipt" in compact, (
  "RED Task 14: missing public certified-v3 startup receipt"
)
assert re.search(r"public::[a-z0-9_,]*validate_rt_dg_hybrid_v3_startup",compact), (
  "RED Task 14: missing public post-redistribution v3 startup validator"
)
assert re.search(
  r"subroutineinitialize_rt_dg_hybrid_from_checkpoint\([^)]*xctype,tolerances,state,ok,message\)",
  compact,
), "RED Task 14: RT initializer does not accept user orbital/density/electron/symmetry tolerances"
for token in (
  "payload%certified_basis%initial_occupied_amplitudes",
  "payload%certified_basis%certified_eigenvalues",
  "payload%rt_space%metric_rows",
  "payload%rt_space%hamiltonian_rows",
  "payload%rt_space%basis_values",
  "payload%rt_space%representation",
  "callvalidate_rt_dg_hybrid_v3_startup",
):
  assert token in compact, f"RED Task 14: initializer omits certified v3 field {token}"
for forbidden in (
  "coefficient_buffer=payload%coefficients",
  "source=payload%eigenvalues",
  "source=payload%basis_values",
):
  assert forbidden not in compact, f"RED Task 14: initializer falls back to legacy construction payload: {forbidden}"
assert "rt_dg_hybrid_ground_state_checkpoint_version" in initialization_source, (
  "RED Task 14: SALMON_DG_GS001 probe does not use the public v3 version"
)
assert "version/=2" not in compact and "version==2" not in compact, (
  "RED Task 14: SALMON_DG_GS001 probe still accepts complete ground-state version 2"
)
assert "read_rt_dg_hybrid_ground_state_checkpoint_coalesced" in initializer, (
  "RED Task 14: initializer does not authenticate the payload after rank redistribution"
)
for forbidden in ("allocate(metric(n,n)","position_rows(3,n,n)","metric_graph(n,n)","operator_graph(n,n)"):
  assert forbidden not in initializer, f"RT initialization replicates a dense construction object: {forbidden}"
assert "callapply_construction_symmetry" not in compact and \
  "callcompute_construction_target_closure" in compact, (
  "RED Task 14: construction-space symmetry closure still broadcasts one row per operation"
)
builder=compact.split("subroutinebuild_certified_rt_state",1)[1].split(
  "endsubroutinebuild_certified_rt_state",1
)[0]
assert "callbuild_rt_dg_hybrid_structural_graph" in builder
assert "entry_tolerance" not in builder and "epsilon(1d0)*entry_scale" not in structural_source, (
  "Hybrid RT structural support must not drop exact basis or fixed-operator pairs by value threshold"
)
for token in ("basis_values(i,p)/=(0d0,0d0)","pair_key", "merge_key_sets"):
  assert token in structural_source.replace(" ",""), f"missing exact structural support token: {token}"
assert "compute_construction_projection" not in builder and \
  "startup_projected_position,startup_projected_basis" in initializer_compact, (
  "RED Task 14: validated construction projections are recomputed while building the RT state"
)

state_declaration=compact.split("type,public::s_rt_dg_hybrid_state",1)[1].split(
  "endtypes_rt_dg_hybrid_state",1
)[0]
assert "certified_rank" in state_declaration, (
  "RED Task 14: Hybrid RT state does not name its certified coefficient rank"
)
for routine in ("reconstruct_rt_dg_hybrid_density","update_rt_dg_hybrid_density"):
  body=density_compact.split(f"subroutine{routine}",1)[1].split(f"endsubroutine{routine}",1)[0]
  assert "callvalidate_certified_rt_state" in body, (
    f"RED Task 14: {routine} does not collectively reject construction-rank state extents"
  )
reconstruction=density_compact.split("subroutinereconstruct_rt_dg_hybrid_density",1)[1].split(
  "endsubroutinereconstruct_rt_dg_hybrid_density",1
)[0]
assert "global_coefficients(state%certified_rank,state%noccupied)" not in reconstruction, (
  "RED: Hybrid density reconstruction still replicates an O(rank*occupied) coefficient matrix"
)
assert "orbital_values(state%noccupied,max_point_count)" in reconstruction, (
  "RED: Hybrid density reconstruction does not reduce grid-local orbital amplitudes"
)

continuation=main_compact.split("subroutinerun_dg_hybrid_continuation_rt()",1)[1].split(
  "endsubroutinerun_dg_hybrid_continuation_rt",1
)[0]
assert "xc_func%xctype,[dg_dc_gs_final_orbital_tolerance,dg_dc_gs_final_density_tolerance," \
  "dg_dc_gs_electron_count_tolerance,dg_ow_symmetry_tolerance],hybrid_state,ok,message)" in continuation, (
  "RED Task 14: main does not pass the four user startup tolerances"
)
for call in (
  "callinitialize_rt_dg_hybrid_stationarity(nproc_group_global,hybrid_state%certified_rank",
  "callevaluate_rt_dg_hybrid_stationarity(nproc_group_global,hybrid_state%certified_rank",
  "callpropagate_rt_dg_hybrid_length_gauge(nproc_group_global,hybrid_state%certified_rank",
):
  assert call in continuation, f"RED Task 14: stale production API call: {call}"
step_tail=continuation.split("dostep=1,nt",1)[1]
outer_loop_end="enddoif(update_count/=nt+1)"
assert outer_loop_end in step_tail, "RED Task 14: cannot prove the certified RT timestep loop boundary"
step_body=step_tail.split(outer_loop_end,1)[0]
positions=[step_body.find(token) for token in (
  "electric_field=", "callpropagate_rt_dg_hybrid_length_gauge", "callreconstruct_rt_dg_hybrid_density",
  "callupdate_rt_dg_hybrid_density", "callevaluate_rt_dg_hybrid_stationarity",
)]
assert min(positions)>=0 and positions==sorted(positions), (
  "RED Task 14: density/stationarity is not evaluated after each propagated state"
)
orbital_tail=step_body.split("doorbital=1,hybrid_state%noccupied",1)[1]
assert "callpropagate_rt_dg_hybrid_length_gauge" in orbital_tail and \
  "enddocallreconstruct_rt_dg_hybrid_density" in orbital_tail, (
  "RED Task 14: density reconstruction can occur before all occupied orbitals are propagated"
)
assert "zero_field_run=zero_field_run.and." in continuation, (
  "RED Task 14: stationarity can be re-enabled after a nonzero pulse"
)
assert "if(step==1.and.trim(ae_shape1)=='impulse')then" in continuation and \
  "electric_field=-vector_potential_samples(:,step)/dt" in continuation, (
  "RED Task 14: Hybrid response route drops SALMON's t=0 impulse field"
)
assert "previous_polarization(3,hybrid_state%noccupied)" in continuation, (
  "RED Task 14: length-gauge branch history is shared between occupied orbitals"
)
assert "owner=mod(hybrid_state%global_count-row,nproc)" not in main_compact, (
  "RED Task 14: local-potential projection still uses reverse construction-era row ownership"
)
physical_invariants=main_compact.split("subroutineevaluate_hybrid_rt_physical_invariants",1)[1].split(
  "endsubroutineevaluate_hybrid_rt_physical_invariants",1
)[0]
assert "electron_count=sum(hybrid_state%occupations)" not in physical_invariants and \
  "sum(hybrid_state%density*hybrid_state%grid_weights)" in physical_invariants, (
  "RED Task 14: stationarity electron count is fixed by occupations instead of reconstructed charge"
)
local_projection=main_compact.split("subroutineproject_salmon_local_rows",1)[1].split(
  "endsubroutineproject_salmon_local_rows",1
)[0]
assert "row_offsets,column_ids" in local_projection and "local_values(:)" in local_projection, (
  "Hybrid local-potential callback does not consume the owned frozen CSR graph"
)
for forbidden in ("allocate(projected_local(hybrid_state%certified_rank,hybrid_state%certified_rank))",
                  "mpi_allreduce(mpi_in_place,projected_local", "doj=1,hybrid_state%certified_rank"):
  assert forbidden not in local_projection, f"Hybrid local-potential projection remains dense: {forbidden}"
assert "callproject_rt_dg_hybrid_sparse_edges" in local_projection, (
  "SALMON physical callback bypasses the tested sparse projection kernel"
)
assert "mpi_reduce_scatter" in projection_source, (
  "Hybrid local-potential projection does not reduce sparse edge contributions to row owners"
)
density_update=compact_density=density_compact.split("subroutineupdate_rt_dg_hybrid_density",1)[1].split(
  "endsubroutineupdate_rt_dg_hybrid_density",1)[0]
for forbidden in ("new_local(size(state%owned_row_ids),state%certified_rank)",
                  "new_h(size(state%owned_row_ids),state%certified_rank)",
                  "size(state%operators%column_ids)/=nowned*r",
                  "state%operators%row_offsets(nowned+1)/=nowned*r+1"):
  assert forbidden not in density_compact, f"Hybrid RT density update requires dense rows: {forbidden}"
assert "project_local(state%owned_row_ids,state%operators%row_offsets," \
  "state%operators%column_ids,state%grid_ids,density,new_local" in density_update, (
  "Hybrid RT density update does not request only frozen operator CSR values"
)

assert "subroutine initialization_rt_dg_hybrid" in rt_environment_source
hybrid_branch=main_rt_code.split("if(yn_rt_dg_hybrid_continuation=='y')then",1)[1].split("endif",1)[0]
assert "call initialization_rt_dg_hybrid" in hybrid_branch
assert "hybrid_basis_only" not in main_rt_source and "hybrid_basis_only" not in rt_environment_source
hybrid_initializer=rt_environment_source.split("subroutine initialization_rt_dg_hybrid",1)[1].split(
  "end subroutine initialization_rt_dg_hybrid",1)[0]
for forbidden in ("spsi_in","spsi_out","tpsi"):
  assert forbidden not in hybrid_initializer, f"hybrid initializer still exposes {forbidden}"

if os.environ.get("SALMON_LAPACK_LIBS"):
  libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip()
  libs=[f"-L{prefix}/lib","-lopenblas"]
else:
  libs=["-llapack","-lblas"]

def run_mpi(exe,nrank,path,mode,env):
  run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe),str(path),mode],
    capture_output=True,text=True,env=env,timeout=90)
  assert run.returncode==0,(nrank,mode,run.stdout,run.stderr)
  assert f"PASS hybrid RT v3 initialization contract on {nrank} ranks mode={mode}" in run.stdout
  return run

def replace_all_i64(source,target,old,new,expected_count):
  payload=bytearray(source.read_bytes());old_bytes=struct.pack("=q",old);new_bytes=struct.pack("=q",new)
  positions=[];start=0
  while (position:=payload.find(old_bytes,start))>=0:
    positions.append(position);payload[position:position+8]=new_bytes;start=position+8
  assert len(positions)==expected_count,(old,positions)
  target.write_bytes(payload)

with tempfile.TemporaryDirectory(prefix="hybrid-rt-v3-init-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_rt_v3_init"
  sources=["src/common/dg_hybrid_sparse_metric.f90","src/common/dg_hybrid_sparse_operators.f90",
    "src/rt/dg/rt_dg_hybrid_checkpoint.f90","src/rt/dg/rt_dg_hybrid_structural_graph.f90",
    "src/rt/dg/rt_dg_hybrid_sparse_projection.f90","src/rt/dg/rt_dg_hybrid_initialization.f90",
    "src/rt/dg/rt_dg_hybrid_density_update.f90",
    "tests/dg/test_rt_dg_hybrid_initialization_mpi.f90"]
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",*[str(root/s) for s in sources],
    *libs,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")

  for nrank in (1,2,4):
    run_mpi(exe,nrank,build/f"roundtrip-{nrank}.chk","roundtrip",env)

  negative_zero=build/"negative-zero.chk"
  run_mpi(exe,1,negative_zero,"write_negative_zero",env)
  run_mpi(exe,2,negative_zero,"read_only",env)
  for mode in ("reject_named_fingerprint","reject_catalog_fingerprint","reject_catalog_semantics",
      "reject_top_fingerprint","reject_receipt_fingerprint","reject_receipt_semantics",
      "reject_handoff_fingerprint","reject_position_convention_semantics","reject_basis_projection",
      "reject_component_projection","reject_orbital_physics","reject_metric_physics",
      "reject_unitarity_physics","reject_density_physics","reject_electron_physics",
      "reject_target_closure","reject_energy_covariance","reject_projector_covariance",
      "reject_vector_covariance","reject_tensor_covariance","density_tolerance_boundary",
      "embedding_projector_tolerance_boundary"):
    run_mpi(exe,2,build/f"{mode}.chk",mode,env)

  four_rank=build/"v3-four-rank.chk"
  write4=run_mpi(exe,4,four_rank,"write_only",env)
  accepted=four_rank.read_bytes()
  assert struct.unpack_from("=i",accepted,16)[0]==3,"fixture is not SALMON_DG_GS001 version 3"
  for nrank in (2,1):
    run_mpi(exe,nrank,four_rank,"read_only",env)
    run_mpi(exe,nrank,four_rank,"tamper_receipt_after_redistribution",env)
    run_mpi(exe,nrank,four_rank,"tamper_provenance_after_redistribution",env)

  two_rank=build/"v3-two-rank.chk"
  run_mpi(exe,2,two_rank,"write_only",env)
  run_mpi(exe,4,two_rank,"read_only",env)

  for version in (2,4):
    invalid_version=build/f"invalid-v{version}.chk";invalid_bytes=bytearray(accepted)
    struct.pack_into("=i",invalid_bytes,16,version);invalid_version.write_bytes(invalid_bytes)
    run_mpi(exe,2,invalid_version,"reject_only",env)

print("PASS hybrid RT certified-v3 initialization on 1, 2, and 4 ranks")
