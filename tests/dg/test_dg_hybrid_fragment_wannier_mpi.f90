#include "config.h"
module dg_hybrid_fragment_wannier_test_stubs
  implicit none
  integer,parameter::history_limit=8
  integer::setup_calls=0,run_calls=0,expected_fragment_id=0
  integer::fail_run_fragment=0
  logical::override_test_centers=.false.
  character(1024)::setup_seed_history(history_limit)='',run_seed_history(history_limit)=''
  integer::run_band_count_history(history_limit)=0
  logical::run_zero_auxiliary_energies(history_limit)=.false.
  logical::setup_saw_dmn(history_limit)=.false.
  logical::setup_saw_site_true(history_limit)=.false.
  logical::setup_saw_site_false(history_limit)=.false.
  logical::setup_saw_symmetrize(history_limit)=.false.
  logical::setup_saw_foreign_fragment(history_limit)=.false.
contains
  subroutine reset_w90_stub_state
    setup_calls=0;run_calls=0;expected_fragment_id=0;fail_run_fragment=0
    override_test_centers=.false.
    setup_seed_history='';run_seed_history='';setup_saw_dmn=.false.
    run_band_count_history=0;run_zero_auxiliary_energies=.false.
    setup_saw_site_true=.false.;setup_saw_site_false=.false.
    setup_saw_symmetrize=.false.;setup_saw_foreign_fragment=.false.
  end subroutine reset_w90_stub_state

  function fragment_token(fragment_id)result(token)
    integer,intent(in)::fragment_id
    character(15)::token
    write(token,'("fragment-",i6.6)')fragment_id
  end function fragment_token

  logical function run_must_fail(seed)
    character(*),intent(in)::seed
    run_must_fail=fail_run_fragment>0.and.&
      index(seed,trim(fragment_token(fail_run_fragment)))>0
  end function run_must_fail
end module dg_hybrid_fragment_wannier_test_stubs

#ifndef DG_W90_STUBS_ONLY
program test_dg_hybrid_fragment_wannier_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_w90,only:apply_dg_w90_gamma_transform
  use dg_hybrid_fragment_wannier,only:s_dg_hybrid_fragment_wannier_cache,&
    build_dg_hybrid_fragment_wannier,export_dg_hybrid_fragment_coordinates,pack_dg_hybrid_fragment_dc_seed,&
    build_dg_hybrid_fragment_wannier_from_dc_seed,map_dg_hybrid_fragment_dc_grid,&
    redistribute_dg_hybrid_fragment_wannier_columns
  use dg_hybrid_fragment_wannier_test_stubs
  use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,&
    initialize_dg_hybrid_fragment_subspace,s_dg_hybrid_fragment_epoch_budget,advance_dg_hybrid_fragment_epoch
  use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,freeze_dg_hybrid_variational_payload
  use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows
  use dg_hybrid_divided_operator,only:extract_dg_hybrid_fragment_self_block,&
    dg_hybrid_fragment_directory_fingerprint,freeze_dg_hybrid_single_owner_payload
  use dg_hybrid_fragment_preconditioner,only:s_dg_hybrid_preconditioner_key,s_dg_hybrid_fragment_preconditioner,&
    prepare_dg_hybrid_fragment_preconditioner,apply_dg_hybrid_fragment_preconditioner
  use dg_hybrid_fragment_solver,only:measure_dg_hybrid_fragment_core_norms
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_local_fragment_basis,&
    s_dg_hybrid_projection_factorization_receipt
  implicit none
  integer,parameter::global_ngrid=8,buffer_global_point=8,projector_global_point=7
  integer,parameter::first_generation=7
  real(real64),parameter::metric_tolerance=1d-12,localization_tolerance=1d-10
  character(*),parameter::artifact_root='fragment-wannier-artifacts'
  integer::comm_total,comm_fragment,total_rank,total_size,fragment_rank,fragment_size
  integer::fragment_id,nseed,ncandidate,nlocal,i,ierr,setup_before,run_before,local_slot
  integer::scf_construction_calls,scf_iterations
  integer(int64),allocatable::grid_ids(:)
  integer(int64),allocatable::reversed_grid_ids(:)
  real(real64),allocatable::grid_weights(:),fractional_coordinates(:,:)
  real(real64),allocatable::reversed_grid_weights(:),reversed_fractional_coordinates(:,:)
  real(real64),allocatable::physical_energies(:),physical_occupations(:)
  real(real64)::real_lattice(3,3),reciprocal_lattice(3,3),atoms_cart(3,1)
  real(real64)::saved_energy,saved_occupation
  complex(real64),allocatable::dc_seed_values(:,:),buffer_candidates(:,:),&
    projector_candidates(:,:)
  complex(real64),allocatable::swapped_seed_values(:,:),reversed_seed_values(:,:),&
    reversed_buffer_candidates(:,:),reversed_projector_candidates(:,:),&
    rank_deficient_seed_values(:,:),duplicate_projector_candidates(:,:)
  complex(real64)::saved_seed_value,saved_buffer_value,saved_projector_value
  character(2)::atom_symbols(1)
  character(512)::message
  logical::ok,cache_unchanged
  type(s_dg_hybrid_fragment_wannier_cache)::cache,first_snapshot,new_cache,failed_cache,&
    rank_deficient_cache,metric_null_cache,mixed_cache,mixed_snapshot,split_cache,&
    corrupted_cache
  complex(real64),allocatable::integration_h(:,:),integration_s(:,:)
  integer(int64),allocatable::integration_rows(:)
  type(s_dg_hybrid_fragment_preconditioner)::integration_preconditioner
  type(s_dg_hybrid_preconditioner_key)::integration_key

  call MPI_Init(ierr);comm_total=MPI_COMM_WORLD
  call MPI_Comm_rank(comm_total,total_rank,ierr)
  call MPI_Comm_size(comm_total,total_size,ierr)
  call require_total(any(total_size==[2,4,8]),'fixture requires exactly 2, 4, or 8 total ranks')
  fragment_id=mod(total_rank,2)+1
  call MPI_Comm_split(comm_total,fragment_id,total_rank,comm_fragment,ierr)
  call MPI_Comm_rank(comm_fragment,fragment_rank,ierr)
  call MPI_Comm_size(comm_fragment,fragment_size,ierr)
  call require_total(fragment_size==total_size/2,'two fragment communicators are not balanced')

  nseed=fragment_id+1;ncandidate=nseed+2
  nlocal=count([(mod(i-1,fragment_size)==fragment_rank,i=1,global_ngrid)])
  allocate(grid_ids(nlocal),grid_weights(nlocal),fractional_coordinates(3,nlocal))
  allocate(physical_energies(nseed),physical_occupations(nseed))
  allocate(dc_seed_values(nseed,nlocal),buffer_candidates(1,nlocal),&
    projector_candidates(1,nlocal))
  call fill_fixture
  call test_dc_tensor_packing
  call test_orbital_distributed_packing
  call test_dc_physical_grid_mapping
  call require_unique_fragment_rows
  call reset_w90_stub_state;expected_fragment_id=fragment_id
  cache%valid=.false.

  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  call require_total(ok,'initial fragment Wannier build failed: '//trim(message))
  call validate_successful_cache(cache,first_generation)
  call require_stub_call(1,first_generation)
  call require_two_distinct_fragment_seeds(1)
  first_snapshot=cache
  call test_coordinate_export

  setup_before=setup_calls;run_before=run_calls
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  call require_total(ok,'valid identical fragment cache was not reusable: '//trim(message))
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'identical fragment cache re-entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'identical cache reuse changed the published payload')

  corrupted_cache=first_snapshot
  if(total_rank==0)corrupted_cache%centers_fractional(1,1)=&
    corrupted_cache%centers_fractional(1,1)+0.125d0
  call require_corrupted_cache_rejected(corrupted_cache,'changed center on one rank')
  corrupted_cache=first_snapshot
  if(total_rank==0)deallocate(corrupted_cache%centers_fractional)
  call require_corrupted_cache_rejected(corrupted_cache,'missing centers on one rank')
  corrupted_cache=first_snapshot
  if(total_rank==0)then
    deallocate(corrupted_cache%centers_fractional)
    allocate(corrupted_cache%centers_fractional(2,ncandidate));corrupted_cache%centers_fractional=0d0
  endif
  call require_corrupted_cache_rejected(corrupted_cache,'wrong center extent on one rank')
  corrupted_cache=first_snapshot
  if(total_rank==0)corrupted_cache%centers_fractional(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call require_corrupted_cache_rejected(corrupted_cache,'nonfinite center on one rank')

  corrupted_cache=first_snapshot
  corrupted_cache%wannier_transform(1,1)=corrupted_cache%wannier_transform(1,1)+&
    cmplx(0.125d0,0d0,real64)
  call require_corrupted_cache_rejected(corrupted_cache,&
    'finite Wannier transform corruption')

  corrupted_cache=first_snapshot
  corrupted_cache%candidate_compression(1,1)=corrupted_cache%candidate_compression(1,1)+&
    cmplx(0d0,0.125d0,real64)
  call require_corrupted_cache_rejected(corrupted_cache,&
    'finite candidate-compression corruption')

  corrupted_cache=first_snapshot
  corrupted_cache%dc_seed_coefficients_in_wannier(1,1)=&
    corrupted_cache%dc_seed_coefficients_in_wannier(1,1)+cmplx(0.125d0,0d0,real64)
  call require_corrupted_cache_rejected(corrupted_cache,&
    'finite physical-seed coefficient corruption')

  corrupted_cache=first_snapshot
  corrupted_cache%physical_dc_seed_energies(1)=&
    corrupted_cache%physical_dc_seed_energies(1)+0.125d0
  call require_corrupted_cache_rejected(corrupted_cache,&
    'finite cached physical-energy corruption')

  corrupted_cache=first_snapshot
  corrupted_cache%physical_dc_seed_occupations(1)=&
    corrupted_cache%physical_dc_seed_occupations(1)-0.125d0
  call require_corrupted_cache_rejected(corrupted_cache,&
    'finite cached physical-occupation corruption')

  corrupted_cache=first_snapshot
  corrupted_cache%receipt%seed_reconstruction_defect=ieee_value(0d0,ieee_quiet_nan)
  call require_corrupted_cache_rejected(corrupted_cache,&
    'nonfinite seed-reconstruction receipt corruption')

  corrupted_cache=first_snapshot
  if(fragment_rank==0)corrupted_cache%wannier_values(1,1)=&
    corrupted_cache%wannier_values(1,1)+cmplx(0d0,0.125d0,real64)
  call require_corrupted_cache_rejected(corrupted_cache,&
    'fragment-root Wannier payload corruption')

  if(fragment_size>1)then
    corrupted_cache=first_snapshot
    if(fragment_rank==fragment_size-1)corrupted_cache%wannier_values(1,1)=&
      corrupted_cache%wannier_values(1,1)+cmplx(0.125d0,0d0,real64)
    call require_corrupted_cache_rejected(corrupted_cache,&
      'non-root rank-local Wannier payload corruption')
  endif

  swapped_seed_values=dc_seed_values
  call swap_seed_grid_rows(swapped_seed_values,1,2)
  call invoke_builder(first_generation,swapped_seed_values,buffer_candidates,cache,ok,message)
  call require_total(.not.ok,'row-associated DC seed fingerprint collision was accepted')
  call require_same_message(message,'row-associated seed collision did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'seed fingerprint')>0,&
    'row-associated seed collision message does not identify the seed fingerprint')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'row-associated seed collision entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'row-associated seed collision changed the valid cache')

  reversed_grid_ids=grid_ids(nlocal:1:-1)
  reversed_grid_weights=grid_weights(nlocal:1:-1)
  reversed_fractional_coordinates=fractional_coordinates(:,nlocal:1:-1)
  reversed_seed_values=dc_seed_values(:,nlocal:1:-1)
  reversed_buffer_candidates=buffer_candidates(:,nlocal:1:-1)
  reversed_projector_candidates=projector_candidates(:,nlocal:1:-1)
  call invoke_builder_layout(comm_total,comm_fragment,first_generation,reversed_grid_ids,&
    reversed_grid_weights,reversed_fractional_coordinates,reversed_seed_values,&
    reversed_buffer_candidates,reversed_projector_candidates,cache,ok,message)
  call require_total(.not.ok,'same physical rows in a different local order reused the cache')
  call require_same_message(message,'local row layout mismatch did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'local row layout')>0,&
    'local row layout mismatch message does not identify its cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'local row layout mismatch entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'local row layout mismatch changed the valid cache')

  saved_energy=physical_energies(1);physical_energies(1)=saved_energy+0.125d0
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  physical_energies(1)=saved_energy
  call require_total(.not.ok,'same-generation stale seed fingerprint was accepted')
  call require_same_message(message,'stale seed mismatch did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'seed fingerprint')>0,&
    'stale seed failure message does not identify its cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'stale seed mismatch entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'stale seed mismatch partially replaced the valid cache')

  local_slot=local_position(1)
  if(local_slot>0)then
    saved_seed_value=dc_seed_values(1,local_slot)
    dc_seed_values(1,local_slot)=saved_seed_value+cmplx(0d0,0.125d0,real64)
  endif
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  if(local_slot>0)dc_seed_values(1,local_slot)=saved_seed_value
  call require_total(.not.ok,'same-generation changed DC seed values were accepted')
  call require_same_message(message,'changed DC seed values did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'seed fingerprint')>0,&
    'changed DC seed value message does not identify its seed-key cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'changed DC seed values entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'changed DC seed values partially replaced the valid cache')

  saved_occupation=physical_occupations(1)
  physical_occupations(1)=saved_occupation-0.125d0
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  physical_occupations(1)=saved_occupation
  call require_total(.not.ok,'same-generation changed DC seed occupation was accepted')
  call require_same_message(message,'changed DC seed occupation did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'seed fingerprint')>0,&
    'changed DC seed occupation message does not identify its seed-key cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'changed DC seed occupation entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'changed DC seed occupation partially replaced the valid cache')

  local_slot=local_position(buffer_global_point)
  if(local_slot>0)then
    saved_buffer_value=buffer_candidates(1,local_slot)
    buffer_candidates(1,local_slot)=saved_buffer_value+cmplx(0.25d0,0d0,real64)
  endif
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  if(local_slot>0)buffer_candidates(1,local_slot)=saved_buffer_value
  call require_total(.not.ok,'same-generation stale basis fingerprint was accepted')
  call require_same_message(message,'stale basis mismatch did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'basis fingerprint')>0,&
    'stale basis failure message does not identify its cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'stale basis mismatch entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'stale basis mismatch partially replaced the valid cache')

  local_slot=local_position(projector_global_point)
  if(local_slot>0)then
    saved_projector_value=projector_candidates(1,local_slot)
    projector_candidates(1,local_slot)=saved_projector_value+cmplx(0d0,0.25d0,real64)
  endif
  call invoke_builder(first_generation,dc_seed_values,buffer_candidates,cache,ok,message)
  if(local_slot>0)projector_candidates(1,local_slot)=saved_projector_value
  call require_total(.not.ok,'same-generation stale projector fingerprint was accepted')
  call require_same_message(message,'stale projector mismatch did not fail collectively')
  call require_total(index(message,'stale')>0.and.index(message,'basis fingerprint')>0,&
    'stale projector failure message does not identify its basis-key cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'stale projector mismatch entered Wannier90')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'stale projector mismatch partially replaced the valid cache')

  new_cache%valid=.false.
  call invoke_builder(first_generation+1,dc_seed_values,buffer_candidates,new_cache,ok,message)
  call require_total(ok,'new-generation fragment Wannier build failed: '//trim(message))
  call validate_successful_cache(new_cache,first_generation+1)
  call require_stub_call(2,first_generation+1)
  call require_two_distinct_fragment_seeds(2)
  call require_total(new_cache%receipt%seed_fingerprint/=&
      first_snapshot%receipt%seed_fingerprint.and.&
      new_cache%receipt%basis_fingerprint/=first_snapshot%receipt%basis_fingerprint,&
    'basis generation is absent from fragment fingerprints')
  call require_total(same_cache_payload(cache,first_snapshot),&
    'new generation overwrote the prior generation cache')

  setup_before=setup_calls;run_before=run_calls;scf_construction_calls=-1
  call mock_scf_loop(4,new_cache,scf_iterations,scf_construction_calls)
  call require_total(scf_iterations==4.and.scf_construction_calls==0,&
    'construction localization was invoked inside the mock SCF loop')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'mock SCF loop re-entered Wannier90')

  failed_cache%valid=.false.;fail_run_fragment=2
  call invoke_builder(first_generation+2,dc_seed_values,buffer_candidates,failed_cache,ok,message)
  fail_run_fragment=0
  call require_total(.not.ok,'one-fragment Wannier90 failure was not collective')
  call require_same_message(message,'one-fragment Wannier90 failure message differs by total rank')
  call require_total(index(message,'Wannier90')>0.and.index(message,'fragment')>0,&
    'Wannier90 failure message does not identify the fragment failure')
  call require_total(cache_is_unpublished(failed_cache),&
    'failed fragment build published a partial cache')
  call require_stub_call(3,first_generation+2)
  call require_two_distinct_fragment_seeds(3)

  rank_deficient_seed_values=dc_seed_values
  rank_deficient_seed_values(2,:)=rank_deficient_seed_values(1,:)
  rank_deficient_cache%valid=.false.;setup_before=setup_calls;run_before=run_calls
  call invoke_builder(first_generation+3,rank_deficient_seed_values,buffer_candidates,&
    rank_deficient_cache,ok,message)
  call require_total(.not.ok,'rank-deficient physical DC seeds were accepted')
  call require_same_message(message,'physical seed rank failure was not collective')
  call require_total(index(message,'physical seed')>0.and.index(message,'rank')>0,&
    'physical seed rank failure message does not identify its cause')
  call require_total(setup_calls==setup_before.and.run_calls==run_before,&
    'rank-deficient physical seeds entered Wannier90')
  call require_total(cache_is_unpublished(rank_deficient_cache),&
    'rank-deficient physical seeds published a partial cache')

  duplicate_projector_candidates=buffer_candidates
  metric_null_cache%valid=.false.;setup_before=setup_calls;run_before=run_calls
  call invoke_builder_layout(comm_total,comm_fragment,first_generation+4,grid_ids,&
    grid_weights,fractional_coordinates,dc_seed_values,buffer_candidates,&
    duplicate_projector_candidates,metric_null_cache,ok,message)
  call require_total(ok,'metric-null auxiliary candidate build failed: '//trim(message))
  call validate_metric_null_cache(metric_null_cache,first_generation+4)
  call require_total(setup_calls==setup_before+merge(1,0,fragment_rank==0).and.&
    run_calls==run_before+merge(1,0,fragment_rank==0),&
    'metric-null candidate did not call Wannier90 exactly once on each fragment root')
  call require_stub_call(4,first_generation+4,ncandidate-1)
  call require_two_distinct_fragment_seeds(4)

  mixed_cache%valid=.false.
  if(fragment_id==1)then
    call invoke_builder_layout(comm_fragment,comm_fragment,first_generation+5,grid_ids,&
      grid_weights,fractional_coordinates,dc_seed_values,buffer_candidates,&
      projector_candidates,mixed_cache,ok,message)
    call require_total(ok,'fragment-local cache prebuild failed: '//trim(message),comm_fragment)
  endif
  call MPI_Barrier(comm_total,ierr)
  mixed_snapshot=mixed_cache;setup_before=setup_calls;run_before=run_calls
  call invoke_builder(first_generation+5,dc_seed_values,buffer_candidates,mixed_cache,ok,message)
  call require_total(ok,'mixed cache-hit/fresh-build construction failed: '//trim(message))
  call require_total(setup_calls==setup_before+merge(1,0,fragment_id==2.and.fragment_rank==0).and.&
    run_calls==run_before+merge(1,0,fragment_id==2.and.fragment_rank==0),&
    'mixed cache-hit/fresh-build path called Wannier90 on the wrong fragment')
  cache_unchanged=.true.
  if(fragment_id==1)cache_unchanged=same_cache_payload(mixed_cache,mixed_snapshot)
  call require_total(cache_unchanged,'cache-hit fragment changed during mixed fragment build')
  call validate_successful_cache(mixed_cache,first_generation+5)
  call require_stub_call(5,first_generation+5)
  call require_two_distinct_fragment_seeds(5)

  if(total_size>=4)then
    split_cache%valid=.false.;setup_before=setup_calls;run_before=run_calls
    call invoke_builder_layout(comm_total,MPI_COMM_SELF,first_generation+6,grid_ids,&
      grid_weights,fractional_coordinates,dc_seed_values,buffer_candidates,&
      projector_candidates,split_cache,ok,message)
    call require_total(.not.ok,'duplicate fragment communicator roots were accepted')
    call require_same_message(message,'fragment communicator root failure was not collective')
    call require_total(index(message,'fragment communicator')>0.and.index(message,'root')>0,&
      'duplicate fragment communicator root message does not identify its cause')
    call require_total(setup_calls==setup_before.and.run_calls==run_before,&
      'duplicate fragment communicator roots entered Wannier90')
    call require_total(cache_is_unpublished(split_cache),&
      'duplicate fragment communicator roots published a partial cache')
  endif

  call test_dc_construction_entry
  if(total_rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment Wannier on ',total_size,' ranks'
  call MPI_Comm_free(comm_fragment,ierr)
  call MPI_Finalize(ierr)
contains

  subroutine test_dc_physical_grid_mapping
    integer(int64),allocatable::local_ids(:),physical(:)
    logical,allocatable::core(:)
    integer::mapping(4,3),a,p,x,y,z,scenario,local_count,core_count,total_core
    logical::correct
    ! DC ordering is core, positive buffer, then wrapped negative buffer.
    ! In particular x=1 is core, whereas x=3 is not.
    mapping=0;mapping(:,1)=[5,6,1,4];mapping(1:2,2)=[3,4];mapping(1:2,3)=[1,2]
    local_count=count([(mod(a-1,min(fragment_size,2))==fragment_rank,a=1,16)])
    allocate(local_ids(local_count));p=0
    do a=16,1,-1
      if(mod(a-1,min(fragment_size,2))/=fragment_rank)cycle
      p=p+1;local_ids(p)=a
    enddo
    call map_dg_hybrid_fragment_dc_grid(comm_fragment,[4,2,2],[2,2,2],[6,4,4],mapping,&
      local_ids,physical,core,ok,message)
    call require_total(ok,'DC physical-grid mapping failed: '//trim(message))
    correct=.true.
    do p=1,local_count
      a=int(local_ids(p))-1;x=mod(a,4)+1;y=mod(a/4,2)+1;z=a/8+1
      correct=correct.and.physical(p)==mapping(x,1)+6*((mapping(y,2)-1)+4*(mapping(z,3)-1))
      correct=correct.and.(core(p).eqv.(x<=2))
    enddo
    core_count=count(core)
    call MPI_Allreduce(core_count,total_core,1,MPI_INTEGER,MPI_SUM,comm_fragment,ierr)
    call require_total(correct.and.total_core==8,'DC mapping assumed centered core or reordered WF rows')
    do scenario=1,3
      if(fragment_rank==0)then
        if(scenario==1)mapping(1,1)=0
        if(scenario==2)local_ids(1)=17_int64
        if(scenario==3)mapping(2,1)=mapping(1,1)
      endif
      call map_dg_hybrid_fragment_dc_grid(comm_fragment,[4,2,2],[2,2,2],[6,4,4],mapping,&
        local_ids,physical,core,ok,message)
      call require_total(.not.ok.and..not.allocated(physical).and..not.allocated(core),&
        'invalid DC mapping published physical ownership')
      mapping(:,1)=[5,6,1,4]
      if(scenario==2.and.fragment_rank==0)local_ids(1)=int(16-mod(15,min(fragment_size,2)),int64)
    enddo
  end subroutine test_dc_physical_grid_mapping

  subroutine test_dc_construction_entry
    real(real64),allocatable::tensor(:,:,:,:,:,:,:),esp(:,:,:),occ(:,:,:)
    complex(real64),allocatable::buffer(:,:),projector(:,:),expected(:,:),reconstructed(:,:)
    integer(int64),allocatable::ids(:),physical_ids(:)
    logical,allocatable::core_mask(:)
    integer::physical_map(8,3)
    real(real64)::density(12),reference_density(12)
    type(s_dg_hybrid_fragment_wannier_cache)::built,snapshot
    integer::a,b,first,last,points,pass,saved_setup,saved_run
    first=1+fragment_rank*nseed/fragment_size;last=(fragment_rank+1)*nseed/fragment_size
    allocate(tensor(0:9,0:2,0:2,1,first:last,1,1),esp(nseed,1,1),occ(nseed,1,1))
    tensor=ieee_value(0d0,ieee_quiet_nan)
    do b=first,last;do a=1,8
      tensor(a,1,1,1,b,1,1)=merge(1d0,0d0,a==b)
    enddo;enddo
    esp(:,1,1)=physical_energies;occ(:,1,1)=physical_occupations
    points=merge(8,0,fragment_rank==0)
    allocate(ids(points),buffer(1,points),projector(1,points),expected(nseed,points))
    buffer=0d0;projector=0d0;expected=0d0
    do a=1,points
      ids(a)=a
      buffer(1,a)=merge(1d0,0d0,a==8);projector(1,a)=merge(1d0,0d0,a==7)
      do b=1,nseed;expected(b,a)=merge(1d0,0d0,a==b);enddo
    enddo
    saved_setup=setup_calls;saved_run=run_calls
    physical_map=0;physical_map(:,1)=[9,10,11,12,1,2,7,8]
    physical_map(1,2:3)=1
    do pass=1,4
      if(pass==3.and.total_rank==0)esp(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
      if(pass==4.and.total_rank==0)ids(1)=99_int64
      call build_dg_hybrid_fragment_wannier_from_dc_seed(comm_total,comm_fragment,comm_fragment,&
        fragment_id,13,artifact_root,[8,1,1],[1,1,1],[8,1,1],tensor,esp,occ,1d0,&
        ids,buffer,projector,metric_tolerance,real_lattice,reciprocal_lattice,atom_symbols,atoms_cart,&
        20,localization_tolerance,10000000_int64,built,ok,message)
      if(pass<=2)then
        call require_total(ok,'direct DC construction failed: '//trim(message))
        if(pass==1)snapshot=built
        reconstructed=matmul(transpose(built%wannier_values),built%dc_seed_coefficients_in_wannier)
        call require_total(all(abs(reconstructed-transpose(expected))<1d-10),&
          'direct DC construction lost original physical seeds')
        call map_dg_hybrid_fragment_dc_grid(comm_fragment,[8,1,1],[4,1,1],[12,1,1],physical_map,&
          built%local_grid_ids,physical_ids,core_mask,ok,message)
        call require_total(ok,'constructed WF physical mapping failed: '//trim(message))
        if(pass==1)call test_column_export(built,physical_ids,core_mask)
        if(pass==1.and.total_size==2)call test_dc_cache_to_projected_basis(built)
        density=0d0;reference_density=0d0
        do a=1,size(physical_ids)
          if(.not.core_mask(a))cycle
          density(int(physical_ids(a)))=sum(physical_occupations*abs(reconstructed(a,:))**2)
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,density,12,MPI_DOUBLE_PRECISION,MPI_SUM,comm_fragment,ierr)
        do b=1,nseed;reference_density(physical_map(b,1))=physical_occupations(b);enddo
        call require_total(all(abs(density-reference_density)<1d-10).and.&
          abs(sum(density)-sum(physical_occupations))<1d-10,&
          'DC-to-WF-to-physical-core mapping changed density or electron count')
        call require_total(setup_calls==saved_setup+merge(1,0,fragment_rank==0).and.&
          run_calls==saved_run+merge(1,0,fragment_rank==0),'DC entry reran construction on cache reuse')
      else
        call require_total(.not.ok,'invalid direct DC input accepted')
        call require_same_message(message,'DC input failure did not propagate across fragments')
      endif
      call require_total(same_cache_payload(built,snapshot),'DC entry mutated the published cache')
      call require_total(setup_calls==saved_setup+merge(1,0,fragment_rank==0).and.&
        run_calls==saved_run+merge(1,0,fragment_rank==0),'rejected DC input entered Wannier90')
      esp(:,1,1)=physical_energies
    enddo
  end subroutine test_dc_construction_entry

  subroutine test_dc_cache_to_projected_basis(built)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::built
    type(s_dg_hybrid_basis_catalog)::pw_catalog
    type(s_dg_hybrid_fragment_basis)::basis
    type(s_dg_hybrid_projection_factorization_receipt)::receipt
    integer(int64),allocatable::physical(:),core_ids(:)
    logical,allocatable::core(:)
    integer,allocatable::owners(:)
    complex(real64),allocatable::q(:,:),coeff(:,:),padded(:,:),reconstructed(:,:),expected(:,:)
    real(real64)::coordinates(3,8),windows(2,8),gv(3,1),density(12),reference(12)
    integer::mapping(8,3),ranks(2),nw,a,b,n,setup_saved,run_saved
    integer(int64)::workspace,fingerprint,wf_fingerprint
    ! Production topology: two fragments, one rank each. All raw DC buffer
    ! points are preserved and wrapped onto the shared twelve-point physical cell.
    setup_saved=setup_calls;run_saved=run_calls
    call export_dg_hybrid_fragment_coordinates(comm_fragment,fragment_id,13,&
      built%receipt%seed_fingerprint,built%receipt%basis_fingerprint,built%local_grid_ids,&
      built%local_row_layout_fingerprint,built,q,coeff,ok,message)
    call require_total(ok,'cache coordinate validation failed before PW projection: '//trim(message))
    mapping=1
    mapping(:,1)=[(1+modulo(a-1+6*(fragment_id-1),12),a=1,8)]
    call map_dg_hybrid_fragment_dc_grid(comm_fragment,[8,1,1],[6,1,1],[12,1,1],mapping,&
      built%local_grid_ids,physical,core,ok,message)
    call require_total(ok,'cache physical mapping failed before PW projection: '//trim(message))
    core_ids=pack(physical,core);nw=built%receipt%retained_rank
    call MPI_Allgather(nw,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm_total,ierr)
    allocate(owners(sum(ranks)))
    owners(:ranks(1))=1;owners(ranks(1)+1:)=2
    windows=0d0;coordinates=0d0;gv=0d0
    do a=1,8
      coordinates(1,a)=real(physical(a)-1_int64,real64)
      b=1+int((physical(a)-1_int64)/6_int64);windows(b,a)=1d0
    enddo
    allocate(pw_catalog%packets(2));pw_catalog%valid=.true.
    pw_catalog%packet_fingerprint=101_int64;pw_catalog%catalog_fingerprint=103_int64
    do a=1,2
      pw_catalog%packets(a)%fragment_id=a;pw_catalog%packets(a)%owner_rank=a-1
      pw_catalog%packets(a)%star_id=1
      allocate(pw_catalog%packets(a)%g_indices(1),source=[1])
    enddo
    call MPI_Allreduce(built%receipt%distributed_wannier_fingerprint,wf_fingerprint,1,&
      MPI_INTEGER8,MPI_BXOR,comm_total,ierr)
    call build_dg_hybrid_projected_local_fragment_basis(comm_total,12,2,fragment_id,core_ids,&
      [1d0,1d0,1d0,1d0,1d0,1d0],coordinates(:,1:6),windows(:,1:6),physical,built%wannier_values,&
      coordinates,windows,pw_catalog,gv,owners,1,metric_tolerance,wf_fingerprint,&
      basis,workspace,fingerprint,ok,message,basis_generation=13,projection_receipt=receipt)
    call require_total(ok,'validated DC cache to projected basis failed: '//trim(message))
    call require_total(basis%generation==13.and.receipt%valid.and.&
      count(basis%sector==1)==nw.and.all(basis%buffer_point_ids==physical),&
      'DC cache handoff lost WF columns, generation or buffer points')
    n=size(basis%global_ids)
    call require_total(n==nw+1.and.sum(abs(basis%buffer_values(:,nw+1))**2)>1d-4,&
      'DC cache integration fixture failed to exercise a nonzero PW complement')
    allocate(padded(n,nseed),expected(8,nseed));padded=0d0;expected=0d0
    padded(:nw,:)=coeff
    reconstructed=matmul(basis%buffer_values,padded)
    do a=1,8;do b=1,nseed
      if(built%local_grid_ids(a)==int(b,int64))expected(a,b)=1d0
    enddo;enddo
    call require_total(maxval(abs(reconstructed-expected))<1d-10,&
      'WF plus PW seed coefficients no longer reconstruct the original DC orbitals')
    density=0d0;reference=0d0
    do a=1,8
      if(.not.core(a))cycle
      density(physical(a))=sum(physical_occupations*abs(reconstructed(a,:))**2)
      reference(physical(a))=sum(physical_occupations*abs(expected(a,:))**2)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,density,12,MPI_DOUBLE_PRECISION,MPI_SUM,comm_total,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,reference,12,MPI_DOUBLE_PRECISION,MPI_SUM,comm_total,ierr)
    call require_total(maxval(abs(density-reference))<1d-10.and.abs(sum(density)-sum(reference))<1d-10,&
      'DC cache to WF plus PW changed physical-core density or electron count')
    call test_projected_basis_local_updates(built,basis,q,padded,core)
    call require_total(setup_calls==setup_saved.and.run_calls==run_saved,&
      'DC cache to projected basis unexpectedly reran Wannier90')
  end subroutine test_dc_cache_to_projected_basis

  subroutine test_projected_basis_local_updates(built,basis,q_wf,seed_coefficients,core)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::built
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    complex(real64),intent(in)::q_wf(:,:),seed_coefficients(:,:)
    logical,intent(in)::core(:)
    type(s_dg_hybrid_fixed_payload)::payload,rejected_payload
    type(s_dg_hybrid_fragment_basis)::bad_basis
    type(s_dg_hybrid_fragment_subspace_state)::state
    type(s_dg_hybrid_fragment_epoch_budget)::budget
    integer,allocatable::owners(:),slots(:),generations(:),selected(:),published_owners(:),&
      published_fragments(:),published_slots(:),published_generations(:)
    complex(real64),allocatable::full_basis(:,:),h_rows(:,:),s_rows(:,:),zero_rows(:,:),q(:,:),&
      physical_states(:,:),orthogonality(:,:)
    complex(real64)::physical_h(12,12)
    complex(real64),allocatable::core_values(:,:),core_gradients(:,:,:),unused_kinetic(:,:),reference_s(:,:)
    real(real64)::volume_diagnostics(4)
    character(64)::audit_mode
    real(real64),allocatable::energies(:),core_norms(:)
    real(real64)::residual,core_error
    integer::counts(2),nb,n,nw,a,b,j,pass,steps,remaining,total_steps
    integer(int64)::directory_fp,workspace,receipt,preconditioner_fp,published_directory_fp
    logical::converged,advanced
    character(256)::reason
    n=size(basis%global_ids);nw=size(q_wf,1)
    call MPI_Allgather(n,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm_total,ierr)
    nb=sum(counts)
    allocate(full_basis(12,nb),owners(nb),slots(nb),generations(nb),q(n,n))
    full_basis=0d0;owners=0;slots=0;generations=13;q=0d0
    do j=1,n
      owners(basis%global_ids(j))=fragment_id;slots(basis%global_ids(j))=j
      do a=1,size(core)
        full_basis(basis%buffer_point_ids(a),basis%global_ids(j))=basis%buffer_values(a,j)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,full_basis,size(full_basis),MPI_DOUBLE_COMPLEX,MPI_SUM,comm_total,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,owners,nb,MPI_INTEGER,MPI_SUM,comm_total,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,slots,nb,MPI_INTEGER,MPI_SUM,comm_total,ierr)
    ! Small independent physical-space oracle, not a production SIPG operator:
    ! an off-diagonal Hermitian H makes the DC seed a nonstationary warm start.
    physical_h=0d0
    do a=1,12
      physical_h(a,a)=0.2d0*a
      if(a==12)cycle
      physical_h(a,a+1)=cmplx(0.03d0,0.01d0,real64)
      physical_h(a+1,a)=conjg(physical_h(a,a+1))
    enddo
    h_rows=matmul(conjg(transpose(full_basis(:,basis%global_ids))),matmul(physical_h,full_basis))
    s_rows=matmul(conjg(transpose(full_basis(:,basis%global_ids))),full_basis)
    call get_command_argument(1,audit_mode)
    if(trim(audit_mode)=='--core-metric-audit')then
      ! Change only S relative to the existing diagnostic H: the production
      ! volume assembler integrates the owner's unique core, not its buffer.
      ! This is a handoff audit, not a complete SIPG/material calculation.
      core_values=transpose(full_basis(pack(basis%buffer_point_ids,core),:))
      allocate(core_gradients(3,nb,count(core)));core_gradients=0d0
      call assemble_dg_hybrid_broken_volume_rows(comm_total,nb,basis%global_ids,owners,&
        pack(basis%buffer_point_ids,core),[(fragment_id,a=1,count(core))],&
        [(1d0,a=1,count(core))],core_values,core_gradients,[(1d0,a=1,count(core))],&
        unused_kinetic,s_rows,volume_diagnostics,ok,message)
      call require_total(ok,'production core-only metric assembly failed: '//trim(message))
      reference_s=matmul(conjg(transpose(basis%buffer_values(pack([(a,a=1,size(core))],core),:))),&
        basis%buffer_values(pack([(a,a=1,size(core))],core),:))
      call require_total(maxval(abs(s_rows(:,basis%global_ids)-reference_s))<1d-12,&
        'production metric differs from independent unique-core Gram matrix')
    endif
    allocate(zero_rows(n,nb));zero_rows=0d0
    directory_fp=dg_hybrid_fragment_directory_fingerprint(owners,slots,generations,501_int64)
    call freeze_dg_hybrid_single_owner_payload(comm_total,2,basis,s_rows,h_rows,zero_rows,&
      zero_rows,501_int64,503_int64,505_int64,payload,published_owners,published_fragments,&
      published_slots,published_generations,published_directory_fp,ok,message)
    call require_total(ok,'DC integration fixed payload failed: '//trim(message))
    call require_total(all(published_owners==owners-1).and.all(published_fragments==owners).and.&
      all(published_slots==slots).and.all(published_generations==generations).and.&
      published_directory_fp==directory_fp.and.payload%basis_directory_fingerprint==directory_fp,&
      'production payload directory differs from independent fragment inventory')
    bad_basis=basis
    if(total_rank==0)bad_basis%global_ids(1)=0_int64
    call freeze_dg_hybrid_single_owner_payload(comm_total,2,bad_basis,s_rows,h_rows,zero_rows,&
      zero_rows,501_int64,503_int64,505_int64,rejected_payload,published_owners,published_fragments,&
      published_slots,published_generations,published_directory_fp,ok,message)
    call require_total(.not.ok.and..not.rejected_payload%frozen.and..not.allocated(published_owners).and.&
      published_directory_fp==0_int64,'invalid basis ID published a production payload or directory')
    integration_rows=[(int(a,int64),a=1,n)]
    call extract_dg_hybrid_fragment_self_block(comm_fragment,fragment_id,integration_rows,owners,payload,&
      zero_rows,integration_h,integration_s,ok,message,basis_local_slot=slots,basis_generation=generations,&
      fragment_catalog_fingerprint=501_int64,fragment_directory_fingerprint=directory_fp)
    call require_total(ok,'DC integration H/S extraction failed: '//trim(message))
    call require_total(maxval(abs(integration_h-h_rows(:,basis%global_ids)))<1d-12.and.&
      maxval(abs(integration_s-s_rows(:,basis%global_ids)))<1d-12,&
      'extracted fragment H/S differs from physical-space oracle')
    orthogonality=integration_s
    do a=1,n;orthogonality(a,a)=orthogonality(a,a)-1d0;enddo
    call require_total(maxval(abs(orthogonality))>1d-3,&
      'DC integration fixture accidentally uses an identity local metric')
    q(:nw,:nw)=q_wf
    do a=nw+1,n;q(a,a)=1d0;enddo
    if(trim(audit_mode)=='--core-metric-audit')then
      reference_s=matmul(conjg(transpose(q)),matmul(integration_s,q))
      write(*,'(a,i0,a,i0,a,*(es12.4,1x))')'CORE-METRIC-AUDIT fragment=',fragment_id,&
        ' columns=',n,' reference norms=',[(real(reference_s(a,a),real64),a=1,n)]
    endif
    integration_key=s_dg_hybrid_preconditioner_key(fragment_id,13,1,501_int64,503_int64,&
      payload%fingerprint,built%receipt%transform_fingerprint)
    call prepare_dg_hybrid_fragment_preconditioner(comm_fragment,n,integration_rows,q,integration_h,&
      integration_s,integration_key,1d-10,integration_preconditioner,preconditioner_fp,ok,message)
    call require_total(ok,'DC integration fixed-frame preconditioner failed: '//trim(message))
    call initialize_dg_hybrid_fragment_subspace(comm_fragment,n,integration_rows,fragment_id,13,&
      501_int64,503_int64,seed_coefficients,built%physical_dc_seed_energies,built%physical_dc_seed_occupations,&
      0,1d-10,1d-8,1d-10,integration_apply_s,state,selected,ok,message)
    call require_total(ok,'DC integration occupied subspace initialization failed: '//trim(message))
    call require_total(state%state_count==nseed.and.state%state_count<n,&
      'DC integration initialized the full local basis instead of the occupied seed inventory')
    allocate(energies(state%state_count));total_steps=0
    do pass=1,2
      call advance_dg_hybrid_fragment_epoch(comm_fragment,n,integration_rows,fragment_id,13,501_int64,503_int64,&
        1,3,integration_apply_h,integration_apply_s,integration_apply_preconditioner,1d-12,1d-10,2d0,&
        budget,state,energies,steps,remaining,residual,converged,advanced,reason,workspace,receipt,ok,message)
      call require_total(ok,'DC integration bounded fragment update failed: '//trim(message))
      total_steps=total_steps+steps
      call require_total(total_steps<=3.and.remaining==3-total_steps,&
        'repeated fragment update restarted the three-step density-epoch budget')
      physical_states=matmul(basis%buffer_values,state%vectors)
      orthogonality=matmul(conjg(transpose(state%vectors)),matmul(integration_s,state%vectors))
      do a=1,state%state_count;orthogonality(a,a)=orthogonality(a,a)-1d0;enddo
      call require_total(maxval(abs(orthogonality))<1d-9,'bounded update lost S orthogonality')
      call measure_dg_hybrid_fragment_core_norms(comm_fragment,basis,state%vectors,core,&
        [(1d0,a=1,size(core))],core_norms,ok,message)
      call require_total(ok,'DC integration current core weights failed: '//trim(message))
      core_error=0d0
      do b=1,state%state_count
        core_error=max(core_error,abs(core_norms(b)-sum(abs(physical_states(:,b))**2,mask=core)))
      enddo
      call require_total(core_error<1d-10,'current core weights disagree with reconstructed bounded-update orbitals')
    enddo
    call require_total(total_steps>0,'DC integration did not exercise a fragment CG update')
  end subroutine test_projected_basis_local_updates

  subroutine integration_apply_h(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    output=matmul(integration_h,input);valid=.true.
  end subroutine integration_apply_h

  subroutine integration_apply_s(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    output=matmul(integration_s,input);valid=.true.
  end subroutine integration_apply_s

  subroutine integration_apply_preconditioner(input,shifts,output,valid)
    complex(real64),intent(in)::input(:,:)
    real(real64),intent(in)::shifts(:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    complex(real64),allocatable::result(:,:)
    integer(int64)::receipt
    character(256)::detail
    call apply_dg_hybrid_fragment_preconditioner(comm_fragment,integration_rows,integration_key,&
      integration_preconditioner,shifts,input,result,receipt,valid,detail)
    if(valid)output=result
  end subroutine integration_apply_preconditioner

  subroutine test_column_export(built,physical_ids,core_mask)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::built
    integer(int64),intent(in)::physical_ids(:)
    logical,intent(in)::core_mask(:)
    integer(int64),allocatable::columns(:),mapped(:),bad_ids(:)
    logical,allocatable::mask(:)
    complex(real64),allocatable::values(:,:),first_values(:,:),seed(:,:)
    type(s_dg_hybrid_fragment_subspace_state)::initial
    type(s_dg_hybrid_fragment_wannier_cache)::corrupted
    real(real64)::initial_occupations(nseed)
    integer,allocatable::selected(:)
    integer::nw,ncols,a,b,p,width
    logical::correct
    nw=built%receipt%retained_rank
    ncols=count([(mod(a-1,min(fragment_size,2))==fragment_rank,a=1,nw)])
    allocate(columns(ncols));p=0
    do a=nw,1,-1
      if(mod(a-1,min(fragment_size,2))/=fragment_rank)cycle
      p=p+1;columns(p)=a
    enddo
    do width=1,3,2
      call redistribute_dg_hybrid_fragment_wannier_columns(comm_fragment,built,columns,physical_ids,&
        core_mask,width,mapped,mask,values,ok,message)
      call require_total(ok,'WF coefficient-column redistribution failed: '//trim(message))
      call require_total(all(shape(values)==[8,ncols]).and.all(mapped==[9,10,11,12,1,2,7,8]).and.&
        all(mask.eqv.[.true.,.true.,.true.,.true.,.false.,.false.,.false.,.false.]),&
        'column redistribution changed core/buffer row mapping')
      if(width==1)first_values=values
      call require_total(all(abs(first_values-values)<1d-14),'column redistribution depends on tile size')
      seed=matmul(values,built%dc_seed_coefficients_in_wannier(int(columns),:))
      call MPI_Allreduce(MPI_IN_PLACE,seed,size(seed),MPI_DOUBLE_COMPLEX,MPI_SUM,comm_fragment,ierr)
      correct=.true.
      do b=1,nseed;do a=1,8
        correct=correct.and.abs(seed(a,b)-merge(1d0,0d0,a==b))<1d-10
      enddo;enddo
      call require_total(correct,'distributed coefficient columns fail to reconstruct original DC orbitals')
    enddo
    initial_occupations=0d0;initial_occupations(1)=2d0
    call initialize_dg_hybrid_fragment_subspace(comm_fragment,nw,columns,fragment_id,13,&
      built%receipt%basis_fingerprint,901_int64,built%dc_seed_coefficients_in_wannier(int(columns),:),&
      physical_energies,initial_occupations,1,1d-8,1d-10,1d-10,fixture_identity_metric,&
      initial,selected,ok,message)
    call require_total(ok,'redistributed WF columns failed CG initialization: '//trim(message))
    call require_total(initial%state_count==2.and.all(selected==[1,2]),&
      'redistributed WF CG initializer changed occupied-plus-guard selection')
    seed=matmul(values,initial%vectors)
    call MPI_Allreduce(MPI_IN_PLACE,seed,size(seed),MPI_DOUBLE_COMPLEX,MPI_SUM,comm_fragment,ierr)
    correct=.true.
    do b=1,2;do a=1,8
      correct=correct.and.abs(seed(a,b)-merge(1d0,0d0,a==b))<1d-10
    enddo;enddo
    call require_total(correct,'redistributed CG initial state changed original physical orbitals')
    corrupted=built
    if(fragment_rank==0)corrupted%wannier_values(1,1)=corrupted%wannier_values(1,1)+0.1d0
    call redistribute_dg_hybrid_fragment_wannier_columns(comm_fragment,corrupted,columns,physical_ids,&
      core_mask,2,mapped,mask,values,ok,message)
    call require_total(.not.ok.and..not.allocated(values),'corrupted WF cache accepted by column export')
    bad_ids=physical_ids
    if(fragment_rank==0)bad_ids(1)=0_int64
    call redistribute_dg_hybrid_fragment_wannier_columns(comm_fragment,built,columns,bad_ids,&
      core_mask,2,mapped,mask,values,ok,message)
    call require_total(.not.ok.and..not.allocated(values).and..not.allocated(mapped),&
      'invalid physical mapping published coefficient columns')
    if(fragment_rank==0)columns(1)=columns(size(columns))
    call redistribute_dg_hybrid_fragment_wannier_columns(comm_fragment,built,columns,physical_ids,&
      core_mask,2,mapped,mask,values,ok,message)
    call require_total(.not.ok.and..not.allocated(values),'duplicate WF column ownership accepted')
  end subroutine test_column_export

  subroutine test_orbital_distributed_packing
    real(real64),allocatable::tensor(:,:,:,:,:,:,:),esp(:,:,:),occ(:,:,:)
    real(real64),allocatable::w(:),e(:),o(:),f(:,:)
    complex(real64),allocatable::v(:,:)
    integer(int64),allocatable::ids(:)
    integer::orbcomm,orb_rank,orb_size,spatial_group,lo(3),hi(3),first,last,a,b,scenario
    logical::correct
    ! 8 total ranks: two spatial slabs times two orbital ranks per fragment.
    ! 4 total ranks: orbital partition only; 2 total ranks: serial degeneration.
    spatial_group=0
    if(fragment_size==4)spatial_group=fragment_rank/2
    call MPI_Comm_split(comm_fragment,spatial_group,fragment_rank,orbcomm,ierr)
    call MPI_Comm_rank(orbcomm,orb_rank,ierr);call MPI_Comm_size(orbcomm,orb_size,ierr)
    lo=[1+4*spatial_group,1,1];hi=[8,1,1]
    if(fragment_size==4)hi(1)=4+4*spatial_group
    first=1+orb_rank*nseed/orb_size;last=(orb_rank+1)*nseed/orb_size
    allocate(tensor(-1:10,0:2,0:2,1,first:last,1,1),esp(nseed,1,1),occ(nseed,1,1))
    tensor=ieee_value(0d0,ieee_quiet_nan)
    do b=first,last;do a=lo(1),hi(1);tensor(a,1,1,1,b,1,1)=100d0*b+a;enddo;enddo
    esp(:,1,1)=physical_energies;occ(:,1,1)=physical_occupations
    call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],lo,hi,tensor,esp,occ,0.25d0,&
      ids,w,v,e,o,f,ok,message,orbital_comm=orbcomm)
    call require_total(ok,'orbital-distributed packing failed: '//trim(message))
    correct=size(ids)==merge(hi(1)-lo(1)+1,0,orb_rank==0)
    do a=1,size(ids)
      correct=correct.and.ids(a)==lo(1)+a-1
      do b=1,nseed;correct=correct.and.v(b,a)==100d0*b+ids(a);enddo
    enddo
    call require_total(correct,'orbital-distributed packing lost orbital order or duplicated grid')
    do scenario=1,4
      if(scenario==1.and.fragment_rank==0)hi(1)=hi(1)-1
      if(scenario==2.and.fragment_rank==0)tensor(lo(1),1,1,1,first,1,1)=ieee_value(0d0,ieee_quiet_nan)
      if(scenario==3.and.fragment_rank==0)then
        deallocate(tensor);allocate(tensor(-1:10,0:2,0:2,1,1:0,1,1))
      endif
      if(scenario==4.and.fragment_rank==0)then
        deallocate(tensor);allocate(tensor(-1:10,0:2,0:2,1,0:nseed,1,1));tensor=0d0
      endif
      call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],lo,hi,tensor,esp,occ,0.25d0,&
        ids,w,v,e,o,f,ok,message,orbital_comm=orbcomm)
      call require_total(.not.ok.and..not.allocated(v),'invalid orbital partition accepted')
      if(scenario==1.and.fragment_rank==0)hi(1)=hi(1)+1
      if(scenario==2.and.fragment_rank==0)tensor(lo(1),1,1,1,first,1,1)=100d0*first+lo(1)
    enddo
    call MPI_Comm_free(orbcomm,ierr)
    ! A single seed distributed across all fragment ranks leaves empty orbital owners.
    deallocate(tensor,esp,occ)
    first=1+fragment_rank/fragment_size;last=(fragment_rank+1)/fragment_size
    allocate(tensor(1:8,1,1,1,first:last,1,1),esp(1,1,1),occ(1,1,1))
    tensor=2d0;esp=0d0;occ=1d0
    call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],[1,1,1],[8,1,1],tensor,esp,occ,0.25d0,&
      ids,w,v,e,o,f,ok,message,orbital_comm=comm_fragment)
    call require_total(ok,'empty orbital owner rejected: '//trim(message))
    call require_total(size(ids)==merge(8,0,fragment_rank==0).and.all(v==2d0),&
      'empty orbital owner packing changed values or point ownership')
    call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],[1,1,1],[8,1,1],tensor,esp,occ,0.25d0,&
      ids,w,v,e,o,f,ok,message,orbital_comm=comm_total)
    call require_total(.not.ok.and..not.allocated(v).and.index(message,'subgroup')>0,&
      'cross-fragment orbital communicator not rejected before orbital collectives')
  end subroutine test_orbital_distributed_packing

  subroutine test_dc_tensor_packing
    real(real64),allocatable::tensor(:,:,:,:,:,:,:),esp(:,:,:),occ(:,:,:)
    real(real64),allocatable::packed_weights(:),packed_energies(:),packed_occupations(:),fractional(:,:)
    complex(real64),allocatable::packed(:,:)
    integer(int64),allocatable::packed_ids(:)
    integer::lo(3),hi(3),a,b,scenario,active_size,c,d,p
    real(real64)::test_hvol
    logical::correct
    active_size=min(fragment_size,2)
    lo=[1+fragment_rank*8/active_size,1,1];hi=[(fragment_rank+1)*8/active_size,1,1]
    if(fragment_rank>=active_size)then;lo=[1,1,1];hi=[0,1,1];endif
    allocate(tensor(-1:10,0:2,0:2,1,1:nseed,1,1),esp(nseed,1,1),occ(nseed,1,1))
    tensor=ieee_value(0d0,ieee_quiet_nan)
    do a=lo(1),hi(1)
      do b=1,nseed;tensor(a,1,1,1,b,1,1)=100d0*b+a;enddo
    enddo
    esp(:,1,1)=physical_energies;occ(:,1,1)=physical_occupations
    call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],lo,hi,tensor,esp,occ,0.25d0,&
      packed_ids,packed_weights,packed,packed_energies,packed_occupations,fractional,ok,message)
    call require_total(ok,'DC tensor packing failed: '//trim(message))
    correct=size(packed_ids)==max(0,hi(1)-lo(1)+1).and.all(packed_weights==0.25d0)
    do a=1,size(packed_ids)
      correct=correct.and.packed_ids(a)==lo(1)+a-1
      do b=1,nseed;correct=correct.and.packed(b,a)==100d0*b+packed_ids(a);enddo
      correct=correct.and.abs(fractional(1,a)-real(packed_ids(a)-1,real64)/8d0)<1d-14
    enddo
    call require_total(correct.and.all(packed_energies==physical_energies).and.&
      all(packed_occupations==physical_occupations),'DC tensor packing mixed halo/grid/orbital indices')
    do scenario=1,7
      test_hvol=0.25d0
      if(scenario==1.and.fragment_rank==0)hi(1)=hi(1)-1
      if(scenario==2.and.fragment_rank==0)esp(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
      if(scenario==3.and.fragment_rank==0)tensor(lo(1),1,1,1,1,1,1)=ieee_value(0d0,ieee_quiet_nan)
      if(scenario==4.and.fragment_rank==0)lo(1)=0
      if(scenario==5.and.fragment_rank==0)then
        deallocate(tensor);allocate(tensor(-1:10,0:2,0:2,1,2:nseed,1,1));tensor=0d0
      endif
      if(scenario==6.and.fragment_rank==0)then
        deallocate(occ)
      endif
      if(scenario==7.and.fragment_rank==0)test_hvol=ieee_value(0d0,ieee_quiet_nan)
      call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[8,1,1],lo,hi,tensor,esp,occ,test_hvol,&
        packed_ids,packed_weights,packed,packed_energies,packed_occupations,fractional,ok,message)
      call require_total(.not.ok.and..not.allocated(packed).and..not.allocated(packed_ids),&
        'invalid/incomplete DC tensor published construction inputs')
      if(scenario==1.and.fragment_rank==0)hi(1)=hi(1)+1
      if(scenario==3.and.fragment_rank==0)tensor(lo(1),1,1,1,1,1,1)=101d0
      if(scenario==4.and.fragment_rank==0)lo(1)=1
      if(scenario==5.and.fragment_rank==0)then
        deallocate(tensor);allocate(tensor(-1:10,0:2,0:2,1,1:nseed,1,1));tensor=0d0
      endif
      if(scenario==6.and.fragment_rank==0)then
        allocate(occ(nseed,1,1));occ(:,1,1)=physical_occupations
      endif
      esp(:,1,1)=physical_energies
    enddo
    if(.not.allocated(occ))allocate(occ(nseed,1,1))
    occ(:,1,1)=physical_occupations
    deallocate(tensor);allocate(tensor(0:3,0:3,0:3,1,1:nseed,1,1))
    lo=[1+fragment_rank*2/active_size,1,1];hi=[(fragment_rank+1)*2/active_size,2,2]
    if(fragment_rank>=active_size)then;lo=[1,1,1];hi=[0,2,2];endif
    tensor=ieee_value(0d0,ieee_quiet_nan)
    do d=1,2;do c=1,2;do a=lo(1),hi(1);do b=1,nseed
      tensor(a,c,d,1,b,1,1)=1000d0*b+a+10*c+100*d
    enddo;enddo;enddo;enddo
    call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[2,2,2],lo,hi,tensor,esp,occ,0.25d0,&
      packed_ids,packed_weights,packed,packed_energies,packed_occupations,fractional,ok,message)
    call require_total(ok,'three-dimensional DC tensor packing failed: '//trim(message))
    p=0;correct=.true.
    do d=1,2;do c=1,2;do a=lo(1),hi(1)
      p=p+1;correct=correct.and.packed_ids(p)==a+2*((c-1)+2*(d-1))
      correct=correct.and.all(fractional(:,p)==real([a-1,c-1,d-1],real64)/2d0)
      do b=1,nseed
        correct=correct.and.packed(b,p)==1000d0*b+a+10*c+100*d
      enddo
    enddo;enddo;enddo
    call require_total(correct,'three-dimensional DC tensor axis ordering changed')
    if(fragment_size>1)then
      if(fragment_rank==1)esp(1,1,1)=esp(1,1,1)+0.5d0
      call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[2,2,2],lo,hi,tensor,esp,occ,0.25d0,&
        packed_ids,packed_weights,packed,packed_energies,packed_occupations,fractional,ok,message)
      call require_total(.not.ok.and..not.allocated(packed),'rank-disagreeing DC spectra accepted')
      esp(:,1,1)=physical_energies
      if(fragment_rank==1)then
        lo(1)=1;hi(1)=1;tensor=0d0
      endif
      call pack_dg_hybrid_fragment_dc_seed(comm_fragment,[2,2,2],lo,hi,tensor,esp,occ,0.25d0,&
        packed_ids,packed_weights,packed,packed_energies,packed_occupations,fractional,ok,message)
      call require_total(.not.ok.and..not.allocated(packed),'duplicate owned DC grid accepted')
    endif
  end subroutine test_dc_tensor_packing

  subroutine test_coordinate_export
    complex(real64),allocatable::q(:,:),seed(:,:),raw(:,:),reference(:,:)
    type(s_dg_hybrid_fragment_wannier_cache)::bad
    integer::case_id,saved_setup,saved_run
    integer(int64)::expected_seed,expected_basis,expected_layout
    integer(int64),allocatable::expected_ids(:)
    saved_setup=setup_calls;saved_run=run_calls
    allocate(raw(ncandidate,nlocal))
    raw(1:nseed,:)=dc_seed_values
    raw(nseed+1:nseed+1,:)=buffer_candidates
    raw(nseed+2:nseed+2,:)=projector_candidates
    reference=matmul(transpose(raw),cache%candidate_compression)
    call export_dg_hybrid_fragment_coordinates(comm_fragment,fragment_id,first_generation,&
      cache%receipt%seed_fingerprint,cache%receipt%basis_fingerprint,grid_ids,&
      cache%local_row_layout_fingerprint,cache,q,seed,ok,message)
    call require_total(ok,'coordinate export failed: '//trim(message))
    call require_total(maxval(abs(matmul(transpose(cache%wannier_values),q)-reference))<1d-10,&
      'exported Q does not restore the pre-Wannier fixed frame')
    call require_total(maxval(abs(matmul(transpose(cache%wannier_values),seed)-&
      transpose(dc_seed_values)))<1d-10,'exported seed coefficients do not reconstruct DC orbitals')
    call test_seed_handoff(seed)
    call require_total(same_cache_payload(cache,first_snapshot),'coordinate export mutated cache')
    do case_id=1,5
      bad=cache
      select case(case_id)
      case(1)
        if(fragment_rank==0)bad%wannier_transform(1,1)=bad%wannier_transform(1,1)+0.1d0
      case(2)
        bad%receipt%basis_generation=first_generation+1
      case(3)
        bad%fragment_comm_size=fragment_size+1
      case(4)
        if(fragment_rank==0)deallocate(bad%local_grid_ids)
      case(5)
        bad%valid=.false.
      end select
      call export_dg_hybrid_fragment_coordinates(comm_fragment,fragment_id,first_generation,&
        cache%receipt%seed_fingerprint,cache%receipt%basis_fingerprint,grid_ids,&
        cache%local_row_layout_fingerprint,bad,q,seed,ok,message)
      call require_total(.not.ok,'coordinate export accepted invalid cache')
      call require_total(.not.allocated(q).and..not.allocated(seed),&
        'failed coordinate export published partial coordinates')
    enddo
    do case_id=1,4
      expected_seed=cache%receipt%seed_fingerprint
      expected_basis=cache%receipt%basis_fingerprint
      expected_layout=cache%local_row_layout_fingerprint
      expected_ids=grid_ids
      if(fragment_rank==0)then
        select case(case_id)
        case(1);expected_seed=ieor(expected_seed,1_int64)
        case(2);expected_basis=ieor(expected_basis,1_int64)
        case(3);expected_layout=ieor(expected_layout,1_int64)
        case(4);expected_ids(1)=expected_ids(1)+1_int64
        end select
      endif
      call export_dg_hybrid_fragment_coordinates(comm_fragment,fragment_id,first_generation,&
        expected_seed,expected_basis,expected_ids,expected_layout,cache,q,seed,ok,message)
      call require_total(.not.ok,'coordinate export ignored caller provenance/layout mismatch')
      call require_total(.not.allocated(q).and..not.allocated(seed),&
        'mismatched coordinate export published outputs')
    enddo
    call require_total(setup_calls==saved_setup.and.run_calls==saved_run,&
      'coordinate export called Wannier90 again')
  end subroutine test_coordinate_export

  subroutine test_seed_handoff(seed)
    complex(real64),intent(in)::seed(:,:)
    type(s_dg_hybrid_fragment_subspace_state)::initial
    complex(real64),allocatable::rows(:,:),gathered(:,:)
    integer(int64),allocatable::coefficient_ids(:)
    integer,allocatable::selected(:)
    real(real64)::occupations(nseed)
    integer::nwf,ntotal,nrows,a,b,p,code
    nwf=size(seed,1);ntotal=nwf+2
    nrows=count([(mod(a-1,min(fragment_size,2))==fragment_rank,a=1,ntotal)])
    allocate(coefficient_ids(nrows),rows(nrows,nseed));rows=0d0;p=0
    do a=ntotal,1,-1
      if(mod(a-1,min(fragment_size,2))/=fragment_rank)cycle
      p=p+1;coefficient_ids(p)=a
      if(a<=nwf)rows(p,:)=seed(a,:)
    enddo
    occupations=0d0;occupations(1)=2d0
    call initialize_dg_hybrid_fragment_subspace(comm_fragment,ntotal,coefficient_ids,&
      fragment_id,first_generation,cache%receipt%basis_fingerprint,901_int64,&
      rows,physical_energies,occupations,1,1d-8,1d-10,1d-10,fixture_identity_metric,&
      initial,selected,ok,message)
    call require_total(ok,'saved Wannier seed to CG handoff failed: '//trim(message))
    call require_total(initial%state_count==2.and.all(selected==[1,2]),'saved seed handoff selected named WFs')
    allocate(gathered(ntotal,2));gathered=0d0
    do a=1,nrows
      b=int(coefficient_ids(a));gathered(b,:)=initial%vectors(a,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,gathered,ntotal*2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_fragment,code)
    call require_total(code==MPI_SUCCESS,'initial CG coefficient gathering failed')
    call require_total(maxval(abs(matmul(transpose(cache%wannier_values),gathered(:nwf,:))-&
      transpose(dc_seed_values(selected,:))))<1d-10,'initial CG state changed the physical DC seed orbitals')
    call require_total(all(gathered(nwf+1:,:)==0d0),'initial DC seeds acquired PW components')
  end subroutine test_seed_handoff

  subroutine fixture_identity_metric(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    ! This fixture's constructed WFs and appended independent PW axes are orthonormal.
    output=input;valid=.true.
  end subroutine fixture_identity_metric
  subroutine fill_fixture
    real(real64)::angle
    integer::point,state,position
    grid_weights=1d0
    fractional_coordinates=0d0
    dc_seed_values=(0d0,0d0);buffer_candidates=(0d0,0d0)
    projector_candidates=(0d0,0d0);position=0
    do point=1,global_ngrid
      if(mod(point-1,fragment_size)/=fragment_rank)cycle
      position=position+1
      grid_ids(position)=int(1000*fragment_id+point,int64)
      fractional_coordinates(1,position)=real(point-1,real64)/real(global_ngrid,real64)
      do state=1,nseed
        if(point/=state)cycle
        angle=0.1d0*real(10*fragment_id+state,real64)
        dc_seed_values(state,position)=cmplx(cos(angle),sin(angle),real64)
      enddo
      if(point==buffer_global_point)&
        buffer_candidates(1,position)=cmplx(0.8d0,0.6d0,real64)
      if(point==projector_global_point)&
        projector_candidates(1,position)=cmplx(0.6d0,-0.8d0,real64)
    enddo
    do state=1,nseed
      physical_energies(state)=-2d0+0.25d0*real(state,real64)+0.1d0*fragment_id
      physical_occupations(state)=2d0-0.25d0*real(state-1,real64)
    enddo
    real_lattice=0d0;reciprocal_lattice=0d0
    do i=1,3
      real_lattice(i,i)=8d0
      reciprocal_lattice(i,i)=2d0*acos(-1d0)/8d0
    enddo
    atom_symbols(1)='H ';atoms_cart=0d0
  end subroutine fill_fixture

  integer function local_position(global_point)result(position)
    integer,intent(in)::global_point
    integer::point
    position=0
    do point=1,global_point
      if(mod(point-1,fragment_size)==fragment_rank)position=position+1
    enddo
    if(mod(global_point-1,fragment_size)/=fragment_rank)position=0
  end function local_position

  subroutine require_unique_fragment_rows
    integer::local_owner_count(global_ngrid),owner_count(global_ngrid),position,point
    logical::ordered
    local_owner_count=0;ordered=.true.
    do position=1,nlocal
      point=int(grid_ids(position)-int(1000*fragment_id,int64))
      if(point<1.or.point>global_ngrid)then
        ordered=.false.
      else
        local_owner_count(point)=local_owner_count(point)+1
        if(position>1)ordered=ordered.and.grid_ids(position)>grid_ids(position-1)
      endif
    enddo
    call MPI_Allreduce(local_owner_count,owner_count,global_ngrid,MPI_INTEGER,MPI_SUM,&
      comm_fragment,ierr)
    call require_total(ordered.and.all(owner_count==1),&
      'fixture grid rows are not uniquely distributed inside each fragment communicator')
  end subroutine require_unique_fragment_rows

  subroutine invoke_builder(generation,seed_values,buffer_values,target_cache,build_ok,build_message)
    integer,intent(in)::generation
    complex(real64),intent(in)::seed_values(:,:),buffer_values(:,:)
    type(s_dg_hybrid_fragment_wannier_cache),intent(inout)::target_cache
    logical,intent(out)::build_ok
    character(*),intent(out)::build_message
    call invoke_builder_layout(comm_total,comm_fragment,generation,grid_ids,grid_weights,&
      fractional_coordinates,seed_values,buffer_values,projector_candidates,target_cache,&
      build_ok,build_message)
  end subroutine invoke_builder

  subroutine invoke_builder_layout(total_communicator,fragment_communicator,generation,&
      layout_grid_ids,layout_grid_weights,layout_fractional,seed_values,buffer_values,&
      projector_values,target_cache,build_ok,build_message)
    integer,intent(in)::total_communicator,fragment_communicator,generation
    integer(int64),intent(in)::layout_grid_ids(:)
    real(real64),intent(in)::layout_grid_weights(:),layout_fractional(:,:)
    complex(real64),intent(in)::seed_values(:,:),buffer_values(:,:),projector_values(:,:)
    type(s_dg_hybrid_fragment_wannier_cache),intent(inout)::target_cache
    logical,intent(out)::build_ok
    character(*),intent(out)::build_message
    call build_dg_hybrid_fragment_wannier(comm_total=total_communicator,&
      comm_fragment=fragment_communicator,fragment_id=fragment_id,&
      basis_generation=generation,seed_directory=artifact_root,&
      grid_ids=layout_grid_ids,grid_weights=layout_grid_weights,dc_seed_values=seed_values,&
      dc_seed_energies=physical_energies,dc_seed_occupations=physical_occupations,&
      buffer_candidate_values=buffer_values,&
      projector_candidate_values=projector_values,&
      metric_tolerance=metric_tolerance,fragment_real_lattice=real_lattice,&
      fragment_reciprocal_lattice=reciprocal_lattice,atom_symbols=atom_symbols,&
      atoms_cart=atoms_cart,fractional_coordinates=layout_fractional,&
      num_iter=20,localization_tolerance=localization_tolerance,&
      coordinator_byte_limit=10000000_int64,cache=target_cache,&
      ok=build_ok,message=build_message)
  end subroutine invoke_builder_layout

  subroutine require_corrupted_cache_rejected(target_cache,description)
    type(s_dg_hybrid_fragment_wannier_cache),intent(inout)::target_cache
    character(*),intent(in)::description
    integer::saved_setup_calls,saved_run_calls
    saved_setup_calls=setup_calls;saved_run_calls=run_calls
    call invoke_builder(first_generation,dc_seed_values,buffer_candidates,target_cache,ok,message)
    call require_total(.not.ok,trim(description)//' was accepted for cache reuse')
    call require_same_message(message,trim(description)//' did not fail collectively')
    call require_total(index(message,'cache')>0.and.&
      (index(message,'payload')>0.or.index(message,'integrity')>0),&
      trim(description)//' message does not identify cache integrity')
    call require_total(setup_calls==saved_setup_calls.and.run_calls==saved_run_calls,&
      trim(description)//' entered Wannier90')
  end subroutine require_corrupted_cache_rejected

  subroutine swap_seed_grid_rows(seed_values,left_point,right_point)
    complex(real64),intent(inout)::seed_values(:,:)
    integer,intent(in)::left_point,right_point
    complex(real64)::left_values(size(seed_values,1)),right_values(size(seed_values,1))
    integer::position
    left_values=(0d0,0d0);right_values=(0d0,0d0)
    position=local_position(left_point)
    if(position>0)left_values=seed_values(:,position)
    position=local_position(right_point)
    if(position>0)right_values=seed_values(:,position)
    call MPI_Allreduce(MPI_IN_PLACE,left_values,size(left_values),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm_fragment,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,right_values,size(right_values),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm_fragment,ierr)
    position=local_position(left_point)
    if(position>0)seed_values(:,position)=right_values
    position=local_position(right_point)
    if(position>0)seed_values(:,position)=left_values
  end subroutine swap_seed_grid_rows

  subroutine validate_metric_null_cache(target_cache,generation)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::target_cache
    integer,intent(in)::generation
    complex(real64),allocatable::reconstructed(:,:)
    real(real64)::gram_defect,buffer_defect
    call require_total(target_cache%valid,'metric-null candidate cache is invalid')
    call require_total(target_cache%receipt%basis_generation==generation.and.&
      target_cache%receipt%candidate_rank==ncandidate.and.&
      target_cache%receipt%retained_rank==ncandidate-1,&
      'metric-null candidate did not remove exactly one auxiliary null mode')
    call require_total(all(shape(target_cache%wannier_values)==[ncandidate-1,nlocal]).and.&
      all(shape(target_cache%candidate_compression)==[ncandidate,ncandidate-1]).and.&
      all(shape(target_cache%wannier_transform)==[ncandidate-1,ncandidate-1]).and.&
      all(shape(target_cache%dc_seed_coefficients_in_wannier)==[ncandidate-1,nseed]),&
      'metric-null candidate cache has the wrong state-major layout')
    allocate(reconstructed(nseed,nlocal))
    reconstructed=matmul(transpose(target_cache%dc_seed_coefficients_in_wannier),&
      target_cache%wannier_values)
    call require_total(maxval(abs(reconstructed-dc_seed_values))<1d-11.and.&
      target_cache%receipt%seed_reconstruction_defect<1d-11,&
      'metric-null candidate build did not preserve all physical DC seeds')
    call compute_metric_span_defect(buffer_candidates(1,:),target_cache%wannier_values,buffer_defect)
    call compute_metric_gram_defect(target_cache%wannier_values,gram_defect)
    call require_total(buffer_defect<1d-11.and.gram_defect<1d-11,&
      'metric-null auxiliary candidate span or metric certificate failed')
  end subroutine validate_metric_null_cache

  subroutine validate_successful_cache(target_cache,generation)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::target_cache
    integer,intent(in)::generation
    complex(real64),allocatable::reconstructed(:,:),global_seed(:,:),global_reconstructed(:,:),&
      projector_before(:,:),projector_after(:,:)
    real(real64),allocatable::density_before(:),density_after(:)
    real(real64)::buffer_defect,projector_defect,gram_defect,buffer_center
    integer::point,position
    call require_total(target_cache%valid,'successful build did not mark its cache valid')
    call require_total(target_cache%receipt%fragment_id==fragment_id.and.&
      target_cache%receipt%basis_generation==generation,'fragment/generation receipt mismatch')
    call require_total(target_cache%receipt%candidate_rank==ncandidate.and.&
      target_cache%receipt%retained_rank==ncandidate,&
      'an independent local candidate, including the neighbor-buffer column, was dropped')
    call require_total(target_cache%receipt%setup_count==1.and.&
      target_cache%receipt%run_count==1,'per-generation Wannier90 receipt count is not one')
    call require_total(target_cache%receipt%seed_fingerprint/=0_int64.and.&
      target_cache%receipt%basis_fingerprint/=0_int64.and.&
      target_cache%receipt%transform_fingerprint/=0_int64.and.&
      target_cache%receipt%replicated_payload_fingerprint/=0_int64.and.&
      target_cache%receipt%distributed_wannier_fingerprint/=0_int64,&
      'fragment Wannier receipt has an empty fingerprint')
    call require_total(allocated(target_cache%wannier_values).and.&
      allocated(target_cache%candidate_compression).and.&
      allocated(target_cache%wannier_transform).and.&
      allocated(target_cache%dc_seed_coefficients_in_wannier).and.&
      allocated(target_cache%physical_dc_seed_energies).and.&
      allocated(target_cache%physical_dc_seed_occupations),&
      'successful fragment cache is incomplete')
    call require_total(all(shape(target_cache%wannier_values)==[ncandidate,nlocal]).and.&
      all(shape(target_cache%candidate_compression)==[ncandidate,ncandidate]).and.&
      all(shape(target_cache%wannier_transform)==[ncandidate,ncandidate]).and.&
      all(shape(target_cache%dc_seed_coefficients_in_wannier)==[ncandidate,nseed]),&
      'fragment cache array layout is not state-major')
    call require_total(ncandidate>size(physical_energies),&
      'fixture does not separate retained candidate count from physical seed energies')
    call require_total(bitwise_real_equal(target_cache%physical_dc_seed_energies,physical_energies).and.&
      bitwise_real_equal(target_cache%physical_dc_seed_occupations,physical_occupations),&
      'auxiliary Wannier90 labels replaced physical energies or occupations')
    call require_complex_same_fragment(target_cache%candidate_compression,&
      'candidate compression differs bitwise within a fragment communicator')
    call require_complex_same_fragment(target_cache%wannier_transform,&
      'Wannier transform differs bitwise within a fragment communicator')
    call require_complex_same_fragment(target_cache%dc_seed_coefficients_in_wannier,&
      'seed-to-WF coefficients differ bitwise within a fragment communicator')
    call require_real_same_fragment(target_cache%physical_dc_seed_energies,&
      'physical energies differ bitwise within a fragment communicator')
    call require_real_same_fragment(target_cache%physical_dc_seed_occupations,&
      'physical occupations differ bitwise within a fragment communicator')
    call require_receipt_same_fragment(target_cache)
    call validate_fragment_local_canonical_gauge(target_cache)
    allocate(reconstructed(nseed,nlocal))
    reconstructed=matmul(transpose(target_cache%dc_seed_coefficients_in_wannier),&
      target_cache%wannier_values)
    call require_total(maxval(abs(reconstructed-dc_seed_values))<1d-11.and.&
      target_cache%receipt%seed_reconstruction_defect<1d-11,&
      'returned seed-to-WF coefficients do not reconstruct every DC seed')
    allocate(global_seed(nseed,global_ngrid),global_reconstructed(nseed,global_ngrid))
    global_seed=(0d0,0d0);global_reconstructed=(0d0,0d0)
    do position=1,nlocal
      point=int(grid_ids(position)-int(1000*fragment_id,int64))
      global_seed(:,point)=dc_seed_values(:,position)
      global_reconstructed(:,point)=reconstructed(:,position)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_seed,size(global_seed),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm_fragment,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,global_reconstructed,size(global_reconstructed),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm_fragment,ierr)
    allocate(projector_before(global_ngrid,global_ngrid),&
      projector_after(global_ngrid,global_ngrid))
    projector_before=matmul(conjg(transpose(global_seed)),global_seed)
    projector_after=matmul(conjg(transpose(global_reconstructed)),global_reconstructed)
    call require_total(maxval(abs(projector_before-projector_after))<1d-11,&
      'occupied DC seed projector changed under fragment localization')
    allocate(density_before(global_ngrid),density_after(global_ngrid))
    do point=1,global_ngrid
      density_before(point)=sum(physical_occupations*abs(global_seed(:,point))**2)
      density_after(point)=sum(physical_occupations*abs(global_reconstructed(:,point))**2)
    enddo
    call require_total(maxval(abs(density_before-density_after))<1d-11,&
      'occupation-weighted DC seed density changed under fragment localization')
    call compute_metric_span_defect(buffer_candidates(1,:),target_cache%wannier_values,buffer_defect)
    call compute_metric_span_defect(projector_candidates(1,:),target_cache%wannier_values,projector_defect)
    call compute_metric_gram_defect(target_cache%wannier_values,gram_defect)
    buffer_center=0d0;position=local_position(buffer_global_point)
    if(position>0)buffer_center=fractional_coordinates(1,position)
    call MPI_Allreduce(MPI_IN_PLACE,buffer_center,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
      comm_fragment,ierr)
    call require_total(abs(buffer_center-0.875d0)<1d-15.and.buffer_defect<1d-11,&
      'neighbor-buffer-centered candidate was not retained in the fragment span')
    call require_total(projector_defect<1d-11.and.gram_defect<1d-11,&
      'projector candidate or fragment Wannier metric certificate failed')
  end subroutine validate_successful_cache

  subroutine validate_fragment_local_canonical_gauge(target_cache)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::target_cache
    complex(real64),allocatable::raw_candidates(:,:),expected_values(:,:),expected_transform(:,:)
    real(real64),allocatable::expected_centers(:,:)
    real(real64)::sine,scale,gauge_marker
    integer::state,setup_count_before,run_count_before
    logical::expected_ok
    character(512)::expected_message
    allocate(raw_candidates(ncandidate,nlocal),expected_values(ncandidate,nlocal))
    raw_candidates(1:nseed,:)=dc_seed_values
    raw_candidates(nseed+1,:)=buffer_candidates(1,:)
    raw_candidates(nseed+2,:)=projector_candidates(1,:)
    expected_values=matmul(transpose(target_cache%candidate_compression),raw_candidates)
    allocate(expected_transform(ncandidate,ncandidate),expected_centers(3,ncandidate))
    expected_transform=(0d0,0d0);expected_centers=0d0
    expected_centers(2,:)=-epsilon(1d0)/4d0
    do state=1,ncandidate
      expected_transform(state,state)=(1d0,0d0)
      expected_centers(1,state)=1d0+0.1d0*real(ncandidate-state,real64)
    enddo
    sine=merge(0.6d0,-0.6d0,fragment_id==1)
    expected_transform(1,1)=0.8d0;expected_transform(2,1)=sine
    expected_transform(1,2)=-sine;expected_transform(2,2)=0.8d0
    expected_centers(1,ncandidate)=-0.125d0
    setup_count_before=setup_calls;run_count_before=run_calls
    call apply_dg_w90_gamma_transform(comm_fragment,grid_ids,expected_values,&
      transform=expected_transform,centers=expected_centers,&
      tolerance=localization_tolerance,ok=expected_ok,message=expected_message)
    call require_total(expected_ok,&
      'fragment-local canonical gauge oracle failed: '//trim(expected_message))
    call require_total(allocated(target_cache%centers_fractional),'cache did not retain WF centers')
    call require_total(all(shape(target_cache%centers_fractional)==[3,ncandidate]),'wrong center shape')
    expected_centers=modulo(expected_centers,1d0)
    where(expected_centers==1d0)expected_centers=0d0
    call require_total(maxval(abs(target_cache%centers_fractional-expected_centers))<1d-13,&
      'cached centers do not match final reordered transform')
    call require_total(all(target_cache%centers_fractional>=0d0).and.&
      all(target_cache%centers_fractional<1d0),'cached centers are not periodically wrapped')
    call require_total(setup_calls==setup_count_before.and.run_calls==run_count_before,&
      'fragment-local gauge oracle entered Wannier90 setup or run')
    scale=max(1d0,maxval(abs(expected_transform)),maxval(abs(target_cache%wannier_transform)))
    call require_total(maxval(abs(target_cache%wannier_transform-expected_transform))<=&
      256d0*epsilon(1d0)*scale,&
      'cache transform contains a cross-fragment gauge alignment')
    scale=max(1d0,maxval(abs(expected_values)),maxval(abs(target_cache%wannier_values)))
    call require_total(maxval(abs(target_cache%wannier_values-expected_values))<=&
      256d0*epsilon(1d0)*scale,&
      'local cache Wannier values contain a cross-fragment gauge alignment')
    gauge_marker=real(target_cache%wannier_transform(2,ncandidate-1)*&
      conjg(target_cache%wannier_transform(1,ncandidate-1)),real64)
    call require_total(merge(gauge_marker>0d0,gauge_marker<0d0,fragment_id==1),&
      'fragment-specific signed stub gauge was collapsed across fragments')
  end subroutine validate_fragment_local_canonical_gauge

  subroutine compute_metric_span_defect(vector,basis,defect)
    complex(real64),intent(in)::vector(:),basis(:,:)
    real(real64),intent(out)::defect
    complex(real64),allocatable::coefficients(:),reconstruction(:)
    real(real64)::norm2
    integer::state
    allocate(coefficients(size(basis,1)),reconstruction(size(vector)))
    do state=1,size(basis,1)
      coefficients(state)=sum(grid_weights*conjg(basis(state,:))*vector)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,coefficients,size(coefficients),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm_fragment,ierr)
    reconstruction=matmul(coefficients,basis)
    norm2=sum(grid_weights*abs(vector-reconstruction)**2)
    call MPI_Allreduce(MPI_IN_PLACE,norm2,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm_fragment,ierr)
    defect=sqrt(max(0d0,norm2))
  end subroutine compute_metric_span_defect

  subroutine compute_metric_gram_defect(basis,defect)
    complex(real64),intent(in)::basis(:,:)
    real(real64),intent(out)::defect
    complex(real64),allocatable::gram(:,:)
    integer::left,right
    allocate(gram(size(basis,1),size(basis,1)));gram=(0d0,0d0)
    do right=1,size(basis,1);do left=1,size(basis,1)
      gram(left,right)=sum(grid_weights*conjg(basis(left,:))*basis(right,:))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,size(gram),MPI_DOUBLE_COMPLEX,MPI_SUM,&
      comm_fragment,ierr)
    do left=1,size(basis,1);gram(left,left)=gram(left,left)-1d0;enddo
    defect=maxval(abs(gram))
  end subroutine compute_metric_gram_defect

  subroutine require_stub_call(call_index,generation,expected_band_count)
    integer,intent(in)::call_index,generation
    integer,intent(in),optional::expected_band_count
    character(1024)::expected_prefix
    integer::band_count
    logical::condition
    band_count=ncandidate;if(present(expected_band_count))band_count=expected_band_count
    condition=setup_calls==merge(call_index,0,fragment_rank==0).and.&
      run_calls==merge(call_index,0,fragment_rank==0)
    if(fragment_rank==0)then
      write(expected_prefix,'(a,"/fragment-",i6.6,"/generation-",i8.8,"/")')&
        artifact_root,fragment_id,generation
      condition=condition.and.index(setup_seed_history(call_index),trim(expected_prefix))==1
      condition=condition.and.setup_seed_history(call_index)==run_seed_history(call_index)
      condition=condition.and.run_band_count_history(call_index)==band_count
      condition=condition.and.run_zero_auxiliary_energies(call_index)
      condition=condition.and..not.setup_saw_dmn(call_index)
      condition=condition.and.setup_saw_site_false(call_index)
      condition=condition.and..not.setup_saw_site_true(call_index)
      condition=condition.and..not.setup_saw_symmetrize(call_index)
      condition=condition.and..not.setup_saw_foreign_fragment(call_index)
    endif
    call require_total(condition,&
      'Wannier90 root count, namespace, or unconstrained fragment contract failed')
  end subroutine require_stub_call

  subroutine require_two_distinct_fragment_seeds(call_index)
    integer,intent(in)::call_index
    character(1024)::local_seed
    character(1024),allocatable::all_seeds(:)
    logical::condition
    integer::rank_index,nonempty
    local_seed='';if(fragment_rank==0)local_seed=setup_seed_history(call_index)
    allocate(all_seeds(total_size));all_seeds=''
    call MPI_Gather(local_seed,len(local_seed),MPI_CHARACTER,all_seeds,len(local_seed),&
      MPI_CHARACTER,0,comm_total,ierr)
    condition=.true.
    if(total_rank==0)then
      nonempty=0
      do rank_index=1,total_size
        if(len_trim(all_seeds(rank_index))>0)nonempty=nonempty+1
      enddo
      condition=nonempty==2.and.trim(all_seeds(1))/=trim(all_seeds(2))
    endif
    call require_total(condition,'fragment artifact seeds are not exactly two distinct namespaces')
  end subroutine require_two_distinct_fragment_seeds

  subroutine require_complex_same_fragment(values,text)
    complex(real64),intent(in)::values(:,:)
    character(*),intent(in)::text
    integer(int64),allocatable::bits(:),root_bits(:)
    allocate(bits(2*size(values)),root_bits(2*size(values)))
    bits=transfer(values,bits);root_bits=bits
    call MPI_Bcast(root_bits,size(root_bits),MPI_INTEGER8,0,comm_fragment,ierr)
    call require_total(all(bits==root_bits),text)
  end subroutine require_complex_same_fragment

  subroutine require_real_same_fragment(values,text)
    real(real64),intent(in)::values(:)
    character(*),intent(in)::text
    integer(int64),allocatable::bits(:),root_bits(:)
    allocate(bits(size(values)),root_bits(size(values)))
    bits=transfer(values,bits);root_bits=bits
    call MPI_Bcast(root_bits,size(root_bits),MPI_INTEGER8,0,comm_fragment,ierr)
    call require_total(all(bits==root_bits),text)
  end subroutine require_real_same_fragment

  subroutine require_receipt_same_fragment(target_cache)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::target_cache
    integer::ints(6),root_ints(6)
    integer(int64)::fingerprints(6),root_fingerprints(6)
    ints=[target_cache%receipt%fragment_id,target_cache%receipt%basis_generation,&
      target_cache%receipt%candidate_rank,target_cache%receipt%retained_rank,&
      target_cache%receipt%setup_count,target_cache%receipt%run_count]
    fingerprints=[target_cache%receipt%seed_fingerprint,target_cache%receipt%basis_fingerprint,&
      target_cache%receipt%transform_fingerprint,&
      target_cache%receipt%replicated_payload_fingerprint,&
      target_cache%receipt%distributed_wannier_fingerprint,&
      transfer(target_cache%receipt%seed_reconstruction_defect,0_int64)]
    root_ints=ints;root_fingerprints=fingerprints
    call MPI_Bcast(root_ints,size(root_ints),MPI_INTEGER,0,comm_fragment,ierr)
    call MPI_Bcast(root_fingerprints,size(root_fingerprints),MPI_INTEGER8,0,comm_fragment,ierr)
    call require_total(all(ints==root_ints).and.all(fingerprints==root_fingerprints),&
      'fragment receipt differs bitwise within its communicator')
  end subroutine require_receipt_same_fragment

  subroutine require_same_message(value,text)
    character(*),intent(in)::value,text
    character(len(value))::root_value
    root_value=value
    call MPI_Bcast(root_value,len(root_value),MPI_CHARACTER,0,comm_total,ierr)
    call require_total(value==root_value,text)
  end subroutine require_same_message

  subroutine mock_scf_loop(iteration_limit,readonly_cache,iterations,construction_calls)
    integer,intent(in)::iteration_limit
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::readonly_cache
    integer,intent(out)::iterations,construction_calls
    real(real64)::physical_energy_sum
    iterations=0;construction_calls=0;physical_energy_sum=0d0
    do while(iterations<iteration_limit)
      iterations=iterations+1
      physical_energy_sum=physical_energy_sum+sum(readonly_cache%physical_dc_seed_energies)
    enddo
    if(.not.readonly_cache%valid.or.physical_energy_sum>=huge(1d0))construction_calls=-1
  end subroutine mock_scf_loop

  logical function cache_is_unpublished(target_cache)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::target_cache
    cache_is_unpublished=.not.target_cache%valid.and.&
      .not.allocated(target_cache%local_grid_ids).and.&
      .not.allocated(target_cache%wannier_values).and.&
      .not.allocated(target_cache%candidate_compression).and.&
      .not.allocated(target_cache%wannier_transform).and.&
      .not.allocated(target_cache%dc_seed_coefficients_in_wannier).and.&
      .not.allocated(target_cache%physical_dc_seed_energies).and.&
      .not.allocated(target_cache%physical_dc_seed_occupations)
  end function cache_is_unpublished

  logical function same_cache_payload(left,right)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::left,right
    same_cache_payload=left%valid.eqv.right%valid
    same_cache_payload=same_cache_payload.and.&
      left%receipt%fragment_id==right%receipt%fragment_id.and.&
      left%receipt%basis_generation==right%receipt%basis_generation.and.&
      left%receipt%candidate_rank==right%receipt%candidate_rank.and.&
      left%receipt%retained_rank==right%receipt%retained_rank.and.&
      left%receipt%setup_count==right%receipt%setup_count.and.&
      left%receipt%run_count==right%receipt%run_count.and.&
      left%receipt%seed_fingerprint==right%receipt%seed_fingerprint.and.&
      left%receipt%basis_fingerprint==right%receipt%basis_fingerprint.and.&
      left%receipt%transform_fingerprint==right%receipt%transform_fingerprint.and.&
      left%receipt%replicated_payload_fingerprint==&
      right%receipt%replicated_payload_fingerprint.and.&
      left%receipt%distributed_wannier_fingerprint==&
      right%receipt%distributed_wannier_fingerprint.and.&
      transfer(left%receipt%seed_reconstruction_defect,0_int64)==&
      transfer(right%receipt%seed_reconstruction_defect,0_int64)
    same_cache_payload=same_cache_payload.and.&
      left%fragment_comm_rank==right%fragment_comm_rank.and.&
      left%fragment_comm_size==right%fragment_comm_size.and.&
      left%local_row_layout_fingerprint==right%local_row_layout_fingerprint.and.&
      all(left%local_grid_ids==right%local_grid_ids)
    same_cache_payload=same_cache_payload.and.&
      bitwise_real_equal(reshape(left%centers_fractional,[size(left%centers_fractional)]),&
        reshape(right%centers_fractional,[size(right%centers_fractional)])).and.&
      bitwise_complex_equal(left%wannier_values,right%wannier_values).and.&
      bitwise_complex_equal(left%candidate_compression,right%candidate_compression).and.&
      bitwise_complex_equal(left%wannier_transform,right%wannier_transform).and.&
      bitwise_complex_equal(left%dc_seed_coefficients_in_wannier,&
        right%dc_seed_coefficients_in_wannier).and.&
      bitwise_real_equal(left%physical_dc_seed_energies,right%physical_dc_seed_energies).and.&
      bitwise_real_equal(left%physical_dc_seed_occupations,right%physical_dc_seed_occupations)
  end function same_cache_payload

  logical function bitwise_complex_equal(left,right)
    complex(real64),intent(in)::left(:,:),right(:,:)
    integer(int64),allocatable::left_bits(:),right_bits(:)
    bitwise_complex_equal=.false.
    if(any(shape(left)/=shape(right)))return
    allocate(left_bits(2*size(left)),right_bits(2*size(right)))
    left_bits=transfer(left,left_bits);right_bits=transfer(right,right_bits)
    bitwise_complex_equal=all(left_bits==right_bits)
  end function bitwise_complex_equal

  logical function bitwise_real_equal(left,right)
    real(real64),intent(in)::left(:),right(:)
    integer(int64),allocatable::left_bits(:),right_bits(:)
    bitwise_real_equal=.false.
    if(size(left)/=size(right))return
    allocate(left_bits(size(left)),right_bits(size(right)))
    left_bits=transfer(left,left_bits);right_bits=transfer(right,right_bits)
    bitwise_real_equal=all(left_bits==right_bits)
  end function bitwise_real_equal

  subroutine require_total(condition,text,communicator)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer,intent(in),optional::communicator
    logical::global_condition
    integer::active_communicator,active_rank
    active_communicator=comm_total;if(present(communicator))active_communicator=communicator
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,active_communicator,ierr)
    if(.not.global_condition)then
      call MPI_Comm_rank(active_communicator,active_rank,ierr)
      if(active_rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(active_communicator,1,ierr)
    endif
  end subroutine require_total
end program test_dg_hybrid_fragment_wannier_mpi
#endif

#ifdef W90_TEST_STUBS
subroutine wannier_setup(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
    recip_lattice_loc,kpt_latt_loc,num_bands_tot,num_atoms_loc,atom_symbols_loc,&
    atoms_cart_loc,gamma_only_loc,spinors_loc,nntot_loc,nnlist_loc,nncell_loc,&
    num_bands_loc,num_wann_loc,proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,&
    proj_z_loc,proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)
  use dg_hybrid_fragment_wannier_test_stubs
  implicit none
  integer,parameter::num_nnmax=12
  character(*),intent(in)::seed_name
  integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_tot,num_atoms_loc
  real(8),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),&
    kpt_latt_loc(3,num_kpts_loc),atoms_cart_loc(3,num_atoms_loc)
  character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
  logical,intent(in)::gamma_only_loc,spinors_loc
  integer,intent(out)::nntot_loc,nnlist_loc(num_kpts_loc,num_nnmax),&
    nncell_loc(3,num_kpts_loc,num_nnmax),num_bands_loc,num_wann_loc
  real(8),intent(out)::proj_site_loc(3,num_bands_tot),proj_z_loc(3,num_bands_tot),&
    proj_x_loc(3,num_bands_tot),proj_zona_loc(num_bands_tot)
  integer,intent(out)::proj_l_loc(num_bands_tot),proj_m_loc(num_bands_tot),&
    proj_radial_loc(num_bands_tot),exclude_bands_loc(num_bands_tot),proj_s_loc(num_bands_tot)
  real(8),intent(out)::proj_s_qaxis_loc(3,num_bands_tot)
  integer::unit,io
  character(512)::line
  setup_calls=setup_calls+1
  if(setup_calls<=history_limit)then
    setup_seed_history(setup_calls)=seed_name
    inquire(file=trim(seed_name)//'.dmn',exist=setup_saw_dmn(setup_calls))
    setup_saw_foreign_fragment(setup_calls)=expected_fragment_id>0.and.&
      index(seed_name,trim(fragment_token(expected_fragment_id)))==0
    open(newunit=unit,file=trim(seed_name)//'.win',status='old',action='read',iostat=io)
    if(io==0)then
      do
        read(unit,'(a)',iostat=io)line
        if(io/=0)exit
        if(index(adjustl(line),'site_symmetry = .true.')==1)&
          setup_saw_site_true(setup_calls)=.true.
        if(index(adjustl(line),'site_symmetry = .false.')==1)&
          setup_saw_site_false(setup_calls)=.true.
        if(index(adjustl(line),'symmetrize_eps')==1)&
          setup_saw_symmetrize(setup_calls)=.true.
      enddo
      close(unit)
    endif
  endif
  nntot_loc=1;nnlist_loc=1;nncell_loc=0
  num_bands_loc=num_bands_tot;num_wann_loc=num_bands_tot
  proj_site_loc=0d0;proj_l_loc=0;proj_m_loc=0;proj_radial_loc=0
  proj_z_loc=0d0;proj_x_loc=0d0;proj_zona_loc=0d0;exclude_bands_loc=0
  proj_s_loc=0;proj_s_qaxis_loc=0d0
end subroutine wannier_setup

subroutine wannier_run(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
    recip_lattice_loc,kpt_latt_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc,&
    atom_symbols_loc,atoms_cart_loc,gamma_only_loc,m_matrix_loc,a_matrix_loc,&
    eigenvalues_loc,u_matrix_loc,u_matrix_opt_loc,lwindow_loc,wann_centres_loc,&
    wann_spreads_loc,spread_loc)
  use dg_hybrid_fragment_wannier_test_stubs
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  character(*),intent(in)::seed_name
  integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc
  real(8),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),&
    kpt_latt_loc(3,num_kpts_loc),atoms_cart_loc(3,num_atoms_loc)
  character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
  logical,intent(in)::gamma_only_loc
  complex(8),intent(in)::m_matrix_loc(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc),&
    a_matrix_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
  real(8),intent(in)::eigenvalues_loc(num_bands_loc,num_kpts_loc)
  complex(8),intent(out)::u_matrix_loc(num_wann_loc,num_wann_loc,num_kpts_loc),&
    u_matrix_opt_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
  logical,intent(out)::lwindow_loc(num_bands_loc,num_kpts_loc)
  real(8),intent(out)::wann_centres_loc(3,num_wann_loc),wann_spreads_loc(num_wann_loc),spread_loc(3)
  integer::i,k,unit,io
  real(8)::sine
  run_calls=run_calls+1
  if(run_calls<=history_limit)then
    run_seed_history(run_calls)=seed_name
    run_band_count_history(run_calls)=num_bands_loc
    run_zero_auxiliary_energies(run_calls)=num_kpts_loc==1.and.&
      num_bands_loc==num_wann_loc.and.&
      all(ieee_is_finite(eigenvalues_loc)).and.all(eigenvalues_loc==0d0)
  endif
  u_matrix_loc=(0d0,0d0);u_matrix_opt_loc=(0d0,0d0)
  do k=1,num_kpts_loc;do i=1,num_wann_loc
    u_matrix_loc(i,i,k)=(1d0,0d0);u_matrix_opt_loc(i,i,k)=(1d0,0d0)
  enddo;enddo
  if(num_wann_loc>=2)then
    sine=merge(0.6d0,-0.6d0,expected_fragment_id==1)
    do k=1,num_kpts_loc
      u_matrix_loc(1,1,k)=0.8d0;u_matrix_loc(2,1,k)=sine
      u_matrix_loc(1,2,k)=-sine;u_matrix_loc(2,2,k)=0.8d0
    enddo
  endif
  lwindow_loc=.true.;wann_centres_loc=0d0;wann_spreads_loc=0d0;spread_loc=0d0
  wann_centres_loc(2,:)=-epsilon(1d0)/4d0*8d0*0.52917721067d0
  do i=1,num_wann_loc
    wann_centres_loc(1,i)=(1d0+0.1d0*real(num_wann_loc-i,8))*8d0*0.52917721067d0
  enddo
  wann_centres_loc(1,num_wann_loc)=-0.125d0*8d0*0.52917721067d0
  if(override_test_centers)wann_centres_loc(1,:)=0.9375d0*8d0*0.52917721067d0
  open(newunit=unit,file=trim(seed_name)//'.wout',status='replace',action='write',iostat=io)
  if(io==0)then
    if(run_must_fail(seed_name))then
      write(unit,'(a)')'Maximum number of Wannier iterations reached without convergence'
    else
      write(unit,'(a)')'      1  -0.100E-13  0.0  1.0  0.0 <-- CONV'
      write(unit,'(a)')'             <<< Wannierisation convergence criteria satisfied >>>'
      write(unit,'(a)')' Final State'
      write(unit,'(a)')' All done: wannier90 exiting'
    endif
    close(unit)
  endif
end subroutine wannier_run
#endif
