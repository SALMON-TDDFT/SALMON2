#include "config.h"
program test_rt_dg_hybrid_initialization_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_hybrid_checkpoint,only:s_rt_dg_hybrid_ground_state_payload,&
    write_rt_dg_hybrid_ground_state_checkpoint,read_rt_dg_hybrid_ground_state_checkpoint_coalesced,&
    fingerprint_rt_dg_hybrid_component,rt_dg_hybrid_ground_state_checkpoint_version,&
    rt_dg_hybrid_energy_window_explicit
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,s_rt_dg_hybrid_v3_startup_receipt,&
    initialize_rt_dg_hybrid_from_checkpoint,validate_rt_dg_hybrid_v3_startup,&
    fingerprint_rt_dg_hybrid_scope,stamp_rt_dg_hybrid_v3_fingerprints
  use rt_dg_hybrid_density_update,only:reconstruct_rt_dg_hybrid_density,update_rt_dg_hybrid_density
  implicit none
  integer,parameter::construction_rank=3,certified_rank=2,occupied_rank=1,operation_count=2,grid_count=2
  integer(int64),parameter::cell_wrapped_position_fingerprint=int(z'43454C4C57524150',int64)
  real(real64),parameter::startup_tolerances(4)=[1d-11,1d-11,1d-11,1d-11]
  integer::comm,rank,nproc,ierr,smoke_grid_count
  integer(int64)::fingerprint
  logical::ok
  character(256)::message,path,mode,grid_count_argument
  type(s_rt_dg_hybrid_ground_state_payload)::payload
  type(s_rt_dg_hybrid_state)::state
  logical::force_callback_failure=.false.
  real(real64)::reference_density_total=0d0
  real(real64)::reference_shift=0d0
  complex(real64)::sparse_local_reference(3,3)=(0d0,0d0)

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call require(rt_dg_hybrid_ground_state_checkpoint_version==3,&
    'Hybrid RT initialization requires complete ground-state checkpoint version 3')
  call get_command_argument(1,path);call get_command_argument(2,mode)
  if(len_trim(mode)==0)mode='roundtrip'
  smoke_grid_count=grid_count
  call get_command_argument(3,grid_count_argument)
  if(len_trim(grid_count_argument)>0)read(grid_count_argument,*)smoke_grid_count

  select case(trim(mode))
  case('write_only')
    call build_v3_payload(payload)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'valid v3 fixture write failed: '//trim(message))
    if(rank==0)write(*,'(a,i0)')'HYBRID_RT_V3_FINGERPRINT=',fingerprint
  case('write_production')
    call build_v3_payload(payload,smoke_grid_count)
    call make_production_smoke_payload(payload)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'production v3 fixture write failed: '//trim(message))
    if(rank==0)write(*,'(a,i0)')'HYBRID_RT_V3_FINGERPRINT=',fingerprint
  case('roundtrip')
    call build_v3_payload(payload)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'valid v3 fixture write failed: '//trim(message))
    call verify_v3_startup(trim(path))
  case('write_negative_zero')
    call build_v3_payload(payload)
    call inject_negative_zero(payload)
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'negative-zero named fingerprint stamping failed: '//trim(message))
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'negative-zero v3 fixture write failed: '//trim(message))
  case('reject_named_fingerprint')
    call build_v3_payload(payload)
    payload%rt_space%basis_fingerprint=payload%rt_space%basis_fingerprint+1_int64
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_catalog_fingerprint')
    call build_v3_payload(payload)
    payload%construction_catalog%ids_fingerprint=payload%construction_catalog%ids_fingerprint+1_int64
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_catalog_semantics')
    call build_v3_payload(payload)
    payload%construction_catalog%ids=payload%construction_catalog%ids+100_int64
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'catalog-semantic named fingerprint stamping failed: '//trim(message))
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_top_fingerprint')
    call build_v3_payload(payload)
    payload%energy_fingerprint=payload%energy_fingerprint+1_int64
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_receipt_fingerprint')
    call build_v3_payload(payload)
    payload%energy_window%fingerprint=payload%energy_window%fingerprint+1_int64
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_receipt_semantics')
    call build_v3_payload(payload)
    payload%symmetry_receipt%scalar_covariance_defect=1d-13
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'receipt-semantic named fingerprint stamping failed: '//trim(message))
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_handoff_fingerprint')
    call build_v3_payload(payload)
    payload%handoff_receipts%fingerprint=payload%handoff_receipts%fingerprint+1_int64
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_position_convention_semantics')
    call build_v3_payload(payload)
    payload%position_convention_fingerprint=ieor(payload%position_convention_fingerprint,1_int64)
    payload%handoff_receipts%position_fingerprint=payload%position_convention_fingerprint
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'position-convention named fingerprint stamping failed: '//trim(message))
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_basis_projection')
    call build_v3_payload(payload)
    call perturb_rt_basis_nullspace(payload)
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'basis-tamper named fingerprint stamping failed: '//trim(message))
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_component_projection')
    call build_v3_payload(payload)
    call perturb_rt_fixed_components(payload)
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'component-tamper named fingerprint stamping failed: '//trim(message))
    call write_and_expect_startup_rejection(payload,trim(path))
  case('reject_orbital_physics','reject_metric_physics','reject_unitarity_physics',&
      'reject_density_physics','reject_electron_physics','reject_target_closure',&
      'reject_energy_covariance','reject_projector_covariance','reject_vector_covariance',&
      'reject_tensor_covariance')
    call build_v3_payload(payload)
    call perturb_startup_invariant(payload,trim(mode))
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    call require(ok,'physical-negative named fingerprint stamping failed: '//trim(message))
    call write_and_expect_invariant_rejection(payload,trim(path),trim(mode))
  case('density_tolerance_boundary')
    call build_v3_payload(payload)
    call exercise_density_tolerance_boundary(payload,trim(path))
  case('embedding_projector_tolerance_boundary')
    call build_v3_payload(payload)
    call exercise_embedding_projector_tolerance_boundary(payload,trim(path))
  case('read_only')
    call verify_v3_startup(trim(path))
  case('reject_only')
    call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(path),'tddft_response',.true.,1,.false.,.false.,.false.,&
      .false.,.false.,[1,0,0],startup_tolerances,state,ok,message)
    call require(.not.ok,'invalid, legacy, or tampered checkpoint was accepted by RT startup')
  case('tamper_receipt_after_redistribution','tamper_provenance_after_redistribution')
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'coalesced v3 fixture read failed before tamper: '//trim(message))
    if(trim(mode)=='tamper_receipt_after_redistribution')then
      payload%symmetry_receipt%maximum_physical_defect=&
        payload%symmetry_receipt%maximum_physical_defect+1d-4
    else
      payload%construction_catalog%provenance_fingerprint=&
        payload%construction_catalog%provenance_fingerprint+1_int64
    endif
    call expect_validator_rejection(payload,fingerprint)
  case default
    call require(.false.,'unknown Hybrid RT initialization test mode: '//trim(mode))
  end select

  if(rank==0)write(*,'(a,i0,a,a)')'PASS hybrid RT v3 initialization contract on ',nproc,' ranks mode=',trim(mode)
  call MPI_Finalize(ierr)
contains
  subroutine verify_v3_startup(checkpoint_path)
    character(*),intent(in)::checkpoint_path
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    integer::i,row,column,saved_certified_rank
    integer(int64)::structure_before,value_before
    logical::extent_ok,coefficient_ok,density_ok,position_ok,diagonal_found
    real(real64)::symmetry_defect,expected_position
    complex(real64),allocatable::hamiltonian_before(:)
    real(real64),allocatable::density_for_update(:)
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,checkpoint_path,payload,fingerprint,ok,message)
    call require(ok,'valid v3 coalesced read failed: '//trim(message))
    call validate_rt_dg_hybrid_v3_startup(comm,payload,fingerprint,startup_tolerances,receipt,ok,message)
    call require(ok,'valid v3 startup receipt failed: '//trim(message))
    call require(receipt%valid.and.receipt%certified_rank==certified_rank,&
      'startup receipt does not identify the certified RT rank')
    call require(receipt%payload_fingerprint==fingerprint,&
      'startup receipt lost the authenticated payload identity')
    call require(ieee_is_finite(receipt%orbital_residual).and.&
      receipt%orbital_residual<=startup_tolerances(1),'startup orbital residual was not recomputed')
    call require(max(receipt%metric_defect,receipt%embedding_defect,receipt%unitarity_defect,&
      receipt%basis_defect)<=startup_tolerances(1),&
      'startup certified metric/embedding/unitarity/basis projection was not recomputed')
    call require(receipt%density_defect<=startup_tolerances(2),&
      'startup density was not reconstructed in the localized RT basis')
    call require(receipt%electron_defect<=startup_tolerances(3),&
      'startup named electron count was not recomputed')
    symmetry_defect=max(receipt%target_closure_defect,receipt%energy_covariance_defect,&
      receipt%projector_defect,receipt%operator_component_defect,&
      receipt%fixed_operator_covariance_defect)
    call require(symmetry_defect<=startup_tolerances(4),&
      'startup certified symmetry and fixed-operator covariance was not recomputed')
    call require(any_rank(reshape(abs(payload%rt_space%tensor_operator_rows(:,:,1,2,1))>0.5d0,&
      [size(payload%rt_space%tensor_operator_rows(:,:,1,2,1))])),&
      'tensor-covariance fixture lacks a non-scalar odd Cartesian component')

    call initialize_rt_dg_hybrid_from_checkpoint(comm,checkpoint_path,'tddft_response',.true.,1,&
      .false.,.false.,.false.,.false.,.false.,[1,0,0],startup_tolerances,state,ok,message)
    call require(ok,'valid certified v3 RT initialization failed: '//trim(message))
    call require(state%valid.and.state%initial_invariants_valid,'certified v3 state lacks startup receipt')
    call require(payload%global_count>payload%rt_space%rank,&
      'fixture does not distinguish construction and certified ranks')
    call require(state%certified_rank==certified_rank.and.state%global_count==certified_rank.and.&
      state%noccupied==occupied_rank,&
      'RT state retained construction-only coefficient directions')
    call require(state%metric%global_count==certified_rank.and.state%operators%global_count==certified_rank,&
      'RT sparse operators retained construction rank')
    call require(size(state%coefficients,1)==size(state%owned_row_ids).and.&
      size(state%coefficients,2)==occupied_rank,'RT coefficient extent is not certified-rank owned')
    call require(size(state%kinetic_rows,2)==certified_rank.and.size(state%nonlocal_rows,2)==certified_rank.and.&
      size(state%local_rows)==size(state%operators%column_ids).and.size(state%sipg_rows,2)==certified_rank,&
      'RT fixed operators retain construction columns')
    call require(size(state%basis_values,1)==certified_rank.and.&
      size(state%basis_values,2)==size(state%grid_ids),'RT basis values are not the localized certified basis')
    call require(size(state%eigenvalues)==certified_rank,'RT eigenvalue extent is not certified rank')
    call require(all(state%eigenvalues==payload%certified_basis%certified_eigenvalues),&
      'RT eigenvalues fell back to the legacy occupied-only array')
    extent_ok=.true.
    if(size(state%metric%column_ids)>0)extent_ok=maxval(state%metric%column_ids)<=certified_rank
    call require(extent_ok,'RT metric graph references a construction-only column')
    extent_ok=.true.
    if(size(state%operators%column_ids)>0)extent_ok=maxval(state%operators%column_ids)<=certified_rank
    call require(extent_ok,'RT operator graph references a construction-only column')
    extent_ok=.true.;coefficient_ok=.true.
    do i=1,size(state%owned_row_ids)
      row=int(state%owned_row_ids(i))
      if(row<1.or.row>certified_rank)extent_ok=.false.
      if(maxval(abs(state%coefficients(i,:)-payload%certified_basis%initial_occupied_amplitudes(row,:)))>&
        startup_tolerances(1))coefficient_ok=.false.
    enddo
    call require(extent_ok,'RT owns a construction-only coefficient row')
    call require(coefficient_ok,'RT initial coefficient was not reconstructed from stored U_rt amplitudes')
    call require(state%payload_fingerprint==fingerprint,'RT state lost the authenticated v3 fingerprint')
    position_ok=.true.
    do i=1,size(state%owned_row_ids)
      row=int(state%owned_row_ids(i))
      diagonal_found=.false.
      do column=state%operators%row_offsets(i),state%operators%row_offsets(i+1)-1
        expected_position=0d0
        if(state%operators%column_ids(column)==row)then
          diagonal_found=.true.;expected_position=0.25d0+merge(1d0,-1d0,row==1)
        endif
        if(abs(state%operators%position_values(1,column)-expected_position)>startup_tolerances(4))position_ok=.false.
        if(any(abs(state%operators%position_values(2:3,column))>startup_tolerances(4)))position_ok=.false.
      enddo
      if(.not.diagonal_found)position_ok=.false.
    enddo
    call require(position_ok,'construction position was not projected as B_rt^H X B_rt')

    call reconstruct_rt_dg_hybrid_density(comm,state,ok,message)
    call require(ok,'certified RT density reconstruction failed: '//trim(message))
    call MPI_Allreduce(sum(state%density),reference_density_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call require(ierr==MPI_SUCCESS,'certified RT reference density reduction failed')
    density_ok=all(abs(state%density-payload%rt_space%density)<=startup_tolerances(2))
    call require(density_ok,'RT density was not reconstructed from certified localized basis amplitudes')
    saved_certified_rank=state%certified_rank
    if(rank==0)state%certified_rank=construction_rank
    call reconstruct_rt_dg_hybrid_density(comm,state,ok,message)
    call require(.not.ok,'density reconstruction accepted a rank-disagreeing construction extent')
    if(rank==0)state%certified_rank=saved_certified_rank

    structure_before=state%operator_structure_fingerprint
    value_before=state%operator_value_fingerprint
    allocate(hamiltonian_before,source=state%operators%hamiltonian_values)
    allocate(density_for_update,source=state%density)
    reference_shift=0.25d0
    call update_rt_dg_hybrid_density(comm,state,density_for_update,project_density_local,ok,message,&
      establish_fixed_density_reference=.true.)
    call require(ok,'certified RT local-potential update failed: '//trim(message))
    call require(all(state%operators%hamiltonian_values==hamiltonian_before),&
      'fixed-density reference did not preserve exact checkpoint H_rt(0)')
    call require(state%fixed_density_reference_valid,'fixed-density reference was not established on every rank')
    call require(any_rank([state%reference_refresh_defect>0d0]),&
      'material t=0 potential refresh did not establish a measured reference correction')
    state%density=state%density+0.125d0
    density_for_update=state%density
    call update_rt_dg_hybrid_density(comm,state,density_for_update,project_density_local,ok,message)
    call require(ok,'perturbed certified RT local-potential update failed: '//trim(message))
    call require(any_rank(state%operators%hamiltonian_values/=hamiltonian_before),&
      'certified density perturbation did not change H_rt')
    call require(state%operator_structure_fingerprint==structure_before.and.&
      state%operators%fingerprint==structure_before.and.state%operator_value_fingerprint/=value_before,&
      'certified density update rebuilt the graph or retained a stale value fingerprint')
    force_callback_failure=.true.
    call update_rt_dg_hybrid_density(comm,state,density_for_update,project_density_local,ok,message)
    call require(.not.ok,'one-rank local-potential callback failure was not rejected collectively')
    force_callback_failure=.false.
    call exercise_sparse_density_update
  end subroutine verify_v3_startup

  subroutine exercise_sparse_density_update
    type(s_rt_dg_hybrid_state)::sparse_state
    integer::owned_count,i,row,edge,nnz,local_degrees(3),global_degrees(3)
    complex(real64)::expected
    logical::sparse_values_ok
    real(real64)::sparse_density(1)
    owned_count=0
    do row=1,3;if(mod(row-1,nproc)==rank)owned_count=owned_count+1;enddo
    allocate(sparse_state%owned_row_ids(owned_count),sparse_state%coefficients(owned_count,1),&
      sparse_state%kinetic_rows(owned_count,3),sparse_state%nonlocal_rows(owned_count,3),&
      sparse_state%sipg_rows(owned_count,3))
    i=0
    do row=1,3
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;sparse_state%owned_row_ids(i)=row
    enddo
    sparse_state%coefficients=(0d0,0d0);sparse_state%kinetic_rows=(0d0,0d0)
    sparse_state%nonlocal_rows=(0d0,0d0);sparse_state%sipg_rows=(0d0,0d0)
    if(owned_count>0)then
      do i=1,owned_count
        row=int(sparse_state%owned_row_ids(i))
        if(row==1)then
          sparse_state%kinetic_rows(i,1)=1d0;sparse_state%kinetic_rows(i,2)=0.2d0
        else if(row==2)then
          sparse_state%kinetic_rows(i,1)=0.2d0;sparse_state%kinetic_rows(i,2)=2d0
        endif
      enddo
    endif
    allocate(sparse_state%metric%owned_row_ids(owned_count),sparse_state%metric%row_offsets(owned_count+1),&
      sparse_state%metric%column_ids(owned_count),sparse_state%metric%values(owned_count),&
      sparse_state%metric%active_rows(3),sparse_state%metric%packet_ids(3))
    sparse_state%metric%owned_row_ids=sparse_state%owned_row_ids;sparse_state%metric%row_offsets(1)=1
    do i=1,owned_count
      sparse_state%metric%column_ids(i)=int(sparse_state%owned_row_ids(i));sparse_state%metric%values(i)=1d0
      sparse_state%metric%row_offsets(i+1)=i+1
    enddo
    sparse_state%metric%active_rows=.true.;sparse_state%metric%packet_ids=1
    sparse_state%metric%global_count=3;sparse_state%metric%numerical_rank=3
    sparse_state%metric%max_row_nnz=1;sparse_state%metric%fingerprint=901_int64
    sparse_state%metric%valid=.true.
    nnz=0
    do i=1,owned_count
      row=int(sparse_state%owned_row_ids(i));if(row<3)nnz=nnz+2
    enddo
    allocate(sparse_state%operators%owned_row_ids(owned_count),&
      sparse_state%operators%row_offsets(owned_count+1),sparse_state%operators%column_ids(nnz),&
      sparse_state%operators%metric_values(nnz),sparse_state%operators%hamiltonian_values(nnz),&
      sparse_state%operators%position_values(3,nnz),sparse_state%local_rows(nnz))
    sparse_state%operators%owned_row_ids=sparse_state%owned_row_ids
    sparse_state%operators%row_offsets(1)=1;edge=0
    do i=1,owned_count
      row=int(sparse_state%owned_row_ids(i))
      if(row<3)then
        edge=edge+1;sparse_state%operators%column_ids(edge)=1
        edge=edge+1;sparse_state%operators%column_ids(edge)=2
      endif
      sparse_state%operators%row_offsets(i+1)=edge+1
    enddo
    sparse_state%operators%metric_values=(0d0,0d0);sparse_state%operators%position_values=(0d0,0d0)
    sparse_state%operators%hamiltonian_values=(0d0,0d0);sparse_state%local_rows=(0d0,0d0)
    sparse_state%operators%global_count=3;sparse_state%operators%metric_fingerprint=901_int64
    sparse_state%operators%fingerprint=902_int64;sparse_state%operators%valid=.true.
    allocate(sparse_state%grid_ids(1),sparse_state%grid_weights(1),sparse_state%density(1),&
      sparse_state%basis_values(3,1),sparse_state%occupations(1),sparse_state%eigenvalues(3))
    sparse_state%grid_ids=rank+1;sparse_state%grid_weights=1d0;sparse_state%density=1d0
    sparse_state%basis_values=(0d0,0d0);sparse_state%occupations=1d0;sparse_state%eigenvalues=0d0
    sparse_state%certified_rank=3;sparse_state%global_count=3;sparse_state%noccupied=1
    sparse_state%operator_structure_fingerprint=902_int64;sparse_state%valid=.true.
    sparse_local_reference=(0d0,0d0)
    sparse_local_reference(1,1)=0.3d0;sparse_local_reference(1,2)=cmplx(0.1d0,0.05d0,real64)
    sparse_local_reference(2,1)=conjg(sparse_local_reference(1,2));sparse_local_reference(2,2)=0.4d0
    sparse_density=1d0
    call update_rt_dg_hybrid_density(comm,sparse_state,sparse_density,project_sparse_density_local,ok,message)
    call require(ok,'unequal sparse density update failed: '//trim(message))
    local_degrees=0;sparse_values_ok=.true.
    do i=1,owned_count
      row=int(sparse_state%owned_row_ids(i))
      local_degrees(row)=sparse_state%operators%row_offsets(i+1)-sparse_state%operators%row_offsets(i)
      do edge=sparse_state%operators%row_offsets(i),sparse_state%operators%row_offsets(i+1)-1
        expected=sparse_state%kinetic_rows(i,sparse_state%operators%column_ids(edge))+&
          sparse_local_reference(row,sparse_state%operators%column_ids(edge))
        if(abs(sparse_state%operators%hamiltonian_values(edge)-expected)>=1d-14)sparse_values_ok=.false.
      enddo
    enddo
    call require(sparse_values_ok,'sparse density update differs from dense reference')
    call MPI_Allreduce(local_degrees,global_degrees,3,MPI_INTEGER,MPI_SUM,comm,ierr)
    call require(ierr==MPI_SUCCESS.and.all(global_degrees==[2,2,0]),&
      'unequal sparse/zero-row graph was not preserved')
    call require(size(sparse_state%local_rows)==nnz,'sparse local potential retained full columns')
  end subroutine exercise_sparse_density_update

  subroutine expect_validator_rejection(tampered,expected_fingerprint)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::tampered
    integer(int64),intent(in)::expected_fingerprint
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    call validate_rt_dg_hybrid_v3_startup(comm,tampered,expected_fingerprint,startup_tolerances,receipt,ok,message)
    call require(.not.ok,'receipt/provenance tamper survived the post-redistribution startup validator')
  end subroutine expect_validator_rejection

  subroutine write_and_expect_startup_rejection(invalid_payload,checkpoint_path)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::invalid_payload
    character(*),intent(in)::checkpoint_path
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,checkpoint_path,invalid_payload,fingerprint,ok,message)
    call require(ok,'freshly authenticated invalid v3 fixture write failed: '//trim(message))
    call initialize_rt_dg_hybrid_from_checkpoint(comm,checkpoint_path,'tddft_response',.true.,1,.false.,.false.,.false.,&
      .false.,.false.,[1,0,0],startup_tolerances,state,ok,message)
    call require(.not.ok,'freshly authenticated physically invalid v3 payload was accepted')
  end subroutine write_and_expect_startup_rejection

  subroutine write_and_expect_invariant_rejection(invalid_payload,checkpoint_path,invariant)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::invalid_payload
    character(*),intent(in)::checkpoint_path,invariant
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    type(s_rt_dg_hybrid_ground_state_payload)::authenticated_payload
    logical::detected
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,checkpoint_path,invalid_payload,fingerprint,ok,message)
    call require(ok,'freshly authenticated invariant fixture write failed: '//trim(message))
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,checkpoint_path,authenticated_payload,&
      fingerprint,ok,message)
    call require(ok,'invariant fixture authentication read failed: '//trim(message))
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,startup_tolerances,&
      receipt,ok,message)
    detected=.false.
    select case(trim(invariant))
    case('reject_orbital_physics');detected=receipt%orbital_residual>startup_tolerances(1)
    case('reject_metric_physics');detected=receipt%metric_defect>startup_tolerances(1)
    case('reject_unitarity_physics');detected=receipt%unitarity_defect>startup_tolerances(1)
    case('reject_density_physics');detected=receipt%density_defect>startup_tolerances(2)
    case('reject_electron_physics');detected=receipt%electron_defect>startup_tolerances(3)
    case('reject_target_closure');detected=receipt%target_closure_defect>startup_tolerances(4)
    case('reject_energy_covariance');detected=receipt%energy_covariance_defect>startup_tolerances(4)
    case('reject_projector_covariance');detected=receipt%projector_defect>startup_tolerances(4)
    case('reject_vector_covariance','reject_tensor_covariance')
      detected=receipt%fixed_operator_covariance_defect>startup_tolerances(4)
    end select
    call require(.not.ok.and.index(message,'startup invariant tolerance exceeded')>0.and.detected,&
      'startup invariant mutation was not physically reported: '//trim(invariant))
    call initialize_rt_dg_hybrid_from_checkpoint(comm,checkpoint_path,'tddft_response',.true.,1,&
      .false.,.false.,.false.,.false.,.false.,[1,0,0],startup_tolerances,state,ok,message)
    call require(.not.ok,'authenticated startup invariant mutation was accepted: '//trim(invariant))
  end subroutine write_and_expect_invariant_rejection

  subroutine perturb_startup_invariant(p,invariant)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    character(*),intent(in)::invariant
    integer::i
    select case(trim(invariant))
    case('reject_orbital_physics')
      p%certified_basis%certified_eigenvalues(1)=p%certified_basis%certified_eigenvalues(1)+1d-4
      p%energy_window%e_homo=p%certified_basis%certified_eigenvalues(1)
      p%energy_window%requested_cutoff=p%energy_window%e_homo+p%energy_window%window_size
      p%energy_window%extension_energy=max(0d0,p%energy_window%certified_cutoff-&
        p%energy_window%requested_cutoff)
    case('reject_metric_physics')
      do i=1,size(p%row_ids)
        if(p%row_ids(i)==1_int64)p%metric_rows(i,1)=p%metric_rows(i,1)+(1d-4,0d0)
      enddo
    case('reject_unitarity_physics')
      ! B*A=C_occ is the physical embedding contract.  Perturb only the
      ! auxiliary localization transform so startup unitarity remains an
      ! independently exercised rejecting gate.
      do i=1,size(p%certified_basis%transformation_row_ids)
        if(p%certified_basis%transformation_row_ids(i)==1_int64)&
          p%certified_basis%u_rt(i,1)=p%certified_basis%u_rt(i,1)+(1d-4,0d0)
      enddo
    case('reject_density_physics')
      do i=1,size(p%grid_ids)
        if(p%grid_ids(i)==1_int64)then
          p%density(i)=p%density(i)+1d-4;p%rt_space%density(i)=p%rt_space%density(i)+1d-4
        endif
      enddo
    case('reject_electron_physics')
      do i=1,size(p%grid_ids)
        if(p%grid_ids(i)==1_int64)p%grid_weights(i)=p%grid_weights(i)+1d-4
      enddo
    case('reject_target_closure')
      p%symmetry_representation(1,1,2)=p%symmetry_representation(1,1,2)+(1d-4,0d0)
    case('reject_energy_covariance')
      do i=1,size(p%rt_space%row_ids)
        if(p%rt_space%row_ids(i)==1_int64)then
          p%rt_space%kinetic_rows(i,1)=p%rt_space%kinetic_rows(i,1)+(1d-4,0d0)
          p%rt_space%hamiltonian_rows(i,1)=p%rt_space%hamiltonian_rows(i,1)+(1d-4,0d0)
          p%rt_space%scalar_operator_rows(i,1,1)=p%rt_space%kinetic_rows(i,1)
          p%rt_space%scalar_operator_rows(i,1,5)=p%rt_space%hamiltonian_rows(i,1)
        endif
      enddo
    case('reject_projector_covariance')
      p%certified_basis%u_rt=(0d0,0d0)
      do i=1,size(p%certified_basis%transformation_row_ids)
        p%certified_basis%u_rt(i,int(p%certified_basis%transformation_row_ids(i)))=(1d0,0d0)
      enddo
      p%certified_basis%b_rt=p%certified_basis%c_cert
      p%certified_basis%initial_occupied_amplitudes(:,1)=[(1d0,0d0),(0d0,0d0)]
    case('reject_vector_covariance')
      do i=1,size(p%rt_space%row_ids)
        if(p%rt_space%row_ids(i)==1_int64)p%rt_space%vector_operator_rows(i,1,1,1)=&
          p%rt_space%vector_operator_rows(i,1,1,1)+(1d-4,0d0)
      enddo
    case('reject_tensor_covariance')
      do i=1,size(p%rt_space%row_ids)
        if(p%rt_space%row_ids(i)==1_int64)then
          p%rt_space%tensor_operator_rows(i,1,1,2,1)=&
            p%rt_space%tensor_operator_rows(i,1,1,2,1)+(1d-4,0d0)
          p%rt_space%tensor_operator_rows(i,1,2,1,1)=&
            p%rt_space%tensor_operator_rows(i,1,2,1,1)+(1d-4,0d0)
        endif
      enddo
    case default
      call require(.false.,'unknown startup invariant mutation: '//trim(invariant))
    end select
  end subroutine perturb_startup_invariant

  subroutine exercise_density_tolerance_boundary(p,checkpoint_path)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    character(*),intent(in)::checkpoint_path
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    type(s_rt_dg_hybrid_ground_state_payload)::authenticated_payload
    real(real64)::relaxed(4),strict(4)
    integer::i
    do i=1,size(p%grid_ids)
      if(p%grid_ids(i)==1_int64)then
        p%density(i)=p%density(i)+5d-8;p%rt_space%density(i)=p%rt_space%density(i)+5d-8
      endif
    enddo
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,p,ok,message)
    call require(ok,'density-boundary named fingerprint stamping failed: '//trim(message))
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,checkpoint_path,p,fingerprint,ok,message)
    call require(ok,'density-boundary fixture write failed: '//trim(message))
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,checkpoint_path,authenticated_payload,&
      fingerprint,ok,message)
    call require(ok,'density-boundary fixture authentication read failed: '//trim(message))
    relaxed=startup_tolerances;relaxed(2)=1d-7
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,relaxed,receipt,ok,message)
    call require(ok.and.receipt%density_defect>1d-9.and.receipt%density_defect<relaxed(2),&
      'user-relaxed density tolerance did not accept the measured defect')
    strict=relaxed;strict(2)=1d-9
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,strict,receipt,ok,message)
    call require(.not.ok.and.receipt%density_defect>strict(2),&
      'user-strict density tolerance did not reject the same measured defect')
  end subroutine exercise_density_tolerance_boundary

  subroutine exercise_embedding_projector_tolerance_boundary(p,checkpoint_path)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    character(*),intent(in)::checkpoint_path
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    type(s_rt_dg_hybrid_ground_state_payload)::authenticated_payload
    real(real64)::relaxed(4),embedding_strict(4),projector_strict(4)
    p%certified_basis%initial_occupied_amplitudes(1,1)=&
      p%certified_basis%initial_occupied_amplitudes(1,1)+(5d-12,0d0)
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,p,ok,message)
    call require(ok,'embedding-boundary named fingerprint stamping failed: '//trim(message))
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,checkpoint_path,p,fingerprint,ok,message)
    call require(ok,'embedding-boundary fixture write failed: '//trim(message))
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,checkpoint_path,authenticated_payload,&
      fingerprint,ok,message)
    call require(ok,'embedding-boundary fixture authentication read failed: '//trim(message))
    relaxed=1d-11
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,relaxed,receipt,ok,message)
    call require(ok.and.receipt%embedding_defect>1d-12.and.receipt%projector_defect>1d-12,&
      'relaxed startup tolerances did not expose the sub-threshold embedding/projector defects')
    embedding_strict=relaxed;embedding_strict(1)=1d-12
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,embedding_strict,receipt,ok,message)
    call require(.not.ok.and.receipt%embedding_defect>embedding_strict(1).and.&
      receipt%projector_defect<embedding_strict(4),&
      'orbital tolerance did not independently gate the measured embedding defect')
    projector_strict=relaxed;projector_strict(4)=1d-12
    call validate_rt_dg_hybrid_v3_startup(comm,authenticated_payload,fingerprint,projector_strict,receipt,ok,message)
    call require(.not.ok.and.receipt%projector_defect>projector_strict(4).and.&
      receipt%embedding_defect<projector_strict(1),&
      'symmetry tolerance did not independently gate the measured projector defect')
  end subroutine exercise_embedding_projector_tolerance_boundary

  subroutine inject_negative_zero(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer::i
    real(real64)::negative_zero
    negative_zero=transfer(int(z'8000000000000000',int64),negative_zero)
    do i=1,size(p%rt_space%row_ids)
      if(p%rt_space%row_ids(i)==1_int64)then
        p%rt_space%vector_operator_rows(i,2,3,1)=cmplx(negative_zero,0d0,real64)
      endif
    enddo
  end subroutine inject_negative_zero

  subroutine perturb_rt_basis_nullspace(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer::i
    do i=1,size(p%grid_ids)
      if(p%grid_ids(i)==1_int64)then
        p%rt_space%basis_values(:,i)=p%rt_space%basis_values(:,i)+[cmplx(0.1d0,0d0,real64),&
          cmplx(-0.1d0,0d0,real64)]
      endif
    enddo
  end subroutine perturb_rt_basis_nullspace

  subroutine perturb_rt_fixed_components(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer::i,row
    do i=1,size(p%rt_space%row_ids)
      row=int(p%rt_space%row_ids(i))
      p%rt_space%kinetic_rows(i,row)=p%rt_space%kinetic_rows(i,row)+(0.1d0,0d0)
      p%rt_space%local_rows(i,row)=p%rt_space%local_rows(i,row)-(0.1d0,0d0)
      p%rt_space%scalar_operator_rows(i,row,1)=p%rt_space%kinetic_rows(i,row)
      p%rt_space%scalar_operator_rows(i,row,3)=p%rt_space%local_rows(i,row)
    enddo
  end subroutine perturb_rt_fixed_components

  subroutine build_v3_payload(p,requested_grid_count)
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::p
    integer,intent(in),optional::requested_grid_count
    integer::row,point,i,nrow,npoint,nrtrow,total_grid_count
    real(real64)::a
    complex(real64)::full_u(certified_rank,certified_rank),construction_h(construction_rank,construction_rank),&
      rt_h(certified_rank,certified_rank),construction_rep(construction_rank,construction_rank,operation_count)
    p%valid=.true.;p%final_refresh_complete=.true.;p%analysis_complete=.true.;p%identity_only=.false.
    total_grid_count=grid_count;if(present(requested_grid_count))total_grid_count=requested_grid_count
    p%global_count=construction_rank;p%global_grid_count=total_grid_count;p%noccupied=occupied_rank
    p%operation_count=operation_count;p%nonidentity_operation_count=1
    p%catalog_fingerprint=101_int64;p%state_fingerprint=102_int64;p%metric_fingerprint=103_int64
    p%operator_structure_fingerprint=104_int64;p%operator_value_fingerprint=105_int64
    p%basis_fingerprint=106_int64;p%face_fingerprint=107_int64;p%dc_seed_fingerprint=108_int64
    p%continuation_fingerprint=109_int64;p%analysis_fingerprint=110_int64
    p%selection_fingerprint=111_int64;p%pseudopotential_fingerprint=112_int64
    p%energy_fingerprint=113_int64
    p%position_convention_fingerprint=cell_wrapped_position_fingerprint

    nrow=count([(mod(row-1,nproc)==rank,row=1,construction_rank)])
    allocate(p%row_ids(nrow),p%metric_rows(nrow,construction_rank),p%kinetic_rows(nrow,construction_rank),&
      p%nonlocal_rows(nrow,construction_rank),p%local_rows(nrow,construction_rank),&
      p%sipg_rows(nrow,construction_rank),p%hamiltonian_rows(nrow,construction_rank),&
      p%coefficients(nrow,occupied_rank),p%position_rows(3,nrow,construction_rank),&
      p%metric_row_offsets(nrow+1),p%metric_column_ids(nrow*construction_rank),&
      p%operator_row_offsets(nrow+1),p%operator_column_ids(nrow*construction_rank))
    construction_h=(0d0,0d0)
    construction_h(1,1)=(-0.5d0,0d0);construction_h(2,2)=(0.25d0,0d0);construction_h(3,3)=(1d0,0d0)
    p%metric_rows=(0d0,0d0);p%kinetic_rows=(0d0,0d0);p%nonlocal_rows=(0d0,0d0)
    p%local_rows=(0d0,0d0);p%sipg_rows=(0d0,0d0);p%hamiltonian_rows=(0d0,0d0)
    p%coefficients=(0d0,0d0);p%position_rows=(0d0,0d0)
    p%metric_row_offsets(1)=1;p%operator_row_offsets(1)=1;i=0
    do row=1,construction_rank
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;p%row_ids(i)=row;p%metric_rows(i,row)=(1d0,0d0)
      p%hamiltonian_rows(i,:)=construction_h(row,:)
      p%kinetic_rows(i,:)=0.5d0*construction_h(row,:);p%nonlocal_rows(i,:)=0.25d0*construction_h(row,:)
      p%local_rows(i,:)=0.125d0*construction_h(row,:);p%sipg_rows(i,:)=0.125d0*construction_h(row,:)
      if(row==1)p%coefficients(i,1)=(1d0,0d0)
      p%position_rows(1,i,row)=(0.25d0,0d0)
      if(row==1)p%position_rows(1,i,2)=(1d0,0d0)
      if(row==2)p%position_rows(1,i,1)=(1d0,0d0)
      p%metric_column_ids((i-1)*construction_rank+1:i*construction_rank)=[(point,point=1,construction_rank)]
      p%operator_column_ids((i-1)*construction_rank+1:i*construction_rank)=[(point,point=1,construction_rank)]
      p%metric_row_offsets(i+1)=i*construction_rank+1;p%operator_row_offsets(i+1)=i*construction_rank+1
    enddo
    construction_rep=(0d0,0d0)
    do row=1,construction_rank
      construction_rep(row,row,1)=(1d0,0d0)
      construction_rep(row,row,2)=merge((1d0,0d0),(-1d0,0d0),row/=2)
    enddo
    allocate(p%symmetry_representation,source=construction_rep)

    npoint=count([(mod(point-1,nproc)==rank,point=1,total_grid_count)])
    allocate(p%grid_ids(npoint),p%grid_weights(npoint),p%partition_ids(npoint),&
      p%basis_values(construction_rank,npoint),p%density(npoint))
    p%grid_weights=1d0;p%partition_ids=1;p%basis_values=(0d0,0d0);i=0
    do point=1,total_grid_count
      if(mod(point-1,nproc)/=rank)cycle
      i=i+1;p%grid_ids(i)=point;p%density(i)=merge(1d0,0d0,point==1)
      p%basis_values(1,i)=merge((1d0,0d0),(0d0,0d0),point==1)
      p%basis_values(2,i)=merge((0d0,0d0),(1d0,0d0),point==1)
      p%basis_values(3,i)=cmplx(0.1d0*point,0.05d0,real64)
    enddo
    call build_common_metadata(p)

    a=1d0/sqrt(2d0);full_u=reshape([cmplx(a,0d0,real64),cmplx(a,0d0,real64),&
      cmplx(a,0d0,real64),cmplx(-a,0d0,real64)],[certified_rank,certified_rank])
    nrtrow=count([(mod(row-1,nproc)==rank,row=1,certified_rank)])
    call build_certified_basis(p,full_u,nrow,nrtrow)
    rt_h=matmul(conjg(transpose(full_u)),matmul(construction_h(1:certified_rank,1:certified_rank),full_u))
    call build_rt_space(p,full_u,rt_h,npoint,nrtrow)
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%kinetic_rows,p%kinetic_fingerprint,ok)
    call require(ok,'construction kinetic fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%nonlocal_rows,p%nonlocal_fingerprint,ok)
    call require(ok,'construction nonlocal fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%local_rows,p%local_fingerprint,ok)
    call require(ok,'construction local fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%sipg_rows,p%sipg_fingerprint,ok)
    call require(ok,'construction SIPG fingerprint failed')
    call build_named_receipts(p)
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,p,ok,message)
    call require(ok,'certified v3 named fingerprint stamping failed: '//trim(message))
  end subroutine build_v3_payload

  subroutine make_production_smoke_payload(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    p%kinetic_rows=(0d0,0d0);p%nonlocal_rows=(0d0,0d0);p%local_rows=(0d0,0d0)
    p%sipg_rows=(0d0,0d0);p%hamiltonian_rows=(0d0,0d0);p%basis_values=(0d0,0d0);p%density=0d0
    p%occupations=0d0;p%eigenvalues=0d0
    p%certified_basis%certified_eigenvalues=0d0;p%certified_basis%occupations=0d0
    p%rt_space%kinetic_rows=(0d0,0d0);p%rt_space%nonlocal_rows=(0d0,0d0)
    p%rt_space%local_rows=(0d0,0d0);p%rt_space%sipg_rows=(0d0,0d0)
    p%rt_space%hamiltonian_rows=(0d0,0d0);p%rt_space%basis_values=(0d0,0d0);p%rt_space%density=0d0
    p%rt_space%scalar_operator_rows=(0d0,0d0)
    p%electron_count%expected_count=0d0;p%electron_count%actual_count=0d0
    p%electron_count%defect=0d0;p%electron_count%omitted_tail=0d0
    p%energy_window%e_homo=0d0;p%energy_window%requested_cutoff=p%energy_window%window_size
    p%energy_window%certified_cutoff=0d0;p%energy_window%extension_energy=0d0
    p%energy_receipt=0d0
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%kinetic_rows,p%kinetic_fingerprint,ok)
    call require(ok,'production construction kinetic fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%nonlocal_rows,p%nonlocal_fingerprint,ok)
    call require(ok,'production construction nonlocal fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%local_rows,p%local_fingerprint,ok)
    call require(ok,'production construction local fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%sipg_rows,p%sipg_fingerprint,ok)
    call require(ok,'production construction SIPG fingerprint failed')
    call stamp_rt_dg_hybrid_v3_fingerprints(comm,p,ok,message)
    call require(ok,'production v3 named fingerprint stamping failed: '//trim(message))
  end subroutine make_production_smoke_payload

  subroutine build_common_metadata(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    allocate(p%requested_ids(2),p%effective_ids(3),p%added_ids(1),p%closure_parent(1),&
      p%closure_reason(1),p%closure_action(1),p%scope_selectors(8),p%xc_types(3))
    p%requested_ids=[11,12];p%effective_ids=[11,12,13];p%added_ids=[13]
    p%closure_parent=[11];p%closure_reason=[1];p%closure_action=[2]
    p%scope_selectors=[1,1,1,0,0,0,0,0];p%xc_types=[1,0,0]
    p%scope_fingerprint=fingerprint_rt_dg_hybrid_scope(p%scope_selectors,p%xc_types)
    allocate(p%occupations(1),p%eigenvalues(1),p%continuation_receipt(1),&
      p%pseudopotential_receipt(1),p%energy_receipt(7))
    p%occupations=[1d0];p%eigenvalues=[-0.5d0];p%continuation_receipt=0d0
    p%pseudopotential_receipt=0d0;p%energy_receipt=0d0
    allocate(p%face_ids(0),p%face_point_ids(0),p%face_metadata(8,0),p%face_offsets(1),&
      p%face_weight_offsets(1),p%face_basis_offsets(1),p%face_value_offsets(1),p%face_observable_offsets(1),&
      p%face_basis_ids(0),p%face_normals(3,0),p%face_weights(0),p%face_values(1,0),&
      p%interface_observables(3,0),p%nonlocal_ids(0),p%nonlocal_owner(0),p%nonlocal_values(1,0))
    p%face_offsets=1;p%face_weight_offsets=1;p%face_basis_offsets=1
    p%face_value_offsets=1;p%face_observable_offsets=1
    p%construction_catalog%valid=.true.;p%construction_catalog%global_count=construction_rank
    allocate(p%construction_catalog%ids(3),p%construction_catalog%generations(3),&
      p%construction_catalog%ordering(3),p%construction_catalog%ownership(3))
    p%construction_catalog%ids=[11_int64,12_int64,13_int64]
    p%construction_catalog%generations=[1,1,2];p%construction_catalog%ordering=[2,3,1]
    p%construction_catalog%ownership=[1,1,1]
    p%construction_catalog%ids_fingerprint=201_int64;p%construction_catalog%generation_fingerprint=202_int64
    p%construction_catalog%ordering_fingerprint=203_int64;p%construction_catalog%ownership_fingerprint=204_int64
    p%construction_catalog%provenance_fingerprint=2105_int64
    p%construction_catalog%catalog_fingerprint=p%catalog_fingerprint
  end subroutine build_common_metadata

  subroutine build_certified_basis(p,u,nrow,nrtrow)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    complex(real64),intent(in)::u(certified_rank,certified_rank)
    integer,intent(in)::nrow,nrtrow
    integer::row,i
    p%certified_basis%valid=.true.;p%certified_basis%localization_converged=.true.
    p%certified_basis%localization_symmetry_constrained=.false.
    p%certified_basis%construction_count=construction_rank;p%certified_basis%certified_count=certified_rank
    p%certified_basis%occupied_count=occupied_rank;p%certified_basis%localization_iterations=3
    allocate(p%certified_basis%construction_row_ids(nrow),p%certified_basis%transformation_row_ids(nrtrow),&
      p%certified_basis%c_cert(nrow,certified_rank),p%certified_basis%u_rt(nrtrow,certified_rank),&
      p%certified_basis%b_rt(nrow,certified_rank),&
      p%certified_basis%initial_occupied_amplitudes(certified_rank,occupied_rank),&
      p%certified_basis%certified_eigenvalues(certified_rank),p%certified_basis%occupations(occupied_rank),&
      p%certified_basis%centers(3,certified_rank),p%certified_basis%spreads_before(certified_rank),&
      p%certified_basis%spreads_after(certified_rank))
    p%certified_basis%c_cert=(0d0,0d0);p%certified_basis%b_rt=(0d0,0d0);i=0
    do row=1,construction_rank
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;p%certified_basis%construction_row_ids(i)=row
      if(row<=certified_rank)p%certified_basis%c_cert(i,row)=(1d0,0d0)
      p%certified_basis%b_rt(i,:)=matmul(p%certified_basis%c_cert(i,:),u)
    enddo
    i=0
    do row=1,certified_rank
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;p%certified_basis%transformation_row_ids(i)=row;p%certified_basis%u_rt(i,:)=u(row,:)
    enddo
    p%certified_basis%initial_occupied_amplitudes=conjg(transpose(u(1:occupied_rank,:)))
    p%certified_basis%certified_eigenvalues=[-0.5d0,0.25d0];p%certified_basis%occupations=[1d0]
    p%certified_basis%centers=0d0;p%certified_basis%spreads_before=[1d0,1d0]
    p%certified_basis%spreads_after=[0.8d0,0.8d0];p%certified_basis%spread_before_total=2d0
    p%certified_basis%spread_after_total=1.6d0;p%certified_basis%spread_improvement=0.4d0
    p%certified_basis%transform_unitarity_defect=0d0;p%certified_basis%certified_metric_defect=0d0
    p%certified_basis%rt_metric_defect=0d0;p%certified_basis%embedding_defect=0d0
    p%certified_basis%projector_invariance_defect=0d0
    p%certified_basis%target_symmetry_defect_before=0d0;p%certified_basis%target_symmetry_defect_after=0d0
    p%certified_basis%energy_symmetry_defect_before=0d0;p%certified_basis%energy_symmetry_defect_after=0d0
    p%certified_basis%symmetry_defect_invariance=0d0;p%certified_basis%scalar_covariance_defect=0d0
    p%certified_basis%vector_covariance_defect=0d0;p%certified_basis%tensor_covariance_defect=0d0
    p%certified_basis%c_cert_fingerprint=301_int64;p%certified_basis%u_rt_fingerprint=302_int64
    p%certified_basis%b_rt_fingerprint=303_int64;p%certified_basis%initial_state_fingerprint=304_int64
    p%certified_basis%transformation_fingerprint=305_int64
    p%certified_basis%operator_fingerprint=306_int64;p%certified_basis%fingerprint=307_int64
  end subroutine build_certified_basis

  subroutine build_rt_space(p,u,h,npoint,nrtrow)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    complex(real64),intent(in)::u(certified_rank,certified_rank),h(certified_rank,certified_rank)
    integer,intent(in)::npoint,nrtrow
    integer::row,i,point
    real(real64)::a
    a=1d0/sqrt(2d0);p%rt_space%valid=.true.;p%rt_space%rank=certified_rank
    p%rt_space%operation_count=operation_count;p%rt_space%scalar_count=5
    p%rt_space%vector_count=1;p%rt_space%tensor_count=1
    allocate(p%rt_space%row_ids(nrtrow),p%rt_space%row_owner_keys(certified_rank),&
      p%rt_space%grid_owner_keys(npoint),p%rt_space%metric_rows(nrtrow,certified_rank),&
      p%rt_space%kinetic_rows(nrtrow,certified_rank),p%rt_space%nonlocal_rows(nrtrow,certified_rank),&
      p%rt_space%local_rows(nrtrow,certified_rank),p%rt_space%sipg_rows(nrtrow,certified_rank),&
      p%rt_space%hamiltonian_rows(nrtrow,certified_rank),&
      p%rt_space%representation(certified_rank,certified_rank,operation_count),&
      p%rt_space%cartesian_rotations(3,3,operation_count),p%rt_space%scalar_operator_rows(nrtrow,certified_rank,5),&
      p%rt_space%vector_operator_rows(nrtrow,certified_rank,3,1),&
      p%rt_space%tensor_operator_rows(nrtrow,certified_rank,3,3,1),&
      p%rt_space%basis_values(certified_rank,npoint),p%rt_space%density(npoint))
    do row=1,certified_rank;p%rt_space%row_owner_keys(row)=mod(row-1,nproc)+1;enddo
    p%rt_space%grid_owner_keys=rank+1;p%rt_space%metric_rows=(0d0,0d0)
    p%rt_space%kinetic_rows=(0d0,0d0);p%rt_space%nonlocal_rows=(0d0,0d0)
    p%rt_space%local_rows=(0d0,0d0);p%rt_space%sipg_rows=(0d0,0d0);p%rt_space%hamiltonian_rows=(0d0,0d0)
    p%rt_space%scalar_operator_rows=(0d0,0d0);p%rt_space%vector_operator_rows=(0d0,0d0)
    p%rt_space%tensor_operator_rows=(0d0,0d0);i=0
    do row=1,certified_rank
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;p%rt_space%row_ids(i)=row;p%rt_space%metric_rows(i,row)=(1d0,0d0)
      p%rt_space%hamiltonian_rows(i,:)=h(row,:);p%rt_space%kinetic_rows(i,:)=0.5d0*h(row,:)
      p%rt_space%nonlocal_rows(i,:)=0.25d0*h(row,:);p%rt_space%local_rows(i,:)=0.125d0*h(row,:)
      p%rt_space%sipg_rows(i,:)=0.125d0*h(row,:)
      p%rt_space%scalar_operator_rows(i,:,1)=p%rt_space%kinetic_rows(i,:)
      p%rt_space%scalar_operator_rows(i,:,2)=p%rt_space%nonlocal_rows(i,:)
      p%rt_space%scalar_operator_rows(i,:,3)=p%rt_space%local_rows(i,:)
      p%rt_space%scalar_operator_rows(i,:,4)=p%rt_space%sipg_rows(i,:)
      p%rt_space%scalar_operator_rows(i,:,5)=p%rt_space%hamiltonian_rows(i,:)
      p%rt_space%vector_operator_rows(i,row,2,1)=(1d0,0d0)
      p%rt_space%tensor_operator_rows(i,row,1,1,1)=(1d0,0d0)
      p%rt_space%tensor_operator_rows(i,row,1,2,1)=merge((1d0,0d0),(-1d0,0d0),row==1)
      p%rt_space%tensor_operator_rows(i,row,2,1,1)=p%rt_space%tensor_operator_rows(i,row,1,2,1)
      if(row==1)p%rt_space%vector_operator_rows(i,row,1,1)=(1d0,0d0)
      if(row==2)p%rt_space%vector_operator_rows(i,row,1,1)=(-1d0,0d0)
    enddo
    p%rt_space%representation=(0d0,0d0)
    p%rt_space%representation(1,1,1)=(1d0,0d0);p%rt_space%representation(2,2,1)=(1d0,0d0)
    p%rt_space%representation(1,2,2)=(1d0,0d0);p%rt_space%representation(2,1,2)=(1d0,0d0)
    p%rt_space%cartesian_rotations=0d0
    do row=1,3
      p%rt_space%cartesian_rotations(row,row,1)=1d0
      p%rt_space%cartesian_rotations(row,row,2)=merge(-1d0,1d0,row==1)
    enddo
    do i=1,npoint
      point=int(p%grid_ids(i));p%rt_space%density(i)=merge(1d0,0d0,point==1)
      p%rt_space%basis_values(:,i)=merge([cmplx(a,0d0,real64),cmplx(a,0d0,real64)],&
        [cmplx(a,0d0,real64),cmplx(-a,0d0,real64)],point==1)
    enddo
    p%rt_space%metric_fingerprint=501_int64;p%rt_space%kinetic_fingerprint=502_int64
    p%rt_space%nonlocal_fingerprint=503_int64;p%rt_space%local_fingerprint=504_int64
    p%rt_space%sipg_fingerprint=505_int64;p%rt_space%hamiltonian_fingerprint=506_int64
    p%rt_space%basis_fingerprint=507_int64;p%rt_space%density_fingerprint=508_int64
    p%rt_space%ownership_fingerprint=509_int64;p%rt_space%scalar_fingerprint=510_int64
    p%rt_space%vector_fingerprint=511_int64;p%rt_space%tensor_fingerprint=512_int64
    p%rt_space%representation_fingerprint=513_int64;p%rt_space%fingerprint=514_int64
  end subroutine build_rt_space

  subroutine build_named_receipts(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    p%electron_count%valid=.true.;p%electron_count%expected_count=1d0;p%electron_count%actual_count=1d0
    p%electron_count%tolerance=startup_tolerances(3);p%electron_count%defect=0d0
    p%electron_count%omitted_tail=0d0;p%electron_count%chemical_potential=-0.1d0
    p%electron_count%fingerprint=401_int64
    p%energy_window%valid=.true.;p%energy_window%compatibility_dynamic_rank=.false.
    p%energy_window%proof_state_present=.true.;p%energy_window%mode=rt_dg_hybrid_energy_window_explicit
    p%energy_window%construction_rank=construction_rank;p%energy_window%solved_rank=construction_rank
    p%energy_window%occupied_rank=occupied_rank;p%energy_window%requested_rank=occupied_rank
    p%energy_window%certified_rank=certified_rank;p%energy_window%extension_states=1
    p%energy_window%boundary_cluster_rank=certified_rank;p%energy_window%proof_status=1
    p%energy_window%window_size=0.4d0;p%energy_window%e_homo=-0.5d0
    p%energy_window%requested_cutoff=p%energy_window%e_homo+p%energy_window%window_size
    p%energy_window%certified_cutoff=p%certified_basis%certified_eigenvalues(certified_rank)
    p%energy_window%extension_energy=max(0d0,p%energy_window%certified_cutoff-&
      p%energy_window%requested_cutoff);p%energy_window%proof_energy=1d0
    p%energy_window%fingerprint=601_int64
    p%symmetry_receipt%valid=.true.;p%symmetry_receipt%worst_operation=2
    p%symmetry_receipt%occupied_subspace_defect=0d0;p%symmetry_receipt%occupied_projector_defect=0d0
    p%symmetry_receipt%target_subspace_defect=0d0;p%symmetry_receipt%target_energy_defect=0d0
    p%symmetry_receipt%density_defect=0d0;p%symmetry_receipt%scalar_covariance_defect=0d0
    p%symmetry_receipt%vector_covariance_defect=0d0;p%symmetry_receipt%tensor_covariance_defect=0d0
    p%symmetry_receipt%final_basis_defect=0d0;p%symmetry_receipt%worst_operation_defect=0d0
    p%symmetry_receipt%maximum_physical_defect=0d0;p%symmetry_receipt%fingerprint=2601_int64
    p%handoff_receipts%valid=.true.;p%handoff_receipts%position_fingerprint=p%position_convention_fingerprint
    p%handoff_receipts%nonlocal_fingerprint=p%nonlocal_fingerprint
    p%handoff_receipts%face_fingerprint=p%face_fingerprint
    p%handoff_receipts%pseudopotential_fingerprint=p%pseudopotential_fingerprint
    p%handoff_receipts%transformation_fingerprint=p%certified_basis%transformation_fingerprint
    p%handoff_receipts%fingerprint=701_int64
  end subroutine build_named_receipts

  subroutine project_density_local(row_ids,row_offsets,column_ids,grid_ids,density,local_values,callback_ok,callback_message)
    integer(int64),intent(in)::row_ids(:),grid_ids(:)
    integer,intent(in)::row_offsets(:),column_ids(:)
    real(real64),intent(in)::density(:)
    complex(real64),intent(out)::local_values(:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    integer::q,edge,row_position,reduction_ierr
    real(real64)::local_sum,global_sum
    local_sum=sum(density)
    call MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,reduction_ierr)
    local_values=(0d0,0d0);callback_ok=reduction_ierr==MPI_SUCCESS;callback_message=''
    if(.not.callback_ok)then;callback_message='fixture density reduction failed';return;endif
    if(force_callback_failure.and.rank==nproc-1)then
      callback_ok=.false.;callback_message='intentional one-rank callback failure';return
    endif
    do q=1,size(row_ids)
      row_position=findloc(payload%rt_space%row_ids,row_ids(q),dim=1)
      if(row_position<1)then;callback_ok=.false.;callback_message='fixture RT row ownership mismatch';return;endif
      do edge=row_offsets(q),row_offsets(q+1)-1
        local_values(edge)=payload%rt_space%local_rows(row_position,column_ids(edge))
        if(column_ids(edge)==int(row_ids(q)))local_values(edge)=local_values(edge)+&
          cmplx(global_sum-reference_density_total+reference_shift,0d0,real64)
      enddo
    enddo
  end subroutine project_density_local

  subroutine project_sparse_density_local(row_ids,row_offsets,column_ids,grid_ids,density,local_values,&
      callback_ok,callback_message)
    integer(int64),intent(in)::row_ids(:),grid_ids(:)
    integer,intent(in)::row_offsets(:),column_ids(:)
    real(real64),intent(in)::density(:)
    complex(real64),intent(out)::local_values(:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    integer::i,edge,row
    local_values=(0d0,0d0)
    do i=1,size(row_ids)
      row=int(row_ids(i))
      do edge=row_offsets(i),row_offsets(i+1)-1
        local_values(edge)=density(1)*sparse_local_reference(row,column_ids(edge))
      enddo
    enddo
    callback_ok=size(grid_ids)==size(density);callback_message=''
  end subroutine project_sparse_density_local

  logical function any_rank(values)
    logical,intent(in)::values(:)
    integer::local_value,global_value,reduction_ierr
    local_value=merge(1,0,any(values))
    call MPI_Allreduce(local_value,global_value,1,MPI_INTEGER,MPI_MAX,comm,reduction_ierr)
    any_rank=reduction_ierr==MPI_SUCCESS.and.global_value==1
  end function any_rank

  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_rt_dg_hybrid_initialization_mpi
