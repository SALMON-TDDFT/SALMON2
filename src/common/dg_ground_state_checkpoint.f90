!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_ground_state_checkpoint
#ifdef USE_MPI
  use mpi
#endif
  use,intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use,intrinsic :: iso_fortran_env, only: int8,int64,iostat_end
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_state, only: nodal_face_count,dg_nodal_face_slot
  use dg_dc_local_basis_ground_state, only: s_dg_dc_local_basis_production_state
  use dg_dc_ground_state, only: s_dg_dc_gs_result,s_dg_dc_gs_diagnostics,s_dg_dc_gs_controls,DG_DC_ACCEPTED
  implicit none
  private
  character(20),parameter :: payload_magic='DG_DC_GS_CHECKPOINT1'
  character(18),parameter :: manifest_magic='DG_DC_GS_MANIFEST1'
  integer,parameter,public :: dg_ground_state_checkpoint_version=1
  integer,parameter :: maximum_dimension=100000000

  type,public :: s_dg_dc_direct_checkpoint_state
    logical :: ready=.false.,accepted=.false.,unmixed_verified=.false.
    integer :: fragment_id=0,nstate=0,state_start=0,state_count=0,nspin=0,scf_iterations=0
    integer :: core_size(3)=0,full_spatial_shape(3)=0
    integer :: orbital_spatial_lower_bound(3)=0
    integer :: orbital_core_local_origin(3)=0,orbital_core_origin(3)=0,orbital_core_size(3)=0
    integer :: global_size(3)=0,fragment_origin(3)=0,fragment_size(3)=0,face_neighbors(6)=0
    integer :: density_origin(3)=0,density_size(3)=0
    integer :: continuation_rollbacks=0,projector_generation=0
    integer :: projector_retained_rank(6)=0,projector_required_rank(6)=0
    integer(int64) :: geometry_fingerprint=0_int64,operator_fingerprint=0_int64
    integer(int64) :: hamiltonian_operator_fingerprint=0_int64
    integer(int64) :: projector_fingerprint=0_int64
    real(8) :: density_residual=huge(1d0),interface_scale=0d0
    real(8) :: orbital_residual=huge(1d0),orthogonality_defect=huge(1d0)
    real(8) :: hermiticity_defect=huge(1d0),charge_error=huge(1d0),face_balance_defect=huge(1d0)
    real(8) :: density_tolerance=0d0,orbital_tolerance=0d0,orthogonality_tolerance=0d0
    real(8) :: hermiticity_tolerance=0d0,charge_tolerance=0d0,face_balance_tolerance=0d0
    real(8) :: projector_projection_residual(6)=huge(1d0),projector_escape_norm(6)=huge(1d0)
    real(8) :: projector_residual_limit(6)=0d0,projector_escape_limit(6)=0d0
    real(8) :: face_action_norm(6)=huge(1d0),face_pair_balance(3)=huge(1d0)
    real(8),allocatable :: fragment_orbitals(:,:,:,:,:,:,:)
    real(8),allocatable :: occupations(:,:,:)
    real(8),allocatable :: continuation_scale_history(:),continuation_step_history(:)
    real(8),allocatable :: density(:,:,:,:),hartree(:,:,:),vxc(:,:,:,:),vlocal(:,:,:,:)
  end type s_dg_dc_direct_checkpoint_state

  type,public :: s_dg_ground_state_checkpoint
    integer :: schema_version=dg_ground_state_checkpoint_version
    integer(int64) :: publication_generation=0_int64
    character(16) :: provenance=''
    integer :: global_size(3)=0
    integer,allocatable :: fragment_origins(:,:),fragment_sizes(:,:),face_neighbors(:)
    type(s_dg_dc_local_basis_production_state) :: state
    type(s_dg_dc_gs_result) :: result
    type(s_dg_dc_gs_controls) :: controls
    real(8),allocatable :: density(:,:,:,:),hartree(:,:,:),vxc(:,:,:,:),vlocal(:,:,:,:)
  end type

  public :: populate_dg_ground_state_checkpoint,validate_dg_ground_state_checkpoint
  public :: publish_dg_ground_state_checkpoint,read_dg_ground_state_checkpoint
  public :: publish_dg_dc_direct_checkpoint,read_dg_dc_direct_checkpoint

  interface
    integer(c_int) function c_rename(old_path,new_path) bind(C,name='rename')
      import c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
    end function
  end interface
contains
  subroutine populate_dg_ground_state_checkpoint(checkpoint,state,result,controls,density,hartree,vxc,vlocal, &
      global_size,fragment_origins,fragment_sizes,face_neighbors,provenance,ok,message)
    type(s_dg_ground_state_checkpoint),intent(out)::checkpoint
    type(s_dg_dc_local_basis_production_state),intent(in)::state
    type(s_dg_dc_gs_result),intent(in)::result
    type(s_dg_dc_gs_controls),intent(in)::controls
    real(8),intent(in)::density(:,:,:,:),hartree(:,:,:),vxc(:,:,:,:),vlocal(:,:,:,:)
    integer,intent(in)::global_size(3),fragment_origins(:,:),fragment_sizes(:,:),face_neighbors(:)
    character(*),intent(in)::provenance
    logical,intent(out)::ok
    character(*),intent(out)::message
    checkpoint%schema_version=dg_ground_state_checkpoint_version
    checkpoint%provenance=provenance
    checkpoint%global_size=global_size
    checkpoint%fragment_origins=fragment_origins
    checkpoint%fragment_sizes=fragment_sizes
    checkpoint%face_neighbors=face_neighbors
    checkpoint%state=state
    checkpoint%result=result
    checkpoint%controls=controls
    checkpoint%density=density
    checkpoint%hartree=hartree
    checkpoint%vxc=vxc
    checkpoint%vlocal=vlocal
    call validate_dg_ground_state_checkpoint(checkpoint,state%geometry_fingerprint, &
      state%operator_fingerprint,'DG_DC_GS',ok,message)
  end subroutine

  subroutine validate_dg_ground_state_checkpoint(checkpoint,expected_geometry,expected_operator, &
      expected_provenance,ok,message)
    type(s_dg_ground_state_checkpoint),intent(in)::checkpoint
    integer(int64),intent(in)::expected_geometry,expected_operator
    character(*),intent(in)::expected_provenance
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8)::diagnostic_values(10)
    ok=.false.
    if(checkpoint%schema_version/=dg_ground_state_checkpoint_version .or. &
       trim(checkpoint%provenance)/='DG_DC_GS' .or. &
       trim(expected_provenance)/='DG_DC_GS')then
      message='DG ground-state checkpoint schema or provenance mismatch';return
    end if
    if(checkpoint%state%geometry_fingerprint/=expected_geometry .or. &
       checkpoint%state%operator_fingerprint/=expected_operator)then
      message='DG ground-state checkpoint fingerprint mismatch';return
    end if
    if(.not.checkpoint%state%ready .or. checkpoint%state%fragment_id<=0 .or. &
       checkpoint%state%local_basis_count<=0 .or. checkpoint%state%global_basis_count<=0 .or. &
       checkpoint%state%global_band_count<=0 .or. checkpoint%state%raw_basis_count<=0 .or. &
       checkpoint%state%hamiltonian_operator_fingerprint==0_int64 .or. &
       any(checkpoint%state%core_size<=0) .or. &
       any(checkpoint%state%full_spatial_shape<checkpoint%state%core_size) .or. &
       .not.allocated(checkpoint%state%coefficient_rows) .or. &
       .not.allocated(checkpoint%state%full_fragment_basis) .or. &
       .not.allocated(checkpoint%state%basis_transform) .or. &
       .not.allocated(checkpoint%state%eigenvalues) .or. &
       .not.allocated(checkpoint%state%occupations) .or. &
       .not.allocated(checkpoint%state%basis_offsets) .or. &
       .not.allocated(checkpoint%state%fragment_ids))then
      message='DG ground-state checkpoint local-basis state is incomplete';return
    end if
    if(.not.allocated(checkpoint%fragment_origins) .or. .not.allocated(checkpoint%fragment_sizes) .or. &
       .not.allocated(checkpoint%face_neighbors) .or. any(checkpoint%global_size<=0))then
      message='DG ground-state checkpoint physical topology is incomplete';return
    end if
    if(any(shape(checkpoint%fragment_origins)/=shape(checkpoint%fragment_sizes)) .or. &
       size(checkpoint%fragment_origins,1)/=3 .or. size(checkpoint%face_neighbors)/=nodal_face_count .or. &
       checkpoint%state%fragment_id<1 .or. checkpoint%state%fragment_id>size(checkpoint%fragment_origins,2))then
      message='DG ground-state checkpoint physical topology dimensions are invalid';return
    end if
    if(any(checkpoint%fragment_origins<0) .or. any(checkpoint%fragment_sizes<=0) .or. &
       any(checkpoint%fragment_origins+checkpoint%fragment_sizes> &
       spread(checkpoint%global_size,2,size(checkpoint%fragment_origins,2))) .or. &
       any(checkpoint%fragment_sizes(:,checkpoint%state%fragment_id)/=checkpoint%state%core_size) .or. &
       any(checkpoint%face_neighbors<0) .or. &
       any(checkpoint%face_neighbors>=size(checkpoint%fragment_origins,2)) .or. &
       .not.face_neighbors_match_geometry(checkpoint))then
      message='DG ground-state checkpoint physical face ownership is invalid';return
    end if
    if(any(shape(checkpoint%state%coefficient_rows)/= &
       [checkpoint%state%local_basis_count,checkpoint%state%global_band_count]) .or. &
       int(size(checkpoint%state%full_fragment_basis,1),int64)/= &
       safe_product3(checkpoint%state%full_spatial_shape) .or. &
       size(checkpoint%state%full_fragment_basis,2)/=checkpoint%state%local_basis_count .or. &
       size(checkpoint%state%basis_transform,1)/=checkpoint%state%raw_basis_count .or. &
       size(checkpoint%state%basis_transform,2)/=checkpoint%state%local_basis_count .or. &
       size(checkpoint%state%eigenvalues)/=checkpoint%state%global_band_count .or. &
       size(checkpoint%state%occupations)/=checkpoint%state%global_band_count .or. &
       size(checkpoint%state%basis_offsets)/=size(checkpoint%state%fragment_ids)+1 .or. &
       lbound(checkpoint%state%basis_offsets,1)/=0 .or. &
       checkpoint%state%basis_offsets(lbound(checkpoint%state%basis_offsets,1))/=0 .or. &
       checkpoint%state%basis_offsets(ubound(checkpoint%state%basis_offsets,1))/= &
       checkpoint%state%global_basis_count .or. &
       any(checkpoint%state%basis_offsets(1:)<=checkpoint%state%basis_offsets(: &
       ubound(checkpoint%state%basis_offsets,1)-1)) .or. &
       any(checkpoint%state%fragment_ids<1) .or. &
       any(checkpoint%state%fragment_ids>size(checkpoint%fragment_origins,2)) .or. &
       .not.valid_fragment_basis_layout(checkpoint%state) .or. &
       .not.all(ieee_is_finite(real(checkpoint%state%coefficient_rows,8))) .or. &
       .not.all(ieee_is_finite(aimag(checkpoint%state%coefficient_rows))) .or. &
       .not.all(ieee_is_finite(real(checkpoint%state%full_fragment_basis,8))) .or. &
       .not.all(ieee_is_finite(aimag(checkpoint%state%full_fragment_basis))) .or. &
       .not.all(ieee_is_finite(real(checkpoint%state%basis_transform,8))) .or. &
       .not.all(ieee_is_finite(aimag(checkpoint%state%basis_transform))) .or. &
       .not.all(ieee_is_finite(checkpoint%state%eigenvalues)) .or. &
       .not.all(ieee_is_finite(checkpoint%state%occupations)) .or. &
       any(checkpoint%state%occupations<0d0))then
      message='DG ground-state checkpoint local-basis payload is invalid';return
    end if
    if(.not.allocated(checkpoint%result%lambda_history) .or. &
       .not.allocated(checkpoint%result%lambda_steps) .or. checkpoint%result%naccepted_steps<=0 .or. &
       checkpoint%result%naccepted_steps>size(checkpoint%result%lambda_history) .or. &
       checkpoint%result%naccepted_steps>size(checkpoint%result%lambda_steps) .or. &
       checkpoint%result%phase/=DG_DC_ACCEPTED .or. checkpoint%result%lambda/=1d0 .or. &
       .not.checkpoint%result%accepted .or. checkpoint%result%failed .or. &
       .not.checkpoint%result%unmixed_verified)then
      message='DG ground-state checkpoint continuation was not accepted';return
    end if
    if(checkpoint%result%lambda_history(1)/=0d0 .or. checkpoint%result%lambda_steps(1)/=0d0 .or. &
       checkpoint%result%lambda_history(checkpoint%result%naccepted_steps)/=1d0 .or. &
       checkpoint%result%lambda/=checkpoint%result%lambda_history(checkpoint%result%naccepted_steps) .or. &
       checkpoint%result%maximum_interface_scale/=maxval(checkpoint%result%lambda_history( &
       :checkpoint%result%naccepted_steps)) .or. &
       any(checkpoint%result%lambda_history(2:checkpoint%result%naccepted_steps)< &
       checkpoint%result%lambda_history(1:checkpoint%result%naccepted_steps-1)) .or. &
       any(checkpoint%result%lambda_steps(2:checkpoint%result%naccepted_steps)/= &
       checkpoint%result%lambda_history(2:checkpoint%result%naccepted_steps)- &
       checkpoint%result%lambda_history(1:checkpoint%result%naccepted_steps-1)) .or. &
       any(checkpoint%result%lambda_history(:checkpoint%result%naccepted_steps)<0d0) .or. &
       any(checkpoint%result%lambda_history(:checkpoint%result%naccepted_steps)>1d0) .or. &
       any(checkpoint%result%lambda_steps(:checkpoint%result%naccepted_steps)<0d0))then
      message='DG ground-state checkpoint continuation history is invalid';return
    end if
    diagnostic_values=diagnostic_vector(checkpoint%result%final_diagnostics)
    if(.not.valid_checkpoint_controls(checkpoint%controls) .or. &
       .not.checkpoint%result%final_diagnostics%finite .or. &
       .not.checkpoint%result%final_diagnostics%eigensolver_converged .or. &
       .not.all(ieee_is_finite(diagnostic_values)) .or. &
       checkpoint%result%final_diagnostics%interface_scale/=1d0 .or. &
       checkpoint%result%final_diagnostics%orbital_residual<0d0 .or. &
       checkpoint%result%final_diagnostics%density_residual<0d0 .or. &
       checkpoint%result%final_diagnostics%subspace_residual<0d0 .or. &
       checkpoint%result%final_diagnostics%hermiticity_defect<0d0 .or. &
       checkpoint%result%final_diagnostics%orthogonality_defect<0d0 .or. &
       checkpoint%result%final_diagnostics%face_balance_defect<0d0 .or. &
       checkpoint%result%final_diagnostics%projector_overlap<0d0 .or. &
       checkpoint%result%final_diagnostics%projector_overlap>1d0 .or. &
       checkpoint%result%final_diagnostics%orbital_residual>checkpoint%controls%final_orbital_tolerance .or. &
       checkpoint%result%final_diagnostics%density_residual>checkpoint%controls%final_density_tolerance .or. &
       checkpoint%result%final_diagnostics%subspace_residual>checkpoint%controls%subspace_tolerance .or. &
       checkpoint%result%final_diagnostics%hermiticity_defect>checkpoint%controls%hermiticity_tolerance .or. &
       checkpoint%result%final_diagnostics%orthogonality_defect>checkpoint%controls%orthogonality_tolerance .or. &
       checkpoint%result%final_diagnostics%face_balance_defect>checkpoint%controls%face_balance_tolerance .or. &
       checkpoint%result%final_diagnostics%projector_overlap<checkpoint%controls%minimum_projector_overlap .or. &
       checkpoint%result%minimum_projector_overlap<checkpoint%controls%minimum_projector_overlap .or. &
       checkpoint%result%minimum_projector_overlap> &
       checkpoint%result%final_diagnostics%projector_overlap .or. &
       abs(checkpoint%result%final_diagnostics%electron_number- &
       checkpoint%result%final_diagnostics%expected_electron_number)>checkpoint%controls%electron_count_tolerance .or. &
       abs(sum(checkpoint%state%occupations)-checkpoint%result%final_diagnostics%electron_number)> &
       checkpoint%controls%electron_count_tolerance .or. &
       checkpoint%result%final_diagnostics%eigensolver_iterations<0 .or. &
       checkpoint%result%final_diagnostics%eigensolver_iterations> &
       checkpoint%controls%maximum_eigensolver_iterations .or. &
       checkpoint%result%final_diagnostics%updated_potential_epoch/= &
       checkpoint%result%final_diagnostics%hamiltonian_potential_epoch+1_int64)then
      message='DG ground-state checkpoint acceptance diagnostics are invalid';return
    end if
    if(.not.allocated(checkpoint%density) .or. .not.allocated(checkpoint%hartree) .or. &
       .not.allocated(checkpoint%vxc) .or. .not.allocated(checkpoint%vlocal))then
      message='DG ground-state checkpoint density or potentials are incomplete';return
    end if
    if(any(shape(checkpoint%density)/=shape(checkpoint%vxc)) .or. &
       any(shape(checkpoint%density)/=shape(checkpoint%vlocal)) .or. &
       any(shape(checkpoint%hartree)/=shape(checkpoint%density(:,:,:,1))) .or. &
       .not.all(ieee_is_finite(checkpoint%density)) .or. &
       .not.all(ieee_is_finite(checkpoint%hartree)) .or. &
       .not.all(ieee_is_finite(checkpoint%vxc)) .or. &
       .not.all(ieee_is_finite(checkpoint%vlocal)))then
      message='DG ground-state checkpoint payload is nonfinite or inconsistent';return
    end if
    ok=.true.;message=''
  end subroutine

  subroutine publish_dg_ground_state_checkpoint(root,checkpoint,communicator,ok,message)
    character(*),intent(in)::root
    type(s_dg_ground_state_checkpoint),intent(inout)::checkpoint
    integer,intent(in)::communicator
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,nproc,ierr,ierr_rank,ierr_size,info,manifest_info
    integer(int64)::local_checksum,generation,common_hash
    integer(int64),allocatable::common_hashes(:)
    integer(int64),allocatable::checksums(:)
    logical::local_ok,path_exists
    character(32)::rank_text
    character(:),allocatable::path
#ifdef USE_MPI
    rank=-1;nproc=-1
    call MPI_Comm_rank(communicator,rank,ierr_rank)
    call MPI_Comm_size(communicator,nproc,ierr_size)
    local_ok=ierr_rank==MPI_SUCCESS.and.ierr_size==MPI_SUCCESS.and.rank>=0.and.nproc>0
#else
    rank=0;nproc=1;local_ok=communicator>=0
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='DG ground-state checkpoint communicator failure';return;end if
    call validate_dg_ground_state_checkpoint(checkpoint,checkpoint%state%geometry_fingerprint, &
      checkpoint%state%operator_fingerprint,'DG_DC_GS',local_ok,message)
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    generation=0_int64
    if(rank==0)then
      call system_clock(count=generation)
      generation=max(generation,checkpoint%publication_generation+1_int64)
    end if
#ifdef USE_MPI
    call MPI_Bcast(generation,1,MPI_INTEGER8,0,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS.and.generation>0_int64
#else
    local_ok=generation>0_int64
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='DG ground-state checkpoint generation failed collectively';return;end if
    do
      path=generation_rank_path(root,generation,rank)
      inquire(file=path,exist=path_exists)
      call collective_and(.not.path_exists,communicator,ok)
      if(ok)exit
      if(rank==0)generation=generation+1_int64
#ifdef USE_MPI
      call MPI_Bcast(generation,1,MPI_INTEGER8,0,communicator,ierr)
      local_ok=ierr==MPI_SUCCESS
#else
      local_ok=.true.
#endif
      call collective_and(local_ok,communicator,ok)
      if(.not.ok)then;message='DG ground-state checkpoint generation retry failed';return;end if
    end do
    checkpoint%publication_generation=generation
    common_hash=common_checkpoint_fingerprint(checkpoint)
    allocate(common_hashes(nproc))
#ifdef USE_MPI
    call MPI_Allgather(common_hash,1,MPI_INTEGER8,common_hashes,1,MPI_INTEGER8,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS.and.all(common_hashes==common_hashes(1))
#else
    common_hashes(1)=common_hash;local_ok=.true.
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='DG ground-state checkpoint rank metadata disagree collectively';return;end if
    write(rank_text,'(i0)')rank
    path=generation_rank_path(root,generation,rank)
    call write_rank_payload(path,checkpoint,info)
    local_ok=info==0
    if(local_ok)call verify_rank_payload(path,checkpoint%state%geometry_fingerprint, &
      checkpoint%state%operator_fingerprint,local_checksum,info)
    local_ok=local_ok.and.info==0
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='DG ground-state checkpoint rank publication/read-back failed collectively';return;end if
    allocate(checksums(nproc))
#ifdef USE_MPI
    call MPI_Allgather(local_checksum,1,MPI_INTEGER8,checksums,1,MPI_INTEGER8,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS
#else
    checksums(1)=local_checksum;local_ok=.true.
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='DG ground-state checkpoint checksum gather failed';return;end if
    manifest_info=0
    if(rank==0)call write_manifest(trim(root)//'.manifest',generation,nproc, &
      checkpoint%state%geometry_fingerprint,checkpoint%state%operator_fingerprint,checksums,manifest_info)
#ifdef USE_MPI
    call MPI_Bcast(manifest_info,1,MPI_INTEGER,0,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS.and.manifest_info==0
#else
    local_ok=manifest_info==0
#endif
    call collective_and(local_ok,communicator,ok)
    if(ok)then
      message=''
    else
      message='DG ground-state checkpoint manifest publication failed collectively'
    end if
  end subroutine

  subroutine read_dg_ground_state_checkpoint(root,expected_geometry,expected_operator,expected_provenance, &
      communicator,checkpoint,reusable,ok,message)
    character(*),intent(in)::root,expected_provenance
    integer(int64),intent(in)::expected_geometry,expected_operator
    integer,intent(in)::communicator
    type(s_dg_ground_state_checkpoint),intent(out)::checkpoint
    logical,intent(out)::reusable,ok
    character(*),intent(out)::message
    integer::rank,nproc,ierr,ierr_rank,ierr_size,info
    integer(int64)::checksum,generation,manifest_geometry,manifest_operator
    integer(int64),allocatable::checksums(:)
    logical::local_ok,manifest_agreement
    character(32)::rank_text
    character(:),allocatable::path
#ifdef USE_MPI
    rank=-1;nproc=-1
    call MPI_Comm_rank(communicator,rank,ierr_rank)
    call MPI_Comm_size(communicator,nproc,ierr_size)
    local_ok=ierr_rank==MPI_SUCCESS.and.ierr_size==MPI_SUCCESS.and.rank>=0.and.nproc>0
#else
    rank=0;nproc=1;local_ok=communicator>=0
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;reusable=.false.;message='DG ground-state checkpoint communicator failure';return;end if
    call read_manifest(trim(root)//'.manifest',nproc,generation,manifest_geometry,manifest_operator,checksums,info)
    local_ok=info==0
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;reusable=.false.;message='DG ground-state checkpoint manifest is absent or incomplete';return;end if
    call manifest_agrees_collectively(generation,manifest_geometry,manifest_operator,checksums, &
      communicator,manifest_agreement)
    if(.not.manifest_agreement)then
      reusable=.false.;ok=.true.;message='DG ground-state checkpoint manifest generation disagrees collectively'
      return
    endif
    write(rank_text,'(i0)')rank
    path=generation_rank_path(root,generation,rank)
    call read_rank_payload(path,checkpoint,checksum,info)
    local_ok=info==0.and.checksum==checksums(rank+1)
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;reusable=.false.;message='DG ground-state checkpoint payload is corrupt collectively';return;end if
    checkpoint%publication_generation=generation
    local_ok=manifest_geometry==expected_geometry.and.manifest_operator==expected_operator
    if(local_ok)call validate_dg_ground_state_checkpoint(checkpoint,expected_geometry,expected_operator, &
      expected_provenance,local_ok,message)
    reusable=local_ok
    call collective_and(reusable,communicator,reusable)
    ok=.true.
    if(.not.reusable)message='DG ground-state checkpoint is stale or mislabeled'
  end subroutine

  subroutine write_rank_payload(path,checkpoint,info)
    character(*),intent(in)::path
    type(s_dg_ground_state_checkpoint),intent(in)::checkpoint
    integer,intent(out)::info
    character(:),allocatable::temporary
    integer::unit,ios,write_ios,flush_ios,close_ios,d(21)
    integer(int64)::checksum
    temporary=trim(path)//'.tmp'
    call checkpoint_dimensions(checkpoint,d)
    checksum=checkpoint_checksum(checkpoint)
    info=1
    open(newunit=unit,file=temporary,access='stream',form='unformatted',status='replace',action='write',iostat=ios)
    if(ios/=0)return
    write_ios=0;flush_ios=0;close_ios=0
    write(unit,iostat=write_ios)payload_magic,checkpoint%schema_version,checkpoint%publication_generation, &
      checkpoint%provenance,checkpoint%global_size, &
      checkpoint%state%fragment_id,checkpoint%state%local_basis_count, &
      checkpoint%state%global_basis_count,checkpoint%state%global_band_count,checkpoint%state%raw_basis_count, &
      checkpoint%state%core_size,checkpoint%state%full_spatial_shape,checkpoint%state%scf_iterations, &
      checkpoint%state%geometry_fingerprint,checkpoint%state%operator_fingerprint, &
      checkpoint%state%hamiltonian_operator_fingerprint, &
      checkpoint%state%density_residual,checkpoint%state%interface_scale,checkpoint%state%ready, &
      checkpoint%result%phase,checkpoint%result%naccepted_steps, &
      checkpoint%result%nrollbacks,checkpoint%result%mixing_reset_count,checkpoint%result%total_scf_iterations, &
      checkpoint%result%lambda,checkpoint%result%maximum_interface_scale, &
      checkpoint%result%minimum_projector_overlap,checkpoint%result%accepted,checkpoint%result%failed, &
      checkpoint%result%unmixed_verified,checkpoint%result%final_diagnostics,checkpoint%controls,d,checksum
    if(write_ios==0)write(unit,iostat=write_ios)checkpoint%fragment_origins,checkpoint%fragment_sizes, &
      checkpoint%face_neighbors,checkpoint%state%coefficient_rows, &
      checkpoint%state%full_fragment_basis,checkpoint%state%basis_transform, &
      checkpoint%state%eigenvalues,checkpoint%state%occupations,checkpoint%state%basis_offsets, &
      checkpoint%state%fragment_ids,checkpoint%result%lambda_history, &
      checkpoint%result%lambda_steps,checkpoint%density,checkpoint%hartree,checkpoint%vxc,checkpoint%vlocal
    if(write_ios==0)flush(unit,iostat=flush_ios)
    close(unit,iostat=close_ios)
    if(write_ios/=0.or.flush_ios/=0.or.close_ios/=0)return
    if(c_rename(trim(temporary)//c_null_char,trim(path)//c_null_char)/=0_c_int)return
    info=0
  end subroutine

  subroutine read_rank_payload(path,checkpoint,checksum,info)
    character(*),intent(in)::path
    type(s_dg_ground_state_checkpoint),intent(out)::checkpoint
    integer(int64),intent(out)::checksum
    integer,intent(out)::info
    character(20)::magic
    integer::unit,ios,d(21)
    integer(int8)::trailing
    info=1;checksum=0_int64
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,checkpoint%schema_version,checkpoint%publication_generation, &
      checkpoint%provenance,checkpoint%global_size, &
      checkpoint%state%fragment_id,checkpoint%state%local_basis_count, &
      checkpoint%state%global_basis_count,checkpoint%state%global_band_count,checkpoint%state%raw_basis_count, &
      checkpoint%state%core_size,checkpoint%state%full_spatial_shape,checkpoint%state%scf_iterations, &
      checkpoint%state%geometry_fingerprint,checkpoint%state%operator_fingerprint, &
      checkpoint%state%hamiltonian_operator_fingerprint, &
      checkpoint%state%density_residual,checkpoint%state%interface_scale,checkpoint%state%ready, &
      checkpoint%result%phase,checkpoint%result%naccepted_steps, &
      checkpoint%result%nrollbacks,checkpoint%result%mixing_reset_count,checkpoint%result%total_scf_iterations, &
      checkpoint%result%lambda,checkpoint%result%maximum_interface_scale, &
      checkpoint%result%minimum_projector_overlap,checkpoint%result%accepted,checkpoint%result%failed, &
      checkpoint%result%unmixed_verified,checkpoint%result%final_diagnostics,checkpoint%controls,d,checksum
    if(ios/=0.or.magic/=payload_magic.or..not.valid_payload_dimensions(checkpoint,d))then;close(unit);return;end if
    call allocate_payload(checkpoint,d,ios)
    if(ios/=0)then;close(unit);return;end if
    read(unit,iostat=ios)checkpoint%fragment_origins,checkpoint%fragment_sizes,checkpoint%face_neighbors, &
      checkpoint%state%coefficient_rows,checkpoint%state%full_fragment_basis, &
      checkpoint%state%basis_transform,checkpoint%state%eigenvalues,checkpoint%state%occupations, &
      checkpoint%state%basis_offsets,checkpoint%state%fragment_ids, &
      checkpoint%result%lambda_history,checkpoint%result%lambda_steps, &
      checkpoint%density,checkpoint%hartree,checkpoint%vxc,checkpoint%vlocal
    if(ios/=0)then;close(unit);return;end if
    read(unit,iostat=ios)trailing;close(unit)
    if(ios/=iostat_end.or.checkpoint_checksum(checkpoint)/=checksum)return
    info=0
  end subroutine

  subroutine verify_rank_payload(path,geometry,operator,checksum,info)
    character(*),intent(in)::path
    integer(int64),intent(in)::geometry,operator
    integer(int64),intent(out)::checksum
    integer,intent(out)::info
    type(s_dg_ground_state_checkpoint)::candidate
    logical::ok
    character(256)::message
    call read_rank_payload(path,candidate,checksum,info)
    if(info/=0)return
    call validate_dg_ground_state_checkpoint(candidate,geometry,operator,'DG_DC_GS',ok,message)
    if(.not.ok)info=1
  end subroutine

  subroutine checkpoint_dimensions(c,d)
    type(s_dg_ground_state_checkpoint),intent(in)::c
    integer,intent(out)::d(21)
    d=[size(c%fragment_origins,2),size(c%face_neighbors),size(c%result%lambda_history), &
      size(c%result%lambda_steps),shape(c%density),shape(c%hartree), &
      shape(c%state%coefficient_rows),shape(c%state%full_fragment_basis), &
      shape(c%state%basis_transform),size(c%state%eigenvalues),size(c%state%occupations), &
      size(c%state%basis_offsets),size(c%state%fragment_ids)]
  end subroutine

  subroutine allocate_payload(c,d,ios)
    type(s_dg_ground_state_checkpoint),intent(inout)::c
    integer,intent(in)::d(21)
    integer,intent(out)::ios
    allocate(c%fragment_origins(3,d(1)),c%fragment_sizes(3,d(1)),c%face_neighbors(d(2)), &
      c%result%lambda_history(d(3)),c%result%lambda_steps(d(4)), &
      c%density(d(5),d(6),d(7),d(8)),c%hartree(d(9),d(10),d(11)), &
      c%vxc(d(5),d(6),d(7),d(8)),c%vlocal(d(5),d(6),d(7),d(8)), &
      c%state%coefficient_rows(d(12),d(13)),c%state%full_fragment_basis(d(14),d(15)), &
      c%state%basis_transform(d(16),d(17)),c%state%eigenvalues(d(18)), &
      c%state%occupations(d(19)),c%state%basis_offsets(0:d(20)-1),c%state%fragment_ids(d(21)),stat=ios)
  end subroutine

  function checkpoint_checksum(c)result(checksum)
    type(s_dg_ground_state_checkpoint),intent(in)::c
    integer(int64)::checksum
    integer(int64),allocatable::words(:)
    checksum=int(z'6A09E667F3BCC909',int64)
    checksum=hash_words(checksum,[int(c%schema_version,int64),c%publication_generation, &
      int(c%global_size(1),int64),int(c%global_size(2),int64),int(c%global_size(3),int64), &
      int(c%state%fragment_id,int64),int(c%state%local_basis_count,int64), &
      int(c%state%global_basis_count,int64),int(c%state%global_band_count,int64), &
      int(c%state%raw_basis_count,int64), &
      int(c%state%core_size,int64),int(c%state%full_spatial_shape,int64), &
      int(c%state%scf_iterations,int64),c%state%geometry_fingerprint,c%state%operator_fingerprint, &
      c%state%hamiltonian_operator_fingerprint, &
      transfer(c%state%density_residual,0_int64),transfer(c%state%interface_scale,0_int64), &
      int(merge(1,0,c%state%ready),int64), &
      int(c%result%phase,int64),int(c%result%naccepted_steps,int64),int(c%result%nrollbacks,int64), &
      int(c%result%mixing_reset_count,int64),int(c%result%total_scf_iterations,int64), &
      transfer(c%result%lambda,0_int64),transfer(c%result%maximum_interface_scale,0_int64), &
      transfer(c%result%minimum_projector_overlap,0_int64),int(merge(1,0,c%result%accepted),int64), &
      int(merge(1,0,c%result%failed),int64),int(merge(1,0,c%result%unmixed_verified),int64)])
    words=transfer(c%provenance,[0_int64],storage_size(c%provenance)/64)
    checksum=hash_words(checksum,words)
    checksum=hash_words(checksum,diagnostic_words(c%result%final_diagnostics))
    checksum=hash_words(checksum,control_words(c%controls))
    block
      integer::dimensions(21)
      call checkpoint_dimensions(c,dimensions)
      checksum=hash_words(checksum,int(dimensions,int64))
    end block
    checksum=hash_words(checksum,reshape(int(c%fragment_origins,int64),[size(c%fragment_origins)]))
    checksum=hash_words(checksum,reshape(int(c%fragment_sizes,int64),[size(c%fragment_sizes)]))
    checksum=hash_words(checksum,int(c%face_neighbors,int64))
    words=transfer(c%state%coefficient_rows,[0_int64],2*size(c%state%coefficient_rows))
    checksum=hash_words(checksum,words)
    words=transfer(c%state%full_fragment_basis,[0_int64],2*size(c%state%full_fragment_basis))
    checksum=hash_words(checksum,words)
    words=transfer(c%state%basis_transform,[0_int64],2*size(c%state%basis_transform))
    checksum=hash_words(checksum,words)
    words=transfer(c%state%eigenvalues,[0_int64],size(c%state%eigenvalues));checksum=hash_words(checksum,words)
    words=transfer(c%state%occupations,[0_int64],size(c%state%occupations));checksum=hash_words(checksum,words)
    checksum=hash_words(checksum,int(c%state%basis_offsets,int64))
    checksum=hash_words(checksum,int(c%state%fragment_ids,int64))
    words=transfer(c%result%lambda_history,[0_int64],size(c%result%lambda_history));checksum=hash_words(checksum,words)
    words=transfer(c%result%lambda_steps,[0_int64],size(c%result%lambda_steps));checksum=hash_words(checksum,words)
    words=transfer(c%density,[0_int64],size(c%density));checksum=hash_words(checksum,words)
    words=transfer(c%hartree,[0_int64],size(c%hartree));checksum=hash_words(checksum,words)
    words=transfer(c%vxc,[0_int64],size(c%vxc));checksum=hash_words(checksum,words)
    words=transfer(c%vlocal,[0_int64],size(c%vlocal));checksum=hash_words(checksum,words)
    if(checksum==0_int64)checksum=1_int64
  end function

  function common_checkpoint_fingerprint(c)result(hash)
    type(s_dg_ground_state_checkpoint),intent(in)::c
    integer(int64)::hash
    integer(int64),allocatable::words(:)
    hash=int(z'510E527FADE682D1',int64)
    hash=hash_words(hash,[int(c%schema_version,int64),c%publication_generation, &
      int(c%global_size,int64),int(c%state%global_basis_count,int64), &
      int(c%state%global_band_count,int64), &
      c%state%geometry_fingerprint,c%state%operator_fingerprint, &
      c%state%hamiltonian_operator_fingerprint, &
      transfer(c%state%density_residual,0_int64),transfer(c%state%interface_scale,0_int64), &
      int(merge(1,0,c%state%ready),int64), &
      int(c%result%phase,int64),int(c%result%naccepted_steps,int64),int(c%result%nrollbacks,int64), &
      int(c%result%mixing_reset_count,int64),int(c%result%total_scf_iterations,int64), &
      transfer(c%result%lambda,0_int64),transfer(c%result%maximum_interface_scale,0_int64), &
      transfer(c%result%minimum_projector_overlap,0_int64),int(merge(1,0,c%result%accepted),int64), &
      int(merge(1,0,c%result%failed),int64),int(merge(1,0,c%result%unmixed_verified),int64)])
    hash=hash_words(hash,reshape(int(c%fragment_origins,int64),[size(c%fragment_origins)]))
    hash=hash_words(hash,reshape(int(c%fragment_sizes,int64),[size(c%fragment_sizes)]))
    hash=hash_words(hash,int(c%state%basis_offsets,int64))
    hash=hash_words(hash,int(c%state%fragment_ids,int64))
    words=transfer(c%state%eigenvalues,[0_int64],size(c%state%eigenvalues));hash=hash_words(hash,words)
    words=transfer(c%state%occupations,[0_int64],size(c%state%occupations));hash=hash_words(hash,words)
    words=transfer(c%provenance,[0_int64],storage_size(c%provenance)/64);hash=hash_words(hash,words)
    hash=hash_words(hash,diagnostic_words(c%result%final_diagnostics))
    hash=hash_words(hash,control_words(c%controls))
    words=transfer(c%result%lambda_history,[0_int64],size(c%result%lambda_history));hash=hash_words(hash,words)
    words=transfer(c%result%lambda_steps,[0_int64],size(c%result%lambda_steps));hash=hash_words(hash,words)
    if(hash==0_int64)hash=1_int64
  end function

  integer(int64) function hash_words(seed,words)result(hash)
    integer(int64),intent(in)::seed,words(:)
    integer::i
    hash=seed
    do i=1,size(words)
      hash=ieor(ishftc(hash,13),words(i))
      hash=ieor(hash,ishft(hash,-7))
    end do
  end function

  pure function diagnostic_words(d)result(words)
    type(s_dg_dc_gs_diagnostics),intent(in)::d
    integer(int64)::words(15)
    words(1:10)=transfer(diagnostic_vector(d),words(1:10))
    words(11:15)=[int(d%eigensolver_iterations,int64),d%hamiltonian_potential_epoch, &
      d%updated_potential_epoch,int(merge(1,0,d%eigensolver_converged),int64), &
      int(merge(1,0,d%finite),int64)]
  end function

  pure function control_words(c)result(words)
    type(s_dg_dc_gs_controls),intent(in)::c
    integer(int64)::words(18)
    real(8)::values(15)
    values=[c%intermediate_orbital_tolerance,c%intermediate_density_tolerance, &
      c%final_orbital_tolerance,c%final_density_tolerance,c%subspace_tolerance, &
      c%initial_lambda_step,c%minimum_lambda_step,c%maximum_lambda_step,c%allowed_residual_growth, &
      c%density_mix_rate,c%hermiticity_tolerance,c%orthogonality_tolerance,c%face_balance_tolerance, &
      c%electron_count_tolerance,c%minimum_projector_overlap]
    words(1:15)=transfer(values,words(1:15))
    words(16:18)=[int(c%maximum_scf_iterations,int64),int(c%maximum_eigensolver_iterations,int64), &
      int(c%maximum_rollbacks,int64)]
  end function

  logical function valid_payload_dimensions(c,d)result(ok)
    type(s_dg_ground_state_checkpoint),intent(in)::c
    integer,intent(in)::d(21)
    real(8)::elements,bytes
    ok=all(d>0).and.all(d<=maximum_dimension).and.c%state%fragment_id>=1.and.c%state%fragment_id<=d(1)
    ok=ok.and.d(3)==d(4).and.d(5)==d(9).and.d(6)==d(10).and.d(7)==d(11)
    ok=ok.and.d(12)==c%state%local_basis_count.and.d(13)==c%state%global_band_count
    ok=ok.and.int(d(14),int64)==safe_product3(c%state%full_spatial_shape)
    ok=ok.and.d(15)==c%state%local_basis_count.and.d(16)==c%state%raw_basis_count
    ok=ok.and.d(17)==c%state%local_basis_count
    ok=ok.and.d(18)==c%state%global_band_count.and.d(19)==c%state%global_band_count
    ok=ok.and.d(20)==d(21)+1.and.d(8)>0.and.d(1)>0.and.d(2)==nodal_face_count
    if(.not.ok)return
    elements=2d0*(real(d(12),8)*real(d(13),8)+real(d(14),8)*real(d(15),8)+ &
      real(d(16),8)*real(d(17),8))+ &
      real(d(18)+d(19)+d(20)+d(21),8)+3d0*product(real(d(5:8),8))+ &
      product(real(d(9:11),8))+6d0*real(d(1),8)+real(d(2)+d(3)+d(4),8)
    bytes=8d0*elements
    ok=ieee_is_finite(bytes).and.bytes>0d0.and.bytes<=8d0*1024d0**3
  end function

  pure integer(int64) function safe_product3(values)result(product_value)
    integer,intent(in)::values(3)
    integer::i
    product_value=1_int64
    do i=1,3
      if(values(i)<=0 .or. product_value>huge(product_value)/int(values(i),int64))then
        product_value=-1_int64
        return
      end if
      product_value=product_value*int(values(i),int64)
    end do
  end function

  logical function face_neighbors_match_geometry(c)result(ok)
    type(s_dg_ground_state_checkpoint),intent(in)::c
    integer::axis,side,iface,fragment,candidate,other_axis,matches
    integer::boundary,current_origin,current_end,candidate_origin,candidate_end
    logical::same_transverse_extent
    fragment=c%state%fragment_id
    ok=.true.
    do axis=1,3
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        matches=0
        current_origin=c%fragment_origins(axis,fragment)
        current_end=current_origin+c%fragment_sizes(axis,fragment)
        if(side<0)then
          boundary=modulo(current_origin,c%global_size(axis))
        else
          boundary=modulo(current_end,c%global_size(axis))
        end if
        do candidate=1,size(c%fragment_origins,2)
          candidate_origin=c%fragment_origins(axis,candidate)
          candidate_end=candidate_origin+c%fragment_sizes(axis,candidate)
          if(side<0)then
            if(modulo(candidate_end,c%global_size(axis))/=boundary)cycle
          else
            if(modulo(candidate_origin,c%global_size(axis))/=boundary)cycle
          end if
          same_transverse_extent=.true.
          do other_axis=1,3
            if(other_axis==axis)cycle
            same_transverse_extent=same_transverse_extent .and. &
              c%fragment_origins(other_axis,candidate)==c%fragment_origins(other_axis,fragment) .and. &
              c%fragment_sizes(other_axis,candidate)==c%fragment_sizes(other_axis,fragment)
          end do
          if(.not.same_transverse_extent)cycle
          matches=matches+1
          if(c%face_neighbors(iface)/=candidate-1)ok=.false.
        end do
        if(matches/=1)ok=.false.
      end do
    end do
  end function

  logical function valid_fragment_basis_layout(state)result(ok)
    type(s_dg_dc_local_basis_production_state),intent(in)::state
    logical,allocatable::seen(:)
    integer::i,position
    allocate(seen(size(state%fragment_ids)))
    seen=.false.
    position=0
    ok=.true.
    do i=1,size(state%fragment_ids)
      if(state%fragment_ids(i)<1 .or. state%fragment_ids(i)>size(seen))then
        ok=.false.;return
      end if
      if(seen(state%fragment_ids(i)))then
        ok=.false.;return
      end if
      seen(state%fragment_ids(i))=.true.
      if(state%fragment_ids(i)==state%fragment_id)position=i
    end do
    ok=all(seen) .and. position>0
    if(ok)ok=state%basis_offsets(position)-state%basis_offsets(position-1)==state%local_basis_count
  end function

  logical function valid_checkpoint_controls(c)result(ok)
    type(s_dg_dc_gs_controls),intent(in)::c
    integer(int64)::words(18)
    real(8)::values(15)
    words=control_words(c)
    values=transfer(words(1:15),values)
    ok=all(ieee_is_finite(values)).and.all(values>0d0).and.c%density_mix_rate<=1d0.and. &
      c%minimum_lambda_step<=c%initial_lambda_step.and.c%initial_lambda_step<=c%maximum_lambda_step.and. &
      c%maximum_lambda_step<=1d0.and.c%minimum_projector_overlap<=1d0.and. &
      c%maximum_scf_iterations>0.and.c%maximum_eigensolver_iterations>0.and.c%maximum_rollbacks>=0
  end function

  function generation_rank_path(root,generation,rank)result(path)
    character(*),intent(in)::root
    integer(int64),intent(in)::generation
    integer,intent(in)::rank
    character(:),allocatable::path
    character(64)::generation_text,rank_text
    write(generation_text,'(i0)')generation
    write(rank_text,'(i0)')rank
    path=trim(root)//'.g'//trim(generation_text)//'.rank'//trim(rank_text)//'.dg_gs'
  end function

  pure function diagnostic_vector(d)result(values)
    type(s_dg_dc_gs_diagnostics),intent(in)::d
    real(8)::values(10)
    values=[d%orbital_residual,d%density_residual,d%subspace_residual,d%projector_overlap, &
      d%hermiticity_defect,d%orthogonality_defect,d%face_balance_defect,d%electron_number, &
      d%expected_electron_number,d%interface_scale]
  end function

  subroutine write_manifest(path,generation,nrank,geometry,operator,checksums,info)
    character(*),intent(in)::path
    integer,intent(in)::nrank
    integer(int64),intent(in)::generation,geometry,operator,checksums(:)
    integer,intent(out)::info
    character(:),allocatable::temporary
    integer::unit,ios,write_ios,flush_ios,close_ios
    temporary=trim(path)//'.tmp';info=1
    open(newunit=unit,file=temporary,access='stream',form='unformatted',status='replace',action='write',iostat=ios)
    if(ios/=0)return
    write_ios=0;flush_ios=0;close_ios=0
    write(unit,iostat=write_ios)manifest_magic,dg_ground_state_checkpoint_version,generation,nrank, &
      geometry,operator,checksums
    if(write_ios==0)flush(unit,iostat=flush_ios)
    close(unit,iostat=close_ios)
    if(write_ios/=0.or.flush_ios/=0.or.close_ios/=0)return
    if(c_rename(trim(temporary)//c_null_char,trim(path)//c_null_char)/=0_c_int)return
    info=0
  end subroutine

  subroutine read_manifest(path,expected_nrank,generation,geometry,operator,checksums,info)
    character(*),intent(in)::path
    integer,intent(in)::expected_nrank
    integer(int64),intent(out)::generation,geometry,operator
    integer(int64),allocatable,intent(out)::checksums(:)
    integer,intent(out)::info
    character(18)::magic
    integer::unit,ios,version,nrank
    integer(int8)::trailing
    info=1
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,version,generation,nrank,geometry,operator
    if(ios/=0.or.magic/=manifest_magic.or.version/=dg_ground_state_checkpoint_version.or. &
       generation<=0_int64.or.nrank/=expected_nrank.or.nrank<=0)then;close(unit);return;end if
    allocate(checksums(nrank),stat=ios)
    if(ios/=0)then;close(unit);return;end if
    read(unit,iostat=ios)checksums
    if(ios/=0)then;close(unit);return;end if
    read(unit,iostat=ios)trailing;close(unit)
    if(ios/=iostat_end.or.any(checksums==0_int64))return
    info=0
  end subroutine

  subroutine manifest_agrees_collectively(generation,geometry,operator,checksums,communicator,agrees)
    integer(int64),intent(in)::generation,geometry,operator,checksums(:)
    integer,intent(in)::communicator
    logical,intent(out)::agrees
    integer(int64)::local_values(4),minimum_values(4),maximum_values(4)
    integer::ierr
    local_values=[generation,geometry,operator, &
      hash_words(int(z'13198A2E03707344',int64),checksums)]
#ifdef USE_MPI
    call MPI_Allreduce(local_values,minimum_values,4,MPI_INTEGER8,MPI_MIN,communicator,ierr)
    agrees=ierr==MPI_SUCCESS
    call MPI_Allreduce(local_values,maximum_values,4,MPI_INTEGER8,MPI_MAX,communicator,ierr)
    agrees=agrees.and.ierr==MPI_SUCCESS.and.all(minimum_values==maximum_values)
#else
    minimum_values=local_values;maximum_values=local_values
    agrees=communicator>=0
#endif
  end subroutine manifest_agrees_collectively

  subroutine collective_and(local_value,communicator,global_value)
    logical,intent(in)::local_value
    integer,intent(in)::communicator
    logical,intent(out)::global_value
#ifdef USE_MPI
    integer::ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS)global_value=.false.
#else
    global_value=local_value.and.communicator>=0
#endif
  end subroutine

  subroutine publish_dg_dc_direct_checkpoint(root,state,communicator,ok,message)
    character(*),intent(in)::root
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    integer,intent(in)::communicator
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(20),parameter::magic='DG_DC_DIRECT_STATE1 '
    character(64)::rank_text
    character(:),allocatable::path,temporary
    integer::rank,nproc,ierr,unit,ios,close_ios,rename_info,manifest_info
    integer(int64)::metadata_hash,minimum_hash,maximum_hash,checksum,generation
    integer(int64),allocatable::checksums(:)
    logical::local_ok,path_exists,unit_open
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    local_ok=ierr==MPI_SUCCESS
    call MPI_Comm_size(communicator,nproc,ierr)
    local_ok=local_ok.and.ierr==MPI_SUCCESS
#else
    rank=0;nproc=1;local_ok=communicator>=0
#endif
    local_ok=local_ok.and.valid_direct_checkpoint_state(state)
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint state is invalid collectively';return;endif
    call validate_direct_checkpoint_topology(state,communicator,local_ok)
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint topology is invalid collectively';return;endif
    metadata_hash=direct_checkpoint_metadata_hash(state)
#ifdef USE_MPI
    call MPI_Allreduce(metadata_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS
    call MPI_Allreduce(metadata_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,communicator,ierr)
    local_ok=local_ok.and.ierr==MPI_SUCCESS
#else
    minimum_hash=metadata_hash;maximum_hash=metadata_hash;local_ok=.true.
#endif
    local_ok=local_ok.and.minimum_hash==maximum_hash
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint rank metadata disagree collectively';return;endif
    generation=0_int64
    if(rank==0)call system_clock(count=generation)
#ifdef USE_MPI
    call MPI_Bcast(generation,1,MPI_INTEGER8,0,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS.and.generation>0_int64
#else
    local_ok=generation>0_int64
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint generation failed collectively';return;endif
    do
      path=direct_generation_rank_path(root,generation,rank)
      inquire(file=path,exist=path_exists)
      call collective_and(.not.path_exists,communicator,ok)
      if(ok)exit
      if(rank==0)generation=generation+1_int64
#ifdef USE_MPI
      call MPI_Bcast(generation,1,MPI_INTEGER8,0,communicator,ierr)
      local_ok=ierr==MPI_SUCCESS
#else
      local_ok=.true.
#endif
      call collective_and(local_ok,communicator,ok)
      if(.not.ok)then;message='direct DG checkpoint generation retry failed collectively';return;endif
    enddo
    write(rank_text,'(i0)')rank
    path=direct_generation_rank_path(root,generation,rank)
    temporary=path//'.tmp'
    checksum=direct_checkpoint_checksum(state)
    ios=0
    open(newunit=unit,file=temporary,access='stream',form='unformatted',status='replace',action='write',iostat=ios)
    unit_open=ios==0
    if(ios==0)write(unit,iostat=ios)magic,dg_ground_state_checkpoint_version,nproc,state%fragment_id, &
      state%nstate,state%state_start,state%state_count,state%nspin,state%scf_iterations, &
      state%accepted,state%unmixed_verified,state%core_size,state%full_spatial_shape, &
      state%orbital_spatial_lower_bound,state%orbital_core_local_origin, &
      state%orbital_core_origin,state%orbital_core_size, &
      state%global_size,state%fragment_origin,state%fragment_size,state%face_neighbors, &
      state%density_origin,state%density_size,state%continuation_rollbacks,state%projector_generation, &
      state%projector_retained_rank,state%projector_required_rank, &
      state%geometry_fingerprint,state%operator_fingerprint,state%hamiltonian_operator_fingerprint, &
      state%projector_fingerprint, &
      state%density_residual,state%interface_scale,state%orbital_residual,state%orthogonality_defect, &
      state%hermiticity_defect,state%charge_error,state%face_balance_defect, &
      state%density_tolerance,state%orbital_tolerance,state%orthogonality_tolerance, &
      state%hermiticity_tolerance,state%charge_tolerance,state%face_balance_tolerance, &
      state%projector_projection_residual,state%projector_escape_norm,state%face_action_norm, &
      state%projector_residual_limit,state%projector_escape_limit, &
      state%face_pair_balance, &
      shape(state%fragment_orbitals),shape(state%occupations),size(state%continuation_scale_history), &
      shape(state%density),checksum
    if(ios==0)write(unit,iostat=ios)state%fragment_orbitals,state%occupations, &
      state%continuation_scale_history,state%continuation_step_history,state%density,state%hartree,state%vxc,state%vlocal
    if(ios==0)flush(unit,iostat=ios)
    close_ios=0
    if(unit_open)close(unit,iostat=close_ios)
    if(ios/=0.or.close_ios/=0)then
      local_ok=.false.
    else
      rename_info=c_rename(trim(temporary)//c_null_char,trim(path)//c_null_char)
      local_ok=rename_info==0
    endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint publication failed collectively';return;endif
    allocate(checksums(nproc))
#ifdef USE_MPI
    call MPI_Allgather(checksum,1,MPI_INTEGER8,checksums,1,MPI_INTEGER8,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS
#else
    checksums(1)=checksum;local_ok=.true.
#endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then;message='direct DG checkpoint checksum gather failed collectively';return;endif
    manifest_info=0
    if(rank==0)call write_manifest(trim(root)//'.manifest',generation,nproc,state%geometry_fingerprint, &
      state%operator_fingerprint,checksums,manifest_info)
#ifdef USE_MPI
    call MPI_Bcast(manifest_info,1,MPI_INTEGER,0,communicator,ierr)
    local_ok=ierr==MPI_SUCCESS.and.manifest_info==0
#else
    local_ok=manifest_info==0
#endif
    call collective_and(local_ok,communicator,ok)
    if(ok)then;message='';else;message='direct DG checkpoint manifest publication failed collectively';endif
  end subroutine publish_dg_dc_direct_checkpoint

  subroutine read_dg_dc_direct_checkpoint(root,expected_geometry,expected_operator,communicator, &
      state,reusable,ok,message)
    character(*),intent(in)::root
    integer(int64),intent(in)::expected_geometry,expected_operator
    integer,intent(in)::communicator
    type(s_dg_dc_direct_checkpoint_state),intent(out)::state
    logical,intent(out)::reusable,ok
    character(*),intent(out)::message
    character(20),parameter::expected_magic='DG_DC_DIRECT_STATE1 '
    character(20)::magic
    character(64)::rank_text
    character(:),allocatable::path
    integer::rank,nproc,file_nproc,version,unit,ios,ierr,manifest_info,inquire_ios
    integer::orbital_shape(7),occupation_shape(3),field_shape(4),history_size
    integer(int64)::checksum,generation,manifest_geometry,manifest_operator
    integer(int64),allocatable::checksums(:)
    integer(int64)::file_position,file_size
    logical::local_ok,manifest_agreement,unit_open
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    local_ok=ierr==MPI_SUCCESS
    call MPI_Comm_size(communicator,nproc,ierr)
    local_ok=local_ok.and.ierr==MPI_SUCCESS
#else
    rank=0;nproc=1;local_ok=communicator>=0
#endif
    generation=0_int64
    call read_manifest(trim(root)//'.manifest',nproc,generation,manifest_geometry,manifest_operator, &
      checksums,manifest_info)
    local_ok=local_ok.and.manifest_info==0
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then
      reusable=.false.;ok=.true.;message='direct DG checkpoint manifest is absent or incomplete collectively'
      return
    endif
    call manifest_agrees_collectively(generation,manifest_geometry,manifest_operator,checksums, &
      communicator,manifest_agreement)
    if(.not.manifest_agreement)then
      reusable=.false.;ok=.true.;message='direct DG checkpoint manifest generation disagrees collectively'
      return
    endif
    write(rank_text,'(i0)')rank
    path=direct_generation_rank_path(root,generation,rank)
    ios=1
    if(local_ok)open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    unit_open=ios==0
    local_ok=local_ok.and.ios==0
    if(local_ok)read(unit,iostat=ios)magic,version,file_nproc,state%fragment_id,state%nstate, &
      state%state_start,state%state_count,state%nspin,state%scf_iterations,state%accepted, &
      state%unmixed_verified,state%core_size,state%full_spatial_shape,state%orbital_spatial_lower_bound, &
      state%orbital_core_local_origin,state%orbital_core_origin,state%orbital_core_size, &
      state%global_size,state%fragment_origin, &
      state%fragment_size,state%face_neighbors,state%density_origin,state%density_size, &
      state%continuation_rollbacks,state%projector_generation,state%projector_retained_rank, &
      state%projector_required_rank, &
      state%geometry_fingerprint, &
      state%operator_fingerprint,state%hamiltonian_operator_fingerprint, &
      state%projector_fingerprint, &
      state%density_residual,state%interface_scale,state%orbital_residual,state%orthogonality_defect, &
      state%hermiticity_defect,state%charge_error,state%face_balance_defect, &
      state%density_tolerance,state%orbital_tolerance,state%orthogonality_tolerance, &
      state%hermiticity_tolerance,state%charge_tolerance,state%face_balance_tolerance, &
      state%projector_projection_residual,state%projector_escape_norm,state%face_action_norm, &
      state%projector_residual_limit,state%projector_escape_limit, &
      state%face_pair_balance,orbital_shape,occupation_shape, &
      history_size,field_shape,checksum
    local_ok=local_ok.and.ios==0.and.magic==expected_magic.and.version==dg_ground_state_checkpoint_version.and. &
      file_nproc==nproc.and.all(orbital_shape>0).and.all(occupation_shape>0).and. &
      all(orbital_shape<=maximum_dimension).and.all(occupation_shape<=maximum_dimension).and. &
      safe_shape_product(orbital_shape)<=int(maximum_dimension,int64).and. &
      safe_shape_product(occupation_shape)<=int(maximum_dimension,int64).and. &
      history_size>1.and.history_size<=maximum_dimension.and.all(field_shape>0).and. &
      all(field_shape<=maximum_dimension).and.safe_shape_product(field_shape)<=int(maximum_dimension,int64)
    if(local_ok)then
      allocate(state%fragment_orbitals(orbital_shape(1),orbital_shape(2),orbital_shape(3),orbital_shape(4), &
        orbital_shape(5),orbital_shape(6),orbital_shape(7)),state%occupations(occupation_shape(1), &
        occupation_shape(2),occupation_shape(3)),state%continuation_scale_history(history_size), &
        state%continuation_step_history(history_size),state%density(field_shape(1),field_shape(2), &
        field_shape(3),field_shape(4)),state%hartree(field_shape(1),field_shape(2),field_shape(3)), &
        state%vxc(field_shape(1),field_shape(2),field_shape(3),field_shape(4)), &
        state%vlocal(field_shape(1),field_shape(2),field_shape(3),field_shape(4)),stat=ios)
      local_ok=ios==0
    endif
    if(local_ok)read(unit,iostat=ios)state%fragment_orbitals,state%occupations, &
      state%continuation_scale_history,state%continuation_step_history,state%density,state%hartree,state%vxc,state%vlocal
    local_ok=local_ok.and.ios==0
    if(local_ok)then
      inquire(unit,pos=file_position,size=file_size,iostat=inquire_ios)
      local_ok=inquire_ios==0.and.file_position==file_size+1_int64
    endif
    if(local_ok)state%ready=.true.
    if(local_ok)local_ok=valid_direct_checkpoint_state(state).and. &
      state%geometry_fingerprint==expected_geometry.and.state%operator_fingerprint==expected_operator.and. &
      direct_checkpoint_checksum(state)==checksum.and.checksum==checksums(rank+1).and. &
      manifest_geometry==expected_geometry.and.manifest_operator==expected_operator
    if(unit_open)close(unit)
    call collective_and(local_ok,communicator,reusable)
    if(reusable)then
      call validate_direct_checkpoint_topology(state,communicator,local_ok)
      call collective_and(local_ok,communicator,reusable)
    endif
    ok=.true.
    if(reusable)then
      message=''
    else
      message='direct DG checkpoint is corrupt, stale, or incomplete collectively'
    endif
  end subroutine read_dg_dc_direct_checkpoint

  function direct_generation_rank_path(root,generation,rank)result(path)
    character(*),intent(in)::root
    integer(int64),intent(in)::generation
    integer,intent(in)::rank
    character(:),allocatable::path
    character(64)::generation_text,rank_text
    write(generation_text,'(i0)')generation
    write(rank_text,'(i0)')rank
    path=trim(root)//'.g'//trim(generation_text)//'.rank'//trim(rank_text)//'.dg_direct'
  end function direct_generation_rank_path

  logical function valid_direct_checkpoint_state(state)result(valid)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    valid=state%ready.and.state%accepted.and.state%unmixed_verified.and. &
      state%fragment_id>0.and.state%nstate>0.and.state%state_start>0.and.state%state_count>0.and. &
      state%state_start+state%state_count-1<=state%nstate.and.state%nspin>0.and. &
      state%scf_iterations>0.and.all(state%core_size>0).and. &
      all(state%full_spatial_shape>0).and.state%geometry_fingerprint/=0_int64.and. &
      all(state%orbital_core_origin>=0).and.all(state%orbital_core_size>=0).and. &
      state%operator_fingerprint/=0_int64.and.state%hamiltonian_operator_fingerprint/=0_int64.and. &
      state%projector_fingerprint/=0_int64.and.state%projector_generation>0.and. &
      all(state%projector_required_rank>0).and. &
      all(state%projector_retained_rank>=state%projector_required_rank).and. &
      state%interface_scale==1d0.and.state%continuation_rollbacks>=0.and.all(state%global_size>0).and. &
      all(state%fragment_origin>=0).and.all(state%fragment_size==state%core_size).and. &
      all(state%fragment_origin+state%fragment_size<=state%global_size).and.all(state%face_neighbors>0).and. &
      all(state%density_origin>=0).and.all(state%density_size>0).and. &
      all(state%density_origin+state%density_size<=state%global_size).and. &
      all(ieee_is_finite([state%density_residual,state%orbital_residual,state%orthogonality_defect, &
      state%hermiticity_defect,state%charge_error,state%face_balance_defect])).and. &
      all([state%density_residual,state%orbital_residual,state%orthogonality_defect, &
      state%hermiticity_defect,state%charge_error,state%face_balance_defect]>=0d0).and. &
      all(ieee_is_finite([state%density_tolerance,state%orbital_tolerance,state%orthogonality_tolerance, &
      state%hermiticity_tolerance,state%charge_tolerance,state%face_balance_tolerance])).and. &
      all([state%density_tolerance,state%orbital_tolerance,state%orthogonality_tolerance, &
      state%hermiticity_tolerance,state%charge_tolerance,state%face_balance_tolerance]>0d0).and. &
      state%density_residual<=state%density_tolerance.and. &
      state%orbital_residual<=state%orbital_tolerance.and. &
      state%orthogonality_defect<=state%orthogonality_tolerance.and. &
      state%hermiticity_defect<=state%hermiticity_tolerance.and. &
      state%charge_error<=state%charge_tolerance.and. &
      state%face_balance_defect<=state%face_balance_tolerance.and. &
      all(ieee_is_finite(state%projector_projection_residual)).and. &
      all(state%projector_projection_residual>=0d0).and. &
      all(ieee_is_finite(state%projector_escape_norm)).and.all(state%projector_escape_norm>=0d0).and. &
      all(ieee_is_finite(state%projector_residual_limit)).and.all(state%projector_residual_limit>0d0).and. &
      all(ieee_is_finite(state%projector_escape_limit)).and.all(state%projector_escape_limit>0d0).and. &
      all(state%projector_projection_residual<=state%projector_residual_limit).and. &
      all(state%projector_escape_norm<=state%projector_escape_limit).and. &
      all(ieee_is_finite(state%face_action_norm)).and.all(state%face_action_norm>=0d0).and. &
      all(ieee_is_finite(state%face_pair_balance)).and.all(state%face_pair_balance>=0d0).and. &
      all(state%face_pair_balance<=state%face_balance_tolerance).and. &
      allocated(state%fragment_orbitals).and.allocated(state%occupations).and. &
      allocated(state%continuation_scale_history).and.allocated(state%continuation_step_history).and. &
      allocated(state%density).and.allocated(state%hartree).and.allocated(state%vxc).and.allocated(state%vlocal)
    if(.not.valid)return
    if(all(state%orbital_core_size>0))valid= &
      all(state%orbital_spatial_lower_bound<=state%orbital_core_local_origin).and. &
      all(state%orbital_core_local_origin+state%orbital_core_size<= &
        state%orbital_spatial_lower_bound+state%full_spatial_shape)
    if(.not.valid)return
    valid=all([size(state%fragment_orbitals,1),size(state%fragment_orbitals,2), &
      size(state%fragment_orbitals,3)]==state%full_spatial_shape).and. &
      size(state%fragment_orbitals,4)==state%nspin.and. &
      size(state%fragment_orbitals,5)==state%state_count.and. &
      size(state%fragment_orbitals,6)==1.and.size(state%fragment_orbitals,7)==1.and. &
      size(state%occupations,1)==state%state_count.and.size(state%occupations,2)==1.and. &
      size(state%occupations,3)==state%nspin.and. &
      size(state%continuation_scale_history)==size(state%continuation_step_history).and. &
      size(state%continuation_scale_history)>1.and.state%continuation_scale_history(1)==0d0.and. &
      state%continuation_scale_history(size(state%continuation_scale_history))==1d0.and. &
      state%continuation_step_history(1)==0d0.and. &
      all(state%continuation_scale_history(2:)>=state%continuation_scale_history(: &
      size(state%continuation_scale_history)-1)).and. &
      all(state%continuation_step_history(2:)==state%continuation_scale_history(2:)- &
      state%continuation_scale_history(:size(state%continuation_scale_history)-1)).and. &
      all(shape(state%density)==[state%density_size,state%nspin]).and. &
      all(shape(state%hartree)==state%density_size).and.all(shape(state%vxc)==shape(state%density)).and. &
      all(shape(state%vlocal)==shape(state%density)).and. &
      all(ieee_is_finite(state%fragment_orbitals)).and.all(ieee_is_finite(state%occupations)).and. &
      all(state%occupations>=0d0).and.all(ieee_is_finite(state%continuation_scale_history)).and. &
      all(ieee_is_finite(state%continuation_step_history)).and.all(ieee_is_finite(state%density)).and. &
      all(ieee_is_finite(state%hartree)).and.all(ieee_is_finite(state%vxc)).and.all(ieee_is_finite(state%vlocal))
  end function valid_direct_checkpoint_state

  subroutine validate_direct_checkpoint_topology(state,communicator,valid)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    integer,intent(in)::communicator
    logical,intent(out)::valid
    integer::local_metadata(39),nproc,nfragment,rank,ierr,i,j,k,iface,axis,component,neighbor_column,opposite
    integer::state_count_sum
    integer,allocatable::metadata(:,:)
    integer(int64)::fragment_volume,density_volume,orbital_core_volume,global_volume
    integer(int64)::local_field_hash,local_occupation_hash
    integer(int64),allocatable::field_hashes(:),occupation_hashes(:)
    logical::found
#ifdef USE_MPI
    call MPI_Comm_size(communicator,nproc,ierr)
    valid=ierr==MPI_SUCCESS.and.nproc>0
    call MPI_Comm_rank(communicator,rank,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
#else
    nproc=1;rank=0;valid=communicator>=0
#endif
    if(.not.valid)return
    local_metadata=[state%fragment_id,state%global_size,state%fragment_origin,state%fragment_size, &
      state%density_origin,state%density_size,state%face_neighbors,state%state_start,state%state_count, &
      state%orbital_spatial_lower_bound,state%full_spatial_shape,state%orbital_core_local_origin, &
      state%orbital_core_origin,state%orbital_core_size]
    allocate(metadata(39,nproc),field_hashes(nproc),occupation_hashes(nproc))
    local_field_hash=direct_checkpoint_field_hash(state)
    local_occupation_hash=direct_checkpoint_occupation_hash(state)
#ifdef USE_MPI
    call MPI_Allgather(local_metadata,39,MPI_INTEGER,metadata,39,MPI_INTEGER,communicator,ierr)
    valid=ierr==MPI_SUCCESS
    call MPI_Allgather(local_field_hash,1,MPI_INTEGER8,field_hashes,1,MPI_INTEGER8,communicator,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
    call MPI_Allgather(local_occupation_hash,1,MPI_INTEGER8,occupation_hashes,1,MPI_INTEGER8,communicator,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
#else
    metadata(:,1)=local_metadata
    field_hashes(1)=local_field_hash
    occupation_hashes(1)=local_occupation_hash
#endif
    if(.not.valid)return
    valid=all(metadata(2:4,:)==spread(metadata(2:4,1),2,nproc))
    nfragment=maxval(metadata(1,:))
    global_volume=product(int(metadata(2:4,1),int64))
    fragment_volume=0_int64;density_volume=0_int64;orbital_core_volume=0_int64
    do i=1,nproc
      valid=valid.and.metadata(1,i)>=1.and.metadata(1,i)<=nfragment
      valid=valid.and.all(metadata(5:7,i)>=0).and.all(metadata(8:10,i)>0).and. &
        all(metadata(5:7,i)+metadata(8:10,i)<=metadata(2:4,i))
      valid=valid.and.all(metadata(11:13,i)>=0).and.all(metadata(14:16,i)>0).and. &
        all(metadata(11:13,i)+metadata(14:16,i)<=metadata(2:4,i))
      valid=valid.and.all(metadata(17:22,i)>=1).and.all(metadata(17:22,i)<=nfragment)
      if(findloc(metadata(1,:),metadata(1,i),dim=1)==i) &
        fragment_volume=fragment_volume+product(int(metadata(8:10,i),int64))
      valid=valid.and.metadata(23,i)>0.and.metadata(24,i)>0.and. &
        metadata(23,i)+metadata(24,i)-1<=state%nstate.and.all(metadata(28:30,i)>0)
      valid=valid.and.all(metadata(34:36,i)>=metadata(5:7,i)).and. &
        all(metadata(34:36,i)+metadata(37:39,i)<=metadata(5:7,i)+metadata(8:10,i)).and. &
        all(metadata(37:39,i)>=0)
      if(all(metadata(37:39,i)>0))valid=valid.and. &
        all(metadata(25:27,i)<=metadata(31:33,i)).and. &
        all(metadata(31:33,i)+metadata(37:39,i)<=metadata(25:27,i)+metadata(28:30,i))
      found=.false.
      do j=1,i-1
        if(all(metadata(11:16,j)==metadata(11:16,i)))found=.true.
      enddo
      if(.not.found)density_volume=density_volume+product(int(metadata(14:16,i),int64))
      found=.false.
      do j=1,i-1
        if(metadata(1,j)==metadata(1,i).and.all(metadata(34:39,j)==metadata(34:39,i)))found=.true.
      enddo
      if(.not.found)orbital_core_volume=orbital_core_volume+product(int(metadata(37:39,i),int64))
      do j=i+1,nproc
        if(metadata(1,i)==metadata(1,j))then
          valid=valid.and.all(metadata(5:10,i)==metadata(5:10,j)).and. &
            all(metadata(17:22,i)==metadata(17:22,j))
        else
          valid=valid.and..not.integer_boxes_overlap(metadata(5:7,i),metadata(8:10,i), &
            metadata(5:7,j),metadata(8:10,j))
        endif
        if(integer_boxes_overlap(metadata(11:13,i),metadata(14:16,i), &
          metadata(11:13,j),metadata(14:16,j)))then
          valid=valid.and.all(metadata(11:16,i)==metadata(11:16,j))
          if(all(metadata(11:16,i)==metadata(11:16,j))) &
            valid=valid.and.field_hashes(i)==field_hashes(j)
        endif
        if(metadata(1,i)==metadata(1,j).and.all(metadata(23:24,i)==metadata(23:24,j))) &
          valid=valid.and.occupation_hashes(i)==occupation_hashes(j)
        if(all(metadata(37:39,i)>0).and.all(metadata(37:39,j)>0))then
          if(metadata(1,i)==metadata(1,j).and.all(metadata(23:24,i)==metadata(23:24,j)).and. &
             all(metadata(34:39,i)==metadata(34:39,j)))valid=.false.
          if(integer_boxes_overlap(metadata(34:36,i),metadata(37:39,i), &
            metadata(34:36,j),metadata(37:39,j))) &
            valid=valid.and.metadata(1,i)==metadata(1,j).and.all(metadata(34:39,i)==metadata(34:39,j))
        endif
      enddo
      state_count_sum=0
      do j=1,nproc
        if(metadata(1,j)==metadata(1,i).and.all(metadata(25:39,j)==metadata(25:39,i)))then
          state_count_sum=state_count_sum+metadata(24,j)
          if(j/=i)valid=valid.and. &
            (metadata(23,i)+metadata(24,i)<=metadata(23,j).or. &
             metadata(23,j)+metadata(24,j)<=metadata(23,i))
        endif
      enddo
      valid=valid.and.state_count_sum==state%nstate
      do j=1,nproc
        if(metadata(1,j)/=metadata(1,i).or. &
           all(metadata(23:24,j)==metadata(23:24,i)))cycle
        found=.false.
        do k=1,nproc
          if(metadata(1,k)==metadata(1,i).and.all(metadata(23:24,k)==metadata(23:24,j)).and. &
             all(metadata(25:39,k)==metadata(25:39,i)))found=.true.
        enddo
        valid=valid.and.found
      enddo
      do iface=1,6
        neighbor_column=findloc(metadata(1,:),metadata(16+iface,i),dim=1)
        valid=valid.and.neighbor_column>0
        if(neighbor_column>0)then
          opposite=merge(iface+1,iface-1,mod(iface,2)==1)
          valid=valid.and.metadata(16+opposite,neighbor_column)==metadata(1,i)
          axis=(iface+1)/2
          do component=1,3
            if(component==axis)cycle
            valid=valid.and.metadata(4+component,i)==metadata(4+component,neighbor_column).and. &
              metadata(7+component,i)==metadata(7+component,neighbor_column)
          enddo
          if(mod(iface,2)==1)then
            if(metadata(4+axis,i)==0)then
              valid=valid.and.metadata(4+axis,neighbor_column)+metadata(7+axis,neighbor_column)== &
                metadata(1+axis,i)
            else
              valid=valid.and.metadata(4+axis,neighbor_column)+metadata(7+axis,neighbor_column)== &
                metadata(4+axis,i)
            endif
          else
            if(metadata(4+axis,i)+metadata(7+axis,i)==metadata(1+axis,i))then
              valid=valid.and.metadata(4+axis,neighbor_column)==0
            else
              valid=valid.and.metadata(4+axis,i)+metadata(7+axis,i)==metadata(4+axis,neighbor_column)
            endif
          endif
        endif
      enddo
    enddo
    do i=1,nfragment
      valid=valid.and.count(metadata(1,:)==i)>0
    enddo
    valid=valid.and.fragment_volume==global_volume.and.density_volume==global_volume.and. &
      orbital_core_volume==global_volume
  end subroutine validate_direct_checkpoint_topology

  pure logical function integer_boxes_overlap(origin_a,size_a,origin_b,size_b)result(overlap)
    integer,intent(in)::origin_a(3),size_a(3),origin_b(3),size_b(3)
    overlap=all(origin_a<origin_b+size_b).and.all(origin_b<origin_a+size_a)
  end function integer_boxes_overlap

  pure integer(int64) function safe_shape_product(dimensions)result(value)
    integer,intent(in)::dimensions(:)
    integer::i
    value=1_int64
    do i=1,size(dimensions)
      if(dimensions(i)<=0.or.value>huge(value)/int(dimensions(i),int64))then
        value=huge(value)
        return
      endif
      value=value*int(dimensions(i),int64)
    enddo
  end function safe_shape_product

  integer(int64) function direct_checkpoint_metadata_hash(state)result(hash)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    hash=hash_words(int(z'243F6A8885A308D3',int64),[int(state%nstate,int64),int(state%nspin,int64), &
      int(state%scf_iterations,int64),int(state%global_size,int64),state%geometry_fingerprint, &
      state%operator_fingerprint,state%hamiltonian_operator_fingerprint,state%projector_fingerprint, &
      int(state%projector_generation,int64),int(state%continuation_rollbacks,int64), &
      transfer(state%interface_scale,0_int64),transfer(state%density_residual,0_int64), &
      transfer(state%charge_error,0_int64)])
    hash=hash_words(hash,transfer(state%continuation_scale_history,[0_int64], &
      size(state%continuation_scale_history)))
    hash=hash_words(hash,transfer(state%continuation_step_history,[0_int64], &
      size(state%continuation_step_history)))
    hash=hash_words(hash,[transfer(state%density_tolerance,0_int64),transfer(state%orbital_tolerance,0_int64), &
      transfer(state%orthogonality_tolerance,0_int64),transfer(state%hermiticity_tolerance,0_int64), &
      transfer(state%charge_tolerance,0_int64),transfer(state%face_balance_tolerance,0_int64)])
  end function direct_checkpoint_metadata_hash

  integer(int64) function direct_checkpoint_checksum(state)result(checksum)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    integer(int64),allocatable::words(:)
    checksum=direct_checkpoint_metadata_hash(state)
    checksum=hash_words(checksum,[int(state%fragment_id,int64),int(state%scf_iterations,int64), &
      int(state%state_start,int64),int(state%state_count,int64), &
      int(state%orbital_spatial_lower_bound,int64),int(state%full_spatial_shape,int64), &
      int(state%orbital_core_local_origin,int64),int(state%orbital_core_origin,int64), &
      int(state%orbital_core_size,int64), &
      int(state%fragment_origin,int64),int(state%face_neighbors,int64),int(state%continuation_rollbacks,int64), &
      int(state%density_origin,int64), &
      int(state%projector_retained_rank,int64),int(state%projector_required_rank,int64), &
      transfer(state%density_residual,0_int64),transfer(state%orbital_residual,0_int64), &
      transfer(state%orthogonality_defect,0_int64),transfer(state%hermiticity_defect,0_int64), &
      transfer(state%charge_error,0_int64),transfer(state%face_balance_defect,0_int64), &
      transfer(state%density_tolerance,0_int64),transfer(state%orbital_tolerance,0_int64), &
      transfer(state%orthogonality_tolerance,0_int64),transfer(state%hermiticity_tolerance,0_int64), &
      transfer(state%charge_tolerance,0_int64),transfer(state%face_balance_tolerance,0_int64), &
      transfer(state%projector_projection_residual,0_int64), &
      transfer(state%projector_escape_norm,0_int64),transfer(state%projector_residual_limit,0_int64), &
      transfer(state%projector_escape_limit,0_int64),transfer(state%face_action_norm,0_int64), &
      transfer(state%face_pair_balance,0_int64)])
    words=transfer(state%fragment_orbitals,[0_int64],size(state%fragment_orbitals))
    checksum=hash_words(checksum,words)
    words=transfer(state%occupations,[0_int64],size(state%occupations))
    checksum=hash_words(checksum,words)
    words=transfer(state%continuation_scale_history,[0_int64],size(state%continuation_scale_history))
    checksum=hash_words(checksum,words)
    words=transfer(state%continuation_step_history,[0_int64],size(state%continuation_step_history))
    checksum=hash_words(checksum,words)
    words=transfer(state%density,[0_int64],size(state%density));checksum=hash_words(checksum,words)
    words=transfer(state%hartree,[0_int64],size(state%hartree));checksum=hash_words(checksum,words)
    words=transfer(state%vxc,[0_int64],size(state%vxc));checksum=hash_words(checksum,words)
    words=transfer(state%vlocal,[0_int64],size(state%vlocal));checksum=hash_words(checksum,words)
    if(checksum==0_int64)checksum=1_int64
  end function direct_checkpoint_checksum

  integer(int64) function direct_checkpoint_field_hash(state)result(hash)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    integer(int64),allocatable::words(:)
    hash=int(z'13198A2E03707344',int64)
    words=transfer(state%density,[0_int64],size(state%density));hash=hash_words(hash,words)
    words=transfer(state%hartree,[0_int64],size(state%hartree));hash=hash_words(hash,words)
    words=transfer(state%vxc,[0_int64],size(state%vxc));hash=hash_words(hash,words)
    words=transfer(state%vlocal,[0_int64],size(state%vlocal));hash=hash_words(hash,words)
    if(hash==0_int64)hash=1_int64
  end function direct_checkpoint_field_hash

  integer(int64) function direct_checkpoint_occupation_hash(state)result(hash)
    type(s_dg_dc_direct_checkpoint_state),intent(in)::state
    integer(int64),allocatable::words(:)
    hash=int(z'A4093822299F31D0',int64)
    words=transfer(state%occupations,[0_int64],size(state%occupations))
    hash=hash_words(hash,words)
    if(hash==0_int64)hash=1_int64
  end function direct_checkpoint_occupation_hash
end module dg_ground_state_checkpoint
