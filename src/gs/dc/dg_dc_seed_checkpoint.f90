#include "config.h"
module dg_dc_seed_checkpoint
  use iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  integer,parameter,public::DG_DC_SEED_ABSENT=0
  integer,parameter,public::DG_DC_SEED_VALID=1
  integer,parameter,public::DG_DC_SEED_INVALID=2
  integer,parameter::seed_version=1
  integer(int64),parameter::max_seed_elements=1099511627776_int64
  integer,parameter::max_seed_directory_length=400
  integer,parameter::stream_chunk_elements=65536
  character(32),parameter::manifest_magic='SALMON_DG_DC_SEED_MANIFEST_V1'
  character(32),parameter::shard_magic='SALMON_DG_DC_SEED_SHARD_V1'
  character(32),parameter::present_magic='SALMON_DG_DC_SEED_PRESENT_V1'
  character(32),parameter::pending_magic='SALMON_DG_DC_SEED_PENDING_V1'

  type,public::s_dg_dc_seed_contract
    integer::version=1,mpi_size=0,rank=0,fragment_id=0
    integer::rwf_bounds(14)=0,rho_bounds(6)=0,vloc_bounds(6)=0
    integer(int64)::immutable_fingerprint=0_int64,ownership_fingerprint=0_int64
  end type

  type,public::s_dg_dc_seed_payload
    real(8),allocatable::rwf(:,:,:,:,:,:,:)
    real(8),allocatable::rho_tot(:,:,:),vloc_tot(:,:,:)
    real(8),allocatable::esp(:,:,:),rocc(:,:,:)
    real(8)::mu=0d0,residual=huge(0d0)
    integer::iteration=0
  end type

  type::s_dg_dc_seed_manifest
    integer::version=0,mpi_size=0
    integer(int64)::publication_id=0_int64,ordered_digest=0_int64
    integer(int64)::manifest_size=0_int64,manifest_digest=0_int64
    real(8)::density_weight=0d0,expected_electrons=0d0
    real(8)::write_electron_tolerance=0d0,write_threshold=0d0
    integer,allocatable::fragment_ids(:),rwf_bounds(:,:),rho_bounds(:,:),vloc_bounds(:,:)
    integer,allocatable::esp_bounds(:,:),rocc_bounds(:,:)
    integer(int64),allocatable::immutable_fingerprints(:),ownership_fingerprints(:)
    integer(int64),allocatable::shard_sizes(:),shard_digests(:)
  end type

  public::write_dg_dc_seed,read_dg_dc_seed,probe_dg_dc_seed
  public::resolve_dg_dc_seed_mode
  public::build_dg_dc_seed_contract
  public::restore_dg_dc_seed_payload

contains

  subroutine resolve_dg_dc_seed_mode(mode,status,run_scf,load_seed,publish_seed,&
      fatal,scf_skipped,ok,message)
    character(*),intent(in)::mode
    integer,intent(in)::status
    logical,intent(out)::run_scf,load_seed,publish_seed,fatal,scf_skipped,ok
    character(*),intent(out)::message
    character(len(mode))::normalized_mode

    normalized_mode=lower_ascii(trim(adjustl(mode)))
    run_scf=.false.;load_seed=.false.;publish_seed=.false.;fatal=.false.
    scf_skipped=.false.;ok=.true.;message=''
    select case(trim(normalized_mode))
    case('off')
      run_scf=.true.
    case('write')
      run_scf=.true.;publish_seed=.true.
    case('read')
      select case(status)
      case(DG_DC_SEED_VALID)
        load_seed=.true.;scf_skipped=.true.
      case(DG_DC_SEED_ABSENT)
        fatal=.true.;message='requested DG DC seed is absent'
      case default
        fatal=.true.;message='requested DG DC seed is invalid'
      end select
    case('auto')
      select case(status)
      case(DG_DC_SEED_VALID)
        load_seed=.true.;scf_skipped=.true.
      case(DG_DC_SEED_ABSENT)
        run_scf=.true.;publish_seed=.true.
      case default
        fatal=.true.;message='automatic DG DC seed is invalid'
      end select
    case default
      fatal=.true.;message='unknown DG DC seed mode'
    end select
    ok=.not.fatal
  end subroutine resolve_dg_dc_seed_mode

  subroutine build_dg_dc_seed_contract(comm,fragment_id,rwf_bounds,rho_bounds,&
      vloc_bounds,immutable_inputs,ownership_map,contract,ok,message)
    integer,intent(in)::comm,fragment_id
    integer,intent(in)::rwf_bounds(14),rho_bounds(6),vloc_bounds(6)
    integer(int64),intent(in)::immutable_inputs(:),ownership_map(:)
    type(s_dg_dc_seed_contract),intent(out)::contract
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_bad,global_bad,i
    integer(int64)::local_immutable,local_ownership
    integer(int64),allocatable::immutable_by_rank(:),ownership_by_rank(:)

    contract=s_dg_dc_seed_contract();ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(fragment_id<1.or.size(immutable_inputs)<1.or.size(ownership_map)<1)local_bad=1
    if(.not.valid_bounds(rwf_bounds,7).or..not.valid_bounds(rho_bounds,3).or.&
       .not.valid_bounds(vloc_bounds,3))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid distributed DG DC seed contract inputs'
      return
    endif

    local_immutable=int(z'510E527FADE682D1',int64)
    call hash_integer8_vector(local_immutable,immutable_inputs)
    local_ownership=int(z'9B05688C2B3E6C1F',int64)
    call hash_integer(local_ownership,rank)
    call hash_integer(local_ownership,fragment_id)
    call hash_integer_vector(local_ownership,rwf_bounds)
    call hash_integer_vector(local_ownership,rho_bounds)
    call hash_integer_vector(local_ownership,vloc_bounds)
    call hash_integer8_vector(local_ownership,ownership_map)
    local_bad=0
    allocate(immutable_by_rank(nproc),ownership_by_rank(nproc),stat=local_bad)
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot allocate distributed DG DC seed fingerprints'
      if(allocated(immutable_by_rank))deallocate(immutable_by_rank)
      if(allocated(ownership_by_rank))deallocate(ownership_by_rank)
      return
    endif
    call MPI_Allgather(local_immutable,1,MPI_INTEGER8,immutable_by_rank,1,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(local_ownership,1,MPI_INTEGER8,ownership_by_rank,1,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot assemble distributed DG DC seed fingerprints'
      deallocate(immutable_by_rank,ownership_by_rank)
      return
    endif

    contract%version=seed_version;contract%mpi_size=nproc;contract%rank=rank
    contract%fragment_id=fragment_id
    contract%rwf_bounds=rwf_bounds;contract%rho_bounds=rho_bounds
    contract%vloc_bounds=vloc_bounds
    contract%immutable_fingerprint=int(z'1F83D9ABFB41BD6B',int64)
    contract%ownership_fingerprint=int(z'5BE0CD19137E2179',int64)
    call hash_integer(contract%immutable_fingerprint,nproc)
    call hash_integer(contract%ownership_fingerprint,nproc)
    do i=1,nproc
      call hash_integer(contract%immutable_fingerprint,i-1)
      call hash_integer8(contract%immutable_fingerprint,immutable_by_rank(i))
      call hash_integer(contract%ownership_fingerprint,i-1)
      call hash_integer8(contract%ownership_fingerprint,ownership_by_rank(i))
    enddo
    deallocate(immutable_by_rank,ownership_by_rank)
    ok=.true.
#else
    contract=s_dg_dc_seed_contract();ok=.false.
    message='DG DC seed contracts require MPI'
#endif
  end subroutine build_dg_dc_seed_contract

  subroutine restore_dg_dc_seed_payload(payload,rwf,rho_owned,vloc_owned,esp,rocc,&
      mu,residual,iteration,ok,message)
    type(s_dg_dc_seed_payload),intent(in)::payload
    real(8),allocatable,intent(inout)::rwf(:,:,:,:,:,:,:)
    real(8),allocatable,intent(inout)::rho_owned(:,:,:),vloc_owned(:,:,:)
    real(8),allocatable,intent(inout)::esp(:,:,:),rocc(:,:,:)
    real(8),intent(out)::mu,residual
    integer,intent(out)::iteration
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::allocation_status

    ok=.false.;message='';mu=0d0;residual=huge(0d0);iteration=0
    if(.not.allocated(payload%rwf).or..not.allocated(payload%rho_tot).or.&
       .not.allocated(payload%vloc_tot).or..not.allocated(payload%esp).or.&
       .not.allocated(payload%rocc))then
      message='incomplete DG DC seed payload'
      return
    endif
    if(.not.all(ieee_is_finite(payload%rwf)).or.&
       .not.all(ieee_is_finite(payload%rho_tot)).or.&
       .not.all(ieee_is_finite(payload%vloc_tot)).or.&
       .not.all(ieee_is_finite(payload%esp)).or.&
       .not.all(ieee_is_finite(payload%rocc)).or.&
       .not.ieee_is_finite(payload%mu).or..not.ieee_is_finite(payload%residual))then
      message='non-finite DG DC seed payload'
      return
    endif
    if(allocated(rwf))deallocate(rwf)
    if(allocated(rho_owned))deallocate(rho_owned)
    if(allocated(vloc_owned))deallocate(vloc_owned)
    if(allocated(esp))deallocate(esp)
    if(allocated(rocc))deallocate(rocc)
    allocation_status=0
    allocate(rwf(lbound(payload%rwf,1):ubound(payload%rwf,1),&
      lbound(payload%rwf,2):ubound(payload%rwf,2),lbound(payload%rwf,3):ubound(payload%rwf,3),&
      lbound(payload%rwf,4):ubound(payload%rwf,4),lbound(payload%rwf,5):ubound(payload%rwf,5),&
      lbound(payload%rwf,6):ubound(payload%rwf,6),lbound(payload%rwf,7):ubound(payload%rwf,7)),&
      rho_owned(lbound(payload%rho_tot,1):ubound(payload%rho_tot,1),&
      lbound(payload%rho_tot,2):ubound(payload%rho_tot,2),&
      lbound(payload%rho_tot,3):ubound(payload%rho_tot,3)),&
      vloc_owned(lbound(payload%vloc_tot,1):ubound(payload%vloc_tot,1),&
      lbound(payload%vloc_tot,2):ubound(payload%vloc_tot,2),&
      lbound(payload%vloc_tot,3):ubound(payload%vloc_tot,3)),&
      esp(lbound(payload%esp,1):ubound(payload%esp,1),lbound(payload%esp,2):ubound(payload%esp,2),&
      lbound(payload%esp,3):ubound(payload%esp,3)),&
      rocc(lbound(payload%rocc,1):ubound(payload%rocc,1),&
      lbound(payload%rocc,2):ubound(payload%rocc,2),&
      lbound(payload%rocc,3):ubound(payload%rocc,3)),stat=allocation_status)
    if(allocation_status/=0)then
      message='cannot allocate restored DG DC seed payload'
      return
    endif
    rwf=payload%rwf;rho_owned=payload%rho_tot;vloc_owned=payload%vloc_tot
    esp=payload%esp;rocc=payload%rocc
    mu=payload%mu;residual=payload%residual;iteration=payload%iteration
    ok=.true.
  end subroutine restore_dg_dc_seed_payload

  subroutine write_dg_dc_seed(comm,directory,contract,payload,density_weight,&
      expected_electrons,electron_tolerance,current_threshold,publication_id,ok,message,&
      failure_injection_rank)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_dc_seed_contract),intent(in)::contract
    type(s_dg_dc_seed_payload),intent(in)::payload
    real(8),intent(in)::density_weight,expected_electrons,electron_tolerance,current_threshold
    integer(int64),intent(out)::publication_id
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::failure_injection_rank
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_bad,global_bad,ios,injected
    integer::esp_bounds(6),rocc_bounds(6)
    integer(int64)::shard_size,shard_digest,digest_min,digest_max
    real(8)::local_density_sum,global_electrons,file_mu,file_residual
    integer::file_iteration
    character(512)::present_path,pending_path,manifest_path,manifest_temporary
    character(512)::shard,shard_temporary,reservation
    type(s_dg_dc_seed_manifest)::manifest
    logical::file_ok,present_exists

    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,nproc,ierr)
    publication_id=0_int64;ok=.false.;message=''
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call validate_current_contract(contract,nproc,rank,local_bad)
    call validate_directory(comm,directory,local_bad)
    call validate_parameters(comm,density_weight,expected_electrons,electron_tolerance,&
      current_threshold,local_bad)
    call validate_payload_local(contract,payload,current_threshold,local_bad,local_density_sum)
    call validate_collective_provenance(comm,payload%mu,payload%residual,&
      payload%iteration,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_density_sum,global_electrons,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    global_electrons=global_electrons*density_weight
    if(.not.ieee_is_finite(global_electrons).or.&
       abs(global_electrons-expected_electrons)>electron_tolerance)global_bad=1
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid or unconverged DG DC seed payload'
      return
    endif

    esp_bounds=[lbound(payload%esp),ubound(payload%esp)]
    rocc_bounds=[lbound(payload%rocc),ubound(payload%rocc)]
    present_path=seed_path(directory,'dg_dc_seed.present')
    pending_path=seed_path(directory,'dg_dc_seed.pending')
    manifest_path=seed_path(directory,'dg_dc_seed.manifest')
    ios=0
    if(rank==0)then
      inquire(file=trim(present_path),exist=present_exists)
      if(.not.present_exists)call write_marker(present_path,present_magic,0_int64,ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot establish DG DC seed publication marker'
      return
    endif
    if(rank==0)then
      call choose_publication_id(directory,contract%immutable_fingerprint,publication_id,reservation,ios)
      if(ios==0)call write_marker(pending_path,pending_magic,publication_id,ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(publication_id,1,MPI_INTEGER8,0,comm,ierr)
    if(ios/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot begin DG DC seed publication'
      return
    endif
    call MPI_Barrier(comm,ierr)

    shard=shard_name(directory,publication_id,rank)
    shard_temporary=trim(shard)//'.temporary'
    shard_digest=digest_shard(contract,payload,publication_id,esp_bounds,rocc_bounds)
    call write_shard_file(shard_temporary,contract,payload,publication_id,esp_bounds,&
      rocc_bounds,shard_digest,shard_size,file_ok)
    local_density_sum=0d0
    file_mu=payload%mu;file_residual=payload%residual;file_iteration=payload%iteration
    if(file_ok)call validate_shard_file_stream(shard_temporary,contract,publication_id,&
      esp_bounds,rocc_bounds,shard_size,shard_digest,current_threshold,&
      local_density_sum,file_mu,file_residual,file_iteration,file_ok)
    local_bad=merge(0,1,file_ok)
    call validate_collective_provenance(comm,file_mu,file_residual,file_iteration,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_density_sum,global_electrons,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    global_electrons=global_electrons*density_weight
    if(.not.ieee_is_finite(global_electrons).or.&
       abs(global_electrons-expected_electrons)>electron_tolerance)global_bad=1
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot write and validate DG DC seed rank shard'
      return
    endif
    call rename(trim(shard_temporary),trim(shard),ios)
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot atomically publish DG DC seed rank shard'
      return
    endif

    injected=0
    if(present(failure_injection_rank))then
      if(rank==failure_injection_rank)injected=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,injected,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(injected/=0.or.ierr/=MPI_SUCCESS)then
      message='injected interruption before DG DC seed manifest publication'
      return
    endif

    call allocate_manifest(manifest,nproc)
    manifest%version=seed_version;manifest%mpi_size=nproc
    manifest%publication_id=publication_id
    manifest%density_weight=density_weight
    manifest%expected_electrons=expected_electrons
    manifest%write_electron_tolerance=electron_tolerance
    manifest%write_threshold=current_threshold
    call MPI_Allgather(contract%fragment_id,1,MPI_INTEGER,manifest%fragment_ids,1,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgather(contract%rwf_bounds,14,MPI_INTEGER,manifest%rwf_bounds,14,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgather(contract%rho_bounds,6,MPI_INTEGER,manifest%rho_bounds,6,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgather(contract%vloc_bounds,6,MPI_INTEGER,manifest%vloc_bounds,6,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgather(esp_bounds,6,MPI_INTEGER,manifest%esp_bounds,6,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(rocc_bounds,6,MPI_INTEGER,manifest%rocc_bounds,6,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(contract%immutable_fingerprint,1,MPI_INTEGER8,&
      manifest%immutable_fingerprints,1,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(contract%ownership_fingerprint,1,MPI_INTEGER8,&
      manifest%ownership_fingerprints,1,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(shard_size,1,MPI_INTEGER8,manifest%shard_sizes,1,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(shard_digest,1,MPI_INTEGER8,manifest%shard_digests,1,MPI_INTEGER8,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    manifest%ordered_digest=digest_ordered_manifest(manifest)
    call MPI_Allreduce(manifest%ordered_digest,digest_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(manifest%ordered_digest,digest_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.digest_min/=digest_max)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='rank-inconsistent DG DC seed manifest metadata'
      return
    endif

    write(manifest_temporary,'(a,".temporary.",z16.16)')trim(manifest_path),publication_id
    ios=0
    if(rank==0)then
      call write_manifest_file(manifest_temporary,manifest,file_ok)
      if(.not.file_ok)ios=1
      if(ios==0)call rename(trim(manifest_temporary),trim(manifest_path),ios)
      if(ios==0)call remove_file_if_present(pending_path)
      if(ios==0)call remove_file_if_present(reservation)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0.or.ierr/=MPI_SUCCESS)then
      message='cannot atomically publish DG DC seed manifest'
      return
    endif
    call MPI_Barrier(comm,ierr)
    ok=ierr==MPI_SUCCESS
    if(.not.ok)message='DG DC seed publication barrier failed'
    call clear_manifest(manifest)
#else
    publication_id=0_int64;ok=.false.;message='DG DC seed checkpoints require MPI'
#endif
  end subroutine write_dg_dc_seed

  subroutine probe_dg_dc_seed(comm,directory,contract,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,status,publication_id,message)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_dc_seed_contract),intent(in)::contract
    real(8),intent(in)::density_weight,expected_electrons,electron_tolerance,current_threshold
    integer,intent(out)::status
    integer(int64),intent(out)::publication_id
    character(*),intent(out)::message
    call inspect_dg_dc_seed(comm,directory,contract,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,status,publication_id,message)
  end subroutine probe_dg_dc_seed

  subroutine read_dg_dc_seed(comm,directory,contract,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,payload,publication_id,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_dc_seed_contract),intent(in)::contract
    real(8),intent(in)::density_weight,expected_electrons,electron_tolerance,current_threshold
    type(s_dg_dc_seed_payload),intent(out)::payload
    integer(int64),intent(out)::publication_id
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::status
    call inspect_dg_dc_seed(comm,directory,contract,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,status,publication_id,message,payload)
    ok=status==DG_DC_SEED_VALID
    if(.not.ok)call clear_payload(payload)
  end subroutine read_dg_dc_seed

  subroutine inspect_dg_dc_seed(comm,directory,contract,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,status,publication_id,message,payload)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_dc_seed_contract),intent(in)::contract
    real(8),intent(in)::density_weight,expected_electrons,electron_tolerance,current_threshold
    integer,intent(out)::status
    integer(int64),intent(out)::publication_id
    character(*),intent(out)::message
    type(s_dg_dc_seed_payload),intent(out),optional::payload
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_bad,global_bad,artifact_state,file_iteration
    integer::local_contract_mismatch,global_contract_mismatch,stored_nproc
    real(8)::local_density_sum,global_electrons,file_mu,file_residual
    type(s_dg_dc_seed_manifest)::manifest
    logical::file_ok
    character(512)::shard

    if(present(payload))call clear_payload(payload)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    status=DG_DC_SEED_INVALID;publication_id=0_int64;message=''
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call validate_current_contract(contract,nproc,rank,local_bad)
    call validate_directory(comm,directory,local_bad)
    call validate_parameters(comm,density_weight,expected_electrons,electron_tolerance,&
      current_threshold,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid current DG DC seed contract or policy'
      return
    endif

    call classify_seed_artifacts(comm,directory,artifact_state)
    if(artifact_state==DG_DC_SEED_ABSENT)then
      status=DG_DC_SEED_ABSENT;message='no DG DC seed artifacts'
      return
    else if(artifact_state==DG_DC_SEED_INVALID)then
      message='status=invalid cause=publication_state detail=incomplete_or_interrupted'
      return
    endif

    call read_manifest_rank_count_collective(comm,&
      seed_path(directory,'dg_dc_seed.manifest'),stored_nproc,file_ok)
    if(.not.file_ok)then
      message='status=invalid cause=manifest_integrity detail=unreadable_header'
      return
    else if(stored_nproc/=nproc)then
      message='status=invalid cause=mpi_rank_count detail=manifest_runtime_mismatch'
      return
    endif

    call read_manifest_collective(comm,seed_path(directory,'dg_dc_seed.manifest'),nproc,manifest,file_ok)
    local_contract_mismatch=0
    if(file_ok)then
      if(manifest%version/=seed_version.or.manifest%mpi_size/=nproc.or.&
         manifest%publication_id==0_int64)local_contract_mismatch=max(local_contract_mismatch,1)
      if(manifest%density_weight/=density_weight.or.&
         manifest%expected_electrons/=expected_electrons)&
        local_contract_mismatch=max(local_contract_mismatch,2)
      if(manifest%immutable_fingerprints(rank+1)/=contract%immutable_fingerprint)&
        local_contract_mismatch=max(local_contract_mismatch,3)
      if(any(manifest%rwf_bounds(:,rank+1)/=contract%rwf_bounds).or.&
         any(manifest%rho_bounds(:,rank+1)/=contract%rho_bounds).or.&
         any(manifest%vloc_bounds(:,rank+1)/=contract%vloc_bounds))&
        local_contract_mismatch=max(local_contract_mismatch,4)
      if(manifest%ownership_fingerprints(rank+1)/=contract%ownership_fingerprint)&
        local_contract_mismatch=max(local_contract_mismatch,5)
      if(manifest%fragment_ids(rank+1)/=contract%fragment_id)&
        local_contract_mismatch=max(local_contract_mismatch,6)
    else
      local_contract_mismatch=1
    endif
    call MPI_Allreduce(local_contract_mismatch,global_contract_mismatch,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_contract_mismatch/=0.or.ierr/=MPI_SUCCESS)then
      select case(global_contract_mismatch)
      case(2)
        message='status=invalid cause=seed_policy detail=density_or_electron_mismatch'
      case(3)
        message='status=invalid cause=immutable_inputs detail=cell_atom_fragment_physics_operator_or_convergence_mismatch'
      case(4)
        message='status=invalid cause=local_bounds detail=rank_local_array_mismatch'
      case(5)
        message='status=invalid cause=ownership_mapping detail=rank_owned_grid_mismatch'
      case(6)
        message='status=invalid cause=rank_fragment_mapping detail=rank_fragment_mismatch'
      case default
        message='status=invalid cause=manifest_integrity detail=corrupt_or_incompatible'
      end select
      call clear_manifest(manifest)
      return
    endif

    publication_id=manifest%publication_id
    shard=shard_name(directory,publication_id,rank)
    local_density_sum=0d0
    file_mu=0d0;file_residual=huge(0d0);file_iteration=-1
    if(present(payload))then
      call read_shard_file(shard,contract,publication_id,manifest%esp_bounds(:,rank+1),&
        manifest%rocc_bounds(:,rank+1),manifest%shard_sizes(rank+1),&
        manifest%shard_digests(rank+1),payload,file_ok)
      local_bad=merge(0,1,file_ok)
      if(file_ok)call validate_payload_local(contract,payload,current_threshold,&
        local_bad,local_density_sum)
      if(file_ok)then
        file_mu=payload%mu;file_residual=payload%residual;file_iteration=payload%iteration
      endif
    else
      call validate_shard_file_stream(shard,contract,publication_id,&
        manifest%esp_bounds(:,rank+1),manifest%rocc_bounds(:,rank+1),&
        manifest%shard_sizes(rank+1),manifest%shard_digests(rank+1),current_threshold,&
        local_density_sum,file_mu,file_residual,file_iteration,file_ok)
      local_bad=merge(0,1,file_ok)
    endif
    call validate_collective_provenance(comm,file_mu,file_residual,file_iteration,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_density_sum,global_electrons,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    global_electrons=global_electrons*density_weight
    if(.not.ieee_is_finite(global_electrons).or.&
       abs(global_electrons-expected_electrons)>electron_tolerance)global_bad=1
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='status=invalid cause=shard_payload detail=missing_corrupt_unconverged_or_wrong_electrons'
      if(present(payload))call clear_payload(payload)
      call clear_manifest(manifest)
      status=DG_DC_SEED_INVALID
      return
    endif
    status=DG_DC_SEED_VALID;message='valid committed DG DC seed'
    call clear_manifest(manifest)
#else
    if(present(payload))call clear_payload(payload)
    status=DG_DC_SEED_INVALID;publication_id=0_int64
    message='DG DC seed checkpoints require MPI'
#endif
  end subroutine inspect_dg_dc_seed

#ifdef USE_MPI
  subroutine validate_parameters(comm,density_weight,expected_electrons,electron_tolerance,&
      current_threshold,local_bad)
    integer,intent(in)::comm
    real(8),intent(in)::density_weight,expected_electrons,electron_tolerance,current_threshold
    integer,intent(inout)::local_bad
    real(8)::values(4),minimum(4),maximum(4)
    integer::ierr
    values=[density_weight,expected_electrons,electron_tolerance,current_threshold]
    if(.not.all(ieee_is_finite(values)).or.density_weight<=0d0.or.&
       expected_electrons<0d0.or.electron_tolerance<0d0.or.current_threshold<=0d0)local_bad=1
    call MPI_Allreduce(values,minimum,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(values,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum))local_bad=1
  end subroutine validate_parameters

  subroutine validate_collective_provenance(comm,mu,residual,iteration,local_bad)
    integer,intent(in)::comm,iteration
    real(8),intent(in)::mu,residual
    integer,intent(inout)::local_bad
    real(8)::values(2),minimum(2),maximum(2)
    integer::iteration_minimum,iteration_maximum,ierr

    values=[mu,residual]
    call MPI_Allreduce(values,minimum,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
      return
    endif
    call MPI_Allreduce(values,maximum,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
      return
    endif
    call MPI_Allreduce(iteration,iteration_minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
      return
    endif
    call MPI_Allreduce(iteration,iteration_maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum).or.&
       iteration_minimum/=iteration_maximum)local_bad=1
  end subroutine validate_collective_provenance

  subroutine validate_directory(comm,directory,local_bad)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    integer,intent(inout)::local_bad
    character(512)::root_directory
    integer::rank,ierr
    root_directory=''
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(rank==0.and.len_trim(directory)<=max_seed_directory_length)&
      root_directory=trim(directory)
    call MPI_Bcast(root_directory,len(root_directory),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(len_trim(directory)<1.or.len_trim(directory)>max_seed_directory_length.or.&
       trim(directory)/=trim(root_directory))local_bad=1
  end subroutine validate_directory

  subroutine validate_current_contract(contract,nproc,rank,local_bad)
    type(s_dg_dc_seed_contract),intent(in)::contract
    integer,intent(in)::nproc,rank
    integer,intent(inout)::local_bad
    if(contract%version/=seed_version.or.contract%mpi_size/=nproc.or.contract%rank/=rank.or.&
       contract%fragment_id<1.or.contract%immutable_fingerprint==0_int64.or.&
       contract%ownership_fingerprint==0_int64)local_bad=1
    if(.not.valid_bounds(contract%rwf_bounds,7).or.&
       .not.valid_bounds(contract%rho_bounds,3).or.&
       .not.valid_bounds(contract%vloc_bounds,3))local_bad=1
  end subroutine validate_current_contract

  subroutine validate_payload_local(contract,payload,current_threshold,local_bad,local_density_sum)
    type(s_dg_dc_seed_contract),intent(in)::contract
    type(s_dg_dc_seed_payload),intent(in)::payload
    real(8),intent(in)::current_threshold
    integer,intent(inout)::local_bad
    real(8),intent(out)::local_density_sum
    local_density_sum=0d0
    if(.not.allocated(payload%rwf).or..not.allocated(payload%rho_tot).or.&
       .not.allocated(payload%vloc_tot).or..not.allocated(payload%esp).or.&
       .not.allocated(payload%rocc))then
      local_bad=1;return
    endif
    if(any([lbound(payload%rwf),ubound(payload%rwf)]/=contract%rwf_bounds).or.&
       any([lbound(payload%rho_tot),ubound(payload%rho_tot)]/=contract%rho_bounds).or.&
       any([lbound(payload%vloc_tot),ubound(payload%vloc_tot)]/=contract%vloc_bounds))local_bad=1
    if(any(lbound(payload%esp)/=lbound(payload%rocc)).or.&
       any(ubound(payload%esp)/=ubound(payload%rocc)).or.size(payload%esp)<1)local_bad=1
    if(.not.all(ieee_is_finite(payload%rwf)).or.&
       .not.all(ieee_is_finite(payload%rho_tot)).or.&
       .not.all(ieee_is_finite(payload%vloc_tot)).or.&
       .not.all(ieee_is_finite(payload%esp)).or.&
       .not.all(ieee_is_finite(payload%rocc)).or..not.ieee_is_finite(payload%mu).or.&
       .not.ieee_is_finite(payload%residual))local_bad=1
    if(payload%residual<0d0.or..not.(payload%residual<current_threshold).or.&
       payload%iteration<0) local_bad=1
    local_density_sum=sum(payload%rho_tot)
    if(.not.ieee_is_finite(local_density_sum))local_bad=1
  end subroutine validate_payload_local

  subroutine classify_seed_artifacts(comm,directory,state)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    integer,intent(out)::state
    integer::rank,ierr,flags(3)
    logical::exists
    ! The persistent present marker is created before reservations, temporary
    ! shards, or manifests, so every artifact this writer can orphan is indexed
    ! by one of these fixed names.  Once a manifest exists it alone is the commit
    ! point; a stale pending marker must not hide the previous committed version.
    call MPI_Comm_rank(comm,rank,ierr);flags=0
    if(rank==0)then
      inquire(file=trim(seed_path(directory,'dg_dc_seed.present')),exist=exists)
      flags(1)=merge(1,0,exists)
      inquire(file=trim(seed_path(directory,'dg_dc_seed.pending')),exist=exists)
      flags(2)=merge(1,0,exists)
      inquire(file=trim(seed_path(directory,'dg_dc_seed.manifest')),exist=exists)
      flags(3)=merge(1,0,exists)
    endif
    call MPI_Bcast(flags,3,MPI_INTEGER,0,comm,ierr)
    if(flags(3)==1)then
      state=DG_DC_SEED_VALID
    else if(flags(1)==0.and.flags(2)==0)then
      state=DG_DC_SEED_ABSENT
    else
      state=DG_DC_SEED_INVALID
    endif
  end subroutine classify_seed_artifacts

  subroutine choose_publication_id(directory,immutable_fingerprint,publication_id,reservation,ios)
    character(*),intent(in)::directory
    integer(int64),intent(in)::immutable_fingerprint
    integer(int64),intent(out)::publication_id
    character(*),intent(out)::reservation
    integer,intent(out)::ios
    integer::unit,attempt,close_ios,inquire_ios
    integer(int64)::clock
    logical::shard_exists
    ios=1;publication_id=0_int64;reservation=''
    call system_clock(clock)
    do attempt=0,1023
      publication_id=ieor(ishftc(clock,mod(attempt,63)),immutable_fingerprint)
      publication_id=ieor(publication_id,int(attempt+1,int64))
      if(publication_id==0_int64)publication_id=int(attempt+1,int64)
      write(reservation,'(a,"/dg_dc_seed.",z16.16,".reservation")')trim(directory),publication_id
      open(newunit=unit,file=trim(reservation),status='new',access='stream',&
        form='unformatted',action='write',iostat=ios)
      if(ios==0)then
        write(unit,iostat=ios)present_magic,publication_id
        close_ios=0;close(unit,iostat=close_ios)
        if(ios==0.and.close_ios/=0)ios=close_ios
        if(ios==0)then
          inquire(file=trim(shard_name(directory,publication_id,0)),exist=shard_exists,&
            iostat=inquire_ios)
          if(inquire_ios==0.and..not.shard_exists)return
          call remove_file_if_present(reservation)
          ios=1
        endif
      endif
      clock=ieor(ishftc(clock,7),int(attempt+1,int64))
    enddo
  end subroutine choose_publication_id

  subroutine write_marker(filename,magic,publication_id,ios)
    character(*),intent(in)::filename,magic
    integer(int64),intent(in)::publication_id
    integer,intent(out)::ios
    integer::unit,flush_ios,close_ios
    open(newunit=unit,file=trim(filename),status='replace',access='stream',&
      form='unformatted',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)magic,publication_id
    flush_ios=0;close_ios=0
    if(ios==0)flush(unit,iostat=flush_ios)
    close(unit,iostat=close_ios)
    if(ios==0.and.flush_ios/=0)ios=flush_ios
    if(ios==0.and.close_ios/=0)ios=close_ios
  end subroutine write_marker

  subroutine write_shard_file(filename,contract,payload,publication_id,esp_bounds,&
      rocc_bounds,shard_digest,shard_size,ok)
    character(*),intent(in)::filename
    type(s_dg_dc_seed_contract),intent(in)::contract
    type(s_dg_dc_seed_payload),intent(in)::payload
    integer(int64),intent(in)::publication_id,shard_digest
    integer,intent(in)::esp_bounds(6),rocc_bounds(6)
    integer(int64),intent(out)::shard_size
    logical,intent(out)::ok
    integer::unit,ios,flush_ios,close_ios
    shard_size=0_int64;ok=.false.
    open(newunit=unit,file=trim(filename),status='replace',access='stream',&
      form='unformatted',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)shard_magic,contract%version,contract%mpi_size,&
      contract%rank,contract%fragment_id,publication_id,contract%immutable_fingerprint,&
      contract%ownership_fingerprint,shard_digest,contract%rwf_bounds,contract%rho_bounds,&
      contract%vloc_bounds,esp_bounds,rocc_bounds,payload%mu,payload%residual,payload%iteration
    if(ios==0)write(unit,iostat=ios)payload%rwf,payload%rho_tot,payload%vloc_tot,&
      payload%esp,payload%rocc
    flush_ios=0;close_ios=0
    if(ios==0)flush(unit,iostat=flush_ios)
    close(unit,iostat=close_ios)
    if(ios==0.and.flush_ios/=0)ios=flush_ios
    if(ios==0.and.close_ios/=0)ios=close_ios
    if(ios==0)inquire(file=trim(filename),size=shard_size,iostat=ios)
    ok=ios==0.and.shard_size>0_int64
  end subroutine write_shard_file

  subroutine validate_shard_file_stream(filename,contract,publication_id,&
      expected_esp_bounds,expected_rocc_bounds,expected_size,expected_digest,&
      current_threshold,local_density_sum,file_mu,file_residual,file_iteration,ok)
    character(*),intent(in)::filename
    type(s_dg_dc_seed_contract),intent(in)::contract
    integer(int64),intent(in)::publication_id,expected_size,expected_digest
    integer,intent(in)::expected_esp_bounds(6),expected_rocc_bounds(6)
    real(8),intent(in)::current_threshold
    real(8),intent(out)::local_density_sum,file_mu,file_residual
    integer,intent(out)::file_iteration
    logical,intent(out)::ok
    character(32)::magic
    integer::unit,ios,close_ios,file_version,file_nproc,file_rank,file_fragment,iteration
    integer::rwf_bounds(14),rho_bounds(6),vloc_bounds(6),esp_bounds(6),rocc_bounds(6)
    integer(int64)::file_publication,immutable_fingerprint,ownership_fingerprint
    integer(int64)::stored_digest,file_size,position,computed_digest
    real(8)::mu,residual,discarded_sum

    ok=.false.;local_density_sum=0d0;discarded_sum=0d0
    file_mu=0d0;file_residual=huge(0d0);file_iteration=-1
    open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,file_version,file_nproc,file_rank,file_fragment,file_publication,&
      immutable_fingerprint,ownership_fingerprint,stored_digest,rwf_bounds,rho_bounds,&
      vloc_bounds,esp_bounds,rocc_bounds,mu,residual,iteration
    if(ios==0)then
      file_mu=mu;file_residual=residual;file_iteration=iteration
      if(magic/=shard_magic.or.file_version/=contract%version.or.&
         file_nproc/=contract%mpi_size.or.file_rank/=contract%rank.or.&
         file_fragment/=contract%fragment_id.or.file_publication/=publication_id.or.&
         immutable_fingerprint/=contract%immutable_fingerprint.or.&
         ownership_fingerprint/=contract%ownership_fingerprint.or.&
         stored_digest/=expected_digest.or.any(rwf_bounds/=contract%rwf_bounds).or.&
         any(rho_bounds/=contract%rho_bounds).or.any(vloc_bounds/=contract%vloc_bounds).or.&
         any(esp_bounds/=expected_esp_bounds).or.any(rocc_bounds/=expected_rocc_bounds))ios=1
    endif
    if(ios==0)then
      if(.not.valid_bounds(esp_bounds,3).or..not.valid_bounds(rocc_bounds,3).or.&
         any(esp_bounds/=rocc_bounds))ios=1
      if(.not.ieee_is_finite(mu).or..not.ieee_is_finite(residual).or.residual<0d0.or.&
         .not.(residual<current_threshold).or.iteration<0)ios=1
    endif
    if(ios==0)then
      call initialize_shard_digest(computed_digest,contract,publication_id,esp_bounds,&
        rocc_bounds,mu,residual,iteration)
      call hash_real_stream(unit,computed_digest,rwf_bounds,7,.false.,discarded_sum,ios)
      if(ios==0)call hash_real_stream(unit,computed_digest,rho_bounds,3,.true.,&
        local_density_sum,ios)
      if(ios==0)call hash_real_stream(unit,computed_digest,vloc_bounds,3,.false.,discarded_sum,ios)
      if(ios==0)call hash_real_stream(unit,computed_digest,esp_bounds,3,.false.,discarded_sum,ios)
      if(ios==0)call hash_real_stream(unit,computed_digest,rocc_bounds,3,.false.,discarded_sum,ios)
    endif
    position=0_int64
    if(ios==0)inquire(unit=unit,pos=position,iostat=ios)
    close_ios=0;close(unit,iostat=close_ios)
    if(ios==0.and.close_ios/=0)ios=close_ios
    file_size=0_int64
    if(ios==0)inquire(file=trim(filename),size=file_size,iostat=ios)
    if(ios==0)then
      if(file_size/=expected_size.or.position-1_int64/=file_size.or.&
         computed_digest/=stored_digest.or..not.ieee_is_finite(local_density_sum))ios=1
    endif
    ok=ios==0
    if(.not.ok)then
      local_density_sum=0d0
      file_mu=0d0;file_residual=huge(0d0);file_iteration=-1
    endif
  end subroutine validate_shard_file_stream

  subroutine hash_real_stream(unit,hash,bounds,ndim,accumulate,value_sum,ios)
    integer,intent(in)::unit,ndim
    integer(int64),intent(inout)::hash
    integer,intent(in)::bounds(2*ndim)
    logical,intent(in)::accumulate
    real(8),intent(inout)::value_sum
    integer,intent(inout)::ios
    real(8)::values(stream_chunk_elements)
    integer(int64)::remaining
    integer::count,index,dimension,normalized_bounds(2*ndim)
    if(ios/=0)return
    if(.not.valid_bounds(bounds,ndim))then
      ios=1;return
    endif
    normalized_bounds(1:ndim)=1
    do dimension=1,ndim
      normalized_bounds(ndim+dimension)=bounds(ndim+dimension)-bounds(dimension)+1
    enddo
    call hash_integer_vector(hash,normalized_bounds)
    remaining=element_count(bounds,ndim)
    do while(remaining>0_int64)
      count=int(min(remaining,int(stream_chunk_elements,int64)))
      read(unit,iostat=ios)values(1:count)
      if(ios/=0)return
      if(.not.all(ieee_is_finite(values(1:count))))then
        ios=1;return
      endif
      if(accumulate)value_sum=value_sum+sum(values(1:count))
      do index=1,count
        call hash_real(hash,values(index))
      enddo
      remaining=remaining-int(count,int64)
    enddo
  end subroutine hash_real_stream

  subroutine read_shard_file(filename,contract,publication_id,expected_esp_bounds,&
      expected_rocc_bounds,expected_size,expected_digest,payload,ok)
    character(*),intent(in)::filename
    type(s_dg_dc_seed_contract),intent(in)::contract
    integer(int64),intent(in)::publication_id,expected_size,expected_digest
    integer,intent(in)::expected_esp_bounds(6),expected_rocc_bounds(6)
    type(s_dg_dc_seed_payload),intent(out)::payload
    logical,intent(out)::ok
    character(32)::magic
    integer::unit,ios,file_version,file_nproc,file_rank,file_fragment,iteration,allocation_status
    integer::rwf_bounds(14),rho_bounds(6),vloc_bounds(6),esp_bounds(6),rocc_bounds(6)
    integer(int64)::file_publication,immutable_fingerprint,ownership_fingerprint
    integer(int64)::stored_digest,file_size,position,computed_digest
    real(8)::mu,residual
    call clear_payload(payload);ok=.false.
    open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,file_version,file_nproc,file_rank,file_fragment,file_publication,&
      immutable_fingerprint,ownership_fingerprint,stored_digest,rwf_bounds,rho_bounds,&
      vloc_bounds,esp_bounds,rocc_bounds,mu,residual,iteration
    if(ios==0)then
      if(magic/=shard_magic.or.file_version/=contract%version.or.&
         file_nproc/=contract%mpi_size.or.file_rank/=contract%rank.or.&
         file_fragment/=contract%fragment_id.or.file_publication/=publication_id.or.&
         immutable_fingerprint/=contract%immutable_fingerprint.or.&
         ownership_fingerprint/=contract%ownership_fingerprint.or.&
         stored_digest/=expected_digest.or.any(rwf_bounds/=contract%rwf_bounds).or.&
         any(rho_bounds/=contract%rho_bounds).or.any(vloc_bounds/=contract%vloc_bounds).or.&
         any(esp_bounds/=expected_esp_bounds).or.any(rocc_bounds/=expected_rocc_bounds))ios=1
    endif
    if(ios==0)then
      if(.not.valid_bounds(esp_bounds,3).or..not.valid_bounds(rocc_bounds,3))ios=1
    endif
    allocation_status=0
    if(ios==0)allocate(payload%rwf(&
      rwf_bounds(1):rwf_bounds(8),rwf_bounds(2):rwf_bounds(9),&
      rwf_bounds(3):rwf_bounds(10),rwf_bounds(4):rwf_bounds(11),&
      rwf_bounds(5):rwf_bounds(12),rwf_bounds(6):rwf_bounds(13),&
      rwf_bounds(7):rwf_bounds(14)),stat=allocation_status)
    if(ios==0.and.allocation_status==0)allocate(payload%rho_tot(&
      rho_bounds(1):rho_bounds(4),rho_bounds(2):rho_bounds(5),&
      rho_bounds(3):rho_bounds(6)),stat=allocation_status)
    if(ios==0.and.allocation_status==0)allocate(payload%vloc_tot(&
      vloc_bounds(1):vloc_bounds(4),vloc_bounds(2):vloc_bounds(5),&
      vloc_bounds(3):vloc_bounds(6)),stat=allocation_status)
    if(ios==0.and.allocation_status==0)allocate(payload%esp(&
      esp_bounds(1):esp_bounds(4),esp_bounds(2):esp_bounds(5),&
      esp_bounds(3):esp_bounds(6)),stat=allocation_status)
    if(ios==0.and.allocation_status==0)allocate(payload%rocc(&
      rocc_bounds(1):rocc_bounds(4),rocc_bounds(2):rocc_bounds(5),&
      rocc_bounds(3):rocc_bounds(6)),stat=allocation_status)
    if(allocation_status/=0)ios=1
    if(ios==0)read(unit,iostat=ios)payload%rwf,payload%rho_tot,payload%vloc_tot,&
      payload%esp,payload%rocc
    position=0_int64
    if(ios==0)inquire(unit=unit,pos=position,iostat=ios)
    close(unit)
    file_size=0_int64
    if(ios==0)inquire(file=trim(filename),size=file_size,iostat=ios)
    if(ios==0.and.(file_size/=expected_size.or.position-1_int64/=file_size))ios=1
    payload%mu=mu;payload%residual=residual;payload%iteration=iteration
    if(ios==0)then
      computed_digest=digest_shard(contract,payload,publication_id,esp_bounds,rocc_bounds)
      if(computed_digest/=stored_digest)ios=1
    endif
    ok=ios==0
    if(.not.ok)call clear_payload(payload)
  end subroutine read_shard_file

  subroutine allocate_manifest(manifest,nproc)
    type(s_dg_dc_seed_manifest),intent(inout)::manifest
    integer,intent(in)::nproc
    call clear_manifest(manifest)
    allocate(manifest%fragment_ids(nproc),manifest%rwf_bounds(14,nproc),&
      manifest%rho_bounds(6,nproc),manifest%vloc_bounds(6,nproc),&
      manifest%esp_bounds(6,nproc),manifest%rocc_bounds(6,nproc),&
      manifest%immutable_fingerprints(nproc),manifest%ownership_fingerprints(nproc),&
      manifest%shard_sizes(nproc),manifest%shard_digests(nproc))
  end subroutine allocate_manifest

  subroutine write_manifest_file(filename,manifest,ok)
    character(*),intent(in)::filename
    type(s_dg_dc_seed_manifest),intent(inout)::manifest
    logical,intent(out)::ok
    integer::ios
    integer(int64)::file_size
    type(s_dg_dc_seed_manifest)::verified
    manifest%manifest_size=0_int64;manifest%manifest_digest=0_int64
    call write_manifest_raw(filename,manifest,ios)
    file_size=0_int64
    if(ios==0)inquire(file=trim(filename),size=file_size,iostat=ios)
    manifest%manifest_size=file_size
    manifest%manifest_digest=digest_manifest(manifest)
    if(ios==0)call write_manifest_raw(filename,manifest,ios)
    if(ios==0)call read_manifest_file(filename,manifest%mpi_size,verified,ok)
    if(ios/=0)ok=.false.
    if(ok)ok=verified%manifest_digest==manifest%manifest_digest.and.&
      verified%ordered_digest==manifest%ordered_digest
    call clear_manifest(verified)
  end subroutine write_manifest_file

  subroutine write_manifest_raw(filename,manifest,ios)
    character(*),intent(in)::filename
    type(s_dg_dc_seed_manifest),intent(in)::manifest
    integer,intent(out)::ios
    integer::unit,flush_ios,close_ios
    open(newunit=unit,file=trim(filename),status='replace',access='stream',&
      form='unformatted',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)manifest_magic,manifest%version,manifest%mpi_size,&
      manifest%publication_id,manifest%ordered_digest,manifest%manifest_size,&
      manifest%manifest_digest,manifest%density_weight,manifest%expected_electrons,&
      manifest%write_electron_tolerance,manifest%write_threshold
    if(ios==0)write(unit,iostat=ios)manifest%fragment_ids,manifest%rwf_bounds,&
      manifest%rho_bounds,manifest%vloc_bounds,manifest%esp_bounds,manifest%rocc_bounds,&
      manifest%immutable_fingerprints,manifest%ownership_fingerprints,&
      manifest%shard_sizes,manifest%shard_digests
    flush_ios=0;close_ios=0
    if(ios==0)flush(unit,iostat=flush_ios)
    close(unit,iostat=close_ios)
    if(ios==0.and.flush_ios/=0)ios=flush_ios
    if(ios==0.and.close_ios/=0)ios=close_ios
  end subroutine write_manifest_raw

  subroutine read_manifest_rank_count_collective(comm,filename,stored_nproc,ok)
    integer,intent(in)::comm
    character(*),intent(in)::filename
    integer,intent(out)::stored_nproc
    logical,intent(out)::ok
    character(32)::magic
    integer::rank,ierr,unit,ios,file_version,root_ok

    stored_nproc=0;ok=.false.;root_ok=0;magic='';file_version=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(rank==0)then
      open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',&
        action='read',iostat=ios)
      if(ios==0)then
        read(unit,iostat=ios)magic,file_version,stored_nproc
        close(unit)
      endif
      if(ios==0.and.magic==manifest_magic.and.file_version==seed_version.and.stored_nproc>0)&
        root_ok=1
    endif
    call MPI_Bcast(root_ok,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.root_ok/=1)return
    call MPI_Bcast(stored_nproc,1,MPI_INTEGER,0,comm,ierr)
    ok=ierr==MPI_SUCCESS
  end subroutine read_manifest_rank_count_collective

  subroutine read_manifest_collective(comm,filename,expected_nproc,manifest,ok)
    integer,intent(in)::comm,expected_nproc
    character(*),intent(in)::filename
    type(s_dg_dc_seed_manifest),intent(out)::manifest
    logical,intent(out)::ok
    integer::rank,ierr,root_ok,ios,local_bad,global_bad
    integer(int64)::metadata8(4)
    real(8)::metadata_real(4)

    call clear_manifest(manifest);ok=.false.;root_ok=0
    call MPI_Comm_rank(comm,rank,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    if(rank==0)then
      call read_manifest_file(filename,expected_nproc,manifest,ok)
      root_ok=merge(1,0,ok)
    endif
    call MPI_Bcast(root_ok,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(root_ok==0)then
      call clear_manifest(manifest)
      ok=.false.
      return
    endif

    ios=0
    if(rank/=0)call allocate_manifest_arrays_only(manifest,expected_nproc,ios)
    if(ios/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      call clear_manifest(manifest)
      ok=.false.
      return
    endif

    if(rank==0)then
      metadata8=[manifest%publication_id,manifest%ordered_digest,&
        manifest%manifest_size,manifest%manifest_digest]
      metadata_real=[manifest%density_weight,manifest%expected_electrons,&
        manifest%write_electron_tolerance,manifest%write_threshold]
    endif
    call MPI_Bcast(manifest%version,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%mpi_size,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(metadata8,4,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(metadata_real,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%fragment_ids,expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%rwf_bounds,14*expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%rho_bounds,6*expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%vloc_bounds,6*expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%esp_bounds,6*expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%rocc_bounds,6*expected_nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%immutable_fingerprints,expected_nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%ownership_fingerprints,expected_nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%shard_sizes,expected_nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(manifest%shard_digests,expected_nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      call clear_manifest(manifest)
      ok=.false.
      return
    endif
    if(rank/=0)then
      manifest%publication_id=metadata8(1)
      manifest%ordered_digest=metadata8(2)
      manifest%manifest_size=metadata8(3)
      manifest%manifest_digest=metadata8(4)
      manifest%density_weight=metadata_real(1)
      manifest%expected_electrons=metadata_real(2)
      manifest%write_electron_tolerance=metadata_real(3)
      manifest%write_threshold=metadata_real(4)
    endif
    ok=.true.
  end subroutine read_manifest_collective

  subroutine read_manifest_file(filename,expected_nproc,manifest,ok)
    character(*),intent(in)::filename
    integer,intent(in)::expected_nproc
    type(s_dg_dc_seed_manifest),intent(out)::manifest
    logical,intent(out)::ok
    character(32)::magic
    integer::unit,ios,file_version,file_nproc
    integer(int64)::position,file_size
    call clear_manifest(manifest);ok=.false.
    open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,file_version,file_nproc,manifest%publication_id,&
      manifest%ordered_digest,manifest%manifest_size,manifest%manifest_digest,&
      manifest%density_weight,manifest%expected_electrons,&
      manifest%write_electron_tolerance,manifest%write_threshold
    if(ios==0.and.(magic/=manifest_magic.or.file_version/=seed_version.or.&
       file_nproc/=expected_nproc))ios=1
    if(ios==0)then
      manifest%version=file_version;manifest%mpi_size=file_nproc
      call allocate_manifest_arrays_only(manifest,file_nproc,ios)
    endif
    if(ios==0)read(unit,iostat=ios)manifest%fragment_ids,manifest%rwf_bounds,&
      manifest%rho_bounds,manifest%vloc_bounds,manifest%esp_bounds,manifest%rocc_bounds,&
      manifest%immutable_fingerprints,manifest%ownership_fingerprints,&
      manifest%shard_sizes,manifest%shard_digests
    position=0_int64
    if(ios==0)inquire(unit=unit,pos=position,iostat=ios)
    close(unit)
    file_size=0_int64
    if(ios==0)inquire(file=trim(filename),size=file_size,iostat=ios)
    if(ios==0.and.(manifest%manifest_size/=file_size.or.position-1_int64/=file_size))ios=1
    if(ios==0)then
      if(manifest%manifest_digest/=digest_manifest(manifest))ios=1
    endif
    if(ios==0)then
      if(manifest%ordered_digest/=digest_ordered_manifest(manifest))ios=1
    endif
    ok=ios==0
    if(.not.ok)call clear_manifest(manifest)
  end subroutine read_manifest_file

  subroutine allocate_manifest_arrays_only(manifest,nproc,ios)
    type(s_dg_dc_seed_manifest),intent(inout)::manifest
    integer,intent(in)::nproc
    integer,intent(inout)::ios
    integer::allocation_status
    allocation_status=0
    allocate(manifest%fragment_ids(nproc),manifest%rwf_bounds(14,nproc),&
      manifest%rho_bounds(6,nproc),manifest%vloc_bounds(6,nproc),&
      manifest%esp_bounds(6,nproc),manifest%rocc_bounds(6,nproc),&
      manifest%immutable_fingerprints(nproc),manifest%ownership_fingerprints(nproc),&
      manifest%shard_sizes(nproc),manifest%shard_digests(nproc),stat=allocation_status)
    if(allocation_status/=0)ios=1
  end subroutine allocate_manifest_arrays_only

  logical function same_manifest_contract(manifest,rank,contract)
    type(s_dg_dc_seed_manifest),intent(in)::manifest
    integer,intent(in)::rank
    type(s_dg_dc_seed_contract),intent(in)::contract
    integer::column
    column=rank+1
    same_manifest_contract=manifest%fragment_ids(column)==contract%fragment_id.and.&
      all(manifest%rwf_bounds(:,column)==contract%rwf_bounds).and.&
      all(manifest%rho_bounds(:,column)==contract%rho_bounds).and.&
      all(manifest%vloc_bounds(:,column)==contract%vloc_bounds).and.&
      manifest%immutable_fingerprints(column)==contract%immutable_fingerprint.and.&
      manifest%ownership_fingerprints(column)==contract%ownership_fingerprint
  end function same_manifest_contract

  integer(int64) function digest_ordered_manifest(manifest)result(hash)
    type(s_dg_dc_seed_manifest),intent(in)::manifest
    integer::rank
    hash=int(z'6A09E667F3BCC909',int64)
    call hash_integer(hash,manifest%version);call hash_integer(hash,manifest%mpi_size)
    call hash_integer8(hash,manifest%publication_id)
    do rank=1,manifest%mpi_size
      call hash_integer(hash,rank-1);call hash_integer(hash,manifest%fragment_ids(rank))
      call hash_integer_vector(hash,manifest%rwf_bounds(:,rank))
      call hash_integer_vector(hash,manifest%rho_bounds(:,rank))
      call hash_integer_vector(hash,manifest%vloc_bounds(:,rank))
      call hash_integer_vector(hash,manifest%esp_bounds(:,rank))
      call hash_integer_vector(hash,manifest%rocc_bounds(:,rank))
      call hash_integer8(hash,manifest%immutable_fingerprints(rank))
      call hash_integer8(hash,manifest%ownership_fingerprints(rank))
      call hash_integer8(hash,manifest%shard_sizes(rank))
      call hash_integer8(hash,manifest%shard_digests(rank))
    enddo
  end function digest_ordered_manifest

  integer(int64) function digest_manifest(manifest)result(hash)
    type(s_dg_dc_seed_manifest),intent(in)::manifest
    hash=int(z'BB67AE8584CAA73B',int64)
    call hash_integer(hash,manifest%version);call hash_integer(hash,manifest%mpi_size)
    call hash_integer8(hash,manifest%publication_id);call hash_integer8(hash,manifest%ordered_digest)
    call hash_integer8(hash,manifest%manifest_size)
    call hash_real(hash,manifest%density_weight);call hash_real(hash,manifest%expected_electrons)
    call hash_real(hash,manifest%write_electron_tolerance);call hash_real(hash,manifest%write_threshold)
    call hash_integer_vector(hash,manifest%fragment_ids)
    call hash_integer_matrix(hash,manifest%rwf_bounds);call hash_integer_matrix(hash,manifest%rho_bounds)
    call hash_integer_matrix(hash,manifest%vloc_bounds);call hash_integer_matrix(hash,manifest%esp_bounds)
    call hash_integer_matrix(hash,manifest%rocc_bounds)
    call hash_integer8_vector(hash,manifest%immutable_fingerprints)
    call hash_integer8_vector(hash,manifest%ownership_fingerprints)
    call hash_integer8_vector(hash,manifest%shard_sizes)
    call hash_integer8_vector(hash,manifest%shard_digests)
  end function digest_manifest

  integer(int64) function digest_shard(contract,payload,publication_id,esp_bounds,rocc_bounds)result(hash)
    type(s_dg_dc_seed_contract),intent(in)::contract
    type(s_dg_dc_seed_payload),intent(in)::payload
    integer(int64),intent(in)::publication_id
    integer,intent(in)::esp_bounds(6),rocc_bounds(6)
    call initialize_shard_digest(hash,contract,publication_id,esp_bounds,rocc_bounds,&
      payload%mu,payload%residual,payload%iteration)
    call hash_real_rank7(hash,payload%rwf)
    call hash_real_rank3(hash,payload%rho_tot);call hash_real_rank3(hash,payload%vloc_tot)
    call hash_real_rank3(hash,payload%esp);call hash_real_rank3(hash,payload%rocc)
  end function digest_shard

  subroutine initialize_shard_digest(hash,contract,publication_id,esp_bounds,rocc_bounds,&
      mu,residual,iteration)
    integer(int64),intent(out)::hash
    type(s_dg_dc_seed_contract),intent(in)::contract
    integer(int64),intent(in)::publication_id
    integer,intent(in)::esp_bounds(6),rocc_bounds(6),iteration
    real(8),intent(in)::mu,residual
    hash=int(z'3C6EF372FE94F82B',int64)
    call hash_integer(hash,contract%version);call hash_integer(hash,contract%mpi_size)
    call hash_integer(hash,contract%rank);call hash_integer(hash,contract%fragment_id)
    call hash_integer8(hash,publication_id)
    call hash_integer8(hash,contract%immutable_fingerprint)
    call hash_integer8(hash,contract%ownership_fingerprint)
    call hash_integer_vector(hash,contract%rwf_bounds)
    call hash_integer_vector(hash,contract%rho_bounds)
    call hash_integer_vector(hash,contract%vloc_bounds)
    call hash_integer_vector(hash,esp_bounds);call hash_integer_vector(hash,rocc_bounds)
    call hash_real(hash,mu);call hash_real(hash,residual)
    call hash_integer(hash,iteration)
  end subroutine initialize_shard_digest

  subroutine hash_integer_vector(hash,values)
    integer(int64),intent(inout)::hash
    integer,intent(in)::values(:)
    integer::i
    call hash_integer(hash,size(values))
    do i=1,size(values);call hash_integer(hash,values(i));enddo
  end subroutine hash_integer_vector

  subroutine hash_integer_matrix(hash,values)
    integer(int64),intent(inout)::hash
    integer,intent(in)::values(:,:)
    integer::i,j
    call hash_integer(hash,size(values,1));call hash_integer(hash,size(values,2))
    do j=1,size(values,2);do i=1,size(values,1)
      call hash_integer(hash,values(i,j))
    enddo;enddo
  end subroutine hash_integer_matrix

  subroutine hash_integer8_vector(hash,values)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::values(:)
    integer::i
    call hash_integer(hash,size(values))
    do i=1,size(values);call hash_integer8(hash,values(i));enddo
  end subroutine hash_integer8_vector

  subroutine hash_real_rank3(hash,values)
    integer(int64),intent(inout)::hash
    real(8),intent(in)::values(:,:,:)
    integer::i,j,k
    call hash_integer_vector(hash,[lbound(values),ubound(values)])
    do k=lbound(values,3),ubound(values,3)
      do j=lbound(values,2),ubound(values,2)
        do i=lbound(values,1),ubound(values,1);call hash_real(hash,values(i,j,k));enddo
      enddo
    enddo
  end subroutine hash_real_rank3

  subroutine hash_real_rank7(hash,values)
    integer(int64),intent(inout)::hash
    real(8),intent(in)::values(:,:,:,:,:,:,:)
    integer::i1,i2,i3,i4,i5,i6,i7
    call hash_integer_vector(hash,[lbound(values),ubound(values)])
    do i7=lbound(values,7),ubound(values,7);do i6=lbound(values,6),ubound(values,6)
      do i5=lbound(values,5),ubound(values,5);do i4=lbound(values,4),ubound(values,4)
        do i3=lbound(values,3),ubound(values,3);do i2=lbound(values,2),ubound(values,2)
          do i1=lbound(values,1),ubound(values,1);call hash_real(hash,values(i1,i2,i3,i4,i5,i6,i7));enddo
        enddo;enddo
      enddo;enddo
    enddo;enddo
  end subroutine hash_real_rank7

  subroutine hash_integer(hash,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::value
    call hash_integer8(hash,int(value,int64))
  end subroutine hash_integer

  subroutine hash_integer8(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    integer::ibyte
    do ibyte=0,7;call hash_byte(hash,int(ibits(value,8*ibyte,8)));enddo
    if(hash==0_int64)hash=1_int64
  end subroutine hash_integer8

  subroutine hash_real(hash,value)
    integer(int64),intent(inout)::hash
    real(8),intent(in)::value
    integer(int64)::bits
    bits=transfer(value,bits);call hash_integer8(hash,bits)
  end subroutine hash_real

  subroutine hash_byte(hash,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::value
    integer(int64),parameter::polynomial=int(z'C96C5795D7870F42',int64)
    integer::bit
    hash=ieor(hash,int(iand(value,255),int64))
    do bit=1,8
      if(btest(hash,0))then
        hash=ieor(shiftr(hash,1),polynomial)
      else
        hash=shiftr(hash,1)
      endif
    enddo
  end subroutine hash_byte

  logical function valid_bounds(bounds,ndim)
    integer,intent(in)::bounds(:),ndim
    integer::dimension
    integer(int64)::count,extent
    valid_bounds=size(bounds)==2*ndim
    if(.not.valid_bounds)return
    count=1_int64
    do dimension=1,ndim
      if(bounds(ndim+dimension)<bounds(dimension))then;valid_bounds=.false.;return;endif
      extent=int(bounds(ndim+dimension),int64)-int(bounds(dimension),int64)+1_int64
      if(extent<1_int64.or.count>max_seed_elements/extent)then;valid_bounds=.false.;return;endif
      count=count*extent
    enddo
  end function valid_bounds

  integer(int64) function element_count(bounds,ndim)result(count)
    integer,intent(in)::bounds(:),ndim
    integer::dimension
    count=1_int64
    do dimension=1,ndim
      count=count*(int(bounds(ndim+dimension),int64)-int(bounds(dimension),int64)+1_int64)
    enddo
  end function element_count
#endif

  pure function lower_ascii(value)result(lowered)
    character(*),intent(in)::value
    character(len(value))::lowered
    integer::i,code
    lowered=value
    do i=1,len(value)
      code=iachar(value(i:i))
      if(code>=iachar('A').and.code<=iachar('Z'))lowered(i:i)=achar(code+32)
    enddo
  end function lower_ascii

  function seed_path(directory,name)result(path)
    character(*),intent(in)::directory,name
    character(512)::path
    integer::length
    length=len_trim(directory)
    if(length>0.and.directory(length:length)=='/')then
      path=trim(directory)//trim(name)
    else
      path=trim(directory)//'/'//trim(name)
    endif
  end function seed_path

  function shard_name(directory,publication_id,rank)result(path)
    character(*),intent(in)::directory
    integer(int64),intent(in)::publication_id
    integer,intent(in)::rank
    character(512)::path
    write(path,'(a,"/dg_dc_seed.",z16.16,".rank",i8.8,".shard")')&
      trim(directory),publication_id,rank
  end function shard_name

  subroutine remove_file_if_present(filename)
    character(*),intent(in)::filename
    integer::unit,ios
    open(newunit=unit,file=trim(filename),status='old',iostat=ios)
    if(ios==0)close(unit,status='delete')
  end subroutine remove_file_if_present

  subroutine clear_payload(payload)
    type(s_dg_dc_seed_payload),intent(inout)::payload
    if(allocated(payload%rwf))deallocate(payload%rwf)
    if(allocated(payload%rho_tot))deallocate(payload%rho_tot)
    if(allocated(payload%vloc_tot))deallocate(payload%vloc_tot)
    if(allocated(payload%esp))deallocate(payload%esp)
    if(allocated(payload%rocc))deallocate(payload%rocc)
    payload%mu=0d0;payload%residual=huge(0d0);payload%iteration=0
  end subroutine clear_payload

  subroutine clear_manifest(manifest)
    type(s_dg_dc_seed_manifest),intent(inout)::manifest
    if(allocated(manifest%fragment_ids))deallocate(manifest%fragment_ids)
    if(allocated(manifest%rwf_bounds))deallocate(manifest%rwf_bounds)
    if(allocated(manifest%rho_bounds))deallocate(manifest%rho_bounds)
    if(allocated(manifest%vloc_bounds))deallocate(manifest%vloc_bounds)
    if(allocated(manifest%esp_bounds))deallocate(manifest%esp_bounds)
    if(allocated(manifest%rocc_bounds))deallocate(manifest%rocc_bounds)
    if(allocated(manifest%immutable_fingerprints))deallocate(manifest%immutable_fingerprints)
    if(allocated(manifest%ownership_fingerprints))deallocate(manifest%ownership_fingerprints)
    if(allocated(manifest%shard_sizes))deallocate(manifest%shard_sizes)
    if(allocated(manifest%shard_digests))deallocate(manifest%shard_digests)
  end subroutine clear_manifest
end module dg_dc_seed_checkpoint
