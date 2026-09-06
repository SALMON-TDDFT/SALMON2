#include "config.h"
module dg_fragment_wf_checkpoint
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter,public::DG_FRAGMENT_WF_ABSENT=0
  integer,parameter,public::DG_FRAGMENT_WF_VALID=1
  integer,parameter,public::DG_FRAGMENT_WF_INVALID=2
  integer,parameter::format_version=1,message_length=512
  integer(int64),parameter::max_payload_elements=1099511627776_int64
  integer(int64),parameter::max_payload_bytes=1099511627776_int64
  character(32),parameter::manifest_magic='SALMON_DG_FRAGMENT_WF_MANIFEST1'
  character(32),parameter::shard_magic='SALMON_DG_FRAGMENT_WF_SHARD_V1'
  character(32),parameter::footer_magic='SALMON_DG_FRAGMENT_WF_COMPLETE1'

  type,public::s_dg_fragment_wf_contract
    integer::version=1,mpi_size=0,rank=-1,fragment_id=0,basis_generation=0
    integer::gauge_algorithm_version=0,candidate_rank=0,retained_rank=0
    integer::local_row_count=0,seed_count=0
    character(16)::gauge_mode=''
    integer(int64)::mapping_fingerprint=0_int64
    integer(int64)::dc_seed_publication_id=0_int64,dc_seed_fingerprint=0_int64
    integer(int64)::grid_fingerprint=0_int64,cell_fingerprint=0_int64
    integer(int64)::pseudopotential_fingerprint=0_int64
    integer(int64)::fragment_geometry_fingerprint=0_int64,boundary_fingerprint=0_int64
    integer(int64)::inventory_fingerprint=0_int64,ordering_fingerprint=0_int64
    integer(int64)::selection_fingerprint=0_int64,gauge_fingerprint=0_int64
    integer(int64)::local_layout_fingerprint=0_int64
  end type s_dg_fragment_wf_contract

  type,public::s_dg_fragment_wf_payload
    integer(int64),allocatable::local_grid_ids(:),selected_state_ids(:)
    complex(real64),allocatable::wannier_values(:,:),candidate_compression(:,:)
    complex(real64),allocatable::wannier_transform(:,:),dc_seed_coefficients(:,:)
    real(real64),allocatable::centers_fractional(:,:),dc_seed_energies(:),dc_seed_occupations(:)
    real(real64)::seed_reconstruction_defect=huge(1d0)
  end type s_dg_fragment_wf_payload

  type::s_manifest
    integer::mpi_size=0
    integer(int64)::publication_id=0_int64,mapping_fingerprint=0_int64,digest=0_int64
    integer,allocatable::fragment_ids(:)
    integer(int64),allocatable::contract_hashes(:),payload_hashes(:),shard_sizes(:)
  end type s_manifest

  public::write_dg_fragment_wf_checkpoint,read_dg_fragment_wf_checkpoint
  public::probe_dg_fragment_wf_checkpoint,decide_dg_fragment_wf_restart
contains
  subroutine decide_dg_fragment_wf_restart(comm,mode,status,reuse,regenerate,publish,fatal,ok,message)
    integer,intent(in)::comm,status
    character(*),intent(in)::mode
    logical,intent(out)::reuse,regenerate,publish,fatal,ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_code,minimum_code,maximum_code,minimum_status,maximum_status
    character(16)::normalized
    reuse=.false.;regenerate=.false.;publish=.false.;fatal=.false.;ok=.false.;message=''
    normalized=lower_ascii(trim(adjustl(mode)))
    select case(trim(normalized))
    case('off');local_code=1
    case('write');local_code=2
    case('read');local_code=3
    case('auto');local_code=4
    case default;local_code=0
    end select
    call MPI_Allreduce(local_code,minimum_code,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_code,maximum_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(status,minimum_status,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(status,maximum_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_code/=maximum_code.or.minimum_status/=maximum_status)then
      fatal=.true.;message='fragment-WF checkpoint modes or hit states disagree across ranks';return
    endif
    if(status<DG_FRAGMENT_WF_ABSENT.or.status>DG_FRAGMENT_WF_INVALID)then
      fatal=.true.;message='invalid fragment-WF checkpoint hit state';return
    endif
    if(local_code==0)then;fatal=.true.;message='unknown fragment-WF checkpoint mode';return;endif
    select case(local_code)
    case(1);regenerate=.true.
    case(2);regenerate=.true.;publish=.true.
    case(3)
      if(status==DG_FRAGMENT_WF_VALID)then
        reuse=.true.
      else
        fatal=.true.;message='strict fragment-WF checkpoint read missed or rejected the generation'
      endif
    case(4)
      if(status==DG_FRAGMENT_WF_VALID)then
        reuse=.true.
      else
        regenerate=.true.;publish=.true.
      endif
    end select
    ok=.not.fatal
#else
    reuse=.false.;regenerate=.false.;publish=.false.;fatal=.true.;ok=.false.
    message='fragment-WF checkpoint policy requires MPI'
#endif
  end subroutine decide_dg_fragment_wf_restart

  subroutine write_dg_fragment_wf_checkpoint(comm,directory,contract,payload,publication_id,&
      ok,message,failure_injection_rank)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_fragment_wf_contract),intent(in)::contract
    type(s_dg_fragment_wf_payload),intent(in)::payload
    integer(int64),intent(out)::publication_id
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::failure_injection_rank
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_code,global_code,ios,injected,attempt,collision
    integer(int64)::local_contract_hash,local_payload_hash,shard_size,clock_count
    integer,allocatable::fragment_ids(:)
    integer(int64),allocatable::contract_hashes(:),payload_hashes(:),shard_sizes(:)
    character(message_length)::shard,temporary,manifest_path,manifest_temporary,local_message
    type(s_dg_fragment_wf_payload)::verified
    type(s_manifest)::manifest
    logical::file_ok,exists
    ok=.false.;message='';publication_id=0_int64;local_code=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_code=90
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_code=90
    call validate_collective_request(comm,directory,contract,nproc,rank,local_code)
    if(contract%mpi_size/=nproc)then;local_code=90;endif
    call validate_payload(contract,payload,local_code)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      message='invalid distributed fragment-WF checkpoint publication';return
    endif
    local_contract_hash=contract_hash(contract);local_payload_hash=payload_hash(payload)
    if(rank==0)then
      call system_clock(clock_count)
      publication_id=mix_hash(mix_hash(clock_count,contract%mapping_fingerprint),int(nproc,int64))
      if(publication_id==0_int64)publication_id=1_int64
    endif
    collision=1
    do attempt=0,1023
      call MPI_Bcast(publication_id,1,MPI_INTEGER8,0,comm,ierr)
      call make_shard_name(directory,publication_id,rank,shard)
      inquire(file=trim(shard),exist=exists);collision=merge(1,0,exists)
      call MPI_Allreduce(MPI_IN_PLACE,collision,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)collision=1
      if(collision==0)exit
      if(rank==0)publication_id=mix_hash(publication_id,int(attempt+1,int64))
      if(rank==0.and.publication_id==0_int64)publication_id=int(attempt+2,int64)
    enddo
    if(collision/=0)then;message='cannot reserve a collision-free fragment-WF publication ID';return;endif
    temporary=trim(shard)//'.temporary'
    call write_shard_file(temporary,contract,payload,publication_id,local_payload_hash,ios)
    file_ok=ios==0
    if(file_ok)then
      inquire(file=trim(temporary),size=shard_size,iostat=ios);file_ok=ios==0.and.shard_size>0_int64
    endif
    if(file_ok)then
      call read_shard_file(temporary,contract,publication_id,shard_size,verified,local_code,local_message)
      file_ok=local_code==0.and.same_payload(payload,verified)
    endif
    local_code=merge(0,80,file_ok)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      message='cannot write and validate fragment-WF rank payload';call remove_file(temporary);return
    endif
    call rename(trim(temporary),trim(shard),ios)
    local_code=merge(0,80,ios==0)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      message='cannot atomically publish fragment-WF rank payload';call remove_file(temporary);return
    endif
    injected=0;if(present(failure_injection_rank))then
      if(rank==failure_injection_rank)injected=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,injected,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.injected/=0)then
      message='injected failure before fragment-WF manifest commit';return
    endif
    allocate(fragment_ids(nproc),contract_hashes(nproc),payload_hashes(nproc),shard_sizes(nproc),stat=ios)
    call MPI_Allreduce(ios,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then;message='cannot allocate fragment-WF manifest';return;endif
    call MPI_Allgather(contract%fragment_id,1,MPI_INTEGER,fragment_ids,1,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(local_contract_hash,1,MPI_INTEGER8,contract_hashes,1,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(local_payload_hash,1,MPI_INTEGER8,payload_hashes,1,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(shard_size,1,MPI_INTEGER8,shard_sizes,1,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='cannot assemble fragment-WF manifest';return;endif
    manifest%mpi_size=nproc;manifest%publication_id=publication_id
    manifest%mapping_fingerprint=contract%mapping_fingerprint
    manifest%fragment_ids=fragment_ids;manifest%contract_hashes=contract_hashes
    manifest%payload_hashes=payload_hashes;manifest%shard_sizes=shard_sizes
    manifest%digest=manifest_hash(manifest)
    call make_manifest_name(directory,manifest_path)
    write(manifest_temporary,'(a,".temporary.",z16.16)')trim(manifest_path),publication_id
    ios=0;if(rank==0)then
      call write_manifest_file(manifest_temporary,manifest,ios)
      if(ios==0)call rename(trim(manifest_temporary),trim(manifest_path),ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then
      message='cannot commit fragment-WF manifest';return
    endif
    ok=.true.;message=''
#else
    publication_id=0_int64;ok=.false.;message='fragment-WF checkpoints require MPI'
#endif
  end subroutine write_dg_fragment_wf_checkpoint

  subroutine probe_dg_fragment_wf_checkpoint(comm,directory,expected,status,publication_id,message)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_fragment_wf_contract),intent(in)::expected
    integer,intent(out)::status
    integer(int64),intent(out)::publication_id
    character(*),intent(out)::message
    type(s_dg_fragment_wf_payload)::discarded
    call load_checkpoint(comm,directory,expected,discarded,status,publication_id,message)
  end subroutine probe_dg_fragment_wf_checkpoint

  subroutine read_dg_fragment_wf_checkpoint(comm,directory,expected,payload,publication_id,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_fragment_wf_contract),intent(in)::expected
    type(s_dg_fragment_wf_payload),intent(out)::payload
    integer(int64),intent(out)::publication_id
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::status
    call load_checkpoint(comm,directory,expected,payload,status,publication_id,message)
    ok=status==DG_FRAGMENT_WF_VALID
    if(.not.ok)payload=s_dg_fragment_wf_payload()
  end subroutine read_dg_fragment_wf_checkpoint

  subroutine load_checkpoint(comm,directory,expected,payload,status,publication_id,message)
    integer,intent(in)::comm
    character(*),intent(in)::directory
    type(s_dg_fragment_wf_contract),intent(in)::expected
    type(s_dg_fragment_wf_payload),intent(out)::payload
    integer,intent(out)::status
    integer(int64),intent(out)::publication_id
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_code,global_code,ios,stored_nproc
    integer(int64)::file_contract_hash,file_payload_hash
    character(message_length)::manifest_path,shard,local_message
    type(s_manifest)::manifest
    type(s_dg_fragment_wf_contract)::stored
    logical::exists
    payload=s_dg_fragment_wf_payload();status=DG_FRAGMENT_WF_INVALID
    publication_id=0_int64;message='';local_code=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_code=90
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_code=90
    call validate_collective_request(comm,directory,expected,nproc,rank,local_code)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call code_message(max(global_code,90),message);return
    endif
    call make_manifest_name(directory,manifest_path);exists=.false.
    if(rank==0)inquire(file=trim(manifest_path),exist=exists)
    call MPI_Bcast(exists,1,MPI_LOGICAL,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call code_message(90,message);return;endif
    if(.not.exists)then;status=DG_FRAGMENT_WF_ABSENT;message='cause=absent_manifest';return;endif
    ios=0;if(rank==0)call read_manifest_file(manifest_path,manifest,local_code)
    call MPI_Bcast(local_code,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.local_code/=0)then
      call code_message(max(local_code,13),message);return
    endif
    stored_nproc=0;if(rank==0)stored_nproc=manifest%mpi_size
    call MPI_Bcast(stored_nproc,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call code_message(90,message);return;endif
    if(stored_nproc/=nproc)then;call code_message(1,message);return;endif
    ios=0
    if(rank/=0)call allocate_manifest(manifest,nproc,ios)
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then;call code_message(90,message);return;endif
    call MPI_Bcast(manifest%mpi_size,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(manifest%publication_id,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(manifest%mapping_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(manifest%digest,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(manifest%fragment_ids,nproc,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(manifest%contract_hashes,nproc,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(manifest%payload_hashes,nproc,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(manifest%shard_sizes,nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call code_message(90,message);return;endif
    publication_id=manifest%publication_id
    local_code=0
    if(manifest%mpi_size/=nproc)then
      local_code=1
    else if(manifest%mapping_fingerprint/=expected%mapping_fingerprint)then
      local_code=2
    else if(manifest%fragment_ids(rank+1)/=expected%fragment_id)then
      local_code=2
    else if(manifest%digest/=manifest_hash(manifest))then
      local_code=17
    endif
    if(local_code==0)then
      call make_shard_name(directory,publication_id,rank,shard)
      inquire(file=trim(shard),exist=exists)
      if(.not.exists)then
        local_code=14
      else
        call read_shard_file(shard,expected,publication_id,manifest%shard_sizes(rank+1),&
          payload,local_code,local_message,stored_contract=stored,&
          stored_contract_hash=file_contract_hash,stored_payload_hash=file_payload_hash)
        if(local_code==0.and.file_contract_hash/=manifest%contract_hashes(rank+1))local_code=17
        if(local_code==0.and.file_payload_hash/=manifest%payload_hashes(rank+1))local_code=16
      endif
    endif
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      payload=s_dg_fragment_wf_payload();call code_message(max(global_code,90*merge(1,0,ierr/=MPI_SUCCESS)),message)
      return
    endif
    status=DG_FRAGMENT_WF_VALID;message=''
#else
    payload=s_dg_fragment_wf_payload();status=DG_FRAGMENT_WF_INVALID
    publication_id=0_int64;message='fragment-WF checkpoints require MPI'
#endif
  end subroutine load_checkpoint

#ifdef USE_MPI
  subroutine validate_collective_request(comm,directory,contract,nproc,rank,code)
    integer,intent(in)::comm,nproc,rank
    character(*),intent(in)::directory
    type(s_dg_fragment_wf_contract),intent(in)::contract
    integer,intent(inout)::code
    integer::ierr,allocation_status,i
    integer(int64)::path_hash,minimum_hash,maximum_hash,mapping_min,mapping_max
    integer,allocatable::fragments(:),counts(:)
    if(len_trim(directory)<1.or.len_trim(directory)>400)code=max(code,90)
    if(.not.valid_contract(contract,rank))code=max(code,90)
    path_hash=hash_character(trim(directory))
    call MPI_Allreduce(path_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(path_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(contract%mapping_fingerprint,mapping_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(contract%mapping_fingerprint,mapping_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash.or.mapping_min/=mapping_max)code=max(code,90)
    allocate(fragments(nproc),counts(nproc),stat=allocation_status)
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_status/=0)then;code=max(code,90);return;endif
    call MPI_Allgather(contract%fragment_id,1,MPI_INTEGER,fragments,1,MPI_INTEGER,comm,ierr)
    counts=0
    if(ierr==MPI_SUCCESS.and.all(fragments>=1).and.all(fragments<=nproc))then
      do i=1,nproc;counts(fragments(i))=counts(fragments(i))+1;enddo
      if(any(counts/=1))code=max(code,90)
    else
      code=max(code,90)
    endif
  end subroutine validate_collective_request

  logical function valid_contract(value,rank)
    type(s_dg_fragment_wf_contract),intent(in)::value
    integer,intent(in)::rank
    valid_contract=value%version==format_version.and.value%mpi_size>0.and.value%rank==rank.and.&
      value%fragment_id>0.and.value%basis_generation>0.and.value%gauge_algorithm_version>0.and.&
      value%candidate_rank>=value%retained_rank.and.value%retained_rank>0.and.&
      value%local_row_count>0.and.value%seed_count>0.and.len_trim(value%gauge_mode)>0.and.&
      value%mapping_fingerprint/=0_int64.and.value%dc_seed_publication_id/=0_int64.and.&
      value%dc_seed_fingerprint/=0_int64.and.value%grid_fingerprint/=0_int64.and.&
      value%cell_fingerprint/=0_int64.and.value%pseudopotential_fingerprint/=0_int64.and.&
      value%fragment_geometry_fingerprint/=0_int64.and.value%boundary_fingerprint/=0_int64.and.&
      value%inventory_fingerprint/=0_int64.and.value%ordering_fingerprint/=0_int64.and.&
      value%selection_fingerprint/=0_int64.and.value%gauge_fingerprint/=0_int64.and.&
      value%local_layout_fingerprint/=0_int64
    if(valid_contract)valid_contract=valid_payload_extents(value)
  end function valid_contract

  logical function valid_payload_extents(value)result(valid)
    type(s_dg_fragment_wf_contract),intent(in)::value
    integer(int64)::candidate,retained,rows,seeds
    candidate=int(value%candidate_rank,int64);retained=int(value%retained_rank,int64)
    rows=int(value%local_row_count,int64);seeds=int(value%seed_count,int64);valid=.true.
    call extent_fits(retained,rows,valid);call extent_fits(candidate,retained,valid)
    call extent_fits(retained,retained,valid);call extent_fits(3_int64,retained,valid)
    call extent_fits(retained,seeds,valid)
  end function valid_payload_extents

  subroutine extent_fits(left,right,valid)
    integer(int64),intent(in)::left,right
    logical,intent(inout)::valid
    if(.not.valid)return
    if(left<1_int64.or.right<1_int64.or.left>max_payload_elements/right)then
      valid=.false.;return
    endif
    if(left*right>max_payload_elements)valid=.false.
  end subroutine extent_fits

  subroutine validate_payload(contract,payload,code)
    type(s_dg_fragment_wf_contract),intent(in)::contract
    type(s_dg_fragment_wf_payload),intent(in)::payload
    integer,intent(inout)::code
    logical::valid
    valid=allocated(payload%local_grid_ids).and.allocated(payload%selected_state_ids).and.&
      allocated(payload%wannier_values).and.allocated(payload%candidate_compression).and.&
      allocated(payload%wannier_transform).and.allocated(payload%centers_fractional).and.&
      allocated(payload%dc_seed_coefficients).and.allocated(payload%dc_seed_energies).and.&
      allocated(payload%dc_seed_occupations)
    if(.not.valid)then;code=max(code,90);return;endif
    valid=size(payload%local_grid_ids)==contract%local_row_count.and.&
      size(payload%selected_state_ids)==contract%retained_rank.and.&
      all(shape(payload%wannier_values)==[contract%retained_rank,contract%local_row_count]).and.&
      all(shape(payload%candidate_compression)==[contract%candidate_rank,contract%retained_rank]).and.&
      all(shape(payload%wannier_transform)==[contract%retained_rank,contract%retained_rank]).and.&
      all(shape(payload%centers_fractional)==[3,contract%retained_rank]).and.&
      all(shape(payload%dc_seed_coefficients)==[contract%retained_rank,contract%seed_count]).and.&
      size(payload%dc_seed_energies)==contract%seed_count.and.&
      size(payload%dc_seed_occupations)==contract%seed_count
    valid=valid.and.all(payload%local_grid_ids>0_int64).and.all(payload%selected_state_ids>0_int64).and.&
      finite_complex(payload%wannier_values).and.finite_complex(payload%candidate_compression).and.&
      finite_complex(payload%wannier_transform).and.finite_complex(payload%dc_seed_coefficients).and.&
      all(ieee_is_finite(payload%centers_fractional)).and.&
      all(payload%centers_fractional>=0d0).and.all(payload%centers_fractional<1d0).and.&
      all(ieee_is_finite(payload%dc_seed_energies)).and.&
      all(ieee_is_finite(payload%dc_seed_occupations)).and.&
      ieee_is_finite(payload%seed_reconstruction_defect).and.payload%seed_reconstruction_defect>=0d0
    if(.not.valid)code=max(code,90)
  end subroutine validate_payload

  subroutine write_shard_file(path,contract,payload,publication_id,digest,ios)
    character(*),intent(in)::path
    type(s_dg_fragment_wf_contract),intent(in)::contract
    type(s_dg_fragment_wf_payload),intent(in)::payload
    integer(int64),intent(in)::publication_id,digest
    integer,intent(out)::ios
    integer::unit,flush_ios,close_ios
    ios=0
    open(newunit=unit,file=trim(path),status='replace',access='stream',form='unformatted',action='write',iostat=ios)
    if(ios==0)write(unit,iostat=ios)shard_magic,format_version,&
      storage_size(0)/8,storage_size(0_int64)/8,storage_size(0d0)/8,storage_size((0d0,0d0))/8,&
      publication_id,contract%version,contract%mpi_size,contract%rank,contract%fragment_id,&
      contract%basis_generation,contract%gauge_algorithm_version,contract%candidate_rank,&
      contract%retained_rank,contract%local_row_count,contract%seed_count,contract%gauge_mode,&
      contract%mapping_fingerprint,contract%dc_seed_publication_id,contract%dc_seed_fingerprint,&
      contract%grid_fingerprint,contract%cell_fingerprint,contract%pseudopotential_fingerprint,&
      contract%fragment_geometry_fingerprint,contract%boundary_fingerprint,&
      contract%inventory_fingerprint,contract%ordering_fingerprint,contract%selection_fingerprint,&
      contract%gauge_fingerprint,contract%local_layout_fingerprint
    if(ios==0)write(unit,iostat=ios)payload%seed_reconstruction_defect,payload%local_grid_ids,&
      payload%selected_state_ids,payload%wannier_values,payload%candidate_compression,&
      payload%wannier_transform,payload%centers_fractional,payload%dc_seed_coefficients,&
      payload%dc_seed_energies,payload%dc_seed_occupations,footer_magic,digest
    flush_ios=0;close_ios=0;if(ios==0)flush(unit,iostat=flush_ios);close(unit,iostat=close_ios)
    if(ios==0.and.flush_ios/=0)ios=flush_ios;if(ios==0.and.close_ios/=0)ios=close_ios
  end subroutine write_shard_file

  subroutine read_shard_file(path,expected,publication_id,expected_size,payload,code,message,&
      stored_contract,stored_contract_hash,stored_payload_hash)
    character(*),intent(in)::path
    type(s_dg_fragment_wf_contract),intent(in)::expected
    integer(int64),intent(in)::publication_id,expected_size
    type(s_dg_fragment_wf_payload),intent(out)::payload
    integer,intent(out)::code
    character(*),intent(out)::message
    type(s_dg_fragment_wf_contract),intent(out),optional::stored_contract
    integer(int64),intent(out),optional::stored_contract_hash,stored_payload_hash
    type(s_dg_fragment_wf_contract)::file_contract
    character(32)::magic,footer
    integer::unit,ios,version,integer_bytes,int64_bytes,real_bytes,complex_bytes,allocation_status
    integer(int64)::file_publication,file_size,file_digest,computed_digest
    payload=s_dg_fragment_wf_payload();code=0;message=''
    inquire(file=trim(path),size=file_size,iostat=ios)
    if(ios/=0.or.file_size/=expected_size)then;code=15;message='payload size';return;endif
    open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios==0)read(unit,iostat=ios)magic,version,integer_bytes,int64_bytes,real_bytes,complex_bytes,&
      file_publication,file_contract%version,file_contract%mpi_size,file_contract%rank,&
      file_contract%fragment_id,file_contract%basis_generation,file_contract%gauge_algorithm_version,&
      file_contract%candidate_rank,file_contract%retained_rank,file_contract%local_row_count,&
      file_contract%seed_count,file_contract%gauge_mode,file_contract%mapping_fingerprint,&
      file_contract%dc_seed_publication_id,file_contract%dc_seed_fingerprint,file_contract%grid_fingerprint,&
      file_contract%cell_fingerprint,file_contract%pseudopotential_fingerprint,&
      file_contract%fragment_geometry_fingerprint,file_contract%boundary_fingerprint,&
      file_contract%inventory_fingerprint,file_contract%ordering_fingerprint,&
      file_contract%selection_fingerprint,file_contract%gauge_fingerprint,&
      file_contract%local_layout_fingerprint
    if(ios/=0)then;close(unit);code=15;message='payload header';return;endif
    if(magic/=shard_magic.or.version/=format_version.or.file_contract%version/=format_version.or.&
        integer_bytes/=storage_size(0)/8.or.int64_bytes/=storage_size(0_int64)/8.or.&
        real_bytes/=storage_size(0d0)/8.or.complex_bytes/=storage_size((0d0,0d0))/8)then
      close(unit);code=13;message='payload version';return
    endif
    if(file_publication/=publication_id)then;close(unit);code=17;message='publication';return;endif
    code=compatibility_code(expected,file_contract)
    if(code/=0)then;close(unit);call code_message(code,message);return;endif
    if(.not.payload_bytes_fit(file_contract,expected_size))then
      close(unit);code=15;message='payload byte extent';return
    endif
    allocation_status=0
    allocate(payload%local_grid_ids(file_contract%local_row_count),&
      payload%selected_state_ids(file_contract%retained_rank),&
      payload%wannier_values(file_contract%retained_rank,file_contract%local_row_count),&
      payload%candidate_compression(file_contract%candidate_rank,file_contract%retained_rank),&
      payload%wannier_transform(file_contract%retained_rank,file_contract%retained_rank),&
      payload%centers_fractional(3,file_contract%retained_rank),&
      payload%dc_seed_coefficients(file_contract%retained_rank,file_contract%seed_count),&
      payload%dc_seed_energies(file_contract%seed_count),&
      payload%dc_seed_occupations(file_contract%seed_count),stat=allocation_status)
    if(allocation_status/=0)then;close(unit);code=90;message='allocation';return;endif
    read(unit,iostat=ios)payload%seed_reconstruction_defect,payload%local_grid_ids,&
      payload%selected_state_ids,payload%wannier_values,payload%candidate_compression,&
      payload%wannier_transform,payload%centers_fractional,payload%dc_seed_coefficients,&
      payload%dc_seed_energies,payload%dc_seed_occupations,footer,file_digest
    close(unit)
    if(ios/=0.or.footer/=footer_magic)then;payload=s_dg_fragment_wf_payload();code=15;message='payload footer';return;endif
    call validate_payload(file_contract,payload,code)
    if(code/=0)then;payload=s_dg_fragment_wf_payload();code=15;message='payload validation';return;endif
    computed_digest=payload_hash(payload)
    if(computed_digest/=file_digest)then
      payload=s_dg_fragment_wf_payload();code=16;message='payload hash';return
    endif
    if(present(stored_contract))stored_contract=file_contract
    if(present(stored_contract_hash))stored_contract_hash=contract_hash(file_contract)
    if(present(stored_payload_hash))stored_payload_hash=file_digest
  end subroutine read_shard_file

  logical function payload_bytes_fit(value,file_size)result(valid)
    type(s_dg_fragment_wf_contract),intent(in)::value
    integer(int64),intent(in)::file_size
    integer(int64)::bytes,complex_bytes,real_bytes,int64_bytes
    valid=valid_payload_extents(value);bytes=0_int64
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    real_bytes=int(storage_size(0d0)/8,int64);int64_bytes=int(storage_size(0_int64)/8,int64)
    call add_extent_bytes(bytes,int(value%local_row_count,int64),1_int64,int64_bytes,valid)
    call add_extent_bytes(bytes,int(value%retained_rank,int64),1_int64,int64_bytes,valid)
    call add_extent_bytes(bytes,int(value%retained_rank,int64),int(value%local_row_count,int64),complex_bytes,valid)
    call add_extent_bytes(bytes,int(value%candidate_rank,int64),int(value%retained_rank,int64),complex_bytes,valid)
    call add_extent_bytes(bytes,int(value%retained_rank,int64),int(value%retained_rank,int64),complex_bytes,valid)
    call add_extent_bytes(bytes,3_int64,int(value%retained_rank,int64),real_bytes,valid)
    call add_extent_bytes(bytes,int(value%retained_rank,int64),int(value%seed_count,int64),complex_bytes,valid)
    call add_extent_bytes(bytes,2_int64,int(value%seed_count,int64),real_bytes,valid)
    call add_extent_bytes(bytes,1_int64,1_int64,real_bytes,valid)
    valid=valid.and.bytes<=max_payload_bytes.and.bytes<=file_size
  end function payload_bytes_fit

  subroutine add_extent_bytes(total,left,right,element_bytes,valid)
    integer(int64),intent(inout)::total
    integer(int64),intent(in)::left,right,element_bytes
    logical,intent(inout)::valid
    integer(int64)::elements,addition
    if(.not.valid)return
    if(left<1_int64.or.right<1_int64.or.left>huge(elements)/right)then;valid=.false.;return;endif
    elements=left*right
    if(elements>huge(addition)/element_bytes)then;valid=.false.;return;endif
    addition=elements*element_bytes
    if(total>huge(total)-addition)then;valid=.false.;return;endif
    total=total+addition
  end subroutine add_extent_bytes

  integer function compatibility_code(expected,stored)result(code)
    type(s_dg_fragment_wf_contract),intent(in)::expected,stored
    code=0
    if(stored%mpi_size/=expected%mpi_size)then;code=1;return;endif
    if(stored%rank/=expected%rank.or.stored%fragment_id/=expected%fragment_id.or.&
        stored%mapping_fingerprint/=expected%mapping_fingerprint)then;code=2;return;endif
    if(stored%dc_seed_publication_id/=expected%dc_seed_publication_id.or.&
        stored%dc_seed_fingerprint/=expected%dc_seed_fingerprint)then;code=3;return;endif
    if(stored%grid_fingerprint/=expected%grid_fingerprint.or.&
        stored%local_layout_fingerprint/=expected%local_layout_fingerprint.or.&
        stored%local_row_count/=expected%local_row_count)then;code=4;return;endif
    if(stored%cell_fingerprint/=expected%cell_fingerprint)then;code=5;return;endif
    if(stored%pseudopotential_fingerprint/=expected%pseudopotential_fingerprint)then;code=6;return;endif
    if(stored%fragment_geometry_fingerprint/=expected%fragment_geometry_fingerprint)then;code=7;return;endif
    if(stored%boundary_fingerprint/=expected%boundary_fingerprint)then;code=8;return;endif
    if(stored%inventory_fingerprint/=expected%inventory_fingerprint.or.&
        stored%ordering_fingerprint/=expected%ordering_fingerprint.or.&
        stored%candidate_rank/=expected%candidate_rank.or.stored%retained_rank/=expected%retained_rank.or.&
        stored%seed_count/=expected%seed_count)then;code=9;return;endif
    if(stored%selection_fingerprint/=expected%selection_fingerprint)then;code=10;return;endif
    if(stored%basis_generation/=expected%basis_generation)then;code=11;return;endif
    if(stored%gauge_algorithm_version/=expected%gauge_algorithm_version.or.&
        stored%gauge_fingerprint/=expected%gauge_fingerprint.or.&
        trim(stored%gauge_mode)/=trim(expected%gauge_mode))then;code=12;return;endif
    if(stored%version/=expected%version)code=13
  end function compatibility_code

  subroutine write_manifest_file(path,manifest,ios)
    character(*),intent(in)::path
    type(s_manifest),intent(in)::manifest
    integer,intent(out)::ios
    integer::unit,flush_ios,close_ios
    ios=0;open(newunit=unit,file=trim(path),status='replace',access='stream',form='unformatted',&
      action='write',iostat=ios)
    if(ios==0)write(unit,iostat=ios)manifest_magic,format_version,storage_size(0)/8,&
      storage_size(0_int64)/8,storage_size(0d0)/8,storage_size((0d0,0d0))/8,&
      manifest%mpi_size,manifest%publication_id,manifest%mapping_fingerprint,manifest%digest,&
      manifest%fragment_ids,manifest%contract_hashes,manifest%payload_hashes,manifest%shard_sizes,footer_magic
    flush_ios=0;close_ios=0;if(ios==0)flush(unit,iostat=flush_ios);close(unit,iostat=close_ios)
    if(ios==0.and.flush_ios/=0)ios=flush_ios;if(ios==0.and.close_ios/=0)ios=close_ios
  end subroutine write_manifest_file

  subroutine read_manifest_file(path,manifest,code)
    character(*),intent(in)::path
    type(s_manifest),intent(out)::manifest
    integer,intent(out)::code
    integer::unit,ios,version,integer_bytes,int64_bytes,real_bytes,complex_bytes,file_nproc
    character(32)::magic,footer
    code=0;open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios)
    if(ios/=0)then;code=17;return;endif
    read(unit,iostat=ios)magic,version,integer_bytes,int64_bytes,real_bytes,complex_bytes,file_nproc
    if(ios/=0)then;close(unit);code=17;return;endif
    if(magic/=manifest_magic.or.version/=format_version.or.integer_bytes/=storage_size(0)/8.or.&
        int64_bytes/=storage_size(0_int64)/8.or.real_bytes/=storage_size(0d0)/8.or.&
        complex_bytes/=storage_size((0d0,0d0))/8)then;close(unit);code=13;return;endif
    if(file_nproc<1.or.file_nproc>100000)then;close(unit);code=17;return;endif
    call allocate_manifest(manifest,file_nproc,ios)
    if(ios/=0)then;close(unit);code=90;return;endif
    manifest%mpi_size=file_nproc
    read(unit,iostat=ios)manifest%publication_id,manifest%mapping_fingerprint,manifest%digest,&
      manifest%fragment_ids,manifest%contract_hashes,manifest%payload_hashes,manifest%shard_sizes,footer
    close(unit)
    if(ios/=0.or.footer/=footer_magic.or.manifest%digest/=manifest_hash(manifest))code=17
  end subroutine read_manifest_file

  subroutine allocate_manifest(manifest,nproc,ios)
    type(s_manifest),intent(inout)::manifest
    integer,intent(in)::nproc
    integer,intent(out)::ios
    ios=0;manifest%mpi_size=nproc
    allocate(manifest%fragment_ids(nproc),manifest%contract_hashes(nproc),&
      manifest%payload_hashes(nproc),manifest%shard_sizes(nproc),stat=ios)
  end subroutine allocate_manifest

  integer(int64) function contract_hash(value)result(hash)
    type(s_dg_fragment_wf_contract),intent(in)::value
    integer::i
    hash=int(z'6A09E667F3BCC909',int64)
    call add_integer(hash,value%version);call add_integer(hash,value%mpi_size);call add_integer(hash,value%rank)
    call add_integer(hash,value%fragment_id);call add_integer(hash,value%basis_generation)
    call add_integer(hash,value%gauge_algorithm_version);call add_integer(hash,value%candidate_rank)
    call add_integer(hash,value%retained_rank);call add_integer(hash,value%local_row_count)
    call add_integer(hash,value%seed_count)
    do i=1,len(value%gauge_mode);call add_integer(hash,iachar(value%gauge_mode(i:i)));enddo
    hash=mix_hash(hash,value%mapping_fingerprint);hash=mix_hash(hash,value%dc_seed_publication_id)
    hash=mix_hash(hash,value%dc_seed_fingerprint);hash=mix_hash(hash,value%grid_fingerprint)
    hash=mix_hash(hash,value%cell_fingerprint);hash=mix_hash(hash,value%pseudopotential_fingerprint)
    hash=mix_hash(hash,value%fragment_geometry_fingerprint);hash=mix_hash(hash,value%boundary_fingerprint)
    hash=mix_hash(hash,value%inventory_fingerprint);hash=mix_hash(hash,value%ordering_fingerprint)
    hash=mix_hash(hash,value%selection_fingerprint);hash=mix_hash(hash,value%gauge_fingerprint)
    hash=mix_hash(hash,value%local_layout_fingerprint);if(hash==0_int64)hash=1_int64
  end function contract_hash

  integer(int64) function payload_hash(value)result(hash)
    type(s_dg_fragment_wf_payload),intent(in)::value
    integer::i
    hash=int(z'BB67AE8584CAA73B',int64)
    hash=mix_hash(hash,transfer(value%seed_reconstruction_defect,hash))
    do i=1,size(value%local_grid_ids);hash=mix_hash(hash,value%local_grid_ids(i));enddo
    do i=1,size(value%selected_state_ids);hash=mix_hash(hash,value%selected_state_ids(i));enddo
    call hash_complex_matrix(hash,value%wannier_values);call hash_complex_matrix(hash,value%candidate_compression)
    call hash_complex_matrix(hash,value%wannier_transform);call hash_real_matrix(hash,value%centers_fractional)
    call hash_complex_matrix(hash,value%dc_seed_coefficients);call hash_real_vector(hash,value%dc_seed_energies)
    call hash_real_vector(hash,value%dc_seed_occupations);if(hash==0_int64)hash=2_int64
  end function payload_hash

  integer(int64) function manifest_hash(value)result(hash)
    type(s_manifest),intent(in)::value
    integer::i
    hash=int(z'3C6EF372FE94F82B',int64);call add_integer(hash,value%mpi_size)
    hash=mix_hash(hash,value%publication_id);hash=mix_hash(hash,value%mapping_fingerprint)
    do i=1,value%mpi_size
      call add_integer(hash,value%fragment_ids(i));hash=mix_hash(hash,value%contract_hashes(i))
      hash=mix_hash(hash,value%payload_hashes(i));hash=mix_hash(hash,value%shard_sizes(i))
    enddo
    if(hash==0_int64)hash=3_int64
  end function manifest_hash

  subroutine hash_complex_matrix(hash,values)
    integer(int64),intent(inout)::hash;complex(real64),intent(in)::values(:,:);integer::i,j
    do j=1,size(values,2);do i=1,size(values,1)
      hash=mix_hash(hash,transfer(real(values(i,j),real64),hash))
      hash=mix_hash(hash,transfer(aimag(values(i,j)),hash))
    enddo;enddo
  end subroutine hash_complex_matrix
  subroutine hash_real_matrix(hash,values)
    integer(int64),intent(inout)::hash;real(real64),intent(in)::values(:,:);integer::i,j
    do j=1,size(values,2);do i=1,size(values,1);hash=mix_hash(hash,transfer(values(i,j),hash));enddo;enddo
  end subroutine hash_real_matrix
  subroutine hash_real_vector(hash,values)
    integer(int64),intent(inout)::hash;real(real64),intent(in)::values(:);integer::i
    do i=1,size(values);hash=mix_hash(hash,transfer(values(i),hash));enddo
  end subroutine hash_real_vector
  subroutine add_integer(hash,value)
    integer(int64),intent(inout)::hash;integer,intent(in)::value
    hash=mix_hash(hash,int(value,int64))
  end subroutine add_integer
  integer(int64) function mix_hash(hash,value)result(mixed)
    integer(int64),intent(in)::hash,value
    mixed=ieor(ishftc(hash,11),value);mixed=ieor(mixed,ishftc(value,23))
  end function mix_hash
  integer(int64) function hash_character(value)result(hash)
    character(*),intent(in)::value;integer::i
    hash=int(z'A54FF53A5F1D36F1',int64)
    do i=1,len_trim(value);call add_integer(hash,iachar(value(i:i)));enddo
  end function hash_character
  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex
  logical function same_payload(a,b)
    type(s_dg_fragment_wf_payload),intent(in)::a,b
    same_payload=payload_hash(a)==payload_hash(b)
  end function same_payload
  subroutine make_manifest_name(directory,path)
    character(*),intent(in)::directory;character(*),intent(out)::path
    path=trim(directory)//'/dg_fragment_wf.manifest'
  end subroutine make_manifest_name
  subroutine make_shard_name(directory,id,rank,path)
    character(*),intent(in)::directory;integer(int64),intent(in)::id;integer,intent(in)::rank
    character(*),intent(out)::path
    write(path,'(a,"/dg_fragment_wf.publication-",z16.16,".rank-",i6.6,".bin")')trim(directory),id,rank
  end subroutine make_shard_name
  subroutine remove_file(path)
    character(*),intent(in)::path;integer::unit,ios;logical::exists
    inquire(file=trim(path),exist=exists);if(.not.exists)return
    open(newunit=unit,file=trim(path),status='old',iostat=ios);if(ios==0)close(unit,status='delete')
  end subroutine remove_file
  subroutine code_message(code,message)
    integer,intent(in)::code;character(*),intent(out)::message
    select case(code)
    case(1);message='cause=mpi_rank_count'
    case(2);message='cause=rank_fragment_mapping'
    case(3);message='cause=dc_seed'
    case(4);message='cause=grid'
    case(5);message='cause=cell'
    case(6);message='cause=pseudopotential'
    case(7);message='cause=fragment_geometry'
    case(8);message='cause=boundary'
    case(9);message='cause=inventory_or_ordering'
    case(10);message='cause=selection'
    case(11);message='cause=basis_generation'
    case(12);message='cause=gauge'
    case(13);message='cause=version'
    case(14);message='cause=missing_payload'
    case(15);message='cause=payload_format'
    case(16);message='cause=payload_hash'
    case(17);message='cause=manifest_integrity'
    case default;message='cause=collective_contract'
    end select
  end subroutine code_message
  pure character(len=len(value)) function lower_ascii(value)result(lowered)
    character(*),intent(in)::value;integer::i,code
    lowered=value
    do i=1,len(value);code=iachar(lowered(i:i));if(code>=65.and.code<=90)lowered(i:i)=achar(code+32);enddo
  end function lower_ascii
#endif
end module dg_fragment_wf_checkpoint
