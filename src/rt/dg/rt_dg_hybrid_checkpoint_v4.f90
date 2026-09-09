#include "config.h"
module rt_dg_hybrid_checkpoint_v4
  use,intrinsic::iso_fortran_env,only:int64,real64
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  character(32),parameter::manifest_magic='SALMON_HYBRID_DG_MANIFEST_V4'
  character(32),parameter::shard_magic='SALMON_HYBRID_DG_RANK_SHARD_V4'
  integer,parameter::schema_version=4
  type,public::s_rt_dg_hybrid_v4_shard
    integer::global_count=0,global_grid_count=0,nocc=0,certified_rank=0,fragment_id=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64,&
      operator_structure_fingerprint=0_int64,scope_fingerprint=0_int64,payload_fingerprint=0_int64
    integer(int64),allocatable::row_ids(:),grid_ids(:)
    integer,allocatable::metric_offsets(:),metric_columns(:)
    complex(real64),allocatable::metric_values(:)
    integer,allocatable::operator_offsets(:),operator_columns(:)
    complex(real64),allocatable::operator_values(:),kinetic_values(:),nonlocal_values(:),&
      local_values(:),sipg_values(:),position_values(:,:)
    integer,allocatable::basis_point_offsets(:),basis_support_ids(:)
    complex(real64),allocatable::basis_support_values(:)
    complex(real64),allocatable::initial_occupied_amplitudes(:,:)
    integer,allocatable::scope_selectors(:),xc_types(:)
    real(real64),allocatable::grid_weights(:),density(:),occupations(:),eigenvalues(:),&
      acceptance_receipts(:),pseudopotential_receipt(:),energy_receipt(:)
  end type s_rt_dg_hybrid_v4_shard
  public::write_rt_dg_hybrid_checkpoint_v4,read_rt_dg_hybrid_checkpoint_v4,&
    checked_rt_dg_hybrid_extent_product
contains
  subroutine write_rt_dg_hybrid_checkpoint_v4(comm,prefix,payload,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    type(s_rt_dg_hybrid_v4_shard),intent(in)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,local_bad,global_bad,flush_ios
    integer(int64)::transaction_id,shard_digest,shard_size
    integer(int64),allocatable::shard_sizes(:),shard_digests(:)
    integer,allocatable::fragment_ids(:)
    character(512)::manifest,manifest_tmp,shard,shard_tmp
    character(256)::iomsg
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    ok=.false.;message='';local_bad=validate_local(payload,rank)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid distributed-v4 rank shard payload';return
    endif
    call validate_common(comm,payload,local_bad)
    if(local_bad/=0)then;message='rank-inconsistent distributed-v4 manifest metadata';return;endif
    if(rank==0)then
      call system_clock(count=transaction_id)
      if(transaction_id<=0_int64)transaction_id=1_int64
    endif
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    call shard_name(prefix,transaction_id,rank,shard)
    write(shard_tmp,'(a,".temporary")')trim(shard)
    shard_digest=digest_payload(payload,rank,nproc,transaction_id)
    ios=0;unit=-1
    open(newunit=unit,file=trim(shard_tmp),status='replace',access='stream',form='unformatted',&
      action='write',iostat=ios,iomsg=iomsg)
    if(ios==0)then
      write(unit,iostat=ios,iomsg=iomsg)shard_magic,schema_version,rank,nproc,payload%fragment_id,&
        payload%global_count,payload%global_grid_count,payload%nocc,payload%certified_rank,transaction_id,&
        payload%basis_fingerprint,payload%operator_fingerprint,payload%operator_structure_fingerprint,&
        payload%scope_fingerprint,payload%payload_fingerprint,&
        shard_digest,size(payload%row_ids),size(payload%metric_offsets),size(payload%metric_columns),&
        size(payload%operator_offsets),size(payload%operator_columns),size(payload%basis_point_offsets),&
        size(payload%basis_support_ids),size(payload%initial_occupied_amplitudes,1),&
        size(payload%initial_occupied_amplitudes,2),size(payload%scope_selectors),size(payload%xc_types)
      if(ios==0)write(unit,iostat=ios,iomsg=iomsg)payload%row_ids,payload%metric_offsets,&
        payload%metric_columns,payload%metric_values,payload%operator_offsets,payload%operator_columns,&
        payload%operator_values,payload%kinetic_values,payload%nonlocal_values,payload%local_values,&
        payload%sipg_values,payload%position_values,payload%grid_ids,payload%basis_point_offsets,&
        payload%basis_support_ids,payload%basis_support_values,payload%grid_weights,payload%density,&
        payload%initial_occupied_amplitudes,payload%occupations,payload%eigenvalues,&
        payload%scope_selectors,payload%xc_types,payload%acceptance_receipts,&
        payload%pseudopotential_receipt,payload%energy_receipt
      flush_ios=0
      if(ios==0)flush(unit,iostat=flush_ios)
      if(ios==0.and.flush_ios/=0)ios=flush_ios
    endif
    call close_if_open(unit,ios)
    if(ios==0)then
      inquire(file=trim(shard_tmp),size=shard_size,iostat=ios)
      if(ios==0)call rename(trim(shard_tmp),trim(shard),ios)
    endif
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='cannot atomically publish distributed-v4 rank shard';return
    endif
    allocate(shard_sizes(nproc),shard_digests(nproc),fragment_ids(nproc))
    call MPI_Gather(shard_size,1,MPI_INTEGER8,shard_sizes,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Gather(shard_digest,1,MPI_INTEGER8,shard_digests,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Gather(payload%fragment_id,1,MPI_INTEGER,fragment_ids,1,MPI_INTEGER,0,comm,ierr)
    manifest=trim(prefix)//'.manifest';manifest_tmp=trim(manifest)//'.temporary'
    ios=0;unit=-1
    if(rank==0)then
      open(newunit=unit,file=trim(manifest_tmp),status='replace',access='stream',form='unformatted',&
        action='write',iostat=ios,iomsg=iomsg)
      if(ios==0)write(unit,iostat=ios,iomsg=iomsg)manifest_magic,schema_version,nproc,payload%global_count,&
        payload%global_grid_count,payload%nocc,payload%certified_rank,transaction_id,payload%basis_fingerprint,&
        payload%operator_fingerprint,payload%operator_structure_fingerprint,payload%scope_fingerprint,&
        payload%payload_fingerprint,&
        shard_sizes,shard_digests,fragment_ids
      flush_ios=0
      if(ios==0)flush(unit,iostat=flush_ios)
      if(ios==0.and.flush_ios/=0)ios=flush_ios
      call close_if_open(unit,ios)
      if(ios==0)call rename(trim(manifest_tmp),trim(manifest),ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then
      message='cannot atomically publish distributed-v4 manifest';return
    endif
    ok=.true.
#else
    ok=.false.;message='distributed-v4 checkpoint requires MPI'
#endif
  end subroutine write_rt_dg_hybrid_checkpoint_v4

  subroutine read_rt_dg_hybrid_checkpoint_v4(comm,prefix,payload,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    type(s_rt_dg_hybrid_v4_shard),intent(out)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,file_nproc,version,global_count,global_grid_count,nocc,certified_rank,&
      shard_rank,shard_nproc,&
      shard_global_count,shard_global_grid_count,shard_nocc,shard_certified_rank,&
      fragment_id,nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,&
      nscope,nxc,local_bad,global_bad
    integer::failure_kind,global_failure_kind
    integer(int64)::transaction_id,basis_fingerprint,operator_fingerprint,operator_structure_fingerprint,&
      scope_fingerprint,payload_fingerprint,stored_digest,actual_digest,actual_size,manifest_size,&
      shard_transaction_id,shard_basis_fingerprint,shard_operator_fingerprint,&
      shard_operator_structure_fingerprint,shard_scope_fingerprint,shard_payload_fingerprint
    integer(int64),allocatable::shard_sizes(:),shard_digests(:)
    integer,allocatable::fragment_ids(:)
    character(32)::magic
    character(512)::manifest,shard
    character(256)::iomsg
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    ok=.false.;message='';payload=s_rt_dg_hybrid_v4_shard();manifest=trim(prefix)//'.manifest';ios=0;unit=-1
    allocate(shard_sizes(nproc),shard_digests(nproc),fragment_ids(nproc))
    if(rank==0)then
      inquire(file=trim(manifest),size=manifest_size,iostat=ios)
      if(ios==0.and.manifest_size<104_int64+20_int64*int(nproc,int64))ios=1
      if(ios==0)open(newunit=unit,file=trim(manifest),status='old',access='stream',form='unformatted',&
        action='read',iostat=ios,iomsg=iomsg)
      if(ios==0)read(unit,iostat=ios,iomsg=iomsg)magic,version,file_nproc,global_count,global_grid_count,nocc,&
        certified_rank,transaction_id,basis_fingerprint,operator_fingerprint,operator_structure_fingerprint,&
        scope_fingerprint,payload_fingerprint,shard_sizes,shard_digests,fragment_ids
      call close_if_open(unit,ios)
      if(ios==0.and.(magic/=manifest_magic.or.version/=schema_version.or.file_nproc/=nproc))ios=1
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then
      message='distributed-v4 manifest missing, corrupt, or MPI rank mapping changed';return
    endif
    call MPI_Bcast(global_count,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(global_grid_count,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(nocc,1,MPI_INTEGER,0,comm,ierr);call MPI_Bcast(certified_rank,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(basis_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(operator_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(operator_structure_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(scope_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(payload_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(shard_sizes,nproc,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(shard_digests,nproc,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(fragment_ids,nproc,MPI_INTEGER,0,comm,ierr)
    call shard_name(prefix,transaction_id,rank,shard);ios=0;unit=-1;failure_kind=0
    open(newunit=unit,file=trim(shard),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios,iomsg=iomsg)
    if(ios==0)then
      inquire(unit=unit,size=actual_size)
      if(actual_size<164_int64)then;ios=1;failure_kind=1;endif
      if(ios==0)read(unit,iostat=ios,iomsg=iomsg)magic,version,shard_rank,shard_nproc,fragment_id,shard_global_count,&
        shard_global_grid_count,shard_nocc,shard_certified_rank,shard_transaction_id,shard_basis_fingerprint,&
        shard_operator_fingerprint,shard_operator_structure_fingerprint,shard_scope_fingerprint,&
        shard_payload_fingerprint,stored_digest,nrow,&
        nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc
    endif
    if(ios/=0.and.failure_kind==0)failure_kind=1
    if(ios==0)then
      if(magic/=shard_magic.or.version/=schema_version.or.shard_rank/=rank.or.shard_nproc/=nproc.or.&
        fragment_id/=fragment_ids(rank+1).or.actual_size/=shard_sizes(rank+1).or.&
        stored_digest/=shard_digests(rank+1).or.shard_global_count/=global_count.or.&
        shard_global_grid_count/=global_grid_count.or.shard_nocc/=nocc.or.&
        shard_certified_rank/=certified_rank.or.shard_transaction_id/=transaction_id.or.&
        shard_basis_fingerprint/=basis_fingerprint.or.shard_operator_fingerprint/=operator_fingerprint.or.&
        shard_operator_structure_fingerprint/=operator_structure_fingerprint.or.&
        shard_scope_fingerprint/=scope_fingerprint.or.shard_payload_fingerprint/=payload_fingerprint)then
        ios=1;failure_kind=2
      endif
    endif
    if(ios==0.and..not.valid_read_dimensions(actual_size,global_count,global_grid_count,nocc,certified_rank,&
      nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc))then
      ios=1;failure_kind=3
    endif
    if(ios==0)then
      allocate(payload%row_ids(nrow),payload%metric_offsets(nmetric_offsets),payload%metric_columns(nmetric),&
        payload%metric_values(nmetric),payload%operator_offsets(noperator_offsets),&
        payload%operator_columns(noperator),payload%operator_values(noperator),payload%kinetic_values(noperator),&
        payload%nonlocal_values(noperator),payload%local_values(noperator),payload%sipg_values(noperator),&
        payload%position_values(3,noperator),payload%grid_ids(npoint_offsets-1),&
        payload%basis_point_offsets(npoint_offsets),payload%basis_support_ids(nsupport),&
        payload%basis_support_values(nsupport),payload%grid_weights(npoint_offsets-1),payload%density(npoint_offsets-1),&
        payload%initial_occupied_amplitudes(ncoeff1,ncoeff2),payload%occupations(nocc),payload%eigenvalues(nocc),&
        payload%scope_selectors(nscope),payload%xc_types(nxc),payload%acceptance_receipts(8),&
        payload%pseudopotential_receipt(6),payload%energy_receipt(7),stat=local_bad)
      if(local_bad/=0)ios=1
    endif
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)payload%row_ids,payload%metric_offsets,&
      payload%metric_columns,payload%metric_values,payload%operator_offsets,payload%operator_columns,&
      payload%operator_values,payload%kinetic_values,payload%nonlocal_values,payload%local_values,&
      payload%sipg_values,payload%position_values,payload%grid_ids,payload%basis_point_offsets,&
      payload%basis_support_ids,payload%basis_support_values,payload%grid_weights,payload%density,&
      payload%initial_occupied_amplitudes,payload%occupations,payload%eigenvalues,&
      payload%scope_selectors,payload%xc_types,payload%acceptance_receipts,&
      payload%pseudopotential_receipt,payload%energy_receipt
    call close_if_open(unit,ios)
    payload%global_count=global_count;payload%global_grid_count=global_grid_count
    payload%nocc=nocc;payload%certified_rank=certified_rank;payload%fragment_id=fragment_id
    payload%basis_fingerprint=basis_fingerprint;payload%operator_fingerprint=operator_fingerprint
    payload%operator_structure_fingerprint=operator_structure_fingerprint
    payload%scope_fingerprint=scope_fingerprint;payload%payload_fingerprint=payload_fingerprint
    if(ios==0)then
      actual_digest=digest_payload(payload,rank,nproc,transaction_id)
      if(actual_digest/=stored_digest)ios=1
    endif
    if(ios/=0.and.failure_kind==0)failure_kind=4
    call MPI_Allreduce(failure_kind,global_failure_kind,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure_kind/=0)then
      select case(global_failure_kind)
      case(1);message='distributed-v4 rank shard is truncated or has an invalid fixed header'
      case(2);message='distributed-v4 rank shard disagrees with manifest common metadata'
      case(3);message='distributed-v4 rank shard has negative, overflowing, or invalid dimensions'
      case default;message='distributed-v4 rank shard is partial, stale, or corrupt'
      end select
      return
    endif
    local_bad=validate_local(payload,rank)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid distributed-v4 rank shard payload';return
    endif
    call validate_common(comm,payload,local_bad)
    if(local_bad/=0)then;message='rank-inconsistent distributed-v4 shard common metadata';return;endif
    ok=.true.
#else
    ok=.false.;message='distributed-v4 checkpoint requires MPI'
#endif
  end subroutine read_rt_dg_hybrid_checkpoint_v4

  logical function valid_read_dimensions(file_size,global_count,global_grid_count,nocc,certified_rank,&
      nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc)
    integer(int64),intent(in)::file_size
    integer,intent(in)::global_count,global_grid_count,nocc,certified_rank,nrow,nmetric_offsets,nmetric,&
      noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc
    integer(int64)::expected_size
    logical::extent_ok
    valid_read_dimensions=.false.
    if(global_count<1.or.global_grid_count<1.or.nocc<1.or.certified_rank<1.or.&
      certified_rank>global_count.or.nrow<0.or.nmetric<0.or.noperator<0.or.npoint_offsets<1.or.&
      nsupport<0.or.ncoeff1<0.or.ncoeff2<0.or.nscope<0.or.nxc<0)return
    if(nrow==huge(0).or.nmetric_offsets/=nrow+1.or.noperator_offsets/=nrow+1.or.&
      ncoeff1/=nrow.or.ncoeff2/=nocc)return
    call expected_shard_extent(nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,&
      nsupport,ncoeff1,ncoeff2,nocc,nscope,nxc,expected_size,extent_ok)
    if(.not.extent_ok.or.expected_size/=file_size)return
    valid_read_dimensions=.true.
  end function valid_read_dimensions

  subroutine expected_shard_extent(nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,&
      nsupport,ncoeff1,ncoeff2,nocc,nscope,nxc,extent,ok)
    integer,intent(in)::nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,&
      nsupport,ncoeff1,ncoeff2,nocc,nscope,nxc
    integer(int64),intent(out)::extent
    logical,intent(out)::ok
    integer(int64),parameter::integer_bytes=int(storage_size(0)/8,int64),&
      int64_bytes=int(storage_size(0_int64)/8,int64),real_bytes=int(storage_size(0d0)/8,int64),&
      complex_bytes=int(storage_size(cmplx(0d0,0d0,real64))/8,int64),&
      character_bytes=int(storage_size('a')/8,int64)
    integer(int64)::npoint,ncoefficient
    ok=.false.;extent=0_int64
    if(any([nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,&
      ncoeff1,ncoeff2,nocc,nscope,nxc]<0))return
    npoint=int(npoint_offsets,int64)-1_int64
    if(npoint<0_int64)return
    call checked_product(int(ncoeff1,int64),int(ncoeff2,int64),ncoefficient,ok)
    if(.not.ok)return
    extent=32_int64*character_bytes+19_int64*integer_bytes+7_int64*int64_bytes
    ok=.true.
    call add_extent(extent,int(nrow,int64),int64_bytes,ok)
    call add_extent(extent,int(nmetric_offsets,int64),integer_bytes,ok)
    call add_extent(extent,int(nmetric,int64),integer_bytes,ok)
    call add_extent(extent,int(nmetric,int64),complex_bytes,ok)
    call add_extent(extent,int(noperator_offsets,int64),integer_bytes,ok)
    call add_extent(extent,int(noperator,int64),integer_bytes,ok)
    call add_extent(extent,8_int64*int(noperator,int64),complex_bytes,ok)
    call add_extent(extent,npoint,int64_bytes,ok)
    call add_extent(extent,int(npoint_offsets,int64),integer_bytes,ok)
    call add_extent(extent,int(nsupport,int64),integer_bytes,ok)
    call add_extent(extent,int(nsupport,int64),complex_bytes,ok)
    call add_extent(extent,2_int64*npoint,real_bytes,ok)
    call add_extent(extent,ncoefficient,complex_bytes,ok)
    call add_extent(extent,2_int64*int(nocc,int64),real_bytes,ok)
    call add_extent(extent,int(nscope,int64)+int(nxc,int64),integer_bytes,ok)
    call add_extent(extent,21_int64,real_bytes,ok)
  end subroutine expected_shard_extent

  subroutine checked_rt_dg_hybrid_extent_product(left,right,value,ok)
    integer(int64),intent(in)::left,right
    integer(int64),intent(out)::value
    logical,intent(out)::ok
    call checked_product(left,right,value,ok)
  end subroutine checked_rt_dg_hybrid_extent_product

  subroutine checked_product(left,right,value,ok)
    integer(int64),intent(in)::left,right
    integer(int64),intent(out)::value
    logical,intent(out)::ok
    value=0_int64;ok=left>=0_int64.and.right>=0_int64
    if(.not.ok)return
    if(left/=0_int64.and.right>huge(value)/left)then;ok=.false.;return;endif
    value=left*right
  end subroutine checked_product

  subroutine add_extent(total,count,element_bytes,ok)
    integer(int64),intent(inout)::total
    integer(int64),intent(in)::count,element_bytes
    logical,intent(inout)::ok
    integer(int64)::increment
    if(.not.ok)return
    call checked_product(count,element_bytes,increment,ok)
    if(.not.ok)return
    if(increment>huge(total)-total)then;ok=.false.;return;endif
    total=total+increment
  end subroutine add_extent

  subroutine close_if_open(unit,status)
    integer,intent(inout)::unit,status
    integer::close_status,inquire_status
    logical::opened
    if(unit==-1)return
    inquire(unit=unit,opened=opened,iostat=inquire_status)
    if(inquire_status/=0)then
      if(status==0)status=inquire_status
    else if(opened)then
      close(unit,iostat=close_status)
      if(status==0.and.close_status/=0)status=close_status
    endif
    unit=-1
  end subroutine close_if_open

  integer function validate_local(payload,rank) result(bad)
    type(s_rt_dg_hybrid_v4_shard),intent(in)::payload
    integer,intent(in)::rank
    bad=0
    if(payload%global_count<1.or.payload%global_grid_count<1.or.payload%nocc<1.or.&
      payload%certified_rank<1.or.payload%certified_rank>payload%global_count.or.payload%fragment_id/=rank+1)bad=1
    if(.not.allocated(payload%row_ids).or..not.allocated(payload%metric_offsets).or.&
       .not.allocated(payload%metric_columns).or..not.allocated(payload%metric_values).or.&
       .not.allocated(payload%operator_offsets).or..not.allocated(payload%operator_columns).or.&
       .not.allocated(payload%operator_values).or..not.allocated(payload%kinetic_values).or.&
       .not.allocated(payload%nonlocal_values).or..not.allocated(payload%local_values).or.&
       .not.allocated(payload%sipg_values).or..not.allocated(payload%position_values).or.&
       .not.allocated(payload%grid_ids).or..not.allocated(payload%basis_point_offsets).or.&
       .not.allocated(payload%basis_support_ids).or..not.allocated(payload%basis_support_values).or.&
       .not.allocated(payload%grid_weights).or..not.allocated(payload%density).or.&
       .not.allocated(payload%initial_occupied_amplitudes).or..not.allocated(payload%occupations).or.&
       .not.allocated(payload%eigenvalues).or..not.allocated(payload%scope_selectors).or.&
       .not.allocated(payload%xc_types).or..not.allocated(payload%acceptance_receipts).or.&
       .not.allocated(payload%pseudopotential_receipt).or..not.allocated(payload%energy_receipt))then;bad=1;return;endif
    if(size(payload%metric_offsets)/=size(payload%row_ids)+1.or.&
       size(payload%operator_offsets)/=size(payload%row_ids)+1.or.&
       size(payload%metric_columns)/=size(payload%metric_values).or.&
       size(payload%operator_columns)/=size(payload%operator_values).or.&
       size(payload%kinetic_values)/=size(payload%operator_values).or.&
       size(payload%nonlocal_values)/=size(payload%operator_values).or.&
       size(payload%local_values)/=size(payload%operator_values).or.&
       size(payload%sipg_values)/=size(payload%operator_values).or.&
       any(shape(payload%position_values)/=[3,size(payload%operator_values)]).or.&
       size(payload%basis_point_offsets)<1.or.&
       size(payload%basis_support_ids)/=size(payload%basis_support_values).or.&
       size(payload%grid_ids)/=size(payload%basis_point_offsets)-1.or.&
       size(payload%grid_weights)/=size(payload%grid_ids).or.size(payload%density)/=size(payload%grid_ids).or.&
       size(payload%initial_occupied_amplitudes,1)/=size(payload%row_ids).or.&
       size(payload%initial_occupied_amplitudes,2)/=payload%nocc.or.&
       size(payload%occupations)/=payload%nocc.or.size(payload%eigenvalues)/=payload%nocc.or.&
       size(payload%acceptance_receipts)/=8.or.size(payload%pseudopotential_receipt)/=6.or.&
       size(payload%energy_receipt)/=7)bad=1
    if(bad/=0)return
    if(payload%metric_offsets(1)/=1.or.payload%metric_offsets(size(payload%metric_offsets))/=&
       size(payload%metric_columns)+1.or.payload%operator_offsets(1)/=1.or.&
       payload%operator_offsets(size(payload%operator_offsets))/=size(payload%operator_columns)+1.or.&
       payload%basis_point_offsets(1)/=1.or.payload%basis_point_offsets(size(payload%basis_point_offsets))/=&
       size(payload%basis_support_ids)+1)bad=1
    if(any(payload%row_ids<1_int64).or.any(payload%row_ids>int(payload%global_count,int64)).or.&
       any(payload%metric_columns<1).or.any(payload%metric_columns>payload%global_count).or.&
       any(payload%operator_columns<1).or.any(payload%operator_columns>payload%global_count).or.&
       any(payload%basis_support_ids<1).or.any(payload%basis_support_ids>payload%global_count))bad=1
  end function validate_local

#ifdef USE_MPI
  subroutine validate_common(comm,payload,bad)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_v4_shard),intent(in)::payload
    integer,intent(out)::bad
    integer::ierr,imin,imax
    integer(int64)::lmin,lmax
    bad=0
    call MPI_Allreduce(payload%global_count,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%global_count,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.imin/=imax)bad=1
    call MPI_Allreduce(payload%global_grid_count,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload%global_grid_count,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.imin/=imax)bad=1
    call MPI_Allreduce(payload%nocc,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload%nocc,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.imin/=imax)bad=1
    call MPI_Allreduce(payload%certified_rank,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload%certified_rank,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.imin/=imax)bad=1
    call MPI_Allreduce(payload%basis_fingerprint,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload%basis_fingerprint,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lmin/=lmax)bad=1
    call MPI_Allreduce(payload%operator_fingerprint,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload%operator_fingerprint,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lmin/=lmax)bad=1
    call agree_int64(payload%operator_structure_fingerprint)
    call agree_int64(payload%scope_fingerprint)
    call agree_int64(payload%payload_fingerprint)
  contains
    subroutine agree_int64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.lmin/=lmax)bad=1
    end subroutine agree_int64
  end subroutine validate_common
#endif

  integer(int64) function digest_payload(payload,rank,nproc,transaction_id) result(hash)
    type(s_rt_dg_hybrid_v4_shard),intent(in)::payload
    integer,intent(in)::rank,nproc
    integer(int64),intent(in)::transaction_id
    integer::i,j
    integer(int64)::bits(2)
    hash=1469598103934665603_int64
    call mix(hash,int(rank,int64));call mix(hash,int(nproc,int64));call mix(hash,transaction_id)
    call mix(hash,int(payload%fragment_id,int64));call mix(hash,int(payload%global_count,int64))
    call mix(hash,int(payload%global_grid_count,int64));call mix(hash,int(payload%nocc,int64))
    call mix(hash,int(payload%certified_rank,int64));call mix(hash,payload%basis_fingerprint)
    call mix(hash,payload%operator_fingerprint)
    call mix(hash,payload%operator_structure_fingerprint);call mix(hash,payload%scope_fingerprint)
    call mix(hash,payload%payload_fingerprint)
    do i=1,size(payload%row_ids);call mix(hash,payload%row_ids(i));enddo
    do i=1,size(payload%metric_offsets);call mix(hash,int(payload%metric_offsets(i),int64));enddo
    do i=1,size(payload%metric_columns);call mix(hash,int(payload%metric_columns(i),int64));enddo
    do i=1,size(payload%metric_values);bits=transfer(payload%metric_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%operator_offsets);call mix(hash,int(payload%operator_offsets(i),int64));enddo
    do i=1,size(payload%operator_columns);call mix(hash,int(payload%operator_columns(i),int64));enddo
    do i=1,size(payload%operator_values);bits=transfer(payload%operator_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%kinetic_values);bits=transfer(payload%kinetic_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%nonlocal_values);bits=transfer(payload%nonlocal_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%local_values);bits=transfer(payload%local_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%sipg_values);bits=transfer(payload%sipg_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do j=1,size(payload%position_values,2);do i=1,3
      bits=transfer(payload%position_values(i,j),bits);call mix(hash,bits(1));call mix(hash,bits(2))
    enddo;enddo
    do i=1,size(payload%grid_ids);call mix(hash,payload%grid_ids(i));enddo
    do i=1,size(payload%basis_point_offsets);call mix(hash,int(payload%basis_point_offsets(i),int64));enddo
    do i=1,size(payload%basis_support_ids);call mix(hash,int(payload%basis_support_ids(i),int64));enddo
    do i=1,size(payload%basis_support_values);bits=transfer(payload%basis_support_values(i),bits);call mix(hash,bits(1));call mix(hash,bits(2));enddo
    do i=1,size(payload%grid_weights);call mix(hash,transfer(payload%grid_weights(i),hash));enddo
    do i=1,size(payload%density);call mix(hash,transfer(payload%density(i),hash));enddo
    do j=1,size(payload%initial_occupied_amplitudes,2);do i=1,size(payload%initial_occupied_amplitudes,1)
      bits=transfer(payload%initial_occupied_amplitudes(i,j),bits);call mix(hash,bits(1));call mix(hash,bits(2))
    enddo;enddo
    do i=1,size(payload%occupations);call mix(hash,transfer(payload%occupations(i),hash));enddo
    do i=1,size(payload%eigenvalues);call mix(hash,transfer(payload%eigenvalues(i),hash));enddo
    do i=1,size(payload%scope_selectors);call mix(hash,int(payload%scope_selectors(i),int64));enddo
    do i=1,size(payload%xc_types);call mix(hash,int(payload%xc_types(i),int64));enddo
    do i=1,size(payload%acceptance_receipts);call mix(hash,transfer(payload%acceptance_receipts(i),hash));enddo
    do i=1,size(payload%pseudopotential_receipt);call mix(hash,transfer(payload%pseudopotential_receipt(i),hash));enddo
    do i=1,size(payload%energy_receipt);call mix(hash,transfer(payload%energy_receipt(i),hash));enddo
  end function digest_payload

  subroutine mix(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    hash=ieor(hash,value);hash=hash*1099511628211_int64
  end subroutine mix

  subroutine shard_name(prefix,transaction_id,rank,name)
    character(*),intent(in)::prefix
    integer(int64),intent(in)::transaction_id
    integer,intent(in)::rank
    character(*),intent(out)::name
    write(name,'(a,".v4.",z16.16,".rank",i6.6,".shard")')trim(prefix),transaction_id,rank
  end subroutine shard_name
end module rt_dg_hybrid_checkpoint_v4
