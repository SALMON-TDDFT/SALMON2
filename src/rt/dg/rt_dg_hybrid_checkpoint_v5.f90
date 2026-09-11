#include "config.h"
module rt_dg_hybrid_checkpoint_v5
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_portable_sha256,only:s_dg_sha256_context,dg_sha256_schema,dg_sha256_init,&
    dg_sha256_update_int64,dg_sha256_final
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  character(32),parameter::manifest_magic='SALMON_HYBRID_DG_MANIFEST_V5'
  character(32),parameter::shard_magic='SALMON_HYBRID_DG_RANK_SHARD_V5'
  integer,parameter::schema_version=5
  type,public::s_rt_dg_hybrid_v5_publication_authorization
    logical::valid=.false.
    integer::checkpoint_version=0,published_rank=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64
  end type s_rt_dg_hybrid_v5_publication_authorization
  public::collective_rt_dg_hybrid_publication_precondition,&
    collective_rt_dg_hybrid_publication_mapping_precondition,publish_rt_dg_hybrid_checkpoint_v5
  type,public::s_rt_dg_hybrid_v5_shard
    integer::global_count=0,global_grid_count=0,nocc=0,certified_rank=0,fragment_id=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64,&
      operator_structure_fingerprint=0_int64,scope_fingerprint=0_int64,payload_fingerprint=0_int64
    integer(int64)::system_fingerprint(4)=0_int64,pseudopotential_digest(4)=0_int64
    integer(int64)::pseudopotential_fingerprint=0_int64
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
  end type s_rt_dg_hybrid_v5_shard
  public::write_rt_dg_hybrid_checkpoint_v5,read_rt_dg_hybrid_checkpoint_v5,&
    checked_rt_dg_hybrid_extent_product
contains
  subroutine collective_rt_dg_hybrid_publication_precondition(comm,local_valid,local_n,local_nocc,ok,message)
    integer,intent(in)::comm,local_n,local_nocc
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_bad,global_bad,local_signature(2),minimum_signature(2),maximum_signature(2)
    local_bad=merge(0,1,local_valid);local_signature=[local_n,local_nocc]
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v5 publication validity reduction failed';return;endif
    call MPI_Allreduce(local_signature,minimum_signature,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v5 publication minimum reduction failed';return;endif
    call MPI_Allreduce(local_signature,maximum_signature,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v5 publication maximum reduction failed';return;endif
    ok=global_bad==0.and.all(minimum_signature==maximum_signature)
    if(ok)then;message='';else;message='terminal divided v5 publication collective precondition failed';endif
#else
    ok=.false.;message='terminal divided v5 publication precondition requires MPI'
#endif
  end subroutine collective_rt_dg_hybrid_publication_precondition

  subroutine collective_rt_dg_hybrid_publication_mapping_precondition(comm,global_count,row_ids,row_owner,&
      occupied_row_ids,local_valid,ok,message)
    integer,intent(in)::comm,global_count,row_owner(:)
    integer(int64),intent(in)::row_ids(:),occupied_row_ids(:)
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,i,ierr,local_bad,global_bad,allocation_status
    integer,allocatable::row_counts(:),owner_minimum(:),owner_maximum(:)
    ok=.false.;message='terminal divided v5 publication stage=pre-gather global row mapping failed'
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(local_valid)then;local_bad=0;else;local_bad=1;endif
    if(global_count<1.or.size(row_owner)/=global_count.or.size(occupied_row_ids)/=size(row_ids))local_bad=1
    if(local_bad==0)then
      if(any(row_owner<0).or.any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or.&
        any(occupied_row_ids<1_int64).or.any(occupied_row_ids>int(global_count,int64)))local_bad=1
    endif
    if(local_bad==0.and.size(row_ids)>0)then
      if(any(row_owner(int(row_ids))/=rank).or.any(occupied_row_ids/=row_ids))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)return
    allocate(row_counts(global_count),owner_minimum(global_count),&
      owner_maximum(global_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)return
    row_counts=0
    do i=1,size(row_ids);row_counts(int(row_ids(i)))=row_counts(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(row_owner,owner_minimum,global_count,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(row_owner,owner_maximum,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,row_counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    ok=global_bad==0.and.all(row_counts==1).and.all(owner_minimum==owner_maximum)
    if(ok)message=''
#else
    ok=.false.;message='terminal divided v5 publication mapping precondition requires MPI'
#endif
  end subroutine collective_rt_dg_hybrid_publication_mapping_precondition

  subroutine publish_rt_dg_hybrid_checkpoint_v5(comm,path,global_count,noccupied,row_ids,row_owner,&
      occupied_row_ids,payload,authorization,local_valid,ok,message)
    integer,intent(in)::comm,global_count,noccupied,row_owner(:)
    character(*),intent(in)::path
    integer(int64),intent(in)::row_ids(:),occupied_row_ids(:)
    type(s_rt_dg_hybrid_v5_shard),intent(in)::payload
    type(s_rt_dg_hybrid_v5_publication_authorization),intent(in)::authorization
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(512)::detail

    call collective_rt_dg_hybrid_publication_precondition(comm,authorization%valid.and.&
      authorization%checkpoint_version==5.and.authorization%published_rank==global_count.and.&
      authorization%basis_fingerprint==payload%basis_fingerprint.and.&
      authorization%operator_fingerprint==payload%operator_fingerprint,global_count,noccupied,ok,detail)
    if(.not.ok)then;message='distributed-v5 endpoint authorization failed: '//trim(detail);return;endif
    call collective_rt_dg_hybrid_publication_precondition(comm,local_valid,global_count,noccupied,ok,detail)
    if(.not.ok)then;message='distributed-v5 endpoint precondition failed: '//trim(detail);return;endif
    call collective_rt_dg_hybrid_publication_mapping_precondition(comm,global_count,row_ids,row_owner,&
      occupied_row_ids,local_valid,ok,detail)
    if(.not.ok)then;message='distributed-v5 endpoint row mapping failed: '//trim(detail);return;endif
    call write_rt_dg_hybrid_checkpoint_v5(comm,path,payload,ok,detail)
    if(.not.ok)then;message='distributed-v5 endpoint publication failed: '//trim(detail);return;endif
    message=''
  end subroutine publish_rt_dg_hybrid_checkpoint_v5
  subroutine write_rt_dg_hybrid_checkpoint_v5(comm,prefix,payload,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    type(s_rt_dg_hybrid_v5_shard),intent(in)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,local_bad,global_bad,flush_ios,allocation_status
    integer(int64)::transaction_id,shard_digest(4),shard_size
    integer(int64),allocatable::shard_sizes(:),shard_digests(:,:)
    integer,allocatable::fragment_ids(:)
    character(512)::manifest,manifest_tmp,shard,shard_tmp
    character(256)::iomsg
    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 writer communicator rank failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 writer communicator size failed';return;endif
    local_bad=validate_local(payload,rank)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 payload validation reduction failed';return;endif
    if(global_bad/=0)then
      message='invalid distributed-v5 rank shard payload';return
    endif
    call validate_common(comm,payload,local_bad)
    if(local_bad/=0)then;message='rank-inconsistent distributed-v5 manifest metadata';return;endif
    if(rank==0)then
      call system_clock(count=transaction_id)
      if(transaction_id<=0_int64)transaction_id=1_int64
    endif
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 transaction broadcast failed';return;endif
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
        payload%scope_fingerprint,payload%payload_fingerprint,payload%system_fingerprint,&
        payload%pseudopotential_fingerprint,payload%pseudopotential_digest,&
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
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard publication reduction failed';return;endif
    if(global_bad/=0)then
      message='cannot atomically publish distributed-v5 rank shard';return
    endif
    allocate(shard_sizes(nproc),shard_digests(4,nproc),fragment_ids(nproc),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 manifest allocation reduction failed';return;endif
    if(global_bad/=0)then;message='cannot allocate distributed-v5 manifest gathers';return;endif
    call MPI_Gather(shard_size,1,MPI_INTEGER8,shard_sizes,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard-size gather failed';return;endif
    call MPI_Gather(shard_digest,4,MPI_INTEGER8,shard_digests,4,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard-digest gather failed';return;endif
    call MPI_Gather(payload%fragment_id,1,MPI_INTEGER,fragment_ids,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 fragment-map gather failed';return;endif
    manifest=trim(prefix)//'.manifest';manifest_tmp=trim(manifest)//'.temporary'
    ios=0;unit=-1
    if(rank==0)then
      open(newunit=unit,file=trim(manifest_tmp),status='replace',access='stream',form='unformatted',&
        action='write',iostat=ios,iomsg=iomsg)
      if(ios==0)write(unit,iostat=ios,iomsg=iomsg)manifest_magic,schema_version,nproc,payload%global_count,&
        payload%global_grid_count,payload%nocc,payload%certified_rank,transaction_id,payload%basis_fingerprint,&
        payload%operator_fingerprint,payload%operator_structure_fingerprint,payload%scope_fingerprint,&
        payload%payload_fingerprint,payload%system_fingerprint,payload%pseudopotential_fingerprint,&
        payload%pseudopotential_digest,&
        shard_sizes,shard_digests,fragment_ids
      flush_ios=0
      if(ios==0)flush(unit,iostat=flush_ios)
      if(ios==0.and.flush_ios/=0)ios=flush_ios
      call close_if_open(unit,ios)
      if(ios==0)call rename(trim(manifest_tmp),trim(manifest),ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 manifest status broadcast failed';return;endif
    if(ios/=0)then
      message='cannot atomically publish distributed-v5 manifest';return
    endif
    ok=.true.
#else
    ok=.false.;message='distributed-v5 checkpoint requires MPI'
#endif
  end subroutine write_rt_dg_hybrid_checkpoint_v5

  subroutine read_rt_dg_hybrid_checkpoint_v5(comm,prefix,payload,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    type(s_rt_dg_hybrid_v5_shard),intent(out)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,file_nproc,version,global_count,global_grid_count,nocc,certified_rank,&
      shard_rank,shard_nproc,&
      shard_global_count,shard_global_grid_count,shard_nocc,shard_certified_rank,&
      fragment_id,nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,&
      nscope,nxc,local_bad,global_bad,allocation_status,manifest_failure_kind
    integer::failure_kind,global_failure_kind
    integer(int64)::transaction_id,basis_fingerprint,operator_fingerprint,operator_structure_fingerprint,&
      scope_fingerprint,payload_fingerprint,system_fingerprint(4),pseudopotential_fingerprint,&
      pseudopotential_digest(4),stored_digest(4),actual_digest(4),actual_size,manifest_size,&
      shard_transaction_id,shard_basis_fingerprint,shard_operator_fingerprint,&
      shard_operator_structure_fingerprint,shard_scope_fingerprint,shard_payload_fingerprint,&
      shard_system_fingerprint(4),shard_pseudopotential_fingerprint,shard_pseudopotential_digest(4)
    integer(int64),allocatable::shard_sizes(:),shard_digests(:,:)
    integer,allocatable::fragment_ids(:)
    character(32)::magic
    character(512)::manifest,shard
    character(256)::iomsg
    ok=.false.;message='';payload=s_rt_dg_hybrid_v5_shard()
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 reader communicator rank failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 reader communicator size failed';return;endif
    manifest=trim(prefix)//'.manifest';ios=0;unit=-1;manifest_failure_kind=0
    allocate(shard_sizes(nproc),shard_digests(4,nproc),fragment_ids(nproc),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 reader allocation reduction failed';return;endif
    if(global_bad/=0)then;message='cannot allocate distributed-v5 manifest metadata';return;endif
    if(rank==0)then
      inquire(file=trim(manifest),size=manifest_size,iostat=ios)
      if(ios==0)open(newunit=unit,file=trim(manifest),status='old',access='stream',form='unformatted',&
        action='read',iostat=ios,iomsg=iomsg)
      if(ios==0)read(unit,pos=1,iostat=ios,iomsg=iomsg)magic,version,file_nproc
      if(ios==0.and.magic=='SALMON_HYBRID_DG_MANIFEST_V4'.and.version==4)then
        ios=1;manifest_failure_kind=2
      elseif(ios==0.and.(magic/=manifest_magic.or.version/=schema_version))then
        ios=1;manifest_failure_kind=3
      elseif(ios==0.and.file_nproc/=nproc)then
        ios=1;manifest_failure_kind=4
      elseif(ios==0.and.manifest_size/=176_int64+44_int64*int(nproc,int64))then
        ios=1;manifest_failure_kind=3
      endif
      if(ios==0)read(unit,pos=1,iostat=ios,iomsg=iomsg)magic,version,file_nproc,global_count,global_grid_count,nocc,&
        certified_rank,transaction_id,basis_fingerprint,operator_fingerprint,operator_structure_fingerprint,&
        scope_fingerprint,payload_fingerprint,system_fingerprint,pseudopotential_fingerprint,pseudopotential_digest,&
        shard_sizes,shard_digests,fragment_ids
      call close_if_open(unit,ios)
      if(ios==0.and.file_nproc/=nproc)ios=1
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 manifest status broadcast failed';return;endif
    call MPI_Bcast(manifest_failure_kind,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 version status broadcast failed';return;endif
    if(ios/=0)then
      if(manifest_failure_kind==2)then
        message='unsupported distributed checkpoint schema v4; regenerate authenticated v5'
      elseif(manifest_failure_kind==3)then
        message='distributed-v5 manifest magic/schema/extent is corrupt'
      elseif(manifest_failure_kind==4)then
        message='distributed-v5 manifest MPI rank mapping changed'
      else
        message='distributed-v5 manifest missing, corrupt, or MPI rank mapping changed'
      endif
      return
    endif
    call MPI_Bcast(global_count,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 global-count broadcast failed';return;endif
    call MPI_Bcast(global_grid_count,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 grid-count broadcast failed';return;endif
    call MPI_Bcast(nocc,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 occupation-count broadcast failed';return;endif
    call MPI_Bcast(certified_rank,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 certified-rank broadcast failed';return;endif
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 transaction broadcast failed';return;endif
    call MPI_Bcast(basis_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 basis fingerprint broadcast failed';return;endif
    call MPI_Bcast(operator_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 operator fingerprint broadcast failed';return;endif
    call MPI_Bcast(operator_structure_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 structure fingerprint broadcast failed';return;endif
    call MPI_Bcast(scope_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 scope fingerprint broadcast failed';return;endif
    call MPI_Bcast(payload_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 payload fingerprint broadcast failed';return;endif
    call MPI_Bcast(system_fingerprint,4,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 system fingerprint broadcast failed';return;endif
    call MPI_Bcast(pseudopotential_fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 pseudopotential fingerprint broadcast failed';return;endif
    call MPI_Bcast(pseudopotential_digest,4,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 pseudopotential digest broadcast failed';return;endif
    call MPI_Bcast(shard_sizes,nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard-size broadcast failed';return;endif
    call MPI_Bcast(shard_digests,4*nproc,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard-digest broadcast failed';return;endif
    call MPI_Bcast(fragment_ids,nproc,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 fragment-map broadcast failed';return;endif
    call shard_name(prefix,transaction_id,rank,shard);ios=0;unit=-1;failure_kind=0
    open(newunit=unit,file=trim(shard),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios,iomsg=iomsg)
    if(ios==0)then
      inquire(unit=unit,size=actual_size)
      if(actual_size<260_int64)then;ios=1;failure_kind=1;endif
      if(ios==0)read(unit,iostat=ios,iomsg=iomsg)magic,version,shard_rank,shard_nproc,fragment_id,shard_global_count,&
        shard_global_grid_count,shard_nocc,shard_certified_rank,shard_transaction_id,shard_basis_fingerprint,&
        shard_operator_fingerprint,shard_operator_structure_fingerprint,shard_scope_fingerprint,&
        shard_payload_fingerprint,shard_system_fingerprint,shard_pseudopotential_fingerprint,&
        shard_pseudopotential_digest,stored_digest,nrow,&
        nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc
    endif
    if(ios/=0.and.failure_kind==0)failure_kind=1
    if(ios==0)then
      if(magic/=shard_magic.or.version/=schema_version.or.shard_rank/=rank.or.shard_nproc/=nproc.or.&
        fragment_id/=fragment_ids(rank+1).or.actual_size/=shard_sizes(rank+1).or.&
        any(stored_digest/=shard_digests(:,rank+1)).or.shard_global_count/=global_count.or.&
        shard_global_grid_count/=global_grid_count.or.shard_nocc/=nocc.or.&
        shard_certified_rank/=certified_rank.or.shard_transaction_id/=transaction_id.or.&
        shard_basis_fingerprint/=basis_fingerprint.or.shard_operator_fingerprint/=operator_fingerprint.or.&
        shard_operator_structure_fingerprint/=operator_structure_fingerprint.or.&
        shard_scope_fingerprint/=scope_fingerprint.or.shard_payload_fingerprint/=payload_fingerprint.or.&
        any(shard_system_fingerprint/=system_fingerprint).or.&
        shard_pseudopotential_fingerprint/=pseudopotential_fingerprint.or.&
        any(shard_pseudopotential_digest/=pseudopotential_digest))then
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
    payload%system_fingerprint=system_fingerprint
    payload%pseudopotential_fingerprint=pseudopotential_fingerprint
    payload%pseudopotential_digest=pseudopotential_digest
    if(ios==0)then
      actual_digest=digest_payload(payload,rank,nproc,transaction_id)
      if(any(actual_digest/=stored_digest))ios=1
    endif
    if(ios/=0.and.failure_kind==0)failure_kind=4
    call MPI_Allreduce(failure_kind,global_failure_kind,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard failure reduction failed';return;endif
    if(global_failure_kind/=0)then
      select case(global_failure_kind)
      case(1);message='distributed-v5 rank shard is truncated or has an invalid fixed header'
      case(2);message='distributed-v5 rank shard disagrees with manifest common metadata'
      case(3);message='distributed-v5 rank shard has negative, overflowing, or invalid dimensions'
      case default;message='distributed-v5 rank shard is partial, stale, or corrupt'
      end select
      return
    endif
    local_bad=validate_local(payload,rank)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed-v5 shard validation reduction failed';return;endif
    if(global_bad/=0)then
      message='invalid distributed-v5 rank shard payload';return
    endif
    call validate_common(comm,payload,local_bad)
    if(local_bad/=0)then;message='rank-inconsistent distributed-v5 shard common metadata';return;endif
    ok=.true.
#else
    ok=.false.;message='distributed-v5 checkpoint requires MPI'
#endif
  end subroutine read_rt_dg_hybrid_checkpoint_v5

  logical function valid_read_dimensions(file_size,global_count,global_grid_count,nocc,certified_rank,&
      nrow,nmetric_offsets,nmetric,noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc)
    integer(int64),intent(in)::file_size
    integer,intent(in)::global_count,global_grid_count,nocc,certified_rank,nrow,nmetric_offsets,nmetric,&
      noperator_offsets,noperator,npoint_offsets,nsupport,ncoeff1,ncoeff2,nscope,nxc
    integer(int64)::expected_size
    logical::extent_ok
    valid_read_dimensions=.false.
    if(global_count<1.or.global_grid_count<1.or.nocc<1.or.certified_rank<nocc.or.&
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
    extent=32_int64*character_bytes+19_int64*integer_bytes+19_int64*int64_bytes
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
    type(s_rt_dg_hybrid_v5_shard),intent(in)::payload
    integer,intent(in)::rank
    bad=0
    if(payload%global_count<1.or.payload%global_grid_count<1.or.payload%nocc<1.or.&
      payload%certified_rank<payload%nocc.or.payload%certified_rank>payload%global_count.or.payload%fragment_id/=rank+1.or.&
      all(payload%system_fingerprint==0_int64).or.payload%pseudopotential_fingerprint==0_int64.or.&
      all(payload%pseudopotential_digest==0_int64))bad=1
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
    type(s_rt_dg_hybrid_v5_shard),intent(in)::payload
    integer,intent(out)::bad
    integer::ierr,imin,imax
    integer(int64)::lmin,lmax
    bad=0
    call MPI_Allreduce(payload%global_count,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%global_count,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(imin/=imax)bad=1
    if(bad/=0)return
    call MPI_Allreduce(payload%global_grid_count,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%global_grid_count,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(imin/=imax)bad=1
    if(bad/=0)return
    call MPI_Allreduce(payload%nocc,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%nocc,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(imin/=imax)bad=1
    if(bad/=0)return
    call MPI_Allreduce(payload%certified_rank,imin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%certified_rank,imax,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(imin/=imax)bad=1
    if(bad/=0)return
    call MPI_Allreduce(payload%basis_fingerprint,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%basis_fingerprint,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(lmin/=lmax)bad=1
    if(bad/=0)return
    call MPI_Allreduce(payload%operator_fingerprint,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(payload%operator_fingerprint,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    if(lmin/=lmax)bad=1
    if(bad/=0)return
    call agree_int64(payload%operator_structure_fingerprint)
    if(bad/=0)return
    call agree_int64(payload%scope_fingerprint)
    if(bad/=0)return
    call agree_int64(payload%payload_fingerprint)
    if(bad/=0)return
    call agree_digest(payload%system_fingerprint)
    if(bad/=0)return
    call agree_int64(payload%pseudopotential_fingerprint)
    if(bad/=0)return
    call agree_digest(payload%pseudopotential_digest)
  contains
    subroutine agree_int64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,lmin,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,lmax,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      if(lmin/=lmax)bad=1
    end subroutine agree_int64
    subroutine agree_digest(value)
      integer(int64),intent(in)::value(4)
      integer::word
      do word=1,4
        call agree_int64(value(word));if(bad/=0)return
      enddo
    end subroutine agree_digest
  end subroutine validate_common
#endif

  function digest_payload(payload,rank,nproc,transaction_id) result(digest)
    type(s_rt_dg_hybrid_v5_shard),intent(in)::payload
    integer,intent(in)::rank,nproc
    integer(int64),intent(in)::transaction_id
    integer::i,j;integer(int64)::bits(2),digest(4);type(s_dg_sha256_context)::hash
    call dg_sha256_init(hash);call mix(dg_sha256_schema)
    call mix(int(rank,int64));call mix(int(nproc,int64));call mix(transaction_id)
    call mix(int(payload%fragment_id,int64));call mix(int(payload%global_count,int64))
    call mix(int(payload%global_grid_count,int64));call mix(int(payload%nocc,int64))
    call mix(int(payload%certified_rank,int64));call mix(payload%basis_fingerprint)
    call mix(payload%operator_fingerprint);call mix(payload%operator_structure_fingerprint)
    call mix(payload%scope_fingerprint);call mix(payload%payload_fingerprint)
    do i=1,4;call mix(payload%system_fingerprint(i));enddo
    call mix(payload%pseudopotential_fingerprint)
    do i=1,4;call mix(payload%pseudopotential_digest(i));enddo
    do i=1,size(payload%row_ids);call mix(payload%row_ids(i));enddo
    do i=1,size(payload%metric_offsets);call mix(int(payload%metric_offsets(i),int64));enddo
    do i=1,size(payload%metric_columns);call mix(int(payload%metric_columns(i),int64));enddo
    do i=1,size(payload%metric_values);bits=transfer(payload%metric_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%operator_offsets);call mix(int(payload%operator_offsets(i),int64));enddo
    do i=1,size(payload%operator_columns);call mix(int(payload%operator_columns(i),int64));enddo
    do i=1,size(payload%operator_values);bits=transfer(payload%operator_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%kinetic_values);bits=transfer(payload%kinetic_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%nonlocal_values);bits=transfer(payload%nonlocal_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%local_values);bits=transfer(payload%local_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%sipg_values);bits=transfer(payload%sipg_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do j=1,size(payload%position_values,2);do i=1,3
      bits=transfer(payload%position_values(i,j),bits);call mix(bits(1));call mix(bits(2))
    enddo;enddo
    do i=1,size(payload%grid_ids);call mix(payload%grid_ids(i));enddo
    do i=1,size(payload%basis_point_offsets);call mix(int(payload%basis_point_offsets(i),int64));enddo
    do i=1,size(payload%basis_support_ids);call mix(int(payload%basis_support_ids(i),int64));enddo
    do i=1,size(payload%basis_support_values);bits=transfer(payload%basis_support_values(i),bits);call mix(bits(1));call mix(bits(2));enddo
    do i=1,size(payload%grid_weights);call mix(transfer(payload%grid_weights(i),bits(1)));enddo
    do i=1,size(payload%density);call mix(transfer(payload%density(i),bits(1)));enddo
    do j=1,size(payload%initial_occupied_amplitudes,2);do i=1,size(payload%initial_occupied_amplitudes,1)
      bits=transfer(payload%initial_occupied_amplitudes(i,j),bits);call mix(bits(1));call mix(bits(2))
    enddo;enddo
    do i=1,size(payload%occupations);call mix(transfer(payload%occupations(i),bits(1)));enddo
    do i=1,size(payload%eigenvalues);call mix(transfer(payload%eigenvalues(i),bits(1)));enddo
    do i=1,size(payload%scope_selectors);call mix(int(payload%scope_selectors(i),int64));enddo
    do i=1,size(payload%xc_types);call mix(int(payload%xc_types(i),int64));enddo
    do i=1,size(payload%acceptance_receipts);call mix(transfer(payload%acceptance_receipts(i),bits(1)));enddo
    do i=1,size(payload%pseudopotential_receipt);call mix(transfer(payload%pseudopotential_receipt(i),bits(1)));enddo
    do i=1,size(payload%energy_receipt);call mix(transfer(payload%energy_receipt(i),bits(1)));enddo
    call dg_sha256_final(hash,digest)
  contains
    subroutine mix(value);integer(int64),intent(in)::value;call dg_sha256_update_int64(hash,value);end subroutine
  end function digest_payload

  subroutine shard_name(prefix,transaction_id,rank,name)
    character(*),intent(in)::prefix
    integer(int64),intent(in)::transaction_id
    integer,intent(in)::rank
    character(*),intent(out)::name
    write(name,'(a,".v5.",z16.16,".rank",i6.6,".shard")')trim(prefix),transaction_id,rank
  end subroutine shard_name
end module rt_dg_hybrid_checkpoint_v5
