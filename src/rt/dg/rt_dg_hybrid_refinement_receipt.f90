#include "config.h"
module rt_dg_hybrid_refinement_receipt
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_portable_sha256,only:s_dg_sha256_context,dg_sha256_init,dg_sha256_update_int64,&
    dg_sha256_update_real64,dg_sha256_update_logical,dg_sha256_update_character,dg_sha256_final
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  character(32),parameter::manifest_magic='SALMON_DG_REFINEMENT_MANIFEST_V1'
  character(32),parameter::shard_magic='SALMON_DG_REFINEMENT_SHARD_V1'
  integer,parameter::receipt_version=1
  type,public::s_rt_dg_hybrid_refinement_receipt
    integer::version=receipt_version,fragment_id=0,total_solve_count=0,additional_refinement_count=0
    integer(int64)::v5_publication_fingerprint=0_int64,digest(4)=0_int64
    real(real64)::density_change=huge(0d0),energy_change=huge(0d0)
    logical::converged=.false.,exhausted=.false.
    character(64)::exit_reason=''
  end type
  public::write_rt_dg_hybrid_refinement_receipt,read_rt_dg_hybrid_refinement_receipt
  interface
    function c_rename(old_path,new_path) bind(C,name="rename") result(status)
      import::c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
      integer(c_int)::status
    end function
  end interface
contains
  subroutine write_rt_dg_hybrid_refinement_receipt(comm,prefix,receipt,ok,message)
    integer,intent(in)::comm;character(*),intent(in)::prefix
    type(s_rt_dg_hybrid_refinement_receipt),intent(in)::receipt
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,local_bad,global_bad
    integer(int64)::transaction,local_digest(4),common(4),minimum_common(4),maximum_common(4),shard_size
    integer(int64),allocatable::sizes(:),digests(:,:)
    integer,allocatable::fragments(:)
    character(512)::manifest,temporary,shard
    ok=.false.;message='';call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,valid_receipt(receipt,rank,nproc))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid refinement receipt';return;endif
    common=common_receipt_digest(receipt)
    call MPI_Allreduce(common,minimum_common,4,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(common,maximum_common,4,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_common/=maximum_common))then
      message='rank-disagreeing refinement receipt';return
    endif
    if(rank==0)then;call system_clock(count=transaction);if(transaction<=0_int64)transaction=1_int64;endif
    call MPI_Bcast(transaction,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='refinement transaction broadcast failed';return;endif
    call shard_name(prefix,transaction,rank,shard);temporary=trim(shard)//'.temporary'
    local_digest=shard_receipt_digest(receipt,rank,nproc,transaction)
    ios=0;unit=-1
    open(newunit=unit,file=trim(temporary),status='replace',access='stream',form='unformatted',action='write',iostat=ios)
    if(ios==0)write(unit,iostat=ios)shard_magic,receipt_version,rank,nproc,receipt%fragment_id,&
      receipt%total_solve_count,receipt%additional_refinement_count,transaction,&
      receipt%v5_publication_fingerprint,receipt%density_change,receipt%energy_change,&
      receipt%converged,receipt%exhausted,receipt%exit_reason,local_digest
    if(unit/=-1)then;if(ios==0)flush(unit,iostat=ios);close(unit);endif
    if(ios==0)then;inquire(file=trim(temporary),size=shard_size,iostat=ios);endif
    if(ios==0)call atomic_rename(trim(temporary),trim(shard),ios)
    local_bad=merge(0,1,ios==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot atomically publish refinement shard';return;endif
    allocate(sizes(nproc),digests(4,nproc),fragments(nproc))
    call MPI_Gather(shard_size,1,MPI_INTEGER8,sizes,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Gather(local_digest,4,MPI_INTEGER8,digests,4,MPI_INTEGER8,0,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Gather(receipt%fragment_id,1,MPI_INTEGER,fragments,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='refinement manifest gather failed';return;endif
    manifest=trim(prefix)//'.manifest';temporary=trim(manifest)//'.temporary';ios=0;unit=-1
    if(rank==0)then
      open(newunit=unit,file=trim(temporary),status='replace',access='stream',form='unformatted',action='write',iostat=ios)
      if(ios==0)write(unit,iostat=ios)manifest_magic,receipt_version,nproc,transaction,&
        receipt%v5_publication_fingerprint,common,sizes,digests,fragments
      if(unit/=-1)then;if(ios==0)flush(unit,iostat=ios);close(unit);endif
      if(ios==0)call atomic_rename(trim(temporary),trim(manifest),ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then;message='cannot atomically publish refinement manifest';return;endif
    ok=.true.;message=''
#else
    ok=.false.;message='refinement receipt requires MPI'
#endif
  end subroutine

  subroutine read_rt_dg_hybrid_refinement_receipt(comm,prefix,expected_v5_fingerprint,&
      receipt,present,ok,message)
    integer,intent(in)::comm;character(*),intent(in)::prefix;integer(int64),intent(in)::expected_v5_fingerprint
    type(s_rt_dg_hybrid_refinement_receipt),intent(out)::receipt
    logical,intent(out)::present,ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,ios,unit,version,file_nproc,shard_rank,shard_nproc,fragment,total,additional
    integer::local_bad,global_bad
    integer(int64)::transaction,v5,common(4),stored(4),actual(4),actual_size,shard_transaction,shard_v5
    integer(int64),allocatable::sizes(:),digests(:,:)
    integer,allocatable::fragments(:)
    real(real64)::density_change,energy_change
    logical::exists,converged,exhausted
    character(32)::magic
    character(64)::reason
    character(512)::manifest,shard
    receipt=s_rt_dg_hybrid_refinement_receipt();present=.false.;ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    manifest=trim(prefix)//'.manifest';exists=.false.;if(rank==0)inquire(file=trim(manifest),exist=exists)
    call MPI_Bcast(exists,1,MPI_LOGICAL,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='refinement manifest presence broadcast failed';return;endif
    if(.not.exists)then;ok=.true.;return;endif
    present=.true.;allocate(sizes(nproc),digests(4,nproc),fragments(nproc));ios=0;unit=-1
    if(rank==0)then
      open(newunit=unit,file=trim(manifest),status='old',access='stream',form='unformatted',action='read',iostat=ios)
      if(ios==0)read(unit,iostat=ios)magic,version,file_nproc,transaction,v5,common,sizes,digests,fragments
      if(unit/=-1)close(unit)
      if(ios==0.and.(magic/=manifest_magic.or.version/=receipt_version.or.file_nproc/=nproc))ios=1
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then;message='corrupt refinement manifest';return;endif
    call MPI_Bcast(transaction,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(v5,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(common,4,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(sizes,nproc,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(digests,4*nproc,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(fragments,nproc,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(v5/=expected_v5_fingerprint)then;message='refinement receipt v5 binding mismatch';return;endif
    call shard_name(prefix,transaction,rank,shard);ios=0;unit=-1;actual_size=-1_int64
    magic='';version=0;shard_rank=-1;shard_nproc=-1;fragment=0;total=0;additional=0
    shard_transaction=0_int64;shard_v5=0_int64;density_change=huge(0d0);energy_change=huge(0d0)
    converged=.false.;exhausted=.false.;reason='';stored=0_int64
    inquire(file=trim(shard),size=actual_size,iostat=ios)
    if(ios==0)then;if(actual_size/=sizes(rank+1))ios=1;endif
    if(ios==0)open(newunit=unit,file=trim(shard),status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios==0)read(unit,iostat=ios)magic,version,shard_rank,shard_nproc,fragment,total,additional,&
      shard_transaction,shard_v5,density_change,energy_change,converged,exhausted,reason,stored
    if(unit/=-1)close(unit)
    local_bad=merge(0,1,ios==0.and.magic==shard_magic.and.version==receipt_version.and.&
      shard_rank==rank.and.shard_nproc==nproc.and.fragment==fragments(rank+1).and.&
      fragment==rank+1.and.shard_transaction==transaction.and.shard_v5==v5.and.&
      all(stored==digests(:,rank+1)))
    if(local_bad==0)then
      receipt%version=version;receipt%fragment_id=fragment;receipt%total_solve_count=total
      receipt%additional_refinement_count=additional;receipt%v5_publication_fingerprint=shard_v5
      receipt%density_change=density_change;receipt%energy_change=energy_change
      receipt%converged=converged;receipt%exhausted=exhausted;receipt%exit_reason=reason;receipt%digest=stored
      actual=shard_receipt_digest(receipt,rank,nproc,transaction)
      if(any(actual/=stored).or.any(common_receipt_digest(receipt)/=common))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='corrupt or rank-inconsistent refinement shard';return;endif
    ok=.true.;message=''
#else
    receipt=s_rt_dg_hybrid_refinement_receipt();present=.false.;ok=.false.;message='refinement receipt requires MPI'
#endif
  end subroutine

  logical function valid_receipt(receipt,rank,nproc)
    type(s_rt_dg_hybrid_refinement_receipt),intent(in)::receipt;integer,intent(in)::rank,nproc
    valid_receipt=receipt%version==receipt_version.and.receipt%fragment_id==rank+1.and.nproc>=1.and.&
      receipt%v5_publication_fingerprint/=0_int64.and.receipt%total_solve_count>=1.and.&
      receipt%total_solve_count<=4.and.receipt%additional_refinement_count==receipt%total_solve_count-1.and.&
      ieee_is_finite(receipt%density_change).and.receipt%density_change>=0d0.and.&
      ieee_is_finite(receipt%energy_change).and.receipt%energy_change>=0d0.and.&
      receipt%converged.neqv.receipt%exhausted.and.len_trim(receipt%exit_reason)>0
  end function
  function common_receipt_digest(receipt)result(digest)
    type(s_rt_dg_hybrid_refinement_receipt),intent(in)::receipt;integer(int64)::digest(4)
    type(s_dg_sha256_context)::hash
    call dg_sha256_init(hash);call dg_sha256_update_character(hash,'SALMON-DG-REFINEMENT-COMMON-v1')
    call update_common(hash,receipt);call dg_sha256_final(hash,digest)
  end function
  function shard_receipt_digest(receipt,rank,nproc,transaction)result(digest)
    type(s_rt_dg_hybrid_refinement_receipt),intent(in)::receipt;integer,intent(in)::rank,nproc
    integer(int64),intent(in)::transaction;integer(int64)::digest(4);type(s_dg_sha256_context)::hash
    call dg_sha256_init(hash);call dg_sha256_update_character(hash,'SALMON-DG-REFINEMENT-SHARD-v1')
    call dg_sha256_update_int64(hash,transaction);call dg_sha256_update_int64(hash,int(rank,int64))
    call dg_sha256_update_int64(hash,int(nproc,int64));call dg_sha256_update_int64(hash,int(receipt%fragment_id,int64))
    call update_common(hash,receipt);call dg_sha256_final(hash,digest)
  end function
  subroutine update_common(hash,receipt)
    type(s_dg_sha256_context),intent(inout)::hash;type(s_rt_dg_hybrid_refinement_receipt),intent(in)::receipt
    call dg_sha256_update_int64(hash,int(receipt%version,int64))
    call dg_sha256_update_int64(hash,receipt%v5_publication_fingerprint)
    call dg_sha256_update_int64(hash,int(receipt%total_solve_count,int64))
    call dg_sha256_update_int64(hash,int(receipt%additional_refinement_count,int64))
    call dg_sha256_update_real64(hash,receipt%density_change);call dg_sha256_update_real64(hash,receipt%energy_change)
    call dg_sha256_update_logical(hash,receipt%converged);call dg_sha256_update_logical(hash,receipt%exhausted)
    call dg_sha256_update_character(hash,receipt%exit_reason)
  end subroutine
  subroutine shard_name(prefix,transaction,rank,name)
    character(*),intent(in)::prefix;integer(int64),intent(in)::transaction;integer,intent(in)::rank
    character(*),intent(out)::name
    write(name,'(a,".transaction-",i0,".rank-",i8.8)')trim(prefix),transaction,rank
  end subroutine
  subroutine atomic_rename(old_path,new_path,status)
    character(*),intent(in)::old_path,new_path;integer,intent(out)::status
    character(c_char),allocatable::old_c(:),new_c(:);integer::i
    allocate(old_c(len_trim(old_path)+1),new_c(len_trim(new_path)+1))
    do i=1,len_trim(old_path);old_c(i)=old_path(i:i);enddo;old_c(size(old_c))=c_null_char
    do i=1,len_trim(new_path);new_c(i)=new_path(i:i);enddo;new_c(size(new_c))=c_null_char
    status=int(c_rename(old_c,new_c))
  end subroutine
end module
