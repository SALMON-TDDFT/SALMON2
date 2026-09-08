#include "config.h"
module dg_hybrid_wannier_selection
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::select_dg_hybrid_wannier_blocks
contains
  subroutine select_dg_hybrid_wannier_blocks(comm,nblock,block_ids,conjugate_blocks,localization,&
      threshold,accepted_blocks,rejected_blocks,complement_rank,fingerprint,ok,message)
    integer,intent(in)::comm,nblock
    integer,intent(in)::block_ids(:),conjugate_blocks(:)
    real(real64),intent(in)::localization(:),threshold
    integer,allocatable,intent(out)::accepted_blocks(:),rejected_blocks(:),complement_rank(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,b,norb,naccepted,nrejected,ierr,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer
    integer,allocatable::member_count(:)
    logical,allocatable::eligible(:),selected(:)
    integer(int64)::local_bits,minimum_bits,maximum_bits
    ok=.false.;message='';fingerprint=0_int64;local_bad=0;norb=size(block_ids)
    call agree_integer(nblock,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid Wannier block count';return
    endif
    call agree_integer(norb,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid Wannier orbital count';return
    endif
    local_bits=transfer(threshold,local_bits)
    call agree_int64(local_bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent hybrid Wannier selection threshold';return
    endif
    if(nblock<1.or.norb<1.or.size(localization)/=norb.or.size(conjugate_blocks)/=nblock)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid Wannier selection shape';return
    endif
    do i=1,norb
      call agree_integer(block_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid Wannier block membership';return
      endif
      local_bits=transfer(localization(i),local_bits)
      call agree_int64(local_bits,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
        message='inconsistent hybrid Wannier localization receipt';return
      endif
    enddo
    do b=1,nblock
      call agree_integer(conjugate_blocks(b),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid Wannier conjugate catalog';return
      endif
    enddo
    if(.not.ieee_is_finite(threshold))local_bad=1
    if(.not.all(ieee_is_finite(localization)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='nonfinite hybrid Wannier selection receipt';return
    endif
    if(threshold<=0d0)local_bad=1
    if(any(block_ids<1).or.any(block_ids>nblock))local_bad=1
    if(any(conjugate_blocks<1).or.any(conjugate_blocks>nblock))local_bad=1
    if(any(localization<0d0))local_bad=1
    if(local_bad==0)then
      do b=1,nblock
        if(conjugate_blocks(conjugate_blocks(b))/=b)local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid Wannier selection contract';return
    endif
    allocate(member_count(nblock),eligible(nblock),selected(nblock),complement_rank(nblock),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid Wannier selection workspace';return
    endif
    member_count=0;eligible=.true.
    do i=1,norb
      b=block_ids(i);member_count(b)=member_count(b)+1
      if(localization(i)>threshold)eligible(b)=.false.
    enddo
    if(any(member_count==0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid Wannier catalog has an empty symmetry block';return
    endif
    do b=1,nblock
      selected(b)=eligible(b).and.eligible(conjugate_blocks(b))
    enddo
    naccepted=count(selected);nrejected=nblock-naccepted
    allocate(accepted_blocks(naccepted),rejected_blocks(nrejected),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid Wannier selection outputs';return
    endif
    naccepted=0;nrejected=0;complement_rank=0
    do b=1,nblock
      if(selected(b))then
        naccepted=naccepted+1;accepted_blocks(naccepted)=b
      else
        nrejected=nrejected+1;rejected_blocks(nrejected)=b;complement_rank(b)=member_count(b)
      endif
    enddo
    fingerprint=ieor(int(z'6A09E667F3BCC909',int64),transfer(threshold,local_bits))
    do b=1,nblock
      fingerprint=ishftc(fingerprint,7)
      fingerprint=ieor(fingerprint,int(b,int64))
      fingerprint=ieor(fingerprint,ishft(int(conjugate_blocks(b),int64),16))
      fingerprint=ieor(fingerprint,ishft(int(member_count(b),int64),32))
      if(selected(b))fingerprint=not(fingerprint)
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='hybrid Wannier selection requires MPI';fingerprint=0_int64
    allocate(accepted_blocks(0),rejected_blocks(0),complement_rank(0))
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup()
      if(allocated(member_count))deallocate(member_count)
      if(allocated(eligible))deallocate(eligible)
      if(allocated(selected))deallocate(selected)
      if(allocated(accepted_blocks))deallocate(accepted_blocks)
      if(allocated(rejected_blocks))deallocate(rejected_blocks)
      if(allocated(complement_rank))deallocate(complement_rank)
    end subroutine cleanup
#endif
  end subroutine select_dg_hybrid_wannier_blocks

#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer

  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value
    integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module dg_hybrid_wannier_selection
