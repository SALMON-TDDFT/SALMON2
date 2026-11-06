! A distributed FFT plan must use one common layout on its parent communicator.
module exx_fft_consensus
  implicit none
  private
  public :: agree_fft_layout
contains
  subroutine agree_fft_layout(times,contiguous,status,comm,automatic)
    use communication,only:comm_summation,comm_get_max,comm_get_groupinfo
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    real(8),intent(inout)::times(2)
    logical,intent(inout)::contiguous
    logical,intent(in)::automatic
    integer,intent(inout)::status
    integer,intent(in)::comm
    real(8)::worst(2)
    integer::bad,total,i,nauto,ncontiguous,rank,np
    bad=merge(1,0,status/=0)
    do i=1,2
      if(.not.ieee_is_finite(times(i)))bad=1
      if(times(i)<0d0)bad=1
    enddo
    call comm_summation(bad,total,comm)
    status=merge(1,0,total/=0)
    if(status/=0)return
    call comm_get_groupinfo(comm,rank,np)
    call comm_summation(merge(1,0,automatic),nauto,comm)
    if(nauto/=0.and.nauto/=np)then
      status=1;return
    endif
    if(nauto==0)then
      call comm_summation(merge(1,0,contiguous),ncontiguous,comm)
      if(ncontiguous/=0.and.ncontiguous/=np)status=1
      return
    endif
    ! Optimize the slowest rank, preserving the existing noise floor and 5% margin.
    call comm_get_max(times,worst,2,comm)
    times=worst
    contiguous=minval(times)>1d-5.and.times(2)<.95d0*times(1)
  end subroutine
end module
