! Copyright 2026 SALMON developers. Apache-2.0.
! Portable real64 little-endian Zarr bytes-codec chunk interchange.
! This is not a SALMON native-restart parser or a JSON/manifest implementation.
module salmonpy_gs_stream
  use iso_fortran_env, only: real64,int8,int16,int64,file_storage_size
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_support_datatype
  implicit none
  private
  public :: gs_read_chunk,gs_write_chunk
contains
  logical function supported_host()
    integer(int8) :: bytes(2)
    bytes=transfer(1_int16,bytes)
    supported_host=(bytes(1)==1 .and. storage_size(0._real64)==64 .and. file_storage_size==8 &
                    .and. ieee_support_datatype(0._real64))
  end function

  subroutine gs_read_chunk(filename,values,status)
    character(*),intent(in) :: filename
    real(real64),intent(inout) :: values(:)
    integer,intent(out) :: status
    real(real64),allocatable :: trial(:)
    integer(int8),allocatable :: raw(:)
    integer :: unit,ios,err
    integer(int64) :: bytes,n,i
    status=1
    if(.not.supported_host())return
    n=size(values,kind=int64)
    if(n>huge(n)/8)return
    inquire(file=filename,size=bytes,iostat=ios)
    if(ios/=0)return
    if(bytes/=8*n)return
    allocate(trial(n),raw(8*n),stat=err)
    if(err/=0)then
      status=3;return
    end if
    open(newunit=unit,file=filename,status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)raw
    close(unit,iostat=err)
    if(ios/=0 .or. err/=0)return
    do i=1,n
      trial(i)=transfer(raw(8*(i-1)+1:8*i),trial(i))
      if(.not.ieee_is_finite(trial(i)))then
        status=4;return
      end if
    end do
    values=trial;status=0
  end subroutine

  subroutine gs_write_chunk(filename,values,status)
    character(*),intent(in) :: filename
    real(real64),intent(in) :: values(:)
    integer,intent(out) :: status
    integer :: unit,ios,err
    integer(int64) :: i,n
    integer(int8),allocatable :: raw(:)
    status=1
    if(.not.supported_host())return
    n=size(values,kind=int64)
    if(n>huge(n)/8)return
    do i=1,n
      if(.not.ieee_is_finite(values(i)))then
        status=4;return
      end if
    end do
    allocate(raw(8*n),stat=err)
    if(err/=0)then
      status=3;return
    end if
    do i=1,n
      raw(8*(i-1)+1:8*i)=transfer(values(i),raw(1:8))
    end do
    open(newunit=unit,file=filename,status='new',access='stream',form='unformatted',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)raw
    if(ios/=0)then
      close(unit,status='delete',iostat=err);return
    end if
    close(unit,iostat=ios)
    if(ios/=0)return
    status=0
  end subroutine
end module
