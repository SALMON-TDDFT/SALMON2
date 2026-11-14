! Synthetic Fortran/Python chunk interchange probe; no SALMON runtime or GS.
program test_gs_stream
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  use salmonpy_gs_stream
  implicit none
  real(real64) :: a(72),before(72),bad(72)
  character(4096) :: input,output
  integer :: status,i
  call get_command_argument(1,input);call get_command_argument(2,output)
  a=-77;before=a
  call gs_read_chunk(trim(input),a,status)
  if(status/=0)stop 1
  do i=1,72
    if(a(i)/=real(i,real64)/8)stop 2
  end do
  call gs_write_chunk(trim(output),a,status)
  if(status/=0)stop 3
  call gs_write_chunk(trim(output),a,status)
  if(status==0)stop 4
  call gs_read_chunk(trim(input)//'.wrong',before,status)
  if(status==0 .or. any(before/=-77._real64))stop 5
  call gs_read_chunk(trim(input)//'.nan',before,status)
  if(status/=4 .or. any(before/=-77._real64))stop 6
  call gs_read_chunk(trim(input)//'.inf',before,status)
  if(status/=4 .or. any(before/=-77._real64))stop 7
  bad=a;bad(72)=ieee_value(0._real64,ieee_quiet_nan)
  call gs_write_chunk(trim(output)//'.nan',bad,status)
  if(status/=4)stop 8
  bad(72)=ieee_value(0._real64,ieee_positive_inf)
  call gs_write_chunk(trim(output)//'.inf',bad,status)
  if(status/=4)stop 9
  print *, 'Fortran real64 stream parity and nonfinite rejection PASS'
end program
