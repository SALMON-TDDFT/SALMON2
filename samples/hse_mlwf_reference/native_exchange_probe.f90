program native_exchange_probe
  use iso_fortran_env, only: real64
  use hse_exchange, only: hse_kernel, hse_kernel_init, hse_kernel_apply, hse_kernel_destroy
  use hse_ace
  implicit none
  type(hse_kernel) :: kernel
  type(hse_ace_state) :: ace
  integer :: n,mesh,no,nt,nk,ng,ierr,unit,rank,nproc
  real(real64) :: h,omega,t0,t1
  real(real64), allocatable :: k(:,:)
  complex(real64), allocatable :: source(:,:,:),target(:,:,:),action(:,:,:),w(:,:,:),interpolated(:,:,:)
  character(1024) :: input,output,arg
  call get_command_argument(1,input);call get_command_argument(2,output)
  call get_command_argument(3,arg);read(arg,*)rank
  call get_command_argument(4,arg);read(arg,*)nproc
  open(newunit=unit,file=trim(input),access='stream',form='unformatted',status='old')
  read(unit)n,mesh,no,nt;read(unit)h,omega
  nk=mesh**3;ng=n**3
  allocate(k(3,nk),source(ng,no,nk),target(ng,nt,nk),action(ng,nt,nk))
  read(unit)k,source,target;close(unit)
  call cpu_time(t0)
  call hse_kernel_init(kernel,n,mesh,h,k,omega,4,ierr)
  if(ierr/=0)error stop 'kernel initialization failed'
  call hse_kernel_apply(kernel,source,target,action,rank,nproc,ierr)
  if(ierr/=0)error stop 'kernel action failed'
  call cpu_time(t1)
  open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  write(unit)action;close(unit)
  print *, 'CPU seconds',t1-t0
  if(nproc==1)then
    allocate(w(ng,no,nk),interpolated(ng,no,nk))
    call hse_kernel_apply(kernel,source,source,w,0,1,ierr)
    if(ierr/=0)error stop 'occupied action failed'
    call hse_ace_build(ace,source,w,h**3,ierr)
    if(ierr/=0)error stop 'ACE construction failed'
    call hse_ace_apply(ace,target,action,ierr)
    if(ierr/=0)error stop 'ACE apply failed'
    open(newunit=unit,file=trim(output)//'.ace',access='stream',form='unformatted',status='replace')
    write(unit)action;close(unit)
    call hse_ace_apply(ace,source,interpolated,ierr)
    if(ierr/=0)error stop 'ACE interpolation failed'
    open(newunit=unit,file=trim(output)//'.occupied',access='stream',form='unformatted',status='replace')
    write(unit)interpolated;close(unit)
  endif
  call hse_kernel_destroy(kernel)
end program
