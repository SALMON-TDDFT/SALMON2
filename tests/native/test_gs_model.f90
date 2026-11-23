program test_gs_model
  use salmonpy_gs_model, only: s_gs_model
  use math_constants, only: set_math_constants
  use parallelization, only: setup_parallel,end_parallel,nproc_id_global
  use mpi, only: MPI_Init_thread,MPI_THREAD_FUNNELED,MPI_COMM_WORLD,MPI_Finalize
  use communication, only: comm_attach,comm_child_counts
  use sym_sub, only: read_sw_symmetry,swap_symmetry_switch
  use timer, only: timer_initialize
  use salmon_global, only: xc,num_rgrid,nelec,calc_mode
  implicit none
  type(s_gs_model) :: a,b,empty
  character(4096) :: pa,pb
  character(:),allocatable :: first,third
  integer :: status,u,n,provided,mpi_status,live,slots,slots_first,i
  logical :: sym,dirs(3)
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,mpi_status)
  if(mpi_status/=0)stop 11
  call comm_attach(MPI_COMM_WORLD,mpi_status)
  if(mpi_status/=0)stop 12
  call setup_parallel
  call set_math_constants
  call timer_initialize
  call get_command_argument(1,pa)
  call get_command_argument(2,pb)
  call empty%run(status)
  if(status/=1)stop 1
  call a%load(trim(pa),status)
  if(status/=0)stop 2
  call b%load(trim(pb),status)
  if(status/=0)stop 3
  call a%load(trim(pa)//'.unsupported',status)
  if(status/=2)stop 10
  ! A rejected replacement must leave the existing model usable.
  ! The run must use A's saved values, not B or poisoned shared settings.
  xc='invalid';num_rgrid=-1;nelec=-1;calc_mode='bad'
  call read_sw_symmetry('yyy')
  call a%run(status)
  if(status/=0)stop 4
  if(xc/='invalid'.or.any(num_rgrid/=-1).or.nelec/=-1.or.calc_mode/='bad')stop 5
  call comm_child_counts(live,slots_first)
  if(live/=0)stop 13
  sym=.false.;dirs=.false.
  call swap_symmetry_switch(sym,dirs)
  if(.not.sym.or..not.all(dirs))stop 14
  call swap_symmetry_switch(sym,dirs)
  call save_status('a-first-status.dat')
  if(nproc_id_global==0)then
    inquire(file='gs_status.dat',size=n)
    allocate(character(n)::first,third)
    open(newunit=u,file='gs_status.dat',access='stream',form='unformatted',status='old')
    read(u)first
    close(u)
  endif
  call b%run(status)
  if(status/=0)stop 6
  call comm_child_counts(live,slots)
  if(live/=0)stop 17
  slots_first=max(slots_first,slots)
  call save_status('b-status.dat')
  call b%clear()
  call a%run(status)
  if(status/=0)stop 7
  if(nproc_id_global==0)then
    open(newunit=u,file='gs_status.dat',access='stream',form='unformatted',status='old')
    read(u)third
    close(u)
    if(first/=third)stop 8
    write(*,*)'GS_MODEL_A_B_A_EXACT'
  endif
  call save_status('a-final-status.dat')
  do i=1,6
    call a%run(status)
    if(status/=0)stop 15
    call comm_child_counts(live,slots)
    if(nproc_id_global==0)write(*,*)'GS_RESOURCE_COUNTS',i,live,slots,slots_first
    if(live/=0.or.slots/=slots_first)stop 16
  enddo
  call a%clear()
  call a%run(status)
  if(status/=1)stop 9
  call end_parallel
  call MPI_Finalize(mpi_status)
contains
  subroutine save_status(filename)
    character(*),intent(in)::filename
    character(:),allocatable::bytes
    integer::length,source,dest
    if(nproc_id_global/=0)return
    inquire(file='gs_status.dat',size=length)
    allocate(character(length)::bytes)
    open(newunit=source,file='gs_status.dat',access='stream',form='unformatted',status='old')
    read(source)bytes
    close(source)
    open(newunit=dest,file=filename,access='stream',form='unformatted',status='new')
    write(dest)bytes
    close(dest)
  end subroutine
end program
