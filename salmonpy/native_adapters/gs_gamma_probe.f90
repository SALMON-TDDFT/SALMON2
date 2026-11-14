! Copyright 2026 SALMON developers. Apache-2.0.
! Minimal scaffolding for compiling UNCHANGED external native wfn readers/writers.
! This is a synthetic I/O probe, not a SALMON Hamiltonian or a solver context.
module structures
  implicit none
  type s_rgrid
    integer :: is(3)=1,ie(3)=8
  end type
  type s_dft_system
    integer :: nspin=1,nk=1,no=2
  end type
  type s_parallel_info
    integer :: io_s=1,io_e=2,ik_s=1,ik_e=1,im_s=1,im_e=1
  end type
  type s_orbital
    real(8),allocatable :: rwf(:,:,:,:,:,:,:)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:)
  end type
contains
  subroutine allocate_orbital_real()
  end subroutine
  subroutine deallocate_orbital()
  end subroutine
end module
module communication
contains
  logical function comm_is_root(i)
    integer :: i
    comm_is_root=i==0
  end function
  subroutine comm_summation()
  end subroutine
  subroutine comm_bcast()
  end subroutine
end module
module filesystem
contains
  subroutine create_directory(name)
    character(*),intent(in) :: name
    error stop 'sliced file layout outside synthetic probe scope'
  end subroutine
end module
module salmon_global
  character(8) :: method_wf_distributor='single'
  integer :: nblock_wf_distribute=1
end module

program probe
  use structures
  implicit none
  type(s_rgrid) :: lg,mg
  type(s_dft_system) :: system
  type(s_parallel_info) :: info
  type(s_orbital) :: psi
  character(16) :: mode
  character(2048) :: path
  integer :: u,nk,nb,iter,nprocs,b,x,y,z,i
  logical :: real_orbital
  real(8) :: occ(2),r,v
  complex(8) :: expected
  call get_command_argument(1,mode);call get_command_argument(2,path)
  allocate(psi%zwf(8,8,8,1,2,1,1))
  if(trim(mode)=='write')then
    do b=1,2
    do z=1,8
    do y=1,8
    do x=1,8
      i=x+8*(y-1)+64*(z-1)
      psi%zwf(x,y,z,1,b,1,1)=cmplx(dble(i+1000*b)/1024.d0,dble(i-1000*b)/2048.d0,8)
    end do
    end do
    end do
    end do
    real_orbital=.false.
    open(newunit=u,file=trim(path)//'info.bin',status='new',form='unformatted')
    write(u)1
    write(u)2
    write(u)0
    write(u)1
    write(u)real_orbital
    close(u)
    occ=[2.d0,0.d0]
    open(newunit=u,file=trim(path)//'occupation.bin',status='new',form='unformatted')
    write(u)occ
    close(u)
    call write_wavefunction(trim(path),lg,mg,system,info,psi,.false.)
  else if(trim(mode)=='read')then
    open(newunit=u,file=trim(path)//'info.bin',status='old',form='unformatted')
    read(u)nk
    read(u)nb
    read(u)iter
    read(u)nprocs
    read(u)real_orbital
    close(u)
    if(nk/=1 .or. nb/=2 .or. iter/=0 .or. nprocs/=1 .or. real_orbital)stop 1
    open(newunit=u,file=trim(path)//'occupation.bin',status='old',form='unformatted')
    read(u)occ
    close(u)
    if(occ(1)/=2.d0 .or. occ(2)/=0.d0)stop 2
    call read_wavefunction(trim(path),lg,mg,system,info,psi,nk,nb,real_orbital,.false.)
    do b=1,2
    do z=1,8
    do y=1,8
    do x=1,8
      i=x+8*(y-1)+64*(z-1)
      expected=cmplx(dble(i+1000*b)/1024.d0,dble(i-1000*b)/2048.d0,8)
      if(psi%zwf(x,y,z,1,b,1,1)/=expected)stop 3
    end do
    end do
    end do
    end do
  else
    stop 4
  end if
  print *, 'Unchanged native wavefunction '//trim(mode)//' probe PASS'
end program
