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
  character(32) :: mode
  character(2048) :: path,arg
  integer :: u,nk,nb,iter,nprocs,b,x,y,z,i,k,expected_nk
  logical :: real_orbital
  real(8),allocatable :: occ(:,:)
  real(8) :: r,v
  complex(8) :: expected
  call get_command_argument(1,mode);call get_command_argument(2,path)
  nk=1
  call get_command_argument(3,arg)
  if(len_trim(arg)>0)read(arg,*)nk
  if(nk<1)stop 5
  expected_nk=nk;info%ik_e=nk;system%nk=nk
  allocate(occ(2,nk))
  allocate(psi%zwf(8,8,8,1,2,nk,1))
  real_orbital=index(mode,'real')>0
  if(real_orbital)allocate(psi%rwf(8,8,8,1,2,nk,1))
  if(index(mode,'write')==1)then
    do k=1,nk
    do b=1,2
    do z=1,8
    do y=1,8
    do x=1,8
      i=x+8*(y-1)+64*(z-1)
      psi%zwf(x,y,z,1,b,k,1)=cmplx(dble(i+1000*b+10000*(k-1))/1024.d0,dble(i-1000*b-5000*(k-1))/2048.d0,8)
      if(index(mode,'norm')>0)psi%zwf(x,y,z,1,b,k,1)=normalized(b,x,real_orbital,k)
      if(real_orbital)psi%rwf(x,y,z,1,b,k,1)=real(psi%zwf(x,y,z,1,b,k,1),8)
    end do
    end do
    end do
    end do
    end do
    open(newunit=u,file=trim(path)//'info.bin',status='new',form='unformatted')
    write(u)nk
    write(u)2
    iter=0
    if(index(mode,'gs')>0)iter=1
    write(u)iter
    write(u)1
    write(u)real_orbital
    close(u)
    occ(1,:)=2.d0;occ(2,:)=0.d0
    if(index(mode,'norm')>0)occ=2.d0
    open(newunit=u,file=trim(path)//'occupation.bin',status='new',form='unformatted')
    write(u)occ
    close(u)
    call write_wavefunction(trim(path),lg,mg,system,info,psi,.false.)
  else if(index(mode,'read')==1)then
    open(newunit=u,file=trim(path)//'info.bin',status='old',form='unformatted')
    read(u)nk
    read(u)nb
    read(u)iter
    read(u)nprocs
    read(u)real_orbital
    close(u)
    if(nk/=expected_nk .or. nb/=2 .or. nprocs/=1)stop 1
    if(real_orbital.neqv.(index(mode,'real')>0))stop 1
    open(newunit=u,file=trim(path)//'occupation.bin',status='old',form='unformatted')
    read(u)occ
    close(u)
    if(any(occ(1,:)/=2.d0))stop 2
    if(index(mode,'weighted')>0)then
      if(nk/=2)stop 2
      if(occ(2,1)/=0.d0.or.occ(2,2)/=2.d0)stop 2
    else if(index(mode,'norm')>0)then
      if(any(occ(2,:)/=2.d0))stop 2
    else
      if(any(occ(2,:)/=0.d0))stop 2
    end if
    call read_wavefunction(trim(path),lg,mg,system,info,psi,nk,nb,real_orbital,.false.)
    do k=1,nk
    do b=1,2
    do z=1,8
    do y=1,8
    do x=1,8
      i=x+8*(y-1)+64*(z-1)
      expected=cmplx(dble(i+1000*b+10000*(k-1))/1024.d0,dble(i-1000*b-5000*(k-1))/2048.d0,8)
      if(index(mode,'norm')>0)expected=normalized(b,x,real_orbital.or.index(mode,'cos')>0,k)
      if(real_orbital)then
        if(psi%rwf(x,y,z,1,b,k,1)/=real(expected,8))stop 3
      else
        if(psi%zwf(x,y,z,1,b,k,1)/=expected)stop 3
      end if
    end do
    end do
    end do
    end do
    end do
  else
    stop 4
  end if
  print *, 'Unchanged native wavefunction '//trim(mode)//' probe PASS'
contains
  complex(8) function normalized(b,x,cosine,k) result(v)
    integer,intent(in) :: b,x,k
    logical,intent(in) :: cosine
    real(8) :: phase
    phase=2.d0*acos(-1.d0)*dble(x-1)/8.d0*dble(1+mod(k-1,3))
    if(b==1)then
      v=cmplx(1.d0/sqrt(512.d0),0.d0,8)
    else if(cosine)then
      v=cmplx(sqrt(2.d0/512.d0)*cos(phase),0.d0,8)
    else
      v=cmplx(cos(phase),sin(phase),8)/sqrt(512.d0)
    end if
  end function
end program
