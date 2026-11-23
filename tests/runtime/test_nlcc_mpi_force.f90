#include "config.h"
! Actual MPI spatial reduction with an additional replicated orbital group.
program test_nlcc_mpi_force
  use mpi
  use communication, only: comm_summation
  use structures, only: s_pp_info,s_dft_system,s_rgrid,s_scalar
  use salmon_pp, only: calc_nlcc_force
  use salmon_global, only: kion,iperiodic
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter :: n=6,nr=500,na=2
  type(s_pp_info) :: pp
  type(s_dft_system) :: sys
  type(s_rgrid) :: full,local
  type(s_scalar) :: vxc(1)
  real(8) :: reference(3,na),partial(3,na),reduced(3,na),radius,error
  integer :: rank,np,ierr,comm_r,i,j,k,status,global_status,rejected,all_rejected,format_,formats
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  if(np/=2 .and. np/=4)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  ! Two spatial ranks in each orbital replica, just as icomm_r fixes KO.
  call MPI_Comm_split(MPI_COMM_WORLD,rank/2,rank,comm_r,ierr)
  pp%flag_nlcc=.true.;pp%nrmax=nr
  allocate(pp%rad(nr,1),pp%rho_nlcc_tbl(nr,1),pp%tau_nlcc_tbl(nr,1))
  do i=1,nr
    radius=(i-1)*.012d0
    pp%rad(i,1)=radius;pp%rho_nlcc_tbl(i,1)=exp(-(radius/.65d0)**2);pp%tau_nlcc_tbl(i,1)=0d0
  end do
  allocate(kion(na));kion=1;iperiodic=3
  sys%nion=na;sys%Nspin=1;sys%hgs=.71d0;sys%Hvol=.71d0**3
  sys%rmatrix_a=0d0
  do i=1,3
    sys%rmatrix_a(i,i)=1d0
  end do
  sys%rmatrix_a(2,1)=.19d0
  sys%primitive_a=sys%rmatrix_a*(n*.71d0)
  allocate(sys%Rion(3,na))
  sys%Rion(:,1)=[.314d0,.217d0,.173d0];sys%Rion(:,2)=[1.213d0,.873d0,2.197d0]
  full%is=1;full%ie=n;full%num=n;local=full
  local%is(1)=1+mod(rank,2)*(n/2);local%ie(1)=local%is(1)+n/2-1;local%num=local%ie-local%is+1
  allocate(vxc(1)%f(n,n,n))
  do k=1,n
  do j=1,n
  do i=1,n
    vxc(1)%f(i,j,k)=.3d0+.07d0*sin(2*acos(-1d0)*(i-1)/n)+.05d0*cos(2*acos(-1d0)*(j+k-2)/n)
  end do
  end do
  end do
  formats=1
#ifdef SALMONPY_FAST_CORE
  allocate(sys%kion(na));sys%kion=1;pp%input_config%iperiodic=iperiodic
  allocate(pp%full_nlcc_support(1),pp%mr(1));pp%full_nlcc_support=.false.;pp%mr=nr-2
  formats=2
#endif
  do format_=1,formats
#ifdef SALMONPY_FAST_CORE
    if(format_==2)then
      pp%full_nlcc_support=.true.
      do i=1,nr
        radius=(i-1)*.04d0;pp%rad(i,1)=radius
        pp%rho_nlcc_tbl(i,1)=exp(-(radius/3.8d0)**2)
      enddo
    endif
#endif
  call calc_nlcc_force(pp,sys,full,vxc,reference,status)
  if(status/=0)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  call calc_nlcc_force(pp,sys,local,vxc,partial,status)
  call comm_summation(status,global_status,MPI_COMM_WORLD)
  if(global_status/=0)call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  call comm_summation(partial,reduced,3*na,comm_r)
  error=maxval(abs(reference-reduced))
  if(rank==0)write(*,*)'MPI_NLCC_FORCE_COMPARE',format_,np,error
  if(error>1d-13)call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  ! A single bad spatial/KO rank must make every replica reject collectively.
  if(rank==0)vxc(1)%f(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
  call calc_nlcc_force(pp,sys,local,vxc,partial,status)
  call comm_summation(status,global_status,MPI_COMM_WORLD)
  rejected=merge(1,0,global_status/=0)
  call comm_summation(rejected,all_rejected,MPI_COMM_WORLD)
  if(all_rejected/=np)call MPI_Abort(MPI_COMM_WORLD,5,ierr)
  if(rank==0)write(*,*)'MPI_NLCC_FORCE_PASS',np,error,all_rejected
  if(rank==0)vxc(1)%f(1,1,1)=.35d0
  enddo
  call MPI_Comm_free(comm_r,ierr)
  call MPI_Finalize(ierr)
end program
