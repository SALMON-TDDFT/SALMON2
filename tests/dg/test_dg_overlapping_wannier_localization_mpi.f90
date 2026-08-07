#include "config.h"
program test_dg_overlapping_wannier_localization_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_localization,only:evaluate_dg_periodic_localization
  implicit none
  complex(8)::values(2,4),phases(3,4),shifted_phases(3,4),moment(3,2),shifted_moment(3,2)
  real(8)::weights(4),norm(2),shifted_norm(2),spread,shifted_spread,nan_value
  logical::ok
  character(256)::message
  integer::ierr,rank,nproc,point

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  weights=1d0;values=(0d0,0d0);phases=(1d0,0d0)
  phases(1,:)=[(1d0,0d0),(0d0,1d0),(-1d0,0d0),(0d0,-1d0)]
  values(1,1)=1d0
  values(2,1)=1d0/sqrt(2d0);values(2,3)=1d0/sqrt(2d0)
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(ok,'valid periodic localization payload')
  call require(maxval(abs(norm-1d0))<1d-14,'Wannier norms')
  call require(abs(moment(1,1)-1d0)<1d-14,'delta-localized phase moment')
  call require(abs(moment(1,2))<1d-14,'opposite-site delocalized phase moment')
  call require(abs(spread-1d0)<1d-14,'bounded periodic spread distinguishes localization')

  do point=1,4
    shifted_phases(1,point)=phases(1,point)*exp(cmplx(0d0,0.37d0,8))
    shifted_phases(2,point)=phases(2,point)*exp(cmplx(0d0,-0.21d0,8))
    shifted_phases(3,point)=phases(3,point)*exp(cmplx(0d0,0.13d0,8))
  end do
  call evaluate_dg_periodic_localization(values,weights,shifted_phases,shifted_norm,&
    shifted_moment,shifted_spread,ok,message)
  call require(ok.and.abs(shifted_spread-spread)<1d-14,'periodic spread is origin invariant')

  nan_value=ieee_value(0d0,ieee_quiet_nan);weights(2)=nan_value
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'finite')>0,'nonfinite weights rejected')
  weights=1d0;values(1,:)=0d0
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'norm')>0,'zero-norm Wannier rejected')
  values=(0d0,0d0);values(1,1)=1d0;values(2,2)=1d0;phases(1,1)=2d0
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'phase')>0,'non-unit periodic phase rejected')

  if(rank==0)write(*,'(a,i0,a)')'PASS buffer-local periodic Wannier spread on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_failure/=0)error stop label
  end subroutine require
end program test_dg_overlapping_wannier_localization_mpi
