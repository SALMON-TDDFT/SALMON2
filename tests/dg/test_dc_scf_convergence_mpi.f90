#include "config.h"
program test_dc_scf_convergence_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dc_scf_convergence,only:reduce_dc_density_convergence
  implicit none
  integer::comm,rank,nproc,ierr,mode_index
  real(real64)::local_abs,local_square,global_abs,global_square,hvol,nelec,value,expected
  logical::ok
  character(256)::message
  character(12),parameter::modes(3)=[character(12)::'rho_dne','norm_rho','norm_rho_dng']

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  local_abs=real(rank+1,real64)
  local_square=real((rank+1)*(rank+1),real64)
  call MPI_Allreduce(local_abs,global_abs,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call MPI_Allreduce(local_square,global_square,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  hvol=0.25d0;nelec=2.5d0

  do mode_index=1,size(modes)
    call reduce_dc_density_convergence(comm,trim(modes(mode_index)),local_abs,local_square,&
      hvol,nelec,4*nproc,value,ok,message)
    call require(ok,trim(message))
    select case(trim(modes(mode_index)))
    case('rho_dne');expected=hvol*global_abs/nelec
    case('norm_rho');expected=global_square
    case('norm_rho_dng');expected=global_square/real(4*nproc,real64)
    end select
    call require(abs(value-expected)<=32d0*epsilon(1d0)*max(1d0,abs(expected)),&
      trim(modes(mode_index))//' does not match conventional DC')
  enddo

  call reduce_dc_density_convergence(comm,'maximum_density',local_abs,local_square,&
    hvol,nelec,4*nproc,value,ok,message)
  call require(.not.ok,'unsupported convergence mode was accepted')

  call reduce_dc_density_convergence(comm,'rho_dne',&
    transfer(int(z'7ff8000000000000',int64),0d0),local_square,&
    hvol,nelec,4*nproc,value,ok,message)
  call require(.not.ok,'non-finite local accumulator was accepted')

  call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
    0d0,nelec,4*nproc,value,ok,message)
  call require(.not.ok,'non-positive grid volume was accepted')

  call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
    hvol,0d0,4*nproc,value,ok,message)
  call require(.not.ok,'non-positive electron target was accepted')

  call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
    hvol,nelec,0,value,ok,message)
  call require(.not.ok,'non-positive global grid count was accepted')

  if(nproc>1)then
    hvol=merge(0.25d0,0.5d0,rank==0)
    call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
      hvol,nelec,4*nproc,value,ok,message)
    call require(.not.ok,'rank-disagreeing cell volumes were accepted')
    hvol=0.25d0

    nelec=merge(2.5d0,3d0,rank==0)
    call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
      hvol,nelec,4*nproc,value,ok,message)
    call require(.not.ok,'rank-disagreeing electron counts were accepted')
    nelec=2.5d0

    call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
      hvol,nelec,4*nproc+merge(0,1,rank==0),value,ok,message)
    call require(.not.ok,'rank-disagreeing grid counts were accepted')

    mode_index=merge(1,2,rank==0)
    call reduce_dc_density_convergence(comm,trim(modes(mode_index)),local_abs,local_square,&
      hvol,nelec,4*nproc,value,ok,message)
    call require(.not.ok,'rank-disagreeing convergence modes were accepted')
  endif

  if(rank==0)write(*,'(a,i0,a)')'PASS DC SCF convergence on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dc_scf_convergence_mpi
