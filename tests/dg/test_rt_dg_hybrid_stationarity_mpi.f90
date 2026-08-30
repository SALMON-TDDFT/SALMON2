#include "config.h"
program test_rt_dg_hybrid_stationarity_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use mpi
  use rt_dg_hybrid_stationarity,only:s_rt_dg_hybrid_stationarity_reference,&
    s_rt_dg_hybrid_stationarity_receipt,initialize_rt_dg_hybrid_stationarity,&
    evaluate_rt_dg_hybrid_stationarity
  use dg_hybrid_total_energy,only:evaluate_dg_hybrid_fixed_energy
  implicit none
  type(s_rt_dg_hybrid_stationarity_reference)::reference
  type(s_rt_dg_hybrid_stationarity_receipt)::receipt
  complex(real64)::coefficients(2,2),s_coefficients(2,2),rotated(2,2)
  complex(real64)::kinetic(2,2),sipg(2,2),nonlocal(2,2)
  real(real64)::density(2),occupations(2),rotation_scale,kinetic_energy,nonlocal_energy
  integer(int64)::row_ids(2)
  integer::comm,rank,nproc,ierr
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  row_ids=[1_int64,2_int64];density=[0.75d0,1.25d0];occupations=[1d0,1d0]
  coefficients=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[2,2])
  s_coefficients=coefficients
  kinetic=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(2d0,0d0)],[2,2])
  sipg=reshape([(0.25d0,0d0),(0d0,0d0),(0d0,0d0),(0.5d0,0d0)],[2,2])
  nonlocal=reshape([(-0.1d0,0d0),(0d0,0d0),(0d0,0d0),(-0.2d0,0d0)],[2,2])
  call evaluate_dg_hybrid_fixed_energy(MPI_COMM_SELF,row_ids,coefficients,occupations,kinetic,sipg,&
    nonlocal,kinetic_energy,nonlocal_energy,ok,message)
  if(.not.ok)error stop trim(message)
  if(abs(kinetic_energy-3.75d0)>1d-12)error stop 'complete SIPG kinetic energy was not counted exactly once'
  if(abs(nonlocal_energy+0.3d0)>1d-12)error stop 'nonlocal projector energy mismatch'
  call initialize_rt_dg_hybrid_stationarity(comm,row_ids,density,-3d0,coefficients,&
    s_coefficients,occupations,2d0,0d0,reference,ok,message)
  if(.not.ok)error stop trim(message)

  rotation_scale=1d0/sqrt(2d0)
  rotated(:,1)=rotation_scale*(coefficients(:,1)+coefficients(:,2))
  rotated(:,2)=rotation_scale*(-coefficients(:,1)+coefficients(:,2))
  call evaluate_rt_dg_hybrid_stationarity(comm,reference,density,-3d0,rotated,rotated,&
    2d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  if(.not.ok.or..not.receipt%accepted)error stop 'degenerate occupied rotation was not stationary'
  if(maxval([receipt%density_drift,receipt%energy_drift,receipt%projector_drift,&
      receipt%electron_drift,receipt%hamiltonian_residual])>1d-12)&
    error stop 'stationary receipt is nonzero'

  rotated(:,1)=coefficients(:,1);rotated(:,2)=(0d0,0d0)
  call evaluate_rt_dg_hybrid_stationarity(comm,reference,density,-3d0,rotated,rotated,&
    1d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  if(ok.or.receipt%accepted)error stop 'changed occupied subspace was accepted'

  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid RT stationarity on ',nproc,' ranks'
  call MPI_Finalize(ierr)
end program test_rt_dg_hybrid_stationarity_mpi
