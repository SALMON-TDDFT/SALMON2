#include "config.h"
program test_rt_dg_hybrid_stationarity_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use mpi
  use rt_dg_hybrid_stationarity,only:s_rt_dg_hybrid_stationarity_reference,&
    s_rt_dg_hybrid_stationarity_receipt,initialize_rt_dg_hybrid_stationarity,&
    evaluate_rt_dg_hybrid_stationarity
  use dg_hybrid_total_energy,only:evaluate_dg_hybrid_fixed_energy
  implicit none
  integer,parameter::certified_rank=3,noccupied=2,grid_count=2
  type(s_rt_dg_hybrid_stationarity_reference)::reference,rejected_reference
  type(s_rt_dg_hybrid_stationarity_receipt)::receipt
  complex(real64),allocatable::coefficients(:,:),s_coefficients(:,:),changed(:,:),&
    extra_coefficients(:,:),extra_s_coefficients(:,:)
  complex(real64)::energy_coefficients(2,2),kinetic(2,2),sipg(2,2),nonlocal(2,2)
  real(real64),allocatable::density(:),changed_density(:)
  real(real64)::occupations(noccupied),kinetic_energy,nonlocal_energy,rotation_scale
  integer(int64),allocatable::row_ids(:),extra_row_ids(:)
  integer::comm,rank,nproc,ierr,nowned,extra_nowned,row,position,extra_owner,npoint,point
  logical::ok,weighted_projector_ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  energy_coefficients=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[2,2])
  kinetic=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(2d0,0d0)],[2,2])
  sipg=reshape([(0.25d0,0d0),(0d0,0d0),(0d0,0d0),(0.5d0,0d0)],[2,2])
  nonlocal=reshape([(-0.1d0,0d0),(0d0,0d0),(0d0,0d0),(-0.2d0,0d0)],[2,2])
  call evaluate_dg_hybrid_fixed_energy(MPI_COMM_SELF,[1_int64,2_int64],energy_coefficients,[1d0,1d0],&
    kinetic,sipg,nonlocal,kinetic_energy,nonlocal_energy,ok,message)
  call require(ok,trim(message))
  call require(abs(kinetic_energy-3.75d0)<=1d-12,&
    'complete SIPG kinetic energy was not counted exactly once')
  call require(abs(nonlocal_energy+0.3d0)<=1d-12,'nonlocal projector energy mismatch')

  ! Exercise the production row-owned contract instead of replicating every
  ! coefficient row on every rank.
  nowned=count([(mod(row-1,nproc)==rank,row=1,certified_rank)])
  npoint=count([(mod(point-1,nproc)==rank,point=1,grid_count)])
  allocate(row_ids(nowned),coefficients(nowned,noccupied),s_coefficients(nowned,noccupied),&
    changed(nowned,noccupied),density(npoint))
  coefficients=(0d0,0d0);position=0
  rotation_scale=1d0/sqrt(2d0)
  do row=1,certified_rank
    if(mod(row-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=int(row,int64)
    select case(row)
    case(1,2)
      coefficients(position,1)=cmplx(rotation_scale,0d0,real64)
    case(3)
      coefficients(position,2)=(1d0,0d0)
    end select
  enddo
  s_coefficients=coefficients;occupations=[1d0,0.25d0]
  density=1d0/real(grid_count,real64)

  call initialize_rt_dg_hybrid_stationarity(comm,certified_rank-1,row_ids,density,-3d0,&
    coefficients,s_coefficients,occupations,1.25d0,0d0,rejected_reference,ok,message)
  call require(.not.ok,'certified-rank mismatch was accepted during stationarity initialization')

  ! An otherwise-zero construction-only direction must not enter the RT state.
  extra_owner=mod(certified_rank,nproc)
  extra_nowned=nowned+merge(1,0,rank==extra_owner)
  allocate(extra_row_ids(extra_nowned),extra_coefficients(extra_nowned,noccupied),&
    extra_s_coefficients(extra_nowned,noccupied))
  extra_row_ids(1:nowned)=row_ids
  extra_coefficients(1:nowned,:)=coefficients
  extra_s_coefficients(1:nowned,:)=s_coefficients
  if(rank==extra_owner)then
    extra_row_ids(extra_nowned)=int(certified_rank+1,int64)
    extra_coefficients(extra_nowned,:)=(0d0,0d0)
    extra_s_coefficients(extra_nowned,:)=(0d0,0d0)
  endif
  call initialize_rt_dg_hybrid_stationarity(comm,certified_rank,extra_row_ids,density,-3d0,&
    extra_coefficients,extra_s_coefficients,occupations,1.25d0,0d0,rejected_reference,ok,message)
  call require(.not.ok,&
    'construction-only coefficient direction was accepted by stationarity initialization')

  call initialize_rt_dg_hybrid_stationarity(comm,certified_rank,row_ids,density,-3d0,&
    coefficients,s_coefficients,occupations,1.25d0,0d0,reference,ok,message)
  call require(ok,trim(message))
  weighted_projector_ok=.true.
  do position=1,nowned
    if(row_ids(position)==3_int64)&
      weighted_projector_ok=abs(reference%projector(position,3)-(0.25d0,0d0))<=1d-12
  enddo
  call require(weighted_projector_ok,'fractional occupation was omitted from the stationarity projector')
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-3d0,&
    coefficients,s_coefficients,1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(ok.and.receipt%accepted,'zero-field certified-rank state was not stationary')
  call require(maxval([receipt%density_drift,receipt%energy_drift,receipt%projector_drift,&
    receipt%electron_drift,receipt%hamiltonian_residual])<=1d-12,&
    'zero-field density/energy/projector/charge/H receipt is nonzero')

  changed=coefficients
  changed(:,1)=changed(:,1)*exp(cmplx(0d0,0.37d0,real64))
  changed(:,2)=changed(:,2)*exp(cmplx(0d0,-0.21d0,real64))
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-3d0,changed,changed,&
    1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(ok.and.receipt%accepted.and.receipt%projector_drift<=1d-12,&
    'phase-equivalent occupied orbitals were not recognized as stationary')

  allocate(changed_density,source=density)
  if(rank==0.and.size(changed_density)>0)changed_density(1)=changed_density(1)+1d-3
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,changed_density,-3d0,&
    coefficients,s_coefficients,1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and.receipt%density_drift>1d-12,'independent density drift was accepted')
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-2.9d0,&
    coefficients,s_coefficients,1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and.receipt%energy_drift>1d-12,'independent energy drift was accepted')
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-3d0,&
    coefficients,s_coefficients,1.3d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and.receipt%electron_drift>1d-12,'independent electron drift was accepted')
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-3d0,&
    coefficients,s_coefficients,1.25d0,1d-3,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and.receipt%hamiltonian_residual>1d-12,&
    'independent Hamiltonian residual was accepted')

  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank+1,reference,density,-3d0,&
    coefficients,s_coefficients,1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and..not.receipt%accepted,&
    'stationarity evaluation accepted a certified-rank mismatch')

  ! This sign change preserves every local diagonal projector block but changes
  ! an off-rank block of the distributed occupied projector.
  changed=coefficients
  do position=1,nowned
    if(row_ids(position)==2_int64)changed(position,1)=-changed(position,1)
  enddo
  call evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,-3d0,changed,changed,&
    1.25d0,0d0,[1d-12,1d-12,1d-12,1d-12,1d-12],receipt,ok,message)
  call require(.not.ok.and..not.receipt%accepted,&
    'changed distributed occupied subspace was accepted as stationary')

  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid RT stationarity on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_rt_dg_hybrid_stationarity_mpi
