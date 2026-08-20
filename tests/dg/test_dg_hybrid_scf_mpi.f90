#include "config.h"
program test_dg_hybrid_scf_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_scf,only:run_dg_hybrid_self_consistent_ground_state
  implicit none
  integer,parameter::nglobal=8
  integer::comm,rank,nproc,ierr,nlocal,p,position,iterations,reset_count
  integer(int64),allocatable::point_ids(:)
  real(real64),allocatable::density(:),converged_density(:),target(:),output_density(:)
  real(real64)::density_residual,energy_residual,eigensystem_residual,electron_count_defect,symmetry_defect
  integer(int64)::fingerprint,reference_fingerprint
  logical::ok,oscillatory
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(point_ids(nlocal),density(nlocal),target(nlocal),output_density(nlocal));position=0
  do p=nglobal,1,-1
    if(mod(p-1,nproc)/=rank)cycle
    position=position+1;point_ids(position)=p;target(position)=0.5d0+0.02d0*p
  enddo
  density=0.2d0;oscillatory=.false.;reset_count=0
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    40,1d-9,1d-10,1d-10,1d-10,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(maxval(abs(converged_density-target))<2d-8.and.iterations>1.and.iterations<=40,&
    'hybrid SCF did not converge to the nonlinear fixed point')
  oscillatory=.true.;density=0.2d0;reset_count=0
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    60,1d-8,1d-9,1d-10,1d-10,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok.and.reset_count>0,'oscillating hybrid SCF did not reset/reduce Pulay history')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_SCF ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid SCF on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine update_potential(input_density,callback_ok)
    real(real64),intent(in)::input_density(:);logical,intent(out)::callback_ok
    callback_ok=all(input_density>=0d0)
  end subroutine update_potential
  subroutine assemble_hamiltonian(iteration,callback_ok)
    integer,intent(in)::iteration;logical,intent(out)::callback_ok
    callback_ok=iteration>0
  end subroutine assemble_hamiltonian
  subroutine solve_occupied_states(iteration,total_energy,residual,electron_defect,group_defect,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(out)::total_energy,residual,electron_defect,group_defect
    logical,intent(out)::callback_ok
    total_energy=sum((density-target)**2);call MPI_Allreduce(MPI_IN_PLACE,total_energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    residual=min(1d-12,1d-13*iteration);electron_defect=0d0;group_defect=0d0;callback_ok=.true.
  end subroutine solve_occupied_states
  subroutine reconstruct_density(reconstructed,callback_ok)
    real(real64),intent(out)::reconstructed(:);logical,intent(out)::callback_ok
    if(oscillatory)then;reconstructed=target-1.2d0*(density-target)
    else;reconstructed=target+0.25d0*(density-target);endif
    output_density=reconstructed;callback_ok=.true.
  end subroutine reconstruct_density
  subroutine pulay_mix(iteration,input_density,reconstructed,reset_history,reduce_rate,mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(in)::input_density(:),reconstructed(:)
    logical,intent(in)::reset_history,reduce_rate
    real(real64),intent(out)::mixed_density(:);logical,intent(out)::callback_ok
    real(real64)::beta
    beta=0.7d0;if(reduce_rate)beta=0.2d0
    if(reset_history)reset_count=reset_count+1
    mixed_density=input_density+beta*(reconstructed-input_density);density=mixed_density;callback_ok=.true.
  end subroutine pulay_mix
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_scf_mpi
