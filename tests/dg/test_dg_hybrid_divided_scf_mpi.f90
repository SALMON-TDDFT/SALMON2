#include "config.h"
program test_dg_hybrid_divided_scf_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_divided_scf,only:run_dg_hybrid_divided_scf
  implicit none
  integer,parameter::nglobal=8
  integer::comm,rank,nproc,ierr,nlocal,p,position,iterations,reference_iterations
  integer::callback_trace,reference_trace
  integer(int64),allocatable::core_ids(:)
  real(real64),allocatable::density(:),converged(:),reference(:),target(:)
  real(real64)::convergence_value,reference_value,electron_count
  logical::ok,empty_pw_control
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(core_ids(nlocal),density(nlocal),target(nlocal));position=0
  do p=1,nglobal
    if(mod(p-1,nproc)/=rank)cycle
    position=position+1;core_ids(position)=p;target(position)=0.1d0*p
  enddo

  density=0d0;callback_trace=0;empty_pw_control=.false.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,1d0,3.6d0,'norm_rho_dng',1d-11,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,80,&
    converged,iterations,convergence_value,ok,message)
  call require(ok,trim(message));call require(maxval(abs(converged-target))<1d-10,'divided density mismatch')
  electron_count=sum(converged);call MPI_Allreduce(MPI_IN_PLACE,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call require(abs(electron_count-3.6d0)<1d-10,'divided electron count mismatch')
  reference=converged;reference_iterations=iterations;reference_value=convergence_value;reference_trace=callback_trace

  density=0d0;callback_trace=0;empty_pw_control=.true.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,1d0,3.6d0,'norm_rho_dng',1d-11,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,80,&
    converged,iterations,convergence_value,ok,message)
  call require(ok.and.all(converged==reference),'empty-PW control changed the density path')
  call require(iterations==reference_iterations.and.convergence_value==reference_value,&
    'empty-PW control changed convergence inputs or iterations')
  call require(callback_trace==reference_trace,'empty-PW control changed callback ordering')

  if(nproc>1)then
    callback_trace=0
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,&
      merge(1d0,2d0,rank==0),3.6d0,'rho_dne',1d-8,&
      update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
      converged,iterations,convergence_value,ok,message)
    call require(.not.ok,'rank-disagreeing divided SCF controls were accepted')
    call require(callback_trace==0,'callbacks ran before divided SCF control agreement')

    callback_trace=0
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,1d0,3.6d0,&
      'rho_dne',1d-8,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,merge(0,2,rank==0),converged,iterations,convergence_value,ok,message)
    call require(.not.ok,'rank-local invalid divided SCF control was accepted')
    call require(callback_trace==0,'callbacks ran after invalid divided SCF controls')
  endif

  ! Duplicate core ownership must be rejected collectively.
  if(nproc==1)then
    core_ids(2)=1_int64
  elseif(nlocal>0)then
    core_ids(1)=1_int64
  endif
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,1d0,3.6d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    converged,iterations,convergence_value,ok,message)
  call require(.not.ok,'duplicate core ownership was accepted')

  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_DIVIDED_SCF ranks=',nproc,' iterations=',reference_iterations
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid divided SCF on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine update_total_potential(input_density,callback_ok)
    real(real64),intent(in)::input_density(:);logical,intent(out)::callback_ok
    callback_trace=modulo(callback_trace*5+1,1000003);callback_ok=all(input_density>=0d0)
  end subroutine update_total_potential
  subroutine solve_fragments(iteration,callback_ok)
    integer,intent(in)::iteration;logical,intent(out)::callback_ok
    callback_trace=modulo(callback_trace*5+2,1000003);callback_ok=iteration>0
  end subroutine solve_fragments
  subroutine assemble_core_density(output_density,electron_count_value,callback_ok)
    real(real64),intent(out)::output_density(:),electron_count_value;logical,intent(out)::callback_ok
    callback_trace=modulo(callback_trace*5+3,1000003);output_density=target
    electron_count_value=sum(output_density)
    call MPI_Allreduce(MPI_IN_PLACE,electron_count_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    callback_ok=.true.
  end subroutine assemble_core_density
  subroutine mix_dc_density(iteration,input_density,output_density,mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(in)::input_density(:),output_density(:);real(real64),intent(out)::mixed_density(:)
    logical,intent(out)::callback_ok
    callback_trace=modulo(callback_trace*5+4,1000003);mixed_density=0.5d0*(input_density+output_density)
    callback_ok=iteration>0.and.(empty_pw_control.or..not.empty_pw_control)
  end subroutine mix_dc_density
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_divided_scf_mpi
