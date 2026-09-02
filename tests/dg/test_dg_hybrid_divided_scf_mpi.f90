#include "config.h"
program test_dg_hybrid_divided_scf_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_divided_scf,only:run_dg_hybrid_divided_scf
  implicit none
  integer,parameter::nglobal=8
  integer::comm,rank,nproc,ierr,nlocal,p,position,iterations,reference_iterations
  integer::callback_trace,reference_trace,potential_calls,solve_calls,density_calls,mix_calls
  integer(int64),allocatable::core_ids(:)
  real(real64),allocatable::core_weights(:),density(:),converged(:),reference(:),target(:),&
    last_potential_density(:)
  real(real64)::convergence_value,reference_value,electron_count,electron_defect
  real(real64)::callback_electron_offset,density_offset,consistent_density_offset,mixed_offset
  logical::ok,empty_pw_control,bad_density_nan,fail_terminal_potential
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(core_ids(nlocal),core_weights(nlocal),density(nlocal),target(nlocal),&
    last_potential_density(nlocal));position=0
  do p=1,nglobal
    if(mod(p-1,nproc)/=rank)cycle
    position=position+1;core_ids(position)=p;target(position)=0.1d0*p
  enddo
  core_weights=0.25d0

  call reset_case;density=0.45d0;empty_pw_control=.false.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'norm_rho_dng',1d-11,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,80,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(ok,trim(message));call require(maxval(abs(converged-target))<1d-10,'divided density mismatch')
  electron_count=sum(core_weights*converged)
  call MPI_Allreduce(MPI_IN_PLACE,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call require(abs(electron_count-0.9d0)<1d-10,'divided electron count mismatch')
  call require(electron_defect<=1d-12,'accepted divided density has an electron defect')
  call require(potential_calls==iterations+1,'terminal potential refresh was not called exactly once')
  call require(solve_calls==iterations.and.density_calls==iterations.and.mix_calls==iterations-1,&
    'successful divided SCF callback counts are inconsistent')
  call require(all(last_potential_density==converged),'terminal potential does not use converged density')
  call require(callback_trace==expected_success_trace(iterations),'successful callback order changed')
  reference=converged;reference_iterations=iterations;reference_value=convergence_value;reference_trace=callback_trace

  call reset_case;density=0.45d0;empty_pw_control=.true.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'norm_rho_dng',1d-11,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,80,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(ok.and.all(converged==reference),'empty-PW control changed the density path')
  call require(iterations==reference_iterations.and.convergence_value==reference_value,&
    'empty-PW control changed convergence inputs or iterations')
  call require(callback_trace==reference_trace,'empty-PW control changed callback ordering')

  call reset_case;density=0.45d0
  if(rank==0)density(1)=density(1)+0.5d0*1d-12/core_weights(1)
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',10d0,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(ok.and.allocated(converged),'within-tolerance initial electron defect was rejected')
  call require(potential_calls>0.and.solve_calls>0,'accepted initial density did not enter SCF')

  call reset_case;density=0.45d0
  if(rank==0)density(1)=density(1)+1.5d0*1d-12/core_weights(1)
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',10d0,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'invalid initial electron count was published')
  call require(potential_calls==0.and.solve_calls==0.and.density_calls==0.and.mix_calls==0,&
    'invalid initial electron count reached an SCF callback')
  call require(electron_defect>1d-12,'invalid initial electron defect was not reported')
  call require(convergence_value==huge(1d0),'convergence ran before initial electron gate')

  call reset_case;density=0.45d0;callback_electron_offset=1d-4
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'bad callback electron count was published')
  call require(mix_calls==0.and.electron_defect>1d-12,'bad callback electron count reached the mixer')
  call require(convergence_value==huge(1d0),'convergence ran before callback electron gate')

  call reset_case;density=0.45d0;density_offset=1d-4
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'independently invalid density was published')
  call require(mix_calls==0.and.electron_defect>1d-12,'independent electron gate ran after the mixer')
  call require(convergence_value==huge(1d0),'convergence ran before independent electron gate')

  call reset_case;density=0.45d0;consistent_density_offset=1d-4
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'target-inconsistent callback density was published')
  call require(mix_calls==0.and.electron_defect>1d-12,'target electron gate ran after the mixer')
  call require(convergence_value==huge(1d0),'convergence ran before target electron gate')

  call reset_case;density=0.45d0;bad_density_nan=.true.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'non-finite callback density was published')
  call require(mix_calls==0,'non-finite callback density reached the mixer')
  call require(convergence_value==huge(1d0),'convergence ran on a non-finite callback density')

  call reset_case;density=0.45d0;mixed_offset=1d-4
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'invalid mixed density was published')
  call require(mix_calls==1.and.potential_calls==1,'invalid mixed density reached another potential epoch')

  call reset_case;density=0.45d0;fail_terminal_potential=.true.
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'norm_rho_dng',1d-11,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,80,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok.and..not.allocated(converged),'failed terminal potential refresh was published')
  call require(potential_calls==iterations+1,'terminal potential failure was not observed at final epoch')
  call require(mix_calls==iterations-1,'terminal potential failure performed an extra density mix')

  if(nproc>1)then
    call reset_case;density=0.45d0
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,&
      merge(0.5d0,0.75d0,rank==0),'rho_dne',1d-8,&
      update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
      core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
    call require(.not.ok,'rank-disagreeing divided SCF controls were accepted')
    call require(callback_trace==0,'callbacks ran before divided SCF control agreement')

    call reset_case
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,&
      'rho_dne',1d-8,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,merge(0,2,rank==0),core_weights,0.9d0,1d-12,&
      converged,iterations,convergence_value,electron_defect,ok,message)
    call require(.not.ok,'rank-local invalid divided SCF control was accepted')
    call require(callback_trace==0,'callbacks ran after invalid divided SCF controls')

    call reset_case;density=0.45d0
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,&
      'rho_dne',1d-8,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,2,core_weights,merge(0.9d0,1d0,rank==0),1d-12,&
      converged,iterations,convergence_value,electron_defect,ok,message)
    call require(.not.ok,'rank-disagreeing expected electron count was accepted')
    call require(callback_trace==0,'callbacks ran before expected electron-count agreement')

    call reset_case
    call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,&
      'rho_dne',1d-8,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,2,core_weights,0.9d0,merge(1d-12,2d-12,rank==0),&
      converged,iterations,convergence_value,electron_defect,ok,message)
    call require(.not.ok,'rank-disagreeing electron tolerance was accepted')
    call require(callback_trace==0,'callbacks ran before electron-tolerance agreement')
  endif

  call reset_case
  if(rank==0)core_weights(1)=0d0
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok,'non-positive core weight was accepted')
  call require(callback_trace==0,'callbacks ran with an invalid core weight')
  if(rank==0)core_weights(1)=0.25d0

  ! Duplicate core ownership must be rejected collectively.
  if(nproc==1)then
    core_ids(2)=1_int64
  elseif(nlocal>0)then
    core_ids(1)=1_int64
  endif
  call reset_case
  call run_dg_hybrid_divided_scf(comm,nglobal,core_ids,density,0.5d0,'rho_dne',1d-8,&
    update_total_potential,solve_fragments,assemble_core_density,mix_dc_density,2,&
    core_weights,0.9d0,1d-12,converged,iterations,convergence_value,electron_defect,ok,message)
  call require(.not.ok,'duplicate core ownership was accepted')

  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_DIVIDED_SCF ranks=',nproc,' iterations=',reference_iterations
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid divided SCF on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine update_total_potential(input_density,callback_ok)
    real(real64),intent(in)::input_density(:);logical,intent(out)::callback_ok
    potential_calls=potential_calls+1;callback_trace=modulo(callback_trace*5+1,1000003)
    last_potential_density=input_density
    callback_ok=all(ieee_is_finite(input_density))
    if(fail_terminal_potential.and.rank==0.and.all(input_density==target))callback_ok=.false.
  end subroutine update_total_potential
  subroutine solve_fragments(iteration,callback_ok)
    integer,intent(in)::iteration;logical,intent(out)::callback_ok
    solve_calls=solve_calls+1;callback_trace=modulo(callback_trace*5+2,1000003);callback_ok=iteration>0
  end subroutine solve_fragments
  subroutine assemble_core_density(output_density,electron_count_value,callback_ok)
    real(real64),intent(out)::output_density(:),electron_count_value;logical,intent(out)::callback_ok
    density_calls=density_calls+1;callback_trace=modulo(callback_trace*5+3,1000003);output_density=target
    if(consistent_density_offset/=0d0.and.rank==0)&
      output_density(1)=output_density(1)+consistent_density_offset
    electron_count_value=sum(core_weights*output_density)
    call MPI_Allreduce(MPI_IN_PLACE,electron_count_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(density_offset/=0d0.and.rank==0)output_density(1)=output_density(1)+density_offset
    if(bad_density_nan.and.rank==0)&
      output_density(1)=transfer(int(z'7ff8000000000000',int64),0d0)
    if(rank==0)electron_count_value=electron_count_value+callback_electron_offset
    callback_ok=.true.
  end subroutine assemble_core_density
  subroutine mix_dc_density(iteration,input_density,output_density,mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(in)::input_density(:),output_density(:);real(real64),intent(out)::mixed_density(:)
    logical,intent(out)::callback_ok
    mix_calls=mix_calls+1;callback_trace=modulo(callback_trace*5+4,1000003)
    mixed_density=0.5d0*(input_density+output_density)
    if(mixed_offset/=0d0.and.rank==0)mixed_density(1)=mixed_density(1)+mixed_offset
    callback_ok=iteration>0.and.(empty_pw_control.or..not.empty_pw_control)
  end subroutine mix_dc_density
  subroutine reset_case
    callback_trace=0;potential_calls=0;solve_calls=0;density_calls=0;mix_calls=0
    callback_electron_offset=0d0;density_offset=0d0;consistent_density_offset=0d0;mixed_offset=0d0
    bad_density_nan=.false.;fail_terminal_potential=.false.;empty_pw_control=.false.
    last_potential_density=huge(1d0)
  end subroutine reset_case
  integer function expected_success_trace(iteration_count)result(trace)
    integer,intent(in)::iteration_count
    integer::iteration
    trace=0
    do iteration=1,iteration_count-1
      trace=modulo(trace*5+1,1000003);trace=modulo(trace*5+2,1000003)
      trace=modulo(trace*5+3,1000003);trace=modulo(trace*5+4,1000003)
    enddo
    trace=modulo(trace*5+1,1000003);trace=modulo(trace*5+2,1000003)
    trace=modulo(trace*5+3,1000003);trace=modulo(trace*5+1,1000003)
  end function expected_success_trace
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_divided_scf_mpi
