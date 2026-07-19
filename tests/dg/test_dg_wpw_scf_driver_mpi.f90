program test_dg_wpw_scf_driver_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_scf_driver,only:run_dg_wpw_scf_command_loop,broadcast_dg_wpw_occupied_coefficients
  use dg_wpw_scf_driver,only:s_dg_wpw_scf_iteration_state,initialize_dg_wpw_scf_iteration_state,&
    advance_dg_wpw_scf_iteration,evaluate_dg_wpw_scf_convergence,verify_dg_wpw_scf_fixed_point
  implicit none
  integer::ierr,rank,info,algebra_calls,map_calls,final_command
  complex(8)::cw(2,1),cp(1,1)
  type(s_dg_wpw_scf_iteration_state)::state
  type(s_dg_wpw_scf_iteration_state)::variable_state
  real(8),allocatable::variable_density(:)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  allocate(variable_density(merge(2,1,rank==0)));variable_density=0.25d0
  call initialize_dg_wpw_scf_iteration_state(variable_state,MPI_COMM_WORLD,0,2,1,1,variable_density,info)
  if(info/=0.or.size(variable_state%density)/=size(variable_density))call MPI_Abort(MPI_COMM_WORLD,34,ierr)
  algebra_calls=0;map_calls=0
  cw=(0d0,0d0);cp=(0d0,0d0)
  if(rank==0)then;cw(:,1)=[(1d0,0.2d0),(2d0,-0.1d0)];cp(1,1)=(3d0,0.4d0);endif
  call broadcast_dg_wpw_occupied_coefficients(MPI_COMM_WORLD,0,2,1,1,cw,cp,info)
  if(info/=0.or.maxval(abs(cw(:,1)-[(1d0,0.2d0),(2d0,-0.1d0)]))>1d-14.or.&
    abs(cp(1,1)-(3d0,0.4d0))>1d-14)call MPI_Abort(MPI_COMM_WORLD,9,ierr)
  if(rank==0)cw(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call broadcast_dg_wpw_occupied_coefficients(MPI_COMM_WORLD,0,2,1,1,cw,cp,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,8,ierr)
  cw=(0d0,0d0);cp=(0d0,0d0)
  call broadcast_dg_wpw_occupied_coefficients(MPI_COMM_WORLD,0,merge(2,1,rank==0),1,1,cw,cp,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,7,ierr)
  call initialize_dg_wpw_scf_iteration_state(state,MPI_COMM_WORLD,0,2,1,1,[0.2d0,0.4d0],info)
  if(info/=0)call MPI_Abort(MPI_COMM_WORLD,6,ierr)
  call advance_dg_wpw_scf_iteration(state,algebra_coefficients,potential_candidate,info)
  if(info/=0.or.state%iteration/=1.or.maxval(abs(state%density-[0.3d0,0.5d0]))>1d-14.or.&
    abs(state%energy-0.8d0)>1d-14.or.maxval(abs(state%q_old_occ(:,1)-[(1d0,0d0),(1d0,0d0),&
    (2d0,0d0)]))>1d-14.or.maxval(abs(state%previous_density-[0.2d0,0.4d0]))>1d-14.or.&
    abs(state%previous_energy)>1d-14.or.abs(state%gap-1.5d0)>1d-14.or.&
    abs(state%generalized_residual-1d-9)>1d-20.or.abs(state%metric_orthonormality-2d-10)>1d-20.or.&
    abs(state%projector_residual-3d-9)>1d-20.or.state%converged)call MPI_Abort(MPI_COMM_WORLD,5,ierr)
  call evaluate_dg_wpw_scf_convergence(state,1d0,1d-8,0d0,info)
  if(info/=0.or.state%converged)call MPI_Abort(MPI_COMM_WORLD,51,ierr)
  call advance_dg_wpw_scf_iteration(state,algebra_coefficients,potential_candidate_failure,info)
  if(info==0.or.state%iteration/=1.or.maxval(abs(state%density-[0.3d0,0.5d0]))>1d-14.or.&
    abs(state%energy-0.8d0)>1d-14.or.maxval(abs(state%cw-(1d0,0d0)))>1d-14.or.&
    maxval(abs(state%cp-(2d0,0d0)))>1d-14.or.maxval(abs(state%previous_density-[0.2d0,0.4d0]))>1d-14.or.&
    abs(state%previous_energy)>1d-14)call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  call advance_dg_wpw_scf_iteration(state,algebra_coefficients,potential_energy_mismatch,info)
  if(info==0.or.state%iteration/=1.or.abs(state%energy-0.8d0)>1d-14)call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  call advance_dg_wpw_scf_iteration(state,algebra_coefficients,potential_fixed_point,info)
  if(info/=0.or.state%iteration/=2)call MPI_Abort(MPI_COMM_WORLD,31,ierr)
  call verify_dg_wpw_scf_fixed_point(state,1d0,1d-8,potential_fixed_point,info)
  if(info/=0.or..not.state%converged)call MPI_Abort(MPI_COMM_WORLD,32,ierr)
  call verify_dg_wpw_scf_fixed_point(state,1d0,1d-8,potential_fixed_point_energy_failure,info)
  if(info==0.or..not.state%converged)call MPI_Abort(MPI_COMM_WORLD,321,ierr)
  call evaluate_dg_wpw_scf_convergence(state,1d0,1d-8,dble(rank)*1d-9,info)
  if(info==0.or..not.state%converged)call MPI_Abort(MPI_COMM_WORLD,33,ierr)
  call run_dg_wpw_scf_command_loop(MPI_COMM_WORLD,0,4,algebra_step,potential_step,info,final_command)
  if(info/=0.or.final_command/=3.or.algebra_calls/=merge(3,0,rank==0).or.map_calls/=3)call MPI_Abort(MPI_COMM_WORLD,10,ierr)
  algebra_calls=0;map_calls=0
  call run_dg_wpw_scf_command_loop(MPI_COMM_WORLD,0,4,algebra_step,potential_failure,info,final_command)
  if(info==0.or.final_command/=4.or.algebra_calls/=merge(2,0,rank==0).or.map_calls/=2)&
    call MPI_Abort(MPI_COMM_WORLD,11,ierr)
  if(rank==0)print '(a)','PASS all fragment ranks follow SCF potential commands and terminate on failure'
  call MPI_Finalize(ierr)
contains
  subroutine algebra_step(iter,converged,step_info)
    integer,intent(in)::iter;logical,intent(out)::converged;integer,intent(out)::step_info
    algebra_calls=algebra_calls+1;converged=iter>=3;step_info=0
  end subroutine
  subroutine potential_step(iter,step_info)
    integer,intent(in)::iter;integer,intent(out)::step_info
    map_calls=map_calls+1;step_info=0
  end subroutine
  subroutine potential_failure(iter,step_info)
    integer,intent(in)::iter;integer,intent(out)::step_info
    map_calls=map_calls+1;step_info=merge(1,0,rank==1.and.iter==2)
  end subroutine
  subroutine algebra_coefficients(iter,cw_out,cp_out,gap,residual,orth,projector,converged,step_info)
    integer,intent(in)::iter;complex(8),intent(out)::cw_out(:,:),cp_out(:,:)
    real(8),intent(out)::gap,residual,orth,projector
    logical,intent(out)::converged;integer,intent(out)::step_info
    cw_out=cmplx(dble(iter),0d0,8);cp_out=cmplx(2d0*dble(iter),0d0,8)
    gap=1.5d0;residual=1d-9;orth=2d-10;projector=3d-9
    converged=iter>=2;step_info=0
  end subroutine
  subroutine potential_candidate(iter,cw_in,cp_in,density_in,density_out,pot_res,energy,charge,step_info)
    integer,intent(in)::iter;complex(8),intent(in)::cw_in(:,:),cp_in(:,:)
    real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::step_info
    density_out=density_in+0.1d0;pot_res=0d0;energy=sum(density_out);charge=0d0;step_info=0
  end subroutine
  subroutine potential_candidate_failure(iter,cw_in,cp_in,density_in,density_out,pot_res,energy,charge,step_info)
    integer,intent(in)::iter;complex(8),intent(in)::cw_in(:,:),cp_in(:,:)
    real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::step_info
    density_out=density_in+10d0;pot_res=0d0;energy=sum(density_out);charge=0d0;step_info=merge(1,0,rank==1)
  end subroutine
  subroutine potential_energy_mismatch(iter,cw_in,cp_in,density_in,density_out,pot_res,energy,charge,step_info)
    integer,intent(in)::iter;complex(8),intent(in)::cw_in(:,:),cp_in(:,:)
    real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::step_info
    density_out=density_in;pot_res=0d0;energy=dble(rank);charge=0d0;step_info=0
  end subroutine
  subroutine potential_fixed_point(iter,cw_in,cp_in,density_in,density_out,pot_res,energy,charge,step_info)
    integer,intent(in)::iter;complex(8),intent(in)::cw_in(:,:),cp_in(:,:)
    real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::step_info
    density_out=density_in;pot_res=0d0;energy=0.8d0;charge=0d0;step_info=0
  end subroutine
  subroutine potential_fixed_point_energy_failure(iter,cw_in,cp_in,density_in,density_out,pot_res,energy,charge,step_info)
    integer,intent(in)::iter;complex(8),intent(in)::cw_in(:,:),cp_in(:,:)
    real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::step_info
    density_out=density_in;pot_res=0d0;energy=state%energy+1d-3;charge=0d0;step_info=0
  end subroutine
end program
