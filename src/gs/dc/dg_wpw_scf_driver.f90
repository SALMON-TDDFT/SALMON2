module dg_wpw_scf_driver
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter,public::DG_WPW_SCF_ALGEBRA=1,DG_WPW_SCF_POTENTIAL=2
  integer,parameter,public::DG_WPW_SCF_CONVERGED=3,DG_WPW_SCF_FAILED=4
  type,public::s_dg_wpw_scf_iteration_state
    integer::comm=MPI_COMM_NULL,root=0,nw=0,np=0,nocc=0,iteration=0
    complex(8),allocatable::cw(:,:),cp(:,:),q_old_occ(:,:)
    real(8),allocatable::density(:),previous_density(:)
    real(8)::energy=0d0,previous_energy=0d0
    real(8)::gap=0d0,generalized_residual=huge(1d0),metric_orthonormality=huge(1d0),&
      projector_residual=huge(1d0)
    real(8)::density_residual=huge(1d0),potential_residual=huge(1d0),&
      energy_residual=huge(1d0),charge_error=huge(1d0),fixed_point_residual=huge(1d0)
    logical::algebra_ready=.false.,converged=.false.,valid=.false.
  end type
  abstract interface
    subroutine dg_wpw_algebra_step(iteration,converged,info)
      integer,intent(in)::iteration;logical,intent(out)::converged;integer,intent(out)::info
    end subroutine
    subroutine dg_wpw_potential_step(iteration,info)
      integer,intent(in)::iteration;integer,intent(out)::info
    end subroutine
    subroutine dg_wpw_coefficient_algebra(iteration,cw,cp,gap,residual,orth,projector,converged,info)
      integer,intent(in)::iteration;complex(8),intent(out)::cw(:,:),cp(:,:)
      real(8),intent(out)::gap,residual,orth,projector
      logical,intent(out)::converged;integer,intent(out)::info
    end subroutine
    subroutine dg_wpw_candidate_potential(iteration,cw,cp,density_in,density_out,potential_residual,&
        energy,charge_error,info)
      integer,intent(in)::iteration;complex(8),intent(in)::cw(:,:),cp(:,:)
      real(8),intent(in)::density_in(:);real(8),intent(out)::density_out(:),potential_residual,energy,charge_error
      integer,intent(out)::info
    end subroutine
  end interface
  public::run_dg_wpw_scf_command_loop
  public::broadcast_dg_wpw_occupied_coefficients
  public::initialize_dg_wpw_scf_iteration_state,advance_dg_wpw_scf_iteration
  public::evaluate_dg_wpw_scf_convergence
  public::verify_dg_wpw_scf_fixed_point
contains
  subroutine initialize_dg_wpw_scf_iteration_state(state,comm,root,nw,np,nocc,density,info)
    type(s_dg_wpw_scf_iteration_state),intent(out)::state
    integer,intent(in)::comm,root,nw,np,nocc;real(8),intent(in)::density(:);integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,dims(3),dmin(3),dmax(3)
    info=1;call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    dims=[nw,np,nocc]
    call MPI_Allreduce(dims,dmin,3,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(dims,dmax,3,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,root>=0.and.root<nrank.and.nw>=0.and.np>=0.and.nocc>0.and.size(density)>0.and.&
      all(dmin==dmax).and.all(ieee_is_finite(density)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    state%comm=comm;state%root=root;state%nw=nw;state%np=np;state%nocc=nocc
    allocate(state%cw(nw,nocc),state%cp(np,nocc),state%q_old_occ(nw+np,nocc))
    state%cw=0;state%cp=0;state%q_old_occ=0
    state%density=density;state%previous_density=density
    state%energy=0d0;state%previous_energy=0d0;state%iteration=0;state%valid=.true.;info=0
  end subroutine

  subroutine advance_dg_wpw_scf_iteration(state,algebra,potential,info)
    type(s_dg_wpw_scf_iteration_state),intent(inout)::state
    procedure(dg_wpw_coefficient_algebra)::algebra
    procedure(dg_wpw_candidate_potential)::potential
    integer,intent(out)::info
    complex(8),allocatable::cw_candidate(:,:),cp_candidate(:,:)
    real(8),allocatable::density_candidate(:)
    real(8)::energy_candidate,energy_min,energy_max,metrics(4),potential_candidate,charge_candidate
    real(8)::local_norms(2),global_norms(2),scalars(3),scalar_min(3),scalar_max(3)
    logical::converged_candidate
    integer::rank,ierr,step_info,local_bad,global_bad,next_iteration
    info=1;if(.not.state%valid)return
    call MPI_Comm_rank(state%comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    next_iteration=state%iteration+1
    allocate(cw_candidate(state%nw,state%nocc),cp_candidate(state%np,state%nocc),&
      density_candidate(size(state%density)))
    cw_candidate=0;cp_candidate=0;step_info=0;converged_candidate=.false.
    metrics=0d0
    if(rank==state%root)call algebra(next_iteration,cw_candidate,cp_candidate,metrics(1),metrics(2),&
      metrics(3),metrics(4),converged_candidate,step_info)
    if(rank==state%root.and.step_info==0)then
      if(.not.all(ieee_is_finite(metrics)).or.any(metrics<0d0))step_info=1
    endif
    local_bad=merge(1,0,rank==state%root.and.step_info/=0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Bcast(converged_candidate,1,MPI_LOGICAL,state%root,state%comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(metrics,4,MPI_DOUBLE_PRECISION,state%root,state%comm,ierr);if(ierr/=MPI_SUCCESS)return
    call broadcast_dg_wpw_occupied_coefficients(state%comm,state%root,state%nw,state%np,state%nocc,&
      cw_candidate,cp_candidate,step_info)
    if(step_info/=0)return
    call potential(next_iteration,cw_candidate,cp_candidate,state%density,density_candidate,potential_candidate,&
      energy_candidate,charge_candidate,step_info)
    local_bad=merge(0,1,step_info==0.and.all(ieee_is_finite(density_candidate)).and.&
      all(ieee_is_finite([potential_candidate,energy_candidate,charge_candidate])))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    scalars=[potential_candidate,energy_candidate,charge_candidate]
    call MPI_Allreduce(scalars,scalar_min,3,MPI_DOUBLE_PRECISION,MPI_MIN,state%comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(scalars,scalar_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,state%comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(abs(scalar_max-scalar_min)>100d0*epsilon(1d0)*max(1d0,abs(scalar_max))))return
    local_norms=[sum((density_candidate-state%density)**2),sum(density_candidate**2)]
    call MPI_Allreduce(local_norms,global_norms,2,MPI_DOUBLE_PRECISION,MPI_SUM,state%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    state%previous_density=state%density;state%previous_energy=state%energy
    state%cw=cw_candidate;state%cp=cp_candidate
    state%q_old_occ(1:state%nw,:)=cw_candidate
    state%q_old_occ(state%nw+1:state%nw+state%np,:)=cp_candidate
    state%density=density_candidate;state%energy=energy_candidate
    state%density_residual=sqrt(max(0d0,global_norms(1)))/max(1d-30,sqrt(max(0d0,global_norms(2))))
    state%potential_residual=potential_candidate;state%energy_residual=abs(energy_candidate-state%previous_energy)
    state%charge_error=charge_candidate
    state%gap=metrics(1);state%generalized_residual=metrics(2)
    state%metric_orthonormality=metrics(3);state%projector_residual=metrics(4)
    state%iteration=next_iteration;state%algebra_ready=converged_candidate
    state%converged=.false.;info=0
  end subroutine

  subroutine evaluate_dg_wpw_scf_convergence(state,gap_threshold,tolerance,fixed_point_residual,info)
    type(s_dg_wpw_scf_iteration_state),intent(inout)::state
    real(8),intent(in)::gap_threshold,tolerance,fixed_point_residual
    integer,intent(out)::info
    integer::local_bad,global_bad,ierr
    real(8)::values(3),values_min(3),values_max(3)
    info=1
    local_bad=merge(0,1,state%valid.and.gap_threshold>0d0.and.tolerance>0d0.and.&
      ieee_is_finite(fixed_point_residual).and.fixed_point_residual>=0d0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    values=[gap_threshold,tolerance,fixed_point_residual]
    call MPI_Allreduce(values,values_min,3,MPI_DOUBLE_PRECISION,MPI_MIN,state%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(values,values_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(values_min/=values_max))return
    state%fixed_point_residual=fixed_point_residual
    state%converged=state%iteration>1.and.state%algebra_ready.and.state%gap>=gap_threshold.and.&
      max(state%density_residual,state%potential_residual)<tolerance.and.state%energy_residual<tolerance.and.&
      state%projector_residual<tolerance.and.state%generalized_residual<tolerance.and.&
      state%metric_orthonormality<tolerance.and.abs(state%charge_error)<tolerance.and.&
      state%fixed_point_residual<tolerance
    info=0
  end subroutine

  subroutine verify_dg_wpw_scf_fixed_point(state,gap_threshold,tolerance,potential,info)
    type(s_dg_wpw_scf_iteration_state),intent(inout)::state
    real(8),intent(in)::gap_threshold,tolerance
    procedure(dg_wpw_candidate_potential)::potential
    integer,intent(out)::info
    type(s_dg_wpw_scf_iteration_state)::candidate
    real(8),allocatable::density_check(:)
    real(8)::potential_residual,energy_check,charge_error,local_norms(2),global_norms(2),fixed_residual
    real(8)::scalars(3),scalar_min(3),scalar_max(3)
    integer::step_info,local_bad,global_bad,ierr
    info=1
    if(.not.state%valid)return
    allocate(density_check(size(state%density)))
    call potential(state%iteration+1,state%cw,state%cp,state%density,density_check,potential_residual,&
      energy_check,charge_error,step_info)
    local_bad=merge(0,1,step_info==0.and.all(ieee_is_finite(density_check)).and.&
      all(ieee_is_finite([potential_residual,energy_check,charge_error])))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    scalars=[potential_residual,energy_check,charge_error]
    call MPI_Allreduce(scalars,scalar_min,3,MPI_DOUBLE_PRECISION,MPI_MIN,state%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(scalars,scalar_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(abs(scalar_max-scalar_min)>100d0*epsilon(1d0)*max(1d0,abs(scalar_max))))return
    local_norms=[sum((density_check-state%density)**2),sum(state%density**2)]
    call MPI_Allreduce(local_norms,global_norms,2,MPI_DOUBLE_PRECISION,MPI_SUM,state%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    fixed_residual=max(sqrt(max(0d0,global_norms(1)))/max(1d-30,sqrt(max(0d0,global_norms(2)))),&
      abs(potential_residual),abs(energy_check-state%energy),abs(charge_error))
    candidate=state
    call evaluate_dg_wpw_scf_convergence(candidate,gap_threshold,tolerance,fixed_residual,step_info)
    if(step_info/=0.or..not.candidate%converged)return
    state=candidate;info=0
  end subroutine
  subroutine broadcast_dg_wpw_occupied_coefficients(comm,root,nw,np,nocc,cw,cp,info)
    integer,intent(in)::comm,root,nw,np,nocc
    complex(8),intent(inout)::cw(nw,nocc),cp(np,nocc)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,dims(3),dims_min(3),dims_max(3)
    info=1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    dims=[nw,np,nocc]
    call MPI_Allreduce(dims,dims_min,3,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(dims,dims_max,3,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,root>=0.and.root<nrank.and.nw>=0.and.np>=0.and.nocc>0.and.&
      all(dims_min==dims_max))
    if(rank==root.and.local_bad==0)then
      if(.not.finite_complex(cw).or..not.finite_complex(cp))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Bcast(cw,nw*nocc,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(cp,np*nocc,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.finite_complex(cw).or..not.finite_complex(cp))return
    info=0
  end subroutine

  logical function finite_complex(a)result(ok)
    complex(8),intent(in)::a(:,:)
    ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
  subroutine run_dg_wpw_scf_command_loop(comm,root,max_iter,algebra,potential,info,final_command)
    integer,intent(in)::comm,root,max_iter
    procedure(dg_wpw_algebra_step)::algebra
    procedure(dg_wpw_potential_step)::potential
    integer,intent(out)::info
    integer,intent(out),optional::final_command
    integer::rank,nrank,ierr,iteration,command,step_info,local_bad,global_bad
    logical::converged
    info=1;if(present(final_command))final_command=DG_WPW_SCF_FAILED
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS.or.root<0.or.root>=nrank.or.max_iter<1)return
    do iteration=1,max_iter
      command=DG_WPW_SCF_ALGEBRA
      call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)return
      step_info=0;converged=.false.
      if(rank==root)call algebra(iteration,converged,step_info)
      local_bad=merge(1,0,rank==root.and.step_info/=0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Bcast(converged,1,MPI_LOGICAL,root,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then
        command=DG_WPW_SCF_FAILED;call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr)
        if(present(final_command))final_command=command
        info=10;return
      endif
      command=DG_WPW_SCF_POTENTIAL
      call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call potential(iteration,step_info)
      local_bad=merge(1,0,step_info/=0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then
        command=DG_WPW_SCF_FAILED;call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr)
        if(present(final_command))final_command=command
        info=20;return
      endif
      if(converged)then
        command=DG_WPW_SCF_CONVERGED
        call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr)
        if(present(final_command))final_command=command
        if(ierr==MPI_SUCCESS)info=0
        return
      endif
    enddo
    command=DG_WPW_SCF_FAILED
    call MPI_Bcast(command,1,MPI_INTEGER,root,comm,ierr)
    if(present(final_command))final_command=command
    if(ierr==MPI_SUCCESS)info=30
  end subroutine
end module
