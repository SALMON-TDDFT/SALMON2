program test_dg_dc_ground_state_mpi
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use dg_nodal_state
  use dg_dc_ground_state
  implicit none

  type(s_dg_nodal_common_state) :: nodal
  type(s_dg_dc_gs_controls) :: controls
  type(s_dg_dc_gs_result) :: result
  real(8), allocatable :: density(:,:,:,:)
  integer(8) :: owner(2,1,1)
  integer :: ierr,rank,nproc,mode,forced_failure_count
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'
  owner(:,1,1)=[int(2*rank+1,8),int(2*rank+2,8)]
  call initialize_dg_nodal_common_state(nodal,rank+1,[2,1,1],[4,3,3],[1,1,1],1,2,1,owner, &
    101_8,202_8,MPI_COMM_WORLD,ok,message)
  call require(ok,'nodal initialization')
  nodal%psi_core=(0d0,0d0)
  nodal%psi_core(1,1,1,1,1)=1d0
  nodal%psi_core(2,1,1,2,1)=1d0
  nodal%occupations(:,1)=[2d0,0d0]
  allocate(density(2,1,1,1)); density=0.5d0

  controls=default_dg_dc_gs_controls()
  controls%initial_lambda_step=0.25d0
  controls%minimum_lambda_step=0.0625d0
  controls%maximum_lambda_step=0.5d0
  controls%maximum_scf_iterations=50
  controls%maximum_rollbacks=4
  controls%allowed_residual_growth=1.5d0
  controls%intermediate_orbital_tolerance=1d-5
  controls%intermediate_density_tolerance=1d-5
  controls%final_orbital_tolerance=1d-7
  controls%final_density_tolerance=1d-7
  controls%subspace_tolerance=1d-7
  controls%hermiticity_tolerance=1d-12
  controls%orthogonality_tolerance=1d-12
  controls%face_balance_tolerance=1d-12
  controls%electron_count_tolerance=1d-12

  mode=0
  forced_failure_count=0
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  if(rank==0 .and. .not.ok) write(*,'(a,1x,a,3(a,i0),a,f8.4)') 'baseline failure:',trim(message), &
    ' phase=',result%phase,' accepted=',result%naccepted_steps,' rollback=',result%nrollbacks,' lambda=',result%lambda
  call require(ok .and. result%accepted,'known nonlinear fixed point acceptance')
  call require(result%lambda==1d0 .and. result%unmixed_verified,'final lambda and unmixed fixed point')
  call require(result%lambda_history(1)==0d0,'lambda zero handoff')
  call require(result%naccepted_steps>=3,'intermediate convergence before lambda advance')
  call require(maxval(result%lambda_steps)>controls%initial_lambda_step,'adaptive step growth')
  call require(result%maximum_interface_scale==1d0,'complete Hermitian interface scaling')
  call require(result%minimum_projector_overlap>0.99d0,'occupied projector tracking not state labels')
  call require(result%mixing_reset_count==result%naccepted_steps,'density mixing reset per lambda')

  call reset_fixture(nodal,density)
  mode=1
  forced_failure_count=0
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(ok .and. result%nrollbacks>=1,'rollback and step halving recovery')

  call reset_fixture(nodal,density)
  mode=2
  forced_failure_count=0
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(.not.ok .and. result%failed,'minimum-step failure')

  call reset_fixture(nodal,density)
  mode=3
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(.not.ok,'unmixed fixed-point rejection')

  call reset_fixture(nodal,density)
  mode=4
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(.not.ok,'stale fingerprint collective failure')

  call reset_fixture(nodal,density)
  mode=5
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(.not.ok,'nonfinite collective failure')

  call reset_fixture(nodal,density)
  mode=6
  call run_dg_dc_ground_state(nodal,density,controls,toy_step,MPI_COMM_WORLD,result,ok,message)
  call require(.not.ok,'rank-local callback failure propagates collectively')

  if(rank==0) print '(a)','PASS adaptive self-consistent DG continuation'
  call MPI_Finalize(ierr)

contains

  subroutine toy_step(state,density_arg,lambda,density_mix,reset_mixing,unmixed,diagnostics,communicator, &
                      step_ok,step_message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density_arg(:,:,:,:)
    real(8), intent(in) :: lambda,density_mix
    logical, intent(in) :: reset_mixing,unmixed
    type(s_dg_dc_gs_diagnostics), intent(out) :: diagnostics
    integer, intent(in) :: communicator
    logical, intent(out) :: step_ok
    character(*), intent(out) :: step_message
    real(8) :: target,current,new_density,residual
    complex(8) :: orbital_swap(2,1,1)

    target=0.5d0+0.25d0*lambda+0.1d0*density_arg(1,1,1,1)**2
    current=density_arg(1,1,1,1)
    new_density=target
    if(.not.unmixed) new_density=(1d0-density_mix)*current+density_mix*target
    density_arg=new_density
    residual=abs(new_density-current)
    diagnostics%orbital_residual=residual
    diagnostics%density_residual=residual
    diagnostics%subspace_residual=0.25d0*residual
    diagnostics%projector_overlap=0d0
    diagnostics%hermiticity_defect=0d0
    diagnostics%orthogonality_defect=0d0
    diagnostics%face_balance_defect=0d0
    diagnostics%electron_number=2d0
    diagnostics%expected_electron_number=2d0
    diagnostics%interface_scale=lambda
    diagnostics%eigensolver_iterations=2
    diagnostics%eigensolver_converged=.true.
    diagnostics%finite=.true.
    if(lambda>0.3d0) then
      if(state%occupations(1,1)>0d0) then
        orbital_swap=state%psi_core(:,:,:,1,1)
        state%psi_core(:,:,:,1,1)=state%psi_core(:,:,:,2,1)
        state%psi_core(:,:,:,2,1)=orbital_swap
        state%occupations(:,1)=[0d0,2d0]
      end if
    end if
    if(mode==1 .and. lambda>=0.5d0 .and. forced_failure_count<controls%maximum_scf_iterations) then
      diagnostics%orbital_residual=1d0
      forced_failure_count=forced_failure_count+1
    end if
    if(mode==2 .and. lambda>0d0) diagnostics%orbital_residual=1d0
    if(mode==3 .and. unmixed) diagnostics%density_residual=1d-2
    if(mode==4 .and. rank==1 .and. lambda>0d0) state%operator_fingerprint=999_8
    if(mode==5 .and. rank==1 .and. lambda>0d0) &
      diagnostics%orbital_residual=ieee_value(0d0,ieee_quiet_nan)
    if(mode==6 .and. rank==1 .and. lambda>0d0) step_ok=.false.
    if(.not.(mode==6 .and. rank==1 .and. lambda>0d0)) step_ok=.true.
    step_message=''
  end subroutine toy_step

  subroutine reset_fixture(state,rho)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: rho(:,:,:,:)
    state%operator_fingerprint=202_8
    state%psi_core=(0d0,0d0)
    state%psi_core(1,1,1,1,1)=1d0
    state%psi_core(2,1,1,2,1)=1d0
    state%occupations(:,1)=[2d0,0d0]
    rho=0.5d0
  end subroutine reset_fixture

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: global
    integer :: mpi_error
    call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if(.not.global) then
      if(rank==0) write(*,'(a,1x,a)') 'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    end if
  end subroutine require
end program test_dg_dc_ground_state_mpi
