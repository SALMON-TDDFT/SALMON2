program test_dg_dc_ground_state_adapter_mpi
  use mpi
  use dg_nodal_state
  use dg_dc_ground_state_adapter
  use dg_dc_ground_state, only: s_dg_dc_gs_diagnostics
  implicit none
  type(s_dg_nodal_common_state) :: state
  type(s_dg_dc_gs_diagnostics) :: production_diagnostics
  integer(8) :: owner(2,1,1)
  complex(8) :: hvolume(2,1,1,3,1),hface(2,1,1,3,1),hpsi(2,1,1,3,1)
  complex(8), allocatable :: hface_global(:,:,:,:,:)
  complex(8) :: constant_psi(3,1,1,1,1),zero_local(3,1,1,1,1),interior_volume(3,1,1,1,1)
  complex(8) :: probe_x(3,1,1),probe_y(3,1,1),probe_hx(3,1,1),probe_hy(3,1,1)
  real(8) :: local_potential(3,1,1,1),laplacian_coefficients(2,3)
  real(8) :: raw_density(2,1,1,1),accepted_density(2,1,1,1),mixed_density(2,1,1,1)
  real(8) :: accepted_potential(2,1,1,1),trial_potential(2,1,1,1)
  complex(8) :: lift_core(3,2,1,1,1),test_core(3,2,1,1,1)
  complex(8) :: face_value(2,1,1,1),face_normal(2,1,1,1),discrete_pairing,lift_pairing
  complex(8) :: neighbor_value(2,1,1,1),neighbor_normal(2,1,1,1)
  complex(8) :: sipg_value(2,1,1,1),sipg_normal(2,1,1,1)
  complex(8) :: exchanged_value(2,1,1,1),exchanged_normal(2,1,1,1)
  real(8) :: derivative_weights(3)
  real(8) :: fd_weights(5),xcoord
  integer :: degree,node
  integer :: production_sequence
  integer :: fragment_origins(3,2),fragment_sizes(3,2),neighbor_fragment
  real(8) :: face_hermiticity,face_balance
  real(8) :: measured_balance
  real(8) :: estimated_bytes
  integer :: ierr,rank,nproc
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'
  owner(:,1,1)=[int(2*rank+1,8),int(2*rank+2,8)]
  call initialize_dg_nodal_common_state(state,rank+1,[2,1,1],[4,3,3],[1,1,1],1,3,1,owner, &
    101_8,202_8,MPI_COMM_WORLD,ok,message)
  call require(ok,'state initialization')
  state%psi_core=(0d0,0d0)
  state%psi_core(1,1,1,1,1)=1d0
  state%psi_core(2,1,1,2,1)=1d0
  state%psi_core(:,:,:,3,1)=100d0 ! MPI-shape padding: must never contribute.
  state%occupations(:,1)=[2d0,1d0,0d0]

  hvolume=(2d0,0d0)
  hface=(4d0,0d0)
  call compose_dg_dc_hamiltonian(hvolume,hface,0.25d0,hpsi,ok,message)
  call require(ok .and. all(abs(hpsi-(3d0,0d0))<1d-14),'only SIPG action is lambda scaled')

  call reconstruct_dg_dc_core_density(state,2,raw_density,ok,message)
  call require(ok,'valid-candidate density reconstruction')
  call require(abs(sum(raw_density)-3d0)<1d-14,'padding candidate excluded from density')

  accepted_density=1d0
  accepted_potential=7d0
  trial_potential=9d0
  call mix_dg_dc_density_transaction(accepted_density,raw_density,0.25d0,.false.,mixed_density,ok,message)
  call require(ok .and. maxval(abs(mixed_density-(0.75d0*accepted_density+0.25d0*raw_density)))<1d-14, &
    'simple transactional density mix')
  call rollback_dg_dc_density_potential(accepted_density,accepted_potential,mixed_density,trial_potential)
  call require(all(mixed_density==accepted_density) .and. all(trial_potential==accepted_potential), &
    'density and potential rollback together')

  lift_core=(0d0,0d0)
  test_core=(0d0,0d0)
  test_core(:,1,1,1,1)=[(1d0,0.5d0),(2d0,-0.25d0),(3d0,0.75d0)]
  test_core(:,2,1,1,1)=[(-1d0,0.25d0),(0.5d0,0.5d0),(2d0,-1d0)]
  derivative_weights=[-1.5d0,2d0,-0.5d0]
  face_value(:,1,1,1)=[(0.25d0,0.5d0),(-0.75d0,0.125d0)]
  face_normal(:,1,1,1)=[(0.5d0,-0.25d0),(0.125d0,0.75d0)]
  call lift_dg_dc_sipg_face(1,-1,derivative_weights,face_value,face_normal,lift_core,ok,message)
  call require(ok,'adjoint SIPG face lift')
  discrete_pairing=conjg(test_core(1,1,1,1,1))*face_value(1,1,1,1)+ &
    conjg(sum(derivative_weights*test_core(:,1,1,1,1)))*face_normal(1,1,1,1)+ &
    conjg(test_core(1,2,1,1,1))*face_value(2,1,1,1)+ &
    conjg(sum(derivative_weights*test_core(:,2,1,1,1)))*face_normal(2,1,1,1)
  lift_pairing=sum(conjg(test_core)*lift_core)
  call require(abs(discrete_pairing-lift_pairing)<1d-13,'face lift is exact derivative adjoint')

  neighbor_value=-face_value
  neighbor_normal=0.5d0*face_normal
  call evaluate_dg_dc_local_sipg_face(face_value,face_normal,neighbor_value,neighbor_normal, &
    0.5d0,0.25d0,8d0,1,local_face_value=sipg_value,local_face_normal=sipg_normal,ok=ok,message=message)
  call require(ok .and. all(abs(sipg_value)>0d0),'local SIPG face evaluation')
  face_value=cmplx(dble(rank+1),0d0,8)
  face_normal=cmplx(dble(10*(rank+1)),0d0,8)
  call exchange_dg_dc_face_traces_mpi(face_value,face_normal,1-rank,41,MPI_COMM_WORLD, &
    exchanged_value,exchanged_normal,ok,message)
  call require(ok .and. all(real(exchanged_value,8)==dble(2-rank)) .and. &
    all(real(exchanged_normal,8)==dble(10*(2-rank))),'paired face trace exchange')
  fragment_origins=reshape([0,0,0,3,0,0],shape(fragment_origins))
  fragment_sizes=reshape([3,2,1,3,2,1],shape(fragment_sizes))
  call find_dg_dc_periodic_face_neighbor(rank+1,1,1,fragment_origins,fragment_sizes,[6,2,1], &
    neighbor_fragment,ok,message)
  call require(ok .and. neighbor_fragment==2-rank,'periodic physical face neighbor')

  call build_dg_dc_one_sided_derivative_weights(0.25d0,-1,fd_weights,ok,message)
  call require(ok,'one-sided derivative weights')
  do degree=0,4
    xcoord=0d0
    do node=1,5
      xcoord=xcoord+fd_weights(node)*(0.25d0*dble(node-1))**degree
    end do
    call require(abs(xcoord-merge(-1d0,0d0,degree==1))<1d-11,'outward one-sided polynomial exactness')
  end do

  call expand_dg_dc_global_candidate_axis(state,2,MPI_COMM_WORLD,ok,message)
  call require(ok .and. state%nstate==4,'global state axis contains every fragment candidate')
  if(rank==0) then
    call require(any(abs(state%psi_core(:,:,:,1:2,1))>0d0) .and. &
      all(abs(state%psi_core(:,:,:,3:4,1))==0d0),'rank zero owns only its global candidate columns')
  else
    call require(all(abs(state%psi_core(:,:,:,1:2,1))==0d0) .and. &
      any(abs(state%psi_core(:,:,:,3:4,1))>0d0),'rank one owns only its global candidate columns')
  end if
  fragment_origins=reshape([0,0,0,2,0,0],shape(fragment_origins))
  fragment_sizes=reshape([2,1,1,2,1,1],shape(fragment_sizes))
  state%psi_core=(0d0,0d0)
  state%psi_core(:,:,:,1,1)=(1d0,0d0)
  allocate(hface_global,mold=state%psi_core)
  if(rank==0) state%fragment=2
  call initialize_dg_dc_physical_faces(state,fragment_origins,fragment_sizes,[4,1,1],MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-local topology mismatch fails collectively')
  state%fragment=rank+1
  call initialize_dg_dc_physical_faces(state,fragment_origins,fragment_sizes,[4,1,1],MPI_COMM_WORLD,ok,message)
  call require(ok,'physical face topology initialization')
  call apply_dg_dc_sipg_operator_mpi(state,fragment_origins,fragment_sizes,[4,1,1], &
    [0.5d0,1d0,1d0],8d0,MPI_COMM_WORLD,hface_global,face_hermiticity,face_balance,ok,message)
  call require(ok .and. maxval(abs(hface_global))<1d-13,'constant-state SIPG null action')
  call measure_dg_dc_face_action_balance(sipg_value,sipg_normal,-sipg_value,-sipg_normal,measured_balance,ok)
  call require(ok .and. measured_balance<1d-14,'paired value and normal action balance')
  call measure_dg_dc_face_action_balance(sipg_value,sipg_normal,-sipg_value,-sipg_normal+(1d-3,0d0), &
    measured_balance,ok)
  call require(ok .and. measured_balance>1d-4,'normal-action imbalance is detected')
  constant_psi=(1d0,0d0)
  local_potential=0d0
  laplacian_coefficients=0d0
  laplacian_coefficients(1,1)=1d0
  call build_dg_dc_interior_volume_action(constant_psi,local_potential,-2d0,laplacian_coefficients, &
    zero_local,interior_volume,ok,message)
  call require(ok .and. maxval(abs(interior_volume))<1d-13,'constant-state interior kinetic null action')
  call require(abs(zero_local(1,1,1,1,1))>0d0 .and. abs(zero_local(3,1,1,1,1))>0d0, &
    'zero-extension boundary artifact is isolated for subtraction')
  probe_x=(0d0,0d0); probe_y=(0d0,0d0); probe_hx=(0d0,0d0); probe_hy=(0d0,0d0)
  probe_x(1,1,1)=(1d0,0d0)
  probe_y(2,1,1)=(1d0,0d0)
  probe_hy(1,1,1)=(1d0,0d0)
  call measure_dg_dc_cross_hermiticity_mpi(probe_x,probe_hx,probe_y,probe_hy,MPI_COMM_WORLD, &
    measured_balance,ok)
  call require(ok .and. measured_balance>0.5d0,'cross-bilinear non-Hermiticity is detected')
  call validate_dg_dc_candidate_memory([100,100,100],10000,1,1024d0*1024d0,estimated_bytes,ok)
  call require(.not.ok .and. estimated_bytes>1024d0*1024d0,'bounded candidate workspace rejection')
  production_sequence=0
  call execute_dg_dc_production_iteration(state,raw_density,0.5d0,0.25d0,.true.,.false., &
    mock_restore,mock_solve,mock_update,MPI_COMM_WORLD,production_diagnostics,ok,message)
  call require(ok .and. production_sequence==123,'production restore-solve-update call order')

  if(rank==0) print '(a)','PASS DG-DC production adapter'
  call release_dg_nodal_common_state(state)
  call MPI_Finalize(ierr)

contains
  subroutine mock_restore(state_arg,density_arg,communicator,callback_ok,callback_message)
    type(s_dg_nodal_common_state), intent(inout) :: state_arg
    real(8), intent(inout) :: density_arg(:,:,:,:)
    integer, intent(in) :: communicator
    logical, intent(out) :: callback_ok
    character(*), intent(out) :: callback_message
    production_sequence=10*production_sequence+1
    callback_ok=communicator==MPI_COMM_WORLD
    callback_message=''
  end subroutine mock_restore

  subroutine mock_solve(state_arg,lambda,communicator,diagnostics,callback_ok,callback_message)
    type(s_dg_nodal_common_state), intent(inout) :: state_arg
    real(8), intent(in) :: lambda
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
    logical, intent(out) :: callback_ok
    character(*), intent(out) :: callback_message
    production_sequence=10*production_sequence+2
    diagnostics%orbital_residual=0d0
    diagnostics%eigensolver_converged=.true.
    callback_ok=lambda==0.5d0 .and. communicator==MPI_COMM_WORLD
    callback_message=''
  end subroutine mock_solve

  subroutine mock_update(state_arg,density_arg,density_mix,unmixed,communicator,diagnostics, &
      callback_ok,callback_message)
    type(s_dg_nodal_common_state), intent(inout) :: state_arg
    real(8), intent(inout) :: density_arg(:,:,:,:)
    real(8), intent(in) :: density_mix
    logical, intent(in) :: unmixed
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
    logical, intent(out) :: callback_ok
    character(*), intent(out) :: callback_message
    production_sequence=10*production_sequence+3
    diagnostics%density_residual=0d0
    callback_ok=density_mix==0.25d0 .and. .not.unmixed .and. communicator==MPI_COMM_WORLD
    callback_message=''
  end subroutine mock_update

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
end program test_dg_dc_ground_state_adapter_mpi
