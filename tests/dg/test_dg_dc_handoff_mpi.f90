program test_dg_dc_handoff_mpi
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use dg_nodal_state, only: s_dg_nodal_common_state
  use dg_dc_handoff
  implicit none

  type(s_dg_dc_handoff_state) :: handoff, snapshot
  type(s_dg_nodal_common_state) :: nodal
  complex(8) :: candidates(5,3,3,2,1), updated(3,1,1,2,1)
  integer(8) :: owner(3,1,1),snapshot_owner(3)
  real(8) :: density(3),potential(3,1)
  real(8) :: fragment_orbitals(5,3,3,1,2,1,1),restored_orbitals(5,3,3,1,2,1,1)
  integer :: ierr, rank, nproc
  logical :: ok, accept,advance,rollback,complete
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if (nproc /= 2) error stop 'test requires two MPI ranks'

  call initialize_dg_dc_handoff(handoff,.true.,3,1.0d-3,1,1.0d-8,ok,message)
  call require_collective(ok,'valid handoff controls')
  call evaluate_dg_dc_handoff(handoff,2,5.0d-4,4.0d0,4.0d0,.true.,.false.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'minimum iteration gate')
  call evaluate_dg_dc_handoff(handoff,3,2.0d-3,4.0d0,4.0d0,.true.,.false.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'residual gate')
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4,4.1d0,4.0d0,.true.,.false.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'charge gate')
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4,4.0d0,4.0d0,rank==0,.false.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'collective fragment-solve gate')
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4,4.0d0,4.0d0,.true.,.true.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'pre-handoff LCFO rejection')
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4,4.0d0,4.0d0,.true.,.false.,.true.,.false., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'pre-handoff Wannier rejection')
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4,4.0d0,4.0d0,.true.,.false.,.false.,.true., &
                              MPI_COMM_WORLD,accept,message)
  call require_collective(.not.accept,'pre-handoff window truncation rejection')

  call check_threshold(1.0d-2)
  call check_threshold(1.0d-3)
  call check_threshold(1.0d-4)
  call initialize_dg_dc_handoff(handoff,.true.,3,1.0d-2,1,1.0d-8,ok,message)
  call evaluate_dg_dc_handoff(handoff,3,5.0d-4+dble(rank)*1.0d-5,4.0d0,4.0d0,.true., &
                              .false.,.false.,.false.,MPI_COMM_WORLD,accept,message)
  call require_collective(accept,'collective maximum residual below threshold')
  call require_collective(abs(handoff%residual-5.1d-4)<1d-14,'canonical maximum residual')

  owner(:,1,1)=[int(3*rank+1,8),int(3*rank+2,8),int(3*rank+3,8)]
  candidates=(0.0d0,0.0d0)
  candidates(2,2,2,1,1)=1.0d0
  candidates(3,2,2,2,1)=1.0d0
  call materialize_dg_dc_candidates(handoff,nodal,candidates,rank+1,[3,1,1],[1,1,1],owner, &
                                    2,2,101_8,202_8,MPI_COMM_WORLD,ok,message)
  call require_collective(ok,'complete candidate materialization')
  call require_collective(handoff%materialized .and. handoff%candidate_count==2,'candidate completeness diagnostic')
  call require_collective(nodal%fragment==rank+1,'physical fragment identity')
  call require_collective(handoff%metric_rank==2 .and. handoff%minimum_metric_pivot>0d0, &
                          'metric-rank diagnostics')
  call require_collective(all(nodal%box_size==[5,3,3]) .and. all(nodal%buffer==[1,1,1]), &
                          'local-periodic core plus buffer metadata')
  call require_collective(all(nodal%core_owner==owner),'explicit unique core ownership')

  snapshot=handoff
  candidates(:,:,:,2,1)=candidates(:,:,:,1,1)
  call materialize_dg_dc_candidates(handoff,nodal,candidates,rank+1,[3,1,1],[1,1,1],owner, &
                                    2,2,303_8,404_8,MPI_COMM_WORLD,ok,message)
  call require_collective(.not.ok,'rank-deficient metric rejection')
  call require_collective(handoff%geometry_fingerprint==snapshot%geometry_fingerprint, &
                          'materialization failure is transactional')
  call require_collective(all(nodal%psi_core==snapshot%dc_reference_core), &
                          'nodal materialization failure is transactional')

  density=[1d0,2d0,1d0]; potential(:,1)=[-1d0,-2d0,-1d0]
  snapshot_owner=[int(3*rank+1,8),int(3*rank+2,8),int(3*rank+3,8)]
  call preserve_dg_dc_density_potential(handoff,density,potential,snapshot_owner,MPI_COMM_WORLD,ok,message)
  call require_collective(ok .and. handoff%density_potential_preserved .and. &
                          all(handoff%density_snapshot==density) .and. &
                          all(handoff%potential_snapshot==potential) .and. &
                          all(handoff%snapshot_owner==snapshot_owner),'distributed density and potential preservation')
  if(rank==1) snapshot_owner(1)=3_8
  call preserve_dg_dc_density_potential(handoff,density,potential,snapshot_owner,MPI_COMM_WORLD,ok,message)
  call require_collective(.not.ok,'snapshot duplicate ownership and matching hole rejection')
  call mark_dg_dc_mixing_discarded(handoff)
  call require_collective(handoff%mixing_history_discarded,'mixing history discard transition')
  fragment_orbitals=dble(rank+1)
  fragment_orbitals(:,:,:,:,2,:,:)=2d0*dble(rank+1)
  call preserve_dg_dc_fragment_orbitals(handoff,fragment_orbitals,MPI_COMM_WORLD,ok,message)
  fragment_orbitals=-1d0
  call restore_dg_dc_fragment_orbitals(handoff,restored_orbitals,MPI_COMM_WORLD,ok,message)
  call require_collective(ok .and. all(restored_orbitals(:,:,:,:,1,:,:)==dble(rank+1)) .and. &
                          all(restored_orbitals(:,:,:,:,2,:,:)==2d0*dble(rank+1)), &
                          'rank-local fragment orbital rollback snapshot')

  updated=nodal%psi_core
  updated(3,1,1,1,1)=updated(3,1,1,1,1)+cmplx(0.25d0,0.0d0,8)
  call update_dg_dc_nodal_state(handoff,nodal,updated,MPI_COMM_WORLD,ok,message)
  call require_collective(ok .and. handoff%span_escape_norm>0.0d0, &
                          'post-handoff update may leave original DC coefficient span')

  call initialize_dg_dc_handoff(handoff,.true.,1,1d-2,1,1d-8,ok,message)
  call evaluate_dg_dc_handoff(handoff,1,1d-3,4d0,4d0,.true.,.false.,.false.,.false., &
                              MPI_COMM_WORLD,accept,message)
  candidates=(0d0,0d0)
  candidates(2,2,2,1,1)=1d0
  if(rank==1) candidates(3,2,2,2,1)=1d0
  call materialize_dg_dc_candidates(handoff,nodal,candidates,rank+1,[3,1,1],[1,1,1],owner, &
                                    rank+1,rank+1,505_8,606_8,MPI_COMM_WORLD,ok,message)
  call require_collective(ok .and. nodal%nstate==2 .and. handoff%candidate_count==rank+1, &
                          'variable local candidate counts use a common padded state axis')

  call initialize_dg_dc_direct_continuation(handoff,0.25d0,0.125d0,0.5d0,4d0,3,ok,message)
  call require_collective(ok .and. abs(handoff%interface_scale-0.25d0)<1d-14, &
                          'direct continuation starts from a bounded nonzero scale')
  call evaluate_dg_dc_direct_continuation(handoff,1d-6,1d-5,.false.,MPI_COMM_WORLD, &
                                          advance,rollback,complete,message)
  call require_collective(advance .and. .not.rollback .and. .not.complete .and. &
                          abs(handoff%accepted_interface_scale-0.25d0)<1d-14 .and. &
                          handoff%interface_scale>handoff%accepted_interface_scale .and. &
                          all(handoff%continuation_scale_history==[0d0,0.25d0]) .and. &
                          all(handoff%continuation_step_history==[0d0,0.25d0]), &
                          'converged DC stage advances the DG scale')
  call evaluate_dg_dc_direct_continuation(handoff,1d2,1d-5,.false.,MPI_COMM_WORLD, &
                                          advance,rollback,complete,message)
  call require_collective(rollback .and. .not.advance .and. &
                          handoff%interface_scale==handoff%accepted_interface_scale .and. &
                          abs(handoff%accepted_stage_residual-1d-6)<1d-14 .and. &
                          abs(handoff%trial_stage_first_residual-1d2)<1d-12, &
                          'abnormal first residual rolls back against the accepted baseline')
  call initialize_dg_dc_direct_continuation(handoff,1d0,0.5d0,1d0,4d0,3,ok,message)
  call evaluate_dg_dc_direct_continuation(handoff,1d-6,1d-5,rank==0,MPI_COMM_WORLD, &
                                          advance,rollback,complete,message)
  call require_collective(.not.complete .and. advance .and. .not.rollback .and. &
                          .not.handoff%direct_ground_state_complete, &
                          'full-scale completion requires collective unmixed agreement')

  if (rank==0) print '(a)','PASS early DC-to-nodal handoff'
  call MPI_Finalize(ierr)

contains

  subroutine check_threshold(tolerance)
    real(8), intent(in) :: tolerance
    call initialize_dg_dc_handoff(handoff,.true.,3,tolerance,1,1.0d-8,ok,message)
    call require_collective(ok,'threshold controls')
    call evaluate_dg_dc_handoff(handoff,3,0.5d0*tolerance,4.0d0,4.0d0,.true.,.false.,.false.,.false., &
                                MPI_COMM_WORLD,accept,message)
    call require_collective(accept,'configured threshold acceptance')
  end subroutine check_threshold

  subroutine require_collective(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: all_condition
    integer :: mpi_error
    call MPI_Allreduce(condition,all_condition,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if (.not.all_condition) then
      if (rank==0) write(*,'(a,1x,a)') 'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    end if
  end subroutine require_collective
end program test_dg_dc_handoff_mpi
