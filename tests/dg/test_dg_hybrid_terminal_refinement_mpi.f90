program test_dg_hybrid_terminal_refinement_mpi
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use mpi
  use dg_hybrid_terminal_refinement, only: &
    s_dg_hybrid_terminal_refinement_controls, &
    s_dg_hybrid_terminal_refinement_state, &
    s_dg_hybrid_terminal_refinement_receipt, &
    initialize_dg_hybrid_terminal_refinement, &
    observe_dg_hybrid_terminal_refinement
  implicit none

  integer :: ierr, rank, nproc
  type(s_dg_hybrid_terminal_refinement_controls) :: controls
  type(s_dg_hybrid_terminal_refinement_state) :: state
  type(s_dg_hybrid_terminal_refinement_receipt) :: receipt
  logical :: ok, request_another
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  controls%maximum_additional_solves=3
  controls%density_tolerance=1.0e-6_real64
  controls%energy_tolerance=1.0e-7_real64

  call run_converged_sequence(1)
  call run_converged_sequence(2)
  call run_converged_sequence(3)
  call run_converged_sequence(4)
  call run_exhausted_sequence
  call run_nonfinite_rejection
  call run_rank_disagreement_rejection

  if(rank==0) write(*,'(a,i0,a)') 'PASS terminal LCFO refinement policy on ',nproc,' ranks'
  call MPI_Finalize(ierr)

contains

  subroutine run_converged_sequence(target_count)
    integer,intent(in) :: target_count
    integer :: sample

    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,controls,state,ok,message)
    call require(ok,'valid policy initialization: '//trim(message))
    do sample=1,target_count
      if(sample<target_count)then
        call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,1.0e-3_real64,&
          1.0e-4_real64,.true.,request_another,receipt,ok,message)
        call require(ok.and.request_another,'unconverged sample must request another solve')
      else
        call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,1.0e-8_real64,&
          1.0e-9_real64,.true.,request_another,receipt,ok,message)
        call require(ok.and..not.request_another,'converged sample requested another solve')
      endif
    enddo
    call require(receipt%valid.and.receipt%converged.and..not.receipt%exhausted,&
      'converged receipt flags')
    call require(.not.receipt%publish_last_valid,'converged receipt requested warning publication')
    call require(receipt%total_solve_count==target_count,'wrong converged total solve count')
    call require(receipt%additional_refinement_count==target_count-1,&
      'wrong converged additional solve count')
    call require(receipt%density_change==1.0e-8_real64.and.&
      receipt%energy_change==1.0e-9_real64,'wrong converged terminal metrics')
    call require(receipt%fingerprint/=0_int64,'converged receipt fingerprint is zero')
  end subroutine run_converged_sequence

  subroutine run_exhausted_sequence
    integer :: sample

    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,controls,state,ok,message)
    call require(ok,'exhaustion policy initialization')
    do sample=1,4
      call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,1.0e-3_real64,&
        1.0e-4_real64,.true.,request_another,receipt,ok,message)
      call require(ok,'finite failed sample was rejected')
      if(sample<4) call require(request_another,'policy stopped before three additional solves')
    enddo
    call require(.not.request_another,'policy requested a fifth total solve')
    call require(receipt%valid.and..not.receipt%converged.and.receipt%exhausted,&
      'exhausted receipt flags')
    call require(receipt%publish_last_valid,'exhausted receipt did not retain last valid state')
    call require(receipt%total_solve_count==4.and.receipt%additional_refinement_count==3,&
      'exhausted solve counts')
    call require(receipt%fingerprint/=0_int64,'exhausted receipt fingerprint is zero')
    call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,1.0e-3_real64,&
      1.0e-4_real64,.true.,request_another,receipt,ok,message)
    call require(.not.ok,'sample after terminal exhaustion was accepted')
  end subroutine run_exhausted_sequence

  subroutine run_nonfinite_rejection
    real(real64) :: nan_value
    nan_value=ieee_value(0.0_real64,ieee_quiet_nan)
    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,controls,state,ok,message)
    call require(ok,'nonfinite test initialization')
    call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,nan_value,&
      0.0_real64,.true.,request_another,receipt,ok,message)
    call require(.not.ok,'nonfinite density metric was accepted')

    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,controls,state,ok,message)
    call require(ok,'invalid-state test initialization')
    call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,0.0_real64,&
      0.0_real64,.false.,request_another,receipt,ok,message)
    call require(.not.ok,'invalid LCFO state was accepted')
  end subroutine run_nonfinite_rejection

  subroutine run_rank_disagreement_rejection
    type(s_dg_hybrid_terminal_refinement_controls) :: disagreeing
    real(real64) :: density_sample

    disagreeing=controls
    if(nproc>1.and.rank==nproc-1) disagreeing%maximum_additional_solves=2
    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,disagreeing,state,ok,message)
    if(nproc>1)then
      call require(.not.ok,'rank-disagreeing controls were accepted')
    else
      call require(ok,'single-rank controls were rejected')
    endif

    call initialize_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,controls,state,ok,message)
    call require(ok,'sample disagreement test initialization')
    density_sample=1.0e-8_real64
    if(nproc>1.and.rank==nproc-1) density_sample=2.0e-8_real64
    call observe_dg_hybrid_terminal_refinement(MPI_COMM_WORLD,state,density_sample,&
      1.0e-9_real64,.true.,request_another,receipt,ok,message)
    if(nproc>1)then
      call require(.not.ok,'rank-disagreeing metrics were accepted')
    else
      call require(ok,'single-rank sample was rejected')
    endif
  end subroutine run_rank_disagreement_rejection

  subroutine require(condition,why)
    logical,intent(in) :: condition
    character(*),intent(in) :: why
    if(.not.condition)then
      write(0,'(a,i0,2a)') 'rank ',rank,': ',trim(why)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_terminal_refinement_mpi
