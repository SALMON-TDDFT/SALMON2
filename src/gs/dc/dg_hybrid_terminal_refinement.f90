#include "config.h"
module dg_hybrid_terminal_refinement
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public :: s_dg_hybrid_terminal_refinement_controls
    integer :: maximum_additional_solves=3
    real(real64) :: density_tolerance=0.0_real64
    real(real64) :: energy_tolerance=0.0_real64
  end type s_dg_hybrid_terminal_refinement_controls

  type,public :: s_dg_hybrid_terminal_refinement_state
    logical :: valid=.false.,terminal=.false.
    integer :: total_solve_count=0
    type(s_dg_hybrid_terminal_refinement_controls) :: controls
  end type s_dg_hybrid_terminal_refinement_state

  type,public :: s_dg_hybrid_terminal_refinement_receipt
    logical :: valid=.false.,converged=.false.,exhausted=.false.
    logical :: publish_last_valid=.false.
    integer :: total_solve_count=0,additional_refinement_count=0
    real(real64) :: density_change=huge(0.0_real64)
    real(real64) :: energy_change=huge(0.0_real64)
    integer(int64) :: fingerprint=0_int64
  end type s_dg_hybrid_terminal_refinement_receipt

  public :: initialize_dg_hybrid_terminal_refinement
  public :: observe_dg_hybrid_terminal_refinement

contains

  subroutine initialize_dg_hybrid_terminal_refinement(comm,controls,state,ok,message)
    integer,intent(in) :: comm
    type(s_dg_hybrid_terminal_refinement_controls),intent(in) :: controls
    type(s_dg_hybrid_terminal_refinement_state),intent(out) :: state
    logical,intent(out) :: ok
    character(*),intent(out) :: message
#ifdef USE_MPI
    integer :: ierr,local_bad,global_bad,minimum_count,maximum_count
    integer(int64) :: local_bits(2),minimum_bits(2),maximum_bits(2)

    state=s_dg_hybrid_terminal_refinement_state();ok=.false.;message=''
    local_bad=merge(0,1,controls%maximum_additional_solves>=0.and.&
      controls%maximum_additional_solves<=3.and.&
      valid_positive(controls%density_tolerance).and.valid_positive(controls%energy_tolerance))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid terminal LCFO refinement controls';return
    endif
    call MPI_Allreduce(controls%maximum_additional_solves,minimum_count,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(controls%maximum_additional_solves,maximum_count,1,&
      MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bits=[transfer(controls%density_tolerance,0_int64),transfer(controls%energy_tolerance,0_int64)]
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_bits,minimum_bits,2,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_bits,maximum_bits,2,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='terminal LCFO refinement control agreement failed';return
    endif
    if(minimum_count/=maximum_count.or.any(minimum_bits/=maximum_bits))then
      message='rank-disagreeing terminal LCFO refinement controls';return
    endif
    state%valid=.true.;state%terminal=.false.;state%total_solve_count=0
    state%controls=controls;ok=.true.;message=''
#else
    state=s_dg_hybrid_terminal_refinement_state();ok=.false.
    message='terminal LCFO refinement requires MPI'
#endif
  end subroutine initialize_dg_hybrid_terminal_refinement

  subroutine observe_dg_hybrid_terminal_refinement(comm,state,density_change,energy_change,&
      finite_state,request_another,receipt,ok,message)
    integer,intent(in) :: comm
    type(s_dg_hybrid_terminal_refinement_state),intent(inout) :: state
    real(real64),intent(in) :: density_change,energy_change
    logical,intent(in) :: finite_state
    logical,intent(out) :: request_another
    type(s_dg_hybrid_terminal_refinement_receipt),intent(out) :: receipt
    logical,intent(out) :: ok
    character(*),intent(out) :: message
#ifdef USE_MPI
    integer :: ierr,local_bad,global_bad,local_finite,minimum_finite,maximum_finite
    integer :: minimum_count,maximum_count
    integer(int64) :: local_bits(2),minimum_bits(2),maximum_bits(2)
    logical :: converged

    request_another=.false.;receipt=s_dg_hybrid_terminal_refinement_receipt()
    ok=.false.;message=''
    local_bad=merge(0,1,state%valid.and..not.state%terminal.and.finite_state.and.&
      state%total_solve_count>=0.and.&
      state%total_solve_count<1+state%controls%maximum_additional_solves.and.&
      valid_nonnegative(density_change).and.valid_nonnegative(energy_change))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid terminal LCFO refinement sample';return
    endif
    call MPI_Allreduce(state%total_solve_count,minimum_count,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%total_solve_count,maximum_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_finite=merge(1,0,finite_state)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_finite,minimum_finite,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_finite,maximum_finite,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bits=[transfer(density_change,0_int64),transfer(energy_change,0_int64)]
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_bits,minimum_bits,2,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_bits,maximum_bits,2,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='terminal LCFO refinement sample agreement failed';return
    endif
    if(minimum_count/=maximum_count.or.minimum_finite/=maximum_finite.or.&
       any(minimum_bits/=maximum_bits))then
      message='rank-disagreeing terminal LCFO refinement sample';return
    endif

    state%total_solve_count=state%total_solve_count+1
    converged=density_change<=state%controls%density_tolerance.and.&
      energy_change<=state%controls%energy_tolerance
    receipt%valid=.true.;receipt%converged=converged
    receipt%total_solve_count=state%total_solve_count
    receipt%additional_refinement_count=max(0,state%total_solve_count-1)
    receipt%density_change=density_change;receipt%energy_change=energy_change
    if(converged)then
      state%terminal=.true.
    else if(state%total_solve_count<1+state%controls%maximum_additional_solves)then
      request_another=.true.
    else
      state%terminal=.true.;receipt%exhausted=.true.;receipt%publish_last_valid=.true.
    endif
    receipt%fingerprint=receipt_fingerprint(state%controls,receipt)
    ok=.true.;message=''
#else
    request_another=.false.;receipt=s_dg_hybrid_terminal_refinement_receipt()
    ok=.false.;message='terminal LCFO refinement requires MPI'
#endif
  end subroutine observe_dg_hybrid_terminal_refinement

  pure integer(int64) function receipt_fingerprint(controls,receipt) result(hash)
    type(s_dg_hybrid_terminal_refinement_controls),intent(in) :: controls
    type(s_dg_hybrid_terminal_refinement_receipt),intent(in) :: receipt
    integer(int64) :: density_bits,energy_bits,tolerance_bits

    density_bits=transfer(receipt%density_change,density_bits)
    energy_bits=transfer(receipt%energy_change,energy_bits)
    tolerance_bits=ieor(transfer(controls%density_tolerance,tolerance_bits),&
      ishftc(transfer(controls%energy_tolerance,tolerance_bits),11))
    hash=ieor(int(receipt%total_solve_count,int64),ishft(int(receipt%additional_refinement_count,int64),8))
    hash=ieor(hash,ishftc(density_bits,17));hash=ieor(hash,ishftc(energy_bits,31))
    hash=ieor(hash,ishftc(tolerance_bits,43))
    hash=ieor(hash,ishft(int(merge(1,0,receipt%converged),int64),56))
    hash=ieor(hash,ishft(int(merge(1,0,receipt%exhausted),int64),57))
    hash=ieor(hash,ishft(int(merge(1,0,receipt%publish_last_valid),int64),58))
    if(hash==0_int64)hash=int(z'524546494E454D54',int64)
  end function receipt_fingerprint

  pure logical function valid_positive(value)
    real(real64),intent(in) :: value
    if(.not.ieee_is_finite(value))then
      valid_positive=.false.
    else
      valid_positive=value>0.0_real64
    endif
  end function valid_positive

  pure logical function valid_nonnegative(value)
    real(real64),intent(in) :: value
    if(.not.ieee_is_finite(value))then
      valid_nonnegative=.false.
    else
      valid_nonnegative=value>=0.0_real64
    endif
  end function valid_nonnegative
end module dg_hybrid_terminal_refinement
