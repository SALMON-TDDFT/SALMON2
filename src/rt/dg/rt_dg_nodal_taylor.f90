!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Reference Taylor propagator for the core-owned nodal DG route.
!=======================================================================
module rt_dg_nodal_taylor
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use rt_dg_nodal_halo, only: refresh_nodal_single_fragment_halos
  use rt_dg_nodal_mpi, only: exchange_nodal_face_halos
  use rt_dg_nodal_hamiltonian, only: apply_nodal_local_hamiltonian
  implicit none
  private
  public :: propagate_nodal_taylor, propagate_nodal_taylor_mpi, propagate_nodal_taylor_action
  public :: nodal_taylor_hamiltonian_action

  abstract interface
    subroutine nodal_taylor_hamiltonian_action(state,hpsi)
      import :: s_dg_nodal_state
      type(s_dg_nodal_state), intent(inout) :: state
      complex(8), intent(out) :: hpsi(:,:,:,:,:)
    end subroutine nodal_taylor_hamiltonian_action
  end interface

contains

  subroutine propagate_nodal_taylor(state, vlocal, kinetic_center, kinetic_offset, &
                                    gradient_offset, avec, dt, expansion_order)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: vlocal(:,:,:,:)
    real(8), intent(in) :: kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:)
    real(8), intent(in) :: avec(3), dt
    integer, intent(in) :: expansion_order
    complex(8), allocatable :: initial(:,:,:,:,:), term(:,:,:,:,:)
    complex(8), allocatable :: accumulated(:,:,:,:,:), hterm(:,:,:,:,:)
    integer :: iorder

    if (.not. state%ground_state_ready) stop 'nodal DG Taylor requires a verified DG ground state'
    if (state%dg_ground_state_residual /= state%dg_ground_state_residual) &
      stop 'nodal DG Taylor received an invalid DG ground-state residual'
    if (dt <= 0.0d0 .or. dt > 0.02d0) stop 'nodal DG Taylor requires 0 < dt <= 0.02 au'
    if (expansion_order < 1) stop 'nodal DG Taylor expansion order must be positive'
    allocate(initial, mold=state%psi_core)
    allocate(term, mold=state%psi_core)
    allocate(accumulated, mold=state%psi_core)
    allocate(hterm, mold=state%psi_core)
    initial = state%psi_core
    term = initial
    accumulated = initial

    do iorder = 1, expansion_order
      state%psi_core = term
      call refresh_nodal_single_fragment_halos(state)
      call apply_nodal_local_hamiltonian(state, vlocal, kinetic_center, kinetic_offset, &
                                         gradient_offset, avec, hterm)
      term = -(0.0d0,1.0d0) * dt * hterm / dble(iorder)
      accumulated = accumulated + term
    end do
    state%psi_core = accumulated
    deallocate(initial, term, accumulated, hterm)
  end subroutine propagate_nodal_taylor

  subroutine propagate_nodal_taylor_mpi(state, vlocal, kinetic_center, kinetic_offset, &
                                        gradient_offset, avec, dt, expansion_order, communicator)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: vlocal(:,:,:,:)
    real(8), intent(in) :: kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:)
    real(8), intent(in) :: avec(3), dt
    integer, intent(in) :: expansion_order, communicator
    complex(8), allocatable :: term(:,:,:,:,:), accumulated(:,:,:,:,:), hterm(:,:,:,:,:)
    integer :: iorder

    if (.not. state%ground_state_ready) stop 'nodal DG Taylor requires a verified DG ground state'
    if (state%dg_ground_state_residual /= state%dg_ground_state_residual) &
      stop 'nodal DG Taylor received an invalid DG ground-state residual'
    if (dt <= 0.0d0 .or. dt > 0.02d0) stop 'nodal DG Taylor requires 0 < dt <= 0.02 au'
    if (expansion_order < 1) stop 'nodal DG Taylor expansion order must be positive'
    allocate(term, mold=state%psi_core)
    allocate(accumulated, mold=state%psi_core)
    allocate(hterm, mold=state%psi_core)
    term = state%psi_core
    accumulated = state%psi_core

    do iorder = 1, expansion_order
      state%psi_core = term
      call exchange_nodal_face_halos(state, communicator)
      call apply_nodal_local_hamiltonian(state, vlocal, kinetic_center, kinetic_offset, &
                                         gradient_offset, avec, hterm)
      term = -(0.0d0,1.0d0) * dt * hterm / dble(iorder)
      accumulated = accumulated + term
    end do
    state%psi_core = accumulated
    deallocate(term, accumulated, hterm)
  end subroutine propagate_nodal_taylor_mpi

  subroutine propagate_nodal_taylor_action(state,apply_hamiltonian,dt,expansion_order)
    type(s_dg_nodal_state), intent(inout) :: state
    procedure(nodal_taylor_hamiltonian_action) :: apply_hamiltonian
    real(8), intent(in) :: dt
    integer, intent(in) :: expansion_order
    complex(8), allocatable :: term(:,:,:,:,:), accumulated(:,:,:,:,:), hterm(:,:,:,:,:)
    integer :: iorder

    if (.not. state%ground_state_ready) stop 'nodal DG Taylor requires a verified DG ground state'
    if (state%dg_ground_state_residual /= state%dg_ground_state_residual) &
      stop 'nodal DG Taylor received an invalid DG ground-state residual'
    if (dt <= 0.0d0 .or. dt > 0.02d0) stop 'nodal DG Taylor requires 0 < dt <= 0.02 au'
    if (expansion_order < 1) stop 'nodal DG Taylor expansion order must be positive'
    allocate(term,mold=state%psi_core)
    allocate(accumulated,mold=state%psi_core)
    allocate(hterm,mold=state%psi_core)
    term=state%psi_core
    accumulated=state%psi_core
    do iorder=1,expansion_order
      state%psi_core=term
      call apply_hamiltonian(state,hterm)
      term=-(0.0d0,1.0d0)*dt*hterm/dble(iorder)
      accumulated=accumulated+term
    end do
    state%psi_core=accumulated
    deallocate(term,accumulated,hterm)
  end subroutine propagate_nodal_taylor_action

end module rt_dg_nodal_taylor
