!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Stationarity verification for nodal DG initial orbitals.
!=======================================================================
module rt_dg_nodal_ground_state
#ifdef USE_MPI
  use mpi
#endif
  use rt_dg_nodal_types, only: s_dg_nodal_state, accept_nodal_dg_ground_state
  use rt_dg_nodal_mpi, only: exchange_nodal_face_halos
  use rt_dg_nodal_hamiltonian, only: apply_nodal_local_hamiltonian
  implicit none
  private
  public :: verify_nodal_dg_eigenstate_mpi

contains

  subroutine verify_nodal_dg_eigenstate_mpi(state, eigenvalues, vlocal, kinetic_center, &
                                             kinetic_offset, gradient_offset, avec, communicator, &
                                             tolerance, max_residual)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: eigenvalues(:,:)
    real(8), intent(in) :: vlocal(:,:,:,:), kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:), avec(3)
    integer, intent(in) :: communicator
    real(8), intent(in) :: tolerance
    real(8), intent(out) :: max_residual
    complex(8), allocatable :: hpsi(:,:,:,:,:), residual_field(:,:,:)
    real(8), allocatable :: residual_local(:,:), residual_global(:,:)
    real(8), allocatable :: norm_local(:,:), norm_global(:,:)
    integer :: istate, ispin
#ifdef USE_MPI
    integer :: ierr
#endif

    if (any(shape(eigenvalues) /= [state%nstate,state%nspin])) &
      stop 'nodal DG: eigenvalue shape mismatch in stationarity check'
    allocate(hpsi, mold=state%psi_core)
    allocate(residual_local(state%nstate,state%nspin), residual_global(state%nstate,state%nspin))
    allocate(norm_local(state%nstate,state%nspin), norm_global(state%nstate,state%nspin))
    residual_local = 0.0d0
    norm_local = 0.0d0

    call exchange_nodal_face_halos(state, communicator)
    call apply_nodal_local_hamiltonian(state, vlocal, kinetic_center, kinetic_offset, &
                                       gradient_offset, avec, hpsi)
    do ispin = 1, state%nspin
      do istate = 1, state%nstate
        allocate(residual_field(state%core_size(1),state%core_size(2),state%core_size(3)))
        ! The residual is hpsi - eigenvalues * psi for exactly the same
        ! Hamiltonian action that will be used by real-time propagation.
        residual_field = hpsi(:,:,:,istate,ispin) - eigenvalues(istate,ispin) * &
                         state%psi_core(:,:,:,istate,ispin)
        residual_local(istate,ispin) = sum(abs(residual_field)**2)
        norm_local(istate,ispin) = sum(abs(state%psi_core(:,:,:,istate,ispin))**2)
        deallocate(residual_field)
      end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(residual_local, residual_global, size(residual_local), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, communicator, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: residual MPI reduction failed'
    call MPI_Allreduce(norm_local, norm_global, size(norm_local), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, communicator, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: norm MPI reduction failed'
#else
    residual_global = residual_local
    norm_global = norm_local
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
    if (any(norm_global <= 0.0d0)) stop 'nodal DG: zero-norm orbital in stationarity check'
    max_residual = maxval(sqrt(residual_global / norm_global))
    call accept_nodal_dg_ground_state(state, max_residual, tolerance)
    deallocate(hpsi, residual_local, residual_global, norm_local, norm_global)
  end subroutine verify_nodal_dg_eigenstate_mpi

end module rt_dg_nodal_ground_state
