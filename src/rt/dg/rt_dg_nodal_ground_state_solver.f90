!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Matrix-free occupied-subspace relaxation for nodal DG initial states.
!=======================================================================
module rt_dg_nodal_ground_state_solver
#ifdef USE_MPI
  use mpi
#endif
  use rt_dg_nodal_types, only: s_dg_nodal_state, accept_nodal_dg_ground_state
  use dg_nodal_interfaces, only: nodal_hamiltonian_action => nodal_complete_h_action
  use dg_nodal_interfaces, only: nodal_subspace_rotation
  use rt_dg_nodal_mpi, only: exchange_nodal_face_halos
  use rt_dg_nodal_hamiltonian, only: apply_nodal_local_hamiltonian
  use rt_dg_nodal_ground_state, only: verify_nodal_dg_eigenstate_mpi
  implicit none
  private
  public :: relax_nodal_dg_ground_state_mpi, relax_nodal_ground_state_action_mpi
  public :: orthonormalize_nodal_states_mpi, nodal_hamiltonian_action

contains

  subroutine relax_nodal_dg_ground_state_mpi(state, vlocal, kinetic_center, kinetic_offset, &
                                              gradient_offset, relaxation_step, max_iteration, &
                                              tolerance, communicator, eigenvalues, max_residual, niteration)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: vlocal(:,:,:,:), kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:)
    real(8), intent(in) :: relaxation_step, tolerance
    integer, intent(in) :: max_iteration, communicator
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin), max_residual
    integer, intent(out) :: niteration
    complex(8), allocatable :: hpsi(:,:,:,:,:), residual(:,:,:)
    complex(8) :: expectation_local, expectation_global
    real(8) :: norm_local, norm_global, residual_local, residual_global
    integer :: iteration, istate, ispin

    if (relaxation_step <= 0.0d0) stop 'nodal DG: ground-state relaxation step must be positive'
    if (max_iteration < 1 .or. tolerance <= 0.0d0) stop 'nodal DG: invalid ground-state solver control'
    state%ground_state_ready = .false.
    state%dg_ground_state_residual = huge(1.0d0)
    allocate(hpsi, mold=state%psi_core)
    call orthonormalize_nodal_states_mpi(state, communicator)

    do iteration = 1, max_iteration
      call exchange_nodal_face_halos(state, communicator)
      call apply_nodal_local_hamiltonian(state, vlocal, kinetic_center, kinetic_offset, &
                                         gradient_offset, [0.0d0,0.0d0,0.0d0], hpsi)
      max_residual = 0.0d0
      do ispin = 1, state%nspin
        do istate = 1, state%nstate
          expectation_local = sum(conjg(state%psi_core(:,:,:,istate,ispin)) * hpsi(:,:,:,istate,ispin))
          norm_local = sum(abs(state%psi_core(:,:,:,istate,ispin))**2)
          call global_complex_sum(expectation_local, expectation_global, communicator)
          call global_real_sum(norm_local, norm_global, communicator)
          if (norm_global <= 0.0d0) stop 'nodal DG: zero norm during ground-state relaxation'
          eigenvalues(istate,ispin) = real(expectation_global,8) / norm_global
          allocate(residual(state%core_size(1),state%core_size(2),state%core_size(3)))
          residual = hpsi(:,:,:,istate,ispin) - eigenvalues(istate,ispin) * &
                     state%psi_core(:,:,:,istate,ispin)
          residual_local = sum(abs(residual)**2)
          call global_real_sum(residual_local, residual_global, communicator)
          max_residual = max(max_residual, sqrt(residual_global / norm_global))
          if (max_residual > tolerance) then
            state%psi_core(:,:,:,istate,ispin) = state%psi_core(:,:,:,istate,ispin) - &
                                                  relaxation_step * residual
          end if
          deallocate(residual)
        end do
      end do
      if (max_residual <= tolerance) exit
      call orthonormalize_nodal_states_mpi(state, communicator)
    end do
    niteration = iteration
    if (max_residual > tolerance) stop 'nodal DG: ground-state relaxation did not converge'
    call verify_nodal_dg_eigenstate_mpi(state, eigenvalues, vlocal, kinetic_center, kinetic_offset, &
                                        gradient_offset, [0.0d0,0.0d0,0.0d0], communicator, &
                                        tolerance, max_residual)
    deallocate(hpsi)
  end subroutine relax_nodal_dg_ground_state_mpi

  subroutine relax_nodal_ground_state_action_mpi(state,apply_hamiltonian,relaxation_step,max_iteration, &
                                                  tolerance,communicator,eigenvalues,max_residual,niteration, &
                                                  rotate_subspace)
    type(s_dg_nodal_state), intent(inout) :: state
    procedure(nodal_hamiltonian_action) :: apply_hamiltonian
    real(8), intent(in) :: relaxation_step, tolerance
    integer, intent(in) :: max_iteration, communicator
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin), max_residual
    integer, intent(out) :: niteration
    procedure(nodal_subspace_rotation), optional :: rotate_subspace
    complex(8), allocatable :: hpsi(:,:,:,:,:), residual(:,:,:)
    complex(8) :: expectation_local, expectation_global
    real(8) :: norm_local, norm_global, residual_local, residual_global
    integer :: iteration, istate, ispin, myrank

    if (relaxation_step <= 0.0d0) stop 'nodal DG: callback GS relaxation step must be positive'
    if (max_iteration < 1 .or. tolerance <= 0.0d0) stop 'nodal DG: invalid callback GS control'
    state%ground_state_ready = .false.
    state%dg_ground_state_residual = huge(1.0d0)
    allocate(hpsi,mold=state%psi_core)
    myrank=0
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,myrank,istate)
#endif
    call orthonormalize_nodal_states_mpi(state,communicator)
    if (myrank == 0) then
      write(*,'(1x,a)') '[DG-NODAL-GS-PHASE] initial orthonormalization complete'
      flush(6)
    end if

    do iteration=1,max_iteration
      call apply_hamiltonian(state,hpsi)
      if (present(rotate_subspace)) call rotate_subspace(state,hpsi,eigenvalues,communicator)
      if (myrank == 0 .and. iteration == 1) then
        write(*,'(1x,a)') '[DG-NODAL-GS-PHASE] first complete-H action complete'
        flush(6)
      end if
      max_residual=0.0d0
      do ispin=1,state%nspin
        do istate=1,state%nstate
          expectation_local=sum(conjg(state%psi_core(:,:,:,istate,ispin))*hpsi(:,:,:,istate,ispin))
          norm_local=sum(abs(state%psi_core(:,:,:,istate,ispin))**2)
          call global_complex_sum(expectation_local,expectation_global,communicator)
          call global_real_sum(norm_local,norm_global,communicator)
          if (norm_global <= 0.0d0) stop 'nodal DG: zero norm in callback GS solver'
          eigenvalues(istate,ispin)=real(expectation_global,8)/norm_global
          allocate(residual(state%core_size(1),state%core_size(2),state%core_size(3)))
          residual=hpsi(:,:,:,istate,ispin)-eigenvalues(istate,ispin)*state%psi_core(:,:,:,istate,ispin)
          residual_local=sum(abs(residual)**2)
          call global_real_sum(residual_local,residual_global,communicator)
          max_residual=max(max_residual,sqrt(residual_global/norm_global))
          deallocate(residual)
        end do
      end do
      if (myrank == 0 .and. (iteration == 1 .or. mod(iteration,10) == 0 .or. max_residual <= tolerance)) then
        write(*,'(1x,a,i0,3(a,1pe13.5))') '[DG-NODAL-GS-ITER] iter=',iteration, &
          ' residual=',max_residual,' eigen_min=',minval(eigenvalues),' eigen_max=',maxval(eigenvalues)
        flush(6)
      end if
      if (max_residual <= tolerance) exit
      do ispin=1,state%nspin
        do istate=1,state%nstate
          state%psi_core(:,:,:,istate,ispin)=state%psi_core(:,:,:,istate,ispin)-relaxation_step* &
            (hpsi(:,:,:,istate,ispin)-eigenvalues(istate,ispin)*state%psi_core(:,:,:,istate,ispin))
        end do
      end do
      call orthonormalize_nodal_states_mpi(state,communicator)
    end do
    niteration=iteration
    if (max_residual > tolerance) stop 'nodal DG: callback ground-state relaxation did not converge'

    ! Re-evaluate the accepted state with the exact callback used above.
    call orthonormalize_nodal_states_mpi(state,communicator)
    call apply_hamiltonian(state,hpsi)
    if (present(rotate_subspace)) call rotate_subspace(state,hpsi,eigenvalues,communicator)
    max_residual=0.0d0
    do ispin=1,state%nspin
      do istate=1,state%nstate
        expectation_local=sum(conjg(state%psi_core(:,:,:,istate,ispin))*hpsi(:,:,:,istate,ispin))
        norm_local=sum(abs(state%psi_core(:,:,:,istate,ispin))**2)
        call global_complex_sum(expectation_local,expectation_global,communicator)
        call global_real_sum(norm_local,norm_global,communicator)
        eigenvalues(istate,ispin)=real(expectation_global,8)/norm_global
        allocate(residual(state%core_size(1),state%core_size(2),state%core_size(3)))
        residual=hpsi(:,:,:,istate,ispin)-eigenvalues(istate,ispin)*state%psi_core(:,:,:,istate,ispin)
        residual_local=sum(abs(residual)**2)
        call global_real_sum(residual_local,residual_global,communicator)
        max_residual=max(max_residual,sqrt(residual_global/norm_global))
        deallocate(residual)
      end do
    end do
    call accept_nodal_dg_ground_state(state,max_residual,tolerance)
    deallocate(hpsi)
  end subroutine relax_nodal_ground_state_action_mpi

  subroutine orthonormalize_nodal_states_mpi(state, communicator)
    type(s_dg_nodal_state), intent(inout) :: state
    integer, intent(in) :: communicator
    complex(8), allocatable :: overlap_local(:,:), overlap_global(:,:)
    integer :: ispin,npoint,info
#ifdef USE_MPI
    integer :: ierr
#endif

    allocate(overlap_local(state%nstate,state%nstate),overlap_global(state%nstate,state%nstate))
    npoint=size(state%psi_core,1)*size(state%psi_core,2)*size(state%psi_core,3)
    do ispin = 1, state%nspin
      call zgemm('C','N',state%nstate,state%nstate,npoint,(1.0d0,0.0d0), &
        state%psi_core(1,1,1,1,ispin),npoint,state%psi_core(1,1,1,1,ispin),npoint, &
        (0.0d0,0.0d0),overlap_local(1,1),state%nstate)
#ifdef USE_MPI
      call MPI_Allreduce(overlap_local,overlap_global,state%nstate*state%nstate,MPI_DOUBLE_COMPLEX, &
        MPI_SUM,communicator,ierr)
      if(ierr/=MPI_SUCCESS) stop 'nodal DG: overlap-matrix reduction failed'
#else
      overlap_global=overlap_local
      if(communicator<0) stop 'nodal DG: invalid serial communicator'
#endif
      call zpotrf('U',state%nstate,overlap_global,state%nstate,info)
      if(info/=0) stop 'nodal DG: linearly dependent ground-state trial orbitals'
      call ztrsm('R','U','N','N',npoint,state%nstate,(1.0d0,0.0d0),overlap_global,state%nstate, &
        state%psi_core(1,1,1,1,ispin),npoint)
    end do
    deallocate(overlap_local,overlap_global)
  end subroutine orthonormalize_nodal_states_mpi

  subroutine global_complex_array_sum(local_value,global_value,communicator)
    complex(8), intent(in) :: local_value(:)
    complex(8), intent(out) :: global_value(:)
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,size(local_value),MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: complex array reduction failed'
#else
    global_value=local_value
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
  end subroutine global_complex_array_sum

  subroutine global_real_sum(local_value, global_value, communicator)
    real(8), intent(in) :: local_value
    real(8), intent(out) :: global_value
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: real reduction failed'
#else
    global_value = local_value
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
  end subroutine global_real_sum

  subroutine global_complex_sum(local_value, global_value, communicator)
    complex(8), intent(in) :: local_value
    complex(8), intent(out) :: global_value
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, communicator, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: complex reduction failed'
#else
    global_value = local_value
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
  end subroutine global_complex_sum

end module rt_dg_nodal_ground_state_solver
