!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Ground-state relaxation using the complete SALMON nodal Hamiltonian.
!=======================================================================
module rt_dg_nodal_salmon_ground_state
  use structures, only: s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use rt_dg_nodal_cg, only: solve_nodal_ground_state_cg_mpi
  use rt_dg_nodal_salmon_hamiltonian, only: apply_nodal_salmon_hamiltonian_mpi
  use rt_dg_nodal_rayleigh_ritz, only: rayleigh_ritz_nodal_subspace_mpi
  implicit none
  private
  public :: relax_nodal_salmon_ground_state_mpi

contains

  subroutine relax_nodal_salmon_ground_state_mpi(state,dg_frag,ppg,ik,vlocal,kinetic_center, &
                                                  kinetic_offset,gradient_offset,relaxation_step, &
                                                  max_iteration,tolerance,communicator,eigenvalues, &
                                                  max_residual,niteration)
    type(s_dg_nodal_state), intent(inout) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik, max_iteration, communicator
    real(8), intent(in) :: vlocal(:,:,:,:), kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:)
    real(8), intent(in) :: relaxation_step, tolerance
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin), max_residual
    integer, intent(out) :: niteration

    call solve_nodal_ground_state_cg_mpi(state,apply_complete_h,max_iteration,tolerance,communicator, &
                                         eigenvalues,max_residual,niteration,rotate_complete_h_subspace)

  contains

    subroutine apply_complete_h(state_work,hpsi)
      type(s_dg_nodal_state), intent(inout) :: state_work
      complex(8), intent(out) :: hpsi(:,:,:,:,:)
      call apply_nodal_salmon_hamiltonian_mpi(state_work,dg_frag,ppg,ik,vlocal,kinetic_center, &
                                               kinetic_offset,gradient_offset,[0.0d0,0.0d0,0.0d0], &
                                               communicator,hpsi)
    end subroutine apply_complete_h

    subroutine rotate_complete_h_subspace(state_work,hpsi,eigenvalues_work,communicator_work)
      type(s_dg_nodal_state), intent(inout) :: state_work
      complex(8), intent(inout) :: hpsi(:,:,:,:,:)
      real(8), intent(out) :: eigenvalues_work(state_work%nstate,state_work%nspin)
      integer, intent(in) :: communicator_work
      call rayleigh_ritz_nodal_subspace_mpi(state_work,hpsi,eigenvalues_work,communicator_work)
    end subroutine rotate_complete_h_subspace

  end subroutine relax_nodal_salmon_ground_state_mpi

end module rt_dg_nodal_salmon_ground_state
