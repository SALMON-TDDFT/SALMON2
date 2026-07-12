!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Taylor propagation using the complete SALMON nodal Hamiltonian.
!=======================================================================
module rt_dg_nodal_salmon_taylor
  use structures, only: s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use rt_dg_nodal_taylor, only: propagate_nodal_taylor_action
  use rt_dg_nodal_salmon_hamiltonian, only: apply_nodal_salmon_hamiltonian_mpi
  implicit none
  private
  public :: propagate_nodal_salmon_taylor

contains

  subroutine propagate_nodal_salmon_taylor(state,dg_frag,ppg,ik,vlocal,kinetic_center,kinetic_offset, &
                                            gradient_offset,avec,dt,expansion_order,communicator)
    type(s_dg_nodal_state), intent(inout) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik, expansion_order, communicator
    real(8), intent(in) :: vlocal(:,:,:,:), kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:), avec(3), dt

    call propagate_nodal_taylor_action(state,apply_complete_h,dt,expansion_order)

  contains

    subroutine apply_complete_h(state_work,hpsi)
      type(s_dg_nodal_state), intent(inout) :: state_work
      complex(8), intent(out) :: hpsi(:,:,:,:,:)
      call apply_nodal_salmon_hamiltonian_mpi(state_work,dg_frag,ppg,ik,vlocal,kinetic_center, &
                                               kinetic_offset,gradient_offset,avec,communicator,hpsi)
    end subroutine apply_complete_h

  end subroutine propagate_nodal_salmon_taylor

end module rt_dg_nodal_salmon_taylor
