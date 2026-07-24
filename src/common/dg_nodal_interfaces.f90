!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
module dg_nodal_interfaces
  use dg_nodal_state, only: s_dg_nodal_common_state
  implicit none
  private
  public :: nodal_complete_h_action, nodal_density_builder
  public :: nodal_subspace_rotation, nodal_collective_validator
  public :: nodal_scaled_complete_h_action,nodal_density_potential_update

  abstract interface
    subroutine nodal_complete_h_action(state,hpsi)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(inout) :: state
      complex(8), intent(out) :: hpsi(:,:,:,:,:)
    end subroutine nodal_complete_h_action

    subroutine nodal_density_builder(state,density,electron_number,communicator)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(in) :: state
      real(8), intent(out) :: density(:,:,:,:)
      real(8), intent(out) :: electron_number
      integer, intent(in) :: communicator
    end subroutine nodal_density_builder

    subroutine nodal_subspace_rotation(state,hpsi,eigenvalues,communicator)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(inout) :: state
      complex(8), intent(inout) :: hpsi(:,:,:,:,:)
      real(8), intent(out) :: eigenvalues(state%nstate,state%nspin)
      integer, intent(in) :: communicator
    end subroutine nodal_subspace_rotation

    subroutine nodal_collective_validator(state,communicator,ok,message)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(in) :: state
      integer, intent(in) :: communicator
      logical, intent(out) :: ok
      character(*), intent(out) :: message
    end subroutine nodal_collective_validator

    subroutine nodal_scaled_complete_h_action(state,lambda,hpsi)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(in) :: lambda
      complex(8), intent(out) :: hpsi(:,:,:,:,:)
    end subroutine nodal_scaled_complete_h_action

    subroutine nodal_density_potential_update(state,density,mix_rate,communicator,density_residual,electron_number,ok)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(inout) :: density(:,:,:,:)
      real(8), intent(in) :: mix_rate
      integer, intent(in) :: communicator
      real(8), intent(out) :: density_residual,electron_number
      logical, intent(out) :: ok
    end subroutine nodal_density_potential_update
  end interface
end module dg_nodal_interfaces
