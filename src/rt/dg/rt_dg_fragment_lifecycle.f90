! DG fragment lifecycle helpers shared by non-SOI and SOI RT modules.
#include "config.h"
module rt_dg_fragment_lifecycle
  use rt_dg_fragment_types, only: s_dg_fragment_rt, invalidate_coef_exchange_cache
  implicit none
  private
  public :: init_rk_coefficients, reset_basis_dependent_operator_storage

contains

  subroutine init_rk_coefficients(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    select case(dg_frag%time_integrator)
    case(1)
      dg_frag%rk_stages = 3
      allocate(dg_frag%rk_alpha(3), dg_frag%rk_beta(3), dg_frag%rk_gamma(3))
      dg_frag%rk_alpha = [1.0d0, 0.75d0, 1.0d0/3.0d0]
      dg_frag%rk_beta  = [0.0d0, 0.25d0, 2.0d0/3.0d0]
      dg_frag%rk_gamma = [1.0d0, 0.25d0, 2.0d0/3.0d0]

    case(3)
      dg_frag%rk_stages = 4
      allocate(dg_frag%rk_alpha(4), dg_frag%rk_beta(4), dg_frag%rk_gamma(4))
      dg_frag%rk_alpha = [0.5d0, 0.5d0, 1.0d0, 0.0d0]
      dg_frag%rk_beta  = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
      dg_frag%rk_gamma = [1.0d0/6.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/6.0d0]

    case default
      dg_frag%rk_stages = 0
    end select
  end subroutine init_rk_coefficients

  subroutine reset_basis_dependent_operator_storage(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: i

    call invalidate_coef_exchange_cache(dg_frag)

    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_c)) deallocate(dg_frag%H_mat_c)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    if (allocated(dg_frag%H_ref_Vh_buffer)) deallocate(dg_frag%H_ref_Vh_buffer)
    if (allocated(dg_frag%H_ref_Vxc_buffer)) deallocate(dg_frag%H_ref_Vxc_buffer)
    dg_frag%H_delta_reference_valid = .false.
    if (allocated(dg_frag%H_mat_blocks)) then
      do i = 1, size(dg_frag%H_mat_blocks)
        if (allocated(dg_frag%H_mat_blocks(i)%val)) deallocate(dg_frag%H_mat_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_blocks)
    end if
    if (allocated(dg_frag%H_mat_core_blocks)) then
      do i = 1, size(dg_frag%H_mat_core_blocks)
        if (allocated(dg_frag%H_mat_core_blocks(i)%val)) deallocate(dg_frag%H_mat_core_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_core_blocks)
    end if
    if (allocated(dg_frag%H_mat_kinetic_blocks)) then
      do i = 1, size(dg_frag%H_mat_kinetic_blocks)
        if (allocated(dg_frag%H_mat_kinetic_blocks(i)%val)) deallocate(dg_frag%H_mat_kinetic_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_kinetic_blocks)
    end if
    if (allocated(dg_frag%H_nl_blocks)) then
      do i = 1, size(dg_frag%H_nl_blocks)
        if (allocated(dg_frag%H_nl_blocks(i)%val)) deallocate(dg_frag%H_nl_blocks(i)%val)
      end do
      deallocate(dg_frag%H_nl_blocks)
    end if
    if (allocated(dg_frag%H_block_map)) deallocate(dg_frag%H_block_map)
    if (allocated(dg_frag%H_nl_block_map)) deallocate(dg_frag%H_nl_block_map)
    if (allocated(dg_frag%H_local_block_ids)) deallocate(dg_frag%H_local_block_ids)
    if (allocated(dg_frag%H_nl_local_block_ids)) deallocate(dg_frag%H_nl_local_block_ids)
    if (allocated(dg_frag%H_local_rows)) deallocate(dg_frag%H_local_rows)
    if (allocated(dg_frag%nl_projector_overlap)) deallocate(dg_frag%nl_projector_overlap)
    if (allocated(dg_frag%nl_projector_r_overlap)) deallocate(dg_frag%nl_projector_r_overlap)
    if (allocated(dg_frag%nl_projector_overlap_halo)) deallocate(dg_frag%nl_projector_overlap_halo)
    if (allocated(dg_frag%nl_projector_r_overlap_halo)) deallocate(dg_frag%nl_projector_r_overlap_halo)
    dg_frag%n_H_blocks = 0
    dg_frag%n_H_nl_blocks = 0
    dg_frag%H_local_block_ids_valid = .false.
    dg_frag%nl_projector_cache_nlma = 0
    dg_frag%has_nl_projector_cache = .false.

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    if (allocated(dg_frag%S_mat_blocks)) then
      do i = 1, size(dg_frag%S_mat_blocks)
        if (allocated(dg_frag%S_mat_blocks(i)%val)) deallocate(dg_frag%S_mat_blocks(i)%val)
      end do
      deallocate(dg_frag%S_mat_blocks)
    end if
    if (allocated(dg_frag%S_mat_prop_blocks)) then
      do i = 1, size(dg_frag%S_mat_prop_blocks)
        if (allocated(dg_frag%S_mat_prop_blocks(i)%val)) deallocate(dg_frag%S_mat_prop_blocks(i)%val)
      end do
      deallocate(dg_frag%S_mat_prop_blocks)
    end if
    if (allocated(dg_frag%S_block_map)) deallocate(dg_frag%S_block_map)
    dg_frag%n_S_blocks = 0
    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .false.

    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    if (allocated(dg_frag%momentum_blocks)) then
      do i = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(i)%val)) deallocate(dg_frag%momentum_blocks(i)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%runtime_neighbor_frag_count)) deallocate(dg_frag%runtime_neighbor_frag_count)
    if (allocated(dg_frag%runtime_neighbor_frag_ids)) deallocate(dg_frag%runtime_neighbor_frag_ids)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_frag_count)) deallocate(dg_frag%momentum_neighbor_frag_count)
    if (allocated(dg_frag%momentum_neighbor_frag_ids)) deallocate(dg_frag%momentum_neighbor_frag_ids)
    dg_frag%n_momentum_blocks = 0
    dg_frag%momentum_blocks_include_dg_flux = .false.
    if (allocated(dg_frag%gradient_basis_cache)) deallocate(dg_frag%gradient_basis_cache)
    dg_frag%gradient_basis_cache_valid = .false.
    if (allocated(dg_frag%nl_projector_overlap)) deallocate(dg_frag%nl_projector_overlap)
    if (allocated(dg_frag%nl_projector_r_overlap)) deallocate(dg_frag%nl_projector_r_overlap)
    if (allocated(dg_frag%nl_projector_overlap_halo)) deallocate(dg_frag%nl_projector_overlap_halo)
    if (allocated(dg_frag%nl_projector_r_overlap_halo)) deallocate(dg_frag%nl_projector_r_overlap_halo)
    dg_frag%has_nl_projector_cache = .false.
    dg_frag%nl_projector_cache_nlma = 0
    dg_frag%Ac_nl_projector_cache = 0.0d0

    if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
    if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
    dg_frag%density_phi_block_size = 0
    dg_frag%density_phi_block_cache_valid = .false.
    if (allocated(dg_frag%coef_ref_all)) deallocate(dg_frag%coef_ref_all)
    dg_frag%coef_ref_ready = .false.
    dg_frag%has_nl_cache = .false.
  end subroutine reset_basis_dependent_operator_storage

end module rt_dg_fragment_lifecycle
