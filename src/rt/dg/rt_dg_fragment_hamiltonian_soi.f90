!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

#include "config.h"
module rt_dg_fragment_hamiltonian_soi
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info, halo_info
  implicit none

  private
  public :: calculate_hamiltonian_matrix, calculate_momentum_matrix, calculate_overlap_matrix
  public :: find_matrix_block, init_matrix_blocks, reduce_matrix_blocks
  public :: build_total_potential_grid, build_hpsi_for_basis, integrate_basis_with_field

contains

!=======================================================================
  ! Calculate Hamiltonian matrix in fragment basis
  !=======================================================================
  !=======================================================================
  ! Calculate initial Hamiltonian matrix from basis functions
  !
  ! Includes halo (ghost cell) exchange for accurate boundary treatment.
  ! System boundaries use PERIODIC boundary conditions (full system is periodic).
  ! Fragment boundaries are handled via MPI communication between neighboring fragments.
  !=======================================================================
  subroutine calculate_hamiltonian_matrix(dg_frag, system, lg, mg, stencil, &
                                         Vh, Vxc, Vpsl, pp, ppg)
    use structures
    use communication, only: comm_is_root, comm_summation
    use parallelization, only: nproc_size_global
    use rt_dg_fragment_ops, only: symmetrize_real_matrix_blocks, rebuild_local_h_block_ids
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    
    integer :: ifrag, ispin, io, jo, i_local, nbf, ig_i, ig_j, iblk
    real(8) :: hvol, integral_re, integral_im, gval, phi_re, phi_im
    complex(8) :: integral_t, integral_h
    real(8) :: max_p
    integer :: is(3), ie(3)
    complex(8), allocatable :: T_phi(:,:,:)   ! Kinetic energy operator applied to basis
    complex(8), allocatable :: H_phi(:,:,:)   ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j>
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    if (.not. dg_frag%has_real_space_basis) then
      return
    end if

    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Preparing Hamiltonian Matrix ==="
    end if
    
    ! Step 1: Calculate momentum matrix elements (transition moments)
    ! Required for velocity gauge A·p coupling
    if (.not. allocated(dg_frag%momentum_blocks)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Calculating momentum matrix elements (p_ij)..."
        write(*,*) "        Using 4th-order finite difference stencil"
      end if
      call calculate_momentum_matrix(dg_frag, system, mg, stencil)
      call calculate_overlap_matrix(dg_frag, system, mg)
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "        Momentum matrix calculated (for A·p coupling)"
        write(*,*) "        Overlap matrix S calculated (for generalized propagation)"
      end if
    else
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Momentum matrix already available"
      end if
      if (.not. allocated(dg_frag%S_mat_blocks)) call calculate_overlap_matrix(dg_frag, system, mg)
    end if
    
    ! Step 2: Allocate rank-local Hamiltonian blocks
    write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_c)) deallocate(dg_frag%H_mat_c)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                            diagonal_only=.false.)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                            diagonal_only=.false.)
    call rebuild_local_h_block_ids(dg_frag)
    
    ! Halo exchange removed: stencil operations use local phi_frag with fragment PBC buffer only.
    
    hvol = system%hvol
    is = mg%is
    ie = mg%ie
    
    ! Allocate work arrays
    allocate(T_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(H_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
      
      ! Loop over fragments assigned to this rank
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        
        ! Calculate Hamiltonian matrix elements for this fragment
        ! H_ij = <φ_i | T + V | φ_j> = T_ij + V_ij
        nbf = dg_frag%n_basis(ifrag, ispin)
        do jo = 1, nbf
          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi, &
                                    ispin_comp_opt=ispin)

          ! Calculate matrix elements with all φ_i
          iblk = 0
          if (allocated(dg_frag%H_block_map)) iblk = dg_frag%H_block_map(ifrag, ifrag)
          if (iblk <= 0) cycle
          !$omp parallel do private(io, integral_t, integral_h, ig_i)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            ! Kinetic energy matrix element: T_ij = ∫ φ_i (T|φ_j>) dr
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, integral_t, &
                                            ispin_comp_opt=ispin)

            ! Store kinetic part.
            ! The work domain is clipped to mg%is:mg%ie, so each rank can carry
            ! a distinct local-box contribution when nproc_rgrid > 1.
            dg_frag%H_mat_kinetic_blocks(iblk)%val(io, jo, ispin) = real(integral_t, kind=8)
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, integral_h, &
                                            ispin_comp_opt=ispin)
            dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = real(integral_h, kind=8)

          end do
          !$omp end parallel do

        end do  ! jo loop
          
        
      end do  ! ifrag loop
      
    end do  ! ispin loop
    
    call add_dg_surface_flux_blocks_soi(dg_frag, system, mg, stencil, &
                                        dg_frag%H_mat_blocks, dg_frag%H_mat_kinetic_blocks, &
                                        dg_frag%H_block_map)

    ! Global Hamiltonian aggregation via fragment blocks to avoid a monolithic dense allreduce.
    ! In orbital mode each subgroup rank owns a row slice of the same fragment
    ! row, so reducing over the global communicator would duplicate or mismatch
    ! row-local neighbor blocks.
    if (.not. dg_frag%parallel_mode_orbital) then
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-soi", dg_frag%icomm)
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, "hmat-kinetic-soi", dg_frag%icomm)
    end if
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_blocks)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks)
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Kinetic and potential terms computed"
    end if
    
    ! Step 3: Non-local pseudopotential is handled in time evolution
    ! with vector potential A(t), so it is not added to H_mat here.
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [3/3] Non-local PP handled in time evolution (A-dependent)"
    end if

    deallocate(T_phi, H_phi, V_total)
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "=== Hamiltonian Matrix Ready ==="
      write(*,*)
    end if
    
  end subroutine calculate_hamiltonian_matrix

  integer function map_global_to_phi_box_coord_ham_soi(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord_ham_soi

  integer function map_global_to_fragment_phi_box_coord_ham_soi(dg_frag, ifrag, axis, ig, lb, ub) result(iloc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, ig, lb, ub
    integer :: core_lo, ndom, ig_wrap

    core_lo = dg_frag%ixyz_frag(axis, ifrag)
    ndom = dg_frag%nxyz_domain(axis, ifrag)
    if (ndom <= 0) then
      iloc = 0
      return
    end if
    ig_wrap = core_lo + modulo(ig - core_lo, ndom)
    ig_wrap = modulo(ig_wrap - 1, dg_frag%lgnum_total(axis)) + 1
    iloc = map_global_to_phi_box_coord_ham_soi(ig_wrap, lb, ub, dg_frag%lgnum_total(axis))
  end function map_global_to_fragment_phi_box_coord_ham_soi

  subroutine enforce_fragment_periodic_buffer_for_state_ham_soi(dg_frag, ifrag, i_local, jo)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo

    integer :: ix, iy, iz
    integer :: sx, sy, sz
    integer :: lb1, ub1, lb2, ub2, lb3, ub3

    lb1 = lbound(dg_frag%phi_frag, 1)
    ub1 = ubound(dg_frag%phi_frag, 1)
    lb2 = lbound(dg_frag%phi_frag, 2)
    ub2 = ubound(dg_frag%phi_frag, 2)
    lb3 = lbound(dg_frag%phi_frag, 3)
    ub3 = ubound(dg_frag%phi_frag, 3)

    do iz = lb3, ub3
      sz = map_global_to_fragment_phi_box_coord_ham_soi(dg_frag, ifrag, 3, iz, lb3, ub3)
      if (sz == 0) cycle
      do iy = lb2, ub2
        sy = map_global_to_fragment_phi_box_coord_ham_soi(dg_frag, ifrag, 2, iy, lb2, ub2)
        if (sy == 0) cycle
        do ix = lb1, ub1
          sx = map_global_to_fragment_phi_box_coord_ham_soi(dg_frag, ifrag, 1, ix, lb1, ub1)
          if (sx == 0) cycle
          if (allocated(dg_frag%phi_frag_c)) then
            dg_frag%phi_frag_c(ix, iy, iz, jo, i_local) = dg_frag%phi_frag_c(sx, sy, sz, jo, i_local)
          end if
          dg_frag%phi_frag(ix, iy, iz, jo, i_local) = dg_frag%phi_frag(sx, sy, sz, jo, i_local)
        end do
      end do
    end do
  end subroutine enforce_fragment_periodic_buffer_for_state_ham_soi

  !=======================================================================
  ! Build total local potential on the given grid:
  !   V_total = Vpsl + Vh + Vxc(ispin)
  !=======================================================================
  subroutine build_total_potential_grid(grid, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz
    integer :: grid_x_lo, grid_x_hi, grid_y_lo, grid_y_hi, grid_z_lo, grid_z_hi

    grid_x_lo = grid%is(1)
    grid_x_hi = grid%ie(1)
    grid_y_lo = grid%is(2)
    grid_y_hi = grid%ie(2)
    grid_z_lo = grid%is(3)
    grid_z_hi = grid%ie(3)
!$omp parallel do collapse(2) private(ix,iy,iz)
    do iz = grid_z_lo, grid_z_hi
      do iy = grid_y_lo, grid_y_hi
!$omp simd
        do ix = grid_x_lo, grid_x_hi
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_total_potential_grid

  subroutine get_halo_block_point_indices(halo, ix_buf, iy_buf, iz_buf, send_idx, recv_idx)
    implicit none
    type(halo_info), intent(in) :: halo
    integer, intent(in) :: ix_buf, iy_buf, iz_buf
    integer, intent(out) :: send_idx(3), recv_idx(3)

    send_idx(1) = halo%send_lo(1) + ix_buf - 1
    send_idx(2) = halo%send_lo(2) + iy_buf - 1
    send_idx(3) = halo%send_lo(3) + iz_buf - 1
    recv_idx(1) = halo%recv_lo(1) + ix_buf - 1
    recv_idx(2) = halo%recv_lo(2) + iy_buf - 1
    recv_idx(3) = halo%recv_lo(3) + iz_buf - 1
  end subroutine get_halo_block_point_indices

  !=======================================================================
  ! Build T|phi_j> and H|phi_j>=T|phi_j>+V|phi_j> for one fragment/basis state
  !=======================================================================
  subroutine build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi, ispin_comp_opt)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    integer, intent(in), optional :: ispin_comp_opt
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: gx, gy, gz, bx, by, bz, ispin_comp
    integer :: iorg(3), ndom(3)
    integer :: gx_lo, gx_hi, gy_lo, gy_hi, gz_lo, gz_hi
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3

    ispin_comp = 1
    if (present(ispin_comp_opt)) ispin_comp = ispin_comp_opt

    call apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi, ispin_comp_opt=ispin_comp)
    H_phi(:, :, :) = T_phi(:, :, :)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    gx_lo = max(iorg(1), mg%is(1))
    gx_hi = min(iorg(1) + ndom(1) - 1, mg%ie(1))
    gy_lo = max(iorg(2), mg%is(2))
    gy_hi = min(iorg(2) + ndom(2) - 1, mg%ie(2))
    gz_lo = max(iorg(3), mg%is(3))
    gz_hi = min(iorg(3) + ndom(3) - 1, mg%ie(3))
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    if (gx_lo > gx_hi .or. gy_lo > gy_hi .or. gz_lo > gz_hi) return
!$omp parallel do private(gz, gy, gx, bx, by, bz) schedule(static)
    do gz = gz_lo, gz_hi
      bz = map_global_to_phi_box_coord_ham_soi(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
      if (bz == 0) cycle
      do gy = gy_lo, gy_hi
        by = map_global_to_phi_box_coord_ham_soi(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
        if (by == 0) cycle
!$omp simd
        do gx = gx_lo, gx_hi
          bx = map_global_to_phi_box_coord_ham_soi(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
          if (bx == 0) cycle
          if (allocated(dg_frag%phi_frag_spinor_c) .and. ispin_comp >= 1 .and. &
              ispin_comp <= size(dg_frag%phi_frag_spinor_c, 4)) then
            H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * &
              dg_frag%phi_frag_spinor_c(bx, by, bz, ispin_comp, jo, i_local)
          else if (allocated(dg_frag%phi_frag_c)) then
            H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + &
              V_total(gx, gy, gz) * dg_frag%phi_frag_c(bx, by, bz, jo, i_local)
          else
            H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * &
              cmplx(dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
          end if
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_hpsi_for_basis

  !=======================================================================
  ! Integrate one bra basis function against a real-space field
  !   integral = <phi_io | field>
  !=======================================================================
  subroutine integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, field, hvol, integral, ispin_comp_opt)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    integer, intent(in), optional :: ispin_comp_opt
    type(s_rgrid), intent(in) :: mg
    complex(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral
    integer :: gx, gy, gz, bx, by, bz, ispin_comp
    integer :: iorg(3), ndom(3)
    integer :: gx_lo, gx_hi, gy_lo, gy_hi, gz_lo, gz_hi
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3

    ispin_comp = 1
    if (present(ispin_comp_opt)) ispin_comp = ispin_comp_opt

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    gx_lo = max(iorg(1), mg%is(1))
    gx_hi = min(iorg(1) + ndom(1) - 1, mg%ie(1))
    gy_lo = max(iorg(2), mg%is(2))
    gy_hi = min(iorg(2) + ndom(2) - 1, mg%ie(2))
    gz_lo = max(iorg(3), mg%is(3))
    gz_hi = min(iorg(3) + ndom(3) - 1, mg%ie(3))
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    integral = (0.0d0, 0.0d0)
    if (gx_lo > gx_hi .or. gy_lo > gy_hi .or. gz_lo > gz_hi) return
    do gz = gz_lo, gz_hi
      bz = map_global_to_phi_box_coord_ham_soi(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
      if (bz == 0) cycle
      do gy = gy_lo, gy_hi
        by = map_global_to_phi_box_coord_ham_soi(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
        if (by == 0) cycle
        do gx = gx_lo, gx_hi
          bx = map_global_to_phi_box_coord_ham_soi(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
          if (bx == 0) cycle
          if (allocated(dg_frag%phi_frag_spinor_c) .and. ispin_comp >= 1 .and. &
              ispin_comp <= size(dg_frag%phi_frag_spinor_c, 4)) then
            integral = integral + conjg(dg_frag%phi_frag_spinor_c(bx, by, bz, ispin_comp, io, i_local)) * &
                       field(gx, gy, gz) * hvol
          else if (allocated(dg_frag%phi_frag_c)) then
            integral = integral + conjg(dg_frag%phi_frag_c(bx, by, bz, io, i_local)) * field(gx, gy, gz) * hvol
          else
            integral = integral + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end if
        end do
      end do
    end do
  end subroutine integrate_basis_with_field
  
  !=======================================================================
  ! Apply kinetic energy operator to a single basis function
  !
  ! T|φ> = -∇²/2 |φ> = -0.5 * Laplacian(φ)
  !
  ! Uses 4th-order finite difference stencil (requires ±4 grid points).
  ! With halo exchange, computation is valid over entire domain (1:nx, 1:ny, 1:nz).
  !
  ! System boundaries: PERIODIC (full system is periodic)
  ! Fragment boundaries: Halo exchange provides neighbor data via MPI
  !=======================================================================
  subroutine apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi, ispin_comp_opt)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer,                intent(in) :: i_local, jo
    integer,                intent(in), optional :: ispin_comp_opt
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    complex(8),             intent(out) :: T_phi(mg%is(1):mg%ie(1), &
                                                 mg%is(2):mg%ie(2), &
                                                 mg%is(3):mg%ie(3))
    
    integer :: gx, gy, gz, ifrag, ispin_comp
    integer :: ix0, iy0, iz0
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    integer :: lgx, lgy, lgz
    complex(8) :: v
    real(8) :: lap0
    real(8) :: lapt(4,3)
    integer :: is(3), ie(3), iorg(3), ndom(3)
    integer :: gx_lo, gx_hi, gy_lo, gy_hi, gz_lo, gz_hi
    logical :: use_complex_phi
    
    ! Extract stencil coefficients
    ispin_comp = 1
    if (present(ispin_comp_opt)) ispin_comp = ispin_comp_opt

    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    is = mg%is
    ie = mg%ie
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
      write(*,*) "[FATAL] SOI kinetic invalid i_local: rank=", dg_frag%id, &
        " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
      stop 1
    end if
    if (jo < 1 .or. jo > size(dg_frag%phi_frag, 4)) then
      write(*,*) "[FATAL] SOI kinetic invalid jo: rank=", dg_frag%id, &
        " jo=", jo, " phi_frag_dim4=", size(dg_frag%phi_frag, 4)
      stop 1
    end if
    ifrag = dg_frag%ifrag_start + i_local - 1
    if (ifrag < dg_frag%ifrag_start .or. ifrag > dg_frag%ifrag_end) then
      write(*,*) "[FATAL] SOI kinetic invalid ifrag from i_local: rank=", dg_frag%id, &
        " ifrag=", ifrag, " i_local=", i_local, " ifrag_start/end=", dg_frag%ifrag_start, dg_frag%ifrag_end
      stop 1
    end if
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    if (any(ndom <= 0)) then
      write(*,*) "[FATAL] SOI kinetic non-positive domain size: rank=", dg_frag%id, &
        " ifrag=", ifrag, " ndom=", ndom
      stop 1
    end if
    if (.not. allocated(dg_frag%phi_frag_has_seed_buffer) .or. &
        .not. dg_frag%phi_frag_has_seed_buffer(i_local)) then
      call enforce_fragment_periodic_buffer_for_state_ham_soi(dg_frag, ifrag, i_local, jo)
    end if
    
    gx_lo = max(iorg(1), mg%is(1))
    gx_hi = min(iorg(1) + ndom(1) - 1, mg%ie(1))
    gy_lo = max(iorg(2), mg%is(2))
    gy_hi = min(iorg(2) + ndom(2) - 1, mg%ie(2))
    gz_lo = max(iorg(3), mg%is(3))
    gz_hi = min(iorg(3) + ndom(3) - 1, mg%ie(3))
    lgx = dg_frag%lgnum_total(1)
    lgy = dg_frag%lgnum_total(2)
    lgz = dg_frag%lgnum_total(3)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    T_phi = (0.0d0, 0.0d0)
    use_complex_phi = allocated(dg_frag%phi_frag_c)
    
    if (gx_lo > gx_hi .or. gy_lo > gy_hi .or. gz_lo > gz_hi) return
    if (allocated(dg_frag%phi_frag_spinor_c) .and. ispin_comp >= 1 .and. &
        ispin_comp <= size(dg_frag%phi_frag_spinor_c, 4)) then
!$omp parallel do private(gz, gy, gx, v, ix0, iy0, iz0) schedule(static)
      do gz = gz_lo, gz_hi
        iz0 = map_global_to_phi_box_coord_ham_soi(gz, p_lb3, p_ub3, lgz)
        do gy = gy_lo, gy_hi
          iy0 = map_global_to_phi_box_coord_ham_soi(gy, p_lb2, p_ub2, lgy)
!$omp simd private(v)
          do gx = gx_lo, gx_hi
            ix0 = map_global_to_phi_box_coord_ham_soi(gx, p_lb1, p_ub1, lgx)

            v = lapt(1,1) * (dg_frag%phi_frag_spinor_c(ix0 + 1, iy0, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0 - 1, iy0, iz0, ispin_comp, jo, i_local)) + &
                lapt(2,1) * (dg_frag%phi_frag_spinor_c(ix0 + 2, iy0, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0 - 2, iy0, iz0, ispin_comp, jo, i_local)) + &
                lapt(3,1) * (dg_frag%phi_frag_spinor_c(ix0 + 3, iy0, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0 - 3, iy0, iz0, ispin_comp, jo, i_local)) + &
                lapt(4,1) * (dg_frag%phi_frag_spinor_c(ix0 + 4, iy0, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0 - 4, iy0, iz0, ispin_comp, jo, i_local))
            v = v + &
                lapt(1,2) * (dg_frag%phi_frag_spinor_c(ix0, iy0 + 1, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0 - 1, iz0, ispin_comp, jo, i_local)) + &
                lapt(2,2) * (dg_frag%phi_frag_spinor_c(ix0, iy0 + 2, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0 - 2, iz0, ispin_comp, jo, i_local)) + &
                lapt(3,2) * (dg_frag%phi_frag_spinor_c(ix0, iy0 + 3, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0 - 3, iz0, ispin_comp, jo, i_local)) + &
                lapt(4,2) * (dg_frag%phi_frag_spinor_c(ix0, iy0 + 4, iz0, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0 - 4, iz0, ispin_comp, jo, i_local))
            v = v + &
                lapt(1,3) * (dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 + 1, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 - 1, ispin_comp, jo, i_local)) + &
                lapt(2,3) * (dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 + 2, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 - 2, ispin_comp, jo, i_local)) + &
                lapt(3,3) * (dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 + 3, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 - 3, ispin_comp, jo, i_local)) + &
                lapt(4,3) * (dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 + 4, ispin_comp, jo, i_local) + &
                             dg_frag%phi_frag_spinor_c(ix0, iy0, iz0 - 4, ispin_comp, jo, i_local))
            T_phi(gx, gy, gz) = lap0 * dg_frag%phi_frag_spinor_c(ix0, iy0, iz0, ispin_comp, jo, i_local) - 0.5d0 * v
          end do
        end do
      end do
!$omp end parallel do
    else if (use_complex_phi) then
!$omp parallel do private(gz, gy, gx, v, ix0, iy0, iz0) schedule(static)
      do gz = gz_lo, gz_hi
        iz0 = map_global_to_phi_box_coord_ham_soi(gz, p_lb3, p_ub3, lgz)
        do gy = gy_lo, gy_hi
          iy0 = map_global_to_phi_box_coord_ham_soi(gy, p_lb2, p_ub2, lgy)
!$omp simd private(v)
          do gx = gx_lo, gx_hi
            ix0 = map_global_to_phi_box_coord_ham_soi(gx, p_lb1, p_ub1, lgx)

            v = lapt(1,1) * (dg_frag%phi_frag_c(ix0 + 1, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0 - 1, iy0, iz0, jo, i_local)) + &
                lapt(2,1) * (dg_frag%phi_frag_c(ix0 + 2, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0 - 2, iy0, iz0, jo, i_local)) + &
                lapt(3,1) * (dg_frag%phi_frag_c(ix0 + 3, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0 - 3, iy0, iz0, jo, i_local)) + &
                lapt(4,1) * (dg_frag%phi_frag_c(ix0 + 4, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0 - 4, iy0, iz0, jo, i_local))
            v = v + &
                lapt(1,2) * (dg_frag%phi_frag_c(ix0, iy0 + 1, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0 - 1, iz0, jo, i_local)) + &
                lapt(2,2) * (dg_frag%phi_frag_c(ix0, iy0 + 2, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0 - 2, iz0, jo, i_local)) + &
                lapt(3,2) * (dg_frag%phi_frag_c(ix0, iy0 + 3, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0 - 3, iz0, jo, i_local)) + &
                lapt(4,2) * (dg_frag%phi_frag_c(ix0, iy0 + 4, iz0, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0 - 4, iz0, jo, i_local))
            v = v + &
                lapt(1,3) * (dg_frag%phi_frag_c(ix0, iy0, iz0 + 1, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0, iz0 - 1, jo, i_local)) + &
                lapt(2,3) * (dg_frag%phi_frag_c(ix0, iy0, iz0 + 2, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0, iz0 - 2, jo, i_local)) + &
                lapt(3,3) * (dg_frag%phi_frag_c(ix0, iy0, iz0 + 3, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0, iz0 - 3, jo, i_local)) + &
                lapt(4,3) * (dg_frag%phi_frag_c(ix0, iy0, iz0 + 4, jo, i_local) + &
                             dg_frag%phi_frag_c(ix0, iy0, iz0 - 4, jo, i_local))
            T_phi(gx, gy, gz) = lap0 * dg_frag%phi_frag_c(ix0, iy0, iz0, jo, i_local) - 0.5d0 * v
          end do
        end do
      end do
!$omp end parallel do
    else
!$omp parallel do private(gz, gy, gx, v, ix0, iy0, iz0) schedule(static)
      do gz = gz_lo, gz_hi
        iz0 = map_global_to_phi_box_coord_ham_soi(gz, p_lb3, p_ub3, lgz)
        do gy = gy_lo, gy_hi
          iy0 = map_global_to_phi_box_coord_ham_soi(gy, p_lb2, p_ub2, lgy)
!$omp simd private(v)
          do gx = gx_lo, gx_hi
            ix0 = map_global_to_phi_box_coord_ham_soi(gx, p_lb1, p_ub1, lgx)

            v = cmplx(lapt(1,1) * (dg_frag%phi_frag(ix0 + 1, iy0, iz0, jo, i_local) + &
                                   dg_frag%phi_frag(ix0 - 1, iy0, iz0, jo, i_local)) + &
                      lapt(2,1) * (dg_frag%phi_frag(ix0 + 2, iy0, iz0, jo, i_local) + &
                                   dg_frag%phi_frag(ix0 - 2, iy0, iz0, jo, i_local)) + &
                      lapt(3,1) * (dg_frag%phi_frag(ix0 + 3, iy0, iz0, jo, i_local) + &
                                   dg_frag%phi_frag(ix0 - 3, iy0, iz0, jo, i_local)) + &
                      lapt(4,1) * (dg_frag%phi_frag(ix0 + 4, iy0, iz0, jo, i_local) + &
                                   dg_frag%phi_frag(ix0 - 4, iy0, iz0, jo, i_local)), 0.0d0, kind=8)
            v = v + cmplx(lapt(1,2) * (dg_frag%phi_frag(ix0, iy0 + 1, iz0, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0 - 1, iz0, jo, i_local)) + &
                          lapt(2,2) * (dg_frag%phi_frag(ix0, iy0 + 2, iz0, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0 - 2, iz0, jo, i_local)) + &
                          lapt(3,2) * (dg_frag%phi_frag(ix0, iy0 + 3, iz0, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0 - 3, iz0, jo, i_local)) + &
                          lapt(4,2) * (dg_frag%phi_frag(ix0, iy0 + 4, iz0, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0 - 4, iz0, jo, i_local)), 0.0d0, kind=8)
            v = v + cmplx(lapt(1,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 1, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0, iz0 - 1, jo, i_local)) + &
                          lapt(2,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 2, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0, iz0 - 2, jo, i_local)) + &
                          lapt(3,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 3, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0, iz0 - 3, jo, i_local)) + &
                          lapt(4,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 4, jo, i_local) + &
                                       dg_frag%phi_frag(ix0, iy0, iz0 - 4, jo, i_local)), 0.0d0, kind=8)
            T_phi(gx, gy, gz) = cmplx(lap0 * dg_frag%phi_frag(ix0, iy0, iz0, jo, i_local), 0.0d0, kind=8) - 0.5d0 * v
          end do
        end do
      end do
!$omp end parallel do
    end if
    
  end subroutine apply_kinetic_to_basis

  subroutine get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_rgrid), intent(in) :: mg
    integer, intent(out) :: loc_s(3), loc_e(3)
    logical, intent(out) :: has_overlap

    integer :: iorg(3), ndom(3), g_s(3), g_e(3), ov_s(3), ov_e(3)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    if (dg_frag%parallel_mode_orbital) then
      loc_s(:) = 1
      loc_e(:) = ndom(:)
      has_overlap = all(ndom(:) > 0)
      return
    end if
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))

    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) then
      loc_s(:) = 1
      loc_e(:) = 0
      return
    end if

    loc_s(:) = ov_s(:) - iorg(:) + 1
    loc_e(:) = ov_e(:) - iorg(:) + 1
  end subroutine get_fragment_owned_range

  subroutine get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ndom(3)
    integer, intent(out) :: loc_s(3), loc_e(3)

    integer :: ipx, ipy, ipz, coords(3), nsize

    ipx = max(1, nproc_rgrid(1))
    ipy = max(1, nproc_rgrid(2))
    ipz = max(1, nproc_rgrid(3))

    if (dg_frag%id_frag < 0 .or. dg_frag%id_frag >= ipx * ipy * ipz) then
      stop "DG-Fragment RT: invalid fragment-local MPI rank in get_fragment_local_range"
    end if

    coords(1) = mod(dg_frag%id_frag, ipx)
    coords(2) = mod(dg_frag%id_frag / ipx, ipy)
    coords(3) = dg_frag%id_frag / max(1, ipx * ipy)

    nsize = (ndom(1) + ipx - 1) / ipx
    loc_s(1) = 1 + nsize * coords(1)
    loc_e(1) = min(ndom(1), loc_s(1) + nsize - 1)

    nsize = (ndom(2) + ipy - 1) / ipy
    loc_s(2) = 1 + nsize * coords(2)
    loc_e(2) = min(ndom(2), loc_s(2) + nsize - 1)

    nsize = (ndom(3) + ipz - 1) / ipz
    loc_s(3) = 1 + nsize * coords(3)
    loc_e(3) = min(ndom(3), loc_s(3) + nsize - 1)
  end subroutine get_fragment_local_range

  !=======================================================================
  ! Calculate momentum matrix elements in fragment basis (velocity gauge)
  !=======================================================================
  subroutine calculate_momentum_matrix(dg_frag, system, mg, stencil)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir, iblk
    integer :: lx, ly, lz, iorg(3), loc_s(3), loc_e(3), phi_loc_s(3), phi_loc_e(3)
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: nrow
    integer :: gidx(3)
    real(8) :: hvol, max_p, momentum_gb, momentum_block_gb
    complex(8) :: bra_val, grad_vec(3), integral_vec(3)
    logical :: has_overlap
    
    if (.not. dg_frag%has_real_space_basis) return

    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Computing transition moments: <φ_i|∇|φ_j>"
    end if
    
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    if (allocated(dg_frag%momentum_mat_c)) deallocate(dg_frag%momentum_mat_c)
    call init_momentum_blocks(dg_frag, diagonal_only=(.not. dg_frag%parallel_mode_orbital))
    momentum_gb = real(3_8 * int(dg_frag%n_mat_max, kind=8) * int(dg_frag%n_mat_max, kind=8) * &
      int(dg_frag%nspin, kind=8) * 16_8, 8) / 1.0d9
    momentum_block_gb = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      momentum_block_gb = momentum_block_gb + real(3_8 * int(dg_frag%momentum_blocks(iblk)%nrow_max, kind=8) * &
        int(dg_frag%momentum_blocks(iblk)%ncol_max, kind=8) * int(dg_frag%nspin, kind=8) * 16_8, 8) / 1.0d9
    end do
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,f10.3,a)') "        n_mat_max=", dg_frag%n_mat_max, &
        " dense complex momentum_mat GB=", momentum_gb, " (not allocated)"
      write(*,'(1x,a,i0,a,f10.6,a)') "        momentum_blocks=", dg_frag%n_momentum_blocks, &
        " complex allocated GB=", momentum_block_gb, " per rank"
      flush(6)
    end if
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    hvol = system%hvol
    
    ! Loop over spin
    do ispin = 1, system%nspin
      ! Loop over local fragments
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ! In orbital mode every rank integrates the full fragment box but only
        ! stores its owned orbital rows.  Legacy mode keeps the parent-grid box.
        call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
        if (.not. has_overlap) cycle
        phi_loc_s(:) = iorg(:) + loc_s(:) - 1
        phi_loc_e(:) = iorg(:) + loc_e(:) - 1
        if (phi_loc_s(1) < phi_lb1 .or. phi_loc_e(1) > phi_ub1 .or. &
            phi_loc_s(2) < phi_lb2 .or. phi_loc_e(2) > phi_ub2 .or. &
            phi_loc_s(3) < phi_lb3 .or. phi_loc_e(3) > phi_ub3) then
          write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x))') &
            "DG-Fragment RT SOI momentum phi_frag local range out of bounds: phi_loc_s=", &
            phi_loc_s(1), phi_loc_s(2), phi_loc_s(3), "phi_loc_e=", phi_loc_e(1), phi_loc_e(2), phi_loc_e(3), &
            "lb=", phi_lb1, phi_lb2, phi_lb3, "ub=", phi_ub1, phi_ub2, phi_ub3
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
            "        soi momentum context: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, " ifrag=", ifrag, &
            " i_local=", i_local, " mg_is=", mg%is(1), mg%is(2), mg%is(3), " mg_ie=", mg%ie(1), mg%ie(2), mg%ie(3), &
            " loc_s=", loc_s(1), loc_s(2), loc_s(3), " loc_e=", loc_e(1), loc_e(2), loc_e(3)
          call flush(6)
          stop "DG-Fragment RT SOI: momentum phi_frag local range out of bounds"
        end if
        
        iblk = find_momentum_block(dg_frag, ifrag, ifrag)
        if (iblk <= 0) cycle
        nrow = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
        if (nrow <= 0) cycle
        do jo = 1, nrow
          do io = 1, nrow
            if (dg_frag%parallel_mode_orbital) then
              if (.not. basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, io)) cycle
            end if
            integral_vec(:) = (0.0d0, 0.0d0)
            do lz = phi_loc_s(3), phi_loc_e(3)
              do ly = phi_loc_s(2), phi_loc_e(2)
                do lx = phi_loc_s(1), phi_loc_e(1)
                  gidx = [lx, ly, lz]
                  bra_val = conjg(phi_local_value_ham_soi(dg_frag, i_local, ispin, io, gidx))
                  call phi_local_grad_ham_soi(dg_frag, i_local, stencil, ispin, jo, gidx, grad_vec)
                  integral_vec(:) = integral_vec(:) + bra_val * grad_vec(:) * hvol
                end do
              end do
            end do
            do idir = 1, 3
              dg_frag%momentum_blocks_c(iblk)%val(idir, io, jo, ispin) = integral_vec(idir)
              dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) = real(integral_vec(idir), kind=8)
            end do
          end do
        end do

      end do  ! ifrag
    end do  ! ispin

    call add_dg_surface_momentum_blocks_soi(dg_frag, system)

    max_p = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      max_p = max(max_p, maxval(abs(dg_frag%momentum_blocks(iblk)%val)))
      if (allocated(dg_frag%momentum_blocks_c)) then
        max_p = max(max_p, maxval(abs(dg_frag%momentum_blocks_c(iblk)%val)))
      end if
    end do
    if (comm_is_root(dg_frag%id)) then
      write(*,'(a,1pe12.4)') "        Max |momentum_mat|: ", max_p
      write(*,'(a,i0,a)') "        Momentum blocks: ", dg_frag%n_momentum_blocks, &
                          " (fragment-local storage)"
    end if

  end subroutine calculate_momentum_matrix

  integer function find_matrix_block(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block

  logical function matrix_row_is_locally_owned(dg_frag, row_idx, ispin) result(is_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_idx, ispin

    is_local = .false.
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > size(dg_frag%coef_owner, 2)) return
    if (row_idx < 1 .or. row_idx > size(dg_frag%coef_owner, 1)) return
    is_local = (dg_frag%coef_owner(row_idx, ispin) == dg_frag%id)
  end function matrix_row_is_locally_owned

  subroutine init_matrix_blocks(dg_frag, blocks, block_map, n_blocks, diagonal_only)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    logical, intent(in), optional :: diagonal_only
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max
    logical :: diagonal_blocks_only, local_fragment_only, include_pair

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    ! Keep SOI volume Galerkin operators fragment-local as in the non-SOI path.
    ! Spinor structure changes matrix values, not the block topology.
    diagonal_blocks_only = .true.
    if (present(diagonal_only)) diagonal_blocks_only = diagonal_only
    local_fragment_only = dg_frag%parallel_mode_orbital

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (include_pair) n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (.not. include_pair) cycle
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk
        blocks(iblk)%ifrag_row = ifrag_row
        blocks(iblk)%ifrag_col = ifrag_col
        blocks(iblk)%nrow_max = nrow_max
        blocks(iblk)%ncol_max = ncol_max
        allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_matrix_blocks


  logical function is_momentum_neighbor_axis(lg, s1, n1, s2, n2) result(ok)
    implicit none
    integer, intent(in) :: lg, s1, n1, s2, n2
    integer :: e1, e2, s1_next, s2_next

    e1 = s1 + n1 - 1
    e2 = s2 + n2 - 1
    s1_next = modulo(e1, lg) + 1
    s2_next = modulo(e2, lg) + 1
    ok = ((s1 == s2) .and. (n1 == n2)) .or. (s1 == s2_next) .or. (s2 == s1_next)
  end function is_momentum_neighbor_axis

  logical function is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col) result(is_pair)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col
    integer :: axis, n_same, n_adjacent
    logical :: same_axis, adjacent_axis

    is_pair = .false.
    if (ifrag_row == ifrag_col) then
      is_pair = .true.
      return
    end if
    n_same = 0
    n_adjacent = 0
    do axis = 1, 3
      same_axis = (dg_frag%ixyz_frag(axis, ifrag_row) == dg_frag%ixyz_frag(axis, ifrag_col) .and. &
                   dg_frag%nxyz_domain(axis, ifrag_row) == dg_frag%nxyz_domain(axis, ifrag_col))
      adjacent_axis = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col)) .and. .not. same_axis
      if (same_axis) n_same = n_same + 1
      if (adjacent_axis) n_adjacent = n_adjacent + 1
    end do

    is_pair = (n_same == 2 .and. n_adjacent == 1)
  end function is_momentum_neighbor_pair

  subroutine ensure_momentum_neighbor_pair_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col

    if (allocated(dg_frag%momentum_neighbor_pair_cache)) return
    allocate(dg_frag%momentum_neighbor_pair_cache(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_neighbor_pair_cache(:, :) = .false.
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (ifrag_row == ifrag_col) then
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = .true.
        else
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = &
            is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
      end do
    end do
  end subroutine ensure_momentum_neighbor_pair_cache

  integer function find_momentum_block(dg_frag, ifrag_row, ifrag_col) result(iblk)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (.not. allocated(dg_frag%momentum_block_map)) return
    if (ifrag_row < 1 .or. ifrag_row > size(dg_frag%momentum_block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(dg_frag%momentum_block_map, 2)) return
    iblk = dg_frag%momentum_block_map(ifrag_row, ifrag_col)
  end function find_momentum_block

  subroutine init_momentum_blocks(dg_frag, diagonal_only)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, intent(in), optional :: diagonal_only
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max
    logical :: diagonal_blocks_only, local_fragment_only, include_pair

    if (allocated(dg_frag%momentum_blocks)) then
      do iblk = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(iblk)%val)) deallocate(dg_frag%momentum_blocks(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
    end if
    if (allocated(dg_frag%momentum_blocks_c)) then
      do iblk = 1, size(dg_frag%momentum_blocks_c)
        if (allocated(dg_frag%momentum_blocks_c(iblk)%val)) deallocate(dg_frag%momentum_blocks_c(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks_c)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    diagonal_blocks_only = .true.
    if (present(diagonal_only)) diagonal_blocks_only = diagonal_only
    local_fragment_only = dg_frag%parallel_mode_orbital

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (include_pair) nblk = nblk + 1
      end do
    end do

    dg_frag%n_momentum_blocks = nblk
    if (nblk <= 0) return
    allocate(dg_frag%momentum_blocks(nblk))
    allocate(dg_frag%momentum_blocks_c(nblk))
    allocate(dg_frag%momentum_block_map(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (.not. include_pair) cycle
        iblk = iblk + 1
        dg_frag%momentum_block_map(ifrag_row, ifrag_col) = iblk
        dg_frag%momentum_blocks(iblk)%ifrag_row = ifrag_row
        dg_frag%momentum_blocks(iblk)%ifrag_col = ifrag_col
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        dg_frag%momentum_blocks(iblk)%nrow_max = nrow_max
        dg_frag%momentum_blocks(iblk)%ncol_max = ncol_max
        allocate(dg_frag%momentum_blocks(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        dg_frag%momentum_blocks(iblk)%val = 0.0d0
        dg_frag%momentum_blocks_c(iblk)%ifrag_row = ifrag_row
        dg_frag%momentum_blocks_c(iblk)%ifrag_col = ifrag_col
        dg_frag%momentum_blocks_c(iblk)%nrow_max = nrow_max
        dg_frag%momentum_blocks_c(iblk)%ncol_max = ncol_max
        allocate(dg_frag%momentum_blocks_c(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        dg_frag%momentum_blocks_c(iblk)%val = (0.0d0, 0.0d0)
      end do
    end do
  end subroutine init_momentum_blocks

  integer function face_neighbor_fragment_ham_soi(dg_frag, ifrag, axis, side) result(jfrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: coords(3), rem

    jfrag = 0
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
    if (axis < 1 .or. axis > 3) return
    if (dg_frag%num_fragment(axis) <= 1) return

    coords(1) = (ifrag - 1) / (dg_frag%num_fragment(2) * dg_frag%num_fragment(3)) + 1
    rem = modulo(ifrag - 1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3))
    coords(2) = rem / dg_frag%num_fragment(3) + 1
    coords(3) = modulo(rem, dg_frag%num_fragment(3)) + 1
    coords(axis) = modulo(coords(axis) - 1 + side + dg_frag%num_fragment(axis), &
                          dg_frag%num_fragment(axis)) + 1
    jfrag = ((coords(1) - 1) * dg_frag%num_fragment(2) + coords(2) - 1) * &
            dg_frag%num_fragment(3) + coords(3)
  end function face_neighbor_fragment_ham_soi

  subroutine read_fragment_buffer_basis_box_ham_soi(dg_frag, ifrag, phi_box, box_lo, box_hi)
    use filesystem, only: get_filehandle
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    complex(8), allocatable, intent(out) :: phi_box(:,:,:,:,:)
    integer, intent(out) :: box_lo(3), box_hi(3)

    character(32), parameter :: bdir_frag = './data_dcdft/fragments/'
    character(32), parameter :: binfile_bfb = 'basis_functions_buffer_soi.bin'
    integer, parameter :: basis_buffer_magic = -22022213
    integer, parameter :: basis_buffer_version = 2
    character(256) :: filename
    integer :: iunit, magic_file, version_file
    integer :: nxyz_domain(3), nxyz_buffer_file(3), nxyz_box(3)
    integer :: nspin_file, nstate_frag_file, n_basis_keep, nspin_keep
    integer :: ispin_file, istate, iaxis, ixg, iyg, izg
    integer :: ix_box, iy_box, iz_box, ix_src, iy_src, iz_src
    integer, allocatable :: n_basis_frag(:)
    complex(8), allocatable :: phi_read(:,:,:)
    logical :: has_file

    if (allocated(phi_box)) deallocate(phi_box)
    box_lo(:) = 0
    box_hi(:) = -1
    write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bfb
    inquire(file=filename, exist=has_file)
    if (.not. has_file) then
      write(*,'(1x,a,i0,a,a)') '[FATAL] missing neighbor DG-SOI buffer basis at ifrag=', ifrag, &
        ' file=', trim(filename)
      stop 'DG-Fragment RT SOI: missing neighbor basis buffer file'
    end if

    iunit = get_filehandle()
    open(iunit, file=filename, form='unformatted', access='stream', status='old')
    read(iunit) magic_file, version_file
    if (magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
      write(*,'(1x,a,i0,4(a,i0))') '[FATAL] invalid neighbor SOI basis buffer header at ifrag=', &
        ifrag, ' magic=', magic_file, ' expected_magic=', basis_buffer_magic, &
        ' version=', version_file, ' expected_version=', basis_buffer_version
      stop 'DG-Fragment RT SOI: invalid neighbor basis buffer file'
    end if
    read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
    do iaxis = 1, 3
      if (dg_frag%num_fragment(iaxis) > 1 .and. nxyz_buffer_file(iaxis) < dg_frag%nxyz_buffer(iaxis)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] neighbor SOI buffer too small at ifrag=', &
          ifrag, ' axis=', iaxis, ' seed_buffer=', nxyz_buffer_file(iaxis), &
          ' required=', dg_frag%nxyz_buffer(iaxis)
        stop 'DG-Fragment RT SOI: insufficient neighbor basis buffer'
      end if
    end do

    allocate(n_basis_frag(max(1, nspin_file)))
    read(iunit) n_basis_frag(1:max(1, nspin_file))
    n_basis_keep = 0
    if (nspin_file >= 1) then
      n_basis_keep = min(maxval(dg_frag%n_basis(ifrag, 1:dg_frag%nspin)), nstate_frag_file, dg_frag%nstate_frag)
    end if
    nspin_keep = min(max(1, nspin_file), max(1, dg_frag%nspin))
    nxyz_box(:) = nxyz_domain(:) + 2 * nxyz_buffer_file(:)
    box_lo(:) = dg_frag%ixyz_frag(:, ifrag) - nxyz_buffer_file(:)
    box_hi(:) = dg_frag%ixyz_frag(:, ifrag) + nxyz_domain(:) - 1 + nxyz_buffer_file(:)
    allocate(phi_box(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3), nspin_keep, max(1, n_basis_keep)))
    allocate(phi_read(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3)))
    phi_box = (0.0d0, 0.0d0)

    do ispin_file = 1, nspin_file
      do istate = 1, nstate_frag_file
        read(iunit) phi_read(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3))
        if (ispin_file > nspin_keep .or. istate > n_basis_keep) cycle
        do izg = box_lo(3), box_hi(3)
          iz_box = izg - (dg_frag%ixyz_frag(3, ifrag) - nxyz_buffer_file(3)) + 1
          if (iz_box < 1 .or. iz_box > nxyz_box(3)) then
            iz_src = nxyz_buffer_file(3) + modulo(izg - dg_frag%ixyz_frag(3, ifrag), nxyz_domain(3)) + 1
          else
            iz_src = iz_box
          end if
          do iyg = box_lo(2), box_hi(2)
            iy_box = iyg - (dg_frag%ixyz_frag(2, ifrag) - nxyz_buffer_file(2)) + 1
            if (iy_box < 1 .or. iy_box > nxyz_box(2)) then
              iy_src = nxyz_buffer_file(2) + modulo(iyg - dg_frag%ixyz_frag(2, ifrag), nxyz_domain(2)) + 1
            else
              iy_src = iy_box
            end if
            do ixg = box_lo(1), box_hi(1)
              ix_box = ixg - (dg_frag%ixyz_frag(1, ifrag) - nxyz_buffer_file(1)) + 1
              if (ix_box < 1 .or. ix_box > nxyz_box(1)) then
                ix_src = nxyz_buffer_file(1) + modulo(ixg - dg_frag%ixyz_frag(1, ifrag), nxyz_domain(1)) + 1
              else
                ix_src = ix_box
              end if
              phi_box(ixg - box_lo(1) + 1, iyg - box_lo(2) + 1, izg - box_lo(3) + 1, ispin_file, istate) = &
                phi_read(ix_src, iy_src, iz_src)
            end do
          end do
        end do
      end do
    end do
    close(iunit)
    deallocate(phi_read, n_basis_frag)
  end subroutine read_fragment_buffer_basis_box_ham_soi

  complex(8) function phi_box_value_ham_soi(phi_box, box_lo, box_hi, lgtot, ispin_comp, istate, gidx) result(val)
    implicit none
    complex(8), intent(in) :: phi_box(:,:,:,:,:)
    integer, intent(in) :: box_lo(3), box_hi(3), lgtot(3), ispin_comp, istate, gidx(3)
    integer :: ix, iy, iz

    val = (0.0d0, 0.0d0)
    if (ispin_comp < 1 .or. ispin_comp > size(phi_box, 4)) return
    if (istate < 1 .or. istate > size(phi_box, 5)) return
    ix = map_global_to_phi_box_coord_ham_soi(gidx(1), box_lo(1), box_hi(1), lgtot(1))
    iy = map_global_to_phi_box_coord_ham_soi(gidx(2), box_lo(2), box_hi(2), lgtot(2))
    iz = map_global_to_phi_box_coord_ham_soi(gidx(3), box_lo(3), box_hi(3), lgtot(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    val = phi_box(ix - box_lo(1) + 1, iy - box_lo(2) + 1, iz - box_lo(3) + 1, ispin_comp, istate)
  end function phi_box_value_ham_soi

  complex(8) function phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, gidx) result(val)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, ispin_comp, istate, gidx(3)
    integer :: ix, iy, iz

    val = (0.0d0, 0.0d0)
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) return
    if (istate < 1 .or. istate > size(dg_frag%phi_frag, 4)) return
    ix = map_global_to_phi_box_coord_ham_soi(gidx(1), lbound(dg_frag%phi_frag, 1), &
                                            ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
    iy = map_global_to_phi_box_coord_ham_soi(gidx(2), lbound(dg_frag%phi_frag, 2), &
                                            ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
    iz = map_global_to_phi_box_coord_ham_soi(gidx(3), lbound(dg_frag%phi_frag, 3), &
                                            ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    if (allocated(dg_frag%phi_frag_spinor_c) .and. ispin_comp >= 1 .and. &
        ispin_comp <= size(dg_frag%phi_frag_spinor_c, 4)) then
      val = dg_frag%phi_frag_spinor_c(ix, iy, iz, ispin_comp, istate, i_local)
    else if (allocated(dg_frag%phi_frag_c)) then
      val = dg_frag%phi_frag_c(ix, iy, iz, istate, i_local)
    else
      val = cmplx(dg_frag%phi_frag(ix, iy, iz, istate, i_local), 0.0d0, kind=8)
    end if
  end function phi_local_value_ham_soi

  complex(8) function phi_box_dn_ham_soi(phi_box, box_lo, box_hi, lgtot, stencil, ispin_comp, istate, gidx, axis, side) result(dn)
    use structures, only: s_stencil
    implicit none
    complex(8), intent(in) :: phi_box(:,:,:,:,:)
    integer, intent(in) :: box_lo(3), box_hi(3), lgtot(3), ispin_comp, istate, gidx(3), axis, side
    type(s_stencil), intent(in) :: stencil
    integer :: ix, iy, iz, ixc, iyc, izc
    real(8) :: nabt(4,3)
    complex(8) :: grad_axis

    dn = (0.0d0, 0.0d0)
    if (ispin_comp < 1 .or. ispin_comp > size(phi_box, 4)) return
    if (istate < 1 .or. istate > size(phi_box, 5)) return
    ix = map_global_to_phi_box_coord_ham_soi(gidx(1), box_lo(1), box_hi(1), lgtot(1))
    iy = map_global_to_phi_box_coord_ham_soi(gidx(2), box_lo(2), box_hi(2), lgtot(2))
    iz = map_global_to_phi_box_coord_ham_soi(gidx(3), box_lo(3), box_hi(3), lgtot(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    ixc = ix - box_lo(1) + 1
    iyc = iy - box_lo(2) + 1
    izc = iz - box_lo(3) + 1
    nabt = stencil%coef_nab
    select case(axis)
    case(1)
      if (ixc - 4 < 1 .or. ixc + 4 > size(phi_box, 1)) return
      grad_axis = nabt(1,1) * (phi_box(ixc + 1, iyc, izc, ispin_comp, istate) - phi_box(ixc - 1, iyc, izc, ispin_comp, istate)) + &
                  nabt(2,1) * (phi_box(ixc + 2, iyc, izc, ispin_comp, istate) - phi_box(ixc - 2, iyc, izc, ispin_comp, istate)) + &
                  nabt(3,1) * (phi_box(ixc + 3, iyc, izc, ispin_comp, istate) - phi_box(ixc - 3, iyc, izc, ispin_comp, istate)) + &
                  nabt(4,1) * (phi_box(ixc + 4, iyc, izc, ispin_comp, istate) - phi_box(ixc - 4, iyc, izc, ispin_comp, istate))
    case(2)
      if (iyc - 4 < 1 .or. iyc + 4 > size(phi_box, 2)) return
      grad_axis = nabt(1,2) * (phi_box(ixc, iyc + 1, izc, ispin_comp, istate) - phi_box(ixc, iyc - 1, izc, ispin_comp, istate)) + &
                  nabt(2,2) * (phi_box(ixc, iyc + 2, izc, ispin_comp, istate) - phi_box(ixc, iyc - 2, izc, ispin_comp, istate)) + &
                  nabt(3,2) * (phi_box(ixc, iyc + 3, izc, ispin_comp, istate) - phi_box(ixc, iyc - 3, izc, ispin_comp, istate)) + &
                  nabt(4,2) * (phi_box(ixc, iyc + 4, izc, ispin_comp, istate) - phi_box(ixc, iyc - 4, izc, ispin_comp, istate))
    case(3)
      if (izc - 4 < 1 .or. izc + 4 > size(phi_box, 3)) return
      grad_axis = nabt(1,3) * (phi_box(ixc, iyc, izc + 1, ispin_comp, istate) - phi_box(ixc, iyc, izc - 1, ispin_comp, istate)) + &
                  nabt(2,3) * (phi_box(ixc, iyc, izc + 2, ispin_comp, istate) - phi_box(ixc, iyc, izc - 2, ispin_comp, istate)) + &
                  nabt(3,3) * (phi_box(ixc, iyc, izc + 3, ispin_comp, istate) - phi_box(ixc, iyc, izc - 3, ispin_comp, istate)) + &
                  nabt(4,3) * (phi_box(ixc, iyc, izc + 4, ispin_comp, istate) - phi_box(ixc, iyc, izc - 4, ispin_comp, istate))
    case default
      return
    end select
    dn = real(side, 8) * grad_axis
  end function phi_box_dn_ham_soi

  complex(8) function phi_local_dn_ham_soi(dg_frag, i_local, stencil, ispin_comp, istate, gidx, axis, side) result(dn)
    use structures, only: s_stencil
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, ispin_comp, istate, gidx(3), axis, side
    type(s_stencil), intent(in) :: stencil
    integer :: ix, iy, iz, ixc, iyc, izc
    real(8) :: nabt(4,3)
    complex(8) :: grad_axis

    dn = (0.0d0, 0.0d0)
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) return
    if (istate < 1 .or. istate > size(dg_frag%phi_frag, 4)) return
    ix = map_global_to_phi_box_coord_ham_soi(gidx(1), lbound(dg_frag%phi_frag, 1), &
                                            ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
    iy = map_global_to_phi_box_coord_ham_soi(gidx(2), lbound(dg_frag%phi_frag, 2), &
                                            ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
    iz = map_global_to_phi_box_coord_ham_soi(gidx(3), lbound(dg_frag%phi_frag, 3), &
                                            ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    ixc = ix
    iyc = iy
    izc = iz
    nabt = stencil%coef_nab
    select case(axis)
    case(1)
      if (ixc - 4 < lbound(dg_frag%phi_frag, 1) .or. ixc + 4 > ubound(dg_frag%phi_frag, 1)) return
      grad_axis = nabt(1,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc + 1, iyc, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc - 1, iyc, izc])) + &
                  nabt(2,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc + 2, iyc, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc - 2, iyc, izc])) + &
                  nabt(3,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc + 3, iyc, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc - 3, iyc, izc])) + &
                  nabt(4,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc + 4, iyc, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc - 4, iyc, izc]))
    case(2)
      if (iyc - 4 < lbound(dg_frag%phi_frag, 2) .or. iyc + 4 > ubound(dg_frag%phi_frag, 2)) return
      grad_axis = nabt(1,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc + 1, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc - 1, izc])) + &
                  nabt(2,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc + 2, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc - 2, izc])) + &
                  nabt(3,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc + 3, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc - 3, izc])) + &
                  nabt(4,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc + 4, izc]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc - 4, izc]))
    case(3)
      if (izc - 4 < lbound(dg_frag%phi_frag, 3) .or. izc + 4 > ubound(dg_frag%phi_frag, 3)) return
      grad_axis = nabt(1,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc + 1]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc - 1])) + &
                  nabt(2,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc + 2]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc - 2])) + &
                  nabt(3,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc + 3]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc - 3])) + &
                  nabt(4,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc + 4]) - &
                               phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [ixc, iyc, izc - 4]))
    case default
      return
    end select
    dn = real(side, 8) * grad_axis
  end function phi_local_dn_ham_soi

  subroutine phi_local_grad_ham_soi(dg_frag, i_local, stencil, ispin_comp, istate, gidx, grad)
    use structures, only: s_stencil
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, ispin_comp, istate, gidx(3)
    type(s_stencil), intent(in) :: stencil
    complex(8), intent(out) :: grad(3)
    real(8) :: nabt(4,3)

    nabt = stencil%coef_nab
    grad(:) = (0.0d0, 0.0d0)
    grad(1) = nabt(1,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) + 1, gidx(2), gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) - 1, gidx(2), gidx(3)])) + &
              nabt(2,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) + 2, gidx(2), gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) - 2, gidx(2), gidx(3)])) + &
              nabt(3,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) + 3, gidx(2), gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) - 3, gidx(2), gidx(3)])) + &
              nabt(4,1) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) + 4, gidx(2), gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1) - 4, gidx(2), gidx(3)]))
    grad(2) = nabt(1,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) + 1, gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) - 1, gidx(3)])) + &
              nabt(2,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) + 2, gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) - 2, gidx(3)])) + &
              nabt(3,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) + 3, gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) - 3, gidx(3)])) + &
              nabt(4,2) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) + 4, gidx(3)]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2) - 4, gidx(3)]))
    grad(3) = nabt(1,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) + 1]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) - 1])) + &
              nabt(2,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) + 2]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) - 2])) + &
              nabt(3,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) + 3]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) - 3])) + &
              nabt(4,3) * (phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) + 4]) - &
                            phi_local_value_ham_soi(dg_frag, i_local, ispin_comp, istate, [gidx(1), gidx(2), gidx(3) - 4]))
  end subroutine phi_local_grad_ham_soi

  logical function basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, ibasis) result(is_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ispin, ibasis
    integer :: row_idx

    is_local = .false.
    if (ibasis < 1 .or. ibasis > size(dg_frag%index_basis, 1)) return
    if (ifrag < 1 .or. ifrag > size(dg_frag%index_basis, 2)) return
    if (ispin < 1 .or. ispin > size(dg_frag%index_basis, 3)) return
    row_idx = dg_frag%index_basis(ibasis, ifrag, ispin)
    is_local = matrix_row_is_locally_owned(dg_frag, row_idx, ispin)
  end function basis_row_is_locally_owned_soi

  subroutine add_dg_surface_flux_blocks_soi(dg_frag, system, mg, stencil, H_blocks, T_blocks, block_map)
    use structures, only: s_dft_system, s_rgrid, s_stencil
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    type(matrix_block_info), intent(inout) :: H_blocks(:)
    type(matrix_block_info), intent(inout) :: T_blocks(:)
    integer, intent(in) :: block_map(:, :)

    real(8), parameter :: surface_penalty_factor = 10.0d0
    integer :: ifrag, jfrag, i_local, axis, side, ispin, io, jo
    integer :: iblk_self, iblk_cross, nrow, ncol, ix, iy, iz
    integer :: g_row(3), g_col(3), loop_lo(3), loop_hi(3)
    integer :: col_lo(3), col_hi(3)
    real(8) :: area_weight, alpha
    complex(8) :: u_l, u_r, v_l, dnu_l, dnu_r, dnv_l
    complex(8) :: term_self, term_cross
    complex(8), allocatable :: phi_col(:,:,:,:,:)

    if (.not. allocated(dg_frag%phi_frag)) return
    if (.not. allocated(dg_frag%phi_frag_spinor_c)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. dg_frag%parallel_mode_orbital) return

    ! SOI uses the same local DG volume operator as non-SOI.  Only the
    ! interface term from Eq. (4) is stored in off-diagonal face-neighbor
    ! blocks; the spinor basis enters through complex conjugated bra values.
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) cycle

      do axis = 1, 3
        if (dg_frag%num_fragment(axis) <= 1) cycle
        if (system%hgs(axis) <= 0.0d0) cycle
        area_weight = system%hvol / system%hgs(axis)
        alpha = surface_penalty_factor / system%hgs(axis)
        do side = -1, 1, 2
          jfrag = face_neighbor_fragment_ham_soi(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          iblk_self = find_matrix_block(block_map, ifrag, ifrag)
          iblk_cross = find_matrix_block(block_map, ifrag, jfrag)
          if (iblk_self <= 0 .and. iblk_cross <= 0) cycle

          call read_fragment_buffer_basis_box_ham_soi(dg_frag, jfrag, phi_col, col_lo, col_hi)

          loop_lo(:) = dg_frag%ixyz_frag(:, ifrag)
          loop_hi(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          if (side > 0) then
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
            loop_hi(axis) = loop_lo(axis)
          else
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag)
            loop_hi(axis) = loop_lo(axis)
          end if

          do ispin = 1, dg_frag%nspin
            nrow = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
            ncol = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), size(phi_col, 5))
            if (nrow <= 0) cycle

            do iz = loop_lo(3), loop_hi(3)
              do iy = loop_lo(2), loop_hi(2)
                do ix = loop_lo(1), loop_hi(1)
                  g_row = [ix, iy, iz]
                  g_col = g_row
                  if (side > 0) then
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag)
                  else
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag) + dg_frag%nxyz_domain(axis, jfrag) - 1
                  end if

                  if (iblk_self > 0) then
                    do jo = 1, nrow
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, io)) cycle
                        u_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, jo, g_row)
                        dnu_l = phi_local_dn_ham_soi(dg_frag, i_local, stencil, ispin, jo, g_row, axis, side)
                        v_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, io, g_row)
                        dnv_l = phi_local_dn_ham_soi(dg_frag, i_local, stencil, ispin, io, g_row, axis, side)
                        term_self = (-0.25d0 * conjg(v_l) * dnu_l - &
                                     0.25d0 * conjg(dnv_l) * u_l + &
                                     alpha * conjg(v_l) * u_l) * area_weight
                        H_blocks(iblk_self)%val(io, jo, ispin) = &
                          H_blocks(iblk_self)%val(io, jo, ispin) + real(term_self, kind=8)
                        T_blocks(iblk_self)%val(io, jo, ispin) = &
                          T_blocks(iblk_self)%val(io, jo, ispin) + real(term_self, kind=8)
                      end do
                    end do
                  end if

                  if (iblk_cross > 0 .and. ncol > 0) then
                    do jo = 1, ncol
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, io)) cycle
                        u_r = phi_box_value_ham_soi(phi_col, col_lo, col_hi, dg_frag%lgnum_total, &
                                                   ispin, jo, g_col)
                        dnu_r = phi_box_dn_ham_soi(phi_col, col_lo, col_hi, dg_frag%lgnum_total, stencil, &
                                                  ispin, jo, g_col, axis, side)
                        v_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, io, g_row)
                        dnv_l = phi_local_dn_ham_soi(dg_frag, i_local, stencil, ispin, io, g_row, axis, side)
                        term_cross = (-0.25d0 * conjg(v_l) * dnu_r + &
                                      0.25d0 * conjg(dnv_l) * u_r - &
                                      alpha * conjg(v_l) * u_r) * area_weight
                        H_blocks(iblk_cross)%val(io, jo, ispin) = &
                          H_blocks(iblk_cross)%val(io, jo, ispin) + real(term_cross, kind=8)
                        T_blocks(iblk_cross)%val(io, jo, ispin) = &
                          T_blocks(iblk_cross)%val(io, jo, ispin) + real(term_cross, kind=8)
                      end do
                    end do
                  end if
                end do
              end do
            end do
          end do

          if (allocated(phi_col)) deallocate(phi_col)
        end do
      end do
    end do
  end subroutine add_dg_surface_flux_blocks_soi

  subroutine add_dg_surface_momentum_blocks_soi(dg_frag, system)
    use structures, only: s_dft_system
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system

    integer :: ifrag, jfrag, i_local, axis, side, ispin, io, jo
    integer :: iblk_self, iblk_cross, nrow, ncol, ix, iy, iz
    integer :: g_row(3), g_col(3), loop_lo(3), loop_hi(3)
    integer :: col_lo(3), col_hi(3)
    real(8) :: area_weight, normal_sign
    complex(8) :: u_l, u_r, v_l, term_self, term_cross
    complex(8), allocatable :: phi_col(:,:,:,:,:)

    if (.not. allocated(dg_frag%momentum_blocks)) return
    if (.not. allocated(dg_frag%momentum_block_map)) return
    if (.not. allocated(dg_frag%phi_frag)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. dg_frag%parallel_mode_orbital) return

    ! Current/momentum keeps the same Galerkin topology as H: volume
    ! gradients are local to a fragment, and only the normal face correction
    ! connects neighboring fragments.
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) cycle

      do axis = 1, 3
        if (dg_frag%num_fragment(axis) <= 1) cycle
        if (system%hgs(axis) <= 0.0d0) cycle
        area_weight = system%hvol / system%hgs(axis)
        do side = -1, 1, 2
          normal_sign = real(side, kind=8)
          jfrag = face_neighbor_fragment_ham_soi(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          iblk_self = find_momentum_block(dg_frag, ifrag, ifrag)
          iblk_cross = find_momentum_block(dg_frag, ifrag, jfrag)
          if (iblk_self <= 0 .and. iblk_cross <= 0) cycle

          call read_fragment_buffer_basis_box_ham_soi(dg_frag, jfrag, phi_col, col_lo, col_hi)

          loop_lo(:) = dg_frag%ixyz_frag(:, ifrag)
          loop_hi(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          if (side > 0) then
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
            loop_hi(axis) = loop_lo(axis)
          else
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag)
            loop_hi(axis) = loop_lo(axis)
          end if

          do ispin = 1, dg_frag%nspin
            nrow = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
            ncol = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), size(phi_col, 5))
            if (nrow <= 0) cycle

            do iz = loop_lo(3), loop_hi(3)
              do iy = loop_lo(2), loop_hi(2)
                do ix = loop_lo(1), loop_hi(1)
                  g_row = [ix, iy, iz]
                  g_col = g_row
                  if (side > 0) then
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag)
                  else
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag) + dg_frag%nxyz_domain(axis, jfrag) - 1
                  end if

                  if (iblk_self > 0) then
                    do jo = 1, nrow
                      u_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, jo, g_row)
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, io)) cycle
                        v_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, io, g_row)
                        term_self = -0.5d0 * normal_sign * conjg(v_l) * u_l * area_weight
                        dg_frag%momentum_blocks_c(iblk_self)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks_c(iblk_self)%val(axis, io, jo, ispin) + term_self
                        dg_frag%momentum_blocks(iblk_self)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks(iblk_self)%val(axis, io, jo, ispin) + real(term_self, kind=8)
                      end do
                    end do
                  end if

                  if (iblk_cross > 0 .and. ncol > 0) then
                    do jo = 1, ncol
                      u_r = phi_box_value_ham_soi(phi_col, col_lo, col_hi, dg_frag%lgnum_total, &
                                                 ispin, jo, g_col)
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned_soi(dg_frag, ifrag, ispin, io)) cycle
                        v_l = phi_local_value_ham_soi(dg_frag, i_local, ispin, io, g_row)
                        term_cross = 0.5d0 * normal_sign * conjg(v_l) * u_r * area_weight
                        dg_frag%momentum_blocks_c(iblk_cross)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks_c(iblk_cross)%val(axis, io, jo, ispin) + term_cross
                        dg_frag%momentum_blocks(iblk_cross)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks(iblk_cross)%val(axis, io, jo, ispin) + real(term_cross, kind=8)
                      end do
                    end do
                  end if
                end do
              end do
            end do
          end do

          if (allocated(phi_col)) deallocate(phi_col)
        end do
      end do
    end do
  end subroutine add_dg_surface_momentum_blocks_soi

  subroutine init_complex_momentum_blocks(dg_frag, blocks_re, blocks_im, block_map, n_blocks)
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(momentum_block_info), allocatable, intent(inout) :: blocks_re(:), blocks_im(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max

    if (allocated(blocks_re)) then
      do iblk = 1, size(blocks_re)
        if (allocated(blocks_re(iblk)%val)) deallocate(blocks_re(iblk)%val)
      end do
      deallocate(blocks_re)
    end if
    if (allocated(blocks_im)) then
      do iblk = 1, size(blocks_im)
        if (allocated(blocks_im(iblk)%val)) deallocate(blocks_im(iblk)%val)
      end do
      deallocate(blocks_im)
    end if
    if (allocated(block_map)) deallocate(block_map)

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) nblk = nblk + 1
      end do
    end do

    n_blocks = nblk
    if (nblk <= 0) return
    allocate(blocks_re(nblk), blocks_im(nblk))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        block_map(ifrag_row, ifrag_col) = iblk
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))

        blocks_re(iblk)%ifrag_row = ifrag_row
        blocks_re(iblk)%ifrag_col = ifrag_col
        blocks_re(iblk)%nrow_max = nrow_max
        blocks_re(iblk)%ncol_max = ncol_max
        allocate(blocks_re(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        blocks_re(iblk)%val = 0.0d0

        blocks_im(iblk)%ifrag_row = ifrag_row
        blocks_im(iblk)%ifrag_col = ifrag_col
        blocks_im(iblk)%nrow_max = nrow_max
        blocks_im(iblk)%ncol_max = ncol_max
        allocate(blocks_im(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        blocks_im(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_complex_momentum_blocks

  subroutine reduce_complex_momentum_blocks(dg_frag, blocks_re, blocks_im, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(momentum_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, idir, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat
    logical, parameter :: trace_block_reduce = .false.

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = 3 * nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = 3 * nrow * ncol

        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              send_block(offset_flat) = blocks_re(iblk)%val(idir, ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              blocks_re(iblk)%val(idir, ii, jj, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          end do
        end do

        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              send_block(offset_flat) = blocks_im(iblk)%val(idir, ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              blocks_im(iblk)%val(idir, ii, jj, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_complex_momentum_blocks

  subroutine reduce_matrix_blocks(dg_frag, blocks, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat
    logical, parameter :: trace_block_reduce = .false.

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do

        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_blocks

  subroutine init_complex_matrix_blocks(dg_frag, blocks_re, blocks_im, block_map, n_blocks)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks_re(:), blocks_im(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max

    if (allocated(blocks_re)) then
      do iblk = 1, size(blocks_re)
        if (allocated(blocks_re(iblk)%val)) deallocate(blocks_re(iblk)%val)
      end do
      deallocate(blocks_re)
    end if
    if (allocated(blocks_im)) then
      do iblk = 1, size(blocks_im)
        if (allocated(blocks_im(iblk)%val)) deallocate(blocks_im(iblk)%val)
      end do
      deallocate(blocks_im)
    end if
    if (allocated(block_map)) deallocate(block_map)

    n_blocks = dg_frag%n_frag * dg_frag%n_frag
    if (n_blocks <= 0) return

    allocate(blocks_re(n_blocks), blocks_im(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk

        blocks_re(iblk)%ifrag_row = ifrag_row
        blocks_re(iblk)%ifrag_col = ifrag_col
        blocks_re(iblk)%nrow_max = nrow_max
        blocks_re(iblk)%ncol_max = ncol_max
        allocate(blocks_re(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks_re(iblk)%val = 0.0d0

        blocks_im(iblk)%ifrag_row = ifrag_row
        blocks_im(iblk)%ifrag_col = ifrag_col
        blocks_im(iblk)%nrow_max = nrow_max
        blocks_im(iblk)%ncol_max = ncol_max
        allocate(blocks_im(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks_im(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_complex_matrix_blocks



  subroutine reduce_complex_matrix_blocks(dg_frag, blocks_re, blocks_im, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat
    logical, parameter :: trace_block_reduce = .false.

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks_re(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks_re(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks_im(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks_im(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (trace_block_reduce .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_complex_matrix_blocks

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
  subroutine calculate_overlap_matrix(dg_frag, system, mg)
    use structures
    use rt_dg_fragment_types, only: matrix_block_info
    use communication, only: comm_is_root, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg

    integer :: ifrag, i_local, ispin, io, jo, nbf_local, iblk, iblk_t
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo
    integer :: ig_row, ig_col, l(3), d(3), ii, jj, halo_send_idx(3), halo_recv_idx(3)
    integer :: lx, ly, lz, iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: ixg, iyg, izg, bx, by, bz
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: n_eval, lwork, info_eig, n_blocks, icomm_reduce
    logical :: use_complex_halo, use_complex_phi, has_overlap
    real(8) :: hvol, s_min, s_max, cond_est
    complex(8) :: integral, savg
    complex(8) :: cwork_query(1)
    complex(8), allocatable :: S_eval(:,:), eig_work(:)
    real(8), allocatable :: eigvals(:), rwork(:)
    type(matrix_block_info), allocatable :: S_blocks_re(:), S_blocks_im(:)
    integer, allocatable :: S_block_map_local(:,:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    if (allocated(dg_frag%S_mat_c)) deallocate(dg_frag%S_mat_c)
    if (allocated(dg_frag%S_mat_prop_c)) deallocate(dg_frag%S_mat_prop_c)
    call init_matrix_blocks(dg_frag, S_blocks_re, S_block_map_local, n_blocks, diagonal_only=.true.)
    if (allocated(S_blocks_im)) then
      do iblk = 1, size(S_blocks_im)
        if (allocated(S_blocks_im(iblk)%val)) deallocate(S_blocks_im(iblk)%val)
      end do
      deallocate(S_blocks_im)
    end if
    if (allocated(S_blocks_re)) then
      allocate(S_blocks_im(size(S_blocks_re)))
      do iblk = 1, size(S_blocks_re)
        S_blocks_im(iblk)%ifrag_row = S_blocks_re(iblk)%ifrag_row
        S_blocks_im(iblk)%ifrag_col = S_blocks_re(iblk)%ifrag_col
        S_blocks_im(iblk)%nrow_max = S_blocks_re(iblk)%nrow_max
        S_blocks_im(iblk)%ncol_max = S_blocks_re(iblk)%ncol_max
        allocate(S_blocks_im(iblk)%val(S_blocks_re(iblk)%nrow_max, S_blocks_re(iblk)%ncol_max, dg_frag%nspin))
        S_blocks_im(iblk)%val = 0.0d0
      end do
    end if

    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)

    use_complex_phi = allocated(dg_frag%phi_frag_c)

    do ispin = 1, system%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        nbf_local = dg_frag%n_basis(ifrag, ispin)
        if (nbf_local <= 0) cycle
        iblk = find_matrix_block(S_block_map_local, ifrag, ifrag)
        if (iblk <= 0) cycle
        ! S uses the same box ownership as H and momentum.  Each rank
        ! integrates its local real-space box, then the block reduction below
        ! reconstructs the full fragment overlap matrix.
        call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
        if (.not. has_overlap) cycle

        do jo = 1, nbf_local
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle

          do io = 1, nbf_local
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = (0.0d0, 0.0d0)
            do lz = loc_s(3), loc_e(3)
              izg = iorg(3) + lz - 1
              bz = map_global_to_phi_box_coord_ham_soi(izg, phi_lb3, phi_ub3, dg_frag%lgnum_total(3))
              if (bz == 0) cycle
              do ly = loc_s(2), loc_e(2)
                iyg = iorg(2) + ly - 1
                by = map_global_to_phi_box_coord_ham_soi(iyg, phi_lb2, phi_ub2, dg_frag%lgnum_total(2))
                if (by == 0) cycle
                do lx = loc_s(1), loc_e(1)
                  ixg = iorg(1) + lx - 1
                  bx = map_global_to_phi_box_coord_ham_soi(ixg, phi_lb1, phi_ub1, dg_frag%lgnum_total(1))
                  if (bx == 0) cycle
                  if (allocated(dg_frag%phi_frag_spinor_c) .and. ispin >= 1 .and. &
                      ispin <= size(dg_frag%phi_frag_spinor_c, 4)) then
                    integral = integral + conjg(dg_frag%phi_frag_spinor_c(bx, by, bz, ispin, io, i_local)) * &
                               dg_frag%phi_frag_spinor_c(bx, by, bz, ispin, jo, i_local) * hvol
                  else if (use_complex_phi) then
                    integral = integral + conjg(dg_frag%phi_frag_c(bx, by, bz, io, i_local)) * &
                               dg_frag%phi_frag_c(bx, by, bz, jo, i_local) * hvol
                  else
                    integral = integral + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8) * &
                               cmplx(dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8) * hvol
                  end if
                end do
              end do
            end do
            S_blocks_re(iblk)%val(io, jo, ispin) = real(integral, kind=8)
            S_blocks_im(iblk)%val(io, jo, ispin) = aimag(integral)
          end do

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            iblk = find_matrix_block(S_block_map_local, jfrag, ifrag)
            iblk_t = find_matrix_block(S_block_map_local, ifrag, jfrag)
            if (iblk <= 0 .or. iblk_t <= 0) cycle
            l = dg_frag%halo(i_halo)%length
            use_complex_halo = allocated(dg_frag%halo(i_halo)%buf_recv_c) .and. use_complex_phi

            do io = 1, n_basis_halo
              ig_row = dg_frag%index_basis(io, jfrag, ispin)
              if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
              integral = (0.0d0, 0.0d0)
              if (use_complex_halo) then
                do iz = 1, l(3)
                  do iy = 1, l(2)
                    do ix = 1, l(1)
                      call get_halo_block_point_indices(dg_frag%halo(i_halo), ix, iy, iz, halo_send_idx, halo_recv_idx)
                      integral = integral + conjg(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, io, 1)) * &
                                 dg_frag%phi_frag_c(halo_recv_idx(1), halo_recv_idx(2), halo_recv_idx(3), jo, i_local) * hvol
                    end do
                  end do
                end do
              else
                do iz = 1, l(3)
                  do iy = 1, l(2)
                    do ix = 1, l(1)
                      call get_halo_block_point_indices(dg_frag%halo(i_halo), ix, iy, iz, halo_send_idx, halo_recv_idx)
                      integral = integral + &
                        cmplx(dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1), 0.0d0, kind=8) * &
                        cmplx(dg_frag%phi_frag(halo_recv_idx(1), halo_recv_idx(2), halo_recv_idx(3), &
                                               jo, i_local), 0.0d0, kind=8) * hvol
                    end do
                  end do
                end do
              end if
              S_blocks_re(iblk)%val(io, jo, ispin) = &
                S_blocks_re(iblk)%val(io, jo, ispin) + 0.5d0 * real(integral, kind=8)
              S_blocks_im(iblk)%val(io, jo, ispin) = &
                S_blocks_im(iblk)%val(io, jo, ispin) + 0.5d0 * aimag(integral)
              S_blocks_re(iblk_t)%val(jo, io, ispin) = &
                S_blocks_re(iblk_t)%val(jo, io, ispin) + 0.5d0 * real(integral, kind=8)
              S_blocks_im(iblk_t)%val(jo, io, ispin) = &
                S_blocks_im(iblk_t)%val(jo, io, ispin) - 0.5d0 * aimag(integral)
            end do
          end do

        end do
      end do
    end do

    icomm_reduce = dg_frag%icomm
    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) icomm_reduce = dg_frag%icomm_frag

    ! Reduce complex S through compact blocks.  Non-SOI and SOI must share this
    ! ownership rule: box contributions are summed, while orbital ranks retain
    ! the replicated reduced matrix for propagation.
    call reduce_complex_matrix_blocks(dg_frag, S_blocks_re, S_blocks_im, "smat-soi", icomm_reduce)
    if (icomm_reduce == dg_frag%icomm_frag .and. .not. dg_frag%is_frag_root &
        .and. .not. dg_frag%parallel_mode_orbital) then
      ! Non-root fragment ranks keep no authoritative dense overlap copy.
      ! Block application paths below remain valid through S_mat_blocks.
    end if

    do iblk = 1, n_blocks
      ifrag = S_blocks_re(iblk)%ifrag_row
      jfrag = S_blocks_re(iblk)%ifrag_col
      iblk_t = find_matrix_block(S_block_map_local, jfrag, ifrag)
      if (iblk_t <= 0) cycle
      do ispin = 1, dg_frag%nspin
        nbf_local = dg_frag%n_basis(ifrag, ispin)
        n_basis_halo = dg_frag%n_basis(jfrag, ispin)
        if (nbf_local <= 0 .or. n_basis_halo <= 0) cycle
        do jj = 1, n_basis_halo
          do ii = 1, nbf_local
            if (ifrag == jfrag .and. ii == jj) then
              if (S_blocks_re(iblk)%val(ii, jj, ispin) < 1.0d-12) S_blocks_re(iblk)%val(ii, jj, ispin) = 1.0d0
              S_blocks_im(iblk)%val(ii, jj, ispin) = 0.0d0
            else
              savg = 0.5d0 * (cmplx(S_blocks_re(iblk)%val(ii, jj, ispin), &
                                    S_blocks_im(iblk)%val(ii, jj, ispin), kind=8) + &
                              conjg(cmplx(S_blocks_re(iblk_t)%val(jj, ii, ispin), &
                                          S_blocks_im(iblk_t)%val(jj, ii, ispin), kind=8)))
              S_blocks_re(iblk)%val(ii, jj, ispin) = real(savg, kind=8)
              S_blocks_im(iblk)%val(ii, jj, ispin) = aimag(savg)
              S_blocks_re(iblk_t)%val(jj, ii, ispin) = real(savg, kind=8)
              S_blocks_im(iblk_t)%val(jj, ii, ispin) = -aimag(savg)
            end if
          end do
        end do
      end do
    end do

    call init_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
                            diagonal_only=.true.)
    call init_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
                            diagonal_only=.true.)
    call init_complex_overlap_blocks_soi(dg_frag, dg_frag%S_mat_blocks_c)
    call init_complex_overlap_blocks_soi(dg_frag, dg_frag%S_mat_prop_blocks_c)
    do iblk = 1, min(size(S_blocks_re), size(dg_frag%S_mat_blocks), size(dg_frag%S_mat_prop_blocks))
      dg_frag%S_mat_blocks(iblk)%val(:, :, :) = S_blocks_re(iblk)%val(:, :, :)
      dg_frag%S_mat_prop_blocks(iblk)%val(:, :, :) = S_blocks_re(iblk)%val(:, :, :)
      if (allocated(dg_frag%S_mat_blocks_c) .and. iblk <= size(dg_frag%S_mat_blocks_c)) then
        dg_frag%S_mat_blocks_c(iblk)%val(:, :, :) = cmplx(S_blocks_re(iblk)%val(:, :, :), &
                                                          S_blocks_im(iblk)%val(:, :, :), kind=8)
      end if
      if (allocated(dg_frag%S_mat_prop_blocks_c) .and. iblk <= size(dg_frag%S_mat_prop_blocks_c)) then
        dg_frag%S_mat_prop_blocks_c(iblk)%val(:, :, :) = cmplx(S_blocks_re(iblk)%val(:, :, :), &
                                                               S_blocks_im(iblk)%val(:, :, :), kind=8)
      end if
    end do
    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .false.

    if (allocated(S_blocks_re)) then
      do ii = 1, size(S_blocks_re)
        if (allocated(S_blocks_re(ii)%val)) deallocate(S_blocks_re(ii)%val)
      end do
      deallocate(S_blocks_re)
    end if
    if (allocated(S_blocks_im)) then
      do ii = 1, size(S_blocks_im)
        if (allocated(S_blocks_im(ii)%val)) deallocate(S_blocks_im(ii)%val)
      end do
      deallocate(S_blocks_im)
    end if
    if (allocated(S_block_map_local)) deallocate(S_block_map_local)

  end subroutine calculate_overlap_matrix

  subroutine init_complex_overlap_blocks_soi(dg_frag, blocks)
    use rt_dg_fragment_types, only: complex_matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), allocatable, intent(inout) :: blocks(:)

    integer :: iblk, nrow_max, ncol_max

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (.not. allocated(dg_frag%S_mat_blocks)) return
    allocate(blocks(size(dg_frag%S_mat_blocks)))
    do iblk = 1, size(dg_frag%S_mat_blocks)
      blocks(iblk)%ifrag_row = dg_frag%S_mat_blocks(iblk)%ifrag_row
      blocks(iblk)%ifrag_col = dg_frag%S_mat_blocks(iblk)%ifrag_col
      blocks(iblk)%nrow_max = dg_frag%S_mat_blocks(iblk)%nrow_max
      blocks(iblk)%ncol_max = dg_frag%S_mat_blocks(iblk)%ncol_max
      nrow_max = max(1, blocks(iblk)%nrow_max)
      ncol_max = max(1, blocks(iblk)%ncol_max)
      allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
      blocks(iblk)%val(:, :, :) = (0.0d0, 0.0d0)
    end do
  end subroutine init_complex_overlap_blocks_soi

end module rt_dg_fragment_hamiltonian_soi
