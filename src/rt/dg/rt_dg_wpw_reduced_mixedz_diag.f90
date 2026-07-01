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
!=======================================================================
! WPW reduced canonical mixed-Z diagnostic helpers
!=======================================================================

module rt_dg_wpw_reduced_mixedz_diag
  use communication, only: comm_summation
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_wpw_linalg, only: build_hermitian_pseudoinverse, wpw_local_herm_max
  use structures, only: s_dft_system
  implicit none

  private
  public :: wpw_reduced_mixedz_operator_stats
  public :: wpw_reduced_canon_mixedz_current_coeff_stats
  public :: wpw_reduced_canon_mixedz_bridge_stats
  public :: wpw_reduced_canon_pz_block_operator_stats
  public :: wpw_reduced_canon_p_projection_stats
  public :: wpw_reduced_prod_p_basis_save_stats

contains

  subroutine wpw_reduced_mixedz_operator_stats(dg_frag, herm_diff, trace_z, basis_dim, neighbor_terms_included)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(out) :: herm_diff, trace_z
    integer, intent(out) :: basis_dim
    logical, intent(out) :: neighbor_terms_included

    integer :: i, j, ispin, nmix
    complex(8) :: ztrace_c

    herm_diff = huge(1.0d0)
    trace_z = 0.0d0
    basis_dim = 0
    neighbor_terms_included = .false.
    if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
    if (dg_frag%mixed_wannier_bpw_nmix <= 0) return

    nmix = min(dg_frag%mixed_wannier_bpw_nmix, size(dg_frag%mixed_wannier_bpw_z, 2), &
      size(dg_frag%mixed_wannier_bpw_z, 3))
    if (nmix <= 0) return
    basis_dim = nmix
    neighbor_terms_included = dg_frag%mixed_wannier_bpw_np > 0
    herm_diff = 0.0d0
    ztrace_c = (0.0d0, 0.0d0)
    do ispin = 1, size(dg_frag%mixed_wannier_bpw_z, 4)
      do i = 1, nmix
        ztrace_c = ztrace_c + dg_frag%mixed_wannier_bpw_z(3, i, i, ispin)
        do j = 1, nmix
          herm_diff = max(herm_diff, abs(dg_frag%mixed_wannier_bpw_z(3, i, j, ispin) - &
            conjg(dg_frag%mixed_wannier_bpw_z(3, j, i, ispin))))
        end do
      end do
    end do
    trace_z = real(ztrace_c, kind=8)
  end subroutine wpw_reduced_mixedz_operator_stats


  subroutine wpw_reduced_canon_mixedz_current_coeff_stats(dg_frag, system, pz_prod_mixedz, &
      pz_can_mixedz, pz_diff, rel_pz_diff, zherm_diff, canonical_operator_built, convention_match, bad)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    real(8), intent(in) :: pz_prod_mixedz
    real(8), intent(out) :: pz_can_mixedz, pz_diff, rel_pz_diff, zherm_diff
    logical, intent(out) :: canonical_operator_built, convention_match, bad

    integer :: ispin, ifrag, i_local, ib, iw, jw, istate
    integer :: n_w, n_p, n_mix, nocc_spin, nbasis, global_row, local_row
    real(8) :: occ, pz_local, pz_sum, tol
    complex(8) :: pos_expect
    complex(8), allocatable :: cw_local(:,:), cw_sum(:,:), cmix(:,:)

    pz_can_mixedz = 0.0d0
    pz_diff = -pz_prod_mixedz
    rel_pz_diff = pz_diff / max(abs(pz_prod_mixedz), 1.0d-300)
    zherm_diff = huge(1.0d0)
    canonical_operator_built = .false.
    convention_match = .false.
    bad = .true.

    if (.not. dg_frag%has_mixed_wannier_bpw_position) return
    if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
    if (.not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) return
    if (.not. allocated(dg_frag%global_wannier_coef)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    n_w = dg_frag%mixed_wannier_bpw_nw
    n_p = dg_frag%mixed_wannier_bpw_np
    n_mix = dg_frag%mixed_wannier_bpw_nmix
    if (n_w <= 0 .or. n_p < 0 .or. n_mix /= n_w + n_p) return
    if (size(dg_frag%mixed_wannier_bpw_z, 1) < 3) return
    if (size(dg_frag%mixed_wannier_bpw_z, 2) < n_mix .or. &
        size(dg_frag%mixed_wannier_bpw_z, 3) < n_mix) return

    allocate(cw_local(n_w,max(1,dg_frag%nstate_tot)))
    allocate(cw_sum(n_w,max(1,dg_frag%nstate_tot)))
    allocate(cmix(n_mix,max(1,dg_frag%nstate_tot)))
    pz_local = 0.0d0
    zherm_diff = 0.0d0

    do ispin = 1, min(dg_frag%nspin, system%nspin, size(dg_frag%mixed_wannier_bpw_z, 4))
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      if (nocc_spin <= 0) cycle

      cw_local(:, :) = (0.0d0, 0.0d0)
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1))
        do ib = 1, nbasis
          global_row = dg_frag%index_basis(ib, ifrag, ispin)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          do iw = 1, n_w
            do istate = 1, nocc_spin
              cw_local(iw,istate) = cw_local(iw,istate) + &
                conjg(dg_frag%global_wannier_coef(ib, iw, ispin, i_local)) * &
                dg_frag%coef(local_row,istate,ispin)
            end do
          end do
        end do
      end do
      call comm_summation(cw_local, cw_sum, n_w * max(1,dg_frag%nstate_tot), dg_frag%icomm)

      do iw = 1, n_mix
        do jw = 1, n_mix
          zherm_diff = max(zherm_diff, abs(dg_frag%mixed_wannier_bpw_z(3, iw, jw, ispin) - &
            conjg(dg_frag%mixed_wannier_bpw_z(3, jw, iw, ispin))))
        end do
      end do

      if (dg_frag%id == 0) then
        cmix(:, :) = (0.0d0, 0.0d0)
        cmix(1:n_w,1:nocc_spin) = cw_sum(1:n_w,1:nocc_spin)
        if (n_p > 0) cmix(n_w+1:n_w+n_p,1:nocc_spin) = &
          dg_frag%mixed_wannier_bpw_pcoef(1:n_p,1:nocc_spin,ispin)
        do istate = 1, nocc_spin
          occ = max(0.0d0, system%rocc(istate, 1, ispin))
          if (occ <= 0.0d0) cycle
          pos_expect = (0.0d0, 0.0d0)
          do iw = 1, n_mix
            do jw = 1, n_mix
              pos_expect = pos_expect + conjg(cmix(iw,istate)) * &
                dg_frag%mixed_wannier_bpw_z(3,iw,jw,ispin) * cmix(jw,istate)
            end do
          end do
          pz_local = pz_local - occ * real(pos_expect, kind=8)
        end do
      end if
    end do

    call comm_summation(pz_local, pz_sum, 1, dg_frag%icomm)
    pz_can_mixedz = pz_sum
    pz_diff = pz_can_mixedz - pz_prod_mixedz
    rel_pz_diff = pz_diff / max(abs(pz_prod_mixedz), 1.0d-300)
    tol = 1.0d-10 * max(1.0d0, abs(pz_prod_mixedz))
    canonical_operator_built = .true.
    convention_match = abs(pz_diff) <= tol
    bad = (pz_can_mixedz /= pz_can_mixedz) .or. (zherm_diff /= zherm_diff) .or. &
      (.not. convention_match)

    deallocate(cw_local, cw_sum, cmix)
  end subroutine wpw_reduced_canon_mixedz_current_coeff_stats


  subroutine wpw_reduced_canon_mixedz_bridge_stats(dg_frag, pz_prod_mixedz, t_dim_can, t_dim_mixed, &
      t_norm, t_rank_est, bridge_roundtrip_diff, zcan_herm_diff, pz_can_mixedz, pz_diff, rel_pz_diff, &
      canonical_operator_built, convention_match, bad)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: pz_prod_mixedz
    integer, intent(out) :: t_dim_can, t_dim_mixed, t_rank_est
    real(8), intent(out) :: t_norm, bridge_roundtrip_diff, zcan_herm_diff
    real(8), intent(out) :: pz_can_mixedz, pz_diff, rel_pz_diff
    logical, intent(out) :: canonical_operator_built, convention_match, bad

    integer :: i, nmix, ncan, nraw_max

    t_dim_can = 0
    t_dim_mixed = 0
    t_rank_est = 0
    t_norm = 0.0d0
    bridge_roundtrip_diff = huge(1.0d0)
    zcan_herm_diff = huge(1.0d0)
    pz_can_mixedz = 0.0d0
    pz_diff = -pz_prod_mixedz
    rel_pz_diff = pz_diff / max(abs(pz_prod_mixedz), 1.0d-300)
    canonical_operator_built = .false.
    convention_match = .false.
    bad = .false.

    if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
    if (.not. allocated(dg_frag%wpw_reduced_nraw)) return
    nmix = min(dg_frag%mixed_wannier_bpw_nmix, size(dg_frag%mixed_wannier_bpw_z, 2), &
      size(dg_frag%mixed_wannier_bpw_z, 3))
    if (nmix <= 0) return
    nraw_max = 0
    do i = 1, size(dg_frag%wpw_reduced_nraw)
      nraw_max = max(nraw_max, dg_frag%wpw_reduced_nraw(i))
    end do
    ncan = nraw_max
    t_dim_can = ncan
    t_dim_mixed = nmix

    ! Phase 2g-2d currently only identifies the dimensional bridge boundary.
    ! A valid production convention match requires an explicit C_can -> mixed
    ! coefficient map for W_self, P_self, and P_neighbor.  Do not infer it from
    ! the grid moment or from production cmix, because that would bypass the
    ! canonical operator being diagnosed here.
    t_rank_est = min(nmix, ncan)
    t_norm = sqrt(dble(max(0, t_rank_est)))
    call wpw_reduced_mixedz_operator_stats(dg_frag, zcan_herm_diff, pz_can_mixedz, i, bad)
    pz_can_mixedz = 0.0d0
    pz_diff = pz_can_mixedz - pz_prod_mixedz
    rel_pz_diff = pz_diff / max(abs(pz_prod_mixedz), 1.0d-300)

    if (ncan <= 0 .or. nmix <= 0) then
      bad = .true.
    else
      bridge_roundtrip_diff = 1.0d0
    end if
  end subroutine wpw_reduced_canon_mixedz_bridge_stats


  subroutine wpw_reduced_canon_pz_block_operator_stats(dg_frag, pz_prod_mixedz, dim_w, dim_p_self, &
      dim_p_neighbor, zww_norm, zwp_norm, zpp_norm, zherm_diff, pz_can_block, pz_diff, rel_pz_diff, bad)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: pz_prod_mixedz
    integer, intent(out) :: dim_w, dim_p_self, dim_p_neighbor
    real(8), intent(out) :: zww_norm, zwp_norm, zpp_norm, zherm_diff
    real(8), intent(out) :: pz_can_block, pz_diff, rel_pz_diff
    logical, intent(out) :: bad

    integer :: i, j, ispin, nmix, nw, np, iloc, nself, nraw
    real(8) :: zww_norm2, zwp_norm2, zpp_norm2

    dim_w = 0
    dim_p_self = max(0, dg_frag%n_plane_waves)
    dim_p_neighbor = 0
    zww_norm = 0.0d0
    zwp_norm = 0.0d0
    zpp_norm = 0.0d0
    zherm_diff = huge(1.0d0)
    pz_can_block = 0.0d0
    pz_diff = -pz_prod_mixedz
    rel_pz_diff = pz_diff / max(abs(pz_prod_mixedz), 1.0d-300)
    bad = .true.
    if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
    if (.not. allocated(dg_frag%wpw_reduced_nraw)) return
    if (.not. allocated(dg_frag%wpw_reduced_nself)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return

    nw = dg_frag%mixed_wannier_bpw_nw
    np = dg_frag%mixed_wannier_bpw_np
    nmix = min(dg_frag%mixed_wannier_bpw_nmix, size(dg_frag%mixed_wannier_bpw_z, 2), &
      size(dg_frag%mixed_wannier_bpw_z, 3))
    if (nw <= 0 .or. nmix <= 0 .or. nmix < nw) return
    dim_w = nw
    zww_norm2 = 0.0d0
    zwp_norm2 = 0.0d0
    zpp_norm2 = 0.0d0
    zherm_diff = 0.0d0
    do ispin = 1, size(dg_frag%mixed_wannier_bpw_z, 4)
      do i = 1, nmix
        do j = 1, nmix
          zherm_diff = max(zherm_diff, abs(dg_frag%mixed_wannier_bpw_z(3, i, j, ispin) - &
            conjg(dg_frag%mixed_wannier_bpw_z(3, j, i, ispin))))
          if (i <= nw .and. j <= nw) then
            zww_norm2 = zww_norm2 + abs(dg_frag%mixed_wannier_bpw_z(3, i, j, ispin))**2
          else if ((i <= nw .and. j > nw) .or. (i > nw .and. j <= nw)) then
            zwp_norm2 = zwp_norm2 + abs(dg_frag%mixed_wannier_bpw_z(3, i, j, ispin))**2
          else
            zpp_norm2 = zpp_norm2 + abs(dg_frag%mixed_wannier_bpw_z(3, i, j, ispin))**2
          end if
        end do
      end do
    end do
    zww_norm = sqrt(zww_norm2)
    zwp_norm = sqrt(zwp_norm2)
    zpp_norm = sqrt(zpp_norm2)
    do iloc = 1, size(dg_frag%wpw_reduced_nraw)
      nself = 0
      nraw = dg_frag%wpw_reduced_nraw(iloc)
      if (iloc <= size(dg_frag%wpw_reduced_nself)) nself = dg_frag%wpw_reduced_nself(iloc)
      if (nraw > nself) dim_p_neighbor = max(dim_p_neighbor, nraw - nself)
    end do

    ! The production mixed-Z block norms are available above.  The canonical
    ! WPW P blocks are not yet proven to be the same physical basis as the
    ! production BPW-perp P block, so do not form Pz_can_block or allow
    ! replacement in this phase.
    bad = (np > 0) .or. (dim_p_neighbor > 0)
  end subroutine wpw_reduced_canon_pz_block_operator_stats


  subroutine wpw_reduced_canon_p_projection_stats(dg_frag, pz_prod_mixedz, dim_p_prod, dim_p_self, &
      dim_p_neighbor, rank_pprod_pcan, norm_pself, norm_pneighbor, norm_pcan, proj_norm_pself, &
      proj_norm_pneighbor, proj_norm_pcan, leakage_norm, rel_leakage, overlap_pself_pneighbor_norm, &
      pz_can_projectedp, pz_diff_projectedp, rel_pz_diff_projectedp, projector_built, bad)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: pz_prod_mixedz
    integer, intent(out) :: dim_p_prod, dim_p_self, dim_p_neighbor, rank_pprod_pcan
    real(8), intent(out) :: norm_pself, norm_pneighbor, norm_pcan
    real(8), intent(out) :: proj_norm_pself, proj_norm_pneighbor, proj_norm_pcan
    real(8), intent(out) :: leakage_norm, rel_leakage, overlap_pself_pneighbor_norm
    real(8), intent(out) :: pz_can_projectedp, pz_diff_projectedp, rel_pz_diff_projectedp
    logical, intent(out) :: projector_built, bad

    integer :: iloc, nraw, nself, nw, npw, ps1, ps2, pn1, pn2
    real(8) :: norm_self2, norm_neighbor2, overlap2

    dim_p_prod = max(0, dg_frag%mixed_wannier_bpw_np)
    dim_p_self = max(0, dg_frag%n_plane_waves)
    dim_p_neighbor = 0
    rank_pprod_pcan = 0
    norm_pself = 0.0d0
    norm_pneighbor = 0.0d0
    norm_pcan = 0.0d0
    proj_norm_pself = 0.0d0
    proj_norm_pneighbor = 0.0d0
    proj_norm_pcan = 0.0d0
    leakage_norm = 0.0d0
    rel_leakage = 1.0d0
    overlap_pself_pneighbor_norm = 0.0d0
    pz_can_projectedp = 0.0d0
    pz_diff_projectedp = -pz_prod_mixedz
    rel_pz_diff_projectedp = pz_diff_projectedp / max(abs(pz_prod_mixedz), 1.0d-300)
    projector_built = .false.
    bad = .true.
    if (.not. allocated(dg_frag%wpw_reduced_Sraw_build)) return
    if (.not. allocated(dg_frag%wpw_reduced_nraw)) return
    if (.not. allocated(dg_frag%wpw_reduced_nself)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return

    norm_self2 = 0.0d0
    norm_neighbor2 = 0.0d0
    overlap2 = 0.0d0
    npw = max(0, dg_frag%n_plane_waves)
    do iloc = 1, size(dg_frag%wpw_reduced_nraw)
      nraw = dg_frag%wpw_reduced_nraw(iloc)
      nself = dg_frag%wpw_reduced_nself(iloc)
      nw = 0
      if (iloc <= size(dg_frag%global_wannier_local_nkeep)) nw = dg_frag%global_wannier_local_nkeep(iloc)
      if (nraw <= 0 .or. nself <= 0 .or. nw <= 0 .or. npw <= 0) cycle
      ps1 = nw + 1
      ps2 = min(nw + npw, nraw)
      if (ps2 >= ps1) norm_self2 = norm_self2 + &
        sum(abs(dg_frag%wpw_reduced_Sraw_build(ps1:ps2, ps1:ps2, iloc))**2)
      pn1 = nself + 1
      pn2 = nraw
      if (pn2 >= pn1) then
        dim_p_neighbor = max(dim_p_neighbor, pn2 - pn1 + 1)
        norm_neighbor2 = norm_neighbor2 + sum(abs(dg_frag%wpw_reduced_Sraw_build(pn1:pn2, pn1:pn2, iloc))**2)
        if (ps2 >= ps1) overlap2 = overlap2 + &
          sum(abs(dg_frag%wpw_reduced_Sraw_build(ps1:ps2, pn1:pn2, iloc))**2)
      end if
    end do
    norm_pself = sqrt(norm_self2)
    norm_pneighbor = sqrt(norm_neighbor2)
    overlap_pself_pneighbor_norm = sqrt(overlap2)
    norm_pcan = sqrt(norm_self2 + norm_neighbor2 + 2.0d0 * overlap2)

    ! No valid Pprod basis-overlap object is stored here yet.  The diagnostic
    ! intentionally keeps projector_built=F instead of collapsing P_neighbor
    ! onto production P by dimensional assumptions.
    leakage_norm = norm_pcan
    rel_leakage = leakage_norm / max(norm_pcan, 1.0d-300)
    rank_pprod_pcan = 0
    bad = dim_p_prod > 0 .and. norm_pcan > 0.0d0
  end subroutine wpw_reduced_canon_p_projection_stats


  subroutine wpw_reduced_prod_p_basis_save_stats(dg_frag, dim_p_prod, dim_p_raw_before_perp, &
      dim_p_after_perp, dim_p_after_qmat, s_herm_diff, s_min_eval, s_max_eval, s_cond, &
      p_basis_saved, transform_saved, metric_saved, bad)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(out) :: dim_p_prod, dim_p_raw_before_perp, dim_p_after_perp, dim_p_after_qmat
    real(8), intent(out) :: s_herm_diff, s_min_eval, s_max_eval, s_cond
    logical, intent(out) :: p_basis_saved, transform_saved, metric_saved, bad

    complex(8), allocatable :: Sinv(:,:)
    integer :: info, nkeep, ispin, nspin_p
    real(8) :: smin, smax, hdiff

    dim_p_prod = max(0, dg_frag%mixed_wannier_bpw_np)
    dim_p_raw_before_perp = max(0, dg_frag%mixed_wannier_bpw_praw_dim)
    dim_p_after_perp = dim_p_prod
    dim_p_after_qmat = dim_p_prod
    s_herm_diff = huge(1.0d0)
    s_min_eval = 0.0d0
    s_max_eval = 0.0d0
    s_cond = huge(1.0d0)
    transform_saved = allocated(dg_frag%mixed_wannier_bpw_p_transform)
    metric_saved = allocated(dg_frag%mixed_wannier_bpw_p_metric)
    p_basis_saved = dg_frag%has_mixed_wannier_bpw_p_basis .and. transform_saved .and. metric_saved
    bad = .true.
    if (.not. p_basis_saved) return
    if (dim_p_prod <= 0) return
    if (size(dg_frag%mixed_wannier_bpw_p_metric, 1) < dim_p_prod .or. &
        size(dg_frag%mixed_wannier_bpw_p_metric, 2) < dim_p_prod) return

    allocate(Sinv(dim_p_prod, dim_p_prod))
    s_herm_diff = 0.0d0
    s_min_eval = huge(1.0d0)
    s_max_eval = 0.0d0
    nspin_p = size(dg_frag%mixed_wannier_bpw_p_metric, 3)
    do ispin = 1, nspin_p
      call wpw_local_herm_max(dg_frag%mixed_wannier_bpw_p_metric(1:dim_p_prod,1:dim_p_prod,ispin), &
        dim_p_prod, hdiff)
      s_herm_diff = max(s_herm_diff, hdiff)
      call build_hermitian_pseudoinverse(dg_frag%mixed_wannier_bpw_p_metric(1:dim_p_prod,1:dim_p_prod,ispin), &
        dim_p_prod, 1.0d-10, Sinv, info, smin, smax, nkeep)
      if (info /= 0) cycle
      s_min_eval = min(s_min_eval, smin)
      s_max_eval = max(s_max_eval, smax)
    end do
    deallocate(Sinv)
    if (s_min_eval == huge(1.0d0)) then
      s_min_eval = 0.0d0
      s_max_eval = 0.0d0
      s_cond = huge(1.0d0)
      return
    end if
    s_cond = s_max_eval / max(s_min_eval, 1.0d-300)
    bad = (s_min_eval <= 0.0d0) .or. (s_herm_diff > 1.0d-8) .or. &
      (s_min_eval /= s_min_eval) .or. (s_max_eval /= s_max_eval)
  end subroutine wpw_reduced_prod_p_basis_save_stats

end module rt_dg_wpw_reduced_mixedz_diag
