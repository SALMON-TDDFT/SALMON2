! DG fixed-density Flux-SCF initial-state relaxation
#include "config.h"
module rt_dg_flux_scf
  use structures
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, apply_matrix_blocks_batch, &
    apply_complex_matrix_blocks_batch, zero_nonowned_coefficients
  use communication, only: comm_summation
  use eigen_subdiag_sub, only: eigen_zheev
  implicit none
  private
  public :: prepare_fragment_flux_scf_coefficients

contains

  subroutine prepare_fragment_flux_scf_coefficients(dg_frag, system, mg, ppg, did_prepare, coef_delta)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    logical,                intent(out)   :: did_prepare
    real(8), optional,      intent(out)   :: coef_delta

    integer, parameter :: cg_steps_default = 5
    real(8), parameter :: cg_step_default = 1.0d-2
    real(8) :: Ac_zero(3)

    did_prepare = .false.
    if (present(coef_delta)) coef_delta = 0.0d0
    if (.not. dg_frag%fragment_basis_contracted) return
    if (dg_frag%nspin /= 1) return
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) return
    if (.not. dg_frag%parallel_mode_orbital) return
    if (.not. allocated(dg_frag%coef)) return
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%S_mat_blocks)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%coef_owner)) return

    Ac_zero(:) = 0.0d0
    call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_zero)
    call relax_flux_coefficients_fixed_basis(dg_frag, cg_steps_default, cg_step_default, &
                                             did_prepare, coef_delta)
  end subroutine prepare_fragment_flux_scf_coefficients

  subroutine relax_flux_coefficients_fixed_basis(dg_frag, cg_steps, cg_step, did_prepare, coef_delta_out)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: cg_steps
    real(8), intent(in) :: cg_step
    logical, intent(out) :: did_prepare
    real(8), optional, intent(out) :: coef_delta_out

    integer :: ispin, iter, c0, c1, nb, nvalid, nstate_use, ist
    integer, parameter :: batch_size = 64
    real(8) :: step_eff, den, norm_val
    real(8) :: coef_delta, rel_res
    real(8) :: delta_local(2), delta_global(2)
    real(8) :: res_local(2), res_global(2)
    logical :: ortho_ok
    real(8), allocatable :: eig_num_l(:), eig_den_l(:), eig_num_g(:), eig_den_g(:)
    real(8), allocatable :: norm_l(:), norm_g(:)
    complex(8), allocatable :: c_old(:,:), c_new(:,:), hc(:,:), sc(:,:), residual(:,:)
    complex(8), allocatable :: coef_iter_start(:,:,:)

    did_prepare = .false.
    if (present(coef_delta_out)) coef_delta_out = 0.0d0
    nvalid = size(dg_frag%coef, 1)
    nstate_use = min(dg_frag%nstate_tot, size(dg_frag%coef, 2))
    if (nvalid <= 0 .or. nstate_use <= 0) return
    step_eff = max(1.0d-8, min(1.0d0, cg_step))
    coef_delta = 0.0d0

    if (dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') &
        '[DG-FLUX-BASIS] fixed-basis coefficient CG: nstate=', nstate_use, &
        ' steps=', max(1, cg_steps), ' step=', step_eff
      flush(6)
    end if

    do iter = 1, max(1, cg_steps)
      delta_local(:) = 0.0d0
      res_local(:) = 0.0d0
      allocate(coef_iter_start(nvalid, nstate_use, dg_frag%nspin))
      coef_iter_start(:, :, :) = dg_frag%coef(:, 1:nstate_use, :)
      do ispin = 1, dg_frag%nspin
        do c0 = 1, nstate_use, batch_size
          c1 = min(nstate_use, c0 + batch_size - 1)
          nb = c1 - c0 + 1
          allocate(c_old(nvalid, nb), c_new(nvalid, nb), hc(nvalid, nb), sc(nvalid, nb))
          allocate(residual(nvalid, nb))
          allocate(eig_num_l(nb), eig_den_l(nb), eig_num_g(nb), eig_den_g(nb))
          allocate(norm_l(nb), norm_g(nb))

          c_old(:, :) = dg_frag%coef(:, c0:c1, ispin)
          call apply_flux_hamiltonian(dg_frag, ispin, c_old, hc)
          sc(:, :) = (0.0d0, 0.0d0)
          call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, c_old, sc)

          eig_num_l(:) = 0.0d0
          eig_den_l(:) = 0.0d0
          do ist = 1, nb
            eig_num_l(ist) = real(sum(conjg(c_old(:, ist)) * hc(:, ist)), kind=8)
            eig_den_l(ist) = real(sum(conjg(c_old(:, ist)) * sc(:, ist)), kind=8)
          end do
          call comm_summation(eig_num_l, eig_num_g, nb, dg_frag%icomm)
          call comm_summation(eig_den_l, eig_den_g, nb, dg_frag%icomm)

          residual(:, :) = hc(:, :)
          do ist = 1, nb
            den = eig_den_g(ist)
            if (abs(den) < 1.0d-300) den = sign(1.0d-300, den + 1.0d-300)
            if (allocated(dg_frag%esp) .and. c0 + ist - 1 <= size(dg_frag%esp, 1)) &
              dg_frag%esp(c0 + ist - 1, ispin) = eig_num_g(ist) / den
            residual(:, ist) = residual(:, ist) - (eig_num_g(ist) / den) * sc(:, ist)
          end do
          res_local(1) = res_local(1) + sum(abs(residual(:, :))**2)
          res_local(2) = res_local(2) + sum(abs(hc(:, :))**2)

          c_new(:, :) = c_old(:, :) - step_eff * residual(:, :)
          sc(:, :) = (0.0d0, 0.0d0)
          call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, c_new, sc)
          norm_l(:) = 0.0d0
          do ist = 1, nb
            norm_l(ist) = real(sum(conjg(c_new(:, ist)) * sc(:, ist)), kind=8)
          end do
          call comm_summation(norm_l, norm_g, nb, dg_frag%icomm)
          do ist = 1, nb
            norm_val = max(norm_g(ist), 1.0d-300)
            c_new(:, ist) = c_new(:, ist) / sqrt(norm_val)
          end do

          dg_frag%coef(:, c0:c1, ispin) = c_new(:, :)

          deallocate(c_old, c_new, hc, sc, residual)
          deallocate(eig_num_l, eig_den_l, eig_num_g, eig_den_g, norm_l, norm_g)
        end do
      end do

      call zero_nonowned_coefficients(dg_frag)
      do ispin = 1, dg_frag%nspin
        call s_orthonormalize_flux_coef_columns(dg_frag, ispin, nstate_use, ortho_ok)
        if (.not. ortho_ok) then
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0)') '[WARN] DG-FLUX-BASIS-CG S-orthonormalization failed: step=', iter, &
              ' ispin=', ispin
            flush(6)
          end if
          deallocate(coef_iter_start)
          return
        end if
      end do
      call zero_nonowned_coefficients(dg_frag)
      delta_local(1) = sum(abs(dg_frag%coef(:, 1:nstate_use, :) - coef_iter_start(:, :, :))**2)
      delta_local(2) = sum(abs(coef_iter_start(:, :, :))**2)
      deallocate(coef_iter_start)
      delta_global(:) = 0.0d0
      res_global(:) = 0.0d0
      call comm_summation(delta_local, delta_global, 2, dg_frag%icomm)
      call comm_summation(res_local, res_global, 2, dg_frag%icomm)
      coef_delta = sqrt(max(0.0d0, delta_global(1)) / max(delta_global(2), 1.0d-300))
      rel_res = sqrt(max(0.0d0, res_global(1)) / max(res_global(2), 1.0d-300))
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,2(a,1pe13.5))') '[DG-FLUX-BASIS-CG] step=', iter, &
          ' coef_delta=', coef_delta, ' rel_res=', rel_res
        flush(6)
      end if
    end do

    if (allocated(dg_frag%coef_work)) dg_frag%coef_work(:, :, :) = dg_frag%coef(:, :, :)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    if (present(coef_delta_out)) coef_delta_out = coef_delta
    did_prepare = .true.
  end subroutine relax_flux_coefficients_fixed_basis

  subroutine s_orthonormalize_flux_coef_columns(dg_frag, ispin, nstate, success)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin, nstate
    logical, intent(out) :: success

    integer, parameter :: panel_size = 32, trail_chunk = 512
    integer :: ipanel, panel_end, bsz, itrail, trail_end, nts, nrow, i, j
    real(8), parameter :: norm_floor = 1.0d-24
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: q(:,:), sq(:,:), gram_local(:,:), gram_global(:,:)
    complex(8), allocatable :: gram_panel(:,:), evec(:,:), evec_scaled(:,:), inv_sqrt(:,:)
    complex(8), allocatable :: proj_local(:,:), proj_global(:,:)

    success = .false.
    if (nstate <= 0) return
    if (.not. allocated(dg_frag%coef)) return
    nrow = size(dg_frag%coef, 1)
    if (nrow <= 0) return

    allocate(q(nrow, panel_size), sq(nrow, panel_size))

    do ipanel = 1, nstate, panel_size
      panel_end = min(nstate, ipanel + panel_size - 1)
      bsz = panel_end - ipanel + 1

      allocate(gram_local(bsz, bsz), gram_global(bsz, bsz))
      q(:, :) = (0.0d0, 0.0d0)
      sq(:, :) = (0.0d0, 0.0d0)
      q(:, 1:bsz) = dg_frag%coef(:, ipanel:panel_end, ispin)
      call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, q(:, 1:bsz), sq(:, 1:bsz))

      gram_local(:, :) = (0.0d0, 0.0d0)
      do j = 1, bsz
        do i = 1, bsz
          gram_local(i, j) = sum(conjg(q(:, i)) * sq(:, j))
        end do
      end do
      call comm_summation(gram_local, gram_global, bsz * bsz, dg_frag%icomm)
      do j = 1, bsz
        do i = j + 1, bsz
          gram_global(i, j) = 0.5d0 * (gram_global(i, j) + conjg(gram_global(j, i)))
          gram_global(j, i) = conjg(gram_global(i, j))
        end do
        gram_global(j, j) = cmplx(real(gram_global(j, j), kind=8), 0.0d0, kind=8)
      end do

      allocate(gram_panel(bsz, bsz), evec(bsz, bsz), evec_scaled(bsz, bsz), inv_sqrt(bsz, bsz), eval(bsz))
      gram_panel(:, :) = gram_global(1:bsz, 1:bsz)
      call eigen_zheev(gram_panel, eval, evec)
      if (minval(eval(1:bsz)) <= norm_floor) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,1pe13.5)') &
            '[WARN] DG-FLUX-BASIS-CG S-orthonormalization found null subspace: rank=', dg_frag%id, &
            ' ispin=', ispin, ' first_col=', ipanel, ' min_norm=', minval(eval(1:bsz))
          flush(6)
        end if
        deallocate(gram_panel, evec, evec_scaled, inv_sqrt, eval)
        deallocate(q, sq, gram_local, gram_global)
        return
      end if
      evec_scaled(:, :) = evec(:, :)
      do i = 1, bsz
        evec_scaled(:, i) = evec_scaled(:, i) / sqrt(eval(i))
      end do
      inv_sqrt(:, :) = matmul(evec_scaled(:, :), conjg(transpose(evec(:, :))))
      q(:, 1:bsz) = matmul(q(:, 1:bsz), inv_sqrt(:, :))
      dg_frag%coef(:, ipanel:panel_end, ispin) = q(:, 1:bsz)
      deallocate(gram_panel, evec, evec_scaled, inv_sqrt, eval)
      deallocate(gram_local, gram_global)

      sq(:, :) = (0.0d0, 0.0d0)
      call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, q(:, 1:bsz), sq(:, 1:bsz))
      do itrail = panel_end + 1, nstate, trail_chunk
        trail_end = min(nstate, itrail + trail_chunk - 1)
        nts = trail_end - itrail + 1
        allocate(proj_local(bsz, nts), proj_global(bsz, nts))
        proj_local(:, :) = matmul(conjg(transpose(sq(:, 1:bsz))), &
                                  dg_frag%coef(:, itrail:trail_end, ispin))
        call comm_summation(proj_local, proj_global, bsz * nts, dg_frag%icomm)
        dg_frag%coef(:, itrail:trail_end, ispin) = dg_frag%coef(:, itrail:trail_end, ispin) - &
          matmul(q(:, 1:bsz), proj_global)
        deallocate(proj_local, proj_global)
      end do
    end do

    deallocate(q, sq)
    success = .true.
  end subroutine s_orthonormalize_flux_coef_columns

  subroutine apply_flux_hamiltonian(dg_frag, ispin, coef_in, hcoef_out)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: coef_in(:, :)
    complex(8), intent(out) :: hcoef_out(:, :)

    hcoef_out(:, :) = (0.0d0, 0.0d0)
    call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_in, hcoef_out)
    if (allocated(dg_frag%H_nl_blocks)) &
      call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, coef_in, hcoef_out)
  end subroutine apply_flux_hamiltonian

end module rt_dg_flux_scf
