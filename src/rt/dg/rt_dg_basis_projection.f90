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
! Basis projection and stabilization utilities for DG-Fragment RT-TDDFT
!=======================================================================

module rt_dg_basis_projection
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  implicit none


  private
  public :: calculate_new_old_basis_overlap
  public :: stabilize_basis_overlap
  public :: project_wavefunction_to_new_basis
  public :: reorthonormalize_occupied_subspace
  public :: diagonalize_and_update_basis

contains

  !=======================================================================
  ! Calculate overlap matrix between new and old basis functions
  ! S_ji = <phi_new_j|phi_old_i>
  !=======================================================================
  subroutine calculate_new_old_basis_overlap(dg_frag, phi_frag_old)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    real(8), intent(in) :: phi_frag_old(:,:,:,:,:)  ! Real basis from DC-LCFO

    integer :: ifrag, i_local, istate_new, istate_old, ispin
    integer :: ix, iy, iz
    integer :: ix_l, ix_u, iy_l, iy_u, iz_l, iz_u
    real(8) :: overlap_sum
    real(8) :: kahan_c, term, y, t
    real(8) :: hvol

    hvol = dg_frag%hgs(1) * dg_frag%hgs(2) * dg_frag%hgs(3)

    do ispin = 1, dg_frag%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        ! Use only physical domain for overlap evaluation (interior, excluding halo).
        ! BUG FIX: dg_frag%lg%is/ie are global grid bounds (e.g., [1,1,Nx]).
        ! phi_frag is fragment-local with bounds (1-nb:nxyz_domain+nb).
        ! Accessing phi_frag(ix,...) with ix up to Nx caused out-of-bounds errors.
        ! Fix: use fragment interior 1:nxyz_domain as loop bounds.
        ix_l = 1
        ix_u = dg_frag%nxyz_domain(1, ifrag)
        iy_l = 1
        iy_u = dg_frag%nxyz_domain(2, ifrag)
        iz_l = 1
        iz_u = dg_frag%nxyz_domain(3, ifrag)

        do istate_new = 1, dg_frag%nstate_frag
          do istate_old = 1, dg_frag%nstate_frag

            overlap_sum = 0.0d0
            kahan_c = 0.0d0

            do iz = iz_l, iz_u
              do iy = iy_l, iy_u
                do ix = ix_l, ix_u
                  term = dg_frag%phi_frag(ix, iy, iz, istate_new, i_local) * &
                        phi_frag_old(ix, iy, iz, istate_old, i_local) * hvol
                  y = term - kahan_c
                  t = overlap_sum + y
                  kahan_c = (t - overlap_sum) - y
                  overlap_sum = t
                end do
              end do
            end do

            dg_frag%basis_overlap(istate_new, istate_old, ispin) = overlap_sum

          end do
        end do

      end do
    end do

  end subroutine calculate_new_old_basis_overlap

  !=======================================================================
  ! Stabilize basis overlap matrix for projection
  ! 1) Sign alignment per new basis state
  ! 2) Procrustes rotation via SVD: Q = U * V^T
  !=======================================================================
  subroutine stabilize_basis_overlap(dg_frag, system)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: dt
    use phys_constants, only: au_fs
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system

    integer :: ispin, n, n_sel, irow, icol, info, lwork, iunit
    real(8) :: max_abs, sign_val
    real(8) :: sigma_clipped, theta, theta_max, theta_mean, q_minus_i_frob
    real(8) :: proj_residual_fro, a_norm_fro
    real(8) :: time_fs
    real(8), allocatable :: A_sel(:,:), A_orig(:,:), U(:,:), VT(:,:), S(:), work(:)
    real(8) :: work_query(1)
    logical :: file_exists
    external :: dgesvd

    n = dg_frag%nstate_frag
    if (.not. allocated(dg_frag%basis_overlap)) return

    n_sel = min(system%no, n)
    if (n_sel <= 0) return

    allocate(A_sel(n_sel, n_sel))
    allocate(A_orig(n_sel, n_sel))
    allocate(U(n_sel, n_sel))
    allocate(VT(n_sel, n_sel))
    allocate(S(n_sel))

    do ispin = 1, dg_frag%nspin
      A_sel(:, :) = dg_frag%basis_overlap(1:n_sel, 1:n_sel, ispin)

      do irow = 1, n_sel
        max_abs = 0.0d0
        sign_val = 1.0d0
        do icol = 1, n_sel
          if (abs(A_sel(irow, icol)) > max_abs) then
            max_abs = abs(A_sel(irow, icol))
            if (A_sel(irow, icol) < 0.0d0) then
              sign_val = -1.0d0
            else
              sign_val = 1.0d0
            end if
          end if
        end do
        if (sign_val < 0.0d0) then
          A_sel(irow, :) = -A_sel(irow, :)
        end if
      end do
      A_orig(:, :) = A_sel(:, :)

      lwork = -1
      call dgesvd('A', 'A', n_sel, n_sel, A_sel, n_sel, S, U, n_sel, VT, n_sel, work_query, lwork, info)
      if (info == 0) then
        lwork = max(1, int(work_query(1)))
        allocate(work(lwork))
        call dgesvd('A', 'A', n_sel, n_sel, A_sel, n_sel, S, U, n_sel, VT, n_sel, work, lwork, info)
        deallocate(work)
      end if

      if (info == 0) then
        theta_max = 0.0d0
        theta_mean = 0.0d0
        do irow = 1, n_sel
          sigma_clipped = min(1.0d0, max(0.0d0, S(irow)))
          theta = acos(sigma_clipped)
          theta_mean = theta_mean + theta
          if (theta > theta_max) theta_max = theta
        end do
        theta_mean = theta_mean / dble(n_sel)

        A_sel(:, :) = matmul(U, VT)
        q_minus_i_frob = 0.0d0
        proj_residual_fro = 0.0d0
        a_norm_fro = 0.0d0
        do irow = 1, n_sel
          do icol = 1, n_sel
            if (irow == icol) then
              q_minus_i_frob = q_minus_i_frob + (A_sel(irow, icol) - 1.0d0)**2
            else
              q_minus_i_frob = q_minus_i_frob + A_sel(irow, icol)**2
            end if
            proj_residual_fro = proj_residual_fro + (A_orig(irow, icol) - A_sel(irow, icol))**2
            a_norm_fro = a_norm_fro + A_orig(irow, icol)**2
          end do
        end do
        q_minus_i_frob = sqrt(q_minus_i_frob)
        proj_residual_fro = sqrt(proj_residual_fro)
        a_norm_fro = sqrt(a_norm_fro)
        if (a_norm_fro > 1.0d-14) then
          proj_residual_fro = proj_residual_fro / a_norm_fro
        else
          proj_residual_fro = 0.0d0
        end if

        dg_frag%basis_overlap(:, :, ispin) = 0.0d0
        dg_frag%basis_overlap(1:n_sel, 1:n_sel, ispin) = A_sel(:, :)
        do irow = n_sel + 1, n
          dg_frag%basis_overlap(irow, irow, ispin) = 1.0d0
        end do

        if (nproc_id_global == 0) then
          time_fs = dble(dg_frag%last_basis_update_step) * dt * au_fs
          write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
            "[BASIS-METRIC] update=", dg_frag%nbasis_update_count, &
            " step=", dg_frag%last_basis_update_step, &
            " spin=", ispin, " theta_max(rad)=", theta_max, &
            " theta_mean(rad)=", theta_mean, " ||Q-I||_F=", q_minus_i_frob, &
            " proj_res(F)=", proj_residual_fro

          inquire(file='basis_update_metrics.csv', exist=file_exists)
          open(newunit=iunit, file='basis_update_metrics.csv', status='unknown', &
               position='append', action='write')
          if (.not. file_exists) then
            write(iunit,'(a)') 'update_count,step,time_fs,spin,theta_max_rad,theta_mean_rad,q_minus_i_fro,proj_residual_fro'
          end if
          write(iunit,'(i0,a,i0,a,es23.15,a,i0,a,es23.15,a,es23.15,a,es23.15,a,es23.15)') &
            dg_frag%nbasis_update_count, ',', dg_frag%last_basis_update_step, ',', time_fs, ',', ispin, ',', &
            theta_max, ',', theta_mean, ',', q_minus_i_frob, ',', proj_residual_fro
          close(iunit)
        end if
      end if
    end do

    deallocate(A_sel, A_orig, U, VT, S)

  end subroutine stabilize_basis_overlap

  !=======================================================================
  ! Project current wave function onto new basis
  !=======================================================================
  subroutine project_wavefunction_to_new_basis(dg_frag, system)
    use structures
    use communication, only: comm_is_root
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system

    integer :: ispin, io, istate_new, istate_old, nbasis_coef
    integer :: ix, iy, iz, info
    integer :: ix_l, ix_u, iy_l, iy_u, iz_l, iz_u
    complex(8), allocatable :: coef_new(:,:,:)
    real(8), allocatable :: s_new(:,:), trans_mat(:,:), rhs(:,:)
    integer, allocatable :: ipiv(:)
    real(8) :: norm_check
    real(8) :: overlap_sum
    real(8) :: kahan_c, term, y, t
    real(8) :: hvol
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    external :: dgesv

    nbasis_coef = size(dg_frag%coef, 1)
    allocate(coef_new(nbasis_coef, dg_frag%nstate_tot, dg_frag%nspin))
    coef_new = (0.0d0, 0.0d0)
    allocate(s_new(dg_frag%nstate_frag, dg_frag%nstate_frag))
    allocate(trans_mat(dg_frag%nstate_frag, dg_frag%nstate_frag))
    allocate(rhs(dg_frag%nstate_frag, dg_frag%nstate_frag))
    allocate(ipiv(dg_frag%nstate_frag))

    hvol = dg_frag%hgs(1) * dg_frag%hgs(2) * dg_frag%hgs(3)
    ! BUG FIX: dg_frag%lg%is/ie are global grid bounds (e.g., [1,1,Nx]).
    ! phi_frag is fragment-local with bounds (1-nb:nxyz_domain+nb).
    ! Use interior bounds 1:nxyz_domain for the first local fragment (index 1).
    ix_l = 1
    ix_u = dg_frag%nxyz_domain(1, dg_frag%ifrag_start)
    iy_l = 1
    iy_u = dg_frag%nxyz_domain(2, dg_frag%ifrag_start)
    iz_l = 1
    iz_u = dg_frag%nxyz_domain(3, dg_frag%ifrag_start)

    do ispin = 1, dg_frag%nspin
      allocate(coef_frag_all(nbasis_coef, dg_frag%nstate_tot))
      coef_frag_all(:, :) = dg_frag%coef(1:nbasis_coef, 1:dg_frag%nstate_tot, ispin)
      s_new = 0.0d0
!$omp parallel do private(istate_new,istate_old,iz,iy,ix,overlap_sum) schedule(static)
      do istate_old = 1, dg_frag%nstate_frag
        do istate_new = 1, istate_old
          overlap_sum = 0.0d0
          do iz = iz_l, iz_u
            do iy = iy_l, iy_u
              do ix = ix_l, ix_u
                overlap_sum = overlap_sum + &
                  dg_frag%phi_frag(ix, iy, iz, istate_new, 1) * &
                  dg_frag%phi_frag(ix, iy, iz, istate_old, 1) * hvol
              end do
            end do
          end do
          s_new(istate_new, istate_old) = overlap_sum
          if (istate_new /= istate_old) s_new(istate_old, istate_new) = overlap_sum
        end do
      end do
!$omp end parallel do

      trans_mat(:, :) = s_new(:, :)
      rhs(:, :) = dg_frag%basis_overlap(:, :, ispin)
      do istate_new = 1, dg_frag%nstate_frag
        trans_mat(istate_new, istate_new) = trans_mat(istate_new, istate_new) + 1.0d-12
      end do
      call dgesv(dg_frag%nstate_frag, dg_frag%nstate_frag, trans_mat, dg_frag%nstate_frag, &
                 ipiv, rhs, dg_frag%nstate_frag, info)
      if (info /= 0) then
        rhs(:, :) = dg_frag%basis_overlap(:, :, ispin)
      end if

!$omp parallel do collapse(2) private(io,istate_new,istate_old)
      do io = 1, dg_frag%nstate_tot
        do istate_new = 1, dg_frag%nstate_frag
          do istate_old = 1, dg_frag%nstate_frag
            coef_new(istate_new, io, ispin) = coef_new(istate_new, io, ispin) + &
              rhs(istate_new, istate_old) * &
              coef_frag_all(istate_old, io)
          end do
        end do
      end do
!$omp end parallel do
      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
    end do

    dg_frag%coef(:, :, :) = coef_new(:, :, :)
    dg_frag%coef_new(:, :, :) = coef_new(:, :, :)
    call zero_nonowned_coefficients(dg_frag)

    call reorthonormalize_occupied_subspace(dg_frag, system)

    if (comm_is_root(dg_frag%icomm)) then
      do ispin = 1, dg_frag%nspin
        norm_check = 0.0d0
        kahan_c = 0.0d0
        do io = 1, min(5, dg_frag%nstate_tot)
          do istate_new = 1, dg_frag%nstate_frag
            term = abs(coef_new(istate_new, io, ispin))**2
            y = term - kahan_c
            t = norm_check + y
            kahan_c = (t - norm_check) - y
            norm_check = t
          end do
          write(*,'(1x,a,i0,a,i0,a,f10.6)') "  Norm check orbital ", io, &
                " spin ", ispin, ": ", norm_check
          norm_check = 0.0d0
        end do
      end do
    end if

    deallocate(coef_new, s_new, trans_mat, rhs, ipiv)

  end subroutine project_wavefunction_to_new_basis

  !=======================================================================
  ! Re-orthonormalize occupied subspace after projection
  !=======================================================================
  subroutine reorthonormalize_occupied_subspace(dg_frag, system)
    use structures
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system

    integer :: ispin, iocc, jocc, n_occ, n_basis
    complex(8) :: proj
    real(8) :: norm_val
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)

    n_basis = dg_frag%nstate_frag
    n_occ = min(system%no, dg_frag%nstate_tot)
    if (n_occ <= 0) return

    do ispin = 1, dg_frag%nspin
      allocate(coef_frag_all(size(dg_frag%coef,1), n_occ))
      coef_frag_all(:, :) = dg_frag%coef(:, 1:n_occ, ispin)
      do iocc = 1, n_occ
        do jocc = 1, iocc - 1
          proj = sum(conjg(coef_frag_all(:, jocc)) * coef_frag_all(:, iocc))
          coef_frag_all(:, iocc) = coef_frag_all(:, iocc) - proj * coef_frag_all(:, jocc)
        end do

        norm_val = sqrt(sum(abs(coef_frag_all(:, iocc))**2))
        if (norm_val > 1.0d-12) then
          coef_frag_all(:, iocc) = coef_frag_all(:, iocc) / norm_val
        end if
      end do
      dg_frag%coef(:, 1:n_occ, ispin) = coef_frag_all(:, 1:n_occ)
      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
    end do

    dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    call zero_nonowned_coefficients(dg_frag)

  end subroutine reorthonormalize_occupied_subspace

  !=======================================================================
  ! Diagonalize Hamiltonian matrix and rotate coefficients
  !=======================================================================
  subroutine diagonalize_and_update_basis(dg_frag, system)
    use structures
    use communication, only: comm_is_root
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients, copy_matrix_blocks_metric_to_real_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system

    integer :: ispin, n, lda, lwork, info, i, j, io
    real(8), allocatable :: H_work(:,:), eigenvalues_tmp(:), eigenvectors(:,:), work(:)
    complex(8), allocatable :: coef_old(:,:,:), coef_new(:,:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    logical :: use_mixed_basis
    external :: dsyev

    use_mixed_basis = dg_frag%use_plane_wave_basis .and. (dg_frag%n_plane_waves > 0)

    if (use_mixed_basis) then
      if (comm_is_root(dg_frag%icomm)) then
        write(*,*) "ERROR: diagonalize_and_update_basis cannot be used for DG+PW."
        write(*,*) "       Mixed-basis updates require diagonalize_mixed_basis with current potentials."
        write(*,*) "       Refusing raw fragment-only fallback during DG+PW basis update."
      end if
      stop "DG+PW basis update requires mixed-basis diagonalization"
    end if

    do ispin = 1, dg_frag%nspin
      n = dg_frag%nstate_frag
      lda = n

      allocate(H_work(n, n))
      allocate(eigenvalues_tmp(n))
      allocate(eigenvectors(n, n))

      if (allocated(dg_frag%H_mat)) then
        H_work(:, :) = dg_frag%H_mat(:, :, ispin)
      else if (allocated(dg_frag%H_mat_blocks)) then
        H_work(:, :) = 0.0d0
        call copy_matrix_blocks_metric_to_real_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n, H_work)
      else
        write(*,*) "ERROR: Hamiltonian matrix unavailable in diagonalize_and_update_basis"
        stop
      end if

      lwork = -1
      allocate(work(1))
      call DSYEV('V', 'U', n, H_work, lda, eigenvalues_tmp, &
                 work, lwork, info)
      lwork = int(work(1)) + 1
      deallocate(work)
      allocate(work(lwork))

      call DSYEV('V', 'U', n, H_work, lda, eigenvalues_tmp, &
                 work, lwork, info)

      if (info /= 0) then
        write(*,*) "ERROR: Hamiltonian diagonalization failed in diagonalize_and_update_basis, info=", info
        stop
      end if

      dg_frag%eigenvalues(:, ispin) = eigenvalues_tmp(:)
      eigenvectors(:, :) = H_work(:, :)

      allocate(coef_old(n, dg_frag%nstate_tot, 1))
      allocate(coef_new(n, dg_frag%nstate_tot, 1))

      allocate(coef_frag_all(n, dg_frag%nstate_tot))
      coef_frag_all(:, :) = dg_frag%coef(1:n, 1:dg_frag%nstate_tot, ispin)
      coef_old(:, :, 1) = coef_frag_all(:, :)
      coef_new = 0.0d0

      do io = 1, dg_frag%nstate_tot
        do j = 1, n
          do i = 1, n
            coef_new(j, io, 1) = coef_new(j, io, 1) + &
                                 eigenvectors(i, j) * coef_old(i, io, 1)
          end do
        end do
      end do

      dg_frag%coef(:, :, ispin) = coef_new(:, :, 1)
      dg_frag%coef_new(:, :, ispin) = coef_new(:, :, 1)
      call zero_nonowned_coefficients(dg_frag)

      deallocate(H_work, eigenvalues_tmp, eigenvectors, work)
      deallocate(coef_old, coef_new)
      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
      if (comm_is_root(dg_frag%icomm)) then
        write(*,'(1x,a,i0,a)') "    Spin ", ispin, " diagonalized successfully"
      end if

    end do

    if (comm_is_root(dg_frag%icomm)) then
      write(*,*) "  Coefficients transformed to new eigenbasis"
      write(*,*) "  (Note: Basis space not expanded - limitations apply)"
    end if

  end subroutine diagonalize_and_update_basis

end module rt_dg_basis_projection
