module rt_dg_fragment_ops
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  implicit none

  private
  public :: ensure_nonlocal_pp_matrix_A
  public :: calculate_microscopic_current_dg
  public :: build_spatial_A_coupling_matrices
  public :: apply_gradient_to_basis

contains

  !=======================================================================
  ! Build non-local pseudopotential matrix with vector potential A(t)
  ! V_NL(A) = e^{-i A·r} V_NL e^{i A·r}
  !
  ! Uses the SALMON approximation: A is nearly constant within PP cutoff
  !=======================================================================
  subroutine build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, nspin, hvol, Ac_tot, use_micro_A, Ac_micro, H_nl)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    integer,                intent(in) :: nspin
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    logical,                intent(in) :: use_micro_A
    real(8),                intent(in), optional :: Ac_micro(:, :, :, :)
    complex(8),             intent(out) :: H_nl(dg_frag%n_mat_max, dg_frag%n_mat_max, nspin)

    integer :: ifrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz, ig_i, ig_j, nbf
    integer :: iorg(3), ndom(3), lx, ly, lz
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase
    real(8) :: A_local(3)
    complex(8), allocatable :: uVpsi(:,:,:)  ! (nstate_frag, Nlma, nspin)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib

    if (ppg%Nlma == 0) then
      H_nl = (0.0d0, 0.0d0)
      return
    end if

    is = mg%is
    ie = mg%ie
    H_nl = (0.0d0, 0.0d0)

    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      uVpsi = (0.0d0, 0.0d0)

      do ispin = 1, nspin
        !$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, x, y, z, phase, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            ia = ppg%ia_tbl(ilma)
            overlap_i = (0.0d0, 0.0d0)
            do j = 1, ppg%mps(ia)
              ix = ppg%jxyz(1, j, ia)
              iy = ppg%jxyz(2, j, ia)
              iz = ppg%jxyz(3, j, ia)

              if (ix >= is(1) .and. ix <= ie(1) .and. &
                  iy >= is(2) .and. iy <= ie(2) .and. &
                  iz >= is(3) .and. iz <= ie(3)) then
                lx = mod(ix - iorg(1), dg_frag%lgnum_total(1)) + 1
                ly = mod(iy - iorg(2), dg_frag%lgnum_total(2)) + 1
                lz = mod(iz - iorg(3), dg_frag%lgnum_total(3)) + 1
                if (lx < 1 .or. lx > ndom(1)) cycle
                if (ly < 1 .or. ly > ndom(2)) cycle
                if (lz < 1 .or. lz > ndom(3)) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                if (use_micro_A .and. present(Ac_micro)) then
                  A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                else
                  A_local(1:3) = Ac_tot(1:3)
                end if
                phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                overlap_i = overlap_i + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                            ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end if
            end do
            uVpsi(io, ilma, ispin) = overlap_i
          end do
        end do
        !$omp end parallel do

        nbf = dg_frag%n_basis(ifrag, ispin)
        !$omp parallel do collapse(2) private(jo, io, ig_i, ig_j, ilma, nlpp_contrib, overlap_i, overlap_j)
        do jo = 1, nbf
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            nlpp_contrib = (0.0d0, 0.0d0)
            do ilma = 1, ppg%Nlma
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              nlpp_contrib = nlpp_contrib + overlap_i * overlap_j * ppg%rinv_uvu(ilma)
            end do
            H_nl(ig_i, ig_j, ispin) = H_nl(ig_i, ig_j, ispin) + nlpp_contrib
          end do
        end do
        !$omp end parallel do
      end do
    end do

    deallocate(uVpsi)

  end subroutine build_nonlocal_pp_matrix_A

  !=======================================================================
  ! Ensure cached non-local PP matrix for current A(t)
  !=======================================================================
  subroutine ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
    use structures
    use salmon_global, only: ae_shape1, ae_shape2, theory
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(in)    :: Ac_tot(3)

    real(8) :: delta_A
    logical :: reuse_allowed
    logical :: use_micro_A

    if (ppg%Nlma == 0 .or. .not. allocated(ppg%uV)) then
      if (allocated(dg_frag%H_nl_cache)) deallocate(dg_frag%H_nl_cache)
      dg_frag%has_nl_cache = .false.
      return
    end if

    if (.not. allocated(dg_frag%H_nl_cache)) then
      allocate(dg_frag%H_nl_cache(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      dg_frag%has_nl_cache = .false.
    end if

    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))

    reuse_allowed = (.not. use_micro_A) .and. &
                    (trim(ae_shape1) == 'impulse' .or. trim(ae_shape1) == 'none') .and. &
                    (trim(ae_shape2) == 'impulse' .or. trim(ae_shape2) == 'none')

    delta_A = maxval(abs(Ac_tot - dg_frag%Ac_nl_cache))
    if (.not. dg_frag%has_nl_cache .or. (.not. reuse_allowed) .or. delta_A > dg_frag%Ac_nl_cache_tol) then
      if (use_micro_A) then
        call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .true., system%Ac_micro%v, dg_frag%H_nl_cache)
      else
        call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .false., H_nl=dg_frag%H_nl_cache)
      end if
      dg_frag%Ac_nl_cache = Ac_tot
      dg_frag%has_nl_cache = .true.
    end if

  end subroutine ensure_nonlocal_pp_matrix_A  !=======================================================================
  ! Apply gradient operator to a basis function using finite differences
  !=======================================================================
  subroutine apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(out) :: grad_phi(mg%is(1):mg%ie(1), &
                                                    mg%is(2):mg%ie(2), &
                                                    mg%is(3):mg%ie(3), 3)

    integer :: lx, ly, lz, gx, gy, gz, ifrag
    real(8) :: nabt(4,3)
    integer :: is(3), ie(3), iorg(3), ndom(3)

    nabt = stencil%coef_nab
    is = mg%is
    ie = mg%ie
    ifrag = dg_frag%ifrag_start + i_local - 1
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)

    grad_phi = 0.0d0

    !$omp parallel do collapse(3) private(lx, ly, lz, gx, gy, gz) schedule(static)
    do lz = 1, ndom(3)
      do ly = 1, ndom(2)
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          gy = iorg(2) + ly - 1
          gz = iorg(3) + lz - 1
          if (gx < is(1) .or. gx > ie(1)) cycle
          if (gy < is(2) .or. gy > ie(2)) cycle
          if (gz < is(3) .or. gz > ie(3)) cycle

          grad_phi(gx, gy, gz, 1) = &
              nabt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
              nabt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
              nabt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
              nabt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-4, ly, lz, jo, i_local))

          grad_phi(gx, gy, gz, 2) = &
              nabt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
              nabt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
              nabt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
              nabt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-4, lz, jo, i_local))

          grad_phi(gx, gy, gz, 3) = &
              nabt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
              nabt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
              nabt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
              nabt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_gradient_to_basis

  subroutine calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr)
    ! NOTE: DG microscopic current here intentionally excludes non-local PP contribution.
    use structures
    use salmon_global, only: nelec
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_vector),         intent(inout) :: curr

    integer :: ifrag, i_local, ispin, io, istate_frag, jstate_frag
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: nxyz(3), ixyz0(3)
    integer :: nocc_per_spin
    real(8) :: occ_factor
    complex(8) :: coef_pair
    real(8) :: phi_i
    real(8), allocatable :: grad_phi(:,:,:,:)
    real(8), allocatable :: curr_local(:,:,:,:), curr_sum(:,:,:,:)
    real(8), allocatable :: w_local(:,:,:), w_sum(:,:,:)

    curr%v = 0.0d0
    if (.not. allocated(dg_frag%phi_frag)) return

    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    occ_factor = 2.0d0 / real(system%nspin, 8)

    allocate(curr_local(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(curr_sum(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    curr_local = 0.0d0
    curr_sum = 0.0d0
    w_local = 0.0d0
    w_sum = 0.0d0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
            iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
            izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
            if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
            if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
            if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle
            w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
          end do
        end do
      end do

      do ispin = 1, system%nspin
        do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
          ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

          allocate(grad_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), 3))
          call apply_gradient_to_basis(dg_frag, i_local, jstate_frag, mg, stencil, grad_phi)

          do istate_frag = 1, dg_frag%n_basis(ifrag, ispin)
            ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            coef_pair = (0.0d0, 0.0d0)
            do io = 1, nocc_per_spin
              coef_pair = coef_pair + occ_factor * conjg(dg_frag%coef(ig_i, io, ispin)) * dg_frag%coef(ig_j, io, ispin)
            end do
            if (abs(aimag(coef_pair)) < 1.0d-18) cycle

            do iz = 1, nxyz(3)
              do iy = 1, nxyz(2)
                do ix = 1, nxyz(1)
                  ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                  iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                  izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                  if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                  if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                  if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                  phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                  curr_local(1, ixg, iyg, izg) = curr_local(1, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ixg, iyg, izg, 1)
                  curr_local(2, ixg, iyg, izg) = curr_local(2, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ixg, iyg, izg, 2)
                  curr_local(3, ixg, iyg, izg) = curr_local(3, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ixg, iyg, izg, 3)
                end do
              end do
            end do
          end do

          deallocate(grad_phi)
        end do
      end do
    end do

    call comm_summation(curr_local, curr_sum, size(curr_local), dg_frag%icomm)
    call comm_summation(w_local, w_sum, size(w_local), dg_frag%icomm)

    do izg = mg%is(3), mg%ie(3)
      do iyg = mg%is(2), mg%ie(2)
        do ixg = mg%is(1), mg%ie(1)
          if (w_sum(ixg, iyg, izg) > 1.0d-12) then
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg) / w_sum(ixg, iyg, izg)
          else
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg)
          end if
        end do
      end do
    end do

    deallocate(curr_local, curr_sum, w_local, w_sum)

  end subroutine calculate_microscopic_current_dg

  !=======================================================================
  ! Build spatially resolved A·∇ and A^2/2 matrices from Ac_micro(r)
  !=======================================================================
  subroutine build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    integer,                intent(in)    :: ispin
    real(8),                intent(out)   :: Ap_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)
    real(8),                intent(out)   :: A2_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)

    integer :: ifrag, i_local, io, jo
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: nxyz(3), ixyz0(3)
    real(8) :: phi_i, A2val, Ap_int, A2_int
    real(8), allocatable :: grad_phi(:,:,:,:)
    real(8), allocatable :: Ap_flat(:), A2_flat(:), Ap_tmp(:), A2_tmp(:)
    integer :: n

    Ap_mat = 0.0d0
    A2_mat = 0.0d0
    if (.not. allocated(system%Ac_micro%v)) return

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      do jo = 1, dg_frag%n_basis(ifrag, ispin)
        ig_j = dg_frag%index_basis(jo, ifrag, ispin)
        if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

        allocate(grad_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), 3))
        call apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)

        !$omp parallel do private(io, ig_i, Ap_int, A2_int, iz, iy, ix, ixg, iyg, izg, phi_i, A2val) schedule(static)
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig_i = dg_frag%index_basis(io, ifrag, ispin)
          if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

          Ap_int = 0.0d0
          A2_int = 0.0d0
          do iz = 1, nxyz(3)
            do iy = 1, nxyz(2)
              do ix = 1, nxyz(1)
                ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                phi_i = dg_frag%phi_frag(ix, iy, iz, io, i_local)
                Ap_int = Ap_int + phi_i * ( &
                         system%Ac_micro%v(1, ixg, iyg, izg) * grad_phi(ixg, iyg, izg, 1) + &
                         system%Ac_micro%v(2, ixg, iyg, izg) * grad_phi(ixg, iyg, izg, 2) + &
                         system%Ac_micro%v(3, ixg, iyg, izg) * grad_phi(ixg, iyg, izg, 3) ) * system%hvol

                A2val = system%Ac_micro%v(1, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(2, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(3, ixg, iyg, izg)**2
                A2_int = A2_int + 0.5d0 * phi_i * A2val * dg_frag%phi_frag(ix, iy, iz, jo, i_local) * system%hvol
              end do
            end do
          end do

          Ap_mat(ig_i, ig_j) = Ap_mat(ig_i, ig_j) + Ap_int
          A2_mat(ig_i, ig_j) = A2_mat(ig_i, ig_j) + A2_int
        end do
        !$omp end parallel do

        deallocate(grad_phi)
      end do
    end do

    n = dg_frag%n_mat_max * dg_frag%n_mat_max
    allocate(Ap_flat(n), A2_flat(n), Ap_tmp(n), A2_tmp(n))
    Ap_flat = reshape(Ap_mat, [n])
    A2_flat = reshape(A2_mat, [n])
    call comm_summation(Ap_flat, Ap_tmp, n, dg_frag%icomm)
    call comm_summation(A2_flat, A2_tmp, n, dg_frag%icomm)
    Ap_mat = reshape(Ap_tmp, [dg_frag%n_mat_max, dg_frag%n_mat_max])
    A2_mat = reshape(A2_tmp, [dg_frag%n_mat_max, dg_frag%n_mat_max])
    deallocate(Ap_flat, A2_flat, Ap_tmp, A2_tmp)

  end subroutine build_spatial_A_coupling_matrices

end module rt_dg_fragment_ops
