! src/ssbe/test/test_gicov_hex.f90
! Unit test for the gicov_rhs field-term DUAL PATH (non-orthogonal b_matrix
! support, design: gw/gw_design/specs/2026-07-05-gicov-nonorthogonal-bmatrix.md).
!
! The field term of gicov_rhs projects E . grad_k onto the reduced k-axes.
! Legacy (strictly diagonal b_matrix) path: dk = b_ii/N_i, weights = Efield --
! kept VERBATIM (bit-for-bit for all current campaigns; guarded by the
! unchanged existing suite, not by this file).  General (non-diagonal) path:
! k = sum_i s_i b_i and a_i . b_j = 2 pi delta_ij give
!     E . grad_k rho = sum_i c_i drho/ds_i,  c_i = (E . a_i)/(2 pi),
! with the stencil in the REDUCED coordinate, dk_i = 1/N_i.  This test proves:
!
!   O  (orthogonal equivalence): on a strictly diagonal-b fixture (cubic,
!      b = 2 pi I, a = I) the general path FORCED via the module test hook
!      set_gicov_force_general_field agrees with the legacy path to 1e-12
!      (algebraically identical, ULP-level float difference only).
!   OT (teeth for O): the forced general path really executes and really reads
!      gs%a_matrix: rescaling a_matrix -> 3*a_matrix must scale the FIELD term
!      (drho(E) - drho(0), isolating it from energy/dephasing) by exactly 3.
!      If the hook were dead (O trivially true), OT would fail.
!   HX (hexagonal analytic): 60-degree hexagonal 2D lattice (a1 = (1,0,0),
!      a2 = (1/2, sqrt3/2, 0), b_matrix NON-diagonal => general path
!      auto-selected, no forcing), U == identity, delta_omega == 0 (so the RHS
!      is the pure field term).  rho is a known smooth lattice-periodic
!      superposition of plane waves e^{i R.k}, R = m1 a1 + m2 a2 (|m| <= 1 --
!      resolved by the full 4-shell stencil at N=12), for which the analytic
!      directional derivative is E . grad_k rho = i (E . R) rho_mode, computed
!      here PURELY in Cartesian form (E . R with R = m1 a1 + m2 a2; no
!      b_matrix, no 1/2pi -- independent of the c_i formula under test).
!      Checked to stencil accuracy (measured ~1e-5 rel at N=12, gated 1e-4)
!      for three fields: E || x (which on a hexagonal cell must activate BOTH
!      in-plane reduced axes -- a wrongly-skipped axis fails at O(1)), a
!      generic in-plane E, and an E with a z-component (singleton third axis:
!      Dq_3 == 0 must match the analytic zero -- no corruption).
!   P/H (hex sanity): the hex field term is trace-preserving over the periodic
!      grid (stencil telescopes) and Hermitian (real c_i, Hermitian rho).
!
! Both fixtures share nb=4, num_kgrid=(12,12,1), nk=144: the gicov halo plan
! (module state, built once per program) is size-compatible across the parts
! (serial run: no partners, all k local).
!
! BUILD (standalone; same pattern as test_gicov_rhs.f90 -- links the salmon
! build tree's objects minus main.f90.o):
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_gicov_hex.f90 -o <scratch_dir>/test_gicov_hex.o
!   gfortran <flags> $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' \
!     ! -name 'main.f90.o') <scratch_dir>/test_gicov_hex.o -o <scratch_dir>/test_gicov_hex \
!     -framework Accelerate -lm -ldl
!   <scratch_dir>/test_gicov_hex
!
program test_gicov_hex
  use gs_info_ssbe,      only: s_sbe_gs_info
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                prepare_qnm, gicov_rhs, q_ij_from_rho, &
                                set_gicov_force_general_field
  use salmon_global,     only: epdir_re1, am_s, num_kgrid, t_2, sbe_lg_degen, &
                                sbe_lg_diag, yn_sbe_gw_collision, sbe_deph_mode
  implicit none

  integer, parameter :: nb = 4, n1 = 12, n2 = 12, n3 = 1, nk = n1 * n2 * n3
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: two_pi = 6.28318530717958647692d0
  ! HX plane-wave modes (m1, m2): |m| <= 1 so the 4-shell stencil at N=12 is
  ! deep in its convergence regime (per-axis rel. error ~ theta^8/630 ~ 9e-6).
  integer, parameter :: nmode = 3
  integer, parameter :: mvec(2, nmode) = reshape((/ 1, 0,  0, 1,  1, -1 /), (/2, nmode/))
  integer :: nfail

  nfail = 0
  call set_globals()
  call test_orthogonal_equivalence(nfail)
  call test_hexagonal_analytic(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_hex checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4
    num_kgrid(1) = n1; num_kgrid(2) = n2; num_kgrid(3) = n3
    t_2 = 10d0
    sbe_lg_diag = 0
    sbe_lg_degen = 'gicov'
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
  end subroutine set_globals

  !======================= assert helper (ssbe style) =========================
  subroutine check_true(cond, label, nfail)
    implicit none
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (cond) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a)') "FAIL  " // label
      nfail = nfail + 1
    end if
  end subroutine check_true

  logical function is_finite(x) result(ok)
    implicit none
    real(8), intent(in) :: x
    ok = (x == x) .and. (abs(x) <= huge(1d0))
  end function is_finite

  ! reduced coordinates (s1, s2) of grid point ik; ordering matches
  ! build_ik_neighbor: ik = i3*N2*N1 + i2*N1 + i1 + 1 (i1 fastest, 0-based).
  subroutine s_of_ik(ik, s1, s2)
    implicit none
    integer, intent(in) :: ik
    real(8), intent(out) :: s1, s2
    integer :: i1, i2
    i1 = mod(ik - 1, n1)
    i2 = mod((ik - 1) / n1, n2)
    s1 = dble(i1) / dble(n1)
    s2 = dble(i2) / dble(n2)
  end subroutine s_of_ik

  ! deterministic generic complex mode amplitudes A_mu(ib, jb)
  function amp(mu, ib, jb) result(a)
    implicit none
    integer, intent(in) :: mu, ib, jb
    complex(8) :: a
    a = (0.05d0 + 0.02d0 * dble(ib) + 0.013d0 * dble(jb)) &
      & * exp(zi_ * (0.3d0 * dble(ib) + 0.7d0 * dble(jb) + 1.1d0 * dble(mu)))
  end function amp

  !======================= shared gs allocation ===============================
  subroutine alloc_gs(gs)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    gs%nk = nk; gs%nb = nb; gs%ne = 6
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_mod_matrix(nb, nb, 3, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))
    allocate(gs%u_transport(nb, nb, 3, nk))
    allocate(gs%block_id(nb, nk))
    allocate(gs%bvec(3, 3))
    gs%nbvec = 3
    gs%bvec(:, 1) = (/ 1, 0, 0 /)
    gs%bvec(:, 2) = (/ 0, 1, 0 /)
    gs%bvec(:, 3) = (/ 0, 0, 1 /)
    gs%p_mod_matrix = (0d0, 0d0)
    gs%d_matrix     = (0d0, 0d0)
  end subroutine alloc_gs

  ! inject physical rho into sbe%qnm via the q_ij_from_rho bridge
  subroutine set_qnm_from_rho(sbe, rho)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    complex(8), intent(in) :: rho(nb, nb, nk)
    integer :: ik, ib, jb
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, nb
        do jb = 1, nb
          sbe%qnm(ib, jb, ik) = q_ij_from_rho(sbe, rho(ib, jb, ik), ib, jb, ik)
        end do
      end do
    end do
  end subroutine set_qnm_from_rho

  !===========================================================================
  ! Part O/OT: orthogonal (cubic b = 2 pi I) fixture with a nontrivial block
  ! Wilson transport on BOTH in-plane axes and a degenerate block {2,3} --
  ! the general path is FORCED via set_gicov_force_general_field and must agree
  ! with the legacy path to 1e-12.
  !===========================================================================
  subroutine build_gs_diag(gs)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8) :: eigen(nb), t, ang, c, s, phz
    integer :: ik, ib, jb, axis

    call alloc_gs(gs)

    eigen(1) = 0.30d0
    eigen(2) = 0.90d0          ! block {2,3}: exactly degenerate
    eigen(3) = 0.90d0
    eigen(4) = 1.50d0

    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi
    gs%b_matrix(2, 2) = two_pi
    gs%b_matrix(3, 3) = two_pi
    gs%a_matrix = 0d0          ! a_i . b_j = 2 pi delta_ij  (cubic, a = 1)
    gs%a_matrix(1, 1) = 1d0
    gs%a_matrix(2, 2) = 1d0
    gs%a_matrix(3, 3) = 1d0

    gs%u_transport = (0d0, 0d0)
    do ik = 1, nk
      t = dble(ik - 1)
      gs%eigen(:, ik) = eigen(:)
      gs%block_id(:, ik) = (/ 1, 2, 2, 3 /)
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0
      gs%occup(3, ik) = 2d0; gs%occup(4, ik) = 0d0
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do
      ! block Wilson transport: 2x2 U(2) on {2,3} on axes 1 and 2 (both have
      ! live stencils on the 12x12 grid), identity on the singleton axis 3
      do axis = 1, 2
        gs%u_transport(1, 1, axis, ik) = (1d0, 0d0)
        gs%u_transport(4, 4, axis, ik) = (1d0, 0d0)
        ang = 0.30d0 + 0.15d0 * t + 0.4d0 * dble(axis)
        c = cos(ang); s = sin(ang)
        phz = 0.10d0 + 0.05d0 * t - 0.3d0 * dble(axis)
        gs%u_transport(2, 2, axis, ik) = c * exp( zi_ * phz)
        gs%u_transport(2, 3, axis, ik) = s
        gs%u_transport(3, 2, axis, ik) = -s
        gs%u_transport(3, 3, axis, ik) = c * exp(-zi_ * phz)
      end do
      do ib = 1, nb
        gs%u_transport(ib, ib, 3, ik) = (1d0, 0d0)
      end do
    end do
  end subroutine build_gs_diag

  ! smooth Hermitian test density on the 2D reduced grid
  subroutine build_rho_diag(rho)
    implicit none
    complex(8), intent(out) :: rho(nb, nb, nk)
    real(8) :: s1, s2, th1, th2, phz
    integer :: ik, ib, jb
    rho = (0d0, 0d0)
    do ik = 1, nk
      call s_of_ik(ik, s1, s2)
      th1 = two_pi * s1;  th2 = two_pi * s2
      rho(1, 1, ik) = dcmplx(0.80d0 + 0.15d0 * cos(th1) + 0.05d0 * sin(th2), 0d0)
      rho(2, 2, ik) = dcmplx(0.60d0 + 0.10d0 * sin(th1 + 0.5d0) + 0.04d0 * cos(th2 + 0.3d0), 0d0)
      rho(3, 3, ik) = dcmplx(0.55d0 + 0.10d0 * cos(th1 + 1.0d0) - 0.06d0 * sin(th2 + 0.7d0), 0d0)
      rho(4, 4, ik) = dcmplx(0.20d0 + 0.05d0 * sin(th1 + 2.0d0) + 0.03d0 * cos(th2 + 1.1d0), 0d0)
      do ib = 1, nb
        do jb = ib + 1, nb
          phz = 0.4d0 * dble(ib) + 0.7d0 * dble(jb) + th1 + 0.6d0 * th2
          rho(ib, jb, ik) = (0.10d0 + 0.03d0 * dble(ib + jb)) * exp(zi_ * phz)
          rho(jb, ib, ik) = conjg(rho(ib, jb, ik))
        end do
      end do
    end do
  end subroutine build_rho_diag

  subroutine test_orthogonal_equivalence(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    complex(8) :: rho(nb, nb, nk)
    complex(8) :: drho_leg(nb, nb, nk), drho_gen(nb, nb, nk)
    complex(8) :: drho0(nb, nb, nk), drho3(nb, nb, nk)
    real(8) :: E(3), E0(3), eqerr, mag, scerr, fmag
    complex(8) :: f_leg, f_3
    integer :: icomm, ik, ib, jb

    icomm = 0
    E  = (/ 0.11d0, 0.07d0, 0d0 /)   ! both in-plane axes field-active
    E0 = 0d0

    call build_gs_diag(gs)
    call build_rho_diag(rho)
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call set_qnm_from_rho(sbe, rho)

    ! ---- O: legacy vs FORCED general path on the same diagonal-b fixture ----
    call set_gicov_force_general_field(.false.)
    call gicov_rhs(sbe, gs, E, drho_leg, icomm)
    call set_gicov_force_general_field(.true.)
    call gicov_rhs(sbe, gs, E, drho_gen, icomm)

    eqerr = 0d0;  mag = 0d0
    do ik = 1, nk
      do jb = 1, nb
        do ib = 1, nb
          eqerr = max(eqerr, abs(drho_gen(ib, jb, ik) - drho_leg(ib, jb, ik)))
          mag   = max(mag,   abs(drho_leg(ib, jb, ik)))
        end do
      end do
    end do
    call check_true(is_finite(eqerr) .and. eqerr < 1d-12, &
      "O: forced general path == legacy path on diagonal b_matrix (1e-12)", nfail)
    call check_true(mag > 1d-6, &
      "O: fixture RHS is nonzero (equivalence check not vacuous)", nfail)

    ! ---- OT (teeth): the forced general path really reads gs%a_matrix -------
    ! field term = drho(E) - drho(0); scaling a_matrix by 3 scales every c_i
    ! by 3, so the general-path field term must scale by exactly 3.
    call set_gicov_force_general_field(.false.)
    call gicov_rhs(sbe, gs, E0, drho0, icomm)         ! energy+dephasing only
    call set_gicov_force_general_field(.true.)
    gs%a_matrix = 3d0 * gs%a_matrix
    call gicov_rhs(sbe, gs, E, drho3, icomm)
    gs%a_matrix = gs%a_matrix / 3d0
    call set_gicov_force_general_field(.false.)

    scerr = 0d0;  fmag = 0d0
    do ik = 1, nk
      do jb = 1, nb
        do ib = 1, nb
          f_leg = drho_leg(ib, jb, ik) - drho0(ib, jb, ik)
          f_3   = drho3(ib, jb, ik)   - drho0(ib, jb, ik)
          scerr = max(scerr, abs(f_3 - 3d0 * f_leg))
          fmag  = max(fmag,  abs(f_leg))
        end do
      end do
    end do
    call check_true(fmag > 1d-6, &
      "OT: field term is nonzero (teeth check not vacuous)", nfail)
    call check_true(is_finite(scerr) .and. scerr < 1d-12, &
      "OT: a_matrix -> 3a scales the forced-general field term by exactly 3 (hook + a_matrix really used)", nfail)

    write(*, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') &
      "      O residuals  eq=", eqerr, "  |drho|max=", mag, &
      "  OT scale=", scerr, "  |field|max=", fmag
  end subroutine test_orthogonal_equivalence

  !===========================================================================
  ! Part HX: 60-degree hexagonal 2D lattice, U == I, delta_omega == 0.
  ! drho from gicov_rhs (general path AUTO-selected off the non-diagonal
  ! b_matrix) must equal the ANALYTIC directional derivative
  !     E . grad_k rho = sum_mode i (E . R_mode) * (mode term),
  ! R_mode = m1 a1 + m2 a2 (pure Cartesian; independent of the c_i formula).
  !===========================================================================
  subroutine build_gs_hex(gs, aerr)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8), intent(out) :: aerr        ! max |a_i . b_j - 2 pi delta_ij|
    real(8) :: a(3, 3), b(3, 3), cx(3), vol, dev
    integer :: ik, ib, i, j

    call alloc_gs(gs)

    ! 60-degree hexagonal, lattice constant 1 (c = 1)
    a(:, 1) = (/ 1.0d0,  0.0d0,               0.0d0 /)
    a(:, 2) = (/ 0.5d0,  0.5d0 * sqrt(3.0d0), 0.0d0 /)
    a(:, 3) = (/ 0.0d0,  0.0d0,               1.0d0 /)

    ! reciprocal rows b_i = 2 pi (a_j x a_k)/V, same construction as
    ! calc_lattice_info (gs_info_ssbe.f90)
    cx = cross(a(:, 2), a(:, 3));  vol = dot_product(cx, a(:, 1))
    b(1, :) = two_pi * cx / vol
    b(2, :) = two_pi * cross(a(:, 3), a(:, 1)) / vol
    b(3, :) = two_pi * cross(a(:, 1), a(:, 2)) / vol

    gs%a_matrix = a
    gs%b_matrix = b

    ! fixture self-consistency: a_i . b_j = 2 pi delta_ij
    aerr = 0d0
    do i = 1, 3
      do j = 1, 3
        dev = dot_product(gs%a_matrix(:, i), gs%b_matrix(j, :))
        if (i == j) dev = dev - two_pi
        aerr = max(aerr, abs(dev))
      end do
    end do

    gs%eigen = 1.0d0                    ! all-degenerate: delta_omega == 0,
    gs%delta_omega = 0d0                ! T2 gated off => RHS = pure field term
    gs%block_id = 1
    gs%occup = 0d0
    gs%u_transport = (0d0, 0d0)         ! bare stencil: U == identity
    do ik = 1, nk
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0; gs%occup(3, ik) = 2d0
      do ib = 1, nb
        gs%u_transport(ib, ib, 1, ik) = (1d0, 0d0)
        gs%u_transport(ib, ib, 2, ik) = (1d0, 0d0)
        gs%u_transport(ib, ib, 3, ik) = (1d0, 0d0)
      end do
    end do
  end subroutine build_gs_hex

  function cross(u, v) result(w)
    implicit none
    real(8), intent(in) :: u(3), v(3)
    real(8) :: w(3)
    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)
  end function cross

  ! smooth lattice-periodic Hermitian density: superposition of plane waves
  !   rho(ib,jb,k) = sum_mu [ A_mu(ib,jb) f_mu + conjg(A_mu(jb,ib) f_mu) ],
  !   f_mu(k) = e^{2 pi i (m1 s1 + m2 s2)} = e^{i R_mu . k},  R_mu = m1 a1 + m2 a2
  subroutine build_rho_hex(rho)
    implicit none
    complex(8), intent(out) :: rho(nb, nb, nk)
    complex(8) :: f
    real(8) :: s1, s2
    integer :: ik, ib, jb, mu
    rho = (0d0, 0d0)
    do ik = 1, nk
      call s_of_ik(ik, s1, s2)
      do mu = 1, nmode
        f = exp(zi_ * two_pi * (dble(mvec(1, mu)) * s1 + dble(mvec(2, mu)) * s2))
        do jb = 1, nb
          do ib = 1, nb
            rho(ib, jb, ik) = rho(ib, jb, ik) &
              & + amp(mu, ib, jb) * f + conjg(amp(mu, jb, ib) * f)
          end do
        end do
      end do
    end do
  end subroutine build_rho_hex

  ! analytic E . grad_k rho at grid point ik (Cartesian: no b_matrix, no 2 pi)
  subroutine analytic_dirderiv(gs, E, ik, dref)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: E(3)
    integer, intent(in) :: ik
    complex(8), intent(out) :: dref(nb, nb)
    complex(8) :: f, g
    real(8) :: s1, s2, ER, R(3)
    integer :: ib, jb, mu
    dref = (0d0, 0d0)
    call s_of_ik(ik, s1, s2)
    do mu = 1, nmode
      R = dble(mvec(1, mu)) * gs%a_matrix(:, 1) + dble(mvec(2, mu)) * gs%a_matrix(:, 2)
      ER = E(1) * R(1) + E(2) * R(2) + E(3) * R(3)
      f = exp(zi_ * two_pi * (dble(mvec(1, mu)) * s1 + dble(mvec(2, mu)) * s2))
      g = zi_ * ER * f
      do jb = 1, nb
        do ib = 1, nb
          dref(ib, jb) = dref(ib, jb) + amp(mu, ib, jb) * g + conjg(amp(mu, jb, ib) * g)
        end do
      end do
    end do
  end subroutine analytic_dirderiv

  subroutine test_hexagonal_analytic(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    complex(8) :: rho(nb, nb, nk), drho(nb, nb, nk), dref(nb, nb)
    real(8) :: Ecase(3, 3), E(3), aerr, err, ref, tr, hn, relerr(3)
    integer :: icomm, ic, ik, ib, jb
    character(64) :: tag

    icomm = 0
    Ecase(:, 1) = (/ 0.10d0, 0.00d0, 0.00d0 /)  ! E || x: BOTH reduced axes active
    Ecase(:, 2) = (/ 0.08d0, 0.05d0, 0.00d0 /)  ! generic in-plane
    Ecase(:, 3) = (/ 0.08d0, 0.05d0, 0.03d0 /)  ! + z-component (singleton axis)

    call build_gs_hex(gs, aerr)
    call check_true(aerr < 1d-13, &
      "HX fixture: a_i . b_j == 2 pi delta_ij (hexagonal duality)", nfail)
    call check_true(abs(gs%b_matrix(1, 2)) > 1d0, &
      "HX fixture: b_matrix is genuinely non-diagonal (general path auto-selected)", nfail)

    call build_rho_hex(rho)
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call set_qnm_from_rho(sbe, rho)
    call set_gicov_force_general_field(.false.)   ! auto-selection, not forcing

    do ic = 1, 3
      E = Ecase(:, ic)
      call gicov_rhs(sbe, gs, E, drho, icomm)

      err = 0d0;  ref = 0d0
      tr = 0d0;   hn = 0d0
      do ik = 1, nk
        call analytic_dirderiv(gs, E, ik, dref)
        do jb = 1, nb
          do ib = 1, nb
            err = max(err, abs(drho(ib, jb, ik) - dref(ib, jb)))
            ref = max(ref, abs(dref(ib, jb)))
            hn  = max(hn, abs(drho(ib, jb, ik) - conjg(drho(jb, ib, ik))))
          end do
          tr = tr + dble(drho(jb, jb, ik))
        end do
      end do
      relerr(ic) = err / ref

      write(tag, '(a,i0,a)') "HX case ", ic, ""
      call check_true(ref > 1d-3, trim(tag) // ": analytic reference nonzero (not vacuous)", nfail)
      call check_true(is_finite(relerr(ic)) .and. relerr(ic) < 1d-4, &
        trim(tag) // ": E.grad_k rho matches analytic directional derivative (rel < 1e-4, stencil accuracy)", nfail)
      call check_true(is_finite(tr) .and. abs(tr) < 1d-10, &
        trim(tag) // ": field term trace-preserving over the periodic grid", nfail)
      call check_true(is_finite(hn) .and. hn < 1d-10, &
        trim(tag) // ": drho Hermitian", nfail)
    end do

    write(*, '(a,es10.2,a,es10.2,a,es10.2)') &
      "      HX rel errors  E||x=", relerr(1), "  in-plane=", relerr(2), "  +z=", relerr(3)
  end subroutine test_hexagonal_analytic

end program test_gicov_hex
