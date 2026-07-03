! src/ssbe/test/test_gicov_rhs.f90
! LG-SBE gicov Phase 3, Task 5a: single-evaluation property test of the gicov
! RHS operator gicov_rhs (the physics core drho/dt) + the R-1 qnm<->rho bridge.
!
! UPDATED for X-full (Task 1, plans/2026-07-03-gicov-xfull.md): gicov_rhs no
! longer reads gs%d_matrix at all -- the full-band covariant transport already
! supplies the WHOLE field term (intraband + interband), so the separate
! analytic dipole commutator this test used to exercise is GONE (double-
! counting if it were still added; proven equivalent by test_gicov_xfull.f90's
! Test D). The fixture below now sets d_matrix=0 (matching production: Task 2
! zeroes it in gs_info_ssbe.f90's gicov branch); the former "N: dipole nonzero"
! check is dropped (there is no dipole term left to be non-vacuous) and the
! former "E: direct-assembly" equality no longer includes the dropped -i*cterm
! commutator. P/H/G are unchanged and still pass.
!
! WHY THIS FILE EXISTS
!   Task 5a builds the gicov right-hand side as a callable routine (gicov_rhs)
!   and the rho_ij_from_q / q_ij_from_rho representation bridge, plus the
!   prepare_qnm gicov change (X-full: exp_iphi=1/abs_dnm=0 for ALL off-diagonal
!   pairs).  It does NOT build the integrator (Strang/Taylor4 + AB4
!   conservation gate = Task 5b).  Isolating the RHS lets a RHS/representation/
!   sign bug be told apart from an integrator-stability issue, so this test
!   asserts on ONE evaluation of drho (no time stepping) against the
!   properties a coherent, gauge-covariant generator must satisfy:
!
!     P (trace-preserving): sum_ik sum_n Re(drho(n,n,ik)) = 0.  Protects total
!         population conservation -- the covariant transport telescopes to
!         zero total-trace over the periodic grid, and the energy/dephasing
!         terms are traceless per k.
!     H (Hermitian): drho(:,:,k) Hermitian for every k.  Protects "drho/dt of a
!         Hermitian rho under a Hermitian generator is Hermitian" -- guards a
!         wrong sign or a broken commutator.
!     G (gauge-covariant): under a random per-k BLOCK gauge W(k) applied to
!         u_transport (Wilson-line) and rho, drho_gauged(k) =
!         W(k)^H drho_ungauged(k) W(k).  This is THE gicov property
!         (Approach-B'): the representation bridge and every sign must be
!         right or it fails.  Exact-degenerate block => the energy term is
!         block-constant hence covariant.  Dephasing now runs at a FINITE,
!         physical-scale t_2 (was 1e30 purely to hide it): the Delta-omega gate
!         skips the non-covariant scalar T2 INSIDE the degenerate block {2,3}
!         while keeping it on energy-distinct pairs, so G holds to 1e-9 WITH
!         dephasing active (ungated, the {2,3} coherence breaks G at ~1/t_2).
!     T2G (Delta-omega dephasing gate): at zero field each off-diagonal drho is
!         (coherent energy) + (gated T2) only.  Inside the exactly degenerate
!         block {2,3} both vanish => drho(2,3) = 0 to machine precision (ungated:
!         -rho(2,3)/t_2 ~ 1/t_2); across the energy-distinct pair (1,2) the T2 is
!         ACTIVE and drho(1,2) = (-i*delta_omega - 1/t_2) rho(1,2) exactly.
!     N (nonzero): the covariant contribution is nonzero -- guards against a
!         vacuous all-zero RHS that would satisfy P/H/G trivially.
!
!   Fixture: nb=4, nk=8 (num_kgrid=(8,1,1)); a fixed EXACTLY-degenerate composite
!   block {2,3} at eigen 0.90 isolated by a 0.60 gap from singletons band 1
!   (0.30) and band 4 (1.50) -- the k-dependent smooth Hermitian test density +
!   a nontrivial 2x2 block Wilson transport keep the covariant term active.
!   d_matrix = 0 (X-full). sbe_lg_degen='gicov' and the REAL prepare_qnm is
!   run (exercises the all-off-diagonal exp_iphi=1 change); the test density is
!   injected into sbe%qnm via q_ij_from_rho so gicov_rhs reconstructs it
!   exactly via rho_ij_from_q.
!
! BUILD (already-built ninja tree at build_local/; single-process communication
! dummy).  Links the SAME objects the salmon executable built, minus main.f90.o
! (modelled on test_gifix_propagator.f90):
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/src/ssbe -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_gicov_rhs.f90 -o <scratch_dir>/test_gicov_rhs.o
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o') \
!     <scratch_dir>/test_gicov_rhs.o -o <scratch_dir>/test_gicov_rhs \
!     -framework Accelerate -lm -ldl
!
!   <scratch_dir>/test_gicov_rhs
!
program test_gicov_rhs
  use gs_info_ssbe,          only: s_sbe_gs_info
  use bloch_solver_ssbe,     only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                    prepare_qnm, gicov_rhs, rho_ij_from_q, q_ij_from_rho
  use degenerate_block_ssbe, only: covariant_grad_block, theta_off
  use salmon_global,         only: epdir_re1, am_s, num_kgrid, t_2, sbe_lg_degen, &
                                    sbe_lg_diag, yn_sbe_gw_collision, sbe_deph_mode
  implicit none

  integer, parameter :: nb = 4, nk = 8
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: two_pi = 6.28318530717958647692d0
  ! block partition: band1 -> block 1; bands 2,3 -> block 2; band4 -> block 3
  integer, parameter :: blk(nb) = (/ 1, 2, 2, 3 /)
  integer :: nfail

  nfail = 0
  call set_globals()
  call test_rhs(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_rhs checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 10d0                 ! FINITE dephasing (was 1d30 to hide it): the
                               ! Delta-omega gate must keep G covariant with T2 active
    sbe_lg_diag = 0            ! no diagnostic knockouts
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

  !======================= 2x2 U(2) gauge on the {2,3} block ==================
  ! W = e^{i psi} [[c e^{i chi}, s e^{i xi}], [-s e^{-i xi}, c e^{-i chi}]]
  ! (c=cos phi, s=sin phi) is an exactly-unitary U(2), plus independent phases
  ! on the singleton bands 1 and 4.  Block-diagonal wrt blk() by construction.
  function make_W(ik) result(W)
    implicit none
    integer, intent(in) :: ik
    complex(8) :: W(nb, nb)
    real(8) :: t, phi, chi, xi, psi, a1, a4, c, s
    t   = dble(ik)
    phi = 0.37d0 + 0.13d0 * t
    chi = 0.21d0 + 0.29d0 * t
    xi  = 0.53d0 - 0.17d0 * t
    psi = 0.11d0 + 0.19d0 * t
    a1  = 0.41d0 + 0.23d0 * t     ! singleton band-1 phase
    a4  = 0.61d0 - 0.31d0 * t     ! singleton band-4 phase
    c = cos(phi); s = sin(phi)
    W = (0d0, 0d0)
    W(1, 1) = exp(zi_ * a1)
    W(2, 2) = exp(zi_ * psi) * c * exp( zi_ * chi)
    W(2, 3) = exp(zi_ * psi) * s * exp( zi_ * xi )
    W(3, 2) = exp(zi_ * psi) * (-s) * exp(-zi_ * xi )
    W(3, 3) = exp(zi_ * psi) * c * exp(-zi_ * chi)
    W(4, 4) = exp(zi_ * a4)
  end function make_W

  function hconj(A) result(Ah)
    implicit none
    complex(8), intent(in) :: A(nb, nb)
    complex(8) :: Ah(nb, nb)
    Ah = transpose(conjg(A))
  end function hconj

  ! forward-neighbour k-index along a Cartesian axis for num_kgrid=(nk,1,1)
  integer function knext(ik, axis) result(kn)
    implicit none
    integer, intent(in) :: ik, axis
    if (axis == 1) then
      kn = mod(ik, nk) + 1
    else
      kn = ik                       ! num_kgrid=1 on axes 2,3: neighbour is self
    end if
  end function knext

  !======================= d_matrix (X-full: always zero) =====================
  ! X-full: gicov_rhs no longer reads gs%d_matrix (the full-band covariant
  ! transport supplies the whole field term, incl. the interband dipole) --
  ! zero it, matching the Task 2 gs_info_ssbe.f90 gicov branch's own
  ! d_matrix=0 wiring. Kept as a named subroutine (still called from build_gs
  ! and the G-test's gauged gsg rebuild) so the fixture's shape/call sites are
  ! unchanged; p_mod_matrix is still built (harmless, just unused for d now).
  subroutine build_d_from_p(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    gs%d_matrix = (0d0, 0d0)
  end subroutine build_d_from_p

  !======================= synthetic gicov gs fixture =========================
  ! Exactly-degenerate composite block {2,3}; singletons {1},{4}; 0.60
  ! out-of-block gap (kept for delta_omega/energy-term realism; d_matrix=0
  ! regardless, X-full).  Hermitian p on axis 1 only (epdir=(1,0,0));
  ! nontrivial 2x2 block Wilson transport on axis 1.
  subroutine build_gs(gs)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8) :: eigen(nb), t, ang, c, s, phz
    integer :: ik, ib, jb

    gs%nk = nk; gs%nb = nb; gs%ne = 6
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_mod_matrix(nb, nb, 3, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))
    allocate(gs%u_transport(nb, nb, 3, nk))
    allocate(gs%block_id(nb, nk))
    allocate(gs%bvec(3, 3))

    eigen(1) = 0.30d0
    eigen(2) = 0.90d0          ! block {2,3}: EXACTLY degenerate
    eigen(3) = 0.90d0
    eigen(4) = 1.50d0

    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi
    gs%b_matrix(2, 2) = two_pi
    gs%b_matrix(3, 3) = two_pi

    gs%nbvec = 3
    gs%bvec(:, 1) = (/ 1, 0, 0 /)
    gs%bvec(:, 2) = (/ 0, 1, 0 /)
    gs%bvec(:, 3) = (/ 0, 0, 1 /)

    gs%p_mod_matrix = (0d0, 0d0)
    gs%u_transport  = (0d0, 0d0)

    do ik = 1, nk
      t = dble(ik - 1)
      gs%eigen(:, ik) = eigen(:)
      gs%block_id(:, ik) = blk(:)
      ! equal full occupation on the block {2,3}; band1 full, band4 empty
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0
      gs%occup(3, ik) = 2d0; gs%occup(4, ik) = 0d0
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do

      ! Hermitian momentum on axis 1 (upper triangle set, lower = conjg)
      do ib = 1, nb
        do jb = ib + 1, nb
          phz = 0.3d0 * dble(ib) + 0.5d0 * dble(jb) + 0.2d0 * t
          gs%p_mod_matrix(ib, jb, 1, ik) = (0.4d0 + 0.1d0 * dble(ib + jb)) * exp(zi_ * phz)
          gs%p_mod_matrix(jb, ib, 1, ik) = conjg(gs%p_mod_matrix(ib, jb, 1, ik))
        end do
      end do

      ! block Wilson transport on axis 1: 2x2 U(2) on {2,3}, identity elsewhere
      gs%u_transport(1, 1, 1, ik) = (1d0, 0d0)
      gs%u_transport(4, 4, 1, ik) = (1d0, 0d0)
      ang = 0.30d0 + 0.15d0 * t
      c = cos(ang); s = sin(ang)
      phz = 0.10d0 + 0.05d0 * t
      gs%u_transport(2, 2, 1, ik) = c * exp( zi_ * phz)
      gs%u_transport(2, 3, 1, ik) = s
      gs%u_transport(3, 2, 1, ik) = -s
      gs%u_transport(3, 3, 1, ik) = c * exp(-zi_ * phz)
      ! axes 2,3: identity (num_kgrid=1 there, no transport)
      do ib = 1, nb
        gs%u_transport(ib, ib, 2, ik) = (1d0, 0d0)
        gs%u_transport(ib, ib, 3, ik) = (1d0, 0d0)
      end do
    end do

    call build_d_from_p(gs)
  end subroutine build_gs

  !======================= smooth Hermitian test density ======================
  ! Full nb x nb Hermitian rho with k-dependent diagonal (drives the covariant
  ! gradient) and k-dependent within-block + out-of-block coherences.
  subroutine build_rho(rho)
    implicit none
    complex(8), intent(out) :: rho(nb, nb, nk)
    real(8) :: th, phz
    integer :: ik, ib, jb
    rho = (0d0, 0d0)
    do ik = 1, nk
      th = two_pi * dble(ik - 1) / dble(nk)
      rho(1, 1, ik) = dcmplx(0.80d0 + 0.15d0 * cos(th),        0d0)
      rho(2, 2, ik) = dcmplx(0.60d0 + 0.10d0 * sin(th + 0.5d0), 0d0)
      rho(3, 3, ik) = dcmplx(0.55d0 + 0.10d0 * cos(th + 1.0d0), 0d0)
      rho(4, 4, ik) = dcmplx(0.20d0 + 0.05d0 * sin(th + 2.0d0), 0d0)
      do ib = 1, nb
        do jb = ib + 1, nb
          phz = 0.4d0 * dble(ib) + 0.7d0 * dble(jb) + th
          rho(ib, jb, ik) = (0.10d0 + 0.03d0 * dble(ib + jb)) * exp(zi_ * phz)
          rho(jb, ib, ik) = conjg(rho(ib, jb, ik))
        end do
      end do
    end do
  end subroutine build_rho

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

  !======================= the test ===========================================
  subroutine test_rhs(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs, gsg
    type(s_sbe_bloch_solver) :: sbe, sbeg
    real(8) :: E(3), dk(3), tr, hn, gerr, cov_mag, rt_err, eqerr
    real(8) :: E0(3), gate_in, gate_out_err
    logical :: deph_by_gw
    complex(8) :: rho(nb, nb, nk), rhog(nb, nb, nk)
    complex(8) :: drho(nb, nb, nk), drhog(nb, nb, nk), drho_expected(nb, nb, nk)
    complex(8) :: drho0(nb, nb, nk), tgt12
    complex(8) :: Dq(nb, nb, 3, nk), W(nb, nb), Wh(nb, nb), tgt(nb, nb)
    integer :: icomm, ik, ib, jb, axis, kk

    icomm = 0
    E(1) = 0.1d0; E(2) = 0d0; E(3) = 0d0

    ! ------- ungauged: real prepare_qnm, inject smooth Hermitian rho ----------
    call build_gs(gs)
    call build_rho(rho)
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call set_qnm_from_rho(sbe, rho)

    ! bridge round-trip: rho_ij_from_q o q_ij_from_rho == identity (exercises R-1)
    rt_err = 0d0
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          rt_err = max(rt_err, abs(rho_ij_from_q(sbe, ib, jb, ik) - rho(ib, jb, ik)))
        end do
      end do
    end do
    call check_true(is_finite(rt_err) .and. rt_err < 1d-12, &
      "bridge: rho_ij_from_q(q_ij_from_rho(rho)) == rho (R-1 round-trip)", nfail)

    call gicov_rhs(sbe, gs, E, drho, icomm)

    ! ---- P: trace-preserving (total population conserved by coherent RHS) -----
    tr = 0d0
    do ik = 1, nk
      do ib = 1, nb
        tr = tr + dble(drho(ib, ib, ik))
      end do
    end do
    call check_true(is_finite(tr) .and. abs(tr) < 1d-10, &
      "P: sum_ik sum_n Re(drho(n,n,ik)) = 0 to 1e-10 (trace-preserving)", nfail)

    ! ---- H: Hermitian drho for every k ---------------------------------------
    hn = 0d0
    do ik = 1, nk
      do jb = 1, nb
        do ib = 1, nb
          hn = max(hn, abs(drho(ib, jb, ik) - conjg(drho(jb, ib, ik))))
        end do
      end do
    end do
    call check_true(is_finite(hn) .and. hn < 1d-10, &
      "H: drho(:,:,k) Hermitian for every k to 1e-10", nfail)

    ! ---- N: covariant contribution nonzero (X-full: no separate dipole term) --
    do axis = 1, 3
      dk(axis) = gs%b_matrix(axis, axis) / dble(num_kgrid(axis))
    end do
    call covariant_grad_block(nb, nk, gs%nbvec, gs%bvec, num_kgrid, &
                              gs%u_transport, rho, dk, Dq)
    cov_mag = 0d0
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          cov_mag = max(cov_mag, abs(E(1)*Dq(ib,jb,1,ik) + E(2)*Dq(ib,jb,2,ik) + E(3)*Dq(ib,jb,3,ik)))
        end do
      end do
    end do
    call check_true(cov_mag > 1d-6, &
      "N: covariant transport contribution is nonzero (not vacuous)", nfail)
    call check_true(maxval(abs(gs%d_matrix)) < 1d-300, &
      "N: d_matrix is exactly 0 (X-full: gicov_rhs no longer reads it)", nfail)

    ! ---- E: direct-assembly equality (drho == independently assembled RHS) ----
    ! P/H/G/N above each pass even if the covariant term were dropped from
    ! gicov_rhs (every term is separately traceless/Hermitian/gauge-covariant,
    ! and N only recomputes magnitudes from the fixture, not from gicov_rhs's
    ! own output).  Here drho_expected is independently re-assembled from the
    ! SAME ingredients gicov_rhs uses -- the SAME covariant_grad_block output Dq
    ! (already built above for the N check) and the SAME -i*delta_omega*rho /
    ! -rho/t_2 off-diagonal terms -- and compared elementwise against
    ! gicov_rhs's actual OUTPUT drho.  X-full: NO dipole commutator term (it
    ! was dropped from gicov_rhs -- see file header); dropping the covariant
    ! transport term would change drho but leave drho_expected (computed here,
    ! not from gicov_rhs) unchanged, so this check would then fail.
    deph_by_gw = (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw')
    eqerr = 0d0
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          drho_expected(ib, jb, ik) = &
            (E(1)*Dq(ib,jb,1,ik) + E(2)*Dq(ib,jb,2,ik) + E(3)*Dq(ib,jb,3,ik))
          if (ib /= jb) then
            drho_expected(ib, jb, ik) = drho_expected(ib, jb, ik) &
              - zi_ * gs%delta_omega(ib, jb, ik) * rho(ib, jb, ik)
            ! mirror production: Delta-omega gate on the scalar T2 dephasing
            if (.not. deph_by_gw .and. abs(gs%delta_omega(ib, jb, ik)) > theta_off) then
              drho_expected(ib, jb, ik) = drho_expected(ib, jb, ik) - rho(ib, jb, ik) / t_2
            end if
          end if
          eqerr = max(eqerr, abs(drho(ib, jb, ik) - drho_expected(ib, jb, ik)))
        end do
      end do
    end do
    call check_true(is_finite(eqerr) .and. eqerr < 1d-12, &
      "E: drho == direct-assembly gterm + energy + dephasing (X-full, non-vacuous equality)", nfail)

    ! ---- T2G: Delta-omega gate on the phenomenological dephasing --------------
    ! Zero field => covariant transport vanishes, so each off-diagonal drho is
    ! (coherent energy) + (gated T2) only.  Degenerate block {2,3}: delta_omega=0
    ! AND T2 gated off => drho(2,3) = 0 to machine precision (ungated the code
    ! would leave -rho(2,3)/t_2 ~ 1/t_2 here).  Energy-distinct pair (1,2),
    ! |delta_omega| = 0.60 > theta_off: T2 ACTIVE, drho(1,2) = (-i*dw - 1/t_2)*rho.
    E0 = 0d0
    call gicov_rhs(sbe, gs, E0, drho0, icomm)
    gate_in = 0d0; gate_out_err = 0d0
    do ik = 1, nk
      gate_in = max(gate_in, abs(drho0(2, 3, ik)), abs(drho0(3, 2, ik)))
      tgt12 = (- zi_ * gs%delta_omega(1, 2, ik) - 1d0 / t_2) * rho(1, 2, ik)
      gate_out_err = max(gate_out_err, abs(drho0(1, 2, ik) - tgt12))
    end do
    call check_true(is_finite(gate_in) .and. gate_in < 1d-12, &
      "T2G-in: degenerate {2,3} coherence has NO scalar T2 (gated off) at E=0", nfail)
    call check_true(is_finite(gate_out_err) .and. gate_out_err < 1d-12, &
      "T2G-out: energy-distinct (1,2) keeps T2: drho = (-i*dw - 1/t_2)*rho", nfail)

    ! ---- G: gauge covariance under a random per-k block gauge W(k) ------------
    ! gauged rho, u_transport (Wilson line), p (-> d_matrix); eigen/block_id kept.
    call build_gs(gsg)
    do ik = 1, nk
      W  = make_W(ik)
      Wh = hconj(W)
      rhog(:, :, ik) = matmul(matmul(Wh, rho(:, :, ik)), W)
      ! p -> W^H p W (per axis); d rebuilt from gauged p below
      do axis = 1, 3
        gsg%p_mod_matrix(:, :, axis, ik) = matmul(matmul(Wh, gs%p_mod_matrix(:, :, axis, ik)), W)
      end do
    end do
    ! u_transport(k,axis) -> W(k)^H U(k,axis) W(k+e_axis)   (Wilson-line transform)
    do ik = 1, nk
      Wh = hconj(make_W(ik))
      do axis = 1, 3
        kk = knext(ik, axis)
        gsg%u_transport(:, :, axis, ik) = matmul(matmul(Wh, gs%u_transport(:, :, axis, ik)), make_W(kk))
      end do
    end do
    call build_d_from_p(gsg)

    call init_sbe_bloch_solver(sbeg, gsg, nb, icomm)
    call prepare_qnm(sbeg, gsg, icomm)
    call set_qnm_from_rho(sbeg, rhog)
    call gicov_rhs(sbeg, gsg, E, drhog, icomm)

    gerr = 0d0
    do ik = 1, nk
      W  = make_W(ik)
      Wh = hconj(W)
      tgt = matmul(matmul(Wh, drho(:, :, ik)), W)     ! W^H drho_ungauged W
      do ib = 1, nb
        do jb = 1, nb
          gerr = max(gerr, abs(drhog(ib, jb, ik) - tgt(ib, jb)))
        end do
      end do
    end do
    call check_true(is_finite(gerr) .and. gerr < 1d-9, &
      "G: drho_gauged(k) = W(k)^H drho_ungauged(k) W(k) to 1e-9 (gauge covariant)", nfail)

    write(*, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2,a,es10.2)') &
      "      residuals  P=", abs(tr), "  H=", hn, "  G=", gerr, "  E=", eqerr, "  round-trip=", rt_err
    write(*, '(a,es10.2)') &
      "      nonzero    cov=", cov_mag
    write(*, '(a,es10.2,a,es10.2)') &
      "      T2 gate    in(deg {2,3})=", gate_in, "  out(distinct 1,2)=", gate_out_err
  end subroutine test_rhs

end program test_gicov_rhs
