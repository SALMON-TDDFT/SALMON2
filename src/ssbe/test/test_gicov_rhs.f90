! src/ssbe/test/test_gicov_rhs.f90
! LG-SBE gicov Phase 3, Task 5a: single-evaluation property test of the gicov
! RHS operator gicov_rhs (the physics core drho/dt) + the R-1 qnm<->rho bridge.
!
! WHY THIS FILE EXISTS
!   Task 5a builds the gicov right-hand side as a callable routine (gicov_rhs)
!   and the rho_ij_from_q / q_ij_from_rho representation bridge, plus the
!   prepare_qnm gicov change (same-block exp_iphi=1 / abs_dnm=0).  It does NOT
!   build the integrator (Strang/Taylor4 + AB4 conservation gate = Task 5b).
!   Isolating the RHS lets a RHS/representation/sign bug be told apart from an
!   integrator-stability issue, so this test asserts on ONE evaluation of drho
!   (no time stepping) against the four properties a coherent, gauge-covariant
!   generator must satisfy:
!
!     P (trace-preserving): sum_ik sum_n Re(drho(n,n,ik)) = 0.  Protects total
!         population conservation -- the covariant intraband gradient telescopes
!         to zero total-trace over the periodic grid, and both the dipole
!         commutator and the energy/dephasing terms are traceless per k.
!     H (Hermitian): drho(:,:,k) Hermitian for every k.  Protects "drho/dt of a
!         Hermitian rho under a Hermitian generator is Hermitian" -- guards a
!         wrong sign / a non-Hermitian d_matrix / a broken commutator.
!     G (gauge-covariant): under a random per-k BLOCK gauge W(k) applied to
!         u_transport (Wilson-line), p_mod_matrix (-> d_matrix) and rho,
!         drho_gauged(k) = W(k)^H drho_ungauged(k) W(k).  This is THE gicov
!         property (Approach-B'): the representation bridge and every sign must
!         be right or it fails.  Exact-degenerate block => the energy term is
!         block-constant hence covariant; t_2 -> huge => the (strictly legacy,
!         not block-covariant) dephasing is negligible at the 1e-9 gate.
!     N (nonzero): the covariant contribution AND the dipole contribution are
!         each individually nonzero -- guards against a vacuous all-zero RHS
!         that would satisfy P/H/G trivially.
!
!   Fixture: nb=4, nk=8 (num_kgrid=(8,1,1)); a fixed EXACTLY-degenerate composite
!   block {2,3} at eigen 0.90 isolated by a 0.60 gap from singletons band 1
!   (0.30) and band 4 (1.50) -- the nonzero out-of-block gap keeps the dipole
!   active; the k-dependent smooth Hermitian test density + a nontrivial 2x2
!   block Wilson transport keep the covariant term active.  sbe_lg_degen='gicov'
!   and the REAL prepare_qnm is run (exercises the same-block exp_iphi=1 change);
!   the test density is injected into sbe%qnm via q_ij_from_rho so gicov_rhs
!   reconstructs it exactly via rho_ij_from_q.
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
  use degenerate_block_ssbe, only: covariant_grad_block
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
    t_2 = 1d30                 ! dephasing negligible at the 1e-9 gauge gate
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

  !======================= build d_matrix from p (gicov rule) =================
  ! Mirrors prepare_matrix's gicov branch: same-block -> 0; out-of-block with
  ! |delta_omega|>eps -> i*p/delta_omega; else 0.  d is Hermitian iff p is.
  subroutine build_d_from_p(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), parameter :: eps = 1d-12
    integer :: ik, ib, jb
    gs%d_matrix = (0d0, 0d0)
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          if (blk(ib) == blk(jb)) then
            gs%d_matrix(ib, jb, 1:3, ik) = (0d0, 0d0)
          else if (abs(gs%delta_omega(ib, jb, ik)) > eps) then
            gs%d_matrix(ib, jb, 1:3, ik) = zi_ * gs%p_mod_matrix(ib, jb, 1:3, ik) &
                                         & / gs%delta_omega(ib, jb, ik)
          else
            gs%d_matrix(ib, jb, 1:3, ik) = (0d0, 0d0)
          end if
        end do
      end do
    end do
  end subroutine build_d_from_p

  !======================= synthetic gicov gs fixture =========================
  ! Exactly-degenerate composite block {2,3}; singletons {1},{4}; nonzero
  ! out-of-block gap (0.60) => dipole active.  Hermitian p on axis 1 only
  ! (epdir=(1,0,0)); nontrivial 2x2 block Wilson transport on axis 1.
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
    real(8) :: E(3), dk(3), tr, hn, gerr, cov_mag, dip_mag, rt_err, eqerr
    logical :: deph_by_gw
    complex(8) :: rho(nb, nb, nk), rhog(nb, nb, nk)
    complex(8) :: drho(nb, nb, nk), drhog(nb, nb, nk), drho_expected(nb, nb, nk)
    complex(8) :: Dq(nb, nb, 3, nk), dE(nb, nb), cc, W(nb, nb), Wh(nb, nb), tgt(nb, nb)
    integer :: icomm, ik, ib, jb, lb, axis, kk

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

    ! ---- N: covariant AND dipole contributions each individually nonzero ------
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
    dip_mag = 0d0
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          dE(ib, jb) = E(1)*gs%d_matrix(ib,jb,1,ik) + E(2)*gs%d_matrix(ib,jb,2,ik) + E(3)*gs%d_matrix(ib,jb,3,ik)
        end do
      end do
      do ib = 1, nb
        do jb = 1, nb
          cc = (0d0, 0d0)
          do lb = 1, nb
            cc = cc + dE(ib, lb) * rho(lb, jb, ik) - rho(ib, lb, ik) * dE(lb, jb)
          end do
          dip_mag = max(dip_mag, abs(cc))
        end do
      end do
    end do
    call check_true(cov_mag > 1d-6, &
      "N: covariant intraband contribution is nonzero (not vacuous)", nfail)
    call check_true(dip_mag > 1d-6, &
      "N: out-of-block dipole contribution is nonzero (not vacuous)", nfail)

    ! ---- E: direct-assembly equality (drho == independently assembled RHS) ----
    ! P/H/G/N above each pass even if a WHOLE TERM were dropped from gicov_rhs
    ! (every term is separately traceless/Hermitian/gauge-covariant, and N only
    ! recomputes magnitudes from the fixture, not from gicov_rhs's own output).
    ! Here drho_expected is independently re-assembled from the SAME ingredients
    ! gicov_rhs uses -- the SAME covariant_grad_block output Dq (already built
    ! above for the N check), the SAME -i*sum_a E_a[d_matrix_a,rho] commutator,
    ! and the SAME -i*delta_omega*rho / -rho/t_2 off-diagonal terms -- and
    ! compared elementwise against gicov_rhs's actual OUTPUT drho.  Dropping
    ! EITHER the covariant transport term OR the dipole commutator term from
    ! gicov_rhs would change drho but leave drho_expected (computed here, not
    ! from gicov_rhs) unchanged, so this check would then fail.
    deph_by_gw = (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw')
    eqerr = 0d0
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          dE(ib, jb) = E(1)*gs%d_matrix(ib,jb,1,ik) + E(2)*gs%d_matrix(ib,jb,2,ik) + E(3)*gs%d_matrix(ib,jb,3,ik)
        end do
      end do
      do ib = 1, nb
        do jb = 1, nb
          cc = (0d0, 0d0)
          do lb = 1, nb
            cc = cc + dE(ib, lb) * rho(lb, jb, ik) - rho(ib, lb, ik) * dE(lb, jb)
          end do
          drho_expected(ib, jb, ik) = &
            (E(1)*Dq(ib,jb,1,ik) + E(2)*Dq(ib,jb,2,ik) + E(3)*Dq(ib,jb,3,ik)) - zi_ * cc
          if (ib /= jb) then
            drho_expected(ib, jb, ik) = drho_expected(ib, jb, ik) &
              - zi_ * gs%delta_omega(ib, jb, ik) * rho(ib, jb, ik)
            if (.not. deph_by_gw) then
              drho_expected(ib, jb, ik) = drho_expected(ib, jb, ik) - rho(ib, jb, ik) / t_2
            end if
          end if
          eqerr = max(eqerr, abs(drho(ib, jb, ik) - drho_expected(ib, jb, ik)))
        end do
      end do
    end do
    call check_true(is_finite(eqerr) .and. eqerr < 1d-12, &
      "E: drho == direct-assembly gterm - i*cterm + energy + dephasing (non-vacuous equality)", nfail)

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
    write(*, '(a,es10.2,a,es10.2)') &
      "      nonzero    cov=", cov_mag, "  dip=", dip_mag
  end subroutine test_rhs

end program test_gicov_rhs
