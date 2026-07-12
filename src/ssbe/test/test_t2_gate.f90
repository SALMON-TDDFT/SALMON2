! src/ssbe/test/test_t2_gate.f90
! T2 Delta-omega dephasing gate: runtime-configurable shape (hard 'step' /
! smooth Gaussian-notch 'gauss') + its own namelist parameters.
! Design: DFT/gw/gw_design/plans/2026-07-08-t2-gate-shape.md
!
! WHY THIS FILE EXISTS
!   The gicov length-gauge SBE's phenomenological T2 dephasing term was
!   hard-gated by the SHARED block-partition threshold theta_off (see
!   bloch_solver_ssbe.f90's gicov_rhs, degenerate_block_ssbe.f90).  This adds
!   THREE new &sbe keys (sbe_t2_gate_shape/theta/width) read ONLY by the T2
!   gate -- theta_off/theta_on/gap_margin (block partition / gi-blend /
!   VG-guard) are untouched.  Task 1 proves the plumbing defaults reproduce
!   the current behaviour bit-for-bit (shape='step', theta=2d-3 == the old
!   theta_off literal, width=0).  Task 2 unit-tests the pure weight helper
!   t2_gate_weight for both shapes + the exact-degeneracy floor clamp.
!   Task 3 unit-tests the checker predicate t2_gate_params_ok.  Task 4 proves
!   (a) shape='step' reproduces gicov_rhs's CURRENT hard-gate output bit-for-
!   bit (regression) and (b) shape='gauss' leaves an exactly-degenerate pair
!   undephased and keeps the dephasing increment Hermitian.
!
! BUILD (mirrors test_gicov_rhs.f90 -- already-built ninja tree at
! build_local/; single-process communication dummy).  Links the SAME objects
! the salmon executable built, minus main.f90.o:
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/src/ssbe -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_t2_gate.f90 -o <scratch_dir>/test_t2_gate.o
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o') \
!     <scratch_dir>/test_t2_gate.o -o <scratch_dir>/test_t2_gate \
!     -framework Accelerate -lm -ldl
!
!   cd <scratch_dir> && <scratch_dir>/test_t2_gate
!     (must run from a writable cwd: test_defaults creates/reads .namelist.tmp
!     via the REAL read_input_common, so it genuinely exercises the &sbe
!     defaults block in inputoutput.f90, not a hand-set fixture.)
!
program test_t2_gate
  use salmon_global, only: sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_t2_gate_width
  use inputoutput,   only: read_input_common
  use gs_info_ssbe,  only: t2_gate_weight, sbe_gs_set_replicated_kmap
  use input_checker_sbe, only: t2_gate_params_ok
  implicit none

  ! Task 4 gicov_rhs-integration fixture constants, shared (host association)
  ! by test_gicov_integration and its SIBLING fixture builders below
  ! (build_gicov_gs/build_gicov_rho/set_qnm_from_rho) -- Fortran does not
  ! allow an internal procedure to itself have a CONTAINS block, so the
  ! fixture builders live at the program level, not nested inside
  ! test_gicov_integration.
  integer, parameter :: nb = 4, nk = 8
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: two_pi = 6.28318530717958647692d0
  integer, parameter :: blk(nb) = (/ 1, 2, 2, 3 /)   ! band1 | {2,3} | band4

  integer :: nfail
  nfail = 0

  call test_defaults(nfail)
  call test_gate_weight(nfail)
  call test_gate_params_ok(nfail)
  call test_gicov_integration(nfail)
  call test_step_floor_clamp(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_t2_gate checks passed."
  end if

contains

  !======================= Task 1: namelist defaults ===========================
  ! Runs the REAL read_input_common with an empty .namelist.tmp (no &sbe group
  ! present at all) so this genuinely exercises the inputoutput.f90 defaults
  ! block, not a hand-set fixture (unlike the other ssbe standalone tests,
  ! which never run read_input -- see test_sbe_spinor.f90's fixture note).
  subroutine test_defaults(nfail)
    implicit none
    integer, intent(inout) :: nfail
    integer :: fh

    open(newunit=fh, file='.namelist.tmp', status='replace')
    close(fh)
    call read_input_common()

    call check_true(trim(sbe_t2_gate_shape) == 'step', 'default shape=step', nfail)
    call check_close_r(sbe_t2_gate_theta, 2d-3, 1d-15, 'default theta=2e-3', nfail)
    call check_close_r(sbe_t2_gate_width, 0d0,  1d-15, 'default width=0', nfail)
  end subroutine test_defaults

  !======================= Task 2: t2_gate_weight helper =======================
  ! Property tests for both shapes, per the "Gate weight definition" table
  ! (plan Sec. "Gate weight definition (single source of truth)"): step is a
  ! strict '>' hard gate; gauss is the Gaussian notch W=1-exp(-(dw/w)^2)
  ! (Thuemmler Eq.28, W(0)=0, quadratic onset, saturates to 1). Both shapes
  ! clamp to 0 for |delta_omega| <= floor (exact-degeneracy protection).
  subroutine test_gate_weight(nfail)
    implicit none
    integer, intent(inout) :: nfail

    ! step: strict > theta, exact-degeneracy 0
    call check_close_r(t2_gate_weight( 3d-3, 'step', 2d-3, 0d0, 1d-9), 1d0, 1d-15, &
      'step above', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step below', nfail)
    call check_close_r(t2_gate_weight( 2d-3, 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step at-threshold=skip (strict >)', nfail)
    call check_close_r(t2_gate_weight( 0d0 , 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step at 0 = 0', nfail)

    ! gauss: W(0)=0, symmetric, saturates
    call check_close_r(t2_gate_weight( 0d0 , 'gauss', 0d0, 1d-3, 1d-9), 0d0, 1d-12, &
      'gauss at 0 = 0', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'gauss', 0d0, 1d-3, 1d-9), 1d0 - exp(-1d0), 1d-12, &
      'gauss at w', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'gauss', 0d0, 1d-3, 1d-9), &
                        t2_gate_weight(-1d-3, 'gauss', 0d0, 1d-3, 1d-9), 1d-15, &
      'gauss symmetric', nfail)
    call check_true(t2_gate_weight( 5d-3, 'gauss', 0d0, 1d-3, 1d-9) > 0.99d0, &
      'gauss saturates', nfail)

    ! floor clamp (both shapes): |dw|<=floor -> 0
    call check_close_r(t2_gate_weight( 5d-10, 'gauss', 0d0, 1d-3, 1d-9), 0d0, 1d-15, &
      'gauss floor clamp', nfail)
  end subroutine test_gate_weight

  !======================= Task 3: checker predicate t2_gate_params_ok ========
  ! shape must be 'step' or 'gauss' (case-sensitive: production never
  ! lowercases sbe_deph_mode-style string keys either -- see
  ! inputoutput.f90's "convert lowercase" block, which omits sbe_* strings);
  ! theta/width must be non-negative; gauss additionally requires width>0
  ! (width=0 would divide by zero in t2_gate_weight's gauss branch).
  subroutine test_gate_params_ok(nfail)
    implicit none
    integer, intent(inout) :: nfail

    call check_true(      t2_gate_params_ok('step' , 2d-3, 0d0 ), 'valid step ok', nfail)
    call check_true(      t2_gate_params_ok('gauss', 0d0 , 1d-3), 'valid gauss ok', nfail)
    call check_true(.not. t2_gate_params_ok('bogus', 2d-3, 0d0 ), 'bad shape rejected', nfail)
    call check_true(.not. t2_gate_params_ok('gauss', 0d0 , 0d0 ), 'gauss width=0 rejected', nfail)
    call check_true(.not. t2_gate_params_ok('step' ,-1d0 , 0d0 ), 'negative theta rejected', nfail)
  end subroutine test_gate_params_ok

  !======================= Task 4: gicov_rhs integration =======================
  ! Fixture copied/adapted from test_gicov_rhs.f90 (nb=4,nk=8; block partition
  ! band1|{2,3}|band4; block {2,3} EXACTLY degenerate at eigen=0.90, isolated
  ! by a 0.60 gap from singletons band1 (0.30) and band4 (1.50)).
  !   (a) step regression: with the DEFAULT params (shape='step',
  !       theta=sbe_t2_gate_theta=2d-3, width=0), gicov_rhs's T2 term must be
  !       BIT-IDENTICAL (tol=0) to an independently-assembled reference using
  !       the literal 2d-3 threshold gicov_rhs used before this feature
  !       existed (theta_off's value) -- the byte-identical regression gate.
  !   (b) gauss: the exactly-degenerate pair {2,3} (delta_omega=0, clamped by
  !       the floor) stays undephased at E=0, and the E=0 RHS (coherent
  !       energy + T2 dephasing only, the covariant transport vanishing at
  !       zero field) is Hermitian -- the dephasing increment cannot break
  !       Hermiticity, which requires the per-pair gate weight be symmetric
  !       under (n,m)->(m,n) (W depends only on |delta_omega|).
  subroutine test_gicov_integration(nfail)
    use gs_info_ssbe,          only: s_sbe_gs_info, sbe_gs_set_replicated_kmap
    use bloch_solver_ssbe,     only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                      prepare_qnm, gicov_rhs, q_ij_from_rho
    use degenerate_block_ssbe, only: covariant_grad_block, theta_off, identity_kmap
    use salmon_global,         only: epdir_re1, am_s, num_kgrid, t_2, sbe_lg_degen, &
                                      sbe_lg_diag, yn_sbe_gw_collision, sbe_deph_mode, &
                                      sbe_lg_degen_floor, sbe_t2_gate_shape, &
                                      sbe_t2_gate_theta, sbe_t2_gate_width
    implicit none
    integer, intent(inout) :: nfail

    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe_step, sbe_gauss
    complex(8) :: rho(nb, nb, nk), drho_step(nb, nb, nk), drho_ref(nb, nb, nk)
    complex(8) :: drho_g0(nb, nb, nk), Dq(nb, nb, 3, nk)
    complex(8) :: tgt12
    real(8) :: E(3), E0(3), dk(3), step_err, gate_in, herm_err, w_gauss, gate_dist_err
    integer :: icomm, ik, ib, jb, axis

    icomm = 0
    E(1) = 0.1d0; E(2) = 0d0; E(3) = 0d0
    E0 = 0d0

    ! ---- shared salmon_global fixture (mirrors test_gicov_rhs.f90's set_globals) ----
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 10d0
    sbe_lg_diag = 0
    sbe_lg_degen = 'gicov'
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
    sbe_lg_degen_floor = 1d-9

    call build_gicov_gs(gs)
    call build_gicov_rho(rho)

    ! ================= (a) step: byte-identical regression ====================
    sbe_t2_gate_shape = 'step'
    sbe_t2_gate_theta = 2d-3    ! == theta_off's literal value (untouched elsewhere)
    sbe_t2_gate_width = 0d0

    call init_sbe_bloch_solver(sbe_step, gs, nb, icomm)
    call prepare_qnm(sbe_step, gs, icomm)
    call set_qnm_from_rho(sbe_step, rho)
    call gicov_rhs(sbe_step, gs, E, drho_step, icomm)

    do axis = 1, 3
      dk(axis) = gs%b_matrix(axis, axis) / dble(num_kgrid(axis))
    end do
    call covariant_grad_block(nb, nk, gs%nbvec, gs%bvec, num_kgrid, &
                              gs%u_transport, rho, dk, Dq, 1, nk, (/.true.,.true.,.true./), nk, identity_kmap(nk))
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          drho_ref(ib, jb, ik) = &
            (E(1)*Dq(ib,jb,1,ik) + E(2)*Dq(ib,jb,2,ik) + E(3)*Dq(ib,jb,3,ik))
          if (ib /= jb) then
            drho_ref(ib, jb, ik) = drho_ref(ib, jb, ik) &
              - zi_ * gs%delta_omega(ib, jb, ik) * rho(ib, jb, ik)
            ! the pre-feature hard gate, verbatim: strict '>' against the
            ! theta_off LITERAL (2d-3) -- NOT read from sbe_t2_gate_theta, so
            ! this reference is independent of the production code under test.
            if (abs(gs%delta_omega(ib, jb, ik)) > theta_off) then
              drho_ref(ib, jb, ik) = drho_ref(ib, jb, ik) - rho(ib, jb, ik) / t_2
            end if
          end if
        end do
      end do
    end do
    step_err = maxval(abs(drho_step - drho_ref))
    call check_close_r(step_err, 0d0, 0d0, &
      'step byte-identical to the pre-feature hard gate', nfail)

    ! ================= (b) gauss: degenerate-safe + Hermitian =================
    ! width = 0.60 == the (1,2) energy gap EXACTLY, so W(0.60) = 1-exp(-1) ~=
    ! 0.632 -- meaningfully LESS than step's hard 1 (not saturated). This is
    ! the assertion that actually distinguishes 'gauss' from 'step': before
    ! gicov_rhs branches on shape, it always applies the old unconditional
    ! step gate (W=1 whenever |delta_omega|>theta_off), so gate_dist_err below
    ! would be ~0.368/t_2 (FAIL) until Task 4's gicov_rhs change lands.
    sbe_t2_gate_shape = 'gauss'
    sbe_t2_gate_theta = 0d0     ! unused by the gauss branch
    sbe_t2_gate_width = 0.60d0

    call init_sbe_bloch_solver(sbe_gauss, gs, nb, icomm)
    call prepare_qnm(sbe_gauss, gs, icomm)
    call set_qnm_from_rho(sbe_gauss, rho)
    call gicov_rhs(sbe_gauss, gs, E0, drho_g0, icomm)

    gate_in = 0d0
    herm_err = 0d0
    do ik = 1, nk
      gate_in = max(gate_in, abs(drho_g0(2, 3, ik)), abs(drho_g0(3, 2, ik)))
      do jb = 1, nb
        do ib = 1, nb
          herm_err = max(herm_err, abs(drho_g0(ib, jb, ik) - conjg(drho_g0(jb, ib, ik))))
        end do
      end do
    end do
    call check_true(gate_in < 1d-12, &
      'gauss: degenerate pair {2,3} not dephased at E=0', nfail)
    call check_true(herm_err < 1d-10, &
      'gauss: E=0 RHS (energy+dephasing increment) Hermitian', nfail)

    w_gauss = 1d0 - exp( -(gs%delta_omega(1, 2, 1) / sbe_t2_gate_width)**2 )
    gate_dist_err = 0d0
    do ik = 1, nk
      tgt12 = (- zi_ * gs%delta_omega(1, 2, ik) - w_gauss / t_2) * rho(1, 2, ik)
      gate_dist_err = max(gate_dist_err, abs(drho_g0(1, 2, ik) - tgt12))
    end do
    call check_true(gate_dist_err < 1d-12, &
      'gauss: energy-distinct pair (1,2) gets the smooth non-saturated weight 1-exp(-1)', nfail)
  end subroutine test_gicov_integration

  !======================= step-shape exact-degeneracy floor clamp =============
  ! Spec: |Delta-omega| <= sbe_lg_degen_floor => weight 0 for BOTH shapes. The
  ! step branch's bare `abs(dw) > theta` test does NOT enforce this when
  ! theta < floor -- and theta=0d0 is checker-valid. Override the (1,4)
  ! singleton pair's delta_omega to a tiny 5d-10 (0 < |dw| <= floor=1d-9), set
  ! shape='step', theta=0d0, and assert that at E=0 (transport term vanishes)
  ! the pair carries ONLY the coherent energy term -i*dw*rho, with NO -rho/t_2
  ! dephasing. Pre-fix the bare gate subtracts rho/t_2 (~0.025 here) => FAIL;
  ! with the floor clamp it is skipped => PASS. gicov_rhs reads gs%delta_omega
  ! directly for both terms, so this override is self-consistent (gs%eigen is
  ! not read by the RHS).
  subroutine test_step_floor_clamp(nfail)
    use gs_info_ssbe,      only: s_sbe_gs_info, sbe_gs_set_replicated_kmap
    use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                  prepare_qnm, gicov_rhs
    use salmon_global,     only: epdir_re1, am_s, num_kgrid, t_2, sbe_lg_degen, &
                                  sbe_lg_diag, yn_sbe_gw_collision, sbe_deph_mode, &
                                  sbe_lg_degen_floor, sbe_t2_gate_shape, &
                                  sbe_t2_gate_theta, sbe_t2_gate_width
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    complex(8) :: rho(nb, nb, nk), drho0(nb, nb, nk), tgt14
    real(8) :: E0(3), dw_small, clamp_err
    integer :: icomm, ik

    icomm = 0
    E0 = 0d0
    dw_small = 5d-10        ! 0 < dw_small <= floor (1d-9): must be clamped to 0

    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 10d0
    sbe_lg_diag = 0
    sbe_lg_degen = 'gicov'
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
    sbe_lg_degen_floor = 1d-9

    call build_gicov_gs(gs)
    call build_gicov_rho(rho)

    ! override the (1,4)/(4,1) singleton pair to a within-floor splitting
    do ik = 1, nk
      gs%delta_omega(1, 4, ik) =  dw_small
      gs%delta_omega(4, 1, ik) = -dw_small
    end do

    sbe_t2_gate_shape = 'step'
    sbe_t2_gate_theta = 0d0          ! checker-valid; bare `>theta` would NOT clamp
    sbe_t2_gate_width = 0d0

    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call set_qnm_from_rho(sbe, rho)
    call gicov_rhs(sbe, gs, E0, drho0, icomm)

    ! correctly-clamped (1,4) pair: drho = -i*dw*rho only (no -rho/t_2).
    clamp_err = 0d0
    do ik = 1, nk
      tgt14 = - zi_ * gs%delta_omega(1, 4, ik) * rho(1, 4, ik)
      clamp_err = max(clamp_err, abs(drho0(1, 4, ik) - tgt14))
    end do
    call check_true(clamp_err < 1d-14, &
      'step: 0<|dw|<=floor pair clamped to weight 0 (not dephased)', nfail)
  end subroutine test_step_floor_clamp

  !======================= synthetic gicov gs fixture (copied/adapted from
  !======================= test_gicov_rhs.f90's build_gs) ======================
  ! SIBLING of test_gicov_integration (Fortran forbids an internal procedure
  ! from itself having a CONTAINS block), sharing nb/nk/blk/zi_/two_pi via
  ! program-level host association.
  subroutine build_gicov_gs(gs)
    use gs_info_ssbe, only: s_sbe_gs_info, sbe_gs_set_replicated_kmap
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8) :: eigen(nb), t, ang, c, s, phz
    integer :: ik, ib, jb

      gs%nk = nk; gs%nb = nb; gs%ne = 6

      call sbe_gs_set_replicated_kmap(gs, nk)   ! replicated k layout (kmap = identity)
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
        gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0
        gs%occup(3, ik) = 2d0; gs%occup(4, ik) = 0d0
        do jb = 1, nb
          do ib = 1, nb
            gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
          end do
        end do

        do ib = 1, nb
          do jb = ib + 1, nb
            phz = 0.3d0 * dble(ib) + 0.5d0 * dble(jb) + 0.2d0 * t
            gs%p_mod_matrix(ib, jb, 1, ik) = (0.4d0 + 0.1d0 * dble(ib + jb)) * exp(zi_ * phz)
            gs%p_mod_matrix(jb, ib, 1, ik) = conjg(gs%p_mod_matrix(ib, jb, 1, ik))
          end do
        end do

        gs%u_transport(1, 1, 1, ik) = (1d0, 0d0)
        gs%u_transport(4, 4, 1, ik) = (1d0, 0d0)
        ang = 0.30d0 + 0.15d0 * t
        c = cos(ang); s = sin(ang)
        phz = 0.10d0 + 0.05d0 * t
        gs%u_transport(2, 2, 1, ik) = c * exp( zi_ * phz)
        gs%u_transport(2, 3, 1, ik) = s
        gs%u_transport(3, 2, 1, ik) = -s
        gs%u_transport(3, 3, 1, ik) = c * exp(-zi_ * phz)
        do ib = 1, nb
          gs%u_transport(ib, ib, 2, ik) = (1d0, 0d0)
          gs%u_transport(ib, ib, 3, ik) = (1d0, 0d0)
        end do
      end do

      gs%d_matrix(:, :, :, :) = (0d0, 0d0)   ! X-full: gicov_rhs does not read it
  end subroutine build_gicov_gs

  !======================= smooth Hermitian test density (copied/adapted
  !======================= from test_gicov_rhs.f90's build_rho) ================
  subroutine build_gicov_rho(rho)
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
  end subroutine build_gicov_rho

  ! inject physical rho into sbe%qnm via the q_ij_from_rho bridge (SIBLING of
  ! test_gicov_integration, same reason as build_gicov_gs above)
  subroutine set_qnm_from_rho(sbe, rho)
    use bloch_solver_ssbe, only: s_sbe_bloch_solver, q_ij_from_rho
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

  !======================= assert helpers (ssbe style, copied verbatim
  !======================= from test/test_block_transport.f90) ================
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

  subroutine check_close_r(got, want, tol, label, nfail)
    implicit none
    real(8), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (abs(got - want) <= tol) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,es12.4)') "FAIL  " // label // "  |got-want|=", abs(got - want)
      nfail = nfail + 1
    end if
  end subroutine check_close_r

  subroutine check_int(got, want, label, nfail)
    implicit none
    integer, intent(in) :: got, want
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (got == want) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,2(1x,i0))') "FAIL  " // label // "  got/want=", got, want
      nfail = nfail + 1
    end if
  end subroutine check_int

end program test_t2_gate
