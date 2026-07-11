! src/ssbe/test/test_sbe_spinor.f90
! Unit test (single rank, no MPI) for the SO-SBE spinor wiring
! (spin='noncollinear'): the SBE reader-side occupation convention switches
! from "2 electrons per band, occupied manifold 1..nelec/2" to "1 electron
! per SPINOR band, occupied manifold 1..nelec" (Kramers partners are separate
! band indices), and the gauge-covariant length-gauge path must NOT dephase
! coherences INSIDE a Kramers-degenerate pair (the delta-omega-gated T2 of
! gicov_rhs, commit a3be2780) while still dephasing energy-distinct pairs.
!
! Fixture: nb_file = 6 spinor bands = 3 exact Kramers pairs (1,2), (3,4),
! (5,6) at k-uniform energies e12 < e34 < e56 (pair splittings 0, inter-pair
! gaps >> theta_off); nelec = 4 -> the occupied manifold is the two lower
! pairs = bands 1..4.  prod_dk is the EXACT identity, so the X-full
! Wilson-line transport is U == I and the covariant gradient of a k-uniform
! density vanishes identically -- every nonzero RHS entry is then an exact
! closed-form expression.
!
! What is checked:
!   A  (occupation): spinor init gives gs%spinor/.focc=1/.nvb=nelec;
!      occup = 1 on bands 1..4, 0 above; VG trace = nelec = 4.
!   A2 (contrast): a SPINLESS init of the SAME files gives occup = 2 on
!      bands 1..2 only (same trace, different distribution) -- the branch
!      demonstrably has teeth.
!   B  (band window): lower-cut lo=3 freezes the bottom Kramers pair:
!      spinor nelec_eff = nelec - (lo-1)*1 = 2 (the spinless formula would
!      give 0), occup = 1 on window bands 1..2, trace = 2.
!   C  (Kramers x gated T2, gicov): with E != 0, U == I, rho k-uniform,
!      finite t_2:
!        - intra-pair coherence rho(1,2):   drho == 0 EXACTLY (delta_omega=0
!          -> no energy term, gate -> no T2, transport of k-uniform rho = 0);
!        - cross-pair coherence rho(1,3):   drho == -i*dw*rho - rho/t_2;
!        - diagonal rho(1,1):               drho == 0 (no T2 on diagonal).
!   D  (checker guards): spin='noncollinear' rejects length_gauge without
!      sbe_lg_degen='gicov', theory='maxwell_sbe', a window smaller than the
!      nelec occupied bands, odd nelec (half-occupied Kramers doublet), and
!      a frozen cut through a doublet; accepts LG+gicov and velocity_gauge;
!      the spinless config stays accepted.
!
! The fixture files (sofix_*.data) are written to the CURRENT directory and
! deleted at the end.
!
! BUILD (standalone; links the already-built build_local objects minus main,
! same pattern as test_gicov_window.f90):
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_sbe_spinor.f90 -o <scratch_dir>/test_sbe_spinor.o
!   gfortran <flags> $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' \
!     ! -name 'main.f90.o') <scratch_dir>/test_sbe_spinor.o \
!     -o <scratch_dir>/test_sbe_spinor <blas> -lm -ldl
!   <scratch_dir>/test_sbe_spinor
!
program test_sbe_spinor
  use gs_info_ssbe,      only: s_sbe_gs_info, init_sbe_gs_info
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_trace, &
                               prepare_qnm, gicov_rhs
  use input_checker_sbe, only: check_input_variables_sbe
  use salmon_global
  implicit none

  integer, parameter :: nb_file = 6          ! 3 Kramers pairs of spinor bands
  integer, parameter :: ne_so   = 4          ! spinor electrons -> occupied bands 1..4
  integer, parameter :: nk      = 64         ! 4x4x4 (>= 4 per axis for grad_k)
  integer, parameter :: ndk     = 1
  character(*), parameter :: sysname_fix = 'sofix'
  character(*), parameter :: gsdir   = './'
  real(8), parameter :: a1(3) = (/ 1.1d0, 0d0, 0d0 /)
  real(8), parameter :: a2(3) = (/ 0d0, 1.1d0, 0d0 /)
  real(8), parameter :: a3(3) = (/ 0d0, 0d0, 1.1d0 /)
  ! k-uniform Kramers-pair energies: splitting 0 inside a pair, gaps >> theta_off
  real(8), parameter :: e_pair(3) = (/ -0.30d0, -0.20d0, 0.30d0 /)
  integer :: nfail

  nfail = 0
  call set_globals_common()
  call write_fixture()
  call run_tests(nfail)
  call cleanup_fixture()

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_sbe_spinor checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals_common()
    implicit none
    num_kgrid(1) = 4; num_kgrid(2) = 4; num_kgrid(3) = 4
    sbe_lg_degen_floor = 1d-9
    file_sbe_prod_dk = 'sofix_prod_dk.data'
    yn_spinorbit = 'y'
    ! standalone tests never run the inputoutput default assignment, so every
    ! key the input checker validates must be set here explicitly
    yn_sbe_vnl_exact = 'n'
    file_sbe_vnl_kappa = ''
  end subroutine set_globals_common

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

  pure function e_val(ib) result(v)
    implicit none
    integer, intent(in) :: ib
    real(8) :: v
    v = e_pair((ib + 1) / 2)      ! bands (1,2)->pair1, (3,4)->pair2, (5,6)->pair3
  end function e_val

  pure function ptm_val(ib, jb, ix, ik) result(v)
    implicit none
    integer, intent(in) :: ib, jb, ix, ik
    complex(8) :: v
    v = dcmplx(0.1d0*ib + 0.01d0*jb + 0.001d0*ik + 1d0*ix, &
             & -0.2d0*ib + 0.02d0*jb - 0.002d0*ik + 0.1d0*ix)
  end function ptm_val

  !======================= write the synthetic exports ========================
  subroutine write_fixture()
    implicit none
    integer :: fh, ik, ib, jb, io, jo, ix, ispin
    integer :: jdk1, jdk2, jdk3
    complex(8) :: z(3), zp
    real(8) :: occ

    ! --- SYSNAME_k.data: 5 header lines + nk rows (ik, k1, k2, k3, weight)
    open(newunit=fh, file=gsdir//sysname_fix//'_k.data', status='replace', action='write')
    write(fh, '(a)') "# fixture k-data (test_sbe_spinor)"
    write(fh, '(a)') "# header 2"
    write(fh, '(a)') "# header 3"
    write(fh, '(a)') "# header 4"
    write(fh, '(a)') "# 1:ik 2:k1 3:k2 4:k3 5:weight"
    do ik = 1, nk
      write(fh, '(i10,4(es25.16e3))') ik, 0.01d0*ik, -0.02d0*ik, 0.03d0*ik, 1d0
    end do
    close(fh)

    ! --- SYSNAME_eigen.data: like write_eigen for noncollinear -- TWO
    ! duplicate spin blocks (is=1,2), each nk k-blocks of nb rows; occ = 1
    ! per occupied spinor band.  The reader must consume the FIRST block and
    ! ignore the duplicate.
    open(newunit=fh, file=gsdir//sysname_fix//'_eigen.data', status='replace', action='write')
    write(fh, '(a)') "#esp: single-particle energies (eigen energies)"
    write(fh, '(a)') "#occ: occupation numbers, io: orbital index"
    write(fh, '(a)') "# 1:io, 2:esp[a.u.], 3:occ"
    do ispin = 1, 2
      do ik = 1, nk
        write(fh, '(a,i5,a,i5)') "k=", ik, ",  spin=", ispin
        do ib = 1, nb_file
          occ = 0d0
          if (ib <= ne_so) occ = 1d0
          write(fh, '(i10,2(es25.16e3))') ib, e_val(ib), occ
        end do
      end do
    end do
    close(fh)

    ! --- SYSNAME_tm.data: 3 header lines + p block, 1 header line + rvnl block
    open(newunit=fh, file=gsdir//sysname_fix//'_tm.data', status='replace', action='write')
    write(fh, '(a)') "# fixture tm-data (test_sbe_spinor)"
    write(fh, '(a)') "# header 2"
    write(fh, '(a)') "# 1:ik 2:ib 3:jb 4:Re(px) 5:Im(px) ... (p_tm)"
    do ik = 1, nk
      do ib = 1, nb_file
        do jb = 1, nb_file
          do ix = 1, 3
            z(ix) = ptm_val(ib, jb, ix, ik)
          end do
          write(fh, '(3(i10),6(es25.16e3))') ik, ib, jb, &
            & dble(z(1)), aimag(z(1)), dble(z(2)), aimag(z(2)), dble(z(3)), aimag(z(3))
        end do
      end do
    end do
    write(fh, '(a)') "# rvnl block"
    do ik = 1, nk
      do ib = 1, nb_file
        do jb = 1, nb_file
          write(fh, '(3(i10),6(es25.16e3))') ik, ib, jb, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
        end do
      end do
    end do
    close(fh)

    ! --- prod_dk file: EXACT identity overlaps -> Wilson transport U == I
    open(newunit=fh, file=trim(file_sbe_prod_dk), status='replace', action='write')
    write(fh, '(a,6(i10))') "#", nb_file, nk, num_kgrid(1), num_kgrid(2), num_kgrid(3), ndk
    write(fh, '(a)') "# 1:ik 2:ik1 3:ik2 4:ik3 5:jdk1 6:jdk2 7:jdk3 8:io 9:jo 10:Re 11:Im"
    do ik = 1, nk
      do jdk3 = -ndk, ndk
        do jdk2 = -ndk, ndk
          do jdk1 = -ndk, ndk
            do io = 1, nb_file
              do jo = 1, nb_file
                zp = (0d0, 0d0)
                if (io == jo) zp = (1d0, 0d0)
                write(fh, '(9(i10),2(es25.16e3))') ik, 0, 0, 0, jdk1, jdk2, jdk3, &
                  & io, jo, dble(zp), aimag(zp)
              end do
            end do
          end do
        end do
      end do
    end do
    close(fh)
  end subroutine write_fixture

  subroutine cleanup_fixture()
    implicit none
    integer :: fh
    open(newunit=fh, file=gsdir//sysname_fix//'_k.data', status='old');     close(fh, status='delete')
    open(newunit=fh, file=gsdir//sysname_fix//'_eigen.data', status='old'); close(fh, status='delete')
    open(newunit=fh, file=gsdir//sysname_fix//'_tm.data', status='old');    close(fh, status='delete')
    open(newunit=fh, file=trim(file_sbe_prod_dk), status='old');            close(fh, status='delete')
  end subroutine cleanup_fixture

  !======================= the tests ==========================================
  subroutine run_tests(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs_sl, gs_so, gs_so_w, gs_sl_w, gs_cov
    type(s_sbe_bloch_solver) :: sbe_so, sbe_so_w, sbe_cov
    integer, parameter :: lo = 3            ! freeze the bottom Kramers pair
    integer :: icomm, ik, ib, jb
    real(8) :: err, tr
    complex(8), allocatable :: drho(:, :, :)
    complex(8) :: c12, c13, expct
    real(8) :: efield(3)
    logical :: ok

    icomm = 0

    ! ---- A: spinor occupation bookkeeping ----------------------------------
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    call init_sbe_gs_info(gs_so, sysname_fix, gsdir, nk, nb_file, 1, nb_file, ne_so, &
                        & a1, a2, a3, .false., icomm)
    call check_true(gs_so%spinor .and. gs_so%focc == 1d0 .and. gs_so%nvb == ne_so &
      & .and. gs_so%ne == ne_so, "A spinor init: focc=1, nvb=ne=nelec", nfail)
    err = 0d0
    do ik = 1, nk
      do ib = 1, ne_so
        err = max(err, abs(gs_so%occup(ib, ik) - 1d0))
      end do
      do ib = ne_so + 1, nb_file
        err = max(err, abs(gs_so%occup(ib, ik)))
      end do
    end do
    call check_true(err == 0d0, "A spinor occup: 1 on bands 1..nelec, 0 above (exact)", nfail)

    call init_sbe_bloch_solver(sbe_so, gs_so, nb_file, icomm)
    tr = calc_trace(sbe_so, gs_so, nb_file, icomm)
    call check_true(abs(tr - dble(ne_so)) < 1d-13, "A spinor VG trace = nelec", nfail)

    ! Kramers pairs are exactly degenerate in the fixture spectrum
    err = 0d0
    do ik = 1, nk
      err = max(err, abs(gs_so%delta_omega(1, 2, ik)))
      err = max(err, abs(gs_so%delta_omega(3, 4, ik)))
      err = max(err, abs(gs_so%delta_omega(5, 6, ik)))
    end do
    call check_true(err == 0d0, "A Kramers pairs: delta_omega == 0 inside each pair (exact)", nfail)

    ! ---- A2: spinless CONTRAST on the same files ----------------------------
    spin = 'unpolarized'
    call init_sbe_gs_info(gs_sl, sysname_fix, gsdir, nk, nb_file, 1, nb_file, ne_so, &
                        & a1, a2, a3, .false., icomm)
    ok = (.not. gs_sl%spinor) .and. gs_sl%focc == 2d0 .and. gs_sl%nvb == ne_so/2
    err = 0d0
    do ik = 1, nk
      do ib = 1, ne_so/2
        err = max(err, abs(gs_sl%occup(ib, ik) - 2d0))
      end do
      do ib = ne_so/2 + 1, nb_file
        err = max(err, abs(gs_sl%occup(ib, ik)))
      end do
    end do
    call check_true(ok .and. err == 0d0, &
      & "A2 spinless contrast: occup = 2 on bands 1..nelec/2 only (branch has teeth)", nfail)

    ! ---- B: band window (freeze the bottom Kramers pair) --------------------
    spin = 'noncollinear'
    call init_sbe_gs_info(gs_so_w, sysname_fix, gsdir, nk, nb_file, lo, nb_file, ne_so, &
                        & a1, a2, a3, .false., icomm)
    call check_true(gs_so_w%ne == ne_so - (lo - 1) .and. gs_so_w%nvb == ne_so - (lo - 1), &
      & "B spinor window lo=3: nelec_eff = nelec - (lo-1) = 2 (spinless formula would give 0)", nfail)
    err = 0d0
    do ik = 1, nk
      do ib = 1, gs_so_w%nvb
        err = max(err, abs(gs_so_w%occup(ib, ik) - 1d0))
      end do
      do ib = gs_so_w%nvb + 1, gs_so_w%nb
        err = max(err, abs(gs_so_w%occup(ib, ik)))
      end do
    end do
    call check_true(err == 0d0, "B spinor window: occup = 1 exactly on window bands 1..nelec_eff", nfail)
    call init_sbe_bloch_solver(sbe_so_w, gs_so_w, gs_so_w%nb, icomm)
    tr = calc_trace(sbe_so_w, gs_so_w, gs_so_w%nb, icomm)
    call check_true(abs(tr - dble(gs_so_w%ne)) < 1d-13, "B spinor window trace = nelec_eff", nfail)

    ! spinless window contrast: lo=3 freezes 2 bands x 2 electrons = ALL of nelec
    spin = 'unpolarized'
    call init_sbe_gs_info(gs_sl_w, sysname_fix, gsdir, nk, nb_file, lo, nb_file, ne_so, &
                        & a1, a2, a3, .false., icomm)
    call check_true(gs_sl_w%ne == 0 .and. gs_sl_w%nvb == 0 .and. maxval(abs(gs_sl_w%occup)) == 0d0, &
      & "B spinless window contrast: lo=3 freezes the whole valence manifold (ne_eff = 0)", nfail)

    ! ---- C: Kramers x delta-omega-gated T2 (gicov RHS) ----------------------
    spin = 'noncollinear'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov'
    t_2 = 10d0
    yn_sbe_gw_collision = 'n'
    sbe_lg_diag = 0
    epdir_re1(1:3) = (/ 1d0, 0d0, 0d0 /)
    call init_sbe_gs_info(gs_cov, sysname_fix, gsdir, nk, nb_file, 1, nb_file, ne_so, &
                        & a1, a2, a3, .false., icomm)
    call init_sbe_bloch_solver(sbe_cov, gs_cov, nb_file, icomm)
    call prepare_qnm(sbe_cov, gs_cov, icomm)

    ! k-UNIFORM density: ground-state diagonal + one intra-pair and one
    ! cross-pair coherence at EVERY k (U == I => covariant transport of a
    ! k-uniform rho vanishes identically; prepare_qnm sets exp_iphi = 1 on
    ! every off-diagonal pair in gicov, so qnm IS the physical rho).
    c12 = (0.05d0, 0.02d0)
    c13 = (0.03d0, -0.04d0)
    do ik = sbe_cov%ik_min, sbe_cov%ik_max
      sbe_cov%qnm(:, :, ik) = (0d0, 0d0)
      do ib = 1, nb_file
        sbe_cov%qnm(ib, ib, ik) = gs_cov%occup(ib, ik)
      end do
      sbe_cov%qnm(1, 2, ik) = c12
      sbe_cov%qnm(2, 1, ik) = conjg(c12)
      sbe_cov%qnm(1, 3, ik) = c13
      sbe_cov%qnm(3, 1, ik) = conjg(c13)
    end do

    allocate(drho(nb_file, nb_file, nk))
    efield(1:3) = (/ 0.1d0, 0d0, 0d0 /)
    call gicov_rhs(sbe_cov, gs_cov, efield, drho, icomm)

    ! intra-pair (Kramers) coherence: NO dephasing, NO energy rotation -> 0
    err = 0d0
    do ik = 1, nk
      err = max(err, abs(drho(1, 2, ik)))
      err = max(err, abs(drho(2, 1, ik)))
    end do
    call check_true(err == 0d0, &
      & "C gated T2: drho == 0 EXACTLY for the coherence INSIDE a Kramers pair", nfail)

    ! cross-pair coherence: full closed form -i*dw*rho - rho/t_2
    err = 0d0
    do ik = 1, nk
      expct = -dcmplx(0d0, 1d0) * gs_cov%delta_omega(1, 3, ik) * c13 - c13 / t_2
      err = max(err, abs(drho(1, 3, ik) - expct))
    end do
    call check_true(err < 1d-15, &
      & "C gated T2: drho == -i*dw*rho - rho/t_2 for the ENERGY-DISTINCT pair", nfail)

    ! diagonal: no T2, transport 0 -> 0
    err = 0d0
    do ik = 1, nk
      do ib = 1, nb_file
        err = max(err, abs(drho(ib, ib, ik)))
      end do
    end do
    call check_true(err == 0d0, "C diagonal: drho == 0 (k-uniform, U == I)", nfail)
    deallocate(drho)

    ! ---- D: input-checker guards --------------------------------------------
    call set_checker_defaults()

    spin = 'noncollinear'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'off'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: spinor + length_gauge + sbe_lg_degen='off' is REJECTED", nfail)

    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov'
    call check_true(check_input_variables_sbe(), &
      & "D checker: spinor + length_gauge + sbe_lg_degen='gicov' is accepted", nfail)

    ! ---- gicov_int (integral covariant-Houston transport) checker guards -----
    ! Shares the X-full spinor accept path with gicov (uses_xfull_links), plus
    ! the v1 scope reject group.  The reject cases use an UNPOLARIZED baseline so
    ! the gicov_int rule is the one under test (a noncollinear baseline would let
    ! the separate spinor+maxwell rule fire on the theory='maxwell_sbe' case).
    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov_int'
    yn_sbe_gw_collision = 'n'
    call check_true(check_input_variables_sbe(), &
      & "D checker: spinor + length_gauge + sbe_lg_degen='gicov_int' is accepted", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov_int'
    yn_sbe_gw_collision = 'n'
    call check_true(check_input_variables_sbe(), &
      & "D checker: unpolarized + length_gauge + gicov_int is accepted (baseline)", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'gicov_int'
    yn_sbe_gw_collision = 'n'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: gicov_int + velocity_gauge is REJECTED (requires length_gauge)", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov_int'
    theory = 'maxwell_sbe'
    yn_sbe_gw_collision = 'n'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: gicov_int + theory='maxwell_sbe' is REJECTED (v1 scope)", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov_int'
    yn_sbe_gw_collision = 'y'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: gicov_int + yn_sbe_gw_collision='y' is REJECTED (v1 scope)", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'gicov_int'
    yn_sbe_gw_collision = 'n'
    yn_sbe_gs_current_subtract = 'y'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: gicov_int + yn_sbe_gs_current_subtract='y' is REJECTED (v1 scope)", nfail)

    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    call check_true(check_input_variables_sbe(), &
      & "D checker: spinor + velocity_gauge is accepted", nfail)

    ! spinor: nstate_sbe must reach the SPINOR occupied count nelec (not nelec/2)
    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    nstate_sbe(1) = ne_so - 1                ! < nelec: window misses occupied bands
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: spinor window smaller than the nelec occupied bands is REJECTED", nfail)

    ! spinor + maxwell_sbe is rejected (multiscale bookkeeping unaudited)
    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    theory = 'maxwell_sbe'
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: spinor + theory='maxwell_sbe' is REJECTED", nfail)

    ! Kramers-pair alignment: odd nelec = a half-occupied doublet
    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    nelec = ne_so - 1
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: spinor odd nelec (half-occupied Kramers doublet) is REJECTED", nfail)

    ! Kramers-pair alignment: a frozen cut through the doublet (1,2)
    call set_checker_defaults()
    spin = 'noncollinear'
    gauge_sbe = 'velocity_gauge'
    sbe_lg_degen = 'off'
    nband_sbe_min = 2
    call check_true(.not. check_input_variables_sbe(), &
      & "D checker: spinor frozen cut splitting a Kramers doublet is REJECTED", nfail)

    call set_checker_defaults()
    spin = 'unpolarized'
    gauge_sbe = 'length_gauge'
    sbe_lg_degen = 'off'
    call check_true(check_input_variables_sbe(), &
      & "D checker: the spinless config stays accepted (regression)", nfail)
  end subroutine run_tests

  ! Minimal self-consistent salmon_global state for check_input_variables_sbe
  ! (standalone tests never run read_input, so every global the checker reads
  ! is set explicitly here -- including the maxwell_sbe fdtd knobs, so the
  ! maxwell rejection case raises ONLY the spinor error).
  subroutine set_checker_defaults()
    implicit none
    theory = 'sbe'
    num_sbe = 1
    sysname = 'sofix'
    sysname_sbe(1) = 'default'
    al(1:3) = (/ 1.1d0, 1.1d0, 1.1d0 /)
    al_vec1(1:3) = 0d0; al_vec2(1:3) = 0d0; al_vec3(1:3) = 0d0
    al_sbe(1:3, 1) = 0d0
    al_vec1_sbe(1:3, 1) = 0d0; al_vec2_sbe(1:3, 1) = 0d0; al_vec3_sbe(1:3, 1) = 0d0
    nstate = nb_file
    nelec = ne_so
    nstate_sbe(1) = nb_file
    nk_sbe(1) = 0
    nelec_sbe(1) = 0
    nband_sbe_min = 1
    t_2 = 10d0
    am_s = 4
    sbe_lg_diag = 0
    sbe_lg_degen_floor = 1d-9
    yn_sbe_gs_current_subtract = 'n'
    norder_correction = 0
    yn_spinorbit = 'y'
    ! maxwell_sbe branch inputs (valid values -> no extra raises there)
    nx_m = 1; ny_m = 1; nz_m = 1
    nxvac_m(1:2) = 0; nxvacl_m = 0; nxvacr_m = 0
    hx_m = 1d0; hy_m = 1d0; hz_m = 1d0
    dl_em(1:3) = 0d0
    fdtddim = '1d'
  end subroutine set_checker_defaults

end program test_sbe_spinor
