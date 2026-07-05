! src/ssbe/test/test_gicov_window.f90
! Unit test (single rank, no MPI) for the ssbe band-window selection
! (nband_sbe_min lower cut + nstate_sbe upper cut; init_sbe_gs_info's
! [nb_min, nb_hi] window, READER-side windowing): the exports are read in
! FULL (all nb_file records consumed, preserving file alignment) and only the
! contiguous window [lo : hi] is stored into gs%* (window index w <-> absolute
! band w + lo - 1); the frozen bands 1..lo-1 are inert fully-occupied and
! enter the bookkeeping only through gs%ne = ne - 2*(lo-1).
!
! What is checked (fixture: synthetic full exports, nb_file=6, ne=8 so the
! occupied manifold is bands 1..4; grid 4x4x4, prod_dk ndk=1; NOTE the grid
! must have >= 4 points per axis -- grad_k_array_nb1d_double's kdx wrap table
! mod(j+N-1,N)+1 hits index 0 for N<4 (Fortran mod of negative), a
! PRE-EXISTING out-of-bounds never triggered by production grids):
!
!   A (slicing):     a lo=3 windowed read of the SAME files equals the
!                    hand-sliced lo=1 (ground-truth, existing-behavior) read
!                    EXACTLY (bit-for-bit: same text -> same doubles) for
!                    eigen / p_tm / rvnl_tm / prod_dk / delta_omega /
!                    grad_k_eigen / kpoint / kweight / bvec.
!   A2 (absolute):   windowed entries also equal the closed-form fixture
!                    value at the ABSOLUTE band index (catches a slicing bug
!                    that shifts full and windowed reads the same way).
!   B (occupation):  gs%ne = nelec_eff = ne - 2*(lo-1); occupied window bands
!                    are exactly 1..nelec_eff/2 (absolute lo..ne/2).
!   C (teeth):       lo=3 cuts 2 occupied bands -> the propagated trace drops
!                    by EXACTLY 2*(lo-1) = 4 (8 -> 4 = nelec_eff); n_ex
!                    bookkeeping vs nelec_eff is 0 at rest while the WRONG
!                    absolute-nelec bookkeeping would report 2*(lo-1) fake
!                    holes -- the window selection demonstrably has teeth.
!   D (edge):        lo = ne/2+1 freezes the whole valence manifold ->
!                    nelec_eff = 0, occup all zero, slice still correct.
!   E (dual cut):    window [2, 5] (lower AND upper cut, hi < nb_file) slices
!                    both edges correctly; bands above hi are excluded.
!
! The fixture files (winfix_*.data) are written to the CURRENT directory and
! deleted at the end.
!
! BUILD (standalone; links the already-built build_local objects minus main,
! same pattern as test_gicov_rhs.f90):
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_gicov_window.f90 -o <scratch_dir>/test_gicov_window.o
!   gfortran <flags> $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' \
!     ! -name 'main.f90.o') <scratch_dir>/test_gicov_window.o \
!     -o <scratch_dir>/test_gicov_window -framework Accelerate -lm -ldl
!   <scratch_dir>/test_gicov_window
!
program test_gicov_window
  use gs_info_ssbe,      only: s_sbe_gs_info, init_sbe_gs_info
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_trace
  use salmon_global,     only: gauge_sbe, file_sbe_prod_dk, sbe_lg_degen, &
                               num_kgrid, sbe_lg_degen_floor
  implicit none

  integer, parameter :: nb_file = 6          ! bands in the synthetic exports
  integer, parameter :: ne_abs  = 8          ! electrons -> occupied bands 1..4
  integer, parameter :: lo      = 3          ! window lower edge under test
  integer, parameter :: nk      = 64         ! 4x4x4
  integer, parameter :: ndk     = 1
  integer, parameter :: nbvec   = 27         ! (2*ndk+1)**3
  character(*), parameter :: sysname = 'winfix'
  character(*), parameter :: gsdir   = './'
  real(8), parameter :: a1(3) = (/ 1.1d0, 0d0, 0d0 /)
  real(8), parameter :: a2(3) = (/ 0d0, 1.1d0, 0d0 /)
  real(8), parameter :: a3(3) = (/ 0d0, 0d0, 1.1d0 /)
  integer :: nfail

  nfail = 0
  call set_globals()
  call write_fixture()
  call run_tests(nfail)
  call cleanup_fixture()

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_window checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    num_kgrid(1) = 4; num_kgrid(2) = 4; num_kgrid(3) = 4
    gauge_sbe = 'length_gauge'          ! exercises the grad_k_eigen slice too
    sbe_lg_degen = 'gicov'              ! exercises read_prod_dk_data + transport
    sbe_lg_degen_floor = 1d-9
    file_sbe_prod_dk = 'winfix_prod_dk.data'
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

  !======================= closed-form fixture values =========================
  ! (functions of the ABSOLUTE band index, so A2 can catch a consistent shift)
  pure function e_val(ib, ik) result(v)
    implicit none
    integer, intent(in) :: ib, ik
    real(8) :: v
    v = dble(ib) + 0.01d0 * dble(ik)
  end function e_val

  pure function ptm_val(ib, jb, ix, ik) result(v)
    implicit none
    integer, intent(in) :: ib, jb, ix, ik
    complex(8) :: v
    v = dcmplx(0.1d0*ib + 0.01d0*jb + 0.001d0*ik + 1d0*ix, &
             & -0.2d0*ib + 0.02d0*jb - 0.002d0*ik + 0.1d0*ix)
  end function ptm_val

  pure function rvnl_val(ib, jb, ix, ik) result(v)
    implicit none
    integer, intent(in) :: ib, jb, ix, ik
    complex(8) :: v
    v = dcmplx(0.3d0*ib - 0.03d0*jb + 0.004d0*ik - 0.5d0*ix, &
             & 0.05d0*ib + 0.4d0*jb + 0.006d0*ik + 0.2d0*ix)
  end function rvnl_val

  pure function prod_val(io, jo, iv, ik) result(v)
    implicit none
    integer, intent(in) :: io, jo, iv, ik
    complex(8) :: v
    real(8) :: d
    d = 0d0
    if (io == jo) d = 1d0
    ! near-identity (keeps the polar factorization in build_block_transport
    ! healthy) with every entry distinct in (io, jo, iv, ik)
    v = dcmplx(d + 1d-3*(io + 2*jo) + 1d-4*ik + 1d-5*iv, &
             & 1d-3*(io - jo) + 1d-5*ik + 1d-6*iv)
  end function prod_val

  !======================= write the synthetic exports ========================
  subroutine write_fixture()
    implicit none
    integer :: fh, ik, ib, jb, io, jo, ix, iv
    integer :: jdk1, jdk2, jdk3
    complex(8) :: z(3), zp
    real(8) :: occ

    ! --- SYSNAME_k.data: 5 header lines + nk rows (ik, k1, k2, k3, weight)
    open(newunit=fh, file=gsdir//sysname//'_k.data', status='replace', action='write')
    write(fh, '(a)') "# fixture k-data (test_gicov_window)"
    write(fh, '(a)') "# header 2"
    write(fh, '(a)') "# header 3"
    write(fh, '(a)') "# header 4"
    write(fh, '(a)') "# 1:ik 2:k1 3:k2 4:k3 5:weight"
    do ik = 1, nk
      write(fh, '(i10,4(es25.16e3))') ik, 0.01d0*ik, -0.02d0*ik, 0.03d0*ik, 1d0
    end do
    close(fh)

    ! --- SYSNAME_eigen.data: 3 header lines, per ik: 1 header + nb rows
    open(newunit=fh, file=gsdir//sysname//'_eigen.data', status='replace', action='write')
    write(fh, '(a)') "# fixture eigen-data (test_gicov_window)"
    write(fh, '(a)') "# header 2"
    write(fh, '(a)') "# 1:io, 2:esp[a.u.], 3:occ"
    do ik = 1, nk
      write(fh, '(a,i0)') "# ik = ", ik
      do ib = 1, nb_file
        occ = 0d0
        if (ib <= ne_abs/2) occ = 2d0
        write(fh, '(i10,2(es25.16e3))') ib, e_val(ib, ik), occ
      end do
    end do
    close(fh)

    ! --- SYSNAME_tm.data: 3 header lines + p block, 1 header line + rvnl block
    open(newunit=fh, file=gsdir//sysname//'_tm.data', status='replace', action='write')
    write(fh, '(a)') "# fixture tm-data (test_gicov_window)"
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
          do ix = 1, 3
            z(ix) = rvnl_val(ib, jb, ix, ik)
          end do
          write(fh, '(3(i10),6(es25.16e3))') ik, ib, jb, &
            & dble(z(1)), aimag(z(1)), dble(z(2)), aimag(z(2)), dble(z(3)), aimag(z(3))
        end do
      end do
    end do
    close(fh)

    ! --- prod_dk file: metadata line + legend + nk*nbvec*no*no records
    ! (ik outermost/contiguous, as the writer emits and the reader validates)
    open(newunit=fh, file=trim(file_sbe_prod_dk), status='replace', action='write')
    write(fh, '(a,6(i10))') "#", nb_file, nk, num_kgrid(1), num_kgrid(2), num_kgrid(3), ndk
    write(fh, '(a)') "# 1:ik 2:ik1 3:ik2 4:ik3 5:jdk1 6:jdk2 7:jdk3 8:io 9:jo 10:Re 11:Im"
    do ik = 1, nk
      do jdk3 = -ndk, ndk
        do jdk2 = -ndk, ndk
          do jdk1 = -ndk, ndk
            iv = (jdk3 + ndk)*(2*ndk+1)**2 + (jdk2 + ndk)*(2*ndk+1) + (jdk1 + ndk) + 1
            do io = 1, nb_file
              do jo = 1, nb_file
                zp = prod_val(io, jo, iv, ik)
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
    open(newunit=fh, file=gsdir//sysname//'_k.data', status='old');     close(fh, status='delete')
    open(newunit=fh, file=gsdir//sysname//'_eigen.data', status='old'); close(fh, status='delete')
    open(newunit=fh, file=gsdir//sysname//'_tm.data', status='old');    close(fh, status='delete')
    open(newunit=fh, file=trim(file_sbe_prod_dk), status='old');        close(fh, status='delete')
  end subroutine cleanup_fixture

  !======================= the tests ==========================================
  subroutine run_tests(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs_full, gs_win, gs_edge, gs_dual
    type(s_sbe_bloch_solver) :: sbe_full, sbe_win
    integer :: icomm, ik, ib, jb, iv, ix, w
    integer :: nb_eff, ne_eff, lo_edge
    real(8) :: err, errf, occ_err
    real(8) :: tr_full, tr_win, tr_win_vb

    icomm = 0
    nb_eff = nb_file - (lo - 1)          ! = 4
    ne_eff = ne_abs - 2*(lo - 1)         ! = 4

    ! ground truth: full window [1, nb_file] (existing behavior)
    call init_sbe_gs_info(gs_full, sysname, gsdir, nk, nb_file, 1, nb_file, ne_abs, &
                        & a1, a2, a3, .false., icomm)
    ! windowed read of the SAME files: [lo, nb_file]
    call init_sbe_gs_info(gs_win, sysname, gsdir, nk, nb_file, lo, nb_file, ne_abs, &
                        & a1, a2, a3, .false., icomm)

    call check_true(gs_full%nb == nb_file .and. gs_full%ne == ne_abs, &
      & "lo=1: gs%nb/gs%ne keep the unwindowed values", nfail)
    call check_true(gs_win%nb == nb_eff .and. gs_win%ne == ne_eff, &
      & "lo=3: gs%nb = nb_eff = hi-lo+1 and gs%ne = nelec_eff = nelec-2*(lo-1)", nfail)
    call check_true(size(gs_win%eigen, 1) == nb_eff .and. &
      &             size(gs_win%p_tm_matrix, 1) == nb_eff .and. &
      &             size(gs_win%prod_dk, 1) == nb_eff .and. &
      &             size(gs_win%u_transport, 1) == nb_eff, &
      & "lo=3: gs arrays allocated at nb_eff", nfail)

    ! ---- A: windowed read == hand-sliced full read (EXACT) -----------------
    err = 0d0
    do ik = 1, nk
      do w = 1, nb_eff
        err = max(err, abs(gs_win%eigen(w, ik) - gs_full%eigen(w + lo - 1, ik)))
      end do
    end do
    call check_true(err == 0d0, "A eigen: window == full[lo:hi] (exact)", nfail)

    err = 0d0
    do ik = 1, nk
      do ix = 1, 3
        do jb = 1, nb_eff
          do ib = 1, nb_eff
            err = max(err, abs(gs_win%p_tm_matrix(ib, jb, ix, ik) &
              &              - gs_full%p_tm_matrix(ib + lo - 1, jb + lo - 1, ix, ik)))
            err = max(err, abs(gs_win%rvnl_tm_matrix(ib, jb, ix, ik) &
              &              - gs_full%rvnl_tm_matrix(ib + lo - 1, jb + lo - 1, ix, ik)))
          end do
        end do
      end do
    end do
    call check_true(err == 0d0, "A tm: p_tm/rvnl_tm window == full[lo:hi, lo:hi] (exact)", nfail)

    call check_true(gs_win%nbvec == nbvec .and. gs_full%nbvec == nbvec .and. &
      & all(gs_win%bvec == gs_full%bvec), "A prod_dk: nbvec/bvec identical", nfail)
    err = 0d0
    do ik = 1, nk
      do iv = 1, nbvec
        do jb = 1, nb_eff
          do ib = 1, nb_eff
            err = max(err, abs(gs_win%prod_dk(ib, jb, iv, ik) &
              &              - gs_full%prod_dk(ib + lo - 1, jb + lo - 1, iv, ik)))
          end do
        end do
      end do
    end do
    call check_true(err == 0d0, "A prod_dk: window == full[lo:hi, lo:hi] (exact)", nfail)

    err = 0d0
    do ik = 1, nk
      do jb = 1, nb_eff
        do ib = 1, nb_eff
          err = max(err, abs(gs_win%delta_omega(ib, jb, ik) &
            &              - gs_full%delta_omega(ib + lo - 1, jb + lo - 1, ik)))
        end do
      end do
    end do
    call check_true(err == 0d0, "A delta_omega: consistent with the sliced eigen (exact)", nfail)

    err = 0d0
    do ik = 1, nk
      do ix = 1, 3
        do w = 1, nb_eff
          err = max(err, abs(gs_win%grad_k_eigen(w, ix, ik) &
            &              - gs_full%grad_k_eigen(w + lo - 1, ix, ik)))
        end do
      end do
    end do
    call check_true(err == 0d0, "A grad_k_eigen: per-band gradient commutes with the slice (exact)", nfail)

    call check_true(all(gs_win%kpoint == gs_full%kpoint) .and. &
      &             all(gs_win%kweight == gs_full%kweight), &
      & "A k-data: kpoint/kweight unaffected by the window", nfail)

    ! ---- A2: windowed entries == closed-form value at the ABSOLUTE band ----
    errf = 0d0
    do ik = 1, nk
      do w = 1, nb_eff
        errf = max(errf, abs(gs_win%eigen(w, ik) - e_val(w + lo - 1, ik)))
      end do
      do ix = 1, 3
        do jb = 1, nb_eff
          do ib = 1, nb_eff
            errf = max(errf, abs(gs_win%p_tm_matrix(ib, jb, ix, ik) &
              &                - ptm_val(ib + lo - 1, jb + lo - 1, ix, ik)))
            errf = max(errf, abs(gs_win%rvnl_tm_matrix(ib, jb, ix, ik) &
              &                - rvnl_val(ib + lo - 1, jb + lo - 1, ix, ik)))
          end do
        end do
      end do
      do iv = 1, nbvec
        do jb = 1, nb_eff
          do ib = 1, nb_eff
            errf = max(errf, abs(gs_win%prod_dk(ib, jb, iv, ik) &
              &                - prod_val(ib + lo - 1, jb + lo - 1, iv, ik)))
          end do
        end do
      end do
    end do
    call check_true(errf < 1d-13, "A2 absolute-index check: window[w] == fixture(w+lo-1) " // &
      & "(catches a consistent full+window shift)", nfail)

    ! ---- B: occupation / nelec_eff bookkeeping ------------------------------
    occ_err = 0d0
    do ik = 1, nk
      do ib = 1, ne_abs/2
        occ_err = max(occ_err, abs(gs_full%occup(ib, ik) - 2d0))
      end do
      do ib = ne_abs/2 + 1, nb_file
        occ_err = max(occ_err, abs(gs_full%occup(ib, ik)))
      end do
    end do
    call check_true(occ_err == 0d0, "B lo=1: occup = 2 on bands 1..ne/2, 0 above", nfail)

    occ_err = 0d0
    do ik = 1, nk
      do ib = 1, ne_eff/2                       ! absolute lo..ne/2
        occ_err = max(occ_err, abs(gs_win%occup(ib, ik) - 2d0))
      end do
      do ib = ne_eff/2 + 1, nb_eff              ! absolute ne/2+1..nb_file
        occ_err = max(occ_err, abs(gs_win%occup(ib, ik)))
      end do
    end do
    call check_true(occ_err == 0d0, "B lo=3: occup = 2 exactly on window bands 1..nelec_eff/2", nfail)

    ! ---- C: teeth -- the lower-cut demonstrably changes the trace ----------
    gauge_sbe = 'velocity_gauge'                ! calc_trace reads sbe%rho
    call init_sbe_bloch_solver(sbe_full, gs_full, nb_file, icomm)
    call init_sbe_bloch_solver(sbe_win,  gs_win,  nb_eff,  icomm)
    tr_full   = calc_trace(sbe_full, gs_full, nb_file, icomm)
    tr_win    = calc_trace(sbe_win,  gs_win,  nb_eff,  icomm)
    tr_win_vb = calc_trace(sbe_win,  gs_win,  ne_eff/2, icomm)
    call check_true(abs(tr_full - dble(ne_abs)) < 1d-12, &
      & "C lo=1: propagated trace = nelec = 8", nfail)
    call check_true(abs(tr_win - dble(ne_eff)) < 1d-12, &
      & "C teeth: lo=3 cuts 2 occupied bands -> trace = nelec_eff = 4 " // &
      & "(changes by exactly -2*(lo-1))", nfail)
    call check_true(abs((tr_win - tr_win_vb)) < 1d-12, &
      & "C n_ex bookkeeping: tr_all - tr_vb(nelec_eff/2 bands) = 0 at rest", nfail)
    call check_true(abs(dble(ne_abs) - tr_win_vb - 2d0*(lo - 1)) < 1d-12, &
      & "C teeth: ABSOLUTE-nelec bookkeeping would report 2*(lo-1) fake holes " // &
      & "(nelec_eff is required)", nfail)

    ! ---- D: edge -- freeze the entire valence manifold ----------------------
    gauge_sbe = 'length_gauge'
    lo_edge = ne_abs/2 + 1                      ! = 5: window = absolute 5..6
    call init_sbe_gs_info(gs_edge, sysname, gsdir, nk, nb_file, lo_edge, nb_file, ne_abs, &
                        & a1, a2, a3, .false., icomm)
    err = 0d0
    do ik = 1, nk
      do w = 1, nb_file - (lo_edge - 1)
        err = max(err, abs(gs_edge%eigen(w, ik) - e_val(w + lo_edge - 1, ik)))
      end do
    end do
    call check_true(gs_edge%ne == 0 .and. maxval(abs(gs_edge%occup)) == 0d0 .and. &
      & err < 1d-13, "D edge lo=ne/2+1: nelec_eff = 0, occup all zero, slice still correct", nfail)

    ! ---- E: dual cut -- lower AND upper edge, hi < nb_file -------------------
    call init_sbe_gs_info(gs_dual, sysname, gsdir, nk, nb_file, 2, 5, ne_abs, &
                        & a1, a2, a3, .false., icomm)
    err = 0d0
    do ik = 1, nk
      do w = 1, 4                                       ! absolute bands 2..5
        err = max(err, abs(gs_dual%eigen(w, ik) - gs_full%eigen(w + 1, ik)))
      end do
      do ix = 1, 3
        do jb = 1, 4
          do ib = 1, 4
            err = max(err, abs(gs_dual%p_tm_matrix(ib, jb, ix, ik) &
              &              - gs_full%p_tm_matrix(ib + 1, jb + 1, ix, ik)))
            err = max(err, abs(gs_dual%rvnl_tm_matrix(ib, jb, ix, ik) &
              &              - gs_full%rvnl_tm_matrix(ib + 1, jb + 1, ix, ik)))
          end do
        end do
      end do
      do iv = 1, nbvec
        do jb = 1, 4
          do ib = 1, 4
            err = max(err, abs(gs_dual%prod_dk(ib, jb, iv, ik) &
              &              - gs_full%prod_dk(ib + 1, jb + 1, iv, ik)))
          end do
        end do
      end do
    end do
    occ_err = 0d0
    do ik = 1, nk
      do ib = 1, 3                                      ! occupied: absolute 2..4
        occ_err = max(occ_err, abs(gs_dual%occup(ib, ik) - 2d0))
      end do
      occ_err = max(occ_err, abs(gs_dual%occup(4, ik))) ! absolute 5: empty
    end do
    call check_true(gs_dual%nb == 4 .and. gs_dual%ne == ne_abs - 2 .and. &
      & err == 0d0 .and. occ_err == 0d0, &
      & "E dual cut [2,5]: hi < nb_file sliced exactly on both edges; occup/ne consistent", nfail)
  end subroutine run_tests

end program test_gicov_window
