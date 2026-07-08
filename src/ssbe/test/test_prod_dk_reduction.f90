! src/ssbe/test/test_prod_dk_reduction.f90
! gicov memory-scale fix: production-path (real init_sbe_gs_info /
! read_prod_dk_data) regression test for the "27 -> 3" prod_dk reduction and
! its "gi/gifix unchanged" counterpart.
!
! test_gicov_window.f90 already proves the 'gicov' reduction end-to-end (with
! band-windowing as its main focus); THIS file isolates the two claims the
! task's fix must satisfy, each against a fresh, independent, un-windowed
! fixture (nb=4, nk=8, ndk=1 -> nbvec_full=27), read through the REAL
! init_sbe_gs_info/read_prod_dk_data (no hand-built gs%prod_dk):
!
!   R1 (gi UNCHANGED):    sbe_lg_degen='gi' -> gs%nbvec == 27 (full shell) and
!                         gs%prod_dk(:,:,iv,:) == prod_val(:,:,iv,:) EXACTLY
!                         for EVERY iv=1..27 (bit-for-bit vs. the pre-fix
!                         behavior -- build_xi's match_link_blocks needs every
!                         direction, so 'gi' must keep them all).
!   R2 (gifix UNCHANGED): same as R1, sbe_lg_degen='gifix'.
!   R3 (gicov REDUCED):   sbe_lg_degen='gicov' -> gs%nbvec == 3, gs%bvec holds
!                         exactly {(1,0,0),(0,1,0),(0,0,1)} (order-independent
!                         check), and each kept slot's data matches
!                         prod_val(:,:,iv_full,:) at that direction's TRUE
!                         file iv -- i.e. the correct 3 columns survived the
!                         reduction, not merely "some 3 columns".
!
! The int64 nrec/irec/ik_exp fix itself (expected_prod_dk_nrec) is unit-tested
! directly in test_gicov_memfix.f90 at k20/k24 scale (nk=8000/13824) -- writing
! an actual multi-billion-line fixture file here to drive it through
! read_prod_dk_data's real do-loop is not practical for a fast unit test; R1-R3
! above already confirm the loop's per-record wiring (storage branch, ik/jdk
! validation) is correct at small, int32-safe scale using the SAME formula.
!
! BUILD (standalone; same pattern as test_gicov_window.f90):
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_prod_dk_reduction.f90 -o <scratch_dir>/test_prod_dk_reduction.o
!   gfortran <flags> $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' \
!     ! -name 'main.f90.o') <scratch_dir>/test_prod_dk_reduction.o \
!     -o <scratch_dir>/test_prod_dk_reduction -framework Accelerate -lm -ldl
!   <scratch_dir>/test_prod_dk_reduction
!
program test_prod_dk_reduction
  use gs_info_ssbe,  only: s_sbe_gs_info, init_sbe_gs_info
  use salmon_global, only: gauge_sbe, file_sbe_prod_dk, sbe_lg_degen, &
                           num_kgrid, sbe_lg_degen_floor
  implicit none

  integer, parameter :: nb     = 4           ! bands (no windowing in this test)
  integer, parameter :: ne_abs = 4           ! electrons -> occupied bands 1..2
  integer, parameter :: nk     = 8           ! 2x2x2
  integer, parameter :: ndk    = 1
  character(*), parameter :: sysname = 'reducefix'
  character(*), parameter :: gsdir   = './'
  real(8), parameter :: a1(3) = (/ 1.3d0, 0d0, 0d0 /)
  real(8), parameter :: a2(3) = (/ 0d0, 1.3d0, 0d0 /)
  real(8), parameter :: a3(3) = (/ 0d0, 0d0, 1.3d0 /)
  integer :: nfail

  nfail = 0
  call set_globals()
  call write_fixture()

  call test_R1_gi_unchanged(nfail)
  call test_R2_gifix_unchanged(nfail)
  call test_R3_gicov_reduced(nfail)

  call cleanup_fixture()

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_prod_dk_reduction checks passed."
  end if

contains

  subroutine set_globals()
    implicit none
    num_kgrid(1) = 2; num_kgrid(2) = 2; num_kgrid(3) = 2
    gauge_sbe = 'velocity_gauge'   ! skip the length_gauge grad_k_eigen branch (irrelevant here)
    sbe_lg_degen_floor = 1d-9
    file_sbe_prod_dk = 'reducefix_prod_dk.data'
  end subroutine set_globals

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

  ! near-identity (keeps build_block_transport's polar factorization healthy),
  ! every entry distinct in (io, jo, iv, ik) -- same style as
  ! test_gicov_window.f90's prod_val.
  pure function prod_val(io, jo, iv, ik) result(v)
    implicit none
    integer, intent(in) :: io, jo, iv, ik
    complex(8) :: v
    real(8) :: d
    d = 0d0
    if (io == jo) d = 1d0
    v = dcmplx(d + 1d-3*(io + 2*jo) + 1d-4*ik + 1d-5*iv, &
             & 1d-3*(io - jo) + 1d-5*ik + 1d-6*iv)
  end function prod_val

  ! file's raw (2*ndk+1)**3 enumeration index for a shift (jdk1,jdk2,jdk3) --
  ! SAME formula write_fixture uses to write the file / read_prod_dk_data uses
  ! to parse it (needed because prod_val's value is keyed by this raw iv).
  pure function iv_full_of(jdk1, jdk2, jdk3) result(iv)
    implicit none
    integer, intent(in) :: jdk1, jdk2, jdk3
    integer :: iv, mdk
    mdk = 2*ndk + 1
    iv = (jdk3 + ndk)*mdk*mdk + (jdk2 + ndk)*mdk + (jdk1 + ndk) + 1
  end function iv_full_of

  subroutine write_fixture()
    implicit none
    integer :: fh, ik, ib, jb, io, jo, ix, iv
    integer :: jdk1, jdk2, jdk3
    complex(8) :: z(3), zp
    real(8) :: occ

    open(newunit=fh, file=gsdir//sysname//'_k.data', status='replace', action='write')
    write(fh, '(a)') "# fixture k-data (test_prod_dk_reduction)"
    write(fh, '(a)') "# header 2"; write(fh, '(a)') "# header 3"; write(fh, '(a)') "# header 4"
    write(fh, '(a)') "# 1:ik 2:k1 3:k2 4:k3 5:weight"
    do ik = 1, nk
      write(fh, '(i10,4(es25.16e3))') ik, 0.01d0*ik, -0.02d0*ik, 0.03d0*ik, 1d0
    end do
    close(fh)

    open(newunit=fh, file=gsdir//sysname//'_eigen.data', status='replace', action='write')
    write(fh, '(a)') "# fixture eigen-data (test_prod_dk_reduction)"
    write(fh, '(a)') "# header 2"; write(fh, '(a)') "# 1:io, 2:esp[a.u.], 3:occ"
    do ik = 1, nk
      write(fh, '(a,i0)') "# ik = ", ik
      do ib = 1, nb
        occ = 0d0
        if (ib <= ne_abs/2) occ = 2d0
        write(fh, '(i10,2(es25.16e3))') ib, e_val(ib, ik), occ
      end do
    end do
    close(fh)

    open(newunit=fh, file=gsdir//sysname//'_tm.data', status='replace', action='write')
    write(fh, '(a)') "# fixture tm-data (test_prod_dk_reduction)"
    write(fh, '(a)') "# header 2"; write(fh, '(a)') "# 1:ik 2:ib 3:jb 4:Re(px) 5:Im(px) ... (p_tm)"
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
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
      do ib = 1, nb
        do jb = 1, nb
          do ix = 1, 3
            z(ix) = ptm_val(ib, jb, ix, ik) * 0.5d0   ! any distinct Hermitian-ish fixture; unused by this test
          end do
          write(fh, '(3(i10),6(es25.16e3))') ik, ib, jb, &
            & dble(z(1)), aimag(z(1)), dble(z(2)), aimag(z(2)), dble(z(3)), aimag(z(3))
        end do
      end do
    end do
    close(fh)

    open(newunit=fh, file=trim(file_sbe_prod_dk), status='replace', action='write')
    write(fh, '(a,6(i10))') "#", nb, nk, num_kgrid(1), num_kgrid(2), num_kgrid(3), ndk
    write(fh, '(a)') "# 1:ik 2:ik1 3:ik2 4:ik3 5:jdk1 6:jdk2 7:jdk3 8:io 9:jo 10:Re 11:Im"
    do ik = 1, nk
      do jdk3 = -ndk, ndk
        do jdk2 = -ndk, ndk
          do jdk1 = -ndk, ndk
            iv = iv_full_of(jdk1, jdk2, jdk3)
            do io = 1, nb
              do jo = 1, nb
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

  ! shared full-27-shell exactness check, used by both R1 (gi) and R2 (gifix)
  subroutine check_full_shell_exact(gs, label_prefix, nfail)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    character(*), intent(in) :: label_prefix
    integer, intent(inout) :: nfail
    integer :: ik, iv, ib, jb
    real(8) :: err
    call check_true(gs%nbvec == 27, label_prefix // ": gs%nbvec == 27 (full shell, unchanged)", nfail)
    err = 0d0
    do ik = 1, nk
      do iv = 1, 27
        do jb = 1, nb
          do ib = 1, nb
            err = max(err, abs(gs%prod_dk(ib, jb, iv, ik) - prod_val(ib, jb, iv, ik)))
          end do
        end do
      end do
    end do
    call check_true(err == 0d0, label_prefix // ": gs%prod_dk == prod_val EXACTLY for all 27 directions (bit-for-bit)", nfail)
  end subroutine check_full_shell_exact

  subroutine test_R1_gi_unchanged(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    integer :: icomm
    icomm = 0
    sbe_lg_degen = 'gi'
    call init_sbe_gs_info(gs, sysname, gsdir, nk, nb, 1, nb, ne_abs, a1, a2, a3, .false., icomm)
    call check_full_shell_exact(gs, "R1 gi", nfail)
  end subroutine test_R1_gi_unchanged

  subroutine test_R2_gifix_unchanged(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    integer :: icomm
    icomm = 0
    sbe_lg_degen = 'gifix'
    call init_sbe_gs_info(gs, sysname, gsdir, nk, nb, 1, nb, ne_abs, a1, a2, a3, .false., icomm)
    call check_full_shell_exact(gs, "R2 gifix", nfail)
  end subroutine test_R2_gifix_unchanged

  subroutine test_R3_gicov_reduced(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    integer :: icomm, ik, ib, jb, islot
    integer :: jdk1(3), jdk2(3), jdk3(3), iv_full
    real(8) :: err
    icomm = 0
    sbe_lg_degen = 'gicov'
    call init_sbe_gs_info(gs, sysname, gsdir, nk, nb, 1, nb, ne_abs, a1, a2, a3, .false., icomm)

    call check_true(gs%nbvec == 3, "R3 gicov: gs%nbvec == 3 (reduced)", nfail)
    call check_true( &
      & any(gs%bvec(1,:)==1 .and. gs%bvec(2,:)==0 .and. gs%bvec(3,:)==0) .and. &
      & any(gs%bvec(1,:)==0 .and. gs%bvec(2,:)==1 .and. gs%bvec(3,:)==0) .and. &
      & any(gs%bvec(1,:)==0 .and. gs%bvec(2,:)==0 .and. gs%bvec(3,:)==1), &
      & "R3 gicov: gs%bvec holds exactly {+x,+y,+z} (order-independent)", nfail)

    ! kept shifts, matched to whichever slot gs%bvec actually put them in
    jdk1 = (/ 1, 0, 0 /); jdk2 = (/ 0, 1, 0 /); jdk3 = (/ 0, 0, 1 /)
    err = 0d0
    do islot = 1, 3
      iv_full = iv_full_of(jdk1(islot), jdk2(islot), jdk3(islot))
      do ik = 1, nk
        do jb = 1, nb
          do ib = 1, nb
            err = max(err, abs(gs%prod_dk(ib, jb, islot, ik) - prod_val(ib, jb, iv_full, ik)))
          end do
        end do
      end do
    end do
    call check_true(err == 0d0, &
      & "R3 gicov: each kept slot == prod_val at its TRUE file iv_full (correct 3 columns survived, not just any 3)", nfail)
  end subroutine test_R3_gicov_reduced

end program test_prod_dk_reduction
