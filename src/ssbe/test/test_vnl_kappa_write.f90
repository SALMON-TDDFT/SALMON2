! src/ssbe/test/test_vnl_kappa_write.f90
! I/O speedup for write_sbe_vnl_kappa_data (src/io/write.f90): the per-k-block
! root write loop was replaced by ONE write() per block instead of
! 2*(2*ns+1) small writes (V_s, then W_{1:3,s}, for every ik in the block).
! The on-disk layout (byte order, file size, header) MUST be unchanged --
! read_vnl_kappa_data (gs_info_ssbe.f90) is untouched and still parses the
! file as a header followed by, for every ik, (2*ns+1) pairs of
! [V_s block][W_{1:3,s} block].
!
! This is a STANDALONE test (no MPI, no SALMON derived types): it mirrors,
! byte for byte, the exact two write patterns that exist in write.f90 around
! src/io/write.f90:1757-1778 --
!
!   OLD (before this change):
!     do ik = ikb_s, ikb_e
!       do s = -ns, ns
!         write(fh) buf(1:NB,1:NB,0,s,ik)        ! V_s
!         write(fh) buf(1:NB,1:NB,1:3,s,ik)      ! W_{1:3,s}
!       end do
!     end do
!
!   NEW (after this change):
!     write(fh) buf(1:NB,1:NB,0:3,-ns:ns,ikb_s:ikb_e)
!
! against a synthetic buf(NB,NB,0:3,-ns:ns,1:NK) with a k-block size (nblk)
! that does NOT evenly divide NK, so the last block is a partial block --
! the one shape edge case in the real subroutine's `do ikb_s = 1, NK, nblk`
! loop.  Checks:
!
!   1. File size:  old.bin and new.bin (156-byte header + body) both equal
!      the reader's own size fingerprint formula
!        expect = 156 + NK*(2ns+1)*4*NB^2*16   (gs_info_ssbe.f90:728)
!   2. Byte-for-byte identity: old.bin == new.bin (the actual invariant this
!      change must preserve).
!   3. Negative control: the byte-compare routine used in (2) actually
!      detects a difference when one byte is deliberately flipped (the check
!      cannot pass vacuously just because both files have the same size).
!   4. Round-trip: reading new.bin with reader-style per-(ik,s) sequential
!      reads [V_s][W_{1:3,s}] recovers buf exactly (bit-exact, no arithmetic
!      in the round trip -- a pure copy must come back identical).
!
! BUILD (fully standalone; no SALMON objects needed):
!   gfortran -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     src/ssbe/test/test_vnl_kappa_write.f90 -o test_vnl_kappa_write
!   ./test_vnl_kappa_write
!
program test_vnl_kappa_write
  implicit none

  integer, parameter :: NB = 5, ns = 6, NK = 8, nblk = 3   ! nblk does NOT divide NK
  integer, parameter :: nkap = 2*ns + 1
  complex(8), allocatable :: buf(:,:,:,:,:)   ! (NB,NB,0:3,-ns:ns,1:NK)
  integer :: nfail

  ! deterministic header fixture (values are irrelevant to the invariant
  ! under test; only their round-trip fidelity matters)
  integer, parameter :: h_ndim = 1, h_kgrid(3) = (/2,2,2/), h_nelec = 4, h_natom = 2
  real(8), parameter :: h_h = 0.031d0, h_edir(3) = (/0d0,0d0,1d0/), h_amax = 1.24d0
  real(8), parameter :: h_avec(3,3) = reshape( &
    (/ 10.2d0,0d0,0d0,  0d0,10.2d0,0d0,  0d0,0d0,10.2d0 /), (/3,3/))

  character(len=*), parameter :: f_old = 'vnlk_test_old.bin'
  character(len=*), parameter :: f_new = 'vnlk_test_new.bin'
  character(len=*), parameter :: f_bad = 'vnlk_test_bad.bin'

  nfail = 0
  allocate(buf(NB,NB,0:3,-ns:ns,1:NK))
  call fill_buf(buf)

  call write_old_style(f_old, buf)
  call write_new_style(f_new, buf)

  call check_sizes(nfail)
  call check_bytes_identical(f_old, f_new, nfail)
  call check_negative_control(nfail)
  call check_roundtrip(f_new, buf, nfail)

  call delete_file(f_old)
  call delete_file(f_new)

  if (nfail > 0) then
    write(*,'(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*,'(a)') "All test_vnl_kappa_write checks passed."
  end if

contains

  !======================= fixture ============================================
  subroutine fill_buf(buf)
    implicit none
    complex(8), intent(out) :: buf(NB,NB,0:3,-ns:ns,1:NK)
    integer :: ib, jb, ich, s, ik
    real(8) :: re, im
    do ik = 1, NK
      do s = -ns, ns
        do ich = 0, 3
          do jb = 1, NB
            do ib = 1, NB
              re = sin(0.10d0*ib + 0.20d0*jb + 0.30d0*ich + 0.05d0*s + 0.70d0*ik)
              im = cos(0.15d0*ib - 0.25d0*jb + 0.35d0*ich - 0.06d0*s + 0.50d0*ik)
              buf(ib,jb,ich,s,ik) = dcmplx(re, im)
            end do
          end do
        end do
      end do
    end do
  end subroutine fill_buf

  subroutine write_header(fh)
    implicit none
    integer, intent(in) :: fh
    write(fh) 'SBEVNLK1'
    write(fh) h_ndim, NK, NB, ns
    write(fh) h_kgrid(1), h_kgrid(2), h_kgrid(3)
    write(fh) h_nelec, h_natom
    write(fh) h_h, h_edir(1), h_edir(2), h_edir(3), h_amax
    write(fh) h_avec(1:3,1), h_avec(1:3,2), h_avec(1:3,3)
  end subroutine write_header

  ! ---- OLD pattern: 2*(2ns+1) small writes per k-block ----
  subroutine write_old_style(fname, buf)
    implicit none
    character(*), intent(in) :: fname
    complex(8), intent(in) :: buf(NB,NB,0:3,-ns:ns,1:NK)
    integer :: fh, ikb_s, ikb_e, ik, s
    open(newunit=fh, file=fname, form='unformatted', access='stream', status='replace')
    call write_header(fh)
    do ikb_s = 1, NK, nblk
      ikb_e = min(ikb_s + nblk - 1, NK)
      do ik = ikb_s, ikb_e
        do s = -ns, ns
          write(fh) buf(1:NB,1:NB,0,s,ik)        ! V_s
          write(fh) buf(1:NB,1:NB,1:3,s,ik)      ! W_{1:3,s}
        end do
      end do
    end do
    close(fh)
  end subroutine write_old_style

  ! ---- NEW pattern: 1 write per k-block ----
  subroutine write_new_style(fname, buf)
    implicit none
    character(*), intent(in) :: fname
    complex(8), intent(in) :: buf(NB,NB,0:3,-ns:ns,1:NK)
    integer :: fh, ikb_s, ikb_e
    open(newunit=fh, file=fname, form='unformatted', access='stream', status='replace')
    call write_header(fh)
    do ikb_s = 1, NK, nblk
      ikb_e = min(ikb_s + nblk - 1, NK)
      write(fh) buf(1:NB,1:NB,0:3,-ns:ns,ikb_s:ikb_e)
    end do
    close(fh)
  end subroutine write_new_style

  !======================= checks ==============================================
  subroutine check_sizes(nfail)
    implicit none
    integer, intent(inout) :: nfail
    integer(8) :: fsize_old, fsize_new, expect
    logical :: ok
    inquire(file=f_old, size=fsize_old)
    inquire(file=f_new, size=fsize_new)
    expect = 156_8 + int(NK,8) * int(2*ns+1,8) * 4_8 * int(NB,8)**2 * 16_8
    ok = (fsize_old == expect) .and. (fsize_new == expect)
    write(*,'(a,i0,a,i0,a,i0)') "  file size old/new/expect: ", fsize_old, " / ", fsize_new, &
      & " / ", expect
    call check_true(ok, "file size: old and new both match the reader's size fingerprint", nfail)
  end subroutine check_sizes

  subroutine read_whole_file(fname, bytes)
    implicit none
    character(*), intent(in) :: fname
    integer(1), allocatable, intent(out) :: bytes(:)
    integer :: fh
    integer(8) :: fsize
    inquire(file=fname, size=fsize)
    allocate(bytes(fsize))
    open(newunit=fh, file=fname, form='unformatted', access='stream', status='old', action='read')
    read(fh) bytes
    close(fh)
  end subroutine read_whole_file

  ! returns ndiff = number of differing bytes, ifirst = 1-based offset of the
  ! first difference (0 if none or if sizes differ)
  subroutine compare_bytes(fname_a, fname_b, ndiff, ifirst)
    implicit none
    character(*), intent(in) :: fname_a, fname_b
    integer(8), intent(out) :: ndiff, ifirst
    integer(1), allocatable :: a(:), b(:)
    integer(8) :: i
    call read_whole_file(fname_a, a)
    call read_whole_file(fname_b, b)
    ndiff = 0;  ifirst = 0
    if (size(a) /= size(b)) then
      ndiff = -1;  ifirst = -1
      return
    end if
    do i = 1, size(a, kind=8)
      if (a(i) /= b(i)) then
        ndiff = ndiff + 1
        if (ifirst == 0) ifirst = i
      end if
    end do
  end subroutine compare_bytes

  subroutine check_bytes_identical(fname_a, fname_b, nfail)
    implicit none
    character(*), intent(in) :: fname_a, fname_b
    integer, intent(inout) :: nfail
    integer(8) :: ndiff, ifirst
    call compare_bytes(fname_a, fname_b, ndiff, ifirst)
    write(*,'(a,i0,a,i0)') "  byte compare old vs new: ndiff=", ndiff, "  first_diff_offset=", ifirst
    call check_true(ndiff == 0, &
      & "bytes: single-write-per-block output is byte-for-byte identical to the old nested writes", nfail)
  end subroutine check_bytes_identical

  ! sanity: prove compare_bytes is not vacuously "equal" -- corrupt one byte
  ! well past the header and re-check that a difference IS detected.
  subroutine check_negative_control(nfail)
    implicit none
    integer, intent(inout) :: nfail
    integer(1), allocatable :: bytes(:)
    integer :: fh
    integer(8) :: ndiff, ifirst
    call read_whole_file(f_new, bytes)
    bytes(200) = bytes(200) + 1_1     ! flip a byte inside the body
    open(newunit=fh, file=f_bad, form='unformatted', access='stream', status='replace')
    write(fh) bytes
    close(fh)
    call compare_bytes(f_old, f_bad, ndiff, ifirst)
    write(*,'(a,i0,a,i0)') "  negative control ndiff=", ndiff, "  first_diff_offset=", ifirst
    call check_true(ndiff > 0, &
      & "negative control: byte-compare correctly flags a deliberately corrupted file", nfail)
    call delete_file(f_bad)
  end subroutine check_negative_control

  ! round-trip: read new.bin with the SAME per-(ik,s) sequential access
  ! pattern read_vnl_kappa_data uses (gs_info_ssbe.f90:753-758), and check
  ! every element of buf comes back bit-exact.
  subroutine check_roundtrip(fname, buf, nfail)
    implicit none
    character(*), intent(in) :: fname
    complex(8), intent(in) :: buf(NB,NB,0:3,-ns:ns,1:NK)
    integer, intent(inout) :: nfail
    integer :: fh, ik, s
    integer :: r_ndim, r_nk, r_nb, r_ns, r_kgrid(3), r_nelec, r_natom
    real(8) :: r_h, r_edir(3), r_amax, r_avec(3,3)
    complex(8) :: v0(NB,NB), w13(NB,NB,3)
    logical :: header_ok, body_ok
    real(8) :: maxdiff

    open(newunit=fh, file=fname, form='unformatted', access='stream', status='old', action='read')
    block
      character(8) :: magic
      read(fh) magic
      read(fh) r_ndim, r_nk, r_nb, r_ns
      read(fh) r_kgrid(1), r_kgrid(2), r_kgrid(3)
      read(fh) r_nelec, r_natom
      read(fh) r_h, r_edir(1), r_edir(2), r_edir(3), r_amax
      read(fh) r_avec(1:3,1), r_avec(1:3,2), r_avec(1:3,3)
      header_ok = (magic == 'SBEVNLK1') .and. (r_ndim == h_ndim) .and. (r_nk == NK) .and. &
        & (r_nb == NB) .and. (r_ns == ns) .and. all(r_kgrid == h_kgrid) .and. &
        & (r_nelec == h_nelec) .and. (r_natom == h_natom) .and. &
        & (r_h == h_h) .and. all(r_edir == h_edir) .and. (r_amax == h_amax) .and. &
        & all(r_avec == h_avec)
    end block
    call check_true(header_ok, "round-trip: header fields read back exactly as written", nfail)

    body_ok = .true.
    maxdiff = 0d0
    do ik = 1, NK
      do s = -ns, ns
        read(fh) v0(1:NB,1:NB)
        read(fh) w13(1:NB,1:NB,1:3)
        if (any(v0 /= buf(1:NB,1:NB,0,s,ik))) body_ok = .false.
        if (any(w13 /= buf(1:NB,1:NB,1:3,s,ik))) body_ok = .false.
        maxdiff = max(maxdiff, maxval(abs(v0 - buf(1:NB,1:NB,0,s,ik))))
        maxdiff = max(maxdiff, maxval(abs(w13 - buf(1:NB,1:NB,1:3,s,ik))))
      end do
    end do
    close(fh)
    write(*,'(a,es12.4)') "  round-trip max|read - original| (expect exactly 0): ", maxdiff
    call check_true(body_ok .and. maxdiff == 0d0, &
      & "round-trip: reader-style per-(ik,s) V_s/W reads recover buf bit-exact", nfail)
  end subroutine check_roundtrip

  !======================= helpers ==============================================
  subroutine delete_file(fname)
    implicit none
    character(*), intent(in) :: fname
    integer :: fh
    logical :: ex
    inquire(file=fname, exist=ex)
    if (ex) then
      open(newunit=fh, file=fname, status='old')
      close(fh, status='delete')
    end if
  end subroutine delete_file

  subroutine check_true(cond, label, nfail)
    implicit none
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (cond) then
      write(*,'(a)') "ok    " // label
    else
      write(*,'(a)') "FAIL  " // label
      nfail = nfail + 1
    end if
  end subroutine check_true

end program test_vnl_kappa_write
