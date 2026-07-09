! mpiio_export_ssbe.f90
!
! Collective MPI-IO helper module for the SBE/gicov GS-export I/O paths
! (vnl_kappa / prod_dk / tm.data). Provides generic collective
! write_at_all/read_at_all wrappers, canonical byte-offset arithmetic,
! independent (rank-0) header field writers/readers, a splitmix64-based
! order-independent checksum, a finiteness check, and a fail-closed abort.
!
! #ifdef USE_MPI builds the real MPI-IO implementation (use mpi,
! MPI_File_open/set_view-free at-offset calls). #ifndef USE_MPI builds a
! serial access='stream' fallback with the same disp/nz semantics, so the
! module type-checks under a local (non-MPI) gfortran CI build.
!
! Format constants (magic strings, header sizes, element size) are the
! single source of truth for all SBE/gicov export readers/writers
! (Tasks 2-9 of the sbe-mpiio-export plan consume this module).

#include "config.h"

module mpiio_export_ssbe
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: mpiio_open_write, mpiio_open_read, mpiio_close
  public :: mpiio_write_at_all_z, mpiio_read_at_all_z
  public :: mpiio_wr_c8, mpiio_wr_i, mpiio_wr_d, mpiio_wr_i8
  public :: mpiio_rd_c8, mpiio_rd_i, mpiio_rd_d, mpiio_rd_i8
  public :: mpiio_disp_k, mpiio_csum_local_z, mpiio_csum_reduce, mpiio_finite_z, mpiio_abort
  public :: mpiio_validate_z, mpiio_assert_dim
  public :: MAGIC_VNL, MAGIC_PRODK, MAGIC_TM, MAGIC_CHK
  public :: HDR_VNL, HDR_PRODK, HDR_TM, SBE_FMT_VERSION, SZ_Z

#ifndef USE_MPI
  ! serial fallback so the type name resolves the same way for callers
  ! (kept identical in byte-width to a real MPI file-offset kind)
  integer, parameter :: MPI_OFFSET_KIND = int64
#endif

  ! ---- Format constants (single source of truth; verbatim from the spec) ----
  ! magic strings (all character(len=8), byte-identical on disk)
  character(8), parameter :: MAGIC_VNL   = 'SBEVNLK1'   ! existing (write.f90:1677) -- DO NOT CHANGE
  character(8), parameter :: MAGIC_PRODK = 'SBEPRDK1'   ! NEW binary prod_dk
  character(8), parameter :: MAGIC_TM    = 'SBETMD01'   ! NEW binary tm.data
  character(8), parameter :: MAGIC_CHK   = 'SBEVNCHK'   ! vnl sidecar checksum file

  ! header sizes (bytes)
  integer(MPI_OFFSET_KIND), parameter :: HDR_VNL   = 156_MPI_OFFSET_KIND  ! existing, immutable
  integer(MPI_OFFSET_KIND), parameter :: HDR_PRODK = 36_MPI_OFFSET_KIND  ! magic8+ver/nk/no/ndk/nbvec_full(5 i32)+csum(i64)
  integer(MPI_OFFSET_KIND), parameter :: HDR_TM    = 32_MPI_OFFSET_KIND  ! magic8+ver/nk/nb/nblocks(4 i32)+csum(i64)
  integer, parameter :: SBE_FMT_VERSION = 1
  integer, parameter :: SZ_Z = 16   ! bytes per complex(8)

contains

  ! =====================================================================
  ! open / close
  ! =====================================================================

  subroutine mpiio_open_write(path, comm, fh, ierr)
    implicit none
    character(*), intent(in)  :: path
    integer,      intent(in)  :: comm
    integer,      intent(out) :: fh
    integer,      intent(out) :: ierr
#ifdef USE_MPI
    integer :: minfo
    call mpiio_set_info_write(minfo)
    call MPI_File_open(comm, path, ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), minfo, fh, ierr)
    if (ierr /= MPI_SUCCESS) call mpiio_abort('mpiio_open_write: MPI_File_open failed: '//trim(path), comm)
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND, ierr)
    if (ierr /= MPI_SUCCESS) call mpiio_abort('mpiio_open_write: MPI_File_set_size failed: '//trim(path), comm)
    call MPI_Info_free(minfo, ierr)
#else
    open(newunit=fh, file=trim(path), form='unformatted', access='stream', &
      &  status='replace', action='write', iostat=ierr)
    if (ierr /= 0) call mpiio_abort('mpiio_open_write: open failed: '//trim(path), comm)
#endif
  end subroutine mpiio_open_write

  subroutine mpiio_open_read(path, comm, fh, ierr)
    implicit none
    character(*), intent(in)  :: path
    integer,      intent(in)  :: comm
    integer,      intent(out) :: fh
    integer,      intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if (ierr /= MPI_SUCCESS) call mpiio_abort('mpiio_open_read: MPI_File_open failed: '//trim(path), comm)
#else
    open(newunit=fh, file=trim(path), form='unformatted', access='stream', &
      &  status='old', action='read', iostat=ierr)
    if (ierr /= 0) call mpiio_abort('mpiio_open_read: open failed: '//trim(path), comm)
#endif
  end subroutine mpiio_open_read

  subroutine mpiio_close(fh, ierr)
    implicit none
    integer, intent(inout) :: fh
    integer, intent(out)   :: ierr
#ifdef USE_MPI
    call MPI_File_close(fh, ierr)
#else
    close(fh, iostat=ierr)
#endif
  end subroutine mpiio_close

#ifdef USE_MPI
  ! Advisory MPI_Info hints for the write-open path. Kept minimal here
  ! (spread the shared file across all FEFS OSTs); further collective-
  ! buffering / striping-unit tuning is added by Task 10.
  subroutine mpiio_set_info_write(minfo)
    implicit none
    integer, intent(out) :: minfo
    integer :: ierr
    call MPI_Info_create(minfo, ierr)
    call MPI_Info_set(minfo, 'striping_factor', '-1', ierr)
  end subroutine mpiio_set_info_write
#endif

  ! =====================================================================
  ! collective bulk complex(8) I/O
  ! =====================================================================

  subroutine mpiio_write_at_all_z(fh, disp, buf, nz, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    complex(8),                intent(in)  :: buf(*)
    integer(8),                intent(in)  :: nz
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_write_at_all(fh, disp, buf, int(nz), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
#else
    ierr = 0
    if (nz > 0_8) then
      write(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) buf(1:nz)
    end if
#endif
  end subroutine mpiio_write_at_all_z

  subroutine mpiio_read_at_all_z(fh, disp, buf, nz, comm, tag, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    complex(8),                intent(out) :: buf(*)
    integer(8),                intent(in)  :: nz
    integer,                   intent(in)  :: comm
    character(*),               intent(in)  :: tag
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    integer :: st(MPI_STATUS_SIZE), ngot
    call MPI_File_read_at_all(fh, disp, buf, int(nz), MPI_DOUBLE_COMPLEX, st, ierr)
    if (ierr /= MPI_SUCCESS) call mpiio_abort(trim(tag)//': MPI_File_read_at_all failed', comm)
    call MPI_Get_count(st, MPI_DOUBLE_COMPLEX, ngot, ierr)
    if (ngot /= int(nz)) call mpiio_abort(trim(tag)//': short read', comm)
#else
    ierr = 0
    if (nz > 0_8) then
      read(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) buf(1:nz)
      if (ierr /= 0) call mpiio_abort(trim(tag)//': short read', comm)
    end if
#endif
  end subroutine mpiio_read_at_all_z

  ! =====================================================================
  ! independent (rank-0-only, by caller convention) header field I/O
  ! =====================================================================

  subroutine mpiio_wr_c8(fh, disp, str, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    character(*),               intent(in)  :: str
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_write_at(fh, disp, str, len(str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
#else
    write(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) str
#endif
  end subroutine mpiio_wr_c8

  subroutine mpiio_rd_c8(fh, disp, str, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    character(*),               intent(out) :: str
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_read_at(fh, disp, str, len(str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
#else
    read(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) str
#endif
  end subroutine mpiio_rd_c8

  subroutine mpiio_wr_i(fh, disp, iarr, n, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer,                   intent(in)  :: n
    integer,                   intent(in)  :: iarr(n)
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_write_at(fh, disp, iarr, n, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
    write(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) iarr(1:n)
#endif
  end subroutine mpiio_wr_i

  subroutine mpiio_rd_i(fh, disp, iarr, n, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer,                   intent(in)  :: n
    integer,                   intent(out) :: iarr(n)
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_read_at(fh, disp, iarr, n, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
    read(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) iarr(1:n)
#endif
  end subroutine mpiio_rd_i

  subroutine mpiio_wr_d(fh, disp, darr, n, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer,                   intent(in)  :: n
    real(8),                   intent(in)  :: darr(n)
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_write_at(fh, disp, darr, n, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
#else
    write(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) darr(1:n)
#endif
  end subroutine mpiio_wr_d

  subroutine mpiio_rd_d(fh, disp, darr, n, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer,                   intent(in)  :: n
    real(8),                   intent(out) :: darr(n)
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_read_at(fh, disp, darr, n, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
#else
    read(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) darr(1:n)
#endif
  end subroutine mpiio_rd_d

  subroutine mpiio_wr_i8(fh, disp, i8val, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer(8),                intent(in)  :: i8val
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_write_at(fh, disp, i8val, 1, MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
#else
    write(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) i8val
#endif
  end subroutine mpiio_wr_i8

  subroutine mpiio_rd_i8(fh, disp, i8val, ierr)
    implicit none
    integer,                   intent(in)  :: fh
    integer(MPI_OFFSET_KIND),  intent(in)  :: disp
    integer(8),                intent(out) :: i8val
    integer,                   intent(out) :: ierr
#ifdef USE_MPI
    call MPI_File_read_at(fh, disp, i8val, 1, MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
#else
    read(fh, pos=disp+1_MPI_OFFSET_KIND, iostat=ierr) i8val
#endif
  end subroutine mpiio_rd_i8

  ! =====================================================================
  ! offset arithmetic
  ! =====================================================================

  ! Per-k-chunk byte displacement: hdr + SZ_Z * nz_k * (ik-1), all
  ! arithmetic carried in integer(8) then cast to MPI_OFFSET_KIND.
  pure function mpiio_disp_k(hdr, nz_k, ik) result(disp)
    implicit none
    integer(MPI_OFFSET_KIND), intent(in) :: hdr
    integer(8),                intent(in) :: nz_k
    integer,                   intent(in) :: ik
    integer(MPI_OFFSET_KIND) :: disp
    disp = hdr + int(SZ_Z, 8) * nz_k * (int(ik, 8) - 1_8)
  end function mpiio_disp_k

  ! =====================================================================
  ! splitmix64 checksum (order-independent XOR-fold, index-sensitive)
  ! =====================================================================

  pure function sm64(x) result(z)
    implicit none
    integer(int64), intent(in) :: x
    integer(int64) :: z
    z = x + int(z'9E3779B97F4A7C15', int64)
    z = ieor(z, ishft(z,-30)) * int(z'BF58476D1CE4E5B9', int64)
    z = ieor(z, ishft(z,-27)) * int(z'94D049BB133111EB', int64)
    z = ieor(z, ishft(z,-31))
  end function sm64

  ! Local (per-rank) checksum over nz complex(8) elements whose global
  ! 0-based start index is g0. XOR-fold so per-rank partials combine via
  ! MPI_BXOR regardless of rank count or ordering; the 2*g / 2*g+1 index
  ! term makes it sensitive to truncation/reordering of elements.
  function mpiio_csum_local_z(buf, nz, g0) result(c)
    implicit none
    complex(8), intent(in) :: buf(*)
    integer(8), intent(in) :: nz, g0
    integer(8) :: c
    integer(8) :: i, g
    c = 0_8
    do i = 1_8, nz
      g = g0 + (i - 1_8)
      c = ieor(c, sm64(ieor(transfer(real(buf(i)),  0_int64), 2_8*g)))
      c = ieor(c, sm64(ieor(transfer(aimag(buf(i)), 0_int64), 2_8*g + 1_8)))
    end do
  end function mpiio_csum_local_z

  function mpiio_csum_reduce(clocal, comm) result(cglob)
    implicit none
    integer(8), intent(in) :: clocal
    integer,    intent(in) :: comm
    integer(8) :: cglob
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(clocal, cglob, 1, MPI_INTEGER8, MPI_BXOR, comm, ierr)
    if (ierr /= MPI_SUCCESS) call mpiio_abort('mpiio_csum_reduce: MPI_Allreduce failed', comm)
#else
    cglob = clocal
#endif
  end function mpiio_csum_reduce

  ! =====================================================================
  ! finiteness check
  ! =====================================================================

  function mpiio_finite_z(buf, nz) result(ok)
    implicit none
    complex(8), intent(in) :: buf(*)
    integer(8), intent(in) :: nz
    logical :: ok
    integer(8) :: i
    ok = .true.
    do i = 1_8, nz
      if (.not. (ieee_is_finite(real(buf(i))) .and. ieee_is_finite(aimag(buf(i))))) then
        ok = .false.
        return
      end if
    end do
  end function mpiio_finite_z

  ! =====================================================================
  ! fail-closed validation layer (Task 2)
  !
  ! Combines the finiteness check and the order-independent splitmix64
  ! checksum into one call for readers: layer (3) finiteness is checked
  ! first (a NaN/Inf abort names itself distinctly from a checksum
  ! mismatch), then layer (4) checksum. Layer (1) short-read is already
  ! enforced inside mpiio_read_at_all_z (MPI_Get_count, Task 1); layer (2)
  ! dimension mismatch is mpiio_assert_dim below. Both routines terminate
  ! via mpiio_abort on failure (fail-closed: never returns past a failed
  ! check) and are safe to call from any communicator (comm is forwarded
  ! to mpiio_abort for MPI_Abort scope).
  ! =====================================================================

  subroutine mpiio_validate_z(buf, nz, g0, comm, expect_csum, tag, ierr)
    implicit none
    complex(8),   intent(in)  :: buf(*)
    integer(8),   intent(in)  :: nz
    integer(8),   intent(in)  :: g0
    integer,      intent(in)  :: comm
    integer(8),   intent(in)  :: expect_csum
    character(*), intent(in)  :: tag
    integer,      intent(out) :: ierr
    integer(8) :: csum_local, csum_glob

    ierr = 0

    if (.not. mpiio_finite_z(buf, nz)) then
      call mpiio_abort(trim(tag)//': non-finite element (NaN/Inf)', comm)
    end if

    csum_local = mpiio_csum_local_z(buf, nz, g0)
    csum_glob  = mpiio_csum_reduce(csum_local, comm)
    if (csum_glob /= expect_csum) then
      call mpiio_abort(trim(tag)//': checksum mismatch', comm)
    end if
  end subroutine mpiio_validate_z

  subroutine mpiio_assert_dim(got, want, name, comm)
    implicit none
    integer,      intent(in) :: got
    integer,      intent(in) :: want
    character(*), intent(in) :: name
    integer,      intent(in) :: comm
    character(len=200) :: msg

    if (got /= want) then
      write(msg, '(a,": file ",i0," vs expected ",i0)') trim(name), got, want
      call mpiio_abort(trim(msg), comm)
    end if
  end subroutine mpiio_assert_dim

  ! =====================================================================
  ! fail-closed abort (never fail-open; no dynamic error-stop code)
  ! =====================================================================

  subroutine mpiio_abort(msg, comm)
    implicit none
    character(*), intent(in) :: msg
    integer,      intent(in) :: comm
#ifdef USE_MPI
    integer :: ierr
    print '(a)', 'ERROR(mpiio): '//trim(msg)
    call MPI_Abort(comm, 1, ierr)
    error stop
#else
    print '(a)', 'ERROR(mpiio): '//trim(msg)
    error stop
#endif
  end subroutine mpiio_abort

end module mpiio_export_ssbe
