! Loads the GW per-state inelastic linewidth file (Si_imsigma.data) into an
! (nb,nk) array in atomic units, honoring the file's "# unit_system =" tag.
! Re-exports the core collision routines so callers need only one use.
module sbe_collision_gw
  use sbe_collision_gw_core, only: add_collision_diag, add_collision_offdiag, &
    & add_collision_vg
  implicit none
  private
  public :: load_gw_rate
  public :: add_collision_diag, add_collision_offdiag, add_collision_vg

  ! a.u. of time in fs: a rate r[1/fs] equals r*au_time_fs in 1/au_time.
  real(8), parameter :: au_time_fs = 0.02418884326505d0

contains

  subroutine load_gw_rate(filename, nb, nk, gamma_nk)
    implicit none
    character(*), intent(in) :: filename
    integer,      intent(in) :: nb, nk
    real(8),      intent(out) :: gamma_nk(nb, nk)
    integer :: fh, ios, in, ik
    integer :: nread, nmissing
    real(8) :: cols(8), rate_fac
    character(512) :: line
    logical :: have(nb, nk)

    gamma_nk = 0.0d0

    have     = .false.
    rate_fac = 1.0d0          ! a.u. fallback for tagless/unknown files
    nread    = 0
    nmissing = 0
    open(newunit=fh, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') "ERROR: cannot open GW rate file: "//trim(filename)
      stop 1
    end if
    do
      read(fh, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') then
        if (index(line, 'unit_system') > 0) then
          block
            integer :: eq_pos
            character(64) :: tok
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
              tok = adjustl(trim(line(eq_pos+1:)))
            else
              tok = ''
            end if
            if      (trim(tok) == 'au' .or. trim(tok) == 'a.u.') then
              rate_fac = 1.0d0
            else if (trim(tok) == 'A_eV_fs') then
              rate_fac = au_time_fs
            else
              rate_fac = 1.0d0
              write(*,'(a)') "# load_gw_rate: WARNING unrecognized unit_system tag '"// &
                trim(tok)//"', assuming a.u."
            end if
          end block
        end if
        cycle
      end if
      read(line, *, iostat=ios) cols(1:8)
      if (ios /= 0) cycle
      in = nint(cols(1)); ik = nint(cols(2))
      if (ik > nk .or. ik < 1) then
        write(*,'(a,2i6)') "ERROR: GW rate file k-index out of range: ", ik, nk
        stop 1
      end if
      if (in < 1 .or. in > nb) cycle      ! band outside window: skip silently
      gamma_nk(in, ik) = cols(7) * rate_fac
      have(in, ik) = .true.
      nread = nread + 1
    end do
    close(fh)
    nmissing = count(.not. have)
    if (nmissing > 0) then
      write(*,'(a,i0,a)') "# load_gw_rate: ", nmissing, &
        & " (n,k) entries absent from rate file -> Gamma=0 (no extra scattering)"
    end if
    write(*,'(a,i0,a,es12.4)') "# load_gw_rate: read ", nread, &
      & " rows; rate->au factor = ", rate_fac
  end subroutine load_gw_rate

end module sbe_collision_gw
