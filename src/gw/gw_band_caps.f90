!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module gw_band_caps_sub
  implicit none
  public :: resolve_band_caps

contains

  ! Resolve the 0=all convention for the chi0 (eps) and Sigma_c band caps and
  ! validate them against the occupied band count.  The two caps are INDEPENDENT
  ! (neither bounds the other): the dielectric RPA transition count and the
  ! self-energy intermediate-state count are separate physical sums.  A cap below
  ! the occupied count is a physical error (the static remainder corrects the
  ! high-energy empty-state tail, NOT dropped occupied states).
  subroutine resolve_band_caps(no, nocc, neps_in, nsig_in, neps, nsig)
    integer, intent(in)  :: no       ! total bands available (system%no)
    integer, intent(in)  :: nocc     ! highest occupied band index
    integer, intent(in)  :: neps_in  ! nband_eps input   (0 = all)
    integer, intent(in)  :: nsig_in  ! nband_sigma input (0 = all)
    integer, intent(out) :: neps     ! resolved chi0/eps unoccupied cap
    integer, intent(out) :: nsig     ! resolved Sigma_c intermediate cap

    neps = no
    if (neps_in > 0) neps = min(neps_in, no)
    nsig = no
    if (nsig_in > 0) nsig = min(nsig_in, no)

    if (neps < nocc) stop 'gw: nband_eps must be >= occupied band count'
    if (nsig < nocc) stop 'gw: nband_sigma must be >= occupied band count'
  end subroutine resolve_band_caps

end module gw_band_caps_sub
