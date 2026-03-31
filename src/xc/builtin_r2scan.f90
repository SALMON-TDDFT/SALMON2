module builtin_r2scan

  implicit none

contains

  subroutine exc_cor_r2scan(nl, rho, rho_s, grho_s, lrho_s, tau_s, eexc, vexc)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho(nl)
    real(8), intent(in) :: rho_s(nl)
    real(8), intent(in) :: grho_s(nl,3)
    real(8), intent(in) :: lrho_s(nl)
    real(8), intent(in) :: tau_s(nl)
    real(8), intent(out) :: eexc(nl)
    real(8), intent(out) :: vexc(nl)

    eexc = 0d0
    vexc = 0d0

    stop "builtin r2SCAN kernel not implemented"
  end subroutine exc_cor_r2scan

end module builtin_r2scan
