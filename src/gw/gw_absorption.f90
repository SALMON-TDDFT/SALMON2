!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Macroscopic dielectric function / absorption from the full-frequency
! dielectric matrix (spec-b1).  The macroscopic (local-field-included)
! dielectric is the head of the inverse dielectric matrix:
!
!   eps_M(q->0, w) = 1 / [ eps^{-1}(q->0, w) ]_{G=0, G'=0}
!
! and Im eps_M(w) is the RPA absorption spectrum (to be compared with the
! RT-TDDFT delta-kick Im eps; the difference is the xc kernel, RPA vs ALDA).
! epsinv_w comes from calc_w_freq evaluated at the smallest-|q| proxy for
! q->0 (StageA: finite-q head; StageB will use the velocity matrix element).
!
! No proper nouns appear in this file per the project constraint.
!
module gw_absorption_sub
  implicit none
  private
  public :: calc_absorption

  real(8), parameter :: au_to_eV = 27.211386d0

contains

  ! --------------------------------------------------------------------
  ! eps_M(w) = 1 / epsinv_w(ig0,ig0,iw), and (on root) a text dump of the
  ! spectrum to <sysname>_absorption.data with columns
  !   omega[eV]   Re eps_M   Im eps_M
  ! ig0 is the G=0 index in the dielectric G-list; omega_grid is in a.u.
  ! --------------------------------------------------------------------
  subroutine calc_absorption(ng, ig0, nomega, omega_grid, epsinv_w, eps_macro, &
                             is_root, sysname)
    implicit none
    integer,          intent(in)  :: ng, ig0, nomega
    real(8),          intent(in)  :: omega_grid(nomega)      ! a.u.
    complex(8),       intent(in)  :: epsinv_w(ng,ng,nomega)
    complex(8),       intent(out) :: eps_macro(nomega)
    logical,          intent(in)  :: is_root
    character(*),     intent(in)  :: sysname

    integer :: iw, fh
    complex(8) :: zhead

    do iw = 1, nomega
      zhead = epsinv_w(ig0,ig0,iw)
      if (abs(zhead) < 1.0d-30) then
        eps_macro(iw) = (0.0d0, 0.0d0)
      else
        eps_macro(iw) = (1.0d0, 0.0d0) / zhead
      end if
    end do

    if (is_root) then
      open(newunit=fh, file=trim(sysname)//'_absorption.data', status='replace')
      write(fh,'(A)') "# full-frequency G0W0+RPA macroscopic dielectric (finite-q head proxy)"
      write(fh,'(A)') "# omega[eV]            Re eps_M             Im eps_M"
      do iw = 1, nomega
        write(fh,'(3ES22.12)') omega_grid(iw)*au_to_eV, &
                               dble(eps_macro(iw)), aimag(eps_macro(iw))
      end do
      close(fh)
    end if
  end subroutine calc_absorption

end module gw_absorption_sub
