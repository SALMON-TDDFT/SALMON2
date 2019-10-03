module update_kvector_plusu_sub

  use plusU_global, only: PLUS_U_ON

  implicit none

  private
  public :: update_kvector_plusU
  public :: PLUS_U_ON

contains

  subroutine update_kvector_plusU( ppg, kAc, ik_s, ik_e )
    use structures, only: s_pp_grid
    implicit none
    type(s_pp_grid)    :: ppg
    integer,intent(in) :: ik_s,ik_e
    real(8),intent(in) :: kAc(3,ik_s:ik_e)
    !
    integer :: ilma,iatom,j,ik,Nlma,jangl,ispin
    real(8) :: x,y,z,k(3),kr
    complex(8) :: conjg_ekr

    Nlma = size( ppg%ia_tbl_ao )

    if ( .not.allocated(ppg%zekr_phi_ao) ) then
      allocate( ppg%zekr_phi_ao(ppg%nps_ao,Nlma,ik_s:ik_e) )
      ppg%zekr_phi_ao=(0.0d0,0.0d0)
    end if

    do ik=ik_s,ik_e

      k(:) = kAc(:,ik)

      do ilma=1,Nlma

        iatom = ppg%ia_tbl_ao(ilma)

        do j=1,ppg%mps(iatom)
          x = ppg%rxyz(1,j,iatom)
          y = ppg%rxyz(2,j,iatom)
          z = ppg%rxyz(3,j,iatom)
          kr = k(1)*x + k(2)*y + k(3)*z
          conjg_ekr = dcmplx( cos(kr),-sin(kr) )
          ppg%zekr_phi_ao(j,ilma,ik) = conjg_ekr * ppg%phi_ao(j,ilma)
        end do

      end do !ilma

    end do !ik

  end subroutine update_kvector_plusU

end module update_kvector_plusu_sub
