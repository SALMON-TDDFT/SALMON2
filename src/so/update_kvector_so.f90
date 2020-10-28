module update_kvector_so_sub

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: update_kvector_so
  public :: SPIN_ORBIT_ON

contains

  subroutine update_kvector_so( ppg, kAc, ik_s, ik_e )
    use structures, only: s_pp_grid
    implicit none
    type(s_pp_grid)    :: ppg
    integer,intent(in) :: ik_s,ik_e
    real(8),intent(in) :: kAc(3,ik_s:ik_e)
    !
    integer :: ilma,iatom,j,ik,Nlma,jangl,ispin
    real(8) :: x,y,z,k(3),kr
    complex(8) :: conjg_ekr

    Nlma = size( ppg%ia_tbl_so )

    if ( .not.allocated(ppg%zekr_uV_so) ) then
      allocate(ppg%zekr_uV_so(ppg%nps,Nlma,ik_s:ik_e,2,1))
    end if

!$omp parallel do collapse(2) default(shared) private(ik,k,ilma,iatom,ispin,j,x,y,z,kr,conjg_ekr)
    do ik=ik_s,ik_e
      do ilma=1,Nlma

        k(:) = kAc(:,ik)
        iatom = ppg%ia_tbl_so(ilma)

        do ispin=1,2

          do j=1,ppg%mps(iatom)
            x = ppg%rxyz(1,j,iatom)
            y = ppg%rxyz(2,j,iatom)
            z = ppg%rxyz(3,j,iatom)
            kr = k(1)*x + k(2)*y + k(3)*z
            conjg_ekr = dcmplx( cos(kr),-sin(kr) )
            ppg%zekr_uV_so(j,ilma,ik,ispin,1) = conjg_ekr * ppg%uv_so(j,ilma,ispin,1)
          end do

        end do !ispin

      end do !ilma
    end do !ik

  end subroutine update_kvector_so

end module update_kvector_so_sub
