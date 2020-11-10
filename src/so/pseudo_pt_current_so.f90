!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

module pseudo_pt_current_so

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none
  private
  public :: SPIN_ORBIT_ON
  public :: calc_current_nonlocal_so
  public :: calc_current_nonlocal_rdivided_so

contains

  subroutine calc_current_nonlocal_so( jw, psi, ppg, is_array, ie_array,ik )
    use structures, only: s_pp_grid
    implicit none
    real(8)   ,intent(out) :: jw(3)
    integer   ,intent(in)  :: is_array(3),ie_array(3),ik
    complex(8),intent(in)  :: psi(is_array(1):ie_array(1),   &
                                  is_array(2):ie_array(2),   &
                                  is_array(3):ie_array(3),2)
    type(s_pp_grid),intent(in) :: ppg

    real(8) :: x,y,z
    integer :: Nlma,ilma,ia,j,jx,jy,jz,ispin
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: wrk,wrk0,wrk1,wrk2,wrk3
    complex(8) :: uVpsi0(2),uVpsiX(2),uVpsiY(2),uVpsiZ(2)

    Nlma = size(ppg%ia_tbl_so)

    jw = 0.0d0

!$omp parallel do default(shared) private(ilma,ia,ispin,wrk0,wrk1,wrk2,wrk3,j,x,y,z,jx,jy,jz,wrk) &
!$omp private(uVpsi0,uVpsiX,uVpsiY,uVpsiZ) reduction(+:jw)
    do ilma=1,Nlma

       ia=ppg%ia_tbl_so(ilma)

       do ispin=1,2

          wrk0=zero
          wrk1=zero
          wrk2=zero
          wrk3=zero
          do j=1,ppg%mps(ia)
             x = ppg%rxyz(1,j,ia)
             y = ppg%rxyz(2,j,ia)
             z = ppg%rxyz(3,j,ia)
             jx = ppg%jxyz(1,j,ia)
             jy = ppg%jxyz(2,j,ia)
             jz = ppg%jxyz(3,j,ia)
             wrk  = conjg(ppg%zekr_uV_so(j,ilma,ik,ispin,1))*psi(jx,jy,jz,ispin)
             wrk0 = wrk0 + wrk 
             wrk1 = wrk1 + x*wrk
             wrk2 = wrk2 + y*wrk
             wrk3 = wrk3 + z*wrk
          end do
          uVpsi0(ispin) = wrk0 * ppg%rinv_uvu_so(ilma)
          uVpsiX(ispin) = wrk1
          uVpsiY(ispin) = wrk2
          uVpsiZ(ispin) = wrk3

       end do !ispin

       wrk1 = conjg( uVpsiX(1) )*uVpsi0(1) - uVpsiX(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiX(1) )*uVpsi0(2) - uVpsiX(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiX(2) )*uVpsi0(1) - uVpsiX(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiX(2) )*uVpsi0(2) - uVpsiX(2)*conjg( uVpsi0(2) )
       wrk2 = conjg( uVpsiY(1) )*uVpsi0(1) - uVpsiY(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiY(1) )*uVpsi0(2) - uVpsiY(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiY(2) )*uVpsi0(1) - uVpsiY(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiY(2) )*uVpsi0(2) - uVpsiY(2)*conjg( uVpsi0(2) )
       wrk3 = conjg( uVpsiZ(1) )*uVpsi0(1) - uVpsiZ(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiZ(1) )*uVpsi0(2) - uVpsiZ(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiZ(2) )*uVpsi0(1) - uVpsiZ(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiZ(2) )*uVpsi0(2) - uVpsiZ(2)*conjg( uVpsi0(2) )
       jw(1) = jw(1) + aimag(wrk1)
       jw(2) = jw(2) + aimag(wrk2)
       jw(3) = jw(3) + aimag(wrk3)

    end do !ilma

    jw=jw*0.5d0

    return
  end subroutine calc_current_nonlocal_so


  subroutine calc_current_nonlocal_rdivided_so( jw,psi,ppg,is_array,ie_array,ik,icomm_r )
    use structures, only: s_pp_grid
    use communication, only: comm_summation
    use timer
    implicit none
    real(8)   ,intent(out) :: jw(3)
    integer   ,intent(in)  :: is_array(3),ie_array(3),ik,icomm_r
    complex(8),intent(in)  :: psi(is_array(1):ie_array(1),   &
                                  is_array(2):ie_array(2),   &
                                  is_array(3):ie_array(3),2)
    type(s_pp_grid),intent(in) :: ppg

    real(8) :: x,y,z
    integer :: Nlma,ilma,ia,j,jx,jy,jz,ispin
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: wrk,wrk0,wrk1,wrk2,wrk3
    complex(8) :: uVpsi0(2),uVpsiX(2),uVpsiY(2),uVpsiZ(2)
    complex(8),allocatable :: uVpsibox1(:,:), uVpsibox2(:,:)

    Nlma = size(ppg%ia_tbl_so)

    call timer_begin(LOG_CURRENT_CALC_UVPSI_RDIVIDED)

    allocate( uVpsibox1(2,Nlma) ); uVpsibox1=zero
    allocate( uVpsibox2(2,Nlma) ); uVpsibox2=zero

!$omp parallel do default(shared) private(ilma,ia,ispin,wrk0,j,x,y,z,jx,jy,jz)
    do ilma=1,Nlma
       ia=ppg%ia_tbl_so(ilma)
       do ispin=1,2
          wrk0=zero
          do j=1,ppg%mps(ia)
             x = ppg%rxyz(1,j,ia)
             y = ppg%rxyz(2,j,ia)
             z = ppg%rxyz(3,j,ia)
             jx = ppg%jxyz(1,j,ia)
             jy = ppg%jxyz(2,j,ia)
             jz = ppg%jxyz(3,j,ia)
             wrk0 = wrk0 + conjg(ppg%zekr_uV_so(j,ilma,ik,ispin,1))*psi(jx,jy,jz,ispin)
          end do
          uVpsibox1(ispin,ilma) = wrk0 * ppg%rinv_uvu_so(ilma)
       end do !ispin
    end do !ilma

    call comm_summation(uVpsibox1,uVpsibox2,size(uVpsibox2),icomm_r)

    call timer_end(LOG_CURRENT_CALC_UVPSI_RDIVIDED)

    jw = 0.0d0
!$omp parallel do default(shared) private(ilma,ia,ispin,wrk1,wrk2,wrk3,j,x,y,z,jx,jy,jz,wrk) &
!$omp private(uVpsi0,uVpsiX,uVpsiY,uVpsiZ) reduction(+:jw)
    do ilma=1,Nlma

       ia=ppg%ia_tbl_so(ilma)

       do ispin=1,2

          wrk1=zero
          wrk2=zero
          wrk3=zero
          do j=1,ppg%mps(ia)
             x = ppg%rxyz(1,j,ia)
             y = ppg%rxyz(2,j,ia)
             z = ppg%rxyz(3,j,ia)
             jx = ppg%jxyz(1,j,ia)
             jy = ppg%jxyz(2,j,ia)
             jz = ppg%jxyz(3,j,ia)
             wrk  = conjg(ppg%zekr_uV_so(j,ilma,ik,ispin,1))*psi(jx,jy,jz,ispin)
             wrk1 = wrk1 + x*wrk
             wrk2 = wrk2 + y*wrk
             wrk3 = wrk3 + z*wrk
          end do
          uVpsi0(ispin) = uVpsibox2(ispin,ilma)
          uVpsiX(ispin) = wrk1
          uVpsiY(ispin) = wrk2
          uVpsiZ(ispin) = wrk3

       end do !ispin

       wrk1 = conjg( uVpsiX(1) )*uVpsi0(1) - uVpsiX(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiX(1) )*uVpsi0(2) - uVpsiX(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiX(2) )*uVpsi0(1) - uVpsiX(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiX(2) )*uVpsi0(2) - uVpsiX(2)*conjg( uVpsi0(2) )
       wrk2 = conjg( uVpsiY(1) )*uVpsi0(1) - uVpsiY(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiY(1) )*uVpsi0(2) - uVpsiY(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiY(2) )*uVpsi0(1) - uVpsiY(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiY(2) )*uVpsi0(2) - uVpsiY(2)*conjg( uVpsi0(2) )
       wrk3 = conjg( uVpsiZ(1) )*uVpsi0(1) - uVpsiZ(1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiZ(1) )*uVpsi0(2) - uVpsiZ(1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsiZ(2) )*uVpsi0(1) - uVpsiZ(2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsiZ(2) )*uVpsi0(2) - uVpsiZ(2)*conjg( uVpsi0(2) )
       jw(1) = jw(1) + aimag(wrk1)
       jw(2) = jw(2) + aimag(wrk2)
       jw(3) = jw(3) + aimag(wrk3)

    end do !ilma

    jw=jw*0.5d0

    deallocate( uVpsibox2 )
    deallocate( uVpsibox1 )

    return
  end subroutine calc_current_nonlocal_rdivided_so

end module pseudo_pt_current_so
