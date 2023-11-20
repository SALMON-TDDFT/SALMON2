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

  implicit none
  private
  public :: calc_current_nonlocal_so
  public :: calc_current_nonlocal_rdivided_so
  public :: calc_spin_current_nonlocal

contains

  subroutine calc_current_nonlocal_so( jw, psi, ppg, is_array, ie_array,ik )
    !$acc routine worker
    use structures, only: s_pp_grid
    implicit none
    real(8)   ,intent(out) :: jw(3)
    integer   ,intent(in)  :: is_array(3),ie_array(3),ik
    complex(8),intent(in)  :: psi(is_array(1):ie_array(1),   &
                                  is_array(2):ie_array(2),   &
                                  is_array(3):ie_array(3),2)
    type(s_pp_grid),intent(in) :: ppg

#ifdef USE_OPENACC
    real(8)               :: jw_1, jw_2, jw_3
#endif
    real(8) :: x,y,z
    integer :: Nlma,ilma,ia,j,jx,jy,jz,ispin
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: wrk1,wrk2,wrk3
    complex(8) :: uVpsi0(2),uVpsi_r(3,2)

    Nlma = size(ppg%ia_tbl_so)

    jw = 0.0d0

#ifdef USE_OPENACC
    jw_1 = 0d0
    jw_2 = 0d0
    jw_3 = 0d0
!$acc loop worker private(ilma,ia,wrk1,wrk2,wrk3,j,x,y,z,jx,jy,jz) &
!$acc private(uVpsi0,uVpsi_r) reduction(+:jw_1, jw_2, jw_3)
#else
!$omp parallel do default(shared) private(ilma,ia,wrk1,wrk2,wrk3,j,x,y,z,jx,jy,jz) &
!$omp private(uVpsi0,uVpsi_r) reduction(+:jw)
#endif
    do ilma=1,Nlma

       ia=ppg%ia_tbl_so(ilma)

       uVpsi0(:)=zero
       uVpsi_r(:,:)=zero

       ! To avoid the compiler problem of GPGPU, cases of ispin=1 and ispin=2 are written redundantly.
       ! start 
       do j=1,ppg%mps(ia)
          x = ppg%rxyz(1,j,ia)
          y = ppg%rxyz(2,j,ia)
          z = ppg%rxyz(3,j,ia)
          jx = ppg%jxyz(1,j,ia)
          jy = ppg%jxyz(2,j,ia)
          jz = ppg%jxyz(3,j,ia)
          uVpsi0(1) = uVpsi0(1) + conjg(ppg%zekr_uV_so(j,ilma,ik,1,1))*psi(jx,jy,jz,1)
          uVpsi_r(1,1) = uVpsi_r(1,1) + conjg(ppg%zekr_uV_so(j,ilma,ik,1,1))*x*psi(jx,jy,jz,1)
          uVpsi_r(2,1) = uVpsi_r(2,1) + conjg(ppg%zekr_uV_so(j,ilma,ik,1,1))*y*psi(jx,jy,jz,1)
          uVpsi_r(3,1) = uVpsi_r(3,1) + conjg(ppg%zekr_uV_so(j,ilma,ik,1,1))*z*psi(jx,jy,jz,1)
       end do
       do j=1,ppg%mps(ia)
          x = ppg%rxyz(1,j,ia)
          y = ppg%rxyz(2,j,ia)
          z = ppg%rxyz(3,j,ia)
          jx = ppg%jxyz(1,j,ia)
          jy = ppg%jxyz(2,j,ia)
          jz = ppg%jxyz(3,j,ia)
          uVpsi0(2) = uVpsi0(2) + conjg(ppg%zekr_uV_so(j,ilma,ik,2,1))*psi(jx,jy,jz,2)
          uVpsi_r(1,2) = uVpsi_r(1,2) + conjg(ppg%zekr_uV_so(j,ilma,ik,2,1))*x*psi(jx,jy,jz,2)
          uVpsi_r(2,2) = uVpsi_r(2,2) + conjg(ppg%zekr_uV_so(j,ilma,ik,2,1))*y*psi(jx,jy,jz,2)
          uVpsi_r(3,2) = uVpsi_r(3,2) + conjg(ppg%zekr_uV_so(j,ilma,ik,2,1))*z*psi(jx,jy,jz,2)
       end do
       uVpsi0(1) = uVpsi0(1) * ppg%rinv_uvu_so(ilma)
       uVpsi0(2) = uVpsi0(2) * ppg%rinv_uvu_so(ilma)
       ! end

       wrk1 = conjg( uVpsi_r(1,1) )*uVpsi0(1) - uVpsi_r(1,1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(1,1) )*uVpsi0(2) - uVpsi_r(1,1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsi_r(1,2) )*uVpsi0(1) - uVpsi_r(1,2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(1,2) )*uVpsi0(2) - uVpsi_r(1,2)*conjg( uVpsi0(2) )
       wrk2 = conjg( uVpsi_r(2,1) )*uVpsi0(1) - uVpsi_r(2,1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(2,1) )*uVpsi0(2) - uVpsi_r(2,1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsi_r(2,2) )*uVpsi0(1) - uVpsi_r(2,2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(2,2) )*uVpsi0(2) - uVpsi_r(2,2)*conjg( uVpsi0(2) )
       wrk3 = conjg( uVpsi_r(3,1) )*uVpsi0(1) - uVpsi_r(3,1)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(3,1) )*uVpsi0(2) - uVpsi_r(3,1)*conjg( uVpsi0(2) ) &
            + conjg( uVpsi_r(3,2) )*uVpsi0(1) - uVpsi_r(3,2)*conjg( uVpsi0(1) ) &
            + conjg( uVpsi_r(3,2) )*uVpsi0(2) - uVpsi_r(3,2)*conjg( uVpsi0(2) )
#ifdef USE_OPENACC
       jw_1 = jw_1 + aimag(wrk1)
       jw_2 = jw_2 + aimag(wrk2)
       jw_3 = jw_3 + aimag(wrk3)
#else
       jw(1) = jw(1) + aimag(wrk1)
       jw(2) = jw(2) + aimag(wrk2)
       jw(3) = jw(3) + aimag(wrk3)
#endif

    end do !ilma

#ifdef USE_OPENACC
    jw(1) = jw_1
    jw(2) = jw_2
    jw(3) = jw_3
#endif

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
  
  
  subroutine calc_spin_current_nonlocal( jspin, psi, ppg, is_array, ie_array, ik )
    !$acc routine worker
    use structures, only: s_pp_grid
    implicit none
    real(8)   ,intent(out) :: jspin(3,0:3)
    integer   ,intent(in)  :: is_array(3),ie_array(3),ik
    complex(8),intent(in)  :: psi(is_array(1):ie_array(1),   &
                                  is_array(2):ie_array(2),   &
                                  is_array(3):ie_array(3),2)
    type(s_pp_grid),intent(in) :: ppg
    !
    real(8) :: x,y,z
    integer :: Nlma,ilma,ia,j,jx,jy,jz,ispin
    complex(8),parameter :: zero=(0.0d0,0.0d0),zi=(0.0d0,1.0d0)
    complex(8) :: eu,ed,pu,pd,sig(0:3)
    complex(8) :: uVpsi0,uVpsi_r(3,0:3),jtmp(3,0:3)
    real(8) :: j0_x,j0_y,j0_z, j1_x,j1_y,j1_z, j2_x,j2_y,j2_z, j3_x,j3_y,j3_z

    Nlma = size(ppg%ia_tbl_so)

    j0_x = 0d0
    j0_y = 0d0
    j0_z = 0d0
    j1_x = 0d0
    j1_y = 0d0
    j1_z = 0d0
    j2_x = 0d0
    j2_y = 0d0
    j2_z = 0d0
    j3_x = 0d0
    j3_y = 0d0
    j3_z = 0d0
    
#ifdef USE_OPENACC
!$acc loop worker private(ilma,ia,j,x,y,z,jx,jy,jz,eu,ed,pu,pd) &
!$acc private(uVpsi0,uVpsi_r,sig,jtmp) &
!$acc reduction(+:j0_x,j0_y,j0_z, j1_x,j1_y,j1_z, j2_x,j2_y,j2_z, j3_x,j3_y,j3_z)
#else
!$omp parallel do default(shared) private(ilma,ia,j,x,y,z,jx,jy,jz,eu,ed,pu,pd) &
!$omp private(uVpsi0,uVpsi_r,sig,jtmp) &
!$omp reduction(+:j0_x,j0_y,j0_z, j1_x,j1_y,j1_z, j2_x,j2_y,j2_z, j3_x,j3_y,j3_z)
#endif
    do ilma=1,Nlma

       ia=ppg%ia_tbl_so(ilma)

       uVpsi0 = zero
       uVpsi_r(:,:) = zero

       do j=1,ppg%mps(ia)
          x = ppg%rxyz(1,j,ia)
          y = ppg%rxyz(2,j,ia)
          z = ppg%rxyz(3,j,ia)
          jx = ppg%jxyz(1,j,ia)
          jy = ppg%jxyz(2,j,ia)
          jz = ppg%jxyz(3,j,ia)
          eu = conjg(ppg%zekr_uV_so(j,ilma,ik,1,1))
          ed = conjg(ppg%zekr_uV_so(j,ilma,ik,2,1))
          pu = psi(jx,jy,jz,1)
          pd = psi(jx,jy,jz,2)
          
          sig(0) = eu * pu + ed * pd
          sig(1) = eu * pd + ed * pu
          sig(2) = -zi* ( eu * pd - ed * pu )
          sig(3) = eu * pu - ed * pd
          
          uVpsi0 = uVpsi0 + sig(0)
          uVpsi_r(1,:) = uVpsi_r(1,:) + x * sig(:)
          uVpsi_r(2,:) = uVpsi_r(2,:) + y * sig(:)
          uVpsi_r(3,:) = uVpsi_r(3,:) + z * sig(:)
       end do
       uVpsi0 = uVpsi0 * ppg%rinv_uvu_so(ilma)
       jtmp = conjg( uVpsi_r )*uVpsi0
       
       j0_x = j0_x + 2d0* aimag(jtmp(1,0))
       j0_y = j0_y + 2d0* aimag(jtmp(2,0))
       j0_z = j0_z + 2d0* aimag(jtmp(3,0))
       
       j1_x = j1_x + 2d0* aimag(jtmp(1,1))
       j1_y = j1_y + 2d0* aimag(jtmp(2,1))
       j1_z = j1_z + 2d0* aimag(jtmp(3,1))
       
       j2_x = j2_x + 2d0* aimag(jtmp(1,2))
       j2_y = j2_y + 2d0* aimag(jtmp(2,2))
       j2_z = j2_z + 2d0* aimag(jtmp(3,2))

       j3_x = j3_x + 2d0* aimag(jtmp(1,3))
       j3_y = j3_y + 2d0* aimag(jtmp(2,3))
       j3_z = j3_z + 2d0* aimag(jtmp(3,3))

    end do !ilma

    jspin(1,0) = j0_x
    jspin(2,0) = j0_y
    jspin(3,0) = j0_z
    
    jspin(1,1) = j1_x
    jspin(2,1) = j1_y
    jspin(3,1) = j1_z
    
    jspin(1,2) = j2_x
    jspin(2,2) = j2_y
    jspin(3,2) = j2_z
    
    jspin(1,3) = j3_x
    jspin(2,3) = j3_y
    jspin(3,3) = j3_z
    return
  end subroutine calc_spin_current_nonlocal

end module pseudo_pt_current_so
