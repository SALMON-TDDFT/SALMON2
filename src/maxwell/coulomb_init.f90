!
!  Copyright 2019 SALMON developers
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
!-----------------------------------------------------------------------------------------
subroutine coulomb_init(ng_sta,ng_end,lg_sta,lg_end,hgs,fw)
  use structures,     only: s_fdtd_system, s_fdtd_field
  use salmon_maxwell, only: ls_fdtd_work
  use salmon_global, only: sysname
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use salmon_initialization, only: setbn
  implicit none
  integer,dimension(3),intent(in) :: ng_sta,ng_end,lg_sta,lg_end
  real(8),intent(in) :: hgs(3)
  type(ls_fdtd_work) :: fw
  !
  integer,parameter :: Nd=4
  integer :: lg_num(3),ii,jj
  lg_num = lg_end - lg_sta + 1

  fw%Energy_poynting = 0d0

  call setbn(fw%bnmat)
  do jj=1,3
    do ii=1,4
      fw%coef_nab(ii,jj) = fw%bnmat(ii,4)/hgs(jj)
    end do
  end do

  allocate( fw%vecA(3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
  allocate( fw%vecA_stock(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) )
  allocate( fw%Vh_n      (  ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) )
  allocate( fw%curr1_m   (  ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),3))

!1st element: time step (-1->m-1, 0->m, 1->m+1)
!5th element: components of A vector (1->Ax, 2->Ay, 3->Az)
  allocate( fw%vecA_m(-1:1,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:3) )

  allocate( fw%vecA_boundary_bottom(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
           ,fw%vecA_boundary_bottom_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
           ,fw%vecA_boundary_top(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
           ,fw%vecA_boundary_top_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) )

  allocate(fw%integral_poynting(lg_sta(3):lg_end(3)))
  allocate(fw%Ax_zt(lg_sta(3):lg_end(3)))

  fw%vecA = 0d0
  fw%vecA_stock = 0d0
  fw%vecA_m = 0d0
  fw%vecA_boundary_bottom = 0d0
  fw%vecA_boundary_bottom_old = 0d0
  fw%vecA_boundary_top = 0d0
  fw%vecA_boundary_top_old = 0d0
  fw%integral_poynting = 0d0
  fw%Ax_zt = 0d0

! temporary arrays
  allocate(fw%box(ng_sta(1)-Nd:ng_end(1)+Nd,ng_sta(2)-Nd:ng_end(2)+Nd,ng_sta(3)-Nd:ng_end(3)+Nd) &
        & ,fw%grad_Vh   (3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
        & ,fw%gradient_V(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
        & ,fw%rotation_A(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
   & ,fw%poynting_vector(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
        & ,fw%divergence_A(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
        & ,fw%vbox      (3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
        & ,fw%lgbox1      (lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
        & ,fw%lgbox2      (lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
        & ,fw%integral_poynting_tmp(lg_num(3)),fw%integral_poynting_tmp2(lg_num(3)) &
        & ,fw%vecA_ext    (ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),0:1,1:3) &
        & ,fw%vecA_ext_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),0:1,1:3) )

  if(comm_is_root(nproc_id_global)) then
    fw%file_777 = trim(sysname)//".777"
    open(777,file=fw%file_777)

  ! for spatial distribution of excitation energy
    fw%file_excitation = trim(sysname)//".excitation.data"
    open(555,file=fw%file_excitation)

  ! for the vector potential Ax(z,t)
    fw%file_Ax = trim(sysname)//".Ac.data"
    open(333,file=fw%file_Ax)
  end if

end subroutine coulomb_init
