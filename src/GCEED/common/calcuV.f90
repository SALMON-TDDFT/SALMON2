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
SUBROUTINE calcuV(lg)
use structures,           only: s_rgrid
use salmon_communication, only: comm_is_root
use prep_pp_sub, only: set_nlma,init_lma_tbl,init_uv,set_lma_tbl,calc_uv
use scf_data
use allocate_psl_sub
use prep_pp_so_sub, only: calc_uv_so, SPIN_ORBIT_ON
use prep_pp_plusU_sub, only: calc_uv_plusU, PLUS_U_ON
implicit none
  type(s_rgrid),intent(in) :: lg

  integer :: i,ix,iy,iz
  integer :: nl

  real(8) :: hx,hy,hz
  integer :: lx(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: ly(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: lz(lg%num(1)*lg%num(2)*lg%num(3))

  character(17) :: property
  
  real(8),allocatable :: save_udVtbl_a(:,:,:)
  real(8),allocatable :: save_udVtbl_b(:,:,:)
  real(8),allocatable :: save_udVtbl_c(:,:,:)
  real(8),allocatable :: save_udVtbl_d(:,:,:)

  logical :: flag_use_grad_wf_on_force
  
  real(8) :: rinv_hvol 
  

  property='initial'
  flag_use_grad_wf_on_force=.false. 

  nl=lg%num(1)*lg%num(2)*lg%num(3)

  hx=Hgs(1)
  hy=Hgs(2)
  hz=Hgs(3)

  if(iperiodic==0)then
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      i=(iz-lg%is(3))*lg%num(1)*lg%num(2)+(iy-lg%is(2))*lg%num(1)+ix-lg%is(1)+1
      lx(i)=ix
      ly(i)=iy
      lz(i)=iz
    end do
    end do
    end do
  else if(iperiodic==3)then
    do iz=1,lg%num(3)
    do iy=1,lg%num(2)
    do ix=1,lg%num(1)
      i=(iz-1)*lg%num(1)*lg%num(2)+(iy-1)*lg%num(1)+ix
      lx(i)=ix-1
      ly(i)=iy-1
      lz(i)=iz-1
    end do
    end do
    end do
  end if
  
  call set_nlma(pp,ppg)
  call set_nlma(pp,ppg_all)

  call init_lma_tbl(pp,ppg)
  call init_lma_tbl(pp,ppg_all)

  call init_uv(pp,ppg)
  call init_uv(pp,ppg_all)

  call set_lma_tbl(pp,ppg)
  call set_lma_tbl(pp,ppg_all)

  allocate( save_udVtbl_a(pp%nrmax,0:2*pp%lmax+1,natom) )
  allocate( save_udVtbl_b(pp%nrmax,0:2*pp%lmax+1,natom) )
  allocate( save_udVtbl_c(pp%nrmax,0:2*pp%lmax+1,natom) )
  allocate( save_udVtbl_d(pp%nrmax,0:2*pp%lmax+1,natom) )

  call calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
               lx,ly,lz,nl,hx,hy,hz,  &
               flag_use_grad_wf_on_force,property,Hvol)

  call calc_uv(pp,ppg_all,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
               lx,ly,lz,nl,hx,hy,hz,  &
               flag_use_grad_wf_on_force,property,Hvol)

  if ( SPIN_ORBIT_ON ) then
    call calc_uv_so(pp,ppg,lx,ly,lz,nl,hx,hy,hz,flag_use_grad_wf_on_force,property,Hvol)
  end if

  if ( PLUS_U_ON ) then
    call calc_uv_plusU( pp, ppg, flag_use_grad_wf_on_force, property )
  end if

  rinv_hvol=1.d0/Hvol

return

END SUBROUTINE calcuV
