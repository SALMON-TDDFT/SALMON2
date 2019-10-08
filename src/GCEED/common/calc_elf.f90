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
subroutine calc_elf(lg,mg,ng,system,info,stencil,srho,srg,srg_ng,tpsi,elf)
use structures
use math_constants, only: pi
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use sendrecv_grid, only: update_overlap_complex8,update_overlap_real8
use stencil_sub, only: calc_gradient_psi,calc_gradient_field
implicit none
type(s_rgrid)           ,intent(in) :: lg,mg,ng
type(s_dft_system)      ,intent(in) :: system
type(s_orbital_parallel),intent(in) :: info
type(s_stencil)         ,intent(in) :: stencil
type(s_scalar)          ,intent(in) :: srho
type(s_sendrecv_grid)               :: srg,srg_ng
type(s_orbital)                     :: tpsi
real(8)                             :: elf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
!
integer :: nspin,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
integer :: io,ik,ispin,ix,iy,iz
integer,parameter :: Nd=4
complex(8) :: ztmp

real(8) :: elftau(mg%is(1):mg%ie(1),   &
                  mg%is(2):mg%ie(2),   &
                  mg%is(3):mg%ie(3))
real(8) :: mrelftau(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: curden(mg%is(1):mg%ie(1),   &
                  mg%is(2):mg%ie(2),   &
                  mg%is(3):mg%ie(3))
real(8) :: mrcurden(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
complex(8) :: gradzpsi(3,mg%is_array(1):mg%ie_array(1),   &
                     mg%is_array(2):mg%ie_array(2),   &
                     mg%is_array(3):mg%ie_array(3))
real(8) :: gradrho(3,ng%is(1):ng%ie(1),   &
                     ng%is(2):ng%ie(2),   &
                     ng%is(3):ng%ie(3))
real(8) :: gradrho2(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: elfc(mg%is(1):mg%ie(1),   &
                mg%is(2):mg%ie(2),   &
                mg%is(3):mg%ie(3))
real(8) :: elfcuni(mg%is(1):mg%ie(1),   &
                   mg%is(2):mg%ie(2),   &
                   mg%is(3):mg%ie(3))
real(8) :: rho_half(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: box(ng%is(1)-Nd:ng%ie(1)+Nd,   &
               ng%is(2)-Nd:ng%ie(2)+Nd,   &
               ng%is(3)-Nd:ng%ie(3)+Nd)
real(8) :: matbox_l(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))

if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ calc_elf"

nspin = system%nspin
no = system%no
nk = system%nk
is = mg%is
ie = mg%ie
ik_s = info%ik_s
ik_e = info%ik_e
io_s = info%io_s
io_e = info%io_e

!$OMP parallel do private(iz,iy,ix)
do iz=is(3),ie(3)
do iy=is(2),ie(2)
do ix=is(1),ie(1)
  rho_half(ix,iy,iz)=srho%f(ix,iy,iz)/2.d0
end do
end do
end do
mrelftau=0.d0
mrcurden=0.d0

if(info%if_divide_rspace) then
   call update_overlap_complex8(srg, mg, tpsi%zwf)
end if

ik=1    ! --> do loop (future work)
ispin=1 ! --> do loop (future work)

  do io=1,no

      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,1),gradzpsi,mg%is_array,mg%ie_array,is,ie, &
      & mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)

!$OMP parallel do private(iz,iy,ix)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        ztmp = tpsi%zwf(ix,iy,iz,ispin,io,ik,1)
  
        mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                           +abs(gradzpsi(2,ix,iy,iz))**2      &
                           +abs(gradzpsi(3,ix,iy,iz))**2
  
        mrcurden(ix,iy,iz)=mrcurden(ix,iy,iz)      &
             +( abs(conjg(ztmp)*gradzpsi(1,ix,iy,iz)      &
                  -ztmp*conjg(gradzpsi(1,ix,iy,iz)))**2      &
               +abs(conjg(ztmp)*gradzpsi(2,ix,iy,iz)      &
                  -ztmp*conjg(gradzpsi(2,ix,iy,iz)))**2      &
               +abs(conjg(ztmp)*gradzpsi(3,ix,iy,iz)      &
                  -ztmp*conjg(gradzpsi(3,ix,iy,iz)))**2 )/2.d0
  
      end do
      end do
      end do
      
  end do

  call comm_summation(mrelftau,elftau,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_o)
  call comm_summation(mrcurden,curden,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_o)
  
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    box(ix,iy,iz) = srho%f(ix,iy,iz)
  end do
  end do
  end do
  
  call update_overlap_real8(srg_ng, ng, box)
  call calc_gradient_field(ng%is,ng%ie,stencil%coef_nab,box,gradrho)

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
          +gradrho(2,ix,iy,iz)**2      &
          +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0  &
                                   -curden(ix,iy,iz)/rho_half(ix,iy,iz)
  end do
  end do
  end do

! matbox_l stores ELF
matbox_l=0.d0
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  elfcuni(ix,iy,iz)=3.d0/5.d0*(6.d0*Pi**2)**(2.d0/3.d0)      &
            *rho_half(ix,iy,iz)**(5.d0/3.d0)
  matbox_l(ix,iy,iz)=1.d0/(1.d0+elfc(ix,iy,iz)**2/elfcuni(ix,iy,iz)**2)
end do
end do
end do

call comm_summation(matbox_l,elf,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

end subroutine calc_elf
