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
subroutine subdip(rt,ng,srho,rNe,poisson,Etot,system,pp)
use structures, only: s_rt,s_rgrid,s_scalar,s_poisson,s_dft_system,s_pp_info
use salmon_parallel, only: nproc_group_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use inputoutput, only: au_length_aa, natom, au_time_fs
use scf_data
use allocate_mat_sub
use inputoutput, only: yn_md
implicit none
type(s_rgrid) ,intent(in) :: ng
type(s_scalar),intent(in) :: srho
type(s_rt),intent(inout) :: rt
real(8),intent(in)        :: Etot
real(8),intent(out)       :: rNe
type(s_poisson),intent(in) :: poisson
type(s_dft_system),intent(in) :: system
type(s_pp_info),intent(in) :: pp
integer :: i1,ix,iy,iz, ia
real(8) :: rbox_array(10),  rbox_arrayq(3,3)
real(8) :: rbox_array2(10), rbox_arrayq2(3,3)
real(8) :: rbox1
real(8) :: fact, time, Hvol,Hgs(3)


Hvol   = system%Hvol
Hgs(:) = system%Hgs(:)

!$OMP parallel do
   do i1=1,4
     rbox_array(i1)=0.d0
   end do
   do i1=1,3
     rbox1=0d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
     do iz=ng%is(3),ng%ie(3)
     do iy=ng%is(2),ng%ie(2)
     do ix=ng%is(1),ng%ie(1)
        rbox1 = rbox1 + vecR(i1,ix,iy,iz)*srho%f(ix,iy,iz)
     end do
     end do
     end do
     rbox_array(i1)=rbox1
   end do
   
   rbox1=0d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
   do iz=ng%is(3),ng%ie(3)
   do iy=ng%is(2),ng%ie(2)
   do ix=ng%is(1),ng%ie(1)
      rbox1 = rbox1 + srho%f(ix,iy,iz)
   end do
   end do
   end do
   rbox_array(4)=rbox1
   
   call comm_summation(rbox_array,rbox_array2,4,nproc_group_global)
   call comm_summation(rbox_arrayq,rbox_arrayq2,9,nproc_group_global)

   !(ionic dipole) -- defined as plus charge (ordinary definition))
   rt%Dp_i(:,itt) = 0d0
   do ia=1,natom
      rt%Dp_i(:,itt) = rt%Dp_i(:,itt) + pp%Zps(Kion(ia)) * system%Rion(:,ia)
   enddo

   !(electronic dipole/quadrapole) -- defined as plus charge (opposite definition))
   rt%Dp_e(1:3,itt)  = -rbox_array2(1:3) * Hgs(1:3) * Hvol
   rt%dDp_e(1:3,itt) = rt%Dp_e(1:3,itt) - rt%Dp0_e(1:3)
   do i1=1,3
      rt%Qp_e(1:3,i1,itt) = -rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
   end do

   !(Number of electrons)
   rNe = rbox_array2(4)*Hvol
   rt%rIe(itt) = rNe

   if(comm_is_root(nproc_id_global))then
      time = dble(itt) * dt * au_time_fs
      select case(iperiodic)
      case(0)
         write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
             itt,time, rt%dDp_e(1:3,itt)*au_length_aa, rNe, Etot*2d0*Ry,poisson%iterVh
      case(3)
         write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8)')       &
             itt, time, rt%curr(1:3,itt), rNe, Etot*2d0*Ry
      end select
   end if

   fact=1.d0

   if(ilsda==0)then
      if(rNe.lt.ifMST(1)*2.d0*10.d0*fact.and.rNe.gt.ifMST(1)*2.d0/10.d0*fact)then
         continue
      else
         write(*,*) nproc_id_global,"t=",itt
         write(*,*) nproc_id_global,"rbox1=",rbox1
         write(*,*) nproc_id_global,"Ne=",rNe
         stop
      end if
   else if(ilsda==1)then
      if(rNe.lt.(ifMST(1)+ifMST(2))*10d0*fact.and.rNe.gt.(ifMST(1)+ifMST(2))/10d0*fact)then
         continue
      else
         write(*,*) nproc_id_global,"t=",itt
         write(*,*) nproc_id_global,"rbox1=",rbox1
         write(*,*) nproc_id_global,"Ne=",rNe
         stop
      end if
   end if


end subroutine subdip
