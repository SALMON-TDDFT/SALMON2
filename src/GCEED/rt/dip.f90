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
subroutine subdip(ng,srho,rNe,ifunc)
use structures
use salmon_parallel, only: nproc_group_h, nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use scf_data
use new_world_sub
use allocate_mat_sub
use timer
implicit none
type(s_rgrid) ,intent(in) :: ng
type(s_scalar),intent(in) :: srho
integer :: ifunc
integer :: i1,ix,iy,iz,i2
real(8) :: rNe
real(8) :: rbox_array(10), rbox_arrayq(3, 3)
real(8) :: rbox_array2(10), rbox_arrayq2(3, 3)
real(8) :: rbox1, rbox1q, rbox1q12, rbox1q23, rbox1q31
real(8) :: fact
real(8) :: absr2

call timer_begin(LOG_CALC_DP)

!$OMP parallel do
   do i1=1,4
     rbox_array(i1)=0.d0
   end do
   do i1=1,3
     rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
     do iz=ng%is(3),ng%ie(3)
     do iy=ng%is(2),ng%ie(2)
     do ix=ng%is(1),ng%ie(1)
       rbox1=rbox1+vecR(i1,ix,iy,iz)*srho%f(ix,iy,iz)
     end do
     end do
     end do
     rbox_array(i1)=rbox1
   end do
   
   rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
   do iz=ng%is(3),ng%ie(3)
   do iy=ng%is(2),ng%ie(2)
   do ix=ng%is(1),ng%ie(1)
     rbox1=rbox1+srho%f(ix,iy,iz)
   end do
   end do
   end do
   rbox_array(4)=rbox1
   
   call timer_begin(LOG_ALLREDUCE_DIPOLE)
   call comm_summation(rbox_array,rbox_array2,4,nproc_group_h)
   call comm_summation(rbox_arrayq,rbox_arrayq2,9,nproc_group_h)
   call timer_end(LOG_ALLREDUCE_DIPOLE)

   rNe=rbox_array2(4)*Hvol               ! Number of electrons
   if(ifunc==1)then
     Dp(1:3,itt)=rbox_array2(1:3)*Hgs(1:3)*Hvol-vecDs(1:3)
     do i1=1,3
       Qp(1:3,i1,itt)=rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
     end do
     rIe(itt)=rNe
   else if(ifunc==2)then
     Dp(1:3,itt-1)=rbox_array2(1:3)*Hgs(1:3)*Hvol-vecDs(1:3)
     do i1=1,3
       Qp(1:3,i1,itt-1)=rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
     end do
     rIe(itt-1)=rNe
   end if

  if(comm_is_root(nproc_id_global))then
    select case(iperiodic)
    case(0)
      if(ifunc==1)then
          write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
            itt,dble(itt)*dt*2.41888d-2, (Dp(i1,itt)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh
        tene(itt)=Etot
      else if(ifunc==2)then
          write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
            itt-1,dble(itt-1)*dt*2.41888d-2, (Dp(i1,itt-1)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh
        tene(itt-1)=Etot
      end if
    case(3)
      write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8)')       &
        itt,dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, Etot*2d0*Ry
      if(iflag_md==1) then
        write(16,'(f14.8, 6e16.8, f15.8,f18.8)')       &
        dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3),(curr_ion(i1,itt),i1=1,3),rNe, Etot*2d0*Ry
      else
        write(16,'(f14.8, 3e16.8, f15.8,f18.8)')       &
        dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, Etot*2d0*Ry
      endif
      write(17,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_tot(i1,itt),i1=1,3)
      write(18,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ext(i1,itt),i1=1,3)
      write(19,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ind(i1,itt),i1=1,3)
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
    if(rNe.lt.(ifMST(1)+ifMST(2))*10.d0*fact.and.rNe.gt.(ifMST(1)+ifMST(2))/10.d0*fact)then
      continue
    else
      write(*,*) nproc_id_global,"t=",itt
      write(*,*) nproc_id_global,"rbox1=",rbox1
      write(*,*) nproc_id_global,"Ne=",rNe
      stop
    end if
  end if

call timer_end(LOG_CALC_DP)

end subroutine subdip
