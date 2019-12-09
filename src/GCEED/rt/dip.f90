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
subroutine subdip(itt,rt,lg,ng,srho,rNe,poisson,Etot,system,pp)
use structures, only: s_rt,s_rgrid,s_scalar,s_poisson,s_dft_system,s_pp_info
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use inputoutput, only: au_length_aa, au_energy_ev, natom, au_time_fs
use salmon_global
use inputoutput, only: yn_md
implicit none
integer,intent(in) :: itt
type(s_rgrid) ,intent(in) :: lg
type(s_rgrid) ,intent(in) :: ng
type(s_scalar),intent(in) :: srho
type(s_rt),intent(inout) :: rt
real(8),intent(in)        :: Etot
real(8),intent(out)       :: rNe
type(s_poisson),intent(in) :: poisson
type(s_dft_system),intent(in) :: system
type(s_pp_info),intent(in) :: pp
integer :: ia
real(8) :: rbox_array2(4)
real(8) :: rbox1
real(8) :: time, Hvol,Hgs(3)


Hvol   = system%Hvol
Hgs(:) = system%Hgs(:)

   call dip(lg,ng,srho,rbox_array2)

   !(ionic dipole) -- defined as plus charge (ordinary definition))
   rt%Dp_i(:,itt) = 0d0
   do ia=1,natom
      rt%Dp_i(:,itt) = rt%Dp_i(:,itt) + pp%Zps(Kion(ia)) * system%Rion(:,ia)
   enddo

   !(electronic dipole/quadrapole) -- defined as plus charge (opposite definition))
   rt%Dp_e(1:3,itt)  = -rbox_array2(1:3) * Hgs(1:3) * Hvol
   rt%dDp_e(1:3,itt) = rt%Dp_e(1:3,itt) - rt%Dp0_e(1:3)

   !(Number of electrons)
   rNe = rbox_array2(4)*Hvol
   rt%rIe(itt) = rNe

   if(comm_is_root(nproc_id_global))then
      time = dble(itt) * dt * au_time_fs
      select case(iperiodic)
      case(0)
         write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
             itt,time, rt%dDp_e(1:3,itt)*au_length_aa, rNe, Etot*au_energy_ev,poisson%iterVh
      case(3)
         write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8)')       &
             itt, time, rt%curr(1:3,itt), rNe, Etot*au_energy_ev
      end select
   end if

   if(system%nspin==2.and.sum(nelec_spin(:))>0)then
      if(rNe.lt.dble(sum(nelec_spin(:)))*10.d0.and.rNe.gt.dble(sum(nelec_spin(:)))/10.d0)then
         continue
      else
         write(*,*) nproc_id_global,"t=",itt
         write(*,*) nproc_id_global,"rbox1=",rbox1
         write(*,*) nproc_id_global,"Ne=",rNe
         stop
      end if
   else 
      if(rNe.lt.dble(nelec)*10d0.and.rNe.gt.dble(nelec)/10d0)then
         continue
      else
         write(*,*) nproc_id_global,"t=",itt
         write(*,*) nproc_id_global,"rbox1=",rbox1
         write(*,*) nproc_id_global,"Ne=",rNe
         stop
      end if
   end if


end subroutine subdip

subroutine dip(lg,ng,srho,rbox_array2)
  use structures, only: s_rgrid,s_scalar
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  implicit none
  type(s_rgrid) ,intent(in) :: lg
  type(s_rgrid) ,intent(in) :: ng
  type(s_scalar),intent(in) :: srho
  integer :: ix,iy,iz
  real(8) :: rbox_array(4),rbox_array2(4)
  real(8) :: rbox

  rbox_array=0d0
  
  rbox=0.d0
  select case(mod(lg%num(1),2))
  case(1)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       rbox = rbox + dble(ix) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  case(0)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rbox = rbox + (dble(ix)-0.5d0) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  end select
  rbox_array(1)=rbox

  rbox=0.d0
  select case(mod(lg%num(2),2))
  case(1)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rbox = rbox + dble(iy) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  case(0)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rbox = rbox + (dble(iy)-0.5d0) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  end select
  rbox_array(2)=rbox
 
  rbox=0.d0
  select case(mod(lg%num(3),2))
  case(1)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       rbox = rbox + dble(iz) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  case(0)
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       rbox = rbox + (dble(iz)-0.5d0) * srho%f(ix,iy,iz)
    end do
    end do
    end do
  end select
  rbox_array(3)=rbox

  rbox=0.d0
!$OMP parallel do reduction( + : rbox ) private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
     rbox = rbox + srho%f(ix,iy,iz)
  end do
  end do
  end do
  rbox_array(4)=rbox

  call comm_summation(rbox_array,rbox_array2,4,nproc_group_global)

end subroutine dip
