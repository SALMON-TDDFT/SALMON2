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
subroutine init_k_rd(tk_rd,tksquare,imode,brl)
use scf_data
!$ use omp_lib
implicit none
integer :: jj
real(8) :: tk_rd(3,num_kpoints_rd)
real(8) :: tksquare(num_kpoints_rd)
integer :: imode
real(8),intent(in) :: brl(3,3)
!
integer :: ix,iy,iz
integer :: ik
complex(8) :: vecA(3)
real(8) :: shift_k(3),k(3)

if(iSCFRT==1)then
  vecA=0.d0
else if(iSCFRT==2)then
  if(imode==1)then
    vecA(:)=A_ind(:,itt)
    if(epdir_re1(1)==1.d0)then
      vecA(1)=vecA(1)+A_ext(1,itt)
    else if(epdir_re1(2)==1.d0)then
      vecA(2)=vecA(2)+A_ext(2,itt)
    else if(epdir_re1(3)==1.d0)then
      vecA(3)=vecA(3)+A_ext(3,itt)
    end if
  else if(imode==2.or.imode==3)then
    vecA=0.d0
  else if(imode==4)then
    vecA(:)=A_ind(:,itt+1)
    if(epdir_re1(1)==1.d0)then
      vecA(1)=vecA(1)+A_ext(1,itt+1)
    else if(epdir_re1(2)==1.d0)then
      vecA(2)=vecA(2)+A_ext(2,itt+1)
    else if(epdir_re1(3)==1.d0)then
      vecA(3)=vecA(3)+A_ext(3,itt+1)
    end if
  end if
end if

if(ik_oddeven==1)then
  shift_k(1:3)=0.d0
  do jj=1,3
    if(num_kpoints_3d(jj)==1)then
      shift_k(jj)=0.5d0  ! tk_rd becomes zero
    end if
  end do
else if(ik_oddeven==2)then
  shift_k(1:3)=0.5d0
end if

do ik=1,num_kpoints_rd
  ix=mod(ik-1,num_kpoints_3d(1))+1
  iy=mod((ik-1)/num_kpoints_3d(1),num_kpoints_3d(2))+1
  iz=mod((ik-1)/(num_kpoints_3d(1)*num_kpoints_3d(2)),num_kpoints_3d(3))+1
  k(1) = (dble(ix)-shift_k(1))/dble(num_kpoints_3d(1))-0.5d0
  k(2) = (dble(iy)-shift_k(2))/dble(num_kpoints_3d(2))-0.5d0
  k(3) = (dble(iz)-shift_k(3))/dble(num_kpoints_3d(3))-0.5d0
  tk_rd(1,ik) = k(1)*Brl(1,1) + k(2)*Brl(1,2) + k(3)*Brl(1,3)
  tk_rd(2,ik) = k(1)*Brl(2,1) + k(2)*Brl(2,2) + k(3)*Brl(2,3)
  tk_rd(3,ik) = k(1)*Brl(3,1) + k(2)*Brl(3,2) + k(3)*Brl(3,3)
  tk_rd(:,ik) = tk_rd(:,ik) + vecA
  tksquare(ik)=tk_rd(1,ik)**2+tk_rd(2,ik)**2+tk_rd(3,ik)**2
end do

end subroutine init_k_rd
