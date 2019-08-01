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
subroutine calc_vecAc(vec_Ac,imode)
use scf_data
!$ use omp_lib
implicit none
real(8) :: vec_Ac(3)
integer :: imode
!
complex(8) :: vecA(3)

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
vec_Ac = vecA

end subroutine calc_vecAc
