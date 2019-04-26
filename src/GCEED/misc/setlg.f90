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
subroutine setlg(lg,lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
                 Hgs,Nd,rLsize1,imesh_oddeven,iperiodic,iscfrt)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(out) :: lg
  integer :: iperiodic
  integer :: lg_sta(3),lg_end(3),lg_num(3)
  integer,intent(in) :: iscfrt
  integer :: ista_Mx_ori(3),iend_Mx_ori(3),inum_Mx_ori(3)
  integer :: Nd
  integer :: imesh_oddeven(3)
  integer :: jj
  real(8) :: Hgs(3)
  real(8) :: rLsize1(3)
  real(8),parameter :: epsilon=1.d-10
  integer :: j
  
  
  select case(iperiodic)
    case(0)
      iend_Mx_ori(:)=int((rLsize1(:)+epsilon)/2.d0/Hgs(:))+Nd
      lg_end(:)=int((rLsize1(:)+epsilon)/2.d0/Hgs(:))
      
      do jj=1,3
        select case(imesh_oddeven(jj))
          case(1)
            ista_Mx_ori(jj)=-(int((rLsize1(jj)+epsilon)/2.d0/Hgs(jj))+Nd)
            lg_sta(jj)=-(int((rLsize1(jj)+epsilon)/2.d0/Hgs(jj)))
          case(2)
            ista_Mx_ori(jj)=-(int((rLsize1(jj)+epsilon)/2.d0/Hgs(jj))+Nd)+1
            lg_sta(jj)=-(int((rLsize1(jj)+epsilon)/2.d0/Hgs(jj)))+1
        end select
      end do
      
    case(3)
      ista_Mx_ori(:)=1-Nd
      lg_sta(:)=1
      iend_Mx_ori(:)=int((rLsize1(:)+epsilon)/Hgs(:))+Nd
      lg_end(:)=int((rLsize1(:)+epsilon)/Hgs(:))
  end select

  inum_Mx_ori(:)=iend_Mx_ori(:)-ista_Mx_ori(:)+1
  lg_num(:)=lg_end(:)-lg_sta(:)+1
  
  lg%is(1:3)=lg_sta(1:3)
  lg%ie(1:3)=lg_end(1:3)
  lg%num(1:3)=lg_num(1:3)
  lg%is_overlap(1:3)=lg_sta(1:3)-nd
  lg%ie_overlap(1:3)=lg_end(1:3)+nd
  lg%is_array(1:3)=lg_sta(1:3)-nd
  lg%ie_array(1:3)=lg_end(1:3)+nd

  if(allocated(lg%idx)) deallocate(lg%idx)
  if(allocated(lg%idy)) deallocate(lg%idy)
  if(allocated(lg%idz)) deallocate(lg%idz)
  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))
  do j=lg%is_overlap(1),lg%ie_overlap(1)
    lg%idx(j) = j
  end do
  do j=lg%is_overlap(2),lg%ie_overlap(2)
    lg%idy(j) = j
  end do
  do j=lg%is_overlap(3),lg%ie_overlap(3)
    lg%idz(j) = j
  end do

  lg%nd=4
  
end subroutine setlg
