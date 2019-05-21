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
subroutine setmg(mg,mg_sta,mg_end,mg_num,mg_sta_all,mg_end_all,mg_num_all,  &
                 lg_sta,lg_num,nproc,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential,iscfrt)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(out) :: mg
  integer :: i1,i2,j1,j2,j3
  integer :: ibox
  integer :: isequential
  integer :: nproc,nproc_id_global
  integer :: nproc_Mxin(3)
  integer :: nproc_k
  integer :: nproc_ob
  integer :: mg_sta(3),mg_end(3),mg_num(3)
  integer :: lg_sta(3),lg_num(3)
  integer,intent(in) :: iscfrt
  integer :: mg_sta_all(3,0:nproc-1),mg_end_all(3,0:nproc-1),mg_num_all(3,0:nproc-1)
  integer :: nproc_Mxin_mul
  integer :: j
  integer,parameter :: nd=4
  
  nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
  if(isequential==1)then
    do j3=0,nproc_Mxin(3)-1
    do j2=0,nproc_Mxin(2)-1
    do j1=0,nproc_Mxin(1)-1
      do i2=0,nproc_k-1
      do i1=0,nproc_ob-1
        ibox = nproc_ob*i2 + i1 + nproc_k*nproc_ob*( j1 + nproc_Mxin(1)*j2 + nproc_Mxin(1)*nproc_Mxin(2)*j3 )
        mg_sta_all(1,ibox)=j1*lg_num(1)/nproc_Mxin(1)+lg_sta(1)
        mg_end_all(1,ibox)=(j1+1)*lg_num(1)/nproc_Mxin(1)+lg_sta(1)-1
        mg_sta_all(2,ibox)=j2*lg_num(2)/nproc_Mxin(2)+lg_sta(2)
        mg_end_all(2,ibox)=(j2+1)*lg_num(2)/nproc_Mxin(2)+lg_sta(2)-1
        mg_sta_all(3,ibox)=j3*lg_num(3)/nproc_Mxin(3)+lg_sta(3)
        mg_end_all(3,ibox)=(j3+1)*lg_num(3)/nproc_Mxin(3)+lg_sta(3)-1
      end do
      end do
    end do
    end do
    end do
  else if(isequential==2)then
    do i2=0,nproc_k-1
    do i1=0,nproc_ob-1
      do j3=0,nproc_Mxin(3)-1
      do j2=0,nproc_Mxin(2)-1
      do j1=0,nproc_Mxin(1)-1
        ibox = j1 + nproc_Mxin(1)*j2 + nproc_Mxin(1)*nproc_Mxin(2)*j3 + (nproc_ob*i2 + i1)*nproc_Mxin_mul
        mg_sta_all(1,ibox)=j1*lg_num(1)/nproc_Mxin(1)+lg_sta(1)
        mg_end_all(1,ibox)=(j1+1)*lg_num(1)/nproc_Mxin(1)+lg_sta(1)-1
        mg_sta_all(2,ibox)=j2*lg_num(2)/nproc_Mxin(2)+lg_sta(2)
        mg_end_all(2,ibox)=(j2+1)*lg_num(2)/nproc_Mxin(2)+lg_sta(2)-1
        mg_sta_all(3,ibox)=j3*lg_num(3)/nproc_Mxin(3)+lg_sta(3)
        mg_end_all(3,ibox)=(j3+1)*lg_num(3)/nproc_Mxin(3)+lg_sta(3)-1
      end do
      end do
      end do
    end do
    end do
  end if
  mg_num_all(:,:)=mg_end_all(:,:)-mg_sta_all(:,:)+1
  
  mg_sta(:)=mg_sta_all(:,nproc_id_global)
  mg_end(:)=mg_end_all(:,nproc_id_global)
  mg_num(:)=mg_num_all(:,nproc_id_global)
  
  mg%is(1:3)=mg_sta(1:3)
  mg%ie(1:3)=mg_end(1:3)
  mg%num(1:3)=mg_num(1:3)
  mg%is_overlap(1:3)=mg_sta(1:3)-nd
  mg%ie_overlap(1:3)=mg_end(1:3)+nd
  if(iscfrt==1)then
    mg%is_array(1:3)=mg_sta(1:3)-nd
    mg%ie_array(1:3)=mg_end(1:3)+nd
  else if(iscfrt==2)then
    mg%is_array(1:3)=mg_sta(1:3)-nd
    mg%ie_array(1)=mg_end(1)+nd+1
    mg%ie_array(2:3)=mg_end(2:3)+nd
  end if

  if(allocated(mg%idx)) deallocate(mg%idx)
  if(allocated(mg%idy)) deallocate(mg%idy)
  if(allocated(mg%idz)) deallocate(mg%idz)
  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))
  do j=mg%is_overlap(1),mg%ie_overlap(1)
    mg%idx(j) = j
  end do
  do j=mg%is_overlap(2),mg%ie_overlap(2)
    mg%idy(j) = j
  end do
  do j=mg%is_overlap(3),mg%ie_overlap(3)
    mg%idz(j) = j
  end do

  mg%nd=4
  
end subroutine setmg
