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
subroutine setng(ng,ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s,   &
                 nproc,nproc_id_global,nproc_d_o,nproc_d_g_dm,ista_Mxin,iend_Mxin)
  use inputoutput, only: process_allocation
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(out) :: ng
  integer :: ii,i1,i2,i3,i4
  integer :: ibox
  integer :: nproc,nproc_id_global
  integer :: ng_sta(3),ng_end(3),ng_num(3)
  integer :: ista_Mxin_s(3,0:nproc-1),iend_Mxin_s(3,0:nproc-1),inum_Mxin_s(3,0:nproc-1)
  integer :: ista_Mxin(3,0:nproc-1),iend_Mxin(3,0:nproc-1)
  integer :: nproc_d_o_mul,nproc_d_g_mul_dm
  integer :: nproc_d_o(3),nproc_d_g_dm(3)
  integer :: j
  integer,parameter :: nd=4
  
  nproc_d_o_mul=nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
  nproc_d_g_mul_dm=nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
  
  if(process_allocation=='orbital_sequential')then
    do ii=0,nproc_d_o_mul-1
      do i4=0,nproc/nproc_d_o_mul/nproc_d_g_mul_dm-1
      do i3=0,nproc_d_g_dm(3)-1
      do i2=0,nproc_d_g_dm(2)-1
      do i1=0,nproc_d_g_dm(1)-1
        ibox= i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)       &
                +i4*nproc_d_g_mul_dm   &
                +ii*nproc/nproc_d_o_mul
        ista_Mxin_s(1,ibox)=i1*(iend_Mxin(1,ibox)-ista_Mxin(1,ibox)+1)/nproc_d_g_dm(1)+ista_Mxin(1,ibox)
        iend_Mxin_s(1,ibox)=(i1+1)*(iend_Mxin(1,ibox)-ista_Mxin(1,ibox)+1)/nproc_d_g_dm(1)+ista_Mxin(1,ibox)-1
        ista_Mxin_s(2,ibox)=i2*(iend_Mxin(2,ibox)-ista_Mxin(2,ibox)+1)/nproc_d_g_dm(2)+ista_Mxin(2,ibox)
        iend_Mxin_s(2,ibox)=(i2+1)*(iend_Mxin(2,ibox)-ista_Mxin(2,ibox)+1)/nproc_d_g_dm(2)+ista_Mxin(2,ibox)-1
        ista_Mxin_s(3,ibox)=i3*(iend_Mxin(3,ibox)-ista_Mxin(3,ibox)+1)/nproc_d_g_dm(3)+ista_Mxin(3,ibox)
        iend_Mxin_s(3,ibox)=(i3+1)*(iend_Mxin(3,ibox)-ista_Mxin(3,ibox)+1)/nproc_d_g_dm(3)+ista_Mxin(3,ibox)-1
      end do
      end do
      end do
      end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i4=0,nproc/nproc_d_o_mul/nproc_d_g_mul_dm-1
    do i3=0,nproc_d_g_dm(3)-1
    do i2=0,nproc_d_g_dm(2)-1
    do i1=0,nproc_d_g_dm(1)-1
      do ii=0,nproc_d_o_mul-1
        ibox=ii+(i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2))*nproc_d_o_mul   &
              +i4*nproc_d_o_mul*nproc_d_g_mul_dm
        ista_Mxin_s(1,ibox)=i1*(iend_Mxin(1,ii)-ista_Mxin(1,ii)+1)/nproc_d_g_dm(1)+ista_Mxin(1,ii)
        iend_Mxin_s(1,ibox)=(i1+1)*(iend_Mxin(1,ii)-ista_Mxin(1,ii)+1)/nproc_d_g_dm(1)+ista_Mxin(1,ii)-1
        ista_Mxin_s(2,ibox)=i2*(iend_Mxin(2,ii)-ista_Mxin(2,ii)+1)/nproc_d_g_dm(2)+ista_Mxin(2,ii)
        iend_Mxin_s(2,ibox)=(i2+1)*(iend_Mxin(2,ii)-ista_Mxin(2,ii)+1)/nproc_d_g_dm(2)+ista_Mxin(2,ii)-1
        ista_Mxin_s(3,ibox)=i3*(iend_Mxin(3,ii)-ista_Mxin(3,ii)+1)/nproc_d_g_dm(3)+ista_Mxin(3,ii)
        iend_Mxin_s(3,ibox)=(i3+1)*(iend_Mxin(3,ii)-ista_Mxin(3,ii)+1)/nproc_d_g_dm(3)+ista_Mxin(3,ii)-1
      end do
    end do
    end do
    end do
    end do
  end if
  inum_Mxin_s(1:3,0:nproc-1)=iend_Mxin_s(1:3,0:nproc-1)-ista_Mxin_s(1:3,0:nproc-1)+1
  
  ng_sta(:)=ista_Mxin_s(:,nproc_id_global)
  ng_end(:)=iend_Mxin_s(:,nproc_id_global)
  ng_num(:)=inum_Mxin_s(:,nproc_id_global)

  ng%is(1:3)=ng_sta(1:3)
  ng%ie(1:3)=ng_end(1:3)
  ng%num(1:3)=ng_num(1:3)
  ng%is_overlap(1:3)=ng_sta(1:3)-nd
  ng%ie_overlap(1:3)=ng_end(1:3)+nd
  ng%is_array(1:3)=ng_sta(1:3)-nd
  ng%ie_array(1:3)=ng_end(1:3)+nd

  if(allocated(ng%idx)) deallocate(ng%idx)
  if(allocated(ng%idy)) deallocate(ng%idy)
  if(allocated(ng%idz)) deallocate(ng%idz)
  allocate(ng%idx(ng%is_overlap(1):ng%ie_overlap(1)) &
          ,ng%idy(ng%is_overlap(2):ng%ie_overlap(2)) &
          ,ng%idz(ng%is_overlap(3):ng%ie_overlap(3)))
  do j=ng%is_overlap(1),ng%ie_overlap(1)
    ng%idx(j) = j
  end do
  do j=ng%is_overlap(2),ng%ie_overlap(2)
    ng%idy(j) = j
  end do
  do j=ng%is_overlap(3),ng%ie_overlap(3)
    ng%idz(j) = j
  end do

  ng%nd=4
  
end subroutine setng
