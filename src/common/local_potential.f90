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
MODULE local_potential
  implicit none

CONTAINS

!===================================================================================================================================

subroutine allgatherv_vlocal(ng,mg,info_field,nspin,Vh,Vpsl,Vxc,Vlocal)
  use structures
  use salmon_communication, only: comm_allgatherv
  use salmon_global, only: process_allocation
  use timer
  implicit none
  type(s_rgrid),         intent(in) :: ng,mg
  type(s_field_parallel),intent(in) :: info_field
  integer       ,intent(in) :: nspin
  type(s_scalar),intent(in) :: Vh,Vpsl,Vxc(nspin)
  type(s_scalar)            :: Vlocal(nspin)
  !
  integer :: i,is,myrank,nproc,i1,i2,i3,ix,iy,iz,ibox,ibox2,ibox3,iscnt
  real(8),allocatable :: matbox11(:),matbox12(:)
  integer,allocatable :: ircnt(:)
  integer,allocatable :: idisp(:)
  integer,dimension(3,0:info_field%isize_all-1) :: ista_Mxin_s,iend_Mxin_s,inum_Mxin_s

  myrank = info_field%id_all
  nproc = info_field%isize_all
  ista_mxin_s = ng%is_all
  iend_mxin_s = ng%ie_all
  inum_mxin_s = ng%ie_all - ng%is_all + 1

  allocate(ircnt(0:info_field%ngo_xyz-1))
  allocate(idisp(0:info_field%ngo_xyz-1))
  allocate(matbox11(0:(inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank))-1))
  allocate(matbox12(0:(mg%num(1)*mg%num(2)*mg%num(3))-1))

  iscnt = inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank)
  if(process_allocation=='orbital_sequential')then
    do i=0,info_field%ngo_xyz-1
      ibox=(myrank/info_field%ngo_xyz)*info_field%ngo_xyz+i
      ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
    end do
  else if(process_allocation=='grid_sequential')then
    do i=0,info_field%ngo_xyz-1
      ibox=mod(myrank,info_field%nproc_o)+i*info_field%nproc_o
      ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
    end do
  end if

  idisp(0)=0
  do i=1,info_field%ngo_xyz-1
    idisp(i)=idisp(i-1)+ircnt(i-1)
  end do

  do is=1,nspin
  !$OMP parallel do private(ibox3,ix,iy,iz)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      ibox3=ix-ng%is(1)+(iy-ng%is(2))*inum_Mxin_s(1,myrank)   &
                    +(iz-ng%is(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
      matbox11(ibox3) = Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(is)%f(ix,iy,iz)
    end do
    end do
    end do

    call timer_begin(LOG_ALLGATHERV_TOTAL)
    call comm_allgatherv(matbox11,matbox12,ircnt,idisp,info_field%icomm_v)
    call timer_end(LOG_ALLGATHERV_TOTAL)

    if(process_allocation=='orbital_sequential')then
  !$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
      do i3=0,info_field%ngo(3)-1
      do i2=0,info_field%ngo(2)-1
      do i1=0,info_field%ngo(1)-1
        ibox=(myrank/info_field%ngo_xyz)*info_field%ngo_xyz    &
              +(i1+i2*info_field%ngo(1)+i3*info_field%ngo(1)*info_field%ngo(2))
        ibox2=i1+i2*info_field%ngo(1)+i3*info_field%ngo(1)*info_field%ngo(2)

        call copyVlocal(mg,ista_Mxin_s(:,ibox),iend_Mxin_s(:,ibox), &
        & matbox12(idisp(ibox2):(idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))-1), &
        & Vlocal(is)%f)

      end do
      end do
      end do
    else if(process_allocation=='grid_sequential')then
  !$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
      do i3=0,info_field%ngo(3)-1
      do i2=0,info_field%ngo(2)-1
      do i1=0,info_field%ngo(1)-1
        ibox=mod(myrank,info_field%nproc_o)    &
            +(i1+i2*info_field%ngo(1)+i3*info_field%ngo(1)*info_field%ngo(2))*info_field%nproc_o
        ibox2=i1+i2*info_field%ngo(1)+i3*info_field%ngo(1)*info_field%ngo(2)

        call copyVlocal(mg,ista_Mxin_s(:,ibox),iend_Mxin_s(:,ibox), &
        & matbox12(idisp(ibox2):(idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))-1), &
        & Vlocal(is)%f)

      end do
      end do
      end do
    end if

  end do ! is=1,nspin

  deallocate (ircnt,idisp,matbox11,matbox12)

end subroutine allgatherv_vlocal

subroutine copyVlocal(mg,ista,iend,matbox,V)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: ista(3),iend(3)
  real(8),intent(in) :: matbox(ista(1):iend(1),     &
                               ista(2):iend(2),     &
                               ista(3):iend(3))
  real(8) :: V(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  V(ista(1):iend(1),     &
    ista(2):iend(2),     &
    ista(3):iend(3)) = &
  matbox(ista(1):iend(1),     &
         ista(2):iend(2),     &
         ista(3):iend(3))

  return
end subroutine copyVlocal

END MODULE local_potential

