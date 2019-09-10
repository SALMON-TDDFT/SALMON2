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
MODULE init_communicator
  implicit none

CONTAINS

!=======================================================================
subroutine init_communicator_dft(comm,info,info_field)
  use salmon_global, only: process_allocation,nproc_domain_orbital,nproc_domain_general,nproc_ob,nproc_k,ispin
  use structures, only: s_orbital_parallel, s_field_parallel
  use salmon_communication, only: comm_create_group, comm_get_groupinfo, &
                                  comm_summation
  use misc_routines, only: get_wtime
  use sendrecv_grid, only: s_sendrecv_grid, init_sendrecv_grid !??
  use init_sendrecv_sub!??
  implicit none
  integer,      intent(in) :: comm
  type(s_orbital_parallel) :: info
  type(s_field_parallel)   :: info_field
  !
  integer :: myrank,nproc
  integer :: nproc_d_o(3),nproc_d_g(3),nproc_d_o_mul,nproc_d_g_dm(3),nproc_d_g_mul_dm,nproc_ob_spin(2)
  integer :: imr(3),imrs(3),igroup
  integer :: i1,i2,i3,i4,i5,ix,iy,iz,ixs,iys,izs
  integer :: ibox,icolor,ikey
  integer :: npuy,npuz

  call comm_get_groupinfo(comm, myrank, nproc)

  ! info

  info%icomm_rko  = comm
  info%id_rko = myrank
  info%isize_rko = nproc

  nproc_d_o = nproc_domain_orbital
  nproc_d_g = nproc_domain_general
  nproc_d_o_mul = nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
  nproc_d_g_dm = nproc_d_g/nproc_d_o
  nproc_d_g_mul_dm = nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

  if(ispin==1)then
    nproc_ob_spin(1) = (nproc_ob+1)/2
    nproc_ob_spin(2) = nproc_ob/2
  end if

  !new_world for comm_kgrid
  if(process_allocation=='orbital_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        if(myrank==ibox)then
          icolor=i5*nproc_ob+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
          ikey=i4
        end if
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
        if(myrank==ibox)then
          icolor=i5+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_k
          ikey=i4
        end if
      end do
      end do
    end do
    end do
    end do
  end if

  info%icomm_o = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info%icomm_o, info%id_o, info%isize_o)

  !new_world for comm_korbital
  if(process_allocation=='orbital_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        if(myrank==ibox)then
          icolor=i5*nproc_ob+i4
          ikey=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
        end if
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
        if(myrank==ibox)then
          icolor=i5*nproc_ob+i4
          ikey=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
        end if
      end do
      end do
    end do
    end do
    end do
  end if

  info%icomm_r = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info%icomm_r, info%id_r, info%isize_r)

  !new_world for comm_k
  if(process_allocation=='orbital_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        if(ispin==0)then
          if(myrank==ibox)then
            icolor=i5
            ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob
          end if
        else
          if(i4<nproc_ob_spin(1))then
            if(myrank==ibox)then
              icolor=2*i5+0
              ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(1)
            end if
          else
            if(myrank==ibox)then
              icolor=2*i5+1
              ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(2)
            end if
          end if
        end if
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
        if(ispin==0)then
          if(myrank==ibox)then
            icolor=i5
            ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob
          end if
        else
          if(i4<nproc_ob_spin(1))then
            if(myrank==ibox)then
              icolor=2*i5+0
              ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(1)
            end if
          else
            if(myrank==ibox)then
              icolor=2*i5+1
              ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(2)
            end if
          end if
        end if
      end do
      end do
    end do
    end do
    end do
  end if

  info%icomm_ro = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info%icomm_ro, info%id_ro, info%isize_ro)

  !new_world for comm_grid
  if(process_allocation=='orbital_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        if(ispin==0)then
          if(myrank==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob+i4
          end if
        else
          if(i4<nproc_ob_spin(1))then
            if(myrank==ibox)then
              icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
              ikey=i5*nproc_ob_spin(1)+i4
            end if
          else
            if(myrank==ibox)then
              icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
              ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
            end if
          end if
        end if
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
        if(ispin==0)then
          if(myrank==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob+i4
          end if
        else
          if(i4<nproc_ob_spin(1))then
            if(myrank==ibox)then
              icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
              ikey=i5*nproc_ob_spin(1)+i4
            end if
          else
            if(myrank==ibox)then
              icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
              ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
            end if
          end if
        end if
      end do
      end do
    end do
    end do
    end do
  end if

  info%icomm_ko = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info%icomm_ko, info%id_ko, info%isize_ko)

  !new_world for comm_orbitalgrid
  if(process_allocation=='orbital_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        if(myrank==ibox)then
          icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)+i4*nproc_d_o_mul
          ikey=i5
        end if
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      do i5=0,nproc_k-1
      do i4=0,nproc_ob-1
        ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
        if(myrank==ibox)then
          icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)+i4*nproc_d_o_mul
          ikey=i5
        end if
      end do
      end do
    end do
    end do
    end do
  end if

  info%icomm_k = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info%icomm_k, info%id_k, info%isize_k)

  ! info_field

  info_field%icomm_all = info%icomm_rko
  info_field%id_all = info%id_rko
  info_field%isize_all = info%isize_rko

  if(process_allocation=='orbital_sequential')then
    do iz=0,nproc_d_o(3)-1
    do iy=0,nproc_d_o(2)-1
    do ix=0,nproc_d_o(1)-1
      do i4=0,nproc/nproc_d_o_mul/nproc_d_g_mul_dm-1
      do izs=0,nproc_d_g_dm(3)-1
      do iys=0,nproc_d_g_dm(2)-1
      do ixs=0,nproc_d_g_dm(1)-1
        ibox=ixs+iys*nproc_d_g_dm(1)   &
                +izs*nproc_d_g_dm(1)*nproc_d_g_dm(2)   &
                +i4*nproc_d_g_mul_dm    &
                +ix*nproc/nproc_d_o_mul    &
                +iy*nproc/nproc_d_o_mul*nproc_d_o(1)   &
                +iz*nproc/nproc_d_o_mul*nproc_d_o(1)*nproc_d_o(2)
        if(myrank==ibox)then
          imr(1)=ix
          imr(2)=iy
          imr(3)=iz
          imrs(1)=ixs
          imrs(2)=iys
          imrs(3)=izs
          igroup=i4
        end if
      end do
      end do
      end do
    end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i4=0,nproc/nproc_d_o_mul/nproc_d_g_mul_dm-1
    do izs=0,nproc_d_g_dm(3)-1
    do iys=0,nproc_d_g_dm(2)-1
    do ixs=0,nproc_d_g_dm(1)-1
      do iz=0,nproc_d_o(3)-1
      do iy=0,nproc_d_o(2)-1
      do ix=0,nproc_d_o(1)-1
        ibox=ix+iy*nproc_d_o(1)+iz*nproc_d_o(1)*nproc_d_o(2)  &
               +ixs*nproc_d_o_mul  &
               +iys*nproc_d_o_mul*nproc_d_g_dm(1)  &
               +izs*nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
               +i4*nproc_d_o_mul*nproc_d_g_mul_dm
        if(myrank==ibox)then
          imr(1)=ix
          imr(2)=iy
          imr(3)=iz
          imrs(1)=ixs
          imrs(2)=iys
          imrs(3)=izs
          igroup=i4
        end if
      end do
      end do
      end do
    end do
    end do
    end do
    end do
  end if

  if(process_allocation=='orbital_sequential')then
    icolor=imrs(2)+imrs(3)*nproc_d_g_dm(2)   &
                  +igroup*nproc_d_g_dm(2)*nproc_d_g_dm(3)   &
                  +imr(2)*nproc/nproc_d_o_mul/nproc_d_g_dm(1)   &
                  +imr(3)*nproc/nproc_d_o_mul/nproc_d_g_dm(1)*nproc_d_o(2)
    ikey=imrs(1)+imr(1)*nproc_d_g_dm(1)
  else if(process_allocation=='grid_sequential')then
    icolor=imr(2)+imr(3)*nproc_d_o(2)   &
                 +imrs(2)*nproc_d_o(2)*nproc_d_o(3)   &
                 +imrs(3)*nproc_d_o(2)*nproc_d_o(3)*nproc_d_g_dm(2)  &
                 +myrank/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(2)  &
                 *nproc_d_o(3)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
  !  ikey=imr(1)+imrs(1)*nproc_d_o(1)
    ikey=imrs(1)+imr(1)*nproc_d_g_dm(1)
  end if

  info_field%icomm(1) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm(1), info_field%id(1), info_field%isize(1))

  if(process_allocation=='orbital_sequential')then
    icolor=imrs(1)+imrs(3)*nproc_d_g_dm(1)   &
                  +igroup*nproc_d_g_dm(1)*nproc_d_g_dm(3)   &
                  +imr(1)*nproc/nproc_d_o_mul/nproc_d_g_dm(2)   &
                  +imr(3)*nproc/nproc_d_o_mul/nproc_d_g_dm(2)*nproc_d_o(1)
    ikey=imrs(2)+imr(2)*nproc_d_g_dm(2)
  else if(process_allocation=='grid_sequential')then
    icolor=imr(1)+imr(3)*nproc_d_o(1)   &
                 +imrs(1)*nproc_d_o(1)*nproc_d_o(3)   &
                 +imrs(3)*nproc_d_o(1)*nproc_d_o(3)*nproc_d_g_dm(1)  &
                 +myrank/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(1)*nproc_d_o(3)  &
                 *nproc_d_g_dm(1)*nproc_d_g_dm(3)
  !  ikey=imr(2)+imrs(2)*nproc_d_o(2)
    ikey=imrs(2)+imr(2)*nproc_d_g_dm(2)
  end if

  info_field%icomm(2) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm(2), info_field%id(2), info_field%isize(2))

  if(process_allocation=='orbital_sequential')then
    icolor=imrs(1)+imrs(2)*nproc_d_g_dm(1)   &
                  +igroup*nproc_d_g_dm(1)*nproc_d_g_dm(2)   &
                  +imr(1)*nproc/nproc_d_o_mul/nproc_d_g_dm(3)   &
                  +imr(2)*nproc/nproc_d_o_mul/nproc_d_g_dm(3)*nproc_d_o(1)
    ikey=imrs(3)+imr(3)*nproc_d_g_dm(3)
  else if(process_allocation=='grid_sequential')then
    icolor=imr(1)+imr(2)*nproc_d_o(1)   &
                 +imrs(1)*nproc_d_o(1)*nproc_d_o(2)   &
                 +imrs(2)*nproc_d_o(1)*nproc_d_o(2)*nproc_d_g_dm(1)  &
                 +myrank/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(1)*nproc_d_o(2)  &
                 *nproc_d_g_dm(1)*nproc_d_g_dm(2)
  !  ikey=imr(3)+imrs(3)*nproc_d_o(3)
    ikey=imrs(3)+imr(3)*nproc_d_g_dm(3)
  end if

  info_field%icomm(3) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm(3), info_field%id(3), info_field%isize(3))

! communicators for FFTE routine
  npuy = nproc_d_g_dm(2)*nproc_d_o(2)
  npuz = nproc_d_g_dm(3)*nproc_d_o(3)

  icolor=info_field%id(3)+info_field%id(1)*npuz
  ikey=info_field%id(2)
  info_field%icomm_ffte(2) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(2), info_field%id_ffte(2), info_field%isize_ffte(2))

  icolor=info_field%id(2)+info_field%id(1)*npuy
  ikey=info_field%id(3)
  info_field%icomm_ffte(3) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(3), info_field%id_ffte(3), info_field%isize_ffte(3))

  icolor=info_field%id(2)+info_field%id(3)*npuy
  ikey=info_field%id(1)
  info_field%icomm_ffte(1) = comm_create_group(comm, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(1), info_field%id_ffte(1), info_field%isize_ffte(1))

  info%imr = imr

  info_field%imr = imr
  info_field%imrs = imrs

end subroutine init_communicator_dft

END MODULE init_communicator

