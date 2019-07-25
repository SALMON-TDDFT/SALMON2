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
<<<<<<< HEAD
subroutine exc_cor_ns(ppn, ng, srg_ng)
=======
subroutine exc_cor_ns(nspin,srho_s,ppn,sVxc,E_xc)
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
  use salmon_parallel, only: nproc_group_h
  use salmon_communication, only: comm_summation
  use salmon_xc, only: calc_xc
  use scf_data
  use new_world_sub
  use allocate_mat_sub
<<<<<<< HEAD
  use structures, only: s_pp_nlcc, s_rgrid,s_sendrecv_grid
  use sendrecv_grid, only: update_overlap_real8
=======
  use sendrecvh_sub
  use structures, only: s_pp_nlcc,s_scalar
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
  implicit none
  integer         ,intent(in) :: nspin
  type(s_scalar)  ,intent(in) :: srho_s(nspin)
  type(s_pp_nlcc), intent(in) :: ppn
<<<<<<< HEAD
  type(s_rgrid),intent(in) :: ng
  type(s_sendrecv_grid),intent(inout) :: srg_ng

=======
  type(s_scalar)              :: sVxc(nspin)
  real(8)                     :: E_xc
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0

  integer :: ix,iy,iz,is
  real(8) :: tot_exc
  real(8),allocatable :: rhd(:,:,:), delr(:,:,:,:)
  integer :: iwk_dum

  if(ilsda==0)then
<<<<<<< HEAD
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
    do ix=1,ng%num(1)
      rho_tmp(ix,iy,iz)=rho(ng%is(1)+ix-1,ng%is(2)+iy-1,ng%is(3)+iz-1)
=======
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      rho_tmp(ix,iy,iz)=srho_s(1)%f(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
    end do
    end do
    end do
  else if(ilsda==1)then
    do is=1,2
<<<<<<< HEAD
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
    do ix=1,ng%num(1)
      rho_s_tmp(ix,iy,iz,is)=rho_s(ng%is(1)+ix-1,ng%is(2)+iy-1,ng%is(3)+iz-1,is)
=======
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      rho_s_tmp(ix,iy,iz,is)=srho_s(is)%f(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
    end do
    end do
    end do
    end do
  end if

  if(xc=='pz'.or.xc=='PZ')then
    continue
  else
    allocate (rhd (ng%is_array(1):ng%ie_array(1), &
                   ng%is_array(2):ng%ie_array(2), &
                   ng%is_array(3):ng%ie_array(3)))
    allocate (delr(ng%is(1):ng%ie(1), &
                   ng%is(2):ng%ie(2), &
                   ng%is(3):ng%ie(3),3))
  
!$OMP parallel do private(ix,iy,iz)
<<<<<<< HEAD
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rhd(ix,iy,iz)=dble(rho(ix,iy,iz))
=======
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rhd(ix,iy,iz)=dble(srho_s(1)%f(ix,iy,iz))
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
    enddo
    enddo
    enddo
  
!$omp end parallel do
    iwk_dum=iwk_size
    iwk_size=12
    call make_iwksta_iwkend
    call update_overlap_real8(srg_ng, ng, rhd)

    iwk_size=iwk_dum

!$OMP parallel do private(ix,iy,iz)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       delr(ix,iy,iz,1)= (1.0d0/Hgs(1))*(   &
          bN1*( rhd(ix+1,iy,iz) - rhd(ix-1,iy,iz))  &
        + bN2*( rhd(ix+2,iy,iz) - rhd(ix-2,iy,iz)) &
        + bN3*( rhd(ix+3,iy,iz) - rhd(ix-3,iy,iz)) &
        + bN4*( rhd(ix+4,iy,iz) - rhd(ix-4,iy,iz)))
       delr(ix,iy,iz,2)= (1.0d0/Hgs(2))*( &
          bN1*( rhd(ix,iy+1,iz) - rhd(ix,iy-1,iz)) &
        + bN2*( rhd(ix,iy+2,iz) - rhd(ix,iy-2,iz)) &
        + bN3*( rhd(ix,iy+3,iz) - rhd(ix,iy-3,iz)) &
        + bN4*( rhd(ix,iy+4,iz) - rhd(ix,iy-4,iz)))
       delr(ix,iy,iz,3)= (1.0d0/Hgs(3))*( &
          bN1*( rhd(ix,iy,iz+1) - rhd(ix,iy,iz-1)) &
        + bN2*( rhd(ix,iy,iz+2) - rhd(ix,iy,iz-2)) &
        + bN3*( rhd(ix,iy,iz+3) - rhd(ix,iy,iz-3)) &
        + bN4*( rhd(ix,iy,iz+4) - rhd(ix,iy,iz-4)))
    enddo
    enddo
    enddo
!$omp end parallel do
  end if

  if(xc=='pz'.or.xc=='PZ')then
    if(ilsda==0)then
      call calc_xc(xc_func, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rho_nlcc=ppn%rho_nlcc)
    else if(ilsda==1)then
      call calc_xc(xc_func, rho_s=rho_s_tmp, eexc=eexc_tmp, vxc_s=vxc_s_tmp, rho_nlcc=ppn%rho_nlcc)
    end if
  else
    call calc_xc(xc_func, rho=rho_tmp, grho=delr, eexc=eexc_tmp, vxc=vxc_tmp, rho_nlcc=ppn%rho_nlcc)
  end if

  if(ilsda==0)then
<<<<<<< HEAD
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
    do ix=1,ng%num(1)
      vxc(ng%is(1)+ix-1,ng%is(2)+iy-1,ng%is(3)+iz-1)=vxc_tmp(ix,iy,iz)
=======
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      sVxc(1)%f(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)=vxc_tmp(ix,iy,iz)
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
    end do
    end do
    end do
  else if(ilsda==1)then 
    do is=1,2
<<<<<<< HEAD
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
    do ix=1,ng%num(1)
      vxc_s(ng%is(1)+ix-1,ng%is(2)+iy-1,ng%is(3)+iz-1,is)=vxc_s_tmp(ix,iy,iz,is)
=======
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      sVxc(is)%f(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)=vxc_s_tmp(ix,iy,iz,is)
>>>>>>> 43c8d5496b16d74c85549914cb33f849b8af01f0
    end do
    end do
    end do
    end do
  end if

  tot_exc=0.d0
  do iz=1,ng%num(3)
  do iy=1,ng%num(2)
  do ix=1,ng%num(1)
    tot_exc=tot_exc+eexc_tmp(ix,iy,iz)
  end do
  end do
  end do
  tot_exc=tot_exc*hvol
 
  call comm_summation(tot_exc,E_xc,nproc_group_h)

  return
end subroutine exc_cor_ns
