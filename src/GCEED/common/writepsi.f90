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
!======================================================================
!======================================================================
subroutine writepsi(lg,info)
  use inputoutput, only: au_length_aa
  use structures, only: s_rgrid,s_orbital_parallel
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use writefile3d
  use calc_myob_sub
  use check_corrkob_sub
  use scf_data
  use allocate_mat_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_orbital_parallel),intent(in) :: info
  integer :: iob,ix,iy,iz
  integer :: p0,iob_myob,icheck_corrkob
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit
 
  if(iSCFRT==1)then
    do p0=1,itotMST
      call conv_p0(p0,iob)
      call calc_myob(iob,iob_myob,ilsda,nproc_ob,itotmst,mst)
      call check_corrkob(iob,info,1,icheck_corrkob,ilsda,nproc_ob,k_sta,k_end,mst)
  !OMP parallel do private(iz,iy,ix)
      do iz=lg%is(3),lg%ie(3)
      do iy=lg%is(2),lg%ie(2)
      do ix=lg%is(1),lg%ie(1)
        matbox_l(ix,iy,iz)=0.d0
      end do
      end do
      end do
      if(icheck_corrkob==1)then
  !OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_l(ix,iy,iz)=psi(ix,iy,iz,iob_myob,1)
        end do
        end do
        end do
      end if

      if(format_voxel_data=='avs')then
  !OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_l(ix,iy,iz)=matbox_l(ix,iy,iz)/sqrt(au_length_aa)**3
        end do
        end do
        end do
      end if

      call comm_summation(matbox_l,matbox_l2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

      write(filenum, '(i5)') p0
      suffix = "psi"//trim(adjustl(filenum))
      phys_quantity = "psi"
      if(format_voxel_data=='avs')then
        header_unit = "A**(-3/2)"
        call writeavs(lg,103,suffix,header_unit,matbox_l2,icoo1d)
      else if(format_voxel_data=='cube')then
        call writecube(lg,103,suffix,phys_quantity,matbox_l2,hgs)
      end if
    end do
  end if
  
end subroutine writepsi

