!
!  Copyright 2017 SALMON developers
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
!=======================================================================
module writefield
  implicit none

contains

subroutine writedns(lg,mg,ng,rho,rmat,rmat2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,rho0,itt)
  use inputoutput, only: format3d,au_length_aa
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use writefile3d
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  real(8),intent(in) :: rho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8),intent(out) :: rmat2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  integer,intent(in) :: icoo1d(3,lg%num(1)*lg%num(2)*lg%num(3))
  real(8),intent(in) :: hgs(3)
  integer,intent(in) :: igc_is,igc_ie
  real(8),intent(in) :: gridcoo(igc_is:igc_ie,3)
  integer,intent(in) :: iscfrt
  real(8),intent(in),optional :: rho0(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer,intent(in),optional :: itt
  integer :: ix,iy,iz
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit

  !$OMP parallel do collapse(2) private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    rmat(ix,iy,iz)=0.d0
  end do
  end do
  end do

  !$OMP parallel do collapse(2) private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rmat(ix,iy,iz)=rho(ix,iy,iz)
  end do
  end do
  end do

  if(format3d=='avs')then
    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rmat(ix,iy,iz)=rmat(ix,iy,iz)/(au_length_aa**3)
    end do
    end do
    end do
  end if
  
  call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

  if(iSCFRT==1)then 
    suffix = "dns"
  else if(iSCFRT==2)then
    write(filenum, '(i6.6)') itt
    suffix = "dns_"//adjustl(filenum)
  end if
  phys_quantity = "dns"
  if(format3d=='avs')then
    header_unit='A**(-3)'
    call writeavs(lg,103,suffix,header_unit,rmat2,icoo1d)
  else if(format3d=='cube')then
    call writecube(lg,103,suffix,phys_quantity,rmat2,hgs,igc_is,igc_ie,gridcoo)
  else if(format3d=='vtk')then
    call writevtk(lg,103,suffix,rmat2,hgs,igc_is,igc_ie,gridcoo)
  end if

  if(iscfrt==2)then
    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      rmat(ix,iy,iz)=0.d0
    end do
    end do
    end do

    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rmat(ix,iy,iz)=rho(ix,iy,iz)-rho0(ix,iy,iz)
    end do
    end do
    end do
  
    if(format3d=='avs')then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=rmat(ix,iy,iz)/(au_length_aa**3)
      end do
      end do
      end do
    end if

    call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

    write(filenum, '(i6.6)') itt
    suffix = "dnsdiff_"//adjustl(filenum)
    phys_quantity = "dnsdiff"
    if(format3d=='avs')then
      header_unit='A**(-3)'
      call writeavs(lg,103,suffix,header_unit,rmat2,icoo1d)
    else if(format3d=='cube')then
      call writecube(lg,103,suffix,phys_quantity,rmat2,hgs,igc_is,igc_ie,gridcoo)
    else if(format3d=='vtk')then
      call writevtk(lg,103,suffix,rmat2,hgs,igc_is,igc_ie,gridcoo)
    end if
  end if
 
end subroutine writedns

end module writefield
