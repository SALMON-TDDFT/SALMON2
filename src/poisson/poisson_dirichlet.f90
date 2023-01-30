!
!  Copyright 2019-2023 SALMON developers
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
!===================================================================================================================================

module poisson_dirichlet
  implicit none

contains

subroutine jones(lg,mg,info,system,rho,Vh,poisson)
  use structures
  use communication, only: comm_summation
  use math_constants, only : pi
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_dft_system)     ,intent(in) :: system
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  !
  integer           :: ix,iy,iz
  real(8)           :: dx,dy,dz
  integer           :: lx_s,lx_e
  integer           :: ly_s,ly_e
  integer           :: lz_s,lz_e
  integer           :: mx_s,mx_e
  integer           :: my_s,my_e
  integer           :: mz_s,mz_e
  integer           :: lx,ly,lz
  integer           :: mx,my,mz
  character(1)      :: yn_discrete_green
  real(8),allocatable :: rho_shift(:,:,:)
  real(8),allocatable :: sigma_1(:,:,:)
  real(8),allocatable :: sigma_2(:,:,:)
  real(8),allocatable :: theta_1(:,:,:)
  real(8),allocatable :: theta_2(:,:,:)
  real(8),allocatable :: phi_b(:,:,:)
  real(8),allocatable :: phi_tilde(:,:,:)

  lx=lg%num(1)
  ly=lg%num(2)
  lz=lg%num(3)
  mx=mg%num(1)
  my=mg%num(2)
  mz=mg%num(3)

  lx_s=1
  lx_e=lg%num(1)
  ly_s=1
  ly_e=lg%num(2)
  lz_s=1
  lz_e=lg%num(3)
  mx_s=mg%is(1)-lg%is(1)+1
  mx_e=mg%ie(1)-lg%is(1)+1
  my_s=mg%is(2)-lg%is(2)+1
  my_e=mg%ie(2)-lg%is(2)+1
  mz_s=mg%is(3)-lg%is(3)+1
  mz_e=mg%ie(3)-lg%is(3)+1
 
  allocate(rho_shift(mx_s:mx_e,my_s:my_e,mz_s:mz_e))
  allocate(sigma_1(0:lx+1,0:ly+1,0:lz+1))
  allocate(sigma_2(0:lx+1,0:ly+1,0:lz+1))
  allocate(theta_1(0:lx+1,0:ly+1,0:lz+1))
  allocate(theta_2(0:lx+1,0:ly+1,0:lz+1))
  allocate(phi_b(0:lx+1,0:ly+1,0:lz+1))
  allocate(phi_tilde(mx_s:mx_e,my_s:my_e,mz_s:mz_e))
  
  dx=system%hgs(1)
  dy=system%hgs(2)
  dz=system%hgs(3)

!$OMP parallel do private(iz,iy,ix)
  do iz=mz_s,mz_e
  do iy=my_s,my_e
  do ix=mx_s,mx_e
    rho_shift(ix,iy,iz)=rho%f(ix+lg%is(1)-1,iy+lg%is(2)-1,iz+lg%is(3)-1)
  end do
  end do
  end do

  yn_discrete_green='n'
  call calc_interior_potential(lx_s,lx_e,ly_s,ly_e,lz_s,lz_e,     &
                               mx_s,mx_e,my_s,my_e,mz_s,mz_e,     &
                               dx,dy,dz,rho_shift,phi_tilde,yn_discrete_green,info)


  sigma_1=0.d0
  if(mz_s==lz_s)then
!$OMP parallel do private(iy,ix)
    do iy=my_s,my_e
    do ix=mx_s,mx_e
      sigma_1(ix,iy,0)=phi_tilde(ix,iy,1)/4.d0/pi/dz**2
    end do 
    end do
  end if
  if(mz_e==lz_e)then
!$OMP parallel do private(iy,ix)
    do iy=my_s,my_e
    do ix=mx_s,mx_e
      sigma_1(ix,iy,lz+1)=phi_tilde(ix,iy,lz)/4.d0/pi/dz**2
    end do 
    end do 
  end if
  if(mx_s==lx_s)then
!$OMP parallel do private(iz,iy)
    do iz=mz_s,mz_e
    do iy=my_s,my_e
      sigma_1(0,iy,iz)=phi_tilde(1,iy,iz)/4.d0/pi/dx**2
    end do 
    end do
  end if
  if(mx_e==lx_e)then
!$OMP parallel do private(iz,iy)
    do iz=mz_s,mz_e
    do iy=my_s,my_e
      sigma_1(lx+1,iy,iz)=phi_tilde(lx,iy,iz)/4.d0/pi/dx**2
    end do 
    end do 
  end if
  if(my_s==ly_s)then
!$OMP parallel do private(iz,ix)
    do iz=mz_s,mz_e
    do ix=mx_s,mx_e
      sigma_1(ix,0,iz)=phi_tilde(ix,1,iz)/4.d0/pi/dy**2
    end do 
    end do 
  end if
  if(my_e==ly_e)then
!$OMP parallel do private(iz,ix)
    do iz=mz_s,mz_e
    do ix=mx_s,mx_e
      sigma_1(ix,ly+1,iz)=phi_tilde(ix,ly,iz)/4.d0/pi/dy**2
    end do 
    end do 
  end if

  call comm_summation(sigma_1,sigma_2,(lx+2)*(ly+2)*(lz+2),info%icomm_r)

  theta_1=0.d0
  if(mx_s==lx_s)then
    do iz=mz_s,mz_e
    do iy=my_s,my_e
      call calc_theta(lx,ly,lz,0,iy,iz,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if
  if(mx_e==lx_e)then
    do iz=mz_s,mz_e
    do iy=my_s,my_e
      call calc_theta(lx,ly,lz,lx+1,iy,iz,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if
  if(my_s==ly_s)then
    do iz=mz_s,mz_e
    do ix=mx_s,mx_e
      call calc_theta(lx,ly,lz,ix,0,iz,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if
  if(my_e==ly_e)then
    do iz=mz_s,mz_e
    do ix=mx_s,mx_e
      call calc_theta(lx,ly,lz,ix,ly+1,iz,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if
  if(mz_s==lz_s)then
    do iy=my_s,my_e
    do ix=mx_s,mx_e
      call calc_theta(lx,ly,lz,ix,iy,0,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if
  if(mz_e==lz_e)then
    do iy=my_s,my_e
    do ix=mx_s,mx_e
      call calc_theta(lx,ly,lz,ix,iy,lz+1,dx,dy,dz,poisson%dgf,sigma_2,theta_1)
    end do
    end do
  end if

  call comm_summation(theta_1,theta_2,(lx+2)*(ly+2)*(lz+2),info%icomm_r)

!$OMP parallel do private(iz,iy,ix)
  do iz=0,lz+1
  do iy=0,ly+1
  do ix=0,lx+1
    phi_b(ix,iy,iz)=-theta_2(ix,iy,iz)
  end do
  end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=mz_s,mz_e
  do iy=my_s,my_e
  do ix=mx_s,mx_e
    rho_shift(ix,iy,iz)=rho_shift(ix,iy,iz)-(phi_b(ix-1,iy,iz)+phi_b(ix+1,iy,iz))/dx**2/4.d0/pi &
                                           -(phi_b(ix,iy-1,iz)+phi_b(ix,iy+1,iz))/dy**2/4.d0/pi &
                                           -(phi_b(ix,iy,iz-1)+phi_b(ix,iy,iz+1))/dz**2/4.d0/pi
    ! the term phi_b(ix,iy,iz) is omitted because it is zero in the active cell.
  end do
  end do
  end do

  yn_discrete_green='n'
  call calc_interior_potential(lx_s,lx_e,ly_s,ly_e,lz_s,lz_e,     &
                               mx_s,mx_e,my_s,my_e,mz_s,mz_e,     &
                               dx,dy,dz,rho_shift,phi_tilde,yn_discrete_green,info)
  

!$OMP parallel do private(iz,iy,ix)
  do iz=mz_s,mz_e
  do iy=my_s,my_e
  do ix=mx_s,mx_e
    Vh%f(ix+lg%is(1)-1,iy+lg%is(2)-1,iz+lg%is(3)-1)=-phi_tilde(ix,iy,iz)
  end do
  end do
  end do

  deallocate(rho_shift)
  deallocate(sigma_1)
  deallocate(sigma_2)
  deallocate(theta_1)
  deallocate(theta_2)
  deallocate(phi_b)
  deallocate(phi_tilde)

  return

end subroutine jones

!==================================================================================================

subroutine calc_interior_potential(lx_s,lx_e,ly_s,ly_e,lz_s,lz_e,  &
                                   mx_s,mx_e,my_s,my_e,mz_s,mz_e,  &
                                   dx,dy,dz,rho,phi_tilde,yn_discrete_green,info)
  use structures, only: s_parallel_info
  use communication, only: comm_summation
  use pack_unpack, only: copy_data
  use math_constants, only : pi
  implicit none
  integer, intent(in)               :: lx_s,lx_e
  integer, intent(in)               :: ly_s,ly_e
  integer, intent(in)               :: lz_s,lz_e
  integer, intent(in)               :: mx_s,mx_e
  integer, intent(in)               :: my_s,my_e
  integer, intent(in)               :: mz_s,mz_e
  real(8), intent(in)               :: dx,dy,dz
  real(8), intent(in)               :: rho(mx_s:mx_e,my_s:my_e,mz_s:mz_e)
  real(8), intent(out)              :: phi_tilde(mx_s:mx_e,my_s:my_e,mz_s:mz_e)
  character(1), intent(in)          :: yn_discrete_green
  type(s_parallel_info) ,intent(in) :: info
  integer                           :: lx,ly,lz
  integer                           :: mx,my,mz
  integer                           :: ix,iy,iz
  integer                           :: kx,ky,kz
  real(8)                           :: xil(lx_s:lx_e,lx_s:lx_e)
  real(8)                           :: yjm(ly_s:ly_e,ly_s:ly_e)
  real(8)                           :: zkn(lz_s:lz_e,lz_s:lz_e)
  real(8)                           :: rlambda_x(lx_s:lx_e)
  real(8)                           :: rlambda_y(ly_s:ly_e)
  real(8)                           :: rlambda_z(lz_s:lz_e)
  real(8)                           :: ff1x(lx_s:lx_e,my_s:my_e,mz_s:mz_e)
  real(8)                           :: ff2x(lx_s:lx_e,my_s:my_e,mz_s:mz_e)
  real(8)                           :: ff1y(mx_s:mx_e,ly_s:ly_e,mz_s:mz_e)
  real(8)                           :: ff2y(mx_s:mx_e,ly_s:ly_e,mz_s:mz_e)
  real(8)                           :: ff1z(mx_s:mx_e,my_s:my_e,lz_s:lz_e)
  real(8)                           :: ff2z(mx_s:mx_e,my_s:my_e,lz_s:lz_e)
  real(8)                           :: phi_tilde_lmn(mx_s:mx_e,my_s:my_e,mz_s:mz_e)

  lx=lx_e-lx_s+1
  ly=ly_e-ly_s+1
  lz=lz_e-lz_s+1
  mx=mx_e-mx_s+1
  my=my_e-my_s+1
  mz=mz_e-mz_s+1

!$OMP parallel do private(kx,ix)
  do kx=lx_s,lx_e
    do ix=lx_s,lx_e
      xil(ix,kx)=sin(pi*dble(ix-lx_s+1)*dble(kx-lx_s+1)/dble(lx+1))
    end do
  end do
!$OMP parallel do private(ky,iy)
  do ky=ly_s,ly_e
    do iy=ly_s,ly_e
      yjm(iy,ky)=sin(pi*dble(iy-ly_s+1)*dble(ky-ly_s+1)/dble(ly+1))
    end do
  end do
!$OMP parallel do private(kz,iz)
  do kz=lz_s,lz_e
    do iz=lz_s,lz_e
      zkn(iz,kz)=sin(pi*dble(iz-lz_s+1)*dble(kz-lz_s+1)/dble(lz+1))
    end do
  end do

!$OMP parallel do private(kx)
  do kx=lx_s,lx_e
    rlambda_x(kx)=-(sin(pi*dble(kx-lx_s+1)/(2.d0*dble(lx+1)))*2.d0/dx)**2
  end do
!$OMP parallel do private(ky)
  do ky=ly_s,ly_e
    rlambda_y(ky)=-(sin(pi*dble(ky-ly_s+1)/(2.d0*dble(ly+1)))*2.d0/dy)**2
  end do
!$OMP parallel do private(kz)
  do kz=lz_s,lz_e
    rlambda_z(kz)=-(sin(pi*dble(kz-lz_s+1)/(2.d0*dble(lz+1)))*2.d0/dz)**2
  end do

  ff1x=0.d0  
  ff1y=0.d0  
  ff1z=0.d0

!$OMP parallel do private(iz,iy,ix)
  do iz=mz_s,mz_e
  do iy=my_s,my_e
  do ix=mx_s,mx_e
    ff1z(ix,iy,iz)=rho(ix,iy,iz)
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1z,ff2z)
  else
    call comm_summation(ff1z,ff2z,mx*my*lz,info%icomm_z)
  end if

!$OMP parallel do private(kz,iy,ix)
  do kz = mz_s,mz_e
  do iy = my_s,my_e
  do ix = mx_s,mx_e
    ff1y(ix,iy,kz) = sum(ff2z(ix,iy,:)*zkn(:,kz))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1y,ff2y)
  else
    call comm_summation(ff1y,ff2y,mx*ly*mz,info%icomm_y)
  end if

!$OMP parallel do private(kz,ky,ix)
  do kz = mz_s,mz_e
  do ky = my_s,my_e
  do ix = mx_s,mx_e
    ff1x(ix,ky,kz) = sum(ff2y(ix,:,kz)*yjm(:,ky))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1x,ff2x)
  else
    call comm_summation(ff1x,ff2x,lx*my*mz,info%icomm_x)
  end if

!$OMP parallel do private(kz,ky,kx)
  do kz = mz_s,mz_e
  do ky = my_s,my_e
  do kx = mx_s,mx_e
    ff1x(kx,ky,kz) = sum(ff2x(:,ky,kz)*xil(:,kx))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1x,ff2x)
  else
    call comm_summation(ff1x,ff2x,lx*my*mz,info%icomm_x)
  end if

  phi_tilde_lmn=0.d0
!$OMP parallel do private(kz,ky,kx)
  do kz=mz_s,mz_e
  do ky=my_s,my_e
  do kx=mx_s,mx_e
    phi_tilde_lmn(kx,ky,kz)=4.d0*pi*ff2x(kx,ky,kz)/(rlambda_x(kx)+rlambda_y(ky)+rlambda_z(kz))
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,kx)
  do kz = mz_s,mz_e
  do ky = my_s,my_e
  do kx = mx_s,mx_e
    ff1z(kx,ky,kz)=phi_tilde_lmn(kx,ky,kz)
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1z,ff2z)
  else
    call comm_summation(ff1z,ff2z,mx*my*lz,info%icomm_z)
  end if

!$OMP parallel do private(iz,ky,kx)
  do iz = mz_s,mz_e
  do ky = my_s,my_e
  do kx = mx_s,mx_e
    ff1y(kx,ky,iz) = sum(ff2z(kx,ky,:)*zkn(iz,:))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1y,ff2y)
  else
    call comm_summation(ff1y,ff2y,mx*ly*mz,info%icomm_y)
  end if

!$OMP parallel do private(iz,iy,kx)
  do iz = mz_s,mz_e
  do iy = my_s,my_e
  do kx = mx_s,mx_e
    ff1x(kx,iy,iz) = sum(ff2y(kx,:,iz)*yjm(iy,:))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1x,ff2x)
  else
    call comm_summation(ff1x,ff2x,lx*my*mz,info%icomm_x)
  end if

!$OMP parallel do private(iz,iy,ix)
  do iz = mz_s,mz_e
  do iy = my_s,my_e
  do ix = mx_s,mx_e
    ff1x(ix,iy,iz) = sum(ff2x(:,iy,iz)*xil(ix,:))
  end do
  end do
  end do
  if(yn_discrete_green=='y')then
    call copy_data(ff1x,ff2x)
  else
    call comm_summation(ff1x,ff2x,lx*my*mz,info%icomm_x)
  end if

!$OMP parallel do private(iz,iy,ix)
  do iz=mz_s,mz_e
  do iy=my_s,my_e
  do ix=mx_s,mx_e
    phi_tilde(ix,iy,iz)=8.d0/dble(lx+1)/dble(ly+1)/dble(lz+1)*ff2x(ix,iy,iz)
  end do
  end do
  end do

  return

end subroutine calc_interior_potential

!==================================================================================================

subroutine calc_theta(lx,ly,lz,ix,iy,iz,dx,dy,dz,dgf,sigma,theta)
  implicit none
  integer, intent(in)    :: lx,ly,lz
  integer, intent(in)    :: ix,iy,iz
  real(8), intent(in)    :: dx,dy,dz
  real(8), intent(in)    :: sigma(0:lx+1,0:ly+1,0:lz+1)
  real(8), intent(in)    :: dgf(-lx:lx+2,-ly:ly+2,-lz:lz+2)
  real(8), intent(inout) :: theta(0:lx+1,0:ly+1,0:lz+1)
  integer                :: kx,ky,kz
  real(8)                :: tmp


  tmp=0.d0
!$OMP parallel do reduction(+:tmp) private(ky,kz)
  do ky=1,ly
  do kx=1,lx
    tmp=tmp+dgf(ix-kx+1,iy-ky+1,iz-0+1)*sigma(kx,ky,0)*dx*dy*dz
    tmp=tmp+dgf(ix-kx+1,iy-ky+1,iz-(lz+1)+1)*sigma(kx,ky,lz+1)*dx*dy*dz
  end do
  end do
  theta(ix,iy,iz)=theta(ix,iy,iz)+tmp

  tmp=0.d0
!$OMP parallel do reduction(+:tmp) private(kz,ky)
  do kz=1,lz
  do ky=1,ly
    tmp=tmp+dgf(ix-0+1,iy-ky+1,iz-kz+1)*sigma(0,ky,kz)*dx*dy*dz
    tmp=tmp+dgf(ix-(lx+1)+1,iy-ky+1,iz-kz+1)*sigma(lx+1,ky,kz)*dx*dy*dz
  end do
  end do
  theta(ix,iy,iz)=theta(ix,iy,iz)+tmp

  tmp=0.d0
!$OMP parallel do reduction(+:tmp) private(kz,kx)
  do kz=1,lz
  do kx=1,lx
    tmp=tmp+dgf(ix-kx+1,iy-0+1,iz-kz+1)*sigma(kx,0,kz)*dx*dy*dz
    tmp=tmp+dgf(ix-kx+1,iy-(ly+1)+1,iz-kz+1)*sigma(kx,ly+1,kz)*dx*dy*dz
  end do
  end do
  theta(ix,iy,iz)=theta(ix,iy,iz)+tmp

  return

end subroutine calc_theta

!==================================================================================================

subroutine calc_dgf(lg,mg,system,info,poisson)
  use structures
  use math_constants, only : pi
  use communication, only: comm_summation
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_dft_system)     ,intent(in) :: system
  type(s_parallel_info)  ,intent(in) :: info
  type(s_poisson)                    :: poisson
  !
  integer           :: ix,iy,iz
  real(8)           :: dx,dy,dz
  integer           :: iix,iiy,iiz
  integer           :: lx_s,lx_e
  integer           :: ly_s,ly_e
  integer           :: lz_s,lz_e
  integer           :: mx_s,mx_e
  integer           :: my_s,my_e
  integer           :: mz_s,mz_e
  integer           :: lx,ly,lz
  integer           :: mx,my,mz
  character(1)      :: yn_discrete_green
  real(8),allocatable :: rho_unit(:,:,:) 
  real(8),allocatable :: phi_b_for_dgf(:,:,:) 
  real(8),allocatable :: dgf_partial(:,:,:)
  
  lx=lg%num(1)
  ly=lg%num(2)
  lz=lg%num(3)
  mx=mg%num(1)
  my=mg%num(2)
  mz=mg%num(3)

  lx_s=1
  lx_e=lg%num(1)
  ly_s=1
  ly_e=lg%num(2)
  lz_s=1
  lz_e=lg%num(3)
  mx_s=mg%is(1)-lg%is(1)+1
  mx_e=mg%ie(1)-lg%is(1)+1
  my_s=mg%is(2)-lg%is(2)+1
  my_e=mg%ie(2)-lg%is(2)+1
  mz_s=mg%is(3)-lg%is(3)+1
  mz_e=mg%ie(3)-lg%is(3)+1
  
  allocate(rho_unit(-15:lx+2,-15:ly+2,-15:lz+2))
  allocate(phi_b_for_dgf(-16:lx+3,-16:ly+3,-16:lz+3))
  allocate(dgf_partial(-15:lx+2,-15:ly+2,-15:lz+2))

  dx=system%hgs(1)
  dy=system%hgs(2)
  dz=system%hgs(3)

  rho_unit=0.d0
  rho_unit(1,1,1)=1.d0/dx/dy/dz

  phi_b_for_dgf=0.d0
  do iz=-16,-16
  do iy=-15,ly+2
  do ix=-15,lx+2
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
  do iz=lz+3,lz+3
  do iy=-15,ly+2
  do ix=-15,lx+2
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
  do iz=-15,lz+2
  do iy=-16,-16
  do ix=-15,lx+2
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
  do iz=-15,lz+2
  do iy=ly+3,ly+3
  do ix=-15,lx+2
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
  do iz=-15,lz+2
  do iy=-15,ly+2
  do ix=-16,-16
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
  do iz=-15,lz+1
  do iy=-15,ly+1
  do ix=lx+3,lx+3
    phi_b_for_dgf(ix,iy,iz) = -1.d0/sqrt((dble(ix)-dble(1))**2*dx**2  &
                                        +(dble(iy)-dble(1))**2*dy**2  &
                                        +(dble(iz)-dble(1))**2*dz**2)
  end do
  end do
  end do
 
  do iz=-15,lz+2
  do iy=-15,ly+2
  do ix=-15,lx+2
    rho_unit(ix,iy,iz)=rho_unit(ix,iy,iz)-(phi_b_for_dgf(ix-1,iy,iz)+phi_b_for_dgf(ix+1,iy,iz))/dx**2/4.d0/pi &
                                         -(phi_b_for_dgf(ix,iy-1,iz)+phi_b_for_dgf(ix,iy+1,iz))/dy**2/4.d0/pi &
                                         -(phi_b_for_dgf(ix,iy,iz-1)+phi_b_for_dgf(ix,iy,iz+1))/dz**2/4.d0/pi
    ! the term phi_b_for_dgf(ix,iy,iz) is omitted because it is zero in the active cell.
  end do
  end do
  end do


  yn_discrete_green='y'
  call calc_interior_potential(-15,lx+2,-15,ly+2,-15,lz+2,     &
                               -15,lx+2,-15,ly+2,-15,lz+2,     &
                               dx,dy,dz,rho_unit,dgf_partial,yn_discrete_green,info)

  do iz=-lz,lz+2
  do iy=-ly,ly+2
  do ix=-lx,lx+2
    iix=abs(ix-1)+1
    iiy=abs(iy-1)+1
    iiz=abs(iz-1)+1
    poisson%dgf(ix,iy,iz)=dgf_partial(iix,iiy,iiz)
  end do
  end do
  end do

  deallocate(rho_unit)
  deallocate(phi_b_for_dgf)
  deallocate(dgf_partial)

  return
end subroutine calc_dgf

end module poisson_dirichlet

