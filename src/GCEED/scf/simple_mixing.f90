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
SUBROUTINE simple_mixing(ng,nspin,c1,c2,srho_s)
use structures, only: s_rgrid, s_scalar
use scf_data
implicit none
type(s_rgrid),intent(in) :: ng
integer,intent(in) :: nspin
real(8),intent(in) :: c1,c2
type(s_scalar),intent(inout) :: srho_s(nspin)

integer :: ix,iy,iz

if(ilsda == 0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_out(ix,iy,iz,num_rho_stock)=srho_s(1)%f(ix,iy,iz)
  end do
  end do
  end do
elseif( ilsda==1 )then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_s_out(ix,iy,iz,num_rho_stock,1)=srho_s(1)%f(ix,iy,iz)
    rho_s_out(ix,iy,iz,num_rho_stock,2)=srho_s(2)%f(ix,iy,iz)
  end do
  end do
  end do
end if

!rho = c1*rho + c2*matmul( psi**2, occ )
if(ilsda == 0)then
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    srho_s(1)%f(ix,iy,iz) = c1*rho_in(ix,iy,iz,num_rho_stock) + c2*rho_out(ix,iy,iz,num_rho_stock)
    rho_in(ix,iy,iz,num_rho_stock+1) = srho_s(1)%f(ix,iy,iz)
  end do
  end do
  end do
else if(ilsda == 1)then
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    srho_s(1)%f(ix,iy,iz) = c1*rho_s_in(ix,iy,iz,num_rho_stock,1) + c2*rho_s_out(ix,iy,iz,num_rho_stock,1)
    srho_s(2)%f(ix,iy,iz) = c1*rho_s_in(ix,iy,iz,num_rho_stock,2) + c2*rho_s_out(ix,iy,iz,num_rho_stock,2)
    rho_s_in(ix,iy,iz,num_rho_stock+1,1) = srho_s(1)%f(ix,iy,iz)
    rho_s_in(ix,iy,iz,num_rho_stock+1,2) = srho_s(2)%f(ix,iy,iz)
  end do
  end do
  end do
end if


return

END SUBROUTINE simple_mixing

