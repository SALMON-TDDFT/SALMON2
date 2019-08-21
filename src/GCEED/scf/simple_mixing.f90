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
SUBROUTINE simple_mixing(ng,system,c1,c2,srho_s,mixing)
use structures, only: s_rgrid, s_dft_system, s_scalar, s_mixing
implicit none
type(s_rgrid),intent(in) :: ng
type(s_dft_system),intent(in) :: system
real(8),intent(in) :: c1,c2
type(s_scalar),intent(inout) :: srho_s(system%nspin)
type(s_mixing),intent(inout) :: mixing

integer :: ix,iy,iz

if(system%nspin == 1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_out(mixing%num_rho_stock)%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
  end do
  end do
  end do
elseif(system%nspin == 2)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
    mixing%srho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=srho_s(2)%f(ix,iy,iz)
  end do
  end do
  end do
end if

!rho = c1*rho + c2*matmul( psi**2, occ )
if(system%nspin == 1)then
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    srho_s(1)%f(ix,iy,iz) = c1*mixing%srho_in(mixing%num_rho_stock)%f(ix,iy,iz) &
                            + c2*mixing%srho_out(mixing%num_rho_stock)%f(ix,iy,iz)
    mixing%srho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = srho_s(1)%f(ix,iy,iz)
  end do
  end do
  end do
else if(system%nspin == 2)then
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    srho_s(1)%f(ix,iy,iz) = c1*mixing%srho_s_in(mixing%num_rho_stock,1)%f(ix,iy,iz) &
                            + c2*mixing%srho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)
    srho_s(2)%f(ix,iy,iz) = c1*mixing%srho_s_in(mixing%num_rho_stock,2)%f(ix,iy,iz) &
                            + c2*mixing%srho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)
    mixing%srho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz) = srho_s(1)%f(ix,iy,iz)
    mixing%srho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz) = srho_s(2)%f(ix,iy,iz)
  end do
  end do
  end do
end if


return

END SUBROUTINE simple_mixing

