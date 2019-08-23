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
subroutine copy_density(nspin,ng,srho_s,mixing)
use structures, only: s_rgrid, s_scalar, s_mixing
use scf_data
implicit none
integer       ,intent(in) :: nspin
type(s_rgrid), intent(in) :: ng
type(s_scalar),intent(in) :: srho_s(nspin)
type(s_mixing),intent(inout) :: mixing
integer :: iiter
integer :: is
integer :: ix,iy,iz

if(Miter==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
  end do
  end do
  end do
  if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      mixing%srho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
      mixing%srho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz)=srho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
end if

do iiter=1,mixing%num_rho_stock
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_in(iiter)%f(ix,iy,iz)=mixing%srho_in(iiter+1)%f(ix,iy,iz)
  end do
  end do
  end do
end do
do iiter=1,mixing%num_rho_stock-1
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_out(iiter)%f(ix,iy,iz)=mixing%srho_out(iiter+1)%f(ix,iy,iz)
  end do
  end do
  end do
end do

if(ilsda==1)then
  do iiter=1,mixing%num_rho_stock
    do is=1,2
!$OMP parallel do private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        mixing%srho_s_in(iiter,is)%f(ix,iy,iz)=mixing%srho_s_in(iiter+1,is)%f(ix,iy,iz)
      end do
      end do
      end do
    end do
  end do
  do iiter=1,mixing%num_rho_stock-1
    do is=1,2
!$OMP parallel do private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        mixing%srho_s_out(iiter,is)%f(ix,iy,iz)=mixing%srho_s_out(iiter+1,is)%f(ix,iy,iz)
      end do
      end do
      end do
    end do
  end do
end if

end subroutine copy_density
