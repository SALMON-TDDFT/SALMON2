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
subroutine buffer_broyden_ns(ng,system,srho,srho_s,mst,ifmst,iter)
  use structures, only: s_rgrid,s_dft_system,s_scalar
  use salmon_parallel, only: nproc_group_global
  use broyden_sub
  use scf_data, only: num_rho_stock,rho,rho_in,rho_out,  &
                      rho_s,rho_s_in,rho_s_out
  implicit none
  type(s_rgrid) :: ng
  type(s_dft_system),intent(in) :: system
  type(s_scalar),intent(inout) :: srho(1,1)
  type(s_scalar),intent(inout) :: srho_s(2,1)
  integer,intent(in) :: mst(2),ifmst(2)
  integer,intent(in) :: iter
  integer :: ix,iy,iz,is
  real(8) :: vecr(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))

  if(system%nspin==1)then

    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      vecr(ix,iy,iz)=rho(ix,iy,iz)
    end do
    end do
    end do

    call broyden(vecr,rho_in,rho_out,ng%num(1)*ng%num(2)*ng%num(3),  &
                 iter,num_rho_stock,num_rho_stock,nproc_group_global)
  
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rho(ix,iy,iz)= vecr(ix,iy,iz)
    end do
    end do
    end do

  else if(system%nspin==2)then
    
    do is=1,2
      if(ifmst(is)>=1.and.mst(is)>=1)then
        do iz=ng%is(3),ng%ie(3)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          vecr(ix,iy,iz)=rho_s(ix,iy,iz,is)
        end do
        end do
        end do
  
        call broyden(vecr,rho_s_in(ng%is(1):,ng%is(2):,ng%is(3):,1:,is),  &
                     rho_s_out(ng%is(1):,ng%is(2):,ng%is(3):,1:,is),  &
                     ng%num(1)*ng%num(2)*ng%num(3),  &
                     iter,num_rho_stock,num_rho_stock,nproc_group_global)
  
        do iz=ng%is(3),ng%ie(3)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          rho_s(ix,iy,iz,is)= vecr(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do

    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rho(ix,iy,iz)= rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
    end do
    end do
    end do
  end if

end subroutine buffer_broyden_ns
