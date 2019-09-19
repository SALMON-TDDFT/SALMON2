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
module mixing_sub
  implicit none

contains

!===================================================================================================================================
subroutine simple_mixing(ng,system,c1,c2,srho_s,mixing)
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
  
end subroutine simple_mixing

!===================================================================================================================================

subroutine wrapper_broyden(ng,system,srho_s,mst,ifmst,iter,mixing)
  use structures, only: s_rgrid,s_dft_system,s_scalar,s_mixing
  use salmon_parallel, only: nproc_group_global
  use broyden_sub
  implicit none
  type(s_rgrid) :: ng
  type(s_dft_system),intent(in) :: system
  type(s_scalar),intent(inout) :: srho_s(system%nspin)
  integer,intent(in) :: mst(2),ifmst(2)
  integer,intent(in) :: iter
  type(s_mixing),intent(inout) :: mixing
  integer :: ix,iy,iz,is
  integer :: i
  real(8) :: vecr(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
  real(8) :: vecr_in(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),mixing%num_rho_stock+1)
  real(8) :: vecr_out(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),mixing%num_rho_stock+1)

  if(system%nspin==1)then

    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      vecr(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do

    do i=1,mixing%num_rho_stock+1
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        vecr_in(ix,iy,iz,i)=mixing%srho_in(i)%f(ix,iy,iz)
        vecr_out(ix,iy,iz,i)=mixing%srho_out(i)%f(ix,iy,iz)
      end do
      end do
      end do
    end do

    call broyden(vecr,vecr_in,vecr_out,ng%num(1)*ng%num(2)*ng%num(3),  &
                 iter,mixing%num_rho_stock,mixing%num_rho_stock,nproc_group_global)
  
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      srho_s(1)%f(ix,iy,iz)= vecr(ix,iy,iz)
    end do
    end do
    end do

    do i=1,mixing%num_rho_stock+1
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        mixing%srho_in(i)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
        mixing%srho_out(i)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
      end do
      end do
      end do
    end do

  else if(system%nspin==2)then
    
    do is=1,2
      if(ifmst(is)>=1.and.mst(is)>=1)then
        do iz=ng%is(3),ng%ie(3)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          vecr(ix,iy,iz)=srho_s(is)%f(ix,iy,iz)
        end do
        end do
        end do
  
        do i=1,mixing%num_rho_stock+1
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            vecr_in(ix,iy,iz,i)=mixing%srho_s_in(i,is)%f(ix,iy,iz)
            vecr_out(ix,iy,iz,i)=mixing%srho_s_out(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
        end do

        call broyden(vecr,vecr_in, vecr_out, ng%num(1)*ng%num(2)*ng%num(3),  &
                     iter,mixing%num_rho_stock,mixing%num_rho_stock,nproc_group_global)
  
        do iz=ng%is(3),ng%ie(3)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          srho_s(is)%f(ix,iy,iz)= vecr(ix,iy,iz)
        end do
        end do
        end do

        do i=1,mixing%num_rho_stock+1
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            mixing%srho_s_in(i,is)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
            mixing%srho_s_out(i,is)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
          end do
          end do
          end do
        end do
      end if
    end do

  end if

end subroutine wrapper_broyden

!===================================================================================================================================

subroutine init_mixing(nspin,ng,mixing)
  use structures
  implicit none
  integer      ,intent(in) :: nspin
  type(s_rgrid),intent(in) :: ng
  type(s_mixing)           :: mixing
  !
  integer :: i,j

  allocate(mixing%srho_in(1:mixing%num_rho_stock+1))
  allocate(mixing%srho_out(1:mixing%num_rho_stock+1))
  do i=1,mixing%num_rho_stock+1
    allocate(mixing%srho_in(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
    allocate(mixing%srho_out(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
    mixing%srho_in(i)%f(:,:,:) =0.d0
    mixing%srho_out(i)%f(:,:,:)=0.d0
  end do

  if(nspin==2)then
    allocate(mixing%srho_s_in(1:mixing%num_rho_stock+1,2))
    allocate(mixing%srho_s_out(1:mixing%num_rho_stock+1,2))
    do j=1,2
      do i=1,mixing%num_rho_stock+1
        allocate(mixing%srho_s_in(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
        allocate(mixing%srho_s_out(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
        mixing%srho_s_in(i,j)%f(:,:,:) =0.d0
        mixing%srho_s_out(i,j)%f(:,:,:)=0.d0
      end do
    end do
  end if

end subroutine init_mixing

!===================================================================================================================================

subroutine copy_density(Miter,nspin,ng,srho_s,mixing)
  use structures, only: s_rgrid, s_scalar, s_mixing
  implicit none
  integer       ,intent(in) :: Miter,nspin
  type(s_rgrid), intent(in) :: ng
  type(s_scalar),intent(in) :: srho_s(nspin)
  type(s_mixing),intent(inout) :: mixing
  !
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
    if(nspin==2)then
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

  if(nspin==2)then
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

end module mixing_sub
