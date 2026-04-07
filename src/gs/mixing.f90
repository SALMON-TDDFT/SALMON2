!
!  Copyright 2019-2020 SALMON developers
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
subroutine simple_mixing(mg,system,c1,c2,rho_s,mixing)
  use structures, only: s_rgrid, s_dft_system, s_scalar, s_mixing  
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  real(8),intent(in) :: c1,c2
  type(s_scalar),intent(inout) :: rho_s(system%nspin)
  type(s_mixing),intent(inout) :: mixing
  
  integer :: ix,iy,iz
  
  if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  elseif(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
      mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  !rho = c1*rho + c2*matmul( psi**2, occ )
  if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(1)%f(ix,iy,iz) = c1*mixing%rho_in(mixing%num_rho_stock)%f(ix,iy,iz) &
                              + c2*mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)
      mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else if(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(1)%f(ix,iy,iz) = c1*mixing%rho_s_in(mixing%num_rho_stock,1)%f(ix,iy,iz) &
                              + c2*mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)
      rho_s(2)%f(ix,iy,iz) = c1*mixing%rho_s_in(mixing%num_rho_stock,2)%f(ix,iy,iz) &
                              + c2*mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)
      mixing%rho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz) = rho_s(1)%f(ix,iy,iz)
      mixing%rho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz) = rho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  
  return
  
end subroutine simple_mixing
!===================================================================================================================================
subroutine simple_mixing_potential(mg,system,c1,c2,Vh,Vxc,mixing)
  use structures, only: s_rgrid, s_dft_system, s_scalar, s_mixing  
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  real(8),intent(in) :: c1,c2
  type(s_scalar),intent(inout) :: Vh,Vxc(system%nspin)
  type(s_mixing),intent(inout) :: mixing
  
  integer :: ix,iy,iz
  

  if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%Vh_out(mixing%num_rho_stock)%f(ix,iy,iz)=Vh%f(ix,iy,iz)
      mixing%Vxc_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=Vxc(1)%f(ix,iy,iz)
    end do
    end do
    end do

!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vh%f(ix,iy,iz) = c1*mixing%Vh_in(mixing%num_rho_stock)%f(ix,iy,iz) &
                      +c2*mixing%Vh_out(mixing%num_rho_stock)%f(ix,iy,iz)
      Vxc(1)%f(ix,iy,iz)= &
          c1*mixing%Vxc_in(mixing%num_rho_stock,1)%f(ix,iy,iz) &
         +c2*mixing%Vxc_out(mixing%num_rho_stock,1)%f(ix,iy,iz)

      mixing%Vh_in(mixing%num_rho_stock)%f(ix,iy,iz)=Vh%f(ix,iy,iz)
      mixing%Vxc_in(mixing%num_rho_stock,1)%f(ix,iy,iz)=Vxc(1)%f(ix,iy,iz)
    end do
    end do
    end do

  else if(system%nspin == 2)then  

!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%Vh_out(mixing%num_rho_stock)%f(ix,iy,iz)=Vh%f(ix,iy,iz)
      mixing%Vxc_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=Vxc(1)%f(ix,iy,iz)
      mixing%Vxc_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=Vxc(2)%f(ix,iy,iz)
    end do
    end do
    end do

!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vh%f(ix,iy,iz) = c1*mixing%Vh_in(mixing%num_rho_stock)%f(ix,iy,iz) &
                      +c2*mixing%Vh_out(mixing%num_rho_stock)%f(ix,iy,iz)
      Vxc(1)%f(ix,iy,iz)= &
          c1*mixing%Vxc_in(mixing%num_rho_stock,1)%f(ix,iy,iz) &
         +c2*mixing%Vxc_out(mixing%num_rho_stock,1)%f(ix,iy,iz)
      Vxc(2)%f(ix,iy,iz)= &
          c1*mixing%Vxc_in(mixing%num_rho_stock,2)%f(ix,iy,iz) &
         +c2*mixing%Vxc_out(mixing%num_rho_stock,2)%f(ix,iy,iz)

      mixing%Vh_in(mixing%num_rho_stock)%f(ix,iy,iz)=Vh%f(ix,iy,iz)
      mixing%Vxc_in(mixing%num_rho_stock,1)%f(ix,iy,iz)=Vxc(1)%f(ix,iy,iz)
      mixing%Vxc_in(mixing%num_rho_stock,2)%f(ix,iy,iz)=Vxc(2)%f(ix,iy,iz)
    end do
    end do
    end do

  end if
end subroutine simple_mixing_potential

!===================================================================================================================================

subroutine ensure_scalar_history_storage(mg,num_rho_stock,field_in,field_out)
  use structures, only: s_rgrid, s_scalar, allocate_scalar
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer, intent(in) :: num_rho_stock
  type(s_scalar),allocatable,intent(inout) :: field_in(:), field_out(:)
  integer :: i

  if (allocated(field_in)) return

  allocate(field_in(1:num_rho_stock+1))
  allocate(field_out(1:num_rho_stock+1))
  do i=1,num_rho_stock+1
    call allocate_scalar(mg, field_in(i))
    call allocate_scalar(mg, field_out(i))
    field_in(i)%f = 0d0
    field_out(i)%f = 0d0
  end do
end subroutine ensure_scalar_history_storage

!===================================================================================================================================

subroutine ensure_aux_mixing_storage(mg,mixing)
  use structures, only: s_rgrid, s_mixing, allocate_scalar, allocate_xc_aux_fields
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_mixing),intent(inout) :: mixing

  if (mixing%use_aux_tau) then
    call ensure_scalar_history_storage(mg,mixing%num_rho_stock,mixing%tau_in,mixing%tau_out)
    if (.not. allocated(mixing%tau_work%f)) then
      call allocate_scalar(mg, mixing%tau_work)
      mixing%tau_work%f = 0d0
    end if
  end if

  if (mixing%use_aux_j) then
    call ensure_scalar_history_storage(mg,mixing%num_rho_stock,mixing%jx_in,mixing%jx_out)
    call ensure_scalar_history_storage(mg,mixing%num_rho_stock,mixing%jy_in,mixing%jy_out)
    call ensure_scalar_history_storage(mg,mixing%num_rho_stock,mixing%jz_in,mixing%jz_out)
    if (.not. allocated(mixing%j_work(1)%f)) then
      call allocate_scalar(mg, mixing%j_work(1))
      mixing%j_work(1)%f = 0d0
    end if
    if (.not. allocated(mixing%j_work(2)%f)) then
      call allocate_scalar(mg, mixing%j_work(2))
      mixing%j_work(2)%f = 0d0
    end if
    if (.not. allocated(mixing%j_work(3)%f)) then
      call allocate_scalar(mg, mixing%j_work(3))
      mixing%j_work(3)%f = 0d0
    end if
  end if

  if (mixing%use_aux_tau .or. mixing%use_aux_j) then
    if ((mixing%aux_work%use_tau .neqv. mixing%use_aux_tau) .or. &
        (mixing%aux_work%use_j .neqv. mixing%use_aux_j)) then
      call allocate_xc_aux_fields(mg,mixing%use_aux_tau,mixing%use_aux_j,mixing%aux_work)
    end if
  end if
end subroutine ensure_aux_mixing_storage

!===================================================================================================================================

subroutine copy_scalar_history(mg,num_rho_stock,field_in,field_out)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid), intent(in) :: mg
  integer, intent(in) :: num_rho_stock
  type(s_scalar), intent(inout) :: field_in(1:num_rho_stock+1), field_out(1:num_rho_stock+1)
  integer :: iiter, ix, iy, iz

  do iiter=1,num_rho_stock
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      field_in(iiter)%f(ix,iy,iz)=field_in(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do

  do iiter=1,num_rho_stock-1
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      field_out(iiter)%f(ix,iy,iz)=field_out(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do
end subroutine copy_scalar_history

!===================================================================================================================================

subroutine copy_aux_history(mg,mixing)
  use structures, only: s_rgrid, s_mixing
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_mixing), intent(inout) :: mixing

  if (mixing%use_aux_tau .and. mixing%tau_history_ready) then
    call ensure_aux_mixing_storage(mg, mixing)
    call copy_scalar_history(mg,mixing%num_rho_stock,mixing%tau_in,mixing%tau_out)
  end if

  if (mixing%use_aux_j .and. mixing%j_history_ready) then
    call ensure_aux_mixing_storage(mg, mixing)
    call copy_scalar_history(mg,mixing%num_rho_stock,mixing%jx_in,mixing%jx_out)
    call copy_scalar_history(mg,mixing%num_rho_stock,mixing%jy_in,mixing%jy_out)
    call copy_scalar_history(mg,mixing%num_rho_stock,mixing%jz_in,mixing%jz_out)
  end if
end subroutine copy_aux_history

!===================================================================================================================================

subroutine simple_mixing_scalar(mg,c1,c2,field,field_in,field_out,history_ready,num_rho_stock)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in) :: c1,c2
  type(s_scalar),intent(inout) :: field
  type(s_scalar),intent(inout) :: field_in(1:num_rho_stock+1), field_out(1:num_rho_stock+1)
  logical, intent(inout) :: history_ready
  integer, intent(in) :: num_rho_stock
  integer :: ix, iy, iz

  if (.not. history_ready) then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      field_out(num_rho_stock)%f(ix,iy,iz)=field%f(ix,iy,iz)
      field_in(num_rho_stock+1)%f(ix,iy,iz)=field%f(ix,iy,iz)
    end do
    end do
    end do
    history_ready = .true.
    return
  end if

!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    field_out(num_rho_stock)%f(ix,iy,iz)=field%f(ix,iy,iz)
    field%f(ix,iy,iz)=c1*field_in(num_rho_stock)%f(ix,iy,iz) &
                    + c2*field_out(num_rho_stock)%f(ix,iy,iz)
    field_in(num_rho_stock+1)%f(ix,iy,iz)=field%f(ix,iy,iz)
  end do
  end do
  end do
end subroutine simple_mixing_scalar

!===================================================================================================================================

subroutine simple_mixing_tau(mg,c1,c2,tau,mixing)
  use structures, only: s_rgrid, s_scalar, s_mixing
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in) :: c1,c2
  type(s_scalar),intent(inout) :: tau
  type(s_mixing),intent(inout) :: mixing

  call ensure_aux_mixing_storage(mg, mixing)
  call simple_mixing_scalar(mg,c1,c2,tau,mixing%tau_in,mixing%tau_out,mixing%tau_history_ready,mixing%num_rho_stock)
end subroutine simple_mixing_tau

!===================================================================================================================================

subroutine simple_mixing_j(mg,c1,c2,j,mixing)
  use structures, only: s_rgrid, s_scalar, s_mixing
  implicit none
  type(s_rgrid), intent(in) :: mg
  real(8), intent(in) :: c1, c2
  type(s_scalar), intent(inout) :: j(3)
  type(s_mixing), intent(inout) :: mixing

  call ensure_aux_mixing_storage(mg, mixing)
  call simple_mixing_scalar(mg,c1,c2,j(1),mixing%jx_in,mixing%jx_out,mixing%j_history_ready,mixing%num_rho_stock)
  call simple_mixing_scalar(mg,c1,c2,j(2),mixing%jy_in,mixing%jy_out,mixing%j_history_ready,mixing%num_rho_stock)
  call simple_mixing_scalar(mg,c1,c2,j(3),mixing%jz_in,mixing%jz_out,mixing%j_history_ready,mixing%num_rho_stock)
  mixing%j_history_ready = .true.
end subroutine simple_mixing_j

!===================================================================================================================================

subroutine damp_scalar_output(mg,field_in,field_raw,field_out,mixrate)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_scalar), intent(in) :: field_in, field_raw
  type(s_scalar), intent(inout) :: field_out
  real(8), intent(in) :: mixrate
  integer :: ix, iy, iz

!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    field_out%f(ix,iy,iz) = field_in%f(ix,iy,iz) + mixrate * (field_raw%f(ix,iy,iz) - field_in%f(ix,iy,iz))
  end do
  end do
  end do
end subroutine damp_scalar_output

!===================================================================================================================================

subroutine copy_scalar_output(mg,field,field_out)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_scalar), intent(in) :: field
  type(s_scalar), intent(inout) :: field_out
  integer :: ix, iy, iz

!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    field_out%f(ix,iy,iz) = field%f(ix,iy,iz)
  end do
  end do
  end do
end subroutine copy_scalar_output

!===================================================================================================================================

subroutine ensure_aux_vector_workspace(mixing,vec_size)
  use structures, only: s_mixing
  implicit none
  type(s_mixing), intent(inout) :: mixing
  integer, intent(in) :: vec_size

  if (mixing%aux_vector_size == vec_size .and. allocated(mixing%aux_vec)) return

  if (allocated(mixing%aux_vec)) deallocate(mixing%aux_vec)
  if (allocated(mixing%aux_vec_in)) deallocate(mixing%aux_vec_in)
  if (allocated(mixing%aux_vec_out)) deallocate(mixing%aux_vec_out)
  if (allocated(mixing%aux_vec_x)) deallocate(mixing%aux_vec_x)
  if (allocated(mixing%aux_vec_y)) deallocate(mixing%aux_vec_y)

  allocate(mixing%aux_vec(vec_size))
  allocate(mixing%aux_vec_in(vec_size,mixing%num_rho_stock+1))
  allocate(mixing%aux_vec_out(vec_size,mixing%num_rho_stock+1))
  allocate(mixing%aux_vec_x(vec_size))
  allocate(mixing%aux_vec_y(vec_size))
  mixing%aux_vec = 0d0
  mixing%aux_vec_in = 0d0
  mixing%aux_vec_out = 0d0
  mixing%aux_vec_x = 0d0
  mixing%aux_vec_y = 0d0
  mixing%aux_vector_size = vec_size
end subroutine ensure_aux_vector_workspace

!===================================================================================================================================

subroutine pack_aux_mixing_vector(mg,rho,vec,tau_scale,j_scale,tau,jx,jy,jz)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_scalar), intent(in) :: rho
  real(8), intent(out) :: vec(:)
  real(8), intent(in) :: tau_scale, j_scale
  type(s_scalar), intent(in), optional :: tau, jx, jy, jz
  integer :: ix, iy, iz, i, nxyz, nxy, offset

  nxy = mg%num(1)*mg%num(2)
  nxyz = mg%num(1)*mg%num(2)*mg%num(3)

!$omp parallel do collapse(2) private(ix,iy,iz,i,offset)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    i = (iz-mg%is(3))*nxy + (iy-mg%is(2))*mg%num(1) + (ix-mg%is(1)) + 1
    vec(i) = rho%f(ix,iy,iz)
    offset = nxyz
    if (present(tau)) then
      vec(offset+i) = tau_scale * tau%f(ix,iy,iz)
      offset = offset + nxyz
    end if
    if (present(jx)) then
      vec(offset+i) = j_scale * jx%f(ix,iy,iz)
      offset = offset + nxyz
      vec(offset+i) = j_scale * jy%f(ix,iy,iz)
      offset = offset + nxyz
      vec(offset+i) = j_scale * jz%f(ix,iy,iz)
    end if
  end do
  end do
  end do
!$omp end parallel do
end subroutine pack_aux_mixing_vector

!===================================================================================================================================

subroutine unpack_aux_mixing_vector(mg,vec,tau_scale,j_scale,rho,tau,jx,jy,jz)
  use structures, only: s_rgrid, s_scalar
  implicit none
  type(s_rgrid), intent(in) :: mg
  real(8), intent(in) :: vec(:)
  real(8), intent(in) :: tau_scale, j_scale
  type(s_scalar), intent(inout) :: rho
  type(s_scalar), intent(inout), optional :: tau, jx, jy, jz
  integer :: ix, iy, iz, i, nxyz, nxy, offset

  nxy = mg%num(1)*mg%num(2)
  nxyz = mg%num(1)*mg%num(2)*mg%num(3)

!$omp parallel do collapse(2) private(ix,iy,iz,i,offset)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    i = (iz-mg%is(3))*nxy + (iy-mg%is(2))*mg%num(1) + (ix-mg%is(1)) + 1
    rho%f(ix,iy,iz) = vec(i)
    offset = nxyz
    if (present(tau)) then
      tau%f(ix,iy,iz) = vec(offset+i) / tau_scale
      offset = offset + nxyz
    end if
    if (present(jx)) then
      jx%f(ix,iy,iz) = vec(offset+i) / j_scale
      offset = offset + nxyz
      jy%f(ix,iy,iz) = vec(offset+i) / j_scale
      offset = offset + nxyz
      jz%f(ix,iy,iz) = vec(offset+i) / j_scale
    end if
  end do
  end do
  end do
!$omp end parallel do
end subroutine unpack_aux_mixing_vector

!===================================================================================================================================

subroutine wrapper_broyden(comm,mg,system,rho_s,tau,j,iter,mixing)
  use structures, only: s_rgrid,s_dft_system,s_scalar,s_mixing
  use broyden_sub
  implicit none
  integer,intent(in) :: comm
  type(s_rgrid) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_scalar),intent(inout) :: rho_s(system%nspin)
  type(s_scalar),intent(inout),optional :: tau
  type(s_scalar),intent(inout),optional :: j(3)
  integer,intent(in) :: iter
  type(s_mixing),intent(inout) :: mixing
  integer :: ix,iy,iz,is
  integer :: i
  integer :: nxyz, vec_size
  real(8) :: tau_scale, j_scale
  logical :: use_tau_block, use_j_block, seed_aux_only
  real(8) :: vecr(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: vecr_in(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),mixing%num_rho_stock+1)
  real(8) :: vecr_out(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),mixing%num_rho_stock+1)

#ifdef USE_OPENACC
!$acc data copyin(mixing)
#endif

  use_tau_block = mixing%use_aux_tau .and. present(tau) .and. system%nspin == 1
  use_j_block = mixing%use_aux_j .and. present(j) .and. system%nspin == 1

  if (use_tau_block .or. use_j_block) then
    call ensure_aux_mixing_storage(mg, mixing)
    seed_aux_only = (use_tau_block .and. .not. mixing%tau_history_ready) .or. &
                    (use_j_block .and. .not. mixing%j_history_ready)
    if (use_tau_block .and. .not. mixing%tau_history_ready) then
      call simple_mixing_tau(mg,1d0,0d0,tau,mixing)
    end if
    if (use_j_block .and. .not. mixing%j_history_ready) then
      call simple_mixing_j(mg,1d0,0d0,j,mixing)
    end if
    if (.not. seed_aux_only) then
      nxyz = mg%num(1)*mg%num(2)*mg%num(3)
      tau_scale = 1d0
      j_scale = 1d0
      if (use_tau_block) tau_scale = sqrt(mixing%tau_metric_weight)
      if (use_j_block) j_scale = sqrt(mixing%j_metric_weight)
      vec_size = nxyz
      if (use_tau_block) vec_size = vec_size + nxyz
      if (use_j_block) vec_size = vec_size + 3*nxyz
      call ensure_aux_vector_workspace(mixing,vec_size)

      if (use_tau_block) call copy_scalar_output(mg,tau,mixing%tau_out(mixing%num_rho_stock))
      if (use_j_block) then
        call damp_scalar_output(mg,mixing%jx_in(mixing%num_rho_stock),j(1),mixing%jx_out(mixing%num_rho_stock),mixing%j_mixrate)
        call damp_scalar_output(mg,mixing%jy_in(mixing%num_rho_stock),j(2),mixing%jy_out(mixing%num_rho_stock),mixing%j_mixrate)
        call damp_scalar_output(mg,mixing%jz_in(mixing%num_rho_stock),j(3),mixing%jz_out(mixing%num_rho_stock),mixing%j_mixrate)
      end if

      if (use_tau_block .and. use_j_block) then
        call pack_aux_mixing_vector(mg,rho_s(1),mixing%aux_vec,tau_scale,j_scale, &
             tau=mixing%tau_out(mixing%num_rho_stock),jx=mixing%jx_out(mixing%num_rho_stock), &
             jy=mixing%jy_out(mixing%num_rho_stock),jz=mixing%jz_out(mixing%num_rho_stock))
        do i=1,mixing%num_rho_stock+1
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i),tau_scale,j_scale, &
               tau=mixing%tau_in(i),jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i),tau_scale,j_scale, &
               tau=mixing%tau_out(i),jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        end do
      else if (use_tau_block) then
        call pack_aux_mixing_vector(mg,rho_s(1),mixing%aux_vec,tau_scale,j_scale,tau=mixing%tau_out(mixing%num_rho_stock))
        do i=1,mixing%num_rho_stock+1
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i),tau_scale,j_scale,tau=mixing%tau_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i),tau_scale,j_scale,tau=mixing%tau_out(i))
        end do
      else
        call pack_aux_mixing_vector(mg,rho_s(1),mixing%aux_vec,tau_scale,j_scale, &
             jx=mixing%jx_out(mixing%num_rho_stock),jy=mixing%jy_out(mixing%num_rho_stock),jz=mixing%jz_out(mixing%num_rho_stock))
        do i=1,mixing%num_rho_stock+1
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i),tau_scale,j_scale, &
               jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i),tau_scale,j_scale, &
               jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        end do
      end if

      call broyden(mixing%alpha_mb,mixing%aux_vec,mixing%aux_vec_in,mixing%aux_vec_out,vec_size,iter, &
                   mixing%num_rho_stock,mixing%num_rho_stock,comm,mixing%flag_mix_zero)

      if (use_tau_block .and. use_j_block) then
        call unpack_aux_mixing_vector(mg,mixing%aux_vec,tau_scale,j_scale,rho_s(1),tau=tau,jx=j(1),jy=j(2),jz=j(3))
        do i=1,mixing%num_rho_stock+1
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_in(:,i),tau_scale,j_scale,mixing%rho_in(i), &
               tau=mixing%tau_in(i),jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_out(:,i),tau_scale,j_scale,mixing%rho_out(i), &
               tau=mixing%tau_out(i),jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        end do
      else if (use_tau_block) then
        call unpack_aux_mixing_vector(mg,mixing%aux_vec,tau_scale,j_scale,rho_s(1),tau=tau)
        do i=1,mixing%num_rho_stock+1
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_in(:,i),tau_scale,j_scale,mixing%rho_in(i),tau=mixing%tau_in(i))
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_out(:,i),tau_scale,j_scale,mixing%rho_out(i),tau=mixing%tau_out(i))
        end do
      else
        call unpack_aux_mixing_vector(mg,mixing%aux_vec,tau_scale,j_scale,rho_s(1),jx=j(1),jy=j(2),jz=j(3))
        do i=1,mixing%num_rho_stock+1
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_in(:,i),tau_scale,j_scale,mixing%rho_in(i), &
               jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call unpack_aux_mixing_vector(mg,mixing%aux_vec_out(:,i),tau_scale,j_scale,mixing%rho_out(i), &
               jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        end do
      end if
      return
    end if
  end if

  if(system%nspin==1)then

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$omp parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       vecr(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do

#ifdef USE_OPENACC
!$acc parallel loop private(i,iz,iy,ix) collapse(3)
#else
!$omp parallel do private(i,iz,iy,ix) collapse(3)
#endif
    do i=1,mixing%num_rho_stock+1
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          vecr_in(ix,iy,iz,i) =mixing%rho_in(i)%f(ix,iy,iz)
          vecr_out(ix,iy,iz,i)=mixing%rho_out(i)%f(ix,iy,iz)
       end do
       end do
       end do
    end do

    call broyden(mixing%alpha_mb,vecr,vecr_in,vecr_out,mg%num(1)*mg%num(2)*mg%num(3),iter,    &
                 mixing%num_rho_stock,mixing%num_rho_stock,comm,&
                 mixing%flag_mix_zero)

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$omp parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       rho_s(1)%f(ix,iy,iz)= vecr(ix,iy,iz)
    end do
    end do
    end do

#ifdef USE_OPENACC
!$acc parallel loop private(i,iz,iy,ix) collapse(3)
#else
!$omp parallel do private(i,iz,iy,ix) collapse(3)
#endif
    do i=1,mixing%num_rho_stock+1
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          mixing%rho_in(i)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
          mixing%rho_out(i)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
       end do
       end do
       end do
    end do

  else if(system%nspin==2)then
    
    do is=1,2
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         vecr(ix,iy,iz)=rho_s(is)%f(ix,iy,iz)
      end do
      end do
      end do
  
!$omp parallel do private(i,iz,iy,ix) collapse(3)
      do i=1,mixing%num_rho_stock+1
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           vecr_in(ix,iy,iz,i)=mixing%rho_s_in(i,is)%f(ix,iy,iz)
           vecr_out(ix,iy,iz,i)=mixing%rho_s_out(i,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do

      call broyden(mixing%alpha_mb,vecr,vecr_in, vecr_out, mg%num(1)*mg%num(2)*mg%num(3),iter,  &
                   mixing%num_rho_stock,mixing%num_rho_stock,comm,&
                   mixing%flag_mix_zero )

!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         rho_s(is)%f(ix,iy,iz)= vecr(ix,iy,iz)
      end do
      end do
      end do

!$omp parallel do private(i,iz,iy,ix) collapse(3)
      do i=1,mixing%num_rho_stock+1
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           mixing%rho_s_in(i,is)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
           mixing%rho_s_out(i,is)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
        end do
        end do
        end do
      end do
    end do

  end if

#ifdef USE_OPENACC
!$acc end data
#endif

end subroutine wrapper_broyden

!===================================================================================================================================

subroutine pulay(mg,info,system,rho_s,tau,j,iter,mixing)
  use salmon_global, only: nmemory_p
  use structures, only: s_rgrid,s_parallel_info,s_dft_system,s_scalar,s_mixing,allocate_scalar,deallocate_scalar
  use communication, only: comm_summation
  implicit none
  type(s_rgrid),intent(in)            :: mg
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in)       :: system
  type(s_scalar),intent(inout)        :: rho_s(system%nspin)
  type(s_scalar),intent(inout),optional :: tau
  type(s_scalar),intent(inout),optional :: j(3)
  integer,intent(in)                  :: iter
  type(s_mixing),intent(inout)        :: mixing
  integer :: nsize
  integer, allocatable :: ipiv(:)
  type(s_scalar) :: x,y
  real(8), allocatable :: b1(:)
  real(8), allocatable :: a1(:,:)
  real(8), allocatable :: a0(:,:)
  integer :: i,k,i0,j0
  integer :: is,ix,iy,iz
  integer :: ierr
  real(8) :: ss
  real(8) :: rc
  integer :: nxyz, ivec, vec_size, offset
  real(8) :: tau_scale, j_scale
  logical :: use_tau_block, use_j_block

  use_tau_block = mixing%use_aux_tau .and. present(tau) .and. system%nspin == 1
  use_j_block = mixing%use_aux_j .and. present(j) .and. system%nspin == 1

  if(iter==1.or.nmemory_p==1)then

    call simple_mixing(mg,system,1.d0-mixing%beta_p,mixing%beta_p,rho_s,mixing)
    if (use_tau_block) then
      call simple_mixing_tau(mg,1.d0-mixing%beta_p,mixing%beta_p,tau,mixing)
    end if
    if (use_j_block) then
      call simple_mixing_j(mg,1.d0-mixing%j_mixrate,mixing%j_mixrate,j,mixing)
    end if

  else
!pulay mixing

    if (use_tau_block .or. use_j_block) then
      call ensure_aux_mixing_storage(mg, mixing)

!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
      end do
      end do
      end do

      if (use_tau_block) call copy_scalar_output(mg,tau,mixing%tau_out(mixing%num_rho_stock))
      if (use_j_block) then
        call damp_scalar_output(mg,mixing%jx_in(mixing%num_rho_stock),j(1),mixing%jx_out(mixing%num_rho_stock),mixing%j_mixrate)
        call damp_scalar_output(mg,mixing%jy_in(mixing%num_rho_stock),j(2),mixing%jy_out(mixing%num_rho_stock),mixing%j_mixrate)
        call damp_scalar_output(mg,mixing%jz_in(mixing%num_rho_stock),j(3),mixing%jz_out(mixing%num_rho_stock),mixing%j_mixrate)
      end if

      if ((use_tau_block .and. .not. mixing%tau_history_ready) .or. (use_j_block .and. .not. mixing%j_history_ready)) then
        if (use_tau_block) call simple_mixing_tau(mg,1.d0-mixing%beta_p,mixing%beta_p,tau,mixing)
        if (use_j_block) call simple_mixing_j(mg,1.d0-mixing%j_mixrate,mixing%j_mixrate,j,mixing)
        return
      end if

      if(iter>=nmemory_p)then
        nsize=nmemory_p
      else
        nsize=iter
      end if

      nxyz = mg%num(1)*mg%num(2)*mg%num(3)
      tau_scale = 1d0
      j_scale = 1d0
      if (use_tau_block) tau_scale = sqrt(mixing%tau_metric_weight)
      if (use_j_block) j_scale = sqrt(mixing%j_metric_weight)
      vec_size = nxyz
      if (use_tau_block) vec_size = vec_size + nxyz
      if (use_j_block) vec_size = vec_size + 3*nxyz
      call ensure_aux_vector_workspace(mixing,vec_size)

      allocate(ipiv(nsize))
      allocate(b1(nsize))
      allocate(a1(nsize,nsize))
      allocate(a0(nsize,nsize))

      b1(:) = 0.d0
      a1(:,:) = 0.d0
      a0(:,:) = 0.d0
      mixing%aux_vec_x(:) = 0.d0
      mixing%aux_vec_y(:) = 0.d0

      do i0=1,nsize
        i=mixing%num_rho_stock-nsize+i0
        if (use_tau_block .and. use_j_block) then
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i0),tau_scale,j_scale, &
               tau=mixing%tau_in(i),jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i0),tau_scale,j_scale, &
               tau=mixing%tau_out(i),jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        else if (use_tau_block) then
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i0),tau_scale,j_scale,tau=mixing%tau_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i0),tau_scale,j_scale,tau=mixing%tau_out(i))
        else
          call pack_aux_mixing_vector(mg,mixing%rho_in(i),mixing%aux_vec_in(:,i0),tau_scale,j_scale, &
               jx=mixing%jx_in(i),jy=mixing%jy_in(i),jz=mixing%jz_in(i))
          call pack_aux_mixing_vector(mg,mixing%rho_out(i),mixing%aux_vec_out(:,i0),tau_scale,j_scale, &
               jx=mixing%jx_out(i),jy=mixing%jy_out(i),jz=mixing%jz_out(i))
        end if
      end do

      do j0=1,nsize
      do i0=j0,nsize
        ss = dot_product(mixing%aux_vec_out(:,i0)-mixing%aux_vec_in(:,i0), mixing%aux_vec_out(:,j0)-mixing%aux_vec_in(:,j0))
        a0(i0,j0)=ss
        a0(j0,i0)=ss
      end do
      end do

      call comm_summation(a0,a1,nsize*nsize,info%icomm_r)
      b1(1:nsize) = 1.d0
      call dgesv(nsize,1,a1,nsize,ipiv,b1,nsize,ierr)

      rc=1.d0/sum(b1(1:nsize))
      b1(1:nsize)=rc*b1(1:nsize)

      do i0=1,nsize
        mixing%aux_vec_x(:) = mixing%aux_vec_x(:) + b1(i0) * mixing%aux_vec_in(:,i0)
        mixing%aux_vec_y(:) = mixing%aux_vec_y(:) + b1(i0) * mixing%aux_vec_out(:,i0)
      end do

      mixing%aux_vec(:) = mixing%aux_vec_x(:) + mixing%beta_p * (mixing%aux_vec_y(:) - mixing%aux_vec_x(:))

      ivec = 0
!$omp parallel do private(iz,iy,ix,ivec) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        ivec = (iz-mg%is(3))*mg%num(1)*mg%num(2) + (iy-mg%is(2))*mg%num(1) + (ix-mg%is(1)) + 1
        mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = max(1.d-20, mixing%aux_vec(ivec))
        rho_s(1)%f(ix,iy,iz) = mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
        offset = nxyz
        if (use_tau_block) then
          mixing%tau_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = mixing%aux_vec(offset+ivec) / tau_scale
          tau%f(ix,iy,iz) = mixing%tau_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
          offset = offset + nxyz
        end if
        if (use_j_block) then
          mixing%jx_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = mixing%aux_vec(offset+ivec) / j_scale
          j(1)%f(ix,iy,iz) = mixing%jx_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
          offset = offset + nxyz
          mixing%jy_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = mixing%aux_vec(offset+ivec) / j_scale
          j(2)%f(ix,iy,iz) = mixing%jy_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
          offset = offset + nxyz
          mixing%jz_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = mixing%aux_vec(offset+ivec) / j_scale
          j(3)%f(ix,iy,iz) = mixing%jz_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
        end if
      end do
      end do
      end do

      deallocate(a0,a1,b1,ipiv)
      return
    end if

    if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
      end do
      end do
      end do
    elseif(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
        mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
      end do
      end do
      end do
    end if
  
  !rho = c1*rho + c2*matmul( psi**2, occ )

    if(iter>=nmemory_p)then
      nsize=nmemory_p
    else
      nsize=iter
    end if
  
    allocate( ipiv(nsize) )
    allocate( b1(nsize) )
    allocate( a1(nsize,nsize) )
    allocate( a0(nsize,nsize) )
    call allocate_scalar(mg,x)
    call allocate_scalar(mg,y)
  
  
    b1(:)   = 0.d0
    a1(:,:) = 0.d0
    a0(:,:) = 0.d0
  
    do j0=1 ,nsize
    do i0=j0,nsize
      i=mixing%num_rho_stock-nsize+i0
      k=mixing%num_rho_stock-nsize+j0
      ss=0.d0
      if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2) reduction(+:ss)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          ss=ss+(mixing%rho_out(i)%f(ix,iy,iz)-mixing%rho_in(i)%f(ix,iy,iz))* &
                (mixing%rho_out(k)%f(ix,iy,iz)-mixing%rho_in(k)%f(ix,iy,iz))
        end do
        end do
        end do
      else
!$omp parallel do private(ix,iy,iz) collapse(3) reduction(+:ss)
        do is=1,system%nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ss=ss+(mixing%rho_s_out(i,is)%f(ix,iy,iz)-mixing%rho_s_in(i,is)%f(ix,iy,iz))* &
                  (mixing%rho_s_out(k,is)%f(ix,iy,iz)-mixing%rho_s_in(k,is)%f(ix,iy,iz))
          end do
          end do
          end do
        end do
      end if
      a0(i0,j0)=ss
      a0(j0,i0)=ss
    end do
    end do

    call comm_summation(a0,a1,nsize*nsize,info%icomm_r)

    b1(1:nsize) = 1.d0
  
    call dgesv(nsize,1,a1,nsize,ipiv,b1,nsize,ierr)

    rc=1.d0/sum( b1(1:nsize) )
    b1(1:nsize)=rc*b1(1:nsize)
 
    do is=1,system%nspin
!$omp parallel do private(ix,iy,iz) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        x%f(ix,iy,iz)=0.d0
        y%f(ix,iy,iz)=0.d0
      end do
      end do
      end do
  
      do i0=1,nsize
        i=mixing%num_rho_stock-nsize+i0
        if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            x%f(ix,iy,iz)=x%f(ix,iy,iz)+b1(i0)*mixing%rho_in(i)%f(ix,iy,iz)
            y%f(ix,iy,iz)=y%f(ix,iy,iz)+b1(i0)*mixing%rho_out(i)%f(ix,iy,iz)
          end do
          end do
          end do
        else
!$omp parallel do private(ix,iy,iz) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            x%f(ix,iy,iz)=x%f(ix,iy,iz)+b1(i0)*mixing%rho_s_in(i,is)%f(ix,iy,iz)
            y%f(ix,iy,iz)=y%f(ix,iy,iz)+b1(i0)*mixing%rho_s_out(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do
  
      if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = max(1.d-20,  &
                                                                   x%f(ix,iy,iz) + mixing%beta_p*( y%f(ix,iy,iz)-x%f(ix,iy,iz) ))
          rho_s(is)%f(ix,iy,iz) = mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
        end do
        end do
        end do
      else
!$omp parallel do private(ix,iy,iz) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_in(mixing%num_rho_stock+1,is)%f(ix,iy,iz) = max(1.d-20,  &
                                                                     x%f(ix,iy,iz) + mixing%beta_p*( y%f(ix,iy,iz)-x%f(ix,iy,iz) ))
          rho_s(is)%f(ix,iy,iz) = mixing%rho_s_in(mixing%num_rho_stock+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end if
  
    end do

    deallocate( a0,a1,b1,ipiv )
    call deallocate_scalar(x)
    call deallocate_scalar(y)

  end if

end subroutine

!===================================================================================================================================

subroutine init_mixing(nspin,mg,mixing)
  use salmon_global, only: mixrate,alpha_mb,beta_p,yn_aux_mixing,tau_mixrate,tau_metric_weight,j_mixrate,j_metric_weight
  use structures
  implicit none
  integer      ,intent(in) :: nspin
  type(s_rgrid),intent(in) :: mg
  type(s_mixing)           :: mixing
  !
  integer :: i,j

  mixing%mixrate=mixrate
  mixing%alpha_mb=alpha_mb
  mixing%beta_p=beta_p
  mixing%use_aux_mixing = (yn_aux_mixing == 'y')
  mixing%use_aux_tau = mixing%use_aux_mixing .and. (tau_mixrate > 0d0)
  mixing%use_aux_j = mixing%use_aux_mixing .and. (j_mixrate > 0d0)
  mixing%tau_history_ready = .false.
  mixing%j_history_ready = .false.
  mixing%tau_mixrate = tau_mixrate
  mixing%tau_metric_weight = tau_metric_weight
  mixing%j_mixrate = j_mixrate
  mixing%j_metric_weight = j_metric_weight
  mixing%aux_vector_size = 0
  mixing%convergence_value_prev=1.d10

  allocate(mixing%rho_in(1:mixing%num_rho_stock+1))
  allocate(mixing%rho_out(1:mixing%num_rho_stock+1))
  do i=1,mixing%num_rho_stock+1
    allocate(mixing%rho_in(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(mixing%rho_out(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    mixing%rho_in(i)%f(:,:,:) =0.d0
    mixing%rho_out(i)%f(:,:,:)=0.d0
  end do

  if(nspin==2)then
    allocate(mixing%rho_s_in(1:mixing%num_rho_stock+1,2))
    allocate(mixing%rho_s_out(1:mixing%num_rho_stock+1,2))
    do j=1,2
      do i=1,mixing%num_rho_stock+1
        allocate(mixing%rho_s_in(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
        allocate(mixing%rho_s_out(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
        mixing%rho_s_in(i,j)%f(:,:,:) =0.d0
        mixing%rho_s_out(i,j)%f(:,:,:)=0.d0
      end do
    end do
  end if

! allocate arrays for potential mixing
  allocate(mixing%Vh_in(1:mixing%num_rho_stock+1))
  allocate(mixing%Vh_out(1:mixing%num_rho_stock+1))
  allocate(mixing%Vxc_in(1:mixing%num_rho_stock+1,nspin))
  allocate(mixing%Vxc_out(1:mixing%num_rho_stock+1,nspin))

  do i=1,mixing%num_rho_stock+1
    allocate(mixing%Vh_in(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(mixing%Vh_out(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    mixing%Vh_in(i)%f(:,:,:)=0d0
    mixing%Vh_out(i)%f(:,:,:)=0d0
    do j = 1,nspin
      allocate(mixing%Vxc_in(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(mixing%Vxc_out(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      mixing%Vxc_in(i,j)%f(:,:,:)=0d0
      mixing%Vxc_out(i,j)%f(:,:,:)=0d0
    end do
  end do


  mixing%flag_mix_zero=.false.

end subroutine init_mixing

!===================================================================================================================================

subroutine copy_density(Miter,nspin,mg,rho_s,mixing)
  use structures, only: s_rgrid, s_scalar, s_mixing
  implicit none
  integer       ,intent(in) :: Miter,nspin
  type(s_rgrid), intent(in) :: mg
  type(s_scalar),intent(in) :: rho_s(nspin)
  type(s_mixing),intent(inout) :: mixing
  !
  integer :: iiter
  integer :: is
  integer :: ix,iy,iz

#ifdef USE_OPENACC
!$acc data copyin(mixing)
#endif

  if(Miter==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
    if(nspin==2)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
        mixing%rho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
      end do
      end do
      end do
    end if
  end if

  do iiter=1,mixing%num_rho_stock
#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$OMP parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_in(iiter)%f(ix,iy,iz)=mixing%rho_in(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do
  do iiter=1,mixing%num_rho_stock-1
#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$OMP parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_out(iiter)%f(ix,iy,iz)=mixing%rho_out(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do

  if(nspin==2)then
    do iiter=1,mixing%num_rho_stock
      do is=1,2
!$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_in(iiter,is)%f(ix,iy,iz)=mixing%rho_s_in(iiter+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
    do iiter=1,mixing%num_rho_stock-1
      do is=1,2
!$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_out(iiter,is)%f(ix,iy,iz)=mixing%rho_s_out(iiter+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
  end if

#ifdef USE_OPENACC
!$acc end data
#endif

end subroutine copy_density

!===================================================================================================================================
subroutine check_mixing_half(Miter,convergence_value,mixing)
  use salmon_global, only: method_mixing,update_mixing_ratio
  use structures, only: s_mixing
  use parallelization, only: nproc_id_global, nproc_group_global
  use communication, only: comm_is_root, comm_bcast
  implicit none
  integer, intent(in) :: Miter
  real(8), intent(in) :: convergence_value
  type(s_mixing), intent(inout) :: mixing
  integer :: icheck
  
  if(comm_is_root(nproc_id_global)) then
    if(convergence_value > update_mixing_ratio * mixing%convergence_value_prev)then
      icheck=1
    else
      icheck=0
    end if
  end if

  call comm_bcast(icheck,nproc_group_global)

  if(icheck==1)then
    select case(method_mixing)
    case('simple')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" mixrate decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%mixrate, mixing%mixrate*0.5d0, Miter
      end if
      mixing%mixrate=mixing%mixrate*0.5d0
      if (mixing%use_aux_tau) mixing%tau_mixrate = mixing%tau_mixrate*0.5d0
      if (mixing%use_aux_j) mixing%j_mixrate = mixing%j_mixrate*0.5d0
    case('simple_potential')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" mixrate decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%mixrate, mixing%mixrate*0.5d0, Miter
      end if
      mixing%mixrate=mixing%mixrate*0.5d0
      if (mixing%use_aux_tau) mixing%tau_mixrate = mixing%tau_mixrate*0.5d0
      if (mixing%use_aux_j) mixing%j_mixrate = mixing%j_mixrate*0.5d0
    case('broyden')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" alpha_mb decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%alpha_mb, mixing%alpha_mb*0.5d0, Miter 
      end if
      mixing%alpha_mb=mixing%alpha_mb*0.5d0
    case('pulay')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" beta_p decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%beta_p, mixing%beta_p*0.5d0, Miter 
      end if
      mixing%beta_p=mixing%beta_p*0.5d0
    end select
  end if

  mixing%convergence_value_prev=convergence_value

end subroutine check_mixing_half

subroutine reset_mixing_rate(mixing)
  use salmon_global, only: mixrate, alpha_mb, beta_p, tau_mixrate, j_mixrate
  use structures, only: s_mixing
  implicit none
  type(s_mixing), intent(inout) :: mixing

  mixing%mixrate  = mixrate
  mixing%alpha_mb = alpha_mb
  mixing%beta_p   = beta_p
  mixing%tau_mixrate = tau_mixrate
  mixing%j_mixrate = j_mixrate
  mixing%convergence_value_prev = 1.d10

end subroutine reset_mixing_rate


!===================================================================================================================================
end module mixing_sub
