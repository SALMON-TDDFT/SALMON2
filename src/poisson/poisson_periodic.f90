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
!===================================================================================================================================

module poisson_periodic
  implicit none

  character(16), save :: poisson_contract_context = 'UNSET'
  logical, save :: dumped_dc_ft = .false.
  logical, save :: dumped_dg_ft = .false.
  logical, save :: dumped_dc_ffte = .false.
  logical, save :: dumped_dg_ffte = .false.

contains

subroutine set_poisson_contract_context(context)
  implicit none
  character(*), intent(in) :: context
  poisson_contract_context = 'UNSET'
  poisson_contract_context = adjustl(trim(context))
end subroutine set_poisson_contract_context

subroutine poisson_ft(lg,mg,info,fg,rho,Vh,poisson,kernel_mode,omega)
  use structures
  use communication, only: comm_summation
  use math_constants, only : pi
  use nvtx
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  character(*), optional ,intent(in) :: kernel_mode
  real(8),      optional ,intent(in) :: omega
  !
  integer :: ix,iy,iz,kx,ky,kz
  logical :: use_hse_sr
  logical :: enable_contract_trace
  real(8) :: omega_eff, g2, sr_factor
  call nvtxStartRange('poission_ft', __LINE__)

  call resolve_poisson_kernel_mode(kernel_mode, omega, use_hse_sr, omega_eff)
  enable_contract_trace = poisson_contract_trace_enabled()

  if (enable_contract_trace) then
    call print_poisson_solver_contract('POISSON-FT-IN', lg, mg, info, rho)
  end if

#ifdef USE_OPENACC
!$acc data copyin(poisson, fg, rho, mg, lg, vh)
#endif

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    if (iz < mg%is(3) .or. iz > mg%ie(3)) then
      poisson%ff1z(ix,iy,iz) = 0d0
    else
      poisson%ff1z(ix,iy,iz) = dcmplx(rho%f(ix,iy,iz))
    end if
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1z = 0d0
!$omp end workshare

!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    poisson%ff1z(ix,iy,iz) = dcmplx(rho%f(ix,iy,iz))
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(kz,iy,ix)
  do kz = mg%is(3),mg%ie(3)
  do iy = lg%is(2),lg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    if (iy < mg%is(2) .or. iy > mg%ie(2)) then
      poisson%ff1y(ix,iy,kz) = 0d0
    else
      poisson%ff1y(ix,iy,kz) = sum(fg%egzc(kz,:)*poisson%ff2z(ix,iy,:))
    endif
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1y = 0d0
!$omp end workshare

!$OMP parallel do private(kz,iy,ix)
  do kz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    poisson%ff1y(ix,iy,kz) = sum(fg%egzc(kz,:)*poisson%ff2z(ix,iy,:))
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(kz,ky,ix)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do ix = lg%is(1),lg%ie(1)
    if (ix < mg%is(1) .or. ix > mg%ie(1)) then
      poisson%ff1x(ix,ky,kz) = 0d0
    else
      poisson%ff1x(ix,ky,kz) = sum(fg%egyc(ky,:)*poisson%ff2y(ix,:,kz))
    endif
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,ix)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    poisson%ff1x(ix,ky,kz) = sum(fg%egyc(ky,:)*poisson%ff2y(ix,:,kz))
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop private(kz,ky,kx)
#else
!$OMP parallel do private(kz,ky,kx)
#endif
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1x(kx,ky,kz) = sum(fg%egxc(kx,:)*poisson%ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end kernels
#endif

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(kz,ky,kx)
  do kz = lg%is(3),lg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    if (kz < mg%is(3) .or. kz > mg%ie(3)) then
      poisson%ff1z(kx,ky,kz) = 0.d0
    else
      poisson%zrhoG_ele(kx,ky,kz) = poisson%ff2x(kx,ky,kz)
      if (use_hse_sr) then
        g2 = fg%vec_G(1,kx,ky,kz)**2 + fg%vec_G(2,kx,ky,kz)**2 + fg%vec_G(3,kx,ky,kz)**2
        sr_factor = 1.d0 - exp(-g2/(4.d0*omega_eff*omega_eff))
        poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*sr_factor*poisson%ff2x(kx,ky,kz)
      else
        poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*poisson%ff2x(kx,ky,kz)
      end if
    end if
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1z = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,kx)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%zrhoG_ele(kx,ky,kz) = poisson%ff2x(kx,ky,kz)
    if (use_hse_sr) then
      g2 = fg%vec_G(1,kx,ky,kz)**2 + fg%vec_G(2,kx,ky,kz)**2 + fg%vec_G(3,kx,ky,kz)**2
      sr_factor = 1.d0 - exp(-g2/(4.d0*omega_eff*omega_eff))
      poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*sr_factor*poisson%ff2x(kx,ky,kz)
    else
      poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*poisson%ff2x(kx,ky,kz)
    end if
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

  if (enable_contract_trace) then
    call print_poisson_kernel_contract('POISSON-FT-KERNEL', mg, info, fg%coef, poisson%zrhoG_ele, poisson%ff2x, poisson%ff1z)
  end if

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(iz,ky,kx)
  do iz = mg%is(3),mg%ie(3)
  do ky = lg%is(2),lg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    if (ky < mg%is(2) .or. ky > mg%ie(2)) then
      poisson%ff1y(kx,ky,iz) = 0.d0
    else
      poisson%ff1y(kx,ky,iz) = sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
    end if
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1y = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,ky,kx)
  do iz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1y(kx,ky,iz) = sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(iz,iy,kx)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do kx = lg%is(1),lg%ie(1)
    if (kx < mg%is(1) .or. kx > mg%ie(1)) then
      poisson%ff1x(kx,iy,iz) = 0.d0
    else
      poisson%ff1x(kx,iy,iz) = sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
    endif
  end do
  end do
  end do
!$acc end kernels
#else
!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,iy,kx)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1x(kx,iy,iz) = sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
  end do
  end do
  end do
#endif

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

#ifdef USE_OPENACC
!$acc kernels
!$acc loop private(iz,iy,ix)
#else
!$OMP parallel do private(iz,iy,ix)
#endif
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    Vh%f(ix,iy,iz) = sum(fg%egx(:,ix)*poisson%ff2x(:,iy,iz))
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end kernels
#endif

#ifdef USE_OPENACC
!$acc end data
#endif
  if (enable_contract_trace) then
    call print_poisson_solver_contract('POISSON-FT-OUT', lg, mg, info, Vh)
  end if
  call nvtxEndRange
  return
end subroutine poisson_ft

subroutine poisson_ft_hse_sr(lg,mg,info,fg,rho,Vout,poisson,omega)
  use structures
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vout
  type(s_poisson)                    :: poisson
  real(8)               ,intent(in)  :: omega
  call poisson_ft(lg,mg,info,fg,rho,Vout,poisson,kernel_mode='hse_sr',omega=omega)
end subroutine poisson_ft_hse_sr

!===================================================================================================================================

subroutine poisson_ffte(lg,mg,info,fg,rho,Vh,poisson,kernel_mode,omega)
  use structures
  use communication, only: comm_summation
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  character(*), optional ,intent(in) :: kernel_mode
  real(8),      optional ,intent(in) :: omega
  !
  integer :: ix,iy,iz
  integer :: iiy,iiz,iix
  real(8) :: inv_lgnum3, g2, sr_factor, omega_eff
  logical :: use_hse_sr
  logical :: enable_contract_trace

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))
  call resolve_poisson_kernel_mode(kernel_mode, omega, use_hse_sr, omega_eff)
  enable_contract_trace = poisson_contract_trace_enabled()

  if (enable_contract_trace) then
    call print_poisson_solver_contract('POISSON-FFTE-IN', lg, mg, info, rho)
  end if

  poisson%b_ffte=0.d0
!$OMP parallel do private(iiz,iiy,ix) collapse(2)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    poisson%b_ffte(mg%is(1):mg%ie(1),iy,iz) = dcmplx(rho%f(mg%is(1):mg%ie(1),iiy,iiz))
  end do
  end do
  call comm_summation(poisson%b_ffte,poisson%a_ffte,size(poisson%a_ffte),info%icomm_x)

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                    info%isize_y,info%isize_z,-1, &
                    info%icomm_y,info%icomm_z)

  poisson%zrhoG_ele=0d0
!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,iiy,iiz,iix,g2,sr_factor) &
!$omp             shared(mg,lg,poisson,inv_lgnum3,fg,use_hse_sr,omega_eff)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    do ix=1,mg%num(1)
      iix=ix+mg%is(1)-1
      poisson%zrhoG_ele(iix,iiy,iiz) = poisson%b_ffte(iix,iy,iz)*inv_lgnum3
    end do
    do ix=1,lg%num(1)
      if (use_hse_sr) then
        g2 = fg%vec_G(1,ix,iiy,iiz)**2 + fg%vec_G(2,ix,iiy,iiz)**2 + fg%vec_G(3,ix,iiy,iiz)**2
        sr_factor = 1.d0 - exp(-g2/(4.d0*omega_eff*omega_eff))
        poisson%b_ffte(ix,iy,iz) = poisson%b_ffte(ix,iy,iz) * fg%coef(ix,iiy,iiz) * sr_factor
      else
        poisson%b_ffte(ix,iy,iz) = poisson%b_ffte(ix,iy,iz) * fg%coef(ix,iiy,iiz)
      end if
    end do
  end do
  end do
!$omp end parallel do

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(1),lg%num(2),lg%num(3), &
                    info%isize_y,info%isize_z,1, &
                    info%icomm_y,info%icomm_z)

!$OMP parallel do private(iiz,iiy) collapse(2)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    Vh%f(mg%is(1):mg%ie(1),iiy,iiz) = poisson%a_ffte(mg%is(1):mg%ie(1),iy,iz)
  end do
  end do

  if (enable_contract_trace) then
    call print_poisson_kernel_contract('POISSON-FFTE-KERNEL', mg, info, fg%coef, poisson%zrhoG_ele, poisson%b_ffte, poisson%b_ffte)
    call print_poisson_solver_contract('POISSON-FFTE-OUT', lg, mg, info, Vh)
  end if

  return
end subroutine poisson_ffte

subroutine poisson_ffte_hse_sr(lg,mg,info,fg,rho,Vout,poisson,omega)
  use structures
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vout
  type(s_poisson)                    :: poisson
  real(8)               ,intent(in)  :: omega
  call poisson_ffte(lg,mg,info,fg,rho,Vout,poisson,kernel_mode='hse_sr',omega=omega)
end subroutine poisson_ffte_hse_sr


#ifdef USE_FFTW
!===================================================================================================================================
subroutine poisson_fftw(lg,mg,info,fg,rho,Vh,poisson,kernel_mode,omega)
  use structures
  use communication, only: comm_summation
  use math_constants, only : pi
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  character(*), optional ,intent(in) :: kernel_mode
  real(8),      optional ,intent(in) :: omega
  !
  type(C_PTR) :: cdata_1, cdata_2
  type(C_PTR) :: plan_forward, plan_backward
  integer :: ix,iy,iz
  real(8) :: inv_lgnum3, g2, sr_factor, omega_eff
  logical :: use_hse_sr
  integer(C_INTPTR_T) :: alloc_local, local_N, local_k_offset
  complex(C_DOUBLE_COMPLEX), pointer :: fdata_1(:,:,:), fdata_2(:,:,:)
  integer(C_INTPTR_T) :: i, j, k
  integer(C_INTPTR_T) :: L, M, N

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))
  call resolve_poisson_kernel_mode(kernel_mode, omega, use_hse_sr, omega_eff)

  L = lg%num(1)
  M = lg%num(2)
  N = lg%num(3)

!$OMP parallel do private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    poisson%fftw1(ix,iy,iz)=dcmplx(rho%f(ix,iy,iz))
  end do
  end do
  end do
  
  call comm_summation(poisson%fftw1,poisson%fftw2,lg%num(1)*lg%num(2)*mg%num(3),info%icomm_xy)

!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N, M, L, &
     &                info%icomm_z, local_N, local_k_offset)
  cdata_1 = fftw_alloc_complex(alloc_local)
  cdata_2 = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cdata_1, fdata_1, [L, M, local_N])
  call c_f_pointer(cdata_2, fdata_2, [L, M, local_N])

!   create MPI plan for forward DFT (note dimension reversal)
  plan_forward = fftw_mpi_plan_dft_3d(N, M, L, fdata_1, fdata_2, &
     &           info%icomm_z, FFTW_FORWARD, FFTW_MEASURE)
!   create MPI plan for backward DFT (note dimension reversal)
  plan_backward = fftw_mpi_plan_dft_3d(N, M, L, fdata_2, fdata_1, &
     &            info%icomm_z, FFTW_BACKWARD, FFTW_MEASURE)

!$OMP parallel do private(i,j,k) collapse(2)
  do k = 1, local_N
  do j = 1, M
  do i = 1, L
    fdata_1(i, j, k) = poisson%fftw2(i, j, k+mg%is(3)-1)
  end do
  end do
  end do

  call fftw_mpi_execute_dft(plan_forward, fdata_1, fdata_2)

!$OMP parallel do private(i,j,k) collapse(2)
  do k = 1, local_N
  do j = 1, M
  do i = 1, L
    poisson%fftw2(i,j,k+mg%is(3)-1) = fdata_2(i,j,k)
  end do
  end do
  end do

  poisson%zrhoG_ele=0d0
!$OMP parallel do private(ix,iy,iz) collapse(2)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    poisson%zrhoG_ele(ix,iy,iz) = poisson%fftw2(ix,iy,iz)*inv_lgnum3
  end do
  end do
  end do
!$OMP parallel do private(ix,iy,iz) collapse(2)
  do iz = mg%is(3),mg%ie(3)
  do iy = lg%is(2),lg%ie(2)
  do ix = lg%is(1),lg%ie(1)
    if (use_hse_sr) then
      g2 = fg%vec_G(1,ix,iy,iz)**2 + fg%vec_G(2,ix,iy,iz)**2 + fg%vec_G(3,ix,iy,iz)**2
      sr_factor = 1.d0 - exp(-g2/(4.d0*omega_eff*omega_eff))
      poisson%fftw2(ix,iy,iz) = poisson%fftw2(ix,iy,iz)*fg%coef(ix,iy,iz)*sr_factor*inv_lgnum3
    else
      poisson%fftw2(ix,iy,iz) = poisson%fftw2(ix,iy,iz)*fg%coef(ix,iy,iz)*inv_lgnum3
    end if
  end do
  end do
  end do

!$OMP parallel do private(i,j,k) collapse(2)
  do k = 1, local_N
  do j = 1, M
  do i = 1, L
    fdata_2(i,j,k) = poisson%fftw2(i,j,k+mg%is(3)-1)
  end do
  end do
  end do

  call fftw_mpi_execute_dft(plan_backward, fdata_2, fdata_1)

!$OMP parallel do private(i,j,k) collapse(2)
  do k = 1, local_N
  do j = 1, M
  do i = 1, L
    poisson%fftw2(i,j,k+mg%is(3)-1) = fdata_1(i,j,k)
  end do
  end do
  end do
  
!$OMP parallel do private(ix,iy,iz) collapse(2)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    Vh%f(ix,iy,iz) = poisson%fftw2(ix,iy,iz)
  end do
  end do
  end do

  call fftw_destroy_plan(plan_forward)
  call fftw_destroy_plan(plan_backward)
  call fftw_free(cdata_1)
  call fftw_free(cdata_2)
 
end subroutine poisson_fftw

subroutine poisson_fftw_hse_sr(lg,mg,info,fg,rho,Vout,poisson,omega)
  use structures
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vout
  type(s_poisson)                    :: poisson
  real(8)               ,intent(in)  :: omega
  call poisson_fftw(lg,mg,info,fg,rho,Vout,poisson,kernel_mode='hse_sr',omega=omega)
end subroutine poisson_fftw_hse_sr

#endif

subroutine resolve_poisson_kernel_mode(kernel_mode, omega, use_hse_sr, omega_eff)
  implicit none
  character(*), optional, intent(in) :: kernel_mode
  real(8),      optional, intent(in) :: omega
  logical,                 intent(out) :: use_hse_sr
  real(8),                intent(out) :: omega_eff
  character(64) :: env_mode, env_omega
  integer :: env_status, ios
  real(8) :: env_omega_val

  ! Priority rule:
  ! explicit arguments (direct call or wrapper-provided) > environment variables > defaults
  use_hse_sr = .false.
  omega_eff = 0.11d0

  env_mode = ''
  call get_environment_variable('SALMON_POISSON_KERNEL_MODE', env_mode, status=env_status)
  if (env_status == 0) then
    if (trim(adjustl(env_mode)) == 'hse_sr') use_hse_sr = .true.
  end if

  env_omega = ''
  call get_environment_variable('SALMON_POISSON_KERNEL_OMEGA', env_omega, status=env_status)
  if (env_status == 0) then
    read(env_omega,*,iostat=ios) env_omega_val
    if (ios == 0) omega_eff = env_omega_val
  end if

  if (present(kernel_mode)) then
    if (trim(adjustl(kernel_mode)) == 'hse_sr') then
      use_hse_sr = .true.
    else
      use_hse_sr = .false.
    end if
  end if
  if (present(omega)) omega_eff = omega

  if (use_hse_sr .and. omega_eff <= 0.d0) then
    stop 'poisson_periodic: omega must be > 0 for hse_sr kernel'
  end if
end subroutine resolve_poisson_kernel_mode

logical function poisson_contract_trace_enabled()
  implicit none
  character(16) :: env_trace
  integer :: env_status
  poisson_contract_trace_enabled = .false.
  env_trace = ''
  call get_environment_variable('SALMON_POISSON_CONTRACT_TRACE', env_trace, status=env_status)
  if (env_status /= 0) return
  select case(trim(adjustl(env_trace)))
  case('1','y','Y','yes','YES','true','TRUE','on','ON')
    poisson_contract_trace_enabled = .true.
  end select
end function poisson_contract_trace_enabled

subroutine print_poisson_solver_contract(tag, lg, mg, info, scalar)
  use structures, only: s_rgrid, s_parallel_info, s_scalar
  use communication, only: comm_summation, comm_get_min, comm_get_max, comm_is_root
  use parallelization, only: nproc_id_global
  implicit none
  character(*),          intent(in) :: tag
  type(s_rgrid),         intent(in) :: lg, mg
  type(s_parallel_info), intent(in) :: info
  type(s_scalar),        intent(in) :: scalar
  real(8) :: lsum, gsum, lsum2, gsum2
  real(8) :: lmin, gmin, lmax, gmax
  real(8) :: min_in(1), min_out(1), max_in(1), max_out(1)
  integer :: npts_global

  lsum = sum(scalar%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  lsum2 = sum(scalar%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))**2)
  lmin = minval(scalar%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  lmax = maxval(scalar%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))

  call comm_summation(lsum, gsum, info%icomm_r)
  call comm_summation(lsum2, gsum2, info%icomm_r)
  min_in(1) = lmin
  max_in(1) = lmax
  call comm_get_min(min_in, min_out, 1, info%icomm_r)
  call comm_get_max(max_in, max_out, 1, info%icomm_r)
  gmin = min_out(1)
  gmax = max_out(1)

  if (.not. comm_is_root(nproc_id_global)) return

  npts_global = lg%num(1) * lg%num(2) * lg%num(3)
  write(*,'(1x,a,1x,a,1x,a,3(i0,1x),a,3(i0,1x),a,es23.15,a,es23.15,a,es23.15,a,es23.15,a,es23.15)') &
    '[POISSON-CONTRACT]', trim(tag), 'lg_num=', lg%num(1), lg%num(2), lg%num(3), &
    ' mg_is=', mg%is(1), mg%is(2), mg%is(3), &
    ' sum=', gsum, ' mean=', gsum / dble(npts_global), ' l2=', sqrt(gsum2), ' min=', gmin, ' max=', gmax
  flush(6)
end subroutine print_poisson_solver_contract

subroutine print_poisson_kernel_contract(tag, mg, info, coef, zrhoG, kernel_in, kernel_out)
  use structures, only: s_rgrid, s_parallel_info
  use communication, only: comm_summation, comm_is_root
  use parallelization, only: nproc_id_global
  implicit none
  character(*),          intent(in) :: tag
  type(s_rgrid),         intent(in) :: mg
  type(s_parallel_info), intent(in) :: info
  real(8),               intent(in) :: coef(:,:,:)
  complex(8),            intent(in) :: zrhoG(:,:,:)
  complex(8),            intent(in) :: kernel_in(:,:,:)
  complex(8),            intent(in) :: kernel_out(:,:,:)
  real(8) :: lnorm_rhog, gnorm_rhog
  real(8) :: lnorm_vhg, gnorm_vhg
  complex(8) :: g0_rhog, g0_kernel_in, g0_kernel_out
  real(8) :: g0_coef

  lnorm_rhog = sum(abs(zrhoG(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))**2)
  lnorm_vhg = sum(abs(kernel_out(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))**2)
  call comm_summation(lnorm_rhog, gnorm_rhog, info%icomm_r)
  call comm_summation(lnorm_vhg, gnorm_vhg, info%icomm_r)

  g0_rhog = (0.d0, 0.d0)
  g0_kernel_in = (0.d0, 0.d0)
  g0_kernel_out = (0.d0, 0.d0)
  g0_coef = 0.d0
  if (1 >= lbound(zrhoG,1) .and. 1 <= ubound(zrhoG,1) .and. &
      1 >= lbound(zrhoG,2) .and. 1 <= ubound(zrhoG,2) .and. &
      1 >= lbound(zrhoG,3) .and. 1 <= ubound(zrhoG,3)) then
    g0_rhog = zrhoG(1,1,1)
  end if
  if (1 >= lbound(kernel_in,1) .and. 1 <= ubound(kernel_in,1) .and. &
      1 >= lbound(kernel_in,2) .and. 1 <= ubound(kernel_in,2) .and. &
      1 >= lbound(kernel_in,3) .and. 1 <= ubound(kernel_in,3)) then
    g0_kernel_in = kernel_in(1,1,1)
  end if
  if (1 >= lbound(kernel_out,1) .and. 1 <= ubound(kernel_out,1) .and. &
      1 >= lbound(kernel_out,2) .and. 1 <= ubound(kernel_out,2) .and. &
      1 >= lbound(kernel_out,3) .and. 1 <= ubound(kernel_out,3)) then
    g0_kernel_out = kernel_out(1,1,1)
  end if
  if (1 >= lbound(coef,1) .and. 1 <= ubound(coef,1) .and. &
      1 >= lbound(coef,2) .and. 1 <= ubound(coef,2) .and. &
      1 >= lbound(coef,3) .and. 1 <= ubound(coef,3)) then
    g0_coef = coef(1,1,1)
  end if

  if (.not. comm_is_root(nproc_id_global)) return
  write(*,'(1x,a,1x,a,1x,a,es23.15,a,2(es23.15,1x),a,2(es23.15,1x),a,es23.15,a,es23.15,a,es23.15)') &
    '[POISSON-CONTRACT]', trim(tag), 'g0_coef=', g0_coef, &
    ' g0_rhoG(re,im)=', real(g0_rhog), aimag(g0_rhog), &
    ' g0_kernel_out(re,im)=', real(g0_kernel_out), aimag(g0_kernel_out), &
    ' ||rhoG||_2=', sqrt(gnorm_rhog), ' ||kernel_out||_2=', sqrt(gnorm_vhg), ' g0_kernel_in_re=', real(g0_kernel_in)
  flush(6)

  call maybe_dump_poisson_spectral(tag, mg, info, zrhoG, kernel_out)
end subroutine print_poisson_kernel_contract

logical function poisson_spectral_dump_enabled()
  implicit none
  character(16) :: env_dump
  integer :: env_status
  poisson_spectral_dump_enabled = .false.
  env_dump = ''
  call get_environment_variable('SALMON_POISSON_SPECTRAL_DUMP', env_dump, status=env_status)
  if (env_status /= 0) return
  select case(trim(adjustl(env_dump)))
  case('1','y','Y','yes','YES','true','TRUE','on','ON')
    poisson_spectral_dump_enabled = .true.
  end select
end function poisson_spectral_dump_enabled

subroutine maybe_dump_poisson_spectral(tag, mg, info, zrhoG, kernel_out)
  use structures, only: s_rgrid, s_parallel_info
  use communication, only: comm_summation, comm_is_root
  use parallelization, only: nproc_id_global
  implicit none
  character(*),          intent(in) :: tag
  type(s_rgrid),         intent(in) :: mg
  type(s_parallel_info), intent(in) :: info
  complex(8),            intent(in) :: zrhoG(:,:,:)
  complex(8),            intent(in) :: kernel_out(:,:,:)
  complex(8), allocatable :: local_rhog(:,:,:), full_rhog(:,:,:)
  complex(8), allocatable :: local_vhg(:,:,:), full_vhg(:,:,:)
  integer :: ixs, ixe, iys, iye, izs, ize
  integer :: ix0, ix1, iy0, iy1, iz0, iz1
  integer :: ix, iy, iz, nall
  integer :: unit
  logical :: do_dump
  character(8) :: solver
  character(16) :: ctx
  character(256) :: file_rhog, file_vhg

  if (.not. poisson_spectral_dump_enabled()) return

  ctx = adjustl(trim(poisson_contract_context))
  solver = 'UNKNOWN'
  if (index(tag, 'POISSON-FT-KERNEL') > 0) solver = 'FT'
  if (index(tag, 'POISSON-FFTE-KERNEL') > 0) solver = 'FFTE'

  do_dump = .false.
  select case (trim(ctx))
  case ('DC')
    if (solver == 'FT' .and. .not. dumped_dc_ft) then
      dumped_dc_ft = .true.
      do_dump = .true.
    else if (solver == 'FFTE' .and. .not. dumped_dc_ffte) then
      dumped_dc_ffte = .true.
      do_dump = .true.
    end if
  case ('DG')
    if (solver == 'FT' .and. .not. dumped_dg_ft) then
      dumped_dg_ft = .true.
      do_dump = .true.
    else if (solver == 'FFTE' .and. .not. dumped_dg_ffte) then
      dumped_dg_ffte = .true.
      do_dump = .true.
    end if
  case default
    if (solver == 'FT' .and. .not. dumped_dc_ft) then
      dumped_dc_ft = .true.
      do_dump = .true.
    else if (solver == 'FFTE' .and. .not. dumped_dc_ffte) then
      dumped_dc_ffte = .true.
      do_dump = .true.
    end if
  end select
  if (.not. do_dump) return

  ixs = lbound(zrhoG,1)
  ixe = ubound(zrhoG,1)
  iys = lbound(zrhoG,2)
  iye = ubound(zrhoG,2)
  izs = lbound(zrhoG,3)
  ize = ubound(zrhoG,3)

  allocate(local_rhog(ixs:ixe, iys:iye, izs:ize), full_rhog(ixs:ixe, iys:iye, izs:ize))
  allocate(local_vhg(ixs:ixe, iys:iye, izs:ize), full_vhg(ixs:ixe, iys:iye, izs:ize))
  local_rhog = (0.d0, 0.d0)
  local_vhg = (0.d0, 0.d0)

  ix0 = max(mg%is(1), ixs)
  ix1 = min(mg%ie(1), ixe)
  iy0 = max(mg%is(2), iys)
  iy1 = min(mg%ie(2), iye)
  iz0 = max(mg%is(3), izs)
  iz1 = min(mg%ie(3), ize)
  if (ix0 <= ix1 .and. iy0 <= iy1 .and. iz0 <= iz1) then
    local_rhog(ix0:ix1, iy0:iy1, iz0:iz1) = zrhoG(ix0:ix1, iy0:iy1, iz0:iz1)
    local_vhg(ix0:ix1, iy0:iy1, iz0:iz1) = kernel_out(ix0:ix1, iy0:iy1, iz0:iz1)
  end if

  nall = size(local_rhog)
  call comm_summation(local_rhog, full_rhog, nall, info%icomm_r)
  call comm_summation(local_vhg, full_vhg, nall, info%icomm_r)

  if (comm_is_root(nproc_id_global)) then
    write(file_rhog,'(a,a,a,a)') 'poisson_spectral_', trim(ctx), '_', trim(solver)//'_rhoG.dat'
    write(file_vhg,'(a,a,a,a)') 'poisson_spectral_', trim(ctx), '_', trim(solver)//'_vhG.dat'

    open(newunit=unit, file=trim(file_rhog), status='replace', action='write')
    write(unit,'(a)') '# ix iy iz re im'
    do iz = izs, ize
      do iy = iys, iye
        do ix = ixs, ixe
          write(unit,'(3(i0,1x),2(es23.15,1x))') ix, iy, iz, real(full_rhog(ix,iy,iz)), aimag(full_rhog(ix,iy,iz))
        end do
      end do
    end do
    close(unit)

    open(newunit=unit, file=trim(file_vhg), status='replace', action='write')
    write(unit,'(a)') '# ix iy iz re im'
    do iz = izs, ize
      do iy = iys, iye
        do ix = ixs, ixe
          write(unit,'(3(i0,1x),2(es23.15,1x))') ix, iy, iz, real(full_vhg(ix,iy,iz)), aimag(full_vhg(ix,iy,iz))
        end do
      end do
    end do
    close(unit)
  end if

  deallocate(local_rhog, full_rhog, local_vhg, full_vhg)
end subroutine maybe_dump_poisson_spectral


end module poisson_periodic
