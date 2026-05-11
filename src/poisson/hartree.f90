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
!=======================================================================
!============================ Hartree potential (Solve Poisson equation)
module hartree_sub
  implicit none

contains

!===================================================================================================================================
subroutine hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  use math_constants,only: pi
  use phys_constants,only: au_aa, au_ev
  use inputoutput, only: iperiodic,yn_ffte,yn_put_wall_z_boundary, &
                         method_poisson
#ifdef USE_FFTW
  use inputoutput, only: yn_fftw
#endif
  use structures, only: s_rgrid,s_dft_system,s_parallel_info,s_poisson,  &
                        s_sendrecv_grid,s_stencil,s_scalar,s_reciprocal_grid,  &
                        allocate_scalar, deallocate_scalar
  use communication, only: comm_is_root, comm_summation, comm_get_min, comm_get_max
  use parallelization, only: nproc_id_global
  use poisson_isolated
  use poisson_periodic
  use poisson_dirichlet, only: jones
  use salmon_global, only: hse_omega
  use nvtx
  implicit none
  type(s_rgrid)          ,intent(in)    :: lg
  type(s_rgrid)          ,intent(in)    :: mg
  type(s_parallel_info)  ,intent(in)    :: info
  type(s_dft_system)     ,intent(in)    :: system
  type(s_reciprocal_grid),intent(in)    :: fg
  type(s_poisson)        ,intent(inout) :: poisson
  type(s_sendrecv_grid)  ,intent(inout) :: srg_scalar
  type(s_stencil)        ,intent(in)    :: stencil
  type(s_scalar)         ,intent(in)    :: rho
  type(s_scalar)         ,intent(inout) :: Vh
  character(16) :: env_hse_sr
  character(16) :: env_contract_trace
  logical :: use_hse_sr_hartree
  logical :: enable_contract_trace
  integer :: env_status

  env_hse_sr = ''
  use_hse_sr_hartree = .false.
  enable_contract_trace = .false.
  call get_environment_variable('SALMON_HSE_SR_HARTREE', env_hse_sr, status=env_status)
  if (env_status == 0) then
    select case(trim(adjustl(env_hse_sr)))
    case('1','y','Y','yes','YES','true','TRUE','on','ON')
      use_hse_sr_hartree = .true.
    end select
  end if
  env_contract_trace = ''
  call get_environment_variable('SALMON_POISSON_CONTRACT_TRACE', env_contract_trace, status=env_status)
  if (env_status == 0) then
    select case(trim(adjustl(env_contract_trace)))
    case('1','y','Y','yes','YES','true','TRUE','on','ON')
      enable_contract_trace = .true.
    end select
  end if

  if (enable_contract_trace) then
    call print_scalar_contract_stats('DC-HARTREE-rho-in', lg, mg, info, rho)
  end if

  call set_poisson_contract_context('DC')

  call nvtxStartRange('hartree', __LINE__)
  
  select case(iperiodic)
  case(0)
    select case(method_poisson)
    case('cg')
      call poisson_isolated_cg(lg,mg,info,system,poisson,rho%f,Vh%f,srg_scalar,stencil)
    case('ft')
#ifdef USE_FFTW
      select case(yn_fftw)
      case('n')
#endif
        select case(yn_ffte)
        case('n')
          call poisson_isolated_ft(lg,mg,info,fg,rho,Vh,poisson)
        case('y')
          call poisson_isolated_ffte(lg,mg,info,fg,rho,Vh,poisson)
        end select
#ifdef USE_FFTW
      case('y')
        call poisson_isolated_fftw(lg,mg,info,fg,rho,Vh,poisson)
      end select
#endif
    case('dirichlet')
      call jones(lg,mg,info,system,rho,Vh,poisson)
    end select
  case(3)
#ifdef USE_FFTW
    select case(yn_fftw)
    case('n')
#endif
      select case(yn_ffte)
      case('n')
        if (use_hse_sr_hartree) then
          call poisson_ft_hse_sr(lg,mg,info,fg,rho,Vh,poisson,hse_omega)
        else
          call poisson_ft(lg,mg,info,fg,rho,Vh,poisson)
        end if
      case('y')
        if (use_hse_sr_hartree) then
          call poisson_ffte_hse_sr(lg,mg,info,fg,rho,Vh,poisson,hse_omega)
        else
          call poisson_ffte(lg,mg,info,fg,rho,Vh,poisson)
        end if
      end select
#ifdef USE_FFTW
    case('y')
      if (use_hse_sr_hartree) then
        call poisson_fftw_hse_sr(lg,mg,info,fg,rho,Vh,poisson,hse_omega)
      else
        call poisson_fftw(lg,mg,info,fg,rho,Vh,poisson)
      end if
    end select
#endif
  end select

  if (enable_contract_trace) then
    call print_scalar_contract_stats('DC-HARTREE-vh-out', lg, mg, info, Vh)
  end if

  !potentiall wall at the boundary on z direction
  if(yn_put_wall_z_boundary=='y') call add_potential_wall

  call nvtxEndRange

  contains

    subroutine  add_potential_wall
      use inputoutput, only: wall_height, wall_width
      implicit none
      integer :: ix,iy,iz
      real(8) :: Vwall_z, z,z0

      !$omp parallel do private(iz,iy,ix,z,z0,Vwall_z)
      do iz = mg%is(3),mg%ie(3)
         z  = iz*system%hgs(3)
         z0 = lg%num(3) * system%hgs(3)
         if( z .le. wall_width ) then
            Vwall_z = wall_height * cos((z/wall_width)*pi/2d0)**2
         else if( z .ge. z0-wall_width ) then
            Vwall_z = wall_height * cos(((z0-z)/wall_width)*pi/2d0)**2
         else
            cycle
         endif

         do iy = mg%is(2),mg%ie(2)
         do ix = mg%is(1),mg%ie(1)
            Vh%f(ix,iy,iz) = Vh%f(ix,iy,iz) + Vwall_z
         end do
         end do
      end do
      !$omp end parallel do

    end subroutine add_potential_wall

    subroutine print_scalar_contract_stats(tag, lg_local, mg_local, info_local, scalar)
      implicit none
      character(*),          intent(in) :: tag
      type(s_rgrid),         intent(in) :: lg_local, mg_local
      type(s_parallel_info), intent(in) :: info_local
      type(s_scalar),        intent(in) :: scalar
      real(8) :: lsum, gsum, lsum2, gsum2
      real(8) :: lmin, gmin, lmax, gmax
      real(8) :: min_in(1), min_out(1), max_in(1), max_out(1)
      integer :: npts_global

      lsum = sum(scalar%f(mg_local%is(1):mg_local%ie(1), mg_local%is(2):mg_local%ie(2), mg_local%is(3):mg_local%ie(3)))
      lsum2 = sum(scalar%f(mg_local%is(1):mg_local%ie(1), mg_local%is(2):mg_local%ie(2), mg_local%is(3):mg_local%ie(3))**2)
      lmin = minval(scalar%f(mg_local%is(1):mg_local%ie(1), mg_local%is(2):mg_local%ie(2), mg_local%is(3):mg_local%ie(3)))
      lmax = maxval(scalar%f(mg_local%is(1):mg_local%ie(1), mg_local%is(2):mg_local%ie(2), mg_local%is(3):mg_local%ie(3)))

      call comm_summation(lsum, gsum, info_local%icomm_r)
      call comm_summation(lsum2, gsum2, info_local%icomm_r)
      min_in(1) = lmin
      max_in(1) = lmax
      call comm_get_min(min_in, min_out, 1, info_local%icomm_r)
      call comm_get_max(max_in, max_out, 1, info_local%icomm_r)
      gmin = min_out(1)
      gmax = max_out(1)

      if (.not. comm_is_root(nproc_id_global)) return

      npts_global = lg_local%num(1) * lg_local%num(2) * lg_local%num(3)
      write(*,'(1x,a,1x,a,1x,a,3(i0,1x),a,3(i0,1x),a,es23.15,a,es23.15,a,es23.15,a,es23.15,a,es23.15)') &
        '[POISSON-CONTRACT]', trim(tag), 'lg_num=', lg_local%num(1), lg_local%num(2), lg_local%num(3), &
        ' mg_is=', mg_local%is(1), mg_local%is(2), mg_local%is(3), &
        ' sum=', gsum, ' mean=', gsum / dble(npts_global), ' l2=', sqrt(gsum2), ' min=', gmin, ' max=', gmax
      flush(6)
    end subroutine print_scalar_contract_stats

end subroutine hartree

end module hartree_sub
