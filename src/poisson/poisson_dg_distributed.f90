!
!  Copyright 2019-2026 SALMON developers
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
! Distributed Hartree potential via slab-FFT across the DG-fragment
! MPI communicator (dg_frag%icomm).
!
! Context:
!   In DG-fragment RT-TDDFT, lg == mg (no spatial grid decomposition).
!   All ranks independently hold the full electron density rho%f.
!   Without distribution, every rank runs the identical O(N*Nz) DFT
!   redundantly.  Here we distribute kz-slabs across dg_frag%isize
!   ranks so each rank does only O(N*Nz/isize) work, then one
!   MPI_Allreduce assembles the full Vh%f on all ranks.
!
! Prerequisites:
!   - rho%f must be identical on every rank on entry (guaranteed by
!     calculate_density_from_fragments + comm_summation in reconstruct).
!   - poisson%ff1*/ff2* work arrays are NOT used here; local temporaries
!     are allocated and freed instead to keep the standard poisson state
!     intact for other callers.
!   - fg%egzc(kz,iz), fg%egyc(ky,iy), fg%egxc(kx,ix): forward DFT phases
!     allocated as (lg%is:lg%ie, lg%is:lg%ie) in initialization.f90.
!   - fg%egz(kz,iz), fg%egy(ky,iy), fg%egx(kx,ix): inverse DFT phases.
!   - fg%coef(kx,ky,kz): Coulomb kernel 4*pi/|G|^2, zero at G=0.
!=======================================================================

module poisson_dg_distributed
  use structures, only: s_reciprocal_grid, s_poisson, s_rgrid, s_parallel_info
  implicit none

contains

  subroutine hartree_dg_via_parent_localrho(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    use structures,           only: s_rgrid, s_reciprocal_grid, s_poisson, s_scalar, s_parallel_info
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use poisson_periodic,     only: poisson_ft, poisson_ffte, poisson_ft_hse_sr, poisson_ffte_hse_sr, set_poisson_contract_context
    use inputoutput,          only: yn_ffte, yn_put_wall_z_boundary, wall_height, wall_width
    use salmon_global,        only: hse_omega
    use math_constants,       only: pi
    use communication,        only: comm_summation, comm_get_min, comm_get_max, comm_is_root
    use parallelization,      only: nproc_id_global
    implicit none

    type(s_parallel_info),  intent(in)    :: info
    type(s_rgrid),           intent(in)    :: lg, mg
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_poisson),         intent(inout) :: poisson
    type(s_dg_fragment_rt),  intent(in)    :: dg_frag
    type(s_scalar),          intent(in)    :: rho
    type(s_scalar),          intent(inout) :: Vh

    character(16) :: env_hse_sr
    character(16) :: env_contract_trace
    logical :: use_hse_sr_hartree
    logical :: enable_contract_trace
    integer :: env_status
    integer :: ix, iy, iz
    real(8) :: z, z0, Vwall_z

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
    if (use_hse_sr_hartree .and. hse_omega <= 0.0d0) then
      stop 'poisson_dg_distributed: hse_omega must be > 0 when SALMON_HSE_SR_HARTREE is enabled'
    end if

    if (enable_contract_trace) then
      call print_scalar_contract_stats('DG-HARTREE-rho-in', lg, mg, info, rho)
    end if

    if (enable_contract_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,4(a,i0))') "        hartree trace: using-parent-localrho-main", &
        " id=", dg_frag%id, " isize=", dg_frag%isize, " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag
      flush(6)
    end if

    call set_poisson_contract_context('DG')

    if (yn_ffte == 'y') then
      if (use_hse_sr_hartree) then
        call poisson_ffte_hse_sr(lg, mg, info, fg, rho, Vh, poisson, hse_omega)
      else
        call poisson_ffte(lg, mg, info, fg, rho, Vh, poisson)
      end if
    else
      if (use_hse_sr_hartree) then
        call poisson_ft_hse_sr(lg, mg, info, fg, rho, Vh, poisson, hse_omega)
      else
        call poisson_ft(lg, mg, info, fg, rho, Vh, poisson)
      end if
    end if

    if (enable_contract_trace) then
      call print_scalar_contract_stats('DG-HARTREE-vh-out-prewall', lg, mg, info, Vh)
    end if

    if (yn_put_wall_z_boundary == 'y') then
      !$omp parallel do private(iz,iy,ix,z,z0,Vwall_z) schedule(static)
      do iz = mg%is(3), mg%ie(3)
        z = iz * dg_frag%hgs(3)
        z0 = lg%num(3) * dg_frag%hgs(3)
        if (z <= wall_width) then
          Vwall_z = wall_height * cos((z / wall_width) * pi / 2.0d0)**2
        else if (z >= z0 - wall_width) then
          Vwall_z = wall_height * cos(((z0 - z) / wall_width) * pi / 2.0d0)**2
        else
          cycle
        end if
        do iy = mg%is(2), mg%ie(2)
          do ix = mg%is(1), mg%ie(1)
            Vh%f(ix, iy, iz) = Vh%f(ix, iy, iz) + Vwall_z
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (enable_contract_trace) then
      call print_scalar_contract_stats('DG-HARTREE-vh-out-postwall', lg, mg, info, Vh)
    end if

  contains

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
  end subroutine hartree_dg_via_parent_localrho

  ! hartree_dg_distributed
  !   Compute Vh%f = Hartree potential from rho%f using a distributed
  !   dimension-by-dimension DFT (same algorithm as poisson_ft, but
  !   kz-work split across dg_frag%icomm).
  subroutine hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    use structures,           only: s_rgrid, s_reciprocal_grid, s_poisson, s_scalar, s_parallel_info
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    implicit none

    type(s_parallel_info),  intent(in)    :: info
    type(s_rgrid),           intent(in)    :: lg, mg
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_poisson),         intent(inout) :: poisson
    type(s_dg_fragment_rt),  intent(in)    :: dg_frag
    type(s_scalar),          intent(in)    :: rho
    type(s_scalar),          intent(inout) :: Vh

    if (lbound(rho%f,1) == mg%is(1) .and. ubound(rho%f,1) == mg%ie(1) .and. &
        lbound(rho%f,2) == mg%is(2) .and. ubound(rho%f,2) == mg%ie(2) .and. &
        lbound(rho%f,3) == mg%is(3) .and. ubound(rho%f,3) == mg%ie(3)) then
      call hartree_dg_via_parent_localrho(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
      return
    end if

    stop 'poisson_dg_distributed: DG Hartree requires rho bounds to match mg bounds'
  end subroutine hartree_dg_distributed

end module poisson_dg_distributed
