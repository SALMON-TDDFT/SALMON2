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
  implicit none

contains

  ! hartree_dg_distributed
  !   Compute Vh%f = Hartree potential from rho%f using a distributed
  !   dimension-by-dimension DFT (same algorithm as poisson_ft, but
  !   kz-work split across dg_frag%icomm).
  !
  !   On exit:
  !     Vh%f              — full Hartree potential, identical on all ranks
  !     poisson%zrhoG_ele — full rho(G), identical on all ranks
  subroutine hartree_dg_distributed(lg, mg, fg, poisson, dg_frag, rho, Vh)
    use structures,           only: s_rgrid, s_reciprocal_grid, s_poisson, s_scalar
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use communication,        only: comm_summation
    use inputoutput,          only: yn_put_wall_z_boundary, wall_height, wall_width
    use math_constants,       only: pi
    use salmon_global,        only: hse_omega
    implicit none

    type(s_rgrid),           intent(in)    :: lg, mg
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_poisson),         intent(inout) :: poisson
    type(s_dg_fragment_rt),  intent(in)    :: dg_frag
    type(s_scalar),          intent(in)    :: rho
    type(s_scalar),          intent(inout) :: Vh

    integer :: Nx, Ny, Nz, N_total
    integer :: ix, iy, iz, kx, ky, kz, kz_loc, nkz_local, nkz_actual, kz_start, kz_end
    real(8) :: inv_N, g2, sr_factor, Vwall_z, z, z0
    logical :: use_hse_sr_hartree
    character(16) :: env_hse_sr
    integer :: env_status

    ! Local work arrays for the kz-slab owned by this rank
    complex(8), allocatable :: ff1(:,:,:)      ! (mg%is(1):ie(1), mg%is(2):ie(2), nkz_local)
    complex(8), allocatable :: ff2(:,:,:)      ! same shape as ff1

    ! Partial contributions on the full real-space / G-space grids
    real(8),    allocatable :: Vh_partial(:,:,:)      ! (mg range)
    complex(8), allocatable :: rhoG_partial(:,:,:)    ! (mg range)

    ! Check SALMON_HSE_SR_HARTREE environment variable (same logic as hartree_sub::hartree)
    env_hse_sr = ''
    use_hse_sr_hartree = .false.
    call get_environment_variable('SALMON_HSE_SR_HARTREE', env_hse_sr, status=env_status)
    if (env_status == 0) then
      select case(trim(adjustl(env_hse_sr)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        use_hse_sr_hartree = .true.
      end select
    end if
    if (use_hse_sr_hartree .and. hse_omega <= 0.0d0) then
      stop 'poisson_dg_distributed: hse_omega must be > 0 when SALMON_HSE_SR_HARTREE is enabled'
    end if

    ! N_total based on mg (the arrays are allocated over mg range).
    ! In DG-fragment lg==mg so this equals lg%num(1)*lg%num(2)*lg%num(3),
    ! but using mg%num is correct for future callers where lg /= mg.
    Nx = mg%num(1)
    Ny = mg%num(2)
    Nz = mg%num(3)
    N_total = Nx * Ny * Nz
    inv_N   = 1.0d0 / dble(lg%num(1) * lg%num(2) * lg%num(3))

    ! -----------------------------------------------------------------------
    ! Distribute kz indices across dg_frag%icomm.
    ! rank r (0-based: dg_frag%id) owns kz in [kz_start : kz_end].
    ! nkz_local  = ceiling(Nz / isize)  — allocation size (may be > actual)
    ! nkz_actual = number of valid kz for this rank (kz_end - kz_start + 1)
    ! Nz here is mg%num(3); in DG-fragment this equals lg%num(3).
    ! -----------------------------------------------------------------------
    nkz_local  = (Nz + dg_frag%isize - 1) / dg_frag%isize
    kz_start   = dg_frag%id * nkz_local + mg%is(3)
    kz_end     = min(kz_start + nkz_local - 1, mg%ie(3))
    nkz_actual = max(0, kz_end - kz_start + 1)

    ! Allocate local slab work arrays and global partial arrays
    allocate(ff1(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), nkz_local))
    allocate(ff2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), nkz_local))
    allocate(Vh_partial   (mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(rhoG_partial (mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))

    ff1          = (0.0d0, 0.0d0)
    ff2          = (0.0d0, 0.0d0)
    Vh_partial   = 0.0d0
    rhoG_partial = (0.0d0, 0.0d0)

    ! =======================================================================
    ! FORWARD 3D DFT: rho(r) -> rho_G for local kz slab
    ! No communication is needed: every rank already holds the full rho%f.
    ! =======================================================================

    ! -- Step 1: z-DFT --
    ! ff1(ix,iy,kz_loc) = sum_iz  egzc(kz, iz) * rho(ix, iy, iz)
    !$omp parallel do private(kz_loc, kz, iy, ix)
    do kz_loc = 1, nkz_actual
      kz = kz_start + kz_loc - 1
      do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          ff1(ix, iy, kz_loc) = sum(fg%egzc(kz, :) * dcmplx(rho%f(ix, iy, :)))
        end do
      end do
    end do
    !$omp end parallel do

    ! -- Step 2: y-DFT --
    ! ff2(ix,ky,kz_loc) = sum_iy  egyc(ky, iy) * ff1(ix, iy, kz_loc)
    !$omp parallel do private(kz_loc, ky, ix)
    do kz_loc = 1, nkz_actual
      do ky = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          ff2(ix, ky, kz_loc) = sum(fg%egyc(ky, :) * ff1(ix, :, kz_loc))
        end do
      end do
    end do
    !$omp end parallel do

    ! -- Step 3: x-DFT + normalize --
    ! ff1(kx,ky,kz_loc) = sum_ix  egxc(kx, ix) * ff2(ix, ky, kz_loc) / N
    !$omp parallel do private(kz_loc, ky, kx)
    do kz_loc = 1, nkz_actual
      do ky = mg%is(2), mg%ie(2)
        do kx = mg%is(1), mg%ie(1)
          ff1(kx, ky, kz_loc) = sum(fg%egxc(kx, :) * ff2(:, ky, kz_loc)) * inv_N
        end do
      end do
    end do
    !$omp end parallel do

    ! -- Save rho(G) for local kz slab (energy / force calculations) --
    !$omp parallel do private(kz_loc, kz, ky, kx)
    do kz_loc = 1, nkz_actual
      kz = kz_start + kz_loc - 1
      do ky = mg%is(2), mg%ie(2)
        do kx = mg%is(1), mg%ie(1)
          rhoG_partial(kx, ky, kz) = ff1(kx, ky, kz_loc)
        end do
      end do
    end do
    !$omp end parallel do

    ! =======================================================================
    ! COULOMB KERNEL: Vh(G) = coef(G) * rho(G)
    ! When SALMON_HSE_SR_HARTREE is set, apply the HSE short-range factor
    ! (1 - exp(-|G|^2 / (4*omega^2))) to match hartree_sub::hartree behaviour.
    ! =======================================================================
    if (use_hse_sr_hartree) then
      !$omp parallel do private(kz_loc, kz, ky, kx, g2, sr_factor)
      do kz_loc = 1, nkz_actual
        kz = kz_start + kz_loc - 1
        do ky = mg%is(2), mg%ie(2)
          do kx = mg%is(1), mg%ie(1)
            g2 = fg%vec_G(1,kx,ky,kz)**2 + fg%vec_G(2,kx,ky,kz)**2 + fg%vec_G(3,kx,ky,kz)**2
            sr_factor = 1.0d0 - exp(-g2 / (4.0d0 * hse_omega * hse_omega))
            ff1(kx, ky, kz_loc) = fg%coef(kx, ky, kz) * sr_factor * ff1(kx, ky, kz_loc)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(kz_loc, kz, ky, kx)
      do kz_loc = 1, nkz_actual
        kz = kz_start + kz_loc - 1
        do ky = mg%is(2), mg%ie(2)
          do kx = mg%is(1), mg%ie(1)
            ff1(kx, ky, kz_loc) = fg%coef(kx, ky, kz) * ff1(kx, ky, kz_loc)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    ! =======================================================================
    ! INVERSE 3D DFT: Vh(G) -> partial Vh(r) contribution from local kz slab
    ! =======================================================================

    ! -- Step 4: inverse x-DFT --
    ! ff2(ix,ky,kz_loc) = sum_kx  egx(kx, ix) * ff1(kx, ky, kz_loc)
    !$omp parallel do private(kz_loc, ky, ix)
    do kz_loc = 1, nkz_actual
      do ky = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          ff2(ix, ky, kz_loc) = sum(fg%egx(:, ix) * ff1(:, ky, kz_loc))
        end do
      end do
    end do
    !$omp end parallel do

    ! -- Step 5: inverse y-DFT --
    ! ff1(ix,iy,kz_loc) = sum_ky  egy(ky, iy) * ff2(ix, ky, kz_loc)
    !$omp parallel do private(kz_loc, iy, ix)
    do kz_loc = 1, nkz_actual
      do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          ff1(ix, iy, kz_loc) = sum(fg%egy(:, iy) * ff2(ix, :, kz_loc))
        end do
      end do
    end do
    !$omp end parallel do

    ! -- Step 6: partial inverse z contribution --
    ! Vh_partial(ix,iy,iz) += sum_{kz in local}  Re[egz(kz,iz) * ff1(ix,iy,kz_loc)]
    ! Vh%f is real; the imaginary parts cancel exactly when all kz slabs are summed.
    !$omp parallel do private(iz, iy, ix, kz_loc, kz) collapse(3)
    do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          do kz_loc = 1, nkz_actual
            kz = kz_start + kz_loc - 1
            Vh_partial(ix, iy, iz) = Vh_partial(ix, iy, iz) + &
                real(fg%egz(kz, iz) * ff1(ix, iy, kz_loc))
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! =======================================================================
    ! ALLREDUCE: sum partial Vh and rho(G) contributions across all ranks.
    ! Each rank owns a disjoint kz slab, so MPI_SUM = MPI_Gatherall.
    ! N_total = mg%num product, matching the allocation of both arrays.
    ! =======================================================================
    call comm_summation(Vh_partial,   Vh%f,               N_total, dg_frag%icomm)
    call comm_summation(rhoG_partial, poisson%zrhoG_ele,  N_total, dg_frag%icomm)

    deallocate(ff1, ff2, Vh_partial, rhoG_partial)

    ! =======================================================================
    ! WALL POTENTIAL at z boundaries (yn_put_wall_z_boundary='y').
    ! Applied after Vh is assembled, matching hartree_sub::hartree behaviour.
    ! dg_frag%hgs(3) is the grid spacing in z.
    ! =======================================================================
    if (yn_put_wall_z_boundary == 'y') then
      z0 = lg%num(3) * dg_frag%hgs(3)
      !$omp parallel do private(iz, iy, ix, z, Vwall_z)
      do iz = mg%is(3), mg%ie(3)
        z = iz * dg_frag%hgs(3)
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

  end subroutine hartree_dg_distributed

end module poisson_dg_distributed
