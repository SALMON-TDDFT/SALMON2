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

  subroutine append_hartree_rank_trace(rank, label, ifrag_group, id_frag, isize_frag, comm_rank, comm_size)
    implicit none
    integer, intent(in) :: rank, ifrag_group, id_frag, isize_frag, comm_rank, comm_size
    character(*), intent(in) :: label
    integer :: iunit
    character(256) :: filename

    write(filename,'(a,i0,a)') 'hartree_rank_trace.', rank, '.log'
    open(newunit=iunit, file=trim(filename), status='unknown', position='append', action='write')
    write(iunit,'(a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') trim(label), &
      'ifrag_group=', ifrag_group, 'id_frag=', id_frag, 'isize_frag=', isize_frag, &
      'comm_rank=', comm_rank, 'comm_size=', comm_size
    close(iunit)
  end subroutine append_hartree_rank_trace

  subroutine allocate_fragment_poisson_periodic(lg, mg, info, use_ffte, poisson_frag)
    use structures, only: s_rgrid, s_parallel_info, s_poisson
    implicit none
    type(s_rgrid), intent(in) :: lg, mg
    type(s_parallel_info), intent(in) :: info
    logical, intent(in) :: use_ffte
    type(s_poisson), intent(inout) :: poisson_frag

    if (use_ffte) then
      allocate(poisson_frag%a_ffte(lg%num(1), mg%num(2), mg%num(3)))
      allocate(poisson_frag%b_ffte(lg%num(1), mg%num(2), mg%num(3)))
      call PZFFT3DV_MOD(poisson_frag%a_ffte, poisson_frag%b_ffte, lg%num(1), lg%num(2), lg%num(3), &
                        info%isize_y, info%isize_z, 0, info%icomm_y, info%icomm_z)
    else
      allocate(poisson_frag%ff1x(lg%is(1):lg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(poisson_frag%ff1y(mg%is(1):mg%ie(1), lg%is(2):lg%ie(2), mg%is(3):mg%ie(3)))
      allocate(poisson_frag%ff1z(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), lg%is(3):lg%ie(3)))
      allocate(poisson_frag%ff2x(lg%is(1):lg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(poisson_frag%ff2y(mg%is(1):mg%ie(1), lg%is(2):lg%ie(2), mg%is(3):mg%ie(3)))
      allocate(poisson_frag%ff2z(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), lg%is(3):lg%ie(3)))
    end if

    allocate(poisson_frag%zrhoG_ele(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  end subroutine allocate_fragment_poisson_periodic

  subroutine deallocate_fragment_poisson_periodic(poisson_frag)
    use structures, only: s_poisson
    implicit none
    type(s_poisson), intent(inout) :: poisson_frag

    if (allocated(poisson_frag%zrhoG_ele)) deallocate(poisson_frag%zrhoG_ele)
    if (allocated(poisson_frag%ff1x)) deallocate(poisson_frag%ff1x)
    if (allocated(poisson_frag%ff1y)) deallocate(poisson_frag%ff1y)
    if (allocated(poisson_frag%ff1z)) deallocate(poisson_frag%ff1z)
    if (allocated(poisson_frag%ff2x)) deallocate(poisson_frag%ff2x)
    if (allocated(poisson_frag%ff2y)) deallocate(poisson_frag%ff2y)
    if (allocated(poisson_frag%ff2z)) deallocate(poisson_frag%ff2z)
    if (allocated(poisson_frag%a_ffte)) deallocate(poisson_frag%a_ffte)
    if (allocated(poisson_frag%b_ffte)) deallocate(poisson_frag%b_ffte)
  end subroutine deallocate_fragment_poisson_periodic

  subroutine hartree_dg_via_poisson_periodic(lg, mg, fg, poisson, dg_frag, rho, Vh)
    use structures, only: s_rgrid, s_reciprocal_grid, s_poisson, s_scalar, s_parallel_info
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_parallel, only: setup_fragment_parallel_grid, finalize_fragment_parallel
    use communication, only: comm_summation, comm_sync_all, comm_get_groupinfo
    use poisson_periodic, only: poisson_ft, poisson_ffte, poisson_ft_hse_sr, poisson_ffte_hse_sr
    use inputoutput, only: yn_ffte, yn_put_wall_z_boundary, wall_height, wall_width
    use salmon_global, only: hse_omega
    use math_constants, only: pi
    implicit none

    type(s_rgrid), intent(in) :: lg, mg
    type(s_reciprocal_grid), intent(in) :: fg
    type(s_poisson), intent(inout) :: poisson
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_scalar), intent(in) :: rho
    type(s_scalar), intent(inout) :: Vh

    type(s_parallel_info) :: info_frag
    type(s_rgrid) :: mg_frag
    type(s_poisson) :: poisson_frag
    type(s_scalar) :: rho_work, rho_reduce, Vh_work, Vh_reduce
    complex(8), allocatable :: rhoG_work(:,:,:), rhoG_reduce(:,:,:)
    character(16) :: env_hse_sr
    logical :: use_hse_sr_hartree
    logical :: use_ffte_local
    integer :: env_status, ix, iy, iz
    integer :: frag_rank_check, frag_size_check
    real(8) :: z, z0, Vwall_z
    real(8) :: rho_probe_in, rho_probe_out

    env_hse_sr = ''
    use_hse_sr_hartree = .false.
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-entry"
      flush(6)
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        hartree trace: bridge-comm-state icomm_frag=", dg_frag%icomm_frag, &
        " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag, " ifrag_group=", dg_frag%ifrag_group
      flush(6)
    end if
    call comm_get_groupinfo(dg_frag%icomm_frag, frag_rank_check, frag_size_check)
    call append_hartree_rank_trace(dg_frag%id, 'hartree-entry', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "        hartree trace: rank-entry rank=", dg_frag%id, &
      " ifrag_group=", dg_frag%ifrag_group, " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag, &
      " comm_rank=", frag_rank_check, " comm_size=", frag_size_check
    flush(6)
    if (frag_rank_check /= dg_frag%id_frag .or. frag_size_check /= dg_frag%isize_frag) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] hartree communicator mismatch rank=", dg_frag%id, &
        " ifrag_group=", dg_frag%ifrag_group, " id_frag=", dg_frag%id_frag, " comm_rank=", frag_rank_check, &
        " comm_size=", frag_size_check
      flush(6)
      stop 'poisson_dg_distributed: dg_frag communicator metadata mismatch before barrier'
    end if
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

    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-before-setup-frag-grid"
      flush(6)
    end if
    call setup_fragment_parallel_grid(dg_frag, lg, info_frag, mg_frag)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,3(a,i0))') "        hartree trace: bridge-after-setup-frag-grid", &
        " mgx=", mg_frag%num(1), " mgy=", mg_frag%num(2), " mgz=", mg_frag%num(3)
      flush(6)
    end if
    use_ffte_local = (yn_ffte == 'y')
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-before-alloc-poisson"
      flush(6)
    end if
    call allocate_fragment_poisson_periodic(lg, mg_frag, info_frag, use_ffte_local, poisson_frag)
    call append_hartree_rank_trace(dg_frag%id, 'after-alloc-poisson', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-after-alloc-poisson"
      flush(6)
    end if

    allocate(rho_work%f(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3)))
    rho_work%f = 0.0d0
    allocate(rho_reduce%f(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3)))
    rho_reduce%f = 0.0d0
    call append_hartree_rank_trace(dg_frag%id, 'after-alloc-rho-work', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    if (lbound(rho%f,1) < lbound(rho_work%f,1) .or. ubound(rho%f,1) > ubound(rho_work%f,1) .or. &
        lbound(rho%f,2) < lbound(rho_work%f,2) .or. ubound(rho%f,2) > ubound(rho_work%f,2) .or. &
        lbound(rho%f,3) < lbound(rho_work%f,3) .or. ubound(rho%f,3) > ubound(rho_work%f,3)) then
      write(*,'(1x,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        [FATAL] rho local bounds outside work rank=", dg_frag%id_frag, &
        " work=", lbound(rho_work%f,1), ":", ubound(rho_work%f,1), ",", lbound(rho_work%f,2), ":", ubound(rho_work%f,2), ",", &
        lbound(rho_work%f,3), ":", ubound(rho_work%f,3), ",", " local=", lbound(rho%f,1), ":", ubound(rho%f,1), ",", &
        lbound(rho%f,2), ":", ubound(rho%f,2), ",", lbound(rho%f,3), ":", ubound(rho%f,3), ","
      flush(6)
      stop 'poisson_dg_distributed: rho local bounds outside work'
    end if
    write(*,'(1x,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        hartree trace: before-rho-allreduce rank=", dg_frag%id_frag, &
      " work=", lbound(rho_work%f,1), ":", ubound(rho_work%f,1), ",", lbound(rho_work%f,2), ":", ubound(rho_work%f,2), ",", &
      lbound(rho_work%f,3), ":", ubound(rho_work%f,3), ",", " local=", lbound(rho%f,1), ":", ubound(rho%f,1), ",", &
      lbound(rho%f,2), ":", ubound(rho%f,2), ",", lbound(rho%f,3), ":", ubound(rho%f,3), ","
    flush(6)
    call append_hartree_rank_trace(dg_frag%id, 'before-rho-copy', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    do iz = lbound(rho%f,3), ubound(rho%f,3)
      do iy = lbound(rho%f,2), ubound(rho%f,2)
        do ix = lbound(rho%f,1), ubound(rho%f,1)
          rho_work%f(ix, iy, iz) = rho%f(ix, iy, iz)
        end do
      end do
    end do
    call append_hartree_rank_trace(dg_frag%id, 'after-rho-copy', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    call append_hartree_rank_trace(dg_frag%id, 'before-rho-barrier', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    write(*,'(1x,a,i0,a,i0,a,i0)') "        hartree trace: before-rho-barrier rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group
    flush(6)
    call comm_sync_all(dg_frag%icomm_frag)
    call append_hartree_rank_trace(dg_frag%id, 'after-rho-barrier', dg_frag%ifrag_group, dg_frag%id_frag, dg_frag%isize_frag, &
      frag_rank_check, frag_size_check)
    write(*,'(1x,a,i0,a,i0,a,i0)') "        hartree trace: after-rho-barrier rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group
    flush(6)
    rho_probe_in = real(dg_frag%id_frag + 1, 8)
    rho_probe_out = -1.0d0
    write(*,'(1x,a,i0,a,1pe12.4)') "        hartree trace: before-rho-probe-allreduce rank=", dg_frag%id_frag, &
      " val=", rho_probe_in
    flush(6)
    call comm_summation(rho_probe_in, rho_probe_out, dg_frag%icomm_frag)
    write(*,'(1x,a,i0,a,1pe12.4)') "        hartree trace: after-rho-probe-allreduce rank=", dg_frag%id_frag, &
      " sum=", rho_probe_out
    flush(6)
    call comm_summation(rho_work%f, rho_reduce%f, size(rho_work%f), dg_frag%icomm_frag)
    write(*,'(1x,a,i0)') "        hartree trace: after-rho-allreduce rank=", dg_frag%id_frag
    flush(6)
    rho_work%f(:, :, :) = rho_reduce%f(:, :, :)

    allocate(Vh_work%f(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3)))
    Vh_work%f = 0.0d0
    allocate(Vh_reduce%f(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3)))
    Vh_reduce%f = 0.0d0
    allocate(rhoG_work(lbound(poisson%zrhoG_ele,1):ubound(poisson%zrhoG_ele,1), &
                       lbound(poisson%zrhoG_ele,2):ubound(poisson%zrhoG_ele,2), &
                       lbound(poisson%zrhoG_ele,3):ubound(poisson%zrhoG_ele,3)))
    rhoG_work = (0.0d0, 0.0d0)
    allocate(rhoG_reduce(lbound(rhoG_work,1):ubound(rhoG_work,1), lbound(rhoG_work,2):ubound(rhoG_work,2), &
                         lbound(rhoG_work,3):ubound(rhoG_work,3)))
    rhoG_reduce = (0.0d0, 0.0d0)

    if (dg_frag%id == 0) then
      if (use_ffte_local) then
        write(*,'(1x,a,a)') "        hartree trace: bridge=", "poisson_ffte"
      else
        write(*,'(1x,a,a)') "        hartree trace: bridge=", "poisson_ft"
      end if
      flush(6)
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-before-poisson-call"
      flush(6)
    end if
    if (use_ffte_local) then
      if (use_hse_sr_hartree) then
        call poisson_ffte_hse_sr(lg, mg_frag, info_frag, fg, rho_work, Vh_work, poisson_frag, hse_omega)
      else
        call poisson_ffte(lg, mg_frag, info_frag, fg, rho_work, Vh_work, poisson_frag)
      end if
    else
      if (use_hse_sr_hartree) then
        call poisson_ft_hse_sr(lg, mg_frag, info_frag, fg, rho_work, Vh_work, poisson_frag, hse_omega)
      else
        call poisson_ft(lg, mg_frag, info_frag, fg, rho_work, Vh_work, poisson_frag)
      end if
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-after-poisson-call"
      flush(6)
    end if

    rhoG_work(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3)) = &
      poisson_frag%zrhoG_ele(mg_frag%is(1):mg_frag%ie(1), mg_frag%is(2):mg_frag%ie(2), mg_frag%is(3):mg_frag%ie(3))
    write(*,'(1x,a,i0,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        hartree trace: before-vh-allreduce rank=", dg_frag%id_frag, &
      " count=", size(Vh_work%f), " work-bounds=", lbound(Vh_work%f,1), ":", ubound(Vh_work%f,1), ",", &
      lbound(Vh_work%f,2), ":", ubound(Vh_work%f,2), ",", lbound(Vh_work%f,3), ":", ubound(Vh_work%f,3), ",", &
      " local-dest=", lbound(Vh%f,1), ":", ubound(Vh%f,1), ",", lbound(Vh%f,2), ":", ubound(Vh%f,2), ",", &
      lbound(Vh%f,3), ":", ubound(Vh%f,3), ","
    flush(6)
    if (lbound(Vh%f,1) < lbound(Vh_work%f,1) .or. ubound(Vh%f,1) > ubound(Vh_work%f,1) .or. &
        lbound(Vh%f,2) < lbound(Vh_work%f,2) .or. ubound(Vh%f,2) > ubound(Vh_work%f,2) .or. &
        lbound(Vh%f,3) < lbound(Vh_work%f,3) .or. ubound(Vh%f,3) > ubound(Vh_work%f,3)) then
      write(*,'(1x,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        [FATAL] Vh local bounds outside work rank=", dg_frag%id_frag, &
        " work=", lbound(Vh_work%f,1), ":", ubound(Vh_work%f,1), ",", lbound(Vh_work%f,2), ":", ubound(Vh_work%f,2), ",", &
        lbound(Vh_work%f,3), ":", ubound(Vh_work%f,3), ",", " local=", lbound(Vh%f,1), ":", ubound(Vh%f,1), ",", &
        lbound(Vh%f,2), ":", ubound(Vh%f,2), ",", lbound(Vh%f,3), ":", ubound(Vh%f,3), ","
      flush(6)
      stop 'poisson_dg_distributed: Vh local bounds outside work'
    end if
    call comm_summation(Vh_work%f, Vh_reduce%f, size(Vh_work%f), dg_frag%icomm_frag)
    Vh%f(lbound(Vh%f,1):ubound(Vh%f,1), lbound(Vh%f,2):ubound(Vh%f,2), lbound(Vh%f,3):ubound(Vh%f,3)) = &
      Vh_reduce%f(lbound(Vh%f,1):ubound(Vh%f,1), lbound(Vh%f,2):ubound(Vh%f,2), lbound(Vh%f,3):ubound(Vh%f,3))
    write(*,'(1x,a,i0)') "        hartree trace: after-vh-allreduce rank=", dg_frag%id_frag
    flush(6)

    write(*,'(1x,a,i0,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        hartree trace: before-rhog-allreduce rank=", dg_frag%id_frag, &
      " count=", size(rhoG_work), " work-bounds=", lbound(rhoG_work,1), ":", ubound(rhoG_work,1), ",", &
      lbound(rhoG_work,2), ":", ubound(rhoG_work,2), ",", lbound(rhoG_work,3), ":", ubound(rhoG_work,3), ",", &
      " local-dest=", lbound(poisson%zrhoG_ele,1), ":", ubound(poisson%zrhoG_ele,1), ",", &
      lbound(poisson%zrhoG_ele,2), ":", ubound(poisson%zrhoG_ele,2), ",", lbound(poisson%zrhoG_ele,3), ":", ubound(poisson%zrhoG_ele,3), ","
    flush(6)
    if (lbound(poisson%zrhoG_ele,1) < lbound(rhoG_work,1) .or. ubound(poisson%zrhoG_ele,1) > ubound(rhoG_work,1) .or. &
        lbound(poisson%zrhoG_ele,2) < lbound(rhoG_work,2) .or. ubound(poisson%zrhoG_ele,2) > ubound(rhoG_work,2) .or. &
        lbound(poisson%zrhoG_ele,3) < lbound(rhoG_work,3) .or. ubound(poisson%zrhoG_ele,3) > ubound(rhoG_work,3)) then
      write(*,'(1x,a,i0,a,3(i0,a,i0,a),a,3(i0,a,i0,a))') "        [FATAL] rhoG local bounds outside work rank=", dg_frag%id_frag, &
        " work=", lbound(rhoG_work,1), ":", ubound(rhoG_work,1), ",", lbound(rhoG_work,2), ":", ubound(rhoG_work,2), ",", &
        lbound(rhoG_work,3), ":", ubound(rhoG_work,3), ",", " local=", lbound(poisson%zrhoG_ele,1), ":", ubound(poisson%zrhoG_ele,1), ",", &
        lbound(poisson%zrhoG_ele,2), ":", ubound(poisson%zrhoG_ele,2), ",", lbound(poisson%zrhoG_ele,3), ":", ubound(poisson%zrhoG_ele,3), ","
      flush(6)
      stop 'poisson_dg_distributed: rhoG local bounds outside work'
    end if
    call comm_summation(rhoG_work, rhoG_reduce, size(rhoG_work), dg_frag%icomm_frag)
    poisson%zrhoG_ele(lbound(poisson%zrhoG_ele,1):ubound(poisson%zrhoG_ele,1), &
                      lbound(poisson%zrhoG_ele,2):ubound(poisson%zrhoG_ele,2), &
                      lbound(poisson%zrhoG_ele,3):ubound(poisson%zrhoG_ele,3)) = &
      rhoG_reduce(lbound(poisson%zrhoG_ele,1):ubound(poisson%zrhoG_ele,1), &
                  lbound(poisson%zrhoG_ele,2):ubound(poisson%zrhoG_ele,2), &
                  lbound(poisson%zrhoG_ele,3):ubound(poisson%zrhoG_ele,3))
    write(*,'(1x,a,i0)') "        hartree trace: after-rhog-allreduce rank=", dg_frag%id_frag
    flush(6)

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

    call deallocate_fragment_poisson_periodic(poisson_frag)
    call finalize_fragment_parallel(info_frag)
    deallocate(rho_work%f)
    deallocate(rho_reduce%f)
    deallocate(Vh_work%f)
    deallocate(Vh_reduce%f)
    deallocate(rhoG_work)
    deallocate(rhoG_reduce)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: bridge-exit"
      flush(6)
    end if
  end subroutine hartree_dg_via_poisson_periodic

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
    real(8) :: t_stage0, t_stage1
    logical :: use_hse_sr_hartree
    logical :: use_fragment_local_poisson
    character(16) :: env_hse_sr
    character(16) :: env_local_poisson
    integer :: env_status

    ! Local work arrays for the kz-slab owned by this rank
    complex(8), allocatable :: ff1(:,:,:)      ! (mg%is(1):ie(1), mg%is(2):ie(2), nkz_local)
    complex(8), allocatable :: ff2(:,:,:)      ! same shape as ff1

    ! Partial contributions on the full real-space / G-space grids
    real(8),    allocatable :: Vh_partial(:,:,:)      ! (mg range)
    complex(8), allocatable :: rhoG_partial(:,:,:)    ! (mg range)

    env_local_poisson = ''
    use_fragment_local_poisson = .false.
    call get_environment_variable('SALMON_USE_FRAGMENT_LOCAL_POISSON', env_local_poisson, status=env_status)
    if (env_status == 0) then
      select case(trim(adjustl(env_local_poisson)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        use_fragment_local_poisson = .true.
      end select
    end if
    ! The direct distributed DFT path assumes every rank sees the full mg/rho box
    ! on the communicator it reduces over. DG fragment RT currently provides
    ! fragment-subgroup/local-grid data instead, so fall back to the bridge path.
    use_fragment_local_poisson = use_fragment_local_poisson .or. dg_frag%icomm_frag /= dg_frag%icomm
    use_fragment_local_poisson = use_fragment_local_poisson .or. dg_frag%isize_frag /= dg_frag%isize
    use_fragment_local_poisson = use_fragment_local_poisson .or. lbound(rho%f,1) /= mg%is(1)
    use_fragment_local_poisson = use_fragment_local_poisson .or. ubound(rho%f,1) /= mg%ie(1)
    use_fragment_local_poisson = use_fragment_local_poisson .or. lbound(rho%f,2) /= mg%is(2)
    use_fragment_local_poisson = use_fragment_local_poisson .or. ubound(rho%f,2) /= mg%ie(2)
    use_fragment_local_poisson = use_fragment_local_poisson .or. lbound(rho%f,3) /= mg%is(3)
    use_fragment_local_poisson = use_fragment_local_poisson .or. ubound(rho%f,3) /= mg%ie(3)
    if (use_fragment_local_poisson) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a,4(a,i0))') "        hartree trace: using-bridge-fallback", &
          " id=", dg_frag%id, " isize=", dg_frag%isize, " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag
        flush(6)
      end if
      call hartree_dg_via_poisson_periodic(lg, mg, fg, poisson, dg_frag, rho, Vh)
      return
    end if

    ! Check SALMON_HSE_SR_HARTREE environment variable (same logic as hartree_sub::hartree)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: entered-subroutine"
      write(*,'(1x,a,4(a,i0))') "        hartree trace: rank-state", &
        " id=", dg_frag%id, " isize=", dg_frag%isize, " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag
      flush(6)
    end if
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
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=after-dims"
      flush(6)
    end if

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
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=after-slab-partition"
      flush(6)
    end if

    ! Allocate local slab work arrays and global partial arrays
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-alloc"
      flush(6)
    end if
    allocate(ff1(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), nkz_local))
    allocate(ff2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), nkz_local))
    allocate(Vh_partial   (mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(rhoG_partial (mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=after-alloc"
      flush(6)
    end if

    ff1          = (0.0d0, 0.0d0)
    ff2          = (0.0d0, 0.0d0)
    Vh_partial   = 0.0d0
    rhoG_partial = (0.0d0, 0.0d0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,3(a,i0),4(a,i0))') "        hartree trace: entry", &
        " Nx=", Nx, " Ny=", Ny, " Nz=", Nz, " nkz_local=", nkz_local, " nkz_actual=", nkz_actual, &
        " kz_start=", kz_start, " kz_end=", kz_end
      flush(6)
    end if

    ! =======================================================================
    ! FORWARD 3D DFT: rho(r) -> rho_G for local kz slab
    ! No communication is needed: every rank already holds the full rho%f.
    ! =======================================================================
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-forward-dft"
      flush(6)
    end if

    ! -- Step 1: z-DFT --
    ! ff1(ix,iy,kz_loc) = sum_iz  egzc(kz, iz) * rho(ix, iy, iz)
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-z-dft"
      flush(6)
    end if
    !$omp parallel do private(kz_loc, kz, iy, ix)
    do kz_loc = 1, nkz_actual
      kz = kz_start + kz_loc - 1
      do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          if (dg_frag%id == 0 .and. kz_loc == 1 .and. iy == mg%is(2) .and. ix == mg%is(1)) then
!$omp critical
            write(*,'(1x,a,3(a,i0),a,6(i0,a),a,6(i0,a),a,4(i0,a),a,l1,a,l1)') "        hartree trace: z-dft first-point", &
              " kz=", kz, " iy=", iy, " ix=", ix, " mg=", mg%is(1), ":", mg%ie(1), ",", mg%is(2), ":", mg%ie(2), ",", &
              mg%is(3), ":", mg%ie(3), ",", " rho=", lbound(rho%f,1), ":", ubound(rho%f,1), ",", lbound(rho%f,2), ":", ubound(rho%f,2), ",", &
              lbound(rho%f,3), ":", ubound(rho%f,3), ",", " egzc=", lbound(fg%egzc,1), ":", ubound(fg%egzc,1), ",", &
              lbound(fg%egzc,2), ":", ubound(fg%egzc,2), ",", " kz_in=", &
              (kz >= lbound(fg%egzc,1) .and. kz <= ubound(fg%egzc,1)), " iz_ok=", &
              (mg%is(3) >= lbound(fg%egzc,2) .and. mg%ie(3) <= ubound(fg%egzc,2))
            flush(6)
!$omp end critical
          end if
          ff1(ix, iy, kz_loc) = sum(fg%egzc(kz, :) * dcmplx(rho%f(ix, iy, :)))
          if (dg_frag%id == 0 .and. kz_loc == 1 .and. iy == mg%is(2) .and. ix == mg%is(1)) then
!$omp critical
            write(*,'(1x,a,2(a,i0),a)') "        hartree trace: z-dft first-point", &
              " kz=", kz, " iy=", iy, " after-sum"
            flush(6)
!$omp end critical
          end if
        end do
      end do
    end do
    !$omp end parallel do
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-z-dft dt=", t_stage1 - t_stage0
      flush(6)
    end if

    ! -- Step 2: y-DFT --
    ! ff2(ix,ky,kz_loc) = sum_iy  egyc(ky, iy) * ff1(ix, iy, kz_loc)
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-y-dft"
      flush(6)
    end if
    !$omp parallel do private(kz_loc, ky, ix)
    do kz_loc = 1, nkz_actual
      do ky = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          ff2(ix, ky, kz_loc) = sum(fg%egyc(ky, :) * ff1(ix, :, kz_loc))
        end do
      end do
    end do
    !$omp end parallel do
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-y-dft dt=", t_stage1 - t_stage0
      flush(6)
    end if

    ! -- Step 3: x-DFT + normalize --
    ! ff1(kx,ky,kz_loc) = sum_ix  egxc(kx, ix) * ff2(ix, ky, kz_loc) / N
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-x-dft"
      flush(6)
    end if
    !$omp parallel do private(kz_loc, ky, kx)
    do kz_loc = 1, nkz_actual
      do ky = mg%is(2), mg%ie(2)
        do kx = mg%is(1), mg%ie(1)
          ff1(kx, ky, kz_loc) = sum(fg%egxc(kx, :) * ff2(:, ky, kz_loc)) * inv_N
        end do
      end do
    end do
    !$omp end parallel do
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-x-dft dt=", t_stage1 - t_stage0
      flush(6)
    end if

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
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-forward-dft dt=", t_stage1 - t_stage0
      flush(6)
    end if

    ! =======================================================================
    ! COULOMB KERNEL: Vh(G) = coef(G) * rho(G)
    ! When SALMON_HSE_SR_HARTREE is set, apply the HSE short-range factor
    ! (1 - exp(-|G|^2 / (4*omega^2))) to match hartree_sub::hartree behaviour.
    ! =======================================================================
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-kernel"
      flush(6)
    end if
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
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-kernel dt=", t_stage1 - t_stage0
      flush(6)
    end if

    ! =======================================================================
    ! INVERSE 3D DFT: Vh(G) -> partial Vh(r) contribution from local kz slab
    ! =======================================================================
    call cpu_time(t_stage0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-inverse-dft"
      flush(6)
    end if

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
    call cpu_time(t_stage1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        hartree trace: stage=after-inverse-dft dt=", t_stage1 - t_stage0
      flush(6)
    end if

    ! =======================================================================
    ! ALLREDUCE: sum partial Vh and rho(G) contributions across all ranks.
    ! Each rank owns a disjoint kz slab, so MPI_SUM = MPI_Gatherall.
    ! N_total = mg%num product, matching the allocation of both arrays.
    ! =======================================================================
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=before-allreduce-vh"
      flush(6)
    end if
    call comm_summation(Vh_partial,   Vh%f,               N_total, dg_frag%icomm)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=after-allreduce-vh"
      write(*,'(1x,a)') "        hartree trace: stage=before-allreduce-rhog"
      flush(6)
    end if
    call comm_summation(rhoG_partial, poisson%zrhoG_ele,  N_total, dg_frag%icomm)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: stage=after-allreduce-rhog"
      flush(6)
    end if

    deallocate(ff1, ff2, Vh_partial, rhoG_partial)

    ! =======================================================================
    ! WALL POTENTIAL at z boundaries (yn_put_wall_z_boundary='y').
    ! Applied after Vh is assembled, matching hartree_sub::hartree behaviour.
    ! dg_frag%hgs(3) is the grid spacing in z.
    ! =======================================================================
    if (yn_put_wall_z_boundary == 'y') then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "        hartree trace: stage=before-wall"
        flush(6)
      end if
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
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "        hartree trace: stage=after-wall"
        flush(6)
      end if
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        hartree trace: exit"
      flush(6)
    end if

  end subroutine hartree_dg_distributed

end module poisson_dg_distributed
