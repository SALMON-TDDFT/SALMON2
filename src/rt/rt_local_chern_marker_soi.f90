#include "config.h"

module rt_local_chern_marker_soi
  use math_constants, only: pi, zi
  use salmon_global, only: yn_spinorbit
  use structures, only: s_rgrid, s_orbital, s_dft_system, s_parallel_info, &
                        allocate_orbital_real, allocate_orbital_complex, deallocate_orbital
  use communication, only: comm_summation, comm_bcast
  use checkpoint_restart_sub, only: read_wavefunction
  use eigen_subdiag_sub, only: eigen_zheev
#ifdef USE_SCALAPACK
  use eigen_scalapack, only: eigen_pzheevd, eigen_pzheevd_red_mem_xlocal
#endif
  implicit none
  private

  public :: compute_local_chern_marker_from_orbital
  public :: compute_local_chern_marker_from_rt_checkpoint

  logical, save :: phase_cache_initialized = .false.
  integer, save :: phase_is(3) = 0, phase_ie(3) = -1
  real(8), save :: phase_b1_cache(3) = 0d0, phase_b2_cache(3) = 0d0
  complex(8), allocatable, save :: phase1_x(:), phase1_y(:), phase1_z(:)
  complex(8), allocatable, save :: phase2_x(:), phase2_y(:), phase2_z(:)

contains

  ! This implementation assumes a gapped system with a well-defined occupied
  ! subspace. It is not intended for metals or smeared/non-idempotent
  ! occupations, because the occupied manifold is constructed from a sharp
  ! occupation cutoff and then Lowdin orthonormalized.
  subroutine compute_local_chern_marker_from_rt_checkpoint(idir, lg, mg, system, info, marker, occ_eps, is_self_checkpoint)
    implicit none
    character(*), intent(in) :: idir
    type(s_rgrid), intent(in) :: lg, mg
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    real(8), intent(out) :: marker(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in), optional :: occ_eps
    logical, intent(in), optional :: is_self_checkpoint

    type(s_orbital) :: spsi
    logical :: is_self

    is_self = .false.
    if (present(is_self_checkpoint)) is_self = is_self_checkpoint

    if (system%if_real_orbital) then
      call allocate_orbital_real(system%nspin, mg, info, spsi)
    else
      call allocate_orbital_complex(system%nspin, mg, info, spsi)
    end if

    call read_wavefunction(idir, lg, mg, system, info, spsi, system%nk, system%no, system%if_real_orbital, is_self)
    call compute_local_chern_marker_from_orbital(mg, system, info, spsi, marker, occ_eps)
    call deallocate_orbital(spsi)
  end subroutine compute_local_chern_marker_from_rt_checkpoint


  subroutine compute_local_chern_marker_from_orbital(mg, system, info, psi, marker, occ_eps)
    implicit none

    type(s_rgrid), intent(in) :: mg
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_orbital), intent(in) :: psi
    real(8), intent(out) :: marker(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in), optional :: occ_eps

    integer :: ik, ispin, ispinor, io, iocc, jocc, ix, iy, iz, im, p, owner, ig
    integer :: nocc, nocc_local, iocc_g, jocc_g, nblk, nloc, nloc_spin
    integer :: jblk_s, jblk_e, nblk_chunk, jloc
    integer, allocatable :: occ_idx(:)
    integer, allocatable :: occ_owner(:), occ_pos_owner(:), local_occ_glob(:), local_occ_io(:)
    integer, allocatable :: owner_blk_s(:), owner_blk_e(:), owner_nblk(:)
    real(8), allocatable :: occ_w(:)
    real(8), allocatable :: local_occ_w(:)
    complex(8), allocatable, target :: zocc(:,:,:,:,:), zt1(:,:,:,:,:), zt2(:,:,:,:,:)
    complex(8), pointer :: zblk(:,:,:,:,:)
    complex(8), pointer :: zocc2d(:,:), zt12d(:,:), zt22d(:,:)
    complex(8), allocatable :: zrhs1(:,:), zrhs2(:,:), ztmp1(:,:), ztmp2(:,:)
    complex(8), allocatable :: g12_blk(:,:), g12_fwd_blk(:,:), g21_blk(:,:), w12(:,:), w21(:,:)
    complex(8), allocatable :: s1_row(:,:), s2_row(:,:), s1_row_sum(:,:), s2_row_sum(:,:)
    complex(8), allocatable :: g12_row(:,:), g12_row_sum(:,:)
    complex(8) :: phase1, phase2, zterm12, zterm21
    real(8) :: eps_occ, rvec(3)
    real(8) :: t_s12_copy, t_s12_pack, t_s12_gemm, t_s12_accum, t0, t1
    real(8) :: t_inv_s1, t_inv_s2, t_zt12
    real(8) :: b1(3), b2(3)
    logical :: use_complex
    real(8), parameter :: eps_ortho = 1.0d-12
    integer, parameter :: s12_blk_chunk = 96
    integer, parameter :: marker_blk_chunk = 64
    real(8), allocatable :: marker_local(:,:,:), marker_sum(:,:,:)
    complex(8), allocatable :: phase1x_l(:), phase1y_l(:), phase1z_l(:)
    complex(8), allocatable :: phase2x_l(:), phase2y_l(:), phase2z_l(:)

    call assert_soi_lcm_ready(psi)

    ! SOI review point: SOI path should treat psi%zwf as the primary orbital representation.
    use_complex = allocated(psi%zwf)
    eps_occ = 1.0d-10
    if (present(occ_eps)) eps_occ = occ_eps

    b1(:) = system%primitive_b(:,1)
    b2(:) = system%primitive_b(:,2)

    marker(:,:,:) = 0.0d0

    do ik = info%ik_s, info%ik_e
      do ispin = 1, 1
        do im = info%im_s, info%im_e

          ! SOI review point: verify that the per-ispin rocc slice still matches the intended
          ! spinor occupied manifold for noncollinear spin-orbit runs.
          call build_occupied_list_soi(ik, ispin, eps_occ, nocc, occ_idx, occ_w)
          if (nocc <= 0) cycle

          call build_occ_distribution_cache(nocc, occ_idx, occ_w, info%id_o, occ_owner, occ_pos_owner, &
            local_occ_glob, local_occ_io, local_occ_w, owner_blk_s, owner_blk_e, owner_nblk, nocc_local)
          allocate(zocc(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, max(1,nocc_local)))
          allocate(marker_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          allocate(marker_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          marker_local(:,:,:) = 0d0
          marker_sum(:,:,:) = 0d0
          nloc = size(marker_local,1) * size(marker_local,2) * size(marker_local,3)
          nloc_spin = nloc * system%nspin
          call copy_occupied_to_temp(ik, ispin, im, nocc, occ_idx, zocc)
          do p = 1, nocc_local
            zocc(:,:,:,:,p) = sqrt(max(0.0d0, local_occ_w(p))) * zocc(:,:,:,:,p)
          end do
          call lowdin_orthonormalize_occupied(nocc, zocc, eps_ortho)

          allocate(zt1(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, max(1,nocc_local)))
          allocate(zt2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, max(1,nocc_local)))
          allocate(s1_row(max(1,nocc_local),nocc), s2_row(max(1,nocc_local),nocc))
          allocate(s1_row_sum(max(1,nocc_local),nocc), s2_row_sum(max(1,nocc_local),nocc))
          allocate(g12_row(max(1,nocc_local),nocc), g12_row_sum(max(1,nocc_local),nocc))

          s1_row(:,:) = (0.0d0, 0.0d0)
          s2_row(:,:) = (0.0d0, 0.0d0)
          t_s12_copy = 0d0
          t_s12_pack = 0d0
          t_s12_gemm = 0d0
          t_s12_accum = 0d0
          t_inv_s1 = 0d0
          t_inv_s2 = 0d0
          t_zt12 = 0d0
          call cpu_time(t0)
          allocate(phase1x_l(mg%is(1):mg%ie(1)))
          allocate(phase1y_l(mg%is(2):mg%ie(2)))
          allocate(phase1z_l(mg%is(3):mg%ie(3)))
          allocate(phase2x_l(mg%is(1):mg%ie(1)))
          allocate(phase2y_l(mg%is(2):mg%ie(2)))
          allocate(phase2z_l(mg%is(3):mg%ie(3)))
          do ix = mg%is(1), mg%ie(1)
            phase1x_l(ix) = exp(-zi * b1(1) * mg%coordinate(ix,1))
            phase2x_l(ix) = exp(-zi * b2(1) * mg%coordinate(ix,1))
          end do
          do iy = mg%is(2), mg%ie(2)
            phase1y_l(iy) = exp(-zi * b1(2) * mg%coordinate(iy,2))
            phase2y_l(iy) = exp(-zi * b2(2) * mg%coordinate(iy,2))
          end do
          do iz = mg%is(3), mg%ie(3)
            phase1z_l(iz) = exp(-zi * b1(3) * mg%coordinate(iz,3))
            phase2z_l(iz) = exp(-zi * b2(3) * mg%coordinate(iz,3))
          end do
          call cpu_time(t1)
          if (info%id_rko == 0) then
            write(*,'(a,1x,es12.4)') 'LCM S1/S2 timing: phase_precompute=', t1 - t0
          end if
          zocc2d(1:nloc_spin,1:max(1,nocc_local)) => zocc(:,:,:,:,1:max(1,nocc_local))
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            call cpu_time(t0)
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
            call cpu_time(t1)
            t_s12_copy = t_s12_copy + (t1 - t0)
            do jblk_s = 1, nblk, s12_blk_chunk
              jblk_e = min(nblk, jblk_s + s12_blk_chunk - 1)
              nblk_chunk = jblk_e - jblk_s + 1
              allocate(zrhs1(max(1,nloc_spin), nblk_chunk))
              allocate(ztmp1(max(1,nocc_local), nblk_chunk))
              call cpu_time(t0)
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,jocc,jloc,phase1) schedule(static)
              do iz = mg%is(3), mg%ie(3)
                do iy = mg%is(2), mg%ie(2)
                  do ix = mg%is(1), mg%ie(1)
                    do ispinor = 1, system%nspin
                      ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                      phase1 = phase1x_l(ix) * phase1y_l(iy) * phase1z_l(iz)
                      do jloc = 1, nblk_chunk
                        jocc = jblk_s + jloc - 1
                        zrhs1(ig,jloc) = phase1 * zblk(ix,iy,iz,ispinor,jocc)
                      end do
                    end do
                  end do
                end do
              end do
              call cpu_time(t1)
              t_s12_pack = t_s12_pack + (t1 - t0)
              call cpu_time(t0)
              call zgemm('C', 'N', nocc_local, nblk_chunk, nloc_spin, cmplx(system%hvol, 0.0d0, kind=8), zocc2d, max(1,nloc_spin), &
                zrhs1, max(1,nloc_spin), (0.0d0, 0.0d0), ztmp1, max(1,nocc_local))
              call cpu_time(t1)
              t_s12_gemm = t_s12_gemm + (t1 - t0)
              call cpu_time(t0)
              do p = 1, nocc_local
                do jloc = 1, nblk_chunk
                  jocc_g = iocc + (jblk_s + jloc - 1) - 1
                  s1_row(p,jocc_g) = s1_row(p,jocc_g) + ztmp1(p,jloc)
                end do
              end do
              call cpu_time(t1)
              t_s12_accum = t_s12_accum + (t1 - t0)
              call cpu_time(t0)
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,jocc,jloc,phase2) schedule(static)
              do iz = mg%is(3), mg%ie(3)
                do iy = mg%is(2), mg%ie(2)
                  do ix = mg%is(1), mg%ie(1)
                    do ispinor = 1, system%nspin
                      ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                      phase2 = phase2x_l(ix) * phase2y_l(iy) * phase2z_l(iz)
                      do jloc = 1, nblk_chunk
                        jocc = jblk_s + jloc - 1
                        zrhs1(ig,jloc) = phase2 * zblk(ix,iy,iz,ispinor,jocc)
                      end do
                    end do
                  end do
                end do
              end do
              call cpu_time(t1)
              t_s12_pack = t_s12_pack + (t1 - t0)
              call cpu_time(t0)
              call zgemm('C', 'N', nocc_local, nblk_chunk, nloc_spin, cmplx(system%hvol, 0.0d0, kind=8), zocc2d, max(1,nloc_spin), &
                zrhs1, max(1,nloc_spin), (0.0d0, 0.0d0), ztmp1, max(1,nocc_local))
              call cpu_time(t1)
              t_s12_gemm = t_s12_gemm + (t1 - t0)
              call cpu_time(t0)
              do p = 1, nocc_local
                do jloc = 1, nblk_chunk
                  jocc_g = iocc + (jblk_s + jloc - 1) - 1
                  s2_row(p,jocc_g) = s2_row(p,jocc_g) + ztmp1(p,jloc)
                end do
              end do
              call cpu_time(t1)
              t_s12_accum = t_s12_accum + (t1 - t0)
              deallocate(ztmp1, zrhs1)
            end do
            deallocate(zblk)
          end do
          call comm_summation(s1_row, s1_row_sum, size(s1_row), info%icomm_r)
          call comm_summation(s2_row, s2_row_sum, size(s2_row), info%icomm_r)
          if (info%id_rko == 0) then
            write(*,'(a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
              'LCM S1/S2 timing: copy=', t_s12_copy, 'pack=', t_s12_pack, 'gemm=', t_s12_gemm, 'accum=', t_s12_accum
          end if
          call cpu_time(t0)
          call build_transposed_inverse_coefficients_rowwise_checked(nocc, occ_idx, s1_row_sum, 'S1')
          call cpu_time(t1)
          t_inv_s1 = t1 - t0
          if (info%id_rko == 0) write(*,'(a,1x,es12.4)') 'LCM timing: inv(S1)=', t_inv_s1
          call cpu_time(t0)
          call build_transposed_inverse_coefficients_rowwise_checked(nocc, occ_idx, s2_row_sum, 'S2')
          call cpu_time(t1)
          t_inv_s2 = t1 - t0
          if (info%id_rko == 0) write(*,'(a,1x,es12.4)') 'LCM timing: inv(S2)=', t_inv_s2

          zt1(:,:,:,:,:) = (0.0d0, 0.0d0)
          zt2(:,:,:,:,:) = (0.0d0, 0.0d0)
          zt12d(1:nloc_spin,1:max(1,nocc_local)) => zt1(:,:,:,:,1:max(1,nocc_local))
          zt22d(1:nloc_spin,1:max(1,nocc_local)) => zt2(:,:,:,:,1:max(1,nocc_local))
          call cpu_time(t0)
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
            allocate(zrhs1(max(1,nloc_spin), nblk))
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,jocc,phase1) schedule(static)
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  do ispinor = 1, system%nspin
                    ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                    phase1 = phase1x_l(ix) * phase1y_l(iy) * phase1z_l(iz)
                    do jocc = 1, nblk
                      zrhs1(ig,jocc) = phase1 * zblk(ix,iy,iz,ispinor,jocc)
                    end do
                  end do
                end do
              end do
            end do
            call zgemm('N', 'T', nloc_spin, nocc_local, nblk, (1.0d0, 0.0d0), zrhs1, max(1,nloc_spin), &
              s1_row_sum(:,iocc:iocc+nblk-1), max(1,nocc_local), (1.0d0, 0.0d0), zt12d, max(1,nloc_spin))
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,jocc,phase2) schedule(static)
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  do ispinor = 1, system%nspin
                    ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                    phase2 = phase2x_l(ix) * phase2y_l(iy) * phase2z_l(iz)
                    do jocc = 1, nblk
                      zrhs1(ig,jocc) = phase2 * zblk(ix,iy,iz,ispinor,jocc)
                    end do
                  end do
                end do
              end do
            end do
            call zgemm('N', 'T', nloc_spin, nocc_local, nblk, (1.0d0, 0.0d0), zrhs1, max(1,nloc_spin), &
              s2_row_sum(:,iocc:iocc+nblk-1), max(1,nocc_local), (1.0d0, 0.0d0), zt22d, max(1,nloc_spin))
            deallocate(zrhs1)
            deallocate(zblk)
          end do
          call cpu_time(t1)
          t_zt12 = t1 - t0
          if (info%id_rko == 0) write(*,'(a,1x,es12.4)') 'LCM timing: ZT1/ZT2=', t_zt12

          deallocate(phase2z_l, phase2y_l, phase2x_l, phase1z_l, phase1y_l, phase1x_l)

          g12_row(:,:) = (0.0d0, 0.0d0)
          zt12d(1:nloc_spin,1:max(1,nocc_local)) => zt1(:,:,:,:,1:max(1,nocc_local))
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zt2, iocc, jocc, zblk)
            allocate(zrhs1(max(1,nloc_spin), nblk), ztmp1(max(1,nocc_local), nblk))
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,jocc) schedule(static)
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  do ispinor = 1, system%nspin
                    ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                    do jocc = 1, nblk
                      zrhs1(ig,jocc) = zblk(ix,iy,iz,ispinor,jocc)
                    end do
                  end do
                end do
              end do
            end do
            call zgemm('C', 'N', nocc_local, nblk, nloc_spin, cmplx(system%hvol, 0.0d0, kind=8), zt12d, max(1,nloc_spin), &
              zrhs1, max(1,nloc_spin), (0.0d0, 0.0d0), ztmp1, max(1,nocc_local))
            do p = 1, nocc_local
              do jocc = 1, nblk
                jocc_g = iocc + jocc - 1
                g12_row(p,jocc_g) = g12_row(p,jocc_g) + ztmp1(p,jocc)
              end do
            end do
            deallocate(ztmp1, zrhs1)
            deallocate(zblk)
          end do
          call comm_summation(g12_row, g12_row_sum, size(g12_row), info%icomm_r)

          marker_local(:,:,:) = 0d0
          allocate(zrhs2(max(1,nloc_spin), max(1,nocc_local)), ztmp2(max(1,nloc_spin), max(1,nocc_local)))
!$omp parallel do collapse(4) private(ix,iy,iz,ig,ispinor,p) schedule(static)
          do iz = mg%is(3), mg%ie(3)
            do iy = mg%is(2), mg%ie(2)
              do ix = mg%is(1), mg%ie(1)
                do ispinor = 1, system%nspin
                  ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                  do p = 1, nocc_local
                    zrhs2(ig,p) = zt1(ix,iy,iz,ispinor,p)
                    ztmp2(ig,p) = zt2(ix,iy,iz,ispinor,p)
                  end do
                end do
              end do
            end do
          end do
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
            do jblk_s = 1, nblk, marker_blk_chunk
              jblk_e = min(nblk, jblk_s + marker_blk_chunk - 1)
              nblk_chunk = jblk_e - jblk_s + 1
              allocate(g12_blk(max(1,nblk_chunk), nocc))
              call get_owned_row_block(nocc, occ_idx, owner, g12_row_sum, iocc + jblk_s - 1, iocc + jblk_e - 1, g12_blk)
              allocate(g12_fwd_blk(max(1,nocc_local), nblk_chunk), g21_blk(max(1,nocc_local), nblk_chunk), &
                w12(max(1,nloc), nblk_chunk), w21(max(1,nloc), nblk_chunk))
              do p = 1, nocc_local
                do jloc = 1, nblk_chunk
                  g12_fwd_blk(p,jloc) = g12_blk(jloc,local_occ_glob(p))
                  g21_blk(p,jloc) = conjg(g12_blk(jloc,local_occ_glob(p)))
                end do
              end do
              call zgemm('N', 'N', nloc_spin, nblk_chunk, nocc_local, (1.0d0, 0.0d0), zrhs2, max(1,nloc_spin), g12_fwd_blk, &
                max(1,nocc_local), (0.0d0, 0.0d0), w12, max(1,nloc_spin))
              call zgemm('N', 'N', nloc_spin, nblk_chunk, nocc_local, (1.0d0, 0.0d0), ztmp2, max(1,nloc_spin), g21_blk, &
                max(1,nocc_local), (0.0d0, 0.0d0), w21, max(1,nloc_spin))
!$omp parallel do collapse(3) private(ix,iy,iz,ig,jocc,jloc,zterm12,zterm21,phase1) schedule(static)
              do iz = mg%is(3), mg%ie(3)
                do iy = mg%is(2), mg%ie(2)
                  do ix = mg%is(1), mg%ie(1)
                    zterm12 = (0.0d0, 0.0d0)
                    zterm21 = (0.0d0, 0.0d0)
                    do jloc = 1, nblk_chunk
                      jocc = jblk_s + jloc - 1
                      do ispinor = 1, system%nspin
                        ig = (((ispinor - 1) * (mg%ie(3) - mg%is(3) + 1) + (iz - mg%is(3))) * (mg%ie(2) - mg%is(2) + 1) + (iy - mg%is(2))) * (mg%ie(1) - mg%is(1) + 1) + (ix - mg%is(1)) + 1
                        phase1 = conjg(zblk(ix,iy,iz,ispinor,jocc))
                        zterm12 = zterm12 + w12(ig,jloc) * phase1
                        zterm21 = zterm21 + w21(ig,jloc) * phase1
                      end do
                    end do
                    marker_local(ix,iy,iz) = marker_local(ix,iy,iz) - aimag(zterm12 - zterm21) * system%wtk(ik) / (2.0d0*pi)
                  end do
                end do
              end do
              deallocate(w21, w12, g21_blk, g12_fwd_blk, g12_blk)
            end do
            deallocate(zblk)
          end do
          deallocate(ztmp2, zrhs2)
          nloc = size(marker_local,1) * size(marker_local,2) * size(marker_local,3)
          call comm_summation(marker_local, marker_sum, nloc, info%icomm_o)
          marker(:,:,:) = marker(:,:,:) + marker_sum(:,:,:)

          deallocate(s1_row, s2_row, s1_row_sum, s2_row_sum, g12_row, g12_row_sum)
          deallocate(zt1, zt2, zocc)
          deallocate(marker_local, marker_sum)
          deallocate(occ_idx, occ_w, occ_owner, occ_pos_owner, local_occ_glob, local_occ_io, local_occ_w, &
            owner_blk_s, owner_blk_e, owner_nblk)

        end do
      end do
    end do

    if (info%isize_k > 1) then
      call reduce_marker_over_k(mg, info, marker)
    end if

  contains

    subroutine build_occupied_list_soi(ik0, ispin0, eps_occ0, nocc0, occ_list0, occ_w0)
      implicit none
      integer, intent(in) :: ik0, ispin0
      real(8), intent(in) :: eps_occ0
      integer, intent(out) :: nocc0
      integer, allocatable, intent(out) :: occ_list0(:)
      real(8), allocatable, intent(out) :: occ_w0(:)

      integer :: io0

      nocc0 = 0
      allocate(occ_list0(1), occ_w0(1))
      do io0 = 1, system%no
        if (system%rocc(io0, ik0, 1) > eps_occ0) nocc0 = nocc0 + 1
      end do
      if (nocc0 <= 0) return

      deallocate(occ_list0, occ_w0)
      allocate(occ_list0(nocc0), occ_w0(nocc0))
      nocc0 = 0
      do io0 = 1, system%no
        if (system%rocc(io0, ik0, 1) > eps_occ0) then
          nocc0 = nocc0 + 1
          occ_list0(nocc0) = io0
          occ_w0(nocc0) = system%rocc(io0, ik0, 1)
        end if
      end do
    end subroutine build_occupied_list_soi

    subroutine assert_soi_lcm_ready(psi0)
      implicit none
      type(s_orbital), intent(in) :: psi0

      if (yn_spinorbit /= 'y') then
        stop 'rt_local_chern_marker_soi requires yn_spinorbit=y'
      end if
      if (.not. allocated(psi0%zwf)) then
        stop 'rt_local_chern_marker_soi requires psi%zwf'
      end if
    end subroutine assert_soi_lcm_ready

    ! SOI review point: make the occupied-orbital copy logic explicit for complex spinor data.
    subroutine copy_occupied_to_temp(ik0, ispin0, im0, nocc0, occ_list, zbuf)
      implicit none
      integer, intent(in) :: ik0, ispin0, im0, nocc0
      integer, intent(in) :: occ_list(nocc0)
      complex(8), intent(out), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :, :) :: zbuf

      integer :: io_g, p, nocc_local0

      zbuf(:,:,:,:,:) = (0.0d0, 0.0d0)
      nocc_local0 = size(zbuf,5)

      if (nocc_local0 <= 0) return

      if (use_complex) then
!$omp parallel do private(p,io_g,ispin,ix,iy,iz) schedule(static)
        do p = 1, nocc_local0
          io_g = local_occ_io(p)
          do ispin = 1, system%nspin
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  zbuf(ix,iy,iz,ispin,p) = psi%zwf(ix,iy,iz,ispin,io_g,ik0,im0)
                end do
              end do
            end do
          end do
        end do
      else
        stop 'rt_local_chern_marker_soi requires psi%zwf'
      end if
    end subroutine copy_occupied_to_temp


    subroutine lowdin_orthonormalize_occupied(nocc0, zbuf, eps_rel)
      use salmon_global, only: yn_scalapack
#ifdef USE_SCALAPACK
      use scalapack_module, only: create_gridmap, init_blacs
#endif
      use structures, only: s_parallel_info
      implicit none
      integer, intent(in) :: nocc0
      complex(8), intent(inout), contiguous, dimension(mg%is(1):, mg%is(2):, mg%is(3):, :, :) :: zbuf
      real(8), intent(in) :: eps_rel

      complex(8), allocatable :: s_local(:,:), umat(:,:), xmat_local(:,:), u_scaled_blk(:,:), umat_blk(:,:)
      real(8), allocatable :: eval(:)
      real(8) :: smax
      integer :: ia, ib, ic, owner, ib_s, ib_e, nblk, ia_g, loc_s, loc_e, nloc_cols, ispin
      integer :: ic_s, ic_e, ncol_chunk
      complex(8), allocatable :: zblk2(:,:,:,:,:)
#ifdef USE_SCALAPACK
      type(s_parallel_info) :: info_sl
#endif
      integer, parameter :: eig_blk_chunk = 64

      allocate(s_local(nocc0,nocc0), eval(nocc0))
      s_local(:,:) = (0.0d0, 0.0d0)

      do owner = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_idx, owner, ib_s, ib_e, nblk)
        if (nblk <= 0) cycle
        allocate(zblk2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
        call get_owned_occ_block(nocc0, occ_idx, owner, zbuf, ib_s, ib_e, zblk2)
        do ia = 1, nocc_local
          ia_g = local_occ_glob(ia)
          do ib = 1, nblk
            do ispin = 1, system%nspin
              do iz = mg%is(3), mg%ie(3)
                do iy = mg%is(2), mg%ie(2)
                  do ix = mg%is(1), mg%ie(1)
                    s_local(ia_g,ib_s + ib - 1) = s_local(ia_g,ib_s + ib - 1) + conjg(zbuf(ix,iy,iz,ispin,ia)) * zblk2(ix,iy,iz,ispin,ib) * system%hvol
                  end do
                end do
              end do
            end do
          end do
        end do
        deallocate(zblk2)
      end do
      call get_owner_occ_block(nocc0, occ_idx, info%id_o, loc_s, loc_e, nloc_cols)
      allocate(xmat_local(nocc0,max(1,nloc_cols)))
      xmat_local(:,:) = (0.0d0, 0.0d0)

      if (yn_scalapack == 'y') then
#ifdef USE_SCALAPACK
        info_sl = info
        if (.not. info_sl%flag_blacs_gridinit) then
          call create_gridmap(info_sl)
          call init_blacs(info_sl, nocc0)
        end if
        allocate(umat(nocc0,nocc0))
        call comm_summation(s_local, umat, nocc0*nocc0, info%icomm_r)
        call comm_summation(umat, s_local, nocc0*nocc0, info%icomm_o)
        umat(:,:) = s_local(:,:)
        deallocate(s_local)
        call eigen_pzheevd_red_mem_xlocal(info_sl, umat, eval, loc_s, loc_e, eps_rel, xmat_local)
        deallocate(umat)
#else
        stop "LCM lowdin: yn_scalapack='y' but ScaLAPACK is not enabled."
#endif
      else
        allocate(umat(nocc0,nocc0))
        call comm_summation(s_local, umat, nocc0*nocc0, info%icomm_r)
        call comm_summation(umat, s_local, nocc0*nocc0, info%icomm_o)
        umat(:,:) = s_local(:,:)
        deallocate(s_local)
        call eigen_zheev(umat, eval, umat)
      end if

      smax = maxval(abs(eval))
      if (smax <= 0.0d0) stop 'lowdin_orthonormalize_occupied: invalid overlap spectrum'
      if (yn_scalapack /= 'y' .and. nloc_cols > 0) then
        do ic_s = 1, nocc0, eig_blk_chunk
          ic_e = min(nocc0, ic_s + eig_blk_chunk - 1)
          ncol_chunk = ic_e - ic_s + 1
          allocate(u_scaled_blk(nocc0,ncol_chunk), umat_blk(max(1,nloc_cols),ncol_chunk))
          u_scaled_blk(:,:) = (0.0d0, 0.0d0)
          do ic = ic_s, ic_e
            umat_blk(1:nloc_cols,ic-ic_s+1) = umat(loc_s:loc_e,ic)
            if (eval(ic) > eps_rel * smax) then
              do ia = 1, nocc0
                u_scaled_blk(ia,ic-ic_s+1) = umat(ia,ic) / sqrt(eval(ic))
              end do
            end if
          end do
          call zgemm('N', 'C', nocc0, nloc_cols, ncol_chunk, (1.0d0, 0.0d0), u_scaled_blk, nocc0, &
            umat_blk, max(1,nloc_cols), (1.0d0, 0.0d0), xmat_local, nocc0)
          deallocate(umat_blk, u_scaled_blk)
        end do
      end if
      call apply_right_transform_occ_inplace(nocc0, zbuf, xmat_local)
      if (allocated(umat)) deallocate(umat)
      if (allocated(s_local)) deallocate(s_local)
      deallocate(xmat_local, eval)
    end subroutine lowdin_orthonormalize_occupied


    subroutine apply_right_transform_occ_inplace(nocc0, zbuf, xmat_local)
      implicit none
      integer, intent(in) :: nocc0
      complex(8), intent(inout), contiguous, target, dimension(mg%is(1):, mg%is(2):, mg%is(3):, :, :) :: zbuf
      complex(8), intent(in) :: xmat_local(:,:)

      complex(8), allocatable, target :: zsrc(:,:,:,:,:), zblk_owner(:,:,:,:,:)
      complex(8), allocatable :: xblk(:,:)
      complex(8), pointer :: zsrc2d(:,:), zbuf2d(:,:), zblk2d(:,:)
      integer :: owner, ia_s, ia_e, nblk, nocc_local0, nloc0

      nocc_local0 = size(zbuf,5)
      nloc0 = (mg%ie(1) - mg%is(1) + 1) * (mg%ie(2) - mg%is(2) + 1) * (mg%ie(3) - mg%is(3) + 1) * system%nspin
      allocate(zsrc(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, max(1,nocc_local0)))
      zsrc(:,:,:,:,:) = zbuf(:,:,:,:,:)
      zbuf(:,:,:,:,:) = (0.0d0, 0.0d0)
      zsrc2d(1:nloc0,1:max(1,nocc_local0)) => zsrc(:,:,:,:,1:max(1,nocc_local0))
      zbuf2d(1:nloc0,1:max(1,nocc_local0)) => zbuf(:,:,:,:,1:max(1,nocc_local0))
      do owner = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_idx, owner, ia_s, ia_e, nblk)
        if (nblk <= 0) cycle
        allocate(zblk_owner(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, nblk))
        call get_owned_occ_block(nocc0, occ_idx, owner, zsrc, ia_s, ia_e, zblk_owner)
        zblk2d(1:nloc0,1:nblk) => zblk_owner(:,:,:,:,1:nblk)
        allocate(xblk(nblk, max(1,nocc_local0)))
        xblk(:,:) = xmat_local(ia_s:ia_e,1:max(1,nocc_local0))
        call zgemm('N', 'N', nloc0, nocc_local0, nblk, (1.0d0, 0.0d0), zblk2d, max(1,nloc0), &
          xblk, max(1,nblk), (1.0d0, 0.0d0), zbuf2d, max(1,nloc0))
        deallocate(xblk)
        deallocate(zblk_owner)
      end do
      deallocate(zsrc)
    end subroutine apply_right_transform_occ_inplace

    subroutine build_occ_distribution_cache(nocc0, occ_list0, occ_w0, owner_id, occ_owner0, occ_pos0, &
      local_glob0, local_io0, local_w0, owner_blk_s0, owner_blk_e0, owner_nblk0, nlocal)
      implicit none
      integer, intent(in) :: nocc0, owner_id
      integer, intent(in) :: occ_list0(nocc0)
      real(8), intent(in) :: occ_w0(nocc0)
      integer, allocatable, intent(out) :: occ_owner0(:), occ_pos0(:)
      integer, allocatable, intent(out) :: local_glob0(:), local_io0(:)
      real(8), allocatable, intent(out) :: local_w0(:)
      integer, allocatable, intent(out) :: owner_blk_s0(:), owner_blk_e0(:), owner_nblk0(:)
      integer, intent(out) :: nlocal
      integer :: p, owner, lidx
      integer, allocatable :: owner_count(:)

      allocate(occ_owner0(nocc0), occ_pos0(nocc0))
      allocate(owner_blk_s0(info%isize_o), owner_blk_e0(info%isize_o), owner_nblk0(info%isize_o))
      allocate(owner_count(info%isize_o))

      owner_count(:) = 0
      owner_blk_s0(:) = 1
      owner_blk_e0(:) = 0
      owner_nblk0(:) = 0

      do p = 1, nocc0
        owner = occ_owner_id(occ_list0(p))
        if (owner < 0) stop 'build_occ_distribution_cache: invalid owner'
        occ_owner0(p) = owner
        owner_count(owner+1) = owner_count(owner+1) + 1
        occ_pos0(p) = owner_count(owner+1)
        owner_nblk0(owner+1) = owner_count(owner+1)
        if (owner_count(owner+1) == 1) owner_blk_s0(owner+1) = p
        owner_blk_e0(owner+1) = p
      end do

      nlocal = owner_nblk0(owner_id+1)
      allocate(local_glob0(max(1,nlocal)), local_io0(max(1,nlocal)), local_w0(max(1,nlocal)))
      local_glob0(:) = 1
      local_io0(:) = 1
      local_w0(:) = 0d0
      if (nlocal > 0) then
        do p = 1, nocc0
          if (occ_owner0(p) == owner_id) then
            lidx = occ_pos0(p)
            local_glob0(lidx) = p
            local_io0(lidx) = occ_list0(p)
            local_w0(lidx) = occ_w0(p)
          end if
        end do
      end if

      deallocate(owner_count)
    end subroutine build_occ_distribution_cache

    integer function occ_owner_id(io_g) result(owner_id)
      implicit none
      integer, intent(in) :: io_g
      owner_id = info%id_o
      do owner_id = 0, info%isize_o - 1
        if (io_g >= info%io_s_all(owner_id) .and. io_g <= info%io_e_all(owner_id)) return
      end do
      owner_id = -1
    end function occ_owner_id

    integer function local_occ_index(nocc0, occ_list0, owner_id, plocal) result(pglob)
      implicit none
      integer, intent(in) :: nocc0, occ_list0(nocc0), owner_id, plocal
      integer :: p, count
      count = 0
      pglob = 0
      do p = 1, nocc0
        if (occ_owner_id(occ_list0(p)) == owner_id) then
          count = count + 1
          if (count == plocal) then
            pglob = p
            return
          end if
        end if
      end do
    end function local_occ_index

    integer function local_occ_position(nocc0, occ_list0, io_g, owner_id) result(plocal)
      implicit none
      integer, intent(in) :: nocc0, occ_list0(nocc0), io_g, owner_id
      integer :: p
      plocal = 0
      do p = 1, nocc0
        if (occ_owner_id(occ_list0(p)) == owner_id) then
          plocal = plocal + 1
          if (occ_list0(p) == io_g) return
        end if
      end do
      plocal = 0
    end function local_occ_position

    real(8) function local_occ_weight(nocc0, occ_list0, occ_w0, owner_id, plocal) result(w)
      implicit none
      integer, intent(in) :: nocc0, occ_list0(nocc0), owner_id, plocal
      real(8), intent(in) :: occ_w0(nocc0)
      integer :: pglob
      pglob = local_occ_index(nocc0, occ_list0, owner_id, plocal)
      if (pglob > 0) then
        w = occ_w0(pglob)
      else
        w = 0d0
      end if
    end function local_occ_weight

    integer function local_occ_global_io(nocc0, occ_list0, owner_id, plocal) result(io_g)
      implicit none
      integer, intent(in) :: nocc0, occ_list0(nocc0), owner_id, plocal
      integer :: pglob
      pglob = local_occ_index(nocc0, occ_list0, owner_id, plocal)
      if (pglob > 0) then
        io_g = occ_list0(pglob)
      else
        io_g = 1
      end if
    end function local_occ_global_io

    subroutine get_owner_occ_block(nocc0, occ_list0, owner_id, blk_s, blk_e, nblk)
      implicit none
      integer, intent(in) :: nocc0, occ_list0(nocc0), owner_id
      integer, intent(out) :: blk_s, blk_e, nblk
      integer :: p
      if (allocated(owner_nblk) .and. allocated(owner_blk_s) .and. allocated(owner_blk_e)) then
        if (size(owner_nblk) == info%isize_o .and. size(owner_blk_s) == info%isize_o .and. size(owner_blk_e) == info%isize_o) then
          nblk = owner_nblk(owner_id+1)
          blk_s = owner_blk_s(owner_id+1)
          blk_e = owner_blk_e(owner_id+1)
          return
        end if
      end if
      blk_s = 1
      blk_e = 0
      nblk = 0
      do p = 1, nocc0
        if (occ_owner_id(occ_list0(p)) == owner_id) then
          if (nblk == 0) blk_s = p
          blk_e = p
          nblk = nblk + 1
        end if
      end do
    end subroutine get_owner_occ_block

    subroutine get_owned_occ_block(nocc0, occ_list0, owner_id, zloc, blk_s, blk_e, zblk_out)
      implicit none
      integer, intent(in) :: nocc0, owner_id, blk_s, blk_e
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(in), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :, :) :: zloc
      complex(8), intent(out) :: zblk_out(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), system%nspin, max(1,blk_e-blk_s+1))
      integer :: nblk0, p, locp
      nblk0 = max(0, blk_e - blk_s + 1)
      zblk_out = (0.0d0, 0.0d0)
      if (nblk0 <= 0) return
      if (info%id_o == owner_id) then
        locp = 0
        do p = blk_s, blk_e
          locp = locp + 1
          if (allocated(occ_pos_owner)) then
            zblk_out(:,:,:,:,locp) = zloc(:,:,:,:,occ_pos_owner(p))
          else
            zblk_out(:,:,:,:,locp) = zloc(:,:,:,:,local_occ_position(nocc0, occ_list0, occ_list0(p), owner_id))
          end if
        end do
      end if
      call comm_bcast(zblk_out(:,:,:,:,1:nblk0), info%icomm_o, owner_id)
    end subroutine get_owned_occ_block

    subroutine assemble_owner_row_matrix_rooted(nocc0, occ_list0, row_local, amat)
      implicit none
      integer, intent(in) :: nocc0
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(in) :: row_local(:,:)
      complex(8), allocatable, intent(inout) :: amat(:,:)
      complex(8), allocatable :: row_blk(:,:)
      integer :: owner_id, blk_s, blk_e, nblk, p

      if (info%id_o == 0) then
        if (.not. allocated(amat)) allocate(amat(nocc0,nocc0))
        amat(:,:) = (0.0d0, 0.0d0)
      else
        if (allocated(amat)) deallocate(amat)
      end if
      do owner_id = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_list0, owner_id, blk_s, blk_e, nblk)
        if (nblk <= 0) cycle
        allocate(row_blk(nblk, nocc0))
        row_blk(:,:) = (0.0d0, 0.0d0)
        if (info%id_o == owner_id) then
          do p = 1, nblk
            row_blk(p,:) = row_local(p,:)
          end do
        end if
        call comm_bcast(row_blk, info%icomm_o, owner_id)
        if (info%id_o == 0) amat(blk_s:blk_e,:) = row_blk(1:nblk,:)
        deallocate(row_blk)
      end do
    end subroutine assemble_owner_row_matrix_rooted

    subroutine get_owned_row_block(nocc0, occ_list0, owner_id, row_local, blk_s, blk_e, row_blk_out)
      implicit none
      integer, intent(in) :: nocc0, owner_id, blk_s, blk_e
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(in) :: row_local(:,:)
      complex(8), intent(out) :: row_blk_out(max(1,blk_e-blk_s+1), nocc0)
      integer :: nblk0, p, p0

      nblk0 = max(0, blk_e - blk_s + 1)
      row_blk_out(:,:) = (0.0d0, 0.0d0)
      if (nblk0 <= 0) return
      if (info%id_o == owner_id) then
        p0 = 1
        if (allocated(owner_blk_s)) p0 = blk_s - owner_blk_s(owner_id+1) + 1
        if (p0 < 1 .or. p0 + nblk0 - 1 > size(row_local,1)) then
          stop 'get_owned_row_block: local row offset out of range'
        end if
        do p = 1, nblk0
          row_blk_out(p,:) = row_local(p0 + p - 1,:)
        end do
      end if
      call comm_bcast(row_blk_out(1:nblk0,:), info%icomm_o, owner_id)
    end subroutine get_owned_row_block

    subroutine get_owned_row_block_selected(nocc0, occ_list0, owner_id, row_local, blk_s, blk_e, col_idx, row_blk_out)
      implicit none
      integer, intent(in) :: nocc0, owner_id, blk_s, blk_e
      integer, intent(in) :: occ_list0(nocc0)
      integer, intent(in) :: col_idx(:)
      complex(8), intent(in) :: row_local(:,:)
      complex(8), intent(out) :: row_blk_out(max(1,blk_e-blk_s+1), max(1,size(col_idx)))
      integer :: nblk0, ncol, p, j, p0

      nblk0 = max(0, blk_e - blk_s + 1)
      ncol = size(col_idx)
      row_blk_out(:,:) = (0.0d0, 0.0d0)
      if (nblk0 <= 0 .or. ncol <= 0) return
      if (info%id_o == owner_id) then
        p0 = 1
        if (allocated(owner_blk_s)) p0 = blk_s - owner_blk_s(owner_id+1) + 1
        if (p0 < 1 .or. p0 + nblk0 - 1 > size(row_local,1)) then
          stop 'get_owned_row_block_selected: local row offset out of range'
        end if
        do p = 1, nblk0
          do j = 1, ncol
            row_blk_out(p,j) = row_local(p0 + p - 1,col_idx(j))
          end do
        end do
      end if
      call comm_bcast(row_blk_out(1:nblk0,1:ncol), info%icomm_o, owner_id)
    end subroutine get_owned_row_block_selected

    subroutine build_transposed_inverse_coefficients_rowwise_checked(nocc0, occ_list0, row_local, label)
      implicit none
      integer, intent(in) :: nocc0
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(inout) :: row_local(:,:)
      character(*), intent(in) :: label
      complex(8), allocatable :: amat(:,:), row_blk(:,:)
      integer :: owner_id, blk_s, blk_e, nblk

      call assemble_owner_row_matrix_rooted(nocc0, occ_list0, row_local, amat)
      if (info%id_o == 0) then
        call build_transposed_inverse_coefficients_checked(amat, label)
      end if

      do owner_id = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_list0, owner_id, blk_s, blk_e, nblk)
        if (nblk <= 0) cycle
        allocate(row_blk(nblk, nocc0))
        row_blk(:,:) = (0.0d0, 0.0d0)
        if (info%id_o == 0) row_blk(:,:) = amat(blk_s:blk_e,:)
        call comm_bcast(row_blk, info%icomm_o, 0)
        if (info%id_o == owner_id) row_local(1:nblk,:) = row_blk(:,:)
        deallocate(row_blk)
      end do

      if (allocated(amat)) deallocate(amat)
    end subroutine build_transposed_inverse_coefficients_rowwise_checked


    subroutine build_transposed_inverse_coefficients_checked(a, label)
      implicit none
      complex(8), intent(inout) :: a(:,:)
      character(*), intent(in) :: label
      integer :: n, nrhs, lwork, info_inv
      integer, allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: a_lu(:,:), rhs(:,:)
      real(8) :: anorm, rcond
      real(8), parameter :: rcond_warn = 1.0d-10
      real(8) :: zlange
      integer :: i

      n = size(a,1)
      if (n <= 0) return

      nrhs = n
      lwork = max(1, n*max(n,64))
      allocate(ipiv(n), work(lwork), rwork(max(1,2*n)), a_lu(n,n), rhs(n,nrhs))
      a_lu(:,:) = a(:,:)
      anorm = zlange('1', n, n, a, n, rwork)

      call zgetrf(n, n, a_lu, n, ipiv, info_inv)
      if (info_inv /= 0) then
        write(*,'(a,1x,a,1x,a,i0)') 'build_transposed_inverse_coefficients_checked:', trim(label), 'zgetrf info=', info_inv
        stop 'build_transposed_inverse_coefficients_checked: LU factorization failed'
      end if

      call zgecon('1', n, a_lu, n, anorm, rcond, work, rwork, info_inv)
      if (info_inv /= 0) then
        write(*,'(a,1x,a,1x,a,i0)') 'build_transposed_inverse_coefficients_checked:', trim(label), 'zgecon info=', info_inv
        stop 'build_transposed_inverse_coefficients_checked: condition estimate failed'
      end if
      if (rcond < rcond_warn) then
        write(*,'(a,1x,a,1x,a,es12.4,1x,a,es12.4)') 'build_transposed_inverse_coefficients_checked:', trim(label), &
          'ill-conditioned matrix: rcond=', rcond, 'anorm=', anorm
      end if

      rhs(:,:) = (0.0d0, 0.0d0)
      do i = 1, n
        rhs(i,i) = (1.0d0, 0.0d0)
      end do

      call zgetrs('T', n, nrhs, a_lu, n, ipiv, rhs, n, info_inv)
      if (info_inv /= 0) then
        write(*,'(a,1x,a,1x,a,i0)') 'build_transposed_inverse_coefficients_checked:', trim(label), 'zgetrs info=', info_inv
        stop 'build_transposed_inverse_coefficients_checked: linear solve failed'
      end if

      a(:,:) = rhs(:,:)

      deallocate(rhs, a_lu, rwork, ipiv, work)
    end subroutine build_transposed_inverse_coefficients_checked


    subroutine reduce_marker_over_k(mg0, info0, lcm)
      implicit none
      type(s_rgrid), intent(in) :: mg0
      type(s_parallel_info), intent(in) :: info0
      real(8), intent(inout) :: lcm(mg0%is(1):mg0%ie(1), mg0%is(2):mg0%ie(2), mg0%is(3):mg0%ie(3))
      real(8), allocatable :: work(:,:,:)
      integer :: nloc

      nloc = size(lcm,1) * size(lcm,2) * size(lcm,3)
      allocate(work(size(lcm,1), size(lcm,2), size(lcm,3)))
      call comm_summation(lcm, work, nloc, info0%icomm_k)
      lcm(:,:,:) = work(:,:,:)
      deallocate(work)
    end subroutine reduce_marker_over_k

  end subroutine compute_local_chern_marker_from_orbital

end module rt_local_chern_marker_soi
