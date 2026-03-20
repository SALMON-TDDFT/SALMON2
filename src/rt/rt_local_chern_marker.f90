#include "config.h"

module rt_local_chern_marker
  use math_constants, only: pi, zi
  use structures, only: s_rgrid, s_orbital, s_dft_system, s_parallel_info, &
                        allocate_orbital_real, allocate_orbital_complex, deallocate_orbital
  use communication, only: comm_summation, comm_bcast
  use checkpoint_restart_sub, only: read_wavefunction
  use eigen_subdiag_sub, only: eigen_zheev
  implicit none
  private

  public :: compute_local_chern_marker_from_orbital
  public :: compute_local_chern_marker_from_rt_checkpoint

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

    integer :: ik, ispin, io, iocc, jocc, ix, iy, iz, im, p, owner, ig
    integer :: nocc, nocc_local, iocc_g, jocc_g, nblk, nloc
    integer, allocatable :: occ_idx(:)
    real(8), allocatable :: occ_w(:)
    complex(8), allocatable :: zocc(:,:,:,:), zt1(:,:,:,:), zt2(:,:,:,:), zblk(:,:,:,:), g12_blk(:,:)
    complex(8), allocatable :: zlhs(:,:), zrhs1(:,:), zrhs2(:,:), ztmp1(:,:), ztmp2(:,:)
    complex(8), allocatable :: s1_row(:,:), s2_row(:,:), s1_row_sum(:,:), s2_row_sum(:,:)
    complex(8), allocatable :: g12_row(:,:), g12_row_sum(:,:)
    real(8), allocatable :: eval_s(:)
    complex(8), allocatable :: u_mat(:,:), x_lowdin(:,:)
    complex(8) :: phase1, phase2, zterm12, zterm21
    real(8) :: eps_occ, rvec(3)
    real(8) :: t_s12_copy, t_s12_pack, t_s12_gemm, t_s12_accum, t0, t1
    real(8) :: b1(3), b2(3)
    logical :: use_complex
    real(8), parameter :: eps_ortho = 1.0d-12
    real(8), allocatable :: marker_local(:,:,:), marker_sum(:,:,:)

    use_complex = allocated(psi%zwf)
    eps_occ = 1.0d-10
    if (present(occ_eps)) eps_occ = occ_eps

    b1(:) = system%primitive_b(:,1)
    b2(:) = system%primitive_b(:,2)

    marker(:,:,:) = 0.0d0

    do ik = info%ik_s, info%ik_e
      do ispin = 1, system%nspin
        do im = info%im_s, info%im_e

          nocc = 0
          do io = 1, system%no
            if (system%rocc(io, ik, ispin) > eps_occ) nocc = nocc + 1
          end do
          if (nocc <= 0) cycle

          allocate(occ_idx(nocc), occ_w(nocc))
          nocc = 0
          do io = 1, system%no
            if (system%rocc(io, ik, ispin) > eps_occ) then
              nocc = nocc + 1
              occ_idx(nocc) = io
              occ_w(nocc) = system%rocc(io, ik, ispin)
            end if
          end do

          call get_local_occ_range(nocc, occ_idx, occ_w, info%id_o, nocc_local)
          allocate(zocc(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,nocc_local)))
          allocate(zt1(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,nocc_local)))
          allocate(zt2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,nocc_local)))
          allocate(s1_row(max(1,nocc_local),nocc), s2_row(max(1,nocc_local),nocc))
          allocate(s1_row_sum(max(1,nocc_local),nocc), s2_row_sum(max(1,nocc_local),nocc))
          allocate(g12_row(max(1,nocc_local),nocc), g12_row_sum(max(1,nocc_local),nocc))
          allocate(eval_s(nocc), u_mat(nocc,nocc), x_lowdin(nocc,nocc))
          allocate(marker_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          allocate(marker_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          marker_local(:,:,:) = 0d0
          marker_sum(:,:,:) = 0d0
          if (info%id_rko == 0) write(*,'(a,4(i0,1x))') 'LCM stage: setup ik ispin im nocc = ', ik, ispin, im, nocc
          call copy_occupied_to_temp(ik, ispin, im, nocc, occ_idx, zocc)
          if (info%id_rko == 0) write(*,*) 'LCM stage: copy_occupied_to_temp done'
          do p = 1, nocc_local
            zocc(:,:,:,p) = sqrt(max(0.0d0, local_occ_weight(nocc, occ_idx, occ_w, info%id_o, p))) * zocc(:,:,:,p)
          end do
          call lowdin_orthonormalize_occupied(nocc, zocc, u_mat, eval_s, x_lowdin, eps_ortho)
          if (info%id_rko == 0) write(*,*) 'LCM stage: lowdin done'

          s1_row(:,:) = (0.0d0, 0.0d0)
          s2_row(:,:) = (0.0d0, 0.0d0)
          t_s12_copy = 0d0
          t_s12_pack = 0d0
          t_s12_gemm = 0d0
          t_s12_accum = 0d0
          allocate(zlhs(max(1,nloc), max(1,nocc_local)))
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            call cpu_time(t0)
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
            call cpu_time(t1)
            t_s12_copy = t_s12_copy + (t1 - t0)
            allocate(zrhs1(max(1,nloc), nblk), zrhs2(max(1,nloc), nblk))
            allocate(ztmp1(max(1,nocc_local), nblk), ztmp2(max(1,nocc_local), nblk))
            call cpu_time(t0)
            ig = 0
            do iz = mg%is(3), mg%ie(3)
              rvec(3) = mg%coordinate(iz,3)
              do iy = mg%is(2), mg%ie(2)
                rvec(2) = mg%coordinate(iy,2)
                do ix = mg%is(1), mg%ie(1)
                  rvec(1) = mg%coordinate(ix,1)
                  phase1 = exp(-zi * dot_product(b1, rvec))
                  phase2 = exp(-zi * dot_product(b2, rvec))
                  ig = ig + 1
                  do p = 1, nocc_local
                    zlhs(ig,p) = zocc(ix,iy,iz,p)
                  end do
                  do jocc = 1, nblk
                    zrhs1(ig,jocc) = phase1 * zblk(ix,iy,iz,jocc)
                    zrhs2(ig,jocc) = phase2 * zblk(ix,iy,iz,jocc)
                  end do
                end do
              end do
            end do
            call cpu_time(t1)
            t_s12_pack = t_s12_pack + (t1 - t0)
            call cpu_time(t0)
            call zgemm('C', 'N', nocc_local, nblk, nloc, cmplx(system%hvol, 0.0d0, kind=8), zlhs, max(1,nloc), &
              zrhs1, max(1,nloc), (1.0d0, 0.0d0), ztmp1, max(1,nocc_local))
            call zgemm('C', 'N', nocc_local, nblk, nloc, cmplx(system%hvol, 0.0d0, kind=8), zlhs, max(1,nloc), &
              zrhs2, max(1,nloc), (1.0d0, 0.0d0), ztmp2, max(1,nocc_local))
            call cpu_time(t1)
            t_s12_gemm = t_s12_gemm + (t1 - t0)
            call cpu_time(t0)
            do p = 1, nocc_local
              do jocc = 1, nblk
                jocc_g = iocc + jocc - 1
                s1_row(p,jocc_g) = s1_row(p,jocc_g) + ztmp1(p,jocc)
                s2_row(p,jocc_g) = s2_row(p,jocc_g) + ztmp2(p,jocc)
              end do
            end do
            call cpu_time(t1)
            t_s12_accum = t_s12_accum + (t1 - t0)
            if (info%id_rko == 0) then
              if (mod(owner + 1, 8) == 0 .or. owner == info%isize_o - 1) then
                write(*,'(a,1x,i0,1x,a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
                  'LCM S1/S2 progress: owner', owner, 'copy=', t_s12_copy, 'pack=', t_s12_pack, 'gemm=', t_s12_gemm, 'accum=', t_s12_accum
              end if
            end if
            deallocate(ztmp1, ztmp2, zrhs1, zrhs2)
            deallocate(zblk)
          end do
          deallocate(zlhs)

          call comm_summation(s1_row, s1_row_sum, size(s1_row), info%icomm_r)
          call comm_summation(s2_row, s2_row_sum, size(s2_row), info%icomm_r)
          if (info%id_rko == 0) then
            write(*,'(a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
              'LCM S1/S2 timing: copy=', t_s12_copy, 'pack=', t_s12_pack, 'gemm=', t_s12_gemm, 'accum=', t_s12_accum
          end if
          if (info%id_rko == 0) write(*,*) 'LCM stage: S1/S2 done'
          call build_transposed_inverse_coefficients_rowwise_checked(nocc, occ_idx, s1_row_sum, 'S1')
          call build_transposed_inverse_coefficients_rowwise_checked(nocc, occ_idx, s2_row_sum, 'S2')

          zt1(:,:,:,:) = (0.0d0, 0.0d0)
          zt2(:,:,:,:) = (0.0d0, 0.0d0)
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
!$omp parallel do private(ix,iy,iz,p,jocc_g,rvec,phase1,phase2,jocc) schedule(static)
            do iz = mg%is(3), mg%ie(3)
              rvec(3) = mg%coordinate(iz,3)
              do iy = mg%is(2), mg%ie(2)
                rvec(2) = mg%coordinate(iy,2)
                do ix = mg%is(1), mg%ie(1)
                  rvec(1) = mg%coordinate(ix,1)
                  phase1 = exp(-zi * dot_product(b1, rvec))
                  phase2 = exp(-zi * dot_product(b2, rvec))
                  do p = 1, nocc_local
                    do jocc = 1, nblk
                      jocc_g = iocc + jocc - 1
                      zt1(ix,iy,iz,p) = zt1(ix,iy,iz,p) + s1_row_sum(p,jocc_g) * phase1 * zblk(ix,iy,iz,jocc)
                      zt2(ix,iy,iz,p) = zt2(ix,iy,iz,p) + s2_row_sum(p,jocc_g) * phase2 * zblk(ix,iy,iz,jocc)
                    end do
                  end do
                end do
              end do
            end do
            deallocate(zblk)
          end do

          g12_row(:,:) = (0.0d0, 0.0d0)
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zt2, iocc, jocc, zblk)
            do p = 1, nocc_local
              iocc_g = local_occ_index(nocc, occ_idx, info%id_o, p)
              do jocc = 1, nblk
                jocc_g = iocc + jocc - 1
                do iz = mg%is(3), mg%ie(3)
                  do iy = mg%is(2), mg%ie(2)
                    do ix = mg%is(1), mg%ie(1)
                      g12_row(p,jocc_g) = g12_row(p,jocc_g) + conjg(zt1(ix,iy,iz,p)) * zblk(ix,iy,iz,jocc) * system%hvol
                    end do
                  end do
                end do
              end do
            end do
            deallocate(zblk)
          end do
          call comm_summation(g12_row, g12_row_sum, size(g12_row), info%icomm_r)
          if (info%id_rko == 0) write(*,*) 'LCM stage: G12 done'

          marker_local(:,:,:) = 0d0
          do owner = 0, info%isize_o - 1
            call get_owner_occ_block(nocc, occ_idx, owner, iocc, jocc, nblk)
            if (nblk <= 0) cycle
            allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
            call get_owned_occ_block(nocc, occ_idx, owner, zocc, iocc, jocc, zblk)
            allocate(g12_blk(max(1,nblk), nocc))
            call get_owned_row_block(nocc, occ_idx, owner, g12_row_sum, iocc, jocc, g12_blk)
!$omp parallel do collapse(3) private(ix,iy,iz,p,jocc,jocc_g,zterm12,zterm21) schedule(static)
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  zterm12 = (0.0d0, 0.0d0)
                  zterm21 = (0.0d0, 0.0d0)
                  do p = 1, nocc_local
                    do jocc = 1, nblk
                      jocc_g = iocc + jocc - 1
                      zterm12 = zterm12 + zt1(ix,iy,iz,p) * g12_row_sum(p,jocc_g) * conjg(zblk(ix,iy,iz,jocc))
                      zterm21 = zterm21 + zt2(ix,iy,iz,p) * conjg(g12_blk(jocc,local_occ_index(nocc, occ_idx, info%id_o, p))) * &
                        conjg(zblk(ix,iy,iz,jocc))
                    end do
                  end do
                  marker_local(ix,iy,iz) = marker_local(ix,iy,iz) - aimag(zterm12 - zterm21) * system%wtk(ik) / (2.0d0*pi)
                end do
              end do
            end do
            deallocate(g12_blk, zblk)
          end do
          nloc = size(marker_local,1) * size(marker_local,2) * size(marker_local,3)
          call comm_summation(marker_local, marker_sum, nloc, info%icomm_o)
          if (info%id_rko == 0) write(*,*) 'LCM stage: marker reduce done'
          marker(:,:,:) = marker(:,:,:) + marker_sum(:,:,:)

          deallocate(eval_s, u_mat, x_lowdin)
          deallocate(s1_row, s2_row, s1_row_sum, s2_row_sum, g12_row, g12_row_sum)
          deallocate(zt1, zt2, zocc)
          deallocate(marker_local, marker_sum)
          deallocate(occ_idx, occ_w)

        end do
      end do
    end do

    if (info%isize_k > 1) then
      call reduce_marker_over_k(mg, info, marker)
    end if

  contains

    subroutine copy_occupied_to_temp(ik0, ispin0, im0, nocc0, occ_list, zbuf)
      implicit none
      integer, intent(in) :: ik0, ispin0, im0, nocc0
      integer, intent(in) :: occ_list(nocc0)
      complex(8), intent(out), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :) :: zbuf

      integer :: io_g, p, nocc_local0

      zbuf(:,:,:,:) = (0.0d0, 0.0d0)
      nocc_local0 = size(zbuf,4)

      if (nocc_local0 <= 0) return

      if (use_complex) then
!$omp parallel do private(p,io_g,ix,iy,iz) schedule(static)
        do p = 1, nocc_local0
          io_g = local_occ_global_io(nocc0, occ_list, info%id_o, p)
          do iz = mg%is(3), mg%ie(3)
            do iy = mg%is(2), mg%ie(2)
              do ix = mg%is(1), mg%ie(1)
                zbuf(ix,iy,iz,p) = psi%zwf(ix,iy,iz,ispin0,io_g,ik0,im0)
              end do
            end do
          end do
        end do
      else
!$omp parallel do private(p,io_g,ix,iy,iz) schedule(static)
        do p = 1, nocc_local0
          io_g = local_occ_global_io(nocc0, occ_list, info%id_o, p)
          do iz = mg%is(3), mg%ie(3)
            do iy = mg%is(2), mg%ie(2)
              do ix = mg%is(1), mg%ie(1)
                zbuf(ix,iy,iz,p) = cmplx(psi%rwf(ix,iy,iz,ispin0,io_g,ik0,im0), 0.0d0, kind=8)
              end do
            end do
          end do
        end do
      end if
    end subroutine copy_occupied_to_temp


    subroutine lowdin_orthonormalize_occupied(nocc0, zbuf, umat, eval, xmat, eps_rel)
      implicit none
      integer, intent(in) :: nocc0
      complex(8), intent(inout), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :) :: zbuf
      complex(8), intent(out) :: umat(nocc0,nocc0), xmat(nocc0,nocc0)
      real(8), intent(out) :: eval(nocc0)
      real(8), intent(in) :: eps_rel

      complex(8), allocatable :: s_local(:,:), s_global(:,:)
      real(8) :: smax
      integer :: ia, ib, ic, owner, ib_s, ib_e, nblk, ia_g
      complex(8), allocatable :: zblk2(:,:,:,:)

      allocate(s_local(nocc0,nocc0), s_global(nocc0,nocc0))
      s_local(:,:) = (0.0d0, 0.0d0)

      do owner = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_idx, owner, ib_s, ib_e, nblk)
        if (nblk <= 0) cycle
        allocate(zblk2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
        call get_owned_occ_block(nocc0, occ_idx, owner, zbuf, ib_s, ib_e, zblk2)
        do ia = 1, size(zbuf,4)
          ia_g = local_occ_index(nocc0, occ_idx, info%id_o, ia)
          do ib = 1, nblk
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  s_local(ia_g,ib_s + ib - 1) = s_local(ia_g,ib_s + ib - 1) + conjg(zbuf(ix,iy,iz,ia)) * zblk2(ix,iy,iz,ib) * system%hvol
                end do
              end do
            end do
          end do
        end do
        deallocate(zblk2)
      end do
      if (info%id_rko == 0) write(*,*) 'LCM lowdin stage: local overlap built'

      call comm_summation(s_local, s_global, nocc0*nocc0, info%icomm_r)
      if (info%id_rko == 0) write(*,*) 'LCM lowdin stage: icomm_r reduction done'
      call comm_summation(s_global, umat, nocc0*nocc0, info%icomm_o)
      if (info%id_rko == 0) write(*,*) 'LCM lowdin stage: icomm_o reduction done'
      deallocate(s_local)
      call eigen_zheev(umat, eval, umat)
      if (info%id_rko == 0) write(*,*) 'LCM lowdin stage: eigen_zheev done'

      smax = maxval(abs(eval))
      if (smax <= 0.0d0) stop 'lowdin_orthonormalize_occupied: invalid overlap spectrum'

      xmat(:,:) = (0.0d0, 0.0d0)
!$omp parallel do collapse(2) private(ia,ib,ic) schedule(static)
      do ia = 1, nocc0
        do ib = 1, nocc0
          do ic = 1, nocc0
            if (eval(ic) > eps_rel * smax) then
              xmat(ia,ib) = xmat(ia,ib) + umat(ia,ic) * conjg(umat(ib,ic)) / sqrt(eval(ic))
            end if
          end do
        end do
      end do

      call apply_right_transform_occ_inplace(nocc0, zbuf, xmat)
      deallocate(s_global)
    end subroutine lowdin_orthonormalize_occupied


    subroutine apply_right_transform_occ_inplace(nocc0, zbuf, xmat)
      implicit none
      integer, intent(in) :: nocc0
      complex(8), intent(inout), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :) :: zbuf
      complex(8), intent(in) :: xmat(nocc0,nocc0)

      complex(8), allocatable :: vin_blk(:), vout(:), zblk2(:,:,:,:), zsrc(:,:,:,:)
      integer :: ia, ib, owner, ia_s, ia_e, nblk, ib_g

      allocate(vout(max(1,size(zbuf,4))))
      allocate(zsrc(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,size(zbuf,4))))
      zsrc(:,:,:,:) = zbuf(:,:,:,:)
      zbuf(:,:,:,:) = (0.0d0, 0.0d0)
      do owner = 0, info%isize_o - 1
        call get_owner_occ_block(nocc0, occ_idx, owner, ia_s, ia_e, nblk)
        if (nblk <= 0) cycle
        allocate(zblk2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nblk))
        call get_owned_occ_block(nocc0, occ_idx, owner, zsrc, ia_s, ia_e, zblk2)
        allocate(vin_blk(nblk))
      do iz = mg%is(3), mg%ie(3)
        do iy = mg%is(2), mg%ie(2)
          do ix = mg%is(1), mg%ie(1)
            do ia = 1, nblk
              vin_blk(ia) = zblk2(ix,iy,iz,ia)
            end do
            do ib = 1, size(zbuf,4)
              ib_g = local_occ_index(nocc0, occ_idx, info%id_o, ib)
              vout(ib) = (0.0d0, 0.0d0)
              do ia = 1, nblk
                vout(ib) = vout(ib) + vin_blk(ia) * xmat(ia_s + ia - 1, ib_g)
              end do
            end do
            do ib = 1, size(zbuf,4)
              zbuf(ix,iy,iz,ib) = zbuf(ix,iy,iz,ib) + vout(ib)
            end do
          end do
        end do
      end do
        deallocate(vin_blk, zblk2)
      end do
      deallocate(zsrc)
      deallocate(vout)
    end subroutine apply_right_transform_occ_inplace

    subroutine get_local_occ_range(nocc0, occ_list0, occ_w0, owner_id, nlocal)
      implicit none
      integer, intent(in) :: nocc0, owner_id
      integer, intent(in) :: occ_list0(nocc0)
      real(8), intent(in) :: occ_w0(nocc0)
      integer, intent(out) :: nlocal
      integer :: p
      nlocal = 0
      do p = 1, nocc0
        if (occ_owner_id(occ_list0(p)) == owner_id) nlocal = nlocal + 1
      end do
    end subroutine get_local_occ_range

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
      complex(8), intent(in), dimension(mg%is(1):, mg%is(2):, mg%is(3):, :) :: zloc
      complex(8), intent(out) :: zblk_out(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,blk_e-blk_s+1))
      integer :: nblk0, p, locp
      nblk0 = max(0, blk_e - blk_s + 1)
      zblk_out = (0.0d0, 0.0d0)
      if (nblk0 <= 0) return
      if (info%id_o == owner_id) then
        locp = 0
        do p = blk_s, blk_e
          locp = locp + 1
          zblk_out(:,:,:,locp) = zloc(:,:,:,local_occ_position(nocc0, occ_list0, occ_list0(p), owner_id))
        end do
      end if
      call comm_bcast(zblk_out(:,:,:,1:nblk0), info%icomm_o, owner_id)
    end subroutine get_owned_occ_block

    subroutine assemble_owner_row_matrix(nocc0, occ_list0, row_local, amat)
      implicit none
      integer, intent(in) :: nocc0
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(in) :: row_local(:,:)
      complex(8), intent(out) :: amat(nocc0,nocc0)
      complex(8), allocatable :: row_blk(:,:)
      integer :: owner_id, blk_s, blk_e, nblk, p

      amat(:,:) = (0.0d0, 0.0d0)
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
        amat(blk_s:blk_e,:) = row_blk(1:nblk,:)
        deallocate(row_blk)
      end do
    end subroutine assemble_owner_row_matrix

    subroutine get_owned_row_block(nocc0, occ_list0, owner_id, row_local, blk_s, blk_e, row_blk_out)
      implicit none
      integer, intent(in) :: nocc0, owner_id, blk_s, blk_e
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(in) :: row_local(:,:)
      complex(8), intent(out) :: row_blk_out(max(1,blk_e-blk_s+1), nocc0)
      integer :: nblk0, p

      nblk0 = max(0, blk_e - blk_s + 1)
      row_blk_out(:,:) = (0.0d0, 0.0d0)
      if (nblk0 <= 0) return
      if (info%id_o == owner_id) then
        do p = 1, nblk0
          row_blk_out(p,:) = row_local(p,:)
        end do
      end if
      call comm_bcast(row_blk_out(1:nblk0,:), info%icomm_o, owner_id)
    end subroutine get_owned_row_block

    subroutine build_transposed_inverse_coefficients_rowwise_checked(nocc0, occ_list0, row_local, label)
      implicit none
      integer, intent(in) :: nocc0
      integer, intent(in) :: occ_list0(nocc0)
      complex(8), intent(inout) :: row_local(:,:)
      character(*), intent(in) :: label
      complex(8), allocatable :: amat(:,:), row_blk(:,:)
      integer :: owner_id, blk_s, blk_e, nblk

      if (info%id_o == 0) then
        allocate(amat(nocc0,nocc0))
        call assemble_owner_row_matrix(nocc0, occ_list0, row_local, amat)
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

end module rt_local_chern_marker
