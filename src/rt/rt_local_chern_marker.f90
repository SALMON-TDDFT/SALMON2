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

    integer :: ik, ispin, io, jo, iocc, jocc, ix, iy, iz, im, p
    integer :: nocc
    integer, allocatable :: occ_idx(:)
    real(8), allocatable :: occ_w(:)
    complex(8), allocatable :: zocc(:,:,:,:), zt1(:,:,:,:), zt2(:,:,:,:)
    complex(8), allocatable :: s1(:,:), s2(:,:), s1_sum(:,:), s2_sum(:,:), g12(:,:)
    real(8), allocatable :: eval_s(:)
    complex(8), allocatable :: u_mat(:,:), x_lowdin(:,:)
    complex(8) :: phase1, phase2, zterm12, zterm21
    real(8) :: eps_occ, rvec(3)
    real(8) :: b1(3), b2(3)
    logical :: use_complex
    real(8), parameter :: eps_ortho = 1.0d-12

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

          allocate(zocc(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc))
          allocate(zt1(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc))
          allocate(zt2(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc))
          allocate(s1(nocc,nocc), s2(nocc,nocc), s1_sum(nocc,nocc), s2_sum(nocc,nocc), g12(nocc,nocc))
          allocate(eval_s(nocc), u_mat(nocc,nocc), x_lowdin(nocc,nocc))
          call copy_occupied_to_temp(ik, ispin, im, nocc, occ_idx, zocc)
          do p = 1, nocc
            zocc(:,:,:,p) = sqrt(max(0.0d0, occ_w(p))) * zocc(:,:,:,p)
          end do
          call lowdin_orthonormalize_occupied(nocc, zocc, u_mat, eval_s, x_lowdin, eps_ortho)

          s1(:,:) = (0.0d0, 0.0d0)
          s2(:,:) = (0.0d0, 0.0d0)

!$omp parallel do collapse(2) private(iocc,jocc,ix,iy,iz,rvec,phase1,phase2) schedule(static)
        do iocc = 1, nocc
          do jocc = 1, nocc
            do iz = mg%is(3), mg%ie(3)
              rvec(3) = mg%coordinate(iz,3)
              do iy = mg%is(2), mg%ie(2)
                rvec(2) = mg%coordinate(iy,2)
                do ix = mg%is(1), mg%ie(1)
                  rvec(1) = mg%coordinate(ix,1)
                  phase1 = exp(-zi * dot_product(b1, rvec))
                  phase2 = exp(-zi * dot_product(b2, rvec))
                  s1(iocc,jocc) = s1(iocc,jocc) + conjg(zocc(ix,iy,iz,iocc)) * phase1 * zocc(ix,iy,iz,jocc) * system%hvol
                  s2(iocc,jocc) = s2(iocc,jocc) + conjg(zocc(ix,iy,iz,iocc)) * phase2 * zocc(ix,iy,iz,jocc) * system%hvol
                end do
              end do
            end do
          end do
        end do

          call comm_summation(s1, s1_sum, nocc*nocc, info%icomm_r)
          call comm_summation(s2, s2_sum, nocc*nocc, info%icomm_r)

          s1(:,:) = s1_sum(:,:)
          s2(:,:) = s2_sum(:,:)
          call invert_complex_matrix_checked(s1, 'S1')
          call invert_complex_matrix_checked(s2, 'S2')

          zt1(:,:,:,:) = (0.0d0, 0.0d0)
          zt2(:,:,:,:) = (0.0d0, 0.0d0)
!$omp parallel do collapse(3) private(ix,iy,iz,iocc,jocc,rvec,phase1,phase2) schedule(static)
        do iz = mg%is(3), mg%ie(3)
          rvec(3) = mg%coordinate(iz,3)
          do iy = mg%is(2), mg%ie(2)
            rvec(2) = mg%coordinate(iy,2)
            do ix = mg%is(1), mg%ie(1)
              rvec(1) = mg%coordinate(ix,1)
              phase1 = exp(-zi * dot_product(b1, rvec))
              phase2 = exp(-zi * dot_product(b2, rvec))
              do iocc = 1, nocc
                do jocc = 1, nocc
                  zt1(ix,iy,iz,iocc) = zt1(ix,iy,iz,iocc) + s1(jocc,iocc) * phase1 * zocc(ix,iy,iz,jocc)
                  zt2(ix,iy,iz,iocc) = zt2(ix,iy,iz,iocc) + s2(jocc,iocc) * phase2 * zocc(ix,iy,iz,jocc)
                end do
              end do
            end do
          end do
        end do

          g12(:,:) = (0.0d0, 0.0d0)
!$omp parallel do collapse(2) private(iocc,jocc,ix,iy,iz) schedule(static)
        do iocc = 1, nocc
          do jocc = 1, nocc
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  g12(iocc,jocc) = g12(iocc,jocc) + conjg(zt1(ix,iy,iz,iocc)) * zt2(ix,iy,iz,jocc) * system%hvol
                end do
              end do
            end do
          end do
        end do
          call comm_summation(g12, s1_sum, nocc*nocc, info%icomm_r)
          g12(:,:) = s1_sum(:,:)

!$omp parallel do collapse(3) private(ix,iy,iz,iocc,jocc,zterm12,zterm21) schedule(static)
        do iz = mg%is(3), mg%ie(3)
          do iy = mg%is(2), mg%ie(2)
            do ix = mg%is(1), mg%ie(1)
              zterm12 = (0.0d0, 0.0d0)
              zterm21 = (0.0d0, 0.0d0)
              do iocc = 1, nocc
                do jocc = 1, nocc
                  zterm12 = zterm12 + zt1(ix,iy,iz,iocc) * g12(iocc,jocc) * conjg(zocc(ix,iy,iz,jocc))
                  zterm21 = zterm21 + zt2(ix,iy,iz,iocc) * conjg(g12(jocc,iocc)) * conjg(zocc(ix,iy,iz,jocc))
                end do
              end do

              marker(ix,iy,iz) = marker(ix,iy,iz) - aimag(zterm12 - zterm21) * system%wtk(ik) / (2.0d0*pi)
            end do
          end do
        end do

          deallocate(eval_s, u_mat, x_lowdin)
          deallocate(s1, s2, s1_sum, s2_sum, g12)
          deallocate(zt1, zt2, zocc)
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
      complex(8), intent(out) :: zbuf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc0)

      complex(8), allocatable :: zblk(:,:,:,:)
      integer :: io_g, p, m, nb, owner, io_s, io_e, nsel, t
      integer, allocatable :: psel(:), glist(:)

      zbuf(:,:,:,:) = (0.0d0, 0.0d0)

      if (.not. info%if_divide_orbit) then
        if (use_complex) then
!$omp parallel do private(p,io_g,ix,iy,iz) schedule(static)
          do p = 1, nocc0
            io_g = occ_list(p)
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
          do p = 1, nocc0
            io_g = occ_list(p)
            do iz = mg%is(3), mg%ie(3)
              do iy = mg%is(2), mg%ie(2)
                do ix = mg%is(1), mg%ie(1)
                  zbuf(ix,iy,iz,p) = cmplx(psi%rwf(ix,iy,iz,ispin0,io_g,ik0,im0), 0.0d0, kind=8)
                end do
              end do
            end do
          end do
        end if
      else
        allocate(zblk(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), max(1,nocc0)))
        allocate(psel(max(1,nocc0)), glist(max(1,nocc0)))

        do m = 0, info%isize_o - 1
          io_s = info%io_s_all(m)
          io_e = info%io_e_all(m)
          nb = info%numo_all(m)
          owner = info%irank_io(io_s)

          nsel = 0
          do p = 1, nocc0
            if (occ_list(p) >= io_s .and. occ_list(p) <= io_e) then
              nsel = nsel + 1
              psel(nsel) = p
              glist(nsel) = occ_list(p)
            end if
          end do

          if (nsel > 0) then
            zblk(:,:,:,1:nsel) = (0.0d0, 0.0d0)
            if (info%id_o == m) then
              if (use_complex) then
!$omp parallel do private(t,io_g,ix,iy,iz) schedule(static)
                do t = 1, nsel
                  io_g = glist(t)
                  do iz = mg%is(3), mg%ie(3)
                    do iy = mg%is(2), mg%ie(2)
                      do ix = mg%is(1), mg%ie(1)
                        zblk(ix,iy,iz,t) = psi%zwf(ix,iy,iz,ispin0,io_g,ik0,im0)
                      end do
                    end do
                  end do
                end do
              else
!$omp parallel do private(t,io_g,ix,iy,iz) schedule(static)
                do t = 1, nsel
                  io_g = glist(t)
                  do iz = mg%is(3), mg%ie(3)
                    do iy = mg%is(2), mg%ie(2)
                      do ix = mg%is(1), mg%ie(1)
                        zblk(ix,iy,iz,t) = cmplx(psi%rwf(ix,iy,iz,ispin0,io_g,ik0,im0), 0.0d0, kind=8)
                      end do
                    end do
                  end do
                end do
              end if
            end if
            call comm_bcast(zblk(:,:,:,1:nsel), info%icomm_o, owner)
!$omp parallel do private(t) schedule(static)
            do t = 1, nsel
              zbuf(:,:,:,psel(t)) = zblk(:,:,:,t)
            end do
          end if
        end do

        deallocate(psel, glist, zblk)
      end if
    end subroutine copy_occupied_to_temp


    subroutine lowdin_orthonormalize_occupied(nocc0, zbuf, umat, eval, xmat, eps_rel)
      implicit none
      integer, intent(in) :: nocc0
      complex(8), intent(inout) :: zbuf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc0)
      complex(8), intent(out) :: umat(nocc0,nocc0), xmat(nocc0,nocc0)
      real(8), intent(out) :: eval(nocc0)
      real(8), intent(in) :: eps_rel

      complex(8), allocatable :: s_local(:,:), s_global(:,:)
      real(8) :: smax
      integer :: ia, ib, ic

      allocate(s_local(nocc0,nocc0), s_global(nocc0,nocc0))
      s_local(:,:) = (0.0d0, 0.0d0)

!$omp parallel do collapse(2) private(ia,ib,ix,iy,iz) schedule(static)
      do ia = 1, nocc0
        do ib = 1, nocc0
          do iz = mg%is(3), mg%ie(3)
            do iy = mg%is(2), mg%ie(2)
              do ix = mg%is(1), mg%ie(1)
                s_local(ia,ib) = s_local(ia,ib) + conjg(zbuf(ix,iy,iz,ia)) * zbuf(ix,iy,iz,ib) * system%hvol
              end do
            end do
          end do
        end do
      end do

      call comm_summation(s_local, s_global, nocc0*nocc0, info%icomm_r)
      deallocate(s_local)

      umat(:,:) = s_global(:,:)
      call eigen_zheev(umat, eval, umat)

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
      complex(8), intent(inout) :: zbuf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), nocc0)
      complex(8), intent(in) :: xmat(nocc0,nocc0)

      complex(8), allocatable :: vin(:), vout(:)
      integer :: ia, ib

      allocate(vin(nocc0), vout(nocc0))
      do iz = mg%is(3), mg%ie(3)
        do iy = mg%is(2), mg%ie(2)
          do ix = mg%is(1), mg%ie(1)
            do ia = 1, nocc0
              vin(ia) = zbuf(ix,iy,iz,ia)
            end do
            do ib = 1, nocc0
              vout(ib) = (0.0d0, 0.0d0)
              do ia = 1, nocc0
                vout(ib) = vout(ib) + vin(ia) * xmat(ia,ib)
              end do
            end do
            do ib = 1, nocc0
              zbuf(ix,iy,iz,ib) = vout(ib)
            end do
          end do
        end do
      end do
      deallocate(vin, vout)
    end subroutine apply_right_transform_occ_inplace


    subroutine invert_complex_matrix_checked(a, label)
      implicit none
      complex(8), intent(inout) :: a(:,:)
      character(*), intent(in) :: label
      integer :: n, lwork, info_inv
      integer, allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)

      n = size(a,1)
      if (n <= 0) return

      lwork = max(1, n*max(n,64))
      allocate(ipiv(n), work(lwork))

      call zgetrf(n, n, a, n, ipiv, info_inv)
      if (info_inv /= 0) then
        write(*,'(a,1x,a,1x,a,i0)') 'invert_complex_matrix_checked:', trim(label), 'zgetrf info=', info_inv
        stop 'invert_complex_matrix_checked: LU factorization failed'
      end if

      call zgetri(n, a, n, ipiv, work, lwork, info_inv)
      if (info_inv /= 0) then
        write(*,'(a,1x,a,1x,a,i0)') 'invert_complex_matrix_checked:', trim(label), 'zgetri info=', info_inv
        stop 'invert_complex_matrix_checked: matrix inversion failed'
      end if

      deallocate(ipiv, work)
    end subroutine invert_complex_matrix_checked


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
