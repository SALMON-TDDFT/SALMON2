  subroutine stabilize_coeff_unitarity(dg_frag, itt)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: itt

    integer :: ispin, io, jo, n, nstab, n_frag, n_pw, n_tot
    complex(8), allocatable :: v(:), Sv(:), u_prev(:)
    complex(8) :: proj
    real(8) :: norm2_v, norm_v, dev_max, dev_post_max
    real(8), parameter :: eps_norm = 1.0d-14
    logical :: use_S

    dev_max = 0.0d0
    dev_post_max = 0.0d0

    do ispin = 1, dg_frag%nspin
      n_frag = dg_frag%n_mat(ispin)
      n_pw = 0
      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
      n_tot = n_frag + n_pw
      n = n_tot
      nstab = min(dg_frag%nstate_tot, n)
      if (n <= 0 .or. nstab <= 0) cycle

      use_S = allocated(dg_frag%S_mat_mixed_prop) .or. allocated(dg_frag%S_mat_prop_c) .or. &
              allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat)
      allocate(v(n), Sv(n), u_prev(n))

      do io = 1, nstab
        v(:) = (0.0d0, 0.0d0)
        v(1:n_frag) = dg_frag%coef(1:n_frag, io, ispin)
        if (n_pw > 0) v(n_frag+1:n_tot) = dg_frag%coef_pw(1:n_pw, io, ispin)

        do jo = 1, io - 1
          u_prev(:) = (0.0d0, 0.0d0)
          u_prev(1:n_frag) = dg_frag%coef(1:n_frag, jo, ispin)
          if (n_pw > 0) u_prev(n_frag+1:n_tot) = dg_frag%coef_pw(1:n_pw, jo, ispin)

          if (use_S) then
            if (n_pw > 0 .and. allocated(dg_frag%S_mat_mixed_prop)) then
              Sv(:) = matmul(dg_frag%S_mat_mixed_prop(1:n, 1:n, ispin), v(:))
            else if (allocated(dg_frag%S_mat_prop_c)) then
              Sv(:) = matmul(dg_frag%S_mat_prop_c(1:n, 1:n, ispin), v(:))
            else if (allocated(dg_frag%S_mat_prop)) then
              Sv(:) = matmul(cmplx(dg_frag%S_mat_prop(1:n, 1:n, ispin), 0.0d0, kind=8), v(:))
            else if (allocated(dg_frag%S_mat_c)) then
              Sv(:) = matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), v(:))
            else
              Sv(:) = matmul(cmplx(dg_frag%S_mat(1:n, 1:n, ispin), 0.0d0, kind=8), v(:))
            end if
            proj = sum(conjg(u_prev(:)) * Sv(:))
          else
            proj = sum(conjg(u_prev(:)) * v(:))
          end if
          v(:) = v(:) - proj * u_prev(:)
        end do

        if (use_S) then
          if (n_pw > 0 .and. allocated(dg_frag%S_mat_mixed_prop)) then
            Sv(:) = matmul(dg_frag%S_mat_mixed_prop(1:n, 1:n, ispin), v(:))
          else if (allocated(dg_frag%S_mat_prop_c)) then
            Sv(:) = matmul(dg_frag%S_mat_prop_c(1:n, 1:n, ispin), v(:))
          else if (allocated(dg_frag%S_mat_prop)) then
            Sv(:) = matmul(cmplx(dg_frag%S_mat_prop(1:n, 1:n, ispin), 0.0d0, kind=8), v(:))
          else if (allocated(dg_frag%S_mat_c)) then
            Sv(:) = matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), v(:))
          else
            Sv(:) = matmul(cmplx(dg_frag%S_mat(1:n, 1:n, ispin), 0.0d0, kind=8), v(:))
          end if
          norm2_v = real(sum(conjg(v(:)) * Sv(:)), kind=8)
        else
          norm2_v = sum(abs(v(:))**2)
        end if

        if (norm2_v < 0.0d0 .and. abs(norm2_v) < 1.0d-12) norm2_v = 0.0d0
        norm_v = sqrt(max(0.0d0, norm2_v))
        dev_max = max(dev_max, abs(norm_v - 1.0d0))

        if (norm_v < eps_norm .or. norm_v /= norm_v) then
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,es12.4,a,l1)') "[WARN] Unitary stabilization failed: spin=", &
              ispin, " state=", io, " itt=", itt, " norm=", norm_v, " use_S=", use_S
          end if
          cycle
        end if

        if (n_pw > 0) then
          dg_frag%coef(1:n_frag, io, ispin) = v(1:n_frag) / norm_v
          dg_frag%coef_pw(1:n_pw, io, ispin) = v(n_frag+1:n_tot) / norm_v
        else
          dg_frag%coef(1:n_frag, io, ispin) = v(1:n_frag) / norm_v
        end if

        ! Post-normalization deviation (actual residual after correction)
        if (use_S) then
          if (n_pw > 0 .and. allocated(dg_frag%S_mat_mixed_prop)) then
            Sv(:) = matmul(dg_frag%S_mat_mixed_prop(1:n, 1:n, ispin), v(:) / norm_v)
          else if (allocated(dg_frag%S_mat_prop_c)) then
            Sv(:) = matmul(dg_frag%S_mat_prop_c(1:n, 1:n, ispin), v(:) / norm_v)
          else if (allocated(dg_frag%S_mat_prop)) then
            Sv(:) = matmul(cmplx(dg_frag%S_mat_prop(1:n, 1:n, ispin), 0.0d0, kind=8), v(:) / norm_v)
          else if (allocated(dg_frag%S_mat_c)) then
            Sv(:) = matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), v(:) / norm_v)
          else
            Sv(:) = matmul(cmplx(dg_frag%S_mat(1:n, 1:n, ispin), 0.0d0, kind=8), v(:) / norm_v)
          end if
          norm2_v = real(sum(conjg(v(:) / norm_v) * Sv(:)), kind=8)
        else
          norm2_v = sum(abs(v(:) / norm_v)**2)
        end if
        dev_post_max = max(dev_post_max, abs(sqrt(max(0.0d0, norm2_v)) - 1.0d0))
      end do

      deallocate(v, Sv, u_prev)
    end do

    if (dg_frag%id == 0 .and. mod(itt, 10) == 0 .and. (dev_max > 1.0d-8 .or. dev_post_max > 1.0d-8)) then
      write(*,'(1x,a,i0,a,es12.4,a,es12.4)') "[INFO] Unitary stabilization at itt=", itt, &
        ", pre max(norm-1)=", dev_max, ", post max(norm-1)=", dev_post_max
    end if

  end subroutine stabilize_coeff_unitarity
