  subroutine stabilize_coeff_unitarity(dg_frag, itt)
    use rt_dg_fragment_ops, only: apply_overlap_operator, gather_full_coef_view, zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: itt

    integer :: ispin, io, jo, n, nstab, n_frag, n_pw, n_tot
    complex(8), allocatable :: v(:), Sv(:), u_prev(:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
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
      if (n_pw > 0) n_frag = dg_frag%n_mat_max
      n_tot = n_frag + n_pw
      n = n_tot
      nstab = min(dg_frag%nstate_tot, n)
      if (n <= 0 .or. nstab <= 0) cycle

      use_S = allocated(dg_frag%S_mat_mixed_prop) .or. allocated(dg_frag%S_mat_frag_pw) .or. allocated(dg_frag%S_mat_prop_c) .or. &
              allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat) .or. &
              allocated(dg_frag%S_mat_prop_blocks) .or. allocated(dg_frag%S_mat_blocks)
      allocate(v(n), Sv(n), u_prev(n))
      call gather_full_coef_view(dg_frag, ispin, n_frag, nstab, coef_frag_all, coef_pw_all)
      do io = 1, nstab
        v(:) = (0.0d0, 0.0d0)
        v(1:n_frag) = coef_frag_all(1:n_frag, io)
        if (n_pw > 0) v(n_frag+1:n_tot) = coef_pw_all(1:n_pw, io)

        do jo = 1, io - 1
          u_prev(:) = (0.0d0, 0.0d0)
          u_prev(1:n_frag) = coef_frag_all(1:n_frag, jo)
          if (n_pw > 0) u_prev(n_frag+1:n_tot) = coef_pw_all(1:n_pw, jo)

          if (use_S) then
            call apply_overlap_operator(dg_frag, ispin, v(:), Sv(:), .true.)
            proj = sum(conjg(u_prev(:)) * Sv(:))
          else
            proj = sum(conjg(u_prev(:)) * v(:))
          end if
          v(:) = v(:) - proj * u_prev(:)
        end do

        if (use_S) then
          call apply_overlap_operator(dg_frag, ispin, v(:), Sv(:), .true.)
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

        coef_frag_all(1:n_frag, io) = v(1:n_frag) / norm_v
        if (n_pw > 0) coef_pw_all(1:n_pw, io) = v(n_frag+1:n_tot) / norm_v
        do jo = 1, n_frag
          if (allocated(dg_frag%coef_owner)) then
            if (dg_frag%coef_owner(jo, ispin) /= dg_frag%id) cycle
          end if
          dg_frag%coef(jo, io, ispin) = coef_frag_all(jo, io)
        end do
        if (n_pw > 0) then
          do jo = 1, n_pw
            if (allocated(dg_frag%coef_pw_owner)) then
              if (dg_frag%coef_pw_owner(jo) /= dg_frag%id) cycle
            end if
            dg_frag%coef_pw(jo, io, ispin) = coef_pw_all(jo, io)
          end do
        end if

        ! Post-normalization deviation (actual residual after correction)
        if (use_S) then
          call apply_overlap_operator(dg_frag, ispin, v(:) / norm_v, Sv(:), .true.)
          norm2_v = real(sum(conjg(v(:) / norm_v) * Sv(:)), kind=8)
        else
          norm2_v = sum(abs(v(:) / norm_v)**2)
        end if
        dev_post_max = max(dev_post_max, abs(sqrt(max(0.0d0, norm2_v)) - 1.0d0))
      end do

      deallocate(v, Sv, u_prev, coef_frag_all, coef_pw_all)
    end do

    call zero_nonowned_coefficients(dg_frag)

    if (dg_frag%id == 0 .and. mod(itt, 10) == 0 .and. (dev_max > 1.0d-8 .or. dev_post_max > 1.0d-8)) then
      write(*,'(1x,a,i0,a,es12.4,a,es12.4)') "[INFO] Unitary stabilization at itt=", itt, &
        ", pre max(norm-1)=", dev_max, ", post max(norm-1)=", dev_post_max
    end if

  end subroutine stabilize_coeff_unitarity
