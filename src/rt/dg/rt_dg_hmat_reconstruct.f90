  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, &
      Vh_buffer, Vxc_buffer, Vpsl_buffer)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use structures
    use communication, only: comm_summation
    use salmon_global, only: theory
    use salmon_global, only: yn_hse
    use rt_dg_fragment_ops, only: rebuild_local_h_block_ids, zero_nonlocal_h_matrix_blocks
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_stencil),        intent(in)    :: stencil
    type(s_scalar),         intent(in)    :: Vh
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    type(s_scalar),         intent(in)    :: Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    real(8),      optional, intent(in)    :: Vh_buffer(:,:,:)
    real(8),      optional, intent(in)    :: Vxc_buffer(:,:,:,:)
    real(8),      optional, intent(in)    :: Vpsl_buffer(:,:,:)

    type(s_rgrid), pointer :: mg
    integer :: ifrag, ispin, io, jo, i_local
    integer :: nbf, nbf_raw, iblk
    integer :: max_nbf_local
    integer :: frag_rank, frag_size
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    real(8) :: hvol, A2val, hij_avg
    real(8) :: t0, t1
    real(8) :: time_local_build, time_subgroup_reduce, time_global_reduce
    real(8) :: time_halo_exchange, time_potential_build, time_post_reduce_cleanup, time_block_hermitize
    complex(8) :: integral_v
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: V_phi(:,:,:)
    real(8), allocatable :: partial_total(:), partial_block(:,:), reduced_block(:,:)
    integer, allocatable :: map_x(:), map_y(:), map_z(:)
    logical :: has_overlap
    logical :: use_block_reconstruct
    logical, save :: debug_static_seed_logged = .false.
    logical, parameter :: enable_reconstruct_timing = .false.
    logical, parameter :: enable_hmat_nan_check = .false.
    real(8), parameter :: hmat_abs_limit = 1.0d12

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. associated(dg_frag%mg)) then
      stop "reconstruct_hamiltonian_matrix requires dg_frag%mg"
    end if
    mg => dg_frag%mg
    time_local_build = 0.0d0
    time_subgroup_reduce = 0.0d0
    time_global_reduce = 0.0d0
    time_halo_exchange = 0.0d0
    time_potential_build = 0.0d0
    time_post_reduce_cleanup = 0.0d0
    time_block_hermitize = 0.0d0

    hvol = system%hvol
    if (hvol /= hvol) then
      write(*,'(a,i0)') "[NaN] hvol in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
      stop "NaN in hvol"
    end if
    if (enable_hmat_nan_check) then
      if (any(Vpsl%f(lbound(Vpsl%f,1):ubound(Vpsl%f,1), lbound(Vpsl%f,2):ubound(Vpsl%f,2), lbound(Vpsl%f,3):ubound(Vpsl%f,3)) /= &
              Vpsl%f(lbound(Vpsl%f,1):ubound(Vpsl%f,1), lbound(Vpsl%f,2):ubound(Vpsl%f,2), lbound(Vpsl%f,3):ubound(Vpsl%f,3)))) then
        write(*,'(a,i0)') "[NaN] Vpsl in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
        stop "NaN in Vpsl"
      end if
      if (any(Vh%f(lbound(Vh%f,1):ubound(Vh%f,1), lbound(Vh%f,2):ubound(Vh%f,2), lbound(Vh%f,3):ubound(Vh%f,3)) /= &
              Vh%f(lbound(Vh%f,1):ubound(Vh%f,1), lbound(Vh%f,2):ubound(Vh%f,2), lbound(Vh%f,3):ubound(Vh%f,3)))) then
        write(*,'(a,i0)') "[NaN] Vh in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
        stop "NaN in Vh"
      end if
      do ispin = 1, system%nspin
        if (any(Vxc(ispin)%f(lbound(Vxc(ispin)%f,1):ubound(Vxc(ispin)%f,1), lbound(Vxc(ispin)%f,2):ubound(Vxc(ispin)%f,2), &
                             lbound(Vxc(ispin)%f,3):ubound(Vxc(ispin)%f,3)) /= &
                Vxc(ispin)%f(lbound(Vxc(ispin)%f,1):ubound(Vxc(ispin)%f,1), lbound(Vxc(ispin)%f,2):ubound(Vxc(ispin)%f,2), &
                             lbound(Vxc(ispin)%f,3):ubound(Vxc(ispin)%f,3)))) then
          write(*,'(a,i0,a,i0)') "[NaN] Vxc in reconstruct_hamiltonian_matrix, rank=", dg_frag%id, " ispin=", ispin
          stop "NaN in Vxc"
        end if
      end do
    end if

    use_block_reconstruct = allocated(dg_frag%H_mat_blocks) .and. allocated(dg_frag%H_mat_kinetic_blocks) .and. &
      allocated(dg_frag%H_block_map)
    if (.not. use_block_reconstruct) then
      stop "reconstruct_hamiltonian_matrix requires block Hamiltonian storage"
    end if
    if (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) then
      stop "reconstruct_hamiltonian_matrix block-only path does not support complex dense Hamiltonian"
    end if
    if (yn_hse == 'y') then
      stop "reconstruct_hamiltonian_matrix block-only path does not support HSE exact exchange"
    end if

    do iblk = 1, size(dg_frag%H_mat_blocks)
      dg_frag%H_mat_blocks(iblk)%val(:, :, :) = 0.0d0
      if (dg_frag%is_frag_root) then
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :)
      end if
    end do

    allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))

    max_nbf_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      nbf_raw = dg_frag%n_basis(ifrag, 1)
      if (system%nspin > 1) nbf_raw = max(nbf_raw, dg_frag%n_basis(ifrag, 2))
      nbf = min(nbf_raw, dg_frag%nstate_frag, size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
      max_nbf_local = max(max_nbf_local, nbf)
    end do
    if (max_nbf_local > 0) then
      allocate(partial_total(max_nbf_local), &
               partial_block(max_nbf_local, max_nbf_local), &
               reduced_block(max_nbf_local, max_nbf_local))
    end if

    frag_size = max(1, dg_frag%isize_frag)
    frag_rank = modulo(dg_frag%id_frag, frag_size)

    do ispin = 1, system%nspin
      call cpu_time(t0)
      if (present(Vh_buffer) .and. present(Vxc_buffer) .and. present(Vpsl_buffer)) then
        call assert_real_hmat_reconstruct_bounded_3d("Vh_buffer", Vh_buffer, &
          hmat_abs_limit, dg_frag%id, 0, ispin, 0)
        call assert_real_hmat_reconstruct_bounded_3d("Vpsl_buffer", Vpsl_buffer, &
          hmat_abs_limit, dg_frag%id, 0, ispin, 0)
        call assert_real_hmat_reconstruct_bounded_3d("Vxc_buffer", Vxc_buffer(:, :, :, ispin), &
          hmat_abs_limit, dg_frag%id, 0, ispin, 0)
        call build_total_potential_grid_with_stage_buffers(mg, dg_frag%lgnum_total, &
          Vh_buffer, Vxc_buffer, Vpsl_buffer, ispin, V_total)
      else if (present(Vh_buffer)) then
        call assert_real_hmat_reconstruct_bounded_3d("Vh_buffer", Vh_buffer, &
          hmat_abs_limit, dg_frag%id, 0, ispin, 0)
        call build_total_potential_grid_with_buffered_hartree(mg, dg_frag, Vh_buffer, Vxc(ispin), Vpsl, V_total)
      else
        call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
      end if
      if (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v)) then
!$omp parallel do collapse(3) private(A2val) schedule(static)
        do i_local = mg%is(3), mg%ie(3)
          do jo = mg%is(2), mg%ie(2)
            do io = mg%is(1), mg%ie(1)
              A2val = system%Ac_micro%v(1, io, jo, i_local)**2 + &
                      system%Ac_micro%v(2, io, jo, i_local)**2 + &
                      system%Ac_micro%v(3, io, jo, i_local)**2
              V_total(io, jo, i_local) = V_total(io, jo, i_local) + 0.5d0 * A2val
            end do
          end do
        end do
!$omp end parallel do
      end if
      call cpu_time(t1)
      time_potential_build = time_potential_build + (t1 - t0)
      call assert_real_hmat_reconstruct_bounded_3d("V_total", V_total, &
        hmat_abs_limit, dg_frag%id, 0, ispin, 0)

      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
        if (nbf_raw < 0) then
          write(*,'(a,i0,a,i0,a,i0,a,i0)') "[FATAL] reconstruct negative n_basis, rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf_raw
          stop "negative n_basis in reconstruct_hamiltonian_matrix"
        end if

        nbf = min(nbf_raw, dg_frag%nstate_frag, size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
        if (nbf <= 0) cycle

        iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if (iblk <= 0) cycle

        partial_block(1:nbf, 1:nbf) = 0.0d0
        call assert_real_hmat_reconstruct_bounded_4d("phi_frag", &
          dg_frag%phi_frag(:, :, :, 1:nbf, i_local), hmat_abs_limit, &
          dg_frag%id, ifrag, ispin, 0)

        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        g_s(:) = iorg(:)
        g_e(:) = iorg(:) + ndom(:) - 1
        ov_s(:) = max(g_s(:), mg%is(:))
        ov_e(:) = min(g_e(:), mg%ie(:))
        has_overlap = all(ov_s(:) <= ov_e(:))
        if (has_overlap) then
          loc_s(:) = ov_s(:) - iorg(:) + 1
          loc_e(:) = ov_e(:) - iorg(:) + 1
          lx_lo = loc_s(1)
          lx_hi = loc_e(1)
          ly_lo = loc_s(2)
          ly_hi = loc_e(2)
          lz_lo = loc_s(3)
          lz_hi = loc_e(3)
          p_lb1 = lbound(dg_frag%phi_frag, 1)
          p_ub1 = ubound(dg_frag%phi_frag, 1)
          p_lb2 = lbound(dg_frag%phi_frag, 2)
          p_ub2 = ubound(dg_frag%phi_frag, 2)
          p_lb3 = lbound(dg_frag%phi_frag, 3)
          p_ub3 = ubound(dg_frag%phi_frag, 3)
          allocate(map_x(ov_s(1):ov_e(1)), map_y(ov_s(2):ov_e(2)), map_z(ov_s(3):ov_e(3)))
          do io = ov_s(1), ov_e(1)
            map_x(io) = map_global_to_phi_box_coord_reconstruct(io, p_lb1, p_ub1, dg_frag%lgnum_total(1))
          end do
          do io = ov_s(2), ov_e(2)
            map_y(io) = map_global_to_phi_box_coord_reconstruct(io, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          end do
          do io = ov_s(3), ov_e(3)
            map_z(io) = map_global_to_phi_box_coord_reconstruct(io, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          end do
        end if

        do jo = frag_rank + 1, nbf, frag_size
          if (.not. has_overlap) cycle

          call cpu_time(t0)
          call build_local_potential_applied_basis_mapped(dg_frag, i_local, jo, mg, V_total, V_phi, &
            lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi, ov_s, ov_e, map_x, map_y, map_z)
          call assert_complex_hmat_reconstruct_bounded_3d("V_phi", V_phi, &
            hmat_abs_limit, dg_frag%id, ifrag, ispin, jo)

          partial_total(1:nbf) = 0.0d0

!$omp parallel do private(io, integral_v) schedule(static)
          do io = 1, nbf
            call integrate_local_basis_with_field_mapped(dg_frag, i_local, io, mg, V_phi, hvol, integral_v, &
              lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi, ov_s, ov_e, map_x, map_y, map_z)
            partial_total(io) = real(integral_v, kind=8)
          end do
!$omp end parallel do
          call cpu_time(t1)
          time_local_build = time_local_build + (t1 - t0)
          call assert_real_hmat_reconstruct_bounded_1d("partial_total", &
            partial_total(1:nbf), hmat_abs_limit, dg_frag%id, ifrag, ispin, jo)

          partial_block(1:nbf, jo) = partial_total(1:nbf)
        end do

        if (allocated(map_x)) deallocate(map_x)
        if (allocated(map_y)) deallocate(map_y)
        if (allocated(map_z)) deallocate(map_z)

        call assert_real_hmat_reconstruct_bounded_2d("partial_block", &
          partial_block(1:nbf, 1:nbf), hmat_abs_limit, dg_frag%id, ifrag, ispin, 0)
        call cpu_time(t0)
        call comm_summation(partial_block(1:nbf, 1:nbf), reduced_block(1:nbf, 1:nbf), nbf * nbf, dg_frag%icomm_frag)
        call cpu_time(t1)
        time_subgroup_reduce = time_subgroup_reduce + (t1 - t0)
        call assert_real_hmat_reconstruct_bounded_2d("reduced_block", &
          reduced_block(1:nbf, 1:nbf), hmat_abs_limit, dg_frag%id, ifrag, ispin, 0)

        if (.not. debug_static_seed_logged .and. dg_frag%is_frag_root .and. ispin == 1 .and. nbf >= 3) then
          write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x),a,3(1pe14.6,1x))') &
            "        reconstruct-diag probe: rank=", dg_frag%id, " ifrag=", ifrag, " seed_t=", &
            dg_frag%H_mat_blocks(iblk)%val(1,1,ispin), &
            dg_frag%H_mat_blocks(iblk)%val(2,2,ispin), &
            dg_frag%H_mat_blocks(iblk)%val(3,3,ispin), &
            " add_block=", reduced_block(1,1), reduced_block(2,2), reduced_block(3,3)
          flush(6)
        end if

        if (dg_frag%is_frag_root) then
          dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) = &
            dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) + reduced_block(1:nbf, 1:nbf)
          call assert_real_hmat_reconstruct_bounded_2d("final_h_block", &
            dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin), &
            hmat_abs_limit, dg_frag%id, ifrag, ispin, 0)
          if (.not. debug_static_seed_logged .and. ispin == 1 .and. nbf >= 3) then
            write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x))') &
              "        reconstruct-diag probe: rank=", dg_frag%id, " ifrag=", ifrag, " final_h=", &
              dg_frag%H_mat_blocks(iblk)%val(1,1,ispin), &
              dg_frag%H_mat_blocks(iblk)%val(2,2,ispin), &
              dg_frag%H_mat_blocks(iblk)%val(3,3,ispin)
            flush(6)
            debug_static_seed_logged = .true.
          end if
        end if
      end do
    end do

    if (.not. allocated(dg_frag%H_mat_blocks) .or. .not. allocated(dg_frag%H_block_map)) then
      call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    end if
    call cpu_time(t0)
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-reconstruct", dg_frag%icomm_frag)
    call cpu_time(t1)
    time_global_reduce = time_global_reduce + (t1 - t0)
    call cpu_time(t0)
    call rebuild_local_h_block_ids(dg_frag)
    call zero_nonlocal_h_matrix_blocks(dg_frag)
    call cpu_time(t1)
    time_post_reduce_cleanup = time_post_reduce_cleanup + (t1 - t0)
    if (enable_reconstruct_timing .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        reconstruct timing: halo=", time_halo_exchange, &
        " potential=", time_potential_build, " local=", time_local_build, &
        " subgroup_reduce=", time_subgroup_reduce, " global_reduce=", time_global_reduce, &
        " post_cleanup=", time_post_reduce_cleanup
      flush(6)
    end if

    call cpu_time(t0)
!$omp parallel do collapse(2) private(ispin,iblk,nbf,jo,io) schedule(static)
    do ispin = 1, system%nspin
      do iblk = 1, size(dg_frag%H_mat_blocks)
        if (dg_frag%H_mat_blocks(iblk)%ifrag_row /= dg_frag%H_mat_blocks(iblk)%ifrag_col) cycle
        nbf = dg_frag%n_basis(dg_frag%H_mat_blocks(iblk)%ifrag_row, ispin)
        do jo = 1, nbf
!$omp simd private(hij_avg)
          do io = jo + 1, nbf
            hij_avg = 0.5d0 * (dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) + &
                               dg_frag%H_mat_blocks(iblk)%val(jo, io, ispin))
            dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = hij_avg
            dg_frag%H_mat_blocks(iblk)%val(jo, io, ispin) = hij_avg
          end do
        end do
      end do
    end do
!$omp end parallel do
    call cpu_time(t1)
    time_block_hermitize = time_block_hermitize + (t1 - t0)

    if (enable_reconstruct_timing .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        reconstruct timing: hermitize=", time_block_hermitize
      flush(6)
    end if

    deallocate(V_total, V_phi)
    if (allocated(partial_total)) deallocate(partial_total)
    if (allocated(partial_block)) deallocate(partial_block)
    if (allocated(reduced_block)) deallocate(reduced_block)

  end subroutine reconstruct_hamiltonian_matrix

  integer function map_global_to_phi_box_coord_reconstruct(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord_reconstruct

  integer function map_global_to_rank_buffer_coord_reconstruct(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    if (lgtot <= 0 .or. ub < lb) then
      iloc = 0
      return
    end if

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_rank_buffer_coord_reconstruct

  subroutine build_total_potential_grid_with_stage_buffers(grid, lgnum_total, &
      Vh_buffer, Vxc_buffer, Vpsl_buffer, ispin, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    integer, intent(in) :: lgnum_total(3)
    real(8), intent(in) :: Vh_buffer(:,:,:)
    real(8), intent(in) :: Vxc_buffer(:,:,:,:)
    real(8), intent(in) :: Vpsl_buffer(:,:,:)
    integer, intent(in) :: ispin
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz
    integer :: bx, by, bz
    integer :: b_lo1, b_lo2, b_lo3, b_hi1, b_hi2, b_hi3

    b_lo1 = lbound(Vh_buffer, 1)
    b_hi1 = ubound(Vh_buffer, 1)
    b_lo2 = lbound(Vh_buffer, 2)
    b_hi2 = ubound(Vh_buffer, 2)
    b_lo3 = lbound(Vh_buffer, 3)
    b_hi3 = ubound(Vh_buffer, 3)

!$omp parallel do private(ix,iy,bx,by,bz) schedule(static)
    do iz = grid%is(3), grid%ie(3)
      bz = map_global_to_rank_buffer_coord_reconstruct(iz, b_lo3, b_hi3, lgnum_total(3))
      do iy = grid%is(2), grid%ie(2)
        by = map_global_to_rank_buffer_coord_reconstruct(iy, b_lo2, b_hi2, lgnum_total(2))
!$omp simd private(bx)
        do ix = grid%is(1), grid%ie(1)
          bx = map_global_to_rank_buffer_coord_reconstruct(ix, b_lo1, b_hi1, lgnum_total(1))
          if (bx == 0 .or. by == 0 .or. bz == 0) then
            V_total(ix, iy, iz) = 0.0d0
          else
            V_total(ix, iy, iz) = Vpsl_buffer(bx, by, bz) + Vh_buffer(bx, by, bz) + &
                                  Vxc_buffer(bx, by, bz, ispin)
          end if
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_total_potential_grid_with_stage_buffers

  subroutine build_local_potential_applied_basis_mapped(dg_frag, i_local, jo, mg, V_total, V_phi, &
      lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi, ov_s, ov_e, map_x, map_y, map_z)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, jo
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer, intent(in) :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer, intent(in) :: ov_s(3), ov_e(3)
    integer, intent(in) :: map_x(ov_s(1):ov_e(1)), map_y(ov_s(2):ov_e(2)), map_z(ov_s(3):ov_e(3))
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    complex(8) :: phi0

    V_phi(:, :, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%phi_frag_c)) then
!$omp parallel do private(lz, ly, lx, gx, gy, gz, bx, by, bz) schedule(static)
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        bz = map_z(gz)
        if (bz == 0) cycle
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          by = map_y(gy)
          if (by == 0) cycle
!$omp simd private(gx, bx, phi0)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_x(gx)
            if (bx == 0) cycle
            phi0 = dg_frag%phi_frag_c(bx, by, bz, jo, i_local)
            V_phi(gx, gy, gz) = V_total(gx, gy, gz) * phi0
          end do
        end do
      end do
!$omp end parallel do
    else
!$omp parallel do private(lz, ly, lx, gx, gy, gz, bx, by, bz) schedule(static)
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        bz = map_z(gz)
        if (bz == 0) cycle
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          by = map_y(gy)
          if (by == 0) cycle
!$omp simd private(gx, bx)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_x(gx)
            if (bx == 0) cycle
            V_phi(gx, gy, gz) = cmplx(V_total(gx, gy, gz) * dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
          end do
        end do
      end do
!$omp end parallel do
    end if
  end subroutine build_local_potential_applied_basis_mapped

  subroutine integrate_local_basis_with_field_mapped(dg_frag, i_local, io, mg, field, hvol, integral, &
      lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi, ov_s, ov_e, map_x, map_y, map_z)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, io
    type(s_rgrid), intent(in) :: mg
    complex(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral
    integer, intent(in) :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer, intent(in) :: ov_s(3), ov_e(3)
    integer, intent(in) :: map_x(ov_s(1):ov_e(1)), map_y(ov_s(2):ov_e(2)), map_z(ov_s(3):ov_e(3))
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    real(8) :: acc_re, acc_im, pr, pi, fr, fi, phi_r

    integral = (0.0d0, 0.0d0)
    if (allocated(dg_frag%phi_frag_c)) then
      acc_re = 0.0d0
      acc_im = 0.0d0
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        bz = map_z(gz)
        if (bz == 0) cycle
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          by = map_y(gy)
          if (by == 0) cycle
!$omp simd reduction(+:acc_re,acc_im) private(gx,bx,pr,pi,fr,fi)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_x(gx)
            if (bx == 0) cycle
            pr = real(dg_frag%phi_frag_c(bx, by, bz, io, i_local), kind=8)
            pi = aimag(dg_frag%phi_frag_c(bx, by, bz, io, i_local))
            fr = real(field(gx, gy, gz), kind=8)
            fi = aimag(field(gx, gy, gz))
            acc_re = acc_re + (pr * fr + pi * fi) * hvol
            acc_im = acc_im + (pr * fi - pi * fr) * hvol
          end do
        end do
      end do
      integral = cmplx(acc_re, acc_im, kind=8)
    else
      acc_re = 0.0d0
      acc_im = 0.0d0
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        bz = map_z(gz)
        if (bz == 0) cycle
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          by = map_y(gy)
          if (by == 0) cycle
!$omp simd reduction(+:acc_re,acc_im) private(gx,bx,phi_r,fr,fi)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_x(gx)
            if (bx == 0) cycle
            phi_r = dg_frag%phi_frag(bx, by, bz, io, i_local)
            fr = real(field(gx, gy, gz), kind=8)
            fi = aimag(field(gx, gy, gz))
            acc_re = acc_re + phi_r * fr * hvol
            acc_im = acc_im + phi_r * fi * hvol
          end do
        end do
      end do
      integral = cmplx(acc_re, acc_im, kind=8)
    end if
  end subroutine integrate_local_basis_with_field_mapped

  subroutine assert_real_hmat_reconstruct_bounded_1d(label, vals, limit, rank, ifrag, ispin, jo)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    real(8), intent(in) :: vals(:)
    real(8), intent(in) :: limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    real(8) :: vmax
    logical :: has_nan, has_big

    has_nan = any(ieee_is_nan(vals))
    has_big = any((.not. ieee_is_finite(vals)) .or. abs(vals) > limit)
    vmax = maxval(abs(vals))
    call stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
  end subroutine assert_real_hmat_reconstruct_bounded_1d

  subroutine assert_real_hmat_reconstruct_bounded_2d(label, vals, limit, rank, ifrag, ispin, jo)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    real(8), intent(in) :: vals(:,:)
    real(8), intent(in) :: limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    real(8) :: vmax
    logical :: has_nan, has_big

    has_nan = any(ieee_is_nan(vals))
    has_big = any((.not. ieee_is_finite(vals)) .or. abs(vals) > limit)
    vmax = maxval(abs(vals))
    call stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
  end subroutine assert_real_hmat_reconstruct_bounded_2d

  subroutine assert_real_hmat_reconstruct_bounded_3d(label, vals, limit, rank, ifrag, ispin, jo)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    real(8), intent(in) :: vals(:,:,:)
    real(8), intent(in) :: limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    real(8) :: vmax
    logical :: has_nan, has_big

    has_nan = any(ieee_is_nan(vals))
    has_big = any((.not. ieee_is_finite(vals)) .or. abs(vals) > limit)
    vmax = maxval(abs(vals))
    call stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
  end subroutine assert_real_hmat_reconstruct_bounded_3d

  subroutine assert_real_hmat_reconstruct_bounded_4d(label, vals, limit, rank, ifrag, ispin, jo)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    real(8), intent(in) :: vals(:,:,:,:)
    real(8), intent(in) :: limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    real(8) :: vmax
    logical :: has_nan, has_big

    has_nan = any(ieee_is_nan(vals))
    has_big = any((.not. ieee_is_finite(vals)) .or. abs(vals) > limit)
    vmax = maxval(abs(vals))
    call stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
  end subroutine assert_real_hmat_reconstruct_bounded_4d

  subroutine assert_complex_hmat_reconstruct_bounded_3d(label, vals, limit, rank, ifrag, ispin, jo)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    complex(8), intent(in) :: vals(:,:,:)
    real(8), intent(in) :: limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    real(8) :: vmax
    logical :: has_nan, has_big

    has_nan = any(ieee_is_nan(real(vals, kind=8))) .or. any(ieee_is_nan(aimag(vals)))
    has_big = any((.not. ieee_is_finite(real(vals, kind=8))) .or. &
                  (.not. ieee_is_finite(aimag(vals))) .or. abs(vals) > limit)
    vmax = maxval(abs(vals))
    call stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
  end subroutine assert_complex_hmat_reconstruct_bounded_3d

  subroutine stop_if_hmat_reconstruct_unbounded(label, vmax, limit, rank, ifrag, ispin, jo, has_nan, has_big)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    character(*), intent(in) :: label
    real(8), intent(in) :: vmax, limit
    integer, intent(in) :: rank, ifrag, ispin, jo
    logical, intent(in) :: has_nan, has_big

    if (has_nan .or. has_big .or. ieee_is_nan(vmax) .or. (.not. ieee_is_finite(vmax)) .or. vmax > limit) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,l1,a,l1,a,es24.16,a,es24.16)') &
        "[FATAL] invalid reduced DG H reconstruct input: label=", trim(label), &
        " rank=", rank, " ifrag=", ifrag, " ispin=", ispin, " jo=", jo, &
        " has_nan=", has_nan, " has_big=", has_big, " max=", vmax, " limit=", limit
      flush(6)
      stop "invalid DG H reconstruct input"
    end if
  end subroutine stop_if_hmat_reconstruct_unbounded
