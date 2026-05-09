  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, Vh_buffer)
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

    type(s_rgrid), pointer :: mg
    integer :: ifrag, ispin, io, jo, i_local
    integer :: nbf, nbf_raw, iblk
    integer :: max_nbf_local
    integer :: jo_s, jo_e
    integer :: loop_ifrag_start, loop_ifrag_end
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
    logical :: is_local_fragment
    logical :: block_is_local_fragment
    logical :: use_block_reconstruct
    logical, parameter :: enable_reconstruct_timing = .false.
    logical, parameter :: enable_hmat_nan_check = .false.

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
      block_is_local_fragment = (dg_frag%H_mat_blocks(iblk)%ifrag_row >= dg_frag%ifrag_start .and. &
                                 dg_frag%H_mat_blocks(iblk)%ifrag_row <= dg_frag%ifrag_end)
      if (dg_frag%is_frag_root .and. (.not. dg_frag%parallel_mode_orbital .or. block_is_local_fragment)) then
        ! In orbital mode each fragment subgroup contributes only its own
        ! block to the final world reduction.  Copying all replicated kinetic
        ! blocks from every subgroup root would multiply T by the number of
        ! fragment roots.
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :)
      end if
    end do

    if (.not. dg_frag%parallel_mode_orbital) then
      allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    end if

    max_nbf_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      nbf_raw = dg_frag%n_basis(ifrag, 1)
      if (system%nspin > 1) nbf_raw = max(nbf_raw, dg_frag%n_basis(ifrag, 2))
      nbf = min(nbf_raw, dg_frag%nstate_frag, size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
      max_nbf_local = max(max_nbf_local, nbf)
    end do
    if (max_nbf_local > 0) then
      allocate(partial_total(max_nbf_local), partial_block(max_nbf_local, max_nbf_local), reduced_block(max_nbf_local, max_nbf_local))
    end if

    frag_size = max(1, dg_frag%isize_frag)
    frag_rank = modulo(dg_frag%id_frag, frag_size)

    do ispin = 1, system%nspin
      call cpu_time(t0)
      if (.not. dg_frag%parallel_mode_orbital) then
        if (present(Vh_buffer)) then
          call build_total_potential_grid_with_buffered_hartree(mg, dg_frag, Vh_buffer, Vxc(ispin), Vpsl, V_total)
        else
          call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
        end if
      end if
      if (.not. dg_frag%parallel_mode_orbital .and. trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v)) then
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

      if (dg_frag%parallel_mode_orbital) then
        loop_ifrag_start = 1
        loop_ifrag_end = dg_frag%n_frag
      else
        loop_ifrag_start = dg_frag%ifrag_start
        loop_ifrag_end = dg_frag%ifrag_end
      end if

      do ifrag = loop_ifrag_start, loop_ifrag_end
        is_local_fragment = (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end)
        i_local = ifrag - dg_frag%ifrag_start + 1
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        if (dg_frag%parallel_mode_orbital) then
          ! All ranks enter this gather in global fragment order.  That lets
          ! the full parent-grid potential be assembled with the world
          ! communicator even though only the owning orbital subgroup will
          ! integrate matrix columns for this fragment.
          allocate(V_total(1:ndom(1), 1:ndom(2), 1:ndom(3)))
          allocate(V_phi(1:ndom(1), 1:ndom(2), 1:ndom(3)))
          call build_reconstruct_fragment_total_potential_grid(dg_frag, ifrag, mg, Vh, Vxc(ispin), Vpsl, V_total)
          if (.not. is_local_fragment) then
            deallocate(V_total, V_phi)
            cycle
          end if
        end if
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
        if (nbf_raw < 0) then
          write(*,'(a,i0,a,i0,a,i0,a,i0)') "[FATAL] reconstruct negative n_basis, rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf_raw
          stop "negative n_basis in reconstruct_hamiltonian_matrix"
        end if

        nbf = min(nbf_raw, dg_frag%nstate_frag, size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
        if (nbf <= 0) then
          if (dg_frag%parallel_mode_orbital) deallocate(V_total, V_phi)
          cycle
        end if

        iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if (iblk <= 0) then
          if (dg_frag%parallel_mode_orbital) deallocate(V_total, V_phi)
          cycle
        end if

        partial_block(1:nbf, 1:nbf) = 0.0d0

        if (dg_frag%parallel_mode_orbital) then
          ! Orbital mode splits matrix construction by ket-side columns.
          ! Each subgroup rank first reconstructs the full fragment potential
          ! box, then integrates only its assigned columns.
          call get_reconstruct_column_range(dg_frag, nbf, jo_s, jo_e)
          do jo = jo_s, jo_e
            call cpu_time(t0)
            call build_fragment_potential_applied_basis_reconstruct(dg_frag, ifrag, i_local, jo, V_total, V_phi)

            partial_total(1:nbf) = 0.0d0
!$omp parallel do private(io, integral_v) schedule(static)
            do io = 1, nbf
              call integrate_fragment_basis_with_field_reconstruct(dg_frag, ifrag, i_local, io, V_phi, hvol, integral_v)
              partial_total(io) = real(integral_v, kind=8)
            end do
!$omp end parallel do
            call cpu_time(t1)
            time_local_build = time_local_build + (t1 - t0)

            partial_block(1:nbf, jo) = partial_total(1:nbf)
          end do
          deallocate(V_total, V_phi)
        else
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

            partial_block(1:nbf, jo) = partial_total(1:nbf)
          end do

          if (allocated(map_x)) deallocate(map_x)
          if (allocated(map_y)) deallocate(map_y)
          if (allocated(map_z)) deallocate(map_z)
        end if

        call cpu_time(t0)
        call comm_summation(partial_block(1:nbf, 1:nbf), reduced_block(1:nbf, 1:nbf), nbf * nbf, dg_frag%icomm_frag)
        call cpu_time(t1)
        time_subgroup_reduce = time_subgroup_reduce + (t1 - t0)

        if (dg_frag%is_frag_root) then
          dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) = dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) + reduced_block(1:nbf, 1:nbf)
        end if
      end do
    end do

    if (.not. allocated(dg_frag%H_mat_blocks) .or. .not. allocated(dg_frag%H_block_map)) then
      call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    end if
    call cpu_time(t0)
    if (dg_frag%parallel_mode_orbital) then
      ! Column-split orbital mode computes one fragment block per subgroup,
      ! but mixed-basis propagation/diagnostics need the rebuilt block set on
      ! every rank.  Finish with a world reduction, not a subgroup-only one.
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-reconstruct", dg_frag%icomm)
    else
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-reconstruct", dg_frag%icomm_frag)
    end if
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
    call trace_reconstruct_matrix_blocks_if_enabled(dg_frag, dg_frag%H_mat_blocks, "HRT")

    if (enable_reconstruct_timing .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        reconstruct timing: hermitize=", time_block_hermitize
      flush(6)
    end if

    if (allocated(V_total)) deallocate(V_total)
    if (allocated(V_phi)) deallocate(V_phi)
    if (allocated(partial_total)) deallocate(partial_total)
    if (allocated(partial_block)) deallocate(partial_block)
    if (allocated(reduced_block)) deallocate(reduced_block)

  end subroutine reconstruct_hamiltonian_matrix

  subroutine trace_reconstruct_matrix_blocks_if_enabled(dg_frag, blocks, label)
    use communication, only: comm_is_root
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    character(*), intent(in) :: label

    logical, save :: initialized = .false.
    logical, save :: enabled = .false.
    character(16) :: env_trace
    integer :: env_status
    integer :: iblk, ispin, nrow, ncol
    real(8) :: frob, maxabs

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_HMAT_BLOCK_TRACE', env_trace, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
        end select
      end if
      initialized = .true.
    end if
    if (.not. enabled) return
    if (.not. comm_is_root(dg_frag%id)) return

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = min(dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin), size(blocks(iblk)%val, 1))
        ncol = min(dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin), size(blocks(iblk)%val, 2))
        if (nrow <= 0 .or. ncol <= 0) cycle
        frob = sqrt(sum(blocks(iblk)%val(1:nrow, 1:ncol, ispin)**2))
        maxabs = maxval(abs(blocks(iblk)%val(1:nrow, 1:ncol, ispin)))
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6)') &
          '[HMAT-BLOCK] label=', trim(label), ' iblk=', iblk, &
          ' row=', blocks(iblk)%ifrag_row, ' col=', blocks(iblk)%ifrag_col, &
          ' ispin=', ispin, ' frob=', frob, ' max=', maxabs
      end do
    end do
    flush(6)
  end subroutine trace_reconstruct_matrix_blocks_if_enabled

  subroutine get_reconstruct_column_range(dg_frag, ncol, col_s, col_e)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ncol
    integer, intent(out) :: col_s, col_e

    integer :: base, extra, rank_in_frag, nworker

    if (ncol <= 0) then
      col_s = 1
      col_e = 0
      return
    end if
    if (.not. dg_frag%parallel_mode_orbital .or. dg_frag%isize_frag <= 1) then
      col_s = 1
      col_e = ncol
      return
    end if

    nworker = max(1, dg_frag%isize_frag)
    rank_in_frag = max(0, min(dg_frag%id_frag, nworker - 1))
    base = ncol / nworker
    extra = mod(ncol, nworker)
    if (rank_in_frag < extra) then
      col_s = rank_in_frag * (base + 1) + 1
      col_e = col_s + base
    else
      col_s = extra * (base + 1) + (rank_in_frag - extra) * base + 1
      col_e = col_s + base - 1
    end if
    if (col_s > ncol) then
      col_s = 1
      col_e = 0
    else
      col_e = min(col_e, ncol)
    end if
  end subroutine get_reconstruct_column_range

  subroutine build_reconstruct_fragment_total_potential_grid(dg_frag, ifrag, mg, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    use communication, only: comm_summation, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_rgrid), intent(in) :: mg
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(inout) :: V_total(:,:,:)

    integer :: lx, ly, lz, gx, gy, gz, gwx, gwy, gwz
    integer :: iorg(3), ndom(3)
    integer :: npt
    real(8) :: vh_val
    real(8), allocatable :: V_reduced(:,:,:)

    V_total(:, :, :) = 0.0d0
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)

!$omp parallel do private(lz,ly,lx,gz,gy,gx,gwz,gwy,gwx,vh_val) schedule(static)
    do lz = 1, ndom(3)
      gz = iorg(3) + lz - 1
      gwz = map_global_to_periodic_box_coord_reconstruct(gz, 1, dg_frag%lgnum_total(3))
      if (gwz < mg%is(3) .or. gwz > mg%ie(3)) cycle
      do ly = 1, ndom(2)
        gy = iorg(2) + ly - 1
        gwy = map_global_to_periodic_box_coord_reconstruct(gy, 1, dg_frag%lgnum_total(2))
        if (gwy < mg%is(2) .or. gwy > mg%ie(2)) cycle
!$omp simd private(gx,gwx,vh_val)
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          gwx = map_global_to_periodic_box_coord_reconstruct(gx, 1, dg_frag%lgnum_total(1))
          if (gwx >= mg%is(1) .and. gwx <= mg%ie(1)) then
            vh_val = Vh%f(gwx, gwy, gwz)
            V_total(lx, ly, lz) = Vpsl%f(gwx, gwy, gwz) + vh_val + Vxc_spin%f(gwx, gwy, gwz)
          end if
        end do
      end do
    end do
!$omp end parallel do

    if (dg_frag%parallel_mode_orbital .and. dg_frag%isize > 1 .and. dg_frag%icomm /= COMM_GROUP_NULL) then
      npt = size(V_total)
      allocate(V_reduced(lbound(V_total,1):ubound(V_total,1), &
                         lbound(V_total,2):ubound(V_total,2), &
                         lbound(V_total,3):ubound(V_total,3)))
      call comm_summation(V_total, V_reduced, npt, dg_frag%icomm)
      V_total(:, :, :) = V_reduced(:, :, :)
      deallocate(V_reduced)
    end if
    call trace_reconstruct_potential_box_if_enabled(dg_frag, ifrag, V_total)
  end subroutine build_reconstruct_fragment_total_potential_grid

  subroutine trace_reconstruct_potential_box_if_enabled(dg_frag, ifrag, V_total)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    real(8), intent(in) :: V_total(:,:,:)

    logical, save :: initialized = .false.
    logical, save :: enabled = .false.
    character(16) :: env_trace
    integer :: env_status
    real(8) :: vsum, vl2, vmin, vmax

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_HMAT_POT_TRACE', env_trace, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
        end select
      end if
      initialized = .true.
    end if
    if (.not. enabled) return
    if (.not. comm_is_root(dg_frag%id)) return

    vsum = sum(V_total)
    vl2 = sqrt(sum(V_total * V_total))
    vmin = minval(V_total)
    vmax = maxval(V_total)
    write(*,'(1x,a,i0,a,i0,4(a,1pe14.6))') &
      '[HMAT-POT] ifrag=', ifrag, ' pfrag=', dg_frag%isize_frag, &
      ' sum=', vsum, ' l2=', vl2, ' min=', vmin, ' max=', vmax
    flush(6)
  end subroutine trace_reconstruct_potential_box_if_enabled

  subroutine build_fragment_potential_applied_basis_reconstruct(dg_frag, ifrag, i_local, jo, V_total, V_phi)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    real(8), intent(in) :: V_total(:,:,:)
    complex(8), intent(out) :: V_phi(:,:,:)

    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: iorg(3), ndom(3)
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    complex(8) :: phi0

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    V_phi(:, :, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%phi_frag_c)) then
!$omp parallel do private(lz,ly,lx,gz,gy,gx,bz,by,bx,phi0) schedule(static)
      do lz = 1, ndom(3)
        gz = iorg(3) + lz - 1
        bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        if (bz == 0) cycle
        do ly = 1, ndom(2)
          gy = iorg(2) + ly - 1
          by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          if (by == 0) cycle
!$omp simd private(lx,gx,bx,phi0)
          do lx = 1, ndom(1)
            gx = iorg(1) + lx - 1
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            if (bx == 0) cycle
            phi0 = dg_frag%phi_frag_c(bx, by, bz, jo, i_local)
            V_phi(lx, ly, lz) = V_total(lx, ly, lz) * phi0
          end do
        end do
      end do
!$omp end parallel do
    else
!$omp parallel do private(lz,ly,lx,gz,gy,gx,bz,by,bx) schedule(static)
      do lz = 1, ndom(3)
        gz = iorg(3) + lz - 1
        bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        if (bz == 0) cycle
        do ly = 1, ndom(2)
          gy = iorg(2) + ly - 1
          by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          if (by == 0) cycle
!$omp simd private(lx,gx,bx)
          do lx = 1, ndom(1)
            gx = iorg(1) + lx - 1
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            if (bx == 0) cycle
            V_phi(lx, ly, lz) = cmplx(V_total(lx, ly, lz) * dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
          end do
        end do
      end do
!$omp end parallel do
    end if
  end subroutine build_fragment_potential_applied_basis_reconstruct

  subroutine integrate_fragment_basis_with_field_reconstruct(dg_frag, ifrag, i_local, io, field, hvol, integral)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    complex(8), intent(in) :: field(:,:,:)
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral

    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: iorg(3), ndom(3)
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    real(8) :: acc_re, acc_im, pr, pi, fr, fi, phi_r

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    acc_re = 0.0d0
    acc_im = 0.0d0
    if (allocated(dg_frag%phi_frag_c)) then
      do lz = 1, ndom(3)
        gz = iorg(3) + lz - 1
        bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        if (bz == 0) cycle
        do ly = 1, ndom(2)
          gy = iorg(2) + ly - 1
          by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          if (by == 0) cycle
!$omp simd reduction(+:acc_re,acc_im) private(lx,gx,bx,pr,pi,fr,fi)
          do lx = 1, ndom(1)
            gx = iorg(1) + lx - 1
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            if (bx == 0) cycle
            pr = real(dg_frag%phi_frag_c(bx, by, bz, io, i_local), kind=8)
            pi = aimag(dg_frag%phi_frag_c(bx, by, bz, io, i_local))
            fr = real(field(lx, ly, lz), kind=8)
            fi = aimag(field(lx, ly, lz))
            acc_re = acc_re + (pr * fr + pi * fi) * hvol
            acc_im = acc_im + (pr * fi - pi * fr) * hvol
          end do
        end do
      end do
    else
      do lz = 1, ndom(3)
        gz = iorg(3) + lz - 1
        bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        if (bz == 0) cycle
        do ly = 1, ndom(2)
          gy = iorg(2) + ly - 1
          by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          if (by == 0) cycle
!$omp simd reduction(+:acc_re,acc_im) private(lx,gx,bx,phi_r,fr,fi)
          do lx = 1, ndom(1)
            gx = iorg(1) + lx - 1
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            if (bx == 0) cycle
            phi_r = dg_frag%phi_frag(bx, by, bz, io, i_local)
            fr = real(field(lx, ly, lz), kind=8)
            fi = aimag(field(lx, ly, lz))
            acc_re = acc_re + phi_r * fr * hvol
            acc_im = acc_im + phi_r * fi * hvol
          end do
        end do
      end do
    end if
    integral = cmplx(acc_re, acc_im, kind=8)
  end subroutine integrate_fragment_basis_with_field_reconstruct

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

  integer function map_global_to_periodic_box_coord_reconstruct(ig, lb, ub) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub
    integer :: extent

    extent = ub - lb + 1
    if (extent <= 0) then
      iloc = lb
      return
    end if
    iloc = modulo(ig - lb, extent) + lb
  end function map_global_to_periodic_box_coord_reconstruct

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
