  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
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

    type(s_rgrid), pointer :: mg
    integer :: ifrag, ispin, io, jo, i_local
    integer :: nbf, nbf_raw, iblk
    real(8) :: hvol, A2val
    real(8) :: t0, t1
    real(8) :: time_local_build, time_subgroup_reduce, time_global_reduce
    real(8) :: time_halo_exchange, time_potential_build, time_post_reduce_cleanup, time_block_hermitize
    complex(8) :: integral_v
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: V_phi(:,:,:)
    real(8), allocatable :: partial_total(:), partial_block(:,:), reduced_block(:,:)
    logical :: use_block_reconstruct
    logical, save :: debug_static_seed_logged = .false.
    logical, parameter :: enable_reconstruct_timing = .false.

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

    do ispin = 1, system%nspin
      call cpu_time(t0)
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
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

        allocate(partial_total(nbf), partial_block(nbf, nbf), reduced_block(nbf, nbf))
        partial_block(:, :) = 0.0d0

        do jo = 1, nbf
          call cpu_time(t0)
          call build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, V_total, V_phi)

          partial_total(:) = 0.0d0

!$omp parallel do private(io, integral_v) schedule(static)
          do io = 1, nbf
            call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, V_phi, hvol, integral_v)
            partial_total(io) = real(integral_v, kind=8)
          end do
!$omp end parallel do
          call cpu_time(t1)
          time_local_build = time_local_build + (t1 - t0)

          partial_block(:, jo) = partial_total(:)
        end do

        call cpu_time(t0)
        call comm_summation(partial_block, reduced_block, nbf * nbf, dg_frag%icomm_frag)
        call cpu_time(t1)
        time_subgroup_reduce = time_subgroup_reduce + (t1 - t0)

        if (.not. debug_static_seed_logged .and. dg_frag%is_frag_root .and. ispin == 1 .and. nbf >= 3) then
          write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x),a,3(1pe14.6,1x))') &
            "        reconstruct-diag probe: rank=", dg_frag%id, " ifrag=", ifrag, " seed_t=", &
            dg_frag%H_mat_blocks(iblk)%val(1,1,ispin), dg_frag%H_mat_blocks(iblk)%val(2,2,ispin), dg_frag%H_mat_blocks(iblk)%val(3,3,ispin), &
            " add_block=", reduced_block(1,1), reduced_block(2,2), reduced_block(3,3)
          flush(6)
        end if

        if (dg_frag%is_frag_root) then
          dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) = dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) + reduced_block(:, :)
          if (.not. debug_static_seed_logged .and. ispin == 1 .and. nbf >= 3) then
            write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x))') &
              "        reconstruct-diag probe: rank=", dg_frag%id, " ifrag=", ifrag, " final_h=", &
              dg_frag%H_mat_blocks(iblk)%val(1,1,ispin), dg_frag%H_mat_blocks(iblk)%val(2,2,ispin), dg_frag%H_mat_blocks(iblk)%val(3,3,ispin)
            flush(6)
            debug_static_seed_logged = .true.
          end if
        end if

        deallocate(partial_total, partial_block, reduced_block)
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
          do io = jo + 1, nbf
            dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = 0.5d0 * &
              (dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) + dg_frag%H_mat_blocks(iblk)%val(jo, io, ispin))
            dg_frag%H_mat_blocks(iblk)%val(jo, io, ispin) = dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin)
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

  subroutine build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, V_total, V_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    complex(8) :: phi0
    logical :: has_overlap

    V_phi(:, :, :) = (0.0d0, 0.0d0)
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))
    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) return
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
    if (allocated(dg_frag%phi_frag_c)) then
!$omp parallel do private(lz, ly, lx, gx, gy, gz, bx, by, bz) schedule(static)
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
!$omp simd private(gx, phi0)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
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
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
!$omp simd private(gx, phi0)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            phi0 = cmplx(dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
            V_phi(gx, gy, gz) = V_total(gx, gy, gz) * phi0
          end do
        end do
      end do
!$omp end parallel do
    end if
  end subroutine build_local_potential_applied_basis

  subroutine integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, field, hvol, integral)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    type(s_rgrid), intent(in) :: mg
    complex(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    logical :: has_overlap

    integral = (0.0d0, 0.0d0)
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))
    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) return
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
    if (allocated(dg_frag%phi_frag_c)) then
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            integral = integral + conjg(dg_frag%phi_frag_c(bx, by, bz, io, i_local)) * field(gx, gy, gz) * hvol
          end do
        end do
      end do
    else
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_reconstruct(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_reconstruct(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_reconstruct(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            integral = integral + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end do
        end do
      end do
    end if
  end subroutine integrate_local_basis_with_field
