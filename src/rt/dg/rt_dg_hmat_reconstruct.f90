  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_summation
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
    integer :: ifrag, ispin, io, jo, i_local, ig_i, ig_j
    integer :: nbf, nbf_raw, iblk, idx_io, idx_jo, valid_basis_count
    integer :: nmat_chk
    real(8) :: hvol
    complex(8) :: integral_t, integral_h
    real(8) :: max_asym_vpsl, max_asym_vh, max_asym_vxc
    real(8), allocatable :: Vpsl_mat(:,:,:), Vh_mat(:,:,:), Vxc_mat(:,:,:)
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: T_phi(:,:,:), H_phi(:,:,:)
    complex(8), allocatable :: Vpsl_phi(:,:,:), Vh_phi(:,:,:), Vxc_phi(:,:,:)
    real(8), allocatable :: partial_total(:), reduced_total(:)
    real(8), allocatable :: partial_vpsl(:), reduced_vpsl(:)
    real(8), allocatable :: partial_vh(:), reduced_vh(:)
    real(8), allocatable :: partial_vxc(:), reduced_vxc(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    logical :: release_dense_h, use_block_reconstruct

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. associated(dg_frag%mg)) then
      stop "reconstruct_hamiltonian_matrix requires dg_frag%mg"
    end if
    mg => dg_frag%mg

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

    release_dense_h = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0) .and. &
      (yn_hse /= 'y')
    use_block_reconstruct = allocated(dg_frag%H_mat_blocks) .and. allocated(dg_frag%H_mat_kinetic_blocks) .and. &
      allocated(dg_frag%H_block_map) .and. release_dense_h

    if ((.not. use_block_reconstruct) .and. (.not. allocated(dg_frag%H_mat))) then
      allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if

    if (use_block_reconstruct) then
      do iblk = 1, size(dg_frag%H_mat_blocks)
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :)
      end do
    else
      dg_frag%H_mat(:, :, :) = dg_frag%H_mat_kinetic(:, :, :)
    end if

    nmat_chk = dg_frag%n_mat_max
    if (.not. use_block_reconstruct) then
      allocate(Vpsl_mat(nmat_chk, nmat_chk, system%nspin))
      allocate(Vh_mat(nmat_chk, nmat_chk, system%nspin))
      allocate(Vxc_mat(nmat_chk, nmat_chk, system%nspin))
      Vpsl_mat = 0.0d0
      Vh_mat = 0.0d0
      Vxc_mat = 0.0d0
      allocate(basis_gid(dg_frag%nstate_frag), valid_basis_ids(dg_frag%nstate_frag))
    end if

    allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    if (.not. use_block_reconstruct) then
      allocate(Vpsl_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(Vh_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(Vxc_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    end if

    call exchange_phi_frag_halo(dg_frag)

    do ispin = 1, system%nspin
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)

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

        if (.not. use_block_reconstruct) then
          valid_basis_count = 0
          do io = 1, nbf
            basis_gid(io) = dg_frag%index_basis(io, ifrag, ispin)
            if (basis_gid(io) < 1 .or. basis_gid(io) > dg_frag%n_mat_max) cycle
            valid_basis_count = valid_basis_count + 1
            valid_basis_ids(valid_basis_count) = io
          end do
        else
          iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        end if

        allocate(partial_total(nbf), reduced_total(nbf))
        if (.not. use_block_reconstruct) then
          allocate(partial_vpsl(nbf), reduced_vpsl(nbf))
          allocate(partial_vh(nbf), reduced_vh(nbf))
          allocate(partial_vxc(nbf), reduced_vxc(nbf))
        end if

        do jo = 1, nbf
          call build_local_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
          if (.not. use_block_reconstruct) then
            call build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, Vpsl%f, Vpsl_phi)
            call build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, Vh%f, Vh_phi)
            call build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, Vxc(ispin)%f, Vxc_phi)
          end if

          partial_total(:) = 0.0d0
          if (.not. use_block_reconstruct) then
            partial_vpsl(:) = 0.0d0
            partial_vh(:) = 0.0d0
            partial_vxc(:) = 0.0d0
          end if

!$omp parallel do private(io, integral_t, integral_h) schedule(static)
          do io = 1, nbf
            call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, integral_t)
            call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, integral_h)
            partial_total(io) = real(integral_h - integral_t, kind=8)
            if (.not. use_block_reconstruct) then
              call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, Vpsl_phi, hvol, integral_h)
              partial_vpsl(io) = real(integral_h, kind=8)
              call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, Vh_phi, hvol, integral_h)
              partial_vh(io) = real(integral_h, kind=8)
              call integrate_local_basis_with_field(dg_frag, ifrag, i_local, io, mg, Vxc_phi, hvol, integral_h)
              partial_vxc(io) = real(integral_h, kind=8)
            end if
          end do
!$omp end parallel do

          call comm_summation(partial_total, reduced_total, nbf, dg_frag%icomm_frag)
          if (.not. use_block_reconstruct) then
            call comm_summation(partial_vpsl, reduced_vpsl, nbf, dg_frag%icomm_frag)
            call comm_summation(partial_vh, reduced_vh, nbf, dg_frag%icomm_frag)
            call comm_summation(partial_vxc, reduced_vxc, nbf, dg_frag%icomm_frag)
          end if

          if (dg_frag%is_frag_root) then
            if (use_block_reconstruct) then
              if (iblk > 0) then
                do io = 1, nbf
                  dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) + reduced_total(io)
                end do
              end if
            else
              ig_j = dg_frag%index_basis(jo, ifrag, ispin)
              if (ig_j >= 1 .and. ig_j <= dg_frag%n_mat_max) then
                do idx_io = 1, valid_basis_count
                  io = valid_basis_ids(idx_io)
                  ig_i = basis_gid(io)
                  dg_frag%H_mat(ig_i, ig_j, ispin) = dg_frag%H_mat(ig_i, ig_j, ispin) + reduced_total(io)
                  Vpsl_mat(ig_i, ig_j, ispin) = Vpsl_mat(ig_i, ig_j, ispin) + reduced_vpsl(io)
                  Vh_mat(ig_i, ig_j, ispin) = Vh_mat(ig_i, ig_j, ispin) + reduced_vh(io)
                  Vxc_mat(ig_i, ig_j, ispin) = Vxc_mat(ig_i, ig_j, ispin) + reduced_vxc(io)
                end do
              end if
            end if
          end if
        end do

        if ((.not. use_block_reconstruct) .and. (yn_hse == 'y')) then
          call add_exact_exchange_hse(dg_frag, system, dg_frag%H_mat(:, :, ispin), ifrag, ispin)
        end if

        deallocate(partial_total, reduced_total)
        if (.not. use_block_reconstruct) then
          deallocate(partial_vpsl, reduced_vpsl, partial_vh, reduced_vh, partial_vxc, reduced_vxc)
        end if
      end do
    end do

    if (.not. allocated(dg_frag%H_mat_blocks) .or. .not. allocated(dg_frag%H_block_map)) then
      call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    end if
    if (.not. use_block_reconstruct) then
      call sync_dense_matrix_to_blocks(dg_frag, dg_frag%H_mat, dg_frag%H_mat_blocks, dg_frag%H_block_map)
    end if
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-reconstruct", dg_frag%icomm)
    call rebuild_local_h_block_ids(dg_frag)
    call zero_nonlocal_h_matrix_blocks(dg_frag)
    if (.not. use_block_reconstruct) then
      call sync_blocks_to_dense_matrix(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%H_mat)
    end if

    if (.not. use_block_reconstruct) then
      do ispin = 1, system%nspin
        max_asym_vpsl = 0.0d0
        max_asym_vh = 0.0d0
        max_asym_vxc = 0.0d0
        !$omp parallel do collapse(2) private(io,jo) reduction(max:max_asym_vpsl,max_asym_vh,max_asym_vxc) schedule(static)
        do jo = 1, nmat_chk
          do io = 1, nmat_chk
            max_asym_vpsl = max(max_asym_vpsl, abs(Vpsl_mat(io, jo, ispin) - Vpsl_mat(jo, io, ispin)))
            max_asym_vh = max(max_asym_vh, abs(Vh_mat(io, jo, ispin) - Vh_mat(jo, io, ispin)))
            max_asym_vxc = max(max_asym_vxc, abs(Vxc_mat(io, jo, ispin) - Vxc_mat(jo, io, ispin)))
          end do
        end do
        !$omp end parallel do
        if (max_asym_vpsl > 1.0d-8 .or. max_asym_vh > 1.0d-8 .or. max_asym_vxc > 1.0d-8) then
          write(*,'(1x,a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4)') "[WARN] term asymmetry: rank=", &
            dg_frag%id, " ispin=", ispin, " Vpsl=", max_asym_vpsl, " Vh=", max_asym_vh, &
            " Vxc=", max_asym_vxc
        end if
      end do
      deallocate(Vpsl_mat, Vh_mat, Vxc_mat)
    end if

    if (use_block_reconstruct) then
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
    else
      block
        integer :: io_sym, jo_sym
!$omp parallel do collapse(2) private(ispin,io_sym,jo_sym) schedule(static)
        do ispin = 1, system%nspin
          do jo_sym = 1, nmat_chk
            do io_sym = jo_sym + 1, nmat_chk
              dg_frag%H_mat(io_sym, jo_sym, ispin) = 0.5d0 * (dg_frag%H_mat(io_sym, jo_sym, ispin) + &
                                                               dg_frag%H_mat(jo_sym, io_sym, ispin))
              dg_frag%H_mat(jo_sym, io_sym, ispin) = dg_frag%H_mat(io_sym, jo_sym, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end block
    end if

    if ((.not. use_block_reconstruct) .and. allocated(dg_frag%H_mat_c)) then
      dg_frag%H_mat_c(:, :, :) = cmplx(dg_frag%H_mat(:, :, :), 0.0d0, kind=8)
    end if
    if (release_dense_h) then
      if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    end if

    deallocate(V_total, T_phi, H_phi)
    if (.not. use_block_reconstruct) then
      deallocate(Vpsl_phi, Vh_phi, Vxc_phi)
      deallocate(basis_gid, valid_basis_ids)
    end if

  end subroutine reconstruct_hamiltonian_matrix

  subroutine build_local_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz
    complex(8) :: phi0, lap_sum
    real(8) :: lap0, lapt(4,3)
    logical :: has_overlap

    T_phi(:, :, :) = (0.0d0, 0.0d0)
    H_phi(:, :, :) = (0.0d0, 0.0d0)
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
    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
!$omp parallel do collapse(2) private(lz, ly, lx, gx, gy, gz) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
!$omp simd private(gx, phi0, lap_sum)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          if (allocated(dg_frag%phi_frag_c)) then
            phi0 = dg_frag%phi_frag_c(lx, ly, lz, jo, i_local)
            lap_sum = lapt(1,1) * (dg_frag%phi_frag_c(lx+1, ly, lz, jo, i_local) + dg_frag%phi_frag_c(lx-1, ly, lz, jo, i_local)) + &
                      lapt(2,1) * (dg_frag%phi_frag_c(lx+2, ly, lz, jo, i_local) + dg_frag%phi_frag_c(lx-2, ly, lz, jo, i_local)) + &
                      lapt(3,1) * (dg_frag%phi_frag_c(lx+3, ly, lz, jo, i_local) + dg_frag%phi_frag_c(lx-3, ly, lz, jo, i_local)) + &
                      lapt(4,1) * (dg_frag%phi_frag_c(lx+4, ly, lz, jo, i_local) + dg_frag%phi_frag_c(lx-4, ly, lz, jo, i_local)) + &
                      lapt(1,2) * (dg_frag%phi_frag_c(lx, ly+1, lz, jo, i_local) + dg_frag%phi_frag_c(lx, ly-1, lz, jo, i_local)) + &
                      lapt(2,2) * (dg_frag%phi_frag_c(lx, ly+2, lz, jo, i_local) + dg_frag%phi_frag_c(lx, ly-2, lz, jo, i_local)) + &
                      lapt(3,2) * (dg_frag%phi_frag_c(lx, ly+3, lz, jo, i_local) + dg_frag%phi_frag_c(lx, ly-3, lz, jo, i_local)) + &
                      lapt(4,2) * (dg_frag%phi_frag_c(lx, ly+4, lz, jo, i_local) + dg_frag%phi_frag_c(lx, ly-4, lz, jo, i_local)) + &
                      lapt(1,3) * (dg_frag%phi_frag_c(lx, ly, lz+1, jo, i_local) + dg_frag%phi_frag_c(lx, ly, lz-1, jo, i_local)) + &
                      lapt(2,3) * (dg_frag%phi_frag_c(lx, ly, lz+2, jo, i_local) + dg_frag%phi_frag_c(lx, ly, lz-2, jo, i_local)) + &
                      lapt(3,3) * (dg_frag%phi_frag_c(lx, ly, lz+3, jo, i_local) + dg_frag%phi_frag_c(lx, ly, lz-3, jo, i_local)) + &
                      lapt(4,3) * (dg_frag%phi_frag_c(lx, ly, lz+4, jo, i_local) + dg_frag%phi_frag_c(lx, ly, lz-4, jo, i_local))
          else
            phi0 = cmplx(dg_frag%phi_frag(lx, ly, lz, jo, i_local), 0.0d0, kind=8)
            lap_sum = lapt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) + dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
                      lapt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) + dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
                      lapt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) + dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
                      lapt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) + dg_frag%phi_frag(lx-4, ly, lz, jo, i_local)) + &
                      lapt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) + dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
                      lapt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) + dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
                      lapt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) + dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
                      lapt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) + dg_frag%phi_frag(lx, ly-4, lz, jo, i_local)) + &
                      lapt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) + dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
                      lapt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) + dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
                      lapt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) + dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
                      lapt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) + dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))
          end if
          T_phi(gx, gy, gz) = lap0 * phi0 - 0.5d0 * lap_sum
          H_phi(gx, gy, gz) = T_phi(gx, gy, gz) + V_total(gx, gy, gz) * phi0
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_local_hpsi_for_basis

  subroutine build_local_potential_applied_basis(dg_frag, ifrag, i_local, jo, mg, field, field_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: field_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz
    logical :: has_overlap

    field_phi(:, :, :) = (0.0d0, 0.0d0)
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
!$omp parallel do collapse(2) private(lz, ly, lx, gx, gy, gz) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
!$omp simd private(gx)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          if (allocated(dg_frag%phi_frag_c)) then
            field_phi(gx, gy, gz) = field(gx, gy, gz) * dg_frag%phi_frag_c(lx, ly, lz, jo, i_local)
          else
            field_phi(gx, gy, gz) = field(gx, gy, gz) * cmplx(dg_frag%phi_frag(lx, ly, lz, jo, i_local), 0.0d0, kind=8)
          end if
        end do
      end do
    end do
!$omp end parallel do
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
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz
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
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          if (allocated(dg_frag%phi_frag_c)) then
            integral = integral + conjg(dg_frag%phi_frag_c(lx, ly, lz, io, i_local)) * field(gx, gy, gz) * hvol
          else
            integral = integral + cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end if
        end do
      end do
    end do
  end subroutine integrate_local_basis_with_field
