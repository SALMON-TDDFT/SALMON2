  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use salmon_global, only: yn_hse
    use rt_dg_fragment_ops, only: rebuild_local_h_block_ids, zero_nonlocal_h_matrix_blocks
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_scalar),         intent(in)    :: Vh
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    type(s_scalar),         intent(in)    :: Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    
    integer :: ifrag, ifrag_pack, ispin, io, jo, ix, iy, iz, idir, i_local, i_local_pack, ig_i, ig_j, nbf, nbf_raw, iblk, ifrag_count
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3), ndom_pack(3)
    real(8) :: hvol, A_squared
    real(8) :: integral_Vh, integral_Vxc, integral_Vpsl  ! Real integrals (real basis)
    integer :: is(3), ie(3)
    real(8), allocatable :: Vpsl_mat(:,:,:), Vh_mat(:,:,:), Vxc_mat(:,:,:)
    real(8), allocatable :: phi_blk(:,:), weighted_phi(:,:), weighted_phi_vpsl(:,:), weighted_phi_vh(:,:), &
                            weighted_phi_vxc(:,:), Vpsl_blk(:), Vh_blk(:), Vxc_blk(:), Vtot_blk(:)
    real(8), allocatable :: total_mat(:,:), Vpsl_local_mat(:,:), Vh_local_mat(:,:), Vxc_local_mat(:,:)
    real(8) :: max_asym_vpsl, max_asym_vh, max_asym_vxc
    integer :: nmat_chk
    integer :: ngrid, ngrid_max, igrid_cache
    logical :: release_dense_h, use_block_reconstruct
    
    ! Check if real-space basis functions are available
    if (.not. dg_frag%has_real_space_basis) then
      return  ! Cannot reconstruct without basis functions
    end if
    
    release_dense_h = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0) .and. &
      (yn_hse /= 'y')
    use_block_reconstruct = allocated(dg_frag%H_mat_blocks) .and. allocated(dg_frag%H_mat_kinetic_blocks) .and. &
      allocated(dg_frag%H_block_map) .and. release_dense_h
    ngrid_max = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(:, ifrag)))
    end do
    ifrag_count = max(1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)

    hvol = system%hvol
    is = dg_frag%lg%is
    ie = dg_frag%lg%ie
    if (hvol /= hvol) then
      write(*,'(a,i0)') "[NaN] hvol in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
      stop "NaN in hvol"
    end if
    if (any(Vpsl%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)) /= Vpsl%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)))) then
      write(*,'(a,i0)') "[NaN] Vpsl in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
      stop "NaN in Vpsl"
    end if
    if (any(Vh%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)) /= Vh%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)))) then
      write(*,'(a,i0)') "[NaN] Vh in reconstruct_hamiltonian_matrix, rank=", dg_frag%id
      stop "NaN in Vh"
    end if
    do ispin = 1, system%nspin
      if (any(Vxc(ispin)%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)) /= Vxc(ispin)%f(is(1):ie(1), is(2):ie(2), is(3):ie(3)))) then
        write(*,'(a,i0,a,i0)') "[NaN] Vxc in reconstruct_hamiltonian_matrix, rank=", dg_frag%id, " ispin=", ispin
        stop "NaN in Vxc"
      end if
    end do
    
    ! Reconstruct Hamiltonian matrix with POTENTIAL-DEPENDENT terms only:
    ! H_mat = T + V_psl + V_H + V_xc
    !
    ! This includes:
    ! 1. Kinetic energy: T (from H_mat_kinetic)
    ! 2. Pseudopotential: V_psl,ij = ∫ φ_i*(r) V_psl(r) φ_j(r) dr
    ! 3. Hartree potential: V_H,ij = ∫ φ_i*(r) V_H(r) φ_j(r) dr
    ! 4. XC potential: V_xc,ij = ∫ φ_i*(r) V_xc(r) φ_j(r) dr
    !
    ! Note: A(t)·p and A²(t)/2 terms are NOT included in H_mat
    !       They are added dynamically in calculate_time_derivative
    !       This avoids double-counting and ensures H_mat remains
    !       the time-independent part of the Hamiltonian
    !
    ! For basis update detection: A temporary H_mat_full is constructed
    ! including A·p and A² terms to capture complete Hamiltonian change
    
    if ((.not. use_block_reconstruct) .and. (.not. allocated(dg_frag%H_mat))) then
      allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if

    ! Step 1: Reset to kinetic parts
    if (use_block_reconstruct) then
      do iblk = 1, size(dg_frag%H_mat_blocks)
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :)
      end do
    else
      dg_frag%H_mat(:, :, :) = dg_frag%H_mat_kinetic(:, :, :)
    end if

    ! Allocate diagnostic matrices for potential terms
    nmat_chk = dg_frag%n_mat_max
    if (.not. use_block_reconstruct) then
      allocate(Vpsl_mat(nmat_chk, nmat_chk, system%nspin))
      allocate(Vh_mat(nmat_chk, nmat_chk, system%nspin))
      allocate(Vxc_mat(nmat_chk, nmat_chk, system%nspin))
      Vpsl_mat = 0.0d0
      Vh_mat = 0.0d0
      Vxc_mat = 0.0d0
    end if
    if (ngrid_max > 0) then
      allocate(phi_blk(ngrid_max, dg_frag%nstate_frag))
      allocate(weighted_phi(ngrid_max, dg_frag%nstate_frag))
      allocate(weighted_phi_vpsl(ngrid_max, dg_frag%nstate_frag))
      allocate(weighted_phi_vh(ngrid_max, dg_frag%nstate_frag))
      allocate(weighted_phi_vxc(ngrid_max, dg_frag%nstate_frag))
      allocate(Vpsl_blk(ngrid_max), Vh_blk(ngrid_max), Vxc_blk(ngrid_max), Vtot_blk(ngrid_max))
      allocate(total_mat(dg_frag%nstate_frag, dg_frag%nstate_frag))
      allocate(Vpsl_local_mat(dg_frag%nstate_frag, dg_frag%nstate_frag))
      allocate(Vh_local_mat(dg_frag%nstate_frag, dg_frag%nstate_frag))
      allocate(Vxc_local_mat(dg_frag%nstate_frag, dg_frag%nstate_frag))
    end if
    
    ! Step 2-4: Add updated potential matrix elements
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = ifrag - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
        write(*,'(a,i0,a,i0,a,i0)') "[FATAL] reconstruct invalid i_local, rank=", dg_frag%id, &
          " ifrag=", ifrag, " i_local=", i_local
        stop "invalid i_local in reconstruct_hamiltonian_matrix"
      end if
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      if (any(ndom <= 0)) then
        write(*,'(a,i0,a,i0,a,3(i0,1x))') "[FATAL] reconstruct non-positive ndom, rank=", dg_frag%id, &
          " ifrag=", ifrag, " ndom=", ndom(1), ndom(2), ndom(3)
        stop "non-positive ndom in reconstruct_hamiltonian_matrix"
      end if
      if (ndom(1) > size(dg_frag%phi_frag, 1) .or. ndom(2) > size(dg_frag%phi_frag, 2) .or. ndom(3) > size(dg_frag%phi_frag, 3)) then
        write(*,'(a,i0,a,i0,a,3(i0,1x),a,3(i0,1x))') "[FATAL] reconstruct ndom exceeds phi_frag bounds, rank=", dg_frag%id, &
          " ifrag=", ifrag, " ndom=", ndom(1), ndom(2), ndom(3), " phi_dims=", &
          size(dg_frag%phi_frag, 1), size(dg_frag%phi_frag, 2), size(dg_frag%phi_frag, 3)
        stop "ndom exceeds phi_frag bounds in reconstruct_hamiltonian_matrix"
      end if
      if (any(dg_frag%phi_frag(1:ndom(1), 1:ndom(2), 1:ndom(3), 1:dg_frag%nstate_frag, i_local) /= &
              dg_frag%phi_frag(1:ndom(1), 1:ndom(2), 1:ndom(3), 1:dg_frag%nstate_frag, i_local))) then
        write(*,'(a,i0,a,i0)') "[NaN] phi_frag in reconstruct_hamiltonian_matrix, rank=", dg_frag%id, " ifrag=", ifrag
        stop "NaN in phi_frag"
      end if
      ngrid = ndom(1) * ndom(2) * ndom(3)
      if ((.not. allocated(dg_frag%density_phi_cache)) .or. (.not. dg_frag%density_phi_cache_valid)) then
        if (.not. allocated(dg_frag%density_phi_cache)) then
          allocate(dg_frag%density_phi_cache(ngrid_max, dg_frag%nstate_frag, ifrag_count))
        end if
        dg_frag%density_phi_cache(:, :, :) = 0.0d0
        do ifrag_pack = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_pack = ifrag_pack - dg_frag%ifrag_start + 1
          ndom_pack(:) = dg_frag%nxyz_domain(:, ifrag_pack)
!$omp parallel do collapse(4) private(lx, ly, lz, io, igrid_cache) schedule(static)
          do lz = 1, ndom_pack(3)
            do ly = 1, ndom_pack(2)
              do lx = 1, ndom_pack(1)
                do io = 1, dg_frag%nstate_frag
                  igrid_cache = lx + ndom_pack(1) * ((ly - 1) + ndom_pack(2) * (lz - 1))
                  dg_frag%density_phi_cache(igrid_cache, io, i_local_pack) = dg_frag%phi_frag(lx, ly, lz, io, i_local_pack)
                end do
              end do
            end do
          end do
!$omp end parallel do
        end do
        dg_frag%density_phi_cache_valid = .true.
      end if
      do ispin = 1, system%nspin
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
        if (.false. .and. dg_frag%id == 0 .and. nbf_raw >= 60) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        reconstruct basis diag: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf_raw
          flush(6)
        end if
        if (nbf_raw < 0) then
          write(*,'(a,i0,a,i0,a,i0,a,i0)') "[FATAL] reconstruct negative n_basis, rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf_raw
          stop "negative n_basis in reconstruct_hamiltonian_matrix"
        end if
        nbf = min(nbf_raw, dg_frag%nstate_frag, size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
        if (nbf_raw > nbf) then
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "[WARN] reconstruct n_basis truncated, rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis_raw=", nbf_raw, " nbf_used=", nbf, &
            " index_dim1=", size(dg_frag%index_basis, 1), " phi_dim4=", size(dg_frag%phi_frag, 4)
        end if
        if (nbf <= 0) cycle
        phi_blk(1:ngrid, 1:nbf) = dg_frag%density_phi_cache(1:ngrid, 1:nbf, i_local)
        igrid_cache = 0
        do lz = 1, ndom(3)
          gz = modulo(iorg(3) + lz - 2, dg_frag%lgnum_total(3)) + 1
          do ly = 1, ndom(2)
            gy = modulo(iorg(2) + ly - 2, dg_frag%lgnum_total(2)) + 1
            do lx = 1, ndom(1)
              gx = modulo(iorg(1) + lx - 2, dg_frag%lgnum_total(1)) + 1
              igrid_cache = igrid_cache + 1
              Vpsl_blk(igrid_cache) = Vpsl%f(gx, gy, gz) * hvol
              Vh_blk(igrid_cache) = Vh%f(gx, gy, gz) * hvol
              Vxc_blk(igrid_cache) = Vxc(ispin)%f(gx, gy, gz) * hvol
              Vtot_blk(igrid_cache) = Vpsl_blk(igrid_cache) + Vh_blk(igrid_cache) + Vxc_blk(igrid_cache)
            end do
          end do
        end do

        do io = 1, nbf
          weighted_phi(1:ngrid, io) = Vtot_blk(1:ngrid) * phi_blk(1:ngrid, io)
        end do
        total_mat(1:nbf, 1:nbf) = matmul(transpose(phi_blk(1:ngrid, 1:nbf)), weighted_phi(1:ngrid, 1:nbf))
        if (.not. use_block_reconstruct) then
!$omp parallel do collapse(2) schedule(static)
          do io = 1, nbf
            do igrid_cache = 1, ngrid
              weighted_phi_vpsl(igrid_cache, io) = Vpsl_blk(igrid_cache) * phi_blk(igrid_cache, io)
              weighted_phi_vh(igrid_cache, io) = Vh_blk(igrid_cache) * phi_blk(igrid_cache, io)
              weighted_phi_vxc(igrid_cache, io) = Vxc_blk(igrid_cache) * phi_blk(igrid_cache, io)
            end do
          end do
!$omp end parallel do
          Vpsl_local_mat(1:nbf, 1:nbf) = matmul(transpose(phi_blk(1:ngrid, 1:nbf)), weighted_phi_vpsl(1:ngrid, 1:nbf))
          Vh_local_mat(1:nbf, 1:nbf) = matmul(transpose(phi_blk(1:ngrid, 1:nbf)), weighted_phi_vh(1:ngrid, 1:nbf))
          Vxc_local_mat(1:nbf, 1:nbf) = matmul(transpose(phi_blk(1:ngrid, 1:nbf)), weighted_phi_vxc(1:ngrid, 1:nbf))
        end if

        if (any(total_mat(1:nbf, 1:nbf) /= total_mat(1:nbf, 1:nbf))) then
          write(*,'(a,i0,a,i0,a,i0)') "[NaN] total potential mat, rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin
          stop "NaN in reconstructed total potential matrix"
        end if

        if (use_block_reconstruct) then
          iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
          if (iblk > 0) then
            dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) = dg_frag%H_mat_blocks(iblk)%val(1:nbf, 1:nbf, ispin) + &
              total_mat(1:nbf, 1:nbf)
          end if
        else
          do jo = 1, nbf
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            do io = 1, nbf
              ig_i = dg_frag%index_basis(io, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
              dg_frag%H_mat(ig_i, ig_j, ispin) = dg_frag%H_mat(ig_i, ig_j, ispin) + total_mat(io, jo)
              Vpsl_mat(ig_i, ig_j, ispin) = Vpsl_mat(ig_i, ig_j, ispin) + Vpsl_local_mat(io, jo)
              Vh_mat(ig_i, ig_j, ispin) = Vh_mat(ig_i, ig_j, ispin) + Vh_local_mat(io, jo)
              Vxc_mat(ig_i, ig_j, ispin) = Vxc_mat(ig_i, ig_j, ispin) + Vxc_local_mat(io, jo)
            end do
          end do
        end if
        
        ! HSE hybrid functional contribution (Plan A: density matrix method)
        if (yn_hse == 'y') then
          call add_exact_exchange_hse(dg_frag, system, dg_frag%H_mat(:, :, ispin), ifrag, ispin)
        end if
        
      end do
    end do
    
    ! Notes on implementation:
    ! 1. H_mat contains ONLY the potential-dependent terms
    !    - This is the time-independent part: H_0 = T + V_psl + V_H + V_xc
    !    - A·p and A²/2 added dynamically in calculate_time_derivative
    ! 2. Assumes basis functions φ_i are real (valid for most systems)
    ! 3. Uses simple rectangular integration (hvol weighting)
    ! 4. Parallelized over matrix elements (OpenMP collapse)
    ! 5. Fragment index mapping: ifrag - ifrag_start + 1 for local storage
    !
    ! For basis update detection (see below):
    ! - A temporary H_mat_full is constructed including A·p and A²
    ! - This captures complete Hamiltonian change H(t) vs H(t-Δt)
    ! - Ensures accurate detection in both weak and strong fields
    !
    ! Future improvements:
    ! - Use grid mapping (jxyz_tot) for proper fragment domain handling
    ! - Handle overlapping fragment regions (weighted average)
    ! - More accurate integration schemes (trapezoid, Simpson)
    ! - Exploit Hermiticity: H_ji = H_ij* (only compute upper triangle)
    ! - MPI communication for cross-fragment matrix elements

    ! Reduce the reconstructed Hamiltonian through fragment blocks to avoid
    ! rebuilding a full dense allreduce path during basis refresh.
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

    ! Per-term asymmetry diagnostics
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

    ! Hermiticity check: max |H - H^T| (H is real)
    if (.not. use_block_reconstruct) then
      block
        integer :: nmat_chk
        real(8) :: max_asym
        nmat_chk = dg_frag%n_mat_max
        do ispin = 1, system%nspin
          max_asym = 0.0d0
          !$omp parallel do collapse(2) private(io,jo) reduction(max:max_asym) schedule(static)
          do jo = 1, nmat_chk
            do io = 1, nmat_chk
              max_asym = max(max_asym, abs(dg_frag%H_mat(io, jo, ispin) - dg_frag%H_mat(jo, io, ispin)))
            end do
          end do
          !$omp end parallel do
          if (max_asym > 1.0d-8) then
            write(*,'(1x,a,i0,a,i0,a,es12.4)') "[WARN] H_mat not Hermitian: rank=", dg_frag%id, &
              " ispin=", ispin, " max|H-H^T|=", max_asym
          end if
        end do
      end block
    end if

    ! Enforce Hermiticity: H <- (H + H^T)/2
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
        integer :: nmat_chk, io, jo
        nmat_chk = dg_frag%n_mat_max
!$omp parallel do collapse(2) private(ispin,io,jo) schedule(static)
        do ispin = 1, system%nspin
          do jo = 1, nmat_chk
            do io = jo + 1, nmat_chk
              dg_frag%H_mat(io, jo, ispin) = 0.5d0 * (dg_frag%H_mat(io, jo, ispin) + dg_frag%H_mat(jo, io, ispin))
              dg_frag%H_mat(jo, io, ispin) = dg_frag%H_mat(io, jo, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end block
    end if

    ! Keep complex Hamiltonian view consistent with reconstructed H_mat.
    if ((.not. use_block_reconstruct) .and. allocated(dg_frag%H_mat_c)) then
      dg_frag%H_mat_c(:, :, :) = cmplx(dg_frag%H_mat(:, :, :), 0.0d0, kind=8)
    end if
    if (release_dense_h) then
      if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    end if
    if (allocated(phi_blk)) deallocate(phi_blk, weighted_phi, weighted_phi_vpsl, weighted_phi_vh, weighted_phi_vxc, &
                                       Vpsl_blk, Vh_blk, Vxc_blk, Vtot_blk)
    if (allocated(total_mat)) deallocate(total_mat, Vpsl_local_mat, Vh_local_mat, Vxc_local_mat)
    
  end subroutine reconstruct_hamiltonian_matrix
