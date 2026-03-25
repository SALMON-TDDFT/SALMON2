  subroutine reconstruct_hamiltonian_matrix(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use salmon_global, only: yn_hse
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_scalar),         intent(in)    :: Vh
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    type(s_scalar),         intent(in)    :: Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    
    integer :: ifrag, ispin, io, jo, ix, iy, iz, idir, i_local, ig_i, ig_j, nbf, nbf_raw
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3)
    real(8) :: hvol, A_squared
    real(8) :: integral_Vh, integral_Vxc, integral_Vpsl  ! Real integrals (real basis)
    integer :: is(3), ie(3)
    real(8), allocatable :: Vpsl_mat(:,:,:), Vh_mat(:,:,:), Vxc_mat(:,:,:)
    real(8) :: max_asym_vpsl, max_asym_vh, max_asym_vxc
    integer :: nmat_chk
    
    ! Check if real-space basis functions are available
    if (.not. dg_frag%has_real_space_basis) then
      return  ! Cannot reconstruct without basis functions
    end if
    
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
    
    ! Step 1: Reset to kinetic parts
    dg_frag%H_mat(:, :, :) = dg_frag%H_mat_kinetic(:, :, :)

    ! Allocate diagnostic matrices for potential terms
    nmat_chk = dg_frag%n_mat_max
    allocate(Vpsl_mat(nmat_chk, nmat_chk, system%nspin))
    allocate(Vh_mat(nmat_chk, nmat_chk, system%nspin))
    allocate(Vxc_mat(nmat_chk, nmat_chk, system%nspin))
    Vpsl_mat = 0.0d0
    Vh_mat = 0.0d0
    Vxc_mat = 0.0d0
    
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
      do ispin = 1, system%nspin
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
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
        !$omp parallel do private(io, jo, ig_i, ig_j, ix, iy, iz, integral_Vh, integral_Vxc, integral_Vpsl) collapse(2)
        do jo = 1, nbf
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            
            ! Calculate pseudopotential matrix element:
            ! V_psl,ij = ∫ φ_i*(r) V_psl(r) φ_j(r) dr
            integral_Vpsl = 0.0d0
            do lz = 1, ndom(3)
              gz = iorg(3) + lz - 1
              do ly = 1, ndom(2)
                gy = iorg(2) + ly - 1
                do lx = 1, ndom(1)
                  gx = iorg(1) + lx - 1
                  integral_Vpsl = integral_Vpsl + &
                    dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                    Vpsl%f(gx, gy, gz) * &
                    dg_frag%phi_frag(lx, ly, lz, jo, i_local) * &
                    hvol
                end do
              end do
            end do
            if (integral_Vpsl /= integral_Vpsl) then
              write(*,'(a,i0,a,i0,a,i0,a,i0)') "[NaN] Vpsl integral, rank=", dg_frag%id, " ifrag=", ifrag, " io=", io, " jo=", jo
              stop "NaN in Vpsl integral"
            end if
            
            ! Calculate Hartree potential matrix element:
            ! V_H,ij = ∫ φ_i*(r) V_H(r) φ_j(r) dr
            integral_Vh = 0.0d0
            do lz = 1, ndom(3)
              gz = iorg(3) + lz - 1
              do ly = 1, ndom(2)
                gy = iorg(2) + ly - 1
                do lx = 1, ndom(1)
                  gx = iorg(1) + lx - 1
                  ! Numerical integration: φ_i * V_H * φ_j * volume_element (real basis)
                  integral_Vh = integral_Vh + &
                    dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                    Vh%f(gx, gy, gz) * &
                    dg_frag%phi_frag(lx, ly, lz, jo, i_local) * &
                    hvol
                end do
              end do
            end do
            if (integral_Vh /= integral_Vh) then
              write(*,'(a,i0,a,i0,a,i0,a,i0)') "[NaN] Vh integral, rank=", dg_frag%id, " ifrag=", ifrag, " io=", io, " jo=", jo
              stop "NaN in Vh integral"
            end if
            
            ! Calculate XC potential matrix element:
            ! V_xc,ij = ∫ φ*_i(r) V_xc(r) φ_j(r) dr
            integral_Vxc = 0.0d0
            do lz = 1, ndom(3)
              gz = iorg(3) + lz - 1
              do ly = 1, ndom(2)
                gy = iorg(2) + ly - 1
                do lx = 1, ndom(1)
                  gx = iorg(1) + lx - 1
                  ! Numerical integration: φ_i * V_xc * φ_j * volume_element (real basis)
                  integral_Vxc = integral_Vxc + &
                    dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                    Vxc(ispin)%f(gx, gy, gz) * &
                    dg_frag%phi_frag(lx, ly, lz, jo, i_local) * &
                    hvol
                end do
              end do
            end do
            if (integral_Vxc /= integral_Vxc) then
              write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') "[NaN] Vxc integral, rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin, " io=", io, " jo=", jo
              stop "NaN in Vxc integral"
            end if
            
            ! Add potential contributions: V_psl + V_H + V_xc
            ! Note: A·p and A²/2 are NOT added here (added in time evolution)
            dg_frag%H_mat(ig_i, ig_j, ispin) = dg_frag%H_mat(ig_i, ig_j, ispin) + &
                                           real(integral_Vpsl + integral_Vh + integral_Vxc, kind=8)
            Vpsl_mat(ig_i, ig_j, ispin) = Vpsl_mat(ig_i, ig_j, ispin) + integral_Vpsl
            Vh_mat(ig_i, ig_j, ispin) = Vh_mat(ig_i, ig_j, ispin) + integral_Vh
            Vxc_mat(ig_i, ig_j, ispin) = Vxc_mat(ig_i, ig_j, ispin) + integral_Vxc
            
          end do
        end do
        !$omp end parallel do
        
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
    call sync_dense_matrix_to_blocks(dg_frag, dg_frag%H_mat, dg_frag%H_mat_blocks, dg_frag%H_block_map)
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-reconstruct", dg_frag%icomm)
    call sync_blocks_to_dense_matrix(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%H_mat)

    ! Per-term asymmetry diagnostics
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

    ! Hermiticity check: max |H - H^T| (H is real)
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

    ! Enforce Hermiticity: H <- (H + H^T)/2
    block
      integer :: nmat_chk, io, jo
      nmat_chk = dg_frag%n_mat_max
      do ispin = 1, system%nspin
        !$omp parallel do private(io,jo) schedule(static)
        do jo = 1, nmat_chk
          do io = jo + 1, nmat_chk
            dg_frag%H_mat(io, jo, ispin) = 0.5d0 * (dg_frag%H_mat(io, jo, ispin) + dg_frag%H_mat(jo, io, ispin))
            dg_frag%H_mat(jo, io, ispin) = dg_frag%H_mat(io, jo, ispin)
          end do
        end do
        !$omp end parallel do
      end do
    end block

    ! Keep complex Hamiltonian view consistent with reconstructed H_mat.
    if (allocated(dg_frag%H_mat_c)) then
      dg_frag%H_mat_c(:, :, :) = cmplx(dg_frag%H_mat(:, :, :), 0.0d0, kind=8)
    end if
    
  end subroutine reconstruct_hamiltonian_matrix
