!=======================================================================
  ! Check Hamiltonian change for assigned fragments and return update flag
  ! Each MPI rank checks its own fragments, then collective decision via Allreduce
  !=======================================================================
  function check_hamiltonian_change_fragments(dg_frag, H_mat_current) result(needs_update)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(in) :: H_mat_current(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin)
    logical :: needs_update
    
    integer :: i, j, ispin
    real(8) :: norm_sq_local, norm_sq_global, diff_re, diff_im
    real(8) :: kahan_c, term, y, t
    
    ! Step 1: Calculate Frobenius norm for this rank's fragments
    ! ||H_new - H_old||_F = sqrt(sum_ij |H_new_ij - H_old_ij|^2)
    norm_sq_local = 0.0d0
    kahan_c = 0.0d0
    do ispin = 1, dg_frag%nspin
      do j = 1, dg_frag%nstate_frag
        do i = 1, dg_frag%nstate_frag
          diff_re = real(H_mat_current(i, j, ispin)) - real(dg_frag%H_mat_old(i, j, ispin))
          diff_im = aimag(H_mat_current(i, j, ispin) - dg_frag%H_mat_old(i, j, ispin))
          term = diff_re**2 + diff_im**2
          y = term - kahan_c
          t = norm_sq_local + y
          kahan_c = (t - norm_sq_local) - y
          norm_sq_local = t
        end do
      end do
    end do
    
    ! Step 2: Sum across all fragments (Allreduce)
    call comm_summation(norm_sq_local, norm_sq_global, dg_frag%icomm)
    dg_frag%hamiltonian_change_norm = sqrt(norm_sq_global)
    
    ! Step 3: Global decision (MPI-size independent)
    ! Use global Frobenius norm against threshold to avoid rank-count dependence.
    needs_update = (dg_frag%hamiltonian_change_norm > dg_frag%basis_update_threshold)
    
    ! Store current Hamiltonian as old for next iteration
    dg_frag%H_mat_old(:, :, :) = H_mat_current(:, :, :)
    
  end function check_hamiltonian_change_fragments

  !=======================================================================
  ! Trigger adaptive basis update via DC-LCFO recalculation
  ! 
  ! Two methods available:
  !   1. DC-LCFO CG solver (RECOMMENDED): Expands basis space
  !   2. Simple diagonalization (FALLBACK): Only rotates basis
  !=======================================================================
  subroutine trigger_basis_update(dg_frag, system, info, itt_update, lg, mg, stencil, srg, &
                                  Vh, Vxc, Vpsl, pp, ppg, tpsi, Ac_tot)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use parallelization, only: nproc_id_global, nproc_size_global
    use salmon_global, only: yn_dc_cg_basis_update
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(inout) :: system
    type(s_parallel_info),  intent(in) :: info
    integer, intent(in) :: itt_update
    type(s_rgrid), intent(in) :: lg, mg
    type(s_stencil), intent(in) :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_scalar), intent(in) :: Vh
    type(s_scalar), intent(in) :: Vxc(system%nspin)
    type(s_scalar), intent(in) :: Vpsl
    type(s_pp_info), intent(in) :: pp
    type(s_pp_grid), intent(in) :: ppg
    type(s_orbital), intent(inout) :: tpsi
    real(8), intent(in) :: Ac_tot(3)
    
    logical :: use_dc_cg_method
    logical :: is_global_root
    real(8) :: e_low, e_high
    integer :: nesp

    is_global_root = (nproc_id_global == 0)
    if (is_global_root) then
      write(*,*)
      write(*,*) "!!! ADAPTIVE BASIS UPDATE TRIGGERED !!!"
      write(*,'(1x,a,f10.6,a)') "  Hamiltonian change: ", &
                                dg_frag%hamiltonian_change_norm, " a.u."
      write(*,'(1x,a,f10.6,a)') "  Threshold: ", &
                                dg_frag%basis_update_threshold, " a.u."
      write(*,*)
    end if

    ! Record update step and increment counter before metric evaluation in projection routines.
    dg_frag%last_basis_update_step = itt_update
    dg_frag%nbasis_update_count = dg_frag%nbasis_update_count + 1
    
    ! Decide which method to use based on user input parameter
    ! yn_dc_cg_basis_update = 'y' : Use DC-CG method (recommended)
    ! yn_dc_cg_basis_update = 'n' : Use simple diagonalization (fallback)
    use_dc_cg_method = (yn_dc_cg_basis_update == 'y')
    
    if (use_dc_cg_method) then
      ! ============================================================
      ! Method 1: DC-LCFO CG solver (RECOMMENDED)
      ! - Solves KS equation with current potentials
      ! - Expands basis space to capture new physics
      ! - Uses ppg (pseudopotential grid) for nonlocal PP
      ! ============================================================
      if (is_global_root) then
        write(*,*) "  Using DC-LCFO CG method (basis expansion)"
      end if
      
      ! Pseudopotential grid now passed from main RT loop (initialized in initialization_rt)
      ! ppg contains mps (number of PP atoms), jxyz (grid-to-atom mapping),
      ! uv (nonlocal projectors), etc.
      
      call update_basis_via_dc_cg(dg_frag, system, info, lg, mg, stencil, srg, ppg, Vh, Vxc, Vpsl, Ac_tot)
      
    else
      ! ============================================================
      ! Method 2: Simple diagonalization (FALLBACK)
      ! - Only rotates basis in existing space
      ! - No basis expansion (limitation)
      ! - Used when DC-CG not available or disabled
      ! ============================================================
      if (is_global_root) then
        write(*,*) "  Using simple diagonalization (fallback)"
        write(*,*) "  WARNING: Basis space not expanded"
        write(*,*) "           For proper updates, enable DC-CG method"
        write(*,*)
      end if
      
      call diagonalize_and_update_basis(dg_frag, system)
      
    end if
    ! ============================================================
    
    e_low = 0.0d0
    e_high = 0.0d0
    if (allocated(dg_frag%esp)) then
      nesp = max(1, min(dg_frag%nstate_tot, size(dg_frag%esp,1)))
      e_low = minval(dg_frag%esp(1:nesp, 1:dg_frag%nspin))
      e_high = maxval(dg_frag%esp(1:nesp, 1:dg_frag%nspin))
    end if

    if (is_global_root) then
      write(*,*) "  Basis updated via diagonalization"
      write(*,'(1x,a,i0)') "  Basis update count: ", dg_frag%nbasis_update_count
      write(*,'(1x,a,f12.6,a)') "  Lowest eigenvalue: ", e_low, " a.u."
      write(*,'(1x,a,f12.6,a)') "  Highest eigenvalue: ", e_high, " a.u."
      write(*,*) "  Continuing time evolution..."
      write(*,*)
    end if
    
  end subroutine trigger_basis_update

  !=======================================================================
  ! Diagonalize mixed basis (fragment + orthogonalized plane waves)
  !=======================================================================
  subroutine diagonalize_fragment_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    use rt_dg_plane_wave, only: compute_fragment_pw_overlap, compute_fragment_pw_hamiltonian, &
                                build_mixed_hamiltonian, assemble_mixed_hamiltonian_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)
    
    integer :: ispin, n_total, n_frag, n_pw, lda, lwork, info, i, j
    complex(8), allocatable :: H_work(:,:), work(:)
    real(8), allocatable :: eigenvalues_tmp(:), rwork(:)
    complex(8), allocatable :: coef_mixed(:,:,:), coef_new(:,:,:)
    complex(8), allocatable :: S_frag_pw(:,:,:)  ! Complex overlap matrix
    complex(8), allocatable :: H_frag_pw(:,:,:)  ! Hamiltonian coupling matrix
    
    n_frag = size(dg_frag%coef, 1)
    n_pw = dg_frag%n_plane_waves
    n_total = n_frag + n_pw
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Diagonalizing mixed basis (Fragment + PW) ==="
      write(*,'(1x,a,i0)') "Fragment basis size: ", n_frag
      write(*,'(1x,a,i0)') "Plane wave basis size: ", n_pw
      write(*,'(1x,a,i0)') "Total mixed basis size: ", n_total
    end if
    
    ! Compute overlap matrix between fragment basis and plane waves
    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
    ! Build mixed Hamiltonian matrix
    ! Note: This is simplified - ideally should use L\u00f6wdin orthogonalization
    call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)
    
    ! Diagonalize for each spin
    allocate(coef_mixed(n_total, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(coef_new(n_frag, dg_frag%nstate_tot, dg_frag%nspin))
    
    do ispin = 1, dg_frag%nspin
      lda = n_total
      
      allocate(H_work(n_total, n_total))
      allocate(eigenvalues_tmp(n_total))
      
      ! Assemble mixed Hamiltonian on demand from FF/FP/PP data.
      call assemble_mixed_hamiltonian_dense(dg_frag, ispin, H_frag_pw, H_work)
      
      ! Query optimal workspace size
      lwork = -1
      allocate(work(1), rwork(max(1, 3*n_total-2)))
      call ZHEEV('V', 'U', n_total, H_work, lda, eigenvalues_tmp, work, lwork, rwork, info)
      lwork = int(real(work(1), kind=8)) + 1
      deallocate(work)
      allocate(work(lwork))
      
      ! Diagonalize mixed Hamiltonian
      call ZHEEV('V', 'U', n_total, H_work, lda, eigenvalues_tmp, work, lwork, rwork, info)
      
      if (info /= 0) then
        write(*,*) "ERROR: Mixed basis diagonalization failed, info=", info
        stop
      end if
      
      ! Store eigenvalues (lowest nstate_tot eigenvalues)
      do i = 1, min(dg_frag%nstate_tot, n_total)
        dg_frag%esp(i, ispin) = eigenvalues_tmp(i)
      end do
      
      ! Eigenvectors are in H_work columns
      ! Transform coefficients from old basis to new eigenbasis
      ! coef_new = U^T * coef_old (project old state onto new eigenstates)
      
      ! For simplicity, initialize with lowest eigenstates
      ! In real implementation, should preserve occupied states
      coef_mixed(:, :, ispin) = (0.0d0, 0.0d0)
      do i = 1, min(dg_frag%nstate_tot, n_total)
        ! Mixed basis eigenvectors (columns of H_work)
        ! Fragment part: H_work(1:n_frag, i)
        ! PW part: H_work(n_frag+1:n_total, i)
        do j = 1, n_frag
          coef_mixed(j, i, ispin) = H_work(j, i)
        end do
      end do
      
      ! Extract fragment part for standard coefficient array
      dg_frag%coef(1:n_frag, 1:dg_frag%nstate_tot, ispin) = &
           coef_mixed(1:n_frag, 1:dg_frag%nstate_tot, ispin)
      
      ! Store plane wave coefficients separately
      if (n_pw > 0) then
        do i = 1, min(dg_frag%nstate_tot, n_total)
          do j = 1, n_pw
            dg_frag%coef_pw(j, i, ispin) = H_work(n_frag+j, i)
          end do
        end do
      end if
      
      deallocate(H_work, eigenvalues_tmp, work, rwork)
    end do
    
    deallocate(coef_mixed, coef_new, S_frag_pw, H_frag_pw)
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "Mixed basis diagonalization complete"
      write(*,'(1x,a,f12.6,a)') "Lowest eigenvalue: ", dg_frag%esp(1, 1), " a.u."
      if (dg_frag%nstate_tot > 1) then
        write(*,'(1x,a,f12.6,a)') "Highest occupied energy: ", &
             dg_frag%esp(min(system%no, dg_frag%nstate_tot), 1), " a.u."
      end if
      write(*,*)
    end if
    
  end subroutine diagonalize_fragment_basis

  !=======================================================================
  ! Update basis via DC-LCFO external calculation
  !
  ! Workflow:
  !   1. Save current potentials (Vh, Vxc, Vpsl) to files
  !   2. Launch external DC-LCFO calculation
  !   3. Read new basis functions
  !   4. Calculate overlap matrix and project wave functions
  !=======================================================================
  subroutine update_basis_via_dc_cg(dg_frag, system, info, lg, mg, stencil, &
                                    srg, ppg, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use parallelization, only: nproc_id_global
    use communication, only: comm_sync_all
    use rt_dg_fragment_parallel, only: setup_fragment_system, finalize_fragment_parallel
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid), intent(in) :: lg, mg
    type(s_stencil), intent(in) :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_pp_grid), intent(in) :: ppg
    type(s_scalar), intent(in) :: Vh, Vxc(system%nspin), Vpsl
    real(8), intent(in) :: Ac_tot(3)

    ! Local variables
    real(8), allocatable :: phi_frag_old(:,:,:,:,:)  ! Old basis for overlap
    complex(8), allocatable :: coef_old(:,:,:), coef_pw_old(:,:,:)
    logical :: basis_functions_changed, overlap_is_valid
    logical :: is_global_root
    integer :: ispin

    is_global_root = (nproc_id_global == 0)
    if (is_global_root) then
      write(*,*)
      write(*,*) "=========================================="
      write(*,*) "UPDATE BASIS (MEMORY-OPTIMIZED)"
      write(*,*) "=========================================="
    end if

    ! ========================================================
    ! Step 1: Save old basis for overlap calculation
    ! ========================================================
    if (allocated(dg_frag%phi_frag)) then
      allocate(phi_frag_old(size(dg_frag%phi_frag,1), &
                           size(dg_frag%phi_frag,2), &
                           size(dg_frag%phi_frag,3), &
                           size(dg_frag%phi_frag,4), &
                           size(dg_frag%phi_frag,5)))
      phi_frag_old = dg_frag%phi_frag
    end if

    allocate(coef_old(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
    coef_old = dg_frag%coef

    if (dg_frag%use_plane_wave_basis) then
      if (allocated(dg_frag%coef_pw)) then
        allocate(coef_pw_old(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
        coef_pw_old = dg_frag%coef_pw
      end if
    end if

    if (is_global_root) then
      write(*,*) "  [1/3] Old basis saved to memory"
    end if

    ! ========================================================
    ! Step 2: Compute updated states from current potentials
    ! ========================================================
    if (is_global_root) then
      write(*,*) "  [2/3] Computing new basis from current potentials..."
    end if

    basis_functions_changed = .false.

    ! ========================================================
    ! Step 3: Basis update and state transfer (no SCF loop)
    ! ========================================================
    if (dg_frag%use_plane_wave_basis) then
      ! Normal RT mixed-basis updates must not require a dense global mixed EVP.
      ! In the current SOI mixed path the real-space fragment basis is not
      ! refreshed here, so the propagated state transfer is the identity in the
      ! existing fragment+PW representation.
      dg_frag%coef(:, :, :) = coef_old(:, :, :)
      dg_frag%coef_new(:, :, :) = coef_old(:, :, :)
      if (allocated(coef_pw_old) .and. allocated(dg_frag%coef_pw)) then
        dg_frag%coef_pw(:, :, :) = coef_pw_old(:, :, :)
      end if
      call zero_nonowned_coefficients(dg_frag)

      if (is_global_root) then
        write(*,'(1x,a)') "  Mixed basis update: skipped dense global re-diagonalization"
        write(*,'(1x,a)') "  Propagated fragment/PW coefficients carried forward in-place"
      end if
    else
      ! Reproduce the intended GS-DC partial workflow in RT:
      !  (1) fragment-wise DC-CG (no SCF),
      !  (2) one-shot full diagonalization in updated basis,
      !  (3) project saved propagated coefficients onto updated basis.
      call update_fragment_basis_via_cg(dg_frag, system, info, mg, stencil, srg, ppg, Vh, Vxc, Vpsl, basis_functions_changed)

      call comm_sync_all(dg_frag%icomm)
      call diagonalize_full_system_dg(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, ppg, basis_functions_changed)

      overlap_is_valid = .true.
      call validate_overlap_matrix(dg_frag, overlap_is_valid)

      if (.not. overlap_is_valid) then
        if (is_global_root) then
          write(*,'(1x,a)') "  [WARN] Updated basis rejected: overlap matrix is ill-conditioned"
          write(*,'(1x,a)') "  [WARN] Restoring previous basis and coefficients"
        end if
        if (allocated(phi_frag_old)) dg_frag%phi_frag = phi_frag_old
        dg_frag%coef = coef_old
        dg_frag%coef_new = coef_old
        call zero_nonowned_coefficients(dg_frag)
	        if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
	        if (allocated(dg_frag%momentum_mat_c)) deallocate(dg_frag%momentum_mat_c)
	        if (allocated(dg_frag%momentum_blocks)) deallocate(dg_frag%momentum_blocks)
	        if (allocated(dg_frag%momentum_blocks_c)) deallocate(dg_frag%momentum_blocks_c)
	        if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
        if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
        if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
        if (allocated(dg_frag%S_mat_c)) deallocate(dg_frag%S_mat_c)
        if (allocated(dg_frag%S_mat_prop_c)) deallocate(dg_frag%S_mat_prop_c)
        dg_frag%has_nl_cache = .false.
        call calculate_momentum_matrix(dg_frag, system, mg, stencil)
        call calculate_overlap_matrix(dg_frag, system, mg)
      else if (allocated(phi_frag_old)) then
        call project_coefficients_fragmentwise(dg_frag, system, phi_frag_old, coef_old)
      end if

      if (.not. allocated(phi_frag_old)) then
        dg_frag%coef = coef_old
        dg_frag%coef_new = coef_old
        call zero_nonowned_coefficients(dg_frag)
        if (is_global_root) then
          write(*,'(1x,a)') "  No local old basis snapshot: restore coefficients and skip projection"
        end if
      end if
    end if

    if (is_global_root) then
      write(*,*) "  [3/3] Wave functions projected to new basis"
      write(*,*) "=========================================="
      write(*,*) "BASIS UPDATE COMPLETE (NO FILE I/O)"
      write(*,*) "=========================================="
      write(*,*)
    end if

    ! Cleanup
    if (allocated(phi_frag_old)) deallocate(phi_frag_old)
    if (allocated(coef_old)) deallocate(coef_old)
    if (allocated(coef_pw_old)) deallocate(coef_pw_old)

  end subroutine update_basis_via_dc_cg

  !=======================================================================
  ! Project old propagated mixed-state coefficients (fragment+PW) onto the
  ! updated mixed basis obtained by diagonalize_mixed_basis_pw.
  !=======================================================================
  subroutine project_coefficients_mixed_state(dg_frag, coef_old, coef_pw_old)
    use rt_dg_fragment_ops, only: apply_matrix_blocks_batch, zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(in) :: coef_old(:,:,:)
    complex(8), intent(in), optional :: coef_pw_old(:,:,:)

    integer :: ispin, n_frag, n_pw, n_tot, nst, ipw_local
    complex(8), allocatable :: U_new(:,:), C_old(:,:), C_new(:,:)
    complex(8), allocatable :: Sm(:,:), tmp(:,:), A(:,:)
    complex(8), parameter :: zone = (1.0d0, 0.0d0), zzero = (0.0d0, 0.0d0)

    n_frag = size(dg_frag%coef, 1)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    nst = dg_frag%nstate_tot
    if (n_tot <= 0 .or. nst <= 0) return

    allocate(U_new(n_tot, nst), C_old(n_tot, nst), C_new(n_tot, nst))
    allocate(Sm(n_tot, n_tot), tmp(n_tot, nst), A(nst, nst))

    do ispin = 1, dg_frag%nspin
      U_new(:, :) = zzero
      C_old(:, :) = zzero
      C_new(:, :) = zzero
      Sm(:, :) = zzero

      U_new(1:n_frag, :) = dg_frag%coef(1:n_frag, 1:nst, ispin)
      if (n_pw > 0) U_new(n_frag+1:n_tot, :) = dg_frag%coef_pw(1:n_pw, 1:nst, ispin)

      C_old(1:n_frag, :) = coef_old(1:n_frag, 1:nst, ispin)
      if (present(coef_pw_old) .and. n_pw > 0) C_old(n_frag+1:n_tot, :) = coef_pw_old(1:n_pw, 1:nst, ispin)

      if (allocated(dg_frag%S_mat_prop_c)) then
        Sm(1:n_frag, 1:n_frag) = dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin)
      else
        Sm(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
        if (n_pw > 0) then
          Sm(n_frag+1:n_tot, n_frag+1:n_tot) = zzero
          do ipw_local = 1, n_pw
            Sm(n_frag+ipw_local, n_frag+ipw_local) = zone
          end do
        end if
      end if
      if (n_pw > 0) then
        Sm(n_frag+1:n_tot, n_frag+1:n_tot) = zzero
        do ipw_local = 1, n_pw
          Sm(n_frag+ipw_local, n_frag+ipw_local) = zone
        end do
      end if

      ! A = U_new^† S C_old, then C_new = U_new A
      tmp(:, :) = zzero
      if (allocated(dg_frag%S_mat_prop_blocks)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_prop_blocks, ispin, C_old(1:n_frag, :), tmp(1:n_frag, :))
      else if (allocated(dg_frag%S_mat_blocks)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, C_old(1:n_frag, :), tmp(1:n_frag, :))
      else
        call zgemm('N', 'N', n_frag, nst, n_frag, zone, Sm(1:n_frag, 1:n_frag), n_tot, C_old(1:n_frag, :), n_tot, zzero, tmp(1:n_frag, :), n_tot)
      end if
      if (n_pw > 0) then
        call zgemm('N', 'N', n_frag, nst, n_pw, zone, Sm(1:n_frag, n_frag+1:n_tot), n_tot, &
                   C_old(n_frag+1:n_tot, :), n_tot, zone, tmp(1:n_frag, :), n_tot)
        tmp(n_frag+1:n_tot, :) = C_old(n_frag+1:n_tot, :)
      end if
      call zgemm('C', 'N', nst, nst, n_tot, zone, U_new, n_tot, tmp, n_tot, zzero, A, nst)
      call zgemm('N', 'N', n_tot, nst, nst, zone, U_new, n_tot, A, nst, zzero, C_new, n_tot)

      dg_frag%coef(1:n_frag, 1:nst, ispin) = C_new(1:n_frag, 1:nst)
      if (n_pw > 0) dg_frag%coef_pw(1:n_pw, 1:nst, ispin) = C_new(n_frag+1:n_tot, 1:nst)
    end do
    call zero_nonowned_coefficients(dg_frag)

    deallocate(U_new, C_old, C_new, Sm, tmp, A)
  end subroutine project_coefficients_mixed_state

  !=======================================================================
  ! Project coefficients fragment by fragment:
  !   c'_j = sum_i <phi'_j|phi_i> c_i  within each fragment subspace
  !=======================================================================
  subroutine project_coefficients_fragmentwise(dg_frag, system, phi_frag_old, coef_old)
    use structures
    use communication, only: comm_summation
    use parallelization, only: nproc_id_global
    use salmon_global, only: dt
    use phys_constants, only: au_fs
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    real(8), intent(in) :: phi_frag_old(:,:,:,:,:)
    complex(8), intent(in) :: coef_old(:,:,:)

    integer :: ifrag, i_local, ispin, io, i, j, ix, iy, iz
    integer :: iidx, jidx, iloc, jloc, nb
    integer :: ndom(3)
    real(8) :: hvol
    real(8), allocatable :: ovlp(:,:), ovlp_orig(:,:), U(:,:), VT(:,:), S(:), work(:)
    complex(8), allocatable :: cnew(:)
    real(8) :: rot_local_sum, rot_global_sum, rot_norm
    real(8) :: time_fs
    integer :: count_local, count_global
    integer :: info, lwork, irow, icol
    real(8) :: max_abs, sign_val
    real(8) :: work_query(1)
    external :: dgesvd

    hvol = system%hvol
    dg_frag%coef(:, :, :) = (0.0d0, 0.0d0)
    rot_local_sum = 0.0d0
    count_local = 0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)

      do ispin = 1, dg_frag%nspin
        nb = dg_frag%n_basis(ifrag, ispin)
        if (nb <= 0) cycle
        allocate(ovlp(nb, nb))
        allocate(ovlp_orig(nb, nb))
        allocate(U(nb, nb))
        allocate(VT(nb, nb))
        allocate(S(nb))
        allocate(cnew(nb))
        ovlp = 0.0d0

        do j = 1, nb
          do i = 1, nb
            do iz = 1, ndom(3)
              do iy = 1, ndom(2)
                do ix = 1, ndom(1)
                  ovlp(j,i) = ovlp(j,i) + &
                    dg_frag%phi_frag(ix,iy,iz,j,i_local) * &
                    phi_frag_old(ix,iy,iz,i,i_local) * hvol
                end do
              end do
            end do
          end do
        end do
        ovlp_orig(:,:) = ovlp(:,:)

        ! Sign alignment + Procrustes rotation (Q=U*VT) to suppress gauge drift
        do irow = 1, nb
          max_abs = 0.0d0
          sign_val = 1.0d0
          do icol = 1, nb
            if (abs(ovlp(irow,icol)) > max_abs) then
              max_abs = abs(ovlp(irow,icol))
              if (ovlp(irow,icol) < 0.0d0) sign_val = -1.0d0
            end if
          end do
          if (sign_val < 0.0d0) ovlp(irow,:) = -ovlp(irow,:)
        end do

        lwork = -1
        call dgesvd('A', 'A', nb, nb, ovlp, nb, S, U, nb, VT, nb, work_query, lwork, info)
        if (info == 0) then
          lwork = max(1, int(work_query(1)))
          allocate(work(lwork))
          call dgesvd('A', 'A', nb, nb, ovlp, nb, S, U, nb, VT, nb, work, lwork, info)
          deallocate(work)
          if (info == 0) then
            ovlp(:,:) = matmul(U, VT)
          else
            ovlp(:,:) = ovlp_orig(:,:)
          end if
        else
          ovlp(:,:) = ovlp_orig(:,:)
        end if

        rot_norm = 0.0d0
        do j = 1, nb
          do i = 1, nb
            if (i == j) then
              rot_norm = rot_norm + (ovlp(j,i) - 1.0d0)**2
            else
              rot_norm = rot_norm + ovlp(j,i)**2
            end if
          end do
        end do
        rot_norm = sqrt(rot_norm / dble(max(1, nb)))
        rot_local_sum = rot_local_sum + rot_norm
        count_local = count_local + 1

        do io = 1, dg_frag%nstate_tot
          cnew(:) = (0.0d0, 0.0d0)
          do j = 1, nb
            do i = 1, nb
              iidx = dg_frag%index_basis(i, ifrag, ispin)
              if (iidx < 1 .or. iidx > dg_frag%n_mat_max) cycle
              if (.not. allocated(dg_frag%coef_global_to_local)) cycle
              iloc = dg_frag%coef_global_to_local(iidx, ispin)
              if (iloc < 1 .or. iloc > size(coef_old,1)) cycle
              cnew(j) = cnew(j) + ovlp(j,i) * coef_old(iloc, io, ispin)
            end do
          end do
          do j = 1, nb
            jidx = dg_frag%index_basis(j, ifrag, ispin)
            if (jidx < 1 .or. jidx > dg_frag%n_mat_max) cycle
            if (.not. allocated(dg_frag%coef_global_to_local)) cycle
            jloc = dg_frag%coef_global_to_local(jidx, ispin)
            if (jloc < 1 .or. jloc > size(dg_frag%coef,1)) cycle
            dg_frag%coef(jloc, io, ispin) = cnew(j)
          end do
        end do

        deallocate(cnew, S, VT, U, ovlp_orig, ovlp)
      end do
    end do

    dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    call zero_nonowned_coefficients(dg_frag)
    call stabilize_coeff_unitarity(dg_frag, dg_frag%last_basis_update_step)
    call comm_summation(rot_local_sum, rot_global_sum, dg_frag%icomm)
    call comm_summation(count_local, count_global, dg_frag%icomm)

    if (nproc_id_global == 0) then
      if (count_global > 0) then
        time_fs = dble(dg_frag%last_basis_update_step) * dt * au_fs
        write(*,'(1x,a,i0,a,i0,a,f10.6,a,1pe12.4)') &
          "[BASIS-ROT] update=", dg_frag%nbasis_update_count, &
          " step=", dg_frag%last_basis_update_step, " time_fs=", time_fs, &
          " mean_frob=", rot_global_sum / dble(count_global)
      end if
      write(*,'(1x,a)') "  Fragment-wise projection complete"
    end if

  end subroutine project_coefficients_fragmentwise

  !=======================================================================
  ! Validate overlap matrix quality after basis update.
  ! Reject nearly singular/indefinite S to avoid RT instability.
  !=======================================================================
  subroutine validate_overlap_matrix(dg_frag, is_valid)
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    logical, intent(out) :: is_valid

    integer :: ispin, n
    integer :: valid_local, valid_global
    real(8) :: eval_min, eval_max, cond_est, herm_err, sigma
    real(8), parameter :: eval_floor = 1.0d-10
    real(8), parameter :: cond_max = 1.0d12
    real(8), parameter :: herm_tol = 1.0d-8
    logical :: ok_local

    is_valid = .true.
    if (.not. allocated(dg_frag%S_mat) .and. .not. allocated(dg_frag%S_mat_c) .and. &
        .not. allocated(dg_frag%S_mat_blocks)) return
    if (.not. allocated(dg_frag%S_mat) .and. .not. allocated(dg_frag%S_mat_c)) return

    do ispin = 1, dg_frag%nspin
      n = dg_frag%n_mat(ispin)
      if (n <= 0) cycle
      if (allocated(dg_frag%S_mat_c)) then
        call estimate_overlap_metrics_complex(ispin, n, &
          lambda_min_est=eval_min, lambda_max_est=eval_max, herm_err=herm_err, sigma=sigma)
      else
        call estimate_overlap_metrics_real(ispin, n, &
          lambda_min_est=eval_min, lambda_max_est=eval_max, herm_err=herm_err, sigma=sigma)
      end if
      if (eval_max <= 0.0d0 .or. sigma <= 0.0d0) then
        is_valid = .false.
      else
        cond_est = eval_max / max(eval_floor, eval_min)
        ok_local = (eval_min > eval_floor) .and. (cond_est < cond_max) .and. &
          (herm_err < herm_tol)
        if (.not. ok_local) is_valid = .false.
      end if
      if (.not. is_valid) exit
    end do

    valid_local = 0
    if (is_valid) valid_local = 1
    call comm_summation(valid_local, valid_global, dg_frag%icomm)
    is_valid = (valid_global == dg_frag%isize)

    if (.not. is_valid .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "[WARN] Overlap validation failed after basis update"
    end if

  contains

    subroutine estimate_overlap_metrics_real(ispin, n, lambda_min_est, lambda_max_est, herm_err, sigma)
      integer, intent(in) :: ispin, n
      real(8), intent(out) :: lambda_min_est, lambda_max_est, herm_err, sigma
      real(8), allocatable :: x(:), y(:)
      real(8) :: norm_y, mu_shift
      integer :: i, j, iter
      integer, parameter :: max_iter = 12

      lambda_min_est = 0.0d0
      lambda_max_est = 0.0d0
      herm_err = 0.0d0
      sigma = 0.0d0

      do i = 1, n
        sigma = max(sigma, dg_frag%S_mat(i, i, ispin) + &
          sum(abs(dg_frag%S_mat(i, 1:n, ispin))) - abs(dg_frag%S_mat(i, i, ispin)))
        do j = i + 1, n
          herm_err = max(herm_err, abs(dg_frag%S_mat(i, j, ispin) - dg_frag%S_mat(j, i, ispin)))
        end do
      end do

      allocate(x(n), y(n))
      x = 1.0d0 / sqrt(dble(n))
      do iter = 1, max_iter
        y = matmul(dg_frag%S_mat(1:n, 1:n, ispin), x)
        norm_y = sqrt(max(dot_product(y, y), 0.0d0))
        if (norm_y <= tiny(1.0d0)) exit
        x = y / norm_y
      end do
      y = matmul(dg_frag%S_mat(1:n, 1:n, ispin), x)
      lambda_max_est = max(0.0d0, dot_product(x, y))

      x = 1.0d0 / sqrt(dble(n))
      do iter = 1, max_iter
        y = sigma * x - matmul(dg_frag%S_mat(1:n, 1:n, ispin), x)
        norm_y = sqrt(max(dot_product(y, y), 0.0d0))
        if (norm_y <= tiny(1.0d0)) exit
        x = y / norm_y
      end do
      y = sigma * x - matmul(dg_frag%S_mat(1:n, 1:n, ispin), x)
      mu_shift = max(0.0d0, dot_product(x, y))
      lambda_min_est = max(0.0d0, sigma - mu_shift)
      deallocate(x, y)
    end subroutine estimate_overlap_metrics_real

    subroutine estimate_overlap_metrics_complex(ispin, n, lambda_min_est, lambda_max_est, herm_err, sigma)
      integer, intent(in) :: ispin, n
      real(8), intent(out) :: lambda_min_est, lambda_max_est, herm_err, sigma
      complex(8), allocatable :: x(:), y(:)
      real(8) :: norm_y, mu_shift
      integer :: i, j, iter
      integer, parameter :: max_iter = 12

      lambda_min_est = 0.0d0
      lambda_max_est = 0.0d0
      herm_err = 0.0d0
      sigma = 0.0d0

      do i = 1, n
        sigma = max(sigma, real(dg_frag%S_mat_c(i, i, ispin), kind=8) + &
          sum(abs(dg_frag%S_mat_c(i, 1:n, ispin))) - abs(dg_frag%S_mat_c(i, i, ispin)))
        do j = i + 1, n
          herm_err = max(herm_err, abs(dg_frag%S_mat_c(i, j, ispin) - &
            conjg(dg_frag%S_mat_c(j, i, ispin))))
        end do
      end do

      allocate(x(n), y(n))
      x = cmplx(1.0d0 / sqrt(dble(n)), 0.0d0, kind=8)
      do iter = 1, max_iter
        y = matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), x)
        norm_y = sqrt(max(real(dot_product(y, y), kind=8), 0.0d0))
        if (norm_y <= tiny(1.0d0)) exit
        x = y / norm_y
      end do
      y = matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), x)
      lambda_max_est = max(0.0d0, real(dot_product(x, y), kind=8))

      x = cmplx(1.0d0 / sqrt(dble(n)), 0.0d0, kind=8)
      do iter = 1, max_iter
        y = sigma * x - matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), x)
        norm_y = sqrt(max(real(dot_product(y, y), kind=8), 0.0d0))
        if (norm_y <= tiny(1.0d0)) exit
        x = y / norm_y
      end do
      y = sigma * x - matmul(dg_frag%S_mat_c(1:n, 1:n, ispin), x)
      mu_shift = max(0.0d0, real(dot_product(x, y), kind=8))
      lambda_min_est = max(0.0d0, sigma - mu_shift)
      deallocate(x, y)
    end subroutine estimate_overlap_metrics_complex
  end subroutine validate_overlap_matrix


  !=======================================================================
  ! Full-system diagonalization in DG basis (DC-LCFO style)
  !=======================================================================
  subroutine diagonalize_full_system_dg(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, ppg, rebuild_rt_operators)
    use structures
    use communication, only: comm_summation, comm_bcast, comm_is_root, comm_isend, comm_irecv, comm_wait_all
    use salmon_global, only: yn_eigenexa
    use eigen_subdiag_sub, only: eigen_dsyev
#ifdef USE_EIGENEXA
    use eigen_libs_mod
#endif
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg, mg
    type(s_stencil), intent(in) :: stencil
    type(s_scalar), intent(in) :: Vh
    type(s_scalar), intent(in) :: Vxc(system%nspin)
    type(s_scalar), intent(in) :: Vpsl
    type(s_pp_grid), intent(in) :: ppg
    logical, intent(in) :: rebuild_rt_operators

    integer :: ifrag, i_local, ispin, io, jo, ix, iy, iz, i_halo, i, j
    integer :: jloc
    integer :: n, ifrag_count, n_basis_local, n_basis_halo, jfrag
    integer :: is(3), ie(3), l(3), d(3), iorg(3), ndom(3), lx, ly, lz, gx, gy, gz
    integer :: halo_send_idx(3), halo_recv_idx(3)
    real(8) :: hvol
    complex(8) :: integral
    complex(8), allocatable :: T_phi(:,:,:)
    complex(8), allocatable :: H_phi(:,:,:)
    real(8), allocatable :: V_total(:,:,:)
    real(8), allocatable :: mat_H_local(:,:,:,:)
    real(8), allocatable :: halo_H_local(:,:,:,:,:)
    complex(8), allocatable :: coef_local(:,:,:,:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "DG full-system diagonalization skipped: index_basis not available"
      end if
      return
    end if
    
    ! Recalculate RT operators only when requested or when not yet available.
    if (rebuild_rt_operators .or. .not. allocated(dg_frag%momentum_blocks) .or. &
        .not. allocated(dg_frag%S_mat_blocks)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  Recalculating momentum matrix for updated basis..."
      end if
	      if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
	      if (allocated(dg_frag%momentum_mat_c)) deallocate(dg_frag%momentum_mat_c)
	      if (allocated(dg_frag%momentum_blocks)) deallocate(dg_frag%momentum_blocks)
	      if (allocated(dg_frag%momentum_blocks_c)) deallocate(dg_frag%momentum_blocks_c)
	      if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
      if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
      if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
      if (allocated(dg_frag%S_mat_c)) deallocate(dg_frag%S_mat_c)
      if (allocated(dg_frag%S_mat_prop_c)) deallocate(dg_frag%S_mat_prop_c)
      call calculate_momentum_matrix(dg_frag, system, mg, stencil)
      call calculate_overlap_matrix(dg_frag, system, mg)
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  Momentum matrix recalculated successfully"
        write(*,*) "  Overlap matrix recalculated successfully"
      end if
    end if

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    allocate(mat_H_local(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin, ifrag_count))
    allocate(halo_H_local(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin, dg_frag%n_halo, ifrag_count))
    mat_H_local = 0.0d0
    halo_H_local = 0.0d0

    is = mg%is
    ie = mg%ie
    hvol = system%hvol

    allocate(T_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(H_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))

    do ispin = 1, system%nspin
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)

      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        n_basis_local = dg_frag%n_basis(ifrag, ispin)
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)

        do jo = 1, n_basis_local
          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

          do io = 1, n_basis_local
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, integral)
            mat_H_local(io, jo, ispin, i_local) = real(integral, kind=8)
          end do

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            l = dg_frag%halo(i_halo)%length
            do io = 1, n_basis_halo
              integral = (0.0d0, 0.0d0)
              do iz = 1, l(3)
                do iy = 1, l(2)
                  do ix = 1, l(1)
                    call get_halo_block_point_indices(dg_frag%halo(i_halo), ix, iy, iz, halo_send_idx, halo_recv_idx)
                    if (allocated(dg_frag%halo(i_halo)%buf_recv_c)) then
                      integral = integral + conjg(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, io, 1)) * &
                                 H_phi(halo_recv_idx(1), halo_recv_idx(2), halo_recv_idx(3)) * hvol
                    else
                      integral = integral + cmplx(dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1), 0.0d0, kind=8) * &
                                 H_phi(halo_recv_idx(1), halo_recv_idx(2), halo_recv_idx(3)) * hvol
                    end if
                  end do
                end do
              end do
              halo_H_local(io, jo, ispin, i_halo, i_local) = real(integral, kind=8)
            end do
          end do

        end do
      end do
    end do

    deallocate(T_phi, H_phi, V_total)

    allocate(coef_local(dg_frag%nstate_frag, dg_frag%nstate_tot, dg_frag%nspin, ifrag_count))
    coef_local = (0.0d0, 0.0d0)

    if (yn_eigenexa == 'y') then
#ifdef USE_EIGENEXA
      call diag_full_eigenexa
#else
      stop "EigenExa does not enabled, please check your build configuration."
#endif
    else
      call diag_full_lapack
    end if

    dg_frag%coef = (0.0d0, 0.0d0)
    do ispin = 1, dg_frag%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        n_basis_local = dg_frag%n_basis(ifrag, ispin)
        do jo = 1, n_basis_local
          j = dg_frag%index_basis(jo, ifrag, ispin)
          if (j < 1 .or. j > dg_frag%n_mat_max) cycle
          if (.not. allocated(dg_frag%coef_global_to_local)) cycle
          jloc = dg_frag%coef_global_to_local(j, ispin)
          if (jloc < 1 .or. jloc > size(dg_frag%coef, 1)) cycle
          dg_frag%coef(jloc, :, ispin) = dg_frag%coef(jloc, :, ispin) + coef_local(jo, :, ispin, i_local)
        end do
      end do
    end do
    dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)

    deallocate(mat_H_local, halo_H_local, coef_local)

  contains

    subroutine diag_full_lapack
      use rt_dg_fragment_ops, only: copy_matrix_blocks_metric_to_real_dense
      implicit none
      integer :: i, j, idx_io, idx_jo, valid_local_count, valid_halo_count
      integer :: local_gid(dg_frag%nstate_frag), halo_gid(dg_frag%nstate_frag)
      integer :: valid_local_ids(dg_frag%nstate_frag), valid_halo_ids(dg_frag%nstate_frag)
      real(8), allocatable :: mat_H(:,:), mat_V(:,:)

      call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                              diagonal_only=.true.)

      do ispin = 1, system%nspin
        n = dg_frag%n_mat(ispin)
        if (n <= 0) cycle
        allocate(mat_V(n, n))
        mat_V = 0.0d0

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          n_basis_local = dg_frag%n_basis(ifrag, ispin)
          valid_local_count = 0
          do io = 1, n_basis_local
            local_gid(io) = dg_frag%index_basis(io, ifrag, ispin)
            if (local_gid(io) < 1 .or. local_gid(io) > n) cycle
            valid_local_count = valid_local_count + 1
            valid_local_ids(valid_local_count) = io
          end do
          do idx_jo = 1, valid_local_count
            jo = valid_local_ids(idx_jo)
            j = local_gid(jo)
            do idx_io = 1, valid_local_count
              io = valid_local_ids(idx_io)
              i = local_gid(io)
              mat_V(i, j) = mat_V(i, j) + mat_H_local(io, jo, ispin, i_local)
            end do
          end do

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            valid_halo_count = 0
            do io = 1, n_basis_halo
              halo_gid(io) = dg_frag%index_basis(io, jfrag, ispin)
              if (halo_gid(io) < 1 .or. halo_gid(io) > n) cycle
              valid_halo_count = valid_halo_count + 1
              valid_halo_ids(valid_halo_count) = io
            end do
            do idx_jo = 1, valid_local_count
              jo = valid_local_ids(idx_jo)
              j = local_gid(jo)
              do idx_io = 1, valid_halo_count
                io = valid_halo_ids(idx_io)
                i = halo_gid(io)
                mat_V(i, j) = mat_V(i, j) + 0.5d0 * halo_H_local(io, jo, ispin, i_halo, i_local)
                mat_V(j, i) = mat_V(j, i) + 0.5d0 * halo_H_local(io, jo, ispin, i_halo, i_local)
              end do
            end do
          end do
        end do
        do idx_jo = 1, dg_frag%n_H_blocks
          ifrag = dg_frag%H_mat_blocks(idx_jo)%ifrag_row
          jfrag = dg_frag%H_mat_blocks(idx_jo)%ifrag_col
          n_basis_local = dg_frag%n_basis(ifrag, ispin)
          n_basis_halo = dg_frag%n_basis(jfrag, ispin)
          do jo = 1, n_basis_halo
            j = dg_frag%index_basis(jo, jfrag, ispin)
            if (j < 1 .or. j > n) cycle
            do io = 1, n_basis_local
              i = dg_frag%index_basis(io, ifrag, ispin)
              if (i < 1 .or. i > n) cycle
              dg_frag%H_mat_blocks(idx_jo)%val(io, jo, ispin) = mat_V(i, j)
            end do
          end do
        end do

        deallocate(mat_V)
      end do

      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-basis-diag-soi", dg_frag%icomm)

      do ispin = 1, system%nspin
        n = dg_frag%n_mat(ispin)
        if (n <= 0) cycle
        allocate(mat_H(n, n), mat_V(n, n))
        mat_H = 0.0d0
        call copy_matrix_blocks_metric_to_real_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n, mat_H)
        call eigen_dsyev(mat_H, dg_frag%esp(1:n, ispin), mat_V)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          n_basis_local = dg_frag%n_basis(ifrag, ispin)
          do i = 1, min(n, dg_frag%nstate_tot)
            do jo = 1, n_basis_local
              j = dg_frag%index_basis(jo, ifrag, ispin)
              if (j < 1 .or. j > n) cycle
              coef_local(jo, i, ispin, i_local) = mat_V(j, i)
            end do
          end do
        end do

        deallocate(mat_H, mat_V)
      end do
    end subroutine diag_full_lapack

#ifdef USE_EIGENEXA
    subroutine diag_full_eigenexa
      implicit none
      integer :: nx, ny, ix_s, ix_e, iy_s, iy_e, ix_loc, iy_loc
      integer :: nnod, x_nnod, y_nnod, inod, x_inod, y_inod
      integer :: ifrag_x, ifrag_y, io_x, io_y
      integer :: jfrag_halo(dg_frag%n_halo)
      integer, allocatable :: io_array(:), ifrag_array(:)
      real(8), allocatable :: h_div(:,:), v_div(:,:), h(:,:,:)
      real(8), allocatable :: v_tmp1(:,:), v_tmp2(:,:)

      do ispin = 1, system%nspin
        n = dg_frag%n_mat(ispin)
        if (n <= 0) cycle

        allocate(io_array(n), ifrag_array(n))
        io_array = 0
        ifrag_array = 0
        do ifrag = 1, dg_frag%n_frag
          do io = 1, dg_frag%n_basis(ifrag, ispin)
            i = dg_frag%index_basis(io, ifrag, ispin)
            if (i < 1 .or. i > n) cycle
            io_array(i) = io
            ifrag_array(i) = ifrag
          end do
        end do

        call eigen_init(dg_frag%icomm)
        call eigen_get_matdims(n, nx, ny)
        call eigen_get_procs(nnod, x_nnod, y_nnod)
        call eigen_get_id(inod, x_inod, y_inod)

        allocate(h_div(nx, ny), v_div(nx, ny))
        allocate(h(dg_frag%nstate_frag, dg_frag%nstate_frag, 0:dg_frag%n_halo))
        allocate(v_tmp1(dg_frag%nstate_frag, dg_frag%nstate_tot))
        allocate(v_tmp2(dg_frag%nstate_frag, dg_frag%nstate_tot))

        ix_s = eigen_loop_start(1, x_nnod, x_inod)
        ix_e = eigen_loop_end(n, x_nnod, x_inod)
        iy_s = eigen_loop_start(1, y_nnod, y_inod)
        iy_e = eigen_loop_end(n, y_nnod, y_inod)

        h_div = 0.0d0
        do ifrag = 1, dg_frag%n_frag
          if (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end) then
            i_local = ifrag - dg_frag%ifrag_start + 1
            h(:, :, 0) = mat_H_local(:, :, ispin, i_local)
            jfrag_halo(:) = -1
            h(:, :, 1:dg_frag%n_halo) = 0.0d0
            do i_halo = 1, dg_frag%n_halo
              if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
              jfrag_halo(i_halo) = dg_frag%halo(i_halo)%ifrag_src
              h(:, :, i_halo) = halo_H_local(:, :, ispin, i_halo, i_local)
            end do
          end if
          call comm_bcast(h, dg_frag%icomm, dg_frag%id_array(ifrag))
          call comm_bcast(jfrag_halo, dg_frag%icomm, dg_frag%id_array(ifrag))

          do iy_loc = iy_s, iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            ifrag_y = ifrag_array(iy)
            io_y = io_array(iy)
            if (ifrag_y <= 0 .or. io_y <= 0) cycle
            do ix_loc = ix_s, ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if (ifrag_x <= 0 .or. io_x <= 0) cycle
              if (ifrag_x == ifrag .and. ifrag_y == ifrag) then
                h_div(ix_loc, iy_loc) = h(io_x, io_y, 0)
              end if
              do i_halo = 1, dg_frag%n_halo
                if (ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag) then
                  h_div(ix_loc, iy_loc) = h_div(ix_loc, iy_loc) + 0.5d0 * h(io_x, io_y, i_halo)
                else if (ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo)) then
                  h_div(ix_loc, iy_loc) = h_div(ix_loc, iy_loc) + 0.5d0 * h(io_y, io_x, i_halo)
                end if
              end do
            end do
          end do
        end do

        call eigen_sx(n, n, h_div, nx, dg_frag%esp(1:n, ispin), v_div, nx)

        do ifrag = 1, dg_frag%n_frag
          v_tmp1 = 0.0d0
          do iy_loc = iy_s, iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            do ix_loc = ix_s, ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if (ifrag_x <= 0 .or. io_x <= 0) cycle
              if (iy <= dg_frag%nstate_tot .and. ifrag_x == ifrag) then
                v_tmp1(io_x, iy) = v_div(ix_loc, iy_loc)
              end if
            end do
          end do
          call comm_summation(v_tmp1, v_tmp2, dg_frag%nstate_frag * dg_frag%nstate_tot, dg_frag%icomm)
          if (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end) then
            i_local = ifrag - dg_frag%ifrag_start + 1
            coef_local(:, :, ispin, i_local) = cmplx(v_tmp2, 0.0d0, kind=8)
          end if
        end do

        deallocate(h_div, v_div, h, v_tmp1, v_tmp2, io_array, ifrag_array)
        call eigen_free()
      end do

    end subroutine diag_full_eigenexa
#endif

  end subroutine diagonalize_full_system_dg

  !=======================================================================
  ! Save checkpoint for adaptive basis update
  !=======================================================================
  subroutine save_adaptive_basis_checkpoint(dg_frag, system, rho, Vh, Vxc)
    use structures
    use communication, only: comm_is_root
    use filesystem, only: create_directory, directory_exists
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: rho, Vh
    type(s_scalar), intent(in) :: Vxc(system%nspin)
    
    character(256) :: dirname, filename
    integer :: iunit, ix, iy, iz, ispin
    real(8) :: total_charge
    
    dirname = './data_for_restart/'
    
    if (comm_is_root(dg_frag%id)) then
      if (.not. directory_exists(dirname)) then
        call create_directory(dirname)
      end if
      
      ! Save update metadata
      filename = trim(dirname) // 'adaptive_basis_update.log'
      open(newunit=iunit, file=trim(filename), status='replace')
      write(iunit,'(a)') "# ============================================="
      write(iunit,'(a)') "# Adaptive Basis Update Checkpoint"
      write(iunit,'(a)') "# ============================================="
      write(iunit,'(a)') ""
      write(iunit,'(a,i0)') "# Update count: ", dg_frag%nbasis_update_count
      write(iunit,'(a,es20.12)') "# Hamiltonian change ||ΔH|| (a.u.): ", &
                                 dg_frag%hamiltonian_change_norm
      write(iunit,'(a,es20.12)') "# Threshold (a.u.): ", dg_frag%basis_update_threshold
      write(iunit,'(a)') ""
      
      ! Calculate total charge for verification
      total_charge = sum(rho%f) * system%hvol
      write(iunit,'(a)') "# Current system state:"
      write(iunit,'(a,es20.12)') "# Total electrons (interpolated): ", total_charge
      write(iunit,'(a,i0)') "# Number of fragments: ", dg_frag%n_frag
      write(iunit,'(a,i0)') "# States per fragment: ", dg_frag%nstate_frag
      write(iunit,'(a)') ""
      
      write(iunit,'(a)') "# Adaptive basis update workflow (AUTOMATED):"
      write(iunit,'(a)') "# ============================================="
      write(iunit,'(a)') "# When ||ΔH|| > threshold, RT calculation:"
      write(iunit,'(a)') "#   1. Saves old basis for overlap calculation"
      write(iunit,'(a)') "#   2. Saves current potentials to ./rt_potentials/"
      write(iunit,'(a)') "#      - Vh.bin, Vxc*.bin, Vpsl.bin"
      write(iunit,'(a)') "#   3. Launches DC-LCFO: mpirun salmon < inputfile_dclcfo_update"
      write(iunit,'(a)') "#   4. DC-LCFO reads saved potentials (no SCF needed)"
      write(iunit,'(a)') "#   5. DC-LCFO expands basis → ./data_dcdft/fragments/"
      write(iunit,'(a)') "#   6. RT reloads new basis functions"
      write(iunit,'(a)') "#   7. Calculates overlap matrix: S_ji = ⟨φ'_j|φ_i⟩"
      write(iunit,'(a)') "#   8. Projects wave functions: c'_j = Σ_i S_ji c_i"
      write(iunit,'(a)') "#   9. Continues time evolution with new basis"
      write(iunit,'(a)') ""
      write(iunit,'(a)') "# Requirements:"
      write(iunit,'(a)') "# - Prepare inputfile_dclcfo_update (see template)"
      write(iunit,'(a)') ""
      
      write(iunit,'(a)') "# Implementation status:"
      write(iunit,'(a)') "# ============================================="
      write(iunit,'(a)') "# ✓ MEMORY-OPTIMIZED (NO FILE I/O):"
      write(iunit,'(a)') "#   ✓ In-place Hamiltonian diagonalization (diagonalize_full_system_dg)"
      write(iunit,'(a)') "#   ✓ Zero file I/O - potentials kept in memory"
      write(iunit,'(a)') "#   ✓ Basis overlap calculation (calculate_new_old_basis_overlap)"
      write(iunit,'(a)') "#   ✓ Wave function projection (project_wavefunction_to_new_basis)"
      write(iunit,'(a)') "#   ✓ Momentum matrix recalculation for new basis"
      write(iunit,'(a)') "#   ✓ Full nonlocal PP support (ppg) in memory"
      write(iunit,'(a)') ""
      
      write(iunit,'(a)') "# Benefits of basis recalculation (if using DC-LCFO):"
      write(iunit,'(a)') "# ============================================="
      write(iunit,'(a)') "# - DC-LCFO starts from near-converged potentials"
      write(iunit,'(a)') "# - Faster convergence in GS calculation (fewer iterations)"
      write(iunit,'(a)') "# - Smooth basis evolution (small overlap changes)"
      write(iunit,'(a)') "# - Real basis from DC-LCFO (no gauge rotation needed)"
      
      close(iunit)
      
      write(*,*) "  Checkpoint saved:", trim(filename)
    end if
    
  end subroutine save_adaptive_basis_checkpoint

  !=======================================================================
  ! Helper routines for DC-LCFO CG integration
  !=======================================================================

  !=======================================================================
  ! Update fragment basis functions with CG minimization using old basis as initial guess
  !=======================================================================
  subroutine update_fragment_basis_via_cg(dg_frag, system, info, mg, stencil, srg, ppg, Vh, Vxc, Vpsl, basis_functions_changed)
    use structures
    use parallelization, only: nproc_id_global
    use Conjugate_Gradient, only: gscg_rwf, gscg_zwf
    use salmon_global, only: ncg, nscf
    use sendrecv_grid, only: init_sendrecv_grid, dealloc_cache
    use communication, only: comm_proc_null
    use rt_dg_fragment_parallel, only: setup_fragment_system, finalize_fragment_parallel
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_pp_grid), intent(in) :: ppg
    type(s_scalar), intent(in) :: Vh, Vxc(system%nspin), Vpsl
    logical, intent(out) :: basis_functions_changed

    type(s_scalar), allocatable :: vlocal_frag(:)
    type(s_dft_system) :: system_frag
    type(s_parallel_info) :: info_frag
    type(s_orbital) :: psi, hpsi, ttpsi
    type(s_cg) :: cg
    integer :: ifrag, i_local, ncg_basis_update
    type(s_sendrecv_grid) :: srg_frag
    type(s_rgrid) :: mg_frag
    integer :: neig(1:2,1:3)

    basis_functions_changed = .false.
    if (nscf > 0) then
      ncg_basis_update = nscf
    else
      ncg_basis_update = max(1, ncg)
    end if

    call setup_fragment_system(dg_frag, system, info, mg, system_frag, info_frag, mg_frag)
    call setup_fragment_potential(dg_frag, mg_frag, Vh, Vxc, Vpsl, vlocal_frag)
    call allocate_fragment_orbitals(dg_frag, system_frag, info_frag, mg_frag, psi, hpsi, ttpsi)
    neig = comm_proc_null
    call init_sendrecv_grid(srg_frag, mg_frag, system_frag%no, info_frag%icomm_r, neig)

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      call initialize_orbitals_from_basis(dg_frag, psi, i_local)

      if (system_frag%if_real_orbital) then
        call gscg_rwf(ncg_basis_update, mg_frag, system_frag, info_frag, stencil, ppg, vlocal_frag, srg_frag, psi, hpsi, ttpsi, cg)
      else
        call gscg_zwf(ncg_basis_update, mg_frag, system_frag, info_frag, stencil, ppg, vlocal_frag, srg_frag, psi, hpsi, ttpsi, cg)
      end if
      ! gscg already performs orthonormalization/subspace handling internally.
      ! Avoid an extra global Gram-Schmidt here to keep basis-update local-safe.

      call extract_basis_from_orbitals(dg_frag, psi, i_local)
    end do

    basis_functions_changed = .true.

    if (allocated(cg%pk%rwf)) deallocate(cg%pk%rwf)
    if (allocated(cg%pk%zwf)) deallocate(cg%pk%zwf)
    if (allocated(cg%pre_gk%rwf)) deallocate(cg%pre_gk%rwf)
    if (allocated(cg%pre_gk%zwf)) deallocate(cg%pre_gk%zwf)
    if (allocated(cg%hpk%rwf)) deallocate(cg%hpk%rwf)
    if (allocated(cg%hpk%zwf)) deallocate(cg%hpk%zwf)

    call deallocate_fragment_orbitals(psi, hpsi, ttpsi)
    call dealloc_cache(srg_frag)
    call finalize_fragment_parallel(info_frag)
    if (allocated(vlocal_frag)) then
      do ifrag = 1, size(vlocal_frag)
        if (allocated(vlocal_frag(ifrag)%f)) deallocate(vlocal_frag(ifrag)%f)
      end do
      deallocate(vlocal_frag)
    end if

    if (nproc_id_global == 0) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "  CG basis-update iterations: ", ncg_basis_update, &
                                      " (nscf=", nscf, ", ncg=", ncg, ")"
      write(*,'(1x,a)') "  Fragment basis updated by CG with old basis as initial guess"
    end if

  end subroutine update_fragment_basis_via_cg

  !=======================================================================
  ! Setup fragment-local potential: V_local = V_ion + Vh + Vxc
  !=======================================================================
  subroutine setup_fragment_potential(dg_frag, mg, Vh, Vxc, Vpsl, vlocal_frag)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid), intent(in) :: mg
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    type(s_scalar), allocatable, intent(inout) :: vlocal_frag(:)

    integer :: ispin

    if (.not. allocated(vlocal_frag)) then
      allocate(vlocal_frag(dg_frag%nspin))
      do ispin = 1, dg_frag%nspin
        allocate(vlocal_frag(ispin)%f(mg%is(1):mg%ie(1), &
                                      mg%is(2):mg%ie(2), &
                                      mg%is(3):mg%ie(3)))
      end do
    end if

    do ispin = 1, dg_frag%nspin
      vlocal_frag(ispin)%f = Vpsl%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)) + &
                             Vh%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)) + &
                             Vxc(ispin)%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    end do

  end subroutine setup_fragment_potential

  !=======================================================================
  ! Setup fragment-local system, communicator topology, and real-space grid
  !=======================================================================
  !=======================================================================
  ! Allocate fragment-local orbitals
  !=======================================================================
  subroutine allocate_fragment_orbitals(dg_frag, system, info, mg, &
                                        psi, hpsi, ttpsi)
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, &
                 allocate_orbital_real, allocate_orbital_complex
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid), intent(in) :: mg
    type(s_orbital), intent(inout) :: psi, hpsi, ttpsi
    
    ! Allocate based on orbital type (real or complex)
    if (system%if_real_orbital) then
      call allocate_orbital_real(system%nspin, mg, info, psi)
      call allocate_orbital_real(system%nspin, mg, info, hpsi)
      call allocate_orbital_real(system%nspin, mg, info, ttpsi)
    else
      call allocate_orbital_complex(system%nspin, mg, info, psi)
      call allocate_orbital_complex(system%nspin, mg, info, hpsi)
      call allocate_orbital_complex(system%nspin, mg, info, ttpsi)
    end if
    
  end subroutine allocate_fragment_orbitals

  !=======================================================================
  ! Initialize orbitals from current basis (use as CG initial guess)
  !=======================================================================
  subroutine initialize_orbitals_from_basis(dg_frag, psi, i_local)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_orbital), intent(inout) :: psi
    integer, intent(in) :: i_local
    
    integer :: istate, jstate, ispin, ix, iy, iz, nstate_use, nb_local, nstate_loc, io_lb, io_idx
    integer :: ifrag, lx, ly, lz, px, py, pz
    integer :: iorg(3), ndom(3)
    integer :: is(3), ie(3)
    
    if (allocated(psi%rwf)) then
      is = [lbound(psi%rwf,1), lbound(psi%rwf,2), lbound(psi%rwf,3)]
      ie = [ubound(psi%rwf,1), ubound(psi%rwf,2), ubound(psi%rwf,3)]
    else if (allocated(psi%zwf)) then
      is = [lbound(psi%zwf,1), lbound(psi%zwf,2), lbound(psi%zwf,3)]
      ie = [ubound(psi%zwf,1), ubound(psi%zwf,2), ubound(psi%zwf,3)]
    else
      return
    end if
    
    ifrag = dg_frag%ifrag_start + i_local - 1
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    nb_local = dg_frag%n_basis(ifrag,1)
    if (allocated(psi%rwf)) then
      io_lb = lbound(psi%rwf,5)
      nstate_use = ubound(psi%rwf,5) - io_lb + 1
    else if (allocated(psi%zwf)) then
      io_lb = lbound(psi%zwf,5)
      nstate_use = ubound(psi%zwf,5) - io_lb + 1
    else
      io_lb = 1
      nstate_use = 0
    end if

    ! Copy basis functions to orbital array
    ! This provides a good initial guess for CG solver
    ! phi_frag is real, so directly copy
    if (allocated(psi%rwf)) then
      ! Real orbitals
      do ispin = 1, dg_frag%nspin
        nstate_loc = min(nb_local, nstate_use)
        do istate = 1, nstate_loc
          io_idx = io_lb + istate - 1
          do iz = is(3), ie(3)
            do iy = is(2), ie(2)
              do ix = is(1), ie(1)
                lx = ix - iorg(1) + 1
                ly = iy - iorg(2) + 1
                lz = iz - iorg(3) + 1
                if (lx >= 1 .and. lx <= ndom(1) .and. &
                    ly >= 1 .and. ly <= ndom(2) .and. &
                    lz >= 1 .and. lz <= ndom(3)) then
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  psi%rwf(ix, iy, iz, ispin, io_idx, 1, 1) = &
                    dg_frag%phi_frag(px, py, pz, istate, i_local)
                else
                  psi%rwf(ix, iy, iz, ispin, io_idx, 1, 1) = 0.0d0
                end if
              end do
            end do
          end do
        end do
      end do
    else if (allocated(psi%zwf)) then
      ! Complex orbitals - convert real basis to complex
      do ispin = 1, dg_frag%nspin
        nstate_loc = min(nb_local, nstate_use)
        do istate = 1, nstate_loc
          io_idx = io_lb + istate - 1
          do iz = is(3), ie(3)
            do iy = is(2), ie(2)
              do ix = is(1), ie(1)
                lx = ix - iorg(1) + 1
                ly = iy - iorg(2) + 1
                lz = iz - iorg(3) + 1
                if (lx >= 1 .and. lx <= ndom(1) .and. &
                    ly >= 1 .and. ly <= ndom(2) .and. &
                    lz >= 1 .and. lz <= ndom(3)) then
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  psi%zwf(ix, iy, iz, ispin, io_idx, 1, 1) = &
                    cmplx(dg_frag%phi_frag(px, py, pz, istate, i_local), 0.0d0, kind=8)
                else
                  psi%zwf(ix, iy, iz, ispin, io_idx, 1, 1) = (0.0d0, 0.0d0)
                end if
              end do
            end do
          end do
        end do
      end do
    end if
    
  end subroutine initialize_orbitals_from_basis

  !=======================================================================
  ! Extract new basis functions from converged CG orbitals
  !=======================================================================
  subroutine extract_basis_from_orbitals(dg_frag, psi, i_local)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_orbital), intent(in) :: psi
    integer, intent(in) :: i_local

    integer :: istate, jstate, ispin, ix, iy, iz, nstate_use, nb_local, nstate_loc, io_lb, io_idx
    integer :: ifrag, lx, ly, lz, px, py, pz
    integer :: iorg(3), ndom(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    real(8) :: normv_local, normv_global, hvol
    complex(8) :: proj_local, proj_global
    integer :: is(3), ie(3)
    real(8), allocatable :: phi_state_sum(:,:,:), phi_state_local(:,:,:)
    complex(8), allocatable :: phi_state_sum_c(:,:,:), phi_state_local_c(:,:,:)

    if (allocated(psi%rwf)) then
      is = [lbound(psi%rwf,1), lbound(psi%rwf,2), lbound(psi%rwf,3)]
      ie = [ubound(psi%rwf,1), ubound(psi%rwf,2), ubound(psi%rwf,3)]
    else if (allocated(psi%zwf)) then
      is = [lbound(psi%zwf,1), lbound(psi%zwf,2), lbound(psi%zwf,3)]
      ie = [ubound(psi%zwf,1), ubound(psi%zwf,2), ubound(psi%zwf,3)]
    else
      return
    end if

    ifrag = dg_frag%ifrag_start + i_local - 1
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    lx_lo = 1
    lx_hi = ndom(1)
    ly_lo = 1
    ly_hi = ndom(2)
    lz_lo = 1
    lz_hi = ndom(3)
    nb_local = dg_frag%n_basis(ifrag,1)
    if (allocated(psi%rwf)) then
      io_lb = lbound(psi%rwf,5)
      nstate_use = ubound(psi%rwf,5) - io_lb + 1
    else if (allocated(psi%zwf)) then
      io_lb = lbound(psi%zwf,5)
      nstate_use = ubound(psi%zwf,5) - io_lb + 1
    else
      io_lb = 1
      nstate_use = 0
    end if

    hvol = dg_frag%hgs(1) * dg_frag%hgs(2) * dg_frag%hgs(3)
    allocate(phi_state_sum(ndom(1), ndom(2), ndom(3)))
    allocate(phi_state_local(ndom(1), ndom(2), ndom(3)))
    if (allocated(dg_frag%phi_frag_c)) then
      allocate(phi_state_sum_c(ndom(1), ndom(2), ndom(3)))
      allocate(phi_state_local_c(ndom(1), ndom(2), ndom(3)))
    end if

    do ispin = 1, dg_frag%nspin
      nstate_loc = min(nb_local, nstate_use)
      do istate = 1, nstate_loc
        io_idx = io_lb + istate - 1
        phi_state_local(:, :, :) = 0.0d0
        if (allocated(dg_frag%phi_frag_c)) phi_state_local_c(:, :, :) = (0.0d0, 0.0d0)
        do lz = lz_lo, lz_hi
          iz = iorg(3) + lz - 1
          do ly = ly_lo, ly_hi
            iy = iorg(2) + ly - 1
            do lx = lx_lo, lx_hi
              ix = iorg(1) + lx - 1
              if (ix < is(1) .or. ix > ie(1) .or. iy < is(2) .or. iy > ie(2) .or. iz < is(3) .or. iz > ie(3)) cycle
              if (allocated(psi%rwf)) then
                if (allocated(dg_frag%phi_frag_c)) then
                  phi_state_local_c(lx, ly, lz) = cmplx(psi%rwf(ix, iy, iz, ispin, io_idx, 1, 1), 0.0d0, kind=8)
                else
                  phi_state_local(lx, ly, lz) = psi%rwf(ix, iy, iz, ispin, io_idx, 1, 1)
                end if
              else if (allocated(psi%zwf)) then
                if (allocated(dg_frag%phi_frag_c)) then
                  phi_state_local_c(lx, ly, lz) = psi%zwf(ix, iy, iz, ispin, io_idx, 1, 1)
                else
                  phi_state_local(lx, ly, lz) = real(psi%zwf(ix, iy, iz, ispin, io_idx, 1, 1), kind=8)
                end if
              end if
            end do
          end do
        end do
        do lz = lz_lo, lz_hi
          do ly = ly_lo, ly_hi
            do lx = lx_lo, lx_hi
              call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
              if (allocated(dg_frag%phi_frag_c)) then
                dg_frag%phi_frag_c(px, py, pz, istate, i_local) = phi_state_local_c(lx, ly, lz)
                dg_frag%phi_frag(px, py, pz, istate, i_local) = real(phi_state_local_c(lx, ly, lz), kind=8)
              else
                dg_frag%phi_frag(px, py, pz, istate, i_local) = phi_state_local(lx, ly, lz)
              end if
            end do
          end do
        end do

        do jstate = 1, istate - 1
          proj_local = (0.0d0, 0.0d0)
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                if (allocated(dg_frag%phi_frag_c)) then
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  proj_local = proj_local + conjg(dg_frag%phi_frag_c(px, py, pz, jstate, i_local)) * &
                                            dg_frag%phi_frag_c(px, py, pz, istate, i_local) * hvol
                else
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  proj_local = proj_local + cmplx(dg_frag%phi_frag(px, py, pz, jstate, i_local), 0.0d0, kind=8) * &
                                            cmplx(dg_frag%phi_frag(px, py, pz, istate, i_local), 0.0d0, kind=8) * hvol
                end if
              end do
            end do
          end do
          call comm_summation(proj_local, proj_global, dg_frag%icomm_frag)
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                if (allocated(dg_frag%phi_frag_c)) then
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  dg_frag%phi_frag_c(px, py, pz, istate, i_local) = &
                    dg_frag%phi_frag_c(px, py, pz, istate, i_local) - &
                    proj_global * dg_frag%phi_frag_c(px, py, pz, jstate, i_local)
                else
                  call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                  dg_frag%phi_frag(px, py, pz, istate, i_local) = &
                    dg_frag%phi_frag(px, py, pz, istate, i_local) - &
                    real(proj_global, kind=8) * dg_frag%phi_frag(px, py, pz, jstate, i_local)
                end if
              end do
            end do
          end do
        end do

        normv_local = 0.0d0
        do lz = lz_lo, lz_hi
          do ly = ly_lo, ly_hi
            do lx = lx_lo, lx_hi
              if (allocated(dg_frag%phi_frag_c)) then
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                normv_local = normv_local + abs(dg_frag%phi_frag_c(px, py, pz, istate, i_local))**2
              else
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                normv_local = normv_local + dg_frag%phi_frag(px, py, pz, istate, i_local)**2
              end if
            end do
          end do
        end do
        call comm_summation(normv_local, normv_global, dg_frag%icomm_frag)
        normv_global = sqrt(max(normv_global * hvol, 1.0d-20))
        if (allocated(dg_frag%phi_frag_c)) then
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                phi_state_local_c(lx, ly, lz) = dg_frag%phi_frag_c(px, py, pz, istate, i_local) / normv_global
              end do
            end do
          end do
          call comm_summation(phi_state_local_c, phi_state_sum_c, ndom(1)*ndom(2)*ndom(3), dg_frag%icomm_frag)
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                dg_frag%phi_frag_c(px, py, pz, istate, i_local) = phi_state_sum_c(lx, ly, lz)
                dg_frag%phi_frag(px, py, pz, istate, i_local) = real(phi_state_sum_c(lx, ly, lz), kind=8)
              end do
            end do
          end do
        else
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                phi_state_local(lx, ly, lz) = dg_frag%phi_frag(px, py, pz, istate, i_local) / normv_global
              end do
            end do
          end do
          call comm_summation(phi_state_local, phi_state_sum, ndom(1)*ndom(2)*ndom(3), dg_frag%icomm_frag)
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              do lx = lx_lo, lx_hi
                call map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
                dg_frag%phi_frag(px, py, pz, istate, i_local) = phi_state_sum(lx, ly, lz)
              end do
            end do
          end do
        end if
      end do
    end do

    if (allocated(phi_state_sum_c)) deallocate(phi_state_sum_c, phi_state_local_c)
    deallocate(phi_state_sum, phi_state_local)

    if (allocated(phi_state_sum_c)) deallocate(phi_state_sum_c)
    deallocate(phi_state_sum)

  end subroutine extract_basis_from_orbitals

  !=======================================================================
  ! Deallocate fragment-local orbitals
  !=======================================================================
  subroutine deallocate_fragment_orbitals(psi, hpsi, ttpsi)
    use structures
    implicit none
    type(s_orbital), intent(inout) :: psi, hpsi, ttpsi
    
    ! Deallocate all orbital arrays
    if (allocated(psi%rwf)) deallocate(psi%rwf)
    if (allocated(psi%zwf)) deallocate(psi%zwf)
    if (allocated(hpsi%rwf)) deallocate(hpsi%rwf)
    if (allocated(hpsi%zwf)) deallocate(hpsi%zwf)
    if (allocated(ttpsi%rwf)) deallocate(ttpsi%rwf)
    if (allocated(ttpsi%zwf)) deallocate(ttpsi%zwf)
    
  end subroutine deallocate_fragment_orbitals
