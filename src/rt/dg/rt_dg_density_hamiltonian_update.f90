  subroutine update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_tot, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use hamiltonian, only: update_vlocal
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use rt_dg_fragment_ops, only: copy_hamiltonian_metric_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt  ! Changed to inout for tpsi0 updates
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: Ac_tot(3)
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_xc_functional),  intent(in)    :: xc_func
    type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_pp_nlcc),        intent(in)    :: ppn
    type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
    type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
    type(s_dft_energy),     intent(inout) :: energy
    logical :: needs_basis_update
    complex(8), allocatable :: H_mat_metric(:,:,:)
    complex(8), allocatable :: H_mat_prev_c(:,:,:)
    real(8), allocatable :: H_mat_prev(:,:,:)
    integer :: n_metric
    logical :: use_hmat_complex
    logical, parameter :: skip_hmat_rebuild = .false.

    ! This implements self-consistent density and Hamiltonian update
    ! Essential for non-perturbative phenomena:
    ! - Photovoltaic effects
    ! - Catalytic reactions under light
    ! - Laser excitation and ionization
    
    ! Check if real-space basis functions are available
    if (.not. dg_frag%has_real_space_basis) then
      if (itt == 1 .and. comm_is_root(nproc_id_global)) then
        write(*,*) "WARNING: Real-space basis functions not available"
        write(*,*) "         Using fixed Hamiltonian approximation"
        write(*,*) "         For self-consistent update, fragment basis functions"
        write(*,*) "         must be read from DC-LCFO data files"
      end if
      return
    end if
    ! Step 1: Calculate electron density from fragment basis coefficients
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
      ! --- NaNチェック: rho ---
      block
        integer :: ix, iy, iz, nan_count
        nan_count = 0
        do ix = mg%is(1), mg%ie(1)
          do iy = mg%is(2), mg%ie(2)
            do iz = mg%is(3), mg%ie(3)
              if (rho%f(ix,iy,iz) /= rho%f(ix,iy,iz)) then
                nan_count = nan_count + 1
              end if
            end do
          end do
        end do
        if (nan_count > 0) then
          write(*,'(1x,a,i0)') '[CHECK] rho NaN count = ', nan_count
        end if
      end block
    if (any(rho%f /= rho%f)) then
      write(*,'(1x,a,i0,a,i0)') "[FATAL] DG-RT invalid density (NaN), rank=", &
                               nproc_id_global, ", itt=", itt
      stop "DG-RT density became NaN"
    end if
    
    ! Step 2: Update Hartree potential from new density
    ! IMPORTANT: Hartree potential is LONG-RANGE (Coulomb interaction)
    !            Must be calculated for the entire system, not per-fragment
    !            Vh(r) = ∫ ρ(r')/|r-r'| dr' includes all fragments
    call hartree_dg_distributed(lg, mg, fg, poisson, dg_frag, rho, Vh)
    if (any(Vh%f /= Vh%f)) then
      write(*,'(1x,a,i0,a,i0)') "[FATAL] DG-RT invalid Hartree (NaN), rank=", nproc_id_global, ", itt=", itt
      stop "DG-RT Hartree became NaN"
    end if
    if (maxval(abs(Vh%f)) > 1.0d120) then
      write(*,'(1x,a,i0,a,i0,a,es12.4)') "[FATAL] DG-RT invalid Hartree, rank=", nproc_id_global, &
        ", itt=", itt, " max=", maxval(abs(Vh%f))
      stop "DG-RT Hartree diverged"
    end if
    
    ! Step 3: Update exchange-correlation potential
    ! IMPORTANT: XC potential is LOCAL/SEMI-LOCAL (LDA/GGA)
    !            Vxc(r) depends only on ρ(r) and ∇ρ(r) at nearby points
    !            Could be calculated per-fragment for efficiency (future optimization)
    !            Current implementation: calculated on full grid for simplicity
    ! Note: For meta-GGA functionals, spsi (wavefunctions) would be needed for τ and j
    !       DG-Fragment RT currently supports LDA/GGA functionals
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    
    ! Step 4: Update DFT+U with explicit on/off branch
    if (dg_frag%use_plusu .and. PLUS_U_ON) then
      if (itt == 1 .and. comm_is_root(nproc_id_global)) then
        write(*,'(1x,a)') "[DG-RT +U] ON branch: updating +U density matrix/potential"
      end if
      call calc_density_matrix_and_energy_plusU(rt%tpsi0, ppg, info, system, energy%E_U)
    else
      if (itt == 1 .and. comm_is_root(nproc_id_global)) then
        write(*,'(1x,a)') "[DG-RT +U] OFF branch: skipping +U update (E_U=0)"
      end if
      energy%E_U = 0.0d0
    end if
    
    ! Step 5: Reconstruct Hamiltonian matrix with updated potentials
    ! H_mat = T + V_psl + V_H + V_xc (potential-dependent terms only)
    ! Note: Ac_tot parameter kept for future use (currently unused)
    if (dg_frag%yn_adaptive_basis) then
      n_metric = min(dg_frag%nstate_frag, dg_frag%n_mat_max)
      allocate(H_mat_prev(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      H_mat_prev(:, :, :) = 0.0d0
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      allocate(H_mat_prev_c(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      H_mat_prev_c(:, :, :) = (0.0d0, 0.0d0)
      if (n_metric > 0) call copy_hamiltonian_metric_to_complex_dense(dg_frag, n_metric, H_mat_prev_c)
      if (n_metric > 0) H_mat_prev(1:n_metric, 1:n_metric, :) = real(H_mat_prev_c(1:n_metric, 1:n_metric, :))
      if (n_metric < dg_frag%nstate_frag .and. itt == 1 .and. comm_is_root(nproc_id_global)) then
        write(*,'(1x,a,i0,a,i0,a)') "[WARN] nstate_frag(", dg_frag%nstate_frag, &
                                     ") > n_mat_max(", dg_frag%n_mat_max, "); truncating adaptive metric block."
      end if
    end if

    if (.not. skip_hmat_rebuild) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    end if
    
    ! Step 6: Check FIELD-FREE Hamiltonian change and trigger adaptive basis update if needed
    if (dg_frag%yn_adaptive_basis) then
      ! Strategy update:
      !   Basis-update metric excludes external-field A contributions.
      !   A-dependent effects are absorbed by fixed PW augmentation in mixed basis.
      allocate(H_mat_metric(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      call copy_hamiltonian_metric_to_complex_dense(dg_frag, n_metric, H_mat_metric)
      
      ! Check field-free Hamiltonian change
      needs_basis_update = check_hamiltonian_change_fragments(dg_frag, H_mat_metric)
      deallocate(H_mat_metric)
      
      if (needs_basis_update) then
        ! User policy: when adaptive basis update is triggered, do not apply
        ! this step's Hamiltonian update to time propagation state.
        if (n_metric > 0) then
          if (allocated(dg_frag%H_mat)) then
            dg_frag%H_mat(1:n_metric, 1:n_metric, :) = H_mat_prev(1:n_metric, 1:n_metric, :)
          end if
          if (use_hmat_complex) then
            dg_frag%H_mat_c(1:n_metric, 1:n_metric, :) = H_mat_prev_c(1:n_metric, 1:n_metric, :)
          end if
        end if

        ! Pass ppg for DC-CG method pseudopotential handling
        if (comm_is_root(nproc_id_global)) then
          write(*,*)
          write(*,'(1x,a,i8)') "!!! Adaptive basis update triggered at step", itt
          write(*,'(1x,a,f10.6,a)') "  Global ||ΔH|| = ", &
                                    dg_frag%hamiltonian_change_norm, " a.u."
          write(*,'(1x,a,f10.6,a)') "  Threshold = ", &
                                    dg_frag%basis_update_threshold, " a.u."
          write(*,*) "  At least one fragment detected significant mean-field change"
        end if
        
        ! Trigger DC-LCFO recalculation and basis rotation
        call trigger_basis_update(dg_frag, system, info, itt, lg, mg, stencil, srg, &
                Vh, Vxc, Vpsl, pp, ppg, rt%tpsi0, Ac_tot)
      else
        ! Log monitoring info every 100 steps
        if (mod(itt, 100) == 0 .and. comm_is_root(nproc_id_global)) then
          write(*,'(1x,a,i8,a,f10.6,a)') "  Step", itt, ": Global ||ΔH|| = ", &
                                         dg_frag%hamiltonian_change_norm, " a.u. (OK)"
        end if
      end if

      if (allocated(H_mat_prev)) deallocate(H_mat_prev)
      if (allocated(H_mat_prev_c)) deallocate(H_mat_prev_c)
    end if
    
  end subroutine update_density_and_hamiltonian
