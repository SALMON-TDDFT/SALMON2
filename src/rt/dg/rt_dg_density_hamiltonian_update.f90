  subroutine update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_tot, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_global, only: yn_dual_rho_vh_only
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use hartree_sub, only: build_hartree_density_from_rho
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
    real(8), allocatable :: rho_buffer(:,:,:), Vh_buffer(:,:,:)
    integer :: n_metric
    logical :: use_hmat_complex
    logical :: use_rank_buffered_potential
    logical :: use_dual_rho_vh
    real(8) :: t_stage0, t_stage1
    real(8) :: time_density, time_hartree, time_xc, time_reconstruct
    logical, parameter :: enable_update_trace = .false.
    type(s_scalar) :: rho_h

    ! This implements self-consistent density and Hamiltonian update
    ! Essential for non-perturbative phenomena:
    ! - Photovoltaic effects
    ! - Catalytic reactions under light
    ! - Laser excitation and ionization
    time_density = 0.0d0
    time_hartree = 0.0d0
    time_xc = 0.0d0
    time_reconstruct = 0.0d0
    
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
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step1-density-begin"
      flush(6)
    end if
    ! Step 1: Calculate electron density from fragment basis coefficients
    call cpu_time(t_stage0)
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    call cpu_time(t_stage1)
    time_density = time_density + (t_stage1 - t_stage0)

    dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (system%nspin > 0) then
      dg_frag%rho_s_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_metric = 1, system%nspin
        dg_frag%rho_s_frag(:, :, :, n_metric) = rho_s(n_metric)%f(:, :, :)
      end do
    end if

    use_dual_rho_vh = (yn_dual_rho_vh_only == 'y')
    if (use_dual_rho_vh) then
      call allocate_scalar(mg, rho_h)
      if (itt == 1 .and. dg_frag%has_seed_rho_h) then
        rho_h%f(:, :, :) = dg_frag%rho_h_frag(:, :, :)
      else
        call build_hartree_density_from_rho(info, rho, rho_h)
        dg_frag%rho_h_frag(:, :, :) = rho_h%f(:, :, :)
      end if
    end if

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    if (use_rank_buffered_potential) then
      allocate(rho_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, rho, rho_buffer)
    end if
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step1-density-end"
      flush(6)
    end if
    
    ! Step 2: Update Hartree potential from new density
    ! IMPORTANT: Hartree potential is LONG-RANGE (Coulomb interaction)
    !            Must be calculated for the entire system, not per-fragment
    !            Vh(r) = ∫ ρ(r')/|r-r'| dr' includes all fragments
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step2-hartree-begin"
      flush(6)
    end if
    call cpu_time(t_stage0)
    if (use_dual_rho_vh) then
      call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho_h, Vh)
    else
      call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    end if
    call cpu_time(t_stage1)
    time_hartree = time_hartree + (t_stage1 - t_stage0)
    dg_frag%Vh_frag(:, :, :) = Vh%f(:, :, :)
    if (use_rank_buffered_potential) then
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vh, Vh_buffer)
    end if
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step2-hartree-end"
      flush(6)
    end if
    if (any(Vh%f /= Vh%f)) then
      write(*,'(1x,a,i0,a,i0)') "[FATAL] DG-RT invalid Hartree (NaN), rank=", nproc_id_global, ", itt=", itt
      stop "DG-RT Hartree became NaN"
    end if
    if (maxval(abs(Vh%f)) > 1.0d120) then
      write(*,'(1x,a,i0,a,i0,a,es12.4)') "[FATAL] DG-RT invalid Hartree, rank=", nproc_id_global, &
        ", itt=", itt, " max=", maxval(abs(Vh%f))
      stop "DG-RT Hartree diverged"
    end if
    if (use_dual_rho_vh) then
      call deallocate_scalar(rho_h)
    end if
    
    ! Step 3: Update exchange-correlation potential
    ! IMPORTANT: XC potential is LOCAL/SEMI-LOCAL (LDA/GGA)
    !            Vxc(r) depends only on ρ(r) and ∇ρ(r) at nearby points
    !            Could be calculated per-fragment for efficiency (future optimization)
    !            Current implementation: calculated on full grid for simplicity
    ! Note: For meta-GGA functionals, spsi (wavefunctions) would be needed for τ and j
    !       DG-Fragment RT currently supports LDA/GGA functionals
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step3-xc-begin"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    call cpu_time(t_stage1)
    time_xc = time_xc + (t_stage1 - t_stage0)
    if (system%nspin > 0) then
      dg_frag%Vxc_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_metric = 1, system%nspin
        dg_frag%Vxc_frag(:, :, :, n_metric) = Vxc(n_metric)%f(:, :, :)
      end do
    end if
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step3-xc-end"
      flush(6)
    end if
    
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

    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step5-reconstruct-begin"
      flush(6)
    end if
    call cpu_time(t_stage0)
    if (use_rank_buffered_potential) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, Vh_buffer)
    else
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    end if
    call cpu_time(t_stage1)
    time_reconstruct = time_reconstruct + (t_stage1 - t_stage0)
    if (enable_update_trace .and. itt == 1) then
      write(*,'(1x,a,i0,a,i0,a)') "        update trace: rank=", nproc_id_global, ", itt=", itt, " stage=step5-reconstruct-end"
      flush(6)
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

    if (allocated(rho_buffer)) deallocate(rho_buffer)
    if (allocated(Vh_buffer)) deallocate(Vh_buffer)

  end subroutine update_density_and_hamiltonian
