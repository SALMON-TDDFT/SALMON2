  subroutine update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_tot, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy, &
                                            skip_hamiltonian_reconstruct, skip_orbital_dependent)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use hamiltonian, only: update_vlocal
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
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
    logical, optional,       intent(in)    :: skip_hamiltonian_reconstruct
    logical, optional,       intent(in)    :: skip_orbital_dependent
    real(8), allocatable :: rho_buffer(:,:,:), Vh_buffer(:,:,:)
    real(8), allocatable :: Vpsl_buffer(:,:,:), Vxc_buffer(:,:,:,:)
    logical :: use_rank_buffered_potential
    logical :: do_hamiltonian_reconstruct
    logical :: allow_orbital_dependent
    integer :: ispin

    do_hamiltonian_reconstruct = .true.
    if (present(skip_hamiltonian_reconstruct)) do_hamiltonian_reconstruct = .not. skip_hamiltonian_reconstruct
    allow_orbital_dependent = .true.
    if (present(skip_orbital_dependent)) allow_orbital_dependent = .not. skip_orbital_dependent

    ! This implements self-consistent density and Hamiltonian update
    ! Essential for non-perturbative phenomena:
    ! - Photovoltaic effects
    ! - Catalytic reactions under light
    ! - Laser excitation and ionization
    
    ! Check if real-space basis functions are available
    if (.not. dg_frag%has_real_space_basis) then
      if (.not. do_hamiltonian_reconstruct) then
        if (comm_is_root(nproc_id_global)) then
          write(*,'(1x,a)') "[FATAL] DG initial density reconstruction requires real-space fragment basis data."
          write(*,'(1x,a)') "[FATAL] Check that DC-LCFO data include fragment basis functions,"
          write(*,'(1x,a)') "[FATAL] not only identity/local coefficients."
        end if
        stop "DG-Fragment RT: missing real-space basis for initial density"
      end if
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

    dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (system%nspin > 0) then
      dg_frag%rho_s_frag(:, :, :, 1:system%nspin) = 0.0d0
      do ispin = 1, system%nspin
        dg_frag%rho_s_frag(:, :, :, ispin) = rho_s(ispin)%f(:, :, :)
      end do
    end if

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    if (use_rank_buffered_potential) then
      allocate(rho_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vpsl_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                           dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                           dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vxc_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3), system%nspin))
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, rho, rho_buffer)
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vpsl, Vpsl_buffer)
    end if
    
    ! Step 2: Update Hartree potential from new density
    ! IMPORTANT: Hartree potential is LONG-RANGE (Coulomb interaction)
    !            Must be calculated for the entire system, not per-fragment
    !            Vh(r) = ∫ ρ(r')/|r-r'| dr' includes all fragments
    call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    dg_frag%Vh_frag(:, :, :) = Vh%f(:, :, :)
    if (use_rank_buffered_potential) then
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vh, Vh_buffer)
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
    ! Step 3: Update exchange-correlation potential
    ! IMPORTANT: XC potential is LOCAL/SEMI-LOCAL (LDA/GGA)
    !            Vxc(r) depends only on ρ(r) and ∇ρ(r) at nearby points
    !            Could be calculated per-fragment for efficiency (future optimization)
    !            Current implementation: calculated on full grid for simplicity
    ! Note: For meta-GGA functionals, spsi (wavefunctions) would be needed for τ and j
    !       DG-Fragment RT currently supports LDA/GGA functionals
    if (.not. allow_orbital_dependent) then
      if (xc_func%use_laplacian .or. &
          xc_func%use_kinetic_energy .or. xc_func%use_current) then
        if (comm_is_root(nproc_id_global)) then
          write(*,'(1x,a)') "[FATAL] DG initial density path cannot evaluate orbital-dependent XC without conventional orbitals."
          write(*,'(1x,a)') "[FATAL] Use an LDA/local XC functional or implement DG tau/current reconstruction."
        end if
        stop "DG-Fragment RT: orbital-dependent XC requires DG reconstruction"
      end if
    end if
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                              info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    if (system%nspin > 0) then
      dg_frag%Vxc_frag(:, :, :, 1:system%nspin) = 0.0d0
      do ispin = 1, system%nspin
        dg_frag%Vxc_frag(:, :, :, ispin) = Vxc(ispin)%f(:, :, :)
        if (use_rank_buffered_potential) then
          call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vxc(ispin), Vxc_buffer(:, :, :, ispin))
        end if
      end do
    end if
    
    ! Step 4: Update DFT+U with explicit on/off branch
    if (dg_frag%use_plusu .and. PLUS_U_ON .and. .not. allow_orbital_dependent) then
      if (comm_is_root(nproc_id_global)) then
        write(*,'(1x,a)') "[FATAL] DG initial density path cannot update DFT+U without conventional orbitals."
        write(*,'(1x,a)') "[FATAL] Implement DG +U density matrix reconstruction before enabling +U in this path."
      end if
      stop "DG-Fragment RT: +U requires DG density matrix reconstruction"
    end if
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

    if (.not. do_hamiltonian_reconstruct) then
      if (allocated(rho_buffer)) deallocate(rho_buffer)
      if (allocated(Vh_buffer)) deallocate(Vh_buffer)
      if (allocated(Vpsl_buffer)) deallocate(Vpsl_buffer)
      if (allocated(Vxc_buffer)) deallocate(Vxc_buffer)
      return
    end if
    
    ! Step 5: Reconstruct Hamiltonian matrix with updated potentials.
    ! Runtime adaptive basis updates are intentionally not supported in DG-Fragment RT.
    if (use_rank_buffered_potential) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, &
                                          Vh_buffer, Vxc_buffer, Vpsl_buffer)
    else
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    end if

    if (allocated(rho_buffer)) deallocate(rho_buffer)
    if (allocated(Vh_buffer)) deallocate(Vh_buffer)
    if (allocated(Vpsl_buffer)) deallocate(Vpsl_buffer)
    if (allocated(Vxc_buffer)) deallocate(Vxc_buffer)

  end subroutine update_density_and_hamiltonian
