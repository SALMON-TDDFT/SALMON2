  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use hartree_sub, only: hartree
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
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
    complex(8), allocatable :: S_frag_pw(:,:,:), H_frag_pw(:,:,:)
    integer :: n_frag, n_pw

    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      call refresh_pw_coef_cache(dg_frag)
    end if
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    call hartree(lg, mg, info, system, fg, poisson, srg_scalar, stencil, rho, Vh)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    if (dg_frag%use_plusu .and. PLUS_U_ON) then
      call calc_density_matrix_and_energy_plusU(rt%tpsi0, ppg, info, system, energy%E_U)
    else
      energy%E_U = 0.0d0
    end if
    call reconstruct_hamiltonian_matrix(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      n_frag = dg_frag%n_mat_max
      n_pw = dg_frag%n_plane_waves
      allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
      allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
      call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
      call build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)
      deallocate(S_frag_pw, H_frag_pw)
    end if

  end subroutine update_density_hamiltonian_stage
