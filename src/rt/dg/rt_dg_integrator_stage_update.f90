  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
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
    integer :: n_frag, n_pw
    real(8) :: t_stage0, t_stage1
    real(8) :: dt_density, dt_hartree, dt_xc, dt_reconstruct, dt_mixed, dt_total
    logical, parameter :: enable_stage_update_trace = .false.
    logical, parameter :: enable_stage_update_progress = .false.
    logical, parameter :: enable_stage_hotspot_probe = .false.
    integer, parameter :: hotspot_probe_stride = 10

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        density-hmat stage trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if

    call cpu_time(t_stage0)
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    call cpu_time(t_stage1)
    dt_density = t_stage1 - t_stage0
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-density dt=", dt_density
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-hartree"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    call cpu_time(t_stage1)
    dt_hartree = t_stage1 - t_stage0
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-hartree dt=", dt_hartree
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-xc"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    call cpu_time(t_stage1)
    dt_xc = t_stage1 - t_stage0
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-xc dt=", dt_xc
      flush(6)
    end if

    if (dg_frag%use_plusu .and. PLUS_U_ON) then
      call calc_density_matrix_and_energy_plusU(rt%tpsi0, ppg, info, system, energy%E_U)
    else
      energy%E_U = 0.0d0
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-reconstruct"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    call cpu_time(t_stage1)
    dt_reconstruct = t_stage1 - t_stage0
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-reconstruct dt=", dt_reconstruct
      flush(6)
    end if

    dt_mixed = 0.0d0
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      n_frag = dg_frag%n_mat_max
      n_pw = dg_frag%n_plane_waves

      if (.not. allocated(dg_frag%S_mat_frag_pw)) then
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      else if (size(dg_frag%S_mat_frag_pw,1) /= n_frag .or. size(dg_frag%S_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%S_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%S_mat_frag_pw)
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      end if

      if (.not. allocated(dg_frag%H_mat_frag_pw)) then
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%H_mat_frag_pw)
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      end if

      call cpu_time(t_stage0)
      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, dg_frag%H_mat_frag_pw)
      call build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, dg_frag%S_mat_frag_pw, dg_frag%H_mat_frag_pw)
      call cpu_time(t_stage1)
      dt_mixed = t_stage1 - t_stage0
      if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        density-hmat stage trace: stage=after-mixed-refresh"
        flush(6)
      end if
    end if

    dt_total = dt_density + dt_hartree + dt_xc + dt_reconstruct + dt_mixed
    if (enable_stage_hotspot_probe .and. dg_frag%id == 0 .and. mod(itt, hotspot_probe_stride) == 0) then
      write(*,'(1x,a,i0,6(a,1pe11.4))') "        stage hotspot: itt=", itt, &
        " total=", dt_total, " dens=", dt_density, " hartree=", dt_hartree, &
        " xc=", dt_xc, " recon=", dt_reconstruct, " mixed=", dt_mixed
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=exit"
      flush(6)
    end if

  end subroutine update_density_hamiltonian_stage
