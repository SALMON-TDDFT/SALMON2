  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_global, only: yn_spinorbit
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use noncollinear_module, only: calc_rho_ud_noncollinear, copy_vxc_mat_noncollinear
    use poisson_dg_distributed, only: hartree_dg_distributed
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use rt_dg_plane_wave, only: compute_fragment_pw_hamiltonian, build_mixed_hamiltonian
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
    complex(8), allocatable :: H_frag_pw(:,:,:)
    real(8) :: t_stage0, t_stage1
    logical :: has_vxc_mat
    integer :: env_len, env_status
    character(len=64) :: env_val
    logical :: enable_density_call_probe
    logical, parameter :: enable_stage_update_trace = .false.
    logical, parameter :: enable_stage_update_progress = .false.

    enable_density_call_probe = .false.
    call get_environment_variable("SALMON_DG_ELECTRON_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_density_call_probe = .true.
      end if
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        density-hmat stage trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if

    call cpu_time(t_stage0)
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    if (allocated(dg_frag%rho_frag)) dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (allocated(dg_frag%rho_s_frag)) then
      dg_frag%rho_s_frag(:, :, :, :) = 0.0d0
      dg_frag%rho_s_frag(:, :, :, 1:system%nspin) = rho_s(1:system%nspin)%f
    end if
    if (allocated(dg_frag%rho_ud_frag)) dg_frag%rho_ud_frag(:, :, :) = (0.0d0, 0.0d0)
    if (yn_spinorbit == 'y' .and. allocated(dg_frag%rho_ud_frag)) then
      call calc_rho_ud_noncollinear(rt%tpsi0, system, info, mg, dg_frag%rho_ud_frag)
      if (abs(dg_frag%rho_scale_factor - 1.0d0) > 1.0d-14) then
        dg_frag%rho_ud_frag(:, :, :) = dg_frag%rho_ud_frag(:, :, :) * dg_frag%rho_scale_factor
      end if
    end if
    call cpu_time(t_stage1)
    if (enable_density_call_probe .and. dg_frag%id == 0 .and. (itt == 1 .or. mod(itt, 10) == 0)) then
      write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6)') "        density-call-probe(stage): itt=", itt, &
        " Ne_raw=", dg_frag%elec_num_raw, " rho_scale=", dg_frag%rho_scale_factor
      flush(6)
    end if
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-density dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-hartree"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    call cpu_time(t_stage1)
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-hartree dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-xc"
      flush(6)
    end if
    call cpu_time(t_stage0)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    if (allocated(dg_frag%Vxc_frag)) then
      dg_frag%Vxc_frag(:, :, :, :) = 0.0d0
      dg_frag%Vxc_frag(:, :, :, 1:system%nspin) = Vxc(1:system%nspin)%f
    end if
    if (allocated(dg_frag%Vxc_mat_frag)) then
      dg_frag%Vxc_mat_frag(:, :, :, :, :) = (0.0d0, 0.0d0)
      if (system%nspin >= 1) dg_frag%Vxc_mat_frag(:, :, :, 1, 1) = cmplx(Vxc(1)%f(:, :, :), 0.0d0, kind=8)
      if (system%nspin >= 2) dg_frag%Vxc_mat_frag(:, :, :, 2, 2) = cmplx(Vxc(2)%f(:, :, :), 0.0d0, kind=8)
      if (yn_spinorbit == 'y' .and. system%nspin == 1) then
        dg_frag%Vxc_mat_frag(:, :, :, 2, 2) = dg_frag%Vxc_mat_frag(:, :, :, 1, 1)
      end if
      if (yn_spinorbit == 'y') then
        call copy_vxc_mat_noncollinear(mg, dg_frag%Vxc_mat_frag, has_vxc_mat)
      end if
    end if
    call cpu_time(t_stage1)
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-xc dt=", t_stage1 - t_stage0
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
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-reconstruct dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      if (.not. dg_frag%mixed_basis_ready) then
        stop "DG+PW runtime requires mixed_basis_ready before stage Hamiltonian update"
      end if
      if (.not. allocated(dg_frag%S_mat_frag_pw)) then
        stop "DG+PW stage update requires precomputed S_mat_frag_pw"
      end if
      allocate(H_frag_pw(dg_frag%n_mat_max, dg_frag%n_plane_waves, dg_frag%nspin))
      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
      if (.not. allocated(dg_frag%H_mat_frag_pw)) then
        allocate(dg_frag%H_mat_frag_pw(dg_frag%n_mat_max, dg_frag%n_plane_waves, dg_frag%nspin))
      else if (size(dg_frag%H_mat_frag_pw, 1) /= dg_frag%n_mat_max .or. &
               size(dg_frag%H_mat_frag_pw, 2) /= dg_frag%n_plane_waves .or. &
               size(dg_frag%H_mat_frag_pw, 3) /= dg_frag%nspin) then
        deallocate(dg_frag%H_mat_frag_pw)
        allocate(dg_frag%H_mat_frag_pw(dg_frag%n_mat_max, dg_frag%n_plane_waves, dg_frag%nspin))
      end if
      dg_frag%H_mat_frag_pw(:, :, :) = H_frag_pw(:, :, :)
      call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, dg_frag%S_mat_frag_pw, H_frag_pw)
      deallocate(H_frag_pw)
      if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        density-hmat stage trace: stage=after-mixed-hmat-refresh"
        flush(6)
      end if
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=exit"
      flush(6)
    end if

  end subroutine update_density_hamiltonian_stage
