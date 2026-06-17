  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: calculate_macroscopic_current_dg, calculate_nonlocal_current_dg, &
                                  calculate_local_wannier_polarization_dg
    use unusedvar_mod, only: salmon_unusedvar
    use salmon_global, only: yn_dg_length_gauge
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    type(s_scalar),         intent(in)    :: Vh, Vpsl
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    type(s_scalar),         intent(in), optional :: rho
    integer :: ispin
    real(8) :: current_raw(3), current_nl_raw(3), polarization_raw(3), polarization_density(3)
    real(8) :: nelec_ref, ne_density
    character(32), save :: pol_trace_env = ''
    logical, save :: pol_trace_initialized = .false.
    logical, save :: enable_pol_trace = .false.
    integer :: env_status, env_len

    call salmon_unusedvar(Vh)
    call salmon_unusedvar(Vxc)
    call salmon_unusedvar(Vpsl)

    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG-Fragment RT: mixed/PW observable route was removed"
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      stop "DG-Fragment RT: non-orbital observable route was removed"
    end if

    call calculate_macroscopic_current_dg(dg_frag, system, mg, stencil, current_raw)
    call calculate_nonlocal_current_dg(dg_frag, system, mg, ppg, rt%Ac_tot(:, itt), current_nl_raw)
    if (system%ngrid > 0) then
      dg_frag%current(:) = current_raw(:) / dble(system%ngrid)
      dg_frag%current_nl(:) = current_nl_raw(:) / dble(system%ngrid)
    else
      dg_frag%current(:) = 0.0d0
      dg_frag%current_nl(:) = 0.0d0
    end if
    dg_frag%current_para(:) = dg_frag%current(:)
    nelec_ref = 0.0d0
    if (allocated(system%rocc)) then
      do ispin = 1, min(dg_frag%nspin, size(system%rocc, 3))
        nelec_ref = nelec_ref + sum(max(0.0d0, system%rocc(1:min(dg_frag%nstate_tot, size(system%rocc, 1)), 1, ispin)))
      end do
    else if (dg_frag%nspin == 1) then
      nelec_ref = 2.0d0 * dble(max(0, dg_frag%nocc_spin(1)))
    else
      nelec_ref = dble(sum(max(0, dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    if (system%ngrid > 0) then
      ne_density = nelec_ref / dble(system%ngrid)
    else
      ne_density = 0.0d0
    end if
    dg_frag%current_dia(:) = ne_density * rt%Ac_tot(:, itt)
    dg_frag%current_total(:) = dg_frag%current_para(:) + dg_frag%current_nl(:) + dg_frag%current_dia(:)
    dg_frag%dipole_lg_raw(:) = 0.0d0
    dg_frag%polarization_lg(:) = 0.0d0
    if (yn_dg_length_gauge == 'y') then
      if (.not. pol_trace_initialized) then
        call get_environment_variable('SALMON_DG_POL_TRACE', pol_trace_env, length=env_len, status=env_status)
        if (env_status == 0) then
          select case(trim(adjustl(pol_trace_env(1:env_len))))
          case('1','y','Y','yes','YES','true','TRUE','on','ON')
            enable_pol_trace = .true.
          end select
        end if
        pol_trace_initialized = .true.
      end if
      call calculate_local_wannier_polarization_dg(dg_frag, system, polarization_raw)
      dg_frag%dipole_lg_raw(:) = polarization_raw(:)
      polarization_density(:) = 0.0d0
      if (system%ngrid > 0) polarization_density(:) = polarization_raw(:) / dble(system%ngrid)
      if (.not. dg_frag%polarization_lg_ref_ready) then
        dg_frag%polarization_lg_ref(:) = polarization_density(:)
        dg_frag%polarization_lg_ref_ready = .true.
      end if
      dg_frag%polarization_lg(:) = polarization_density(:) - dg_frag%polarization_lg_ref(:)
      if (enable_pol_trace .and. itt <= 5 .and. dg_frag%id == 0) then
        write(*,'(1x,a,i0,9(a,1pe13.5))') '[DG-POL-TRACE] itt=', itt, &
          ' raw_x=', polarization_raw(1), ' raw_y=', polarization_raw(2), ' raw_z=', polarization_raw(3), &
          ' den_x=', polarization_density(1), ' den_y=', polarization_density(2), ' den_z=', polarization_density(3), &
          ' dP_x=', dg_frag%polarization_lg(1), ' dP_y=', dg_frag%polarization_lg(2), &
          ' dP_z=', dg_frag%polarization_lg(3)
        flush(6)
      end if
    end if
    rt%curr(:, itt) = dg_frag%current_total(:)
  end subroutine calculate_observables
