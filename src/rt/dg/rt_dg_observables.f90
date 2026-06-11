  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: calculate_macroscopic_current_dg
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
    real(8) :: current_raw(3)
    real(8) :: nelec_ref, ne_density

    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG-Fragment RT: mixed/PW observable route was removed"
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      stop "DG-Fragment RT: non-orbital observable route was removed"
    end if

    call calculate_macroscopic_current_dg(dg_frag, system, mg, stencil, current_raw)
    if (system%ngrid > 0) then
      dg_frag%current(:) = current_raw(:) / dble(system%ngrid)
    else
      dg_frag%current(:) = 0.0d0
    end if
    dg_frag%current_para(:) = dg_frag%current(:)
    dg_frag%current_nl(:) = 0.0d0
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
    rt%curr(:, itt) = dg_frag%current_total(:)
  end subroutine calculate_observables
