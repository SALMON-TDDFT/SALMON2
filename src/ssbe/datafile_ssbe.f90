module datafile_ssbe
    implicit none
contains


subroutine write_sbe_rt_header(fh)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh, '(a)') "# Real time calculation by SBE"
    write(fh, '(a)') "# Ac_ext: External vector potential field"
    write(fh, '(a)') "# E_ext: External electric field"
    write(fh, '(a)') "# Ac_tot: Total vector potential field"
    write(fh, '(a)') "# E_tot: Total electric field"
    write(fh, '(a)') "# Jm: Matter current density (electrons)"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Ac_ext_x", trim(t_unit_ac%name), &
        & 3, "Ac_ext_y", trim(t_unit_ac%name), &
        & 4, "Ac_ext_z", trim(t_unit_ac%name), &
        & 5, "E_ext_x", trim(t_unit_elec%name), &
        & 6, "E_ext_y", trim(t_unit_elec%name), &
        & 7, "E_ext_z", trim(t_unit_elec%name), &
        & 8, "Ac_tot_x", trim(t_unit_ac%name), &
        & 9, "Ac_tot_y", trim(t_unit_ac%name), &
        & 10, "Ac_tot_z", trim(t_unit_ac%name), &
        & 11, "E_tot_x", trim(t_unit_elec%name), &
        & 12, "E_tot_y", trim(t_unit_elec%name), &
        & 13, "E_tot_z", trim(t_unit_elec%name), &
        & 14, "Jm_x", trim(t_unit_current%name), &
        & 15, "Jm_y", trim(t_unit_current%name), &
        & 16, "Jm_z", trim(t_unit_current%name)
    return
end subroutine write_sbe_rt_header



subroutine write_sbe_rt_line(fh, t, Ac_ext, E_ext, Ac_tot, E_tot, Jm)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    real(8), intent(in) :: t, Ac_ext(3), E_ext(3), Ac_tot(3), E_tot(3), Jm(3)
    write(fh, '(F16.8,99(1X,E23.15E3))') &
        & t * t_unit_time%conv, &
        & Ac_ext(1:3) * t_unit_ac%conv, &
        & E_ext(1:3) * t_unit_elec%conv, &
        & Ac_tot(1:3) * t_unit_ac%conv, &
        & E_tot(1:3) * t_unit_elec%conv, &
        & Jm(1:3) * t_unit_current%conv
    return
end subroutine write_sbe_rt_line


subroutine write_sbe_rt_energy_header(fh)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh,'(a)') "# Real time calculation"
    write(fh,'(a)') "# Eall: Total energy"
    write(fh,'(a)') "# Eall0: Initial energy"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Eall-Eall0", trim(t_unit_energy%name)
    return
end subroutine write_sbe_rt_energy_header



subroutine write_sbe_rt_energy_line(fh, t, E_tot, E_tot_delta)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    real(8), intent(in) :: t, E_tot, E_tot_delta
    write(fh, '(F16.8,99(1X,E23.15E3))') &
        & t * t_unit_time%conv,&
        & E_tot * t_unit_energy%conv, &
        & E_tot_delta * t_unit_energy%conv
    return
end subroutine write_sbe_rt_energy_line



subroutine write_sbe_nex_header(fh)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh,'(a)') "# Excitation"
    write(fh,'(a)') "# nelec: Number of excited electrons"
    write(fh,'(a)') "# nhole: Number of excited holes"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "time", trim(t_unit_time%name), &
        & 2, "nelec", trim(t_unit_length%name)//"^-3", &
        & 3, "nhole", trim(t_unit_length%name)//"^-3"
    return
end subroutine write_sbe_nex_header



subroutine write_sbe_nex_line(fh, t, nelec, nhole)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    real(8), intent(in) :: t, nelec, nhole
    write(fh, '(F16.8,99(1X,E23.15E3))') &
        & t * t_unit_time%conv, &
        & nelec * (t_unit_length%conv ** (-3)), &
        & nhole * (t_unit_length%conv ** (-3))
    return
end subroutine write_sbe_nex_line



subroutine write_sbe_wave_header(fh)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh,'(a)') "# Waveform"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_elec%name), &
        & 2, "E_inc_x", trim(t_unit_elec%name), &
        & 3, "E_inc_y", trim(t_unit_elec%name), &
        & 4, "E_inc_z", trim(t_unit_elec%name), &
        & 5, "E_ref_x", trim(t_unit_elec%name), &
        & 6, "E_ref_y", trim(t_unit_elec%name), &
        & 7, "E_ref_z", trim(t_unit_elec%name), &
        & 8, "E_tra_x", trim(t_unit_elec%name), &
        & 9, "E_tra_y", trim(t_unit_elec%name), &
        & 10, "E_tra_z", trim(t_unit_elec%name)
    return
end subroutine write_sbe_wave_header



subroutine write_sbe_wave_line(fh, t, e_inc, e_ref, e_tra)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    real(8), intent(in) :: t, e_inc(3), e_ref(3), e_tra(3)
    write(fh, '(F16.8,99(1X,E23.15E3))') &
        & t * t_unit_time%conv, &
        & e_inc(1:3) * t_unit_elec%conv, &
        & e_ref(1:3) * t_unit_elec%conv, &
        & e_tra(1:3) * t_unit_elec%conv
    return
end subroutine write_sbe_wave_line



subroutine write_sbe_obs_header(fh)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh,'(a)') "# Waveform"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "E_x", trim(t_unit_elec%name), &
        & 3, "E_y", trim(t_unit_elec%name), &
        & 4, "E_z", trim(t_unit_elec%name), &
        & 5, "H_x", "[none]", &
        & 6, "H_y", "[none]", &
        & 7, "H_z", "[none]"
    return
end subroutine write_sbe_obs_header



subroutine write_sbe_obs_line(fh, t, e)
    use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    real(8), intent(in) :: t, e(3)
    write(fh, '(F16.8,99(1X,E23.15E3))') &
        & t * t_unit_time%conv, &
        & e(1:3) * t_unit_elec%conv
    return
end subroutine write_sbe_obs_line


end module datafile_ssbe