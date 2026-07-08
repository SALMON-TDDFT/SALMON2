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



! Band-resolved excited population n_ex(b,t) (yn_sbe_out_occ='y').  One column
! per SBE window band; column b+1 is
!   n_ex(b,t) = sum_k w_k ( Re rho_bb(k,t) - f0_b(k) ) / sum_k w_k,
! f0 = ground-state occupation.  Summing the conduction-band columns reproduces
! the "nelec" column of <sysname>_sbe_nex.data; the valence columns sum to
! -nhole.  Dimensionless (occupation per unit cell), so no unit conversion.
subroutine write_sbe_occ_header(fh, nb)
    use inputoutput, only: t_unit_time
    implicit none
    integer, intent(in) :: fh
    integer, intent(in) :: nb
    write(fh,'(a)') "# SBE band-resolved excited population"
    write(fh,'(a)') "# n_ex(b,t) = sum_k w_k ( Re rho_bb(k,t) - f0_b(k) ) / sum_k w_k"
    write(fh,'(a)') "#   f0 = ground-state occupation; window band index b = 1..nb"
    write(fh,'(a)') "#   sum over conduction bands b reproduces nelec in *_sbe_nex.data"
    write(fh, '("# 1:time[",A,"]  columns 2..",I0,": n_ex_b for b = 1..",I0," [none]")') &
        & trim(t_unit_time%name), nb + 1, nb
    return
end subroutine write_sbe_occ_header



subroutine write_sbe_occ_line(fh, t, nb, nex_b)
    use inputoutput, only: t_unit_time
    implicit none
    integer, intent(in) :: fh
    integer, intent(in) :: nb
    real(8), intent(in) :: t, nex_b(nb)
    ! Unlimited-repeat edit descriptor (F2008): nb (= nstate_sbe window) can
    ! exceed the fixed "99(...)" width used by the other writers.
    write(fh, '(F16.8,*(1X,E23.15E3))') t * t_unit_time%conv, nex_b(1:nb)
    return
end subroutine write_sbe_occ_line



! Energy-resolved nonequilibrium occupation distribution n(eps,t)
! (yn_sbe_out_occ='y').  gnuplot pm3d "block" layout: for each output time one
! block of nbin rows (time, energy, n_occ), blank line between blocks.
!   n(eps,t) = sum_{b,k} w_k Re rho_bb(k,t) g(eps - eps_b(k)) / sum_k w_k,
! g = unit-normalized Gaussian; integral over eps = electrons/cell in the window.
subroutine write_sbe_edist_header(fh)
    use inputoutput, only: t_unit_time, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    write(fh,'(a)') "# SBE energy-resolved nonequilibrium occupation distribution"
    write(fh,'(a)') "# n(eps,t) = sum_{b,k} w_k Re rho_bb(k,t) g(eps-eps_b(k)) / sum_k w_k"
    write(fh,'(a)') "#   g = unit-normalized Gaussian; int n(eps) deps = electrons/cell (window)"
    write(fh,'(a)') "# gnuplot pm3d block format: one time-block of nbin rows, blank line between"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "time",   trim(t_unit_time%name), &
        & 2, "energy", trim(t_unit_energy%name), &
        & 3, "n_occ",  trim(t_unit_energy%name)//"^-1"
    return
end subroutine write_sbe_edist_header



subroutine write_sbe_edist_block(fh, t, nbin, e_lo, de, dist)
    use inputoutput, only: t_unit_time, t_unit_energy
    implicit none
    integer, intent(in) :: fh
    integer, intent(in) :: nbin
    real(8), intent(in) :: t, e_lo, de, dist(nbin)
    integer :: j
    real(8) :: ec
    do j = 1, nbin
        ec = e_lo + dble(j - 1) * de
        write(fh, '(F16.8,2(1X,E23.15E3))') &
            & t * t_unit_time%conv, &
            & ec * t_unit_energy%conv, &
            & dist(j) / t_unit_energy%conv
    end do
    write(fh, '(a)') ""   ! blank line = gnuplot block separator
    return
end subroutine write_sbe_edist_block


end module datafile_ssbe