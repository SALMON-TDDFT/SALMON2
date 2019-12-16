!
!  Copyright 2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!-----------------------------------------------------------------------------------------
subroutine eh_finalize(fs,fw)
  use salmon_global,        only: dt_em,unit_system,yn_periodic,ae_shape1,ae_shape2,e_impulse,sysname, &
                                  nt_em,nenergy,de,base_directory,obs_num_em,obs_samp_em,yn_obs_plane_em
  use inputoutput,          only: utime_from_au,ulength_from_au,uenergy_from_au
  use salmon_parallel,      only: nproc_id_global
  use communication, only: comm_is_root
  use structures,           only: s_fdtd_system
  use salmon_maxwell,       only: ls_fdtd_work
  use math_constants,       only: pi,zi
  implicit none
  type(s_fdtd_system),intent(in)    :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  integer                           :: ii
  complex(8)                        :: z_tmp(3)
  character(128)                    :: save_name
  
  !output linear response(matter dipole pm and current jm are outputted: pm = -dip and jm = -curr)
  if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
    if(yn_periodic=='n') then
      !output time-dependent dipole data
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
        open(fw%ifn,file=save_name)
        write(fw%ifn,'(A)') "# Real time calculation:" 
        write(fw%ifn,'(A)') "# ddm_e: Change of dipole moment (electrons/plus definition)" 
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[a.u.]",                &
                2, "ddm_e_x[a.u.]",             &
                3, "ddm_e_y[a.u.]",             &
                4, "ddm_e_z[a.u.]"
        case('A_eV_fs')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[fs]",                  &
                2, "ddm_e_x[Angstrom]",         &
                3, "ddm_e_y[Angstrom]",         &
                4, "ddm_e_z[Angstrom]"
        end select
        do ii=1,nt_em
          write(fw%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
                fw%time_lr(ii)*utime_from_au,                  &
                fw%dip_lr(ii,1:3)*ulength_from_au
          write(fw%ifn,*)
        end do
        close(fw%ifn)
      end if
      
      !output response data
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,1),fw%fr_lr(:,1),fw%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,2),fw%fr_lr(:,2),fw%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,3),fw%fr_lr(:,3),fw%fi_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
        open(fw%ifn,file=save_name)
        write(fw%ifn,'(A)') "# Fourier-transform spectra:" 
        write(fw%ifn,'(A)') "# alpha: Polarizability" 
        write(fw%ifn,'(A)') "# df/dE: Strength function" 
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Energy[a.u.]",              &
                2, "Re(alpha_x)[a.u.]",         &
                3, "Im(alpha_x)[a.u.]",         &
                4, "Re(alpha_y)[a.u.]",         &
                5, "Im(alpha_y)[a.u.]",         &
                6, "Re(alpha_z)[a.u.]",         &
                7, "Im(alpha_z)[a.u.]",         &
                8, "df_x/dE[none]",             &
                9, "df_y/dE[none]",             &
                10,"df_z/dE[none]"
        case('A_eV_fs')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Energy[eV]",                &
                2, "Re(alpha_x)[Augstrom^2/V]", &
                3, "Im(alpha_x)[Augstrom^2/V]", &
                4, "Re(alpha_y)[Augstrom^2/V]", &
                5, "Im(alpha_y)[Augstrom^2/V]", &
                6, "Re(alpha_z)[Augstrom^2/V]", &
                7, "Im(alpha_z)[Augstrom^2/V]", &
                8, "df_x/dE[none]",             &
                9, "df_y/dE[none]",             &
                10,"df_z/dE[none]"
        end select
        do ii=1,nenergy
          write(fw%ifn,"(F16.8,99(1X,E23.15E3))",advance='no')                      &
                dble(ii)*de*uenergy_from_au,                                        &
                fw%fr_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                fw%fi_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                fw%fr_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                fw%fi_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                fw%fr_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                fw%fi_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fw%uVperm_from_au), &
                2.0d0*dble(ii)*de/pi*fw%fi_lr(ii,1:3)/(-e_impulse)
          write(fw%ifn,*)
        end do
        close(fw%ifn)
      end if
    elseif(yn_periodic=='y') then
      !output time-dependent average matter current density data and average electric field data
      if(comm_is_root(nproc_id_global)) then
        !average matter current density data
        save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
        open(fw%ifn,file=save_name)
        write(fw%ifn,'(A)') "# Real time calculation:" 
        write(fw%ifn,'(A)') "# Jm: Matter current density (electrons)" 
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[a.u.]",                &
                2, "Jm_x[a.u.]",                &
                3, "Jm_y[a.u.]",                &
                4, "Jm_z[a.u.]"
        case('A_eV_fs')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[fs]",                  &
                2, "Jm_x[1/fs*Angstrom^2]",     &
                3, "Jm_y[1/fs*Angstrom^2]",     &
                4, "Jm_z[1/fs*Angstrom^2]"
        end select
        do ii=1,nt_em
          write(fw%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
                fw%time_lr(ii)*utime_from_au,                  &
               -fw%curr_lr(ii,:)*((ulength_from_au/utime_from_au)/ulength_from_au**3)
          write(fw%ifn,*)
        end do
        close(fw%ifn)
        
        !average electric field data
        save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt_e.data'
        open(fw%ifn,file=save_name)
        write(fw%ifn,'(A)') "# Real time calculation:" 
        write(fw%ifn,'(A)') "# E: Averaged electric field in the unit cell" 
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[a.u.]",                &
                2, "E_x[a.u.]",                 &
                3, "E_y[a.u.]",                 &
                4, "E_z[a.u.]"
        case('A_eV_fs')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Time[fs]",                  &
                2, "E_x[V/Angstrom]",           &
                3, "E_y[V/Angstrom]",           &
                4, "E_z[V/Angstrom]"
        end select
        do ii=1,nt_em
          write(fw%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
               (fw%time_lr(ii)+0.5d0*dt_em)*utime_from_au,     &
                fw%e_lr(ii,:)*fw%uVperm_from_au
          write(fw%ifn,*)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 0.5d0*dt_em is introduced to adjust actual time of electric field !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        close(fw%ifn)
      end if
      
      !output response data
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,1),fw%fr_lr(:,1),fw%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,2),fw%fr_lr(:,2),fw%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,3),fw%fr_lr(:,3),fw%fi_lr(:,3))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr+0.5d0*dt_em,fw%e_lr(:,1),fw%er_lr(:,1),fw%ei_lr(:,1))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr+0.5d0*dt_em,fw%e_lr(:,2),fw%er_lr(:,2),fw%ei_lr(:,2))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr+0.5d0*dt_em,fw%e_lr(:,3),fw%er_lr(:,3),fw%ei_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
        open(fw%ifn,file=save_name)
        write(fw%ifn,'(A)') "# Fourier-transform spectra:" 
        write(fw%ifn,'(A)') "# sigma: Conductivity" 
        write(fw%ifn,'(A)') "# eps: Dielectric constant" 
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'("#",99(1X,I0,":",A))') &
                1, "Energy[a.u.]",              &
                2, "Re(sigma_x)[a.u.]",         &
                3, "Im(sigma_x)[a.u.]",         &
                4, "Re(sigma_y)[a.u.]",         &
                5, "Im(sigma_y)[a.u.]",         &
                6, "Re(sigma_z)[a.u.]",         &
                7, "Im(sigma_z)[a.u.]",         &
                8, "Re(eps_x)[none]",           &
                9, "Im(eps_x)[none]",           &
                10,"Re(eps_y)[none]",           &
                11,"Im(eps_y)[none]",           &
                12,"Re(eps_z)[none]",           &
                13,"Im(eps_z)[none]"
        case('A_eV_fs')
          write(fw%ifn,'("#",99(1X,I0,":",A))')    &
                1, "Energy[eV]",                   &
                2, "Re(sigma_x)[1/fs*V*Angstrom]", &
                3, "Im(sigma_x)[1/fs*V*Angstrom]", &
                4, "Re(sigma_y)[1/fs*V*Angstrom]", &
                5, "Im(sigma_y)[1/fs*V*Angstrom]", &
                6, "Re(sigma_z)[1/fs*V*Angstrom]", &
                7, "Im(sigma_z)[1/fs*V*Angstrom]", &
                8, "Re(eps_x)[none]",              &
                9, "Im(eps_x)[none]",              &
                10,"Re(eps_y)[none]",              &
                11,"Im(eps_y)[none]",              &
                12,"Re(eps_z)[none]",              &
                13,"Im(eps_z)[none]"
        end select
        do ii=1,nenergy
          z_tmp(:)=(fw%fr_lr(ii,:)+zi*fw%fi_lr(ii,:)) / (-e_impulse+fw%er_lr(ii,:)+zi*fw%ei_lr(ii,:))
          write(fw%ifn,"(F16.8,99(1X,E23.15E3))",advance='no')                           &
                dble(ii)*de*uenergy_from_au,                                             &
                real( z_tmp(1))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                aimag(z_tmp(1))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                real( z_tmp(2))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                aimag(z_tmp(2))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                real( z_tmp(3))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                aimag(z_tmp(3))*(1.0d0/utime_from_au/fw%uVperm_from_au/ulength_from_au), &
                1.0d0-4.0d0*pi*aimag(z_tmp(1))/(dble(ii)*de),                            &
                4.0d0*pi*real(z_tmp(1))/(dble(ii)*de),                                   &
                1.0d0-4.0d0*pi*aimag(z_tmp(2))/(dble(ii)*de),                            &
                4.0d0*pi*real(z_tmp(2))/(dble(ii)*de),                                   &
                1.0d0-4.0d0*pi*aimag(z_tmp(3))/(dble(ii)*de),                            &
                4.0d0*pi*real(z_tmp(3))/(dble(ii)*de)
          write(fw%ifn,*)
        end do
      end if
    end if
  end if
  
  !observation
  if(obs_num_em>0) then
    if(comm_is_root(nproc_id_global)) then
      !make information file
      open(fw%ifn,file=trim(base_directory)//"/obs0_info.data")
      write(fw%ifn,'(A,A23)')         'unit_system          =',trim(unit_system)
      write(fw%ifn,'(A,A23)')         'yn_periodic          =',yn_periodic
      write(fw%ifn,'(A,E23.15E3)')    'dt_em                =',dt_em*utime_from_au
      write(fw%ifn,'(A,I23)')         'nt_em                =',(fw%iter_end-fw%iter_sta+1)
      write(fw%ifn,'(3(A,E23.15E3))') 'al_em                =',fs%rlsize(1)*ulength_from_au,', ',&
                                                               fs%rlsize(2)*ulength_from_au,', ',&
                                                               fs%rlsize(3)*ulength_from_au
      write(fw%ifn,'(3(A,E23.15E3))') 'dl_em                =',fs%hgs(1)*ulength_from_au,', ',&
                                                               fs%hgs(2)*ulength_from_au,', ',&
                                                               fs%hgs(3)*ulength_from_au
      write(fw%ifn,'(3(A,I23))')      'lg_sta               =',fs%lg%is(1),', ',fs%lg%is(2),', ',fs%lg%is(3)
      write(fw%ifn,'(3(A,I23))')      'lg_end               =',fs%lg%ie(1),', ',fs%lg%ie(2),', ',fs%lg%ie(3)
      write(fw%ifn,'(A,I23)')         'obs_num_em           =',obs_num_em
      write(fw%ifn,'(A,I23)')         'obs_samp_em          =',obs_samp_em
      do ii=1,obs_num_em
        write(fw%ifn,'(A,I3,A,A23)')  'yn_obs_plane_em(',&
                                                      ii,') =',yn_obs_plane_em(ii)
      end do
      write(fw%ifn,'(A,E23.15E3)')    'e_max                =',fw%e_max
      write(fw%ifn,'(A,E23.15E3)')    'h_max                =',fw%h_max
      close(fw%ifn)
    end if
  end if
  
  !deallocate
  deallocate(fw%ex_y,fw%c1_ex_y,fw%c2_ex_y,fw%ex_z,fw%c1_ex_z,fw%c2_ex_z,&
             fw%ey_z,fw%c1_ey_z,fw%c2_ey_z,fw%ey_x,fw%c1_ey_x,fw%c2_ey_x,&
             fw%ez_x,fw%c1_ez_x,fw%c2_ez_x,fw%ez_y,fw%c1_ez_y,fw%c2_ez_y,&
             fw%hx_y,fw%c1_hx_y,fw%c2_hx_y,fw%hx_z,fw%c1_hx_z,fw%c2_hx_z,&
             fw%hy_z,fw%c1_hy_z,fw%c2_hy_z,fw%hy_x,fw%c1_hy_x,fw%c2_hy_x,&
             fw%hz_x,fw%c1_hz_x,fw%c2_hz_x,fw%hz_y,fw%c1_hz_y,fw%c2_hz_y)
  
  !write end
  if(comm_is_root(nproc_id_global)) then
    write(*,*) "-------------------------------------------------------"
    write(*,*) "**************************"
    write(*,*) "FDTD end"
    write(*,*) "**************************"
  end if
  
end subroutine eh_finalize
