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
  use inputoutput,          only: dt_em,utime_from_au,ulength_from_au,uenergy_from_au,unit_system,iperiodic,&
                                  ae_shape1,ae_shape2,e_impulse,sysname,nt_em,nenergy,de, &
                                  directory,iobs_num_em,iobs_samp_em,obs_plane_em
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use structures,           only: s_fdtd_system
  use salmon_maxwell,       only: ls_fdtd_work
  use math_constants,       only: pi
  implicit none
  type(s_fdtd_system),intent(in)    :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  integer                           :: ii
  character(128)                    :: save_name
  
  !output linear response(matter dipole pm and current jm are outputted: pm = -dip and jm = -curr)
  if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
    if(iperiodic==0) then
      !output time-dependent dipole data
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_p.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# time[a.u.], dipoleMoment(x,y,z)[a.u.]" 
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# time[fs], dipoleMoment(x,y,z)[Ang.]" 
        end select
        do ii=1,nt_em
          write(fw%ifn, '(E13.5)',advance="no")     fw%time_lr(ii)*utime_from_au
          write(fw%ifn, '(3E16.6e3)',advance="yes") -fw%dip_lr(ii,:)*ulength_from_au
        end do
        close(fw%ifn)
      end if
      
      !output lr data
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,1),fw%fr_lr(:,1),fw%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,2),fw%fr_lr(:,2),fw%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%dip_lr(:,3),fw%fr_lr(:,3),fw%fi_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_lr.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# energy[a.u.], Re[alpha](x,y,z)[a.u.], Im[alpha](x,y,z)[a.u.], df/dE(x,y,z)[a.u.]"
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# energy[eV], Re[alpha](x,y,z)[Ang.**3], Im[alpha](x,y,z)[Ang.**3], df/dE(x,y,z)[1/eV]"
        end select
        do ii=0,nenergy
          write(fw%ifn, '(E13.5)',advance="no")     dble(ii)*de*uenergy_from_au
          write(fw%ifn, '(3E16.6e3)',advance="no")  fw%fr_lr(ii,:)/(-e_impulse)*(ulength_from_au**3.0d0)
          write(fw%ifn, '(3E16.6e3)',advance="no")  fw%fi_lr(ii,:)/(-e_impulse)*(ulength_from_au**3.0d0)
          write(fw%ifn, '(3E16.6e3)',advance="yes") 2.0d0*dble(ii)*de/pi*fw%fi_lr(ii,:)/(-e_impulse)/uenergy_from_au
        end do
        close(fw%ifn)
      end if
    elseif(iperiodic==3) then
      !output time-dependent dipole data
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_current.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# time[a.u.],  current(x,y,z)[a.u.]" 
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# time[fs],    current(x,y,z)[A/Ang.^2]" 
        end select
        do ii=1,nt_em
          write(fw%ifn, '(E13.5)',advance="no")     fw%time_lr(ii)*utime_from_au
          write(fw%ifn, '(3E16.6e3)',advance="yes") -fw%curr_lr(ii,:)*fw%uAperm_from_au/ulength_from_au
        end do
        close(fw%ifn)
      end if
      
      !output lr data
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,1),fw%fr_lr(:,1),fw%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,2),fw%fr_lr(:,2),fw%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,dt_em,de,fw%time_lr,fw%curr_lr(:,3),fw%fr_lr(:,3),fw%fi_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_lr.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# energy[a.u.], Re[epsilon](x,y,z), Im[epsilon](x,y,z)"
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# energy[eV], Re[epsilon](x,y,z), Im[epsilon](x,y,z)"
        end select
        do ii=1,nenergy
          write(fw%ifn, '(E13.5)',advance="no")     dble(ii)*de*uenergy_from_au
          write(fw%ifn, '(3E16.6e3)',advance="no")  1.0d0-4.0d0*pi*fw%fi_lr(ii,:)/(-e_impulse)/(dble(ii)*de)
          write(fw%ifn, '(3E16.6e3)',advance="yes") 4.0d0*pi*fw%fr_lr(ii,:)/(-e_impulse)/(dble(ii)*de)
        end do
      end if
    end if
  end if
  
  !observation
  if(iobs_num_em>0) then
    if(comm_is_root(nproc_id_global)) then
      !make information file
      open(fw%ifn,file=trim(directory)//"/obs0_info.data")
      write(fw%ifn,'(A,A14)')                      'unit_system       =',trim(unit_system)
      write(fw%ifn,'(A,I14)')                      'iperiodic         =',iperiodic
      write(fw%ifn,'(A,ES14.5)')                   'dt_em             =',dt_em*utime_from_au
      write(fw%ifn,'(A,I14)')                      'nt_em             =',(fw%iter_end-fw%iter_sta+1)
      write(fw%ifn,'(A,ES14.5,A,ES14.5,A,ES14.5)') 'al_em             =',&
            fs%rlsize(1)*ulength_from_au,', ',&
            fs%rlsize(2)*ulength_from_au,', ',&
            fs%rlsize(3)*ulength_from_au
      write(fw%ifn,'(A,ES14.5,A,ES14.5,A,ES14.5)') 'dl_em             =',&
            fs%hgs(1)*ulength_from_au,', ',&
            fs%hgs(2)*ulength_from_au,', ',&
            fs%hgs(3)*ulength_from_au
      write(fw%ifn,'(A,I14,A,I14,A,I14)')          'lg_sta            =',&
            fs%lg%is(1),', ',fs%lg%is(2),', ',fs%lg%is(3)
      write(fw%ifn,'(A,I14,A,I14,A,I14)')          'lg_end            =',&
            fs%lg%ie(1),', ',fs%lg%ie(2),', ',fs%lg%ie(3)
      write(fw%ifn,'(A,I14)')                      'iobs_num_em       =',iobs_num_em
      write(fw%ifn,'(A,I14)')                      'iobs_samp_em      =',iobs_samp_em
      do ii=1,iobs_num_em
        write(fw%ifn,'(A,I3,A,A)')                 'obs_plane_em(',&
                                                                ii,') =             ',obs_plane_em(ii)
      end do
      write(fw%ifn,'(A,ES14.5)')                   'e_max             =',fw%e_max
      write(fw%ifn,'(A,ES14.5)')                   'h_max             =',fw%h_max
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

!=========================================================================================
!= Fourier transformation in eh ==========================================================
subroutine eh_fourier(nt,ne,dt,de,ti,ft,fr,fi)
  use inputoutput,    only: wf_em
  use math_constants, only: zi
  implicit none
  integer,intent(in)   :: nt,ne
  real(8),intent(in)   :: dt,de
  real(8),intent(in)   :: ti(nt),ft(nt)
  real(8),intent(out)  :: fr(0:ne),fi(0:ne)
  integer              :: ie,it
  real(8)              :: ft_wf(nt)
  real(8)              :: hw
  complex(8)           :: zf
  
  !apply window function
  if(wf_em=='y') then
    do it=1,nt
      ft_wf(it)=ft(it)*( 1.0d0 -3.0d0*(ti(it)/maxval(ti(:)))**2.0d0 +2.0d0*(ti(it)/maxval(ti(:)))**3.0d0 )
    end do
  else
    ft_wf(:)=ft(:)
  end if
  
  !Fourier transformation
  do ie=0,ne
    hw=dble(ie)*de; zf=(0.0d0,0.0d0);
!$omp parallel
!$omp do private(it) reduction( + : zf )
    do it=1,nt
      zf=zf+exp(zi*hw*ti(it))*ft_wf(it)
    end do
!$omp end do
!$omp end parallel
    zf=zf*dt; fr(ie)=real(zf,8); fi(ie)=aimag(zf)
  end do
  
end subroutine eh_fourier
