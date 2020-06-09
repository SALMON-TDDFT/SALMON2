!
!  Copyright 2019-2020 SALMON developers
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
!=======================================================================
!=======================================================================

SUBROUTINE time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func,srg,srg_scalar, &
&   pp,ppg,ppn,spsi_in,spsi_out,tpsi,rho,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc,Vpsl,dmat,fg,energy, &
&   ewald,md,ofl,poisson,singlescale)
  use structures
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use density_matrix, only: calc_density, calc_density_matrix, calc_current, calc_current_use_dmat, calc_microscopic_current
  use writefield
  use timer
  use salmon_global
  use taylor_sub
  use const, only: umass
  use sendrecv_grid, only: s_sendrecv_grid
  use hartree_sub, only: hartree
  use Total_Energy, only: calc_Total_Energy_isolated, calc_Total_Energy_periodic, calc_eigen_energy, check_rion_update
  use force_sub, only: calc_force
  use md_sub, only: time_evolution_step_md_part1,time_evolution_step_md_part2, &
                    update_pseudo_rt
  use write_sub
  use hamiltonian, only: update_kvector_nonlocalpt, update_kvector_nonlocalpt_microAc, update_vlocal
  use fdtd_coulomb_gauge, only: fdtd_singlescale,fourier_singlescale
  use sendrecv_grid, only: update_overlap_complex8
  use salmon_xc
  use em_field, only: calcVbox, calc_emfields
  use dip, only: subdip
  use gram_schmidt_orth, only: gram_schmidt
  implicit none
  integer,intent(in)       :: itt
  integer,intent(in)       :: itotNtime
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(inout) :: system
  type(s_rt),intent(inout) :: rt
  type(s_parallel_info),intent(in) :: info
  type(s_stencil),intent(inout) :: stencil
  type(s_xc_functional),intent(in) :: xc_func
  type(s_sendrecv_grid),intent(inout) :: srg,srg_scalar
  type(s_pp_info),intent(inout) :: pp
  type(s_pp_grid) :: ppg
!  type(s_pp_nlcc),intent(in)    :: ppn
  type(s_pp_nlcc),intent(inout)    :: ppn
  type(s_orbital),intent(inout) :: spsi_in,spsi_out
  type(s_orbital),intent(inout) :: tpsi ! temporary wavefunctions
  type(s_scalar), intent(inout) :: rho,rho_s(system%nspin),V_local(system%nspin),Vh,Vxc(system%nspin),Vpsl
  type(s_scalar), intent(inout) :: Vh_stock1,Vh_stock2,Vbox
  type(s_dmatrix),intent(inout) :: dmat
  type(s_poisson),intent(inout) :: poisson
  type(s_singlescale) :: singlescale
  type(s_reciprocal_grid) :: fg
  type(s_dft_energy) :: energy
  type(s_ewald_ion_ion) :: ewald
  type(s_md) :: md
  type(s_ofile) :: ofl

  integer :: ix,iy,iz,iatom,is,nspin,Mit
  integer :: idensity, idiffDensity, ielf
  real(8) :: rNe  !, FionE(3,system%nion)
  real(8) :: curr_e_tmp(3,2), curr_i_tmp(3)  !??curr_e_tmp(3,nspin) ?
  character(100) :: comment_line
  logical :: rion_update,if_use_dmat
  integer :: ihpsieff

  spsi_out%update_zwf_overlap = .false. 
  nspin = system%nspin

  call timer_begin(LOG_CALC_VBOX)
  
  idensity=0
  idiffDensity=1
  ielf=2
  if_use_dmat = .false. !(singlescale%flag_use==.true.) ! .or. if_metaGGA ! (future work)

  ! for calc_total_energy_periodic
  if(yn_md=='y') then
     rion_update = .true.
  else
     rion_update = check_rion_update()
  endif

  if(ae_shape1 == 'impulse')then
    ihpsieff=0
  else 
    ihpsieff=1
  end if
  
  select case(iperiodic)
  case(0)
    if(ae_shape1 /= 'impulse') call calcVbox(mg,lg,itt,system,rt,Vbox)
    if(ihpsieff==1) then
!$OMP parallel do collapse(3) private(is,iz,iy,ix)
      do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        V_local(is)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz) + Vbox%f(ix,iy,iz)
      end do
      end do
      end do
      end do
    end if
  case(3)
    if(.not.singlescale%flag_use) then
      system%vec_Ac(1:3) = rt%Ac_ext(1:3,itt) + rt%Ac_ind(1:3,itt)
      system%vec_E(1:3) = -((rt%Ac_ext(1:3,itt) + rt%Ac_ind(1:3,itt))-(rt%Ac_ext(1:3,itt-1) + rt%Ac_ind(1:3,itt-1)))/dt
      system%vec_Ac_ext(1:3) = rt%Ac_ext(1:3,itt) 
      system%vec_E_ext(1:3) = -(rt%Ac_ext(1:3,itt) - rt%Ac_ext(1:3,itt-1))/dt
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
    end if
  end select

  call timer_end(LOG_CALC_VBOX)

  !(MD:part1 & update of pseudopotential)
  if(yn_md=='y') then
     call time_evolution_step_md_part1(itt,system,md)
     call update_pseudo_rt(itt,info,system,lg,mg,poisson,fg,pp,ppg,ppn,Vpsl)
  endif

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

  if(propagator == 'aetrs')then
    call time_evolution_half_step_etrs
! predictor-corrector for metaGGA
  else if(xc == 'tbmbj') then
    call predictor_corrector
  end if

  if(info%numo.ge.1)then
    ! spsi_in --> spsi_out (tpsi = working array)
    call taylor(mg,system,info,stencil,srg,spsi_in,spsi_out,tpsi,ppg,V_local,rt)
  end if
    
  call timer_end(LOG_CALC_TIME_PROPAGATION)
  
  ! Gram Schmidt orghonormalization
  if((gram_schmidt_interval >= 1) .and. (mod(itt,gram_schmidt_interval) == 0)) then
    call gram_schmidt(system, mg, info, spsi_out)
  end if

  call timer_begin(LOG_CALC_RHO)
  call calc_density(system,rho_s,spsi_out,info,mg)

  if(nspin==1)then
    !$omp workshare
    rho%f = rho_s(1)%f
    !$omp end workshare
  else if(nspin==2)then
    !$omp workshare
    rho%f = rho_s(1)%f + rho_s(2)%f
    !$omp end workshare
  end if
  call timer_end(LOG_CALC_RHO)
  
  if(singlescale%flag_use) then
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, spsi_out%zwf)
    end if
    spsi_out%update_zwf_overlap = .true.
    call calc_microscopic_current(system,mg,stencil,info,spsi_out,rt%j_e)
  end if

  call timer_begin(LOG_CALC_HARTREE)
  if(iperiodic==0 .and. itt/=1)then
    Vh%f = 2.d0*Vh_stock1%f - Vh_stock2%f
    Vh_stock2%f = Vh_stock1%f
  end if
  if(singlescale%flag_use .and. method_singlescale=='1d_fourier' .and. yn_ffte=='y') then
    call fourier_singlescale(lg,mg,info,fg,rho,rt%j_e,Vh,poisson,singlescale)
  else
    call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  end if
  if(iperiodic==0 .and. itt/=1)then
    Vh_stock1%f = Vh%f
  end if
  call timer_end(LOG_CALC_HARTREE)

  call timer_begin(LOG_CALC_EXC_COR)
  call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,ppn,info,spsi_out,stencil,Vxc,energy%E_xc)
  call timer_end(LOG_CALC_EXC_COR)

  call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)

! result

  call timer_begin(LOG_CALC_PROJECTION)
  if(iwrite_projection==1.and.mod(itt,itwproj)==0)then
    call projection(itt,mg,system,info,spsi_out,tpsi) ! tpsi must be GS orbital (future work)
  end if
  call timer_end(LOG_CALC_PROJECTION)

  if(itt==1.or.itt==itotNtime.or.mod(itt,out_rt_energy_step)==0)then
    call timer_begin(LOG_CALC_EIGEN_ENERGY)
    ! tpsi,spsi_in = working arrays
    call calc_eigen_energy(energy,spsi_out,tpsi,spsi_in,system,info,mg,V_local,stencil,srg,ppg)
    call timer_end(LOG_CALC_EIGEN_ENERGY)
  end if

  select case(iperiodic)
  case(0)

    call calc_Total_Energy_isolated(system,info,mg,pp,rho_s,Vh,Vxc,rion_update,energy)

  case(3)

    call timer_begin(LOG_CALC_DENSITY_MATRIX)
    if(if_use_dmat) then
       call calc_density_matrix(system,info,mg,srg,spsi_out,dmat)
       spsi_out%update_zwf_overlap = .true. 
    endif
    call timer_end(LOG_CALC_DENSITY_MATRIX)

    call timer_begin(LOG_CALC_CURRENT)
    if(if_use_dmat) then
       call calc_current_use_dmat(system,mg,stencil,info,spsi_out,ppg,dmat,curr_e_tmp(1:3,1:nspin))
    else
       call calc_current(system,mg,stencil,info,srg,spsi_out,ppg,curr_e_tmp(1:3,1:nspin))
       spsi_out%update_zwf_overlap = .true. 
    end if
    call calc_emfields(itt,nspin,curr_e_tmp(1:3,1:nspin),rt)
    call timer_end(LOG_CALC_CURRENT)

    if(yn_md=='y') then
      call timer_begin(LOG_CALC_CURRENT_ION)
      call calc_current_ion(lg,system,pp,curr_i_tmp)
      call timer_end(LOG_CALC_CURRENT_ION)
    end if

    call timer_begin(LOG_CALC_TOTAL_ENERGY_PERIODIC)
    call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
    call timer_end(LOG_CALC_TOTAL_ENERGY_PERIODIC)

    if(singlescale%flag_use) then
      call timer_begin(LOG_CALC_SINGLESCALE)
      singlescale%E_electron = energy%E_tot
      call fdtd_singlescale(itt,lg,mg,system,info,rho, &
      & Vh,rt%j_e,srg_scalar,system%Ac_micro,system%div_Ac,singlescale)
      call update_kvector_nonlocalpt_microAc(info%ik_s,info%ik_e,system,ppg)
      call timer_end(LOG_CALC_SINGLESCALE)
    end if

  end select

  call timer_begin(LOG_WRITE_ENERGIES)
  call subdip(info%icomm_r,itt,rt,lg,mg,rho,rNe,poisson,energy%E_tot,system,pp)
  call timer_end(LOG_WRITE_ENERGIES)

  !(force)
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then  ! and or rvf flag in future

     call calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,spsi_out,ewald)

     !force on ion directly from field --- should put in calc_force?
!$omp parallel do private(iatom)
     do iatom=1,system%nion
        system%Force(:,iatom)= system%Force(:,iatom) + pp%Zps(Kion(iatom)) * rt%E_tot(:,itt)
     enddo
!$omp end parallel do

  endif

  !(MD: part2)
  if(yn_md=='y') call time_evolution_step_md_part2(system,md)

  call timer_begin(LOG_WRITE_RT_INFOS)
  ! Output 
  !(Export to SYSname_trj.xyz)
  if( (yn_md=='y'.or.yn_out_rvf_rt=='y') .and. mod(itt,out_rvf_rt_step)==0 )then
     write(comment_line,10) itt, itt*dt
10   format("#rt   step=",i8,"   time=",e16.6)
     if(yn_md=='y') write(comment_line,11) trim(comment_line),md%Temperature
11   format(a,"   T=",f12.4)
     if(yn_md=='y' .and. ensemble=="NVT" .and. thermostat=="nose-hoover") &
          &  write(comment_line,12) trim(comment_line), md%xi_nh
12   format(a,"  xi_nh=",e18.10)
     call write_xyz(comment_line,"add","rvf",system)
  endif


!  if( mod(itt,100) == 0 ) then
  if( mod(itt,1) == 0 ) then  !for debug
     !(Export to SYSname_rt.data)
     select case(iperiodic)
     case(0)
        call write_rt_data_0d(itt,ofl,dt,system,rt)
     case(3)
        call write_rt_data_3d(itt,ofl,dt,system,curr_e_tmp,curr_i_tmp)
     end select

     !(Export to SYSname_rt_energy.data)
     call write_rt_energy_data(itt,ofl,dt,energy,md)

  endif

  if(yn_out_dns_rt=='y')then
    if(mod(itt,out_dns_rt_step)==0)then
      call write_dns(lg,mg,system,rho%f,rho%f,itt)
    end if
  end if
  if(yn_out_dns_ac_je=='y' .and. singlescale%flag_use)then
    if(mod(itt,out_dns_ac_je_step)==0)then
      call write_dns_ac_je(info,mg,system,rho%f,rt%j_e,itt,"bin")
    end if
  end if
  if(yn_out_elf_rt=='y')then
    if(mod(itt,out_elf_rt_step)==0)then
      call write_elf(itt,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi_in)
    end if
  end if
  if(yn_out_estatic_rt=='y')then
    if(mod(itt,out_estatic_rt_step)==0)then
      call write_estatic(lg,mg,system,stencil,info,Vh,srg_scalar,itt)
    end if
  end if
  call timer_end(LOG_WRITE_RT_INFOS)

  return
  
contains

  subroutine predictor_corrector
    implicit none
    real(8) :: V_tmp(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin)
    complex(8) :: psi_tmp(mg%is_array(1):mg%ie_array(1), &
                        & mg%is_array(2):mg%ie_array(2), &
                        & mg%is_array(3):mg%ie_array(3), &
                        & 1:system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e)

    !$omp workshare
    V_tmp(:,:,:,1) = V_local(1)%f
    !$omp end workshare
    if(nspin==2) then
      !$omp workshare
      V_tmp(:,:,:,2) = V_local(2)%f
      !$omp end workshare
    end if
    
    !$omp workshare
    psi_tmp = spsi_in%zwf
    !$omp end workshare
    
!  if(functional == 'VS98' .or. functional == 'TPSS')then
!    tmass_t=tmass
!    tjr_t=tjr
!    tjr2_t=tjr2
!  end if
    
    ! spsi_in --> spsi_out (tpsi = working array)
    call taylor(mg,system,info,stencil,srg,spsi_in,spsi_out,tpsi,ppg,V_local,rt)
    
    call calc_density(system,rho_s,spsi_out,info,mg)
    if(nspin==1)then
      !$omp workshare
      rho%f = rho_s(1)%f
      !$omp end workshare
    else if(nspin==2)then
      !$omp workshare
      rho%f = rho_s(1)%f + rho_s(2)%f
      !$omp end workshare
    end if
    call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
    call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,ppn,info,spsi_out,stencil,Vxc,energy%E_xc)
    call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)
    
    !$omp workshare
    V_local(1)%f = 0.5d0* ( V_tmp(:,:,:,1) + V_local(1)%f )
    !$omp end workshare
    if(nspin==2) then
      !$omp workshare
      V_local(2)%f = 0.5d0* ( V_tmp(:,:,:,2) + V_local(2)%f )
      !$omp end workshare
    end if
    
    !$omp workshare
    spsi_in%zwf = psi_tmp
    !$omp end workshare
    
    if(yn_periodic=='y') then
    
      if(trans_longi=="lo")then
        call calc_current(system,mg,stencil,info,srg,spsi_out,ppg,curr_e_tmp(1:3,1:nspin))
        if(nspin==1) then
          rt%curr(1:3,itt) = curr_e_tmp(1:3,1)
        else if(nspin==2) then
          rt%curr(1:3,itt) = curr_e_tmp(1:3,1) + curr_e_tmp(1:3,2)
        end if
        rt%Ac_ind(:,itt+1) = 2d0*rt%Ac_ind(:,itt) -rt%Ac_ind(:,itt-1) -4d0*Pi*rt%curr(:,itt)*dt**2
      else if(trans_longi=="tr")then
        rt%Ac_ind(:,itt+1) = 0d0
      end if
      rt%Ac_tot(:,itt+1) = rt%Ac_ext(:,itt+1) + rt%Ac_ind(:,itt+1)
      
      system%vec_Ac(1:3) = 0.5d0* ( rt%Ac_tot(:,itt) + rt%Ac_tot(:,itt+1) )
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
      
    end if

!  if(functional == 'VS98' .or. functional == 'TPSS')then
!    tmass=0.5d0*(tmass+tmass_t)
!    tjr=0.5d0*(tjr+tjr_t)
!    tjr2=0.5d0*(tjr2+tjr2_t)
!  end if
    
    return
  end subroutine predictor_corrector

  subroutine time_evolution_half_step_etrs

    if(info%numo.ge.1)then
      ! spsi_in --> spsi_out (tpsi = working array)
      call taylor(mg,system,info,stencil,srg,spsi_in,spsi_out,tpsi,ppg,V_local,rt)
    end if

!$OMP parallel do private(is,iz,iy,ix) collapse(3)
    do is=1,nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rt%vloc_t(is)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz)
      rt%vloc_new(is)%f(ix,iy,iz) = 3d0*V_local(is)%f(ix,iy,iz) - 3d0*rt%vloc_old(is,1)%f(ix,iy,iz) + rt%vloc_old(is,2)%f(ix,iy,iz)
      rt%vloc_old(is,2)%f(ix,iy,iz) = rt%vloc_old(is,1)%f(ix,iy,iz)
      rt%vloc_old(is,1)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz)
      V_local(is)%f(ix,iy,iz) = rt%vloc_new(is)%f(ix,iy,iz)
    end do
    end do
    end do
    end do


  select case(iperiodic)
  case(0)
    if(ae_shape1 /= 'impulse') call calcVbox(mg,lg,itt+1,system,rt,Vbox)
    if(ihpsieff==1) then
!$OMP parallel do collapse(3) private(is,iz,iy,ix)
      do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        V_local(is)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz) + Vbox%f(ix,iy,iz)
      end do
      end do
      end do
      end do
    end if
  case(3)
    if(.not.singlescale%flag_use) then
      system%vec_Ac(1:3) = rt%Ac_ext(1:3,itt+1) + 2d0*rt%Ac_ind(1:3,itt) - rt%Ac_ind(1:3,itt-1)
      system%vec_Ac_ext(1:3) = rt%Ac_ext(1:3,itt+1)
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
    end if
  end select

  spsi_in%zwf = spsi_out%zwf

! self-consistent loop will be implemented for ETRS propagator

  end subroutine time_evolution_half_step_etrs

END SUBROUTINE time_evolution_step

subroutine calc_current_ion(lg,system,pp,curr_i)
  use structures
  use salmon_global, only: natom,Kion
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_dft_system) :: system
  type(s_pp_info) :: pp
  integer :: ia
  real(8) :: curr_i(3)

  !AY memo
  !current of ion: defined by positive charge-->pulse sign
  !This is matter current & = electric current.
  !Be carefull, current by electrons is defined by matter current.
  !Then, total electric current = -curr + curr_ion
  !Be carefull, The definition in ARTED & multiscale is different.
  curr_i(:)=0d0
  do ia=1,natom
     curr_i(:) = curr_i(:) + pp%Zps(Kion(ia)) * system%Velocity(:,ia)
    !curr_i(:) = curr_i(:) - pp%Zps(Kion(ia)) * system%Velocity(:,ia)
  enddo
  curr_i(:) = curr_i(:)/(dble(lg%num(1)*lg%num(2)*lg%num(3))*system%Hvol)

end subroutine calc_current_ion

