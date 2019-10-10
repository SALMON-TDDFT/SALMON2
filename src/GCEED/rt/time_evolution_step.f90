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
!=======================================================================
!=======================================================================

SUBROUTINE time_evolution_step(lg,mg,ng,system,info,info_field,stencil,xc_func,srg,srg_ng, &
&   pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,sVh,sVh_stock1,sVh_stock2,sVxc,sVpsl,dmat,fg,energy,md,ofl, &
&   poisson,j_e,singlescale)
  use structures
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use density_matrix, only: calc_density, calc_density_matrix, calc_current, calc_current_use_dmat, calc_microscopic_current
  use writefield
  use timer
  use inputoutput
  use taylor_sub
  use const, only: umass
  use scf_data
  use allocate_mat_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use hartree_sub, only: hartree
  use salmon_Total_Energy, only: calc_Total_Energy_isolated, calc_Total_Energy_periodic, calc_eigen_energy
  use force_sub, only: calc_force_salmon
  use md_sub, only: time_evolution_step_md_part1,time_evolution_step_md_part2, &
                    update_pseudo_rt
  use write_sub
  use hamiltonian, only: update_kvector_nonlocalpt, update_kvector_nonlocalpt_microAc, allgatherv_vlocal
  use fdtd_coulomb_gauge, only: ls_singlescale, fdtd_singlescale
  use salmon_pp, only: calc_nlcc !test hoge
  use salmon_xc
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(inout) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_field_parallel),intent(in) :: info_field
  type(s_stencil),intent(inout) :: stencil
  type(s_xc_functional),intent(in) :: xc_func
  type(s_sendrecv_grid),intent(inout) :: srg,srg_ng
  type(s_pp_info),intent(inout) :: pp
  type(s_pp_grid) :: ppg
!  type(s_pp_nlcc),intent(in)    :: ppn
  type(s_pp_nlcc),intent(inout)    :: ppn !hoge test
  type(s_orbital),intent(inout) :: spsi_in,spsi_out
  type(s_orbital),intent(inout) :: tpsi ! temporary wavefunctions
  type(s_scalar), intent(inout) :: srho,srho_s(system%nspin),V_local(system%nspin),sVh,sVxc(system%nspin),sVpsl
  type(s_scalar), intent(inout) :: sVh_stock1,sVh_stock2
  type(s_dmatrix),intent(inout) :: dmat
  type(s_poisson),intent(inout) :: poisson
  type(s_vector) :: j_e ! microscopic electron number current density
  type(ls_singlescale) :: singlescale
  type(s_reciprocal_grid) :: fg
  type(s_dft_energy) :: energy
  type(s_md) :: md
  type(s_ofile) :: ofl

  integer :: ix,iy,iz,nspin
  integer :: iatom
  integer :: idensity, idiffDensity, ielf
  real(8) :: rNe, FionE(3,MI)
  real(8) :: curr_tmp(3,2)
  integer :: is
  character(100) :: comment_line
  logical :: rion_update,if_use_dmat

  nspin = system%nspin

  call timer_begin(LOG_CALC_VBOX)
  
  idensity=0
  idiffDensity=1
  ielf=2
  if_use_dmat = (use_singlescale=='y') ! .or. if_metaGGA ! (future work)

  ! for calc_total_energy_periodic
  if(yn_md=='y') then
     rion_update = .true.
  else
     rion_update = check_rion_update() .or. (itt == Miter_rt+1)
  endif

  select case(ikind_eext)
    case(0,3,9:12)
      ihpsieff=0
    case(1,2,4,6:8,15)
      ihpsieff=1
  end select
  
  select case(iperiodic)
  case(0)
    if(ikind_eext==1) call calcVbox(lg,itt)
    if(ihpsieff==1) then
!$OMP parallel do collapse(3) private(is,iz,iy,ix)
      do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        V_local(is)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz) + vbox(ix,iy,iz)
      end do
      end do
      end do
      end do
    end if
  case(3)
    if(use_singlescale=='n') then
      system%vec_Ac(1:3) = A_ext(1:3,itt) + A_ind(1:3,itt)
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
    end if
  end select

  call timer_end(LOG_CALC_VBOX)

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

  !(MD:part1 & update of pseudopotential)
  if(yn_md=='y') then
     call time_evolution_step_md_part1(itt,system,md)
     call update_pseudo_rt(itt,info,info_field,system,stencil,lg,mg,ng,poisson,fg,pp,ppg,ppn,sVpsl)
  endif

  if(propagator=='etrs')then
    if(iobnum.ge.1)then
    ! spsi_in --> tpsi, (spsi_out = working array)
      call taylor(mg,system,info,stencil,srg,spsi_in,tpsi,spsi_out,ppg,V_local,zc)
    end if

!$OMP parallel do private(is,iz,iy,ix) collapse(3)
    do is=1,nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      vloc_t(ix,iy,iz,is) = V_local(is)%f(ix,iy,iz)
      vloc_new(ix,iy,iz,is) = 3d0*V_local(is)%f(ix,iy,iz) - 3d0*vloc_old(ix,iy,iz,is,1) + vloc_old(ix,iy,iz,is,2)
      vloc_old(ix,iy,iz,is,2) = vloc_old(ix,iy,iz,is,1)
      vloc_old(ix,iy,iz,is,1) = V_local(is)%f(ix,iy,iz)
      V_local(is)%f(ix,iy,iz) = vloc_new(ix,iy,iz,is)
    end do
    end do
    end do
    end do

    select case(iperiodic)
    case(0)
      if(ikind_eext==1) call calcVbox(lg,itt+1)
      if(ihpsieff==1)then
  !$OMP parallel do collapse(3) private(is,iz,iy,ix)
        do is=1,nspin
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          V_local(is)%f(ix,iy,iz) = V_local(is)%f(ix,iy,iz) + vbox(ix,iy,iz)
        end do
        end do
        end do
        end do
      end if
    case(3)
      if(use_singlescale=='y') then
        stop "etrs mode for single-scale Maxwell-TDDFT is not implemented"
      else
        system%vec_Ac(1:3) = A_ext(1:3,itt+1) + A_ind(1:3,itt+1)
        call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
      end if
    end select

    if(iobnum.ge.1)then
    ! tpsi --> spsi_out (spsi_in = working array)
      call taylor(mg,system,info,stencil,srg,tpsi,spsi_out,spsi_in,ppg,V_local,zc)
    end if

  else 

    if(iobnum.ge.1)then
    ! spsi_in --> spsi_out (tpsi = working array)
      call taylor(mg,system,info,stencil,srg,spsi_in,spsi_out,tpsi,ppg,V_local,zc)
    end if
    
  end if
  call timer_end(LOG_CALC_TIME_PROPAGATION)

  call timer_begin(LOG_CALC_RHO)
  call calc_density(system,srho_s,spsi_out,info,mg)

  if(nspin==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      srho%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else if(nspin==2)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      srho%f(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)+srho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  call timer_end(LOG_CALC_RHO)

  call timer_begin(LOG_CALC_HARTREE)
  if(iperiodic==0 .and. itt/=1)then
    if(mod(itt,2)==1)then
      sVh_stock2%f = 2.d0*sVh_stock1%f - sVh_stock2%f
      sVh%f = sVh_stock2%f
    else
      sVh_stock1%f = 2.d0*sVh_stock2%f - sVh_stock1%f
      sVh%f = sVh_stock1%f
    end if
  end if
  call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
  if(iperiodic==0 .and. itt/=1)then
    if(mod(itt,2)==1)then
      sVh_stock2%f = sVh%f
    else
      sVh_stock1%f = sVh%f
    end if
  end if
  call timer_end(LOG_CALC_HARTREE)

  call timer_begin(LOG_CALC_EXC_COR)
  call exchange_correlation(system,xc_func,ng,srg_ng,srho_s,ppn,info_field%icomm_all,sVxc,energy%E_xc)
  call timer_end(LOG_CALC_EXC_COR)

  call timer_begin(LOG_CALC_VLOCAL) ! FIXME: wrong name
  call allgatherv_vlocal(ng,mg,info_field,system%nspin,sVh,sVpsl,sVxc,V_local)
  call timer_end(LOG_CALC_VLOCAL)

! result

  call timer_begin(LOG_CALC_PROJECTION)
  if(iwrite_projection==1.and.mod(itt,itwproj)==0)then
    call projection(itt,mg,system,info,spsi_out,tpsi) ! tpsi must be GS orbital (future work)
  end if
  call timer_end(LOG_CALC_PROJECTION)

  if(itt==1.or.itt==itotNtime.or.mod(itt,itcalc_ene)==0)then
    call timer_begin(LOG_CALC_EIGEN_ENERGY)
    ! tpsi,spsi_in = working arrays
    call calc_eigen_energy(energy,spsi_out,tpsi,spsi_in,system,info,mg,V_local,stencil,srg,ppg)
    call timer_end(LOG_CALC_EIGEN_ENERGY)
  end if

  select case(iperiodic)
  case(0)

    call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
    Etot = energy%E_tot

  case(3)

    if(if_use_dmat) call calc_density_matrix(system,info,mg,srg,spsi_out,dmat)

    call timer_begin(LOG_CALC_CURRENT)
    if(if_use_dmat) then
      call calc_current_use_dmat(system,mg,stencil,info,spsi_out,ppg,dmat,curr_tmp(1:3,1:nspin))
    else
      call calc_current(system,mg,stencil,info,srg,spsi_out,ppg,curr_tmp(1:3,1:nspin))
    end if
    call calc_emfields(nspin,curr_tmp)
    call timer_end(LOG_CALC_CURRENT)

    if(yn_md=='y') then
      call timer_begin(LOG_CALC_CURRENT_ION)
      call calc_current_ion(lg,system,pp,curr_ion(:,itt))
      call timer_end(LOG_CALC_CURRENT_ION)
    end if

    call timer_begin(LOG_CALC_TOTAL_ENERGY_PERIODIC)
    call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    Etot = energy%E_tot
    call timer_end(LOG_CALC_TOTAL_ENERGY_PERIODIC)

    if(use_singlescale=='y') then
      call calc_microscopic_current(system,mg,stencil,info,spsi_out,dmat,j_e)
      singlescale%E_electron = energy%E_tot
      call fdtd_singlescale(itt,info%icomm_rko,lg,mg,ng,system%hgs,srho, &
      & sVh,j_e,srg_ng,system%Ac_micro,system%div_Ac,singlescale)
      call update_kvector_nonlocalpt_microAc(info%ik_s,info%ik_e,system,ppg)
    end if

  end select

  call timer_begin(LOG_WRITE_ENERGIES)
  call subdip(ng,srho,rNe,poisson)
  call timer_end(LOG_WRITE_ENERGIES)

  call timer_begin(LOG_WRITE_RT_INFOS)

  !(force)
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then  ! and or rvf flag in future

     call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi_out)

     !force on ion directly from field --- should put in calc_force_salmon?
     do iatom=1,MI
        FionE(:,iatom) = pp%Zps(Kion(iatom)) * E_tot(:,itt)
     enddo
     system%Force(:,:) = system%Force(:,:) + FionE(:,:)

  endif

  !(MD: part2)
  if(yn_md=='y') call time_evolution_step_md_part2(system,md)

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
     case(3)
        call write_rt_data_3d(itt,ofl,dt)
     end select

     !(Export to SYSname_rt_energy.data)
     call write_rt_energy_data(itt,ofl,dt,energy,md)

  endif

  if(yn_out_dns_rt=='y')then
    if(mod(itt,out_dns_rt_step)==0)then
      call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,hgs,iscfrt,rho0,itt)
    end if
  end if
  if(yn_out_elf_rt=='y')then
    if(mod(itt,out_elf_rt_step)==0)then
      call write_elf(iscfrt,itt,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi_out)
    end if
  end if
  if(yn_out_estatic_rt=='y')then
    if(mod(itt,out_estatic_rt_step)==0)then
      call calcEstatic(lg, ng, info, sVh, srg_ng)
      call writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,matbox_l,matbox_l2,hgs,itt)
    end if
  end if
  call timer_end(LOG_WRITE_RT_INFOS)

  return

END SUBROUTINE time_evolution_step

subroutine calc_current_ion(lg,system,pp,j_ion)
  use structures
  use salmon_global, only: MI,Kion
  use scf_data, only: Hvol
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_dft_system) :: system
  type(s_pp_info) :: pp
  integer :: ia
  real(8) :: j_ion(3)

  !AY memo
  !current of ion: defined by positive charge-->minus sign
  !This is NOT matter current NOR electric current.... strange definition....
  !This is defined so as to the total electric current = -(curr + curr_ion)
  !Should change this ion current but if you change,
  !please change all part in ARTED, multiscale .....
  j_ion(:)=0d0
  do ia=1,MI
     j_ion(:) = j_ion(:) - pp%Zps(Kion(ia)) * system%Velocity(:,ia)
  enddo
  j_ion(:) = j_ion(:)/(dble(lg%num(1)*lg%num(2)*lg%num(3))*Hvol)

end subroutine calc_current_ion

