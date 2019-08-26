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

SUBROUTINE time_evolution_step(lg,mg,ng,system,info,info_field,stencil,srg,srg_ng,srg_ob_1, &
&   ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,sVh,sVxc,sVpsl,dmat,fg,energy,md,ofl, &
&   poisson_cg,j_e,singlescale)
  use structures
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use density_matrix, only: calc_density, calc_density_matrix, calc_current, calc_current_use_dmat, calc_microscopic_current
  use writefield
  use timer
  use inputoutput
  use taylor_sub
  use const, only: umass
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  use read_pslfile_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use salmon_Total_Energy, only: calc_Total_Energy_isolated, calc_Total_Energy_periodic, calc_eigen_energy
  use force_sub, only: calc_force_salmon
  use md_sub, only: time_evolution_step_md_part1,time_evolution_step_md_part2, &
                    update_pseudo_rt
  use print_sub, only: write_xyz,write_rt_data_3d,write_rt_energy_data
  use hpsi_sub, only: update_kvector_nonlocalpt, update_kvector_nonlocalpt_microAc
  use fdtd_coulomb_gauge, only: ls_singlescale, fdtd_singlescale
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(inout) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_field_parallel),intent(in) :: info_field
  type(s_stencil),intent(inout) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg,srg_ng,srg_ob_1
  type(s_pp_nlcc),intent(in)    :: ppn
  type(s_orbital),intent(inout) :: spsi_in,spsi_out
  type(s_orbital),intent(inout) :: tpsi ! temporary wavefunctions
  type(s_scalar), intent(inout) :: srho,srho_s(system%nspin),V_local(system%nspin),sVh,sVxc(system%nspin),sVpsl
  type(s_dmatrix),intent(inout) :: dmat
  type(s_poisson_cg),intent(in) :: poisson_cg
  type(s_vector) :: j_e ! microscopic electron number current density
  type(ls_singlescale) :: singlescale
  type(s_reciprocal_grid) :: fg
  type(s_dft_energy) :: energy
  type(s_md) :: md
  type(s_ofile) :: ofl

  integer :: ix,iy,iz,nspin
  integer :: iatom,ik
  integer :: idensity, idiffDensity, ielf
  real(8) :: rNe, FionE(3,MI)
  real(8) :: curr_tmp(3,2)
  complex(8),parameter :: zi=(0.d0,1.d0)
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
  rion_update = check_rion_update() .or. (itt == Miter_rt+1)
  
  select case(ikind_eext)
    case(0,3,9:12)
      ihpsieff=0
    case(1,2,4,6:8,15)
      ihpsieff=1
  end select
  
  select case(iperiodic)
  case(0)
    if(ikind_eext==1) call calcVbox(itt)
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
      call calc_vecAc(system%vec_Ac,1)
      do ik=info%ik_s,info%ik_e
        stencil%vec_kAc(:,ik) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
      end do
      call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
    end if
  end select

  call timer_end(LOG_CALC_VBOX)

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

  !(MD:part1 & update of pseudopotential)
  if(iflag_md==1) then
     call time_evolution_step_md_part1(system,md)
     call update_pseudo_rt(itt,info,system,stencil,lg,ng,fg,ppg,ppg_all,ppn)
     sVpsl%f = Vpsl ! future work: remove Vpsl
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
      if(ikind_eext==1) call calcVbox(itt+1)
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
      ! future work: etrs for single-scale Maxwell-TDDFT
!        call calc_density_matrix(nspin,info,mg,srg,tpsi,dmat)
!        call calc_microscopic_current(nspin,mg,stencil,info,tpsi,dmat,j_e)
!        singlescale%E_electron = energy%E_tot
!        call fdtd_singlescale(itt-1,lg,mg,ng,system%hgs,srho,sVh,j_e,srg_ng,system%Ac_micro,system%div_Ac,singlescale)
!        call update_kvector_nonlocalpt_microAc(info%ik_s,info%ik_e,system,ppg)
      else
        call calc_vecAc(system%vec_Ac,4)
        do ik=info%ik_s,info%ik_e
          stencil%vec_kAc(:,ik) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
        end do
        call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
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
  call calc_density(srho_s,spsi_out,info,mg,nspin)

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
      Vh_stock2 = 2.d0*Vh_stock1 - Vh_stock2
      sVh%f = Vh_stock2
    else
      Vh_stock1 = 2.d0*Vh_stock2 - Vh_stock1
      sVh%f = Vh_stock1
    end if
  end if
  call Hartree_ns(lg,mg,ng,info_field,system,poisson_cg,srg_ng,stencil,srho,sVh,fg)
  if(iperiodic==0 .and. itt/=1)then
    if(mod(itt,2)==1)then
      Vh_stock2 = sVh%f
    else
      Vh_stock1 = sVh%f
    end if
  end if
  call timer_end(LOG_CALC_HARTREE)

  call timer_begin(LOG_CALC_EXC_COR)
  if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_d_o_mul*nproc_d_g_mul_dm))then
    call exc_cor_ns(ng, srg_ng, system%nspin, srho_s, ppn, sVxc, energy%E_xc)
  end if
  call timer_end(LOG_CALC_EXC_COR)

  call timer_begin(LOG_CALC_VLOCAL) ! FIXME: wrong name
  call allgatherv_vlocal(ng,info,system%nspin,sVh,sVpsl,sVxc,V_local)
  call timer_end(LOG_CALC_VLOCAL)

! result

  call timer_begin(LOG_CALC_PROJECTION)
  if(iwrite_projection==1.and.mod(itt,itwproj)==0)then
    call projection(mg,info,zpsi_out)
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

    if(if_use_dmat) call calc_density_matrix(nspin,info,mg,srg,spsi_out,dmat)

    call timer_begin(LOG_CALC_CURRENT)
    if(if_use_dmat) then
      call calc_current_use_dmat(nspin,system%ngrid,mg,stencil,info,spsi_out,ppg,dmat,curr_tmp(1:3,1:nspin))
    else
      call calc_current(nspin,system%ngrid,mg,stencil,info,srg,spsi_out,ppg,curr_tmp(1:3,1:nspin))
    end if
    call calc_emfields(nspin,curr_tmp)
    call timer_end(LOG_CALC_CURRENT)

    if(iflag_md==1) then
      call timer_begin(LOG_CALC_CURRENT_ION)
      call calc_current_ion(system,curr_ion(:,itt))
      call timer_end(LOG_CALC_CURRENT_ION)
    end if

    call timer_begin(LOG_CALC_TOTAL_ENERGY_PERIODIC)
    call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    Etot = energy%E_tot
    call timer_end(LOG_CALC_TOTAL_ENERGY_PERIODIC)

    if(use_singlescale=='y') then
      call calc_microscopic_current(nspin,mg,stencil,info,spsi_out,dmat,j_e)
      singlescale%E_electron = energy%E_tot
      call fdtd_singlescale(itt,lg,mg,ng,system%hgs,srho,sVh,j_e,srg_ng,system%Ac_micro,system%div_Ac,singlescale)
      call update_kvector_nonlocalpt_microAc(info%ik_s,info%ik_e,system,ppg)
    end if

  end select

  call timer_begin(LOG_WRITE_ENERGIES)
  call subdip(ng,srho,rNe)
  call timer_end(LOG_WRITE_ENERGIES)

  call timer_begin(LOG_WRITE_RT_INFOS)

  !(force)
  if(icalcforce==1)then  ! and or rvf flag in future

     !(currently does not work)
     if(iperiodic==3)then
        call get_fourier_grid_G_rt(system,lg,ng,fg)
     endif

     call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi_out)

     !force on ion directly from field --- should put in calc_force_salmon?
     do iatom=1,MI
        FionE(:,iatom) = pp%Zps(Kion(iatom)) * E_tot(:,itt)
     enddo
     system%Force(:,:) = system%Force(:,:) + FionE(:,:)

     rforce(:,:) = system%Force(:,:)  !test not necessary

  endif

  !(MD: part2)
  if(iflag_md==1) call time_evolution_step_md_part2(system,md)

  ! Output 
  !(Export to SYSname_trj.xyz)
  if( icalcforce==1 .and. mod(itt,out_rvf_rt_step)==0 )then
     write(comment_line,10) itt, itt*dt
10   format("#rt   step=",i8,"   time=",e16.6)
     if(iflag_md==1) write(comment_line,11) trim(comment_line),md%Temperature
11   format(a,"   T=",f12.4)
     if(iflag_md==1 .and. ensemble=="NVT" .and. thermostat=="nose-hoover") &
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
        call write_rt_data_3d(itt,ofl,iflag_md,dt)
     end select

     !(Export to SYSname_rt_energy.data)
     call write_rt_energy_data(itt,ofl,iflag_md,dt,energy,md)

  endif

  if(yn_out_dns_rt=='y')then
    if(mod(itt,out_dns_rt_step)==0)then
      call writedns(lg,mg,ng,srho%f,matbox_m,matbox_m2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,rho0,itt)
    end if
  end if
  if(yn_out_elf_rt=='y')then
    if(mod(itt,out_elf_rt_step)==0)then
      call calcELF(mg,ng,info,srho,itt,srg_ob_1)
      call writeelf(lg,elf,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,itt)
    end if
  end if
  if(yn_out_estatic_rt=='y')then
    if(mod(itt,out_estatic_rt_step)==0)then
      call calcEstatic(ng, info, sVh, srg_ng)
      call writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,matbox_l,matbox_l2,icoo1d,hgs,igc_is,igc_ie,gridcoo,itt)
    end if
  end if
  call timer_end(LOG_WRITE_RT_INFOS)

  return

END SUBROUTINE time_evolution_step

subroutine get_fourier_grid_G_rt(system,lg,ng,fg)
  use salmon_global, only: nelem
  use structures, only: s_dft_system, s_reciprocal_grid, s_rgrid
  use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global
  use scf_data
  use allocate_psl_sub
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_reciprocal_grid) :: fg

  integer :: jj,ix,iy,iz,n,nn

  if(allocated(fg%Gx))       deallocate(fg%Gx,fg%Gy,fg%Gz)
  if(allocated(fg%zrhoG_ion)) deallocate(fg%zrhoG_ion,fg%zrhoG_ele,fg%zdVG_ion)

  jj = system%ngrid/nproc_size_global
  fg%ig_s = nproc_id_global*jj+1
  fg%ig_e = (nproc_id_global+1)*jj
  if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
  fg%icomm_G = nproc_group_global
  fg%ng = system%ngrid
  allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))
  allocate(fg%zrhoG_ion(fg%ng),fg%zrhoG_ele(fg%ng),fg%zdVG_ion(fg%ng,nelem))
  if(iflag_hartree==2)then
     fg%iGzero = nGzero
     fg%Gx = Gx
     fg%Gy = Gy
     fg%Gz = Gz
     fg%zrhoG_ion = rhoion_G
     fg%zdVG_ion = dVloc_G
  else if(iflag_hartree==4)then
     fg%iGzero = 1
     fg%Gx = 0.d0
     fg%Gy = 0.d0
     fg%Gz = 0.d0
     fg%zrhoG_ion = 0.d0
     fg%zdVG_ion = 0.d0
     do iz=1,lg_num(3)/NPUZ
     do iy=1,lg_num(2)/NPUY
     do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
        n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
        nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
        fg%Gx(nn) = Gx(n)
        fg%Gy(nn) = Gy(n)
        fg%Gz(nn) = Gz(n)
        fg%zrhoG_ion(nn) = rhoion_G(n)
        fg%zdVG_ion(nn,:) = dVloc_G(n,:)
     enddo
     enddo
     enddo
  end if

end subroutine get_fourier_grid_G_rt
