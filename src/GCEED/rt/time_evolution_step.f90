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

SUBROUTINE time_evolution_step(lg,mg,ng,system,nspin,info,stencil,srg,srg_ng, &
&   ppn,spsi_in,spsi_out,tpsi1,tpsi2,fg,energy,force,md,ofl)
  use structures
  use salmon_parallel, only: nproc_id_global, nproc_group_global, nproc_group_h, nproc_group_korbital
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use density_matrix, only: calc_density
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
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(inout) :: system
  integer,intent(in) :: nspin
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil),intent(inout) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg,srg_ng
  type(s_pp_nlcc), intent(in) :: ppn
  type(s_orbital),intent(inout) :: spsi_in,spsi_out
  type(s_orbital),intent(inout) :: tpsi1,tpsi2 ! temporary wavefunctions
  type(s_reciprocal_grid) :: fg
  type(s_force),intent(inout) :: force
  type(s_dft_energy) :: energy
  type(s_scalar) :: sVh
  type(s_scalar),allocatable :: srho(:),V_local(:),sVxc(:)
  type(s_md) :: md
  type(s_ofile) :: ofl

  integer :: ix,iy,iz,i1,mm,jj,jspin,n,nn
  integer :: iob,iatom,iik,ik
  real(8) :: rbox1,rbox1q,rbox1q12,rbox1q23,rbox1q31,rbox1e
  complex(8),allocatable :: cmatbox1(:,:,:),cmatbox2(:,:,:)
  real(8) :: absr2
  
  integer :: idensity, idiffDensity, ielf
  real(8) :: rNe, FionE(3,MI)
  
  complex(8),parameter :: zi=(0.d0,1.d0)
  
  complex(8) :: cbox1,cbox2,cbox3
  integer :: is

  character(100) :: comment_line

  logical :: rion_update

  call timer_begin(LOG_CALC_VBOX)
  
  idensity=0
  idiffDensity=1
  ielf=2

  ! for calc_total_energy_periodic
  rion_update = check_rion_update() .or. (itt == Miter_rt+1)
  
  if(iperiodic==3) call init_k_rd(k_rd,ksquare,1,system%primitive_b)
  
  select case(ikind_eext)
    case(0,3,9:12)
      ihpsieff=0
    case(1,2,4,6:8,15)
      ihpsieff=1
  end select
  
  if(iperiodic==0.and.ikind_eext==1) call calcVbox(itt)
  call timer_end(LOG_CALC_VBOX)
  
  
  call timer_begin(LOG_CALC_TIME_PROPAGATION)
!$OMP parallel do private(ik,iob,is,iz,iy,ix) collapse(5)
  do ik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is_array(3),mg%ie_array(3)
      do iy=mg%is_array(2),mg%ie_array(2)
      do ix=mg%is_array(1),mg%ie_array(1)
        spsi_in%zwf(ix,iy,iz,is,iob,ik,1)=zpsi_in(ix,iy,iz,iob+(is-1)*info%numo,ik)
      end do
      end do
      end do
    end do
  end do
  end do
!$OMP parallel do private(ik,iob,is,iz,iy,ix) collapse(5)
  do ik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is_array(3),mg%ie_array(3)
      do iy=mg%is_array(2),mg%ie_array(2)
      do ix=mg%is_array(1),mg%ie_array(1)
        spsi_out%zwf(ix,iy,iz,is,iob,ik,1)=zpsi_out(ix,iy,iz,iob+(is-1)*info%numo,ik)
      end do
      end do
      end do
    end do
  end do
  end do

  !(MD:part1 & update of pseudopotential)
  if(iflag_md==1) then
     call time_evolution_step_md_part1(system,force,md)
     call update_pseudo_rt(itt,system,stencil,lg,ng,fg,ppg,ppg_all,ppn)
  endif

  if(propagator=='etrs')then
    if(iobnum.ge.1)then
      call taylor(mg,nspin,info,lg_sta,lg_end,stencil,srg,spsi_in,spsi_out,tpsi1,   &
                  ppg,vlocal,vbox,num_kpoints_rd,k_rd,zc,ihpsieff)
    end if

!$OMP parallel do private(is,iz,iy,ix) collapse(3)
    do is=1,nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      vloc_t(ix,iy,iz,is)=vlocal(ix,iy,iz,is)
      vloc_new(ix,iy,iz,is) = 3d0*vlocal(ix,iy,iz,is) - 3d0*vloc_old(ix,iy,iz,is,1) + vloc_old(ix,iy,iz,is,2)
      vloc_old(ix,iy,iz,is,2) = vloc_old(ix,iy,iz,is,1)
      vloc_old(ix,iy,iz,is,1) = vlocal(ix,iy,iz,is)
      vlocal(ix,iy,iz,is) = vloc_new(ix,iy,iz,is)
    end do
    end do
    end do
    end do

    if(iperiodic==0.and.ikind_eext==1) call calcVbox(itt+1)
    if(iperiodic==3) call init_k_rd(k_rd,ksquare,4,system%primitive_b)

    if(iobnum.ge.1)then
      call taylor(mg,nspin,info,lg_sta,lg_end,stencil,srg,spsi_out,spsi_in,tpsi1,   &
                  ppg,vlocal,vbox,num_kpoints_rd,k_rd,zc,ihpsieff)
    end if

  else 

    if(iobnum.ge.1)then
      if(mod(itt,2)==1)then
        call taylor(mg,nspin,info,lg_sta,lg_end,stencil,srg,spsi_in,spsi_out,tpsi1,   &
                    ppg,vlocal,vbox,num_kpoints_rd,k_rd,zc,ihpsieff)
      else
        call taylor(mg,nspin,info,lg_sta,lg_end,stencil,srg,spsi_out,spsi_in,tpsi1,   &
                    ppg,vlocal,vbox,num_kpoints_rd,k_rd,zc,ihpsieff)
      end if
    end if
    
  end if
  call timer_end(LOG_CALC_TIME_PROPAGATION)


!$OMP parallel do private(ik,iob,is,iz,iy,ix) collapse(5)
  do ik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is_array(3),mg%ie_array(3)
      do iy=mg%is_array(2),mg%ie_array(2)
      do ix=mg%is_array(1),mg%ie_array(1)
        zpsi_in(ix,iy,iz,iob+(is-1)*info%numo,ik)=spsi_in%zwf(ix,iy,iz,is,iob,ik,1)
      end do
      end do
      end do
    end do
  end do
  end do
!$OMP parallel do private(ik,iob,iz,iy,ix) collapse(5)
  do ik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is_array(3),mg%ie_array(3)
      do iy=mg%is_array(2),mg%ie_array(2)
      do ix=mg%is_array(1),mg%ie_array(1)
        zpsi_out(ix,iy,iz,iob+(is-1)*info%numo,ik)=spsi_out%zwf(ix,iy,iz,is,iob,ik,1)
      end do
      end do
      end do
    end do
  end do
  end do

  allocate(srho(nspin),V_local(nspin),sVxc(nspin))
  do jspin=1,system%nspin
    allocate(srho(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(V_local(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(sVxc(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  end do
  allocate(sVh%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

  if(iperiodic==0)then
    if(ikind_eext==0.and.itt>=2)then
      do jspin=1,system%nspin
        V_local(jspin)%f = Vlocal(:,:,:,jspin)
      end do
      if(ilsda == 1) then
        do jspin=1,system%nspin
          srho(jspin)%f = rho_s(:,:,:,jspin)
          sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
        end do
      else
        srho(1)%f = rho
        sVxc(1)%f = Vxc
      end if
      sVh%f = Vh
      energy%E_xc = Exc
      if(mod(itt,2)==0.or.propagator=='etrs')then
        call calc_eigen_energy(energy,spsi_in,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
      else
        call calc_eigen_energy(energy,spsi_out,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
      end if
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho,sVh,sVxc)
      Etot = energy%E_tot
      call subdip(rNe,2)
    end if
  end if


  call timer_begin(LOG_CALC_RHO)
  if(mod(itt,2)==0.or.propagator=='etrs')then
    call calc_density(srho,spsi_in,info,mg,nspin)
  else
    call calc_density(srho,spsi_out,info,mg,nspin)
  end if
  
  if(ilsda==0)then  
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho(ix,iy,iz)=srho(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(ix,iy,iz,1)=srho(1)%f(ix,iy,iz)
      rho_s(ix,iy,iz,2)=srho(2)%f(ix,iy,iz)
      rho(ix,iy,iz)=srho(1)%f(ix,iy,iz)+srho(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  call timer_end(LOG_CALC_RHO)


  call timer_begin(LOG_CALC_HARTREE)
  if(itt/=1)then
    if(mod(itt,2)==1)then
!$OMP parallel do private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        Vh_stock2(ix,iy,iz)=2.d0*Vh_stock1(ix,iy,iz)-Vh_stock2(ix,iy,iz)
      end do
      end do
      end do
    else
!$OMP parallel do private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        Vh_stock1(ix,iy,iz)=2.d0*Vh_stock2(ix,iy,iz)-Vh_stock1(ix,iy,iz)
      end do
      end do
      end do
    end if
  end if

  
  call Hartree_ns(lg,mg,ng,system%primitive_b,srg_ng,stencil)
  call timer_end(LOG_CALC_HARTREE)


  call timer_begin(LOG_CALC_EXC_COR)
  if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
    call exc_cor_ns(ppn)
  end if
  call timer_end(LOG_CALC_EXC_COR)


  call timer_begin(LOG_CALC_VLOCAL) ! FIXME: wrong name
  call allgatherv_vlocal
  call timer_end(LOG_CALC_VLOCAL)

  do jspin=1,system%nspin
    V_local(jspin)%f = Vlocal(:,:,:,jspin)
  end do
  if(ilsda == 1) then
    do jspin=1,system%nspin
      sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
    end do
  else
    sVxc(1)%f = Vxc
  end if
  sVh%f = Vh
  energy%E_xc = Exc

! result

  if(iperiodic==0)then
    if(ikind_eext/=0.or.(ikind_eext==0.and.itt==itotNtime))then
  
      ihpsieff=0
      if(mod(itt,2)==0.or.propagator=='etrs')then
        call calc_eigen_energy(energy,spsi_in,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
      else
        call calc_eigen_energy(energy,spsi_out,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
      end if
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho,sVh,sVxc)
      Etot = energy%E_tot
      call subdip(rNe,1)
    end if
  end if

  call timer_begin(LOG_CALC_PROJECTION)
  if(iwrite_projection==1.and.mod(itt,itwproj)==0)then
    if(mod(itt,2)==0.or.propagator=='etrs')then
      call projection(mg,zpsi_in)
    else
      call projection(mg,zpsi_out)
    end if
  end if
  call timer_end(LOG_CALC_PROJECTION)


  call timer_begin(LOG_CALC_QUADRUPOLE) ! FIXME: wrong name
  if(iflag_dip2==1)then
    do jj=1,num_dip2
      do i1=1,3
        rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1=rbox1+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2(i1,jj)=rbox1
      end do
    end do
  
    call comm_summation(rbox_array_dip2,rbox_array2_dip2,4*num_dip2,nproc_group_h)

    do jj=1,num_dip2
      Dp2(1:3,itt,jj)=rbox_array2_dip2(1:3,jj)*Hgs(1:3)*Hvol-vecDs2(1:3,jj)
    end do

!------------QUADRUPOLE-start------------

    if(quadrupole=='y')then
      rho_diff(:,:,:) = rho(:,:,:)-rho0(:,:,:)
      do jj=1,num_dip2
        vecR_tmp(:,:,:,:)=vecR(:,:,:,:)
        vecR_tmp(1,:,:,:)=vecR_tmp(1,:,:,:)-dip2center(jj)/Hgs(1)
        do i1=1,3
          rbox1q=0.d0
 !$OMP parallel do reduction( + : rbox1q ) private(absr2,iz,iy,ix)
          do iz=ng_sta(3),ng_end(3)
          do iy=ng_sta(2),ng_end(2)
          do ix=ng_sta(1),ng_end(1)
            absr2=vecR_tmp(1,ix,iy,iz)**2+vecR_tmp(2,ix,iy,iz)**2+vecR_tmp(3,ix,iy,iz)**2
            rbox1q=rbox1q+(3.d0*vecR_tmp(i1,ix,iy,iz)*vecR_tmp(i1,ix,iy,iz)-absr2)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          end do
          end do
          end do
          rbox_array_dip2q(i1,i1,jj)=rbox1q
        end do
      end do
    
      do jj=1,num_dip2
        rbox1q12=0.d0
        rbox1q23=0.d0
        rbox1q31=0.d0
 !$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1q12=rbox1q12+3.d0*vecR_tmp(1,ix,iy,iz)*vecR_tmp(2,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q23=rbox1q23+3.d0*vecR_tmp(2,ix,iy,iz)*vecR_tmp(3,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q31=rbox1q31+3.d0*vecR_tmp(3,ix,iy,iz)*vecR_tmp(1,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2q(1,2,jj)=rbox1q12 ; rbox_array_dip2q(2,1,jj)=rbox1q12
        rbox_array_dip2q(2,3,jj)=rbox1q23 ; rbox_array_dip2q(3,2,jj)=rbox1q23
        rbox_array_dip2q(3,1,jj)=rbox1q31 ; rbox_array_dip2q(1,3,jj)=rbox1q31
      end do

      call comm_summation(rbox_array_dip2q,rbox_array2_dip2q,9*num_dip2,nproc_group_h)
 
      do jj=1,num_dip2
        do i1=1,3
          Qp2(1:3,i1,itt,jj)=rbox_array2_dip2q(1:3,i1,jj)*Hgs(1:3)**2*Hvol
        end do
      end do

    end if

!------------QUADRUPOLE-end--------------
    if(iflag_intelectron==1)then
    do jj=1,num_dip2
        rbox1e=0.d0
!$OMP parallel do reduction( + : rbox1e ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1e=rbox1e+rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2e(jj)=rbox1e
      end do
  
      call comm_summation(rbox_array_dip2e,rbox_array2_dip2e,num_dip2,nproc_group_h)

      do jj=1,num_dip2
        rIe2(itt,jj)=rbox_array2_dip2e(jj)*Hvol
      end do
    end if
  end if
  call timer_end(LOG_CALC_QUADRUPOLE)


  if(iperiodic==3)then
    call subdip(rNe,1)

    call timer_begin(LOG_WRITE_ENERGIES)
    if(iflag_hartree==2)then
      fg%zrhoG_ele = rhoe_G
    else if(iflag_hartree==4)then
      do iz=1,lg_num(3)/NPUZ
      do iy=1,lg_num(2)/NPUY
      do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
        n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
        nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
        fg%zrhoG_ele(nn) = rhoe_G(n)
      enddo
      enddo
      enddo
    end if
    call timer_end(LOG_WRITE_ENERGIES)

    if(mod(itt,2)==0.or.propagator=='etrs')then
      call calc_current(mg,srg,zpsi_in)
      if(itt==itotNtime.or.mod(itt,itcalc_ene)==0)then
        call timer_begin(LOG_CALC_EIGEN_ENERGY)
        call calc_eigen_energy(energy,spsi_in,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
        call timer_end(LOG_CALC_EIGEN_ENERGY)
      end if
    else
      call calc_current(mg,srg,zpsi_out)
      if(itt==1.or.itt==itotNtime.or.mod(itt,itcalc_ene)==0)then
        call timer_begin(LOG_CALC_EIGEN_ENERGY)
        call calc_eigen_energy(energy,spsi_out,tpsi1,tpsi2,system,info,mg,V_local,stencil,srg,ppg)
        call timer_end(LOG_CALC_EIGEN_ENERGY)
      end if
    end if

    if(iflag_md==1) then
      call timer_begin(LOG_CALC_CURRENT_ION)
      call calc_current_ion(system,curr_ion(:,itt))
      call timer_begin(LOG_CALC_CURRENT_ION)
    end if

    call timer_begin(LOG_CALC_TOTAL_ENERGY_PERIODIC)
    call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    call timer_end(LOG_CALC_TOTAL_ENERGY_PERIODIC)

    call timer_begin(LOG_WRITE_ENERGIES)
    rbox1=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction( + : rbox1 )
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rbox1=rbox1+rho(ix,iy,iz)
    end do
    end do
    end do
    rbox1=rbox1*Hvol
    call comm_summation(rbox1,rNe,nproc_group_h)
  !    write(*,'(1x,i7, 3e16.8, f15.8,f18.8,i5,f16.8)')       &
  !      itt, (curr(i1,itt),i1=1,3), Ne, Etot*2d0*Ry,iterVh,dble(cumnum)
    if(comm_is_root(nproc_id_global))then
      write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8)')       &
        itt,dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, energy%E_tot*2d0*Ry
      if(iflag_md==1) then
        write(16,'(f14.8, 6e16.8, f15.8,f18.8)')       &
        dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3),(curr_ion(i1,itt),i1=1,3),rNe, energy%E_tot*2d0*Ry
      else
        write(16,'(f14.8, 3e16.8, f15.8,f18.8)')       &
        dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, energy%E_tot*2d0*Ry
      endif
      write(17,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_tot(i1,itt),i1=1,3)
      write(18,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ext(i1,itt),i1=1,3)
      write(19,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ind(i1,itt),i1=1,3)
    end if
    call timer_end(LOG_WRITE_ENERGIES)
  end if


  call timer_begin(LOG_WRITE_RT_INFOS)

  !(force)
  if(icalcforce==1)then  ! and or rvf flag in future

     !(currently does not work)
     if(iperiodic==3)then
        call get_fourier_grid_G_rt(system,lg,ng,fg)
     endif

     if(mod(itt,2)==0.or.propagator=='etrs')then
        call calc_force_salmon(force,system,pp,fg,info,mg,stencil,srg,ppg,spsi_in) 
     else
        call calc_force_salmon(force,system,pp,fg,info,mg,stencil,srg,ppg,spsi_out)
     end if

     !force on ion directly from field --- should put in calc_force_salmon?
     do iatom=1,MI
        FionE(:,iatom) = pp%Zps(Kion(iatom)) * E_tot(:,itt)
     enddo
     force%F(:,:) = force%F(:,:) + FionE(:,:)

     rforce(:,:) = force%F(:,:)  !test not necessary

  endif

  !(MD: part2)
  if(iflag_md==1) call time_evolution_step_md_part2(system,force,md)


  if(circular=='y')then
    allocate(cmatbox1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    allocate(cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    
!$OMP parallel do private(iz,iy,ix)
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      cmatbox1(ix,iy,iz)=0.d0
    end do
    end do
    end do
    cbox1=0.d0

    do iik=k_sta,k_end
    do iob=1,iobnum
      if(mod(itt,2)==0.or.propagator=='etrs')then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_in(ix,iy,iz,iob,iik)
        end do
        end do
        end do
      else
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_out(ix,iy,iz,iob,iik)
        end do
        end do
        end do
      end if

      call comm_summation(cmatbox1,cmatbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_korbital)
      cbox3=0.d0
      do iz=lg_sta(3),lg_end(3)
      do ix=1,lg_end(1)
        cbox3=cbox3+(conjg(cmatbox2(ix,0,iz))*(cmatbox2(ix,1,iz)-cmatbox2(ix,-1,iz))/(2.d0*Hgs(2)) &
                   -(conjg(cmatbox2(ix,1,iz))-conjg(cmatbox2(ix,-1,iz)))/(2.d0*Hgs(2))*cmatbox2(ix,0,iz))/Hgs(1)/Hgs(2)
      end do
      end do
      cbox1=cbox1+cbox3

      cbox3=0.d0
      do iz=lg_sta(3),lg_end(3)
      do iy=1,lg_end(2)
        cbox3=cbox3-(conjg(cmatbox2(0,iy,iz))*(cmatbox2(1,iy,iz)-cmatbox2(-1,iy,iz))/(2.d0*Hgs(1)) &
                   +(conjg(cmatbox2(1,iy,iz))-conjg(cmatbox2(-1,iy,iz)))/(2.d0*Hgs(1))*cmatbox2(0,iy,iz))/Hgs(1)/Hgs(2)
      end do
      end do
      cbox1=cbox1+cbox3

    end do
    end do

    call comm_summation(cbox1,cbox2,nproc_group_global)

    cumnum=cumnum+cbox2/zi*dt

    deallocate(cmatbox1,cmatbox2) 
  end if

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
     call write_xyz(comment_line,"add","rvf",system,force)
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

  if(iflag_fourier_omega==1)then
    do mm=1,num_fourier_omega
!$OMP parallel do  private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        zalpha2(ix,iy,iz,mm)=zalpha2(ix,iy,iz,mm)   &
                             +exp(zi*fourier_omega(mm)*(itt*dt))*(rho(ix,iy,iz)-rho0(ix,iy,iz)) & 
                             *(1-3*(itt/itotNtime2)**2+2*(itt/itotNtime2)**3)
      end do
      end do
      end do
    end do
  end if

  if(out_dns_rt=='y')then
    if(mod(itt,out_dns_rt_step)==0)then
      call writedns(lg,mg,ng,rho,matbox_m,matbox_m2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,rho0,itt)
    end if
  end if
  if(out_elf_rt=='y')then
    if(mod(itt,out_elf_rt_step)==0)then
      call calcELF
      call writeelf(lg,elf,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,itt)
    end if
  end if
  if(out_estatic_rt=='y')then
    if(mod(itt,out_estatic_rt_step)==0)then
      call calcEstatic
      call writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,matbox_l,matbox_l2,icoo1d,hgs,igc_is,igc_ie,gridcoo,itt)
    end if
  end if
  call timer_end(LOG_WRITE_RT_INFOS)

  call deallocate_scalar(sVh)
  do jspin=1,nspin
    call deallocate_scalar(srho(jspin))
    call deallocate_scalar(V_local(jspin))
    call deallocate_scalar(sVxc(jspin))
  end do
  deallocate(srho,V_local,sVxc)

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
