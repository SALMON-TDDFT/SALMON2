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
MODULE global_variables_rt
use inputoutput
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
implicit none

integer       :: Ntime
real(8)       :: debye2au   ! [D]  -> [a.u.] 
integer       :: iii
real(8), allocatable :: alpha2(:,:,:,:)

END MODULE global_variables_rt

!=======================================================================

subroutine main_tddft
use structures
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use global_variables_rt
use write_sub, only: write_xyz,write_rt_data_3d,write_rt_energy_data
use code_optimization
use initialization_sub
use input_pp_sub
use prep_pp_sub
use density_matrix, only: calc_density
use writefield
use salmon_pp, only: calc_nlcc
use force_sub, only: calc_force_salmon
use hamiltonian
use md_sub, only: init_md
use fdtd_coulomb_gauge, only: ls_singlescale
use read_write_restart_rt_sub, only: write_checkpoint_rt
use taylor_sub, only: taylor_coe
use checkpoint_restart_sub
use hartree_sub, only: hartree
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_dft_system)  :: system
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_mixing) :: mixing
type(s_scalar) :: sVpsl
type(s_scalar) :: srho,sVh,sVh_stock1,sVh_stock2
type(s_scalar),allocatable :: srho_s(:),V_local(:),sVxc(:)
type(s_dmatrix) :: dmat
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_ng
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_vector)  :: j_e ! microscopic electron number current density
type(ls_singlescale) :: singlescale
type(s_ofile) :: ofile

integer :: iob, i1,i2,i3, iik,jspin
integer :: idensity, idiffDensity
integer :: jj,nn, iene, nntime, ix,iy,iz
real(8) :: rbox_array(10), rbox_array2(10)
real(8),allocatable :: alpha_R(:,:),   alpha_I(:,:) 
real(8),allocatable :: alphaq_R(:,:,:),alphaq_I(:,:,:)
real(8),allocatable :: R1(:,:,:), Sf(:)
real(8),allocatable :: tfourier_integrand(:,:)
complex(8),parameter :: zi=(0.d0,1.d0)
character(10) :: fileLaser
character(100):: comment_line


call timer_begin(LOG_TOTAL)

call timer_begin(LOG_INIT_RT)
call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=2
OC=0

iwdenstep=30 
denplane='xy'
idensum=0
posplane=0.d0

call convert_input_rt(Ntime)
mixing%num_rho_stock=21

call set_filename

if(comm_is_root(nproc_id_global))then
  write(*,*)
  write(*,*) "Total time step      =",Ntime
  write(*,*) "Time step[fs]        =",dt*au_time_fs
  write(*,*) "Field strength[?]    =",Fst
  write(*,*) "Energy range         =",Nenergy
  write(*,*) "Energy resolution[eV]=",dE*au_energy_ev
  write(*,*) "ikind_eext is         ", ikind_eext
  write(*,*) "Step for writing dens=", iwdenstep
  write(*,*) "Plane showing density=", denplane
  write(*,*) "idensum              =", idensum 
  if(idensum==0) write(*,*) "Position of the plane=", posplane
  select case (ikind_eext)
    case(1,6,7,8,15)
      write(*,20) "Laser frequency     =",romega*au_energy_ev,"[eV]"
      write(*,21) "Pulse width of laser=",pulse_T*au_time_fs, "[fs]"
      write(*,22) "Laser intensity     =",rlaser_I,       "[W/cm^2]"
      write(*,22) "tau                 =",tau*au_time_fs, "[fs]"
    case(4,12)
      write(*,23) "Laser frequency     =",romega2(1:2)*au_energy_ev,"[eV]"
      write(*,24) "Pulse width of laser=",pulse_T2(1:2)*au_time_fs, "[fs]"
      write(*,25) "Laser intensity     =",rlaser_I2(1:2),       "[W/cm^2]"
      write(*,21) "delay time          =",delay*au_time_fs,     "[fs]"
      write(*,26) "rcycle              =",rcycle
  end select
20 format(a21,f5.2, a4)
21 format(a21,f16.8,a4)
22 format(a21,e16.8,a8)
23 format(a21,2f5.2, a4)
24 format(a21,2f16.8,a4)
25 format(a21,2e16.8,a8)
26 format(a21,f16.8)

end if

debye2au = 0.393428d0

select case (ikind_eext)
  case(0,10) ; Fst=Fst !/5.14223d1
end select
dE=dE !/2d0/Ry 
dt=dt !*fs2eVinv*2.d0*Ry!a.u. ! 1[fs] = 1.51925 [1/eV]  !2.d0*Ry*1.51925d0

if(idensum==0) posplane=posplane/a_B 

select case (ikind_eext)
  case(1)
    if(I_wcm2_1>=1.d-12)then
      E_amplitude1=sqrt(I_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    if(I_wcm2_2>=1.d-12)then
      E_amplitude2=sqrt(I_wcm2_2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    else
      if(abs(E_amplitude2)<=1.d-12)then
        E_amplitude2=0.d0
      end if
    end if
end select

call timer_end(LOG_INIT_RT)

call timer_begin(LOG_READ_GS_DATA)

! +----------------+
! | initialization |
! +----------------+

call init_dft(nproc_group_global,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)

call init_code_optimization
call old_mesh(lg,mg,ng,system,info) ! future work: remove this line
call allocate_mat(ng) ! future work: remove this line

call allocate_scalar(mg,srho)
call allocate_scalar(mg,sVh)
call allocate_scalar(mg,sVh_stock1)
call allocate_scalar(mg,sVh_stock2)
call allocate_scalar(mg,sVpsl)
allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))
do jspin=1,system%nspin
  call allocate_scalar(mg,srho_s(jspin))
  call allocate_scalar(mg,V_local(jspin))
  call allocate_scalar(mg,sVxc(jspin))
end do
allocate(ppg%Vpsl_atom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),natom))
call read_pslfile(system,pp,ppg)
call init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)

call allocate_orbital_complex(system%nspin,mg,info,spsi_in)
call allocate_orbital_complex(system%nspin,mg,info,spsi_out)
call allocate_orbital_complex(system%nspin,mg,info,tpsi)
call allocate_dmatrix(system%nspin,mg,info,dmat)

call read_bin(lg,mg,ng,system,info,spsi_in,mixing,sVh_stock1,sVh_stock2,miter_rt)
if(read_rt_wfn_k=='n') miter_rt=0

call calc_nlcc(pp, system, mg, ppn)
if (comm_is_root(nproc_id_global)) then
  write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
  write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
end if

call calc_density(system,srho_s,spsi_in,info,mg)
srho%f = 0d0
do jspin=1,system%nspin
   srho%f = srho%f + srho_s(jspin)%f
end do
if(read_rt_wfn_k=='y')then
  sVh%f = 2.d0*sVh_stock1%f - sVh_stock2%f
  sVh_stock2%f = sVh_stock1%f
end if
call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
call exchange_correlation(system,xc_func,ng,srg_ng,srho_s,ppn,info_field%icomm_all,sVxc,energy%E_xc)
call allgatherv_vlocal(ng,mg,info_field,system%nspin,sVh,sVpsl,sVxc,V_local)
if(read_rt_wfn_k=='y')then
  sVh_stock1%f=sVh%f
else if(read_rt_wfn_k=='n')then
  sVh_stock1%f=sVh%f
  sVh_stock2%f=sVh%f
end if

allocate(energy%esp(system%no,system%nk,system%nspin))

call timer_end(LOG_READ_GS_DATA)

! +-------------+
! | old fashion |
! +-------------+

call timer_begin(LOG_READ_RT_DATA)
allocate( Ex_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 
allocate( Ey_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 
allocate( Ez_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
   Ex_static(ix,iy,iz)=0.d0; Ey_static(ix,iy,iz)=0.d0; Ez_static(ix,iy,iz)=0.d0
end do
end do
end do

if(IC_rt==0) then
  allocate( rIe(0:Ntime) )
  allocate( Dp(3,0:Ntime) )
  allocate( Qp(3,3,0:Ntime) )
  allocate( tene(0:Ntime) )
  call initA(Ntime)
  itotNtime=Ntime
  Miter_rt=0
end if
call timer_end(LOG_READ_RT_DATA)


call timer_begin(LOG_INIT_RT)
allocate( alpha_R(3,0:Nenergy), alpha_I(3,0:Nenergy), Sf(3) )
allocate( alphaq_R(3,3,0:Nenergy), alphaq_I(3,3,0:Nenergy) )

ntmg=1
! 'Hartree' parameter

poisson%iterVh = 0        ! Iteration counter


if(comm_is_root(nproc_id_global))then
  write(*, *) 
  write(*, *) "dip2boundary", dip2boundary(1), dip2boundary(2)
  write(*, *) "dip2center", dip2center(1), dip2center(2)
  write(*, *) "dip2boundary[A]", dip2boundary(1)*a_B, dip2boundary(2)*a_B
  write(*, *) "dip2center[A]", dip2center(1)*a_B, dip2center(2)*a_B
  write(*, *) 
end if
call timer_end(LOG_INIT_RT)

!Open output files and print header parts (Please move and put this kinds here!)
!(output file names)
if(comm_is_root(nproc_id_global))then
   write(ofl%file_rt_data,"(2A,'_rt.data')") trim(base_directory),trim(SYSname)
   write(ofl%file_rt_energy_data,"(2A,'_rt_energy.data')") trim(base_directory),trim(SYSname)
endif
call comm_bcast(ofl%file_rt_data,       nproc_group_global)
call comm_bcast(ofl%file_rt_energy_data,nproc_group_global)

!(write header)
if(comm_is_root(nproc_id_global))then
  !(header of standard output)
  select case(iperiodic)
  case(0)
    write(*,'(1x,a10,a11,a48,a15,a18,a10)') &
                "time-step ", "time[fs]",   &
                "Dipole moment(xyz)[A]"     &
               ,"electrons", "Total energy[eV]", "iterVh"
  case(3)
    write(*,'(1x,a10,a11,a48,a15,a18)')   &
                "time-step", "time[fs] ", &
                "Current(xyz)[a.u.]",     &
                "electrons", "Total energy[eV] "
  end select
  write(*,'("#",7("----------"))')

  !(header of SYSname_rt.data)
  select case(iperiodic)
  case(0)!; call write_rt_data_0d --- make in future
  case(3) ; call write_rt_data_3d(-1,ofl,dt)
  end select

  !(header of SYSname_rt_energy.data)
  call write_rt_energy_data(-1,ofl,dt,energy,md)

  !(header in SYSname_proj.data)
  if(iwrite_projection==1)then
    open(41,file=file_Projection)
    write(41,'("#",5X,"time[fs]",4("    projection"))')
    write(41,'("#",13x,4("  orbital",i5))') (iwrite_projection_ob(jj),jj=1,4)
    write(41,'("#",13x,4("        k",i5))') (iwrite_projection_k(jj), jj=1,4)
    write(41,'("#",7("----------"))')
  end if
end if

if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
   call write_xyz(comment_line,"new","rvf",system)
endif

!---------------------------- time-evolution

call timer_begin(LOG_INIT_TIME_PROPAGATION)

if(comm_is_root(nproc_id_global).and.iperiodic==3) then
  open(16,file="current.data")
  open(17,file="Etot.data")
  open(18,file="Eext.data")
  open(19,file="Eind.data")
end if

idensity=0
idiffDensity=1
fileLaser= "laser.out"

allocate( R1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2), &
                               lg%is(3):lg%ie(3)))

!if(ikind_eext.ne.0)then
allocate( Vbox(lg%is(1)-Nd:lg%ie(1)+Nd,lg%is(2)-Nd:lg%ie(2)+Nd, &
                                       lg%is(3)-Nd:lg%ie(3)+Nd))
!endif

allocate(rho0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

! External Field Direction
select case (ikind_eext)
  case(1)
    if(yn_local_field=='y')then
      rlaser_center(1:3)=(rlaserbound_sta(1:3)+rlaserbound_end(1:3))/2.d0
      do jj=1,3
        select case(mod(lg%num(jj),2))
          case(1)
            ilasbound_sta(jj)=nint(rlaserbound_sta(jj)/Hgs(jj))
            ilasbound_end(jj)=nint(rlaserbound_end(jj)/Hgs(jj))
          case(0)
            ilasbound_sta(jj)=nint(rlaserbound_sta(jj)/Hgs(jj)+0.5d0)
            ilasbound_end(jj)=nint(rlaserbound_end(jj)/Hgs(jj)+0.5d0)
        end select
      end do
    else
      rlaser_center(1:3)=0.d0
    end if
end select

select case(mod(lg%num(1),2))
  case(1)
    do i1=lg%is(1),lg%ie(1)
       vecR(1,i1,:,:)=dble(i1)-rlaser_center(1)/Hgs(1)
    end do
  case(0)
    do i1=lg%is(1),lg%ie(1)
       vecR(1,i1,:,:)=dble(i1)-0.5d0-rlaser_center(1)/Hgs(1)
    end do
end select

select case(mod(lg%num(2),2))
  case(1)
    do i2=lg%is(2),lg%ie(2)
       vecR(2,:,i2,:)=dble(i2)-rlaser_center(2)/Hgs(2)
    end do
  case(0)
    do i2=lg%is(2),lg%ie(2)
       vecR(2,:,i2,:)=dble(i2)-0.5d0-rlaser_center(2)/Hgs(2)
    end do
end select

select case(mod(lg%num(3),2))
  case(1)
    do i3=lg%is(3),lg%ie(3)
       vecR(3,:,:,i3)=dble(i3)-rlaser_center(3)/Hgs(3)
    end do
  case(0)
    do i3=lg%is(3),lg%ie(3)
       vecR(3,:,:,i3)=dble(i3)-0.5d0-rlaser_center(3)/Hgs(3)
    end do
end select

!$OMP parallel do collapse(2) private(iz,iy,ix)
do iz=lg%is(3),lg%ie(3)
do iy=lg%is(2),lg%ie(2)
do ix=lg%is(1),lg%ie(1)
   R1(ix,iy,iz)=(epdir_re1(1)*lg%coordinate(ix,1)+   &
                 epdir_re1(2)*lg%coordinate(iy,2)+   &
                 epdir_re1(3)*lg%coordinate(iz,3))
end do
end do
end do

if(num_dipole_source>=1)then
  allocate(vonf_sd(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  allocate(eonf_sd(3,mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  vonf_sd=0.d0
  eonf_sd=0.d0
  call set_vonf_sd
end if

if(IC_rt==0)then
  rbox_array=0.d0
  do i1=1,3
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       rbox_array(i1)=rbox_array(i1)+vecR(i1,ix,iy,iz)*srho%f(ix,iy,iz)
    end do
    end do
    end do
  end do

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
     rbox_array(4)=rbox_array(4)+srho%f(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(rbox_array,rbox_array2,4,nproc_group_global)
  vecDs(1:3)=rbox_array2(1:3)*Hgs(1:3)*Hvol

end if
if(comm_is_root(nproc_id_global))then
  write(*,'(a30)', advance="no") "Static dipole moment(xyz) ="
  write(*,'(3e15.8)') (vecDs(i1)*a_B, i1=1,3)
  write(*,*)
endif

! Initial wave function
if(iperiodic==0 .and. ikind_eext==0 .and. IC_rt==0)then
  do iik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
  do jspin=1,system%nspin
        spsi_in%zwf(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),jspin,iob,iik,1) &
        = exp(zi*Fst*R1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
           mg_sta(3):mg_end(3)))   &
           *  spsi_in%zwf(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),jspin,iob,iik,1)
  end do
  end do
  end do
end if

rIe(0)=rbox_array2(4)*Hvol
Dp(:,0)=0.d0
Qp(:,:,0)=0.d0

rho0 = srho%f

  do itt=0,0
    if(yn_out_dns_rt=='y')then
      call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,hgs,iscfrt,rho0,itt)
    end if
    if(yn_out_elf_rt=='y')then
      call write_elf(iscfrt,itt,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi_in)
    end if
    if(yn_out_estatic_rt=='y')then
      call calcEstatic(lg, ng, info, sVh, srg_ng)
      call writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,matbox_l,matbox_l2,hgs,itt)
    end if
  end do

if(iperiodic==3) call calcAext

if(yn_md=='y') call init_md(system,md)

! single-scale Maxwell-TDDFT
if(use_singlescale=='y') then
  if(comm_is_root(nproc_id_global)) write(*,*) "single-scale Maxwell-TDDFT method"
  call allocate_vector(mg,j_e)
  call allocate_scalar(mg,system%div_Ac)
  call allocate_vector(mg,system%Ac_micro)
end if

!-------------------------------------------------- Time evolution

!(force at initial step)
if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
   call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi_in)

   !open trj file for coordinate, velocity, and force (rvf) in xyz format
   write(comment_line,10) -1, 0.0d0
   if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
        &  write(comment_line,12) trim(comment_line), md%xi_nh
   call write_xyz(comment_line,"add","rvf",system)
10 format("#rt   step=",i8,"   time",e16.6)
12 format(a,"  xi_nh=",e18.10)
end if

call timer_end(LOG_INIT_TIME_PROPAGATION)


call timer_begin(LOG_INIT_RT)
allocate(zc(N_hamil))
call taylor_coe(N_hamil,dt,zc)
call timer_end(LOG_INIT_RT)


call timer_begin(LOG_RT_ITERATION)
TE : do itt=Miter_rt+1,itotNtime

  if(mod(itt,2)==1)then
    call time_evolution_step(lg,mg,ng,system,info,info_field,stencil,xc_func &
     & ,srg,srg_ng,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,sVh,sVh_stock1,sVh_stock2,sVxc,sVpsl,dmat,fg,energy,md,ofl &
     & ,poisson,j_e,singlescale)
  else
    call time_evolution_step(lg,mg,ng,system,info,info_field,stencil,xc_func &
     & ,srg,srg_ng,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,sVh,sVh_stock1,sVh_stock2,sVxc,sVpsl,dmat,fg,energy,md,ofl &
     & ,poisson,j_e,singlescale)
  end if

  if(checkpoint_interval .ge. 1) then
     if(mod(itt,checkpoint_interval)==0) call write_checkpoint_rt(itt,ng,info)
  endif

end do TE
call timer_end(LOG_RT_ITERATION)

close(030) ! laser

if(ikind_eext.ne.0) deallocate (Vbox)
deallocate (R1)


!--------------------------------- end of time-evolution


!------------ Writing part -----------

! write RT: 
call timer_begin(LOG_WRITE_RT_RESULTS)
if(iwrite_external==1)then
  if(comm_is_root(nproc_id_global))then
    open(1,file=file_external)
    if(ikind_eext==1)then
      do nntime=0,itotNtime
        write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
        if(dt*itt <= tau)then
          write(1,'(e16.8)',advance="yes") Fst*sin(romega*dble(nntime)*dt)*sin(Pi*dble(nntime)*dt/pulse_T)**2
        else
          write(1,'(e16.8)',advance="yes") 0.d0
        end if
      end do
    end if
    close(1)
  end if
end if

!
select case(iperiodic)
case(0)

  call Fourier3D(Dp,alpha_R,alpha_I)
  if(comm_is_root(nproc_id_global))then
    open(1,file=file_RT)
    write(1,'(a)') "# time[fs],    dipoleMoment(x,y,z)[A],                        Energy[eV]" 
     do nntime=1,itotNtime
        write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
        write(1,'(3e16.8)',advance="no") (Dp(iii,nntime)*a_B, iii=1,3)
        write(1,'(e16.8)',advance="yes") tene(nntime)*2.d0*Ry
     end do
    close(1)
  
    if(iflag_intelectron==1)then
      open(1,file=file_RT_e)
      write(1,'(a)') "# time[fs],    integrated electron density" 
       do nntime=1,itotNtime
          write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
          write(1,'(e16.8)',advance="yes") rIe(nntime)
       end do
      close(1)
    end if

    ! Alpha
    if(ae_shape1=='impulse')then
      open(1,file=file_alpha_lr)
      write(1,'(a)') "# energy[eV], Re[alpha](x,y,z)[A**3], Im[alpha](x,y,z)[A**3], df/dE(x,y,z)[1/eV]" 
      do iene=0,Nenergy
        Sf(:)=2*iene*dE/(Pi)*alpha_I(:,iene)
        write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
        write(1,'(3e16.8)',advance="no") (alpha_R(iii,iene)*(a_B)**3, iii=1,3)
        write(1,'(3e16.8)',advance="no") (alpha_I(iii,iene)*(a_B)**3, iii=1,3)
        write(1,'(3e16.8)',advance="yes") (Sf(iii)/2d0/Ry, iii=1,3)
      end do
    else
      open(1,file=file_alpha_pulse)
      write(1,'(a)') "# energy[eV], Re[d(w)](x,y,z)[A*fs],  Im[d(w)](x,y,z)[A*fs],  |d(w)|^2(x,y,z)[A**2*fs**2]"
      do iene=0,Nenergy
        write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
        write(1,'(3e16.8)',advance="no") (alpha_R(iii,iene)*(a_B)*(2.d0*Ry*fs2eVinv), iii=1,3)
        write(1,'(3e16.8)',advance="no") (alpha_I(iii,iene)*(a_B)*(2.d0*Ry*fs2eVinv), iii=1,3)
        write(1,'(3e16.8)',advance="yes") ((alpha_R(iii,iene)**2+alpha_I(iii,iene)**2)   &
                                               *(a_B)**2*(2.d0*Ry*fs2eVinv)**2, iii=1,3)
      end do
    end if 
    close(1)

  end if
  
case(3)
  allocate( tfourier_integrand(1:3,0:Ntime) )
  if(trans_longi=="lo")then
    tfourier_integrand(1:3,0:Ntime)=A_ind(1:3,0:Ntime)
  else if(trans_longi=="tr")then
    tfourier_integrand(1:3,0:Ntime)=curr(1:3,0:Ntime)
  end if
  call Fourier3D(tfourier_integrand,alpha_R,alpha_I)
  if(comm_is_root(nproc_id_global))then
    open(1,file=file_alpha_lr)
    write(1,*) "# energy[eV], Re[epsilon](x,y,z), Im[epsilon](x,y,z)" 
    do nn=1,Nenergy
      write(1,'(e13.5)',advance="no") nn*dE*2d0*Ry
!      write(1,'(3e16.8)',advance="no")      &
!           (F*(F+alpha_R(iii,n))/((F+alpha_R(iii,n))**2+alpha_I(iii,n)**2), iii=1,3)
!      write(1,'(3e16.8)',advance="yes")     &
!           (-F*alpha_I(iii,n)/((F+alpha_R(iii,n))**2+alpha_I(iii,n)**2), iii=1,3)
      write(1,'(3e16.8)',advance="no") (alpha_R(iii,nn), iii=1,3)
      write(1,'(3e16.8)',advance="yes") (alpha_I(iii,nn), iii=1,3)
    end do
    close(1)
  end if 
  deallocate( tfourier_integrand )

end select
call timer_end(LOG_WRITE_RT_RESULTS)
call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
  call write_bin(ofile%dir_out_restart,lg,mg,ng,system,info,spsi_out,mixing,sVh_stock1,sVh_stock2,miter)
end if

call deallocate_mat

call finalize_xc(xc_func)

contains

  subroutine init_code_optimization
    implicit none
    integer :: ignum(3)

    call switch_stencil_optimization(mg%num)
    call switch_openmp_parallelization(mg%num)

    if(iperiodic==3 .and. nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)==1) then
       ignum = mg%num
    else
       ignum = mg%num + (nd*2)
    end if
    call set_modulo_tables(ignum)

    if (comm_is_root(nproc_id_global)) then
       call optimization_log(nproc_k, nproc_ob, nproc_d_o, nproc_d_g)
    end if
  end subroutine init_code_optimization

end subroutine main_tddft

!=======================================================================
! Fourier transform for 3D

SUBROUTINE Fourier3D(Dp_t,alpha_R,alpha_I)
use global_variables_rt
implicit none

real(8),intent(IN) :: Dp_t(3,0:Ntime)
real(8),intent(OUT) :: alpha_R(3,0:Nenergy),alpha_I(3,0:Nenergy)
complex(8),parameter   :: zi=(0.d0,1.d0)
complex(8),allocatable :: zalpha(:)
integer :: iene,nntime
real(8) :: t2,hw,TT
allocate(zalpha(3))

! Fourier Transform

TT = dt*itotNtime ! [a.u.]

do iene=0,Nenergy
  hw=iene*dE ; zalpha=(0.d0,0.d0)  ! [a.u.]
  do nntime=1,itotNtime
     t2=nntime*dt ; zalpha(:)=zalpha(:)+exp(zi*hw*t2)*Dp_t(:,nntime) & !hw*t is unitless      
                       *(1-3*(t2/TT)**2+2*(t2/TT)**3)
  end do
  select case(iperiodic)
  case(0)
    if(ikind_eext==0.or.ikind_eext==10)then
      zalpha=zalpha/Fst*dt
    else
      zalpha=zalpha*dt 
    end if
  case(3)
    if(ikind_eext==0.or.ikind_eext==10)then
      zalpha=zalpha/Fst*dt
      if(trans_longi=="tr")then
        if (iene == 0) then ! WARNING: zero divide happens when iene is 0
          zalpha(1:3)=0
        else
          zalpha(1:3)=1.d0+4.d0*Pi*zi*zalpha(1:3)/hw
        end if
      else if(trans_longi=="lo")then
        zalpha(1:3)=1.d0/(1.d0-zi*hw*zalpha(1:3))
      end if
    else
      zalpha=zalpha*dt
    end if
  end select
  alpha_R(:,iene)=real(zalpha(:),8)    ! Real part
  alpha_I(:,iene)=aimag(zalpha(:))      ! Imaginary part
end do

deallocate(zalpha)
END SUBROUTINE Fourier3D

