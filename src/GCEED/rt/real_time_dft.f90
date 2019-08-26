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
use calc_allob_sub
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
use init_sendrecv_sub
use new_world_sub
use read_pslfile_sub
use allocate_psl_sub
use persistent_comm

implicit none

integer       :: Ntime

real(8)       :: debye2au   ! [D]  -> [a.u.] 
integer       :: iii

real(8), allocatable :: alpha2(:,:,:,:)

END MODULE global_variables_rt

!=======================================================================

subroutine Real_Time_DFT
use structures
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc, only: init_xc, finalize_xc
use timer
use global_variables_rt
use print_sub, only: write_xyz,write_rt_data_3d,write_rt_energy_data
use code_optimization
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_dft_system)  :: system
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_poisson_cg) :: poisson_cg
type(s_stencil) :: stencil
type(s_reciprocal_grid) :: fg
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_mixing) :: mixing
real(8),allocatable :: alpha_R(:,:),alpha_I(:,:) 
real(8),allocatable :: alphaq_R(:,:,:),alphaq_I(:,:,:)
real(8),allocatable :: Sf(:)
integer :: jj,nn
integer :: iene,nntime,ix,iy,iz
character(100):: comment_line
integer :: ia,ib
real(8) :: rab
real(8),allocatable :: tfourier_integrand(:,:)

call timer_begin(LOG_TOTAL)

call timer_begin(LOG_INIT_RT)
call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

call check_cep
call check_ae_shape

iSCFRT=2
OC=0
img=1

iwdenstep=30 
denplane='xy'
idensum=0
posplane=0.d0

inumcpu_check=0

call convert_input_rt(Ntime,mixing,poisson_cg)
allocate(system%mass(1:nelem))

call set_filename

if(comm_is_root(nproc_id_global))then
  write(*,*)
  write(*,*) "Total time step      =",Ntime
  write(*,*) "Time step[fs]        =",dt*au_time_fs
  write(*,*) "Field strength[?]    =",Fst
  write(*,*) "Energy range         =",Nenergy
  write(*,*) "Energy resolution[eV]=",dE*au_energy_ev
  write(*,*) "ikind_eext is           ", ikind_eext
  write(*,*) "Step for writing dens=", iwdenstep
  write(*,*) "Plane showing density=", denplane
  write(*,*) "idensum              =", idensum 
  if(idensum==0) write(*,*) "Position of the plane=", posplane
  select case (ikind_eext)
    case(1,6,7,8,15)
      write(*,'(a21,f5.2,a4)') "Laser frequency     =",       &
                           romega*au_energy_ev, "[eV]"
      write(*,'(a21,f16.8,a4)') "Pulse width of laser=",      &
                           pulse_T*au_time_fs,"[fs]"
      write(*,'(a21,e16.8,a8)') "Laser intensity      =",      &
                           rlaser_I, "[W/cm^2]"
      write(*,'(a21,e16.8,a8)') "tau                  =",      &
                           tau*au_time_fs, "[fs]"
    case(4,12)
      write(*,'(a21,2f5.2,a4)') "Laser frequency     =",       &
                          romega2(1)*au_energy_ev &
                          ,romega2(2)*au_energy_ev, "[eV]"
      write(*,'(a21,2f16.8,a4)') "Pulse width of laser=",      &
                          pulse_T2(1)*au_time_fs&
                          ,pulse_T2(2)*au_time_fs,"[fs]"
      write(*,'(a21,2e16.8,a8)') "Laser intensity      =",      &
                          rlaser_I2(1),rlaser_I2(2), "[W/cm^2]"
      write(*,'(a21,f16.8,a4)') "delay time           =",      &
                          delay*au_time_fs, "[fs]"
      write(*,'(a21,f16.8)') "rcycle                =",rcycle
  end select

end if

debye2au = 0.393428d0

select case (ikind_eext)
  case(0,10)
    Fst=Fst !/5.14223d1
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


call timer_begin(LOG_READ_LDA_DATA)
! Read SCF data
call IN_data(lg,mg,ng,info,info_field,system,stencil,mixing)

if(comm_is_root(nproc_id_global))then
  if(iflag_md==1)then
    do jj=1,2
      if(idisnum(jj)>MI) then
        write(*,*) "idisnum is larger than MI"
        stop
      end if
    end do
  end if
end if

if(iperiodic==3 .and. iflag_hartree==4)then
  call prep_poisson_fft(ng)
end if

call read_pslfile(system)
call allocate_psl
call init_ps(ng,system%primitive_a,system%primitive_b,stencil%rmatrix_A,info%icomm_r)

call init_updown(info)
call init_itype
call init_sendrecv_matrix

call allocate_sendrecv
call init_persistent_requests(info)
call init_code_optimization

if(ilsda==0)then
  numspin=1
else if(ilsda==1)then
  numspin=2
end if

if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(ng,poisson_cg)
call set_ig_bound(ng,poisson_cg)

call timer_end(LOG_READ_LDA_DATA)


call timer_begin(LOG_INIT_RT)

Eion=0.d0
do ia=1,MI
do ib=1,ia-1
  rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
           +(Rion(2,ia)-Rion(2,ib))**2      &
           +(Rion(3,ia)-Rion(3,ib))**2)
  Eion=Eion+Zps(Kion(ia))*Zps(Kion(ib))/rab
end do
end do
call timer_end(LOG_INIT_RT)


call timer_begin(LOG_READ_RT_DATA)
allocate(Ex_fast(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(Ec_fast(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho0(ix,iy,iz) = rho(ix,iy,iz)
end do
end do
end do

allocate( Vh0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

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
else if(IC_rt==1) then
  call IN_data_rt(ng,info,Ntime)
end if
call timer_end(LOG_READ_RT_DATA)


call timer_begin(LOG_INIT_RT)
allocate( alpha_R(3,0:Nenergy), & 
                    alpha_I(3,0:Nenergy), Sf(3) )
allocate( alphaq_R(3,3,0:Nenergy), & 
                    alphaq_I(3,3,0:Nenergy) )

ntmg=1
! 'Hartree' parameter

Hconv  = Hconv !/(2d0*Ry)**2d0/a_B**3   ! Convergence criterion
iterVh = 0        ! Iteration counter


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
   write(ofl%file_rt_data,"(2A,'_rt.data')") trim(directory),trim(SYSname)
   write(ofl%file_rt_energy_data,"(2A,'_rt_energy.data')") trim(directory),trim(SYSname)
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
  case(0)
  case(3)
     call write_rt_data_3d(-1,ofl,iflag_md,dt)
  end select

  !(header of SYSname_rt_energy.data)
  call write_rt_energy_data(-1,ofl,iflag_md,dt,energy,md)

  !(header in SYSname_proj.data)
  if(iwrite_projection==1)then
    open(41,file=file_Projection)
    write(41,'("#",5X,"time[fs]",4("    projection"))')
    write(41,'("#",13x,4("  orbital",i5))') (iwrite_projection_ob(jj),jj=1,4)
    write(41,'("#",13x,4("        k",i5))') (iwrite_projection_k(jj), jj=1,4)
    write(41,'("#",7("----------"))')
  end if
end if

if(iflag_md==1 .or. icalcforce==1)then
   call write_xyz(comment_line,"new","rvf",system)
endif


! Go into Time-Evolution
call Time_Evolution(lg,mg,ng,system,info,info_field,stencil,fg,energy,md,ofl,poisson_cg)


call timer_begin(LOG_WRITE_RT_DATA)
if(OC_rt==1) call OUT_data_rt(ng)
call timer_end(LOG_WRITE_RT_DATA)


! Output after time-evolution
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
  if(iflag_indA==1)then
    tfourier_integrand(1:3,0:Ntime)=A_ind(1:3,0:Ntime)
  else if(iflag_indA==0)then
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
end subroutine

END subroutine Real_Time_DFT

!=========%==============================================================

SUBROUTINE Time_Evolution(lg,mg,ng,system,info,info_field,stencil,fg,energy,md,ofl,poisson_cg)
use structures
use salmon_parallel, only: nproc_group_global, nproc_id_global, & 
                           nproc_group_h, nproc_size_global
use salmon_communication, only: comm_is_root, comm_summation
use density_matrix, only: calc_density
use writefield
use timer
use global_variables_rt
use init_sendrecv_sub, only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
use sendrecv_grid, only: init_sendrecv_grid
use salmon_pp, only: calc_nlcc
use calc_iroot_sub
use force_sub, only: calc_force_salmon
use hpsi_sub, only: update_kvector_nonlocalpt
use md_sub, only: init_md
use print_sub, only: write_xyz
use fdtd_coulomb_gauge, only: ls_singlescale
implicit none

type(s_rgrid) :: lg,mg,ng
type(s_dft_system) :: system
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_stencil) :: stencil
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_ob_1,srg_ng
type(s_pp_nlcc) :: ppn
type(s_reciprocal_grid) :: fg
type(s_md) :: md
type(s_dft_energy) :: energy
type(s_ofile) :: ofl
type(s_poisson_cg) :: poisson_cg
type(s_vector) :: j_e ! microscopic electron number current density
type(ls_singlescale) :: singlescale

complex(8),parameter :: zi=(0.d0,1.d0)
integer :: iob,i1,i2,i3,ix,iy,iz,jj,ik,iik,n,nn
integer :: nspin
real(8),allocatable :: R1(:,:,:)
character(10):: fileLaser
integer:: idensity, idiffDensity, ielf
integer :: is,jspin
integer :: neig(1:3, 1:2)
integer :: neig_ng(1:3, 1:2)

real(8)    :: rbox_array(10)
real(8)    :: rbox_array2(10)

character(100) :: comment_line

type(s_scalar) :: srho,sVh,sVpsl
type(s_scalar),allocatable :: srho_s(:),V_local(:),sVxc(:)
type(s_dmatrix) :: dmat

call timer_begin(LOG_INIT_TIME_PROPAGATION)

  if(ispin==0)then
    nspin=1
  else
    nspin=2
  end if

  system%rocc(1:itotMST,1:system%nk,1) = rocc(1:itotMST,1:system%nk)

  allocate(energy%esp(system%no,system%nk,system%nspin))

  info%im_s=1
  info%im_e=1
  info%numm=1
  info%ik_s=k_sta
  info%ik_e=k_end
  info%numk=k_num
  info%io_s=1
  info%io_e=iobnum/nspin
  info%numo=iobnum/nspin

  info%if_divide_rspace = nproc_d_o_mul.ne.1   ! moved just after init_lattice
  info%if_divide_orbit = nproc_ob.ne.1
  info%icomm_rko = nproc_group_global

  allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin,1) &
          ,info%io_tbl(info%io_s:info%io_e), info%jo_tbl(1:system%no) &
          ,info%irank_jo(1:system%no))
  info%jo_tbl(:) = 0 ! info%io_s-1 (initial value)
  do iob=info%io_s,info%io_e
    call calc_allob(iob,jj,itotmst,mst,iobnum)
    info%io_tbl(iob) = jj
    info%jo_tbl(jj) = iob
  end do
  do jj=1, system%no
    call calc_iroot(jj,info%irank_jo(jj),ilsda,nproc_ob,itotmst,mst)
  end do

  do ik=info%ik_s,info%ik_e
    do iob=info%io_s,info%io_e
      do jspin=1,system%nspin
        jj = info%io_tbl(iob)
        info%occ(iob,ik,jspin,1) = system%rocc(jj,ik,1)*system%wtk(ik)
      end do
    end do
  end do

  ! Initialization of s_sendrecv_grid structure (experimental implementation)
  neig(1, 1) = iup_array(1)
  neig(1, 2) = idw_array(1)
  neig(2, 1) = jup_array(1)
  neig(2, 2) = jdw_array(1)
  neig(3, 1) = kup_array(1)
  neig(3, 2) = kdw_array(1)
  call init_sendrecv_grid(srg, mg, iobnum * k_num, info%icomm_r, neig)
  call init_sendrecv_grid(srg_ob_1, mg, 1, info%icomm_r, neig)

  neig_ng(1, 1) = iup_array(2)
  neig_ng(1, 2) = idw_array(2)
  neig_ng(2, 1) = jup_array(2)
  neig_ng(2, 2) = jdw_array(2)
  neig_ng(3, 1) = kup_array(2)
  neig_ng(3, 2) = kdw_array(2)
  call init_sendrecv_grid(srg_ng, ng, 1, &
    & nproc_group_global, neig_ng)

  call allocate_orbital_complex(system%nspin,mg,info,spsi_in)
  call allocate_orbital_complex(system%nspin,mg,info,spsi_out)
  call allocate_orbital_complex(system%nspin,mg,info,tpsi)
  call allocate_dmatrix(system%nspin,mg,info,dmat)

  if(iperiodic==3) then
    allocate(stencil%vec_kAc(3,info%ik_s:info%ik_e))

!????????? get_fourier_grid_G @ real_space_dft.f90
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
  end if

if(comm_is_root(nproc_id_global).and.iperiodic==3) then
  open(16,file="current.data")
  open(17,file="Etot.data")
  open(18,file="Eext.data")
  open(19,file="Eind.data")
end if

cumnum=0.d0

idensity=0
idiffDensity=1
ielf=2
fileLaser= "laser.out"

allocate (R1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3))) 
!if(ikind_eext.ne.0)then
  allocate( Vbox(lg_sta(1)-Nd:lg_end(1)+Nd,lg_sta(2)-Nd:lg_end(2)+Nd, & 
                                           lg_sta(3)-Nd:lg_end(3)+Nd))
!endif

allocate( elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3))) 

allocate(rhobox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
!if(ilsda==1)then
  allocate(rhobox_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2))
!end if

  call allocate_scalar(mg,srho)
  call allocate_scalar(mg,sVh)
  call allocate_scalar(mg,sVpsl)
  allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))
  do jspin=1,system%nspin
    call allocate_scalar(mg,srho_s(jspin))
    call allocate_scalar(mg,V_local(jspin))
    call allocate_scalar(mg,sVxc(jspin))
  end do
  sVpsl%f = Vpsl

  call calc_nlcc(pp, system, mg, ppn)
  if (comm_is_root(nproc_id_global)) then
    write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
    write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
  end if

  call calc_density(srho_s,spsi_in,info,mg,nspin)

  if(ilsda==0)then  
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(ix,iy,iz,1)=srho_s(1)%f(ix,iy,iz)
      rho_s(ix,iy,iz,2)=srho_s(2)%f(ix,iy,iz)
      rho(ix,iy,iz)=srho_s(1)%f(ix,iy,iz)+srho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  
!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho0(ix,iy,iz)=rho(ix,iy,iz)
end do
end do
end do

allocate(zc(N_hamil))

! External Field Direction
select case (ikind_eext)
  case(1)
    if(yn_local_field=='y')then
      rlaser_center(1:3)=(rlaserbound_sta(1:3)+rlaserbound_end(1:3))/2.d0
      do jj=1,3
        select case(imesh_oddeven(jj))
          case(1)
            ilasbound_sta(jj)=nint(rlaserbound_sta(jj)/Hgs(jj))
            ilasbound_end(jj)=nint(rlaserbound_end(jj)/Hgs(jj))
          case(2)
            ilasbound_sta(jj)=nint(rlaserbound_sta(jj)/Hgs(jj)+0.5d0)
            ilasbound_end(jj)=nint(rlaserbound_end(jj)/Hgs(jj)+0.5d0)
        end select
      end do
    else
      rlaser_center(1:3)=0.d0
    end if
end select 

select case(imesh_oddeven(1))
  case(1)
    do i1=lg_sta(1),lg_end(1)
      vecR(1,i1,:,:)=dble(i1)-rlaser_center(1)/Hgs(1)
    end do
  case(2)
    do i1=lg_sta(1),lg_end(1)
      vecR(1,i1,:,:)=dble(i1)-0.5d0-rlaser_center(1)/Hgs(1)
    end do
end select

select case(imesh_oddeven(2))
  case(1)
    do i2=lg_sta(2),lg_end(2)
      vecR(2,:,i2,:)=dble(i2)-rlaser_center(2)/Hgs(2)
    end do
  case(2)
    do i2=lg_sta(2),lg_end(2)
      vecR(2,:,i2,:)=dble(i2)-0.5d0-rlaser_center(2)/Hgs(2)
    end do
end select

select case(imesh_oddeven(3))
  case(1)
    do i3=lg_sta(3),lg_end(3)
      vecR(3,:,:,i3)=dble(i3)-rlaser_center(3)/Hgs(3)
    end do
  case(2)
    do i3=lg_sta(3),lg_end(3)
      vecR(3,:,:,i3)=dble(i3)-0.5d0-rlaser_center(3)/Hgs(3)
    end do
end select

!$OMP parallel do collapse(2) private(iz,iy,ix)
do iz=lg_sta(3),lg_end(3)
do iy=lg_sta(2),lg_end(2)
do ix=lg_sta(1),lg_end(1)
   R1(ix,iy,iz)=(epdir_re1(1)*gridcoo(ix,1)+   &
                 epdir_re1(2)*gridcoo(iy,2)+   &
                 epdir_re1(3)*gridcoo(iz,3))
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
      rbox_array(i1)=rbox_array(i1)+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)
    end do
    end do
    end do
  end do
  
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rbox_array(4)=rbox_array(4)+rho(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(rbox_array,rbox_array2,4,nproc_group_h)
  vecDs(1:3)=rbox_array2(1:3)*Hgs(1:3)*Hvol

end if
if(comm_is_root(nproc_id_global))then
  write(*,'(a30)', advance="no") "Static dipole moment(xyz) ="
  write(*,'(3e15.8)') (vecDs(i1)*a_B, i1=1,3)
  write(*,*)
endif

! Initial wave function
if(iperiodic==0)then
  if(IC_rt==0)then
  if(iobnum.ge.1)then
    do iik=k_sta,k_end
    do iob=1,iobnum
      select case (ikind_eext)
        case(0)
          zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3),iob,iik)  &
          = exp(zi*Fst*R1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3)))   &
             *  zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                     mg_sta(3):mg_end(3),iob,iik) 
      end select 
    end do
    end do
  end if
  end if
end if

!$OMP parallel do private(ik,iob,is,iz,iy,ix) collapse(5)
  do ik=info%ik_s,info%ik_e
  do iob=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        spsi_in%zwf(ix,iy,iz,is,iob,ik,1)=zpsi_in(ix,iy,iz,iob+(is-1)*info%numo,ik)
      end do
      end do
      end do
    end do
  end do
  end do

  do is=1,system%nspin
    V_local(is)%f = Vlocal(:,:,:,is)
  end do

rIe(0)=rbox_array2(4)*Hvol
Dp(:,0)=0.d0
Qp(:,:,0)=0.d0

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Vh0(ix,iy,iz)=Vh(ix,iy,iz)
  sVh%f(ix,iy,iz)=Vh(ix,iy,iz)
end do
end do
end do

  do itt=0,0
    if(yn_out_dns_rt=='y')then
      call writedns(lg,mg,ng,rho,matbox_m,matbox_m2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,rho0,itt)
    end if
    if(yn_out_elf_rt=='y')then
      call calcELF(mg,ng,info,srho,itt,srg_ob_1)
      call writeelf(lg,elf,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt,itt)
    end if
    if(yn_out_estatic_rt=='y')then
      call calcEstatic(ng, info, sVh, srg_ng)
      call writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,matbox_l,matbox_l2,icoo1d,hgs,igc_is,igc_ie,gridcoo,itt)
    end if
  end do

if(iflag_comm_rho==2)then
  allocate(rhobox1_all(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))) 
  allocate(rhobox2_all(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))) 
!$OMP parallel do private(iz,iy,ix)
  do iz=lg_sta(3),lg_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=lg_sta(1),lg_end(1)
    rhobox1_all(ix,iy,iz) = 0.d0
  end do
  end do
  end do
end if

if(iperiodic==3) call calcAext

if(iflag_md==1) call init_md(system,md)

! single-scale Maxwell-TDDFT
if(use_singlescale=='y') then
  if(comm_is_root(nproc_id_global)) write(*,*) "single-scale Maxwell-TDDFT method"
  call allocate_vector(mg,j_e)
  call allocate_scalar(mg,system%div_Ac)
  call allocate_vector(mg,system%Ac_micro)
  do ik=info%ik_s,info%ik_e
    stencil%vec_kAc(:,ik) = system%vec_k(1:3,ik)
  end do
  call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
end if

!-------------------------------------------------- Time evolution

!(force at initial step)
if(iflag_md==1 .or. icalcforce==1)then
   if(iperiodic==3)then
      do ik=info%ik_s,info%ik_e
        stencil%vec_kAc(1:3,ik) = system%vec_k(1:3,ik)
      end do
      call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
      call get_fourier_grid_G_rt(system,lg,ng,fg)
   endif
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
call taylor_coe
call timer_end(LOG_INIT_RT)


call timer_begin(LOG_RT_ITERATION)
TE : do itt=Miter_rt+1,itotNtime
  if(mod(itt,2)==1)then
    call time_evolution_step(lg,mg,ng,system,info,info_field,stencil &
     & ,srg,srg_ng,srg_ob_1,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,sVh,sVxc,sVpsl,dmat,fg,energy,md,ofl &
     & ,poisson_cg,j_e,singlescale)
  else
    call time_evolution_step(lg,mg,ng,system,info,info_field,stencil &
     & ,srg,srg_ng,srg_ob_1,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,sVh,sVxc,sVpsl,dmat,fg,energy,md,ofl &
     & ,poisson_cg,j_e,singlescale)
  end if
end do TE
call timer_end(LOG_RT_ITERATION)

  if(iperiodic==3) deallocate(stencil%vec_kAc)

close(030) ! laser

deallocate (R1)
deallocate (Vlocal)
if(ikind_eext.ne.0)then
  deallocate (Vbox)
endif
END SUBROUTINE Time_Evolution

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
      if(iflag_indA==0)then
        if (iene == 0) then ! WARNING: zero divide happens when iene is 0
          zalpha(1:3)=0
        else
          zalpha(1:3)=1.d0+4.d0*Pi*zi*zalpha(1:3)/hw
        end if
      else if(iflag_indA==1)then
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

