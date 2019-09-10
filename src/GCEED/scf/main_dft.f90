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

MODULE global_variables_scf

use inputoutput
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
use new_world_sub
use change_order_sub
use read_pslfile_sub
use structure_opt_sub
use salmon_total_energy
use hpsi_sub
implicit none

END MODULE global_variables_scf

!=======================================================================

subroutine main_dft
use structures
use salmon_parallel, only: nproc_id_global,nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc, only: init_xc, finalize_xc
use timer
use scf_iteration_sub
use rmmdiis_sub
use density_matrix, only: calc_density
use writefield
use global_variables_scf
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use gram_schmidt_orth, only: gram_schmidt 
use write_sub
use read_gs
use code_optimization
use salmon_initialization
use occupation
use init_poisson_sub
implicit none
integer :: ix,iy,iz,ik,i,j
integer :: iter,iatom,iob,p1,p2,p5,jj,iflag,jspin
real(8) :: sum0,sum1
character(100) :: file_atoms_coo, comment_line
complex(8),allocatable :: zpsi_tmp(:,:,:,:,:)
real(8) :: rNebox1,rNebox2
integer :: itmg,nspin

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_sendrecv_grid) :: srg, srg_ng
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_scalar) :: srho,sVh,sVpsl,rho_old,Vlocal_old
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_cg)  :: cg
type(s_mixing) :: mixing

logical :: rion_update

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
ihpsieff=0
iblacsinit=0

call timer_begin(LOG_TOTAL)


call timer_begin(LOG_INIT_GS)
inumcpu_check=0

call setbN(bnmat)
call setcN(cnmat)

call convert_input_scf(info,info_field,file_atoms_coo,mixing,poisson)

call init_dft(nproc_group_global,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng)

if(stencil%if_orthogonal) then
  if(comm_is_root(nproc_id_global)) write(*,*) "orthogonal cell: using al"
else
  if(comm_is_root(nproc_id_global)) write(*,*) "non-orthogonal cell: using al_vec[1,2,3]"
end if
allocate(system%mass(1:nelem))

call set_filename

k_sta = info%ik_s ! future work: remove this line
k_end = info%ik_e ! future work: remove this line
k_num = info%numk ! future work: remove this line
iobnum = info%numo ! future work: remove this line

if(iflag_opt==1)then
   call structure_opt_ini(MI)
   flag_opt_conv=.false.
   write(comment_line,10) 0
   call write_xyz(comment_line,"new","r  ",system)
10 format("#opt iteration step=",i5)
end if
call timer_end(LOG_INIT_GS)


Structure_Optimization_Iteration : do iopt=1,iter_opt
Multigrid_Iteration : do img=1,ntmg

if(iopt==1)then
  select case( IC )
  case default ! New calculation

    call timer_begin(LOG_INIT_GS)

    Hvol = system%Hvol
    Hgs = system%Hgs
    Miter = 0        ! Miter: Iteration counter set to zero
    itmg=img
    call set_imesh_oddeven(itmg)
    call old_mesh(lg,mg,ng) ! future work: remove this line

  case(1,3) ! Continue the previous calculation

    call read_gs_bin(lg,mg,ng,info,info_field,system,stencil,mixing)

  end select

  select case(iperiodic)
  case(0)
    if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(lg,ng,poisson)
  end select
  call set_ig_bound(lg,ng,poisson)

  call allocate_mat(ng)
  call set_icoo1d(lg)
  call init_code_optimization

  if(ispin==0)then
    nspin=1
  else
    nspin=2
  end if

  allocate(energy%esp(system%no,system%nk,system%nspin))
  allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))

  call allocate_scalar(mg,srho)
  call allocate_scalar(mg,sVh)
  do jspin=1,system%nspin
    call allocate_scalar(mg,srho_s(jspin))
    call allocate_scalar(mg,V_local(jspin))
    call allocate_scalar(mg,sVxc(jspin))
  end do

  select case(iperiodic)
  case(0)
    call allocate_orbital_real(system%nspin,mg,info,spsi)
    call allocate_orbital_real(system%nspin,mg,info,shpsi)
  case(3)
    call allocate_orbital_complex(system%nspin,mg,info,spsi)
    call allocate_orbital_complex(system%nspin,mg,info,shpsi)
    call allocate_orbital_complex(system%nspin,mg,info,sttpsi)
  end select

  if(iperiodic==3)then
    allocate (zpsi_tmp(mg%is_overlap(1):mg%ie_overlap(1) &
    &                 ,mg%is_overlap(2):mg%ie_overlap(2) &
    &                 ,mg%is_overlap(3):mg%ie_overlap(3) &
    &                 ,1:iobnum,k_sta:k_end))
  end if

  if(.not. allocated(Vpsl)) allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
  if(.not. allocated(Vpsl_atom)) allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )
  if(iflag_ps.eq.0)then
    Vpsl=0d0
  else
    call read_pslfile(system)
    call init_ps(lg,mg,ng,system,fg,info_field,poisson,info%icomm_r,sVpsl)
  end if

  if(iperiodic==3) then
    allocate(stencil%vec_kAc(3,info%ik_s:info%ik_e))
    stencil%vec_kAc(:,info%ik_s:info%ik_e) = system%vec_k(:,info%ik_s:info%ik_e)
    call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
  end if

  select case( IC )
  case default ! New calculation

    if(iobnum >= 1)then
      select case(iperiodic)
      case(0)
        allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      case(3)
        allocate( ttpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
        allocate( zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      end select
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      select case(iperiodic)
      case(0)
        allocate( psi_mesh(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:itotMST,1) )
      case(3)
        allocate( zpsi_mesh(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:itotMST,num_kpoints_rd) )
      end select
    end if

    if(read_gs_wfn_k=='n') then
      call init_wf_ns(lg,info,1)
      select case(iperiodic) ! Store to psi/zpsi
      case(0) ; call copy_psi_to_rwf(info,nspin,mg,spsi)
      case(3) ; call copy_zpsi_to_zwf(info,nspin,mg,spsi)
      end select
    else
      if(iperiodic==0) stop "error: read_gs_wfn_k='y' & iperiodic=0"
      call read_wfn(lg,mg,spsi,info,system)
    end if

    call gram_schmidt(system, mg, info, spsi)

    allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    if(ilsda == 1)then
      allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
    end if

    allocate(mixing%srho_in(1:mixing%num_rho_stock+1))
    allocate(mixing%srho_out(1:mixing%num_rho_stock+1))
    do i=1,mixing%num_rho_stock+1
      allocate(mixing%srho_in(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
      allocate(mixing%srho_out(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
      mixing%srho_in(i)%f(:,:,:) =0.d0
      mixing%srho_out(i)%f(:,:,:)=0.d0
    end do

    if(ilsda==1)then
      allocate(mixing%srho_s_in(1:mixing%num_rho_stock+1,2))
      allocate(mixing%srho_s_out(1:mixing%num_rho_stock+1,2))
      do j=1,2
        do i=1,mixing%num_rho_stock+1
          allocate(mixing%srho_s_in(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
          allocate(mixing%srho_s_out(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
          mixing%srho_s_in(i,j)%f(:,:,:) =0.d0
          mixing%srho_s_out(i,j)%f(:,:,:)=0.d0
        end do
      end do
    end if

    if(read_gs_dns_cube == 'n') then
       call calc_density(system,srho_s,spsi,info,mg)
    else
       if(ispin/=0) stop "read_gs_dns_cube=='n' & ispin/=0"
       call read_dns(lg,mg,srho_s(1)%f) ! cube file only
    end if

    srho%f = 0d0
    do jspin=1,nspin
       srho%f = srho%f + srho_s(jspin)%f
    end do
    rho = srho%f

    allocate( Vlocal(mg_sta(1):mg_end(1),  &
                     mg_sta(2):mg_end(2),  &
                     mg_sta(3):mg_end(3),nspin))
    allocate( Vh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    Vh=0.d0

    call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
    Vh = sVh%f

    if(ilsda == 0) then
       allocate( Vxc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    else if(ilsda == 1) then
       allocate( Vxc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
    end if
    allocate( esp(itotMST,num_kpoints_rd) )

    call exc_cor_ns(ng, srg_ng, system%nspin, srho_s, ppn, sVxc, energy%E_xc)

    call allgatherv_vlocal(ng,info,system%nspin,sVh,sVpsl,sVxc,V_local)
    do jspin=1,system%nspin
       Vlocal(:,:,:,jspin) = V_local(jspin)%f
    end do

    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
       call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
    case(3)
       rion_update = .true. ! it's first calculation
       call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    end select
    esp = energy%esp(:,:,1) !++++++++

  case(1,3) ! Continue the previous calculation

    select case(iperiodic) !Store
      case(0) ; call copy_psi_to_rwf(info,nspin,mg,spsi)
      case(3) ; call copy_zpsi_to_zwf(info,nspin,mg,spsi)
    end select
    srho%f = rho
    if(ilsda == 1)then
       srho_s(1)%f = rho_s(:,:,:,1)
       srho_s(2)%f = rho_s(:,:,:,2)
    end if
    do jspin=1,nspin
       V_local(jspin)%f = Vlocal(:,:,:,jspin)
    end do

  end select

  call timer_end(LOG_INIT_GS)

else if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  if(iflag_ps/=0) then
    rion_update = .true.
    call dealloc_init_ps(ppg,ppg_all)
!    call calc_nlcc(pp, system, mg, ppn) !test
    call init_ps(lg,mg,ng,system,fg,info_field,poisson,info%icomm_r,sVpsl)
    if(iperiodic==3) then
       if(.not.allocated(stencil%vec_kAc)) allocate(stencil%vec_kAc(3,info%ik_s:info%ik_e))
       stencil%vec_kAc(:,info%ik_s:info%ik_e) = system%vec_k(:,info%ik_s:info%ik_e)
       call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
    end if

  end if
  call timer_end(LOG_INIT_GS)
end if


if(comm_is_root(nproc_id_global)) then
  write(*,*) '-----------------------------------------------'
  select case(iperiodic)
  case(0) ; write(*,100) Miter,energy%E_tot*2d0*Ry, poisson%iterVh
  case(3) ; write(*,101) Miter,energy%E_tot*2d0*Ry
  end select
100 format(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)
101 format(1x,"iter =",i6,5x,"Total Energy =",f19.8)

  do ik=1,num_kpoints_rd
    if(ik<=3)then
      if(iperiodic==3) write(*,*) "k=",ik
      do p5=1,(itotMST+3)/4
        p1=4*(p5-1)+1
        p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
        write(*,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,ik,1)*2d0*Ry,iob=p1,p2)
      end do
      if(iperiodic==3) write(*,*)
    end if
  end do
end if

!---------------------------------------- Iteration

call timer_begin(LOG_INIT_GS_ITERATION)
iflag=1
poisson%iterVh=1000
sum1=1.0d9

iflag_diisjump=0

allocate(idiis_sd(itotMST))
idiis_sd=0

if(img==1.and.iopt==1) allocate(norm_diff_psi_stock(itotMST,1))
norm_diff_psi_stock=1.0d9

if(img>=2.or.iopt>=2) deallocate(rho_old%f,Vlocal_old%f)
call allocate_scalar(ng,rho_old)
call allocate_scalar(ng,Vlocal_old)

!$OMP parallel do private(iz,iy,ix)
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
   rho_old%f(ix,iy,iz)=srho%f(ix,iy,iz)
   Vlocal_old%f(ix,iy,iz)=V_local(1)%f(ix,iy,iz)
end do
end do
end do

! Setup NLCC term from pseudopotential
call calc_nlcc(pp, system, mg, ppn)

if (comm_is_root(nproc_id_global)) then
  write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
  write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
end if    
call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
DFT_Iteration : do iter=1,iDiter(img)

  if(sum1<threshold) cycle DFT_Iteration

  Miter=Miter+1

  ! for calc_total_energy_periodic
  rion_update = check_rion_update() .or. (iter == 1)

  if(temperature>=0.d0 .and. Miter>iditer_notemperature) then
     call ne2mu(energy,system)
  end if
  rocc(1:itotMST,1:system%nk) = system%rocc(1:itotMST,1:system%nk,1) ! future work: remove this line

  call copy_density(system%nspin,ng,srho_s,mixing)

  if(iscf_order==1)then

    call scf_iteration(lg,mg,ng,system,info,info_field,stencil,srg,srg_ng,spsi,shpsi,srho,srho_s,itotmst,mst, &
                       cg,ppg,V_local,  &
                       iflag_diisjump,energy, &
                       norm_diff_psi_stock, &
                       Miter,iDiterYBCG,   &
                       iflag_subspace_diag,iditer_nosubspace_diag,ifmst,mixing,iter,    &
                       poisson,fg,sVh)

    call timer_begin(LOG_CALC_EXC_COR)
    call exc_cor_ns(ng, srg_ng, system%nspin, srho_s, ppn, sVxc, energy%E_xc)
    call timer_end(LOG_CALC_EXC_COR)

    call allgatherv_vlocal(ng,info,system%nspin,sVh,sVpsl,sVxc,V_local)

    call timer_begin(LOG_CALC_TOTAL_ENERGY)
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0); call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
    case(3); call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    end select
    esp = energy%esp(:,:,1) !++++++++
    call timer_end(LOG_CALC_TOTAL_ENERGY)


    call timer_begin(LOG_CALC_CHANGE_ORDER)
    if(iperiodic==0) call change_order(psi,info)
    call timer_end(LOG_CALC_CHANGE_ORDER)

  end if

  call timer_begin(LOG_WRITE_GS_RESULTS)

  select case(convergence)
    case('rho_dne')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
         sum0 = sum0 + abs(srho%f(ix,iy,iz)-rho_old%f(ix,iy,iz))
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info_field%icomm_all)
      if(ispin==0)then
         sum1 = sum1*Hvol/(dble(ifMST(1))*2.d0)
      else if(ispin==1)then
         sum1 = sum1*Hvol/dble(ifMST(1)+ifMST(2))
      end if
    case('norm_rho','norm_rho_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
         sum0 = sum0 + (srho%f(ix,iy,iz)-rho_old%f(ix,iy,iz))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info_field%icomm_all)
      if(convergence=='norm_rho_dng')then
         sum1 = sum1/dble(lg%num(1)*lg%num(2)*lg%num(3))
      end if
    case('norm_pot','norm_pot_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
         sum0 = sum0 + (V_local(1)%f(ix,iy,iz)-Vlocal_old%f(ix,iy,iz))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info_field%icomm_all)
      if(convergence=='norm_pot_dng')then
         sum1 = sum1/dble(lg%num(1)*lg%num(2)*lg%num(3))
      end if
  end select 

  if(comm_is_root(nproc_id_global)) then
    write(*,*) '-----------------------------------------------'
    select case(iperiodic)
    case(0)
      if(iflag_diisjump == 1) then
         write(*,'("Diisjump occured. Steepest descent was used.")')
      end if
      write(*,100) Miter,energy%E_tot*2d0*Ry, poisson%iterVh
    case(3)
      write(*,101) Miter,energy%E_tot*2d0*Ry
    end select
    do ik=1,num_kpoints_rd
      if(ik<=3)then
        if(iperiodic==3) write(*,*) "k=",ik
        do p5=1,(itotMST+3)/4
          p1=4*(p5-1)+1
          p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
          write(*,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,ik,1)*2d0*Ry,iob=p1,p2)
        end do
        if(iperiodic==3) write(*,*) 
      end if
    end do

    select case(convergence)
      case('rho_dne' )     ; write(*,200) Miter, sum1
      case('norm_rho')     ; write(*,201) Miter, sum1/a_B**6
      case('norm_rho_dng') ; write(*,202) Miter, sum1/a_B**6
      case('norm_pot')     ; write(*,203) Miter, sum1*(2.d0*Ry)**2/a_B**6
      case('norm_pot_dng') ; write(*,204) Miter, sum1*(2.d0*Ry)**2/a_B**6
    end select
200 format("iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec        = ",i6,e15.8)
201 format("iter and ||rho_i(ix)-rho_i-1(ix)||**2              = ",i6,e15.8)
202 format("iter and ||rho_i(ix)-rho_i-1(ix)||**2/(# of grids) = ",i6,e15.8)
203 format("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2              = ",i6,e15.8)
204 format("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2/(# of grids) = ",i6,e15.8)

  end if 
  rNebox1=0.d0 
!$OMP parallel do reduction(+:rNebox1) private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
     rNebox1 = rNebox1 + srho%f(ix,iy,iz)
  end do
  end do
  end do
  call comm_summation(rNebox1,rNebox2,info_field%icomm_all)
  if(comm_is_root(nproc_id_global))then
     write(*,*) "Ne=",rNebox2*Hvol
  end if
  call timer_end(LOG_WRITE_GS_RESULTS)

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
     rho_old%f(ix,iy,iz)    = srho%f(ix,iy,iz)
     Vlocal_old%f(ix,iy,iz) = V_local(1)%f(ix,iy,iz)
  end do
  end do
  end do

end do DFT_Iteration

! for writing GS data
Vpsl = sVpsl%f
if(allocated(Vpsl_atom) .and. allocated(ppg%Vpsl_atom)) &
  Vpsl_atom = ppg%Vpsl_atom
Vh = sVh%f
rho = srho%f
if(ilsda == 1) then
  do jspin=1,system%nspin
     Vxc_s(:,:,:,jspin) = sVxc(jspin)%f
  end do
else
  Vxc = sVxc(1)%f
end if
Exc = energy%E_xc

! Store to psi/zpsi
select case(iperiodic)
  case(0) ; call copy_rwf_to_psi(info,nspin,mg,spsi)
  case(3) ; call copy_zwf_to_zpsi(info,nspin,mg,spsi)
end select

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then
  if(iperiodic==3) then
    call write_wfn(lg,mg,spsi,info,system)
    ! Experimental Implementation of Inner-Product Outputs:
    call write_prod_dk_data(lg, mg, system, info, spsi) 
  else
    write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
  end if
end if

! output transition moment
if(yn_out_tm  == 'y') then
  if(iperiodic==3) then
     call write_k_data(system,stencil)
     call write_tm_data(spsi,system,info,mg,stencil,srg,ppg)
  else
     write(*,*) "error: yn_out_tm='y' & iperiodic=0"
  end if
end if

! force
!if(iflag_opt==1) then
if (iperiodic == 3 .and. iflag_hartree == 4) then
  ! NOTE: calc_force_salmon hangs under this configuration due to ppg%vpsl_atom
  ! does not allocate.
else
   call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi)
   if(comm_is_root(nproc_id_global))then
      write(*,*) "===== force ====="
      do iatom=1,MI
         select case(unit_system)
         case('au','a.u.'); write(*,300)iatom,(system%Force(ix,iatom),ix=1,3)
         case('A_eV_fs'  ); write(*,300)iatom,(system%Force(ix,iatom)*2.d0*Ry/a_B,ix=1,3)
         end select
      end do
300   format(i6,3e16.8)
   end if
end if
!end if

if(iperiodic==3) deallocate(stencil%vec_kAc,ppg%zekr_uV)
deallocate(idiis_sd)
call timer_end(LOG_GS_ITERATION)

call timer_begin(LOG_DEINIT_GS_ITERATION)
if(iflag_opt==1) then
  call structure_opt_check(MI,iopt,flag_opt_conv,system%Force)
  if(.not.flag_opt_conv) call structure_opt(MI,iopt,system)
  !! Rion is old variables to be removed 
  !! but currently it is used in many subroutines.
  Rion(:,:) = system%Rion(:,:) 

  write(comment_line,10) iopt
  call write_xyz(comment_line,"add","r  ",system)

  if(comm_is_root(nproc_id_global))then
    write(*,*) "atomic coordinate"
    do iatom=1,MI
       write(*,20) "'"//trim(AtomName(Kion(iatom)))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
    end do
20  format(a5,3f16.8,i3,a3)
  end if

  if(flag_opt_conv) then
     call structure_opt_fin
     exit Multigrid_Iteration
  end if

else
   select case(iperiodic)
   case(0) ; deallocate(spsi%rwf)
   case(3) ; deallocate(spsi%zwf)
   end select
end if
call timer_end(LOG_DEINIT_GS_ITERATION)


end do Multigrid_Iteration
if(flag_opt_conv)then
  exit Structure_Optimization_Iteration
end if
end do Structure_Optimization_Iteration


!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

call write_band_information
call write_eigen
if(yn_out_psi=='y' ) call write_psi(lg,info)
if(yn_out_dns=='y' ) call write_dns(lg,mg,ng,rho,matbox_m,matbox_m2,icoo1d,hgs,iscfrt)
if(yn_out_dos=='y' ) call calc_dos(info)
if(yn_out_pdos=='y') call calc_pdos(lg,info)
if(yn_out_elf=='y' ) then
  allocate(elf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  call calc_elf(lg,mg,ng,srg,info,srho,0)
  call write_elf(lg,elf,icoo1d,hgs,iscfrt)
  deallocate(elf)
end if
call timer_end(LOG_WRITE_GS_RESULTS)

! write GS data
call timer_begin(LOG_WRITE_GS_DATA)
if( OC==1.or.OC==2.or.OC==3 ) call write_gs_bin(lg,ng,info,mixing)
call timer_end(LOG_WRITE_GS_DATA)

! write GS information
call timer_begin(LOG_WRITE_GS_INFO)
call write_info_data(system,energy)
call timer_end(LOG_WRITE_GS_INFO)

deallocate(Vlocal)
call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)

contains

subroutine copy_psi_to_rwf(info,nspin,mg,spsi)
  use scf_data, only: psi
  implicit none
  integer                 ,intent(in) :: nspin
  type(s_orbital_parallel),intent(in) :: info
  type(s_rgrid)           ,intent(in) :: mg
  type(s_orbital)         ,intent(inout) :: spsi
  integer :: ik,iob,is

  do ik=k_sta,k_end
  do iob=1,info%numo
  do is=1,nspin
     spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob+info%io_s-1,ik,1) = &
         &  psi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
  end do
  end do
  end do
end subroutine copy_psi_to_rwf

subroutine copy_zpsi_to_zwf(info,nspin,mg,spsi)
  use scf_data, only: zpsi
  implicit none
  integer                 ,intent(in) :: nspin
  type(s_orbital_parallel),intent(in) :: info
  type(s_rgrid)           ,intent(in) :: mg
  type(s_orbital)         ,intent(inout) :: spsi
  integer :: ik,iob,is

  do ik=k_sta,k_end
  do iob=1,info%numo
  do is=1,nspin
     spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob+info%io_s-1,ik,1) = &
          & zpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
  end do
  end do
  end do
end subroutine copy_zpsi_to_zwf

subroutine copy_rwf_to_psi(info,nspin,mg,spsi)
  use scf_data, only: psi
  implicit none
  integer                 ,intent(in) :: nspin
  type(s_orbital_parallel),intent(in) :: info
  type(s_rgrid)           ,intent(in) :: mg
  type(s_orbital)         ,intent(in) :: spsi
  integer :: ik,iob,is, ix,iy,iz

  do ik=k_sta,k_end
  do iob=1,info%numo
  do is=1,nspin
     !$OMP parallel do private(iz,iy,ix)
     do iz=mg%is(3),mg%ie(3)
     do iy=mg%is(2),mg%ie(2)
     do ix=mg%is(1),mg%ie(1)
        psi(ix,iy,iz,iob+(is-1)*info%numo,ik) = spsi%rwf(ix,iy,iz,is,iob+info%io_s-1,ik,1)
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine copy_rwf_to_psi

subroutine copy_zwf_to_zpsi(info,nspin,mg,spsi)
  use scf_data, only: zpsi
  implicit none
  integer                 ,intent(in) :: nspin
  type(s_orbital_parallel),intent(in) :: info
  type(s_rgrid)           ,intent(in) :: mg
  type(s_orbital)         ,intent(in) :: spsi
  integer :: ik,iob,is, ix,iy,iz

  do ik=k_sta,k_end
  do iob=1,info%numo
  do is=1,nspin
     !$OMP parallel do private(iz,iy,ix)
     do iz=mg%is(3),mg%ie(3)
     do iy=mg%is(2),mg%ie(2)
     do ix=mg%is(1),mg%ie(1)
        zpsi(ix,iy,iz,iob+(is-1)*info%numo,ik) = spsi%zwf(ix,iy,iz,is,iob+info%io_s-1,ik,1)
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine copy_zwf_to_zpsi

subroutine write_band_information
  implicit none
  integer :: ik
  real(8),dimension(num_kpoints_rd) :: esp_vb_min,esp_vb_max,esp_cb_min,esp_cb_max
  if(comm_is_root(nproc_id_global) .and. itotfMST<itotMST) then
    do ik=1,num_kpoints_rd
      esp_vb_min(ik)=minval(energy%esp(1:itotfMST,ik,1))
      esp_vb_max(ik)=maxval(energy%esp(1:itotfMST,ik,1))
      esp_cb_min(ik)=minval(energy%esp(itotfMST+1:itotMST,ik,1))
      esp_cb_max(ik)=maxval(energy%esp(itotfMST+1:itotMST,ik,1))
    end do
    write(*,*) 'band information-----------------------------------------'
    write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
    write(*,*) 'Top of VB',maxval(esp_vb_max(:))
    write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
    write(*,*) 'Top of CB',maxval(esp_cb_max(:))
    write(*,*) 'Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(*,*) 'Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(*,*) '---------------------------------------------------------'
    write(*,*) 'Bottom of VB[eV]',minval(esp_vb_min(:))*2.0*Ry
    write(*,*) 'Top of VB[eV]',maxval(esp_vb_max(:))*2.0*Ry
    write(*,*) 'Bottom of CB[eV]',minval(esp_cb_min(:))*2.0*Ry
    write(*,*) 'Top of CB[eV]',maxval(esp_cb_max(:))*2.0*Ry
    write(*,*) 'Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*2.0*Ry
    write(*,*) 'BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*2.0*Ry
    write(*,*) '---------------------------------------------------------'
  end if
  return
end subroutine write_band_information

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

end subroutine main_dft


