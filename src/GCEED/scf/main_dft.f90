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
use structure_opt_sub
implicit none

END MODULE global_variables_scf

!=======================================================================

subroutine main_dft
use math_constants, only: pi, zi
use structures
use salmon_parallel, only: nproc_id_global,nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use global_variables_scf
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use input_pp_sub
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use salmon_total_energy
use init_gs, only: init_wf
use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
implicit none
integer :: ix,iy,iz,ik
integer :: iter,iatom,iob,p1,p2,p5,jj,iflag,jspin
real(8) :: sum0,sum1
character(100) :: file_atoms_coo, comment_line
real(8) :: rNebox1,rNebox2
integer :: itmg,nspin

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_process_info) :: pinfo
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_sendrecv_grid) :: srg, srg_ng
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl,rho_old,Vlocal_old
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_ofile)  :: ofile

logical :: rion_update
integer :: iopt,nopt_max
logical :: flag_opt_conv

real(8),allocatable :: esp_old(:,:,:)
real(8) :: tol_esp_diff
integer :: iter_band_kpt, num_band_kpt, nref_band
real(8),allocatable :: band_kpt(:,:)
logical,allocatable :: check_conv_esp(:,:,:)
integer :: iter_band_kpt_end, iter_band_kpt_stride

integer :: i,j, img

if(calc_mode=='DFT_BAND'.and.iperiodic/=3) return

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
ihpsieff=0
iblacsinit=0

call timer_begin(LOG_TOTAL)

call timer_begin(LOG_INIT_GS)

call convert_input_scf(file_atoms_coo)
mixing%num_rho_stock=21

! +----------------+
! | initialization |
! +----------------+

call init_dft(iSCFRT,nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)

call init_code_optimization
call allocate_mat(ng,mg,lg) ! future work: remove this line

allocate( energy%esp(system%no,system%nk,system%nspin) ); energy%esp=0.0d0
allocate( esp_old(system%no,system%nk,system%nspin) ); esp_old=0.0d0

allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))
call allocate_scalar(mg,srho)
call allocate_scalar(mg,sVh)
call allocate_scalar(mg,sVpsl)
do jspin=1,system%nspin
  call allocate_scalar(mg,srho_s(jspin))
  call allocate_scalar(mg,V_local(jspin))
  call allocate_scalar(mg,sVxc(jspin))
end do
allocate(ppg%Vpsl_atom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),natom))
call read_pslfile(system,pp,ppg)
call init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)

select case(iperiodic)
case(0)
  call allocate_orbital_real(system%nspin,mg,info,spsi)
  call allocate_orbital_real(system%nspin,mg,info,shpsi)
case(3)
  call allocate_orbital_complex(system%nspin,mg,info,spsi)
  call allocate_orbital_complex(system%nspin,mg,info,shpsi)
  call allocate_orbital_complex(system%nspin,mg,info,sttpsi)
end select

nspin = system%nspin

if(stencil%if_orthogonal) then
  if(comm_is_root(nproc_id_global)) write(*,*) "orthogonal cell: using al"
else
  if(comm_is_root(nproc_id_global)) write(*,*) "non-orthogonal cell: using al_vec[1,2,3]"
end if

call set_filename

if(yn_opt=='y')then
   call structure_opt_ini(MI)
   flag_opt_conv=.false.
   write(comment_line,10) 0
   call write_xyz(comment_line,"new","r  ",system)
10 format("#opt iteration step=",i5)
end if
call timer_end(LOG_INIT_GS)

if(yn_opt=='y') then ; nopt_max = nopt
else                 ; nopt_max = 1
endif

Structure_Optimization_Iteration : do iopt=1,nopt_max
Multigrid_Iteration : do img=1,ntmg

if(iopt==1)then

  call timer_begin(LOG_INIT_GS)

  call init_mixing(nspin,ng,mixing)

  if (yn_restart == 'y') then
    ! restart from binary
    call restart_gs(lg,mg,ng,system,info,spsi,miter,mixing=mixing)
  else
    ! new calculation
    miter = 0        ! Miter: Iteration counter set to zero
    itmg  = img
    call init_wf(lg,mg,system,info,spsi)
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
  call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
  call exchange_correlation(system,xc_func,ng,srg_ng,srho_s,ppn,info_field%icomm_all,sVxc,energy%E_xc)
  call allgatherv_vlocal(ng,mg,info_field,system%nspin,sVh,sVpsl,sVxc,V_local)

  call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
  select case(iperiodic)
  case(0)
     call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
  case(3)
     rion_update = .true. ! it's first calculation
     call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
  end select

  call timer_end(LOG_INIT_GS)

else if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  rion_update = .true.
  call dealloc_init_ps(ppg)
! call calc_nlcc(pp, system, mg, ppn) !test
  call init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(calc_mode=='DFT_BAND')then
  system%wtk(:) = 0.0d0
  
  call get_band_kpt( band_kpt, nref_band, system )
  
  num_band_kpt = size( band_kpt, 2 )
  !write(*,*) "num_band_kpt=",num_band_kpt
  
  allocate( check_conv_esp(nref_band,system%nk,system%nspin) )
  check_conv_esp=.false.
  
  if ( comm_is_root(nproc_id_global) ) then
  open(100,file='band.dat')
  write(100,*) "Number_of_Bands:",system%no
  write(100,*) "Number_of_kpt_in_each_block:",system%nk
  write(100,*) "Number_of_blocks:",num_band_kpt/system%nk
  end if
end if

if(calc_mode=='DFT_BAND')then
  iter_band_kpt_end = num_band_kpt
  iter_band_kpt_stride = system%nk
else
  iter_band_kpt_end = 1
  iter_band_kpt_stride = 1
end if

Band_Iteration : do iter_band_kpt = 1, iter_band_kpt_end, iter_band_kpt_stride

if(calc_mode=='DFT_BAND')then
  check_conv_esp=.false.
  do ik=1,system%nk
     if ( info%ik_s <= ik .and. ik <= info%ik_e ) then
        system%vec_k(:,ik) = matmul( system%primitive_b, band_kpt(:,iter_band_kpt+ik-1) )
     end if
  end do
  
  if ( comm_is_root(nproc_id_global) ) then
     write(*,'(1x,"iter_band_kpt=",i3," to",i3)') iter_band_kpt, iter_band_kpt+system%nk-1
     write(*,'(1x,3x,2x,a30,2x,a30)') "kpoints","kpoints in Cartesian"
     do ik=iter_band_kpt,iter_band_kpt+system%nk-1
        write(*,'(1x,i3,2x,3f10.5,2x,3f10.5)') ik, band_kpt(:,ik), system%vec_k(:,ik-iter_band_kpt+1)
        write(100,'(1x,i3,2x,3f10.5,2x,3f10.5)') ik, band_kpt(:,ik), system%vec_k(:,ik-iter_band_kpt+1)
     end do
  end if
end if

!---------------------------------------- Iteration

call timer_begin(LOG_INIT_GS_ITERATION)
iflag=1
poisson%iterVh=1000
sum1=1.0d9

iflag_diisjump=0

if(.not.allocated(idiis_sd)) allocate(idiis_sd(itotMST))
idiis_sd=0

if(.not.allocated(norm_diff_psi_stock)) then
  if(img==1.and.iopt==1) allocate(norm_diff_psi_stock(itotMST,1))
end if
norm_diff_psi_stock=1.0d9

if(allocated(rho_old%f)) deallocate(rho_old%f)
if(allocated(Vlocal_old%f)) deallocate(Vlocal_old%f)

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
  if(calc_mode=='DFT_BAND')then
    if(all(check_conv_esp)) cycle DFT_Iteration
  end if

  Miter=Miter+1

  if(calc_mode/='DFT_BAND')then
    ! for calc_total_energy_periodic
    rion_update = check_rion_update() .or. (iter == 1)
  
    if(temperature>=0.d0 .and. Miter>iditer_notemperature) then
       call ne2mu(energy,system)
    end if
  end if

  call copy_density(Miter,system%nspin,ng,srho_s,mixing)

  if(iscf_order==1)then

    call scf_iteration(lg,mg,ng,system,info,info_field,stencil,srg,srg_ng,spsi,shpsi,srho,srho_s,mst, &
                       cg,ppg,V_local,  &
                       Miter,iDiterYBCG,   &
                       iflag_subspace_diag,iditer_nosubspace_diag,ifmst,mixing,iter,    &
                       poisson,fg,sVh,xc_func,ppn,sVxc,energy)

    call allgatherv_vlocal(ng,mg,info_field,system%nspin,sVh,sVpsl,sVxc,V_local)

    call timer_begin(LOG_CALC_TOTAL_ENERGY)
    if( PLUS_U_ON )then
      call calc_density_matrix_and_energy_plusU( spsi, ppg, info, system, energy%E_U )
    end if
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    if(calc_mode/='DFT_BAND')then
      select case(iperiodic)
      case(0); call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
      case(3); call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
      end select
    end if
    call timer_end(LOG_CALC_TOTAL_ENERGY)

    if(calc_mode=='DFT_BAND')then
      tol_esp_diff=1.0d-5
      esp_old=abs(esp_old-energy%esp)
      check_conv_esp(:,:,:)=.false.
      do ispin=1,system%nspin
      do ik=1,system%nk
         i=0
         j=0
         do iob=1,system%no
            if ( esp_old(iob,ik,ispin) <= tol_esp_diff ) then
               i=i+1
               j=max(j,iob)
               if ( iob <= nref_band ) check_conv_esp(iob,ik,ispin)=.true.
            end if
         end do !io
         if ( ispin==1 .and. ik==1 ) then
            write(*,'(/,1x,"ispin","   ik",2x,"converged bands (total, maximum band index)")')
         end if
         write(*,'(1x,2i5,2x,2i5)') ispin,ik,i,j
      end do !ik
      end do !ispin

      esp_old=energy%esp
    end if
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
         sum1 = sum1*system%Hvol/(dble(ifMST(1))*2.d0)
      else if(ispin==1)then
         sum1 = sum1*system%Hvol/dble(ifMST(1)+ifMST(2))
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
100 format(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)
101 format(1x,"iter =",i6,5x,"Total Energy =",f19.8)

    do ik=1,system%nk
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
     write(*,*) "Ne=",rNebox2*system%Hvol
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

if(calc_mode=='DFT_BAND')then
  if ( comm_is_root(nproc_id_global) ) then
  do ik=1,size(energy%esp,2)
  do iob=1,size(energy%esp,1)
    write(100,*) ik,iob,(energy%esp(iob,ik,ispin),ispin=1,system%nspin)
  end do
  end do
  end if
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then   !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(lg,mg,spsi,info,system)
      ! Experimental Implementation of Inner-Product Outputs:
      call write_prod_dk_data(lg, mg, system, info, spsi) 
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg)
   case(0)
     write(*,*) "error: yn_out_tm='y' & iperiodic=0"
  end select
end if

! force
!if(y_opt=='y') then
if(iperiodic == 3 .and. yn_ffte=='y') then
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

deallocate(idiis_sd)
call timer_end(LOG_GS_ITERATION)

end do Band_Iteration

if(calc_mode=='DFT_BAND')then
  if ( comm_is_root(nproc_id_global) ) then
  close(100)
  end if
end if

call timer_begin(LOG_DEINIT_GS_ITERATION)
if(yn_opt=='y') then
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
       write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
    end do
20  format(a5,3f16.8,i3,a3)
  end if

  if(flag_opt_conv) then
     call structure_opt_fin
     exit Multigrid_Iteration
  end if

end if
call timer_end(LOG_DEINIT_GS_ITERATION)


end do Multigrid_Iteration
if(flag_opt_conv)then
  exit Structure_Optimization_Iteration
end if
end do Structure_Optimization_Iteration


!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

! write GS: basic data
call write_band_information(system,energy)
call write_eigen(file_eigen,system,energy)
call write_info_data(Miter,system,energy,pp)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,system%hgs,iscfrt)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(iscfrt,0,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
call write_bin(ofile%dir_out_restart,lg,mg,ng,system,info,spsi,miter,mixing=mixing)
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)

contains

subroutine init_code_optimization
  implicit none
  integer :: ignum(3)

  call switch_stencil_optimization(mg%num)
  call switch_openmp_parallelization(mg%num)

  if(iperiodic==3 .and. product(pinfo%npdomain_orbital)==1) then
    ignum = mg%num
  else
    ignum = mg%num + (nd*2)
  end if
  call set_modulo_tables(ignum)

  if (comm_is_root(nproc_id_global)) then
    call optimization_log(pinfo)
  end if
end subroutine

subroutine read_bandcalc_param( lattice, nref_band, ndiv_segment, kpt, kpt_label )
  implicit none
  character(3),intent(out) :: lattice
  integer,intent(out) :: nref_band
  integer,allocatable,intent(inout) :: ndiv_segment(:)
  real(8),allocatable,intent(inout) :: kpt(:,:)  ! given in reduced coordinates in reciprocal space
  character(1),allocatable,intent(inout) :: kpt_label(:)
  integer,parameter :: unit=100
  integer :: i, num_of_segments, iformat
  if ( comm_is_root(nproc_id_global) ) then
     write(*,'(a50)') repeat("-",24)//"read_bandcalc_param(start)"
  end if
  open(unit,file='bandcalc.dat',status='old')
  read(unit,*) lattice; write(*,*) lattice
  read(unit,*) nref_band
  if ( lattice == "non" ) then
  else
     close(unit)
     if ( comm_is_root(nproc_id_global) ) then
        write(*,'(a50)') repeat("-",23)//"read_bandcalc_param(return)"
     end if
     return
  end if
  read(unit,*) num_of_segments
  allocate( ndiv_segment(num_of_segments) ); ndiv_segment=0
  allocate( kpt(3,num_of_segments+1)      ); kpt=0.0d0
  allocate( kpt_label(num_of_segments+1)  ); kpt_label=""
  read(unit,*) ndiv_segment(:)
  call check_data_format( unit, iformat )
  select case( iformat )
  case( 0 )
     do i=1,num_of_segments+1
        read(unit,*) kpt(1:3,i)
     end do
  case( 1 )
     do i=1,num_of_segments+1
        read(unit,*) kpt_label(i), kpt(1:3,i)
     end do
  end select
  close(unit)
  if ( comm_is_root(nproc_id_global) ) then
     write(*,'(a50)') repeat("-",26)//"read_bandcalc_param(end)"
  end if
end subroutine read_bandcalc_param

subroutine check_data_format( unit, iformat )
  implicit none
  integer,intent(in) :: unit
  integer,intent(out) :: iformat
  character(100) :: ccc
  character(1) :: b(4)
  read(unit,'(a)') ccc
  backspace(unit)
  read(ccc,*,END=9) b
  iformat=1 ! 4 data in one line
  return
9 iformat=0 ! 3 data
end subroutine check_data_format

subroutine get_band_kpt( kpt, nref_band, system )
   use structures, only: s_dft_system
   use salmon_parallel, only: nproc_id_global
   use salmon_communication, only: comm_is_root
   implicit none
   real(8),allocatable,intent(inout) :: kpt(:,:)
   integer,intent(out) :: nref_band ! convergence is checked up to nref_band
   type(s_dft_system),intent(in) :: system 
   real(8) :: G(3),X(3),M(3),R(3),L(3),W(3) ! XYZ coordinates of high-symmetry
   real(8) :: H(3),N(3),P(3),A(3),Q(3)      ! points in the 1st Brillouin zone
   real(8) :: al,cl ! length of the real-space lattice vectors (a- and c-axis)
   real(8) :: dk(3),k0(3),k1(3),pi,c1,c2,c3
   character(3) :: lattice
   integer,allocatable :: ndiv_segment(:)
   real(8),allocatable :: kpt_(:,:)
   character(1),allocatable :: kpt_label(:)
   integer :: nk,nnk,iseg,num_of_segments,i,ik

   if ( comm_is_root(nproc_id_global) ) then
      write(*,'(a60)') repeat("-",41)//"get_band_kpt(start)"
   end if

   pi=acos(-1.0d0)

   call read_bandcalc_param( lattice, nref_band, ndiv_segment, kpt_, kpt_label )

   if ( allocated(ndiv_segment) ) then
      if ( comm_is_root(nproc_id_global) ) then
         write(*,*) "k points are generated from 'bandcalc.dat'"
      end if
      do ik=1,size(kpt_,2)
         k0(:) = matmul( system%primitive_b, kpt_(:,ik) )
         kpt_(:,ik) = k0(:)
      end do
      num_of_segments = size( ndiv_segment )
   else if ( .not.allocated(ndiv_segment) ) then ! set default
      if ( comm_is_root(nproc_id_global) ) then
         write(*,*) "k points are generated by a default setting"
      end if
      select case( lattice )
      case( "sc" , "SC"  ); num_of_segments=5
      case( "fcc", "FCC" ); num_of_segments=5
      case( "bcc", "BCC" ); num_of_segments=5
      case( "hex", "HEX" ); num_of_segments=7
      case default
         write(*,*) "lattice=",lattice
         write(*,*)"default setting is not available for this lattice" 
         stop "stop@get_band_kpt"
      end select
      allocate( ndiv_segment(num_of_segments) ); ndiv_segment=10
      allocate( kpt_(3,num_of_segments+1)     ); kpt_=0.0d0
      allocate( kpt_label(num_of_segments+1)  ); kpt_label=""
      select case( lattice )
      case( "sc" , "SC"  ) ! G -> X -> M -> R -> G -> M  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         X=c1*(/ 0.5d0, 0.0d0, 0.0d0 /)
         M=c1*(/ 0.5d0, 0.5d0, 0.0d0 /)
         R=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=X; kpt_label(2)="X"
         kpt_(:,3)=M; kpt_label(3)="M"
         kpt_(:,4)=R; kpt_label(4)="R"
         kpt_(:,5)=G; kpt_label(5)="G"
         kpt_(:,6)=M; kpt_label(6)="M"
      case( "fcc", "FCC" ) ! G -> X -> W -> G -> L -> X  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))*sqrt(2.0d0)
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         X=c1*(/ 1.0d0, 0.0d0, 0.0d0 /)
         W=c1*(/ 1.0d0, 0.5d0, 0.0d0 /)
         L=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=X; kpt_label(2)="X"
         kpt_(:,3)=W; kpt_label(3)="W"
         kpt_(:,4)=G; kpt_label(4)="G"
         kpt_(:,5)=L; kpt_label(5)="L"
         kpt_(:,6)=X; kpt_label(6)="X"
      case( "bcc", "BCC" ) ! G -> H -> N -> P -> G -> N  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))*2.0d0/sqrt(3.0d0)
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         H=c1*(/ 0.0d0, 1.0d0, 0.0d0 /)
         N=c1*(/ 0.5d0, 0.5d0, 0.0d0 /)
         P=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=H; kpt_label(2)="H"
         kpt_(:,3)=N; kpt_label(3)="N"
         kpt_(:,4)=P; kpt_label(4)="P"
         kpt_(:,5)=G; kpt_label(5)="G"
         kpt_(:,6)=N; kpt_label(6)="N"
      case( "hex", "HEX" ) ! G -> P -> Q -> G -> A -> L -> H -> P  (7 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))
         cl=sqrt(sum(system%primitive_a(:,3)**2))
         c1=2.0d0*pi/al
         c2=1.0d0*pi/al
         c3=1.0d0*pi/cl
         G=(/ 0.0d0, 0.0d0, 0.0d0 /)
         P=c1*(/ 2.0d0/3.0d0, 0.0d0, 0.0d0 /)
         Q=c2*(/ 1.0d0, 1.0d0/sqrt(3.0d0), 0.0d0 /)
         A=c3*(/ 0.0d0, 0.0d0, 1.0d0 /)
         L=c2*(/ 1.0d0, 1.0d0/sqrt(3.0d0), c3/c2 /)
         H=c1*(/ 2.0d0/3.0d0, 0.0d0, c3/c1 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=P; kpt_label(2)="P"
         kpt_(:,3)=Q; kpt_label(3)="Q"
         kpt_(:,4)=G; kpt_label(4)="G"
         kpt_(:,5)=A; kpt_label(5)="A"
         kpt_(:,6)=L; kpt_label(6)="L"
         kpt_(:,7)=H; kpt_label(7)="H"
         kpt_(:,8)=P; kpt_label(8)="P"
      end select
   end if

   nk=system%nk
   nnk=sum( ndiv_segment(1:num_of_segments) )
   if ( mod(nnk,nk) /= 0 ) nnk=nnk-mod(nnk,nk)+nk

   allocate( kpt(3,nnk) ); kpt=0.0d0

   i=0
   do iseg=1,num_of_segments
      k0(:)=kpt_(:,iseg)
      k1(:)=kpt_(:,iseg+1)
      dk(:)=( k1(:) - k0(:) )/ndiv_segment(iseg)
      do ik=0,ndiv_segment(iseg)-1
         i=i+1
         kpt(:,i)=k0(:)+dk(:)*ik
      end do
   end do ! iseg

   if ( i < nnk ) then
      i=i+1
      kpt(:,i)=kpt(:,i-1)+dk(:)
   end if
   if ( i < nnk ) then
      do ik=i+1,nnk
         kpt(:,ik)=kpt(:,ik-1)+dk(:)
      end do
   end if

   if ( comm_is_root(nproc_id_global) ) then
      write(*,*) "Number of computed bands:",nref_band
      write(*,*) "Whole number of bands(system%no):",system%no
      write(*,*) "array size of wf for k points(system%nk):",nk
      write(*,*) "Number of segments:",num_of_segments
      write(*,*) "Total number of k points:",nnk 
      write(*,*) "k points in Cartesian coordinates:"
      do i=1,size(kpt,2)
         write(*,'(1x,i4,3f10.5)') i,kpt(:,i)
      end do
      write(*,'(a60)') repeat("-",43)//"get_band_kpt(end)"
   end if

end subroutine get_band_kpt

end subroutine main_dft
