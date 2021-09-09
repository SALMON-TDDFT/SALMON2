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

subroutine scf_iteration_dft( Miter,rion_update,sum1,  &
                              system,energy,ewald,  &
                              lg,mg,  &
                              info,  &
                              poisson,fg,  &
                              cg,mixing,  &
                              stencil,  &
                              srg,srg_scalar,   &
                              spsi,shpsi,sttpsi,  &
                              rho,rho_jm,rho_s,  &
                              V_local,Vh,Vxc,Vpsl,xc_func,  &
                              pp,ppg,ppn,  &
                              rho_old,Vlocal_old,  &
                              band,ilevel_print )
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global,adjust_elapse_time
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use total_energy
use init_gs, only: init_wf
use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
use noncollinear_module, only: calc_magnetization
implicit none
integer :: ix,iy,iz,ik,is
integer :: ilevel_print !=3:print-all
                        !=2:print-only energy & convergence
                        !=1:print-minimum
                        !=0:no-print
integer :: iter,Miter,iob,p1,p2,p5
real(8) :: sum0,sum1
real(8) :: rNebox1,rNebox2

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,rho_jm,Vh,Vpsl,rho_old,Vlocal_old
!type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_scalar) :: V_local(system%nspin),rho_s(system%nspin),Vxc(system%nspin)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ewald_ion_ion) :: ewald
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_band_dft) :: band

logical :: rion_update, flag_conv
integer :: i,j, icnt_conv_nomix
logical :: is_checkpoint_iter, is_shutdown_time

real(8),allocatable :: esp_old(:,:,:)
real(8) :: ene_gap, magnetization(3)

if(calc_mode=='DFT_BAND') then
   allocate( esp_old(system%no,system%nk,system%nspin) )
   esp_old=0d0
endif

if(nscf_init_mix_zero.gt.1)then
   icnt_conv_nomix = 0
   mixing%flag_mix_zero = .true.
   DFT_NoMix_Iteration : do iter=1,nscf_init_mix_zero

      if(yn_jm=='n') rion_update = check_rion_update() .or. (iter == 1)
      call copy_density(iter,system%nspin,mg,rho_s,mixing)
      call scf_iteration_step(lg,mg,system,info,stencil,  &
                     srg,srg_scalar,spsi,shpsi,rho,rho_jm,rho_s,  &
                     cg,ppg,V_local,  &
                     iter,  &
                     nscf_init_no_diagonal, mixing, iter,  &
                     poisson,fg,Vh,xc_func,ppn,Vxc,energy)
      call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)
      call timer_begin(LOG_CALC_TOTAL_ENERGY)
      call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
      select case(iperiodic)
      case(0); call calc_Total_Energy_isolated(system,info,lg,mg,pp,ppg,fg,poisson,rho_s,Vh,Vxc,rion_update,energy)
      case(3); call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
      end select
      call get_band_gap(system,energy,ene_gap)
      call timer_end(LOG_CALC_TOTAL_ENERGY)
      if(comm_is_root(nproc_id_global)) then
         select case(iperiodic)
         case(0); write(*,300) iter, energy%E_tot*au_energy_ev, ene_gap*au_energy_ev, poisson%iterVh
         case(3)
           if(yn_jm=='n') then
             write(*,301) iter, energy%E_tot*au_energy_ev, ene_gap*au_energy_ev
           else
             write(*,302) iter, ene_gap*au_energy_ev
           end if
         end select
300      format(2x,"no-mixing iter=",i6,5x,"Total Energy=",f19.8,5x,"Gap=",f15.8,5x,"Vh iter=",i4)
301      format(2x,"no-mixing iter=",i6,5x,"Total Energy=",f19.8,5x,"Gap=",f15.8)
302      format(2x,"no-mixing iter=",i6,                         5x,"Gap=",f15.8)
      endif
      !(convergence: energy gap is over specified energy)
      if(ene_gap .ge. conv_gap_mix_zero) then
         icnt_conv_nomix = icnt_conv_nomix + 1
         if(icnt_conv_nomix==5) then
            if(comm_is_root(nproc_id_global)) write(*,*) "  converged no-mixing iteration"
            exit
         endif
      else
         icnt_conv_nomix = 0
      endif

   end do DFT_NoMix_Iteration
   mixing%flag_mix_zero = .false.
endif

flag_conv = .false.
sum1=1d9

!DFT_Iteration : do iter=1,nscf
DFT_Iteration : do iter=Miter+1,nscf

   if( sum1 < threshold ) then
      flag_conv = .true.
      if( ilevel_print.ge.3 .and. comm_is_root(nproc_id_global)) then
         write(*,'(a,i6,a,e15.8)') "  #GS converged at",iter, "  :",sum1
      endif
      exit DFT_Iteration
   endif
   if(calc_mode=='DFT_BAND')then
      if(all(band%check_conv_esp)) then
        if ( comm_is_root(nproc_id_global) ) write(*,*) "cycle!!! : iter=",iter
        cycle DFT_Iteration
      end if
   end if

   Miter=Miter+1

   if(calc_mode/='DFT_BAND')then
      ! for calc_total_energy_periodic
      if(yn_jm=='n') rion_update = check_rion_update() .or. (iter == 1)
  
      if(temperature>=0.d0 .and. Miter>nscf_init_redistribution) then
         call ne2mu(energy,system,ilevel_print)
      end if
   end if
   call copy_density(Miter,system%nspin,mg,rho_s,mixing)
   call scf_iteration_step(lg,mg,system,info,stencil,  &
                     srg,srg_scalar,spsi,shpsi,rho,rho_jm,rho_s,  &
                     cg,ppg,V_local,  &
                     Miter,  &
                     nscf_init_no_diagonal, mixing, iter,  &
                     poisson,fg,Vh,xc_func,ppn,Vxc,energy)
   call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)
   call timer_begin(LOG_CALC_TOTAL_ENERGY)
   if( PLUS_U_ON )then
      call calc_density_matrix_and_energy_plusU( spsi,ppg,info,system,energy%E_U )
   end if
   if(yn_spinorbit=='y') then
     call calc_magnetization(system,mg,info,magnetization)
   end if
   call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
   call get_band_gap(system,energy,ene_gap)
   if(calc_mode/='DFT_BAND')then
      select case(iperiodic)
      case(0); call calc_Total_Energy_isolated(system,info,lg,mg,pp,ppg,fg,poisson,rho_s,Vh,Vxc,rion_update,energy)
      case(3); call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
      end select
   end if
   call timer_end(LOG_CALC_TOTAL_ENERGY)
   if(calc_mode=='DFT_BAND')then
      esp_old=abs(esp_old-energy%esp)
      band%check_conv_esp(:,:,:)=.false.
      do is=1,system%nspin
      do ik=1,system%nk
         i=0
         j=0
         do iob=1,system%no
            if ( esp_old(iob,ik,is) <= tol_esp_diff ) then
               i=i+1
               j=max(j,iob)
               if( iob <= band%nref_band ) band%check_conv_esp(iob,ik,is)=.true.
            end if
         end do !io
         if( ilevel_print.ge.3 ) then
         if( is==1 .and. ik==1 ) then
            write(*,'(/,1x,"ispin","   ik",2x,"converged bands (total, maximum band index)")')
         end if
         write(*,'(1x,2i5,2x,2i5)') is,ik,i,j
         end if
      end do !ik
      end do !is

      esp_old=energy%esp
   end if

   if( calc_mode/='DFT_BAND' )then
   call timer_begin(LOG_WRITE_GS_RESULTS)

   select case(convergence)
   case('rho_dne')
      sum0=0d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3) 
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         sum0 = sum0 + abs(rho%f(ix,iy,iz)-rho_old%f(ix,iy,iz))
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info%icomm_r)
      if(system%nspin==1)then
         sum1 = sum1*system%Hvol/dble(nelec)
      else if(system%nspin==2)then
         if(sum(nelec_spin(:))>0)then
            sum1 = sum1*system%Hvol/dble(sum(nelec_spin(:)))
         else 
            sum1 = sum1*system%Hvol/dble(nelec)
         end if
      end if
   case('norm_rho','norm_rho_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3) 
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         sum0 = sum0 + (rho%f(ix,iy,iz)-rho_old%f(ix,iy,iz))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info%icomm_r)
      if(convergence=='norm_rho_dng')then
         sum1 = sum1/dble(lg%num(1)*lg%num(2)*lg%num(3))
      end if
   case('norm_pot','norm_pot_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3) 
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         sum0 = sum0 + (V_local(1)%f(ix,iy,iz)-Vlocal_old%f(ix,iy,iz))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info%icomm_r)
      if(convergence=='norm_pot_dng')then
         sum1 = sum1/dble(lg%num(1)*lg%num(2)*lg%num(3))
      end if
   end select
   
   if( ilevel_print.ge.3 ) then
   if(comm_is_root(nproc_id_global)) then
      write(*,*) '-----------------------------------------------'
      select case(iperiodic)
      case(0)
         write(*,100) Miter,energy%E_tot*au_energy_ev, ene_gap*au_energy_ev, poisson%iterVh
      case(3)
        if(yn_jm=='n') then
          write(*,101) Miter,energy%E_tot*au_energy_ev, ene_gap*au_energy_ev
        else
          write(*,102) Miter, ene_gap*au_energy_ev
        end if
      end select
100   format(1x,"iter=",i6,5x,"Total Energy=",f19.8,5x,"Gap=",f15.8,5x,"Vh iter=",i4)
101   format(1x,"iter=",i6,5x,"Total Energy=",f19.8,5x,"Gap=",f15.8)
102   format(1x,"iter=",i6,5x,                         "Gap=",f15.8)

      do is=1,system%nspin
         if(system%nspin==2.and.is==1) write(*,*) "for up-spin"
         if(system%nspin==2.and.is==2) write(*,*) "for down-spin"
         do ik=1,system%nk
            if(ik<=3)then
               if(iperiodic==3) write(*,*) "k=",ik
               do p5=1,(system%no+3)/4
                  p1=4*(p5-1)+1
                  p2=4*p5 ; if ( p2 > system%no ) p2=system%no
                  write(*,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,ik,is)*au_energy_ev,iob=p1,p2)
               end do
               if(iperiodic==3) write(*,*) 
            end if
         end do
      end do

      select case(convergence)
      case('rho_dne' )     ; write(*,200) Miter, sum1
      case('norm_rho')     ; write(*,201) Miter, sum1/au_length_aa**6
      case('norm_rho_dng') ; write(*,202) Miter, sum1/au_length_aa**6
      case('norm_pot')     ; write(*,203) Miter, sum1*(au_energy_ev)**2/au_length_aa**6
      case('norm_pot_dng') ; write(*,204) Miter, sum1*(au_energy_ev)**2/au_length_aa**6
      end select
200   format(1x,"iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec        = ",i6,e15.8)
201   format(1x,"iter and ||rho_i(ix)-rho_i-1(ix)||**2              = ",i6,e15.8)
202   format(1x,"iter and ||rho_i(ix)-rho_i-1(ix)||**2/(# of grids) = ",i6,e15.8)
203   format(1x,"iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2             = ",i6,e15.8)
204   format(1x,"iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2/(# of grids)= ",i6,e15.8)

      if(yn_spinorbit=='y') then
        write(*,'(1x,"Magnetization= ",3(e15.8,2x))') magnetization(1:3)
      end if
   end if !comm_is_root

   else if( ilevel_print==2 ) then

   if(comm_is_root(nproc_id_global)) then
      select case(iperiodic)
      case(0)
        write(*,400) Miter,energy%E_tot*au_energy_ev, ene_gap*au_energy_ev, sum1,poisson%iterVh
      case(3)
        if(yn_jm=='n') then
          write(*,401) Miter,energy%E_tot*au_energy_ev, ene_gap*au_energy_ev, sum1
        else
          write(*,402) Miter,                           ene_gap*au_energy_ev, sum1
        end if  
      end select
400   format(5x,"#SCF ",i6,3x,"E(total)=",f19.8,3x,"Gap=",f15.8,3x,"conv[au]=",e15.7,3x,"Vh iter=",i4)
401   format(5x,"#SCF ",i6,3x,"E(total)=",f19.8,3x,"Gap=",f15.8,3x,"conv[au]=",e15.7)
402   format(5x,"#SCF ",i6,3x,                     "Gap=",f15.8,3x,"conv[au]=",e15.7)
   endif

   end if  !ilevel_print

! modification of mixing rate for auto_mixing
   if(yn_auto_mixing=='y')then
     call check_mixing_half(Miter,sum1,mixing)
   end if
 
   rNebox1 = 0d0 
!$OMP parallel do reduction(+:rNebox1) private(iz,iy,ix)
   do iz=mg%is(3),mg%ie(3)
   do iy=mg%is(2),mg%ie(2)
   do ix=mg%is(1),mg%ie(1)
      rNebox1 = rNebox1 + rho%f(ix,iy,iz)
   end do
   end do
   end do
   call comm_summation(rNebox1,rNebox2,info%icomm_r)
   if( ilevel_print.ge.3 ) then
   if(comm_is_root(nproc_id_global))then
      write(*,*) "Ne=",rNebox2*system%Hvol
   end if
   end if
   call timer_end(LOG_WRITE_GS_RESULTS)
!$OMP parallel do private(iz,iy,ix)
   do iz=mg%is(3),mg%ie(3)
   do iy=mg%is(2),mg%ie(2)
   do ix=mg%is(1),mg%ie(1)
      rho_old%f(ix,iy,iz)    = rho%f(ix,iy,iz)
      Vlocal_old%f(ix,iy,iz) = V_local(1)%f(ix,iy,iz)
   end do
   end do
   end do

   end if !calc_mode/=DFT_BAND

   if(theory=='dft' .and. yn_opt=='n')then
     is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(Miter,checkpoint_interval) == 0)
     is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

     if (is_checkpoint_iter .or. is_shutdown_time) then
       if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
         print *, 'shutdown the calculation, iter =', Miter
       end if

       call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing)
       if(comm_is_root(nproc_id_global)) write(*,'(a)')"  checkpoint data is printed"
       call comm_sync_all

       if (is_shutdown_time) then
         exit DFT_Iteration
       end if
     endif
   endif

end do DFT_Iteration

if(calc_mode/='DFT_BAND')then
if(.not.flag_conv) then
   if( ilevel_print.ge.1 .and. comm_is_root(nproc_id_global)) then
      write(*,'(a,e15.8)') "  #GS does not converged :",sum1
   endif
endif
endif

end subroutine scf_iteration_dft

subroutine  get_band_gap(system,energy,gap)
  use structures
  use salmon_global, only: nelec,yn_spinorbit
  use inputoutput, only: au_energy_ev
  use communication, only: comm_is_root
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_dft_energy),intent(in) :: energy
  integer :: ik, index_vbm
  real(8) :: gap
  real(8),dimension(system%nk) :: esp_vb_max, esp_cb_min
  if( yn_spinorbit=='y' )then
     index_vbm=nelec
  else
     index_vbm=nelec/2
  end if
  if( index_vbm >= system%no )then
     gap = 1d99
     return
  endif
  do ik=1,system%nk
     esp_vb_max(ik)=maxval(energy%esp(1:index_vbm,ik,:))
     esp_cb_min(ik)=minval(energy%esp(index_vbm+1:system%no,ik,:))
  end do
  ! 'Fundamental gap' (see write_band_information)
  gap = minval(esp_cb_min(:))-maxval(esp_vb_max(:))
end subroutine get_band_gap
