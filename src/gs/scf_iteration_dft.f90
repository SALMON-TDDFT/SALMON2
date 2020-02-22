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

subroutine scf_iteration_dft( Miter,rion_update,sum1,  &
                              system,energy,ewald,  &
                              lg,mg,ng,  &
                              info,info_field,pinfo,  &
                              poisson,fg,  &
                              cg,mixing,  &
                              stencil,  &
                              srg,srg_ng,   &
                              spsi,shpsi,sttpsi,  &
                              srho,srho_s,  &
                              V_local,sVh,sVxc,sVpsl,xc_func,  &
                              pp,ppg,ppn,  &
                              rho_old,Vlocal_old,  &
                              band,ilevel_print )
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global
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
use input_pp_sub
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use salmon_total_energy
use init_gs, only: init_wf
use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
implicit none
integer :: ix,iy,iz,ik,is
integer :: ilevel_print !=2:print-all, =1:print-minimum, =1:no-print
integer :: iter,Miter,iob,p1,p2,p5
real(8) :: sum0,sum1
real(8) :: rNebox1,rNebox2

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_process_info),intent(in) :: pinfo
type(s_sendrecv_grid) :: srg, srg_ng
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl,rho_old,Vlocal_old
!type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_scalar) :: V_local(system%nspin),srho_s(system%nspin),sVxc(system%nspin)
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
integer :: i,j

real(8),allocatable :: esp_old(:,:,:)
real(8) :: tol_esp_diff

if(calc_mode=='DFT_BAND') then
   allocate( esp_old(system%no,system%nk,system%nspin) )
   esp_old=0d0
endif

flag_conv = .false.
sum1=1d9

!DFT_Iteration : do iter=1,nscf
DFT_Iteration : do iter=Miter+1,nscf

   if(.not. mixing%flag_mix_zero)then
   if( sum1 < threshold ) then
      flag_conv = .true.
      if( ilevel_print.ge.2 .and. comm_is_root(nproc_id_global)) then
         write(*,'(a,i6,a,e15.8)') "  #GS converged at",iter, "  :",sum1
      endif
      exit DFT_Iteration
   endif
   endif
   if(calc_mode=='DFT_BAND')then
      if(all(band%check_conv_esp)) cycle DFT_Iteration
   end if

   Miter=Miter+1

   if(theory=='dft' .and. yn_opt=='n')then
      if(step_initial_mix_zero .ge. Miter) then
         mixing%flag_mix_zero=.true.
      else
         mixing%flag_mix_zero=.false.
      endif
   endif

   if(calc_mode/='DFT_BAND')then
      ! for calc_total_energy_periodic
      rion_update = check_rion_update() .or. (iter == 1)
  
      if(temperature>=0.d0 .and. Miter>iditer_notemperature) then
         call ne2mu(energy,system)
      end if
   end if
   call copy_density(Miter,system%nspin,ng,srho_s,mixing)
   call scf_iteration_step(lg,mg,ng,system,info,info_field,pinfo,stencil,  &
                     srg,srg_ng,spsi,shpsi,srho,srho_s,  &
                     cg,ppg,V_local,  &
                     Miter,  &
                     iditer_nosubspace_diag,mixing,iter,  &
                     poisson,fg,sVh,xc_func,ppn,sVxc,energy)
   call update_vlocal(mg,system%nspin,sVh,sVpsl,sVxc,V_local)
   call timer_begin(LOG_CALC_TOTAL_ENERGY)
   if( PLUS_U_ON )then
      call calc_density_matrix_and_energy_plusU( spsi,ppg,info,system,energy%E_U )
   end if
   call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
   if(calc_mode/='DFT_BAND')then
      select case(iperiodic)
      case(0); call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
      case(3); call calc_Total_Energy_periodic(energy,ewald,system,pp,fg,rion_update)
      end select
   end if
   call timer_end(LOG_CALC_TOTAL_ENERGY)
   if(calc_mode=='DFT_BAND')then
      tol_esp_diff=1.0d-5
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
         if( ilevel_print.ge.2 ) then
         if( is==1 .and. ik==1 ) then
            write(*,'(/,1x,"ispin","   ik",2x,"converged bands (total, maximum band index)")')
         end if
         write(*,'(1x,2i5,2x,2i5)') is,ik,i,j
         end if
      end do !ik
      end do !is

      esp_old=energy%esp
   end if

   call timer_begin(LOG_WRITE_GS_RESULTS)

   select case(convergence)
   case('rho_dne')
      sum0=0d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
         sum0 = sum0 + abs(srho%f(ix,iy,iz)-rho_old%f(ix,iy,iz))
      end do
      end do
      end do
      call comm_summation(sum0,sum1,info_field%icomm_all)
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
  
   if( ilevel_print.ge.2 ) then
   if(comm_is_root(nproc_id_global)) then
      write(*,*) '-----------------------------------------------'
      select case(iperiodic)
      case(0)
         write(*,100) Miter,energy%E_tot*au_energy_ev, poisson%iterVh
      case(3)
         write(*,101) Miter,energy%E_tot*au_energy_ev
      end select
100   format(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)
101   format(1x,"iter =",i6,5x,"Total Energy =",f19.8)

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
200   format("iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec        = ",i6,e15.8)
201   format("iter and ||rho_i(ix)-rho_i-1(ix)||**2              = ",i6,e15.8)
202   format("iter and ||rho_i(ix)-rho_i-1(ix)||**2/(# of grids) = ",i6,e15.8)
203   format("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2             = ",i6,e15.8)
204   format("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2/(# of grids)= ",i6,e15.8)

   end if
   end if
   rNebox1 = 0d0 
!$OMP parallel do reduction(+:rNebox1) private(iz,iy,ix)
   do iz=ng%is(3),ng%ie(3)
   do iy=ng%is(2),ng%ie(2)
   do ix=ng%is(1),ng%ie(1)
      rNebox1 = rNebox1 + srho%f(ix,iy,iz)
   end do
   end do
   end do
   call comm_summation(rNebox1,rNebox2,info_field%icomm_all)
   if( ilevel_print.ge.2 ) then
   if(comm_is_root(nproc_id_global))then
      write(*,*) "Ne=",rNebox2*system%Hvol
   end if
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

   if(theory=='dft' .and. yn_opt=='n')then
   if((checkpoint_interval >= 1) .and. (mod(Miter,checkpoint_interval)==0)) then
      call checkpoint_gs(lg,mg,ng,system,info,spsi,Miter,mixing)
      if(comm_is_root(nproc_id_global)) write(*,'(a)')"  checkpoint data is printed"
      call comm_sync_all
   endif
   endif

end do DFT_Iteration

if(.not.flag_conv) then
   if( ilevel_print.ge.1 .and. comm_is_root(nproc_id_global)) then
      write(*,'(a,e15.8)') "  #GS does not converged :",sum1
   endif
endif



end subroutine scf_iteration_dft
