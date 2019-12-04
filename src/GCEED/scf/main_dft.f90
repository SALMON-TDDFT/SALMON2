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
use inputoutput, only: au_length_aa, au_energy_ev
use salmon_parallel, only: nproc_id_global,nproc_group_global
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
use band_dft_sub
use init_gs, only: init_wf
implicit none
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
real(8) :: sum1
character(100) :: file_atoms_coo, comment_line

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
type(s_band_dft) ::band

logical :: rion_update
logical :: flag_opt_conv
integer :: iopt,nopt_max
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='DFT_BAND'.and.iperiodic/=3) return

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)

call convert_input_scf(file_atoms_coo)

! please move folloings into initialization_dft 
call init_dft(nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)
allocate( srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin) )


call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          pinfo, info, info_field,  &
                          srg, srg_ng,  &
                          srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofile )

call initialization2_dft( Miter, nspin, rion_update,  &
                          system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          info, info_field,   &
                          srg, srg_ng,  &
                          srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )


call timer_end(LOG_INIT_GS)

!---------------------------------------- Opt Iteration

if(yn_opt=='y')then
   call structure_opt_ini(natom)
   flag_opt_conv=.false.
   nopt_max = nopt

   write(comment_line,10) 0
   call write_xyz(comment_line,"new","r  ",system)
10 format("#opt iteration step=",i5)
else
   nopt_max = 1
end if

Structure_Optimization_Iteration : do iopt=1,nopt_max

if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  rion_update = .true.
  call dealloc_init_ps(ppg)
  call init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)
  call calc_nlcc(pp, system, mg, ppn)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(theory=='DFT_BAND')then
   call init_band_dft(system,band)
   iter_band_kpt_end    = band%num_band_kpt
   iter_band_kpt_stride = system%nk
else
   iter_band_kpt_end    = 1
   iter_band_kpt_stride = 1
end if

Band_Iteration : do iter_band_kpt= 1, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='DFT_BAND')then
   call calc_band_write(iter_band_kpt,system,band,info)
end if


call timer_begin(LOG_INIT_GS_ITERATION)

if(allocated(rho_old%f))    deallocate(rho_old%f)
if(allocated(Vlocal_old%f)) deallocate(Vlocal_old%f)
call allocate_scalar(ng,rho_old)
call allocate_scalar(ng,Vlocal_old)

!$OMP parallel do private(iz,iy,ix)
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
   rho_old%f(ix,iy,iz)   = srho%f(ix,iy,iz)
   Vlocal_old%f(ix,iy,iz)= V_local(1)%f(ix,iy,iz)
end do
end do
end do

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        system,energy,  &
                        lg,mg,ng,  &
                        info,info_field,  &
                        poisson,fg,  &
                        cg,mixing,  &
                        stencil,  &
                        srg,srg_ng,   &
                        spsi,shpsi,sttpsi,  &
                        srho,srho_s,  &
                        V_local,sVh,sVxc,sVpsl,xc_func,  &
                        pp,ppg,ppn,  &
                        rho_old,Vlocal_old,  &
                        band, 2 )


if(theory=='DFT_BAND')then
   call write_band(system,energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
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
   call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi)
   if(comm_is_root(nproc_id_global))then
      write(*,*) "===== force ====="
      do iatom=1,natom
         select case(unit_system)
         case('au','a.u.'); write(*,300)iatom,(system%Force(ix,iatom),ix=1,3)
         case('A_eV_fs'  ); write(*,300)iatom,(system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
         end select
      end do
300   format(i6,3e16.8)
   end if

call timer_end(LOG_GS_ITERATION)

end do Band_Iteration


call timer_begin(LOG_DEINIT_GS_ITERATION)
if(yn_opt=='y') then
   call structure_opt_check(natom,iopt,flag_opt_conv,system%Force)
   if(.not.flag_opt_conv) call structure_opt(natom,iopt,system)
   !! Rion is old variables to be removed 
   !! but currently it is used in many subroutines.
   Rion(:,:) = system%Rion(:,:) 

   write(comment_line,10) iopt
   call write_xyz(comment_line,"add","r  ",system)

   if(comm_is_root(nproc_id_global))then
      write(*,*) "atomic coordinate"
      do iatom=1,natom
         write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
      end do
20    format(a5,3f16.8,i3,a3)
   end if

   if(flag_opt_conv) call structure_opt_fin

end if
call timer_end(LOG_DEINIT_GS_ITERATION)


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
if(yn_out_dns =='y') call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,system%hgs)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(0,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
call write_bin(ofile%dir_out_restart,lg,mg,ng,system,info,spsi,Miter,mixing=mixing)
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)


end subroutine main_dft
