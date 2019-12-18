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

subroutine initialization1_dft( system, energy, stencil, fg, poisson,  &
                                lg, mg, ng,  &
                                pinfo, info, info_field,  &
                                srg, srg_ng,  &
                                srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                                spsi, shpsi, sttpsi,  &
                                pp, ppg,  &
                                ofl )
use math_constants, only: pi, zi
use structures
use parallelization, only: nproc_id_global!,nproc_group_global
use communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use salmon_pp, only: calc_nlcc
use timer
use scf_iteration_sub
use writefield
use global_variables_scf
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use input_pp_sub
use prep_pp_sub
use checkpoint_restart_sub
use hamiltonian
use salmon_total_energy
implicit none
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
!type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl !,rho_old,Vlocal_old
!type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_scalar) :: V_local(system%nspin),srho_s(system%nspin),sVxc(system%nspin)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ofile)  :: ofl

integer,parameter :: Nd = 4

integer :: jspin

!call init_dft(nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofl)

call init_code_optimization

allocate( energy%esp(system%no,system%nk,system%nspin) ); energy%esp=0.0d0

!allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))
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
call init_ps(lg,mg,system,info,info_field,fg,poisson,pp,ppg,sVpsl)
call calc_nlcc(pp, system, mg, ppn)  !setup NLCC term from pseudopotential
if(comm_is_root(nproc_id_global)) then
  write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
  write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
end if    

select case(iperiodic)
case(0)
  call allocate_orbital_real(system%nspin,mg,info,spsi)
  call allocate_orbital_real(system%nspin,mg,info,shpsi)
case(3)
  call allocate_orbital_complex(system%nspin,mg,info,spsi)
  call allocate_orbital_complex(system%nspin,mg,info,shpsi)
  call allocate_orbital_complex(system%nspin,mg,info,sttpsi)
end select

if(stencil%if_orthogonal) then
  if(comm_is_root(nproc_id_global)) write(*,*) "orthogonal cell: using al"
else
  if(comm_is_root(nproc_id_global)) write(*,*) "non-orthogonal cell: using al_vec[1,2,3]"
end if


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

end subroutine initialization1_dft


subroutine initialization2_dft( Miter, nspin, rion_update,  &
                                system,energy,stencil,fg,poisson,  &
                                lg,mg,ng,  &
                                info,info_field,  &
                                srg,srg_ng,  &
                                srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                                spsi,shpsi,sttpsi,  &
                                pp,ppg,ppn,  &
                                xc_func,mixing )
use math_constants, only: pi, zi
use structures
use parallelization, only: nproc_id_global
use communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use global_variables_scf
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
type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl
!type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_scalar) :: V_local(system%nspin),srho_s(system%nspin),sVxc(system%nspin)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_mixing) :: mixing

logical :: rion_update
integer :: Miter,jspin, nspin

  select case(method_min)
  case('cg')   ; continue
  case default ; stop 'Specify "cg" for method_min.'
  end select

  select case(method_mixing)
  case ('simple','broyden') ; continue
  case default ; stop 'Specify either "simple" or "broyden" for method_mixing.'
  end select


  nspin = system%nspin

  mixing%num_rho_stock = 21
  call init_mixing(nspin,mg,mixing)

  if (yn_restart == 'y') then
    ! restart from binary
    call restart_gs(lg,mg,ng,system,info,spsi,Miter,mixing=mixing)
  else
    ! new calculation
    Miter = 0        ! Miter: Iteration counter set to zero
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


end subroutine initialization2_dft

!====================================
subroutine initialization_dft_md( Miter, rion_update,  &
                                system,md,energy,stencil,fg,poisson,  &
                                lg,mg,ng,  &
                                info,info_field,  &
                                srg,srg_ng,  &
                                srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                                spsi,shpsi,sttpsi,  &
                                pp,ppg,ppn,  &
                                xc_func,mixing )
  use math_constants, only: pi, zi
  use const, only: hartree2J,kB
  use structures
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_xc
  use timer
  use scf_iteration_sub
  use density_matrix, only: calc_density
  use writefield
  use global_variables_scf
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
  use md_sub, only: init_md
  implicit none
  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_rgrid) :: ng
  type(s_orbital_parallel) :: info
  type(s_field_parallel) :: info_field
  type(s_sendrecv_grid) :: srg, srg_ng
  type(s_orbital) :: spsi,shpsi,sttpsi
  type(s_dft_system) :: system
  type(s_md) :: md
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_xc_functional) :: xc_func
  type(s_scalar) :: srho,sVh,sVpsl,rho_old,Vlocal_old
  type(s_scalar) :: V_local(system%nspin),srho_s(system%nspin),sVxc(system%nspin)
  type(s_reciprocal_grid) :: fg
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_dft_energy) :: energy
  type(s_mixing) :: mixing
  type(s_cg)     :: cg
  type(s_band_dft) ::band
  
  logical :: rion_update
  integer :: Miter,ix,iy,iz
  real(8) :: sum1

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

  !-------------- SCF Iteration ----------------
  !Iteration loop for SCF (DFT_Iteration)
  Miter=0
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
                          band,1 )

  call init_md(system,md)
  call calc_force(system,pp,fg,info,mg,stencil,srg,ppg,spsi)

  md%Uene0 = energy%E_tot
  md%E_tot0= energy%E_tot + md%Tene

  md%Uene  = energy%E_tot
  md%E_tot = energy%E_tot + md%Tene
  md%Htot  = energy%E_tot + md%E_nh  !for NHC

  md%Enh_gkTlns = 0d0
  md%E_nh       = 0d0
  if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
     md%gkT = 3d0*natom * kB/hartree2J * temperature0_ion_k
     md%Qnh = md%gkT * thermostat_tau**2d0
  endif


end subroutine initialization_dft_md
