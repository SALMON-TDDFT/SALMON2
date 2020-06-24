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

#include "config.h"
module initialization_dft
  implicit none
  
contains

subroutine initialization1_dft( system, energy, stencil, fg, poisson,  &
                                lg, mg,  &
                                info,  &
                                srg, srg_scalar,  &
                                rho, rho_s, Vh, V_local, Vpsl, Vxc,  &
                                spsi, shpsi, sttpsi,  &
                                pp, ppg, ppn,  &
                                ofl )
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global!,nproc_group_global
use communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use salmon_pp, only: calc_nlcc, read_pslfile
use timer
use scf_iteration_sub
use writefield
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use prep_pp_sub
use checkpoint_restart_sub
use hamiltonian
use total_energy
implicit none
type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
!type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,Vh,Vpsl !,rho_old,Vlocal_old
!type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_scalar) :: V_local(system%nspin),rho_s(system%nspin),Vxc(system%nspin)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ofile)  :: ofl

integer,parameter :: Nd = 4

integer :: jspin

!call init_dft(nproc_group_global,info,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl)

call init_scalapack(info,system)
call init_eigenexa(info,system)

call init_code_optimization

allocate( energy%esp(system%no,system%nk,system%nspin) ); energy%esp=0.0d0

!allocate(rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin))
call allocate_scalar(mg,rho)
call allocate_scalar(mg,Vh)
call allocate_scalar(mg,Vpsl)
do jspin=1,system%nspin
  call allocate_scalar(mg,rho_s(jspin))
  call allocate_scalar(mg,V_local(jspin))
  call allocate_scalar(mg,Vxc(jspin))
end do
call read_pslfile(system,pp)
call init_ps(lg,mg,system,info,fg,poisson,pp,ppg,Vpsl)
call calc_nlcc(pp, system, mg, ppn)  !setup NLCC term from pseudopotential
if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
  write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
  write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
end if

if(system%if_real_orbital) then
  call allocate_orbital_real(system%nspin,mg,info,spsi)
  call allocate_orbital_real(system%nspin,mg,info,shpsi)
  call allocate_orbital_real(system%nspin,mg,info,sttpsi)
else
  call allocate_orbital_complex(system%nspin,mg,info,spsi)
  call allocate_orbital_complex(system%nspin,mg,info,shpsi)
  call allocate_orbital_complex(system%nspin,mg,info,sttpsi)
end if

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

  if(iperiodic==3 .and. product(info%nprgrid)==1) then
    ignum = mg%num
  else
    ignum = mg%num + (nd*2)
  end if
  call set_modulo_tables(ignum)

  if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
    call optimization_log(info)
  end if
end subroutine

end subroutine initialization1_dft


subroutine initialization2_dft( Miter, nspin, rion_update,  &
                                system,energy,ewald,stencil,fg,poisson,&
                                lg,mg,info,  &
                                srg,srg_scalar,  &
                                rho, rho_s, Vh,V_local, Vpsl, Vxc,  &
                                spsi,shpsi,sttpsi,  &
                                pp,ppg,ppn,  &
                                xc_func,mixing )
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global
use communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
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
use band_dft_sub
use init_gs, only: init_wf
implicit none
type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,Vh,Vpsl
!type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_scalar) :: V_local(system%nspin),rho_s(system%nspin),Vxc(system%nspin)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ewald_ion_ion) :: ewald
type(s_mixing) :: mixing

logical :: rion_update
integer :: Miter,jspin, nspin,i,ix,iy,iz
real(8) :: rNe0,rNe

  select case(method_min)
  case('cg')   ; continue
  case default ; stop 'Specify "cg" for method_min.'
  end select

  select case(method_mixing)
  case ('simple','broyden','pulay') ; continue
  case default ; stop 'Specify any one of "simple" or "broyden" or "pulay" for method_mixing.'
  end select


  nspin = system%nspin
  spsi%update_zwf_overlap = .false.

  mixing%num_rho_stock = 21
  call init_mixing(nspin,mg,mixing)

  ! restart from binary
  if (yn_restart == 'y') then
     call restart_gs(lg,mg,system,info,spsi,Miter,mixing=mixing)
     if(read_gs_restart_data=='wfn') then
        call calc_density(system,rho_s,spsi,info,mg)
     else 
        i = mixing%num_rho_stock+1
        !rho%f      => mg (orbital domain allocation)
        !rho_s(j)%f => mg (orbital domain allocation)
        !mixing%rho_in(i)%f => mg (density domain allocation)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           rho_s(1)%f(ix,iy,iz) = mixing%rho_in(i)%f(ix,iy,iz)
        end do
        end do
        end do
        if(system%nspin==2) then
           do jspin=1,system%nspin
              do iz=mg%is(3),mg%ie(3)
              do iy=mg%is(2),mg%ie(2)
              do ix=mg%is(1),mg%ie(1)
                 rho_s(jspin)%f(ix,iy,iz) = mixing%rho_s_in(i,jspin)%f(ix,iy,iz)
              end do
              end do
              end do
           end do
        endif

        if( read_gs_restart_data=='rho_inout'.or. &
            read_gs_restart_data=='rho'         )then
           call init_wf(lg,mg,system,info,spsi)
        endif
     endif

     if(yn_reset_step_restart=='y' .or. &
        (read_gs_restart_data=='rho'.or.read_gs_restart_data=='rho_inout') ) Miter=0

  else
    ! new calculation
    Miter = 0        ! Miter: Iteration counter set to zero
    call init_wf(lg,mg,system,info,spsi)
    call calc_density(system,rho_s,spsi,info,mg)
  end if

  rho%f = 0d0
  do jspin=1,nspin
     rho%f = rho%f + rho_s(jspin)%f
  end do
  
  if(read_gs_dns_cube == 'y') then
    call read_dns(lg,mg,rho%f) ! cube file only
    rNe0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rNe0 = rNe0 + rho%f(ix,iy,iz) *system%Hvol
    end do
    end do
    end do
    call comm_summation(rNe0,rNe,info%icomm_r)
    rho%f = rho%f *(dble(nelec)/rNe)
    if(system%nspin==1) then
      rho_s(1)%f = rho%f
    else if(system%nspin==2) then
      rho_s(1)%f = rho%f/2d0
      rho_s(2)%f = rho%f/2d0
    end if
  end if

  call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,ppn,info,spsi,stencil,Vxc,energy%E_xc)
  call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)

  select case(iperiodic)
  case(0)
     ewald%yn_bookkeep='n'  !to be input keyword??
  case(3)
     ewald%yn_bookkeep='y'
     call  init_nion_div(system,lg,mg,info)
  end select
  if(ewald%yn_bookkeep=='y') call init_ewald(system,info,ewald)

  call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
  rion_update = .true. ! it's first calculation
  select case(iperiodic)
  case(0)
     call calc_Total_Energy_isolated(system,info,mg,pp,rho_s,Vh,Vxc,rion_update,energy)
  case(3)
     call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
  end select


end subroutine initialization2_dft

!====================================
subroutine initialization_dft_md( Miter, rion_update,  &
                                system,md,energy,ewald,stencil,fg,poisson,&
                                lg,mg,  &
                                info,  &
                                srg,srg_scalar,  &
                                rho, rho_s, Vh,V_local, Vpsl, Vxc,  &
                                spsi,shpsi,sttpsi,  &
                                pp,ppg,ppn,  &
                                xc_func,mixing )
  use math_constants, only: pi, zi
  use const, only: hartree2J,kB
  use structures
  use inputoutput
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_xc
  use timer
  use scf_iteration_sub
  use density_matrix, only: calc_density
  use writefield
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
  use band_dft_sub
  use md_sub, only: init_md
  implicit none
  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_parallel_info) :: info
  type(s_sendrecv_grid) :: srg, srg_scalar
  type(s_orbital) :: spsi,shpsi,sttpsi
  type(s_dft_system) :: system
  type(s_md) :: md
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_xc_functional) :: xc_func
  type(s_scalar) :: rho,Vh,Vpsl,rho_old,Vlocal_old
  type(s_scalar) :: V_local(system%nspin),rho_s(system%nspin),Vxc(system%nspin)
  type(s_reciprocal_grid) :: fg
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_dft_energy) :: energy
  type(s_ewald_ion_ion) :: ewald
  type(s_mixing) :: mixing
  type(s_cg)     :: cg
  type(s_band_dft) ::band

  logical :: rion_update
  integer :: Miter,ix,iy,iz
  real(8) :: sum1

  if(allocated(rho_old%f))    deallocate(rho_old%f)
  if(allocated(Vlocal_old%f)) deallocate(Vlocal_old%f)
  call allocate_scalar(mg,rho_old)
  call allocate_scalar(mg,Vlocal_old)

!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
     rho_old%f(ix,iy,iz)   = rho%f(ix,iy,iz)
     Vlocal_old%f(ix,iy,iz)= V_local(1)%f(ix,iy,iz)
  end do
  end do
  end do

  !-------------- SCF Iteration ----------------
  !Iteration loop for SCF (DFT_Iteration)
  Miter=0
  call scf_iteration_dft( Miter,rion_update,sum1,  &
                          system,energy,ewald,  &
                          lg,mg,  &
                          info,  &
                          poisson,fg,  &
                          cg,mixing,  &
                          stencil,  &
                          srg,srg_scalar,   &
                          spsi,shpsi,sttpsi,  &
                          rho,rho_s,  &
                          V_local,Vh,Vxc,Vpsl,xc_func,  &
                          pp,ppg,ppn,  &
                          rho_old,Vlocal_old,  &
                          band, 1 )

  call init_md(system,md)
  call calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,spsi,ewald)

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

!===================================================================================================================================

subroutine init_scalapack(info,system)
  use salmon_global, only: yn_scalapack
#ifdef USE_SCALAPACK
  use scalapack_module
  use communication, only: comm_is_root,comm_sync_all
#endif
  use structures
  implicit none
  type(s_parallel_info)         :: info
  type(s_dft_system),intent(in) :: system
#ifdef USE_SCALAPACK
  integer :: n
#endif

  if (yn_scalapack == 'y') then
#ifdef USE_SCALAPACK
    call create_gridmap(info)
    n = system%no
    call init_blacs(info,n)
    if (n /= system%no) then
      if (comm_is_root(info%id_rko)) then
        print '(A,I6)', '[FATAL ERROR] nstate must be ',n
      end if
      call comm_sync_all
      stop
    end if
#endif
  end if
end subroutine

!===================================================================================================================================

subroutine init_eigenexa(info,system)
  use salmon_global, only: yn_eigenexa
#ifdef USE_EIGENEXA
  use eigenexa_module, only: init_eigenexa_mod => init_eigenexa
  use communication, only: comm_is_root,comm_sync_all
#endif
  use structures
  implicit none
  type(s_parallel_info)         :: info
  type(s_dft_system),intent(in) :: system

  if (yn_eigenexa == 'y') then
#ifdef USE_EIGENEXA
    call init_eigenexa_mod(info,system%no)
#endif
  end if
end subroutine

end module initialization_dft
