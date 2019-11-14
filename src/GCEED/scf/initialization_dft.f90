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

subroutine initialization_dft( system, energy, stencil, fg, poisson,  &
                               lg, mg, ng,  &
                               pinfo, info, info_field,  &
                               srg, srg_ng,  &
                               srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                               spsi, shpsi, sttpsi,  &
                               pp, ppg,  &
                               ofile,  &
                               nspin, flag_opt_conv, nopt_max )
use math_constants, only: pi, zi
use structures
use salmon_parallel, only: nproc_id_global,nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
!use density_matrix, only: calc_density
use writefield
use global_variables_scf
!use salmon_pp, only: calc_nlcc
!use hartree_sub, only: hartree
!use force_sub
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use input_pp_sub
use prep_pp_sub
!use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use salmon_total_energy
!use init_gs, only: init_wf
!use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
implicit none
integer :: jspin,nspin
!integer :: ix,iy,iz,ik
!integer :: iter,iatom,iob,p1,p2,p5,jj,iflag
!real(8) :: sum0,sum1
character(100) :: file_atoms_coo, comment_line
!real(8) :: rNebox1,rNebox2
!integer :: itmg

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
!type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
!type(s_cg)     :: cg
!type(s_mixing) :: mixing
type(s_ofile)  :: ofile

integer :: nopt_max
logical :: flag_opt_conv

!real(8),allocatable :: esp_old(:,:,:)
!integer :: i,j


!call init_dft(iSCFRT,nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)

call init_code_optimization
call allocate_mat(ng,mg,lg) ! future work: remove this line

allocate( energy%esp(system%no,system%nk,system%nspin) ); energy%esp=0.0d0
!allocate( esp_old(system%no,system%nk,system%nspin) ); esp_old=0.0d0

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
   call structure_opt_ini(natom)  !check later MI->natom
!   call structure_opt_ini(MI)
   flag_opt_conv=.false.
   write(comment_line,10) 0
   call write_xyz(comment_line,"new","r  ",system)
10 format("#opt iteration step=",i5)
end if

if(yn_opt=='y') then ; nopt_max = nopt
else                 ; nopt_max = 1
endif

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

end subroutine initialization_dft
