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

subroutine main_dft2tddft
use structures
use salmon_global, only: ispin,cval,xc,xname,cname,directory_read_data
use salmon_parallel, only: nproc_group_global
use salmon_xc
use timer
use initialization_sub
use checkpoint_restart_sub
implicit none
integer :: Miter
character(100) :: file_atoms_coo

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
type(s_scalar) :: srho,sVh,sVpsl
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ofile)  :: ofl

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)

call convert_input_scf(file_atoms_coo)

! please move folloings into initialization_dft
call init_dft(nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofl)
allocate( srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin) )


call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          pinfo, info, info_field,  &
                          srg, srg_ng,  &
                          srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call read_bin(directory_read_data,lg,mg,ng,system,info,spsi,Miter)

call timer_end(LOG_INIT_GS)

! Redistributed write to use TDDFT calculation.
call timer_begin(LOG_WRITE_GS_DATA)
call write_bin(ofl%dir_out_restart,lg,mg,ng,system,info,spsi,Miter)
call timer_end(LOG_WRITE_GS_DATA)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)

end subroutine main_dft2tddft
