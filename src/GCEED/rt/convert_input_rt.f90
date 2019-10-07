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
subroutine convert_input_rt(Ntime)
use salmon_parallel, only: nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_bcast
use inputoutput
use set_numcpu
use scf_data
implicit none
integer :: Ntime
real(8) :: dip_spacing

ilsda=ispin

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if

if(sum(abs(num_rgrid)) /= 0 .and. sum(abs(dl)) /= 0d0)then
  call err_finalize('Error: [num_rgrid] and [dl] are incompatible input parameters.')
else if(sum(abs(num_rgrid)) /= 0 .and. sum(abs(dl)) == 0d0)then
  Harray(1:3,1)=al(1:3)/dble(num_rgrid(1:3))
else if(sum(abs(num_rgrid)) == 0d0 .and. sum(abs(dl)) /= 0d0)then
  Harray(1:3,1)=dl(1:3)
end if

rLsize(1:3,1)=al(1:3)
iDiter(1) = nscf

if(ispin == 0)then
  MST(1)=nstate
  if(temperature_k>=0.d0)then
    ifMST(1)=nstate
    rNetot=dble(nelec)
  else
    ifMST(1)=nelec/2
  end if
  MST(2)=0
  ifMST(2)=0
else if(ispin == 1)then
  if(nstate /= 0 .and. sum(nstate_spin) ==0)then
    MST(1:2)=nstate
    if(temperature_k>=0.d0)then
      ifMST(1:2)=nstate
      rNetot=dble(nelec)
    else
      ifMST(1)=nelec - nelec/2
      ifMST(2)=nelec/2
    end if
  else if(nstate == 0 .and. sum(nstate_spin) /=0)then
    MST(1)=maxval(nstate_spin(1:2))
    MST(2)=maxval(nstate_spin(1:2))
    if(temperature_k>=0.d0)then
      ifMST(1)=maxval(nstate_spin(1:2))
      ifMST(2)=maxval(nstate_spin(1:2))
      rNetot=dble(nelec_spin(1))+dble(nelec_spin(2))
    else
      ifMST(1:2)=nelec_spin(1:2)
    end if
  else
    write(*,*)"'nstate' or 'nstate_spin' should be spacified in input. "
  end if
else
  write(*,*)"'ispin' should be 0 or 1. "
end if

num_kpoints_3d(1:3)=num_kgrid(1:3)
num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)
allocate(wtk(num_kpoints_rd))
wtk(:)=1.d0/dble(num_kpoints_rd)

if(ilsda==1)then
  nproc_ob_spin(1)=(nproc_ob+1)/2
  nproc_ob_spin(2)=nproc_ob/2
end if

if(ilsda == 0) then
  itotMST=MST(1)
  itotfMST=ifMST(1)
else if(ilsda == 1) then
  itotMST=MST(1)+MST(2)
  itotfMST=ifMST(1)+ifMST(2)
end if

!===== namelist for group_fundamental =====
if(iwrite_projection==1.and.itwproj==-1)then
  write(*,*) "Please specify itwproj when iwrite_projection=1."
  stop
end if

nproc_d_o = nproc_domain_orbital
nproc_d_g = nproc_domain_general

if(nproc_ob==0.and.nproc_d_o(1)==0.and.nproc_d_o(2)==0.and.nproc_d_o(3)==0.and.  &
                   nproc_d_g(1)==0.and.nproc_d_g(2)==0.and.nproc_d_g(3)==0) then
  call set_numcpu_rt(nproc_d_o,nproc_d_g,nproc_d_g_dm)
else
  call check_numcpu(nproc_d_o,nproc_d_g,nproc_d_g_dm)
end if

nproc_d_o_mul=nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
nproc_d_g_mul_dm=nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

!===== namelist for group_propagation =====
Ntime=nt
if(dt<=1.d-10)then
  write(*,*) "please set dt."
  stop
end if
if(Ntime==0)then
  write(*,*) "please set nt."
  stop
end if

!===== namelist for group_hartree =====
if(layout_multipole<=0.or.layout_multipole>=4)then
  stop "layout_multipole must be equal to 1 or 2 or 3."
else if(layout_multipole==3)then
  if(num_multipole_xyz(1)==0.and.num_multipole_xyz(2)==0.and.num_multipole_xyz(3)==0)then
    continue
  else if(num_multipole_xyz(1)<=0.or.num_multipole_xyz(2)<=0.or.num_multipole_xyz(3)<=0)then
    stop "num_multipole_xyz must be largar than 0 when they are not default values."
  end if
  if(num_multipole_xyz(1)==0.and.num_multipole_xyz(2)==0.and.num_multipole_xyz(3)==0)then
    dip_spacing = 8.d0/au_length_aa  ! approximate spacing of multipoles 
    num_multipole_xyz(:)=int((al(:)+dip_spacing)/dip_spacing-1.d-8)
  end if
end if

!===== namelist for group_file =====
if(ic==0)then
  ic=1
end if

if(IC==3.and.num_datafiles_IN/=nproc_size_global)then
  if(comm_is_root(nproc_id_global))then
    write(*,*) "num_datafiles_IN is set to nproc."
  end if
  num_datafiles_IN=nproc_size_global
end if

!===== namelist for group_extfield =====
if(ae_shape1 == 'impulse')then
  ikind_eext = 0
else
  ikind_eext = 1
end if

Fst = e_impulse
romega = omega1
pulse_T = tw1
rlaser_I = I_wcm2_1

!===== namelist for group_others =====

if(comm_is_root(nproc_id_global))then
  if(iwdenoption/=0.and.iwdenoption/=1)then
    write(*,*)  'iwdenoption must be equal to 0 or 1.'
    stop
  end if
end if
if(iwdenoption==0)then
  iwdenstep=0
end if

if(comm_is_root(nproc_id_global))close(fh_namelist)

end subroutine convert_input_rt
