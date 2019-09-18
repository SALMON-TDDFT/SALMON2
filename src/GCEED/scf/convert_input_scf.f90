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
subroutine convert_input_scf(file_atoms_coo)
use salmon_global
use salmon_parallel, only: nproc_group_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_bcast
use set_numcpu
use inputoutput
use scf_data
use new_world_sub
implicit none
integer :: ii  !,iatom
integer :: ibox2
integer :: icheck1,icheck2
character(100) :: file_atoms_coo
real(8) :: dip_spacing

ilsda = ispin

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if
!===== namelist for group_fundamental =====

select case(method_min)
  case('cg','diis','cg-diis')
    continue
  case default
    stop 'Specify either "cg", "diis", or "cg-diis" for method_min.'
end select

select case(method_mixing)
  case ('simple','broyden')
    continue
  case default
    stop 'Specify either "simple" or "broyden" for method_mixing.'
end select
    
Harray(1:3,1:maxntmg)=0.d0
rLsize(1:3,1:maxntmg)=0.d0
iDiter(1:maxntmg)=1000

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

select case(yn_subspace_diagonalization)
case('y')
  iflag_subspace_diag = 1
case('n')
  iflag_subspace_diag = 0
end select

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

allocate( rocc(itotMST,num_kpoints_rd) )

rocc(:,:)=0.d0             
if(ilsda == 0) then
  rocc(:ifMST(1),:num_kpoints_rd) = 2.d0   ! Occupation number
else if(ilsda == 1) then
  rocc(:ifMST(1),:num_kpoints_rd) = 1.d0   ! Occupation number
  rocc(MST(1)+1:MST(1)+ifMST(2),:num_kpoints_rd) = 1.d0   ! Occupation number
end if

!===== namelist for group_parallel =====
if(comm_is_root(nproc_id_global))then
  ibox2=1
  icheck1=0
  icheck2=0
  do ii=1,19
    if(num_datafiles_IN==ibox2) icheck1=1
    if(num_datafiles_OUT==ibox2) icheck2=1
    ibox2=ibox2*2 
  end do
  if(icheck1/=1.or.icheck2/=1)then
    write(*,*) "num_datafiles_IN and num_datafiles_OUT must be equal to nth power of 2. (n: positive integer)"
    stop
  end if
end if

nproc_d_o = nproc_domain_orbital
nproc_d_g = nproc_domain_general

if(nproc_ob==0.and.nproc_d_o(1)==0.and.nproc_d_o(2)==0.and.nproc_d_o(3)==0.and.  &
                   nproc_d_g(1)==0.and.nproc_d_g(2)==0.and.nproc_d_g(3)==0) then
  call set_numcpu_gs(nproc_d_o,nproc_d_g,nproc_d_g_dm)
else
  call check_numcpu(nproc_d_o,nproc_d_g,nproc_d_g_dm)
end if

nproc_d_o_mul=nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
nproc_d_g_mul_dm=nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

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
if(comm_is_root(nproc_id_global))then
  if(IC<0.or.IC>=2)then
    write(*,*) "IC must be equal to 0 or 1."
    stop
  end if
end if

!===== namelist for group_atom =====
iflag_ps=1
MI=natom
MKI=nelem
!iZatom(:)=0
!ipsfileform(:)=1
!ps_format = 'default'
if(file_atom_coor/='none')then
   file_atoms_coo=trim(file_atom_coor)
else
   file_atoms_coo='.atomic_coor.tmp'
end if

!Lmax_ps(:)=-1
!Lloc_ps(:)=-1
if(comm_is_root(nproc_id_global))then
!  read(fh_namelist,NML=group_atom) 
!  rewind(fh_namelist)

  Mlps(:) = Lmax_ps(:)
  Lref(:) = Lloc_ps(:)

!ps format conversion
!  do iatom = 1,MKI
!    select case(ps_format(iatom))
!    case('default')
!    case('KY')        ; ipsfileform(iatom)=n_Yabana_Bertsch_psformat
!    case('ABINIT')    ; ipsfileform(iatom)=n_ABINIT_psformat
!    case('FHI')       ; ipsfileform(iatom)=n_FHI_psformat
!    case('ABINITFHI') !; ipsfileform(iatom)=n_ABINITFHI_psformat
!      write(*,"(A)") "Invalid ps_format. ABINITFHI format is not supported for isolated systems."
!    stop
!    case default
!      write(*,"(A)") "Invalid ps_format."
!    stop
!    end select
!  end do
end if

call comm_bcast(Mlps,nproc_group_global)
call comm_bcast(Lref,nproc_group_global)

if(comm_is_root(nproc_id_global)) write(*,*) "MI =",MI


!===== namelist for group_others =====

!if(iflag_ps==1)then
!  do ii=1,3
!    do iatom=1,MI
!      Rion(ii,iatom)=Rion(ii,iatom)*ulength_to_au
!    end do
!  end do
!end if

nproc_d_o_mul   = nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
nproc_d_g_mul_dm= nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

if(comm_is_root(nproc_id_global))close(fh_namelist)

return

end subroutine convert_input_scf
