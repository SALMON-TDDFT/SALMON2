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
use parallelization, only: nproc_group_global, nproc_id_global
use communication, only: comm_is_root, comm_bcast
use set_numcpu
use inputoutput
implicit none
integer :: ii  !,iatom
integer :: ibox2
integer :: icheck1,icheck2
character(100) :: file_atoms_coo

!===== namelist for group_atom =====
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
!if(comm_is_root(nproc_id_global))then
!  read(fh_namelist,NML=group_atom) 
!  rewind(fh_namelist)

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
!end if

!if(comm_is_root(nproc_id_global)) write(*,*) "MI =",MI


!===== namelist for group_others =====

!if(comm_is_root(nproc_id_global))close(fh_namelist)

return

end subroutine convert_input_scf
