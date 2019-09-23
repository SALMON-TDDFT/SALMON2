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
MODULE new_world_sub

use scf_data

implicit none

integer,allocatable :: icorr_polenum(:)

CONTAINS

!======================================================================

subroutine wrapper_allgatherv_vlocal(ng,mg,info_field) ! --> remove (future works)
  use structures
  use hamiltonian
  implicit none
  type(s_rgrid),         intent(in) :: ng,mg
  type(s_field_parallel),intent(in) :: info_field
  type(s_scalar) :: sVh,sVpsl
  type(s_scalar),allocatable :: sVlocal(:),sVxc(:)
  integer :: nspin,jspin

  nspin = 1
  if(ispin==1) nspin = 2
  allocate(sVlocal(nspin),sVxc(nspin))
  do jspin=1,nspin
    allocate(sVlocal(jspin)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    allocate(sVxc(jspin)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  end do
  allocate(sVh%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  allocate(sVpsl%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  sVpsl%f = Vpsl
  sVh%f = Vh
  if(ilsda == 1) then
    do jspin=1,nspin
      sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
    end do
  else
    sVxc(1)%f = Vxc
  end if

  call allgatherv_vlocal(ng,mg,info_field,nspin,sVh,sVpsl,sVxc,sVlocal)

  do jspin=1,nspin
    Vlocal(:,:,:,jspin) = sVlocal(jspin)%f
    call deallocate_scalar(sVxc(jspin))
    call deallocate_scalar(sVlocal(jspin))
  end do
  call deallocate_scalar(sVh)
  call deallocate_scalar(sVpsl)
  return
end subroutine wrapper_allgatherv_vlocal

END MODULE new_world_sub

