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
module allocate_psl_sub

use scf_data
use read_pslfile_sub

integer :: Nr

real(8), allocatable :: rhopp(:,:)
real(8), allocatable :: uppr(:,:,:)

real(8), allocatable :: rad_psl(:,:)

real(8), allocatable :: ur(:,:)

real(8),allocatable :: rho_core(:,:,:)

real(8),allocatable :: vloctbl(:,:)
real(8),allocatable :: Gx(:),Gy(:),Gz(:)
real(8),allocatable :: dVloc_G(:,:)
real(8),allocatable :: dVloc_G_tmp(:,:)
complex(8),allocatable :: rhoion_G(:),Vion_G(:)
complex(8),allocatable :: rhoion_G_tmp(:)
complex(8),allocatable :: Vion_G_tmp(:)

integer :: nGzero

contains
!==================================================================================================

subroutine allocate_psl(lg)
use structures, only: s_rgrid
use scf_data
implicit none
type(s_rgrid),intent(in) :: lg
real(8) :: r
integer :: Nr0,nprj_u

select case(iperiodic)
case(0)
  r=sqrt((dble(lg%num(1))*Hgs(1))**2  &
        +(dble(lg%num(2))*Hgs(2))**2  &
        +(dble(lg%num(3))*Hgs(3))**2)
case(3)
  r=sqrt((dble(lg%num(1))*Hgs(1))**2  &
        +(dble(lg%num(2))*Hgs(2))**2  &
        +(dble(lg%num(3))*Hgs(3))**2)+maxval(Rps(:))
end select
Nr0=r/rmin_step+1
if(Nr0<maxval(Mr(:)))then
  Nr=maxval(Mr(:))+1
else
  Nr=Nr0
end if

maxMps=int(4.d0/3.d0*Pi*(rmaxRps+4.d0*maxval(Hgs(:)))**3/Hvol)
Mlmps=maxlm

nprj_u=size(upp_f,2)
allocate(rhopp(0:Nr,MKI))
allocate(uppr(0:Nr,0:nprj_u-1,MKI))
allocate(rad_psl(0:Nr,MKI))

allocate(ur(maxMps,Mlmps))

allocate(rho_core(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate(numatom_ps(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

if(iperiodic==3)then
  allocate(vloctbl(Nr,MKI))
  allocate(Gx(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(Gy(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(Gz(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(dVloc_G(lg%num(1)*lg%num(2)*lg%num(3),MKI))
  allocate(dVloc_G_tmp(lg%num(1)*lg%num(2)*lg%num(3),MKI))
  allocate(rhoion_G(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(rhoion_G_tmp(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(Vion_G(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(Vion_G_tmp(lg%num(1)*lg%num(2)*lg%num(3)))
end if

end subroutine allocate_psl
!==================================================================================================
subroutine deallocate_psl

deallocate(rhopp,uppr)
deallocate(ur)
deallocate(rad_psl)

deallocate(rho_core)

deallocate(numatom_ps)

if(iperiodic==3)then
  deallocate(vloctbl)
  deallocate(Gx,Gy,Gz)
  deallocate(dVloc_G,rhoion_G,Vion_G)
  deallocate(dVloc_G_tmp,rhoion_G_tmp,Vion_G_tmp)
end if

end subroutine deallocate_psl
!==================================================================================================
end module allocate_psl_sub
