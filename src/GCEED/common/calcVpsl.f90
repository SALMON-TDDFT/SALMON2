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
SUBROUTINE calcVpsl(lg)
use structures, only: s_rgrid
use salmon_parallel, only: nproc_id_global
use prep_pp_sub, only: bisection
use scf_data
use allocate_psl_sub
implicit none

type(s_rgrid),intent(in)   :: lg
integer :: ix,iy,iz,ak
integer :: j,a,intr
real(8) :: ratio1,ratio2

real(8) :: r

Vpsl=0.d0 

allocate(ppg%Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI))
do a=1,MI
  ak=Kion(a)
  do j=1,3
    if(abs(Rion(j,a))<lg%num(j)*Hgs(j))then
      continue
    else
      write(*,*) "Rion error",nproc_id_global,a,j,Rion(j,a)
    end if
  end do
  do ix=mg_sta(1),mg_end(1)
  do iy=mg_sta(2),mg_end(2)
  do iz=mg_sta(3),mg_end(3)
    r=sqrt( (lg%coordinate(ix,1)-Rion(1,a))**2      &
           +(lg%coordinate(iy,2)-Rion(2,a))**2      &
           +(lg%coordinate(iz,3)-Rion(3,a))**2 )+1.d-50
    call bisection(r,intr,ak,pp%nrmax,pp%rad)
    ratio1=(r-pp%rad(intr,ak))/(pp%rad(intr+1,ak)-pp%rad(intr,ak)) ; ratio2=1.d0-ratio1
    if(intr>0.and.intr<=Nr)then
      continue
    else
      write(*,*) "intr error",nproc_id_global,intr,r
    end if

    Vpsl(ix,iy,iz)=Vpsl(ix,iy,iz)      &
                +ratio1*pp%vpp_f(intr,Lref(ak),ak)      &
                +ratio2*pp%vpp_f(intr-1,Lref(ak),ak)  !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate

    if(icalcforce==1)then
      Vpsl_atom(ix,iy,iz,a)=                    &
                +ratio1*pp%vpp_f(intr,Lref(ak),ak)      &
                +ratio2*pp%vpp_f(intr-1,Lref(ak),ak)
    end if
    ppg%Vpsl_atom(ix,iy,iz,a) = ratio1*pp%vpp_f(intr,Lref(ak),ak) + ratio2*pp%vpp_f(intr-1,Lref(ak),ak)
  end do
  end do
  end do
end do

return

END SUBROUTINE calcVpsl
