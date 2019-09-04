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
SUBROUTINE storevpp
!$ use omp_lib
use scf_data
use allocate_psl_sub
implicit none
integer :: i,ak

!-------------------- Read atomic pseudopotential data

do ak=1,MKI
  select case( ipsfileform(ak) )
    case(1)  
      do i=0,Nr
        rad_psl(i,ak)=dble(i)*step(ak)
      end do
    case(2)
      do i=0,Mr(ak)
        rad_psl(i,ak) = 1.0d2*(dble(i)/dble(Mr(ak))+1.0d-2)**5 - 1.0d-8
      end do
      do i=Mr(ak)+1,Nr
        rad_psl(i,ak) = 1.0d2*(dble(Mr(ak))/dble(Mr(ak))+1.0d-2)**5 - 1.0d-8 + step(ak)*dble(i-Mr(ak))
      end do
    case(3,4)
      do i=0,Mr(ak)
        rad_psl(i,ak) = rad_f(i,ak)
      end do
      do i=Mr(ak)+1,Nr
        rad_psl(i,ak) = rad_f(Mr(ak),ak)+ step(ak)*dble(i-Mr(ak))
      end do
    case default
       do i=0,size(rad_psl,1)-1
         rad_psl(i,ak)=i*step(ak)
       end do
  end select

end do

return

END SUBROUTINE storevpp
