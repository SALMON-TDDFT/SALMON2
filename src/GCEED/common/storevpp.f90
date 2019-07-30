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
integer :: i,ak,L,l0,ll
real(8) :: r

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

  l0=0
  do ll=0,Mlps0(ak)
  do L=l0,l0+pp%nproj(ll,ak)-1
    do i=1,Mr(ak)
      r=rad_psl(i,ak)
      uppr(i,L,ak)=upp_f(i,L,ak)/r**(ll+1)*sqrt((2*ll+1)/(4*Pi))
    end do
    uppr(0,L,ak)=uppr(1,L,ak)
  end do
  l0=L
  end do

  l0=0
  do ll=0,Mlps0(ak)
  do L=l0,l0+pp%nproj(ll,ak)-1
    do i=0,Mr(ak)
      vpp(i,L,ak)=vpp_f(i,L,ak)
    end do
  end do
  l0=L
  end do
  if( Lref(ak) > Mlps0(ak) )then
    do i=0,Mr(ak)
      vpp(i,Lref(ak),ak)=vpp_f(i,Lref(ak),ak)
    end do
  end if

  uppr(Mr(ak)+1:Nr,:,ak)=0.d0

  l0=0
  do ll=0,Mlps0(ak)
  do L=l0,l0+pp%nproj(ll,ak)-1
    do i=Mr(ak)+1,Nr
      r=rad_psl(i,ak) ; vpp(i,L,ak)=-Zps(ak)/r
    end do
  end do
  l0=L
  end do
  if( Lref(ak) > Mlps0(ak) )then
    do i=Mr(ak)+1,Nr
       r=rad_psl(i,ak)
      vpp(i,Lref(ak),ak)=-Zps(ak)/r
    end do
  end if

end do

return

END SUBROUTINE storevpp
