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
!-----------------------------------------------------------------------------------------
subroutine set_coo_em(iperi,Nd,ioe,ista,iend,hgs,coo)
  implicit none
  integer,intent(in)  :: iperi,Nd
  integer,intent(in)  :: ioe(3),ista(3),iend(3)
  real(8),intent(in)  :: hgs(3)
  real(8),intent(out) :: coo(minval(ista(:))-Nd:maxval(iend(:))+Nd,3)
  integer :: ii,ij
  
  do ii=1,3
    select case(iperi)
    case(0)
      select case(ioe(ii))
      case(1)
        do ij=ista(ii)-Nd,iend(ii)+Nd
          coo(ij,ii)=dble(ij)*hgs(ii)
        end do
      case(2)
        do ij=ista(ii)-Nd,iend(ii)+Nd
          coo(ij,ii)=(dble(ij)-0.5d0)*hgs(ii)
        end do
      end select
    case(3)
      do ij=ista(ii)-Nd,iend(ii)+Nd
        coo(ij,ii)=dble(ij-1)*hgs(ii)
      end do
    end select
  end do
  
end subroutine set_coo_em
