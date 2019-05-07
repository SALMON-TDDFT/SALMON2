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
!=======================================================================

subroutine r_hpsi2_buf(tpsi,htpsi,iob,iiik,nn,isub)
  use hpsi2_sub
  implicit none
  real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,   &
                  mg_sta(2)-Nd:mg_end(2)+Nd,   &
                  mg_sta(3)-Nd:mg_end(3)+Nd)
  real(8) :: htpsi(mg_sta(1):mg_end(1),  &
                   mg_sta(2):mg_end(2),      &
                   mg_sta(3):mg_end(3))
  integer :: iob,iiik,nn,isub
  
  iwk_size=2

  call hpsi2(tpsi,htpsi,iob,iiik,nn,isub)


end subroutine r_hpsi2_buf

!=======================================================================

subroutine c_hpsi2_buf(tpsi,htpsi,iob,iiik,nn,isub)
  use hpsi2_sub
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,   &
                     mg_sta(2)-Nd:mg_end(2)+Nd,   &
                     mg_sta(3)-Nd:mg_end(3)+Nd)
  complex(8) :: htpsi(mg_sta(1):mg_end(1),  &
                      mg_sta(2):mg_end(2),      &
                      mg_sta(3):mg_end(3))
  integer :: iob,iiik,nn,isub
  
  iwk_size=2

  call hpsi2(tpsi,htpsi,iob,iiik,nn,isub)


end subroutine c_hpsi2_buf

!=======================================================================
