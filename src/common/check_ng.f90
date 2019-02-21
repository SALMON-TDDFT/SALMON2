!
!  Copyright 2017 SALMON developers
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
subroutine check_ng(ng)
  use inputoutput, only: iperiodic
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid)     :: ng
  integer,parameter :: nd=4

  if(iperiodic==0.and.(ng%num(1)<nd.or.ng%num(2)<nd.or.ng%num(3)<nd))then
    stop "A system is small. Please use less number of processors."
  end if

end subroutine check_ng
