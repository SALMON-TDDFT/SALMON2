!
!  Copyright 2020 SALMON developers
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
module network_address
  implicit none

contains
  ! convert 5-D address to rank
  function get_neighbour_rank(info, idir, idisp) result(irank)
    use salmon_global, only: yn_periodic
    use structures, only: s_parallel_info
    implicit none
    type(s_parallel_info), intent(in) :: info
    integer,               intent(in) :: idir, idisp
    integer :: iaddr(5),irank
    integer :: ishape(5)

    ishape(1:5) = [info%nprgrid(1:3), info%nporbital, info%npk]
    iaddr(1:5)  = info%iaddress(1:5)
    iaddr(idir) = iaddr(idir) + idisp

    if (yn_periodic == 'y') then
      ! periodic boundary
      if (iaddr(idir) <  0)            iaddr(idir) = iaddr(idir) + ishape(idir)
      if (iaddr(idir) >= ishape(idir)) iaddr(idir) = iaddr(idir) - ishape(idir)
    else
      ! dirichlet boundary
      if (iaddr(idir) <  0)            iaddr(idir) = -1
      if (iaddr(idir) >= ishape(idir)) iaddr(idir) = -1
    end if

    if (iaddr(idir) < 0) then
      irank = -1
    else
      irank = info%imap(iaddr(1),iaddr(2),iaddr(3),iaddr(4),iaddr(5))
    end if
  end function get_neighbour_rank

end module
