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
module stencil_sub
complex(8), parameter :: zI=(0.d0,1.d0)

contains

!===================================================================================================================================

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine stencil_R(tpsi,htpsi,ib_array,ie_array &
                    ,V_local,ib,ie &
                    ,idx,idy,idz,lap0,lapt)
  implicit none
  integer,intent(in)  :: ib_array(3),ie_array(3),ib(3),ie(3) &
                        ,idx(ib(1)-4:ie(1)+4),idy(ib(2)-4:ie(2)+4),idz(ib(3)-4:ie(3)+4)
  real(8),intent(in)  :: tpsi(ib_array(1):ie_array(1),ib_array(2):ie_array(2),ib_array(3):ie_array(3)) &
                        ,V_local(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)),lap0,lapt(4,3)
  real(8),intent(out) :: htpsi(ib_array(1):ie_array(1),ib_array(2):ie_array(2),ib_array(3):ie_array(3))
  !
  integer :: iz,iy,ix
  real(8) :: v

!$OMP parallel
!$OMP do private(iz,iy,ix,v)
  do iz=ib(3),ie(3)
  do iy=ib(2),ie(2)
  do ix=ib(1),ie(1)

    v =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
        +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
        +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
        +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
        +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
        +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
        +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4))) + v

    v =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
        +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
        +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
        +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4))) + v

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_R

!===================================================================================================================================

subroutine stencil_C(tpsi,htpsi,ib_array,ie_array &
                    ,V_local,ib,ie &
                    ,idx,idy,idz,lap0,lapt,nabt)
  implicit none
  integer,intent(in)  :: ib_array(3),ie_array(3),ib(3),ie(3) &
                        ,idx(ib(1)-4:ie(1)+4),idy(ib(2)-4:ie(2)+4),idz(ib(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(ib_array(1):ie_array(1),ib_array(2):ie_array(2),ib_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)),lap0,lapt(4,3),nabt(4,3)
  complex(8),intent(out) :: htpsi(ib_array(1):ie_array(1),ib_array(2):ie_array(2),ib_array(3):ie_array(3))
  !
  integer :: iz,iy,ix
  complex(8) :: v,w

!$OMP parallel
!$OMP do private(iz,iy,ix,v,w)
  do iz=ib(3),ie(3)
  do iy=ib(2),ie(2)
  do ix=ib(1),ie(1)

    v =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
        +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
        +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
        +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
        +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
        +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
        +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4))) + v

    v =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
        +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
        +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
        +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4))) + v

    w =  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1))) &
        +nabt(2,1)*(tpsi(DX(2)) - tpsi(DX(-2))) &
        +nabt(3,1)*(tpsi(DX(3)) - tpsi(DX(-3))) &
        +nabt(4,1)*(tpsi(DX(4)) - tpsi(DX(-4)))

    w =  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1))) &
        +nabt(2,2)*(tpsi(DY(2)) - tpsi(DY(-2))) &
        +nabt(3,2)*(tpsi(DY(3)) - tpsi(DY(-3))) &
        +nabt(4,2)*(tpsi(DY(4)) - tpsi(DY(-4))) + w

    w =  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1))) &
        +nabt(2,3)*(tpsi(DZ(2)) - tpsi(DZ(-2))) &
        +nabt(3,3)*(tpsi(DZ(3)) - tpsi(DZ(-3))) &
        +nabt(4,3)*(tpsi(DZ(4)) - tpsi(DZ(-4))) + w

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v - zI * w
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_C

!===================================================================================================================================

end module stencil_sub
