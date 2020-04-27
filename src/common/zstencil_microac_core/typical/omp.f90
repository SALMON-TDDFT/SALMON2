!
!  Copyright 2020-2020 SALMON developers
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

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine zstencil_microAc_typical_omp(is_array,ie_array,is,ie,idx,idy,idz &
                                       ,tpsi,htpsi,V_local,Ac,div_Ac,lap0,lapt,nabt,k &
                                       )
  use math_constants,only : zi
  implicit none

  integer,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
  integer,intent(in) :: idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)

  complex(8),intent(in)  :: tpsi   (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  complex(8),intent(out) :: htpsi  (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8),   intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3))
  real(8),   intent(in)  :: Ac(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
  real(8),   intent(in)  :: div_Ac(is(1):ie(1),is(2):ie(2),is(3):ie(3))
  real(8),   intent(in)  :: lap0
  real(8),   intent(in)  :: lapt(4,3), nabt(4,3)
  real(8),   intent(in)  :: k(3)

  integer    :: ix,iy,iz
  real(8)    :: kAc(3),div
  complex(8) :: w(3),v,psi0

!$OMP parallel
!$OMP do collapse(2) private(iz,iy,ix,w,v,psi0,kAc,div)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    psi0 = tpsi(ix,iy,iz)
    kAc = k + Ac(:,ix,iy,iz)
    div = div_Ac(ix,iy,iz)

  ! laplacian of psi
    v =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
      & +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
      & +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
      & +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
      & +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
      & +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
      & +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4))) + v

    v =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
      & +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
      & +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
      & +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4))) + v

  ! gradient of psi
    w(1) =  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1))) &
         & +nabt(2,1)*(tpsi(DX(2)) - tpsi(DX(-2))) &
         & +nabt(3,1)*(tpsi(DX(3)) - tpsi(DX(-3))) &
         & +nabt(4,1)*(tpsi(DX(4)) - tpsi(DX(-4)))

    w(2) =  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1))) &
         & +nabt(2,2)*(tpsi(DY(2)) - tpsi(DY(-2))) &
         & +nabt(3,2)*(tpsi(DY(3)) - tpsi(DY(-3))) &
         & +nabt(4,2)*(tpsi(DY(4)) - tpsi(DY(-4)))

    w(3) =  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1))) &
         & +nabt(2,3)*(tpsi(DZ(2)) - tpsi(DZ(-2))) &
         & +nabt(3,3)*(tpsi(DZ(3)) - tpsi(DZ(-3))) &
         & +nabt(4,3)*(tpsi(DZ(4)) - tpsi(DZ(-4)))

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )* psi0 - 0.5d0* v           &
                    & + 0.5d0* ( kAc(1)**2 + kAc(2)**2 + kAc(3)**2 ) * psi0   &
                    & - zi* ( kAc(1) * w(1) + kAc(2) * w(2) + kAc(3) * w(3) ) &
                    & - zi*0.5d0* div * psi0 ! deviation from the Coulomb gauge condition (for numerical stability)
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel
end subroutine
