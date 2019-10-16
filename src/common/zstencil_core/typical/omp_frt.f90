!
!  Copyright 2017-2019 SALMON developers
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
!OCL eval_concurrent
subroutine zstencil_typical_omp(is_array,ie_array,is,ie,idx,idy,idz &
                               ,tpsi,htpsi,V_local,lap0,lapt,nabt &
                               )
  implicit none

  integer,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
  integer,intent(in) :: idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)

  complex(8),intent(in)  :: tpsi   (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  complex(8),intent(out) :: htpsi  (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8),   intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3))
  real(8),   intent(in)  :: lap0
  real(8),   intent(in)  :: lapt(12), nabt(12)

  complex(8), parameter :: zI=(0.d0,1.d0)

  integer :: iz,iy,ix
  integer :: endz,endy,endx
  complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,v,w

  endz = ie(3)
  endy = ie(2)
  endx = ie(1)

!$omp parallel do collapse(2) default(none) schedule(runtime) &
!$omp             private(iz,iy,ix,t1,t2,t3,t4,t5,t6,t7,t8,v,w) &
!$omp             firstprivate(lap0,lapt,nabt) &
!$omp             shared(idx,idy,idz,is,endz,endy,endx,htpsi,V_local,tpsi)
  do iz=is(3),endz
  do iy=is(2),endy

#define DX(dt) tpsi(idx(ix+(dt)),iy,iz)
#define DY(dt) tpsi(ix,idy(iy+(dt)),iz)
#define DZ(dt) tpsi(ix,iy,idz(iz+(dt)))

!OCL swp
  do ix=is(1),endx
    htpsi(ix,iy,iz) = (V_local(ix,iy,iz) + lap0) * tpsi(ix,iy,iz)
  end do

!OCL swp
  do ix=is(1),endx
    t8 = DX(-4) ; t7 = DX(-3) ; t6 = DX(-2) ; t5 = DX(-1)
    t1 = DX( 1) ; t2 = DX( 2) ; t3 = DX( 3) ; t4 = DX( 4)

    v=(lapt( 1)*(t1+t5) &
    & +lapt( 2)*(t2+t6) &
    & +lapt( 3)*(t3+t7) &
    & +lapt( 4)*(t4+t8))

    w=(nabt( 1)*(t1-t5) &
    & +nabt( 2)*(t2-t6) &
    & +nabt( 3)*(t3-t7) &
    & +nabt( 4)*(t4-t8))

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do

!OCL swp
  do ix=is(1),endx
    t8 = DY(-4) ; t7 = DY(-3) ; t6 = DY(-2) ; t5 = DY(-1)
    t1 = DY( 1) ; t2 = DY( 2) ; t3 = DY( 3) ; t4 = DY( 4)

    v=(lapt( 5)*(t1+t5) &
    & +lapt( 6)*(t2+t6) &
    & +lapt( 7)*(t3+t7) &
    & +lapt( 8)*(t4+t8))

    w=(nabt( 5)*(t1-t5) &
    & +nabt( 6)*(t2-t6) &
    & +nabt( 7)*(t3-t7) &
    & +nabt( 8)*(t4-t8))

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do

!OCL swp
  do ix=is(1),endx
    t8 = DZ(-4) ; t7 = DZ(-3) ; t6 = DZ(-2) ; t5 = DZ(-1)
    t1 = DZ( 1) ; t2 = DZ( 2) ; t3 = DZ( 3) ; t4 = DZ( 4)

    v=(lapt( 9)*(t1+t5) &
    & +lapt(10)*(t2+t6) &
    & +lapt(11)*(t3+t7) &
    & +lapt(12)*(t4+t8))

    w=(nabt( 9)*(t1-t5) &
    & +nabt(10)*(t2-t6) &
    & +nabt(11)*(t3-t7) &
    & +nabt(12)*(t4-t8))

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do

  end do
  end do
end subroutine
