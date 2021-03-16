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
!OCL eval_concurrent
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
  integer    :: endx,endy,endz
  real(8)    :: kAc1,kAc2,kAc3,div
  complex(8) :: v,w,psi0
  complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8

  endz = ie(3)
  endy = ie(2)
  endx = ie(1)

write(*,'(a, a, a, i0)') "OMP DEBUG STRING" , __FILE__ , ": ",  __LINE__
!$omp parallel do collapse(2) default(none) schedule(runtime) &
!$omp             private(iz,iy,ix,t1,t2,t3,t4,t5,t6,t7,t8,v,w,psi0,kAc1,kAc2,kAc3,div) &
!$omp             firstprivate(lap0,lapt,nabt,k) &
!$omp             shared(idx,idy,idz,is,endz,endy,endx,htpsi,V_local,Ac,div_Ac,tpsi)
  do iz=is(3),endz
  do iy=is(2),endy

#define DX(dt) tpsi(idx(ix+(dt)),iy,iz)
#define DY(dt) tpsi(ix,idy(iy+(dt)),iz)
#define DZ(dt) tpsi(ix,iy,idz(iz+(dt)))

!OCL swp
  do ix=is(1),endx
    psi0   = tpsi(ix,iy,iz)
    kAc1   = Ac(1,ix,iy,iz) + k(1)
    kAc2   = Ac(2,ix,iy,iz) + k(2)
    kAc3   = Ac(3,ix,iy,iz) + k(3)
    div    = div_Ac(ix,iy,iz)
    htpsi(ix,iy,iz) = (V_local(ix,iy,iz) + lap0) * psi0 &
                    + 0.5d0 * (kAc1*kAc1 + kAc2*kAc2 + kAc3*kAc3) * psi0 &
                    - zi * 0.5d0 * div * psi0 ! deviation from the Coulomb gauge condition (for numerical stability)
  end do

!OCL swp
  do ix=is(1),endx
    kAc1 = Ac(1,ix,iy,iz) + k(1)

    t8 = DX(-4) ; t7 = DX(-3) ; t6 = DX(-2) ; t5 = DX(-1)
    t1 = DX( 1) ; t2 = DX( 2) ; t3 = DX( 3) ; t4 = DX( 4)

  ! laplacian of psi
    v =  lapt(1,1)*(t1+t5) &
      & +lapt(2,1)*(t2+t6) &
      & +lapt(3,1)*(t3+t7) &
      & +lapt(4,1)*(t4+t8)

  ! gradient of psi
    w =  nabt(1,1)*(t1-t5) &
      & +nabt(2,1)*(t2-t6) &
      & +nabt(3,1)*(t3-t7) &
      & +nabt(4,1)*(t4-t8)

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v &
                    - zi * (kAc1 * w)
  end do

!OCL swp
!OCL swp_freg_rate(105)
  do ix=is(1),endx
    kAc2 = Ac(2,ix,iy,iz) + k(2)

    t8 = DY(-4) ; t7 = DY(-3) ; t6 = DY(-2) ; t5 = DY(-1)
    t1 = DY( 1) ; t2 = DY( 2) ; t3 = DY( 3) ; t4 = DY( 4)

  ! laplacian of psi
    v =  lapt(1,2)*(t1+t5) &
      & +lapt(2,2)*(t2+t6) &
      & +lapt(3,2)*(t3+t7) &
      & +lapt(4,2)*(t4+t8)

  ! gradient of psi
    w =  nabt(1,2)*(t1-t5) &
      & +nabt(2,2)*(t2-t6) &
      & +nabt(3,2)*(t3-t7) &
      & +nabt(4,2)*(t4-t8)

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v &
                    - zi * (kAc2 * w)
  end do

!OCL swp
!OCL swp_freg_rate(105)
  do ix=is(1),endx
    kAc3 = Ac(3,ix,iy,iz) + k(3)

    t8 = DZ(-4) ; t7 = DZ(-3) ; t6 = DZ(-2) ; t5 = DZ(-1)
    t1 = DZ( 1) ; t2 = DZ( 2) ; t3 = DZ( 3) ; t4 = DZ( 4)

  ! laplacian of psi
    v =  lapt(1,3)*(t1+t5) &
      & +lapt(2,3)*(t2+t6) &
      & +lapt(3,3)*(t3+t7) &
      & +lapt(4,3)*(t4+t8)

  ! gradient of psi
    w =  nabt(1,3)*(t1-t5) &
      & +nabt(2,3)*(t2-t6) &
      & +nabt(3,3)*(t3-t7) &
      & +nabt(4,3)*(t4-t8)

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v &
                    - zi * (kAc3 * w)
  end do

  end do
  end do
!$OMP end parallel do
end subroutine
