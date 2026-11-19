!
!  Copyright 2019-2020 SALMON developers
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

#ifdef __FUJITSU
#else
#define LOOP_BLOCKING
#endif

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

  integer    :: ix,iy,iz
  complex(8) :: v,w
  complex(8) :: t(8)

#ifdef LOOP_BLOCKING
  integer :: bz,bz_e,bz_step
  integer :: by,by_e,by_step
#endif

#ifdef __INTEL_COMPILER
#if defined(__KNC__) || defined(__AVX512F__)
#   define MEM_ALIGN   64
#   define VECTOR_SIZE 4
# else
#   define MEM_ALIGN   32
#   define VECTOR_SIZE 2
# endif

!dir$ assume_aligned V_local:MEM_ALIGN
!dir$ assume_aligned tpsi   :MEM_ALIGN
!dir$ assume_aligned htpsi  :MEM_ALIGN
#endif

#define DX(dt) idx(ix+(dt)),iy,iz
#define DY(dt) ix,idy(iy+(dt)),iz
#define DZ(dt) ix,iy,idz(iz+(dt))

#ifdef LOOP_BLOCKING
  bz_step = 8
  by_step = 8

!$omp parallel do collapse(2) private(ix,iy,iz,v,w,t,bz,by,bz_e,by_e)
  do bz=is(3),ie(3),bz_step
  do by=is(2),ie(2),by_step

    bz_e = min(bz + bz_step - 1, ie(3))
    by_e = min(by + by_step - 1, ie(2))

  do iz=bz,bz_e
  do iy=by,by_e
#else
!$omp parallel do collapse(2) private(ix,iy,iz,v,w,t)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
#endif

!dir$ assume_aligned V_local(is(1),iy,iz)     :MEM_ALIGN
!dir$ assume_aligned tpsi(is_array(1),iy,iz)  :MEM_ALIGN
!dir$ assume_aligned htpsi(is_array(1),iy,iz) :MEM_ALIGN

  do ix=is(1),ie(1)
    t(1) = tpsi(DX( 4))
    t(2) = tpsi(DX( 3))
    t(3) = tpsi(DX( 2))
    t(4) = tpsi(DX( 1))
    t(5) = tpsi(DX(-1))
    t(6) = tpsi(DX(-2))
    t(7) = tpsi(DX(-3))
    t(8) = tpsi(DX(-4))

    v=(lapt(1)*(t(4)+t(5)) &
    & +lapt(2)*(t(3)+t(6)) &
    & +lapt(3)*(t(2)+t(7)) &
    & +lapt(4)*(t(1)+t(8)))
    w=(nabt(1)*(t(4)-t(5)) &
    & +nabt(2)*(t(3)-t(6)) &
    & +nabt(3)*(t(2)-t(7)) &
    & +nabt(4)*(t(1)-t(8)))

    t(1) = tpsi(DY( 1))
    t(2) = tpsi(DY( 2))
    t(3) = tpsi(DY( 3))
    t(4) = tpsi(DY( 4))
    t(5) = tpsi(DY(-1))
    t(6) = tpsi(DY(-2))
    t(7) = tpsi(DY(-3))
    t(8) = tpsi(DY(-4))

    v=(lapt(5)*(t(1)+t(5)) &
    & +lapt(6)*(t(2)+t(6)) &
    & +lapt(7)*(t(3)+t(7)) &
    & +lapt(8)*(t(4)+t(8))) + v
    w=(nabt(5)*(t(1)-t(5)) &
    & +nabt(6)*(t(2)-t(6)) &
    & +nabt(7)*(t(3)-t(7)) &
    & +nabt(8)*(t(4)-t(8))) + w

    htpsi(ix,iy,iz) = V_local(ix,iy,iz)*tpsi(ix,iy,iz) &
                    + lap0*tpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do

#if _OPENMP >= 201307
#ifndef __ARM_FLANG
!$omp simd
#endif
#endif
  do ix=is(1),ie(1)
    t(1) = tpsi(DZ( 1))
    t(2) = tpsi(DZ( 2))
    t(3) = tpsi(DZ( 3))
    t(4) = tpsi(DZ( 4))
    t(5) = tpsi(DZ(-1))
    t(6) = tpsi(DZ(-2))
    t(7) = tpsi(DZ(-3))
    t(8) = tpsi(DZ(-4))

    v=(lapt( 9)*(t(1)+t(5)) &
    & +lapt(10)*(t(2)+t(6)) &
    & +lapt(11)*(t(3)+t(7)) &
    & +lapt(12)*(t(4)+t(8)))
    w=(nabt( 9)*(t(1)-t(5)) &
    & +nabt(10)*(t(2)-t(6)) &
    & +nabt(11)*(t(3)-t(7)) &
    & +nabt(12)*(t(4)-t(8)))

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do
  end do
  end do

#ifdef LOOP_BLOCKING
  end do
  end do
#endif
end subroutine
