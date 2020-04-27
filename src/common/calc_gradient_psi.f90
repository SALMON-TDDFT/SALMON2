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

#include "config.h"

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine calc_gradient_psi(tpsi,gtpsi,is_array,ie_array,is,ie,idx,idy,idz,nabt,matrix_B)
  implicit none
  integer,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                        ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: nabt(4,3),matrix_B(3,3)
  complex(8),intent(out) :: gtpsi(3,is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  !
  integer :: iz,iy,ix
  complex(8) :: w(3)
  real(8) :: Bt(3,3)

  Bt = transpose(matrix_B)

!$OMP parallel do collapse(2) private(iz,iy,ix,w) schedule(runtime)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)

#ifdef __FUJITSU
!OCL swp
  do ix=is(1),ie(1)
    gtpsi(1,ix,iy,iz) = &
         &  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1))) &
         & +nabt(2,1)*(tpsi(DX(2)) - tpsi(DX(-2))) &
         & +nabt(3,1)*(tpsi(DX(3)) - tpsi(DX(-3))) &
         & +nabt(4,1)*(tpsi(DX(4)) - tpsi(DX(-4)))
  end do

!OCL swp
  do ix=is(1),ie(1)
    gtpsi(2,ix,iy,iz) = &
         &  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1))) &
         & +nabt(2,2)*(tpsi(DY(2)) - tpsi(DY(-2))) &
         & +nabt(3,2)*(tpsi(DY(3)) - tpsi(DY(-3))) &
         & +nabt(4,2)*(tpsi(DY(4)) - tpsi(DY(-4)))
  end do

!OCL swp
  do ix=is(1),ie(1)
    gtpsi(3,ix,iy,iz) = &
         &  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1))) &
         & +nabt(2,3)*(tpsi(DZ(2)) - tpsi(DZ(-2))) &
         & +nabt(3,3)*(tpsi(DZ(3)) - tpsi(DZ(-3))) &
         & +nabt(4,3)*(tpsi(DZ(4)) - tpsi(DZ(-4)))
  end do

!OCL swp
  do ix=is(1),ie(1)
    w(:) = gtpsi(:,ix,iy,iz)
    gtpsi(1,ix,iy,iz) = dot_product(matrix_B(:,1),w) ! B^{T} * (nabla) psi
    gtpsi(2,ix,iy,iz) = dot_product(matrix_B(:,2),w) ! B^{T} * (nabla) psi
    gtpsi(3,ix,iy,iz) = dot_product(matrix_B(:,3),w) ! B^{T} * (nabla) psi
  end do

#else
  do ix=is(1),ie(1)
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

    gtpsi(:,ix,iy,iz) = matmul(Bt,w) ! B^{T} * (nabla) psi
  end do
#endif

  end do
  end do
!$OMP end parallel do

  return
end subroutine calc_gradient_psi


