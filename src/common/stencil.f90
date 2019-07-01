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
module stencil_sub
use math_constants,only : zi

contains

!===================================================================================================================================

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine stencil_R(is_array,ie_array,is,ie,idx,idy,idz &
                    ,tpsi,htpsi,V_local,lap0,lapt)
  implicit none
  integer,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                        ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  real(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3)) &
                        ,V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)),lap0,lapt(4,3)
  real(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  !
  integer :: iz,iy,ix
  real(8) :: v

!$OMP parallel
!$OMP do private(iz,iy,ix,v)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)

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

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_R

!===================================================================================================================================

subroutine stencil_C(is_array,ie_array,is,ie,idx,idy,idz &
                    ,tpsi,htpsi,V_local,lap0,lapt,nabt)
#ifdef SALMON_EXPLICIT_VECTORIZATION
  use code_optimization, only: modx,mody,modz,stencil_is_parallelized_by_omp
#endif
  implicit none
  integer,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                        ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)),lap0,lapt(4,3),nabt(4,3)
  complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))

#ifdef SALMON_EXPLICIT_VECTORIZATION
#define idx modx
#define idy mody
#define idz modz
#endif

  if (stencil_is_parallelized_by_omp) then
    call stencil_C_omp(is_array,ie_array,is,ie,idx,idy,idz,tpsi,htpsi,V_local,lap0,lapt,nabt)
  else
    call stencil_C_seq(is_array,ie_array,is,ie,idx,idy,idz,tpsi,htpsi,V_local,lap0,lapt,nabt)
  end if

#ifdef SALMON_EXPLICIT_VECTORIZATION
#undef idx modx
#undef idy mody
#undef idz modz
#endif

  return
end subroutine stencil_C

!===================================================================================================================================

subroutine stencil_nonorthogonal(is_array,ie_array,is,ie,idx,idy,idz,wrk &
                                ,tpsi,htpsi,V_local,lap0,lapt,nabt,Bk,F)
  implicit none
  integer   ,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                           ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)),lap0,lapt(4,3),nabt(4,3)
  real(8)    ,intent(in) :: Bk(3),F(6)
  complex(8)             :: wrk(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),2)
  complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  !
  integer :: ix,iy,iz
  complex(8) :: w(3),v(3)

!$OMP parallel
!$OMP do private(iz,iy,ix,w,v)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)

    v(1) =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
         & +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
         & +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
         & +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v(2) =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
         & +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
         & +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
         & +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4)))

    v(3) =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
         & +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
         & +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
         & +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4)))

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

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )* tpsi(ix,iy,iz) &
                    - 0.5d0* ( F(1)*v(1) + F(2)*v(2) + F(3)*v(3) ) - zI* ( Bk(1)*w(1) + Bk(2)*w(2) + Bk(3)*w(3) )
    wrk(ix,iy,iz,1) = w(1) ! df/dx
    wrk(ix,iy,iz,2) = w(2) ! df/dy
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

!$OMP parallel
!$OMP do private(iz,iy,ix,w)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)

  ! yz: (d/dz) * (df/dy)
    w(1) =  nabt(1,3)*(wrk(DZ(1),2) - wrk(DZ(-1),2)) &
         & +nabt(2,3)*(wrk(DZ(2),2) - wrk(DZ(-2),2)) &
         & +nabt(3,3)*(wrk(DZ(3),2) - wrk(DZ(-3),2)) &
         & +nabt(4,3)*(wrk(DZ(4),2) - wrk(DZ(-4),2))

  ! zx: (d/dz) * (df/dx)
    w(2) =  nabt(1,3)*(wrk(DZ(1),1) - wrk(DZ(-1),1)) &
         & +nabt(2,3)*(wrk(DZ(2),1) - wrk(DZ(-2),1)) &
         & +nabt(3,3)*(wrk(DZ(3),1) - wrk(DZ(-3),1)) &
         & +nabt(4,3)*(wrk(DZ(4),1) - wrk(DZ(-4),1))

  ! xy: (d/dy) * (df/dx)
    w(3) =  nabt(1,2)*(wrk(DY(1),1) - wrk(DY(-1),1)) &
         & +nabt(2,2)*(wrk(DY(2),1) - wrk(DY(-2),1)) &
         & +nabt(3,2)*(wrk(DY(3),1) - wrk(DY(-3),1)) &
         & +nabt(4,2)*(wrk(DY(4),1) - wrk(DY(-4),1))

    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) - 0.5d0* ( F(4)*w(1) + F(5)*w(2) + F(6)*w(3) )

  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_nonorthogonal

!===================================================================================================================================

# define DR(dt) idx(ix+(sx)*(dt)),idy(iy+(sy)*(dt)),idz(iz+(sz)*(dt))

subroutine stencil_nonorthogonal_highsymmetry(is_array,ie_array,is,ie,idx,idy,idz,ndir &
                                             ,tpsi,htpsi,V_local,lap0,lapt,nabt,sign)
  implicit none
  integer   ,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                           ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4),ndir,sign(3,4:ndir)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)),lap0,lapt(4,6),nabt(4,3)
  complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  !
  integer :: iz,iy,ix,sx,sy,sz,idir
  complex(8) :: v,w

!$OMP parallel
!$OMP do private(iz,iy,ix,w,v,idir,sx,sy,sz)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)

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

    do idir=4,ndir
      sx = sign(1,idir)
      sy = sign(2,idir)
      sz = sign(3,idir)
      v = v + lapt(1,idir) * ( tpsi(DR(1))    &
          &                  + tpsi(DR(-1)) ) &
          & + lapt(2,idir) * ( tpsi(DR(2))    &
          &                  + tpsi(DR(-2)) ) &
          & + lapt(3,idir) * ( tpsi(DR(3))    &
          &                  + tpsi(DR(-3)) ) &
          & + lapt(4,idir) * ( tpsi(DR(4))    &
          &                  + tpsi(DR(-4)) )
    end do

    w =  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1))) &
      & +nabt(2,1)*(tpsi(DX(2)) - tpsi(DX(-2))) &
      & +nabt(3,1)*(tpsi(DX(3)) - tpsi(DX(-3))) &
      & +nabt(4,1)*(tpsi(DX(4)) - tpsi(DX(-4)))

    w =  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1))) &
      & +nabt(2,2)*(tpsi(DY(2)) - tpsi(DY(-2))) &
      & +nabt(3,2)*(tpsi(DY(3)) - tpsi(DY(-3))) &
      & +nabt(4,2)*(tpsi(DY(4)) - tpsi(DY(-4))) + w

    w =  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1))) &
      & +nabt(2,3)*(tpsi(DZ(2)) - tpsi(DZ(-2))) &
      & +nabt(3,3)*(tpsi(DZ(3)) - tpsi(DZ(-3))) &
      & +nabt(4,3)*(tpsi(DZ(4)) - tpsi(DZ(-4))) + w

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v - zI * w
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_nonorthogonal_highsymmetry

!===================================================================================================================================

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

!$OMP parallel
!$OMP do private(iz,iy,ix,w)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
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

    gtpsi(:,ix,iy,iz) = matmul(transpose(matrix_B),w) ! B^{T} * (nabla) psi
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine calc_gradient_psi

!===================================================================================================================================

end module stencil_sub
