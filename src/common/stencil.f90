!
!  Copyright 2017-2020 SALMON developers
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

module stencil_sub
use math_constants,only : zi

contains

!===================================================================================================================================

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine dstencil(is_array,ie_array,is,ie,idx,idy,idz &
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
end subroutine dstencil

!===================================================================================================================================

subroutine zstencil(is_array,ie_array,is,ie,idx,idy,idz &
                   ,tpsi,htpsi,V_local,lap0,lapt,nabt)
  use code_optimization, &
&    only: modx,mody,modz,optimized_stencil_is_callable,stencil_is_parallelized_by_omp
  implicit none
  integer,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                        ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)),lap0,lapt(4,3),nabt(4,3)
  complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))

  if (optimized_stencil_is_callable) then
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
    ! optimized version with hand-coding vectorization (AVX-512, SVE...)
    if (stencil_is_parallelized_by_omp) then
      call zstencil_tuned_omp(is_array,ie_array,is,ie,modx,mody,modz,tpsi,htpsi,V_local,lap0,lapt,nabt)
    else
      call zstencil_tuned_seq(is_array,ie_array,is,ie,modx,mody,modz,is,ie,tpsi,htpsi,V_local,lap0,lapt,nabt)
    end if
#else
    stop 'error: explicit vectorization does not support'
#endif
  else
    ! typical version with fortran compiler vectorization
    if (stencil_is_parallelized_by_omp) then
      call zstencil_typical_omp(is_array,ie_array,is,ie,idx,idy,idz,tpsi,htpsi,V_local,lap0,lapt,nabt)
    else
      call zstencil_typical_seq(is_array,ie_array,is,ie,idx,idy,idz,is,ie,tpsi,htpsi,V_local,lap0,lapt,nabt)
    end if
  end if

  return
end subroutine zstencil

!===================================================================================================================================

subroutine zstencil_nonorthogonal(is_array,ie_array,is,ie,idx,idy,idz,wrk &
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
end subroutine zstencil_nonorthogonal

!===================================================================================================================================

# define DR(dt) idx(ix+(sx)*(dt)),idy(iy+(sy)*(dt)),idz(iz+(sz)*(dt))

! (future works)
subroutine zstencil_nonorthogonal_highsymmetry(is_array,ie_array,is,ie,idx,idy,idz,ndir &
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
end subroutine zstencil_nonorthogonal_highsymmetry

!===================================================================================================================================

subroutine zstencil_microAc(is_array,ie_array,is,ie,idx,idy,idz &
                                ,tpsi,htpsi,V_local,Ac,div_Ac,lap0,lapt,nabt,k)
  use code_optimization, only: stencil_is_parallelized_by_omp
  use salmon_global, only: yn_symmetrized_stencil
  implicit none
  integer   ,intent(in)  :: is_array(3),ie_array(3),is(3),ie(3) &
                           ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)  :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                          & ,Ac(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                          & ,div_Ac(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                          & ,lap0,lapt(4,3),nabt(4,3),k(3)
  complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))

  if(yn_symmetrized_stencil=='y') then
    call symmetrized_stencil 
  else
    ! typical version with fortran compiler vectorization
    if (stencil_is_parallelized_by_omp) then
      call zstencil_microAc_typical_omp(is_array,ie_array,is,ie,idx,idy,idz,tpsi,htpsi,V_local,Ac,div_ac,lap0,lapt,nabt,k)
    else
      call zstencil_microAc_typical_seq(is_array,ie_array,is,ie,idx,idy,idz,is,ie,tpsi,htpsi,V_local,Ac,div_ac,lap0,lapt,nabt,k)
    end if
  end if

contains

  subroutine symmetrized_stencil
    implicit none
    integer :: ix,iy,iz,i
    real(8) :: kAc(3),Ac_tmp(1:3,is(1):ie(1),is(2):ie(2),is(3)-4:ie(3)+4)
    complex(8) :: w(3),v,psi0,x
    
# define DTZ(dt) ix,iy,(iz+(dt))
    
!!!!!!!!!!! nproc_rgrid(1:3) must be =1

  !$OMP parallel
 !$OMP do collapse(2) private(iz,iy,ix)
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      Ac_tmp(1:3,ix,iy,iz) = Ac(1:3,ix,iy,iz)
    end do
    end do
    end do
    !$OMP end do
    !$OMP end parallel
    
    do i=1,4
      iz = is(3)
      Ac_tmp(:,:,:,iz-i) = Ac(:,:,:,iz)
      iz = ie(3)
      Ac_tmp(:,:,:,iz+i) = Ac(:,:,:,iz)
    end do

  !$OMP parallel
  !$OMP do collapse(2) private(iz,iy,ix,w,v,psi0,kAc,x)
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      psi0 = tpsi(ix,iy,iz)
      kAc = k + Ac(:,ix,iy,iz)

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
           
    ! divergence of (k+Ac)*psi
      x    =  nabt(1,1)* ( (k(1)+Ac_tmp(1,DX(1) ))*tpsi(DX(1)) - (k(1)+Ac_tmp(1,DX(-1) ))*tpsi(DX(-1)) ) &
           & +nabt(2,1)* ( (k(1)+Ac_tmp(1,DX(2) ))*tpsi(DX(2)) - (k(1)+Ac_tmp(1,DX(-2) ))*tpsi(DX(-2)) ) &
           & +nabt(3,1)* ( (k(1)+Ac_tmp(1,DX(3) ))*tpsi(DX(3)) - (k(1)+Ac_tmp(1,DX(-3) ))*tpsi(DX(-3)) ) &
           & +nabt(4,1)* ( (k(1)+Ac_tmp(1,DX(4) ))*tpsi(DX(4)) - (k(1)+Ac_tmp(1,DX(-4) ))*tpsi(DX(-4)) )
           
      x    =  nabt(1,2)* ( (k(2)+Ac_tmp(2,DY(1) ))*tpsi(DY(1)) - (k(2)+Ac_tmp(2,DY(-1) ))*tpsi(DY(-1)) ) &
           & +nabt(2,2)* ( (k(2)+Ac_tmp(2,DY(2) ))*tpsi(DY(2)) - (k(2)+Ac_tmp(2,DY(-2) ))*tpsi(DY(-2)) ) &
           & +nabt(3,2)* ( (k(2)+Ac_tmp(2,DY(3) ))*tpsi(DY(3)) - (k(2)+Ac_tmp(2,DY(-3) ))*tpsi(DY(-3)) ) &
           & +nabt(4,2)* ( (k(2)+Ac_tmp(2,DY(4) ))*tpsi(DY(4)) - (k(2)+Ac_tmp(2,DY(-4) ))*tpsi(DY(-4)) ) + x
           
      x    =  nabt(1,3)* ( (k(3)+Ac_tmp(3,DTZ(1)))*tpsi(DZ(1)) - (k(3)+Ac_tmp(3,DTZ(-1)))*tpsi(DZ(-1)) ) &
           & +nabt(2,3)* ( (k(3)+Ac_tmp(3,DTZ(2)))*tpsi(DZ(2)) - (k(3)+Ac_tmp(3,DTZ(-2)))*tpsi(DZ(-2)) ) &
           & +nabt(3,3)* ( (k(3)+Ac_tmp(3,DTZ(3)))*tpsi(DZ(3)) - (k(3)+Ac_tmp(3,DTZ(-3)))*tpsi(DZ(-3)) ) &
           & +nabt(4,3)* ( (k(3)+Ac_tmp(3,DTZ(4)))*tpsi(DZ(4)) - (k(3)+Ac_tmp(3,DTZ(-4)))*tpsi(DZ(-4)) ) + x

      htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )* psi0 - 0.5d0* v           &
                      & + 0.5d0* ( kAc(1)**2 + kAc(2)**2 + kAc(3)**2 ) * psi0   &
                      & - 0.5d0*zi* ( ( kAc(1)*w(1) + kAc(2)*w(2) + kAc(3)*w(3) ) + x ) 
    end do
    end do
    end do
  !$OMP end do
  !$OMP end parallel
    
  end subroutine symmetrized_stencil
  
end subroutine zstencil_microAc

!===================================================================================================================================

# define DDX(dt) mg%idx(ix+(dt)),iy,iz
# define DDY(dt) ix,mg%idy(iy+(dt)),iz
# define DDZ(dt) ix,iy,mg%idz(iz+(dt))

subroutine calc_gradient_field(mg,nabt,box,grad)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8)      ,intent(in) :: nabt(4,3)
  real(8)      ,intent(in) :: box(mg%is_array(1):mg%ie_array(1), &
                                & mg%is_array(2):mg%ie_array(2), &
                                & mg%is_array(3):mg%ie_array(3))
  real(8)                  :: grad(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  !
  integer :: ix,iy,iz
  real(8) :: w(3)
!$OMP parallel
!$OMP do private(iz,iy,ix,w)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    w(1) =  nabt(1,1)*(box(DDX(1)) - box(DDX(-1))) &
         & +nabt(2,1)*(box(DDX(2)) - box(DDX(-2))) &
         & +nabt(3,1)*(box(DDX(3)) - box(DDX(-3))) &
         & +nabt(4,1)*(box(DDX(4)) - box(DDX(-4)))
    w(2) =  nabt(1,2)*(box(DDY(1)) - box(DDY(-1))) &
         & +nabt(2,2)*(box(DDY(2)) - box(DDY(-2))) &
         & +nabt(3,2)*(box(DDY(3)) - box(DDY(-3))) &
         & +nabt(4,2)*(box(DDY(4)) - box(DDY(-4)))
    w(3) =  nabt(1,3)*(box(DDZ(1)) - box(DDZ(-1))) &
         & +nabt(2,3)*(box(DDZ(2)) - box(DDZ(-2))) &
         & +nabt(3,3)*(box(DDZ(3)) - box(DDZ(-3))) &
         & +nabt(4,3)*(box(DDZ(4)) - box(DDZ(-4)))
    grad(:,ix,iy,iz) = w
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel
end subroutine calc_gradient_field

subroutine calc_laplacian_field(mg,lapt,lap0,box,lap)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8)      ,intent(in) :: lapt(4,3),lap0
  real(8)      ,intent(in) :: box(mg%is_array(1):mg%ie_array(1), &
                                & mg%is_array(2):mg%ie_array(2), &
                                & mg%is_array(3):mg%ie_array(3))
  real(8)                  :: lap(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  !
  integer :: ix,iy,iz
  real(8) :: v
!$OMP parallel
!$OMP do private(iz,iy,ix,v)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    v = lapt(1,1)*( box(DDX(1)) + box(DDX(-1)) ) &
    & + lapt(2,1)*( box(DDX(2)) + box(DDX(-2)) ) &
    & + lapt(3,1)*( box(DDX(3)) + box(DDX(-3)) ) &
    & + lapt(4,1)*( box(DDX(4)) + box(DDX(-4)) ) &
    & + lapt(1,2)*( box(DDY(1)) + box(DDY(-1)) ) &
    & + lapt(2,2)*( box(DDY(2)) + box(DDY(-2)) ) &
    & + lapt(3,2)*( box(DDY(3)) + box(DDY(-3)) ) &
    & + lapt(4,2)*( box(DDY(4)) + box(DDY(-4)) ) &
    & + lapt(1,3)*( box(DDZ(1)) + box(DDZ(-1)) ) &
    & + lapt(2,3)*( box(DDZ(2)) + box(DDZ(-2)) ) &
    & + lapt(3,3)*( box(DDZ(3)) + box(DDZ(-3)) ) &
    & + lapt(4,3)*( box(DDZ(4)) + box(DDZ(-4)) )
    lap(ix,iy,iz) = v + lap0* box(ix,iy,iz)
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel
end subroutine calc_laplacian_field

end module stencil_sub
