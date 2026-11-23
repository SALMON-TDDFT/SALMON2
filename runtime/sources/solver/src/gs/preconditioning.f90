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
module preconditioning_sub
  implicit none

contains

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine dstencil_preconditioning(is_array,ie_array,is,ie,idx,idy,idz,hgs &
                                   ,tpsi,htpsi,alpha)
  use structures
  implicit none

  integer   ,intent(in)    :: is_array(3),ie_array(3),is(3),ie(3) &
                             ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  real(8)   ,intent(in)    :: hgs(3)
  real(8)   ,intent(in)    :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(inout) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)    :: alpha
  !
  integer :: ix,iy,iz
  real(8) :: d,v
  real(8) :: omega

! low-pass filter precondtioning
! N. Tancogne-Dejean et al., J. Chem. Phys. 152, 124119 (2020). equation (76)
! The actual formula is written in the Octopus code as
! \tilde{\psi} =  \omega D^{-1}(2\psi - \omega D^{-1} (-0.5\Laplacian) \psi),
! where \omega = 2(1-\alpha), and D is the diagonal part of (-0.5\Laplacian).

  omega=2.d0*(1.d0-alpha)
  ! d: diagonal part of -0.5*\Laplacian
  d = -0.5d0* ( (-2.d0)/hgs(1)**2 &
              + (-2.d0)/hgs(2)**2 &
              + (-2.d0)/hgs(3)**2 )
#ifdef USE_OPENACC
!$acc parallel
!$acc loop collapse(2) private(iz,iy,ix,v)
#else
!$OMP parallel
!$OMP do collapse(2) private(iz,iy,ix,v)
#endif
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    ! v: \Laplacian \psi
    v =   (tpsi(DX(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DX(-1)))/hgs(1)**2 &
        + (tpsi(DY(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DY(-1)))/hgs(2)**2 &
        + (tpsi(DZ(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DZ(-1)))/hgs(3)**2
    htpsi(ix,iy,iz)=omega/d*(2.d0*tpsi(ix,iy,iz)-omega/d*(-0.5d0*v))
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$OMP end do
!$OMP end parallel
#endif

end subroutine dstencil_preconditioning

!===================================================================================================================================

subroutine zstencil_preconditioning(is_array,ie_array,is,ie,idx,idy,idz,hgs &
                                   ,tpsi,htpsi,alpha)
  use structures
  implicit none

  integer   ,intent(in)    :: is_array(3),ie_array(3),is(3),ie(3) &
                             ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  real(8)   ,intent(in)    :: hgs(3)
  complex(8),intent(in)    :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  complex(8),intent(inout) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)    :: alpha
  !
  integer :: ix,iy,iz
  real(8) :: d
  complex(8) :: v
  real(8) :: omega

! low-pass filter precondtioning
! N. Tancogne-Dejean et al., J. Chem. Phys. 152, 124119 (2020). equation (76)
! The actual formula is written in the Octopus code as
! \tilde{\psi} =  \omega D^{-1}(2\psi - \omega D^{-1} (-0.5\Laplacian) \psi),
! where \omega = 2(1-\alpha), and D is the diagonal part of (-0.5\Laplacian).

  omega=2.d0*(1.d0-alpha)
  ! d: diagonal part of -0.5*\Laplacian
  d = -0.5d0* ( (-2.d0)/hgs(1)**2 &
              + (-2.d0)/hgs(2)**2 &
              + (-2.d0)/hgs(3)**2 )
#ifdef USE_OPENACC
!$acc parallel
!$acc loop collapse(2) private(iz,iy,ix,v)
#else
!$OMP parallel
!$OMP do collapse(2) private(iz,iy,ix,v)
#endif
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    ! v: \Laplacian \psi
    v =   (tpsi(DX(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DX(-1)))/hgs(1)**2 &
        + (tpsi(DY(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DY(-1)))/hgs(2)**2 &
        + (tpsi(DZ(1))-2.d0*tpsi(ix,iy,iz)+tpsi(DZ(-1)))/hgs(3)**2
    htpsi(ix,iy,iz)=omega/d*(2.d0*tpsi(ix,iy,iz)-omega/d*(-0.5d0*v))
  end do
  end do
  end do

#ifdef USE_OPENACC
!$acc end parallel
#else
!$OMP end do
!$OMP end parallel
#endif

end subroutine zstencil_preconditioning

!===================================================================================================================================

subroutine zstencil_nonorthogonal_preconditioning(is_array,ie_array,is,ie,idx,idy,idz &
                                   ,tpsi,htpsi,lap0,lapt,nabt,F,alpha)
  use structures
  implicit none

  integer   ,intent(in)    :: is_array(3),ie_array(3),is(3),ie(3) &
                             ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  complex(8),intent(in)    :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  complex(8),intent(inout) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8)   ,intent(in)    :: lap0,lapt(1,3),nabt(1,3)
  real(8)   ,intent(in)    :: F(6)
  real(8)   ,intent(in)    :: alpha
  !
  integer :: ix,iy,iz
  real(8) :: omega
  real(8) :: d
  complex(8) :: w(3),v(3)
  complex(8) :: wrk(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),2)

! low-pass filter precondtioning
! N. Tancogne-Dejean et al., J. Chem. Phys. 152, 124119 (2020). equation (76)
! The actual formula is written in the Octopus code as
! \tilde{\psi} =  \omega D^{-1}(2\psi - \omega D^{-1} (-0.5\Laplacian) \psi),
! where \omega = 2(1-\alpha), and D is the diagonal part of (-0.5\Laplacian).

  omega=2.d0*(1.d0-alpha)
  ! d: diagonal part of -0.5*\Laplacian
  d = lap0
#ifdef USE_OPENACC
!$acc parallel
!$acc loop collapse(2) private(iz,iy,ix,v,w)
#else
!$OMP parallel
!$OMP do collapse(2) private(iz,iy,ix,v,w)
#endif
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    v(1) =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1)))
    v(2) =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1)))
    v(3) =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1)))

    w(1) =  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1)))
    w(2) =  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1)))
    w(3) =  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1)))

    htpsi(ix,iy,iz) = lap0 * tpsi(ix,iy,iz) &
                    - 0.5d0* ( F(1)*v(1) + F(2)*v(2) + F(3)*v(3) ) 

    wrk(ix,iy,iz,1) = w(1) ! df/dx
    wrk(ix,iy,iz,2) = w(2) ! df/dy
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$OMP end do
!$OMP end parallel
#endif
    
#ifdef USE_OPENACC
!$acc parallel
!$acc loop collapse(2) private(iz,iy,ix,w)
#else
!$OMP parallel
!$OMP do collapse(2) private(iz,iy,ix,w)
#endif
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
  ! yz: (d/dz) * (df/dy)
    w(1) =  nabt(1,3)*(wrk(DZ(1),2) - wrk(DZ(-1),2))
  ! zx: (d/dz) * (df/dx)
    w(2) =  nabt(1,3)*(wrk(DZ(1),1) - wrk(DZ(-1),1))
  ! xy: (d/dy) * (df/dx)
    w(3) =  nabt(1,2)*(wrk(DY(1),1) - wrk(DY(-1),1))
    htpsi(ix,iy,iz) = htpsi(ix,iy,iz) - 0.5d0* ( F(4)*w(1) + F(5)*w(2) + F(6)*w(3) )
    htpsi(ix,iy,iz) = omega/d*(2.d0*tpsi(ix,iy,iz)-omega/d*htpsi(ix,iy,iz))
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$OMP end do
!$OMP end parallel
#endif
end subroutine zstencil_nonorthogonal_preconditioning

end module preconditioning_sub
