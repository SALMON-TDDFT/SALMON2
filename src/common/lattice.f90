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
module lattice
  implicit none

contains

!===================================================================================================================================

SUBROUTINE init_lattice(system,stencil)
  use math_constants,only : pi
  use structures
  implicit none
  type(s_dft_system) :: system
  type(s_stencil)    :: stencil
  !
  real(8),dimension(3,3) :: A,B,F,wrk
  real(8) :: a1(3),a2(3),a3(3),detA,normA(3),f_uu,f_vv,f_ww,f_uv,f_uw,f_vw

! al = [ a1, a2, a3 ]
  A = system%primitive_a ! primitive lattice vectors
  a1 = A(1:3,1)
  a2 = A(1:3,2)
  a3 = A(1:3,3)
  call calc_inverse(A,wrk,detA)
  system%det_a = detA
  system%Hvol = detA/dble(system%ngrid)
  system%primitive_b = 2d0*pi* transpose(wrk) ! reciprocal primitive lattice vectors
  ! [ b1 b2 b3 ]^{T} = 2*pi* [ a1 a2 a3 ]^{-1}

  normA(1) = sqrt(sum(a1**2))
  normA(2) = sqrt(sum(a2**2))
  normA(3) = sqrt(sum(a3**2))

! cf. A. Natan et al., PRB 78, 075109 (2008).
! A = [ u, v, w ], B = A^{-1}
  A(1:3,1) = a1(1:3) / normA(1) ! u
  A(1:3,2) = a2(1:3) / normA(2) ! v
  A(1:3,3) = a3(1:3) / normA(3) ! w
  call calc_inverse(A,B,detA)

  wrk = transpose(B)
  F = matmul(B,wrk)
  f_uu = F(1,1)
  f_vv = F(2,2)
  f_ww = F(3,3)
  f_uv = F(1,2) + F(2,1)
  f_uw = F(1,3) + F(3,1)
  f_vw = F(2,3) + F(3,2)

  stencil%rmatrix_A = A
  stencil%rmatrix_B = B
  stencil%coef_F(1) = f_uu
  stencil%coef_F(2) = f_vv
  stencil%coef_F(3) = f_ww
  stencil%coef_F(4) = f_vw ! yz
  stencil%coef_F(5) = f_uw ! zx
  stencil%coef_F(6) = f_uv ! xy

  return
end SUBROUTINE init_lattice

SUBROUTINE calc_inverse(a,b,detA) ! b = a^{-1}
  implicit none
  real(8),intent(in) :: a(3,3)
  real(8)            :: b(3,3),detA

  detA=a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3) &
    -a(1,3)*a(2,2)*a(3,1)-a(2,3)*a(3,2)*a(1,1)-a(3,3)*a(1,2)*a(2,1)

  b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
  b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
  b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)

  b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
  b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
  b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)

  b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)

  b=b/detA
  return
end SUBROUTINE calc_inverse

SUBROUTINE init_kvector(num_kgrid,system)
  use structures
  use sym_kvector, only: init_sym_kvector
  implicit none
  integer,intent(in) :: num_kgrid(3)
  type(s_dft_system) :: system
  !
  integer :: ix,iy,iz
  integer :: ik,nk
  real(8) :: shift_k(3),k(3),B(3,3)

  nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
  system%nk = nk
  if ( allocated(system%vec_k) ) deallocate(system%vec_k)
  if ( allocated(system%wtk)   ) deallocate(system%wtk)
  allocate(system%vec_k(3,nk),system%wtk(nk))

  system%wtk  = 1d0/dble(nk)

  B = system%primitive_b
  shift_k(1:3) = 0.5d0
  do ik=1,nk
    ix=mod(ik-1,num_kgrid(1))+1
    iy=mod((ik-1)/num_kgrid(1),num_kgrid(2))+1
    iz=mod((ik-1)/(num_kgrid(1)*num_kgrid(2)),num_kgrid(3))+1
    k(1) = (dble(ix)-shift_k(1))/dble(num_kgrid(1))-0.5d0
    k(2) = (dble(iy)-shift_k(2))/dble(num_kgrid(2))-0.5d0
    k(3) = (dble(iz)-shift_k(3))/dble(num_kgrid(3))-0.5d0
    system%vec_k(1,ik) = k(1)*B(1,1) + k(2)*B(1,2) + k(3)*B(1,3)
    system%vec_k(2,ik) = k(1)*B(2,1) + k(2)*B(2,2) + k(3)*B(2,3)
    system%vec_k(3,ik) = k(1)*B(3,1) + k(2)*B(3,2) + k(3)*B(3,3)
  end do

  call init_sym_kvector( system%vec_k, system%wtk, system%nk, B ) 

  return
end SUBROUTINE init_kvector

end module lattice
