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
!-----------------------------------------------------------------------------------------
module salmon_math
  implicit none

  public :: erf_salmon, &
            erfc_salmon, &
            bessel_j1_salmon, &
            matrix_inverse, &
            xjl, dxjl

  interface matrix_inverse
     module procedure matrix_inverse_double
     module procedure matrix_inverse_complex
  end interface
contains
!--------------------------------------------------------------------------------
!! Error function and its complement are implemented based on the reference
!! W. J. Cody, Math. Comp. 23 631 (1969)
  function erf_salmon(x) result(y)
    !$acc routine seq
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: xa

    xa = abs(x)
    if(xa < 0.46875d0)then
       y = sign( erf_salmon_short(xa) , x )
    else if(xa < 4d0)then
       y = erfc_salmon_mid(xa)
       y = sign(1d0 - y,x)
    else
       y = erfc_salmon_long(xa)
       y = sign(1d0 - y,x)
    end if

  end function erf_salmon
!--------------------------------------------------------------------------------
  function erfc_salmon(x) result(y)
    !$acc routine seq
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: xa

    xa = abs(x)
    if(xa < 0.46875d0)then
       y = 1d0 - erf_salmon(xa)
    else if(xa < 4d0)then
       y = erfc_salmon_mid(xa)
    else
       y = erfc_salmon_long(xa)
    end if

    if(x<0d0) y = 2d0 - y

  end function erfc_salmon
!--------------------------------------------------------------------------------
  function erf_salmon_short(x) result(y) ! 0d0 <=x<0.5d0
    !$acc routine seq
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2
    real(8),parameter :: &
       a(0:4) = (/  &
       & 3.209377589138469472562d3, &
       & 3.774852376853020208137d2, &
       & 1.138641541510501556495d2, &
       & 3.161123743870565596947d0, &
       & 1.857777061846031526730d-1 &
       &/), &
       b(0:4) = (/  &
       & 2.844236833439170622273d3, &
       & 1.282616526077372275645d3, &
       & 2.440246379344441733056d2, &
       & 2.360129095234412093499d1, &
       & 1d0 &
       &/)

    x2 = x**2
    y = x*(a(0) + a(1)*x2 + a(2)*x2**2 + a(3)*x2**3 + a(4)*x2**4 ) &
         /(b(0) + b(1)*x2 + b(2)*x2**2 + b(3)*x2**3 + b(4)*x2**4 )


  end function erf_salmon_short
!--------------------------------------------------------------------------------
  function erfc_salmon_mid(x) result(y) ! 0.46875d0 <x< 4d0
    !$acc routine seq
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2
    real(8),parameter :: &
       p(0:7) = (/  &
       & 3.004592610201616005d2, &
       & 4.519189537118729422d2, &
       & 3.393208167343436870d2, &
       & 1.529892850469404039d2, &
       & 4.316222722205673530d1, &
       & 7.211758250883093659d0, &
       & 5.641955174789739711d-1, &
       &-1.368648573827167067d-7 & 
       &/), &
       q(0:7) = (/ &
       & 3.004592609569832933d2, &
       & 7.909509253278980272d2, &
       & 9.313540948506096211d2, &
       & 6.389802644656311665d2, &
       & 2.775854447439876434d2, &
       & 7.700015293522947295d1, &
       & 1.278272731962942351d1, &
       & 1d0 &
       &/)

    x2 = x**2

    y = (p(0) + p(1)*x + p(2)*x**2 + p(3)*x**3 + p(4)*x**4 & 
         + p(5)*x**5 + p(6)*x**6 + p(7)*x**7)/ &
         (q(0) + q(1)*x + q(2)*x**2 + q(3)*x**3 + q(4)*x**4 & 
         + q(5)*x**5 + q(6)*x**6 + q(7)*x**7)
    y = exp(-x2)*y


  end function erfc_salmon_mid
!--------------------------------------------------------------------------------
  function erfc_salmon_long(x) result(y) ! 4d0 <x
    !$acc routine seq
    use math_constants, only : pi
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2,x2_i
    real(8),parameter :: &
       r(0:4) = (/ &
       & -2.99610707703542174d-3, &
       & -4.94730910623250734d-2, &
       & -2.26956593539686930d-1, &
       & -2.78661308609647788d-1, &
       & -2.23192459734184686d-2  &
       &/), &
       s(0:4) = (/ &
       & 1.06209230528467918d-2, &
       & 1.91308926107829841d-1, &
       & 1.05167510706793207d0,  &
       & 1.98733201817135256d0,  &
       & 1d0 &
       &/)

    x2 = x**2

    x2_i = 1d0/x2
    y = 1d0/sqrt(pi) &
         +x2_i*(r(0) + r(1)*x2_i + r(2)*x2_i**2 + r(3)*x2_i**3 + r(4)*x2_i**4) &
         /(s(0) + s(1)*x2_i + s(2)*x2_i**2 + s(3)*x2_i**3 + s(4)*x2_i**4) 
    y = exp(-x2)/x*y

  end function erfc_salmon_long
  
  
!--------------------------------------------------------------------------------
  !! 1st order Bessel function
  ! The algorithm and numerical table are implemented based on the reference:
  ! Methods and Programs for Mathematical Functions, Prentice-Hall, 1989
  real(8) function bessel_j1_salmon(x) result(y)
    implicit none
    real(8), intent(in) :: x
    real(8), parameter :: &
      &  rp(0:3) = (/ &
      & -8.999712d+08, +4.522283d+11, -7.274942d+13, +3.682957d+15  &
      & /), &
      & rq(0:8) = (/ &
      & +1.000000d+00, +6.208365d+02, +2.569873d+05, +8.351468d+07,  &
      & +2.215116d+10, +4.749141d+12, +7.843696d+14, +8.952223d+16,  &
      & +5.322786d+18 /), &
      & pp(0:6) = (/ &
      & +7.621256d-04, +7.313971d-02, +1.127196d+00, +5.112080d+00,  &
      & +8.424046d+00, +5.214516d+00, +1.000000d+00 /), &
      & pq(0:6) = (/ &
      & +5.713231d-04, +6.884559d-02, +1.105142d+00, +5.073864d+00,  &
      & +8.399856d+00, +5.209828d+00, +1.000000d+00 /), &
      & qp(0:7) = (/ &
      & +5.108626d-02, +4.982139d+00, +7.582383d+01, +3.667796d+02,  &
      & +7.108563d+02, +5.974896d+02, +2.116888d+02, +2.520702d+01  &
      & /), &
      & qq(0:7) = (/ &
      & +1.000000d+00, +7.423733d+01, +1.056449d+03, +4.986411d+03,  &
      & +9.562319d+03, +7.997042d+03, +2.826193d+03, +3.360936d+02  &
      & /)
    real(8), parameter :: z1 = +1.468197d+01
    real(8), parameter :: z2 = +4.921846d+01
    real(8), parameter :: thpio4 =  +2.356194d+00
    real(8), parameter :: sq2opi =  +7.978846d-01
    real(8) :: w, z, p, q, xn
    real(8) :: rpz, rqz, ppz, pqz, qpz, qqz

    if (x < 0) then
      w = -x
    else
      w = +x
    end if

    if( w <= 5.0 ) then
      z = x * x
      rpz = ((rp(0)*z+rp(1))*z+rp(2))*z+rp(3)
      rqz = (((((((rq(0)*z+rq(1))*z+rq(2))*z+rq(3))*z+rq(4))*z &
        & +rq(5))*z+rq(6))*z+rq(7))*z+rq(8)
      w = rpz / rqz
      w = w * x * (z - z1) * (z - z2)
      y = w
      return
    end if

    w = 5.0 / x
    z = w * w
    ppz = (((((pp(0)*z+pp(1))*z+pp(2))*z+pp(3))*z+pp(4))*z+pp(5))*z+pp(6)
    pqz = (((((pq(0)*z+pq(1))*z+pq(2))*z+pq(3))*z+pq(4))*z+pq(5))*z+pq(6)
    qpz = ((((((qp(0)*z+qp(1))*z+qp(2))*z+qp(3))*z &
    & +qp(4))*z+qp(5))*z+qp(6))*z+qp(7)
    qqz = ((((((qq(0)*z+qq(1))*z+qq(2))*z+qq(3))*z &
    & +qq(4))*z+qq(5))*z+qq(6))*z+qq(7)
    p = ppz/pqz
    q = qpz/qqz
    xn = x - thpio4
    p = p * cos(xn) - w * q * sin(xn)
    y = p * sq2opi / sqrt(x)
    return
  end function bessel_j1_salmon
!--------------------------------------------------------------------------------
  subroutine matrix_inverse_double(rmat)
    implicit none
    real(8),intent(inout) :: rmat(:,:)
    integer :: nn
! for lapack
    integer :: lwork, info
    integer, allocatable :: ipiv(:) ! dimension N
    real(8), allocatable :: work(:) ! dimension LWORK  


    nn = size(rmat,dim=1)
    lwork = nn * max(nn, 64)
    
    allocate(ipiv(nn),work(lwork))

    call dgetrf(nn, nn, rmat, nn, ipiv, info)  ! factorize
    call dgetri(nn, rmat, nn, ipiv, work, lwork, info)  ! inverse

    deallocate(ipiv,work)

  end subroutine matrix_inverse_double
!--------------------------------------------------------------------------------
  subroutine matrix_inverse_complex(zmat)
    implicit none
    complex(8),intent(inout) :: zmat(:,:)
    integer :: nn
! for lapack
    integer :: lwork, info
    integer, allocatable :: ipiv(:) ! dimension N
    complex(8), allocatable :: zwork(:) ! dimension LWORK  


    nn = size(zmat,dim=1)
    lwork = nn * max(nn, 64)
    
    allocate(ipiv(nn),zwork(lwork))

    call zgetrf(nn, nn, zmat, nn, ipiv, info)  ! factorize
    call zgetri(nn, zmat, nn, ipiv, zwork, lwork, info)  ! inverse

    deallocate(ipiv,zwork)

  end subroutine matrix_inverse_complex  
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  real(8) function xjl(x,l) result(r_xjl)
    implicit none
!argument
    integer,intent(in) :: l
    real(8),intent(in) :: x
!local variable
    real(8),parameter :: eps=1.0d-1
  
    if (l >= 5) then
      write(*,*) 'xjl function not prepared for l>=5'
      stop
    endif
  
    if (x < eps) then
      select case(l)
      case(-1)
         r_xjl = 1.d0 - x**2/2.d0 + x**4/24.d0
      case(0)
         r_xjl = x - x**3/6.d0 + x**5/120.d0 -x**7/5040.d0 + x**9/362880.d0
      case(1)
         r_xjl = (2.d0/6.d0)*x**2 - (2.d0/60.d0)*x**4 + (1.d0/840.d0)*x**6 - (2.d0/90720.d0)*x**8
      case(2)
         r_xjl = (4.d0/60.d0)*x**3 - (4.d0/840.d0)*x**5 + (2.d0/15120.d0)*x**7 - (2.d0/997920.d0)*x**9
      case(3)
         r_xjl = (8.d0/840.d0)*x**4 - (8.d0/15120.d0)*x**6 + (4.d0/332640.d0)*x**8
      case(4)
         r_xjl = (16.d0/15120.d0)*x**5 - (16.d0/332640.d0)*x**7 + (8.d0/8648640.d0)*x**9
      end select
    else
      select case(l)
      case(-1)
         r_xjl = cos(x)
      case(0)
         r_xjl = sin(x)
      case(1)
         r_xjl = (1.d0/x)*sin(x) - cos(x)
      case(2)
         r_xjl = (3.d0/x**2 - 1.d0)*sin(x) - (3.d0/x)*cos(x)
      case(3)
         r_xjl = (15.d0/x**3 - 6.d0/x)*sin(x) - (15.d0/x**2 - 1.d0)*cos(x)
      case(4)
         r_xjl = (105.d0/x**4 - 45.d0/x**2 + 1.d0)*sin(x) - (105.d0/x**3 - 10.d0/x)*cos(x)
      end select
    end if
  
    return
  end function xjl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  real(8) function dxjl(x,l) result(r_dxjl)
    implicit none
!argument
    integer,intent(in) :: l
    real(8),intent(in) :: x
!local variable
    real(8),parameter :: eps=1.0d-1
  
    if (l >= 5) then
      write(*,*) 'dxjl function not prepared for l>=5'
      stop
    endif
  
    if (x < eps) then
      select case(l)
      case(-1)
         r_dxjl = x - x**3/6.d0 
      case(0)
         r_dxjl = 1.d0 - x**2/2.d0 + x**4/24.d0 - x**6/720.d0 + x**8/40320.d0
      case(1)
         r_dxjl = (2.d0/3.d0)*x**1 - (2.d0/15.d0)*x**3 + (1.d0/140.d0)*x**5 - (2.d0/11340.d0)*x**7
      case(2)
         r_dxjl = (4.d0/20.d0)*x**2 - (4.d0/168.d0)*x**4 + (2.d0/2160.d0)*x**6 - (2.d0/110880.d0)*x**7
      case(3)
         r_dxjl = (8.d0/210.d0)*x**3 - (8.d0/2520.d0)*x**5 + (4.d0/41580.d0)*x**7
      case(4)
         r_dxjl = (16.d0/3024.d0)*x**4 - (16.d0/47520.d0)*x**6 + (8.d0/960960.d0)*x**8
      end select
    else
      select case(l)
      case(-1)
         r_dxjl = -sin(x)
      case(0)
         r_dxjl = cos(x)
      case(1)
         r_dxjl = -(1.d0/x**2 - 1.d0)*sin(x) + (1.d0/x)*cos(x)
      case(2)
         r_dxjl = -(6.d0/x**3 - 3.d0/x)*sin(x) + (6.d0/x**2 - 1.d0)*cos(x)
      case(3)
         r_dxjl = -(45.d0/x**4 - 21.d0/x**2 + 1.d0)*sin(x) + (45.d0/x**3 - 6.d0/x)*cos(x)
      case(4)
         r_dxjl = -(420.d0/x**5 - 195.d0/x**3 + 10.d0/x)*sin(x) + (420.d0/x**4 - 55.d0/x**2 + 1.d0)*cos(x)
      end select
    end if
  
    return
  end function dxjl

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! This function Ylm is related to the real spherical harmonics Ylm0 by
! Ylm=sqrt(4*pi/(2l+1))*r**l*Ylm0 and is a monomial of x,y,z
  Function Ylm(x,y,z,il,im)
    implicit none
    real(8),intent(IN) :: x,y,z
    integer,intent(IN) :: il,im
    integer :: ilm
    real(8) :: Ylm,r2
  
    ilm=il*il+il+1+im
    r2=x*x+y*y+z*z
    select case( ilm )
    case(1)  ; Ylm=1.d0                                     ! ilm=1  (0  0)
    case(2)  ; Ylm=-y                                       ! ilm=2  (1 -1)
    case(3)  ; Ylm=z                                        ! ilm=3  (1  0)
    case(4)  ; Ylm=-x                                       ! ilm=4  (1  1)
    case(5)  ; Ylm=sqrt(3.d0)*x*y                           ! ilm=5  (2 -2)
    case(6)  ; Ylm=-sqrt(3.d0)*y*z                          ! ilm=6  (2 -1)
    case(7)  ; Ylm=(3.d0*z*z-r2)/2.d0                       ! ilm=7  (2  0)
    case(8)  ; Ylm=-sqrt(3.d0)*x*z                          ! ilm=8  (2  1)
    case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)                ! ilm=9  (2  2)
    case(10) ; Ylm=-sqrt(5.d0/8.d0)*y*(3*x*x-y*y)           ! ilm=10 (3 -3)
    case(11) ; Ylm=sqrt(15.d0)*x*y*z                        ! ilm=11 (3 -2)
    case(12) ; Ylm=-sqrt(3.d0/8.d0)*y*(5.d0*z*z-r2)         ! ilm=12 (3 -1)
    case(13) ; Ylm=z*(5.d0*z*z-3*r2)/2.d0                   ! ilm=13 (3  0)
    case(14) ; Ylm=-sqrt(3.d0/8.d0)*x*(5.d0*z*z-r2)         ! ilm=14 (3  1)
    case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)             ! ilm=15 (3  2)
    case(16) ; Ylm=-sqrt(5.d0/8.d0)*x*(x*x-3.d0*y*y)        ! ilm=16 (3  3)
    case(17) ; Ylm=sqrt(35.d0)/2.d0*x*y*(x*x-y*y)           ! ilm=17 (4 -4)
    case(18) ; Ylm=-sqrt(35.d0/8.d0)*y*z*(3.d0*x*x-y*y)     ! ilm=18 (4 -3)
    case(19) ; Ylm=sqrt(5.d0)/2.d0*x*y*(7.d0*z*z-r2)        ! ilm=19 (4 -2)
    case(20) ; Ylm=-sqrt(5.d0/8.d0)*y*z*(7.d0*z*z-3.d0*r2)  ! ilm=20 (4 -1)
    case(21) ; Ylm=(35.d0*z**4-30.d0*z*z*r2+3.d0*r2*r2)/8.d0! ilm=21 (4  0)
    case(22) ; Ylm=-sqrt(5.d0/8.d0)*x*z*(7.d0*z**2-3.d0*r2) ! ilm=22 (4  1)
    case(23) ; Ylm=sqrt(5.d0)/4.d0*(7.d0*z*z-r2)*(x*x-y*y)  ! ilm=23 (4  2)
    case(24) ; Ylm=-sqrt(35.d0/8.d0)*x*z*(x*x-3.d0*y*y)     ! ilm=24 (4  3)
    case(25) ; Ylm=sqrt(35.d0)/8.d0*(x**4+y**4-6.d0*x*x*y*y)! ilm=25 (4  4)
    end select
  
    return
  End Function Ylm
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  Function dYlm(x,y,z,il,im,idir)
    implicit none
    real(8),intent(IN) :: x,y,z
    integer,intent(IN) :: il,im,idir
    integer :: ilm
    real(8) :: dYlm
  
    ilm=il*il+il+1+im
    if ((ilm > 16) .or. (il >= 4)) then
      write(*,*) 'dYlm routine not prepared for il>=4'
      stop
    endif
    if((ilm ==  1 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=1  (0  0)
    if((ilm ==  1 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=1  (0  0)
    if((ilm ==  1 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=1  (0  0)
    if((ilm ==  2 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=2  (1 -1)
    if((ilm ==  2 ).and.(idir == 2)) dYlm=-1.d0                           ! ilm=2  (1 -1)
    if((ilm ==  2 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=2  (1 -1)
    if((ilm ==  3 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=3  (1  0)
    if((ilm ==  3 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=3  (1  0)
    if((ilm ==  3 ).and.(idir == 3)) dYlm=1.d0                            ! ilm=3  (1  0)
    if((ilm ==  4 ).and.(idir == 1)) dYlm=-1.d0                           ! ilm=4  (1  1)
    if((ilm ==  4 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=4  (1  1)
    if((ilm ==  4 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=4  (1  1)
    if((ilm ==  5 ).and.(idir == 1)) dYlm=sqrt(3.d0)*y                    ! ilm=5  (2 -2)
    if((ilm ==  5 ).and.(idir == 2)) dYlm=sqrt(3.d0)*x                    ! ilm=5  (2 -2)
    if((ilm ==  5 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=5  (2 -2)
    if((ilm ==  6 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=6  (2 -1)
    if((ilm ==  6 ).and.(idir == 2)) dYlm=-sqrt(3.d0)*z                   ! ilm=6  (2 -1)
    if((ilm ==  6 ).and.(idir == 3)) dYlm=-sqrt(3.d0)*y                   ! ilm=6  (2 -1)
    if((ilm ==  7 ).and.(idir == 1)) dYlm=-x                              ! ilm=7  (2  0)
    if((ilm ==  7 ).and.(idir == 2)) dYlm=-y                              ! ilm=7  (2  0)
    if((ilm ==  7 ).and.(idir == 3)) dYlm=2.d0*z                          ! ilm=7  (2  0)
    if((ilm ==  8 ).and.(idir == 1)) dYlm=-sqrt(3.d0)*z                   ! ilm=8  (2  1)
    if((ilm ==  8 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=8  (2  1)
    if((ilm ==  8 ).and.(idir == 3)) dYlm=-sqrt(3.d0)*x                   ! ilm=8  (2  1)
    if((ilm ==  9 ).and.(idir == 1)) dYlm=sqrt(3.d0)*x                    ! ilm=9  (2  2)
    if((ilm ==  9 ).and.(idir == 2)) dYlm=-sqrt(3.d0)*y                   ! ilm=9  (2  2)
    if((ilm ==  9 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=9  (2  2)
    if((ilm == 10 ).and.(idir == 1)) dYlm=-3.d0*sqrt(5.d0/2.d0)*x*y                         ! ilm=10 (3 -3)
    if((ilm == 10 ).and.(idir == 2)) dYlm=-3.d0/2.d0*sqrt(5.d0/2.d0)*(x*x-y*y)               ! ilm=10 (3 -3)
    if((ilm == 10 ).and.(idir == 3)) dYlm=0.d0                                              ! ilm=10 (3 -3)
    if((ilm == 11 ).and.(idir == 1)) dYlm=sqrt(15.d0)*y*z                                   ! ilm=11 (3 -2)
    if((ilm == 11 ).and.(idir == 2)) dYlm=sqrt(15.d0)*z*x                                   ! ilm=11 (3 -2)
    if((ilm == 11 ).and.(idir == 3)) dYlm=sqrt(15.d0)*x*y                                   ! ilm=11 (3 -2)
    if((ilm == 12 ).and.(idir == 1)) dYlm=sqrt(3.d0/2.d0)*x*y                               ! ilm=12 (3 -1)
    if((ilm == 12 ).and.(idir == 2)) dYlm=-1.d0/2.d0*sqrt(3.d0/2.d0)*(4.d0*z*z-x*x-3.d0*y*y)! ilm=12 (3 -1)
    if((ilm == 12 ).and.(idir == 3)) dYlm=-4.d0*sqrt(3.d0/2.d0)*y*z                         ! ilm=12 (3 -1)
    if((ilm == 13 ).and.(idir == 1)) dYlm=-3.d0*z*x                                         ! ilm=13 (3  0)
    if((ilm == 13 ).and.(idir == 2)) dYlm=-3.d0*y*z                                         ! ilm=13 (3  0)
    if((ilm == 13 ).and.(idir == 3)) dYlm=3.d0/2.d0*(2.d0*z*z-x*x-y*y)                      ! ilm=13 (3  0)
    if((ilm == 14 ).and.(idir == 1)) dYlm=-1.d0/2.d0*sqrt(3.d0/2.d0)*(4.d0*z*z-3.d0*x*x-y*y)! ilm=14 (3  1)
    if((ilm == 14 ).and.(idir == 2)) dYlm=sqrt(3.d0/2.d0)*x*y                               ! ilm=14 (3  1)
    if((ilm == 14 ).and.(idir == 3)) dYlm=-4.d0*sqrt(3.d0/2.d0)*z*x                         ! ilm=14 (3  1)
    if((ilm == 15 ).and.(idir == 1)) dYlm=sqrt(15.d0)*z*x                                   ! ilm=15 (3  2)
    if((ilm == 15 ).and.(idir == 2)) dYlm=-sqrt(15.d0)*y*z                                  ! ilm=15 (3  2)
    if((ilm == 15 ).and.(idir == 3)) dYlm=1.d0/2.d0*sqrt(15.d0)*(x*x-y*y)                   ! ilm=15 (3  2)
    if((ilm == 16 ).and.(idir == 1)) dYlm=-3.d0/2.d0*sqrt(5.d0/2.d0)*(x*x-y*y)              ! ilm=16 (3  3)
    if((ilm == 16 ).and.(idir == 2)) dYlm=3.d0*sqrt(5.d0/2.d0)*x*y                          ! ilm=16 (3  3)
    if((ilm == 16 ).and.(idir == 3)) dYlm=0.d0                                              ! ilm=16 (3  3)
  
    return
  End Function dYlm

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

  subroutine spline(Np,xn,yn,an,bn,cn,dn)
    implicit none
    integer,intent(in) :: Np
    real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
    real(8),intent(out) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
    integer :: i,Npm2,info
    real(8) :: dxn(0:Np-1),dyn(0:Np-1),u(1:Np-2),v(1:Np-2),Amat(1:Np-2,1:Np-2)
    real(8) :: Amat_t(1:Np-2,1:Np-2)
  ! for lapack
    integer :: LWORK
    integer, allocatable :: IPIV(:) ! dimension N
    real(8), allocatable :: WORK(:) ! dimension LWORK
  ! for check inverse matrix problem
  !  integer :: j,k
  !  real(8) :: Amat_chk(1:Np-2,1:Np-2)
  !  real(8) :: ss

    Npm2 = Np-2
    LWORK = Npm2*Npm2*6
    allocate(IPIV(Npm2),WORK(LWORK))


    do i = 0,Np-2
      dxn(i) = xn(i+1) - xn(i)
      dyn(i) = yn(i+1) - yn(i)
    end do

    do i = 1,Npm2
      v(i) = 6d0*(dyn(i)/dxn(i) - dyn(i-1)/dxn(i-1))
    end do

    Amat = 0d0
    Amat(1,1) = 2d0*(dxn(1) + dxn(0))
    Amat(1,2) = dxn(1)
    do i = 2,Npm2-1
      Amat(i,i+1) = dxn(i)
      Amat(i,i  ) = 2d0*(dxn(i)+dxn(i-1))
      Amat(i,i-1) = dxn(i-1)
    end do
    Amat(Npm2,Npm2  ) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
    Amat(Npm2,Npm2-1) = dxn(Npm2-1)

  ! inverse matrix problem
    Amat_t = Amat


    call DGETRF(Npm2, Npm2, Amat_t, Npm2, IPIV, info)  ! factorize
    call DGETRI(Npm2, Amat_t, Npm2, IPIV, WORK, LWORK, info)  ! inverse

  !  check inverse matrix problem
  !  do i = 1,Npm2
  !    do j = 1,Npm2
  !      ss = 0d0
  !      do k = 1,Npm2
  !        ss = ss + Amat(i,k)*Amat_t(k,j)
  !      end do
  !      Amat_chk(i,j) = ss
  !    end do
  !  end do
  !
  !  do i = 1,Npm2
  !    write(*,'(999e16.6e3)')(Amat_chk(i,j),j=1,Npm2)
  !  end do
  !
  !  stop


    do i = 1,Npm2
      u(i) = sum(Amat_t(i,:)*v(:))
    end do

  ! for b
    bn(0) = 0d0
    bn(1:Np-2) = 0.5d0*u(1:Np-2)
  ! for a
    do i = 0,Npm2-1
      an(i) = (u(i+1) -2d0*bn(i))/(6d0*dxn(i))
    end do
    an(Npm2) = (0d0 -2d0*bn(Npm2))/(6d0*dxn(Npm2))
  ! for d
    dn(0:Npm2) = yn(0:Npm2)
  ! for c
    i=0
    cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*0.d0)/6d0
    do i = 1,Npm2-1
       cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0
    end do
    cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

    return
  end subroutine spline
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

end module salmon_math
!--------------------------------------------------------------------------------
  
