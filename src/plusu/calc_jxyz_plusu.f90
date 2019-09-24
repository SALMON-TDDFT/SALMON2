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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module calc_jxyz_plusu_sub

  implicit none

  private
  public :: calc_jxyz_plusu

contains

  subroutine calc_jxyz_plusu(pp,ppg,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz,al0,matrix_A0)
    use salmon_global,only : natom,kion,rion,iperiodic,yn_domain_parallel
    use structures,only : s_pp_info,s_pp_grid
    implicit none
    type(s_pp_info) :: pp
    type(s_pp_grid) :: ppg
    real(8),intent(in) :: alx,aly,alz
    integer,intent(in) :: nl,ml
    integer,intent(in) :: lx(nl),ly(nl),lz(nl)
    integer,intent(in) :: mx(ml),my(ml),mz(ml)
    real(8),intent(in) :: hx,hy,hz
    real(8),intent(in),optional :: al0(3,3),matrix_A0(3,3)
    integer :: a,i,ik,ix,iy,iz,j
    integer :: nc1,nc2,nc3
    real(8) :: tmpx,tmpy,tmpz
    real(8) :: r,x,y,z,u,v,w
    real(8) :: rshift(3),matrix_a(3,3),rr(3),al(3,3)

    matrix_a = 0.0d0
    matrix_a(1,1) = 1.0d0
    matrix_a(2,2) = 1.0d0
    matrix_a(3,3) = 1.0d0
    if ( present(matrix_A0) ) matrix_a = matrix_A0

    al = 0.0d0
    al(1,1) = alx
    al(2,2) = aly
    al(3,3) = alz
    if ( present(al0) ) al = al0

    if( iperiodic == 0 )then
      nc1=0
      nc2=0
      nc3=0
    else if( iperiodic == 3 )then
      r=maxval( pp%rps_ao )
      nc1=nint(r/alx); if ( nc1*alx < r ) nc1=nc1+1
      nc2=nint(r/aly); if ( nc2*aly < r ) nc2=nc2+1
      nc3=nint(r/alz); if ( nc3*alz < r ) nc3=nc3+1
    end if

    if( iperiodic == 0 )then
      if( mod(lx(nl)-lx(1)+1,2) == 1 )then
        rshift(1)=0.0d0
      else
        rshift(1)=-0.5d0*Hx
      end if
      if( mod(ly(nl)-ly(1)+1,2) == 1 )then
        rshift(2)=0.0d0
      else
        rshift(2)=-0.5d0*Hy
      end if
      if( mod(lz(nl)-lz(1)+1,2) == 1 )then
        rshift(3)=0.0d0
      else
        rshift(3)=-0.5d0*Hz
      end if
    else if( iperiodic == 3 )then 
      if( yn_domain_parallel == 'y' )then
        rshift(1)=-Hx
        rshift(2)=-Hy
        rshift(3)=-Hz
      else
        rshift(:)=0.0d0
      end if
    end if

    do a=1,natom
      ik=kion(a)
      j=0
      do ix=-nc1,nc1
      do iy=-nc2,nc2
      do iz=-nc3,nc3
        rr(1) = ix*al(1,1) + iy*al(1,2) + iz*al(1,3)
        rr(2) = ix*al(2,1) + iy*al(2,2) + iz*al(2,3)
        rr(3) = ix*al(3,1) + iy*al(3,2) + iz*al(3,3)
        tmpx = rion(1,a)+ rr(1)
        tmpy = rion(2,a)+ rr(2)
        tmpz = rion(3,a)+ rr(3)
        do i=1,ml
          u = mx(i)*hx + rshift(1)
          v = my(i)*hy + rshift(2)
          w = mz(i)*hz + rshift(3)
          rr(1) = u*matrix_a(1,1) + v*matrix_a(1,2) + w*matrix_a(1,3)
          rr(2) = u*matrix_a(2,1) + v*matrix_a(2,2) + w*matrix_a(2,3)
          rr(3) = u*matrix_a(3,1) + v*matrix_a(3,2) + w*matrix_a(3,3)
          x = rr(1) - tmpx
          y = rr(2) - tmpy
          z = rr(3) - tmpz
          r = sqrt(x*x+y*y+z*z)
          if( r < pp%rps_ao(ik) )then
            j=j+1
            if( j <= ppg%nps_ao )then
              ppg%jxyz_ao(1,j,a)=mx(i)
              ppg%jxyz_ao(2,j,a)=my(i)
              ppg%jxyz_ao(3,j,a)=mz(i)
              ppg%jxx_ao( j,a)=ix
              ppg%jyy_ao( j,a)=iy
              ppg%jzz_ao( j,a)=iz
              ppg%rxyz_ao(1,j,a) = x
              ppg%rxyz_ao(2,j,a) = y
              ppg%rxyz_ao(3,j,a) = z
            else
            
            end if
          end if
        end do
      end do
      end do
      end do
      ppg%mps_ao(a)=j
    end do


  end subroutine calc_jxyz_plusu

end module calc_jxyz_plusu_sub
