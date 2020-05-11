module calc_nps_plusu_sub

  implicit none
  private
  public :: calc_nps_plusu

contains

  subroutine calc_nps_plusu(pp,ppg,system,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz,al0,matrix_A0)
    use salmon_global,only : natom,kion,iperiodic
    use structures,only : s_pp_info,s_pp_grid,s_dft_system
    implicit none
    type(s_pp_info) :: pp
    type(s_pp_grid) :: ppg
    type(s_dft_system),intent(in) :: system
    real(8),intent(in) :: alx,aly,alz
    integer,intent(in) :: nl,ml
    integer,intent(in) :: lx(nl),ly(nl),lz(nl)
    integer,intent(in) :: mx(ml),my(ml),mz(ml)
    real(8),intent(in) :: hx,hy,hz
    real(8),intent(in),optional :: al0(3,3),matrix_A0(3,3)
    integer :: a,i,ik,ix,iy,iz,j
    integer :: nc1,nc2,nc3
    real(8) :: tmpx,tmpy,tmpz
    real(8) :: x,y,z,r,u,v,w
    real(8) :: rshift(3),matrix_a(3,3),rr(3),al(3,3)
    integer :: mps_tmp(natom)

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
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
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
      tmpx = system%Rion(1,a) + rr(1)
      tmpy = system%Rion(2,a) + rr(2)
      tmpz = system%Rion(3,a) + rr(3)
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
        if ( r < pp%rps_ao(ik)+1.d-12 ) j=j+1
      end do
    end do
    end do
    end do
    mps_tmp(a)=j
    end do

    ppg%nps_ao=maxval(mps_tmp(:))

  end subroutine calc_nps_plusu

end module calc_nps_plusu_sub
