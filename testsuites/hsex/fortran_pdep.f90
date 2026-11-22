program test_pdep
  use hsex_pdep, only: pdep_basis,pdep_inverse,pdep_modes
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  complex(8) :: u(12,12),basis(12,12),x(12),y(12),expected(12)
  real(8) :: spectrum(12),values(12),residual(12),a
  integer :: i,j,status,calls
  integer :: invalid
  invalid=0
  do j=1,12
    do i=1,12
      a=2d0*acos(-1d0)*(i-1)*(j-1)/12d0
      u(i,j)=cmplx(cos(a),sin(a),8)/sqrt(12d0)
    enddo
  enddo
  spectrum=0d0;spectrum(:3)=[4d0,2d0,.5d0]
  do i=1,12
    x(i)=cmplx(sin(dble(i)),cos(2d0*i),8)
  enddo
  call pdep_basis(action,12,3,basis(:,:3),values(:3),residual(:3),calls,status)
  if(status/=0)error stop 'rank-three PDEP rejected'
  if(maxval(residual(:3))>1d-10.or.calls>12)error stop 'low-rank modes not resolved'
  call pdep_inverse(basis(:,:3),values(:3),x,y,status)
  expected=x-matmul(u,(spectrum/(1d0+spectrum))*matmul(conjg(transpose(u)),x))
  if(status/=0.or.maxval(abs(y-expected))>1d-11)error stop 'low-rank inverse differs from analytic W'
  spectrum=[(0.1d0*i,i=1,12)]
  call pdep_basis(action,12,12,basis,values,residual,calls,status)
  if(status/=0.or.maxval(residual)>1d-10)error stop 'full-rank PDEP rejected'
  call pdep_inverse(basis,values,x,y,status)
  expected=matmul(u,matmul(conjg(transpose(u)),x)/(1d0+spectrum))
  if(status/=0.or.maxval(abs(y-expected))>1d-11)error stop 'full-rank inverse differs'
  invalid=1
  call pdep_basis(action,12,3,basis(:,:3),values(:3),residual(:3),calls,status)
  if(status==0)error stop 'NaN operator accepted'
  invalid=2
  call pdep_basis(action,12,3,basis(:,:3),values(:3),residual(:3),calls,status)
  if(status==0)error stop 'infinite operator accepted'
  x(12)=ieee_value(0d0,ieee_positive_inf)
  call pdep_inverse(basis(:,:3),values(:3),x,y,status)
  if(status==0)error stop 'infinite trial accepted'
  block
    complex(8) :: big_u(64,3),big_basis(64,3),big_x(64),big_y(64),big_expected(64)
    real(8) :: big_values(3),big_residual(3)
    do j=1,3
      do i=1,64
        a=2d0*acos(-1d0)*(i-1)*(j-1)/64d0
        big_u(i,j)=cmplx(cos(a),sin(a),8)/8d0
      enddo
    enddo
    do i=1,64
      big_x(i)=cmplx(sin(dble(i)),cos(2d0*i),8)
    enddo
    call pdep_basis(big_action,64,3,big_basis,big_values,big_residual,calls,status)
    if(status/=0.or.calls>=64.or.maxval(big_residual)>1d-10)error stop 'matrix-free low-rank construction failed'
    call pdep_inverse(big_basis,big_values,big_x,big_y,status)
    big_expected=big_x-matmul(big_u,[.8d0,.8d0,1d0/3d0]*matmul(conjg(transpose(big_u)),big_x))
    if(status/=0.or.maxval(abs(big_y-big_expected))>1d-10)error stop 'degenerate low-rank W mismatch'
  end block
  block
    complex(8) :: modes(64,3),warm(64,3)
    real(8) :: vals(3),res(3)
    integer :: cold_calls,warm_calls
    call pdep_modes(big_action,64,3,1d-8,40,modes,vals,res,cold_calls,status)
    if(status/=0.or.maxval(res)>1d-8)error stop 'adaptive PDEP failed'
    warm=modes
    call pdep_modes(big_action,64,3,1d-8,40,modes,vals,res,warm_calls,status,warm)
    if(status/=0.or.maxval(res)>1d-8.or.warm_calls>=cold_calls)error stop 'PDEP warm restart did not reuse subspace'
    print *, 'PDEP_ADAPTIVE_WARM_PASS',cold_calls,warm_calls
  end block
  block
    complex(8) :: modes(6,1),old_mode(6,1)
    real(8) :: vals(1),res(1)
    integer :: count
    call pdep_modes(diagonal_action,6,1,1d-9,40,modes,vals,res,count,status)
    if(status/=0.or.abs(vals(1)-6d0)>1d-9)error stop 'PDEP missed leading diagonal mode'
    old_mode=0d0;old_mode(1,1)=1d0
    call pdep_modes(diagonal_action,6,1,1d-9,40,modes,vals,res,count,status,old_mode)
    if(status/=0.or.abs(vals(1)-6d0)>1d-9)error stop 'PDEP reused an obsolete dominant mode'
    print *, 'PDEP_CHANGED_DOMINANT_MODE_PASS'
  end block
  print *, 'PDEP_ANALYTIC_FULL_LOW_RANK_IEEE_PASS'
contains
  subroutine diagonal_action(v,w,code)
    implicit none
    complex(8),intent(in) :: v(:)
    complex(8),intent(out) :: w(:)
    integer,intent(out) :: code
    integer :: k
    code=0
    do k=1,6
      w(k)=dble(k)*v(k)
    enddo
  end subroutine
  subroutine big_action(v,w,code)
    implicit none
    complex(8),intent(in) :: v(:)
    complex(8),intent(out) :: w(:)
    integer,intent(out) :: code
    complex(8) :: f(64),overlap
    real(8) :: angle,strength
    integer :: mode,site
    w=0d0;code=0
    do mode=1,3
      do site=1,64
        angle=2d0*acos(-1d0)*(site-1)*(mode-1)/64d0
        f(site)=cmplx(cos(angle),sin(angle),8)/8d0
      enddo
      strength=4d0
      if(mode==3)strength=.5d0
      overlap=dot_product(f,v);w=w+strength*f*overlap
    enddo
  end subroutine
  subroutine action(v,w,code)
    implicit none
    complex(8),intent(in) :: v(:)
    complex(8),intent(out) :: w(:)
    integer,intent(out) :: code
    code=0;w=matmul(u,spectrum*matmul(conjg(transpose(u)),v))
    if(invalid==1)w(12)=ieee_value(0d0,ieee_quiet_nan)
    if(invalid==2)w(12)=ieee_value(0d0,ieee_positive_inf)
  end subroutine
end program
