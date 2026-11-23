program test_fdtd_coordinates
  use coordinate_kernel,only:set_coo_em
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_finite
  implicit none
  integer,parameter::nd=1,ista(3)=[-24,-12,-8],iend(3)=[24,12,8]
  real(8),parameter::hgs(3)=[2d0,3d0,4d0]
  real(8)::coo(minval(ista)-nd:maxval(iend)+nd,3),expected
  integer::i,d,mode,poison,ioe(3)
  character::periodic
  do poison=1,2
    do mode=1,3
      periodic='n';ioe=[1,2,1]
      if(mode==2)ioe=[2,1,2]
      if(mode==3)periodic='y'
      if(poison==1)coo=ieee_value(0d0,ieee_quiet_nan)
      if(poison==2)coo=ieee_value(0d0,ieee_positive_inf)
      call set_coo_em(nd,ioe,ista,iend,hgs,coo,periodic)
      do d=1,3
        do i=lbound(coo,1),ubound(coo,1)
          if(.not.ieee_is_finite(coo(i,d)))error stop 'coordinate padding is not initialized'
          expected=0d0
          if(i>=ista(d)-nd.and.i<=iend(d)+nd)then
            if(periodic=='y')then
              expected=dble(i-1)*hgs(d)
            else if(ioe(d)==1)then
              expected=dble(i)*hgs(d)
            else
              expected=(dble(i)-.5d0)*hgs(d)
            endif
          endif
          if(coo(i,d)/=expected)error stop 'valid coordinate changed or padding not zero'
        enddo
      enddo
    enddo
  enddo
  print *,'PASS: coordinate and halo values unchanged; all padding initialized (NaN/Inf poison)'
end program
