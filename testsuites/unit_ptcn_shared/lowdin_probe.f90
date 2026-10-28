program lowdin_probe
 use ptcn_core,only:pt_orthogonalize
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 complex(8)::x(4,2,2),initial(4,2,2),metric(2,2)
 real(8)::errgram,volume,nan,inf
 integer::code,k,j
 volume=.5d0;initial=0
 do k=1,2
  initial(1,1,k)=1;initial(2,1,k)=(.01d0,.02d0)
  initial(1,2,k)=(.03d0,.01d0);initial(2,2,k)=1.01d0
 enddo
 x=initial
 call pt_orthogonalize(x,volume,code)
 if(code/=0)error stop 'Lowdin failed'
 errgram=0
 do k=1,2
  metric=volume*matmul(conjg(transpose(x(:,:,k))),x(:,:,k))
  do j=1,2;metric(j,j)=metric(j,j)-1;enddo
  errgram=max(errgram,maxval(abs(metric)))
 enddo
 if(errgram>1d-12)error stop 'Lowdin did not restore orthogonality'
 nan=ieee_value(0d0,ieee_quiet_nan);inf=ieee_value(0d0,ieee_positive_inf)
 x=initial;x(4,2,2)=cmplx(nan,0d0,8)
 call pt_orthogonalize(x,volume,code)
 if(code==0)error stop 'Lowdin accepted NaN'
 x=initial;x(4,2,2)=cmplx(0d0,inf,8)
 call pt_orthogonalize(x,volume,code)
 if(code==0)error stop 'Lowdin accepted Infinity'
 x=initial
 call pt_orthogonalize(x,nan,code)
 if(code==0)error stop 'Lowdin accepted NaN volume'
 x=initial;x(:,2,:)=x(:,1,:)
 call pt_orthogonalize(x,volume,code)
 if(code==0)error stop 'Lowdin accepted rank deficiency'
 print *,'LOWDIN_SCALAR_MULTI_K_FINITE_PASS',errgram
end program
