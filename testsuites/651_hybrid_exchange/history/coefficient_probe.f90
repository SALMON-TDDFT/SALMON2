program coefficient_probe
 use exx_factor_history
 implicit none
 type(s_factor_history)::h
 complex(8)::x(16,1)
 real(8)::reference(2),error
 integer::stat,i
 x=(1d0,0d0)
 h%count=3;h%last=0;h%teachers=8
 allocate(h%x(16,1,3));h%x=(1d0,0d0)
 h%g(1,:)=[1d0,.999999999999d0];h%g(2,:)=[.999999999999d0,1d0]
 h%b=[1.0000000000001d0,.9999999999999d0]
 reference=[5.00004993529472364e-01_8,4.99994996471027686e-01_8]
 call history_accept(h,x,1d0,8,sumgrid,stat)
 error=maxval(abs(h%coeff-reference))
 print *, 'near-collinear error',error,'coeff',h%coeff
 if(stat/=0.or..not.h%ready.or.error>2d-9)stop 1
 ! Identical equations at a small scale: direct determinant underflows.
 h=s_factor_history();h%count=3;h%last=0;h%teachers=8
 allocate(h%x(16,1,3));h%x=(1d0,0d0)
 h%g=0;h%g(1,1)=1d-200;h%g(2,2)=2d-200
 h%b=[1d-200,2d-200]
 call history_accept(h,x,1d0,8,sumgrid,stat)
 reference=[1d0/(1d0+3d-8),2d0/(2d0+3d-8)]
 error=maxval(abs(h%coeff-reference))
 print *, 'scaled error',error,'ready',h%ready
 if(stat/=0.or..not.h%ready.or.error>1d-13)stop 2
 contains
 subroutine sumgrid(a)
 complex(8),intent(inout)::a(:,:)
 end subroutine
end program
