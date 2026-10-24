program probe
 use exx_stationarity,only:spinor_stationarity
 use,intrinsic::ieee_arithmetic
 implicit none
 complex(8)::p(3,2,1),h(3,2,1)
 real(8)::f(2,1),g,a,c,outside,within
 integer::ierr
 p=0;p(1,1,1)=1;p(2,2,1)=1;f(:,1)=[.8d0,.2d0]
 h=0;h(1,1,1)=2;h(2,2,1)=3
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.max(g,a,c)>1d-14)error stop 'eigenstates'
 h(2,1,1)=.3d0;h(1,2,1)=.3d0
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.abs(c-sqrt(2d0*(.6d0*.3d0)**2/sum(f*f)))>1d-14)error stop 'fractional commutator'
 if(abs(outside)>1d-14.or.abs(within-c)>1d-14)error stop 'inside decomposition'
 f=.5d0
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.c>1d-14)error stop 'equal-occupation gauge'
 h(3,1,1)=.4d0
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.abs(c-.4d0)>1d-14)error stop 'outside subspace'
 if(abs(within)>1d-14.or.abs(outside-c)>1d-14)error stop 'outside decomposition'
 h(1,2,1)=cmplx(.3d0,.1d0,8)
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.abs(a-.1d0)>1d-14)error stop 'nonhermiticity'
 p(1,1,1)=2
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr/=0.or.abs(g-3d0)>1d-14)error stop 'orthogonality'
 h(1,1,1)=cmplx(ieee_value(1d0,ieee_quiet_nan),0d0,8)
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr==0)error stop 'NaN accepted'
 h=0;f(1,1)=ieee_value(1d0,ieee_positive_inf)
 call spinor_stationarity(p,h,f,1d0,g,a,c,ierr,outside=outside,within=within)
 if(ierr==0)error stop 'Inf accepted'
 print *,'PASS spinor stationarity known commutators, gauge, nonhermiticity, norm, NaN/Inf'
end program
