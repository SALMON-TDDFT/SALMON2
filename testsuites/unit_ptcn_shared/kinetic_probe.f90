program kinetic_probe
 use ptcn_kinetic_preconditioner
 use ptcn_rectangular_fft
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 type(ptcn_fft_plan)::p
 complex(8)::a(210,2),b(210,2),den(210)
 real(8)::lap(4,3),nab(4,3),f(6),Bmat(3,3),k(3),q(3),symbol,expected,dvec(3),h(3),lap0,dt,diff
 integer::ix,iy,iz,index,err,j
 h=[.7d0,.8d0,.9d0];lap=0;nab=0
 lap(1,:)=1d0/h**2;nab(1,:)=.5d0/h
 f=[1.2d0,1.3d0,1.4d0,.2d0,.1d0,-.3d0];lap0=sum(f(1:3)/h**2)
 Bmat=0;Bmat(1,:)=[1.0d0,.1d0,.0d0];Bmat(2,:)=[.0d0,1.1d0,.2d0];Bmat(3,:)=[.1d0,.0d0,.9d0]
 k=[.1d0,-.2d0,.3d0];q=2*acos(-1d0)*[1d0/5,2d0/6,-1d0/7];dt=.08d0
 dvec=sin(q)/h
 expected=lap0-sum(f(1:3)*cos(q)/h**2)+.5d0*sum(k*k)+sum(matmul(Bmat,k)*dvec) &
  +.5d0*(f(4)*dvec(2)*dvec(3)+f(5)*dvec(3)*dvec(1)+f(6)*dvec(1)*dvec(2))
 call ptcn_kinetic_symbol(q,lap0,lap,nab,f,Bmat,k,symbol,err)
 if(err/=0.or.abs(symbol-expected)>1d-12)error stop 'nonorthogonal kinetic symbol'
 call ptcn_fft_create(p,[5,6,7],2,err)
 if(err/=0)error stop 'FFT plan'
 index=0
 do iz=0,6;do iy=0,5;do ix=0,4
  index=index+1;a(index,1)=exp((0d0,1d0)*sum(q*[ix,iy,iz]));a(index,2)=(.3d0,.7d0)*a(index,1)
  call ptcn_kinetic_symbol(2*acos(-1d0)*real([ix,iy,iz],8)/[5,6,7],lap0,lap,nab,f,Bmat,k,symbol,err)
  if(err/=0)error stop 'symbol grid'
  den(index)=1+(0d0,.5d0)*dt*symbol
 enddo;enddo;enddo
 call ptcn_fft_apply(p,a,den,b,err)
 diff=maxval(abs(b-a/(1+(0d0,.5d0)*dt*expected)))
 if(err/=0.or.diff>1d-12)error stop 'Fourier inverse plane wave'
 print *,'NONORTHOGONAL_KINETIC_PRECONDITIONER_PASS',diff
 q(1)=ieee_value(0d0,ieee_quiet_nan)
 call ptcn_kinetic_symbol(q,lap0,lap,nab,f,Bmat,k,symbol,err)
 if(err==0)error stop 'NaN not rejected'
 q=0;k(2)=ieee_value(0d0,ieee_positive_inf)
 call ptcn_kinetic_symbol(q,lap0,lap,nab,f,Bmat,k,symbol,err)
 if(err==0)error stop 'Infinity not rejected'
 print *,'SCALAR_FINITE_REJECTION_PASS'
 call ptcn_fft_destroy(p)
end program
