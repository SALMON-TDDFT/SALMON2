! Independent quadrature checks of the representation-independent molecular coupling.
program test_near_field_coupling
 use iso_c_binding
 use salmonpy_near_field_coupling,only:exponential_path_moments,sp_near_couple_samples
 use salmonpy_solid_harmonics,only:near_potentials,near_line_integral,sp_near_fit,sp_near_evaluate
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer,parameter::n=32
 real(8),parameter::arguments(12)=[0d0,1d-12,1d-6,.5d0,2d0,8d0,20d0,21d0,100d0,-1d0,-8d0,-30d0]
 real(8)::nodes(n),weights(n),s(n),w(n),E(3,n),H(3,n),moments(0:9),reference(0:9),a
 real(8)::r(3),ce(3,10),ch(3,10),potential,vector(3),vp,Ac(3),decay(3),geo(4),start(3),finish(3)
 real(8)::t,x,argument,radial,ei(3),hi(3),m1,phase,expected_phase,cross(3),max_error,line_error
 real(8)::ev(3),hv(3)
 integer::i,j,k,segment,segments,status,side
 call gauss_legendre(nodes,weights)
 max_error=0d0
 do i=1,size(arguments)
  a=arguments(i)
  call exponential_path_moments(a,9,moments,status)
  if(status/=0)error stop 'moment evaluation failed'
  reference=0d0;segments=max(1,ceiling(abs(a)/2))
  do segment=1,segments
  do j=1,n
   t=(segment-.5d0+.5d0*nodes(j))/segments
   do k=0,9
    reference(k)=reference(k)+.5d0*weights(j)*t**k*exp(-a*t)/segments
   enddo
  enddo
  enddo
  do k=0,9
   max_error=max(max_error,abs(moments(k)-reference(k))/max(abs(reference(k)),tiny(1d0)))
  enddo
 enddo
 if(max_error>2d-13)error stop 'exponential path moments mismatch'
 ! Coupling of analytic r^2 E/H with an exponential envelope, at an off-grid point.
 r=[.37d0,-.29d0,.83d0];decay=[.1d0,-.2d0,.7d0]
 ev=[.01d0,.02d0,-.03d0];hv=[.003d0,-.004d0,.002d0]
 do j=1,n
  s(j)=.5d0*(nodes(j)+1d0);w(j)=.5d0*weights(j)
  radial=sum((s(j)*r)**2)*exp(-dot_product(decay,s(j)*r))
  E(:,j)=radial*ev;H(:,j)=radial*hv
 enddo
 status=sp_near_couple_samples(n,r,s,w,E,H,potential,vector)
 if(status/=0)error stop 'generic sampled-field coupling failed'
 ce=0d0;ch=0d0
 ! Degree2, l=0 is the last coefficient: r^2*Y00, radius1.
 ce(:,10)=ev*sqrt(4d0*acos(-1d0));ch(:,10)=hv*sqrt(4d0*acos(-1d0))
 call near_potentials(r,2,1d0,ce,ch,vp,Ac,decay,status)
 if(status/=0)error stop 'basis-specific field integral failed'
 if(abs(vp-potential)>1d-15.or.maxval(abs(Ac-vector))>1d-17) &
  error stop 'sampled-field and moment coupling disagree'
 ! Failure preserves outputs; scalar NaN and Inf checks.
 do i=1,2
  if(i==1)E(1,1)=ieee_value(0d0,ieee_quiet_nan)
  if(i==2)E(1,1)=ieee_value(0d0,ieee_positive_inf)
  potential=7d0;vector=8d0
  status=sp_near_couple_samples(n,r,s,w,E,H,potential,vector)
  if(status==0.or.potential/=7d0.or.any(vector/=8d0))error stop 'invalid samples changed outputs'
 enddo
 ! Magnetic phase: compare production line rule to a separate32-point rule and closed-form M1.
 geo=[.17d0,-.23d0,.31d0,1d0]
 start=geo(1:3)+[-.3d0,.4d0,.2d0];finish=geo(1:3)+[.7d0,-.2d0,1.1d0]
 ce=0d0;ch=0d0;ch(:,1)=hv*sqrt(4d0*acos(-1d0));line_error=0d0
 do side=1,3
  decay=[0d0,0d0,.7d0]
  if(side==2)decay(3)=8d0
  if(side==3)decay(3)=16d0
  phase=near_line_integral(start,finish,geo,0,ch,decay)
  expected_phase=0d0
  do j=1,n
   t=.5d0*(nodes(j)+1d0);r=start+t*(finish-start)-geo(1:3)
   argument=dot_product(decay,r)
   m1=(1d0-(1d0+argument)*exp(-argument))/argument**2
   cross=[r(2)*hv(3)-r(3)*hv(2),r(3)*hv(1)-r(1)*hv(3),r(1)*hv(2)-r(2)*hv(1)]
   expected_phase=expected_phase-.5d0*weights(j)*dot_product(cross,finish-start)*m1/137.035999084d0
  enddo
  line_error=max(line_error,abs(phase-expected_phase))
 enddo
 if(line_error>2d-17)error stop 'exponential magnetic projector phase mismatch'
 call production_fit_test()
 print *, 'PASS: generic coupling, exponential moments, invalid samples, magnetic phase'
 print *, 'maximum relative moment error=',max_error,' maximum absolute line-phase error=',line_error
contains
 subroutine production_fit_test()
  integer,parameter::nf=729,nq=257
  real(8)::points(3,nf),queries(3,nq),data(3,nf,2),target(3,nq),targetV(nq)
  real(8)::geometry(7),ce(3,165),ch(3,165),stats(2),Efit(3,nq),Hfit(3,nq),Vfit(nq),Afit(3,nq)
  real(8)::offset(3),exponential,relativeE,relativeV,poly_error,exp_error
  integer::i,j,k,nn,degree,mode,status,seedsize
  integer,allocatable::seed(:)
  geometry=[.17d0,-.23d0,.31d0,8d0,0d0,0d0,0d0]
  nn=0;data=0d0
  do k=0,8
  do j=0,8
  do i=0,8
   nn=nn+1;offset=[real(i,8),real(j,8),real(k,8)]/8d0-.5d0
   points(:,nn)=geometry(1:3)+8d0*offset
   exponential=exp(-8d0*(offset(3)+.5d0))
   data(:,nn,1)=exponential*[sin(8d0*offset(1)),0d0,cos(8d0*offset(1))]
  enddo
  enddo
  enddo
  call random_seed(size=seedsize);allocate(seed(seedsize));seed=2952;call random_seed(put=seed)
  call random_number(queries)
  do i=1,nq
   offset=queries(:,i)-.5d0;queries(:,i)=geometry(1:3)+8d0*offset
   exponential=exp(-8d0*(offset(3)+.5d0))
   target(:,i)=exponential*[sin(8d0*offset(1)),0d0,cos(8d0*offset(1))]
   targetV(i)=exp(-4d0)*(1d0-exp(-8d0*offset(3))*cos(8d0*offset(1)))
  enddo
  do degree=2,8,2
   if(degree==6)cycle
   do mode=0,1
    geometry(5:7)=[0d0,0d0,real(mode,8)]
    ce=0d0;ch=0d0
    status=sp_near_fit(nf,degree,points,data,geometry,0,0,ce,ch,stats)
    if(status/=0)error stop 'production exponential fit failed'
    status=sp_near_evaluate(nq,degree,queries,geometry,ce,ch,Efit,Hfit,Vfit,Afit)
    if(status/=0)error stop 'production exponential evaluation failed'
    relativeE=sqrt(sum((Efit-target)**2)/sum(target**2))
    relativeV=sqrt(sum((Vfit-targetV)**2)/sum(targetV**2))
    print *, 'production degree/mode/E error/V error=',degree,mode,relativeE,relativeV
    if(degree==8.and.mode==0)poly_error=relativeE
    if(degree==8.and.mode==1)then
     exp_error=relativeE
     if(relativeE>.01d0.or.relativeV>.01d0)error stop 'strong-decay off-grid accuracy failed'
    endif
   enddo
  enddo
  if(poly_error<5d0*exp_error)error stop 'strong-decay comparison failed'
 end subroutine
 subroutine gauss_legendre(x,w)
  real(8),intent(out)::x(:),w(:)
  real(8)::z,zold,p,pold,temp,derivative
  integer::i,j,iteration,nn
  nn=size(x)
  do i=1,(nn+1)/2
   z=cos(acos(-1d0)*(i-.25d0)/(nn+.5d0))
   do iteration=1,100
    p=1d0;pold=0d0
    do j=1,nn
     temp=pold;pold=p;p=((2*j-1)*z*pold-(j-1)*temp)/j
    enddo
    derivative=nn*(z*p-pold)/(z*z-1d0);zold=z;z=zold-p/derivative
    if(abs(z-zold)<2d-16)exit
   enddo
   if(iteration>100)error stop 'reference quadrature failed'
   x(i)=-z;x(nn+1-i)=z;w(i)=2d0/((1d0-z*z)*derivative**2);w(nn+1-i)=w(i)
  enddo
 end subroutine
end program
