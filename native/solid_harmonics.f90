! Real orthonormal regular solid harmonics for a source-free electric potential.
! V_e(r)=sum c_lm (|r-origin|/radius)^l Y_lm; physical E=grad V_e.
module salmonpy_solid_harmonics
 use salmonpy_near_field_coupling,only:exponential_path_moments,near_coupling_map
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 integer,parameter::solid_lmax=8,solid_size=80
contains
 subroutine solid_basis(x,lmax,value,gradient)
  real(c_double),intent(in)::x(3)
  integer,intent(in)::lmax
  real(c_double),intent(out)::value(solid_size),gradient(3,solid_size)
  complex(c_double_complex)::p(0:solid_lmax,0:solid_lmax),g(3,0:solid_lmax,0:solid_lmax),z,dp(3)
  real(c_double)::r2,norm,phase
  integer::l,m,k,q
  p=cmplx(0d0,0d0,kind=c_double);g=p(0,0);p(0,0)=1d0
  z=cmplx(x(1),x(2),kind=c_double);dp=[cmplx(1d0,0d0,kind=c_double), &
      cmplx(0d0,1d0,kind=c_double),cmplx(0d0,0d0,kind=c_double)]
  r2=dot_product(x,x);value=0d0;gradient=0d0
  do m=1,lmax
   p(m,m)=-(2*m-1)*z*p(m-1,m-1)
   g(:,m,m)=-(2*m-1)*(dp*p(m-1,m-1)+z*g(:,m-1,m-1))
  enddo
  do m=0,lmax
   if(m<lmax)then
    p(m+1,m)=(2*m+1)*x(3)*p(m,m)
    g(:,m+1,m)=(2*m+1)*x(3)*g(:,m,m)
    g(3,m+1,m)=g(3,m+1,m)+(2*m+1)*p(m,m)
   endif
   do l=m+2,lmax
    p(l,m)=((2*l-1)*x(3)*p(l-1,m)-(l+m-1)*r2*p(l-2,m))/(l-m)
    g(:,l,m)=((2*l-1)*x(3)*g(:,l-1,m)-(l+m-1)*(2d0*x*p(l-2,m)+r2*g(:,l-2,m)))/(l-m)
    g(3,l,m)=g(3,l,m)+(2*l-1)*p(l-1,m)/(l-m)
   enddo
  enddo
  do l=1,lmax
   do m=0,l
    norm=(2*l+1)/(4d0*acos(-1d0))
    do q=l-m+1,l+m
     norm=norm/q
    enddo
    norm=sqrt(norm);k=l*l+l+m
    if(m==0)then
     value(k)=norm*real(p(l,m),c_double);gradient(:,k)=norm*real(g(:,l,m),c_double)
    else
     phase=sqrt(2d0)*norm*(-1d0)**m
     value(k)=phase*real(p(l,m),c_double);gradient(:,k)=phase*real(g(:,l,m),c_double)
     k=l*l+l-m
     value(k)=phase*aimag(p(l,m));gradient(:,k)=phase*aimag(g(:,l,m))
    endif
   enddo
  enddo
 end subroutine
 integer(c_int) function sp_solid_protocol() bind(C)
  sp_solid_protocol=1
 end function
 integer(c_int) function sp_solid_fit(n,lmax,positions,fields,origin,radius,coeff,statistics) bind(C) result(status)
  integer(c_int),value::n,lmax
  real(c_double),intent(in)::positions(3,n),fields(3,n),origin(3)
  real(c_double),value::radius
  real(c_double),intent(inout)::coeff(*),statistics(2)
  real(c_double),allocatable::a(:,:),qmat(:,:),rhs(:)
  real(c_double)::v(solid_size),g(3,solid_size),r(solid_size,solid_size),b(solid_size),c(solid_size)
  real(c_double)::s(2),norm,scale,projection,denominator
  integer::i,j,k,p,pass,allocation_status
  status=1
  if(n<1.or.lmax<1.or.lmax>solid_lmax)return
  p=lmax*(lmax+2)
  if(3*n<p.or..not.ieee_is_finite(radius))return
  if(radius<=0d0)return
  do k=1,3
   if(.not.ieee_is_finite(origin(k)))return
  enddo
  allocate(a(3*n,p),qmat(3*n,p),rhs(3*n),stat=allocation_status)
  status=2
  if(allocation_status/=0)return
  status=1
  do i=1,n
   do k=1,3
    if(.not.ieee_is_finite(positions(k,i)).or..not.ieee_is_finite(fields(k,i)))return
   enddo
   call solid_basis((positions(:,i)-origin)/radius,lmax,v,g)
   a(3*i-2:3*i,:)=g(:,:p)/radius;rhs(3*i-2:3*i)=fields(:,i)
  enddo
  status=4
  do j=1,p
  do i=1,3*n
   if(.not.ieee_is_finite(a(i,j)))return
  enddo
  enddo
  r=0d0;b=0d0;c=0d0;qmat=a;scale=maxval(abs(a))
  if(scale==0d0)return
  ! Twice-reorthogonalized modified Gram-Schmidt; no normal-equation squaring.
  do j=1,p
   do pass=1,2
    do k=1,j-1
     projection=dot_product(qmat(:,k),qmat(:,j))
     r(k,j)=r(k,j)+projection;qmat(:,j)=qmat(:,j)-projection*qmat(:,k)
    enddo
   enddo
   norm=sqrt(dot_product(qmat(:,j),qmat(:,j)))
   status=1
   if(norm<1d-12*scale)return
   status=4
   if(.not.ieee_is_finite(norm))return
   r(j,j)=norm;qmat(:,j)=qmat(:,j)/norm;b(j)=dot_product(qmat(:,j),rhs)
  enddo
  do j=p,1,-1
   c(j)=(b(j)-dot_product(r(j,j+1:p),c(j+1:p)))/r(j,j)
   if(.not.ieee_is_finite(c(j)))return
  enddo
  rhs=matmul(a,c(:p))-rhs
  denominator=sum(fields*fields)
  s(1)=sqrt(sum(rhs*rhs)/max(denominator,tiny(1d0)));s(2)=maxval(abs(rhs))
  do i=1,2
   if(.not.ieee_is_finite(s(i)))return
  enddo
  coeff(1:p)=c(1:p);statistics=s;status=0
 end function
 integer(c_int) function sp_solid_evaluate(n,lmax,positions,origin,radius,coeff,potential,fields) bind(C) result(status)
  integer(c_int),value::n,lmax
  real(c_double),intent(in)::positions(3,n),origin(3),coeff(*)
  real(c_double),value::radius
  real(c_double),intent(inout)::potential(n),fields(3,n)
  real(c_double),allocatable::vp(:),ef(:,:)
  real(c_double)::v(solid_size),g(3,solid_size)
  integer::i,k,p,allocation_status
  status=1
  if(n<1.or.lmax<1.or.lmax>solid_lmax)return
  if(.not.ieee_is_finite(radius))return
  if(radius<=0d0)return
  p=lmax*(lmax+2)
  do i=1,p
   if(.not.ieee_is_finite(coeff(i)))return
  enddo
  do k=1,3
   if(.not.ieee_is_finite(origin(k)))return
   do i=1,n
    if(.not.ieee_is_finite(positions(k,i)))return
   enddo
  enddo
  allocate(vp(n),ef(3,n),stat=allocation_status);status=2
  if(allocation_status/=0)return
  status=4
  do i=1,n
   call solid_basis((positions(:,i)-origin)/radius,lmax,v,g)
   vp(i)=dot_product(v(:p),coeff(1:p));ef(:,i)=matmul(g(:,:p),coeff(1:p))/radius
   if(.not.ieee_is_finite(vp(i)))return
   do k=1,3
    if(.not.ieee_is_finite(ef(k,i)))return
   enddo
  enddo
  potential=vp;fields=ef;status=0
 end function
 ! Vector E/H coefficients include l=0; no longitudinal projection.
 subroutine near_basis(x,lmax,value)
  real(c_double),intent(in)::x(3)
  integer,intent(in)::lmax
  real(c_double),intent(out)::value(165)
  real(c_double)::v(solid_size),g(3,solid_size),radial,r2
  integer::degree,l,k,m,index
  call solid_basis(x,lmax,v,g)
  value=0d0;index=0;r2=dot_product(x,x)
  ! Total degree, then descending angular l=degree,degree-2,... and m=-l..l.
  do degree=0,lmax
   do l=degree,0,-2
    k=(degree-l)/2;radial=r2**k
    if(l==0)then
     index=index+1;value(index)=radial/sqrt(4d0*acos(-1d0))
    else
     do m=-l,l
      index=index+1;value(index)=radial*v(l*l+l+m)
     enddo
    endif
   enddo
  enddo
 end subroutine
 subroutine near_potentials(r,lmax,radius,ce,ch,potential,Ac,decay_vector,status)
  real(c_double),intent(in)::r(3),radius,ce(3,*),ch(3,*)
  integer,intent(in)::lmax
  real(c_double),intent(out)::potential,Ac(3)
  real(c_double),optional,intent(in)::decay_vector(3)
  integer,optional,intent(out)::status
  real(c_double)::v(165),ei(3),hi(3),moments(0:9),argument
  integer::l,first,last,code,k
  logical::weighted
  if(present(status))status=4
  call near_basis(r/radius,lmax,v);ei=0d0;hi=0d0
  weighted=.false.
  if(present(decay_vector))weighted=any(decay_vector/=0d0)
  if(weighted)then
   argument=dot_product(decay_vector,r)
   call exponential_path_moments(argument,lmax+1,moments,code)
   if(code/=0)then
    if(present(status))return
    error stop 'SALMONpy exponential field path is outside numerical range'
   endif
  endif
  do l=0,lmax
   first=l*(l+1)*(l+2)/6+1;last=(l+1)*(l+2)*(l+3)/6
   if(weighted)then
    ei=ei+matmul(ce(:,first:last),v(first:last))*moments(l)
    hi=hi+matmul(ch(:,first:last),v(first:last))*moments(l+1)
   else
    ! Preserve exact polynomial arithmetic when the optional envelope is absent.
    ei=ei+matmul(ce(:,first:last),v(first:last))/(l+1)
    hi=hi+matmul(ch(:,first:last),v(first:last))/(l+2)
   endif
  enddo
  call near_coupling_map(r,ei,hi,potential,Ac)
  code=0
  if(.not.ieee_is_finite(potential))code=4
  do k=1,3
   if(.not.ieee_is_finite(Ac(k)))code=4
  enddo
  if(present(status))then
   status=code
  elseif(code/=0)then
   error stop 'SALMONpy nonfinite near-field coupling'
  endif
 end subroutine
 real(c_double) function near_line_integral(start,finish,geometry,lmax,ch,decay_vector) result(integral)
  real(c_double),intent(in)::start(3),finish(3),geometry(4),ch(3,*)
  integer,intent(in)::lmax
  real(c_double),optional,intent(in)::decay_vector(3)
  real(c_double),parameter::nodes(5)=[-.906179845938663993d0,-.538469310105683091d0,0d0, &
                                      .538469310105683091d0,.906179845938663993d0]
  real(c_double),parameter::weights(5)=[.236926885056189088d0,.478628670499366468d0,.568888888888888889d0, &
                                        .478628670499366468d0,.236926885056189088d0]
  real(c_double),parameter::x16(8)=[.095012509837637440d0,.281603550779258913d0, &
    .458016777657227386d0,.617876244402643748d0,.755404408355003034d0, &
    .865631202387831744d0,.944575023073232576d0,.989400934991649933d0]
  real(c_double),parameter::w16(8)=[.189450610455068496d0,.182603415044923589d0, &
    .169156519395002538d0,.149595988816576732d0,.124628971255533872d0, &
    .095158511682492785d0,.062253523938647893d0,.027152459411754095d0]
  real(c_double)::r(3),Ac(3),potential,zero(3,165),parameter,argument
  integer::k,segment,segments,side
  logical::weighted
  zero=0d0;integral=0d0;weighted=.false.
  if(present(decay_vector))weighted=any(decay_vector/=0d0)
  if(weighted)then
   argument=abs(dot_product(decay_vector,finish-start))
   if(.not.ieee_is_finite(argument)) &
    error stop 'SALMONpy magnetic line is outside supported quadrature range'
   if(argument/2d0>real(huge(segments)-1,c_double)) &
    error stop 'SALMONpy magnetic quadrature segment count overflows'
   segments=max(1,ceiling(argument/2d0))
   do segment=1,segments
   do k=1,8
   do side=-1,1,2
    parameter=(segment-.5d0+.5d0*side*x16(k))/segments
    r=start+parameter*(finish-start)-geometry(1:3)
    call near_potentials(r,lmax,geometry(4),zero,ch,potential,Ac,decay_vector)
    integral=integral+.5d0*w16(k)*dot_product(Ac,finish-start)/segments
   enddo
   enddo
   enddo
  else
   do k=1,5
    r=start+.5d0*(nodes(k)+1d0)*(finish-start)-geometry(1:3)
    call near_potentials(r,lmax,geometry(4),zero,ch,potential,Ac)
    integral=integral+.5d0*weights(k)*dot_product(Ac,finish-start)
   enddo
  endif
 end function
 integer(c_int) function sp_near_protocol() bind(C)
  sp_near_protocol=3
 end function
 integer(c_int) function sp_near_fit(n,lmax,positions,fields,geometry,auto_center,auto_radius,ce,ch,statistics) &
  bind(C) result(status)
  integer(c_int),value::n,lmax,auto_center,auto_radius
  real(c_double),intent(in)::positions(3,n),fields(3,n,2)
  real(c_double),intent(inout)::geometry(7),ce(3,*),ch(3,*),statistics(2)
  real(c_double),allocatable::a(:,:),qmat(:,:),rhs(:,:)
  real(c_double)::rmat(165,165),b(165,6),coeff(165,6),v(165),geo(7),stats(2),norm,scale,projection,denom,err,envelope
  integer::i,j,k,pass,p,which,d,allocation_status
  status=1
  if(n<1.or.lmax<0.or.lmax>solid_lmax)return
  p=(lmax+1)*(lmax+2)*(lmax+3)/6
  if(n<p)return
  do i=1,n
   do d=1,3
    if(.not.ieee_is_finite(positions(d,i)))return
    do which=1,2
     if(.not.ieee_is_finite(fields(d,i,which)))return
    enddo
   enddo
  enddo
  geo=geometry
  if(auto_center==1)then
   geo(1:3)=0d0
   do i=1,n
    geo(1:3)=geo(1:3)+positions(:,i)/n
   enddo
  endif
  do d=1,7
   if(.not.ieee_is_finite(geo(d)))return
  enddo
  if(auto_radius==1)then
   geo(4)=0d0
   do i=1,n
    geo(4)=max(geo(4),sqrt(sum((positions(:,i)-geo(1:3))**2)))
   enddo
   if(geo(4)==0d0.and.lmax==0)geo(4)=1d0
  endif
  if(.not.ieee_is_finite(geo(4)))return
  if(geo(4)<=0d0)return
  allocate(a(n,p),qmat(n,p),rhs(n,6),stat=allocation_status);status=2
  if(allocation_status/=0)return
  status=4
  do i=1,n
   call near_basis((positions(:,i)-geo(1:3))/geo(4),lmax,v)
   do j=1,p
    if(.not.ieee_is_finite(v(j)))return
   enddo
   envelope=-dot_product(geo(5:7),positions(:,i)-geo(1:3))
   if(envelope>log(huge(1d0)))return
   a(i,:)=v(:p)*exp(envelope)
   rhs(i,1:3)=fields(:,i,1);rhs(i,4:6)=fields(:,i,2)
  enddo
  rmat=0d0;b=0d0;coeff=0d0;qmat=a;scale=maxval(abs(a))
  if(scale==0d0)then
   status=1;return
  endif
  do j=1,p
   do pass=1,2
    do k=1,j-1
     projection=dot_product(qmat(:,k),qmat(:,j))
     rmat(k,j)=rmat(k,j)+projection
     qmat(:,j)=qmat(:,j)-projection*qmat(:,k)
    enddo
   enddo
   norm=sqrt(dot_product(qmat(:,j),qmat(:,j)));status=1
   if(norm<1d-12*scale)return
   status=4
   if(.not.ieee_is_finite(norm))return
   rmat(j,j)=norm;qmat(:,j)=qmat(:,j)/norm
   do k=1,6
    b(j,k)=dot_product(qmat(:,j),rhs(:,k))
   enddo
  enddo
  do k=1,6
   do j=p,1,-1
    coeff(j,k)=(b(j,k)-dot_product(rmat(j,j+1:p),coeff(j+1:p,k)))/rmat(j,j)
    if(.not.ieee_is_finite(coeff(j,k)))return
   enddo
  enddo
  do which=1,2
   denom=0d0;err=0d0
   do k=1,3
    denom=denom+sum(rhs(:,3*(which-1)+k)**2)
    rhs(:,3*(which-1)+k)=matmul(a,coeff(:p,3*(which-1)+k))-rhs(:,3*(which-1)+k)
    err=err+sum(rhs(:,3*(which-1)+k)**2)
   enddo
   stats(which)=sqrt(err/max(denom,tiny(1d0)))
   if(.not.ieee_is_finite(stats(which)))return
  enddo
  geometry=geo;ce(:,1:p)=transpose(coeff(:p,1:3));ch(:,1:p)=transpose(coeff(:p,4:6));statistics=stats
  status=0
 end function
 integer(c_int) function sp_near_evaluate(n,lmax,positions,geometry,ce,ch,E,H,potential,Ac) bind(C) result(status)
  integer(c_int),value::n,lmax
  real(c_double),intent(in)::positions(3,n),geometry(7),ce(3,*),ch(3,*)
  real(c_double),intent(inout)::E(3,n),H(3,n),potential(n),Ac(3,n)
  real(c_double),allocatable::out(:,:)
  real(c_double)::v(165),geo(7),envelope
  integer::i,j,k,p,allocation_status
  status=1
  if(n<1.or.lmax<0.or.lmax>solid_lmax)return
  p=(lmax+1)*(lmax+2)*(lmax+3)/6;geo=geometry
  do k=1,7
   if(.not.ieee_is_finite(geo(k)))return
  enddo
  if(geo(4)<=0d0)return
  do i=1,p
   do k=1,3
    if(.not.ieee_is_finite(ce(k,i)).or..not.ieee_is_finite(ch(k,i)))return
   enddo
  enddo
  allocate(out(10,n),stat=allocation_status);status=2
  if(allocation_status/=0)return
  status=1
  do i=1,n
   do k=1,3
    if(.not.ieee_is_finite(positions(k,i)))return
   enddo
   call near_basis((positions(:,i)-geo(1:3))/geo(4),lmax,v)
   envelope=-dot_product(geo(5:7),positions(:,i)-geo(1:3))
   status=4
   if(envelope>log(huge(1d0)))return
   out(1:3,i)=matmul(ce(:,1:p),v(:p))*exp(envelope)
   out(4:6,i)=matmul(ch(:,1:p),v(:p))*exp(envelope)
   call near_potentials(positions(:,i)-geo(1:3),lmax,geo(4),ce,ch,out(7,i),out(8:10,i),geo(5:7),status)
   if(status/=0)return
   status=4
   do j=1,10
    if(.not.ieee_is_finite(out(j,i)))return
   enddo
  enddo
  E=out(1:3,:);H=out(4:6,:);potential=out(7,:);Ac=out(8:10,:);status=0
 end function
end module
