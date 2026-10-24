! Auxiliary Pauli core model built from BOTH occupied Dirac blocks.
! Not a Foldy-Wouthuysen reduction and not the Dirac kinetic-energy density.
module crps_core_aux
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: build_core_aux
contains
 subroutine hermite(x,y,dy,r,s,ds,dds)
  real(8),intent(in)::x(0:),y(0:),dy(0:),r
  real(8),intent(out)::s,ds,dds
  integer::lo,hi,k
  real(8)::h,t,a,b,c,d
  lo=0;hi=ubound(x,1)
  do while(hi-lo>1)
   k=(lo+hi)/2
   if(x(k)>r)then
    hi=k
   else
    lo=k
   endif
  enddo
  h=x(hi)-x(lo);t=(r-x(lo))/h
  a=y(lo);b=h*dy(lo);c=3*(y(hi)-y(lo))-h*(2*dy(lo)+dy(hi));d=2*(y(lo)-y(hi))+h*(dy(lo)+dy(hi))
  s=a+t*(b+t*(c+t*d));ds=(b+t*(2*c+3*t*d))/h;dds=(2*c+6*t*d)/h**2
 end subroutine

 subroutine build_core_aux(n,h,ipc,x,fa,ga,va,eps,occ,amap,la,ka,za,alpha,hc,mc2,rho,nb,ids,ls,lambda,q,amp,damp,tau,cur)
  integer,intent(in)::n,ipc,amap(:),la(:),ka(:)
  real(8),intent(in)::h,x(0:),fa(0:,:),ga(0:,:),va(0:),eps(:),occ(:),za,alpha,hc,mc2,rho(0:n)
  integer,intent(out)::nb,ids(:),ls(:),lambda(:)
  real(8),intent(out)::q(:),amp(0:,:),damp(0:,:),tau(0:n),cur(0:n)
  real(8)::raw(0:ubound(x,1)),dr(0:ubound(x,1)),seed(0:n,size(q)),dseed(0:n,size(q)),density_d(0:n)
  real(8)::pi,rc,r,s,ds,dds,a,b,c,e,t,p,dp,pp,z,fp,gp,vf,scale,den,num,fac,logw,ar
  integer::m,block,l,k,i,j,nx
  logical::tail_started
  real(8)::slope,tail_charge
  integer::tail_first,anchor
  real(8)::decay,affected_charge
  pi=acos(-1d0);nx=ubound(x,1);rc=ipc*h
  if(ipc<2)error stop 'auxiliary core requires nonzero NLCC radius'
  nb=0;ids=0;ls=0;lambda=0;q=0d0;amp=0d0;damp=0d0;tau=0d0;cur=0d0;seed=0d0;dseed=0d0
  do m=1,size(occ)
   if(amap(m)/=0.or.occ(m)==0d0)cycle
   do block=1,2
    nb=nb+1
    if(nb>size(q))error stop 'auxiliary core block capacity'
    ids(nb)=m;q(nb)=occ(m)
    if(block==1)then
     l=la(m);lambda(nb)=-(ka(m)+1)
    else
     l=la(m)-sign(1,ka(m));lambda(nb)=ka(m)-1
    endif
    ! kappa>0: small l=l-1; kappa<0: small l=l+1.
    ls(nb)=l
    do k=0,nx
     r=x(k);vf=va(k)
     fp=-ka(m)*fa(k,m)/r+(za*alpha/r+(2*mc2+eps(m)-vf)/hc)*ga(k,m)
     gp=ka(m)*ga(k,m)/r-(za*alpha/r+(eps(m)-vf)/hc)*fa(k,m)
     if(block==1)then
      raw(k)=fa(k,m)/r;dr(k)=fp/r-fa(k,m)/r**2
     else
      raw(k)=ga(k,m)/r;dr(k)=gp/r-ga(k,m)/r**2
     endif
    enddo
    call hermite(x,raw,dr,rc,s,ds,dds)
    ! s(t)=t**l*(a+b*t**2+c*t**4), C2 matching at t=1.
    e=(rc*ds-l*s)/2;c=(rc**2*dds-l*(l-1)*s-(4*l+2)*e)/8
    b=e-2*c;a=s-b-c
    do i=0,n
     r=i*h
     if(i<=ipc)then
      t=r/rc;p=a+b*t*t+c*t**4;dp=2*b*t+4*c*t**3
      if(l==0)then
       seed(i,nb)=p;dseed(i,nb)=dp/rc
      else
       seed(i,nb)=t**l*p;dseed(i,nb)=(l*t**(l-1)*p+t**l*dp)/rc
      endif
     else
      call hermite(x,raw,dr,r,seed(i,nb),dseed(i,nb),dds)
     endif
    enddo
   enddo
  enddo
  if(nb==0)error stop 'empty auxiliary core'
  density_d=0d0
  do i=2,n-2
   density_d(i)=(rho(i-2)-8*rho(i-1)+8*rho(i+1)-rho(i+2))/(12*h)
  enddo
  density_d(1)=(rho(1)-8*rho(0)+8*rho(2)-rho(3))/(12*h)
  density_d(n-1)=(rho(n)-rho(n-2))/(2*h)
  density_d(n)=(3*rho(n)-4*rho(n-1)+rho(n-2))/(2*h)
  tail_started=.false.;tail_first=-1;tail_charge=0d0
  do i=0,n
   if(.not.ieee_is_finite(rho(i)).or.rho(i)<0d0)error stop 'invalid auxiliary target density'
   if(rho(i)==0d0)cycle
   scale=maxval(abs(seed(i,1:nb)))
   if(scale==0d0.or.tail_started)then
    if(i*h<=10d0)error stop 'positive stable target core without seed support'
    if(.not.tail_started)then
     tail_first=i
     affected_charge=0d0
     do k=i,n
      affected_charge=affected_charge+4*pi*(k*h)**2*rho(k)*h
     enddo
     if(affected_charge>=1d-60)error stop 'nonnegligible core seed-support loss'
     anchor=i-1
     do while(anchor>0)
      if(maxval(abs(seed(anchor,1:nb)))>0d0)exit
      anchor=anchor-1
     enddo
     if(anchor==0)error stop 'auxiliary tail continuation lacks native anchor'
     write(*,*) 'AUX_TAIL_CONTINUATION_FIRST ',i,i*h,rho(i),anchor*h
     write(*,*) 'AUX_TAIL_AFFECTED_CHARGE_RECT ',affected_charge
    endif
    tail_started=.true.
    scale=maxval(abs(seed(anchor,1:nb)))
    do j=1,nb
     ! Bound-Dirac asymptotic F,G ~ exp(-decay*r); angular amplitudes are F/r,G/r.
     decay=sqrt(-eps(ids(j))*(2*mc2+eps(ids(j))))/hc
     seed(i,j)=(seed(anchor,j)/scale)*(anchor*dble(h)/(i*h))*exp(-decay*(i-anchor)*h)
     dseed(i,j)=(-decay-1d0/(i*h))*seed(i,j)
    enddo
    scale=maxval(abs(seed(i,1:nb)))
    tail_charge=tail_charge+4*pi*(i*h)**2*rho(i)*h
   endif
   den=0d0;num=0d0
   do j=1,nb
    z=seed(i,j)/scale;den=den+q(j)*z*z
    num=num+q(j)*z*(dseed(i,j)/scale)
   enddo
   fac=sqrt(4*pi*rho(i)/den)
   logw=.5d0*density_d(i)/rho(i)-num/den
   do j=1,nb
    amp(i,j)=fac*(seed(i,j)/scale)
    damp(i,j)=fac*(dseed(i,j)/scale)+logw*amp(i,j)
    if(.not.ieee_is_finite(amp(i,j)).or..not.ieee_is_finite(damp(i,j)))error stop 'nonfinite auxiliary core primitive'
    l=ls(j)
    if(i==0)then
     ar=0d0
     if(l==1)ar=damp(i,j)
    else
     ar=amp(i,j)/(i*h)
    endif
    tau(i)=tau(i)+q(j)*(damp(i,j)**2+l*(l+1)*ar**2)/(4*pi)
    if(i>0)cur(i)=cur(i)+q(j)*lambda(j)*amp(i,j)**2/(8*pi*i*h)
   enddo
   if(.not.ieee_is_finite(tau(i)).or..not.ieee_is_finite(cur(i)))error stop 'nonfinite auxiliary core closure'
  enddo
  if(tail_started)write(*,*)'AUX_TAIL_CHARGE_UPPER_RECT ',tail_first,tail_charge
 end subroutine
end module
