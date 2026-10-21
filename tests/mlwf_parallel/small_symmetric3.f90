module small_symmetric3
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
private
public::symmetric3_eigen
contains
pure recursive subroutine symmetric3_eigen(input,v,eval,status)
real(8),intent(in)::input(3,3)
real(8),intent(out)::v(3,3),eval(3)
integer,intent(out)::status
real(8)::a(3,3),scale,off,tau,t,c,s,app,aqq,apq,x,y,tmp
integer::i,j,k,p,q,it
status=1;v=0d0;eval=0d0
 do i=1,3;do j=1,3
 if(.not.ieee_is_finite(input(i,j)))return
 enddo;enddo
scale=maxval(abs(input))
do i=1,3;v(i,i)=1d0;enddo
if(scale==0d0)then
 status=0;return
endif
a=input/scale
if(maxval(abs(a-transpose(a)))>64*epsilon(1d0))return
a=(a+transpose(a))*.5d0
 do it=1,64
 p=1;q=2;off=abs(a(1,2))
 if(abs(a(1,3))>off)then;p=1;q=3;off=abs(a(1,3));endif
 if(abs(a(2,3))>off)then;p=2;q=3;off=abs(a(2,3));endif
 if(off<=2*epsilon(1d0))exit
 app=a(p,p);aqq=a(q,q);apq=a(p,q)
 tau=(aqq-app)/(2*apq)
 t=sign(1d0,tau)/(abs(tau)+sqrt(1d0+tau*tau))
 c=1d0/sqrt(1d0+t*t);s=t*c
 a(p,p)=app-t*apq;a(q,q)=aqq+t*apq;a(p,q)=0d0;a(q,p)=0d0
 do k=1,3
 if(k==p.or.k==q)cycle
 x=a(k,p);y=a(k,q)
 a(k,p)=c*x-s*y;a(p,k)=a(k,p)
 a(k,q)=s*x+c*y;a(q,k)=a(k,q)
 enddo
 do k=1,3
 x=v(k,p);y=v(k,q);v(k,p)=c*x-s*y;v(k,q)=s*x+c*y
 enddo
 enddo
if(it>64)return
 do i=1,3;eval(i)=a(i,i)*scale;enddo
 do i=1,2;do j=i+1,3
 if(eval(j)>=eval(i))cycle
 tmp=eval(i);eval(i)=eval(j);eval(j)=tmp
 do k=1,3;tmp=v(k,i);v(k,i)=v(k,j);v(k,j)=tmp;enddo
 enddo;enddo
status=0
end subroutine
end module
