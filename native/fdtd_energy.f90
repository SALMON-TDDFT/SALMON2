! Read-only local Yee energy/flux accounting. No MPI, state or I/O.
module salmonpy_fdtd_energy
use iso_c_binding
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
private
public::sp_fdtd_box_energy,sp_fdtd_box_summary
type real_buffer
 real(c_double),pointer::p(:)=>null()
end type
contains
integer(c_int) function sp_fdtd_box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles, &
 oldptr,newptr,jptr,sourceptr,result) bind(C) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),value::dt,speed
 integer(c_int),value::poles
 type(c_ptr),intent(in)::oldptr(12),newptr(12),jptr(6)
 type(c_ptr),value::sourceptr
 real(c_double),intent(inout)::result(17)
 type(real_buffer)::old(12),new(12),curr(6)
 real(c_double),pointer::source(:)
 integer(c_int64_t)::dims(3),count,index,jindex
 integer::a(3),b(3),i,j,k,d,q,v(3),plus(3)
 real(c_double)::s(17),eo(3),en(3),ho(3),hn(3),hp(3),em(3),curl(3),xp(3),yp(3),zp(3)
 real(c_double)::ldj(3),ext(3),vol,pi,factor
 logical::valid
 status=1
 if(poles<0.or..not.ieee_is_finite(dt).or..not.ieee_is_finite(speed))return
 if(dt<=0d0.or.speed<=0d0)return
 count=1
 do d=1,3
  dims(d)=int(hi(d),c_int64_t)-int(lo(d),c_int64_t)+1
  if(dims(d)<3.or.box(d,1)>box(d,2))return
  if(.not.ieee_is_finite(dx(d)))return
  if(dx(d)<=0d0)return
  if(count>huge(count)/dims(d))return
  count=count*dims(d)
 enddo
 if(poles>0)then
  if(count>huge(count)/poles)return
 endif
 a=max(owned_lo,box(:,1));b=min(owned_hi,box(:,2));s=0d0
 if(any(a>b))then
  result=s;status=0;return
 endif
 if(any(a<lo+1).or.any(b>hi-1))return
 do d=1,12
  if(.not.c_associated(oldptr(d)).or..not.c_associated(newptr(d)))return
  call c_f_pointer(oldptr(d),old(d)%p,[count])
  call c_f_pointer(newptr(d),new(d)%p,[count])
 enddo
 if(poles>0)then
  do d=1,6
   if(.not.c_associated(jptr(d)))return
   call c_f_pointer(jptr(d),curr(d)%p,[count*poles])
  enddo
 endif
 if(c_associated(sourceptr))then
  if(count>huge(count)/3)return
  call c_f_pointer(sourceptr,source,[3*count])
 else
  nullify(source)
 endif
 valid=.true.;vol=product(dx);pi=acos(-1d0);factor=vol/(8d0*pi)
 do k=a(3),b(3)
 do j=a(2),b(2)
 do i=a(1),b(1)
  v=[i,j,k];index=offset(v)
  eo=field(old,0,v);en=field(new,0,v);ho=field(old,6,v);hn=field(new,6,v)
  plus=v;plus(1)=plus(1)+1;xp=field(old,0,plus)
  plus=v;plus(2)=plus(2)+1;yp=field(old,0,plus)
  plus=v;plus(3)=plus(3)+1;zp=field(old,0,plus)
  curl=[(yp(3)-eo(3))/dx(2)-(zp(2)-eo(2))/dx(3), &
        (zp(1)-eo(1))/dx(3)-(xp(3)-eo(3))/dx(1), &
        (xp(2)-eo(2))/dx(1)-(yp(1)-eo(1))/dx(2)]
  hp=ho+speed*dt*curl
  s(1)=s(1)+factor*(dot_product(eo,eo)+dot_product(hp,ho))
  s(2)=s(2)+factor*(dot_product(en,en)+dot_product(ho,hn))
  s(4)=s(4)+factor*(dot_product(eo,eo)+dot_product(.5d0*(hp+ho),.5d0*(hp+ho)))
  s(5)=s(5)+factor*(dot_product(en,en)+dot_product(.5d0*(ho+hn),.5d0*(ho+hn)))
  em=.5d0*eo+.5d0*en;ldj=0d0;ext=0d0
  do q=1,poles
   jindex=index+count*(q-1)
   do d=1,3
    if(.not.ieee_is_finite(curr(d)%p(jindex)).or..not.ieee_is_finite(curr(d+3)%p(jindex)))valid=.false.
    ldj(d)=ldj(d)+.5d0*curr(d)%p(jindex)+.5d0*curr(d+3)%p(jindex)
   enddo
  enddo
  if(associated(source))then
   do d=1,3
    ext(d)=source(3*(index-1)+d)
    if(.not.ieee_is_finite(ext(d)))valid=.false.
   enddo
  endif
  s(7)=s(7)+vol*dt*dot_product(em,ldj)
  s(8)=s(8)+vol*dt*dot_product(em,ext)
  do d=1,3
   if(v(d)==box(d,1))s(8+2*d-1)=s(8+2*d-1)+face(d,-1,v)
   if(v(d)==box(d,2))s(8+2*d)=s(8+2*d)+face(d,1,v)
  enddo
 enddo
 enddo
 enddo
 s(3)=s(2)-s(1);s(6)=s(5)-s(4);s(15)=sum(s(9:14))
 s(16)=s(3)+s(7)+s(8)+s(15);s(17)=s(6)+s(7)+s(8)+s(15)
 status=4
 if(.not.valid)return
 do d=1,17
  if(.not.ieee_is_finite(s(d)))return
 enddo
 result=s;status=0
contains
 integer(c_int64_t) function offset(v) result(n)
  integer,intent(in)::v(3)
  n=int(v(1),c_int64_t)-lo(1)+dims(1)*(int(v(2),c_int64_t)-lo(2)) &
   +dims(1)*dims(2)*(int(v(3),c_int64_t)-lo(3))+1
 end function
 function field(f,q,v) result(x)
  type(real_buffer),intent(in)::f(12)
  integer,intent(in)::q,v(3)
  real(c_double)::x(3)
  integer(c_int64_t)::n
  integer::d
  n=offset(v)
  do d=1,3
   x(d)=f(q+2*d-1)%p(n)+f(q+2*d)%p(n)
   if(.not.ieee_is_finite(x(d)))valid=.false.
  enddo
 end function
 real(c_double) function face(axis,side,v) result(value)
  integer,intent(in)::axis,side,v(3)
  integer::eindex(3),hindex(3)
  real(c_double)::electric(3),magnetic(3),cross(3)
  eindex=v;hindex=v
  if(side>0)eindex(axis)=eindex(axis)+1
  if(side<0)hindex(axis)=hindex(axis)-1
  electric=.5d0*field(old,0,eindex)+.5d0*field(new,0,eindex)
  magnetic=field(old,6,hindex)
  cross=[electric(2)*magnetic(3)-electric(3)*magnetic(2), &
         electric(3)*magnetic(1)-electric(1)*magnetic(3), &
         electric(1)*magnetic(2)-electric(2)*magnetic(1)]
  value=side*speed*dt*vol/dx(axis)/(4d0*pi)*cross(axis)
 end function
end function
! Sum accepted steps through a fresh quantum energy sample. Drift is not removed.
integer(c_int) function sp_fdtd_box_summary(n,steps,quantum,own_work,point_work,result) bind(C) result(status)
 integer(c_int),value::n
 real(c_double),intent(in)::steps(17,n)
 real(c_double),value::quantum,own_work,point_work
 real(c_double),intent(inout)::result(10)
 real(c_double)::s(10)
 integer::i,j
 status=1
 if(n<1)return
 if(.not.ieee_is_finite(quantum).or..not.ieee_is_finite(own_work).or..not.ieee_is_finite(point_work))return
 s=0d0
 do i=1,n
  do j=1,17
   if(.not.ieee_is_finite(steps(j,i)))return
  enddo
  s(1)=s(1)+steps(3,i)
  s(2)=s(2)+steps(7,i)
  s(3)=s(3)+steps(8,i)
  s(4)=s(4)+steps(15,i)
  s(5)=s(5)+steps(16,i)
  s(6)=max(s(6),abs(steps(16,i)))
  s(9)=s(9)+steps(6,i)
 enddo
 s(7)=s(1)+s(2)+s(4)+quantum+own_work
 s(8)=s(3)-point_work
 s(10)=s(9)+s(2)+s(4)+quantum+own_work
 do i=1,10
  if(.not.ieee_is_finite(s(i)))return
 enddo
 result=s;status=0
end function sp_fdtd_box_summary
end module
