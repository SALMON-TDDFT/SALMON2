! Read-only local Yee energy/flux accounting. No MPI, state or I/O.
module salmonpy_fdtd_energy
use iso_c_binding
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
private
public::sp_fdtd_box_energy_conductive,sp_fdtd_cross_energy,sp_fdtd_box_energy,sp_fdtd_box_energy_diagonal,sp_fdtd_box_summary
type real_buffer
 real(c_double),pointer::p(:)=>null()
end type
contains
! K is the symmetric off-diagonal constitutive interpolation operator.
! Caller supplies E_old,E_new,K E_old,K E_new on the same edge layout.
! Returns cross electric energy old/new and outward constitutive exchange.
! This exchange is NOT E cross H flux. No communication or state mutation.
integer(c_int) function sp_fdtd_cross_energy(lo,hi,owned_lo,owned_hi,box,dx,old,new,ko,kn,result) &
 bind(C) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),intent(in)::old(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
 real(c_double),intent(in)::new(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
 real(c_double),intent(in)::ko(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
 real(c_double),intent(in)::kn(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
 real(c_double),intent(inout)::result(3)
 real(c_double)::s(3),scale,e0,e1,k0,k1
 integer::a(3),b(3),d,i,j,k
 status=1;s=0d0
 do d=1,3
  if(hi(d)<lo(d).or.owned_lo(d)<lo(d).or.owned_hi(d)>hi(d))return
  if(box(d,1)>box(d,2).or.owned_lo(d)>owned_hi(d))return
  if(.not.ieee_is_finite(dx(d)))return
  if(dx(d)<=0d0)return
 enddo
 a=max(lo,owned_lo,box(:,1));b=min(hi,owned_hi,box(:,2))
 scale=product(dx)/(8d0*acos(-1d0))
 do k=a(3),b(3)
 do j=a(2),b(2)
 do i=a(1),b(1)
 do d=1,3
  e0=old(d,i,j,k);e1=new(d,i,j,k);k0=ko(d,i,j,k);k1=kn(d,i,j,k)
  if(.not.ieee_is_finite(e0).or..not.ieee_is_finite(e1))return
  if(.not.ieee_is_finite(k0).or..not.ieee_is_finite(k1))return
  s(1)=s(1)+e0*k0;s(2)=s(2)+e1*k1
  s(3)=s(3)+e0*k1-e1*k0
 enddo
 enddo
 enddo
 enddo
 s=scale*s
 do d=1,3
  if(.not.ieee_is_finite(s(d)))return
 enddo
 result=s;status=0
end function sp_fdtd_cross_energy

integer(c_int) function sp_fdtd_box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles, &
 oldptr,newptr,jptr,sourceptr,result) bind(C) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),value::dt,speed
 integer(c_int),value::poles
 type(c_ptr),intent(in)::oldptr(12),newptr(12),jptr(6)
 type(c_ptr),value::sourceptr
 real(c_double),intent(inout)::result(17)
 status=box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles,oldptr,newptr,jptr,sourceptr,result)
end function

! Same17 observables; diagonal positive constitutive mass from actual coefficients.
! coeffptr: c1Ex,c1Ey,c1Ez,cJx,cJy,cJz,cHx_y,cHy_z,cHz_x. Sigma must be zero.
integer(c_int) function sp_fdtd_box_energy_diagonal(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles, &
 oldptr,newptr,jptr,sourceptr,coeffptr,result) bind(C) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),value::dt,speed
 integer(c_int),value::poles
 type(c_ptr),intent(in)::oldptr(12),newptr(12),jptr(6),coeffptr(9)
 type(c_ptr),value::sourceptr
 real(c_double),intent(inout)::result(17)
 status=box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles,oldptr,newptr,jptr,sourceptr,result,coeffptr)
end function

! Actual electric loss minus the implicit LD coefficients gives Ohmic sigma.
! coeffptr(1:9): as diagonal kernel; (10:12): LD c2 for x/y/z, pole-major.
integer(c_int) function sp_fdtd_box_energy_conductive(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles, &
 oldptr,newptr,jptr,sourceptr,coeffptr,result) bind(C) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),value::dt,speed
 integer(c_int),value::poles
 type(c_ptr),intent(in)::oldptr(12),newptr(12),jptr(6),coeffptr(12)
 type(c_ptr),value::sourceptr
 real(c_double),intent(inout)::result(18)
 real(c_double)::base(17),heat
 base=0d0;heat=0d0
 status=box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles,oldptr,newptr,jptr,sourceptr,base,coeffptr,heat)
 if(status/=0)return
 base(16)=base(16)+heat;base(17)=base(17)+heat
 if(.not.ieee_is_finite(base(16)).or..not.ieee_is_finite(base(17)))then
  status=4;return
 endif
 result(1:17)=base;result(18)=heat
end function

integer(c_int) function box_energy(lo,hi,owned_lo,owned_hi,box,dx,dt,speed,poles, &
 oldptr,newptr,jptr,sourceptr,result,coeffptr,ohmic) result(status)
 integer(c_int),intent(in)::lo(3),hi(3),owned_lo(3),owned_hi(3),box(3,2)
 real(c_double),intent(in)::dx(3)
 real(c_double),value::dt,speed
 integer(c_int),value::poles
 type(c_ptr),intent(in)::oldptr(12),newptr(12),jptr(6)
 type(c_ptr),value::sourceptr
 real(c_double),intent(inout)::result(17)
 type(c_ptr),optional,intent(in)::coeffptr(:)
 real(c_double),optional,intent(out)::ohmic
 type(real_buffer)::old(12),new(12),curr(6),coeff(12)
 real(c_double),pointer::source(:)
 integer(c_int64_t)::dims(3),count,index,jindex
 integer::a(3),b(3),i,j,k,d,q,v(3),plus(3)
 real(c_double)::s(17),eo(3),en(3),ho(3),hn(3),hp(3),em(3),curl(3),xp(3),yp(3),zp(3)
 real(c_double)::ldj(3),ext(3),vol,pi,factor,epsilon(3),mu(3),denominator,component
 real(c_double)::heat_value,sigma,loss,ld_loss,c2,roundoff
 integer::derivative(3)
 logical::valid
 status=1;heat_value=0d0
 if(present(ohmic).and..not.present(coeffptr))return
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
  if(present(ohmic))ohmic=0d0
  result=s;status=0;return
 endif
 if(any(a<lo+1).or.any(b>hi-1))return
 do d=1,12
  if(.not.c_associated(oldptr(d)).or..not.c_associated(newptr(d)))return
  call c_f_pointer(oldptr(d),old(d)%p,[count])
  call c_f_pointer(newptr(d),new(d)%p,[count])
 enddo
 if(present(coeffptr))then
  if(size(coeffptr)<9)return
  do d=1,9
   if(.not.c_associated(coeffptr(d)))return
   call c_f_pointer(coeffptr(d),coeff(d)%p,[count])
  enddo
 endif
 if(present(ohmic).and.poles>0)then
  if(size(coeffptr)<12)return
  do d=10,12
   if(.not.c_associated(coeffptr(d)))return
   call c_f_pointer(coeffptr(d),coeff(d)%p,[count*poles])
  enddo
 endif
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
  if(present(coeffptr))then
   ! c4=-4*pi*dt/cJ; epsilon=(c4+c5)/2. The LD implicit term cancels.
   ! H uses the actual staggered coefficient; no material-ID reconstruction.
   derivative=[2,3,1]
   do d=1,3
    component=coeff(d)%p(index);denominator=coeff(d+3)%p(index)
    if(.not.ieee_is_finite(component).or..not.ieee_is_finite(denominator))then
     status=4;return
    endif
    if(denominator>=0d0)return
    epsilon(d)=(-4d0*pi*dt/denominator)*(.5d0+.5d0*component)
    denominator=coeff(d+6)%p(index)
    if(.not.ieee_is_finite(denominator))then
     status=4;return
    endif
    if(denominator>=0d0)return
    mu(d)=-speed*dt/(denominator*dx(derivative(d)))
    if(.not.ieee_is_finite(epsilon(d)).or..not.ieee_is_finite(mu(d)))then
     status=4;return
    endif
    if(epsilon(d)<=0d0.or.mu(d)<=0d0)return
   enddo
   hp=ho+speed*dt*curl/mu
   s(1)=s(1)+factor*(dot_product(epsilon*eo,eo)+dot_product(mu*hp,ho))
   s(2)=s(2)+factor*(dot_product(epsilon*en,en)+dot_product(mu*ho,hn))
   s(4)=s(4)+factor*(dot_product(epsilon*eo,eo)+dot_product(mu*(.5d0*(hp+ho)),.5d0*(hp+ho)))
   s(5)=s(5)+factor*(dot_product(epsilon*en,en)+dot_product(mu*(.5d0*(ho+hn)),.5d0*(ho+hn)))
  else
  hp=ho+speed*dt*curl
  s(1)=s(1)+factor*(dot_product(eo,eo)+dot_product(hp,ho))
  s(2)=s(2)+factor*(dot_product(en,en)+dot_product(ho,hn))
  s(4)=s(4)+factor*(dot_product(eo,eo)+dot_product(.5d0*(hp+ho),.5d0*(hp+ho)))
  s(5)=s(5)+factor*(dot_product(en,en)+dot_product(.5d0*(ho+hn),.5d0*(ho+hn)))
  endif
  em=.5d0*eo+.5d0*en;ldj=0d0;ext=0d0
  if(present(ohmic))then
   do d=1,3
    loss=-(1d0-coeff(d)%p(index))/coeff(d+3)%p(index)
    ld_loss=0d0
    do q=1,poles
     c2=coeff(d+9)%p(index+count*(q-1))
     if(.not.ieee_is_finite(c2))then
      status=4;return
     endif
     if(c2<0d0)return
     ld_loss=ld_loss+c2
    enddo
    sigma=loss-ld_loss
    if(.not.ieee_is_finite(sigma))then
     status=4;return
    endif
    ! Remove only subtraction roundoff, not a physical negative conductivity.
    roundoff=64d0*spacing(1d0)*max(abs(loss),abs(ld_loss))
    if(sigma < -roundoff)return
    if(abs(sigma)<=roundoff)sigma=0d0
    heat_value=heat_value+vol*dt*sigma*em(d)*em(d)
   enddo
  endif
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
 if(.not.ieee_is_finite(heat_value))return
 if(present(ohmic))ohmic=heat_value
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
 real(c_double),intent(in)::steps(19,n)
 real(c_double),value::quantum,own_work,point_work
 real(c_double),intent(inout)::result(12)
 real(c_double)::s(12)
 integer::i,j
 status=1
 if(n<1)return
 if(.not.ieee_is_finite(quantum).or..not.ieee_is_finite(own_work).or..not.ieee_is_finite(point_work))return
 s=0d0
 do i=1,n
  do j=1,19
   if(.not.ieee_is_finite(steps(j,i)))return
  enddo
  s(1)=s(1)+steps(3,i)
  s(2)=s(2)+steps(7,i)
  s(3)=s(3)+steps(8,i)
  s(4)=s(4)+steps(15,i)
  s(11)=s(11)+steps(18,i)
  s(12)=s(12)+steps(19,i)
  s(5)=s(5)+steps(16,i)
  s(6)=max(s(6),abs(steps(16,i)))
  s(9)=s(9)+steps(6,i)
 enddo
 s(7)=s(1)+s(2)+s(4)+s(11)+s(12)+quantum+own_work
 s(8)=s(3)-point_work
 s(10)=s(9)+s(2)+s(4)+s(11)+s(12)+quantum+own_work
 do i=1,12
  if(.not.ieee_is_finite(s(i)))return
 enddo
 result=s;status=0
end function sp_fdtd_box_summary
end module
