! Electric-dipole source mapping and development midpoint/work kernels.
! The oscillator is a regression substitute, not a molecular TDDFT model.
module salmonpy_molecular_feedback
use iso_c_binding
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
contains
integer(c_int) function sp_molecular_protocol() bind(C)
 sp_molecular_protocol=1
end function
integer(c_int) function sp_molecular_midpoint(dt,omega,volume,p,v,total_left,total_right,own_left,own_right, &
                                             p_new,v_new,current,statistics) bind(C) result(status)
 real(c_double),value::dt,omega,volume
 real(c_double),intent(in)::p(3),v(3),total_left(3),total_right(3),own_left(3),own_right(3)
 real(c_double),intent(inout)::p_new(3),v_new(3),current(3),statistics(8)
 real(c_double)::q(3),w(3),j(3),dp(3),drive(3),s(8),factor
 integer::i
 status=1
 if(.not.ieee_is_finite(dt).or..not.ieee_is_finite(omega).or..not.ieee_is_finite(volume))return
 if(dt<=0d0.or.omega<0d0.or.volume<=0d0)return
 do i=1,3
  if(.not.ieee_is_finite(p(i)).or..not.ieee_is_finite(v(i)))return
  if(.not.ieee_is_finite(total_left(i)).or..not.ieee_is_finite(total_right(i)))return
  if(.not.ieee_is_finite(own_left(i)).or..not.ieee_is_finite(own_right(i)))return
 enddo
 status=4
 factor=(omega*dt/2d0)**2
 drive=.5d0*(total_left+total_right-own_left-own_right)
 q=((1d0-factor)*p+dt*v+.5d0*dt*dt*drive)/(1d0+factor)
 dp=q-p;w=2d0*dp/dt-v;j=dp/(volume*dt)
 s(1)=.5d0*(dot_product(v,v)+omega*omega*dot_product(p,p))
 s(2)=.5d0*(dot_product(w,w)+omega*omega*dot_product(q,q))
 s(3)=dot_product(.5d0*(total_left+total_right),dp)
 s(4)=dot_product(.5d0*(own_left+own_right),dp)
 s(5)=dot_product(drive,dp)
 s(6)=s(2)-s(1)-s(5)
 s(7)=s(3)-s(4)-s(5)
 ! The explicit held-left-work negative control exposes the old time split.
 s(8)=dot_product(total_left-own_left,dp)-s(5)
 do i=1,3
  if(.not.ieee_is_finite(q(i)).or..not.ieee_is_finite(w(i)).or..not.ieee_is_finite(j(i)))return
 enddo
 do i=1,8
  if(.not.ieee_is_finite(s(i)))return
 enddo
 p_new=q;v_new=w;current=j;statistics=s
 status=0
end function
integer(c_int) function sp_molecular_work(omega,p,v,q,w,total_left,total_right,own_left,own_right,statistics) &
 bind(C) result(status)
 real(c_double),value::omega
 real(c_double),intent(in)::p(3),v(3),q(3),w(3),total_left(3),total_right(3),own_left(3),own_right(3)
 real(c_double),intent(inout)::statistics(8)
 real(c_double)::drive(3),dp(3),s(8)
 integer::i
 status=1
 if(.not.ieee_is_finite(omega).or.omega<0d0)return
 do i=1,3
  if(.not.ieee_is_finite(p(i)).or..not.ieee_is_finite(v(i)))return
  if(.not.ieee_is_finite(q(i)).or..not.ieee_is_finite(w(i)))return
  if(.not.ieee_is_finite(total_left(i)).or..not.ieee_is_finite(total_right(i)))return
  if(.not.ieee_is_finite(own_left(i)).or..not.ieee_is_finite(own_right(i)))return
 enddo
 drive=.5d0*(total_left+total_right-own_left-own_right);dp=q-p
 s(1)=.5d0*(dot_product(v,v)+omega*omega*dot_product(p,p))
 s(2)=.5d0*(dot_product(w,w)+omega*omega*dot_product(q,q))
 s(3)=dot_product(.5d0*(total_left+total_right),dp)
 s(4)=dot_product(.5d0*(own_left+own_right),dp)
 s(5)=dot_product(drive,dp)
 s(6)=s(2)-s(1)-s(5);s(7)=s(3)-s(4)-s(5)
 s(8)=dot_product(total_left-own_left,dp)-s(5)
 status=4
 do i=1,8
  if(.not.ieee_is_finite(s(i)))return
 enddo
 statistics=s;status=0
end function
integer(c_int) function sp_molecular_summary(n,statistics,result) bind(C) result(status)
 integer(c_int),value::n
 real(c_double),intent(in)::statistics(8,n)
 real(c_double),intent(inout)::result(11)
 real(c_double)::s(11),scale
 integer::i,k
 status=1
 if(n<1)return
 s=0d0
 do k=1,n
 do i=1,8
  if(.not.ieee_is_finite(statistics(i,k)))return
 enddo
 s(1:3)=s(1:3)+statistics(3:5,k)
 s(4)=max(s(4),abs(statistics(6,k)))
 s(5)=s(5)+statistics(7,k);s(6)=s(6)+statistics(8,k)
 s(7)=max(s(7),abs(statistics(7,k)))
 scale=max(abs(statistics(1,k)),abs(statistics(2,k)),abs(statistics(5,k)),tiny(1d0))
 s(9)=max(s(9),abs(statistics(6,k))/scale)
 enddo
 s(8)=statistics(2,n)-statistics(1,1)-s(3)
 scale=max(abs(statistics(2,n)-statistics(1,1)),abs(s(3)),tiny(1d0))
 s(10)=abs(s(8))/scale;s(11)=abs(s(5))/max(abs(s(3)),tiny(1d0))
 status=4
 do i=1,11
  if(.not.ieee_is_finite(s(i)))return
 enddo
 result=s;status=0
end function
integer(c_int) function sp_molecular_select(total,own,selected) bind(C) result(status)
 real(c_double),intent(in)::total(3),own(3)
 real(c_double),intent(inout)::selected(3)
 real(c_double)::s(3)
 integer::i
 status=1
 do i=1,3
  if(.not.ieee_is_finite(total(i)).or..not.ieee_is_finite(own(i)))return
 enddo
 s=total-own;status=4
 do i=1,3
  if(.not.ieee_is_finite(s(i)))return
 enddo
 selected=s;status=0
end function
integer(c_int) function sp_molecular_map(mode,np,nm,dt,weights,volume,p,old,result) bind(C) result(status)
 integer(c_int),value::mode,np,nm
 real(c_double),value::dt
 ! C row-major: weights[point,molecule], p/old[molecule,xyz], result[point,xyz].
 real(c_double),intent(in)::weights(nm,np),volume(np),p(3,nm),old(3,nm)
 real(c_double),intent(inout)::result(3,np)
 real(c_double)::value,increment
 integer::i,m,d
 status=1
 if(np<1.or.nm<1.or.mode<0.or.mode>1)return
 if(.not.ieee_is_finite(dt).or.dt<=0d0)return
 do i=1,np
  if(.not.ieee_is_finite(volume(i)).or.volume(i)<=0d0)return
  do m=1,nm
   if(.not.ieee_is_finite(weights(m,i)).or.weights(m,i)<0d0)return
  enddo
 enddo
 do m=1,nm
 do d=1,3
  if(.not.ieee_is_finite(p(d,m)).or..not.ieee_is_finite(old(d,m)))return
 enddo
 enddo
 status=4
 do i=1,np
 do d=1,3
  value=0d0
  do m=1,nm
   increment=p(d,m)
   if(mode==1)increment=(p(d,m)-old(d,m))/dt
   value=value+weights(m,i)*increment
  enddo
  value=value/volume(i)
  if(.not.ieee_is_finite(value))return
  result(d,i)=value
 enddo
 enddo
 status=0
end function
! Work histories use actually consumed endpoint-held fields; midpoint is
! reported separately and never substituted for the electronic Hamiltonian.
integer(c_int) function sp_molecular_energy_work(n,fields,dipoles,energies,fresh,result) bind(C) result(status)
 integer(c_int),value::n
 real(c_double),intent(in)::fields(3,n),dipoles(3,n),energies(n),fresh(n)
 real(c_double),intent(inout)::result(6,n)
 real(c_double),allocatable::staged(:,:)
 real(c_double)::held,mid,delta(3)
 integer::i,k,ios
 status=1
 if(n<1)return
 do k=1,n
  if(.not.ieee_is_finite(energies(k)).or..not.ieee_is_finite(fresh(k)))return
  if(fresh(k)/=0d0.and.fresh(k)/=1d0)return
  do i=1,3
   if(.not.ieee_is_finite(fields(i,k)).or..not.ieee_is_finite(dipoles(i,k)))return
  enddo
 enddo
 if(fresh(1)/=1d0)return
 status=2
 allocate(staged(6,n),stat=ios)
 if(ios/=0)return
 staged=0d0;held=0d0;mid=0d0
 do k=1,n
  if(k>1)then
   delta=dipoles(:,k)-dipoles(:,k-1)
   held=held+dot_product(fields(:,k),delta)
   mid=mid+dot_product(.5d0*(fields(:,k)+fields(:,k-1)),delta)
  endif
  staged(1,k)=held;staged(2,k)=mid
  if(fresh(k)==1d0)then
   staged(3,k)=energies(k)-energies(1)
   staged(4,k)=staged(3,k)-held;staged(5,k)=staged(3,k)-mid
   staged(6,k)=1d0
  endif
 enddo
 status=4
 do k=1,n
 do i=1,6
  if(.not.ieee_is_finite(staged(i,k)))return
 enddo
 enddo
 result=staged;status=0
end function
integer(c_int) function sp_molecular_interval_field(left,right,selected) bind(C) result(status)
 real(c_double),intent(in)::left(3),right(3)
 real(c_double),intent(inout)::selected(3)
 real(c_double)::staged(3)
 integer::i
 status=1
 do i=1,3
  if(.not.ieee_is_finite(left(i)).or..not.ieee_is_finite(right(i)))return
 enddo
 staged=.5d0*left+.5d0*right
 status=4
 do i=1,3
  if(.not.ieee_is_finite(staged(i)))return
 enddo
 selected=staged;status=0
end function
end module
