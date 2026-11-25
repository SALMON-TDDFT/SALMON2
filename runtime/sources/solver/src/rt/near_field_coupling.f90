! Field-representation-independent molecular coupling in Gaussian atomic units.
! A field provider supplies the electric integral and the s-weighted magnetic integral.
module salmonpy_near_field_coupling
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::exponential_path_moments,near_coupling_map,sp_near_couple_samples
contains
 subroutine exponential_path_moments(argument,degree,moments,status)
  real(c_double),intent(in)::argument
  integer,intent(in)::degree
  real(c_double),intent(out)::moments(0:)
  integer,intent(out)::status
  real(c_double)::term,total,decay,b
  integer::n,k
  status=1
  if(degree<0.or.degree>9.or.size(moments)<degree+1)return
  if(.not.ieee_is_finite(argument))return
  status=4
  if(argument==0d0)then
   do n=0,degree
    moments(n)=1d0/(n+1)
   enddo
  elseif(argument>20d0)then
   ! Upward recurrence is stable here; no small-argument cancellation.
   decay=exp(-argument);moments(0)=(1d0-decay)/argument
   do n=1,degree
    moments(n)=(n*moments(n-1)-decay)/argument
   enddo
  else
   b=-argument
   if(b>log(huge(1d0)))return
   do n=0,degree
    if(argument>0d0)then
     ! Positive series for exp(-a) integral, never subtract nearly equal values.
     term=exp(-argument)/(n+1)
    else
     term=1d0/(n+1)
    endif
    total=term
    do k=1,10000
     if(argument>0d0)then
      term=term*argument/(n+k+1)
     else
      term=term*(b/k)*(real(n+k,c_double)/(n+k+1))
     endif
     total=total+term
     if(.not.ieee_is_finite(total))return
     if(k>abs(argument).and.term<=epsilon(1d0)*total)exit
    enddo
    if(k>10000)return
    moments(n)=total
   enddo
  endif
  do n=0,degree
   if(.not.ieee_is_finite(moments(n)))return
  enddo
  status=0
 end subroutine
 subroutine near_coupling_map(r,electric_integral,magnetic_integral,potential,Ac)
  real(c_double),intent(in)::r(3),electric_integral(3),magnetic_integral(3)
  real(c_double),intent(out)::potential,Ac(3)
  potential=dot_product(r,electric_integral)
  Ac=-[r(2)*magnetic_integral(3)-r(3)*magnetic_integral(2), &
       r(3)*magnetic_integral(1)-r(1)*magnetic_integral(3), &
       r(1)*magnetic_integral(2)-r(2)*magnetic_integral(1)]/137.035999084d0
 end subroutine
 integer(c_int) function sp_near_couple_samples(n,r,parameters,weights,E,H,potential,Ac) bind(C) result(status)
  ! Quadrature samples are at center+s*r, 0<=s<=1; H=B in molecular vacuum.
  ! A future grid interpolator can use this entry without Ylm coefficients.
  integer(c_int),value::n
  real(c_double),intent(in)::r(3),parameters(n),weights(n),E(3,n),H(3,n)
  real(c_double),intent(inout)::potential,Ac(3)
  real(c_double)::ei(3),hi(3),vp,vector(3)
  integer::i,k
  status=1
  if(n<1)return
  do k=1,3
   if(.not.ieee_is_finite(r(k)))return
  enddo
  ei=0d0;hi=0d0
  do i=1,n
   if(.not.ieee_is_finite(parameters(i)).or..not.ieee_is_finite(weights(i)))return
   if(parameters(i)<0d0.or.parameters(i)>1d0.or.weights(i)<0d0)return
   do k=1,3
    if(.not.ieee_is_finite(E(k,i)).or..not.ieee_is_finite(H(k,i)))return
   enddo
   ei=ei+weights(i)*E(:,i);hi=hi+weights(i)*parameters(i)*H(:,i)
  enddo
  call near_coupling_map(r,ei,hi,vp,vector)
  status=4
  if(.not.ieee_is_finite(vp))return
  do k=1,3
   if(.not.ieee_is_finite(vector(k)))return
  enddo
  potential=vp;Ac=vector;status=0
 end function
end module
