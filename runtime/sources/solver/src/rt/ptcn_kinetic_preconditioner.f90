! Fourier symbol of the SALMON nonorthogonal finite-difference kinetic operator.
module ptcn_kinetic_preconditioner
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::ptcn_kinetic_symbol
contains
 subroutine ptcn_kinetic_symbol(q,lap0,lap,nab,f,B,k,value,err)
 real(8),intent(in)::q(3),lap0,lap(4,3),nab(4,3),f(6),B(3,3),k(3)
 real(8),intent(out)::value
 integer,intent(out)::err
 real(8)::grad(3),second(3),bk(3)
 integer::i,j
 err=1;value=0d0
 if(.not.ieee_is_finite(lap0))return
 do i=1,3
  if(.not.ieee_is_finite(q(i)).or..not.ieee_is_finite(k(i)))return
  do j=1,3
   if(.not.ieee_is_finite(B(i,j)))return
  enddo
  do j=1,4
   if(.not.ieee_is_finite(lap(j,i)).or..not.ieee_is_finite(nab(j,i)))return
  enddo
 enddo
 do i=1,6
  if(.not.ieee_is_finite(f(i)))return
 enddo
 grad=0d0;second=0d0
 do i=1,3;do j=1,4
  grad(i)=grad(i)+2*nab(j,i)*sin(j*q(i))
  second(i)=second(i)+2*lap(j,i)*cos(j*q(i))
 enddo;enddo
 bk=matmul(B,k)
 value=lap0+.5d0*sum(k*k)-.5d0*sum(f(1:3)*second)+sum(bk*grad) &
  +.5d0*(f(4)*grad(2)*grad(3)+f(5)*grad(3)*grad(1)+f(6)*grad(1)*grad(2))
 if(.not.ieee_is_finite(value))return
 err=0
 end subroutine
end module
