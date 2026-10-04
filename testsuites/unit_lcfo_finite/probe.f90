program test
use finite_probe
use,intrinsic :: ieee_arithmetic
implicit none
complex(8),allocatable::a(:,:,:,:)
logical::ok
complex(8)::b(4,4),c(4,4,2)
allocate(a(400,400,1,26))
a=(1d0,2d0)
call probe(a,ok)
if(.not.ok) stop 1
call probe(a(1:400:2,:,:,:),ok)
if(.not.ok) stop 2
a(1,1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
call probe(a,ok)
if(ok) stop 3
a(1,1,1,1)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
call probe(a,ok)
if(ok) stop 4
call probe(a(1:0,:,:,:),ok)
if(.not.ok) stop 5
b=(1d0,2d0)
c=(1d0,2d0)
call probe2(b,ok)
if(.not.ok) stop 6
call probe3(c,ok)
if(.not.ok) stop 7
b(4,4)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
c(4,4,2)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
call probe2(b,ok)
if(ok) stop 8
call probe3(c,ok)
if(ok) stop 9
call probe2(b(1:0,:),ok)
if(.not.ok) stop 10
call probe3(c(1:0,:,:),ok)
if(.not.ok) stop 11
print *, 'rank2/3/4 finite/noncontiguous/NaN/Inf/empty PASS'
end program
