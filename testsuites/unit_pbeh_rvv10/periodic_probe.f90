program probe
 use rvv10
 implicit none
 integer,parameter :: n(3)=[5,4,3]
 real(8) :: rho(5,4,3),e(5,4,3),v(5,4,3),ep(5,4,3),em(5,4,3),work(5,4,3)
 real(8) :: h(3)=[.7d0,.8d0,.9d0],coef(4,3),bmat(3,3),delta,original,err
 integer :: i,j,k,status
 bmat=0;do i=1,3
 bmat(i,i)=1
 coef(:,i)=[.8d0,-.2d0,4d0/105,-1d0/280]/h(i)
 enddo
 do k=1,3;do j=1,4;do i=1,5
 rho(i,j,k)=.04d0+.01d0*sin(real(i+3*j+7*k,8))
 enddo;enddo;enddo
 call rvv10_periodic(n,h,coef,bmat,rho,5.3d0,.0093d0,32,e,v,status)
 if(status/=0)error stop 'periodic evaluation'
 delta=1d-7;err=0
 do k=1,3;do j=1,4;do i=1,5
 original=rho(i,j,k);rho(i,j,k)=original+delta
 call rvv10_periodic(n,h,coef,bmat,rho,5.3d0,.0093d0,32,ep,work,status)
 rho(i,j,k)=original-delta
 call rvv10_periodic(n,h,coef,bmat,rho,5.3d0,.0093d0,32,em,work,status)
 rho(i,j,k)=original
 err=max(err,abs(sum(ep-em)/(2*delta)-v(i,j,k)))
 enddo;enddo;enddo
 if(err>2d-8)error stop 'periodic functional derivative'
 print *,err
end program
