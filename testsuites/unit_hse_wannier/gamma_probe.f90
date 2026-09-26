program gamma_probe
 use hse_wannier_gauge
 implicit none
 integer,parameter :: n=8
 complex(8) :: u(n,n,1),q(n,n),raw(n,n,6,1),diag(n,n),grad(n,n,1),overlap(n,n)
 real(8) :: b(3,6),weights(6),spread,variable,gradient,pi,theta,initial
 integer :: i,j,a,iterations,status,neighbors(6,1)
 pi=acos(-1d0);u=0d0;b=0d0;neighbors=1
 do i=1,n
   u(i,i,1)=1d0
   do j=1,n
     q(i,j)=exp(cmplx(0d0,2*pi*(i-1)*(j-1)/n,8))/sqrt(real(n,8))
   enddo
 enddo
 do a=1,3
   b(a,a)=2*pi/merge(164.16d0,10.26d0,a==1);b(a,a+3)=-b(a,a)
   weights(a)=1d0/(2*b(a,a)**2);weights(a+3)=weights(a)
   diag=0d0
   do i=1,n
     theta=.37d0*i+.09d0*a*i*i
     diag(i,i)=exp(cmplx(0d0,theta,8))
   enddo
   raw(:,:,a,1)=matmul(conjg(transpose(q)),matmul(diag,q))
   raw(:,:,a+3,1)=conjg(transpose(raw(:,:,a,1)))
 enddo
 call gauge_minimize_gamma(u,raw,b,weights,200,1d-7,spread,gradient,iterations,status)
 write(*,*) 'Gamma optimizer status/sweeps/spread/gradient',status,iterations,spread,gradient
 if(status/=0.or.gradient>1d-7.or.abs(spread)>1d-7)error stop 'Gamma localization convergence'
 call gauge_functional(u,raw,neighbors,b,weights,spread,variable,grad,status)
 if(status/=0.or.sqrt(sum(abs(grad)**2))>1d-7)error stop 'Independent MV gradient'
 overlap=matmul(conjg(transpose(u(:,:,1))),u(:,:,1))
 do i=1,n;overlap(i,i)=overlap(i,i)-1d0;enddo
 if(maxval(abs(overlap))>1d-12)error stop 'Gamma rotation unitarity'
 write(*,*) 'Gamma MV passed: sweeps, spread, gradient=',iterations,spread,gradient
 raw(1,2,1,1)=raw(1,2,1,1)+cmplx(.003d0,.004d0,8)
 raw(2,1,1,1)=raw(2,1,1,1)+cmplx(-.002d0,.006d0,8)
 raw(:,:,4,1)=conjg(transpose(raw(:,:,1,1)))
 call gauge_functional(u,raw,neighbors,b,weights,initial,variable,grad,status)
 call gauge_minimize_gamma(u,raw,b,weights,200,1d-7,spread,gradient,iterations,status)
 if(status/=0.or.gradient>1d-7.or.spread>initial+1d-10)error stop 'Noncommuting Gamma links'
 write(*,*) 'Noncommuting links passed:',iterations,gradient

end program
