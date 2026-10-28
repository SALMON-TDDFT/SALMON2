program fresh_exchange_probe
 use ptcn_core
 implicit none
 complex(8)::u(2,1,1),x(2,1,1),hu(2,1,1),r(2,1,1),rhs(2,1,1)
 real(8)::res,cache,dt,physical
 integer::builds,apps,err,fresh_count,inner_count
 u=0;u(1,1,1)=1;dt=.01d0;cache=0;fresh_count=0;inner_count=0
 call exx_ptcn_solve(u,dt,1d0,action,prec,norm,x,res,builds,apps,err,use_anderson=.true.)
 if(err/=0.or.res>=1d-10)error stop 'shared solver convergence'
 if(fresh_count<3.or.builds/=fresh_count.or.inner_count/=apps)error stop 'fresh exchange gate skipped'
 hu=0;hu(2,1,1)=.1d0
 call pt_gradient(u,hu,1d0,r,1);rhs=u-(0d0,.5d0)*dt*r
 hu(1,1,1)=.2d0*x(2,1,1);hu(2,1,1)=.2d0*x(1,1,1)
 call pt_gradient(x,hu,1d0,r,1)
 physical=norm(x+(0d0,.5d0)*dt*r-rhs)
 if(physical>=1d-10)error stop 'physical exchange residual differs'
 print *,'FRESH_EXCHANGE_ACCEPTANCE_PASS',physical,builds,apps
contains
 subroutine action(a,b,refresh)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 logical,intent(in)::refresh
 if(refresh)then
  fresh_count=fresh_count+1;cache=.1d0
  if(fresh_count>1)cache=.2d0
 else
  inner_count=inner_count+1
 endif
 b(1,1,1)=cache*a(2,1,1);b(2,1,1)=cache*a(1,1,1)
 end subroutine
 subroutine prec(a,b)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 b=a
 end subroutine
 real(8) function norm(a)
 complex(8),intent(in)::a(:,:,:)
 norm=sqrt(sum(abs(a)**2))
 end function
end program
