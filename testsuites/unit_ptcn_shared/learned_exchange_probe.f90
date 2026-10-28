program learned_exchange_probe
 use ptcn_core
 implicit none
 complex(8)::u(2,1,1),x(2,1,1),reference(2,1,1)
 real(8)::res,cache
 integer::mode,builds,apps,err,fresh_count,inner_count,ordinary_builds
 logical::attempted
 u=0;u(1,1,1)=1
 do mode=0,2
  cache=0;fresh_count=0;inner_count=0;attempted=.false.
  call exx_ptcn_solve(u,.01d0,1d0,action,prec,norm,x,res,builds,apps,err,use_anderson=.true.)
  if(err/=0.or.res>=1d-10)error stop 'physical endpoint did not converge'
  if(builds/=fresh_count)error stop 'fresh action accounting'
  if(mode==0)then
   reference=x;ordinary_builds=builds
  else
   if(maxval(abs(x-reference))>1d-10)error stop 'learned differs from exact solution'
   if(mode==1.and.builds>=ordinary_builds)error stop 'accurate predictor did not save rebuild'
   if(mode==2.and.builds<3)error stop 'bad prediction skipped exact correction'
  endif
  print *,'LEARNED_EXACT_GATE_PASS mode/builds/apps/error',mode,builds,apps,res
 enddo
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
  if(.not.attempted)then
   attempted=.true.
   if(mode==1)cache=.2d0
   if(mode==2)cache=.4d0
  endif
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
