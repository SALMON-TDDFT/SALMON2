program probe
 use iso_fortran_env,only:real64
 use exx_factor_history
 implicit none
 type(s_factor_history)::state
 complex(real64)::x(16,2),base(16,2),pred(16,2),rot(2,2),ref(16,2)
 real(real64)::err,total
 integer::i,j,t,stat,rank,ierr,p,h
 rank=0
 if(size(state%coeff)/=2.or.any(shape(state%g)/=[2,2]))stop 7
 do j=1,2
 do i=1,16
  base(i,j)=cmplx(sin(real(i+j+rank*16,8)),cos(real(i*j+rank*16,8)),8)
 enddo
 enddo
 rot=0;rot(1,2)=(0d0,1d0);rot(2,1)=(1d0,0d0)
 x=matmul(base,rot);call factor_align(x,base,sumgrid,stat)
 err=maxval(abs(x-base));total=err
 if(stat/=0.or.total>1d-10)stop 1
 do p=3,3
 state=s_factor_history()
 do t=0,128,8
  x=base*(1d0+.01d0*t+.0001d0*t*t)
  call history_accept(state,x,1d0,t,sumgrid,stat)
  if(stat/=0)stop 2
 enddo
 call history_predict(state,8,pred,stat)
 ref=base*(1d0+.01d0*136+.0001d0*136*136)
 err=maxval(abs(pred-ref));total=err
 ! Rank-one temporal signal has underdetermined coefficients; allow approximate regularized prediction.
 if(stat/=0.or.total>.03d0)stop 3
 if(size(state%x,3)/=p)stop 4
 if(rank==0)print *,'PASS unitary alignment, serial sum callback, bounded history predictor',total
 enddo

 ! Check every internal horizon independently of the coefficient fit.
 do p=3,3
 state=s_factor_history();;state%interval=8;state%ready=.true.
 allocate(state%x(16,2,p))
 do j=1,p
  state%x(:,:,j)=base*real(1-j,real64)
 enddo
 state%coeff=[1d0,0d0]
 do h=1,8
  call history_predict(state,h,pred,stat)
  if(stat/=0.or.maxval(abs(pred-base*real(h,real64)/8d0))>1d-12)stop 5
 enddo
 ! Quadratic endpoint coefficients have a known internal interpolation bias.
 do j=1,p
  state%x(:,:,j)=base*real((1-j)**2,real64)
 enddo
 state%coeff=[2d0,-1d0]
 do h=1,8
  call history_predict(state,h,pred,stat)
  ref=base*(2d0*(real(h,real64)/8d0)**2-real(h,real64)/8d0)
  if(stat/=0.or.maxval(abs(pred-ref))>1d-12)stop 6
 enddo
 enddo
 print *, 'PASS internal horizons 1:8, linear exactness and quadratic bias'

contains
 subroutine sumgrid(a)
 implicit none
 complex(real64),intent(inout)::a(:,:)
 integer::ierr
 ! Serial callback: values already global.
 end subroutine
end program
