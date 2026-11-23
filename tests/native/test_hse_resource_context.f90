! Synthetic ownership/replay check, not converged-material physics.
program probe
 use hse_resource_context
 use exx_k_exchange,only:exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply,exx_k_kernel_destroy
 use exx_ace,only:s_exx_ace,exx_ace_build,exx_ace_apply,exx_ace_clear
 use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(hse_context),target :: a,b
 complex(8) :: ua(8,2,1),ub(12,2,2),wa(8,2,1),wb(12,2,2),old(8,2,1),bad(8,2,1)
 complex(8) :: fixed(8,2,1),before(8,2,1),other(12,2,2)
 real(8) :: invalid_k(3,2)
 real(8) :: k(3,1),kb(3,2),ha(3),hb(3)
 integer :: status,i
 k=0;ha=.8d0;hb=[.7d0,.9d0,1.1d0]
 kb=0;kb(1,2)=acos(-1d0)/(3*hb(1))
 call orbitals(ua,product(ha));call orbitals(ub,product(hb))
 call hse_open(a,[2,2,2],[1,1,1],ha,k,.11d0,.25d0,status)
 if(status/=0)stop 1
 call hse_open(b,[3,2,2],[2,1,1],hb,kb,.2d0,.4d0,status)
 if(status/=0)stop 2
 call hse_refresh(a,ua,status);if(status/=0)stop 3
 call hse_refresh(b,ub,status);if(status/=0)stop 4
 do i=1,3
  call check(a,[2,2,2],[1,1,1],k,ha,.11d0,.25d0,ua,wa)
  call check(b,[3,2,2],[2,1,1],kb,hb,.2d0,.4d0,ub,wb)
  ! Refresh one material while the other retains accepted factors.
  fixed=ua
  call hse_apply(a,fixed,before,status);if(status/=0)stop 27
  call hse_apply(b,ub,other,status);if(status/=0)stop 28
  ua(:,1,1)=0
  ua(1,1,1)=cos(.2d0*i)/sqrt(product(ha))
  ua(3,1,1)=sin(.2d0*i)/sqrt(product(ha))
  call hse_refresh(a,ua,status);if(status/=0)stop 5
  call hse_apply(a,fixed,wa,status);if(status/=0)stop 29
  if(maxval(abs(wa-before))<1d-8)stop 30
  call hse_apply(b,ub,wb,status);if(status/=0.or.any(wb/=other))stop 31
 enddo
 call hse_apply(a,ua,old,status);if(status/=0)stop 6
 bad=ua;bad(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call hse_refresh(a,bad,status);if(status==0)stop 7
 call hse_refresh(a,ub,status);if(status==0)stop 8
 bad=0
 call hse_refresh(a,bad,status);if(status==0)stop 9
 call hse_apply(a,ua,wa,status);if(status/=0.or.any(wa/=old))stop 10
 wa=7
 call hse_apply(a,bad(:,:1,:),wa,status)
 if(status==0.or.any(wa/=7))stop 11
 call hse_open(a,[2,2,2],[1,1,1],ha,k,.11d0,.25d0,status)
 if(status==0)stop 12
 call hse_apply(a,ua,wa,status);if(status/=0.or.any(wa/=old))stop 13
 call hse_close(a)
 call check(b,[3,2,2],[2,1,1],kb,hb,.2d0,.4d0,ub,wb)
 call hse_close(a)
 call hse_open(a,[2,2,2],[1,1,1],ha,k,.11d0,.25d0,status)
 if(status/=0)stop 14
 call hse_refresh(a,ua,status);if(status/=0)stop 15
 call hse_close(b)
 call check(a,[2,2,2],[1,1,1],k,ha,.11d0,.25d0,ua,wa)
 call hse_close(a);call hse_close(b)
 wa=9;call hse_apply(a,ua,wa,status)
 if(status==0.or.any(wa/=9))stop 16
 call hse_open(a,[0,2,2],[1,1,1],ha,k,.11d0,.25d0,status)
 if(status==0)stop 17
 call hse_close(a)
 invalid_k=0
 call hse_open(a,[3,2,2],[2,1,1],hb,invalid_k,.2d0,.4d0,status)
 if(status==0)stop 32
 call hse_close(a);call hse_close(a)
 call hse_open(a,[2,2,2],[1,1,1],ha,k,.11d0,.25d0,status)
 if(status/=0)stop 33
 call hse_refresh(a,ua,status);if(status/=0)stop 34
 call check(a,[2,2,2],[1,1,1],k,ha,.11d0,.25d0,ua,wa)
 call hse_close(a)
 print *,'HSE resource context: distinct grids/settings, refresh, rejection and both close orders PASS'
contains
 subroutine orbitals(u,dv)
  complex(8),intent(out)::u(:,:,:)
  real(8),intent(in)::dv
  integer::j
  u=0
  do j=1,size(u,2)
   u(j,j,:)=1/sqrt(dv)
  enddo
 end subroutine
 subroutine check(ctx,n,mesh,points,h,omega,fraction,u,found)
  type(hse_context),target,intent(inout)::ctx
  integer,intent(in)::n(3),mesh(3)
  real(8),intent(in)::h(3),points(:,:),omega,fraction
  complex(8),intent(in)::u(:,:,:)
  complex(8),intent(out)::found(:,:,:)
  type(exx_k_kernel),target::reference
  type(s_exx_ace)::ace
  complex(8)::w(size(u,1),size(u,2),size(u,3)),expected(size(u,1),size(u,2),size(u,3))
  complex(8)::target(size(u,1),size(u,2),size(u,3))
  integer::code,ii,jj,kk
  call exx_k_kernel_init(reference,n,mesh,h,points,omega,2,code,fft_layout='strided')
  if(code/=0)stop 20
  call exx_k_kernel_apply(reference,u,u,w,0,1,code);if(code/=0)stop 21
  call exx_ace_build(ace,u,w,product(h),code);if(code/=0)stop 22
  call exx_ace_apply(ace,u,expected,code);if(code/=0)stop 23
  call hse_apply(ctx,u,found,code);if(code/=0)stop 24
  if(maxval(abs(found-fraction*expected))>1d-12)stop 25
  if(maxval(abs(found-fraction*w))>1d-12)stop 26
  do kk=1,size(u,3)
   do jj=1,size(u,2)
    do ii=1,size(u,1)
     target(ii,jj,kk)=cmplx(.01d0*ii,.02d0*jj+.03d0*kk,8)
    enddo
   enddo
  enddo
  call exx_ace_apply(ace,target,expected,code);if(code/=0)stop 35
  call hse_apply(ctx,target,found,code);if(code/=0)stop 36
  if(maxval(abs(found-fraction*expected))>1d-12)stop 37
  call exx_ace_clear(ace);call exx_k_kernel_destroy(reference)
 end subroutine
end program
