program probe
 use exx_k_exchange
 use exx_spinor,only:exx_spinor_apply_distributed,exx_spinor_weight_sources
 use exx_ace
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 type(exx_k_kernel)::op
 type(s_exx_ace)::ace,history,mixed,physical
 complex(8)::u(16,3,1),q(16,3,1),w(16,3,1),ref(16,3,1),wa(16,3,1)
 complex(8)::dm
 real(8)::f(3,1),k(3,1),residual
 integer::i,j,a,s,t,err,delta(3)
 character(30)::mode
 call get_command_argument(1,mode)
 f(:,1)=[1d0,.37d0,0d0]
 select case(trim(mode))
 case('nan');f(2,1)=ieee_value(0d0,ieee_quiet_nan)
 case('inf');f(2,1)=ieee_value(0d0,ieee_positive_inf)
 case('negative');f(2,1)=-.1d0
 case('overfull');f(2,1)=1.1d0
 end select
 do a=1,3;do i=1,16
  u(i,a,1)=cmplx(sin(real(i+3*a,8)),cos(real(2*i-a,8)),8)/5d0
 enddo;enddo
 call exx_spinor_weight_sources(u,f,q,err)
 if(len_trim(mode)>0)then
  if(err==0)error stop 'invalid occupation accepted'
  print *,'PASS invalid occupation rejected';stop
 endif
 if(err/=0)error stop 'valid occupation rejected'
 if(maxval(abs(q(:,3,1)))>0d0)error stop 'empty source nonzero'
 k=0d0
 call exx_k_kernel_init(op,2,1,.7d0,k,.11d0,3,err,fft_layout='strided')
 if(err/=0)stop 1
 call exx_spinor_apply_distributed(op,q,u,w,[1],[1],0,0,transpose_tiles,err)
 if(err/=0)stop 2
 ref=0d0
 do i=1,8;do j=1,8
  delta=modulo(op%point(:,i)-op%point(:,j),op%n)
  do s=1,2;do t=1,2
   dm=sum(f(:,1)*u(i+(s-1)*8,:,1)*conjg(u(j+(t-1)*8,:,1)))
   ref(i+(s-1)*8,:,1)=ref(i+(s-1)*8,:,1)-dm*op%kernel(delta(1),delta(2),delta(3))*u(j+(t-1)*8,:,1)
  enddo;enddo
 enddo;enddo
 if(maxval(abs(ref-w))>1d-11)error stop 'weighted Fock mismatch'
 call exx_ace_build(ace,u,w,.7d0**3,err)
 if(err/=0)error stop 'thermal ACE metric failed'
 call exx_ace_apply(ace,u,wa,err)
 if(err/=0.or.maxval(abs(wa-w))>1d-10)error stop 'thermal ACE interpolation failed'
 f=1d0
 call exx_spinor_weight_sources(u,f,q,err)
 if(err/=0.or.maxval(abs(q-u))>0d0)error stop 'integer limit mismatch'
 call exx_ace_build(history,u,-2d0*u,.7d0**3,err)
 if(err/=0)error stop 'history fixture'
 call exx_ace_build(physical,u,-4d0*u,.7d0**3,err)
 if(err/=0)error stop 'physical fixture'
 call exx_ace_scf_mix(history,physical,u,.25d0,mixed,residual,err)
 if(err/=0.or.abs(residual-.5d0)>1d-10)error stop 'SCF operator residual'
 call exx_ace_apply(mixed,u,wa,err)
 if(err/=0.or.maxval(abs(wa+2.5d0*u))>1d-10)error stop 'SCF operator damping'
 call exx_ace_apply(physical,u,wa,err)
 if(err/=0.or.maxval(abs(wa+4d0*u))>1d-10)error stop 'physical operator mutated'
 call exx_ace_scf_mix(history,physical,u,ieee_value(0d0,ieee_quiet_nan),mixed,residual,err)
 if(err==0)error stop 'NaN mixing accepted'
 call exx_ace_scf_mix(history,physical,u,0d0,mixed,residual,err)
 if(err==0)error stop 'zero mixing accepted'
 print *,'PASS SCF ACE damping, residual and physical operator preservation'
 print *,'PASS weighted spinor Fock, empty source/target, ACE and integer limit'
contains
 subroutine transpose_tiles(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  recv=send
 end subroutine
end program
