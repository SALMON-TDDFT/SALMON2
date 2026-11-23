! Lifecycle/ownership regression using actual pinned fast-PR ACE/history kernels.
program test_state
 use iso_fortran_env,only:int64
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_is_nan,ieee_quiet_nan,ieee_positive_inf,ieee_negative_inf
 use exx_ace,only:s_exx_ace,exx_ace_build,exx_ace_apply
 use exx_factor_history,only:s_factor_history,history_accept,history_predict
 use exx_rt_state,only:s_exx_rt_state,rt_taylor_stage,rt_ptcn_begin,rt_state_clear,rt_state_bytes
 use legacy_lifecycle,only:old_stage=>exx_taylor_stage
 use legacy_lifecycle,only:old_ace=>ace,old_initial=>initial_ace,old_midpoint=>midpoint_ace, &
  old_source=>full_source,old_initial_source=>initial_source,old_active=>taylor_active,old_mid=>taylor_midpoint,propagator
 implicit none
 type(s_exx_rt_state)::a,b,reference_a,reference_b,accounted
 type(s_factor_history)::ha,hb
 complex(8)::u(6,2,1),w(6,2,1),target(6,1,1),out_a(6,1,1),out_b(6,1,1),pred(6,2),expected(6,2)
 complex(8),allocatable::held(:,:,:)
 real(8)::invalid(5),nan
 integer::status,step,k
 integer(int64)::bytes_a,bytes_b,allocated_a
 nan=ieee_value(0d0,ieee_quiet_nan)
 invalid=[nan,ieee_value(0d0,ieee_positive_inf),ieee_value(0d0,ieee_negative_inf),0d0,-1d0]
 u=0;u(1,1,1)=1;u(2,2,1)=1;target=cmplx(.2d0,.1d0,8)
 allocate(a%factor_history(1),b%factor_history(1))
 a%factor_history%interval=1;a%factor_history%n_teach=1
 b%factor_history%interval=1;b%factor_history%n_teach=1
 ha%interval=1;ha%n_teach=1;hb%interval=1;hb%n_teach=1
 call rt_ptcn_begin(a,.001d0,.true.,status)
 if(status/=0)stop 1
 call rt_ptcn_begin(b,.002d0,.true.,status)
 if(status/=0)stop 2
 do step=0,7
  w=-u*(1d0+.02d0*step)
  call exx_ace_build(a%ace,u,w,1d0,status)
  if(status/=0)stop 3
  w=-u*(2d0+.03d0*step)
  call exx_ace_build(b%ace,u,w,1d0,status)
  if(status/=0)stop 4
  reference_a%ace=a%ace;reference_b%ace=b%ace
  call history_accept(a%factor_history(1),a%ace%factors(:,:,1),1d0,step,sumgrid,status)
  if(status/=0)stop 5
  call history_accept(ha,reference_a%ace%factors(:,:,1),1d0,step,sumgrid,status)
  if(status/=0)stop 6
  call history_accept(b%factor_history(1),b%ace%factors(:,:,1),1d0,step,sumgrid,status)
  if(status/=0)stop 7
  call history_accept(hb,reference_b%ace%factors(:,:,1),1d0,step,sumgrid,status)
  if(status/=0)stop 8
  if(any(a%factor_history(1)%x/=ha%x).or.any(b%factor_history(1)%x/=hb%x))stop 9
  call exx_ace_apply(a%ace,target,out_a,status)
  if(status/=0)stop 10
  call exx_ace_apply(reference_a%ace,target,out_b,status)
  if(status/=0.or.any(out_a/=out_b))stop 11
  call exx_ace_apply(b%ace,target,out_a,status)
  if(status/=0)stop 32
  call exx_ace_apply(reference_b%ace,target,out_b,status)
  if(status/=0.or.any(out_a/=out_b))stop 33
  ! Preserve the old Taylor stage0/1/2 behavior, while B remains live.
  old_ace=a%ace;old_midpoint=a%ace;propagator='hse_taylor4'
  a%midpoint_ace=a%ace
  call old_stage(0)
  call rt_taylor_stage(a,0,.false.,status)
  if(status/=0.or.any(a%initial_ace%factors/=old_initial%factors))stop 12
  held=b%ace%factors
  call old_stage(1)
  call rt_taylor_stage(a,1,.false.,status)
  if(status/=0.or.(a%taylor_midpoint.neqv.old_mid))stop 13
  if(any(b%ace%factors/=held))stop 14
  call old_stage(2)
  call rt_taylor_stage(a,2,.false.,status)
  if(status/=0.or.allocated(a%initial_ace%factors).or.allocated(a%midpoint_ace%factors))stop 15
 end do
 if(.not.a%factor_history(1)%ready.or..not.b%factor_history(1)%ready)stop 16
 call history_predict(a%factor_history(1),1,pred,status)
 if(status/=0)stop 17
 call history_predict(ha,1,expected,status)
 if(status/=0.or.any(pred/=expected))stop 18
 ! Full-source Taylor path must copy/clear only that state's source.
 allocate(a%full_source(6,2,1));a%full_source=u
 old_source=u;propagator='hse_taylor4_full'
 call old_stage(0)
 call rt_taylor_stage(a,0,.true.,status)
 if(status/=0.or.any(a%initial_source/=old_initial_source))stop 19
 call rt_taylor_stage(a,99,.true.,status)
 if(status==0.or..not.a%taylor_active.or.any(a%initial_source/=u))stop 20
 call rt_taylor_stage(a,2,.true.,status)
 if(status/=0.or.allocated(a%initial_source))stop 21
 a%ptcn_attempted=.true.;a%ptcn_predicted=.true.
 do k=1,size(invalid)
  call rt_ptcn_begin(a,invalid(k),.true.,status)
  if(status==0.or..not.a%ptcn_attempted.or..not.a%ptcn_predicted.or.a%ptcn_history_dt/=.001d0)stop 22
 end do
 a%ptcn_history_dt=nan
 call rt_ptcn_begin(a,.001d0,.true.,status)
 if(status==0.or..not.a%ptcn_attempted.or..not.ieee_is_nan(a%ptcn_history_dt))stop 34
 a%ptcn_history_dt=.001d0
 call rt_ptcn_begin(a,.003d0,.true.,status)
 if(status==0.or..not.a%ptcn_attempted.or.a%ptcn_history_dt/=.001d0)stop 23
 call rt_ptcn_begin(b,.002d0,.true.,status)
 if(status/=0.or.b%ptcn_attempted.or.b%ptcn_predicted)stop 24
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a<=0)stop 25
 allocated_a=bytes_a
 call rt_state_bytes(b,bytes_b,status)
 if(status/=0.or.bytes_b<=0)stop 26
 call rt_state_clear(a)
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a/=0.or.a%ptcn_history_dt/=0d0)stop 27
 call history_predict(b%factor_history(1),1,pred,status)
 if(status/=0)stop 28
 call history_predict(hb,1,expected,status)
 if(status/=0.or.any(pred/=expected))stop 29
 call rt_state_bytes(b,bytes_a,status)
 if(status/=0.or.bytes_a/=bytes_b)stop 30
 ! Recreate A at a different dt; destroy B first and keep A usable.
 call rt_taylor_stage(a,0,.true.,status)
 if(status==0.or.a%taylor_active)stop 35
 w=-3d0*u
 call exx_ace_build(a%ace,u,w,1d0,status)
 if(status/=0)stop 36
 reference_a%ace=a%ace
 call rt_ptcn_begin(a,.004d0,.true.,status)
 if(status/=0)stop 37
 call rt_state_clear(b)
 call exx_ace_apply(a%ace,target,out_a,status)
 if(status/=0)stop 38
 call exx_ace_apply(reference_a%ace,target,out_b,status)
 if(status/=0.or.any(out_a/=out_b))stop 39
 call rt_state_clear(a)
 call rt_state_bytes(b,bytes_b,status)
 if(status/=0.or.bytes_b/=0)stop 31
 call accounting_check(accounted)
 print *, 'Owned allocated bytes A/B before clear:' ,allocated_a,bytes_a
 print *, 'Owned allocated bytes after clear:',bytes_b
 print *,'Fast-PR HSE RT state: ACE/history/lifecycle isolation and rejected-dt atomicity PASS'
contains
 subroutine accounting_check(state)
  implicit none
  type(s_exx_rt_state),intent(inout)::state
  integer(int64)::actual,expected
  integer::j,code
  call allocate_ace(state%ace);call allocate_ace(state%scf_damped_ace)
  call allocate_ace(state%initial_ace);call allocate_ace(state%midpoint_ace)
  allocate(state%cached_occupation(2,3),state%cached_action(2,3,1),state%factor_occupation(2,3,1))
  allocate(state%full_source(2,3,1),state%initial_source(2,3,1),state%cached_source(2,3,1))
  allocate(state%target_work(2,3,1),state%action_work(2,3,1))
  allocate(state%factor_history(2),state%ptcn_history(3))
  do j=1,2
   allocate(state%factor_history(j)%x(2,3,3))
  enddo
  do j=1,3
   allocate(state%ptcn_history(j)%x(2,3,3))
  enddo
  ! Allocated payload includes the allocated history-array descriptor/fixed
  ! elements, plus their separately allocated factor arrays, not top-level state.
  expected=4_int64*(20_int64*16+14_int64*int(storage_size(0)/8,int64)) &
   +12_int64*8+36_int64*16+5_int64*18*16+5_int64*int(storage_size(ha)/8,int64)
  call rt_state_bytes(state,actual,code)
  if(code/=0.or.actual/=expected)stop 40
  print *, 'All allocated components: expected/actual dynamic bytes:',expected,actual
  call rt_state_clear(state)
  call rt_state_bytes(state,actual,code)
  if(code/=0.or.actual/=0)stop 41
 end subroutine
 subroutine allocate_ace(payload)
  implicit none
  type(s_exx_ace),intent(inout)::payload
  allocate(payload%factors(2,3,1),payload%values(5),payload%metric_factor(3,3))
  allocate(payload%metric_rows(2),payload%metric_cols(3),payload%offset(4),payload%row(5))
  payload%packed=.true.;payload%metric_distributed=.true.;payload%metric_comm=77
 end subroutine
 subroutine sumgrid(matrix)
  implicit none
  complex(8),intent(inout)::matrix(:,:)
 end subroutine
end program
