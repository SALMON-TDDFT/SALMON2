! Serial strict-endpoint collector. Call AFTER corrected-density exact ACE refresh.
! A stage-2 notification alone is insufficient evidence of an exact endpoint ACE.
module exx_surrogate_trace
  use iso_fortran_env,only:real64,int32
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use exx_ace,only:s_exx_ace
  use exx_surrogate_dense,only:surrogate_dense_snapshot,surrogate_project_dense
  use exx_surrogate_model,only:surrogate_fit,surrogate_predict
  implicit none
  private
  type,public::s_surrogate_trace
    integer::epoch=-1,rank_max=0,capacity=0,count=0,rank=0,last_generation=0,last_step=-1
    real(real64)::dt=0,rank_rtol=0,dv=0
    complex(real64),allocatable::q(:,:)
    complex(real64),allocatable,private::b(:,:,:)
    integer,allocatable::steps(:)
    logical,allocatable,private::accepted(:)
  end type
  type,public::s_factor_diagnostic
    integer::unit=0,last_step=0,ng=0,no=0,expected=0
    real(real64)::dv=0,dt=0
    logical::opened=.false.,finished=.false.
  end type
  public::trace_factor_write
  public::trace_shadow_fit,trace_shadow_compare
  public::trace_initialize,trace_endpoint,trace_candidate,trace_teacher,trace_write_diagnostic
contains

 logical function scalar_finite_real(a) result(ok)
  real(8),intent(in)::a(:)
  integer::i
  ok=.false.
  do i=1,size(a)
   if(.not.ieee_is_finite(a(i)))return
  enddo
  ok=.true.
 end function

  ! Bounded raw factor stream for offline response-basis diagnostics only.
  subroutine trace_factor_write(writer,path,ace,step,dt,expected,final,status)
    type(s_factor_diagnostic),intent(inout)::writer
    character(*),intent(in)::path
    type(s_exx_ace),intent(in)::ace
    integer,intent(in)::step,expected
    real(real64),intent(in)::dt
    logical,intent(in)::final
    integer,intent(out)::status
    integer::ios,ng,no,i,j,k
    status=1
    if(writer%finished.or.step/=writer%last_step+1.or.expected<step)return
    if(ace%packed.or.ace%metric_distributed.or..not.allocated(ace%factors))return
    if(size(ace%factors,3)/=1)return
    ng=size(ace%factors,1);no=size(ace%factors,2)
    if(min(ng,no,expected)<1.or..not.scalar_finite_real([dt,ace%dv]).or.min(dt,ace%dv)<=0)return
    ! Scalar IEEE inquiries avoid rank-three array temporaries on Fujitsu.
    do k=1,size(ace%factors,3)
      do j=1,no
        do i=1,ng
          if(.not.ieee_is_finite(real(ace%factors(i,j,k),real64)))return
          if(.not.ieee_is_finite(aimag(ace%factors(i,j,k))))return
        enddo
      enddo
    enddo
    if(final.neqv.(step==expected))return
    if(.not.writer%opened)then
      open(newunit=writer%unit,file=path,status='new',access='stream',form='unformatted',iostat=ios)
      if(ios/=0)return
      writer%opened=.true.;writer%ng=ng;writer%no=no;writer%expected=expected;writer%dv=ace%dv;writer%dt=dt
      write(writer%unit,iostat=ios)'SALMON_FACTOR_DIAG_V1   ',int(z'01020304',int32), &
        int(ng,int32),int(no,int32),int(expected,int32),ace%dv,dt
      if(ios/=0)return
    endif
    if(ng/=writer%ng.or.no/=writer%no.or.expected/=writer%expected)return
    if(abs(dt-writer%dt)>epsilon(dt)*dt.or.abs(ace%dv-writer%dv)>epsilon(ace%dv)*ace%dv)return
    write(writer%unit,iostat=ios)int(step,int32),ace%factors(:,:,1)
    if(ios/=0)return
    writer%last_step=step
    if(final)then
      write(writer%unit,iostat=ios)'COMPLETE'
      if(ios/=0)return
      close(writer%unit,iostat=ios)
      if(ios/=0)return
      writer%finished=.true.;writer%opened=.false.
    endif
    status=0
  end subroutine

  ! Diagnostic-only native fit: uses the fixed prefix, never later test data.
  subroutine trace_shadow_fit(state,train_end,horizon,lambda,coeff,status)
    type(s_surrogate_trace),intent(in)::state
    integer,intent(in)::train_end,horizon
    real(real64),intent(in)::lambda
    real(real64),intent(out)::coeff(3)
    integer,intent(out)::status
    complex(real64),allocatable::x(:,:,:,:),y(:,:,:)
    integer::origin,e,j,n,ne
    status=1;coeff=0
    if(train_end>state%count.or.train_end<5.or.horizon<1)return
    if(.not.all(state%steps(:train_end)==[(j,j=1,train_end)]))return
    ne=train_end-horizon-3;if(ne<1)return
    n=state%rank;allocate(x(n,n,3,ne),y(n,n,ne));e=0
    do origin=4,train_end-horizon
      e=e+1
      do j=1,3
        x(:,:,j,e)=horizon*(state%b(:,:,origin-j+1)-state%b(:,:,origin-j))
      enddo
      y(:,:,e)=state%b(:,:,origin+horizon)-state%b(:,:,origin)
    enddo
    call surrogate_fit(x,y,lambda,coeff,status)
  end subroutine

  subroutine trace_shadow_compare(state,horizon,coeff,pred,truth,frozen,linear,status)
    type(s_surrogate_trace),intent(in)::state
    integer,intent(in)::horizon
    real(real64),intent(in)::coeff(3)
    complex(real64),intent(out)::pred(:,:),truth(:,:),frozen(:,:),linear(:,:)
    integer,intent(out)::status
    integer::origin,j
    real(real64)::times(4)
    status=1;pred=0;truth=0;frozen=0;linear=0
    origin=state%count-horizon
    if(origin<4.or.horizon<1.or.state%count<5)return
    if(any(shape(pred)/=[state%rank,state%rank]).or.any(shape(truth)/=shape(pred)).or. &
       any(shape(frozen)/=shape(pred)).or.any(shape(linear)/=shape(pred)))return
    if(any(state%steps(origin-3:state%count)/=[(j,j=origin-3,state%count)]))return
    times=real(state%steps(origin-3:origin),real64)*state%dt
    call surrogate_predict(state%b(:,:,origin-3:origin),times,state%steps(state%count)*state%dt,coeff,pred,status)
    if(status/=0)return
    truth=state%b(:,:,state%count);frozen=state%b(:,:,origin)
    linear=frozen+horizon*(frozen-state%b(:,:,origin-1))
  end subroutine

  subroutine trace_initialize(state,epoch,dt,rank_max,capacity,rank_rtol,status)
    type(s_surrogate_trace),intent(inout)::state
    integer,intent(in)::epoch,rank_max,capacity
    real(real64),intent(in)::dt,rank_rtol
    integer,intent(out)::status
    type(s_surrogate_trace)::fresh
    status=1
    if(epoch<0.or.rank_max<1.or.capacity<1)return
    if(.not.scalar_finite_real([dt,rank_rtol]).or.dt<=0.or.rank_rtol<=0.or.rank_rtol>=1)return
    fresh%epoch=epoch;fresh%dt=dt;fresh%rank_max=rank_max;fresh%capacity=capacity;fresh%rank_rtol=rank_rtol
    state=fresh;status=0
  end subroutine

  ! Only this explicit accepted-entry API can create a teacher.
  subroutine trace_endpoint(state,ace,step,exact_generation,strict_corrected,accepted,status)
    type(s_surrogate_trace),intent(inout)::state
    type(s_exx_ace),intent(in)::ace
    integer,intent(in)::step,exact_generation
    logical,intent(in)::strict_corrected,accepted
    integer,intent(out)::status
    status=1
    if(.not.accepted)return
    call trace_capture(state,ace,step,exact_generation,strict_corrected,accepted,status)
  end subroutine

  ! RT capture has no acceptance evidence yet. Never label it as a teacher.
  subroutine trace_candidate(state,ace,step,exact_generation,strict_corrected,status)
    type(s_surrogate_trace),intent(inout)::state
    type(s_exx_ace),intent(in)::ace
    integer,intent(in)::step,exact_generation
    logical,intent(in)::strict_corrected
    integer,intent(out)::status
    call trace_capture(state,ace,step,exact_generation,strict_corrected,.false.,status)
  end subroutine

  ! Raw diagnostics, never a certified teacher/model file. External provenance
  ! sealing and accepted-run validation are required before any training use.
  subroutine trace_write_diagnostic(state,path,status)
    type(s_surrogate_trace),intent(in)::state
    character(*),intent(in)::path
    integer,intent(out)::status
    integer::unit,ios,closed,i,j,k
    status=1
    if(state%count<1.or..not.allocated(state%q))return
    open(newunit=unit,file=path,status='new',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,'(a)',iostat=ios)'SALMON_ACE_DIAGNOSTIC_V1'
    if(ios==0)write(unit,*,iostat=ios)state%epoch,size(state%q,1),state%rank,state%count,state%dt,state%dv
    do j=1,state%rank
      do i=1,size(state%q,1)
        if(ios==0)write(unit,'(2(es26.17e3,1x))',iostat=ios)real(state%q(i,j)),aimag(state%q(i,j))
      enddo
    enddo
    do k=1,state%count
      if(ios==0)write(unit,*,iostat=ios)state%steps(k),state%accepted(k)
      do j=1,state%rank
        do i=1,state%rank
          if(ios==0)write(unit,'(2(es26.17e3,1x))',iostat=ios)real(state%b(i,j,k)),aimag(state%b(i,j,k))
        enddo
      enddo
    enddo
    if(ios==0)write(unit,'(a)',iostat=ios)'END_ACE_DIAGNOSTIC'
    close(unit,iostat=closed)
    if(ios==0.and.closed==0)status=0
  end subroutine

  subroutine trace_teacher(state,index,b,status)
    type(s_surrogate_trace),intent(in)::state
    integer,intent(in)::index
    complex(real64),allocatable,intent(out)::b(:,:)
    integer,intent(out)::status
    status=1
    if(index<1.or.index>state%count)return
    if(.not.state%accepted(index))return
    b=state%b(:,:,index);status=0
  end subroutine

  subroutine trace_capture(state,ace,step,exact_generation,strict_corrected,accepted,status)
    type(s_surrogate_trace),intent(inout)::state
    type(s_exx_ace),intent(in)::ace
    integer,intent(in)::step,exact_generation
    logical,intent(in)::strict_corrected,accepted
    integer,intent(out)::status
    complex(real64),allocatable::initial_q(:,:),snapshot(:,:)
    integer::rank
    status=1
    if(state%epoch<0.or..not.strict_corrected.or.step<=state%last_step)return
    if(exact_generation<=state%last_generation)return
    if(state%count>=state%capacity)return
    if(.not.allocated(state%q))then
      call surrogate_dense_snapshot(ace,state%rank_max,state%rank_rtol,initial_q,snapshot,rank,status)
      if(status/=0)return
      state%q=initial_q;state%rank=rank;state%dv=ace%dv
      allocate(state%b(rank,rank,state%capacity),state%steps(state%capacity),state%accepted(state%capacity))
      state%b=0;state%steps=-1;state%accepted=.false.
    else
      ! Same grid-volume convention and fixed coordinates throughout the epoch.
      if(.not.ieee_is_finite(ace%dv))return
      if(abs(ace%dv-state%dv)>epsilon(1d0)*max(abs(state%dv),1d0))return
      allocate(snapshot(state%rank,state%rank))
      call surrogate_project_dense(ace,state%q,snapshot,status)
      if(status/=0)return
    endif
    state%count=state%count+1;state%b(:,:,state%count)=snapshot;state%steps(state%count)=step
    state%accepted(state%count)=accepted
    state%last_step=step
    state%last_generation=exact_generation
    status=0
  end subroutine
end module
