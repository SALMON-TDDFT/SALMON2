module rt_dg_wpw_exp_production
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_exp_state
    logical::valid=.false.
    integer::step=0,max_corrector=0
    real(8)::dt=0d0,corrector_tolerance=0d0,norm_tolerance=0d0
    complex(8),allocatable::coeff(:,:)
  end type
  abstract interface
    subroutine dg_wpw_midpoint_operator(iteration,midpoint_coefficients,h_reduced,info)
      integer,intent(in)::iteration
      complex(8),intent(in)::midpoint_coefficients(:,:)
      complex(8),intent(out)::h_reduced(:,:)
      integer,intent(out)::info
    end subroutine
  end interface
  public::initialize_dg_wpw_exp_state,advance_dg_wpw_midpoint_exp,advance_dg_wpw_fieldoff_exp
  public::advance_dg_wpw_length_gauge_exp
contains
  subroutine initialize_dg_wpw_exp_state(state,coeff,dt,max_corrector,corrector_tolerance,norm_tolerance,info)
    type(s_dg_wpw_exp_state),intent(out)::state
    complex(8),intent(in)::coeff(:,:)
    real(8),intent(in)::dt,corrector_tolerance,norm_tolerance
    integer,intent(in)::max_corrector
    integer,intent(out)::info
    info=1
    if(size(coeff,1)<=0.or.size(coeff,2)<=0.or.dt<=0d0.or.max_corrector<1.or.&
       corrector_tolerance<=0d0.or.norm_tolerance<=0d0.or..not.finite2(coeff))return
    state%coeff=coeff;state%dt=dt;state%max_corrector=max_corrector
    state%corrector_tolerance=corrector_tolerance;state%norm_tolerance=norm_tolerance
    state%valid=.true.;info=0
  end subroutine

  subroutine advance_dg_wpw_length_gauge_exp(state,eigenvalues,position,field,iterations,norm_drift,info)
    type(s_dg_wpw_exp_state),intent(inout)::state
    real(8),intent(in)::eigenvalues(:),field
    complex(8),intent(in)::position(:,:)
    integer,intent(out)::iterations,info
    real(8),intent(out)::norm_drift
    if(size(eigenvalues)/=size(state%coeff,1).or.any(shape(position)/=[size(eigenvalues),size(eigenvalues)]).or.&
      .not.ieee_is_finite(field).or..not.finite2(position))then
      iterations=0;norm_drift=huge(1d0);info=1;return
    endif
    call advance_dg_wpw_midpoint_exp(state,field_operator,iterations,norm_drift,info)
  contains
    subroutine field_operator(iteration,midpoint_coefficients,h_reduced,op_info)
      integer,intent(in)::iteration
      complex(8),intent(in)::midpoint_coefficients(:,:)
      complex(8),intent(out)::h_reduced(:,:)
      integer,intent(out)::op_info
      integer::i
      h_reduced=field*position
      do i=1,size(eigenvalues);h_reduced(i,i)=h_reduced(i,i)+eigenvalues(i);enddo
      op_info=merge(0,1,iteration>0.and.size(midpoint_coefficients,1)==size(eigenvalues))
    end subroutine
  end subroutine

  subroutine advance_dg_wpw_midpoint_exp(state,midpoint_operator,iterations,norm_drift,info)
    type(s_dg_wpw_exp_state),intent(inout)::state
    procedure(dg_wpw_midpoint_operator)::midpoint_operator
    integer,intent(out)::iterations,info
    real(8),intent(out)::norm_drift
    complex(8),allocatable::saved(:,:),previous(:,:),midpoint(:,:),trial(:,:),h(:,:)
    real(8)::saved_norm,trial_norm,delta
    integer::iter,op_info,exp_info,n
    logical::converged
    info=1;iterations=0;norm_drift=huge(1d0)
    if(.not.state%valid)return
    n=size(state%coeff,1);saved=state%coeff;previous=saved
    allocate(midpoint(n,size(saved,2)),trial(n,size(saved,2)),h(n,n))
    saved_norm=sum(abs(saved)**2);converged=.false.
    do iter=1,state%max_corrector
      midpoint=0.5d0*(saved+previous);h=0
      call midpoint_operator(iter,midpoint,h,op_info)
      if(op_info/=0.or..not.finite2(h))return
      if(maxval(abs(h-conjg(transpose(h))))>100d0*epsilon(1d0)*max(1d0,maxval(abs(h))))return
      call apply_hermitian_exponential(h,state%dt,saved,trial,exp_info)
      if(exp_info/=0.or..not.finite2(trial))return
      delta=sqrt(sum(abs(trial-previous)**2))/max(1d-30,sqrt(sum(abs(trial)**2)))
      iterations=iter
      if(delta<=state%corrector_tolerance)then;converged=.true.;exit;endif
      previous=trial
    enddo
    if(.not.converged)return
    trial_norm=sum(abs(trial)**2);norm_drift=trial_norm-saved_norm
    if(abs(norm_drift)>state%norm_tolerance*max(1d0,saved_norm))return
    state%coeff=trial;state%step=state%step+1;info=0
  end subroutine

  subroutine advance_dg_wpw_fieldoff_exp(state,eigenvalues,iterations,norm_drift,info)
    type(s_dg_wpw_exp_state),intent(inout)::state
    real(8),intent(in)::eigenvalues(:)
    integer,intent(out)::iterations,info
    real(8),intent(out)::norm_drift
    if(size(eigenvalues)/=size(state%coeff,1))then
      iterations=0;norm_drift=huge(1d0);info=1;return
    endif
    call advance_dg_wpw_midpoint_exp(state,fieldoff_operator,iterations,norm_drift,info)
  contains
    subroutine fieldoff_operator(iteration,midpoint_coefficients,h_reduced,op_info)
      integer,intent(in)::iteration
      complex(8),intent(in)::midpoint_coefficients(:,:)
      complex(8),intent(out)::h_reduced(:,:)
      integer,intent(out)::op_info
      integer::i
      h_reduced=0
      do i=1,size(eigenvalues);h_reduced(i,i)=eigenvalues(i);enddo
      op_info=merge(0,1,iteration>0.and.size(midpoint_coefficients,1)==size(eigenvalues))
    end subroutine
  end subroutine

  subroutine apply_hermitian_exponential(h,dt,x,y,info)
    complex(8),intent(in)::h(:,:),x(:,:)
    real(8),intent(in)::dt
    complex(8),intent(out)::y(:,:)
    integer,intent(out)::info
    complex(8),allocatable::vectors(:,:),work(:),rotated(:,:)
    real(8),allocatable::eigenvalues(:),rwork(:)
    integer::n,lwork,lapack_info,i
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character(1),intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*)
        real(8),intent(out)::w(*)
        complex(8),intent(inout)::work(*)
        real(8),intent(inout)::rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface
    info=1;n=size(h,1)
    if(size(h,2)/=n.or.size(x,1)/=n.or.any(shape(y)/=shape(x)))return
    vectors=h;allocate(eigenvalues(n),rwork(max(1,3*n-2)))
    lwork=max(1,2*n-1);allocate(work(lwork))
    call zheev('V','U',n,vectors,n,eigenvalues,work,lwork,rwork,lapack_info)
    if(lapack_info/=0)return
    rotated=matmul(conjg(transpose(vectors)),x)
    do i=1,n;rotated(i,:)=exp(cmplx(0d0,-dt*eigenvalues(i),8))*rotated(i,:);enddo
    y=matmul(vectors,rotated);info=0
  end subroutine
  logical function finite2(a)result(ok)
    complex(8),intent(in)::a(:,:)
    ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
end module
