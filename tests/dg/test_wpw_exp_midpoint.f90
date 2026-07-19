program test_wpw_exp_midpoint
  use rt_dg_wpw_exp_production,only:s_dg_wpw_exp_state,initialize_dg_wpw_exp_state,advance_dg_wpw_midpoint_exp,&
    advance_dg_wpw_length_gauge_exp
  implicit none
  type(s_dg_wpw_exp_state)::state
  complex(8)::initial(2,1),expected(2,1),saved(2,1)
  integer::info,iterations
  real(8)::norm_drift
  initial(:,1)=[(1d0,0d0),(0d0,0d0)]
  call initialize_dg_wpw_exp_state(state,initial,0.2d0,12,1d-13,1d-12,info)
  if(info/=0)error stop 1
  call advance_dg_wpw_midpoint_exp(state,static_operator,iterations,norm_drift,info)
  expected(:,1)=[exp(cmplx(0d0,-0.2d0,8)),(0d0,0d0)]
  if(info/=0.or.maxval(abs(state%coeff-expected))>1d-12.or.abs(norm_drift)>1d-12)error stop 2
  saved=state%coeff
  call advance_dg_wpw_midpoint_exp(state,nonlinear_operator,iterations,norm_drift,info)
  if(info/=0.or.iterations<2.or.abs(sum(abs(state%coeff)**2)-1d0)>1d-12)error stop 3
  saved=state%coeff
  call advance_dg_wpw_midpoint_exp(state,never_converges,iterations,norm_drift,info)
  if(info==0.or.maxval(abs(state%coeff-saved))>0d0)error stop 4
  call advance_dg_wpw_midpoint_exp(state,nonhermitian_operator,iterations,norm_drift,info)
  if(info==0.or.maxval(abs(state%coeff-saved))>0d0)error stop 5
  call initialize_dg_wpw_exp_state(state,initial,0.2d0,12,1d-13,1d-12,info)
  call advance_dg_wpw_length_gauge_exp(state,[1d0,2d0],reshape([(0.5d0,0d0),(0d0,0d0),&
    (0d0,0d0),(-0.5d0,0d0)],[2,2]),0.4d0,iterations,norm_drift,info)
  expected(:,1)=[exp(cmplx(0d0,-0.24d0,8)),(0d0,0d0)]
  if(info/=0.or.maxval(abs(state%coeff-expected))>1d-12)error stop 6
  print '(a)','PASS saved-Cn midpoint Exp convergence, unitarity, and rollback'
contains
  subroutine static_operator(iter,mid,h,op_info)
    integer,intent(in)::iter;complex(8),intent(in)::mid(:,:);complex(8),intent(out)::h(:,:);integer,intent(out)::op_info
    h=0;h(1,1)=1d0;h(2,2)=2d0;op_info=0
  end subroutine
  subroutine nonlinear_operator(iter,mid,h,op_info)
    integer,intent(in)::iter;complex(8),intent(in)::mid(:,:);complex(8),intent(out)::h(:,:);integer,intent(out)::op_info
    h=0;h(1,1)=1d0+0.05d0*abs(mid(1,1))**2;h(2,2)=2d0;op_info=0
  end subroutine
  subroutine never_converges(iter,mid,h,op_info)
    integer,intent(in)::iter;complex(8),intent(in)::mid(:,:);complex(8),intent(out)::h(:,:);integer,intent(out)::op_info
    h=0;h(1,1)=merge(1d0,2d0,mod(iter,2)==0);h(2,2)=2d0;op_info=0
  end subroutine
  subroutine nonhermitian_operator(iter,mid,h,op_info)
    integer,intent(in)::iter;complex(8),intent(in)::mid(:,:);complex(8),intent(out)::h(:,:);integer,intent(out)::op_info
    h=0;h(1,2)=1d0;op_info=0
  end subroutine
end program
