program test_frozen_window
  use hsex_frozen_window, only: window_filter
  use hsex_frozen_response, only: frozen_solve,s_frozen_solution,s_frozen_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n=6
  real(8),parameter :: points(9)=[0d0,.1d0,1d0,4d0,7.99d0,8d0,8.01d0,10d0,100d0]
  real(8),parameter :: dawson(9)=[0d0,.09933599239785286115d0,.53807950691276841914d0, &
    .12934800123600511559d0,.06308033315699247654d0,.06300019870755338792d0, &
    .06292026930563153354d0,.05025384718759852803d0,.00500025003750937828d0]
  complex(8) :: value,expected,op(n,n),left(n,n),right(n,n),hl(n,n),hr(n,n),weighted(n,n)
  complex(8) :: u(n,2),v(n,2),c(2,2),extra(2,2,2),z(3),reference(n,n,3),rhs(n,n),actual(n,n),filters(n,n)
  real(8) :: el(n),er(n),times(3),area,pi,delta,error,source_norm,time_error
  integer :: kinds(3),i,j,k,status,capacity
  type(s_frozen_solution) :: solution
  type(s_frozen_stats) :: stats
  pi=acos(-1d0);area=sqrt(pi)
  do i=1,9
    call window_filter(1,area,points(i),value,status)
    expected=cmplx(-2d0*dawson(i),-area*exp(-points(i)**2),8)
    if(status/=0.or.abs(value-expected)>5d-14*max(1d0,abs(expected)))error stop 'Dawson reference mismatch'
    call window_filter(1,area,-points(i),value,status)
    if(status/=0.or.abs(value+conjg(expected))>5d-14*max(1d0,abs(expected)))error stop 'window parity mismatch'
  enddo
  area=20d0
  call window_filter(1,area,0d0,value,status)
  if(status/=0.or.abs(value-cmplx(0d0,-area,8))>1d-13)error stop 'Gaussian finite zero mode'
  call window_filter(2,area,0d0,value,status)
  if(status/=0.or.abs(value-cmplx(0d0,-area,8))>1d-13)error stop 'equal-area zero mode'
  do i=1,4
    delta=.07d0*i
    call window_filter(1,area,delta,value,status)
    call time_reference(1,area,delta,expected)
    if(status/=0.or.abs(value-expected)>2d-10)error stop 'Gaussian time integral mismatch'
    call window_filter(2,area,delta,value,status)
    if(status/=0.or.abs(value-1d0/(cmplx(0d0,1d0/area,8)-delta))>1d-13) &
      error stop 'exponential resolvent mismatch'
    call time_reference(2,area,delta,expected)
    if(abs(value-expected)>2d-9)error stop 'exponential time integral mismatch'
  enddo
  call window_filter(1,-1d0,0d0,value,status)
  if(status==0)error stop 'negative time accepted'
  call window_filter(1,ieee_value(0d0,ieee_quiet_nan),0d0,value,status)
  if(status==0)error stop 'NaN time accepted'
  call window_filter(1,1d0,ieee_value(0d0,ieee_positive_inf),value,status)
  if(status==0)error stop 'infinite energy accepted'
  call window_filter(9,1d0,0d0,value,status)
  if(status==0)error stop 'unknown filter accepted'
  do j=1,n
    do i=1,n
      left(i,j)=exp(cmplx(0d0,2d0*pi*(i-1)*(j-1)/n,8))/sqrt(dble(n))
      right(i,j)=left(i,j)*exp(cmplx(0d0,.17d0*i,8))
    enddo
    el(j)=.13d0*j;er(j)=.11d0*j+.03d0
    weighted(:,j)=left(:,j)*el(j)
  enddo
  hl=matmul(weighted,conjg(transpose(left)))
  do j=1,n
    weighted(:,j)=right(:,j)*er(j)
  enddo
  hr=matmul(weighted,conjg(transpose(right)))
  u=0d0;v=0d0;u(1,1)=1d0;u(3,2)=1d0;v(2,1)=1d0;v(5,2)=cmplx(.7d0,.2d0,8)
  c=reshape([cmplx(1d0,.2d0,8),cmplx(.3d0,-.1d0,8),cmplx(-.1d0,.4d0,8),cmplx(.8d0,0d0,8)],[2,2])
  kinds=[0,1,2];times=[0d0,10d0,10d0];z=cmplx(0d0,.1d0,8)
  op=matmul(u,matmul(c,conjg(transpose(v))));source_norm=sqrt(sum(abs(op)**2))
  rhs=matmul(conjg(transpose(left)),matmul(op,right))
  do k=1,3
    do j=1,n
      do i=1,n
        delta=el(i)-er(j)
        if(k==1)then
          filters(i,j)=1d0/(z(k)-delta)
        else
          call window_filter(kinds(k),times(k),delta,filters(i,j),status)
          if(status/=0)error stop 'spectral reference filter'
        endif
      enddo
    enddo
    reference(:,:,k)=matmul(left,matmul(rhs*filters,conjg(transpose(right))))
  enddo
  do capacity=2,n,4
    call frozen_solve(action_left,action_right,u,c,v,z,capacity,1d-9,solution,stats,status, &
      residual_bound=.true.,window_kind=kinds,window_time=times)
    if(capacity==n.and.status/=0)error stop 'full window solution rejected'
    if(capacity<n.and.status/=2)error stop 'incomplete window space not rejected'
    do k=2,3
      actual=matmul(solution%left,matmul(solution%coefficient(:,:,k),conjg(transpose(solution%right))))
      error=sqrt(sum(abs(actual-reference(:,:,k))**2))
      if(error>solution%residual(k)*solution%error_scale(k)+1d-12)error stop 'window bound underestimated error'
      if(solution%window_kind(k)/=kinds(k))error stop 'missing window error marker'
      if(abs(solution%error_scale(k)-times(k)*source_norm)>1d-12)error stop 'incorrect window normalization'
      if(capacity==n.and.error>1d-10)error stop 'window dense oracle mismatch'
    enddo
    if(capacity==n)then
      if(maxval(abs(solution%coefficient(:,:,1)-solution%coefficient(:,:,3)))>1d-11) &
        error stop 'mixed exponential/resolvent mismatch'
    endif
  enddo
  extra(:,:,1)=.7d0*c;extra(:,:,2)=0d0
  call frozen_solve(action_left,action_right,u,c,v,z,n,1d-9,solution,stats,status, &
    residual_bound=.true.,window_kind=kinds,window_time=times,extra_sources=extra)
  if(status/=0.or.stats%h_actions/=2*n)error stop 'mixed window sources did not share bases'
  if(size(solution%error_scale)/=9)error stop 'mixed error marker source ordering'
  if(maxval(abs(solution%coefficient(:,:,4:6)-.7d0*solution%coefficient(:,:,:3)))>1d-11) &
    error stop 'mixed window source scaling'
  if(maxval(abs(solution%coefficient(:,:,7:9)))>0d0.or.maxval(solution%error_scale(7:9))>0d0) &
    error stop 'zero extra source produced response'
  call frozen_solve(action_left,action_right,u,c,v,z,n,1d-9,solution,stats,status, &
    window_kind=kinds,window_time=times)
  if(status==0)error stop 'uncertified window mode accepted'
  call frozen_solve(action_left,action_right,u,c,v,z,n,1d-9,solution,stats,status, &
    residual_bound=.true.,window_kind=kinds)
  if(status==0)error stop 'partial window arguments accepted'
  times(2)=ieee_value(0d0,ieee_quiet_nan)
  call frozen_solve(action_left,action_right,u,c,v,z,n,1d-9,solution,stats,status, &
    residual_bound=.true.,window_kind=kinds,window_time=times)
  if(status==0)error stop 'nonfinite window solve accepted'
  times=[0d0,10d0,10d0];c=0d0
  call frozen_solve(action_left,action_right,u,c,v,z,n,1d-9,solution,stats,status, &
    residual_bound=.true.,window_kind=kinds,window_time=times)
  if(status/=0.or.stats%h_actions/=0.or.maxval(solution%residual)/=0d0)error stop 'zero window source'
  if(any(solution%error_scale/=0d0))error stop 'zero window source normalization'
  print *,'FROZEN_WINDOW_QUADRATURE_SPECTRAL_BOUND_IEEE_PASS'
contains
  subroutine action_left(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(hl,x);code=0
  end subroutine
  subroutine action_right(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(hr,x);code=0
  end subroutine
  subroutine time_reference(kind,a,d,result)
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: a,d
    complex(8),intent(out) :: result
    integer,parameter :: intervals=200000
    real(8) :: t,dt,weight,limit,sigma
    integer :: j,factor
    sigma=a*sqrt(2d0/pi)
    limit=12d0*sigma
    if(kind==2)limit=40d0*a
    dt=limit/intervals;result=0d0
    do j=0,intervals
      t=j*dt;factor=2
      if(mod(j,2)==1)factor=4
      if(j==0.or.j==intervals)factor=1
      if(kind==1)then
        weight=exp(-t*t/(2d0*sigma*sigma))
      else
        weight=exp(-t/a)
      endif
      result=result+factor*weight*exp(cmplx(0d0,-d*t,8))
    enddo
    result=cmplx(0d0,-1d0,8)*result*dt/3d0
  end subroutine
end program
