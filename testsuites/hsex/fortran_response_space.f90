program test_response_space
  use hsex_response_space, only: response_space_solve,s_response_space_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n=24,nrhs=8
  complex(8) :: target(n,2),rhs(n,nrhs),x(n,nrhs),exact(n,nrhs),z(n),h(n,n)
  complex(8) :: projector(n,n),operator(n,n),dense(n,n),column(n,1),rotation(n,n)
  complex(8) :: defect(n,nrhs),corrected(n,nrhs),difference(nrhs,nrhs),bound(nrhs,nrhs)
  real(8) :: shifts(nrhs),residual(nrhs),tolerance
  integer :: lapack_status
  type(s_response_space_stats) :: stats
  integer :: i,j,code
  target=0d0;target(1,1)=1d0;target(2,2)=1d0;h=0d0;rhs=0d0;exact=0d0
  do i=1,n
    h(i,i)=1d0+.1d0*i
  enddo
  do j=1,nrhs
    shifts(j)=-.2d0*j
    rhs(1,j)=.3d0 ! must be projected away
    do i=3,6
      rhs(i,j)=cmplx(sin(dble(i))*j,cos(dble(i))*(1+mod(j,2)),8)
      exact(i,j)=rhs(i,j)/(real(h(i,i),8)-shifts(j))
    enddo
  enddo
  tolerance=1d-10
  call response_space_solve(apply,target,shifts,rhs,12,tolerance,x,residual,stats,code)
  if(code/=0.or.maxval(abs(x-exact))>1d-10)error stop 'shifted dependent RHS solution'
  if(stats%dimension>4.or.stats%h_actions/=stats%dimension+nrhs)error stop 'response reuse accounting'
  if(maxval(residual)>tolerance)error stop 'true residual'
  print *,'RESPONSE_SPACE_DEPENDENT_SHIFT_PASS',stats%dimension,stats%h_actions
  call response_space_solve(apply,target,shifts,rhs,1,tolerance,x,residual,stats,code)
  if(code==0)error stop 'capacity exhaustion accepted'
  print *,'RESPONSE_SPACE_CAPACITY_REJECTION_PASS'
  do j=1,2
    if(j==1)rhs(n,nrhs)=ieee_value(0d0,ieee_quiet_nan)
    if(j==2)rhs(n,nrhs)=ieee_value(0d0,ieee_positive_inf)
    call response_space_solve(apply,target,shifts,rhs,12,tolerance,x,residual,stats,code)
    if(code==0)error stop 'nonfinite RHS accepted'
  enddo
  print *,'RESPONSE_SPACE_IEEE_REJECTION_PASS'
  ! A complex Hamiltonian whose occupied projector does not commute with H.
  ! Compare Q(H-e)Q inverse with an independent positive-definite dense solve.
  rotation=0d0
  do i=1,n
    rotation(i,i)=1d0
  enddo
  rotation(1,1)=.8d0;rotation(3,3)=.8d0
  rotation(3,1)=(0d0,.6d0);rotation(1,3)=(0d0,.6d0)
  target=rotation(:,:2);h=0d0;projector=0d0
  do i=1,n
    h(i,i)=1d0+.1d0*i;projector(i,i)=1d0
    if(i<n)then
      h(i,i+1)=(.07d0,.04d0);h(i+1,i)=(.07d0,-.04d0)
    endif
  enddo
  projector=projector-matmul(target,conjg(transpose(target)))
  operator=matmul(projector,matmul(h,projector))
  do j=1,nrhs
    do i=1,n
      rhs(i,j)=cmplx(sin(dble(i*j)),cos(dble(i*(j+1))),8)
    enddo
    dense=operator
    do i=1,n
      dense(i,i)=dense(i,i)-shifts(j)
    enddo
    column(:,1)=matmul(projector,rhs(:,j))
    call zposv('L',n,1,dense,n,column,n,lapack_status)
    if(lapack_status/=0)error stop 'independent projected reference failed'
    exact(:,j)=column(:,1)
  enddo
  call response_space_solve(apply,target,shifts,rhs,n-2,tolerance,x,residual,stats,code)
  if(code/=0.or.maxval(abs(x-exact))>1d-10.or.maxval(residual)>tolerance) &
    error stop 'noncommuting complex projector response mismatch'
  print *,'RESPONSE_SPACE_COMPLEX_PROJECTOR_PASS',stats%dimension,maxval(abs(x-exact))
  shifts=-.2d0
  do j=1,nrhs
    dense=operator
    do i=1,n
      dense(i,i)=dense(i,i)-shifts(j)
    enddo
    column(:,1)=matmul(projector,rhs(:,j))
    call zposv('L',n,1,dense,n,column,n,lapack_status)
    if(lapack_status/=0)error stop 'Galerkin identity reference failed'
    exact(:,j)=column(:,1)
  enddo
  call response_space_solve(apply,target,shifts,rhs,1,tolerance,x,residual,stats,code)
  if(code/=2)error stop 'Galerkin identity needs a deliberately truncated space'
  defect=matmul(projector,rhs)-matmul(operator,x)+shifts(1)*x
  do j=1,nrhs
    dense=operator
    do i=1,n
      dense(i,i)=dense(i,i)-shifts(j)
    enddo
    column(:,1)=defect(:,j)
    call zposv('L',n,1,dense,n,column,n,lapack_status)
    if(lapack_status/=0)error stop 'Galerkin identity correction failed'
    corrected(:,j)=column(:,1)
  enddo
  difference=matmul(conjg(transpose(rhs)),exact-x)
  bound=matmul(conjg(transpose(defect)),corrected)
  if(maxval(abs(difference-bound))>1d-10)error stop 'quadratic response error identity'
  print *,'RESPONSE_SPACE_QUADRATIC_OBSERVABLE_ERROR_PASS',maxval(abs(difference-bound))
contains
  subroutine apply(v,y,status)
    implicit none
    complex(8),intent(in) :: v(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    y=matmul(h,v);status=0
  end subroutine
end program
