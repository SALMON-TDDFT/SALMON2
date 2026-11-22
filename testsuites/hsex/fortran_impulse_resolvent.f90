program test_impulse_resolvent
  use hsex_impulse_resolvent, only: impulse_model,impulse_build,impulse_resolve,impulse_trajectory
  use communication, only: comm_init,comm_finalize
  use mpi
  implicit none
  type(impulse_model) :: model,repeated
  complex(8) :: source(1),s(3),value(3),wanted(3)
  real(8) :: residual(3),frequency
  real(8) :: current(0:2000),exact
  real(8) :: original(4,4),rotation(4,4),inverse(4,4)
  integer :: j
  integer :: status,rank,ierr
  call comm_init
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  frequency=.13d0;source=cmplx(0d0,1d0,8)
  s=[cmplx(.002d0,-.08d0,8),cmplx(.002d0,-.13d0,8),cmplx(.002d0,-.18d0,8)]
  call impulse_build(action,observation,source,s,4,1d-10,MPI_COMM_WORLD,rank,model,status)
  if(status/=0)error stop 'oscillator Arnoldi build'
  call impulse_resolve(model,s,value,residual,status)
  wanted=1d0/(s*s+frequency**2)
  if(status/=0.or.maxval(abs(value-wanted))/maxval(abs(wanted))>1d-10)error stop 'oscillator Sternheimer'
  if(maxval(residual)>1d-10)error stop 'oscillator full residual'
  call impulse_trajectory(model,0d0,.03d0,current,status)
  if(status/=0)error stop 'oscillator finite-window companion'
  do j=0,2000
    exact=(1d0-cos(frequency*j*.03d0))/frequency**2
    if(abs(current(j)-exact)>1d-8)error stop 'oscillator time/spectral mismatch'
  enddo
  call impulse_build(action,observation,source,s,8,1d-9,MPI_COMM_WORLD,rank,model,status,precondition=preconditioner)
  if(status/=0)error stop 'preconditioned occupied response'
  call impulse_resolve(model,s,value,residual,status)
  if(status/=0.or.maxval(abs(value-wanted))/maxval(abs(wanted))>1d-9)error stop 'rational response mismatch'
  if(maxval(residual)>1d-9)error stop 'rational true residual'
  allocate(repeated%h(4,4),repeated%g(4));original=0d0;rotation=0d0;inverse=0d0
  do j=1,4
    rotation(j,j)=1d0;inverse(j,j)=1d0
  enddo
  original(1,2)=1d0;original(2,1)=-frequency**2
  original(3,4)=1d0;original(4,3)=-frequency**2
  rotation(3,1)=.3d0;rotation(4,2)=.7d0;inverse(3,1)=-.3d0;inverse(4,2)=-.7d0
  repeated%h=matmul(rotation,matmul(original,inverse));repeated%g=[1d0,.2d0,-.4d0,.3d0]
  repeated%source_norm=1d0;repeated%next_norm=0d0
  call impulse_trajectory(repeated,0d0,.03d0,current,status)
  if(status/=0)error stop 'degenerate nonnormal spectral reconstruction'
  do j=0,2000
    exact=repeated%g(1)*sin(frequency*j*.03d0)/frequency+ &
      dot_product(repeated%g,repeated%h(:,1))*(1d0-cos(frequency*j*.03d0))/frequency**2
    if(abs(current(j)-exact)>1d-8)error stop 'degenerate spectral/time mismatch'
  enddo
  if(rank==0)print *, 'IMPULSE_RESOLVENT_OSCILLATOR_PASS'
  call comm_finalize
contains
  subroutine preconditioner(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    y=cmplx(aimag(x)/frequency**2,-real(x,8),8);status=0
  end subroutine

  subroutine action(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    y=cmplx(aimag(x),-frequency**2*real(x,8),8);status=0
  end subroutine
  subroutine observation(x,value,status)
    implicit none
    complex(8),intent(in) :: x(:)
    real(8),intent(out) :: value
    integer,intent(out) :: status
    value=real(x(1),8);status=0
  end subroutine
end program
