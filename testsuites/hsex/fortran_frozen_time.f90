! Independent orbital-time reference, including evolving-background Hartree RPA.
program test_frozen_time
  use hsex_frozen_response, only: frozen_source,frozen_solve,frozen_trace,frozen_background_rate
  use hsex_frozen_response, only: s_frozen_solution,s_frozen_stats
  implicit none
  integer,parameter :: n=4,ns=2
  complex(8) :: h(n,n),vectors(n,n),psi(n,1),hp(n,1),d(n),z(ns),chi(ns),td(ns),rpa(ns),fine(ns),ws(ns)
  complex(8),allocatable :: u(:,:),v(:,:),c(:,:)
  real(8) :: occ(1),rate,error,drift,g,energies(n)
  type(s_frozen_solution) :: solution
  type(s_frozen_stats) :: statistics
  integer :: a,b,coherent,status
  vectors(:,1)=[1d0,1d0,1d0,1d0]/2d0
  vectors(:,2)=[1d0,-1d0,1d0,-1d0]/2d0
  vectors(:,3)=[1d0,1d0,-1d0,-1d0]/2d0
  vectors(:,4)=[1d0,-1d0,-1d0,1d0]/2d0
  energies=[-.6d0,-.1d0,.5d0,1d0];h=0d0
  do a=1,n
    do b=1,n
      h(:,a)=h(:,a)+energies(b)*vectors(:,b)*conjg(vectors(a,b))
    enddo
  enddo
  d=cmplx([1d0,-.3d0,.7d0,-.5d0],0d0,8);occ=2d0
  z=[cmplx(0d0,.16d0,8),cmplx(0d0,.08d0,8)]
  do coherent=0,1
    psi(:,1)=vectors(:,1)
    if(coherent==1)psi(:,1)=cos(.3d0)*vectors(:,1)+cmplx(0d0,sin(.3d0),8)*vectors(:,3)
    hp=matmul(h,psi)
    call frozen_background_rate(psi,hp,rate,status)
    if(status/=0)error stop 'background diagnostic rejected'
    if(coherent==0.and.rate>1d-13)error stop 'stationary background rate'
    if(coherent==1.and.rate<.1d0)error stop 'coherent background missed'
    call frozen_source(psi,psi,occ,occ,d,u,v,c,status)
    call frozen_solve(action,action,u,c,v,z,n,1d-11,solution,statistics,status)
    if(status/=0)error stop 'time reference solver failed'
    call frozen_trace(solution,d,chi,status)
    call propagate(0d0,.01d0,td,drift)
    call propagate(0d0,.005d0,fine,drift)
    error=maxval(abs(fine-chi))
    if(error>2d-7.or.maxval(abs(fine-td))>2d-7)error stop 'independent orbital impulse mismatch'
    if(drift>1d-8)error stop 'reference norm drift'
    g=.18d0;ws=g/(1d0-g*chi)
    call propagate(g,.005d0,rpa,drift)
    rpa=g+g*g*rpa
    error=maxval(abs(rpa-ws))
    if(coherent==0.and.error>2d-7)error stop 'stationary snapshot RPA mismatch'
    if(coherent==1.and.error<1d-5)error stop 'model fails to distinguish evolving and frozen RPA'
    print '(a,i2,3es17.8)','FROZEN_TIME coherent/chi_error/W_adiabatic_difference/background_rate ', &
      coherent,maxval(abs(fine-chi)),error,rate
  enddo
  print *,'FROZEN_IMPULSE_AND_STATIONARY_RPA_PASS_COHERENT_DIFFERENCE_MEASURED'
contains
  subroutine action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(h,x);code=0
  end subroutine

  subroutine derivative(state,coupling,rhs)
    implicit none
    complex(8),intent(in) :: state(n,2)
    real(8),intent(in) :: coupling
    complex(8),intent(out) :: rhs(n,2)
    real(8) :: density
    density=2d0*occ(1)*real(sum(conjg(state(:,1))*d*state(:,2)),8)
    rhs(:,1)=cmplx(0d0,-1d0,8)*matmul(h,state(:,1))
    rhs(:,2)=cmplx(0d0,-1d0,8)*(matmul(h,state(:,2))+coupling*d*density*state(:,1))
  end subroutine

  subroutine propagate(coupling,dt,values,norm_drift)
    implicit none
    real(8),intent(in) :: coupling,dt
    complex(8),intent(out) :: values(ns)
    real(8),intent(out) :: norm_drift
    complex(8) :: state(n,2),k1(n,2),k2(n,2),k3(n,2),k4(n,2)
    real(8) :: time,density,weight
    integer :: step,steps
    steps=nint(300d0/dt)
    if(mod(steps,2)/=0)error stop 'Simpson grid must be even'
    state(:,1)=psi(:,1);state(:,2)=cmplx(0d0,-1d0,8)*d*psi(:,1);values=0d0
    do step=0,steps
      time=step*dt
      weight=2d0
      if(mod(step,2)==1)weight=4d0
      if(step==0.or.step==steps)weight=1d0
      density=2d0*occ(1)*real(sum(conjg(state(:,1))*d*state(:,2)),8)
      values=values+weight*exp(cmplx(0d0,1d0,8)*z*time)*density
      if(step==steps)exit
      call derivative(state,coupling,k1)
      call derivative(state+dt*k1/2d0,coupling,k2)
      call derivative(state+dt*k2/2d0,coupling,k3)
      call derivative(state+dt*k3,coupling,k4)
      state=state+dt*(k1+2d0*k2+2d0*k3+k4)/6d0
    enddo
    values=values*dt/3d0
    norm_drift=abs(sum(abs(state(:,1))**2)-1d0)
  end subroutine
end program
