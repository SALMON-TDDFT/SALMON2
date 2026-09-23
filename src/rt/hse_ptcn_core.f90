! Parallel-transport Crank-Nicolson with fixed-ACE inner and fresh-EXX outer loops.
module hse_ptcn_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: hse_ptcn_solve
  abstract interface
    subroutine action_callback(a,b,refresh)
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
      logical,intent(in) :: refresh
    end subroutine
    subroutine precondition_callback(a,b)
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
    end subroutine
    real(8) function norm_callback(a)
      complex(8),intent(in) :: a(:,:,:)
    end function
  end interface
contains
  subroutine hse_ptcn_solve(u,dt,dv,action,precondition,norm,x,error,builds,apps,ierr)
    complex(8),intent(in) :: u(:,:,:)
    real(8),intent(in) :: dt,dv
    procedure(action_callback) :: action
    procedure(precondition_callback) :: precondition
    procedure(norm_callback) :: norm
    complex(8),intent(out) :: x(:,:,:)
    real(8),intent(out) :: error
    integer,intent(out) :: builds,apps,ierr
    complex(8),allocatable :: hu(:,:,:),rhs(:,:,:),r(:,:,:),correction(:,:,:)
    real(8) :: scale,errors(150)
    integer :: outer,inner
    ierr=1;builds=0;apps=0;error=huge(1d0)
    if(dt<=0.or.dv<=0.or..not.ieee_is_finite(dt).or..not.ieee_is_finite(dv))return
    if(any(shape(u)/=shape(x)))return
    scale=norm(u);if(scale<=0.or..not.ieee_is_finite(scale))return
    allocate(hu(size(u,1),size(u,2),size(u,3)),rhs(size(u,1),size(u,2),size(u,3)), &
      r(size(u,1),size(u,2),size(u,3)),correction(size(u,1),size(u,2),size(u,3)))
    x=u
    call action(u,hu,.true.);builds=1
    call pt_gradient(u,hu,dv,r)
    rhs=u-(0d0,.5d0)*dt*r
    do outer=1,12
      do inner=1,150
        call action(x,hu,.false.);apps=apps+1
        call pt_gradient(x,hu,dv,r)
        r=x+(0d0,.5d0)*dt*r-rhs
        error=norm(r)/scale;errors(inner)=error
        if(.not.ieee_is_finite(error))return
        if(error<1d-12)exit
        if(inner>=9.and.error<.5d-10)then
          if(minval(errors(inner-3:inner))>.5d0*minval(errors(inner-7:inner-4)))exit
        endif
        call precondition(r,correction)
        x=x-correction
      enddo
      if(inner>150)return
      call action(x,hu,.true.);builds=builds+1
      call pt_gradient(x,hu,dv,r)
      r=x+(0d0,.5d0)*dt*r-rhs
      error=norm(r)/scale
      if(.not.ieee_is_finite(error))return
      if(error<1d-10)then
        ierr=0;return
      endif
    enddo
  end subroutine

  subroutine pt_gradient(u,hu,dv,r)
    complex(8),intent(in) :: u(:,:,:),hu(:,:,:)
    complex(8),intent(out) :: r(:,:,:)
    real(8),intent(in) :: dv
    complex(8) :: projection(size(u,2),size(u,2))
    integer :: ik,ng,no
    complex(8),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    external :: zgemm
    ng=size(u,1);no=size(u,2);r=hu
    do ik=1,size(u,3)
      call zgemm('C','N',no,no,ng,one*dv,u(1,1,ik),ng,hu(1,1,ik),ng,zero,projection(1,1),no)
      call zgemm('N','N',ng,no,no,-one,u(1,1,ik),ng,projection(1,1),no,one,r(1,1,ik),ng)
    enddo
  end subroutine
end module
