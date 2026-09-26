! Copyright 2026 SALMON developers. Licensed under Apache-2.0.
! Fixed orthonormal complex LCFO basis. This module supplies algebra only;
! a self-consistent SALMON Hamiltonian/current adapter is required for TDHSE.
module lcfo_rt_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_cayley_step,lcfo_density,lcfo_midpoint_step
  public :: lcfo_grid_density,lcfo_project_potential
  abstract interface
    subroutine lcfo_hamiltonian_callback(density,time,h,status)
      complex(8),intent(in) :: density(:,:)
      real(8),intent(in) :: time
      complex(8),intent(out) :: h(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  logical function finite_complex(a) result(ok)
    complex(8),intent(in) :: a(:,:)
    ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function

  subroutine lcfo_density(c,occupation,p,status)
    complex(8),intent(in) :: c(:,:)
    real(8),intent(in) :: occupation(:)
    complex(8),intent(out) :: p(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: weighted(:,:)
    integer :: m,n
    status=1;p=(0d0,0d0);m=size(c,1);n=size(c,2)
    if(m<1.or.n<1.or.size(occupation)/=n.or.any(shape(p)/=[m,m]))return
    if(.not.finite_complex(c).or.any(.not.ieee_is_finite(occupation)))return
    if(any(occupation<0d0))return
    allocate(weighted(m,n))
    weighted=c*spread(occupation,1,m)
    p=matmul(weighted,conjg(transpose(c)))
    if(.not.finite_complex(p))then
      p=0;return
    end if
    status=0
  end subroutine

  subroutine lcfo_cayley_step(h,c,dt,next,status)
    complex(8),intent(in) :: h(:,:),c(:,:)
    real(8),intent(in) :: dt
    complex(8),intent(out) :: next(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: lhs(:,:),rhs(:,:)
    integer,allocatable :: pivot(:)
    real(8) :: scale
    integer :: m,n,i,lapack_info
    external zgesv
    status=1;next=0;m=size(c,1);n=size(c,2)
    if(any(shape(next)/=shape(c)))return
    next=c
    if(m<1.or.n<1.or.any(shape(h)/=[m,m]))return
    if(.not.ieee_is_finite(dt).or..not.finite_complex(c).or..not.finite_complex(h))return
    status=2
    scale=max(1d0,maxval(abs(h)))
    if(maxval(abs(h-conjg(transpose(h))))>1d-12*scale)return
    allocate(lhs(m,m),rhs(m,n),pivot(m))
    lhs=cmplx(0d0,dt/2,8)*h
    rhs=c-matmul(lhs,c)
    do i=1,m
      lhs(i,i)=lhs(i,i)+1d0
    end do
    status=1
    if(.not.finite_complex(lhs).or..not.finite_complex(rhs))return
    call zgesv(m,n,lhs,m,pivot,rhs,m,lapack_info)
    status=3
    if(lapack_info/=0.or..not.finite_complex(rhs))return
    next=rhs;status=0
  end subroutine

  subroutine lcfo_midpoint_step(c,occupation,time,dt,hamiltonian,tolerance,maxiter, &
                               next,residual,iterations,status)
    complex(8),intent(in) :: c(:,:)
    real(8),intent(in) :: occupation(:),time,dt,tolerance
    integer,intent(in) :: maxiter
    procedure(lcfo_hamiltonian_callback) :: hamiltonian
    complex(8),intent(out) :: next(:,:)
    real(8),intent(out) :: residual
    integer,intent(out) :: iterations,status
    complex(8),allocatable :: p0(:,:),guess(:,:),updated(:,:),h(:,:),trial(:,:),p1(:,:)
    integer :: m,n,k,ierr
    real(8) :: scale,midtime
    status=1;iterations=0;residual=huge(1d0);next=0
    if(any(shape(next)/=shape(c)))return
    next=c;m=size(c,1);n=size(c,2)
    if(m<1.or.n<1.or.maxiter<1)return
    if(.not.ieee_is_finite(time).or..not.ieee_is_finite(dt))return
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)return
    midtime=time+dt/2
    if(.not.ieee_is_finite(midtime))return
    allocate(p0(m,m),guess(m,m),updated(m,m),h(m,m),p1(m,m),trial(m,n))
    call lcfo_density(c,occupation,p0,ierr)
    if(ierr/=0)return
    guess=p0;scale=max(1d0,sqrt(sum(abs(p0)**2)))
    do k=1,maxiter
      iterations=k
      call hamiltonian(guess,midtime,h,ierr)
      status=5
      if(ierr/=0)return
      call lcfo_cayley_step(h,c,dt,trial,status)
      if(status/=0)return
      call lcfo_density(trial,occupation,p1,status)
      if(status/=0)return
      updated=0.5d0*(p0+p1)
      residual=sqrt(sum(abs(updated-guess)**2))/scale
      status=1
      if(.not.ieee_is_finite(residual))return
      if(residual<=tolerance)then
        next=trial;status=0;return
      end if
      guess=updated
    end do
    ! Failed fixed point: retain input and let caller reduce dt or abort.
    status=4
  end subroutine
  ! A fragment owns disjoint core-grid basis functions. Its coefficient rows
  ! still refer to ALL global occupied states; do not use fragment occupations.
  subroutine lcfo_grid_density(basis,coeff,occupation,rho,status)
    complex(8),intent(in) :: basis(:,:),coeff(:,:)
    real(8),intent(in) :: occupation(:)
    real(8),intent(out) :: rho(:)
    integer,intent(out) :: status
    complex(8),allocatable :: psi(:,:)
    integer :: ng,nb,no,j
    status=1;rho=0d0;ng=size(basis,1);nb=size(basis,2);no=size(coeff,2)
    if(ng<1.or.no<1.or.size(coeff,1)/=nb.or.size(rho)/=ng.or.size(occupation)/=no)return
    if(.not.finite_complex(basis).or..not.finite_complex(coeff))return
    if(any(.not.ieee_is_finite(occupation)).or.any(occupation<0d0))return
    allocate(psi(ng,no));psi=matmul(basis,coeff)
    do j=1,no
      rho=rho+occupation(j)*abs(psi(:,j))**2
    end do
    if(any(.not.ieee_is_finite(rho)))then
      rho=0d0;return
    end if
    status=0
  end subroutine

  subroutine lcfo_project_potential(basis,potential,dv,h,status)
    complex(8),intent(in) :: basis(:,:)
    real(8),intent(in) :: potential(:),dv
    complex(8),intent(out) :: h(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: weighted(:,:)
    integer :: ng,nb
    status=1;h=0d0;ng=size(basis,1);nb=size(basis,2)
    if(ng<1.or.size(potential)/=ng.or.any(shape(h)/=[nb,nb]))return
    if(.not.ieee_is_finite(dv).or.dv<=0d0)return
    if(.not.finite_complex(basis).or.any(.not.ieee_is_finite(potential)))return
    allocate(weighted(ng,nb));weighted=basis*spread(potential,2,nb)
    h=matmul(conjg(transpose(basis)),weighted)*dv
    if(.not.finite_complex(h))then
      h=0;return
    end if
    status=0
  end subroutine
end module
