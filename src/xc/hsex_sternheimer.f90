! Static finite-temperature response for one source band. The caller supplies
! converged target eigenpairs; omitted target bands must be empty and above
! the source energy. The Hamiltonian callback acts on a flattened spinor.
module hsex_sternheimer
  use hsex_screen_math, only: screen_finite
  implicit none
  private
  public :: sternheimer_band

  abstract interface
    subroutine apply_h(x, y, status)
      implicit none
      complex(8), intent(in) :: x(:)
      complex(8), intent(out) :: y(:)
      integer, intent(out) :: status
    end subroutine
  end interface
contains
  subroutine sternheimer_band(apply, target, levels, occupations, energy, occupation, temperature, &
                              vertices, response, status)
    implicit none
    procedure(apply_h) :: apply
    complex(8), intent(in) :: target(:,:), vertices(:,:)
    real(8), intent(in) :: levels(:), occupations(:), energy, occupation, temperature
    complex(8), intent(out) :: response(:,:)
    integer, intent(out) :: status
    complex(8), allocatable :: overlap(:,:), x(:), r(:), p(:), ap(:), rhs(:)
    real(8) :: gap, weight, rr, next_rr, curvature, step, tolerance, midpoint
    integer :: i, j, iteration, ierr

    status=1
    response=0d0
    if (.not.screen_finite([energy,occupation,temperature])) return
    if (temperature<=0d0.or.occupation<0d0.or.occupation>1d0) return
    if (size(target,2)/=size(levels).or.size(occupations)/=size(levels)) return
    if (size(target,1)/=size(vertices,1).or.any(shape(response)/=size(vertices,2))) return
    if (.not.screen_finite(levels).or..not.screen_finite(occupations)) return
    if (any(occupations<0d0).or.any(occupations>1d0)) return
    status=2
    overlap=matmul(conjg(transpose(target)),vertices)
    do i=1,size(levels)
      gap=levels(i)-energy
      if (abs(gap)<=sqrt(epsilon(gap))*temperature) then
        midpoint=0.5d0*(occupation+occupations(i))
        weight=midpoint*(1d0-midpoint)/temperature
      else
        weight=(occupation-occupations(i))/gap
      endif
      ! A second copy of this explicit contribution enters from -q.
      do j=1,size(vertices,2)
        response(:,j)=response(:,j)+0.25d0*weight*conjg(overlap(i,:))*overlap(i,j)
      enddo
    enddo

    if (occupation>0d0.and.size(target,2)<size(target,1)) then
      allocate(x(size(target,1)),r(size(target,1)),p(size(target,1)),ap(size(target,1)),rhs(size(target,1)))
      do j=1,size(vertices,2)
        rhs=vertices(:,j)
        call project(rhs)
        x=0d0
        r=rhs
        p=r
        rr=real(dot_product(r,r),8)
        tolerance=max(1d-10*sqrt(rr),1d-12*norm2(abs(vertices(:,j))))
        do iteration=1,1000
          if (sqrt(rr)<=tolerance) then
            ! Recursive CG residuals can drift from b-Ax. Restart from the
            ! explicitly evaluated residual if the convergence test disagrees.
            call true_residual(ierr)
            if(ierr/=0)return
            if(sqrt(rr)<=tolerance)exit
            p=r
          endif
          status=3
          call apply(p,ap,ierr)
          if (ierr/=0.or..not.screen_finite(ap)) return
          ap=ap-energy*p
          call project(ap)
          curvature=real(dot_product(p,ap),8)
          status=4
          if (curvature<=0d0) return
          step=rr/curvature
          x=x+step*p
          r=r-step*ap
          next_rr=real(dot_product(r,r),8)
          status=5
          if (.not.screen_finite([next_rr])) return
          p=r+(next_rr/rr)*p
          rr=next_rr
        enddo
        call true_residual(ierr)
        if(ierr/=0)return
        if(sqrt(rr)>tolerance)return
        response(:,j)=response(:,j)+0.5d0*occupation*matmul(conjg(transpose(vertices)),x)
      enddo
    endif
    status=7
    if (.not.screen_finite(response)) return
    status=0
  contains
    subroutine true_residual(code)
      implicit none
      integer,intent(out) :: code
      status=3
      call apply(x,ap,code)
      if(code/=0)return
      if(.not.screen_finite(ap))then
        code=1
        return
      endif
      ap=ap-energy*x
      call project(ap)
      r=rhs-ap
      rr=real(dot_product(r,r),8)
      status=6
      if(.not.screen_finite([rr]))code=1
    end subroutine true_residual

    subroutine project(vector)
      implicit none
      complex(8), intent(inout) :: vector(:)
      vector=vector-matmul(target,matmul(conjg(transpose(target)),vector))
    end subroutine
  end subroutine
end module
