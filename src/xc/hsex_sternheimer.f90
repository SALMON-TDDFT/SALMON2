! Static finite-temperature response for one source band. The caller supplies
! converged target eigenpairs; omitted target bands must be empty and above
! the source energy. The Hamiltonian callback acts on a flattened spinor.
module hsex_sternheimer
  implicit none
  private
  public :: sternheimer_band,sternheimer_solve,s_sternheimer_stats
  type :: s_sternheimer_stats
    integer(8) :: rhs=0,iterations=0,h_actions=0
    real(8) :: max_relative_residual=0d0,solve_seconds=0d0
  end type

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
                              vertices, response, status,statistics)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(apply_h) :: apply
    complex(8), intent(in) :: target(:,:), vertices(:,:)
    real(8), intent(in) :: levels(:), occupations(:), energy, occupation, temperature
    complex(8), intent(out) :: response(:,:)
    integer, intent(out) :: status
    type(s_sternheimer_stats),optional,intent(inout) :: statistics
    complex(8), allocatable :: overlap(:,:), x(:)
    real(8) :: gap, weight, midpoint
    integer :: i, j

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
      allocate(x(size(target,1)))
      do j=1,size(vertices,2)
        call sternheimer_solve(apply,target,energy,vertices(:,j),x,status,statistics)
        if(status/=0)return
        response(:,j)=response(:,j)+0.5d0*occupation*matmul(conjg(transpose(vertices)),x)
      enddo
    endif
    status=7
    if (.not.screen_finite(response)) return
    status=0
  end subroutine sternheimer_band

  subroutine sternheimer_solve(apply,target,energy,vertex,x,status,statistics)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(apply_h) :: apply
    complex(8),intent(in) :: target(:,:),vertex(:)
    real(8),intent(in) :: energy
    complex(8),intent(out) :: x(:)
    integer,intent(out) :: status
    type(s_sternheimer_stats),optional,intent(inout) :: statistics
    complex(8) :: r(size(vertex)),p(size(vertex)),ap(size(vertex)),rhs(size(vertex))
    real(8) :: rr,next_rr,curvature,step,tolerance,scale
    integer :: iteration,ierr
    integer(8) :: clock_start,clock_end,clock_rate
    if(present(statistics))call system_clock(clock_start,clock_rate)
    status=1;x=0d0
    if(size(target,1)/=size(vertex).or.size(x)/=size(vertex))return
    if(.not.screen_finite(target).or..not.screen_finite(vertex).or..not.screen_finite([energy]))return
    if(present(statistics))statistics%rhs=statistics%rhs+1
    rhs=vertex;call project(rhs)
    r=rhs;p=r;rr=real(dot_product(r,r),8)
    scale=norm2(abs(vertex))
    tolerance=max(1d-10*sqrt(rr),1d-12*scale)
    do iteration=1,1000
      if(sqrt(rr)<=tolerance)then
        call true_residual(ierr)
        if(ierr/=0)return
        if(sqrt(rr)<=tolerance)exit
        p=r
      endif
      status=3
      call apply(p,ap,ierr)
      if(present(statistics))then
        statistics%iterations=statistics%iterations+1
        statistics%h_actions=statistics%h_actions+1
      endif
      if(ierr/=0.or..not.screen_finite(ap))return
      ap=ap-energy*p;call project(ap)
      curvature=real(dot_product(p,ap),8);status=4
      if(curvature<=0d0)return
      step=rr/curvature;x=x+step*p;r=r-step*ap
      next_rr=real(dot_product(r,r),8);status=5
      if(.not.screen_finite([next_rr]))return
      p=r+(next_rr/rr)*p;rr=next_rr
    enddo
    call true_residual(ierr)
    if(ierr/=0)return
    if(sqrt(rr)>tolerance)return
    if(present(statistics))statistics%max_relative_residual=max(statistics%max_relative_residual, &
      sqrt(rr)/max(scale,tiny(scale)))
    if(present(statistics))then
      call system_clock(clock_end)
      statistics%solve_seconds=statistics%solve_seconds+dble(clock_end-clock_start)/clock_rate
    endif
    status=0
  contains
    subroutine true_residual(code)
      implicit none
      integer,intent(out) :: code
      status=3
      call apply(x,ap,code)
      if(present(statistics))statistics%h_actions=statistics%h_actions+1
      if(code/=0)return
      if(.not.screen_finite(ap))then
        code=1;return
      endif
      ap=ap-energy*x;call project(ap);r=rhs-ap
      rr=real(dot_product(r,r),8);status=6
      if(.not.screen_finite([rr]))code=1
    end subroutine
    subroutine project(vector)
      implicit none
      complex(8),intent(inout) :: vector(:)
      vector=vector-matmul(target,matmul(conjg(transpose(target)),vector))
    end subroutine
  end subroutine sternheimer_solve
end module
