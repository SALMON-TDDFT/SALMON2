! Adaptive periodic support of spatially distributed localized orbitals.
! Each comm_r group must call collectively with identical n, h, fraction and
! source-column ordering/count; lo/m partition the global mesh without overlap.
! Only scalar/small-vector reductions are used, never a source or grid gather.
module exx_adaptive_support
  implicit none
  private
  public :: adaptive_source_mask
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
contains
  subroutine adaptive_source_mask(n,h,lo,m,comm_r,source,fraction,radii,loss,protected,status,fixed_radius)
    use communication, only: comm_summation,comm_get_max
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: n(3),lo(3),m(3),comm_r
    real(8),intent(in) :: h(3),fraction
    real(8),intent(in),optional :: fixed_radius
    complex(8),intent(inout) :: source(:,:)
    real(8),intent(out) :: radii(:),loss(:)
    logical,intent(out) :: protected(:)
    integer,intent(out) :: status
    real(8),allocatable :: weight(:),distance2(:),point(:,:)
    real(8) :: length(3),center(3),delta(3),local_moment(6),moment(6)
    real(8) :: norm,local_norm,kept,lower,upper,middle,pi,angle,tol
    real(8) :: requested_radius
    integer :: bad,no,g,j,x,y,z,axis,iter
    status=1;bad=0;radii=0d0;loss=0d0;protected=.false.
    requested_radius=0d0
    if(present(fixed_radius))requested_radius=fixed_radius
    if(.not.ieee_is_finite(requested_radius).or.requested_radius<0d0)bad=1
    no=size(source,2)
    if(any(n<1).or.any(m<0).or.any(lo<0).or.any(lo+m>n))bad=1
    if(size(source,1)/=product(m))bad=1
    if(size(radii)/=no.or.size(loss)/=no.or.size(protected)/=no)bad=1
    if(any(h<=0d0).or..not.salmon_all_finite(h))bad=1
    if(.not.ieee_is_finite(fraction).or.fraction<=0d0.or.fraction>1d0)bad=1
    if(.not.salmon_all_finite(real(source)).or..not.salmon_all_finite(aimag(source)))bad=1
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    length=n*h
    if(.not.salmon_all_finite(length).or..not.ieee_is_finite(sum(length**2)))bad=1
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    radii=sqrt(sum((length/2d0)**2))
    ! Exact full-support mode: do not alter even subnormal source coefficients.
    if(fraction==1d0.and.requested_radius==0d0)then
      status=0
      return
    endif
    allocate(weight(size(source,1)),distance2(size(source,1)),point(3,size(source,1)))
    pi=acos(-1d0);g=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      g=g+1;point(:,g)=([x,y,z]+lo)*h
    enddo;enddo;enddo
    do j=1,no
      weight=abs(source(:,j))**2
      local_norm=sum(weight)
      if(.not.ieee_is_finite(local_norm))bad=1
      call comm_get_max(bad,comm_r)
      if(bad/=0)return
      call comm_summation(local_norm,norm,comm_r)
      if(.not.ieee_is_finite(norm))return
      if(norm<=tiny(1d0))then
        protected(j)=.true.
        cycle
      endif
      local_moment=0d0
      do g=1,size(weight);do axis=1,3
        angle=2d0*pi*point(axis,g)/length(axis)
        local_moment(axis)=local_moment(axis)+weight(g)*cos(angle)
        local_moment(axis+3)=local_moment(axis+3)+weight(g)*sin(angle)
      enddo;enddo
      call comm_summation(local_moment,moment,6,comm_r)
      moment=moment/norm
      if(any(sqrt(moment(:3)**2+moment(4:)**2)<.1d0))then
        ! A circular center is unreliable for delocalized or multiple-lobe WFs.
        protected(j)=.true.
        cycle
      endif
      center=modulo(atan2(moment(4:),moment(:3))*length/(2d0*pi),length)
      do g=1,size(weight)
        delta=modulo(point(:,g)-center+length/2d0,length)-length/2d0
        distance2(g)=sum(delta**2)
      enddo
      if(requested_radius>0d0)then
        ! Fixed R takes precedence over the diagnostic norm target. Limit the
        ! square to the largest periodic distance to avoid overflow for huge R.
        upper=min(requested_radius,sqrt(sum((length/2d0)**2)))**2
      else
      lower=0d0;upper=sum((length/2d0)**2)
      ! upper always encloses at least fraction of the norm. Bisection in r^2
      ! locates the first retained mesh shell without storing global distances.
      do iter=1,60
        middle=lower+(upper-lower)/2d0
        local_norm=sum(weight,mask=distance2<=middle)
        call comm_summation(local_norm,kept,comm_r)
        if(kept>=fraction*norm)then
          upper=middle
        else
          lower=middle
        endif
      enddo
      endif
      ! Expand by roundoff only, to keep equivalent boundary shells consistent
      ! across MPI decompositions and the subsequent square-root conversion.
      tol=64d0*epsilon(1d0)*max(1d0,upper)
      upper=upper+tol
      local_norm=sum(weight,mask=distance2<=upper)
      call comm_summation(local_norm,kept,comm_r)
      if(requested_radius==0d0.and.kept< fraction*norm-64d0*epsilon(1d0)*norm)then
        protected(j)=.true.
        cycle
      endif
      where(distance2>upper)source(:,j)=(0d0,0d0)
      radii(j)=sqrt(upper)
      if(requested_radius>0d0)radii(j)=requested_radius
      loss(j)=max(0d0,1d0-kept/norm)
    enddo
    status=0
  end subroutine adaptive_source_mask


  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module exx_adaptive_support
