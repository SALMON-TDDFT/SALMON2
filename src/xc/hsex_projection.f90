! Nonthermal projection diagnostics and real-space auxiliary-orbital guesses.
module hsex_projection
  implicit none
  private
  public :: s_hsex_capture,projection_capture,projection_seed
  type :: s_hsex_capture
    real(8) :: missing_norm=0d0,lost_electrons=0d0
    real(8) :: propagated_electrons=0d0,projected_electrons=0d0
    real(8) :: excited_electrons=0d0,holes=0d0
  end type s_hsex_capture
contains
  subroutine projection_capture(overlap,norms,occupations,weights,nv,tolerance,result,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: overlap(:,:,:)
    real(8),intent(in) :: norms(:,:),occupations(:,:),weights(:),tolerance
    integer,intent(in) :: nv
    type(s_hsex_capture),intent(out) :: result
    integer,intent(out) :: status
    real(8) :: represented,missing,population
    integer :: k,i,a
    result=s_hsex_capture();status=1
    if(any(shape(norms)/=[size(overlap,2),size(overlap,3)]))return
    if(any(shape(occupations)/=shape(norms)).or.size(weights)/=size(norms,2))return
    if(nv<1.or.nv>size(overlap,1))return
    if(.not.screen_finite(overlap).or..not.screen_finite(norms))return
    if(.not.screen_finite(occupations).or..not.screen_finite(weights).or..not.screen_finite([tolerance]))return
    if(tolerance<=0d0.or.any(norms<=0d0).or.any(occupations<0d0).or.any(occupations>2d0))return
    if(any(weights<0d0).or.abs(sum(weights)-1d0)>1d-10)return
    status=2
    do k=1,size(norms,2)
      do i=1,size(norms,1)
        represented=sum(abs(overlap(:,i,k))**2)
        if(.not.screen_finite([represented]))return
        missing=norms(i,k)-represented
        if(missing < -tolerance)return ! Nonorthogonal auxiliaries must not be repaired by rescaling.
        if(occupations(i,k)>0d0)result%missing_norm=max(result%missing_norm,max(0d0,missing))
        result%lost_electrons=result%lost_electrons+weights(k)*occupations(i,k)*max(0d0,missing)
        result%propagated_electrons=result%propagated_electrons+weights(k)*occupations(i,k)*norms(i,k)
        result%projected_electrons=result%projected_electrons+weights(k)*occupations(i,k)*represented
      enddo
      do a=1,size(overlap,1)
        population=sum(occupations(:,k)*abs(overlap(a,:,k))**2)
        if(a>nv)then
          result%excited_electrons=result%excited_electrons+weights(k)*population
        else
          result%holes=result%holes+weights(k)*(2d0-population)
        endif
      enddo
    enddo
    if(.not.screen_finite([result%missing_norm,result%lost_electrons,result%propagated_electrons, &
      result%projected_electrons,result%excited_electrons,result%holes]))return
    status=0
  end subroutine projection_capture

  subroutine projection_seed(old,propagated,expanded,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: old(:,:),propagated(:,:)
    complex(8),intent(out) :: expanded(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: candidate(:),gram(:,:)
    real(8) :: length,angle
    integer :: n,have,source,i,j,pass
    status=1;expanded=0d0;n=size(old,1);have=size(old,2)
    if(n<1.or.size(propagated,1)/=n.or.size(expanded,1)/=n)return
    if(size(expanded,2)<have.or.size(expanded,2)>n)return
    if(.not.screen_finite(old).or..not.screen_finite(propagated))return
    gram=matmul(conjg(transpose(old)),old)
    do i=1,have
      gram(i,i)=gram(i,i)-1d0
    enddo
    if(.not.screen_finite(gram))return
    if(maxval(abs(gram))>1d-8)return
    expanded(:,:have)=old
    allocate(candidate(n))
    ! Missing RT components are efficient initial guesses. Deterministic grid
    ! vectors complete the guesses if their rank is insufficient. Neither
    ! source orbitals nor their occupations are modified or renormalized.
    do source=1,size(propagated,2)+n
      if(have==size(expanded,2))exit
      if(source<=size(propagated,2))then
        candidate=propagated(:,source)
      else
        do i=1,n
          angle=2d0*acos(-1d0)*dble(i-1)*(source-size(propagated,2)-1)/n
          candidate(i)=cmplx(cos(angle),sin(angle),8)/sqrt(dble(n))
        enddo
      endif
      do pass=1,2
        do j=1,have
          candidate=candidate-expanded(:,j)*dot_product(expanded(:,j),candidate)
        enddo
      enddo
      length=norm2(abs(candidate))
      if(.not.screen_finite([length]))return
      if(length<=1d-10)cycle
      have=have+1;expanded(:,have)=candidate/length
    enddo
    if(have/=size(expanded,2).or..not.screen_finite(expanded))return
    status=0
  end subroutine projection_seed
end module hsex_projection
