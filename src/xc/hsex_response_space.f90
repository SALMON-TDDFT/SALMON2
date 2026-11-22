! Shared projected response space for exact energy shifts; opt-in diagnostics.
module hsex_response_space
  implicit none
  private
  public :: response_space_solve,s_response_space_stats
  type :: s_response_space_stats
    integer :: dimension=0,h_actions=0,validation_actions=0,iterations=0
    real(8) :: build_seconds=0d0,validation_seconds=0d0
  end type
  abstract interface
    subroutine response_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine response_space_solve(apply,target,energy,vertices,maxrank,tolerance,x,residual,statistics,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(response_action) :: apply
    complex(8),intent(in) :: target(:,:),vertices(:,:)
    real(8),intent(in) :: energy(:),tolerance
    integer,intent(in) :: maxrank
    complex(8),intent(out) :: x(:,:)
    real(8),intent(out) :: residual(:)
    type(s_response_space_stats),intent(out) :: statistics
    integer,intent(out) :: status
    complex(8),allocatable :: q(:,:),aq(:,:),rhs(:,:),r(:,:),t(:,:),coef(:,:),small(:,:),work(:)
    complex(8),allocatable :: candidate(:),action(:),gram(:,:)
    real(8),allocatable :: eigen(:),rwork(:),scale(:),lengths(:)
    real(8) :: length,denominator
    integer :: n,nrhs,capacity,m,applied,i,j,pass,info,pivot,previous,additions
    integer(8) :: started,finished,rate,validated
    statistics=s_response_space_stats();status=1;x=0d0;residual=huge(1d0)
    n=size(vertices,1);nrhs=size(vertices,2)
    if(n<1.or.nrhs<1.or.size(target,1)/=n.or.size(target,2)>=n)return
    if(size(energy)/=nrhs.or.any(shape(x)/=shape(vertices)).or.size(residual)/=nrhs)return
    if(maxrank<1.or..not.screen_finite([tolerance]).or.tolerance<=0d0.or.tolerance>=1d0)return
    if(.not.screen_finite(target).or..not.screen_finite(vertices).or..not.screen_finite(energy))return
    capacity=min(maxrank,n-size(target,2))
    allocate(q(n,capacity),aq(n,capacity),rhs(n,nrhs),r(n,nrhs),t(capacity,capacity), &
      coef(capacity,nrhs),small(capacity,nrhs),work(2*capacity),candidate(n),action(n), &
      gram(size(target,2),size(target,2)),eigen(capacity),rwork(3*capacity),scale(nrhs),lengths(nrhs),stat=info)
    if(info/=0)return
    gram=matmul(conjg(transpose(target)),target)
    do i=1,size(target,2)
      gram(i,i)=gram(i,i)-1d0
    enddo
    if(.not.screen_finite(gram))return
    if(size(gram)>0)then
      if(maxval(abs(gram))>1d-8)return
    endif
    call system_clock(started,rate)
    rhs=vertices-matmul(target,matmul(conjg(transpose(target)),vertices))
    do j=1,nrhs
      scale(j)=max(norm2(abs(rhs(:,j))),.01d0*norm2(abs(vertices(:,j))),tiny(1d0))
      r(:,j)=rhs(:,j)/scale(j)
    enddo
    m=0;applied=0;q=0d0;aq=0d0
    do
      ! Pivoted, twice-orthogonalized residual directions compress dependent
      ! right-hand sides. Reuse stored QHQ actions across every energy shift.
      previous=m;additions=min(16,capacity-m)
      do i=1,additions
        do j=1,nrhs
          lengths(j)=norm2(abs(r(:,j)))
        enddo
        pivot=maxloc(lengths,dim=1)
        if(lengths(pivot)<=max(1d-13,.1d0*tolerance))exit
        candidate=r(:,pivot)
        do pass=1,2
          candidate=candidate-matmul(target,matmul(conjg(transpose(target)),candidate))
          if(m>0)candidate=candidate-matmul(q(:,:m),matmul(conjg(transpose(q(:,:m))),candidate))
        enddo
        length=norm2(abs(candidate))
        if(length<=1d-13)exit
        m=m+1;q(:,m)=candidate/length
        do pass=1,2
          do j=1,nrhs
            r(:,j)=r(:,j)-q(:,m)*dot_product(q(:,m),r(:,j))
          enddo
        enddo
      enddo
      if(m==0)then
        status=0
        exit
      endif
      if(m==previous)then
        status=2
        exit
      endif
      do i=applied+1,m
        call apply(q(:,i),action,info);statistics%h_actions=statistics%h_actions+1
        status=3
        if(info/=0.or..not.screen_finite(action))return
        aq(:,i)=action-matmul(target,matmul(conjg(transpose(target)),action))
      enddo
      applied=m;statistics%dimension=m;statistics%iterations=statistics%iterations+1
      call zgemm('C','N',m,m,n,(1d0,0d0),q,n,aq,n,(0d0,0d0),t,capacity)
      status=4
      if(.not.screen_finite(t(:m,:m)))return
      if(maxval(abs(t(:m,:m)-conjg(transpose(t(:m,:m)))))>1d-8*max(1d0,maxval(abs(t(:m,:m)))))return
      t(:m,:m)=.5d0*(t(:m,:m)+conjg(transpose(t(:m,:m))))
      call zheev('V','L',m,t,capacity,eigen,work,size(work),rwork,info)
      if(info/=0.or..not.screen_finite(eigen(:m)))return
      if(eigen(1)<=maxval(energy))return
      call zgemm('C','N',m,nrhs,n,(1d0,0d0),q,n,rhs,n,(0d0,0d0),coef,capacity)
      call zgemm('C','N',m,nrhs,m,(1d0,0d0),t,capacity,coef,capacity,(0d0,0d0),small,capacity)
      do j=1,nrhs
        do i=1,m
          denominator=eigen(i)-energy(j)
          small(i,j)=small(i,j)/denominator
        enddo
      enddo
      call zgemm('N','N',m,nrhs,m,(1d0,0d0),t,capacity,small,capacity,(0d0,0d0),coef,capacity)
      call zgemm('N','N',n,nrhs,m,(1d0,0d0),q,n,coef,capacity,(0d0,0d0),x,n)
      r=rhs
      call zgemm('N','N',n,nrhs,m,(-1d0,0d0),aq,n,coef,capacity,(1d0,0d0),r,n)
      do j=1,nrhs
        r(:,j)=(r(:,j)+energy(j)*x(:,j))/scale(j)
        residual(j)=norm2(abs(r(:,j)))
      enddo
      status=5
      if(.not.screen_finite(x).or..not.screen_finite(residual))return
      if(maxval(residual)<=tolerance)then
        status=0
        exit
      endif
      if(m==capacity)then
        status=2
        exit
      endif
    enddo
    call system_clock(finished)
    ! Independently apply the original Hamiltonian, rather than accepting only
    ! the Galerkin residual assembled from cached basis actions.
    do j=1,nrhs
      call apply(x(:,j),action,info)
      statistics%h_actions=statistics%h_actions+1
      statistics%validation_actions=statistics%validation_actions+1
      if(info/=0.or..not.screen_finite(action))then
        status=3
        return
      endif
      action=action-matmul(target,matmul(conjg(transpose(target)),action))
      residual(j)=norm2(abs(rhs(:,j)-action+energy(j)*x(:,j)))/scale(j)
    enddo
    if(.not.screen_finite(residual))then
      status=5
    else if(maxval(residual)>tolerance)then
      status=2
    endif
    call system_clock(validated)
    statistics%build_seconds=dble(finished-started)/rate
    statistics%validation_seconds=dble(validated-finished)/rate
  end subroutine
end module
