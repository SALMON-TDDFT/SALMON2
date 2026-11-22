! Matrix-free static polarization eigensubspace, finite-q diagnostic support.
module hsex_pdep
  implicit none
  private
  public :: pdep_basis,pdep_inverse,pdep_modes
  abstract interface
    subroutine polarization(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine mode_consensus(code)
      implicit none
      integer,intent(inout) :: code
    end subroutine
  end interface
contains
  subroutine pdep_basis(apply,n,rank,basis,values,residual,calls,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(polarization) :: apply
    integer,intent(in) :: n,rank
    complex(8),intent(out) :: basis(:,:)
    real(8),intent(out) :: values(:),residual(:)
    integer,intent(out) :: calls,status
    complex(8),allocatable :: q(:,:),bq(:,:),t(:,:),work(:),v(:),r(:,:)
    real(8),allocatable :: ev(:),rwork(:)
    real(8) :: length
    integer :: m,j,i,pass,seed,info
    status=1;calls=0;basis=0d0;values=0d0;residual=0d0
    if(n<1.or.rank<1.or.rank>n)return
    if(any(shape(basis)/=[n,rank]).or.size(values)/=rank.or.size(residual)/=rank)return
    m=min(n,2*rank+8)
    allocate(q(n,m),bq(n,m),t(m,m),work(max(1,2*m)),v(n),r(n,rank),ev(m),rwork(max(1,3*m-2)))
    do i=1,n
      v(i)=cmplx(sin(dble(i)),cos(2d0*i),8)
    enddo
    seed=0
    do j=1,m
      do
        do pass=1,2
          do i=1,j-1
            v=v-q(:,i)*dot_product(q(:,i),v)
          enddo
        enddo
        length=norm2(abs(v))
        if(.not.screen_finite([length]))return
        if(length>1d-12)exit
        seed=seed+1
        if(seed>n)return
        v=0d0;v(seed)=1d0
      enddo
      q(:,j)=v/length
      call apply(q(:,j),bq(:,j),info);calls=calls+1
      status=2
      if(info/=0.or..not.screen_finite(bq(:,j)))return
      v=bq(:,j)
    enddo
    t=matmul(conjg(transpose(q)),bq)
    if(maxval(abs(t-conjg(transpose(t))))>1d-7*max(1d0,maxval(abs(t))))return
    t=.5d0*(t+conjg(transpose(t)))
    call zheev('V','L',m,t,m,ev,work,size(work),rwork,info)
    status=3
    if(info/=0.or..not.screen_finite(ev))return
    if(minval(ev)<-1d-8*max(1d0,maxval(abs(ev))))return
    basis=matmul(q,t(:,m-rank+1:m));values=max(0d0,ev(m-rank+1:m))
    r=matmul(bq,t(:,m-rank+1:m))
    do j=1,rank
      residual(j)=norm2(abs(r(:,j)-basis(:,j)*values(j)))/max(1d0,values(j))
    enddo
    if(.not.screen_finite(basis).or..not.screen_finite(residual))return
    status=0
  end subroutine

  subroutine pdep_modes(apply,n,rank,tolerance,maxiter,basis,values,residual,calls,status,initial,synchronize)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(polarization) :: apply
    integer,intent(in) :: n,rank,maxiter
    real(8),intent(in) :: tolerance
    complex(8),intent(out) :: basis(:,:)
    real(8),intent(out) :: values(:),residual(:)
    integer,intent(out) :: calls,status
    complex(8),optional,intent(in) :: initial(:,:)
    procedure(mode_consensus),optional :: synchronize
    complex(8),allocatable :: q(:,:),aq(:,:),t(:,:),work(:),candidate(:),modes(:,:),actions(:,:),errors(:,:)
    real(8),allocatable :: eig(:),rwork(:)
    integer :: capacity,m,applied,j,i,pass,iteration,info,seed,old_m,seek
    real(8),allocatable :: guard_residual(:)
    real(8) :: length
    status=1;calls=0;basis=0d0;values=0d0;residual=0d0
    if(n<1.or.rank<1.or.rank>n.or.maxiter<1)return
    if(.not.screen_finite([tolerance]).or.tolerance<=0d0)return
    if(any(shape(basis)/=[n,rank]).or.size(values)/=rank.or.size(residual)/=rank)return
    seek=min(n,rank+4)
    capacity=min(n,2*seek+8)
    allocate(q(n,capacity),aq(n,capacity),t(capacity,capacity),work(2*capacity),candidate(n), &
      modes(n,seek),actions(n,seek),errors(n,seek),eig(capacity),rwork(3*capacity),guard_residual(seek),stat=info)
    if(present(synchronize))call synchronize(info)
    if(info/=0)return
    q=0d0;aq=0d0;m=0;applied=0;seed=0
    if(present(initial))then
      if(any(shape(initial)/=[n,rank]).or..not.screen_finite(initial))return
      do j=1,rank
        candidate=initial(:,j);call append()
        if(m/=j)return
      enddo
    endif
    ! Fresh full-support probes are required even with an exact cached space.
    ! Converging guard Ritz modes prevents an invariant but obsolete cache
    ! from masquerading as the dominant current eigensubspace.
    do while(m<seek)
      seed=seed+1
      if(seed>n)return
      do i=1,n
        candidate(i)=cmplx(sin(dble(i*seed)*sqrt(2d0))+cos(dble(i*(seed+1))*sqrt(3d0)),0d0,8)
      enddo
      call append()
    enddo
    do iteration=1,maxiter
      do j=applied+1,m
        call apply(q(:,j),aq(:,j),info);calls=calls+1
        status=2
        if(info/=0.or..not.screen_finite(aq(:,j)))return
      enddo
      applied=m
      t(:m,:m)=matmul(conjg(transpose(q(:,:m))),aq(:,:m))
      if(maxval(abs(t(:m,:m)-conjg(transpose(t(:m,:m)))))> &
        1d-7*max(1d0,maxval(abs(t(:m,:m)))))return
      t(:m,:m)=.5d0*(t(:m,:m)+conjg(transpose(t(:m,:m))))
      call zheev('V','L',m,t,capacity,eig,work,size(work),rwork,info)
      status=3
      if(.not.screen_finite(eig(:m)))info=1
      if(present(synchronize))call synchronize(info)
      if(info/=0)return
      if(minval(eig(:m))< -1d-8*max(1d0,maxval(abs(eig(:m)))))return
      modes=matmul(q(:,:m),t(:m,m-seek+1:m))
      actions=matmul(aq(:,:m),t(:m,m-seek+1:m))
      values=max(0d0,eig(m-rank+1:m));errors=actions
      do j=1,seek
        errors(:,j)=errors(:,j)-eig(m-seek+j)*modes(:,j)
        guard_residual(j)=norm2(abs(errors(:,j)))/max(1d0,abs(eig(m-seek+j)))
      enddo
      residual=guard_residual(seek-rank+1:seek)
      if(.not.screen_finite(guard_residual).or..not.screen_finite(modes))return
      basis=modes(:,seek-rank+1:seek)
      if(maxval(guard_residual)<=tolerance)then
        status=0
        return
      endif
      if(capacity-m<seek)then
        q(:,:seek)=modes;aq(:,:seek)=actions;m=seek;applied=seek
      endif
      old_m=m
      do j=1,seek
        if(m==capacity)exit
        candidate=errors(:,j);call append()
      enddo
      if(m==old_m)then
        do while(m==old_m.and.seed<n.and.m<capacity)
          seed=seed+1;candidate=0d0;candidate(seed)=1d0;call append()
        enddo
      endif
      status=4
      if(m==old_m)return
    enddo
    status=5
  contains
    subroutine append()
      implicit none
      if(m>=capacity)return
      do pass=1,2
        do i=1,m
          candidate=candidate-q(:,i)*dot_product(q(:,i),candidate)
        enddo
      enddo
      length=norm2(abs(candidate))
      if(.not.screen_finite([length]))return
      if(length<=1d-12)return
      m=m+1;q(:,m)=candidate/length
    end subroutine
  end subroutine pdep_modes

  subroutine pdep_inverse(basis,values,x,y,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: basis(:,:),x(:)
    real(8),intent(in) :: values(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    status=1;y=0d0
    if(size(basis,1)/=size(x).or.size(y)/=size(x).or.size(basis,2)/=size(values))return
    if(.not.screen_finite(basis).or..not.screen_finite(values).or..not.screen_finite(x))return
    if(any(values<0d0))return
    y=x-matmul(basis,(values/(1d0+values))*matmul(conjg(transpose(basis)),x))
    if(.not.screen_finite(y))return
    status=0
  end subroutine
end module
