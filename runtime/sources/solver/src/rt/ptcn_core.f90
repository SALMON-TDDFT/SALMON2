! Functional-independent PT-CN: callbacks own frozen-inner/fresh-outer operator semantics.
module ptcn_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: exx_ptcn_solve,pt_orthogonalize,pt_predict_endpoint,pt_gradient,pt_uniform_occupation
  abstract interface
    subroutine mixing_reduce_callback(matrix)
      complex(8),intent(inout)::matrix(:,:)
    end subroutine
    subroutine consensus_callback(status)
      integer,intent(inout)::status
    end subroutine
    subroutine action_callback(a,b,refresh)
      implicit none
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
      logical,intent(in) :: refresh
    end subroutine
    subroutine projection_callback(a,b,dv,r,nocc)
      complex(8),intent(in)::a(:,:,:),b(:,:,:)
      complex(8),intent(out)::r(:,:,:)
      real(8),intent(in)::dv
      integer,intent(in)::nocc
    end subroutine
    subroutine precondition_callback(a,b)
      implicit none
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
    end subroutine
    real(8) function norm_callback(a)
      implicit none
      complex(8),intent(in) :: a(:,:,:)
    end function
  end interface
contains
  logical function pt_uniform_occupation(occ,expected) result(valid)
    real(8),intent(in)::occ(:,:,:),expected
    integer::i,j,k
    valid=.false.
    if(any(shape(occ)<1))return
    if(.not.ieee_is_finite(expected))return
    do k=1,size(occ,3);do j=1,size(occ,2);do i=1,size(occ,1)
      if(.not.ieee_is_finite(occ(i,j,k)))return
      if(abs(occ(i,j,k)-expected)>1d-12)return
    enddo;enddo;enddo
    valid=.true.
  end function
  subroutine pt_predict_endpoint(current,previous,dt,previous_dt,guess,used)
    ! Caller must align equal-occupation gauges and orbital ownership first.
    complex(8),intent(in) :: current(:,:,:),previous(:,:,:)
    real(8),intent(in) :: dt,previous_dt
    complex(8),intent(out) :: guess(:,:,:)
    logical,intent(out) :: used
    real(8) :: ratio,base,change
    used=.false.;guess=0
    if(any(shape(guess)/=shape(current)))return
    guess=current
    if(any(shape(previous)/=shape(current)))return
    if(.not.finite_array(current).or..not.finite_array(previous))return
    if(.not.ieee_is_finite(dt).or..not.ieee_is_finite(previous_dt))return
    if(dt<=0.or.previous_dt<=0)return
    ratio=dt/previous_dt
    if(.not.ieee_is_finite(ratio).or.ratio<.5d0.or.ratio>2d0)return
    base=maxval(abs(current));change=maxval(abs(current-previous))
    if(base<=tiny(1d0).or..not.ieee_is_finite(change))return
    if(change>.25d0*base/ratio)return
    guess=current+ratio*(current-previous)
    if(.not.finite_array(guess))then
      guess=current;return
    endif
    used=.true.
  end subroutine
  logical function finite_array(a) result(ok)
    complex(8),intent(in) :: a(:,:,:)
    integer :: i,j,k
    ok=.false.
    do k=1,size(a,3);do j=1,size(a,2);do i=1,size(a,1)
      if(.not.ieee_is_finite(real(a(i,j,k),8)))return
      if(.not.ieee_is_finite(aimag(a(i,j,k))))return
    enddo;enddo;enddo
    ok=.true.
  end function
  subroutine exx_ptcn_solve(u,dt,dv,action,precondition,norm,x,error,builds,apps,ierr,noccupied,use_anderson,initial_guess,initial_hu,endpoint_hu,projection,consensus,mixing_reduce)
    implicit none
    complex(8),intent(in) :: u(:,:,:)
    real(8),intent(in) :: dt,dv
    procedure(action_callback) :: action
    procedure(precondition_callback) :: precondition
    procedure(norm_callback) :: norm
    procedure(projection_callback),optional :: projection
    procedure(consensus_callback),optional :: consensus
    procedure(mixing_reduce_callback),optional :: mixing_reduce
    complex(8),intent(out) :: x(:,:,:)
    real(8),intent(out) :: error
    integer,intent(out) :: builds,apps,ierr
    integer,intent(in),optional :: noccupied
    logical,intent(in),optional :: use_anderson
    ! initial_hu must be the certified fresh action at exactly u, including fields.
    ! Caller must also install matching ACE and callback time/field state.
    ! Skipping the callback must not leave a first_action/ac_start transition pending.
    ! Caller owns cache identity/invalidation; a wrong RHS can pass final residual.
    ! Optional cache arrays must not alias each other or solver input/output arrays.
    complex(8),intent(in),optional :: initial_guess(:,:,:),initial_hu(:,:,:)
    complex(8),intent(out),optional :: endpoint_hu(:,:,:)
    logical :: mixing
    integer :: history_count,mixing_status
    complex(8),allocatable :: images(:,:,:,:),residuals(:,:,:,:)
    complex(8),allocatable :: hu(:,:,:),rhs(:,:,:),r(:,:,:),correction(:,:,:)
    real(8) :: scale,errors(150)
    integer :: outer,inner,nocc,admission,allocation_status
    ierr=1;builds=0;apps=0;error=huge(1d0)
    admission=0
    if(present(endpoint_hu))then
      endpoint_hu=0
      if(any(shape(endpoint_hu)/=shape(u)))admission=1
    endif
    if(dt<=0.or.dv<=0.or..not.ieee_is_finite(dt).or..not.ieee_is_finite(dv))admission=1
    if(any(shape(u)/=shape(x)))admission=1
    if(.not.finite_array(u))admission=1
    if(present(initial_guess))then
      if(any(shape(initial_guess)/=shape(u)))admission=1
      if(.not.finite_array(initial_guess))admission=1
    endif
    if(present(initial_hu))then
      if(any(shape(initial_hu)/=shape(u)))admission=1
      if(.not.finite_array(initial_hu))admission=1
    endif
    nocc=size(u,2)
    if(present(noccupied))nocc=noccupied
    if(nocc<1)admission=1
    if(.not.present(projection))then
      if(nocc>size(u,2))admission=1
    endif
    if(present(consensus))call consensus(admission)
    if(admission/=0)return
    scale=norm(u)
    admission=0
    if(scale<=0.or..not.ieee_is_finite(scale))admission=1
    if(present(consensus))call consensus(admission)
    if(admission/=0)return
    allocate(hu(size(u,1),size(u,2),size(u,3)),rhs(size(u,1),size(u,2),size(u,3)), &
      r(size(u,1),size(u,2),size(u,3)),correction(size(u,1),size(u,2),size(u,3)),stat=allocation_status)
    admission=abs(allocation_status)
    if(present(consensus))call consensus(admission)
    if(admission/=0)return
    mixing=.false.
    if(present(use_anderson))mixing=use_anderson
    if(mixing)then
      allocate(images(size(u,1),size(u,2),size(u,3),10),residuals(size(u,1),size(u,2),size(u,3),10),stat=allocation_status)
      admission=abs(allocation_status)
      if(present(consensus))call consensus(admission)
      if(admission/=0)return
    endif
    x=u
    if(present(initial_guess))x=initial_guess
    if(present(initial_hu))then
      hu=initial_hu
    else
      call action(u,hu,.true.);builds=1
    endif
    call gradient(u,hu,dv,r,nocc)
    rhs=u-(0d0,.5d0)*dt*r
    do outer=1,12
      history_count=0
      do inner=1,150
        call action(x,hu,.false.);apps=apps+1
        call gradient(x,hu,dv,r,nocc)
        r=x+(0d0,.5d0)*dt*r-rhs
        error=norm(r)/scale;errors(inner)=error
        if(.not.ieee_is_finite(error))return
        ! Inner solves need not be tighter than the independently rebuilt
        ! full residual gate. Reserve a factor-two margin for that check.
        ! Preserve the legacy fixed-point accuracy; Anderson uses the tested margin.
        if(error<1d-12.or.(mixing.and.error<5d-11))exit
        call precondition(r,correction)
        if(mixing)then
          call anderson_update(x,correction,images,residuals,history_count,mixing_status,mixing_reduce)
          if(present(consensus))call consensus(mixing_status)
          if(mixing_status/=0)return
        else
          x=x-correction
        endif
      enddo
      ! A roundoff plateau below the FINAL tolerance still requires a fresh
      ! operator rebuild and residual check; never accept on the inner estimate.
      if(inner>150.and.error>=1d-10)return
      call action(x,hu,.true.);builds=builds+1
      call gradient(x,hu,dv,r,nocc)
      r=x+(0d0,.5d0)*dt*r-rhs
      error=norm(r)/scale
      if(.not.ieee_is_finite(error))return
      if(error<1d-10)then
        if(present(endpoint_hu))endpoint_hu=hu
        ierr=0;return
      endif
    enddo
  contains
    subroutine gradient(a,b,volume,out,count)
      complex(8),intent(in)::a(:,:,:),b(:,:,:)
      complex(8),intent(out)::out(:,:,:)
      real(8),intent(in)::volume
      integer,intent(in)::count
      if(present(projection))then
        call projection(a,b,volume,out,count)
      else
        call pt_gradient(a,b,volume,out,count)
      endif
    end subroutine
  end subroutine

  subroutine anderson_update(x,correction,images,residuals,count,ierr,reduce)
    implicit none
    procedure(mixing_reduce_callback),optional::reduce
    complex(8),intent(inout) :: x(:,:,:),images(:,:,:,:),residuals(:,:,:,:)
    complex(8),intent(in) :: correction(:,:,:)
    integer,intent(inout) :: count
    integer,intent(out) :: ierr
    complex(8) :: packed(min(count+1,size(images,4))**2,size(x,2)*size(x,3))
    complex(8) :: matrix(min(count+1,size(images,4)),min(count+1,size(images,4)))
    complex(8) :: coefficient(min(count+1,size(images,4)),1),total
    real(8) :: scale
    integer :: nh,slot,i,j,io,ik,col
    external :: zposv
    ierr=0;count=count+1;nh=min(count,size(images,4));slot=mod(count-1,size(images,4))+1
    images(:,:,:,slot)=x-.2d0*correction
    if(present(reduce))images(:,:,:,slot)=x-.5d0*correction
    residuals(:,:,:,slot)=correction
    if(present(reduce))then
      ! One common coefficient vector for the coupled global nonlinear state.
      do j=1,nh;do i=1,nh
        matrix(i,j)=sum(conjg(residuals(:,:,:,i))*residuals(:,:,:,j))
      enddo;enddo
      call reduce(matrix)
      do j=1,nh;do i=1,nh
        if(.not.ieee_is_finite(real(matrix(i,j),8)).or..not.ieee_is_finite(aimag(matrix(i,j))))then
          ierr=1;return
        endif
      enddo;enddo
      scale=maxval(abs(matrix))
      if(scale<=tiny(1d0))then
        x=images(:,:,:,slot);return
      endif
      matrix=cmplx(real(matrix,8)/scale,0d0,8)
      do i=1,nh;matrix(i,i)=matrix(i,i)+1d-12;enddo
      coefficient=(1d0,0d0)
      call zposv('U',nh,1,matrix,nh,coefficient,nh,ierr)
      if(ierr/=0)return
      total=sum(coefficient)
      if(.not.ieee_is_finite(real(total,8)).or..not.ieee_is_finite(aimag(total)))then
        ierr=1;return
      endif
      if(abs(total)<=tiny(1d0))then
        ierr=1;return
      endif
      coefficient=coefficient/total
      x=0
      do i=1,nh
        if(.not.ieee_is_finite(real(coefficient(i,1),8)).or..not.ieee_is_finite(aimag(coefficient(i,1))))then
          ierr=1;return
        endif
        x=x+coefficient(i,1)*images(:,:,:,i)
      enddo
      return
    endif
    do ik=1,size(x,3);do io=1,size(x,2)
      col=(ik-1)*size(x,2)+io
      do j=1,nh;do i=1,nh
        packed((j-1)*nh+i,col)=sum(conjg(residuals(:,io,ik,i))*residuals(:,io,ik,j))
      enddo;enddo
    enddo;enddo
    do ik=1,size(x,3);do io=1,size(x,2)
      col=(ik-1)*size(x,2)+io;matrix=reshape(packed(:,col),[nh,nh])
      scale=maxval(abs(matrix))
      if(.not.ieee_is_finite(scale).or.scale<=tiny(1d0))then
        ierr=1;return
      endif
      matrix=cmplx(real(matrix,8)/scale,0d0,8)
      do i=1,nh;matrix(i,i)=matrix(i,i)+1d-12;enddo
      coefficient=(1d0,0d0)
      call zposv('U',nh,1,matrix,nh,coefficient,nh,ierr)
      if(ierr/=0)return
      total=sum(coefficient)
      if(abs(total)<=tiny(1d0))then
        ierr=1;return
      endif
      coefficient=coefficient/total
      do i=1,nh
        if(.not.ieee_is_finite(real(coefficient(i,1),8)).or. &
           .not.ieee_is_finite(aimag(coefficient(i,1))))then
          ierr=1;return
        endif
      enddo
      x(:,io,ik)=(0d0,0d0)
      do i=1,nh;x(:,io,ik)=x(:,io,ik)+coefficient(i,1)*images(:,io,ik,i);enddo
    enddo;enddo
  end subroutine

  subroutine pt_gradient(u,hu,dv,r,nocc)
    implicit none
    complex(8),intent(in) :: u(:,:,:),hu(:,:,:)
    complex(8),intent(out) :: r(:,:,:)
    real(8),intent(in) :: dv
    complex(8) :: projection(size(u,2),size(u,2))
    integer,intent(in) :: nocc
    integer :: ik,ng,no,block,first,last,nb
    complex(8),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    external :: zgemm
    ng=size(u,1);no=size(u,2);r=hu
    do ik=1,size(u,3)
      ! Gauge rotations must commute with the occupation matrix.  Mixing
      ! occupied and empty targets would change the physical density matrix.
      do block=1,2
        first=1;last=nocc
        if(block==2)then
          first=nocc+1;last=no
        endif
        nb=last-first+1
        if(nb==0)cycle
        call zgemm('C','N',nb,nb,ng,one*dv,u(1,first,ik),ng,hu(1,first,ik),ng,zero,projection(1,1),no)
        call zgemm('N','N',ng,nb,nb,-one,u(1,first,ik),ng,projection(1,1),no,one,r(1,first,ik),ng)
      enddo
    enddo
  end subroutine
  subroutine pt_orthogonalize(x,dv,ierr)
    implicit none
    complex(8),intent(inout) :: x(:,:,:)
    real(8),intent(in) :: dv
    integer,intent(out) :: ierr
    complex(8) :: metric(size(x,2),size(x,2)),transform(size(x,2),size(x,2))
    complex(8) :: vectors(size(x,2),size(x,2)),work(max(1,2*size(x,2)))
    real(8) :: eigen(size(x,2)),rwork(max(1,3*size(x,2)-2))
    integer :: ik,j,no
    external :: zheev
    ierr=1
    if(any(shape(x)<1))return
    if(.not.finite_array(x))return
    if(.not.ieee_is_finite(dv))return
    if(dv<=0d0)return
    no=size(x,2);ierr=0
    do ik=1,size(x,3)
      metric=dv*matmul(conjg(transpose(x(:,:,ik))),x(:,:,ik))
      call zheev('V','U',no,metric,no,eigen,work,size(work),rwork,ierr)
      if(ierr/=0)return
      do j=1,no
        if(.not.ieee_is_finite(eigen(j)))then
          ierr=1;return
        endif
      enddo
      if(minval(eigen)<=epsilon(1d0)*maxval(eigen))then
        ierr=1;return
      endif
      vectors=metric
      do j=1,no;vectors(:,j)=vectors(:,j)/sqrt(eigen(j));enddo
      transform=matmul(vectors,conjg(transpose(metric)))
      x(:,:,ik)=matmul(x(:,:,ik),transform)
    enddo
  end subroutine
end module
