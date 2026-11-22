! Frequency-domain resolvent of a real-linear occupied-orbital tangent equation.
module hsex_impulse_resolvent
  implicit none
  private
  public :: impulse_model,impulse_build,impulse_resolve,impulse_trajectory
  type :: impulse_model
    real(8),allocatable :: h(:,:),g(:)
    real(8),allocatable :: residual_factor(:,:)
    real(8) :: source_norm=0d0,next_norm=0d0
    real(8) :: validation_residual=0d0
    integer :: actions=0
  end type
  abstract interface
    subroutine tangent_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine tangent_observation(x,value,status)
      implicit none
      complex(8),intent(in) :: x(:)
      real(8),intent(out) :: value
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine impulse_build(apply,observe,source,shifts,maxrank,tol,icomm,rank,model,status,precondition)
    use communication, only: comm_summation,comm_get_max,comm_bcast,comm_is_root
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(tangent_action) :: apply
    procedure(tangent_observation) :: observe
    procedure(tangent_action),optional :: precondition
    complex(8),intent(in) :: source(:),shifts(:)
    integer,intent(in) :: maxrank,icomm,rank
    real(8),intent(in) :: tol
    type(impulse_model),intent(out) :: model
    integer,intent(out) :: status
    complex(8),allocatable :: basis(:,:),action(:),transfer(:)
    real(8),allocatable :: h(:,:),g(:),local(:),global(:),residual(:)
    real(8) :: squared,norm,value
    integer :: j,pass,code,stop_build,n
    status=1;n=size(source)
    if(present(precondition))then
      call rational_build(apply,observe,precondition,source,shifts,maxrank,tol,icomm,rank,model,status)
      return
    endif
    if(maxrank<1.or.size(shifts)<1.or..not.screen_finite([tol]).or.tol<=0d0)return
    code=0
    if(.not.screen_finite(source).or..not.screen_finite(shifts))code=1
    if(any(real(shifts,8)<=0d0))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    allocate(basis(n,maxrank),action(n),h(maxrank+1,maxrank),g(maxrank),local(maxrank),global(maxrank), &
      transfer(size(shifts)),residual(size(shifts)),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    squared=sum(abs(source)**2);call comm_summation(squared,norm,icomm);norm=sqrt(norm)
    if(norm<=1d-14.or..not.screen_finite([norm]))return
    model%source_norm=norm;basis(:,1)=source/norm;h=0d0;g=0d0
    do j=1,maxrank
      call observe(basis(:,j),value,code)
      if(.not.screen_finite([value]))code=1
      call comm_get_max(code,icomm)
      status=2
      if(code/=0)return
      g(j)=value
      call apply(basis(:,j),action,code)
      if(.not.screen_finite(action))code=1
      call comm_get_max(code,icomm)
      if(code/=0)return
      model%actions=j
      do pass=1,2
        local(:j)=real(matmul(conjg(transpose(basis(:,:j))),action),8)
        call comm_summation(local(:j),global(:j),j,icomm)
        h(:j,j)=h(:j,j)+global(:j)
        action=action-matmul(basis(:,:j),cmplx(global(:j),0d0,8))
      enddo
      squared=sum(abs(action)**2);call comm_summation(squared,norm,icomm);norm=sqrt(norm)
      if(.not.screen_finite([norm]))return
      h(j+1,j)=norm;stop_build=0
      if(mod(j,256)==0.or.j==maxrank.or.norm<1d-12)then
        model%h=h(:j,:j);model%g=g(:j);model%next_norm=norm
        if(comm_is_root(rank))then
          call impulse_resolve(model,shifts,transfer,residual,code)
          if(code/=0)then
            stop_build=-1
          else
            write(*,'(a,i6,2es18.8)')'STERNHEIMER_TANGENT rank/residual/tol ',j,maxval(residual),tol
            if(maxval(residual)<=tol)stop_build=1
          endif
        endif
        call comm_bcast(stop_build,icomm,0)
        if(stop_build<0)return
        if(stop_build==1)then
          status=0;return
        endif
      endif
      if(j<maxrank)then
        if(norm<=1d-12)exit
        basis(:,j+1)=action/norm
      endif
    enddo
    status=3
  end subroutine

  subroutine impulse_resolve(model,shifts,transfer,residual,status,coefficients)
    use hsex_screen_math, only: screen_finite
    implicit none
    type(impulse_model),intent(in) :: model
    complex(8),intent(in) :: shifts(:)
    complex(8),intent(out) :: transfer(:)
    real(8),intent(out) :: residual(:)
    integer,intent(out) :: status
    complex(8),intent(out),optional :: coefficients(:,:)
    real(8),allocatable :: t(:,:),u(:,:),wr(:),wi(:),work(:)
    complex(8),allocatable :: rhs(:),x(:),base(:)
    complex(8) :: aa,bb,cc,dd,det,left,right
    logical,allocatable :: bwork(:)
    logical :: paired
    integer :: n,ns,s,sdim,info,lwork,i
    status=1;transfer=0d0;residual=huge(1d0)
    if(.not.allocated(model%h).or..not.allocated(model%g))return
    n=size(model%g);ns=size(shifts)
    if(n<1.or.any(shape(model%h)/=[n,n]).or.size(transfer)/=ns.or.size(residual)/=ns)return
    if(present(coefficients))then
      if(any(shape(coefficients)/=[n,ns]))return
    endif
    if(allocated(model%residual_factor))then
      if(any(shape(model%residual_factor)/=[n,n]).or..not.screen_finite(model%residual_factor))return
    endif
    if(.not.screen_finite(model%h).or..not.screen_finite(model%g).or..not.screen_finite(shifts))return
    if(.not.screen_finite([model%source_norm,model%next_norm]).or.any(real(shifts,8)<=0d0))return
    allocate(t(n,n),u(n,n),wr(n),wi(n),work(1),rhs(n),x(n),base(n),bwork(n))
    t=model%h
    call dgees('V','N',select_none,n,t,n,sdim,wr,wi,u,n,work,-1,bwork,info)
    status=2
    if(info/=0)return
    lwork=max(1,int(work(1)));deallocate(work);allocate(work(lwork))
    t=model%h
    call dgees('V','N',select_none,n,t,n,sdim,wr,wi,u,n,work,lwork,bwork,info)
    if(info/=0.or..not.screen_finite(t).or..not.screen_finite(u))return
    base=cmplx(u(1,:),0d0,8)
    do s=1,ns
      rhs=base;i=n
      ! Real Schur has 1x1/2x2 diagonal blocks. Preserve its orthogonal basis
      ! and solve the complex shifted quasi-triangular system directly.
      do while(i>=1)
        paired=.false.
        if(i>1)paired=t(i,i-1)/=0d0
        if(paired)then
          left=base(i-1)+sum(t(i-1,i+1:n)*rhs(i+1:n))
          right=base(i)+sum(t(i,i+1:n)*rhs(i+1:n))
          aa=shifts(s)-t(i-1,i-1);bb=-t(i-1,i);cc=-t(i,i-1);dd=shifts(s)-t(i,i)
          det=aa*dd-bb*cc
          if(abs(det)<tiny(1d0))return
          rhs(i-1)=(dd*left-bb*right)/det;rhs(i)=(aa*right-cc*left)/det;i=i-2
        else
          aa=shifts(s)-t(i,i)
          if(abs(aa)<tiny(1d0))return
          rhs(i)=(base(i)+sum(t(i,i+1:n)*rhs(i+1:n)))/aa;i=i-1
        endif
      enddo
      x=matmul(u,rhs)
      if(.not.screen_finite(x))return
      if(present(coefficients))coefficients(:,s)=x
      transfer(s)=model%source_norm*sum(model%g*x)
      rhs=-matmul(model%h,x)+shifts(s)*x;rhs(1)=rhs(1)-1d0
      residual(s)=sqrt(sum(abs(rhs)**2)+(model%next_norm*abs(x(n)))**2)
      if(allocated(model%residual_factor))then
        residual(s)=sqrt(sum(abs(rhs)**2))+sqrt(sum(abs(matmul(model%residual_factor,x))**2))+ &
          64d0*epsilon(1d0)*n*(1d0+sqrt(sum(model%h**2))+abs(shifts(s)))*sqrt(sum(abs(x)**2))
      endif
      residual(s)=max(residual(s),model%validation_residual)
    enddo
    if(.not.screen_finite(transfer).or..not.screen_finite(residual))return
    status=0
  end subroutine

  subroutine rational_build(apply,observe,precondition,source,shifts,maxrank,tol,icomm,rank,model,status)
    use communication, only: comm_summation,comm_get_max,comm_bcast,comm_is_root
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(tangent_action) :: apply,precondition
    procedure(tangent_observation) :: observe
    complex(8),intent(in) :: source(:),shifts(:)
    integer,intent(in) :: maxrank,icomm,rank
    real(8),intent(in) :: tol
    type(impulse_model),intent(out) :: model
    integer,intent(out) :: status
    complex(8),allocatable :: q(:,:),aq(:,:),candidate(:),action(:),small(:,:),coef(:,:),transfer(:)
    complex(8),allocatable :: xr(:),xi(:),ar(:),ai(:),defect(:,:)
    real(8),allocatable :: g(:),dots(:),total_dots(:),residual(:),hlocal(:,:),hglobal(:,:)
    real(8),allocatable :: real_defect(:,:),local_r(:,:),stack(:,:),stack_sum(:,:)
    real(8) :: local_norm,global_norm,value,scale,local_error,global_error,true_worst
    integer :: n,ns,capacity,m,done,j,i,code,np,iteration,selected,previous,pass
    logical :: accepted
    status=1;n=size(source);ns=size(shifts);code=0
    if(maxrank<1.or.ns<1.or..not.screen_finite([tol]).or.tol<=0d0)code=1
    if(.not.screen_finite(source).or..not.screen_finite(shifts))code=1
    if(any(real(shifts,8)<=0d0))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    call comm_summation(1,np,icomm)
    capacity=min(16,maxrank)
    allocate(q(n,capacity),aq(n,capacity),g(capacity),dots(maxrank),total_dots(maxrank),candidate(n), &
      action(n),xr(n),xi(n),ar(n),ai(n),residual(ns),transfer(ns))
    local_norm=sum(abs(source)**2);call comm_summation(local_norm,global_norm,icomm)
    if(.not.screen_finite([global_norm]).or.global_norm<=1d-28)return
    model%source_norm=sqrt(global_norm);m=0;done=0;q=0d0;aq=0d0;g=0d0
    call append(source,accepted)
    call precondition(source,candidate,code)
    if(.not.screen_finite(candidate))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    call append(candidate,accepted)
    do iteration=1,maxrank
      do j=done+1,m
        call apply(q(:,j),aq(:,j),code)
        if(.not.screen_finite(aq(:,j)))code=1
        call comm_get_max(code,icomm)
        if(code/=0)return
        call observe(q(:,j),g(j),code)
        if(.not.screen_finite([g(j)]))code=1
        call comm_get_max(code,icomm)
        if(code/=0)return
      enddo
      done=m;model%actions=m
      allocate(small(m,m),hlocal(m,m),hglobal(m,m),coef(m,ns),defect(n,m),real_defect(2*n,m), &
        local_r(m,m),stack(np*m,m),stack_sum(np*m,m))
      small=0d0
      if(n>0)call zgemm('C','N',m,m,n,(1d0,0d0),q,n,aq,n,(0d0,0d0),small,m)
      hlocal=real(small,8);call comm_summation(hlocal,hglobal,m*m,icomm)
      model%h=hglobal;model%g=g(:m);model%next_norm=0d0
      defect=aq(:,:m)-matmul(q(:,:m),cmplx(hglobal,0d0,8))
      if(.not.screen_finite(defect))code=1
      call comm_get_max(code,icomm)
      if(code/=0)return
      real_defect(:n,:)=real(defect,8);real_defect(n+1:,:)=aimag(defect)
      call real_factor(real_defect,local_r,code)
      call comm_get_max(code,icomm)
      if(code/=0)return
      stack=0d0;stack(rank*m+1:(rank+1)*m,:)=local_r
      call comm_summation(stack,stack_sum,size(stack),icomm)
      model%residual_factor=local_r
      if(comm_is_root(rank))call real_factor(stack_sum,model%residual_factor,code)
      call comm_bcast(code,icomm,0)
      if(code/=0)return
      call comm_bcast(model%residual_factor,icomm,0)
      if(comm_is_root(rank))call impulse_resolve(model,shifts,transfer,residual,code,coefficients=coef)
      call comm_bcast(code,icomm,0)
      if(code/=0)return
      call comm_bcast(coef,icomm,0);call comm_bcast(residual,icomm,0)
      if(comm_is_root(rank))write(*,'(a,i6,2es18.8)')'STERNHEIMER_RATIONAL rank/residual/tol ',m,maxval(residual),tol
      if(maxval(residual)<=tol)then
        ! Check the original callback on the complete response at every
        ! requested frequency, not just its projected Galerkin equation.
        true_worst=0d0
        do i=1,ns
          xr=matmul(q(:,:m),cmplx(real(coef(:,i),8),0d0,8))
          xi=matmul(q(:,:m),cmplx(aimag(coef(:,i)),0d0,8))
          call scaled_action(xr,ar,code)
          if(code/=0)return
          call scaled_action(xi,ai,code)
          if(code/=0)return
          ar=source/model%source_norm-real(shifts(i),8)*xr+aimag(shifts(i))*xi+ar
          ai=-real(shifts(i),8)*xi-aimag(shifts(i))*xr+ai
          local_error=sum(abs(ar)**2)+sum(abs(ai)**2)
          call comm_summation(local_error,global_error,icomm)
          true_worst=max(true_worst,sqrt(global_error))
        enddo
        model%validation_residual=true_worst
        if(comm_is_root(rank))write(*,'(a,es20.10)')'STERNHEIMER_RATIONAL_ORIGINAL_RESIDUAL ',true_worst
        status=4
        if(true_worst>tol)return
        status=0;return
      endif
      selected=maxloc(residual,dim=1);previous=m
      xr=matmul(q(:,:m),cmplx(real(coef(:,selected),8),0d0,8))
      xi=matmul(q(:,:m),cmplx(aimag(coef(:,selected)),0d0,8))
      ar=matmul(aq(:,:m),cmplx(real(coef(:,selected),8),0d0,8))
      ai=matmul(aq(:,:m),cmplx(aimag(coef(:,selected)),0d0,8))
      xr=source/model%source_norm-real(shifts(selected),8)*xr+aimag(shifts(selected))*xi+ar
      ! Recompute the second quadrature using the original projected vectors.
      xi=-real(shifts(selected),8)*xi-aimag(shifts(selected))* &
        matmul(q(:,:m),cmplx(real(coef(:,selected),8),0d0,8))+ai
      call precondition(xr,candidate,code)
      if(.not.screen_finite(candidate))code=1
      call comm_get_max(code,icomm)
      if(code/=0)return
      call append(candidate,accepted)
      call precondition(xi,candidate,code)
      if(.not.screen_finite(candidate))code=1
      call comm_get_max(code,icomm)
      if(code/=0)return
      call append(candidate,accepted)
      deallocate(small,hlocal,hglobal,coef,defect,real_defect,local_r,stack,stack_sum)
      if(m==previous)exit
    enddo
    status=3
  contains
    subroutine append(vector,added)
      implicit none
      complex(8),intent(in) :: vector(:)
      logical,intent(out) :: added
      complex(8),allocatable :: new_q(:,:),new_aq(:,:)
      real(8),allocatable :: new_g(:)
      integer :: new_capacity,ip
      added=.false.
      if(m==maxrank)return
      if(m==capacity)then
        new_capacity=min(maxrank,2*capacity)
        allocate(new_q(n,new_capacity),new_aq(n,new_capacity),new_g(new_capacity))
        new_q=0d0;new_aq=0d0;new_g=0d0
        new_q(:,:m)=q;new_aq(:,:m)=aq;new_g(:m)=g
        call move_alloc(new_q,q);call move_alloc(new_aq,aq);call move_alloc(new_g,g);capacity=new_capacity
      endif
      action=vector
      do ip=1,2
        dots(:m)=real(matmul(conjg(transpose(q(:,:m))),action),8)
        call comm_summation(dots(:m),total_dots(:m),m,icomm)
        action=action-matmul(q(:,:m),cmplx(total_dots(:m),0d0,8))
      enddo
      local_norm=sum(abs(action)**2);call comm_summation(local_norm,global_norm,icomm)
      if(global_norm<1d-26)return
      m=m+1;q(:,m)=action/sqrt(global_norm);added=.true.
    end subroutine
    subroutine scaled_action(vector,output,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: output(:)
      integer,intent(out) :: ierr
      real(8) :: local_size,global_size
      local_size=sum(abs(vector)**2);call comm_summation(local_size,global_size,icomm)
      output=0d0;ierr=0
      if(global_size<=1d-28)return
      global_size=sqrt(global_size)
      call apply(vector/global_size,output,ierr);output=output*global_size
      if(.not.screen_finite(output))ierr=1
      call comm_get_max(ierr,icomm)
    end subroutine
  end subroutine

  subroutine real_factor(x,r,status)
    implicit none
    real(8),intent(in) :: x(:,:)
    real(8),intent(out) :: r(:,:)
    integer,intent(out) :: status
    real(8),allocatable :: a(:,:),tau(:),work(:)
    integer :: rows,cols,i,j,lwork
    rows=size(x,1);cols=size(x,2);r=0d0;status=0
    if(rows==0.or.cols==0)return
    allocate(a(rows,cols),tau(min(rows,cols)),work(1));a=x
    call dgeqrf(rows,cols,a,rows,tau,work,-1,status)
    if(status/=0)return
    lwork=max(1,int(work(1)));deallocate(work);allocate(work(lwork));a=x
    call dgeqrf(rows,cols,a,rows,tau,work,lwork,status)
    if(status/=0)return
    do j=1,cols
      do i=1,min(rows,j)
        r(i,j)=a(i,j)
      enddo
    enddo
  end subroutine

  ! Finite-time window companion of the SAME resolved linear operator. This
  ! small-space trajectory is diagnostic only; physical RT stays on its grid.
  subroutine impulse_trajectory(model,contact,dt,current,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    type(impulse_model),intent(in) :: model
    real(8),intent(in) :: contact,dt
    real(8),intent(out) :: current(0:)
    integer,intent(out) :: status
    real(8),allocatable :: a(:,:),right(:,:),wr(:),wi(:),work(:)
    real(8) :: unused_left(1,1)
    complex(8),allocatable :: weights(:),eigen(:),r(:),matrix(:,:),rhs(:,: )
    integer,allocatable :: pivots(:)
    complex(8) :: check,expected,shift,total,z,increment,term
    integer :: n,i,j,k,lwork,info
    status=1;current=0d0
    if(.not.allocated(model%h).or..not.allocated(model%g))return
    if(.not.screen_finite([contact,dt]).or.dt<=0d0)return
    n=size(model%g)
    allocate(a(n,n),right(n,n),wr(n),wi(n),work(1),weights(n),eigen(n),r(n), &
      matrix(n,n),rhs(n,1),pivots(n))
    a=model%h
    call dgeev('N','V',n,a,n,wr,wi,unused_left,1,right,n,work,-1,info)
    status=2
    if(info/=0)return
    lwork=max(1,int(work(1)));deallocate(work);allocate(work(lwork));a=model%h
    call dgeev('N','V',n,a,n,wr,wi,unused_left,1,right,n,work,lwork,info)
    if(info/=0)return
    i=1
    do while(i<=n)
      eigen(i)=cmplx(wr(i),wi(i),8)
      if(wi(i)>0d0)then
        if(i==n)return
        r=cmplx(right(:,i),right(:,i+1),8)
      else
        r=cmplx(right(:,i),0d0,8)
      endif
      matrix(:,i)=r;weights(i)=sum(model%g*r)
      if(wi(i)>0d0)then
        eigen(i+1)=conjg(eigen(i));matrix(:,i+1)=conjg(r);weights(i+1)=conjg(weights(i));i=i+2
      else
        i=i+1
      endif
    enddo
    ! Degenerate eigenvalues do not define paired left/right eigenvectors.
    ! Obtain source coefficients from the complete right-eigenvector matrix,
    ! rather than dividing by individual left/right overlaps.
    rhs=0d0;rhs(1,1)=model%source_norm
    call zgesv(n,1,matrix,n,pivots,rhs,n,info)
    status=3
    if(info/=0.or..not.screen_finite(rhs))return
    weights=weights*rhs(:,1)
    if(.not.screen_finite(weights).or..not.screen_finite(eigen))return
    ! Validate the partial fractions with a separate direct shifted solve.
    shift=cmplx(.02d0,-.07d0,8);matrix=-cmplx(model%h,0d0,8)
    do i=1,n
      matrix(i,i)=matrix(i,i)+shift
    enddo
    rhs=0d0;rhs(1,1)=model%source_norm
    call zgesv(n,1,matrix,n,pivots,rhs,n,info)
    status=4
    if(info/=0)return
    expected=sum(model%g*rhs(:,1));check=sum(weights/(shift-eigen))
    status=5
    if(abs(check-expected)>1d-7*max(1d-12,abs(expected)))return
    current(0)=contact
    do j=1,ubound(current,1)
      total=cmplx(contact,0d0,8)
      do i=1,n
        z=eigen(i)*j*dt
        status=6
        if(real(z,8)>600d0)return
        if(abs(z)<.1d0)then
          increment=cmplx(j*dt,0d0,8);term=increment
          do k=1,16
            term=term*z/(k+1);increment=increment+term
          enddo
        else
          increment=(exp(z)-1d0)/eigen(i)
        endif
        total=total+weights(i)*increment
      enddo
      status=7
      if(.not.screen_finite([total]))return
      status=8
      if(abs(aimag(total))>1d-7*max(abs(contact),abs(total),1d-12))return
      current(j)=real(total,8)
    enddo
    status=0
  end subroutine

  logical function select_none(real_part,imaginary_part)
    implicit none
    real(8),intent(in) :: real_part,imaginary_part
    select_none=.false.
  end function
end module
