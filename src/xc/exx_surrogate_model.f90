! Shared real coefficients for complex matrix differences. No Python RT dependency.
module exx_surrogate_model
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::surrogate_fit,surrogate_predict,surrogate_difference
contains
  logical function finite_complex(a)result(ok)
    complex(real64),intent(in)::a(:)
    integer::i
    ok=.false.
    do i=1,size(a)
      if(.not.ieee_is_finite(real(a(i),real64)).or..not.ieee_is_finite(aimag(a(i))))return
    enddo
    ok=.true.
  end function

  subroutine surrogate_fit(x,y,lambda,coeff,status)
    ! x(row,col,feature,example), y(row,col,example).
    complex(real64),intent(in)::x(:,:,:,:),y(:,:,:)
    real(real64),intent(in)::lambda
    real(real64),intent(out)::coeff(:)
    integer,intent(out)::status
    real(real64),allocatable::a(:,:),rhs(:,:),work(:)
    real(real64)::query(1),rmax,rmin
    integer::n,k,ne,m,row,i,j,f,e,info,lwork
    external::dgels
    status=1;coeff=0;n=size(x,1);k=size(x,3);ne=size(x,4)
    if(n<1.or.size(x,2)/=n.or.k<1.or.k>3.or.ne<1)return
    if(any(shape(y)/=[n,n,ne]).or.size(coeff)/=k)return
    if(.not.ieee_is_finite(lambda).or.lambda<0)return
    if(.not.finite_complex(reshape(x,[size(x)])))return
    if(.not.finite_complex(reshape(y,[size(y)])))return
    m=2*n*n*ne+k
    allocate(a(m,k),rhs(m,1));a=0;rhs=0;row=0
    do e=1,ne
      do j=1,n
        do i=1,n
          row=row+1
          a(row,:)=real(x(i,j,:,e),real64);rhs(row,1)=real(y(i,j,e),real64)
          row=row+1
          a(row,:)=aimag(x(i,j,:,e));rhs(row,1)=aimag(y(i,j,e))
        enddo
      enddo
    enddo
    do f=1,k
      a(row+f,f)=sqrt(lambda)
    enddo
    ! Augmented least squares QR avoids squaring the condition number.
    call dgels('N',m,k,1,a,m,rhs,m,query,-1,info)
    if(info/=0.or..not.ieee_is_finite(query(1)))return
    lwork=max(1,int(query(1)));allocate(work(lwork))
    call dgels('N',m,k,1,a,m,rhs,m,work,lwork,info)
    if(info/=0)return
    rmax=0;rmin=huge(1d0)
    do f=1,k
      rmax=max(rmax,abs(a(f,f)));rmin=min(rmin,abs(a(f,f)))
    enddo
    if(rmax<=0.or.rmin<=epsilon(1d0)*max(m,k)*rmax)return
    do f=1,k
      if(.not.ieee_is_finite(rhs(f,1)))return
    enddo
    coeff=rhs(:k,1);status=0
  end subroutine

  subroutine surrogate_difference(history,times,target_time,floor,p2,p4,absolute,relative,status)
    ! Diagnostic tangent estimates only. Disagreement is not an error bound.
    complex(real64),intent(in)::history(:,:,:)
    real(real64),intent(in)::times(:),target_time,floor
    complex(real64),intent(out)::p2(:,:),p4(:,:)
    real(real64),intent(out)::absolute,relative
    integer,intent(out)::status
    integer::n,k
    real(real64)::spacing,scale,n2,n4
    status=1;p2=0;p4=0;absolute=0;relative=0
    n=size(history,1)
    if(size(history,3)/=4.or.size(history,2)/=n.or.n<1.or.size(times)/=4)return
    if(any(shape(p2)/=[n,n]).or.any(shape(p4)/=[n,n]))return
    if(.not.all(ieee_is_finite(times)).or..not.ieee_is_finite(floor).or.floor<=0)return
    spacing=times(2)-times(1)
    if(spacing<=0)return
    if(any(abs((times(2:)-times(:3))-spacing)>1d-10*spacing))return
    if(.not.finite_complex(reshape(history,[size(history)])))return
    scale=max(maxval(abs(history)),1d0)
    do k=1,4
      if(maxval(abs(history(:,:,k)-transpose(conjg(history(:,:,k)))))>1d-12*scale)return
    enddo
    call surrogate_predict(history,times,target_time,[1d0,0d0,0d0],p2,status)
    if(status/=0)return
    call surrogate_predict(history,times,target_time,[11d0,-7d0,2d0]/6d0,p4,status)
    if(status/=0)then
      p2=0;return
    endif
    absolute=sqrt(sum(abs(p4-p2)**2))
    n2=sqrt(sum(abs(p2-history(:,:,4))**2));n4=sqrt(sum(abs(p4-history(:,:,4))**2))
    relative=absolute/max(n2,n4,floor)
    if(.not.all(ieee_is_finite([absolute,relative,n2,n4])))then
      p2=0;p4=0;absolute=0;relative=0;status=1
    endif
    ! PSD projection and strict sampled action checks remain separate.
  end subroutine

  subroutine surrogate_predict(history,times,target_time,coeff,b,status)
    complex(real64),intent(in)::history(:,:,:)
    real(real64),intent(in)::times(:),target_time,coeff(:)
    complex(real64),intent(out)::b(:,:)
    integer,intent(out)::status
    integer::p,n,j,newer
    status=1;b=0;p=size(history,3);n=size(history,1)
    if(p<2.or.p>4.or.n<1.or.size(history,2)/=n)return
    if(size(times)/=p.or.size(coeff)/=p-1.or.any(shape(b)/=[n,n]))return
    if(.not.finite_complex(reshape(history,[size(history)])))return
    if(.not.all(ieee_is_finite(times)).or..not.all(ieee_is_finite(coeff)))return
    if(.not.ieee_is_finite(target_time))return
    if(any(times(2:)<=times(:p-1)).or.target_time<times(p))return
    b=history(:,:,p)
    do j=1,p-1
      newer=p-j+1
      b=b+(target_time-times(p))*coeff(j)*(history(:,:,newer)-history(:,:,newer-1))/ &
        (times(newer)-times(newer-1))
    enddo
    if(.not.finite_complex(reshape(b,[size(b)])))then
      b=0;return
    endif
    ! PSD projection and applied-action residual are separate mandatory checks.
    status=0
  end subroutine
end module
