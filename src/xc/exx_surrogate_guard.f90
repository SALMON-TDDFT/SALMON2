! Prediction-only corrections and serial/global weighted action checks.
module exx_surrogate_guard
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::surrogate_psd,surrogate_action_check
contains
  logical function finite(a)result(ok)
    complex(real64),intent(in)::a(:,:)
    integer::i,j
    ok=.false.
    do j=1,size(a,2)
      do i=1,size(a,1)
        if(.not.ieee_is_finite(real(a(i,j),real64)).or..not.ieee_is_finite(aimag(a(i,j))))return
      enddo
    enddo
    ok=.true.
  end function

  subroutine surrogate_psd(b,correction_max,floor,correction,status)
    complex(real64),intent(inout)::b(:,:)
    real(real64),intent(in)::correction_max,floor
    real(real64),intent(out)::correction
    integer,intent(out)::status
    complex(real64),allocatable::vectors(:,:),corrected(:,:),work(:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    integer::n,j,info
    external::zheev
    status=1;correction=huge(1d0);n=size(b,1)
    if(n<1.or.size(b,2)/=n.or..not.finite(b))return
    if(.not.ieee_is_finite(correction_max).or.correction_max<0)return
    if(.not.ieee_is_finite(floor).or.floor<=0)return
    allocate(vectors(n,n),corrected(n,n),eigenvalues(n),work(max(1,2*n-1)),rwork(max(1,3*n-2)))
    vectors=.5d0*(b+conjg(transpose(b)))
    call zheev('V','U',n,vectors,n,eigenvalues,work,size(work),rwork,info)
    if(info/=0)return
    corrected=0
    do j=1,n
      corrected=corrected+max(eigenvalues(j),0d0)* &
        spread(vectors(:,j),2,n)*spread(conjg(vectors(:,j)),1,n)
    enddo
    if(.not.finite(corrected))return
    correction=sqrt(sum(abs(corrected-b)**2))/max(sqrt(sum(abs(b)**2)),floor)
    if(.not.ieee_is_finite(correction).or.correction>correction_max)return
    b=corrected;status=0
  end subroutine

  subroutine surrogate_action_check(q,b,psi,exact,dv,atol,rtol,floor,absolute,relative,accepted,status)
    ! Serial/global assembled rows only. Distributed reduction is not wired yet.
    complex(real64),intent(in)::q(:,:),b(:,:),psi(:,:),exact(:,:)
    real(real64),intent(in)::dv,atol,rtol,floor
    real(real64),intent(out)::absolute,relative
    logical,intent(out)::accepted
    integer,intent(out)::status
    complex(real64),allocatable::predicted(:,:),gram(:,:)
    real(real64)::a,e,norm
    integer::i,j,k
    status=1;accepted=.false.;absolute=huge(1d0);relative=huge(1d0);k=size(q,2)
    if(k<1.or.any(shape(b)/=[k,k]).or.size(psi,1)/=size(q,1))return
    if(any(shape(exact)/=shape(psi)).or.size(psi,2)<1)return
    if(.not.finite(q).or..not.finite(b).or..not.finite(psi).or..not.finite(exact))return
    if(.not.all(ieee_is_finite([dv,atol,rtol,floor])))return
    if(dv<=0.or.atol<0.or.rtol<0.or.floor<=0)return
    gram=matmul(conjg(transpose(q)),q)*dv
    do i=1,k
      gram(i,i)=gram(i,i)-1
    enddo
    if(maxval(abs(gram))>1d-10)return
    predicted=-matmul(q,matmul(b,matmul(conjg(transpose(q)),psi)*dv))
    absolute=0;relative=0
    do j=1,size(psi,2)
      a=sqrt(sum(abs(predicted(:,j)-exact(:,j))**2)*dv)
      norm=sqrt(sum(abs(exact(:,j))**2)*dv);e=a/max(norm,floor)
      absolute=max(absolute,a);relative=max(relative,e)
    enddo
    if(.not.all(ieee_is_finite([absolute,relative])))return
    accepted=absolute<=atol.and.relative<=rtol;status=0
  end subroutine
end module
