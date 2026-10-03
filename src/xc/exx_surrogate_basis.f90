! Storage-neutral projection of the current exact ACE action onto fixed Q.
! Q contains local spatial rows, globally replicated columns, Gamma only.
module exx_surrogate_basis
  use iso_fortran_env,only: real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::surrogate_project,surrogate_build_q,sum_callback
  abstract interface
    subroutine action_callback(target,result,status)
      import real64
      complex(real64),intent(in)::target(:,:)
      complex(real64),intent(out)::result(:,:)
      integer,intent(out)::status
    end subroutine
    subroutine sum_callback(matrix)
      import real64
      complex(real64),intent(inout)::matrix(:,:)
    end subroutine
  end interface
contains
  logical function finite_matrix(a)result(ok)
    complex(real64),intent(in)::a(:,:)
    integer::i,j
    ok=.false.
    do j=1,size(a,2)
      do i=1,size(a,1)
        if(.not.ieee_is_finite(real(a(i,j),real64)))return
        if(.not.ieee_is_finite(aimag(a(i,j))))return
      enddo
    enddo
    ok=.true.
  end function

  subroutine surrogate_build_q(columns,dv,rank_max,rtol,q,rank,status,sum_grid)
    complex(real64),intent(in)::columns(:,:)
    real(real64),intent(in)::dv,rtol
    integer,intent(in)::rank_max
    complex(real64),allocatable,intent(out)::q(:,:)
    integer,intent(out)::rank,status
    procedure(sum_callback),optional::sum_grid
    complex(real64),allocatable::work(:,:),v(:)
    complex(real64)::scalar(1,1),bad(1,1)
    real(real64)::norm,scale
    integer::j,k,pass
    status=1;rank=0;bad=0
    if(rank_max<1.or.size(columns,2)<1)bad=1
    if(.not.ieee_is_finite(dv).or.dv<=0)bad=1
    if(.not.ieee_is_finite(rtol).or.rtol<=0.or.rtol>=1)bad=1
    if(.not.finite_matrix(columns))bad=1
    if(present(sum_grid))call sum_grid(bad)
    if(abs(bad(1,1))>0)return
    allocate(work(size(columns,1),min(rank_max,size(columns,2))),v(size(columns,1)))
    work=0;scale=0
    do j=1,size(columns,2)
      scalar(1,1)=sum(abs(columns(:,j))**2)*dv
      if(present(sum_grid))call sum_grid(scalar)
      if(.not.ieee_is_finite(real(scalar(1,1),real64)))return
      scale=max(scale,sqrt(max(real(scalar(1,1),real64),0d0)))
    enddo
    if(scale<=0)return
    do j=1,size(columns,2)
      v=columns(:,j)
      ! Two-pass modified Gram-Schmidt, using globally reduced inner products.
      do pass=1,2
        do k=1,rank
          scalar(1,1)=sum(conjg(work(:,k))*v)*dv
          if(present(sum_grid))call sum_grid(scalar)
          v=v-work(:,k)*scalar(1,1)
        enddo
      enddo
      scalar(1,1)=sum(abs(v)**2)*dv
      if(present(sum_grid))call sum_grid(scalar)
      if(.not.ieee_is_finite(real(scalar(1,1),real64)))return
      norm=sqrt(max(real(scalar(1,1),real64),0d0))
      if(norm<=rtol*scale)cycle
      ! Do not silently discard independent directions to meet a memory cap.
      if(rank==size(work,2))return
      rank=rank+1;work(:,rank)=v/norm
    enddo
    if(rank==0)return
    allocate(q(size(columns,1),rank));q=work(:,:rank);status=0
  end subroutine

  subroutine surrogate_project(q,dv,apply,b,status,sum_grid)
    complex(real64),intent(in)::q(:,:)
    real(real64),intent(in)::dv
    procedure(action_callback)::apply
    complex(real64),intent(out)::b(:,:)
    integer,intent(out)::status
    procedure(sum_callback),optional::sum_grid
    complex(real64),allocatable::gram(:,:),w(:,:)
    complex(real64)::bad(1,1)
    real(real64)::scale
    integer::i,n,action_status
    status=1;b=0;n=size(q,2);bad=0
    if(n<1.or.any(shape(b)/=[n,n]))bad=1
    if(.not.ieee_is_finite(dv).or.dv<=0)bad=1
    if(.not.finite_matrix(q))bad=1
    if(present(sum_grid))call sum_grid(bad)
    if(abs(bad(1,1))>0)return
    allocate(gram(n,n),w(size(q,1),n))
    gram=matmul(conjg(transpose(q)),q)*dv
    if(present(sum_grid))call sum_grid(gram)
    do i=1,n
      gram(i,i)=gram(i,i)-1
    enddo
    if(maxval(abs(gram))>1d-10)return
    w=0
    call apply(q,w,action_status)
    bad=0
    if(action_status/=0)bad=1
    if(.not.finite_matrix(w))bad=1
    if(present(sum_grid))call sum_grid(bad)
    if(abs(bad(1,1))>0)return
    b=-matmul(conjg(transpose(q)),w)*dv
    if(present(sum_grid))call sum_grid(b)
    if(.not.finite_matrix(b))then
      b=0;return
    endif
    scale=max(maxval(abs(b)),1d0)
    if(maxval(abs(b-conjg(transpose(b))))>1d-10*scale)then
      b=0;return
    endif
    ! Do not symmetrize/PSD-project a rejected exact action.
    ! PSD certification is a separate spectral diagnostic, not done here.
    status=0
  end subroutine
end module
