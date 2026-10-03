! Serial Gamma adapter. Distributed/packed ACE must use orbital_ace_apply instead.
module exx_surrogate_dense
  use iso_c_binding,only:c_double_complex
  use exx_ace,only:s_exx_ace,exx_ace_apply
  use exx_surrogate_basis,only:surrogate_project,surrogate_build_q,sum_callback
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::surrogate_project_dense,surrogate_dense_snapshot,surrogate_dense_basis_error
contains
  ! Relative factor leakage and ACE action error on all factor columns.
  ! These are representation diagnostics, not a strict-Fock sample certificate.
  subroutine surrogate_dense_basis_error(ace,q,eta,action_error,status)
    type(s_exx_ace),intent(in)::ace
    complex(c_double_complex),intent(in)::q(:,:)
    real(8),intent(out)::eta,action_error
    integer,intent(out)::status
    complex(c_double_complex),allocatable::x(:,:),c(:,:),b(:,:),exact(:,:),projected(:,:),gram(:,:)
    real(8)::norm_x,norm_action
    integer::i
    status=1;eta=huge(1d0);action_error=huge(1d0)
    if(ace%packed.or.ace%metric_distributed.or..not.allocated(ace%factors))return
    if(size(ace%factors,3)/=1.or.size(q,1)/=size(ace%factors,1).or.size(q,2)<1)return
    if(.not.ieee_is_finite(ace%dv).or.ace%dv<=0)return
    if(.not.all(ieee_is_finite(real(q))).or..not.all(ieee_is_finite(aimag(q))))return
    x=ace%factors(:,:,1)
    if(.not.all(ieee_is_finite(real(x))).or..not.all(ieee_is_finite(aimag(x))))return
    gram=matmul(conjg(transpose(q)),q)*ace%dv
    do i=1,size(q,2)
      gram(i,i)=gram(i,i)-1
    enddo
    if(maxval(abs(gram))>1d-10)return
    norm_x=sqrt(sum(abs(x)**2));if(norm_x<=tiny(1d0))return
    c=matmul(conjg(transpose(q)),x)*ace%dv
    eta=sqrt(sum(abs(x-matmul(q,c))**2))/norm_x
    b=matmul(c,conjg(transpose(c)))
    exact=-matmul(x,matmul(conjg(transpose(x)),x))*ace%dv
    projected=-matmul(q,matmul(b,c))
    norm_action=sqrt(sum(abs(exact)**2));if(norm_action<=tiny(1d0))return
    action_error=sqrt(sum(abs(exact-projected)**2))/norm_action
    if(.not.all(ieee_is_finite([eta,action_error])))return
    status=0
  end subroutine

  subroutine surrogate_dense_snapshot(ace,rank_max,rtol,q,b,rank,status,sum_grid)
    type(s_exx_ace),intent(in)::ace
    integer,intent(in)::rank_max
    real(8),intent(in)::rtol
    complex(c_double_complex),allocatable,intent(out)::q(:,:),b(:,:)
    integer,intent(out)::rank,status
    procedure(sum_callback),optional::sum_grid
    status=1;rank=0
    if(ace%packed.or.ace%metric_distributed)return
    if(.not.allocated(ace%factors))return
    if(size(ace%factors,3)/=1)return
    call surrogate_build_q(ace%factors(:,:,1),ace%dv,rank_max,rtol,q,rank,status,sum_grid)
    if(status/=0)return
    allocate(b(rank,rank))
    call surrogate_project_dense(ace,q,b,status,sum_grid)
    if(status/=0)then
      deallocate(q,b);rank=0
    endif
  end subroutine

  subroutine surrogate_project_dense(ace,q,b,status,sum_grid)
    type(s_exx_ace),intent(in)::ace
    complex(c_double_complex),intent(in)::q(:,:)
    complex(c_double_complex),intent(out)::b(:,:)
    integer,intent(out)::status
    procedure(sum_callback),optional::sum_grid
    status=1;b=0
    if(ace%packed.or.ace%metric_distributed)return
    if(.not.allocated(ace%factors))return
    if(size(ace%factors,3)/=1.or.size(ace%factors,1)/=size(q,1))return
    call surrogate_project(q,ace%dv,apply_current,b,status,sum_grid)
  contains
    subroutine apply_current(target,result,ierr)
      complex(c_double_complex),intent(in)::target(:,:)
      complex(c_double_complex),intent(out)::result(:,:)
      integer,intent(out)::ierr
      complex(c_double_complex),allocatable::target3(:,:,:),action3(:,:,:)
      allocate(target3(size(target,1),size(target,2),1),action3(size(target,1),size(target,2),1))
      target3(:,:,1)=target
      call exx_ace_apply(ace,target3,action3,ierr,sum_grid)
      result=action3(:,:,1)
    end subroutine
  end subroutine
end module
