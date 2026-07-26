module dg_buffer_window_projector
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public :: s_dg_buffer_projector_diagnostics,build_dg_buffer_window_projector,&
    build_dg_buffer_window_projector_dual

  type :: s_dg_buffer_projector_diagnostics
    integer :: configured_states=0
    integer :: retained_rank=0
    real(8) :: minimum_retained_singular_value=0d0
    real(8) :: projection_residual=huge(1d0)
    real(8) :: escape_norm=huge(1d0)
  end type
  interface
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      character(1),intent(in) :: jobz,uplo
      integer,intent(in) :: n,lda,lwork
      real(8),intent(inout) :: a(lda,*)
      real(8),intent(out) :: w(*),work(*)
      integer,intent(out) :: info
    end subroutine
  end interface
contains
  subroutine build_dg_buffer_window_projector_dual(phi_core,weights,rank_tolerance,dual,&
      diagnostics,ok,message)
    real(8),intent(in) :: phi_core(:,:),weights(:),rank_tolerance
    real(8),intent(out) :: dual(:,:)
    type(s_dg_buffer_projector_diagnostics),intent(out) :: diagnostics
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    real(8),allocatable :: metric(:,:),metric_inverse(:,:),eigenvalues(:),work(:),weighted_phi(:,:)
    real(8) :: basis_scale,weight_scale,maximum_singular_value,work_query(1)
    integer :: npoints,nstates,allocation_status,lapack_info,lwork,index
    ok=.false.;message='';dual=0d0;diagnostics=s_dg_buffer_projector_diagnostics()
    npoints=size(phi_core,1);nstates=size(phi_core,2)
    diagnostics%configured_states=nstates
    if(npoints<1.or.nstates<1.or.size(weights)/=npoints.or.&
       any(shape(dual)/=[nstates,npoints]).or.rank_tolerance<=0d0.or.rank_tolerance>1d0.or.&
       .not.ieee_is_finite(rank_tolerance).or..not.all(ieee_is_finite(phi_core)).or.&
       .not.all(ieee_is_finite(weights)).or.any(weights<0d0).or.maxval(weights)<=0d0)then
      message='buffer projector dual requires valid finite basis, weights, and dimensions';return
    endif
    basis_scale=maxval(abs(phi_core));weight_scale=maxval(weights)
    if(basis_scale<=0d0)then
      message='buffer projector dual overlap metric has no positive retained rank';return
    endif
    allocate(metric(nstates,nstates),metric_inverse(nstates,nstates),eigenvalues(nstates),&
      weighted_phi(npoints,nstates),stat=allocation_status)
    if(allocation_status/=0)then;message='buffer projector dual workspace allocation failed';return;endif
    weighted_phi=(phi_core/basis_scale)*spread(weights/weight_scale,2,nstates)
    metric=matmul(transpose(weighted_phi),phi_core/basis_scale)
    call dsyev('V','U',nstates,metric,nstates,eigenvalues,work_query,-1,lapack_info)
    if(lapack_info/=0.or..not.ieee_is_finite(work_query(1)))then
      message='buffer projector dual eigensolver workspace query failed';return
    endif
    lwork=max(1,ceiling(work_query(1)));allocate(work(lwork),stat=allocation_status)
    if(allocation_status/=0)then;message='buffer projector dual eigensolver allocation failed';return;endif
    call dsyev('V','U',nstates,metric,nstates,eigenvalues,work,lwork,lapack_info)
    if(lapack_info/=0.or..not.all(ieee_is_finite(metric)).or..not.all(ieee_is_finite(eigenvalues)))then
      message='buffer projector dual overlap eigendecomposition failed';return
    endif
    maximum_singular_value=maxval(eigenvalues)
    if(maximum_singular_value<=0d0)then
      message='buffer projector dual overlap metric has no positive retained rank';return
    endif
    diagnostics%retained_rank=count(eigenvalues>=rank_tolerance*maximum_singular_value)
    if(diagnostics%retained_rank<1)then
      message='buffer projector dual rank tolerance removed the complete state window';return
    endif
    diagnostics%minimum_retained_singular_value=restore_metric_scale(minval(eigenvalues,&
      mask=eigenvalues>=rank_tolerance*maximum_singular_value),basis_scale,weight_scale)
    metric_inverse=metric
    do index=1,nstates
      if(eigenvalues(index)>=rank_tolerance*maximum_singular_value)then
        metric_inverse(:,index)=metric_inverse(:,index)/eigenvalues(index)
      else
        metric_inverse(:,index)=0d0
      endif
    enddo
    metric_inverse=matmul(metric_inverse,transpose(metric))
    dual=matmul(metric_inverse,transpose(weighted_phi))/basis_scale
    if(.not.all(ieee_is_finite(dual)))then
      dual=0d0;message='buffer projector dual produced nonfinite output';return
    endif
    diagnostics%projection_residual=0d0;diagnostics%escape_norm=0d0
    ok=.true.
  end subroutine build_dg_buffer_window_projector_dual

  subroutine build_dg_buffer_window_projector(phi_core,buffer_values,weights, &
      rank_tolerance,coefficients,reconstructed,diagnostics,ok,message)
    real(8),intent(in) :: phi_core(:,:),buffer_values(:,:),weights(:),rank_tolerance
    real(8),intent(out) :: coefficients(:,:),reconstructed(:,:)
    type(s_dg_buffer_projector_diagnostics),intent(out) :: diagnostics
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    real(8),allocatable :: metric(:,:),right_hand_side(:,:),eigenvalues(:),work(:)
    real(8),allocatable :: weighted_phi(:,:),eigen_rhs(:,:),normalized_buffer(:,:),normalized_weights(:)
    real(8) :: maximum_singular_value,residual_squared,buffer_norm_squared,work_query(1)
    real(8) :: basis_scale,buffer_scale,weight_scale,coefficient_scale
    integer :: allocation_status,lapack_info,lwork,npoints,nstates,nvectors,i

    ok=.false.
    message=''
    diagnostics=s_dg_buffer_projector_diagnostics()
    coefficients=0d0
    reconstructed=0d0
    npoints=size(phi_core,1)
    nstates=size(phi_core,2)
    nvectors=size(buffer_values,2)
    diagnostics%configured_states=nstates
    if(npoints<=0.or.nstates<=0.or.nvectors<=0)then
      message='buffer projector requires a nonempty point, state, and vector window'
      return
    endif
    if(size(buffer_values,1)/=npoints.or.size(weights)/=npoints.or.&
       size(coefficients,1)/=nstates.or.size(coefficients,2)/=nvectors.or.&
       size(reconstructed,1)/=npoints.or.size(reconstructed,2)/=nvectors)then
      message='buffer projector dimension mismatch'
      return
    endif
    if(.not.ieee_is_finite(rank_tolerance))then
      message='buffer projector rank tolerance must be finite, positive, and at most one'
      return
    endif
    if(rank_tolerance<=0d0.or.rank_tolerance>1d0)then
      message='buffer projector rank tolerance must be finite, positive, and at most one'
      return
    endif
    if(.not.all(ieee_is_finite(phi_core)).or..not.all(ieee_is_finite(buffer_values)).or.&
       .not.all(ieee_is_finite(weights)))then
      message='buffer projector inputs must be finite'
      return
    endif
    if(any(weights<0d0))then
      message='buffer projector point weight must be nonnegative'
      return
    endif
    if(maxval(weights)<=0d0)then
      message='buffer projector requires a positive total point weight'
      return
    endif

    allocate(metric(nstates,nstates),right_hand_side(nstates,nvectors),eigenvalues(nstates),&
      weighted_phi(npoints,nstates),eigen_rhs(nstates,nvectors),normalized_buffer(npoints,nvectors),&
      normalized_weights(npoints),stat=allocation_status)
    if(allocation_status/=0)then
      message='buffer projector workspace allocation failed'
      return
    endif
    basis_scale=maxval(abs(phi_core))
    buffer_scale=maxval(abs(buffer_values))
    weight_scale=maxval(weights)
    if(basis_scale<=0d0.or.weight_scale<=0d0)then
      message='buffer projector overlap metric has no positive retained rank'
      return
    endif
    if(buffer_scale<=0d0)buffer_scale=1d0
    normalized_weights=weights/weight_scale
    normalized_buffer=buffer_values/buffer_scale
    weighted_phi=(phi_core/basis_scale)*spread(normalized_weights,2,nstates)
    metric=matmul(transpose(weighted_phi),phi_core/basis_scale)
    right_hand_side=matmul(transpose(weighted_phi),normalized_buffer)
    if(.not.all(ieee_is_finite(metric)).or..not.all(ieee_is_finite(right_hand_side)))then
      message='buffer projector scaled overlap assembly produced nonfinite values'
      return
    endif
    call dsyev('V','U',nstates,metric,nstates,eigenvalues,work_query,-1,lapack_info)
    if(lapack_info/=0.or..not.ieee_is_finite(work_query(1)))then
      message='buffer projector eigensolver workspace query failed'
      return
    endif
    lwork=max(1,ceiling(work_query(1)))
    allocate(work(lwork),stat=allocation_status)
    if(allocation_status/=0)then
      message='buffer projector eigensolver workspace allocation failed'
      return
    endif
    call dsyev('V','U',nstates,metric,nstates,eigenvalues,work,lwork,lapack_info)
    if(lapack_info/=0.or..not.all(ieee_is_finite(eigenvalues)).or.&
       .not.all(ieee_is_finite(metric)))then
      message='buffer projector overlap eigendecomposition failed'
      return
    endif
    maximum_singular_value=maxval(eigenvalues)
    if(maximum_singular_value<=0d0)then
      message='buffer projector overlap metric has no positive retained rank'
      return
    endif
    diagnostics%retained_rank=count(eigenvalues>=rank_tolerance*maximum_singular_value)
    if(diagnostics%retained_rank<=0)then
      message='buffer projector rank tolerance removed the complete state window'
      return
    endif
    diagnostics%minimum_retained_singular_value=restore_metric_scale(minval(eigenvalues,&
      mask=eigenvalues>=rank_tolerance*maximum_singular_value),basis_scale,weight_scale)
    eigen_rhs=matmul(transpose(metric),right_hand_side)
    do i=1,nstates
      if(eigenvalues(i)>=rank_tolerance*maximum_singular_value)then
        eigen_rhs(i,:)=eigen_rhs(i,:)/eigenvalues(i)
      else
        eigen_rhs(i,:)=0d0
      endif
    enddo
    coefficient_scale=buffer_scale/basis_scale
    if(.not.ieee_is_finite(coefficient_scale))then
      message='buffer projector coefficient scaling overflowed'
      return
    endif
    coefficients=coefficient_scale*matmul(metric,eigen_rhs)
    reconstructed=matmul(phi_core,coefficients)
    residual_squared=sum(spread(normalized_weights,2,nvectors)*&
      ((buffer_values-reconstructed)/buffer_scale)**2)
    buffer_norm_squared=sum(spread(normalized_weights,2,nvectors)*normalized_buffer**2)
    if(buffer_norm_squared<=0d0)then
      diagnostics%projection_residual=sqrt(max(0d0,residual_squared))
    else
      diagnostics%projection_residual=sqrt(max(0d0,residual_squared)/buffer_norm_squared)
    endif
    diagnostics%escape_norm=restore_escape_scale(residual_squared,buffer_scale,weight_scale)
    if(.not.all(ieee_is_finite(coefficients)).or..not.all(ieee_is_finite(reconstructed)).or.&
       .not.ieee_is_finite(diagnostics%projection_residual).or.&
       .not.ieee_is_finite(diagnostics%escape_norm))then
      message='buffer projector produced nonfinite output'
      return
    endif
    ok=.true.
  end subroutine

  real(8) function restore_metric_scale(value,basis_scale,weight_scale)result(restored)
    real(8),intent(in) :: value,basis_scale,weight_scale
    real(8) :: logarithm
    if(value<=0d0)then
      restored=0d0
      return
    endif
    logarithm=log(value)+2d0*log(basis_scale)+log(weight_scale)
    if(logarithm>=log(huge(1d0)))then
      restored=huge(1d0)
    else if(logarithm<=log(tiny(1d0)))then
      restored=0d0
    else
      restored=exp(logarithm)
    endif
  end function

  real(8) function restore_escape_scale(residual_squared,buffer_scale,weight_scale)result(restored)
    real(8),intent(in) :: residual_squared,buffer_scale,weight_scale
    real(8) :: logarithm
    if(residual_squared<=0d0)then
      restored=0d0
      return
    endif
    logarithm=0.5d0*log(residual_squared)+log(buffer_scale)+0.5d0*log(weight_scale)
    if(logarithm>=log(huge(1d0)))then
      restored=huge(1d0)
    else if(logarithm<=log(tiny(1d0)))then
      restored=0d0
    else
      restored=exp(logarithm)
    endif
  end function
end module
