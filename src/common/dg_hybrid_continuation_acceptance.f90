#include "config.h"
module dg_hybrid_continuation_acceptance
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_continuation_residuals,only:validate_cluster_occupations
#ifdef USE_MPI
  use mpi, only: MPI_Allreduce, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_MIN, MPI_SUCCESS, MPI_SUM
#endif
  implicit none
  private

  type,public::s_dg_hybrid_acceptance_result
    logical::valid=.false.,identity_only=.false.
    logical::symmetry_complete=.false.,grid_complete=.false.,face_complete=.false.
    logical::excitation_cutoff_convergence_proven=.false.
    integer::expected_occupied_count=0,checked_occupied_count=0,checked_face_count=0
    integer(int64)::analysis_fingerprint=0_int64,basis_fingerprint=0_int64
    integer(int64)::operator_fingerprint=0_int64,state_fingerprint=0_int64,action_fingerprint=0_int64
    integer(int64)::action_operator_fingerprint=0_int64,face_topology_fingerprint=0_int64
    real(real64)::volume_covariance=huge(1d0),interface_covariance=huge(1d0)
    real(real64)::metric_covariance=huge(1d0),retained_leakage=huge(1d0)
    real(real64)::occupied_covariance=huge(1d0),grid_residual=huge(1d0)
    real(real64)::face_residual(3)=huge(1d0)
  end type s_dg_hybrid_acceptance_result

  public::evaluate_dg_hybrid_acceptance,validate_dg_hybrid_acceptance_receipt
contains
  subroutine validate_dg_hybrid_acceptance_receipt(icomm,receipt,expected_occupied_count,expected_face_count,&
      expected_operator_fingerprint,expected_face_topology_fingerprint,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_acceptance_result),intent(in)::receipt
    integer,intent(in)::expected_occupied_count,expected_face_count
    integer(int64),intent(in)::expected_operator_fingerprint,expected_face_topology_fingerprint
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr,minimum_count,maximum_count,minimum_face_count,maximum_face_count
    integer(int64)::values(7),minimum_values(7),maximum_values(7)
    local_bad=merge(0,1,receipt%valid.and.receipt%symmetry_complete.and.receipt%grid_complete.and.&
      receipt%face_complete.and.receipt%expected_occupied_count>0.and.&
      receipt%expected_occupied_count==expected_occupied_count.and.&
      receipt%checked_occupied_count==expected_occupied_count.and.receipt%checked_face_count==expected_face_count.and.&
      receipt%operator_fingerprint==expected_operator_fingerprint.and.&
      receipt%action_operator_fingerprint==expected_operator_fingerprint.and.&
      receipt%face_topology_fingerprint==expected_face_topology_fingerprint.and.&
      receipt%analysis_fingerprint/=0_int64.and.&
      receipt%basis_fingerprint/=0_int64.and.receipt%operator_fingerprint/=0_int64.and.&
      receipt%state_fingerprint/=0_int64.and.receipt%action_fingerprint/=0_int64)
    if(expected_face_count<0.or.expected_face_topology_fingerprint==0_int64.or.&
        receipt%face_topology_fingerprint==0_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    call MPI_Allreduce(receipt%checked_occupied_count,minimum_count,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(receipt%checked_occupied_count,maximum_count,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    call MPI_Allreduce(receipt%checked_face_count,minimum_face_count,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(receipt%checked_face_count,maximum_face_count,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    values=[receipt%analysis_fingerprint,receipt%basis_fingerprint,receipt%operator_fingerprint,&
      receipt%state_fingerprint,receipt%action_fingerprint,receipt%action_operator_fingerprint,&
      receipt%face_topology_fingerprint]
    call MPI_Allreduce(values,minimum_values,7,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(values,maximum_values,7,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0.and.minimum_count==maximum_count.and.&
      minimum_face_count==maximum_face_count.and.all(minimum_values==maximum_values)
    if(ok)then;message='';else;message='incomplete or rank-disagreeing DG acceptance receipt';endif
#else
    ok=.false.;message='DG acceptance receipt validation requires MPI'
#endif
  end subroutine validate_dg_hybrid_acceptance_receipt

  subroutine evaluate_dg_hybrid_acceptance(icomm,analysis_complete,analysis_fingerprint,expected_operation_count, &
      expected_nonidentity_count,authoritative_identity_only,representation, &
      inverse_representation,h_volume,h_interface,s_metric,q_retained,gamma_occupied,occupations,cluster_ids, &
      face_lambda,h_grid,s_grid_epsilon,face_terms,grid_weights,face_weights,tolerance,result,ok,message)
    integer,intent(in)::icomm
    logical,intent(in)::analysis_complete
    integer(int64),intent(in)::analysis_fingerprint
    integer,intent(in)::expected_operation_count,expected_nonidentity_count
    logical,intent(in)::authoritative_identity_only
    complex(real64),intent(in)::representation(:,:,:),inverse_representation(:,:,:)
    complex(real64),intent(in)::h_volume(:,:),h_interface(:,:),s_metric(:,:),q_retained(:,:),gamma_occupied(:,:)
    real(real64),intent(in)::occupations(:),face_lambda(:),face_terms(:,:),grid_weights(:),face_weights(:),tolerance
    integer,intent(in)::cluster_ids(:)
    complex(real64),intent(in)::h_grid(:),s_grid_epsilon(:)
    type(s_dg_hybrid_acceptance_result),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::identity(:,:),work(:,:),inverse_metric(:,:)
    real(real64)::local_defects(9),global_defects(9),local_grid(3),global_grid(3),scale
    integer::n,noperation,i,local_bad,global_bad,ierr,identity_integer
    integer::minimum_operation_count,maximum_operation_count,minimum_nonidentity_count,maximum_nonidentity_count
    integer::minimum_identity,maximum_identity
    integer(int64)::minimum_fingerprint,maximum_fingerprint
    real(real64)::local_lambda_min,local_lambda_max,global_lambda_min,global_lambda_max
    logical::occupation_ok
    character(256)::occupation_message

    result=s_dg_hybrid_acceptance_result();ok=.false.;message=''
    n=size(h_volume,1);noperation=size(representation,3);local_bad=0
    if(.not.analysis_complete.or.analysis_fingerprint==0_int64.or.n<1.or.noperation<1.or.&
        expected_operation_count/=noperation.or.expected_nonidentity_count/=noperation-1.or.&
        (authoritative_identity_only.neqv.(expected_nonidentity_count==0)).or.tolerance<=0d0.or.&
        .not.ieee_is_finite(tolerance).or.any(shape(h_volume)/=[n,n]).or.any(shape(h_interface)/=[n,n]).or.&
        any(shape(s_metric)/=[n,n]).or.any(shape(q_retained)/=[n,n]).or.any(shape(gamma_occupied)/=[n,n]).or.&
        any(shape(representation)/=[n,n,noperation]).or.any(shape(inverse_representation)/=[n,n,noperation]).or.&
        size(face_lambda)<1.or.size(h_grid)<1.or.size(s_grid_epsilon)/=size(h_grid).or.&
        size(grid_weights)/=size(h_grid).or.size(face_weights)/=size(face_lambda).or.&
        any(shape(face_terms)/=[3,size(face_lambda)]).or.&
        .not.finite_matrix(h_volume).or..not.finite_matrix(h_interface).or..not.finite_matrix(s_metric).or.&
        .not.finite_matrix(q_retained).or..not.finite_matrix(gamma_occupied).or.&
        .not.metric_is_spd(s_metric).or.&
        .not.finite_tensor(representation).or..not.finite_tensor(inverse_representation).or.&
        .not.all(ieee_is_finite(face_lambda)).or..not.all(ieee_is_finite(face_terms)).or.&
        .not.all(ieee_is_finite(grid_weights)).or.any(grid_weights<=0d0).or.&
        .not.all(ieee_is_finite(face_weights)).or.any(face_weights<=0d0))local_bad=1
    call validate_cluster_occupations(occupations,cluster_ids,occupation_ok,occupation_message)
    if(.not.occupation_ok)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    call MPI_Allreduce(expected_operation_count,minimum_operation_count,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(expected_operation_count,maximum_operation_count,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    call MPI_Allreduce(expected_nonidentity_count,minimum_nonidentity_count,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(expected_nonidentity_count,maximum_nonidentity_count,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    call MPI_Allreduce(analysis_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(analysis_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    identity_integer=merge(1,0,authoritative_identity_only)
    call MPI_Allreduce(identity_integer,minimum_identity,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(identity_integer,maximum_identity,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(minimum_operation_count/=maximum_operation_count.or.minimum_nonidentity_count/=maximum_nonidentity_count.or.&
        minimum_fingerprint/=maximum_fingerprint.or.minimum_identity/=maximum_identity)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid or incomplete authoritative DG acceptance payload';return
    endif
    allocate(identity(n,n),work(n,n),inverse_metric(n,n));identity=(0d0,0d0)
    do i=1,n;identity(i,i)=1d0;enddo
    call invert_matrix(s_metric,inverse_metric,occupation_ok)
    local_bad=merge(0,1,occupation_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='DG metric inversion failed collectively';return;endif
    local_defects=0d0
    local_defects(1)=normalized_defect(representation(:,:,1),identity)
    do i=1,noperation
      local_defects(2)=max(local_defects(2),bilinear_defect(&
        matmul(conjg(transpose(representation(:,:,i))),matmul(h_volume,representation(:,:,i))),h_volume,inverse_metric))
      local_defects(3)=max(local_defects(3),bilinear_defect(&
        matmul(conjg(transpose(representation(:,:,i))),matmul(h_interface,representation(:,:,i))),h_interface,inverse_metric))
      local_defects(4)=max(local_defects(4),bilinear_defect(&
        matmul(conjg(transpose(representation(:,:,i))),matmul(s_metric,representation(:,:,i))),s_metric,inverse_metric))
      work=matmul(representation(:,:,i),matmul(q_retained,inverse_representation(:,:,i)))
      local_defects(5)=max(local_defects(5),mixed_defect(work,q_retained,s_metric,inverse_metric))
      work=matmul(identity-q_retained,matmul(representation(:,:,i),q_retained))
      local_defects(6)=max(local_defects(6),mixed_norm(work,s_metric,inverse_metric))
      work=matmul(representation(:,:,i),matmul(gamma_occupied,conjg(transpose(representation(:,:,i)))))
      local_defects(7)=max(local_defects(7),gamma_defect(work,gamma_occupied,s_metric))
      local_defects(8)=max(local_defects(8),normalized_defect(&
        matmul(representation(:,:,i),inverse_representation(:,:,i)),identity))
    enddo
    local_defects(9)=maxval(abs(face_lambda-face_lambda(1)))/max(1d0,abs(face_lambda(1)))
    local_lambda_min=minval(face_lambda);local_lambda_max=maxval(face_lambda)
    call MPI_Allreduce(local_lambda_min,global_lambda_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierr)
    call MPI_Allreduce(local_lambda_max,global_lambda_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierr)
    local_defects(9)=max(local_defects(9),abs(global_lambda_max-global_lambda_min)/max(1d0,abs(global_lambda_min)))
    call MPI_Allreduce(local_defects,global_defects,9,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierr)
    local_grid=[sum(grid_weights*abs(h_grid-s_grid_epsilon)**2),&
      sum(grid_weights*abs(h_grid)**2),sum(grid_weights*abs(s_grid_epsilon)**2)]
    call MPI_Allreduce(local_grid,global_grid,3,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)
    scale=max(1d0,sqrt(global_grid(2)),sqrt(global_grid(3)))
    result%grid_residual=sqrt(global_grid(1))/scale
    do i=1,3
      local_grid(1)=sum(face_weights*face_terms(i,:)**2)
      call MPI_Allreduce(local_grid(1),global_grid(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)
      result%face_residual(i)=sqrt(global_grid(1))
    enddo
    result%analysis_fingerprint=analysis_fingerprint;result%identity_only=authoritative_identity_only
    result%volume_covariance=global_defects(2);result%interface_covariance=global_defects(3)
    result%metric_covariance=global_defects(4);result%retained_leakage=max(global_defects(5),global_defects(6))
    result%occupied_covariance=global_defects(7)
    result%valid=ierr==MPI_SUCCESS.and.maxval(global_defects)<=tolerance.and.&
      result%grid_residual<=tolerance.and.maxval(result%face_residual)<=tolerance
    ok=result%valid
    if(ok)then;message='';else;message='DG symmetry covariance or reconstructed real-space residual gate failed';endif
#else
    result=s_dg_hybrid_acceptance_result();ok=.false.;message='DG acceptance requires MPI'
#endif
  end subroutine evaluate_dg_hybrid_acceptance

  real(real64) function normalized_defect(value,reference) result(defect)
    complex(real64),intent(in)::value(:,:),reference(:,:)
    defect=matrix_norm(value-reference)/max(1d0,matrix_norm(reference))
  end function normalized_defect

  real(real64) function matrix_norm(value) result(norm)
    complex(real64),intent(in)::value(:,:)
    norm=sqrt(sum(abs(value)**2))
  end function matrix_norm

  logical function finite_matrix(value) result(ok)
    complex(real64),intent(in)::value(:,:)
    ok=all(ieee_is_finite(real(value))).and.all(ieee_is_finite(aimag(value)))
  end function finite_matrix

  logical function finite_tensor(value) result(ok)
    complex(real64),intent(in)::value(:,:,:)
    ok=all(ieee_is_finite(real(value))).and.all(ieee_is_finite(aimag(value)))
  end function finite_tensor

  real(real64) function bilinear_defect(value,reference,inverse_metric) result(defect)
    complex(real64),intent(in)::value(:,:),reference(:,:),inverse_metric(:,:)
    defect=bilinear_norm(value-reference,inverse_metric)/max(1d0,bilinear_norm(reference,inverse_metric))
  end function bilinear_defect
  real(real64) function bilinear_norm(value,inverse_metric) result(norm)
    complex(real64),intent(in)::value(:,:),inverse_metric(:,:)
    norm=sqrt(max(0d0,real(sum(conjg(value)*matmul(inverse_metric,matmul(value,inverse_metric))),real64)))
  end function bilinear_norm
  real(real64) function mixed_defect(value,reference,metric,inverse_metric) result(defect)
    complex(real64),intent(in)::value(:,:),reference(:,:),metric(:,:),inverse_metric(:,:)
    defect=mixed_norm(value-reference,metric,inverse_metric)/max(1d0,mixed_norm(reference,metric,inverse_metric))
  end function mixed_defect
  real(real64) function mixed_norm(value,metric,inverse_metric) result(norm)
    complex(real64),intent(in)::value(:,:),metric(:,:),inverse_metric(:,:)
    norm=sqrt(max(0d0,real(sum(conjg(value)*matmul(metric,matmul(value,inverse_metric))),real64)))
  end function mixed_norm
  real(real64) function gamma_defect(value,reference,metric) result(defect)
    complex(real64),intent(in)::value(:,:),reference(:,:),metric(:,:)
    complex(real64)::work(size(value,1),size(value,2))
    work=value-reference
    defect=sqrt(max(0d0,real(sum(conjg(work)*matmul(metric,matmul(work,metric))),real64)))
    defect=defect/max(1d0,sqrt(max(0d0,real(sum(conjg(reference)*matmul(metric,matmul(reference,metric))),real64))))
  end function gamma_defect
  subroutine invert_matrix(value,inverse,ok)
    complex(real64),intent(in)::value(:,:);complex(real64),intent(out)::inverse(:,:);logical,intent(out)::ok
    complex(real64)::work(size(value,1),2*size(value,1)),row(2*size(value,1));integer::i,j,pivot,n
    n=size(value,1);work=(0d0,0d0);work(:,1:n)=value
    do i=1,n;work(i,n+i)=1d0;enddo
    do i=1,n
      pivot=i;do j=i+1,n;if(abs(work(j,i))>abs(work(pivot,i)))pivot=j;enddo
      if(abs(work(pivot,i))<=tiny(1d0))then;ok=.false.;return;endif
      if(pivot/=i)then;row=work(i,:);work(i,:)=work(pivot,:);work(pivot,:)=row;endif
      work(i,:)=work(i,:)/work(i,i)
      do j=1,n;if(j/=i)work(j,:)=work(j,:)-work(j,i)*work(i,:);enddo
    enddo
    inverse=work(:,n+1:2*n);ok=.true.
  end subroutine invert_matrix
  logical function metric_is_spd(value) result(ok)
    complex(real64),intent(in)::value(:,:)
    complex(real64)::factor(size(value,1),size(value,2)),sum_value
    integer::i,j,k,n
    n=size(value,1);ok=.false.
    if(n<1)return
    if(size(value,2)/=n)return
    if(maxval(abs(value-conjg(transpose(value))))>1d-12*max(1d0,matrix_norm(value)))return
    factor=(0d0,0d0)
    do i=1,n
      do j=1,i
        sum_value=value(i,j)
        do k=1,j-1;sum_value=sum_value-factor(i,k)*conjg(factor(j,k));enddo
        if(i==j)then
          if(abs(aimag(sum_value))>1d-12.or.real(sum_value,real64)<=0d0)return
          factor(i,j)=sqrt(real(sum_value,real64))
        else
          factor(i,j)=sum_value/conjg(factor(j,j))
        endif
      enddo
    enddo
    ok=.true.
  end function metric_is_spd
end module dg_hybrid_continuation_acceptance
