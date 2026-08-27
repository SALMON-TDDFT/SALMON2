#include "config.h"
module dg_hybrid_sipg_operator
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public::s_dg_hybrid_sipg_face_operator
    integer::global_face_id=0
    integer::basis_count=0
    complex(real64),allocatable::consistency(:,:)
    complex(real64),allocatable::adjoint_consistency(:,:)
    complex(real64),allocatable::raw_penalty(:,:)
    complex(real64),allocatable::physical_penalty(:,:)
    complex(real64),allocatable::total(:,:)
  end type s_dg_hybrid_sipg_face_operator

  public::assemble_dg_hybrid_sipg_face,scale_dg_hybrid_sipg_faces
contains
  subroutine assemble_dg_hybrid_sipg_face(comm,global_face_id,owner_rank,value_minus,&
      derivative_minus,value_plus,derivative_plus,h_normal,face_weight,penalty_factor,face,ok,message)
    integer,intent(in)::comm,global_face_id,owner_rank
    complex(real64),intent(in)::value_minus(:),derivative_minus(:),value_plus(:),derivative_plus(:)
    real(real64),intent(in)::h_normal,face_weight,penalty_factor
    type(s_dg_hybrid_sipg_face_operator),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,i,j,n,local_bad,global_bad
    complex(real64),allocatable::jump(:),average_derivative(:)

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='SIPG rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='SIPG size query failed';return;endif
    local_bad=merge(0,1,global_face_id>0.and.owner_rank>=0.and.owner_rank<nproc.and.&
      size(value_minus)>0.and.size(value_plus)>0.and.size(derivative_minus)==size(value_minus).and.&
      size(derivative_plus)==size(value_plus).and.h_normal>0d0.and.face_weight>=0d0.and.&
      penalty_factor>0d0.and.ieee_is_finite(h_normal).and.ieee_is_finite(face_weight).and.&
      ieee_is_finite(penalty_factor))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid projected SIPG face';return;endif
    n=size(value_minus)+size(value_plus);face%global_face_id=global_face_id;face%basis_count=n
    allocate(jump(n),average_derivative(n),face%consistency(n,n),face%adjoint_consistency(n,n),&
      face%raw_penalty(n,n),face%physical_penalty(n,n),face%total(n,n))
    jump=[value_minus,-value_plus]
    average_derivative=0.5d0*[derivative_minus,derivative_plus]
    face%consistency=(0d0,0d0);face%adjoint_consistency=(0d0,0d0)
    face%raw_penalty=(0d0,0d0);face%physical_penalty=(0d0,0d0);face%total=(0d0,0d0)
    if(rank==owner_rank)then
      do j=1,n;do i=1,n
        face%consistency(i,j)=-0.5d0*face_weight*conjg(jump(i))*average_derivative(j)
        face%adjoint_consistency(i,j)=-0.5d0*face_weight*conjg(average_derivative(i))*jump(j)
        face%raw_penalty(i,j)=face_weight*(penalty_factor/h_normal)*conjg(jump(i))*jump(j)
        face%physical_penalty(i,j)=0.5d0*face%raw_penalty(i,j)
        face%total(i,j)=face%consistency(i,j)+face%adjoint_consistency(i,j)+face%physical_penalty(i,j)
      enddo;enddo
    endif
    call reduce_matrix(face%consistency,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='SIPG consistency reduction failed';return
    endif
    call reduce_matrix(face%adjoint_consistency,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SIPG adjoint reduction failed';return;endif
    call reduce_matrix(face%raw_penalty,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='SIPG raw penalty reduction failed';return
    endif
    call reduce_matrix(face%physical_penalty,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SIPG physical penalty reduction failed';return;endif
    call reduce_matrix(face%total,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='SIPG total reduction failed';return
    endif
    ok=.true.
#else
    ok=.false.;message='projected SIPG assembly requires MPI'
#endif
  end subroutine assemble_dg_hybrid_sipg_face

  subroutine scale_dg_hybrid_sipg_faces(faces,lambda,face_lambda,scaled,ok,message)
    type(s_dg_hybrid_sipg_face_operator),intent(in)::faces(:)
    real(real64),intent(in)::lambda,face_lambda(:)
    type(s_dg_hybrid_sipg_face_operator),intent(out)::scaled
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,n
    ok=.false.;message=''
    if(size(faces)<1.or.size(face_lambda)/=size(faces).or.lambda<0d0.or.lambda>1d0.or.&
        .not.ieee_is_finite(lambda).or.any(abs(face_lambda-lambda)>1d-15))then
      message='DG interface lambda must be one uniform scalar';return
    endif
    n=faces(1)%basis_count
    if(n<1)then;message='empty SIPG operator';return;endif
    scaled%basis_count=n;scaled%global_face_id=0
    allocate(scaled%consistency(n,n),scaled%adjoint_consistency(n,n),scaled%raw_penalty(n,n),&
      scaled%physical_penalty(n,n),scaled%total(n,n))
    scaled%consistency=(0d0,0d0);scaled%adjoint_consistency=(0d0,0d0)
    scaled%raw_penalty=(0d0,0d0);scaled%physical_penalty=(0d0,0d0);scaled%total=(0d0,0d0)
    do i=1,size(faces)
      if(faces(i)%basis_count/=n)then;message='inconsistent SIPG face basis extent';return;endif
      scaled%consistency=scaled%consistency+lambda*faces(i)%consistency
      scaled%adjoint_consistency=scaled%adjoint_consistency+lambda*faces(i)%adjoint_consistency
      scaled%raw_penalty=scaled%raw_penalty+lambda*faces(i)%raw_penalty
      scaled%physical_penalty=scaled%physical_penalty+lambda*faces(i)%physical_penalty
      scaled%total=scaled%total+lambda*faces(i)%total
    enddo
    ok=.true.
  end subroutine scale_dg_hybrid_sipg_faces

#ifdef USE_MPI
  subroutine reduce_matrix(matrix,comm,ierr)
    complex(real64),intent(inout)::matrix(:,:)
    integer,intent(in)::comm
    integer,intent(out)::ierr
    call MPI_Allreduce(MPI_IN_PLACE,matrix,size(matrix),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine reduce_matrix
#endif
end module dg_hybrid_sipg_operator
