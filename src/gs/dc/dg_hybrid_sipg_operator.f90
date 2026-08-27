#include "config.h"
module dg_hybrid_sipg_operator
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public::s_dg_hybrid_sipg_face_operator
    integer::global_face_id=0
    integer::periodic_shift(3)=0
    integer::basis_count=0
    integer,allocatable::global_basis_ids(:)
    complex(real64),allocatable::consistency(:,:)
    complex(real64),allocatable::adjoint_consistency(:,:)
    complex(real64),allocatable::raw_penalty(:,:)
    complex(real64),allocatable::physical_penalty(:,:)
    complex(real64),allocatable::total(:,:)
  end type s_dg_hybrid_sipg_face_operator

  public::assemble_dg_hybrid_sipg_face,scale_dg_hybrid_sipg_faces
contains
  subroutine assemble_dg_hybrid_sipg_face(comm,global_face_id,owner_rank,periodic_shift,basis_ids_minus,basis_ids_plus,value_minus,&
      derivative_minus,value_plus,derivative_plus,h_normal,face_weight,penalty_factor,face,ok,message)
    integer,intent(in)::comm,global_face_id,owner_rank,periodic_shift(3)
    integer,intent(in)::basis_ids_minus(:),basis_ids_plus(:)
    complex(real64),intent(in)::value_minus(:),derivative_minus(:),value_plus(:),derivative_plus(:)
    real(real64),intent(in)::h_normal,face_weight,penalty_factor
    type(s_dg_hybrid_sipg_face_operator),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,i,j,n,local_bad,global_bad,min_n,max_n
    integer(int64)::local_hash,min_hash,max_hash
    complex(real64),allocatable::jump(:),average_derivative(:)

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='SIPG rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='SIPG size query failed';return;endif
    local_bad=merge(0,1,global_face_id>0.and.owner_rank>=0.and.owner_rank<nproc.and.&
      size(value_minus)>0.and.size(value_plus)>0.and.size(basis_ids_minus)==size(value_minus).and.&
      size(basis_ids_plus)==size(value_plus).and.size(derivative_minus)==size(value_minus).and.&
      size(derivative_plus)==size(value_plus).and.h_normal>0d0.and.face_weight>=0d0.and.&
      penalty_factor>0d0.and.ieee_is_finite(h_normal).and.ieee_is_finite(face_weight).and.&
      ieee_is_finite(penalty_factor).and.all(basis_ids_minus>0).and.all(basis_ids_plus>0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid projected SIPG face';return;endif
    n=size(value_minus)+size(value_plus)
    call MPI_Allreduce(n,min_n,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(n,max_n,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(min_n/=max_n)then;message='rank-disagreeing SIPG trace extent';return;endif
    local_hash=face_input_fingerprint(global_face_id,owner_rank,periodic_shift,basis_ids_minus,basis_ids_plus,value_minus,&
      derivative_minus,value_plus,derivative_plus,h_normal,face_weight,penalty_factor)
    call MPI_Allreduce(local_hash,min_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_hash,max_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(min_hash/=max_hash)then;message='rank-disagreeing SIPG face metadata or traces';return;endif
    face%global_face_id=global_face_id;face%periodic_shift=periodic_shift;face%basis_count=n
    allocate(face%global_basis_ids(n));face%global_basis_ids=[basis_ids_minus,basis_ids_plus]
    do i=1,n;do j=i+1,n
      if(face%global_basis_ids(i)==face%global_basis_ids(j))then;message='duplicate global basis ID on SIPG face';return;endif
    enddo;enddo
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
    integer::i,j,k,n
    integer,allocatable::all_ids(:),positions(:)
    ok=.false.;message=''
    if(size(faces)<1.or.size(face_lambda)/=size(faces).or.lambda<0d0.or.lambda>1d0.or.&
        .not.ieee_is_finite(lambda).or.any(abs(face_lambda-lambda)>1d-15))then
      message='DG interface lambda must be one uniform scalar';return
    endif
    n=sum([(faces(i)%basis_count,i=1,size(faces))])
    if(n<1)then;message='empty SIPG operator';return;endif
    allocate(all_ids(n));k=0
    do i=1,size(faces)
      if(.not.allocated(faces(i)%global_basis_ids).or.size(faces(i)%global_basis_ids)/=faces(i)%basis_count)then
        message='SIPG face lacks global basis IDs';return
      endif
      all_ids(k+1:k+faces(i)%basis_count)=faces(i)%global_basis_ids;k=k+faces(i)%basis_count
    enddo
    call sort_unique(all_ids,n)
    scaled%basis_count=n;scaled%global_face_id=0
    allocate(scaled%global_basis_ids(n));scaled%global_basis_ids=all_ids(:n)
    allocate(scaled%consistency(n,n),scaled%adjoint_consistency(n,n),scaled%raw_penalty(n,n),&
      scaled%physical_penalty(n,n),scaled%total(n,n))
    scaled%consistency=(0d0,0d0);scaled%adjoint_consistency=(0d0,0d0)
    scaled%raw_penalty=(0d0,0d0);scaled%physical_penalty=(0d0,0d0);scaled%total=(0d0,0d0)
    do i=1,size(faces)
      allocate(positions(faces(i)%basis_count))
      do j=1,faces(i)%basis_count
        positions(j)=find_id(faces(i)%global_basis_ids(j),scaled%global_basis_ids)
      enddo
      do j=1,faces(i)%basis_count;do k=1,faces(i)%basis_count
        scaled%consistency(positions(k),positions(j))=scaled%consistency(positions(k),positions(j))+lambda*faces(i)%consistency(k,j)
        scaled%adjoint_consistency(positions(k),positions(j))=scaled%adjoint_consistency(positions(k),positions(j))+&
          lambda*faces(i)%adjoint_consistency(k,j)
        scaled%raw_penalty(positions(k),positions(j))=scaled%raw_penalty(positions(k),positions(j))+lambda*faces(i)%raw_penalty(k,j)
        scaled%physical_penalty(positions(k),positions(j))=scaled%physical_penalty(positions(k),positions(j))+&
          lambda*faces(i)%physical_penalty(k,j)
        scaled%total(positions(k),positions(j))=scaled%total(positions(k),positions(j))+lambda*faces(i)%total(k,j)
      enddo;enddo
      deallocate(positions)
    enddo
    ok=.true.
  end subroutine scale_dg_hybrid_sipg_faces

  subroutine sort_unique(values,count_unique)
    integer,intent(inout)::values(:)
    integer,intent(out)::count_unique
    integer::i,j,key
    do i=2,size(values)
      key=values(i);j=i-1
      do while(j>=1)
        if(values(j)<=key)exit
        values(j+1)=values(j);j=j-1
      enddo
      values(j+1)=key
    enddo
    count_unique=0
    do i=1,size(values)
      if(i==1)then
        count_unique=1;values(1)=values(i)
      else if(values(i)/=values(count_unique))then
        count_unique=count_unique+1;values(count_unique)=values(i)
      endif
    enddo
  end subroutine sort_unique

  integer function find_id(value,values) result(position)
    integer,intent(in)::value,values(:)
    integer::i
    position=0;do i=1,size(values);if(values(i)==value)then;position=i;return;endif;enddo
  end function find_id

#ifdef USE_MPI
  integer(int64) function face_input_fingerprint(face_id,owner,periodic_shift,ids_minus,ids_plus,value_minus,derivative_minus,&
      value_plus,derivative_plus,h_normal,face_weight,penalty_factor) result(hash)
    integer,intent(in)::face_id,owner,periodic_shift(3),ids_minus(:),ids_plus(:)
    complex(real64),intent(in)::value_minus(:),derivative_minus(:),value_plus(:),derivative_plus(:)
    real(real64),intent(in)::h_normal,face_weight,penalty_factor
    integer::i
    hash=int(z'6A09E667F3BCC909',int64);call mix(int(face_id,int64));call mix(int(owner,int64))
    do i=1,3;call mix(int(periodic_shift(i),int64));enddo
    do i=1,size(ids_minus)
      call mix(int(ids_minus(i),int64));call mix_complex(value_minus(i));call mix_complex(derivative_minus(i))
    enddo
    do i=1,size(ids_plus)
      call mix(int(ids_plus(i),int64));call mix_complex(value_plus(i));call mix_complex(derivative_plus(i))
    enddo
    call mix(transfer(h_normal,hash));call mix(transfer(face_weight,hash));call mix(transfer(penalty_factor,hash))
  contains
    subroutine mix(value)
      integer(int64),intent(in)::value
      hash=ieor(ishftc(hash,11),value)
    end subroutine mix
    subroutine mix_complex(value)
      complex(real64),intent(in)::value
      call mix(transfer(real(value,real64),hash));call mix(transfer(aimag(value),hash))
    end subroutine mix_complex
  end function face_input_fingerprint

  subroutine reduce_matrix(matrix,comm,ierr)
    complex(real64),intent(inout)::matrix(:,:)
    integer,intent(in)::comm
    integer,intent(out)::ierr
    call MPI_Allreduce(MPI_IN_PLACE,matrix,size(matrix),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine reduce_matrix
#endif
end module dg_hybrid_sipg_operator
