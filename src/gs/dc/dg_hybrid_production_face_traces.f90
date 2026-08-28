#include "config.h"
module dg_hybrid_production_face_traces
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sipg_operator,only:s_dg_hybrid_sipg_face_operator,assemble_dg_hybrid_sipg_face
#ifdef USE_MPI
  use mpi,only:MPI_Allreduce,MPI_Comm_size,MPI_INTEGER,MPI_INTEGER8,MPI_MAX,MPI_MIN,MPI_SUCCESS
#endif
  implicit none
  private

  type,public::s_dg_hybrid_production_face_trace
    logical::frozen=.false.
    integer::global_face_id=0,owner_rank=-1,minus_fragment=0,plus_fragment=0
    integer::periodic_shift(3)=0
    real(real64)::canonical_normal(3)=0d0,h_normal=0d0
    integer(int64),allocatable::point_ids(:)
    real(real64),allocatable::weights(:)
    integer,allocatable::basis_ids_minus(:),basis_ids_plus(:)
    complex(real64),allocatable::value_minus(:,:),value_plus(:,:)
    complex(real64),allocatable::derivative_minus(:,:),derivative_plus(:,:)
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_production_face_trace

  public::build_dg_hybrid_production_face_trace,assemble_dg_hybrid_production_face
contains
  subroutine build_dg_hybrid_production_face_trace(icomm,face_id,owner,fragment_minus,fragment_plus,periodic_shift,&
      normal,h_normal,point_ids_minus,point_ids_plus,weights,basis_ids_minus,basis_ids_plus,value_minus,&
      outward_minus,value_plus,outward_plus,effective_ids,group_action,face,ok,message)
    integer,intent(in)::icomm,face_id,owner,fragment_minus,fragment_plus,periodic_shift(3)
    integer(int64),intent(in)::point_ids_minus(:),point_ids_plus(:)
    integer,intent(in)::basis_ids_minus(:),basis_ids_plus(:),effective_ids(:),group_action(:,:)
    real(real64),intent(in)::normal(3),h_normal,weights(:)
    complex(real64),intent(in)::value_minus(:,:),outward_minus(:,:),value_plus(:,:),outward_plus(:,:)
    type(s_dg_hybrid_production_face_trace),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,nproc,ierr,local_bad,global_bad
    integer(int64)::local_hash,minimum_hash,maximum_hash

    ok=.false.;message=''
    call MPI_Comm_size(icomm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production face communicator query failed';return;endif
    local_bad=merge(0,1,face_id>0.and.owner>=0.and.owner<nproc.and.fragment_minus>0.and.fragment_plus>0.and.&
      fragment_minus/=fragment_plus.and.size(point_ids_minus)>0.and.size(point_ids_minus)==size(point_ids_plus).and.&
      size(weights)==size(point_ids_minus).and.all(point_ids_minus==point_ids_plus).and.h_normal>0d0.and.&
      abs(sqrt(sum(normal**2))-1d0)<1d-12.and.all(ieee_is_finite(normal)).and.ieee_is_finite(h_normal).and.&
      all(ieee_is_finite(weights)).and.all(weights>0d0).and.&
      size(value_minus,1)==size(weights).and.size(value_minus,2)==size(basis_ids_minus).and.&
      all(shape(outward_minus)==shape(value_minus)).and.size(value_plus,1)==size(weights).and.&
      size(value_plus,2)==size(basis_ids_plus).and.all(shape(outward_plus)==shape(value_plus)).and.&
      all(ieee_is_finite(real(value_minus))).and.all(ieee_is_finite(aimag(value_minus))).and.&
      all(ieee_is_finite(real(outward_minus))).and.all(ieee_is_finite(aimag(outward_minus))).and.&
      all(ieee_is_finite(real(value_plus))).and.all(ieee_is_finite(aimag(value_plus))).and.&
      all(ieee_is_finite(real(outward_plus))).and.all(ieee_is_finite(aimag(outward_plus))).and.&
      valid_closed_action(effective_ids,group_action).and.&
      all([(count(effective_ids==basis_ids_minus(i))==1,i=1,size(basis_ids_minus))]).and.&
      all([(count(effective_ids==basis_ids_plus(i))==1,i=1,size(basis_ids_plus))]))
    if(local_bad==0)then
      do i=1,size(basis_ids_minus)
        if(count(basis_ids_minus==basis_ids_minus(i))/=1.or.any(basis_ids_plus==basis_ids_minus(i)))local_bad=1
      enddo
      do i=1,size(basis_ids_plus)
        if(count(basis_ids_plus==basis_ids_plus(i))/=1)local_bad=1
      enddo
      do i=1,size(point_ids_minus)
        if(point_ids_minus(i)<=0_int64.or.count(point_ids_minus==point_ids_minus(i))/=1)local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid or incomplete production face trace';return;endif
    local_hash=trace_fingerprint(face_id,owner,fragment_minus,fragment_plus,periodic_shift,normal,h_normal,&
      point_ids_minus,weights,basis_ids_minus,basis_ids_plus,value_minus,outward_minus,value_plus,outward_plus,&
      effective_ids,group_action)
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing production face topology or trace';return
    endif
    face%global_face_id=face_id;face%owner_rank=owner
    face%minus_fragment=fragment_minus;face%plus_fragment=fragment_plus
    face%periodic_shift=periodic_shift;face%canonical_normal=normal;face%h_normal=h_normal
    allocate(face%point_ids(size(point_ids_minus)),face%weights(size(weights)),&
      face%basis_ids_minus(size(basis_ids_minus)),face%basis_ids_plus(size(basis_ids_plus)),&
      face%value_minus(size(value_minus,1),size(value_minus,2)),&
      face%derivative_minus(size(outward_minus,1),size(outward_minus,2)),&
      face%value_plus(size(value_plus,1),size(value_plus,2)),&
      face%derivative_plus(size(outward_plus,1),size(outward_plus,2)))
    face%point_ids=point_ids_minus;face%weights=weights
    face%basis_ids_minus=basis_ids_minus;face%basis_ids_plus=basis_ids_plus
    face%value_minus=value_minus;face%derivative_minus=outward_minus
    face%value_plus=value_plus;face%derivative_plus=-outward_plus
    face%fingerprint=local_hash;if(face%fingerprint==0_int64)face%fingerprint=1_int64
    face%frozen=.true.;ok=.true.
#else
    ok=.false.;message='production face traces require MPI'
#endif
  end subroutine build_dg_hybrid_production_face_trace

  subroutine assemble_dg_hybrid_production_face(icomm,trace,penalty_factor,face,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_production_face_trace),intent(in)::trace
    real(real64),intent(in)::penalty_factor
    type(s_dg_hybrid_sipg_face_operator),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_sipg_face_operator)::point_face
    integer::point
    ok=.false.;message=''
    if(.not.trace%frozen.or.trace%fingerprint==0_int64.or..not.allocated(trace%weights))then
      message='production face trace is not frozen';return
    endif
    do point=1,size(trace%weights)
      call assemble_dg_hybrid_sipg_face(icomm,trace%global_face_id,trace%owner_rank,trace%periodic_shift,&
        trace%basis_ids_minus,trace%basis_ids_plus,trace%value_minus(point,:),trace%derivative_minus(point,:),&
        trace%value_plus(point,:),trace%derivative_plus(point,:),trace%h_normal,trace%weights(point),&
        penalty_factor,point_face,ok,message)
      if(.not.ok)return
      if(point==1)then
        face=point_face
      else
        face%consistency=face%consistency+point_face%consistency
        face%adjoint_consistency=face%adjoint_consistency+point_face%adjoint_consistency
        face%raw_penalty=face%raw_penalty+point_face%raw_penalty
        face%physical_penalty=face%physical_penalty+point_face%physical_penalty
        face%total=face%total+point_face%total
      endif
    enddo
    ok=.true.
  end subroutine assemble_dg_hybrid_production_face

  logical function valid_closed_action(ids,action) result(valid)
    integer,intent(in)::ids(:),action(:,:)
    integer::i,j
    valid=size(ids)>0.and.size(action,1)==size(ids).and.size(action,2)>0
    if(.not.valid)return
    valid=all(ids>0).and.all(action>=1).and.all(action<=size(ids)).and.&
      all(action(:,1)==[(i,i=1,size(ids))])
    if(.not.valid)return
    do i=1,size(ids)
      if(count(ids==ids(i))/=1)then;valid=.false.;return;endif
    enddo
    do j=1,size(action,2)
      do i=1,size(ids)
        if(count(action(:,j)==i)/=1)then;valid=.false.;return;endif
      enddo
    enddo
  end function valid_closed_action

  integer(int64) function trace_fingerprint(face_id,owner,fragment_minus,fragment_plus,shift,normal,h_normal,&
      point_ids,weights,ids_minus,ids_plus,value_minus,outward_minus,value_plus,outward_plus,effective_ids,action)&
      result(hash)
    integer,intent(in)::face_id,owner,fragment_minus,fragment_plus,shift(3),ids_minus(:),ids_plus(:),&
      effective_ids(:),action(:,:)
    integer(int64),intent(in)::point_ids(:)
    real(real64),intent(in)::normal(3),h_normal,weights(:)
    complex(real64),intent(in)::value_minus(:,:),outward_minus(:,:),value_plus(:,:),outward_plus(:,:)
    integer::i,j
    hash=int(z'510E527FADE682D1',int64)
    call mix(int(face_id,int64));call mix(int(owner,int64));call mix(int(fragment_minus,int64))
    call mix(int(fragment_plus,int64))
    do i=1,3;call mix(int(shift(i),int64));call mix(transfer(normal(i),hash));enddo
    call mix(transfer(h_normal,hash))
    do i=1,size(point_ids);call mix(point_ids(i));call mix(transfer(weights(i),hash));enddo
    do i=1,size(ids_minus);call mix(int(ids_minus(i),int64));enddo
    do i=1,size(ids_plus);call mix(int(ids_plus(i),int64));enddo
    do j=1,size(value_minus,2);do i=1,size(value_minus,1)
      call mix_complex(value_minus(i,j));call mix_complex(outward_minus(i,j))
    enddo;enddo
    do j=1,size(value_plus,2);do i=1,size(value_plus,1)
      call mix_complex(value_plus(i,j));call mix_complex(outward_plus(i,j))
    enddo;enddo
    do i=1,size(effective_ids);call mix(int(effective_ids(i),int64));enddo
    do j=1,size(action,2);do i=1,size(action,1);call mix(int(action(i,j),int64));enddo;enddo
  contains
    subroutine mix(value)
      integer(int64),intent(in)::value
      hash=ieor(ishftc(hash,11),value);hash=ieor(hash,ishftc(hash,17))
    end subroutine mix
    subroutine mix_complex(value)
      complex(real64),intent(in)::value
      call mix(transfer(real(value,real64),hash));call mix(transfer(aimag(value),hash))
    end subroutine mix_complex
  end function trace_fingerprint
end module dg_hybrid_production_face_traces
