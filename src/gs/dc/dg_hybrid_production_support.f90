#include "config.h"
module dg_hybrid_production_support
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use mpi
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_production_face_traces,only:s_dg_hybrid_production_face_trace,&
    validate_dg_hybrid_production_face_trace_local
  use dg_hybrid_fragment_admission,only:s_dg_hybrid_support_operator,prepare_dg_hybrid_support_operator
  implicit none
  private

  type,public::s_dg_hybrid_production_support_receipt
    logical::valid=.false.
    integer::fragment_id=0,generation=0
    integer::boundary_rows=0,derivative_rows=0,projector_rows=0
    real(real64)::maximum_face_value_defect=huge(1d0)
    real(real64)::maximum_face_derivative_defect=huge(1d0)
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_production_support_receipt

  public::prepare_dg_hybrid_production_support
contains
  subroutine prepare_dg_hybrid_production_support(comm,fragment_id,generation,global_size,coef_nab,basis,&
      faces,expected_face_count,expected_face_fingerprint,projector_offsets,projector_grid_ids,&
      projector_values,projector_weights,operators,&
      operator_fingerprints,receipt,ok,message)
    integer,intent(in)::comm,fragment_id,generation,global_size(3),expected_face_count,projector_offsets(:)
    integer(int64),intent(in)::expected_face_fingerprint
    real(real64),intent(in)::coef_nab(:,:),projector_weights(:)
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    type(s_dg_hybrid_production_face_trace),target,intent(in)::faces(:)
    integer(int64),intent(in)::projector_grid_ids(:)
    complex(real64),intent(in)::projector_values(:)
    type(s_dg_hybrid_support_operator),intent(out)::operators(3)
    integer(int64),intent(out)::operator_fingerprints(3)
    type(s_dg_hybrid_production_support_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_support_operator)::work(3)
    integer(int64),allocatable::boundary_required(:),boundary_samples(:),boundary_points(:),&
      derivative_required(:),derivative_samples(:),derivative_points(:),projector_required(:),projector_samples(:)
    integer,allocatable::boundary_row_offsets(:),derivative_row_offsets(:)
    complex(real64),allocatable::boundary_coefficients(:),derivative_coefficients(:)
    real(real64),allocatable::boundary_weights(:),derivative_weights(:)
    integer::rank,nproc,ierr,local_bad,global_bad,face,row,nface,axis,normal_sign,offset,k,column,slot,&
      point_position(3),plus_position(3),minus_position(3),nz,nprojector
    integer,allocatable::owners(:)
    integer(int64)::work_fingerprints(3),minimum_face_fingerprint,maximum_face_fingerprint,grid_count,&
      recomputed_face_fingerprint
    integer::minimum_face_count,maximum_face_count
    integer(int64),pointer::side_points(:)
    integer,pointer::side_basis_ids(:)
    complex(real64),pointer::side_values(:,:),side_derivatives(:,:)
    complex(real64)::sample
    real(real64)::value_defect,derivative_defect,global_value_defect,global_derivative_defect
    logical::prepared,face_valid
    character(512)::why

    ok=.false.;message='';operator_fingerprints=0_int64;work_fingerprints=0_int64
    receipt=s_dg_hybrid_production_support_receipt()
    minimum_face_count=0;maximum_face_count=-1
    minimum_face_fingerprint=0_int64;maximum_face_fingerprint=1_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=0;grid_count=0_int64
    if(fragment_id<1.or.fragment_id>nproc.or.generation<1.or.any(global_size<1).or.&
        size(coef_nab,1)<1.or.size(coef_nab,2)/=3.or..not.all(ieee_is_finite(coef_nab)).or.&
        basis%fragment_id/=fragment_id.or.basis%generation/=generation.or.&
        .not.allocated(basis%global_ids).or..not.allocated(basis%buffer_point_ids).or.&
        .not.allocated(basis%buffer_values))local_bad=1
    if(local_bad==0)then
      call checked_grid_count(global_size,grid_count,prepared)
      if(.not.prepared)local_bad=1
    endif
    if(local_bad==0)then
      if(any(shape(basis%buffer_values)/=[size(basis%buffer_point_ids),size(basis%global_ids)]).or.&
          any(basis%buffer_point_ids<1_int64).or.&
          any(basis%buffer_point_ids>grid_count))local_bad=1
    endif
    nprojector=max(0,size(projector_offsets)-1)
    if(size(projector_weights)/=nprojector.or.size(projector_grid_ids)/=size(projector_values).or.&
        .not.all(ieee_is_finite(projector_weights)).or.any(projector_weights<=0d0).or.&
        .not.all(ieee_is_finite(real(projector_values))).or.&
        .not.all(ieee_is_finite(aimag(projector_values))))local_bad=1
    if(size(projector_offsets)<1)local_bad=1
    if(local_bad==0)then
      if(projector_offsets(1)/=1.or.projector_offsets(nprojector+1)/=size(projector_values)+1.or.&
          any(projector_offsets(2:)<=projector_offsets(:nprojector)).or.any(projector_grid_ids<1_int64).or.&
          any(projector_grid_ids>grid_count))local_bad=1
    endif
    if(expected_face_count<1.or.expected_face_fingerprint==0_int64.or.size(faces)/=expected_face_count)local_bad=1
    if(local_bad==0)then
      recomputed_face_fingerprint=int(z'6A09E667F3BCC909',int64)
      recomputed_face_fingerprint=ieor(ishftc(recomputed_face_fingerprint,11),int(size(faces),int64))
      do face=1,size(faces)
        if(faces(face)%production_inventory_count/=expected_face_count.or.&
            faces(face)%production_inventory_fingerprint/=expected_face_fingerprint.or.&
            faces(face)%production_inventory_slot_fingerprint==0_int64.or.faces(face)%global_face_id/=face)&
          local_bad=1
        do k=1,face-1
          if(faces(k)%production_inventory_slot_fingerprint==faces(face)%production_inventory_slot_fingerprint)&
            local_bad=1
        enddo
        recomputed_face_fingerprint=ieor(ishftc(recomputed_face_fingerprint,11),&
          faces(face)%production_inventory_slot_fingerprint)
      enddo
      if(recomputed_face_fingerprint/=expected_face_fingerprint)local_bad=1
    endif
    allocate(owners(nproc));owners=0
    if(fragment_id>=1.and.fragment_id<=nproc)owners(fragment_id)=1
    call MPI_Allreduce(MPI_IN_PLACE,owners,nproc,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(owners/=1))local_bad=1
    call MPI_Allreduce(expected_face_count,minimum_face_count,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(expected_face_count,maximum_face_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(expected_face_fingerprint,minimum_face_fingerprint,1,&
      MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(expected_face_fingerprint,maximum_face_fingerprint,1,&
      MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_face_count/=maximum_face_count.or.&
        minimum_face_fingerprint/=maximum_face_fingerprint)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid production support or rank-fragment contract';return
    endif

    nface=0
    do face=1,size(faces)
      if(faces(face)%minus_fragment/=fragment_id.and.faces(face)%plus_fragment/=fragment_id)then
        if(faces(face)%frozen)local_bad=1
        cycle
      endif
      if(.not.faces(face)%frozen)then;local_bad=1;cycle;endif
      call validate_dg_hybrid_production_face_trace_local(faces(face),face_valid)
      if(.not.face_valid)then;local_bad=1;cycle;endif
      if(size(faces(face)%weights)>huge(nface)-nface)then;local_bad=1;cycle;endif
      nface=nface+size(faces(face)%weights)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid frozen production face support inventory';return;endif
    if(nface>0)then
      if(nface>huge(nface)/2)then
        local_bad=1
      else if(size(coef_nab,1)>huge(nface)/(2*nface))then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='production derivative support size overflow';return;endif
    allocate(boundary_required(nface),boundary_samples(nface),boundary_points(nface),&
      boundary_row_offsets(nface+1),boundary_coefficients(nface),boundary_weights(nface),&
      derivative_required(nface),derivative_samples(nface),derivative_row_offsets(nface+1),&
      derivative_points(2*size(coef_nab,1)*nface),&
      derivative_coefficients(2*size(coef_nab,1)*nface),derivative_weights(nface))
    row=0;nz=0;boundary_row_offsets(1)=1;derivative_row_offsets(1)=1
    value_defect=0d0;derivative_defect=0d0;local_bad=0
    do face=1,size(faces)
      if(.not.faces(face)%frozen)cycle
      if(faces(face)%minus_fragment==fragment_id)then
        side_points=>faces(face)%point_ids_minus;side_basis_ids=>faces(face)%basis_ids_minus
        side_values=>faces(face)%value_minus;side_derivatives=>faces(face)%derivative_minus
      else if(faces(face)%plus_fragment==fragment_id)then
        side_points=>faces(face)%point_ids_plus;side_basis_ids=>faces(face)%basis_ids_plus
        side_values=>faces(face)%value_plus;side_derivatives=>faces(face)%derivative_plus
      else
        cycle
      endif
      axis=maxloc(abs(faces(face)%canonical_normal),dim=1)
      normal_sign=nint(faces(face)%canonical_normal(axis))
      if(abs(faces(face)%canonical_normal(axis)-real(normal_sign,real64))>1d-12.or.&
          count(abs(faces(face)%canonical_normal)>1d-12)/=1.or.&
          size(side_points)/=size(faces(face)%weights).or.size(side_basis_ids)/=size(basis%global_ids).or.&
          any(side_basis_ids/=int(basis%global_ids)))then;local_bad=1;cycle;endif
      do k=1,size(side_points)
        row=row+1
        boundary_required(row)=int(row,int64);boundary_samples(row)=boundary_required(row)
        boundary_points(row)=side_points(k);boundary_coefficients(row)=(1d0,0d0)
        boundary_weights(row)=faces(face)%weights(k);boundary_row_offsets(row+1)=row+1
        derivative_required(row)=int(row,int64);derivative_samples(row)=derivative_required(row)
        derivative_weights(row)=faces(face)%weights(k)
        call grid_position(side_points(k),global_size,point_position)
        do offset=1,size(coef_nab,1)
          plus_position=point_position;minus_position=point_position
          plus_position(axis)=modulo(point_position(axis)+offset,global_size(axis))
          minus_position(axis)=modulo(point_position(axis)-offset,global_size(axis))
          nz=nz+1;derivative_points(nz)=physical_grid_id(plus_position,global_size)
          derivative_coefficients(nz)=real(normal_sign,real64)*coef_nab(offset,axis)
          nz=nz+1;derivative_points(nz)=physical_grid_id(minus_position,global_size)
          derivative_coefficients(nz)=-real(normal_sign,real64)*coef_nab(offset,axis)
        enddo
        derivative_row_offsets(row+1)=nz+1
        do column=1,size(basis%global_ids)
          slot=findloc(basis%buffer_point_ids,side_points(k),dim=1)
          if(slot<=0)then;local_bad=1;cycle;endif
          value_defect=max(value_defect,abs(basis%buffer_values(slot,column)-side_values(k,column)))
          sample=(0d0,0d0)
          do offset=derivative_row_offsets(row),derivative_row_offsets(row+1)-1
            slot=findloc(basis%buffer_point_ids,derivative_points(offset),dim=1)
            if(slot<=0)then;local_bad=1;cycle;endif
            sample=sample+derivative_coefficients(offset)*basis%buffer_values(slot,column)
          enddo
          derivative_defect=max(derivative_defect,abs(sample-side_derivatives(k,column)))
        enddo
      enddo
    enddo
    if(row/=nface.or.nz/=size(derivative_points))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(value_defect,global_value_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(derivative_defect,global_derivative_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid production face support inventory';return;endif
    if(global_value_defect>1d-12)then;message='production boundary manifest does not reproduce frozen trace';return;endif
    if(global_derivative_defect>1d-12)then;message='production derivative manifest does not reproduce frozen trace';return;endif

    allocate(projector_required(nprojector),projector_samples(nprojector))
    projector_required=[(int(k,int64),k=1,nprojector)];projector_samples=projector_required
    call prepare_dg_hybrid_support_operator(comm,fragment_id,generation,1,boundary_required,boundary_samples,&
      boundary_row_offsets,boundary_points,boundary_coefficients,boundary_weights,work(1),&
      work_fingerprints(1),prepared,why)
    if(.not.prepared)then;message='production boundary support: '//trim(why);return;endif
    call prepare_dg_hybrid_support_operator(comm,fragment_id,generation,2,derivative_required,derivative_samples,&
      derivative_row_offsets,derivative_points,derivative_coefficients,derivative_weights,work(2),&
      work_fingerprints(2),prepared,why)
    if(.not.prepared)then;message='production derivative support: '//trim(why);return;endif
    call prepare_dg_hybrid_support_operator(comm,fragment_id,generation,3,projector_required,projector_samples,&
      projector_offsets,projector_grid_ids,projector_values,projector_weights,work(3),&
      work_fingerprints(3),prepared,why)
    if(.not.prepared)then;message='production projector support: '//trim(why);return;endif
    operators=work;operator_fingerprints=work_fingerprints
    receipt%fragment_id=fragment_id;receipt%generation=generation
    receipt%boundary_rows=nface;receipt%derivative_rows=nface;receipt%projector_rows=nprojector
    receipt%maximum_face_value_defect=global_value_defect
    receipt%maximum_face_derivative_defect=global_derivative_defect
    receipt%fingerprint=work_fingerprints(1)
    receipt%fingerprint=ieor(ishftc(receipt%fingerprint,13),work_fingerprints(2))
    receipt%fingerprint=ieor(ishftc(receipt%fingerprint,13),work_fingerprints(3))
    if(receipt%fingerprint==0_int64)receipt%fingerprint=1_int64
    receipt%valid=.true.;ok=.true.;message=''
  end subroutine prepare_dg_hybrid_production_support

  pure integer(int64) function physical_grid_id(position,global_size) result(identifier)
    integer,intent(in)::position(3),global_size(3)
    identifier=1_int64+int(position(1),int64)+int(global_size(1),int64)*&
      (int(position(2),int64)+int(global_size(2),int64)*int(position(3),int64))
  end function physical_grid_id

  pure subroutine grid_position(identifier,global_size,position)
    integer(int64),intent(in)::identifier
    integer,intent(in)::global_size(3)
    integer,intent(out)::position(3)
    integer(int64)::offset,plane
    offset=identifier-1_int64;plane=int(global_size(1),int64)*int(global_size(2),int64)
    position(3)=int(offset/plane);offset=modulo(offset,plane)
    position(2)=int(offset/int(global_size(1),int64));position(1)=int(modulo(offset,int(global_size(1),int64)))
  end subroutine grid_position

  pure subroutine checked_grid_count(global_size,count,ok)
    integer,intent(in)::global_size(3)
    integer(int64),intent(out)::count
    logical,intent(out)::ok
    integer::axis
    count=1_int64;ok=all(global_size>0)
    if(.not.ok)return
    do axis=1,3
      if(int(global_size(axis),int64)>huge(count)/count)then;count=0_int64;ok=.false.;return;endif
      count=count*int(global_size(axis),int64)
    enddo
  end subroutine checked_grid_count
end module dg_hybrid_production_support
