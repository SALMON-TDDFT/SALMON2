#include "config.h"
module dg_hybrid_divided_operator
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,s_dg_hybrid_variational_iterate,&
    compose_dg_hybrid_variational_hamiltonian,verify_dg_hybrid_variational_payload_rows
  use dg_hybrid_wannier_complement,only:compute_dg_hybrid_union_to_complete_binding
  implicit none
  private
  public::extract_dg_hybrid_fragment_self_block,compose_dg_hybrid_complete_rows
  public::dg_hybrid_fragment_directory_fingerprint
  interface
    subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
      import::real64
      character(1),intent(in)::jobz,uplo
      integer,intent(in)::n,lda,lwork
      complex(real64),intent(inout)::a(lda,*),work(*)
      real(real64),intent(out)::w(*),rwork(*)
      integer,intent(out)::info
    end subroutine zheev
  end interface
contains
  subroutine extract_dg_hybrid_fragment_self_block(comm_fragment,fragment_id,row_ids,basis_fragment,&
      fixed_payload,local_rows,hff,sff,ok,message,basis_local_slot,basis_generation,&
      fragment_catalog_fingerprint,fragment_directory_fingerprint)
    integer,intent(in)::comm_fragment,fragment_id,basis_fragment(:)
    integer(int64),intent(in)::row_ids(:)
    type(s_dg_hybrid_fixed_payload),intent(in)::fixed_payload
    complex(real64),intent(in)::local_rows(:,:)
    complex(real64),allocatable,intent(out)::hff(:,:),sff(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::basis_local_slot(:),basis_generation(:)
    integer(int64),intent(in),optional::fragment_catalog_fingerprint,fragment_directory_fingerprint
    integer::rank,ierr,nunion,nfragment,a,p,source_position,owner,output_position,allocation_status
    integer,allocatable::fragment_slots(:),source_counts(:),output_owners(:)
    complex(real64),allocatable::local_h(:),local_s(:),reduced_h(:),reduced_s(:),&
      working_hff(:,:),working_sff(:,:)
    logical::stage_ok

    ok=.false.;message=''
    call MPI_Comm_rank(comm_fragment,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment operator rank query failed';return;endif
    call validate_payload_contract(comm_fragment,fixed_payload,local_rows,stage_ok,message)
    if(.not.stage_ok)return
    call collective_gate(comm_fragment,present(basis_local_slot).and.present(basis_generation).and.&
      present(fragment_catalog_fingerprint).and.present(fragment_directory_fingerprint),&
      'fragment catalog provenance is required',stage_ok,message)
    if(.not.stage_ok)return
    nunion=fixed_payload%global_basis_count
    call validate_fragment_contract(comm_fragment,fragment_id,basis_fragment,basis_local_slot,&
      basis_generation,fragment_catalog_fingerprint,fragment_directory_fingerprint,&
      fixed_payload%basis_fingerprint,fixed_payload%basis_directory_fingerprint,nunion,stage_ok,message)
    if(.not.stage_ok)return
    nfragment=count(basis_fragment==fragment_id)
    allocate(fragment_slots(nfragment),source_counts(nfragment),stat=allocation_status)
    call collective_gate(comm_fragment,allocation_status==0,'fragment slot allocation failed',stage_ok,message)
    if(.not.stage_ok)return
    fragment_slots=0
    do p=1,nunion
      if(basis_fragment(p)==fragment_id)fragment_slots(basis_local_slot(p))=p
    enddo

    source_counts=0
    do p=1,size(fixed_payload%row_ids)
      source_position=findloc(fragment_slots,int(fixed_payload%row_ids(p)),dim=1)
      if(source_position>0)source_counts(source_position)=source_counts(source_position)+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,source_counts,nfragment,MPI_INTEGER,MPI_SUM,comm_fragment,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment source-row ownership reduction failed';return;endif
    call collective_gate(comm_fragment,all(source_counts==1),&
      'fragment source rows are not owned exactly once',stage_ok,message)
    if(.not.stage_ok)return
    call build_output_directory(comm_fragment,nfragment,row_ids,output_owners,stage_ok,message)
    if(.not.stage_ok)return

    allocate(local_h(nfragment),local_s(nfragment),reduced_h(nfragment),reduced_s(nfragment),&
      working_hff(size(row_ids),nfragment),working_sff(size(row_ids),nfragment),stat=allocation_status)
    call collective_gate(comm_fragment,allocation_status==0,&
      'fragment self-block workspace allocation failed',stage_ok,message)
    if(.not.stage_ok)return
    working_hff=(0d0,0d0);working_sff=(0d0,0d0)
    do a=1,nfragment
      local_h=(0d0,0d0);local_s=(0d0,0d0)
      source_position=findloc(fixed_payload%row_ids,int(fragment_slots(a),int64),dim=1)
      if(source_position>0)then
        local_h=fixed_payload%kinetic_rows(source_position,fragment_slots)+&
          fixed_payload%nonlocal_rows(source_position,fragment_slots)+&
          fixed_payload%interface_rows(source_position,fragment_slots)+&
          local_rows(source_position,fragment_slots)
        local_s=fixed_payload%metric_rows(source_position,fragment_slots)
      endif
      owner=output_owners(a)-1
      call MPI_Reduce(local_h,reduced_h,nfragment,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm_fragment,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment Hamiltonian row transfer failed';return;endif
      call MPI_Reduce(local_s,reduced_s,nfragment,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm_fragment,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment metric row transfer failed';return;endif
      if(rank==owner)then
        output_position=findloc(row_ids,int(a,int64),dim=1)
        if(output_position<1)then;message='fragment output row owner is inconsistent';return;endif
        working_hff(output_position,:)=reduced_h
        working_sff(output_position,:)=reduced_s
      endif
    enddo
    call collective_gate(comm_fragment,finite_matrix(working_hff).and.finite_matrix(working_sff),&
      'fragment self block contains nonfinite values',stage_ok,message)
    if(.not.stage_ok)return
    call move_alloc(working_hff,hff);call move_alloc(working_sff,sff)
    ok=.true.;message=''
  end subroutine extract_dg_hybrid_fragment_self_block

  subroutine compose_dg_hybrid_complete_rows(comm,row_ids,fixed_payload,local_rows,union_to_complete,&
      hamiltonian_rows,metric_rows,operator_fingerprint,ok,message,metric_tolerance,&
      complete_map_fingerprint,complete_map_rank,complete_transform_binding_fingerprint)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    type(s_dg_hybrid_fixed_payload),intent(in)::fixed_payload
    complex(real64),intent(in)::local_rows(:,:),union_to_complete(:,:)
    complex(real64),allocatable,intent(out)::hamiltonian_rows(:,:),metric_rows(:,:)
    integer(int64),intent(out)::operator_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),intent(in),optional::metric_tolerance
    integer(int64),intent(in),optional::complete_map_fingerprint,complete_transform_binding_fingerprint
    integer,intent(in),optional::complete_map_rank
    type(s_dg_hybrid_variational_iterate)::iterate
    integer::rank,ierr,nunion,ncomplete,a,b,p,source_row,owner,output_position,allocation_status,&
      transform_shape(2),minimum_shape(2),maximum_shape(2)
    integer,allocatable::output_owners(:),source_counts(:)
    integer(int64)::matrix_elements,local_hash,global_hash,minimum_map_fingerprint,maximum_map_fingerprint,&
      minimum_binding_fingerprint,maximum_binding_fingerprint,recomputed_binding_fingerprint
    real(real64)::effective_metric_tolerance,minimum_tolerance,maximum_tolerance
    integer::minimum_complete_map_rank,maximum_complete_map_rank
    complex(real64),allocatable::reference_transform(:,:),right_h(:,:),right_s(:,:),local_h(:),local_s(:),&
      reduced_h(:),reduced_s(:),working_h(:,:),working_s(:,:)
    logical::stage_ok

    ok=.false.;message='';operator_fingerprint=0_int64
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='complete operator rank query failed';return;endif
    call validate_payload_contract(comm,fixed_payload,local_rows,stage_ok,message)
    if(.not.stage_ok)return
    call collective_gate(comm,present(metric_tolerance).and.present(complete_map_fingerprint).and.&
      present(complete_map_rank).and.present(complete_transform_binding_fingerprint),&
      'terminal map provenance, binding, rank, and metric tolerance are required',stage_ok,message)
    if(.not.stage_ok)return
    nunion=fixed_payload%global_basis_count
    transform_shape=shape(union_to_complete)
    call MPI_Allreduce(transform_shape,minimum_shape,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(transform_shape,maximum_shape,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal transform shape agreement failed';return;endif
    call collective_gate(comm,all(minimum_shape==maximum_shape),&
      'terminal transform shape differs between ranks',stage_ok,message)
    if(.not.stage_ok)return
    ncomplete=transform_shape(2)
    effective_metric_tolerance=metric_tolerance
    call MPI_Allreduce(effective_metric_tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(effective_metric_tolerance,maximum_tolerance,1,&
      MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_map_fingerprint,minimum_map_fingerprint,1,&
      MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_map_fingerprint,maximum_map_fingerprint,1,&
      MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_transform_binding_fingerprint,minimum_binding_fingerprint,1,&
      MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_transform_binding_fingerprint,maximum_binding_fingerprint,1,&
      MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_map_rank,minimum_complete_map_rank,1,&
      MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(complete_map_rank,maximum_complete_map_rank,1,&
      MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal map provenance agreement failed';return;endif
    call collective_gate(comm,minimum_complete_map_rank==maximum_complete_map_rank.and.&
      complete_map_rank==ncomplete,&
      'terminal transform rank disagrees with certified complete map rank',stage_ok,message)
    if(.not.stage_ok)return
    matrix_elements=int(size(union_to_complete,1),int64)*int(ncomplete,int64)
    call collective_gate(comm,size(union_to_complete,1)==nunion.and.ncomplete>=1.and.ncomplete<=nunion.and.&
      matrix_elements<=int(huge(0),int64).and.finite_matrix(union_to_complete).and.&
      ieee_is_finite(effective_metric_tolerance).and.effective_metric_tolerance>0d0.and.&
      minimum_tolerance==maximum_tolerance.and.minimum_map_fingerprint==maximum_map_fingerprint.and.&
      complete_map_fingerprint/=0_int64.and.&
      minimum_binding_fingerprint==maximum_binding_fingerprint.and.&
      complete_transform_binding_fingerprint/=0_int64,&
      'invalid union-to-complete transform',stage_ok,message)
    if(.not.stage_ok)return
    allocate(reference_transform(nunion,ncomplete),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'union-to-complete agreement workspace allocation failed',stage_ok,message)
    if(.not.stage_ok)return
    if(rank==0)reference_transform=union_to_complete
    call MPI_Bcast(reference_transform,int(matrix_elements),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='union-to-complete agreement broadcast failed';return;endif
    call collective_gate(comm,all(reference_transform==union_to_complete),&
      'union-to-complete transform differs between ranks',stage_ok,message)
    if(.not.stage_ok)return
    call compute_dg_hybrid_union_to_complete_binding(comm,union_to_complete,complete_map_fingerprint,&
      effective_metric_tolerance,recomputed_binding_fingerprint,stage_ok,message)
    if(.not.stage_ok)return
    call collective_gate(comm,recomputed_binding_fingerprint==complete_transform_binding_fingerprint,&
      'terminal transform binding fingerprint disagrees with the certified complete map',stage_ok,message)
    if(.not.stage_ok)return
    call validate_source_ownership(comm,nunion,fixed_payload%row_ids,source_counts,stage_ok,message)
    if(.not.stage_ok)return
    call validate_terminal_transform(comm,fixed_payload,union_to_complete,effective_metric_tolerance,&
      stage_ok,message)
    if(.not.stage_ok)return
    call build_output_directory(comm,ncomplete,row_ids,output_owners,stage_ok,message)
    if(.not.stage_ok)return
    call compose_dg_hybrid_variational_hamiltonian(comm,fixed_payload,local_rows,1d0,1,iterate,&
      stage_ok,message)
    if(.not.stage_ok)return

    allocate(right_h(size(fixed_payload%row_ids),ncomplete),&
      right_s(size(fixed_payload%row_ids),ncomplete),local_h(ncomplete),local_s(ncomplete),&
      reduced_h(ncomplete),reduced_s(ncomplete),working_h(size(row_ids),ncomplete),&
      working_s(size(row_ids),ncomplete),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'complete operator workspace allocation failed',stage_ok,message)
    if(.not.stage_ok)return
    right_h=matmul(iterate%hamiltonian_rows,union_to_complete)
    right_s=matmul(fixed_payload%metric_rows,union_to_complete)
    working_h=(0d0,0d0);working_s=(0d0,0d0)
    do a=1,ncomplete
      local_h=(0d0,0d0);local_s=(0d0,0d0)
      do p=1,size(fixed_payload%row_ids)
        source_row=int(fixed_payload%row_ids(p))
        local_h=local_h+conjg(union_to_complete(source_row,a))*right_h(p,:)
        local_s=local_s+conjg(union_to_complete(source_row,a))*right_s(p,:)
      enddo
      owner=output_owners(a)-1
      call MPI_Reduce(local_h,reduced_h,ncomplete,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='complete Hamiltonian row reduction failed';return;endif
      call MPI_Reduce(local_s,reduced_s,ncomplete,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='complete metric row reduction failed';return;endif
      if(rank==owner)then
        output_position=findloc(row_ids,int(a,int64),dim=1)
        if(output_position<1)then;message='complete output row owner is inconsistent';return;endif
        working_h(output_position,:)=reduced_h;working_s(output_position,:)=reduced_s
      endif
    enddo
    call collective_gate(comm,finite_matrix(working_h).and.finite_matrix(working_s),&
      'complete operator contains nonfinite values',stage_ok,message)
    if(.not.stage_ok)return

    local_hash=0_int64
    do p=1,size(fixed_payload%row_ids)
      do b=1,nunion
        call hash_complex_element(local_hash,int(fixed_payload%row_ids(p)),b,1,&
          iterate%hamiltonian_rows(p,b))
        call hash_complex_element(local_hash,int(fixed_payload%row_ids(p)),b,2,&
          fixed_payload%metric_rows(p,b))
      enddo
    enddo
    if(rank==0)then
      call hash_integer_element(local_hash,1,nunion)
      call hash_integer_element(local_hash,2,ncomplete)
      call hash_integer64_element(local_hash,3,fixed_payload%basis_fingerprint)
      call hash_integer64_element(local_hash,4,fixed_payload%metric_fingerprint)
      call hash_integer64_element(local_hash,5,fixed_payload%interface_fingerprint)
      call hash_integer64_element(local_hash,6,complete_map_fingerprint)
      call hash_integer64_element(local_hash,7,complete_transform_binding_fingerprint)
      call hash_integer_element(local_hash,8,complete_map_rank)
      do a=1,nunion;do b=1,ncomplete
        call hash_complex_element(local_hash,a,b,3,union_to_complete(a,b))
      enddo;enddo
    endif
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='complete operator fingerprint reduction failed';return;endif
    if(global_hash==0_int64)global_hash=int(z'243F6A8885A308D3',int64)
    call move_alloc(working_h,hamiltonian_rows);call move_alloc(working_s,metric_rows)
    operator_fingerprint=global_hash;ok=.true.;message=''
  end subroutine compose_dg_hybrid_complete_rows

  subroutine validate_payload_contract(comm,payload,local_rows,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    complex(real64),intent(in)::local_rows(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::metadata(2),minimum_metadata(2),maximum_metadata(2),ierr
    integer(int64)::fingerprints(4),minimum_fingerprints(4),maximum_fingerprints(4)
    logical::local_ok
    local_ok=payload%frozen.and.allocated(payload%row_ids).and.allocated(payload%metric_rows).and.&
      allocated(payload%kinetic_rows).and.allocated(payload%nonlocal_rows).and.&
      allocated(payload%interface_rows)
    call collective_gate(comm,local_ok,'fixed variational payload is not frozen and allocated',ok,message)
    if(.not.ok)return
    metadata=[payload%global_basis_count,size(payload%row_ids)]
    call MPI_Allreduce(metadata,minimum_metadata,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(metadata,maximum_metadata,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fixed payload metadata agreement failed';return;endif
    local_ok=minimum_metadata(1)==maximum_metadata(1).and.minimum_metadata(1)>=1.and.&
      size(payload%metric_rows,1)==size(payload%row_ids).and.&
      size(payload%metric_rows,2)==payload%global_basis_count.and.&
      all(shape(payload%kinetic_rows)==shape(payload%metric_rows)).and.&
      all(shape(payload%nonlocal_rows)==shape(payload%metric_rows)).and.&
      all(shape(payload%interface_rows)==shape(payload%metric_rows)).and.&
      all(shape(local_rows)==shape(payload%metric_rows)).and.&
      all(payload%row_ids>=1_int64).and.all(payload%row_ids<=int(payload%global_basis_count,int64)).and.&
      finite_matrix(payload%metric_rows).and.finite_matrix(payload%kinetic_rows).and.&
      finite_matrix(payload%nonlocal_rows).and.finite_matrix(payload%interface_rows).and.finite_matrix(local_rows)
    call collective_gate(comm,local_ok,'invalid fixed variational payload rows',ok,message)
    if(.not.ok)return
    fingerprints=[payload%basis_fingerprint,payload%metric_fingerprint,&
      payload%interface_fingerprint,payload%fingerprint]
    call MPI_Allreduce(fingerprints,minimum_fingerprints,4,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(fingerprints,maximum_fingerprints,4,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fixed payload fingerprint agreement failed';return;endif
    call collective_gate(comm,all(minimum_fingerprints==maximum_fingerprints).and.&
      all(fingerprints/=0_int64),'invalid or rank-disagreeing fixed payload fingerprint',ok,message)
    if(.not.ok)return
    call verify_dg_hybrid_variational_payload_rows(comm,payload,ok,message)
  end subroutine validate_payload_contract

  subroutine validate_fragment_contract(comm,fragment_id,basis_fragment,basis_local_slot,basis_generation,&
      fragment_catalog_fingerprint,fragment_directory_fingerprint,expected_catalog_fingerprint,&
      expected_directory_fingerprint,nunion,ok,message)
    integer,intent(in)::comm,fragment_id,basis_fragment(:),basis_local_slot(:),basis_generation(:),nunion
    integer(int64),intent(in)::fragment_catalog_fingerprint,fragment_directory_fingerprint,&
      expected_catalog_fingerprint,expected_directory_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,ierr,metadata(2),minimum_metadata(2),maximum_metadata(2),allocation_status,p,nfragment
    integer,allocatable::reference(:,:),slot_counts(:)
    integer(int64)::fingerprints(2),minimum_fingerprints(2),maximum_fingerprints(2),computed_directory_fingerprint
    logical::local_ok
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fragment contract rank query failed';return;endif
    metadata=[fragment_id,size(basis_fragment)]
    call MPI_Allreduce(metadata,minimum_metadata,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(metadata,maximum_metadata,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fragment contract agreement failed';return;endif
    local_ok=all(minimum_metadata==maximum_metadata).and.fragment_id>=1.and.&
      size(basis_fragment)==nunion.and.size(basis_local_slot)==nunion.and.size(basis_generation)==nunion.and.&
      all(basis_fragment>=1).and.all(basis_local_slot>=1).and.all(basis_generation>=1).and.&
      count(basis_fragment==fragment_id)>=1.and.fragment_catalog_fingerprint/=0_int64.and.&
      fragment_directory_fingerprint/=0_int64.and.fragment_catalog_fingerprint==expected_catalog_fingerprint.and.&
      fragment_directory_fingerprint==expected_directory_fingerprint
    call collective_gate(comm,local_ok,'invalid fragment basis directory',ok,message)
    if(.not.ok)return
    fingerprints=[fragment_catalog_fingerprint,fragment_directory_fingerprint]
    call MPI_Allreduce(fingerprints,minimum_fingerprints,2,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(fingerprints,maximum_fingerprints,2,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_fingerprints/=maximum_fingerprints))then
      ok=.false.;message='fragment catalog fingerprint differs between ranks';return
    endif
    nfragment=count(basis_fragment==fragment_id)
    allocate(reference(3,nunion),slot_counts(nfragment),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'fragment directory agreement allocation failed',ok,message)
    if(.not.ok)return
    if(rank==0)then
      reference(1,:)=basis_fragment;reference(2,:)=basis_local_slot;reference(3,:)=basis_generation
    endif
    call MPI_Bcast(reference,3*nunion,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fragment directory agreement broadcast failed';return;endif
    call collective_gate(comm,all(reference(1,:)==basis_fragment).and.&
      all(reference(2,:)==basis_local_slot).and.all(reference(3,:)==basis_generation),&
      'fragment basis directory differs between ranks',ok,message)
    if(.not.ok)return
    slot_counts=0
    do p=1,nunion
      if(basis_fragment(p)/=fragment_id)cycle
      if(basis_local_slot(p)<=nfragment)slot_counts(basis_local_slot(p))=slot_counts(basis_local_slot(p))+1
    enddo
    call collective_gate(comm,all(slot_counts==1).and.&
      minval(pack(basis_generation,basis_fragment==fragment_id))==&
      maxval(pack(basis_generation,basis_fragment==fragment_id)),&
      'fragment local slots or generation are inconsistent',ok,message)
    if(.not.ok)return
    computed_directory_fingerprint=dg_hybrid_fragment_directory_fingerprint(&
      basis_fragment,basis_local_slot,basis_generation,fragment_catalog_fingerprint)
    call collective_gate(comm,computed_directory_fingerprint==fragment_directory_fingerprint,&
      'fragment directory fingerprint changed',ok,message)
  end subroutine validate_fragment_contract

  pure integer(int64) function dg_hybrid_fragment_directory_fingerprint(&
      basis_fragment,basis_local_slot,basis_generation,fragment_catalog_fingerprint)result(fingerprint)
    integer,intent(in)::basis_fragment(:),basis_local_slot(:),basis_generation(:)
    integer(int64),intent(in)::fragment_catalog_fingerprint
    integer::p
    fingerprint=0_int64
    if(size(basis_fragment)/=size(basis_local_slot).or.size(basis_fragment)/=size(basis_generation).or.&
        size(basis_fragment)<1.or.fragment_catalog_fingerprint==0_int64)return
    fingerprint=int(z'6A09E667F3BCC909',int64)
    call mix_hash_word(fingerprint,fragment_catalog_fingerprint,1_int64)
    call mix_hash_word(fingerprint,int(size(basis_fragment),int64),2_int64)
    do p=1,size(basis_fragment)
      call mix_hash_word(fingerprint,int(p,int64),3_int64)
      call mix_hash_word(fingerprint,int(basis_fragment(p),int64),4_int64)
      call mix_hash_word(fingerprint,int(basis_local_slot(p),int64),5_int64)
      call mix_hash_word(fingerprint,int(basis_generation(p),int64),6_int64)
    enddo
    if(fingerprint==0_int64)fingerprint=1543_int64
  end function dg_hybrid_fragment_directory_fingerprint

  subroutine validate_terminal_transform(comm,payload,transform,metric_tolerance,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::metric_tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nunion,ncomplete,p,a,b,source_row,ierr,allocation_status
    complex(real64),allocatable::column_gram(:,:),right_metric(:,:),&
      local_retained_metric(:,:),retained_metric(:,:)
    complex(real64)::expected
    real(real64)::column_defect,map_tolerance
    logical::local_ok,identity_map

    nunion=size(transform,1);ncomplete=size(transform,2)
    allocate(column_gram(ncomplete,ncomplete),right_metric(size(payload%row_ids),ncomplete),&
      local_retained_metric(ncomplete,ncomplete),retained_metric(ncomplete,ncomplete),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'terminal transform validation allocation failed',ok,message)
    if(.not.ok)return
    map_tolerance=max(100d0*metric_tolerance,&
      100d0*epsilon(1d0)*real(max(nunion,ncomplete),real64))
    column_gram=matmul(conjg(transpose(transform)),transform);column_defect=0d0
    do a=1,ncomplete;do b=1,ncomplete
      expected=merge(cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),a==b)
      column_defect=max(column_defect,abs(column_gram(a,b)-expected))
    enddo;enddo
    call collective_gate(comm,ieee_is_finite(column_defect).and.column_defect<=map_tolerance,&
      'terminal transform columns are not orthonormal',ok,message)
    if(.not.ok)return
    identity_map=.true.
    if(ncomplete==nunion)then
      do a=1,nunion;do b=1,nunion
        expected=merge(cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),a==b)
        identity_map=identity_map.and.transform(a,b)==expected
      enddo;enddo
    endif
    call collective_gate(comm,ncomplete<nunion.or.identity_map,&
      'full-rank terminal transform must be bitwise identity',ok,message)
    if(.not.ok)return
    call certify_terminal_metric_range(comm,payload,transform,metric_tolerance,map_tolerance,ok,message)
    if(.not.ok)return
    right_metric=matmul(payload%metric_rows,transform)
    local_retained_metric=(0d0,0d0)
    do p=1,size(payload%row_ids)
      source_row=int(payload%row_ids(p))
      do a=1,ncomplete;do b=1,ncomplete
        local_retained_metric(a,b)=local_retained_metric(a,b)+&
          conjg(transform(source_row,a))*right_metric(p,b)
      enddo;enddo
    enddo
    call MPI_Allreduce(local_retained_metric,retained_metric,ncomplete*ncomplete,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='retained terminal metric reduction failed';return;endif
    local_ok=positive_definite_hermitian(retained_metric,metric_tolerance,map_tolerance)
    call collective_gate(comm,local_ok,'terminal transform retained metric is not positive rank',ok,message)
  end subroutine validate_terminal_transform

  subroutine certify_terminal_metric_range(comm,payload,transform,metric_tolerance,map_tolerance,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::metric_tolerance,map_tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::local_row(:),reduced_row(:),metric_vectors(:,:),work(:),&
      retained_projector(:,:),transform_projector(:,:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::scale,cutoff,negative_limit,roundoff_floor,hermitian_defect,projector_defect
    integer::rank,ierr,nunion,ncomplete,a,b,k,source_position,root_extent,retained_rank,&
      allocation_status,info
    logical::local_ok

    ok=.false.;message='';retained_rank=0;projector_defect=0d0
    nunion=size(transform,1);ncomplete=size(transform,2)
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal metric certification rank query failed';return;endif
    allocate(local_row(nunion),reduced_row(nunion),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'terminal metric row workspace allocation failed',ok,message)
    if(.not.ok)return
    root_extent=merge(nunion,0,rank==0)
    allocate(metric_vectors(root_extent,root_extent),eigenvalues(root_extent),&
      work(max(1,2*root_extent-1)),rwork(max(1,3*root_extent-2)),&
      retained_projector(root_extent,root_extent),transform_projector(root_extent,root_extent),&
      stat=allocation_status)
    call collective_gate(comm,allocation_status==0,&
      'terminal metric eigensystem workspace allocation failed',ok,message)
    if(.not.ok)return
    do a=1,nunion
      local_row=(0d0,0d0)
      source_position=findloc(payload%row_ids,int(a,int64),dim=1)
      if(source_position>0)local_row=payload%metric_rows(source_position,:)
      call MPI_Reduce(local_row,reduced_row,nunion,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal metric row collection failed';return;endif
      if(rank==0)metric_vectors(a,:)=reduced_row
    enddo
    local_ok=.true.;info=0
    if(rank==0)then
      scale=max(1d0,maxval(abs(metric_vectors)))
      hermitian_defect=maxval(abs(metric_vectors-conjg(transpose(metric_vectors))))
      local_ok=finite_matrix(metric_vectors).and.&
        hermitian_defect<=100d0*epsilon(1d0)*real(nunion,real64)*scale
      if(local_ok)then
        metric_vectors=0.5d0*(metric_vectors+conjg(transpose(metric_vectors)))
        call zheev('V','U',nunion,metric_vectors,nunion,eigenvalues,work,size(work),rwork,info)
        local_ok=info==0.and.all(ieee_is_finite(eigenvalues))
      endif
      if(local_ok)then
        scale=max(1d0,maxval(abs(eigenvalues)))
        roundoff_floor=64d0*epsilon(1d0)*scale*real(max(1,nunion),real64)
        cutoff=max(metric_tolerance*scale,roundoff_floor);negative_limit=roundoff_floor
        local_ok=minval(eigenvalues)>=-negative_limit
        if(local_ok)local_ok=.not.any(abs(eigenvalues-cutoff)<=16d0*roundoff_floor)
        retained_rank=count(eigenvalues>cutoff)
      endif
    endif
    call collective_gate(comm,local_ok,'terminal metric eigensystem certification failed',ok,message)
    if(.not.ok)return
    call MPI_Bcast(retained_rank,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal metric rank broadcast failed';return;endif
    call collective_gate(comm,retained_rank==ncomplete,&
      'terminal transform rank disagrees with metric cutoff rank',ok,message)
    if(.not.ok)return
    local_ok=.true.
    if(rank==0)then
      retained_projector=(0d0,0d0)
      do k=1,nunion
        if(eigenvalues(k)<=cutoff)cycle
        do b=1,nunion;do a=1,nunion
          retained_projector(a,b)=retained_projector(a,b)+&
            metric_vectors(a,k)*conjg(metric_vectors(b,k))
        enddo;enddo
      enddo
      transform_projector=matmul(transform,conjg(transpose(transform)))
      projector_defect=maxval(abs(transform_projector-retained_projector))
      local_ok=ieee_is_finite(projector_defect).and.projector_defect<=map_tolerance
    endif
    call collective_gate(comm,local_ok,&
      'terminal transform span disagrees with retained metric range',ok,message)
  end subroutine certify_terminal_metric_range

  logical function positive_definite_hermitian(matrix,rank_tolerance,hermitian_tolerance)result(positive)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),intent(in)::rank_tolerance,hermitian_tolerance
    complex(real64),allocatable::factor(:,:)
    complex(real64)::pivot_value,value
    real(real64)::scale,pivot_threshold
    integer::j,k,q,allocation_status
    positive=.false.
    if(size(matrix,1)/=size(matrix,2).or..not.finite_matrix(matrix))return
    scale=max(1d0,maxval(abs(matrix)))
    if(maxval(abs(matrix-conjg(transpose(matrix))))>hermitian_tolerance*scale)return
    allocate(factor(size(matrix,1),size(matrix,2)),stat=allocation_status)
    if(allocation_status/=0)return
    factor=(0d0,0d0);pivot_threshold=max(rank_tolerance,&
      100d0*epsilon(1d0)*real(max(1,size(matrix,1)),real64))*scale
    do j=1,size(matrix,1)
      pivot_value=matrix(j,j)
      do q=1,j-1;pivot_value=pivot_value-factor(j,q)*conjg(factor(j,q));enddo
      if(abs(aimag(pivot_value))>hermitian_tolerance*scale.or.real(pivot_value,real64)<=pivot_threshold)return
      factor(j,j)=sqrt(real(pivot_value,real64))
      do k=j+1,size(matrix,1)
        value=matrix(k,j)
        do q=1,j-1;value=value-factor(k,q)*conjg(factor(j,q));enddo
        factor(k,j)=value/factor(j,j)
      enddo
    enddo
    positive=.true.
  end function positive_definite_hermitian

  subroutine validate_source_ownership(comm,global_count,row_ids,counts,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    integer,allocatable,intent(out)::counts(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::p,ierr,allocation_status
    logical::local_ok
    allocate(counts(global_count),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,'source-row ownership allocation failed',ok,message)
    if(.not.ok)return
    counts=0;local_ok=all(row_ids>=1_int64).and.all(row_ids<=int(global_count,int64))
    call collective_gate(comm,local_ok,'invalid source-row ID',ok,message)
    if(.not.ok)return
    do p=1,size(row_ids);counts(int(row_ids(p)))=counts(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='source-row ownership reduction failed';return;endif
    call collective_gate(comm,all(counts==1),'source rows are not owned exactly once',ok,message)
  end subroutine validate_source_ownership

  subroutine build_output_directory(comm,global_count,row_ids,owners,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    integer,allocatable,intent(out)::owners(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,p,ierr,allocation_status
    integer,allocatable::counts(:)
    logical::local_ok
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='output-row rank query failed';return;endif
    allocate(owners(global_count),counts(global_count),stat=allocation_status)
    call collective_gate(comm,allocation_status==0,'output rows directory allocation failed',ok,message)
    if(.not.ok)return
    owners=0;counts=0
    local_ok=all(row_ids>=1_int64).and.all(row_ids<=int(global_count,int64))
    call collective_gate(comm,local_ok,'invalid output rows',ok,message)
    if(.not.ok)return
    do p=1,size(row_ids)
      counts(int(row_ids(p)))=counts(int(row_ids(p)))+1
      owners(int(row_ids(p)))=rank+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='output rows directory reduction failed';return;endif
    call collective_gate(comm,all(counts==1).and.all(owners>=1),&
      'output rows are not owned exactly once',ok,message)
  end subroutine build_output_directory

  subroutine collective_gate(comm,local_ok,label,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    character(*),intent(in)::label
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,local_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0;message=''
    if(ierr/=MPI_SUCCESS)message=trim(label)//' reduction failed'
    if(ierr==MPI_SUCCESS.and.global_bad/=0)message=label
  end subroutine collective_gate

  subroutine hash_complex_element(hash,row,column,component,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::row,column,component
    complex(real64),intent(in)::value
    integer(int64)::bits,element_hash
    element_hash=int(z'BB67AE8584CAA73B',int64)
    call mix_hash_word(element_hash,int(row,int64),1_int64)
    call mix_hash_word(element_hash,int(column,int64),2_int64)
    call mix_hash_word(element_hash,int(component,int64),3_int64)
    bits=transfer(real(value,real64),bits);call mix_hash_word(element_hash,bits,4_int64)
    bits=transfer(aimag(value),bits);call mix_hash_word(element_hash,bits,5_int64)
    hash=ieor(hash,element_hash)
  end subroutine hash_complex_element

  subroutine hash_integer_element(hash,key,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::key,value
    call hash_integer64_element(hash,key,int(value,int64))
  end subroutine hash_integer_element

  subroutine hash_integer64_element(hash,key,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::key
    integer(int64),intent(in)::value
    hash=ieor(hash,ishftc(ieor(value,int(z'9E3779B97F4A7C15',int64)),mod(23*key,63)))
  end subroutine hash_integer64_element

  pure subroutine mix_hash_word(hash,word,tag)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::word,tag
    integer::shift
    shift=1+int(modulo(ieor(word,ishftc(tag,7)),63_int64))
    hash=ieor(ishftc(hash,shift),word)
    hash=ieor(hash,ishftc(tag,modulo(shift+23,64)))
    hash=ieor(hash,int(z'13198A2E03707344',int64))
  end subroutine mix_hash_word

  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values,real64))).and.&
      all(ieee_is_finite(aimag(values)))
  end function finite_matrix
end module dg_hybrid_divided_operator
