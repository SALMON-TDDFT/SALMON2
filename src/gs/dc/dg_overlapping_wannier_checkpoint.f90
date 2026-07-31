#include "config.h"
module dg_overlapping_wannier_checkpoint
  use iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  character(32),parameter::manifest_magic='SALMON_OW_GS_CHECKPOINT_V3'
  character(32),parameter::shard_magic='SALMON_OW_GS_RANK_SHARD_V3'
  type,public::s_dg_overlapping_wannier_checkpoint
    integer::basis_generation=0,geometry_generation=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64
    integer(int64)::hamiltonian_fingerprint=0_int64,observable_fingerprint=0_int64
    integer(int64)::publication_id=0_int64
    character(64)::field_coupling_convention=''
    integer,allocatable::center_owner(:),tail_center(:),tail_generation(:),tail_offsets(:)
    integer(int64),allocatable::tail_physical_ids(:),core_physical_ids(:),overlap_row_ids(:)
    complex(8),allocatable::overlap(:,:),hamiltonian0(:,:),position(:,:,:),velocity(:,:,:),coefficients(:,:)
    real(8),allocatable::occupations(:),density(:)
    real(8)::density_residual=huge(1d0),coefficient_residual=huge(1d0),charge_error=huge(1d0)
    real(8)::unmixed_density_residual=huge(1d0),orthogonality_defect=huge(1d0)
    real(8)::metric_condition=huge(1d0)
    real(8)::density_tolerance=0d0,coefficient_tolerance=0d0,orthogonality_tolerance=0d0
    real(8)::charge_tolerance=0d0,condition_limit=0d0
    real(8)::symmetry_closure_residual=huge(1d0),symmetry_tolerance=0d0
    logical::accepted=.false.
  end type
  public::write_dg_overlapping_wannier_checkpoint,read_dg_overlapping_wannier_checkpoint
  public::compute_dg_overlapping_wannier_matrix_fingerprints
contains
  subroutine write_dg_overlapping_wannier_checkpoint(comm,prefix,checkpoint,ok,message,failure_injection_rank)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    type(s_dg_overlapping_wannier_checkpoint),intent(inout)::checkpoint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::failure_injection_rank
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,ios,flush_ios,close_ios,local_bad,global_bad,dims(4),&
      local_rejection_code,global_rejection_code
    integer(int64)::transaction_id,previous_publication_id,shard_size,manifest_size,shard_digest,manifest_digest,&
      computed_hamiltonian_fingerprint,computed_observable_fingerprint,digest_position,payload_position,&
      integer_bytes,integer8_bytes,file_bytes
    character(512)::manifest,manifest_temporary,shard,shard_temporary,reservation
    character(256)::iomsg
    logical::matrix_ok
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    ok=.false.;message='';local_bad=0
    previous_publication_id=checkpoint%publication_id
    call validate_checkpoint(checkpoint,nproc,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad==0)then
      call compute_dg_overlapping_wannier_matrix_fingerprints(comm,checkpoint%overlap_row_ids,&
        checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity,&
        computed_hamiltonian_fingerprint,computed_observable_fingerprint,matrix_ok)
      if(.not.matrix_ok.or.computed_hamiltonian_fingerprint/=checkpoint%hamiltonian_fingerprint.or.&
         computed_observable_fingerprint/=checkpoint%observable_fingerprint)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_rejection_code=checkpoint_rejection_code(checkpoint,nproc)
    if(local_bad/=0.and.local_rejection_code==0)local_rejection_code=ibset(local_rejection_code,15)
    call MPI_Allreduce(local_rejection_code,global_rejection_code,1,MPI_INTEGER,MPI_BOR,comm,ierr)
    if(global_bad/=0.or.global_rejection_code/=0)then
      write(message,'(a,i0,a,l1,3(a,es12.4))')&
        'invalid or unaccepted overlapping-Wannier checkpoint code=',global_rejection_code,&
        ' accepted=',checkpoint%accepted,' density=',checkpoint%density_residual,&
        ' unmixed=',checkpoint%unmixed_density_residual,' tolerance=',checkpoint%density_tolerance
      return
    endif
    call validate_replicated_payload(comm,checkpoint,local_bad)
    if(local_bad/=0)then;message='rank-inconsistent overlapping-Wannier checkpoint manifest payload';return;endif
    call validate_shard_ownership(comm,checkpoint,local_bad)
    if(local_bad/=0)then;message='incomplete overlapping-Wannier checkpoint shard ownership';return;endif
    manifest=trim(prefix)//'.manifest'
    if(rank==0)call choose_transaction_id(prefix,checkpoint%basis_generation,&
      checkpoint%operator_fingerprint,transaction_id,reservation,ios)
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0)then;message='cannot reserve overlapping-Wannier publication transaction';return;endif
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    write(manifest_temporary,'(a,".temporary.",z16.16)')trim(manifest),transaction_id
    checkpoint%publication_id=transaction_id
    shard_digest=digest_shard(checkpoint,rank,transaction_id)
    shard_size=expected_shard_size(checkpoint)
    call versioned_shard_name(prefix,checkpoint%basis_generation,checkpoint%operator_fingerprint,&
      transaction_id,rank,shard)
    shard_temporary=trim(shard)//'.temporary'
    open(newunit=unit,file=trim(shard_temporary),status='replace',access='stream',form='unformatted',&
      action='write',iostat=ios,iomsg=iomsg)
    if(ios==0)then
      write(unit,iostat=ios,iomsg=iomsg)shard_magic,rank,nproc,checkpoint%basis_generation,&
        checkpoint%geometry_generation,checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,transaction_id,&
        shard_size,shard_digest,size(checkpoint%tail_center),size(checkpoint%tail_physical_ids),&
        size(checkpoint%core_physical_ids),size(checkpoint%overlap_row_ids)
      if(ios==0)write(unit,iostat=ios,iomsg=iomsg)checkpoint%tail_center,checkpoint%tail_generation,&
        checkpoint%tail_offsets,checkpoint%tail_physical_ids,checkpoint%core_physical_ids,checkpoint%density,&
        checkpoint%overlap_row_ids,checkpoint%overlap,checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity
      flush_ios=0;close_ios=0
      if(ios==0)flush(unit,iostat=flush_ios)
      close(unit,iostat=close_ios)
      if(ios==0.and.flush_ios/=0)ios=flush_ios
      if(ios==0.and.close_ios/=0)ios=close_ios
    endif
    if(ios==0)then
      integer_bytes=int(storage_size(0)/8,int64)
      integer8_bytes=int(storage_size(0_int64)/8,int64)
      digest_position=int(storage_size(shard_magic)/8,int64)+4_int64*integer_bytes+&
        4_int64*integer8_bytes+1_int64
      payload_position=int(storage_size(shard_magic)/8,int64)+8_int64*integer_bytes+&
        5_int64*integer8_bytes+1_int64
      if(ios==0)open(newunit=unit,file=trim(shard_temporary),status='old',access='stream',form='unformatted',&
        action='readwrite',iostat=ios,iomsg=iomsg)
      if(ios==0)then
        inquire(unit=unit,size=file_bytes)
        if(file_bytes/=shard_size)ios=1
        if(ios==0)then
          call digest_serialized_shard(unit,payload_position,checkpoint,rank,transaction_id,shard_digest,ios,iomsg)
        endif
        if(ios==0)then
          write(unit,pos=digest_position,iostat=ios,iomsg=iomsg)shard_digest
        endif
        flush_ios=0;close_ios=0
        if(ios==0)flush(unit,iostat=flush_ios)
        close(unit,iostat=close_ios)
        if(ios==0.and.flush_ios/=0)ios=flush_ios
        if(ios==0.and.close_ios/=0)ios=close_ios
      endif
    endif
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='cannot write overlapping-Wannier rank shard'
      call abort_publication(rank,manifest_temporary,shard_temporary,shard,reservation,checkpoint,previous_publication_id)
      return
    endif
    call rename(trim(shard_temporary),trim(shard),ios)
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='cannot atomically publish overlapping-Wannier rank shard'
      call abort_publication(rank,manifest_temporary,shard_temporary,shard,reservation,checkpoint,previous_publication_id)
      return
    endif
    local_bad=0
    if(present(failure_injection_rank))then
      if(rank==failure_injection_rank)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='injected failure before overlapping-Wannier manifest commit'
      call abort_publication(rank,manifest_temporary,shard_temporary,shard,reservation,checkpoint,previous_publication_id)
      return
    endif
    call MPI_Barrier(comm,ierr)
    if(rank==0)then
      dims=[size(checkpoint%center_owner),size(checkpoint%coefficients,1),size(checkpoint%coefficients,2),&
        size(checkpoint%occupations)]
      manifest_digest=digest_manifest(checkpoint,transaction_id)
      manifest_size=expected_manifest_size(checkpoint)
      open(newunit=unit,file=trim(manifest_temporary),status='replace',access='stream',form='unformatted',&
        action='write',iostat=ios,iomsg=iomsg)
      if(ios==0)then
        write(unit,iostat=ios,iomsg=iomsg)manifest_magic,nproc,checkpoint%basis_generation,&
          checkpoint%geometry_generation,checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,&
          checkpoint%hamiltonian_fingerprint,checkpoint%observable_fingerprint,&
          checkpoint%field_coupling_convention,transaction_id,manifest_size,manifest_digest,dims,&
          checkpoint%density_residual,&
          checkpoint%coefficient_residual,checkpoint%charge_error,&
          checkpoint%unmixed_density_residual,checkpoint%orthogonality_defect,checkpoint%metric_condition,&
          checkpoint%density_tolerance,checkpoint%coefficient_tolerance,checkpoint%orthogonality_tolerance,&
          checkpoint%charge_tolerance,checkpoint%condition_limit,&
          checkpoint%symmetry_closure_residual,checkpoint%symmetry_tolerance,&
          checkpoint%accepted
        if(ios==0)write(unit,iostat=ios,iomsg=iomsg)checkpoint%center_owner,&
          checkpoint%coefficients,checkpoint%occupations
        flush_ios=0;close_ios=0
        if(ios==0)flush(unit,iostat=flush_ios)
        close(unit,iostat=close_ios)
        if(ios==0.and.flush_ios/=0)ios=flush_ios
        if(ios==0.and.close_ios/=0)ios=close_ios
      endif
      if(ios==0)call rename(trim(manifest_temporary),trim(manifest),ios)
      call remove_file_if_present(reservation)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0)then
      message='cannot publish overlapping-Wannier manifest'
      call abort_publication(rank,manifest_temporary,shard_temporary,shard,reservation,checkpoint,previous_publication_id)
      return
    endif
    ok=.true.
#else
    ok=.false.;message='overlapping-Wannier checkpoint requires MPI'
#endif
  end subroutine

  subroutine read_dg_overlapping_wannier_checkpoint(comm,prefix,expected_basis_generation,&
      expected_geometry_generation,expected_basis_fingerprint,expected_operator_fingerprint,&
      expected_acceptance_gates,checkpoint,reusable,ok,message)
    integer,intent(in)::comm,expected_basis_generation,expected_geometry_generation
    integer(int64),intent(in)::expected_basis_fingerprint,expected_operator_fingerprint
    real(8),intent(in)::expected_acceptance_gates(6)
    character(*),intent(in)::prefix
    type(s_dg_overlapping_wannier_checkpoint),intent(out)::checkpoint
    logical,intent(out)::reusable,ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,ios,file_nproc,dims(4),local_bad,global_bad,shard_rank,shard_nproc,&
      gate_bad,allocation_status
    integer::integer_metadata(4)
    integer::expected_integers(2),expected_integer_min(2),expected_integer_max(2)
    integer(int64)::fingerprint_metadata(4),count64,file_bytes,transaction_id,shard_transaction_id,&
      stored_size,stored_digest,&
      expected_fingerprints(2),&
      expected_fingerprint_min(2),expected_fingerprint_max(2)
    real(8)::acceptance_metadata(13),gate_min(6),gate_max(6)
    integer::nbasis_gen,ngeometry_gen,ntail_record,ntail_id,ncore,noverlap_row
    integer(int64)::basis_fp,operator_fp,computed_hamiltonian_fingerprint,computed_observable_fingerprint
    character(32)::magic
    character(512)::manifest,shard
    character(256)::iomsg
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    checkpoint=s_dg_overlapping_wannier_checkpoint();reusable=.false.;ok=.false.;message=''
    expected_integers=[expected_basis_generation,expected_geometry_generation]
    expected_fingerprints=[expected_basis_fingerprint,expected_operator_fingerprint]
    call MPI_Allreduce(expected_integers,expected_integer_min,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_integers,expected_integer_max,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_fingerprints,expected_fingerprint_min,2,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_fingerprints,expected_fingerprint_max,2,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_acceptance_gates,gate_min,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_acceptance_gates,gate_max,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    gate_bad=merge(0,1,all(ieee_is_finite(expected_acceptance_gates)).and.&
      minval(expected_acceptance_gates)>0d0)
    call MPI_Allreduce(MPI_IN_PLACE,gate_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(any(expected_integer_min/=expected_integer_max).or.&
       any(expected_fingerprint_min/=expected_fingerprint_max).or.any(gate_min/=gate_max).or.&
       gate_bad/=0)then
      message='rank-inconsistent overlapping-Wannier checkpoint read provenance';return
    endif
    manifest=trim(prefix)//'.manifest';ios=0;dims=0;file_nproc=0
    if(rank==0)then
      open(newunit=unit,file=trim(manifest),status='old',access='stream',form='unformatted',&
        action='read',iostat=ios,iomsg=iomsg)
      if(ios==0)then
        inquire(unit=unit,size=file_bytes)
        if(file_bytes<128)ios=1
        if(ios==0)read(unit,iostat=ios,iomsg=iomsg)magic,file_nproc,checkpoint%basis_generation,&
          checkpoint%geometry_generation,checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,&
          checkpoint%hamiltonian_fingerprint,checkpoint%observable_fingerprint,&
          checkpoint%field_coupling_convention,transaction_id,stored_size,stored_digest,dims,&
          checkpoint%density_residual,&
          checkpoint%coefficient_residual,checkpoint%charge_error,&
          checkpoint%unmixed_density_residual,checkpoint%orthogonality_defect,checkpoint%metric_condition,&
          checkpoint%density_tolerance,checkpoint%coefficient_tolerance,checkpoint%orthogonality_tolerance,&
          checkpoint%charge_tolerance,checkpoint%condition_limit,&
          checkpoint%symmetry_closure_residual,checkpoint%symmetry_tolerance,&
          checkpoint%accepted
        if(ios==0)then
          if(magic/=manifest_magic.or.file_nproc/=nproc.or.dims(1)<=0.or.dims(2)<=0.or.dims(3)<=0.or.&
             dims(4)/=dims(3).or.any(dims>1000000))ios=1
        endif
        if(ios==0.and.(file_bytes/=stored_size.or.stored_size<128_int64))ios=1
        if(ios==0)then
          count64=int(dims(2),int64)*int(dims(3),int64)
          if(count64>100000000_int64)ios=1
        endif
        if(ios==0)then
          allocate(checkpoint%center_owner(dims(1)),checkpoint%coefficients(dims(2),dims(3)),&
            checkpoint%occupations(dims(4)),stat=allocation_status)
          if(allocation_status/=0)ios=1
          if(ios==0)read(unit,iostat=ios,iomsg=iomsg)checkpoint%center_owner,checkpoint%coefficients,&
            checkpoint%occupations
          if(ios==0)then
            if(stored_size/=expected_manifest_size(checkpoint).or.&
               stored_digest/=digest_manifest(checkpoint,transaction_id))ios=1
          endif
        else
          ios=1
        endif
        close(unit)
      endif
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0)then;message='overlapping-Wannier manifest is missing or belongs to another route';return;endif
    call MPI_Bcast(dims,4,MPI_INTEGER,0,comm,ierr)
    allocation_status=0
    if(rank/=0)allocate(checkpoint%center_owner(dims(1)),checkpoint%coefficients(dims(2),dims(3)),&
      checkpoint%occupations(dims(4)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='cannot allocate overlapping-Wannier manifest payload collectively';return
    endif
    if(rank==0)then
      integer_metadata=[checkpoint%basis_generation,checkpoint%geometry_generation,file_nproc,0]
      fingerprint_metadata=[checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,&
        checkpoint%hamiltonian_fingerprint,checkpoint%observable_fingerprint]
      acceptance_metadata(1:6)=[checkpoint%density_residual,checkpoint%coefficient_residual,checkpoint%charge_error,&
        checkpoint%unmixed_density_residual,checkpoint%orthogonality_defect,checkpoint%metric_condition]
      acceptance_metadata(7:11)=[checkpoint%density_tolerance,checkpoint%coefficient_tolerance,&
        checkpoint%orthogonality_tolerance,checkpoint%charge_tolerance,checkpoint%condition_limit]
      acceptance_metadata(12:13)=[checkpoint%symmetry_closure_residual,checkpoint%symmetry_tolerance]
    endif
    call MPI_Bcast(integer_metadata,4,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(fingerprint_metadata,4,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(checkpoint%field_coupling_convention,len(checkpoint%field_coupling_convention),MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(transaction_id,1,MPI_INTEGER8,0,comm,ierr)
    checkpoint%publication_id=transaction_id
    call MPI_Bcast(acceptance_metadata,13,MPI_DOUBLE_PRECISION,0,comm,ierr)
    checkpoint%basis_generation=integer_metadata(1);checkpoint%geometry_generation=integer_metadata(2)
    checkpoint%basis_fingerprint=fingerprint_metadata(1);checkpoint%operator_fingerprint=fingerprint_metadata(2)
    checkpoint%hamiltonian_fingerprint=fingerprint_metadata(3)
    checkpoint%observable_fingerprint=fingerprint_metadata(4)
    checkpoint%density_residual=acceptance_metadata(1);checkpoint%coefficient_residual=acceptance_metadata(2)
    checkpoint%charge_error=acceptance_metadata(3);checkpoint%unmixed_density_residual=acceptance_metadata(4)
    checkpoint%orthogonality_defect=acceptance_metadata(5);checkpoint%metric_condition=acceptance_metadata(6)
    checkpoint%density_tolerance=acceptance_metadata(7);checkpoint%coefficient_tolerance=acceptance_metadata(8)
    checkpoint%orthogonality_tolerance=acceptance_metadata(9);checkpoint%charge_tolerance=acceptance_metadata(10)
    checkpoint%condition_limit=acceptance_metadata(11)
    checkpoint%symmetry_closure_residual=acceptance_metadata(12);checkpoint%symmetry_tolerance=acceptance_metadata(13)
    call MPI_Bcast(checkpoint%accepted,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(checkpoint%center_owner,dims(1),MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(checkpoint%coefficients,dims(2)*dims(3),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(checkpoint%occupations,dims(4),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call versioned_shard_name(prefix,checkpoint%basis_generation,checkpoint%operator_fingerprint,&
      transaction_id,rank,shard)
    open(newunit=unit,file=trim(shard),status='old',access='stream',form='unformatted',&
      action='read',iostat=ios,iomsg=iomsg)
    if(ios==0)then
      inquire(unit=unit,size=file_bytes)
      if(file_bytes<72)ios=1
      if(ios==0)read(unit,iostat=ios,iomsg=iomsg)magic,shard_rank,shard_nproc,nbasis_gen,ngeometry_gen,basis_fp,&
        operator_fp,shard_transaction_id,stored_size,stored_digest,ntail_record,ntail_id,ncore,noverlap_row
      if(ios==0.and.magic==shard_magic.and.shard_rank==rank.and.shard_nproc==nproc.and.&
          ntail_record>=0.and.ntail_id>=0.and.ncore>=0.and.ntail_record<=1000000.and.&
          ntail_id<=10000000.and.ncore<=10000000.and.noverlap_row>=0.and.noverlap_row<=dims(2))then
        allocate(checkpoint%tail_center(ntail_record),checkpoint%tail_generation(ntail_record),&
          checkpoint%tail_offsets(ntail_record+1),checkpoint%tail_physical_ids(ntail_id),&
          checkpoint%core_physical_ids(ncore),checkpoint%density(ncore),&
          checkpoint%overlap_row_ids(noverlap_row),checkpoint%overlap(noverlap_row,dims(2)),&
          checkpoint%hamiltonian0(noverlap_row,dims(2)),checkpoint%position(3,noverlap_row,dims(2)),&
          checkpoint%velocity(3,noverlap_row,dims(2)),&
          stat=allocation_status)
        if(allocation_status/=0)ios=1
        if(ios==0)read(unit,iostat=ios,iomsg=iomsg)checkpoint%tail_center,checkpoint%tail_generation,&
          checkpoint%tail_offsets,checkpoint%tail_physical_ids,checkpoint%core_physical_ids,checkpoint%density,&
          checkpoint%overlap_row_ids,checkpoint%overlap,checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity
        if(ios==0)then
          if(file_bytes/=stored_size.or.stored_size/=expected_shard_size(checkpoint).or.&
             stored_digest/=digest_shard(checkpoint,rank,shard_transaction_id))ios=1
        endif
      else
        ios=1
      endif
      close(unit)
    endif
    local_bad=merge(0,1,ios==0)
    if(local_bad==0)then
      if(nbasis_gen/=checkpoint%basis_generation.or.ngeometry_gen/=checkpoint%geometry_generation.or.&
         basis_fp/=checkpoint%basis_fingerprint.or.operator_fp/=checkpoint%operator_fingerprint.or.&
         shard_transaction_id/=transaction_id)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='overlapping-Wannier checkpoint shard set is partial or stale';return;endif
    call validate_checkpoint(checkpoint,nproc,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='overlapping-Wannier checkpoint payload is invalid';return;endif
    call validate_shard_ownership(comm,checkpoint,local_bad)
    if(local_bad/=0)then;message='overlapping-Wannier checkpoint ownership is incomplete';return;endif
    call compute_dg_overlapping_wannier_matrix_fingerprints(comm,checkpoint%overlap_row_ids,&
      checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity,&
      computed_hamiltonian_fingerprint,computed_observable_fingerprint,reusable)
    if(.not.reusable.or.computed_hamiltonian_fingerprint/=checkpoint%hamiltonian_fingerprint.or.&
       computed_observable_fingerprint/=checkpoint%observable_fingerprint)then
      reusable=.false.;message='overlapping-Wannier matrix fingerprint mismatch';return
    endif
    reusable=(expected_basis_generation==0.or.checkpoint%basis_generation==expected_basis_generation).and.&
      (expected_geometry_generation==0.or.checkpoint%geometry_generation==expected_geometry_generation).and.&
      (expected_basis_fingerprint==0_int64.or.checkpoint%basis_fingerprint==expected_basis_fingerprint).and.&
      (expected_operator_fingerprint==0_int64.or.checkpoint%operator_fingerprint==expected_operator_fingerprint).and.&
      checkpoint%density_residual<=expected_acceptance_gates(1).and.&
      checkpoint%unmixed_density_residual<=expected_acceptance_gates(1).and.&
      checkpoint%coefficient_residual<=expected_acceptance_gates(2).and.&
      checkpoint%orthogonality_defect<=expected_acceptance_gates(3).and.&
      abs(checkpoint%charge_error)<=expected_acceptance_gates(4).and.&
      checkpoint%metric_condition<=expected_acceptance_gates(5).and.&
      checkpoint%symmetry_closure_residual<=expected_acceptance_gates(6)
    ok=.true.
    if(.not.reusable)message='overlapping-Wannier checkpoint provenance mismatch'
#else
    checkpoint=s_dg_overlapping_wannier_checkpoint();reusable=.false.;ok=.false.
    message='overlapping-Wannier checkpoint requires MPI'
#endif
  end subroutine

  subroutine compute_dg_overlapping_wannier_matrix_fingerprints(comm,row_ids,hamiltonian,position,&
      velocity,hamiltonian_fingerprint,observable_fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(8),intent(in)::hamiltonian(:,:),position(:,:,:),velocity(:,:,:)
    integer(int64),intent(out)::hamiltonian_fingerprint,observable_fingerprint
    logical,intent(out)::ok
#ifdef USE_MPI
    integer::i,j,axis,ierr,local_bad,global_bad
    integer(int64)::local_h,local_o
    local_bad=0
    if(size(hamiltonian,1)/=size(row_ids).or.&
       any(shape(position)/=[3,size(row_ids),size(hamiltonian,2)]).or.&
       any(shape(velocity)/=[3,size(row_ids),size(hamiltonian,2)]).or.&
       .not.finite_complex(hamiltonian).or..not.finite_complex3(position).or.&
       .not.finite_complex3(velocity))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      hamiltonian_fingerprint=0_int64;observable_fingerprint=0_int64;ok=.false.;return
    endif
    local_h=0_int64;local_o=0_int64
    do j=1,size(hamiltonian,2);do i=1,size(row_ids)
      local_h=ieor(local_h,matrix_entry_fingerprint(hamiltonian(i,j),row_ids(i),j,0,1))
      do axis=1,3
        local_o=ieor(local_o,matrix_entry_fingerprint(position(axis,i,j),row_ids(i),j,axis,2))
        local_o=ieor(local_o,matrix_entry_fingerprint(velocity(axis,i,j),row_ids(i),j,axis,3))
      enddo
    enddo;enddo
    call MPI_Allreduce(local_h,hamiltonian_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    call MPI_Allreduce(local_o,observable_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(hamiltonian_fingerprint==0_int64)hamiltonian_fingerprint=1_int64
    if(observable_fingerprint==0_int64)observable_fingerprint=1_int64
    ok=.true.
#else
    hamiltonian_fingerprint=0_int64;observable_fingerprint=0_int64;ok=.false.
#endif
  end subroutine

  integer(int64) function matrix_entry_fingerprint(value,row,column,axis,kind)
    complex(8),intent(in)::value
    integer(int64),intent(in)::row
    integer,intent(in)::column,axis,kind
    integer(int64)::real_bits,imaginary_bits,key
    real_bits=transfer(real(value,8),real_bits);imaginary_bits=transfer(aimag(value),imaginary_bits)
    key=ieor(row,ishft(int(column,int64),17))
    key=ieor(key,ishft(int(axis,int64),37));key=ieor(key,ishft(int(kind,int64),43))
    matrix_entry_fingerprint=ieor(ishftc(real_bits,modulo(int(row),63)),&
      ishftc(imaginary_bits,modulo(column+7*axis+11*kind,63)))
    matrix_entry_fingerprint=ieor(matrix_entry_fingerprint,ishftc(key,23))
  end function

#ifdef USE_MPI
  integer(int64) function expected_manifest_size(checkpoint)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer(int64)::integer_bytes,integer8_bytes,real_bytes,complex_bytes,logical_bytes
    integer_bytes=storage_size(0)/8;integer8_bytes=storage_size(0_int64)/8
    real_bytes=storage_size(0d0)/8;complex_bytes=storage_size((0d0,0d0))/8
    logical_bytes=storage_size(.false.)/8
    expected_manifest_size=int(storage_size(manifest_magic)/8,int64)+&
      integer_bytes*int(7+size(checkpoint%center_owner),int64)+integer8_bytes*7_int64+&
      int(storage_size(checkpoint%field_coupling_convention)/8,int64)+&
      real_bytes*int(13+size(checkpoint%occupations),int64)+logical_bytes+&
      complex_bytes*int(size(checkpoint%coefficients),int64)
  end function

  integer(int64) function expected_shard_size(checkpoint)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer(int64)::integer_bytes,integer8_bytes,real_bytes,complex_bytes
    integer_bytes=storage_size(0)/8;integer8_bytes=storage_size(0_int64)/8;real_bytes=storage_size(0d0)/8
    complex_bytes=storage_size((0d0,0d0))/8
    expected_shard_size=int(storage_size(shard_magic)/8,int64)+integer_bytes*&
      int(9+3*size(checkpoint%tail_center),int64)+integer8_bytes*&
      int(5+size(checkpoint%tail_physical_ids)+size(checkpoint%core_physical_ids)+&
      size(checkpoint%overlap_row_ids),int64)+real_bytes*int(size(checkpoint%density),int64)+&
      complex_bytes*int(size(checkpoint%overlap)+size(checkpoint%hamiltonian0)+&
      size(checkpoint%position)+size(checkpoint%velocity),int64)
  end function

  integer(int64) function digest_manifest(checkpoint,transaction_id)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer(int64),intent(in)::transaction_id
    integer::i,j,position
    integer(int64)::hash
    hash=int(z'BB67AE8584CAA73B',int64);position=1
    call digest_integer(hash,int(checkpoint%basis_generation,int64),position)
    call digest_integer(hash,int(checkpoint%geometry_generation,int64),position)
    call digest_integer(hash,checkpoint%basis_fingerprint,position)
    call digest_integer(hash,checkpoint%operator_fingerprint,position)
    call digest_integer(hash,checkpoint%hamiltonian_fingerprint,position)
    call digest_integer(hash,checkpoint%observable_fingerprint,position)
    do i=1,len(checkpoint%field_coupling_convention)
      call digest_integer(hash,int(iachar(checkpoint%field_coupling_convention(i:i)),int64),position)
    enddo
    call digest_integer(hash,transaction_id,position)
    do i=1,size(checkpoint%center_owner);call digest_integer(hash,int(checkpoint%center_owner(i),int64),position);enddo
    do j=1,size(checkpoint%coefficients,2);do i=1,size(checkpoint%coefficients,1)
      call digest_complex(hash,checkpoint%coefficients(i,j),position)
    enddo;enddo
    do i=1,size(checkpoint%occupations);call digest_real(hash,checkpoint%occupations(i),position);enddo
    call digest_real(hash,checkpoint%density_residual,position)
    call digest_real(hash,checkpoint%coefficient_residual,position)
    call digest_real(hash,checkpoint%charge_error,position)
    call digest_real(hash,checkpoint%unmixed_density_residual,position)
    call digest_real(hash,checkpoint%orthogonality_defect,position)
    call digest_real(hash,checkpoint%metric_condition,position)
    call digest_real(hash,checkpoint%density_tolerance,position)
    call digest_real(hash,checkpoint%coefficient_tolerance,position)
    call digest_real(hash,checkpoint%orthogonality_tolerance,position)
    call digest_real(hash,checkpoint%charge_tolerance,position)
    call digest_real(hash,checkpoint%condition_limit,position)
    call digest_real(hash,checkpoint%symmetry_closure_residual,position)
    call digest_real(hash,checkpoint%symmetry_tolerance,position)
    call digest_integer(hash,int(merge(1,0,checkpoint%accepted),int64),position)
    digest_manifest=hash
  end function

  integer(int64) function digest_shard(checkpoint,rank,transaction_id)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer,intent(in)::rank
    integer(int64),intent(in)::transaction_id
    integer::i,j,position
    integer(int64)::hash
    hash=int(z'3C6EF372FE94F82B',int64);position=1
    call digest_integer(hash,int(rank,int64),position)
    call digest_integer(hash,int(checkpoint%basis_generation,int64),position)
    call digest_integer(hash,int(checkpoint%geometry_generation,int64),position)
    call digest_integer(hash,checkpoint%basis_fingerprint,position)
    call digest_integer(hash,checkpoint%operator_fingerprint,position)
    call digest_integer(hash,transaction_id,position)
    do i=1,size(checkpoint%tail_center);call digest_integer(hash,int(checkpoint%tail_center(i),int64),position);enddo
    do i=1,size(checkpoint%tail_generation);call digest_integer(hash,int(checkpoint%tail_generation(i),int64),position);enddo
    do i=1,size(checkpoint%tail_offsets);call digest_integer(hash,int(checkpoint%tail_offsets(i),int64),position);enddo
    do i=1,size(checkpoint%tail_physical_ids);call digest_integer(hash,checkpoint%tail_physical_ids(i),position);enddo
    do i=1,size(checkpoint%core_physical_ids);call digest_integer(hash,checkpoint%core_physical_ids(i),position);enddo
    do i=1,size(checkpoint%density);call digest_real(hash,checkpoint%density(i),position);enddo
    do i=1,size(checkpoint%overlap_row_ids)
      call digest_integer(hash,checkpoint%overlap_row_ids(i),position)
    enddo
    do j=1,size(checkpoint%overlap,2);do i=1,size(checkpoint%overlap,1)
      call digest_complex(hash,checkpoint%overlap(i,j),position)
    enddo;enddo
    do j=1,size(checkpoint%hamiltonian0,2);do i=1,size(checkpoint%hamiltonian0,1)
      call digest_complex(hash,checkpoint%hamiltonian0(i,j),position)
    enddo;enddo
    do j=1,size(checkpoint%position,3);do i=1,size(checkpoint%position,2)
      call digest_complex(hash,checkpoint%position(1,i,j),position)
      call digest_complex(hash,checkpoint%position(2,i,j),position)
      call digest_complex(hash,checkpoint%position(3,i,j),position)
      call digest_complex(hash,checkpoint%velocity(1,i,j),position)
      call digest_complex(hash,checkpoint%velocity(2,i,j),position)
      call digest_complex(hash,checkpoint%velocity(3,i,j),position)
    enddo;enddo
    digest_shard=hash
  end function

  subroutine digest_serialized_shard(unit,payload_position,checkpoint,rank,transaction_id,hash,ios,iomsg)
    integer,intent(in)::unit,rank
    integer(int64),intent(in)::payload_position,transaction_id
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer(int64),intent(out)::hash
    integer,intent(out)::ios
    character(*),intent(out)::iomsg
    integer,parameter::digest_block_size=4096
    integer::i,j,axis,digest_index,integer_value,first,nblock,k,total
    integer(int64)::integer8_value,position_start,velocity_start,complex_bytes,integer_bytes,integer8_bytes,&
      real_bytes,stream_position,element_offset
    real(8)::real_value
    complex(8)::complex_value
    integer::integer_block(digest_block_size)
    integer(int64)::integer8_block(digest_block_size)
    real(8)::real_block(digest_block_size)
    complex(8)::complex_block(digest_block_size),position_block(3,digest_block_size),&
      velocity_block(3,digest_block_size)
    hash=int(z'3C6EF372FE94F82B',int64);digest_index=1;ios=0;iomsg=''
    integer_bytes=int(storage_size(integer_value)/8,int64)
    integer8_bytes=int(storage_size(integer8_value)/8,int64)
    real_bytes=int(storage_size(real_value)/8,int64)
    complex_bytes=int(storage_size(complex_value)/8,int64)
    stream_position=payload_position
    call digest_integer(hash,int(rank,int64),digest_index)
    call digest_integer(hash,int(checkpoint%basis_generation,int64),digest_index)
    call digest_integer(hash,int(checkpoint%geometry_generation,int64),digest_index)
    call digest_integer(hash,checkpoint%basis_fingerprint,digest_index)
    call digest_integer(hash,checkpoint%operator_fingerprint,digest_index)
    call digest_integer(hash,transaction_id,digest_index)
    total=size(checkpoint%tail_center)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,int(integer_block(k),int64),digest_index)
      enddo
    enddo
    total=size(checkpoint%tail_generation)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,int(integer_block(k),int64),digest_index)
      enddo
    enddo
    total=size(checkpoint%tail_offsets)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,int(integer_block(k),int64),digest_index)
      enddo
    enddo
    total=size(checkpoint%tail_physical_ids)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer8_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer8_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,integer8_block(k),digest_index)
      enddo
    enddo
    total=size(checkpoint%core_physical_ids)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer8_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer8_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,integer8_block(k),digest_index)
      enddo
    enddo
    total=size(checkpoint%density)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)real_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+real_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_real(hash,real_block(k),digest_index)
      enddo
    enddo
    total=size(checkpoint%overlap_row_ids)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)integer8_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+integer8_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_integer(hash,integer8_block(k),digest_index)
      enddo
    enddo
    total=size(checkpoint%overlap)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)complex_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+complex_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_complex(hash,complex_block(k),digest_index)
      enddo
    enddo
    total=size(checkpoint%hamiltonian0)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      read(unit,pos=stream_position,iostat=ios,iomsg=iomsg)complex_block(1:nblock)
      if(ios/=0)return
      stream_position=stream_position+complex_bytes*int(nblock,int64)
      do k=1,nblock
        call digest_complex(hash,complex_block(k),digest_index)
      enddo
    enddo
    position_start=stream_position
    velocity_start=position_start+complex_bytes*int(size(checkpoint%position),int64)
    total=size(checkpoint%position,2)*size(checkpoint%position,3)
    do first=1,total,digest_block_size
      nblock=min(digest_block_size,total-first+1)
      element_offset=3_int64*int(first-1,int64)
      read(unit,pos=position_start+complex_bytes*element_offset,iostat=ios,iomsg=iomsg)&
        position_block(:,1:nblock)
      if(ios/=0)return
      read(unit,pos=velocity_start+complex_bytes*element_offset,iostat=ios,iomsg=iomsg)&
        velocity_block(:,1:nblock)
      if(ios/=0)return
      do k=1,nblock
        do axis=1,3
          call digest_complex(hash,position_block(axis,k),digest_index)
        enddo
        do axis=1,3
          call digest_complex(hash,velocity_block(axis,k),digest_index)
        enddo
      enddo
    enddo
  end subroutine

  subroutine digest_integer(hash,value,position)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    integer,intent(inout)::position
    hash=ieor(ishftc(hash,13),ishftc(value,mod(position,63)));position=position+1
  end subroutine

  subroutine digest_real(hash,value,position)
    integer(int64),intent(inout)::hash
    real(8),intent(in)::value
    integer,intent(inout)::position
    integer(int64)::bits
    bits=transfer(value,bits);call digest_integer(hash,bits,position)
  end subroutine

  subroutine digest_complex(hash,value,position)
    integer(int64),intent(inout)::hash
    complex(8),intent(in)::value
    integer,intent(inout)::position
    call digest_real(hash,real(value,8),position);call digest_real(hash,aimag(value),position)
  end subroutine

  subroutine validate_shard_ownership(comm,checkpoint,bad)
    integer,intent(in)::comm
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer,intent(out)::bad
    integer::nproc,rank,ierr,i,j,total,local_count,expected_tail_count,global_bad,row_count,row_total,&
      tail_first,tail_last
    integer(int64)::total64
    integer,allocatable::counts(:),owners(:),row_owners(:),tail_owners(:)
    logical::tail_offsets_safe
    local_count=size(checkpoint%core_physical_ids)
    call MPI_Comm_size(comm,nproc,ierr);call MPI_Comm_rank(comm,rank,ierr)
    allocate(counts(nproc))
    call MPI_Allgather(local_count,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total64=0_int64;bad=0
    do i=1,nproc
      if(counts(i)<0.or.total64>10000000_int64-int(counts(i),int64))then
        bad=1
      else
        total64=total64+int(counts(i),int64)
      endif
    enddo
    if(total64<1_int64.or.size(checkpoint%density)/=local_count)bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)return
    total=int(total64)
    allocate(owners(max(1,total)));owners=0
    do i=1,local_count
      if(checkpoint%core_physical_ids(i)<1_int64.or.checkpoint%core_physical_ids(i)>int(total,int64))then
        bad=1
      else
        owners(int(checkpoint%core_physical_ids(i)))=owners(int(checkpoint%core_physical_ids(i)))+1
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owners,max(1,total),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(total>0)then
      if(any(owners(1:total)/=1))bad=1
    endif
    expected_tail_count=count(checkpoint%center_owner==rank)
    tail_offsets_safe=size(checkpoint%tail_offsets)==expected_tail_count+1.and.&
      size(checkpoint%tail_offsets)==size(checkpoint%tail_center)+1
    if(size(checkpoint%tail_center)/=expected_tail_count.or.&
       size(checkpoint%tail_generation)/=expected_tail_count.or.&
       size(checkpoint%tail_offsets)/=expected_tail_count+1)bad=1
    if(tail_offsets_safe)then
      tail_offsets_safe=checkpoint%tail_offsets(1)==1.and.&
        checkpoint%tail_offsets(size(checkpoint%tail_offsets))==size(checkpoint%tail_physical_ids)+1.and.&
        all(checkpoint%tail_offsets>=1).and.&
        all(checkpoint%tail_offsets<=size(checkpoint%tail_physical_ids)+1).and.&
        all(checkpoint%tail_offsets(2:)>=checkpoint%tail_offsets(:size(checkpoint%tail_offsets)-1))
      if(.not.tail_offsets_safe)bad=1
    endif
    do i=1,size(checkpoint%tail_center)
      if(checkpoint%tail_center(i)<1.or.checkpoint%tail_center(i)>size(checkpoint%center_owner))then
        bad=1
      else if(checkpoint%center_owner(checkpoint%tail_center(i))/=rank)then
        bad=1
      endif
      if(count(checkpoint%tail_center==checkpoint%tail_center(i))/=1)bad=1
      if(.not.tail_offsets_safe)cycle
      tail_first=checkpoint%tail_offsets(i);tail_last=checkpoint%tail_offsets(i+1)-1
      if(tail_first>tail_last.or.tail_last-tail_first+1<total)then
        bad=1
      else
        allocate(tail_owners(total));tail_owners=0
        do j=tail_first,tail_last
          if(checkpoint%tail_physical_ids(j)<1_int64.or.&
             checkpoint%tail_physical_ids(j)>int(total,int64))then
            bad=1
          else
            tail_owners(int(checkpoint%tail_physical_ids(j)))=&
              tail_owners(int(checkpoint%tail_physical_ids(j)))+1
          endif
        enddo
        if(any(tail_owners<1))bad=1
        deallocate(tail_owners)
      endif
    enddo
    if(any(checkpoint%tail_generation/=checkpoint%basis_generation).or.&
       any(checkpoint%tail_physical_ids<=0_int64))bad=1
    row_count=size(checkpoint%overlap_row_ids)
    if(row_count/=expected_tail_count.or.size(checkpoint%overlap,1)/=row_count.or.&
        size(checkpoint%overlap,2)/=size(checkpoint%center_owner))bad=1
    if(any(shape(checkpoint%hamiltonian0)/=[row_count,size(checkpoint%center_owner)]).or.&
       any(shape(checkpoint%position)/=[3,row_count,size(checkpoint%center_owner)]).or.&
       any(shape(checkpoint%velocity)/=[3,row_count,size(checkpoint%center_owner)]))bad=1
    call MPI_Allreduce(row_count,row_total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(row_total/=size(checkpoint%center_owner))bad=1
    allocate(row_owners(max(1,size(checkpoint%center_owner))));row_owners=0
    do i=1,row_count
      if(checkpoint%overlap_row_ids(i)<1_int64.or.&
          checkpoint%overlap_row_ids(i)>int(size(checkpoint%center_owner),int64))then
        bad=1
      else
        row_owners(int(checkpoint%overlap_row_ids(i)))=&
          row_owners(int(checkpoint%overlap_row_ids(i)))+1
        if(checkpoint%center_owner(int(checkpoint%overlap_row_ids(i)))/=rank)bad=1
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_owners,max(1,size(checkpoint%center_owner)),&
      MPI_INTEGER,MPI_SUM,comm,ierr)
    if(size(checkpoint%center_owner)>0)then
      if(any(row_owners(1:size(checkpoint%center_owner))/=1))bad=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine

  subroutine validate_replicated_payload(comm,checkpoint,bad)
    integer,intent(in)::comm
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer,intent(out)::bad
    integer::ierr,dims(4),dims_min(4),dims_max(4),integer_values(3),integer_min(3),integer_max(3)
    integer(int64)::fingerprints(4),fingerprints_min(4),fingerprints_max(4)
    character(64)::convention_reference
    integer,allocatable::owner_reference(:)
    complex(8),allocatable::coefficients_reference(:,:)
    real(8),allocatable::occupations_reference(:)
    real(8)::local_defect,global_defect,metrics(13),metrics_min(13),metrics_max(13)
    dims=[size(checkpoint%center_owner),size(checkpoint%coefficients,1),size(checkpoint%coefficients,2),&
      size(checkpoint%occupations)]
    integer_values=[checkpoint%basis_generation,checkpoint%geometry_generation,merge(1,0,checkpoint%accepted)]
    fingerprints=[checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,&
      checkpoint%hamiltonian_fingerprint,checkpoint%observable_fingerprint]
    metrics=[checkpoint%density_residual,checkpoint%coefficient_residual,checkpoint%charge_error,&
      checkpoint%unmixed_density_residual,checkpoint%orthogonality_defect,checkpoint%metric_condition,&
      checkpoint%density_tolerance,checkpoint%coefficient_tolerance,checkpoint%orthogonality_tolerance,&
      checkpoint%charge_tolerance,checkpoint%condition_limit,checkpoint%symmetry_closure_residual,&
      checkpoint%symmetry_tolerance]
    call MPI_Allreduce(dims,dims_min,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(dims,dims_max,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(integer_values,integer_min,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(integer_values,integer_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(fingerprints,fingerprints_min,4,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(fingerprints,fingerprints_max,4,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(metrics,metrics_min,13,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(metrics,metrics_max,13,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    bad=0
    if(any(dims_min/=dims_max).or.any(integer_min/=integer_max).or.any(fingerprints_min/=fingerprints_max).or.&
       any(metrics_min/=metrics_max))bad=1
    if(bad/=0)return
    convention_reference=checkpoint%field_coupling_convention
    call MPI_Bcast(convention_reference,len(convention_reference),MPI_CHARACTER,0,comm,ierr)
    if(checkpoint%field_coupling_convention/=convention_reference)then;bad=1;return;endif
    allocate(owner_reference(dims(1)),coefficients_reference(dims(2),dims(3)),&
      occupations_reference(dims(4)))
    owner_reference=checkpoint%center_owner
    coefficients_reference=checkpoint%coefficients;occupations_reference=checkpoint%occupations
    call MPI_Bcast(owner_reference,dims(1),MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(coefficients_reference,dims(2)*dims(3),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(occupations_reference,dims(4),MPI_DOUBLE_PRECISION,0,comm,ierr)
    local_defect=max(maxval(abs(checkpoint%coefficients-coefficients_reference)),&
      maxval(abs(checkpoint%occupations-occupations_reference)))
    if(any(checkpoint%center_owner/=owner_reference))local_defect=huge(1d0)
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    bad=merge(0,1,global_defect==0d0)
  end subroutine

  subroutine validate_checkpoint(checkpoint,nproc,bad)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer,intent(in)::nproc
    integer,intent(out)::bad
    integer(int64)::matrix_count
    bad=0
    if(.not.checkpoint%accepted.or.checkpoint%basis_generation<=0.or.checkpoint%geometry_generation<=0)bad=ibset(bad,0)
    if(checkpoint%basis_fingerprint==0_int64.or.checkpoint%operator_fingerprint==0_int64.or.&
       checkpoint%hamiltonian_fingerprint==0_int64.or.checkpoint%observable_fingerprint==0_int64.or.&
       trim(checkpoint%field_coupling_convention)/='cell_wrapped_length_velocity')bad=ibset(bad,1)
    if(.not.allocated(checkpoint%center_owner).or..not.allocated(checkpoint%tail_center).or.&
       .not.allocated(checkpoint%tail_generation).or..not.allocated(checkpoint%tail_offsets).or.&
       .not.allocated(checkpoint%tail_physical_ids).or..not.allocated(checkpoint%core_physical_ids).or.&
       .not.allocated(checkpoint%overlap_row_ids).or..not.allocated(checkpoint%overlap).or.&
       .not.allocated(checkpoint%hamiltonian0).or..not.allocated(checkpoint%position).or.&
       .not.allocated(checkpoint%velocity).or.&
       .not.allocated(checkpoint%coefficients).or.&
       .not.allocated(checkpoint%occupations).or..not.allocated(checkpoint%density))then
      bad=ibset(bad,2);return
    endif
    if(size(checkpoint%center_owner)<1.or.any(checkpoint%center_owner<0).or.any(checkpoint%center_owner>=nproc))&
      bad=ibset(bad,3)
    if(size(checkpoint%center_owner)/=size(checkpoint%overlap,2).or.&
       size(checkpoint%overlap_row_ids)/=size(checkpoint%overlap,1))bad=ibset(bad,4)
    if(any(shape(checkpoint%hamiltonian0)/=shape(checkpoint%overlap)).or.&
       any(shape(checkpoint%position)/=[3,size(checkpoint%overlap,1),size(checkpoint%overlap,2)]).or.&
       any(shape(checkpoint%velocity)/=[3,size(checkpoint%overlap,1),size(checkpoint%overlap,2)]))bad=ibset(bad,4)
    if(size(checkpoint%coefficients,1)/=size(checkpoint%center_owner).or.&
       size(checkpoint%coefficients,2)/=size(checkpoint%occupations))bad=ibset(bad,5)
    matrix_count=8_int64*int(size(checkpoint%overlap,1),int64)*int(size(checkpoint%overlap,2),int64)+&
      int(size(checkpoint%coefficients,1),int64)*int(size(checkpoint%coefficients,2),int64)
    if(matrix_count>100000000_int64)bad=ibset(bad,6)
    if(any(checkpoint%core_physical_ids<=0_int64).or.any(checkpoint%density<0d0).or.any(checkpoint%occupations<0d0))&
      bad=ibset(bad,7)
    if(size(checkpoint%core_physical_ids)>10000000.or.size(checkpoint%tail_physical_ids)>10000000.or.&
       size(checkpoint%tail_center)>1000000)bad=ibset(bad,8)
    if(.not.all(ieee_is_finite(checkpoint%density)).or..not.all(ieee_is_finite(checkpoint%occupations)).or.&
       .not.all(ieee_is_finite([checkpoint%density_residual,checkpoint%coefficient_residual,&
       checkpoint%charge_error,checkpoint%unmixed_density_residual,checkpoint%orthogonality_defect,&
       checkpoint%metric_condition,checkpoint%density_tolerance,checkpoint%coefficient_tolerance,&
       checkpoint%orthogonality_tolerance,checkpoint%charge_tolerance,checkpoint%condition_limit])))bad=ibset(bad,9)
    if(.not.all(ieee_is_finite([checkpoint%symmetry_closure_residual,checkpoint%symmetry_tolerance])))&
      bad=ibset(bad,10)
    if(min(checkpoint%density_residual,checkpoint%coefficient_residual,checkpoint%unmixed_density_residual,&
       checkpoint%orthogonality_defect)<0d0.or.checkpoint%metric_condition<1d0)bad=ibset(bad,11)
    if(min(checkpoint%density_tolerance,checkpoint%coefficient_tolerance,checkpoint%orthogonality_tolerance,&
       checkpoint%charge_tolerance,checkpoint%condition_limit)<=0d0)bad=ibset(bad,12)
    if(checkpoint%density_residual>checkpoint%density_tolerance.or.&
       checkpoint%unmixed_density_residual>checkpoint%density_tolerance.or.&
       checkpoint%coefficient_residual>checkpoint%coefficient_tolerance.or.&
       checkpoint%orthogonality_defect>checkpoint%orthogonality_tolerance.or.&
       abs(checkpoint%charge_error)>checkpoint%charge_tolerance.or.&
       checkpoint%metric_condition>checkpoint%condition_limit)bad=ibset(bad,13)
    if(checkpoint%symmetry_tolerance<=0d0.or.checkpoint%symmetry_closure_residual<0d0.or.&
       checkpoint%symmetry_closure_residual>checkpoint%symmetry_tolerance)bad=ibset(bad,14)
    if(.not.finite_complex(checkpoint%overlap).or..not.finite_complex(checkpoint%hamiltonian0).or.&
       .not.finite_complex3(checkpoint%position).or..not.finite_complex3(checkpoint%velocity).or.&
       .not.finite_complex(checkpoint%coefficients))bad=ibset(bad,15)
  end subroutine

  integer function checkpoint_rejection_code(checkpoint,nproc) result(code)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    integer,intent(in)::nproc
    code=0
    if(.not.checkpoint%accepted)code=ibset(code,0)
    if(checkpoint%basis_generation<=0.or.checkpoint%geometry_generation<=0)code=ibset(code,1)
    if(checkpoint%basis_fingerprint==0_int64.or.checkpoint%operator_fingerprint==0_int64.or.&
       checkpoint%hamiltonian_fingerprint==0_int64.or.checkpoint%observable_fingerprint==0_int64.or.&
       trim(checkpoint%field_coupling_convention)/='cell_wrapped_length_velocity')code=ibset(code,2)
    if(.not.allocated(checkpoint%center_owner).or..not.allocated(checkpoint%tail_center).or.&
       .not.allocated(checkpoint%tail_generation).or..not.allocated(checkpoint%tail_offsets).or.&
       .not.allocated(checkpoint%tail_physical_ids).or..not.allocated(checkpoint%core_physical_ids).or.&
       .not.allocated(checkpoint%overlap_row_ids).or..not.allocated(checkpoint%overlap).or.&
       .not.allocated(checkpoint%hamiltonian0).or..not.allocated(checkpoint%position).or.&
       .not.allocated(checkpoint%velocity).or.&
       .not.allocated(checkpoint%coefficients).or.&
       .not.allocated(checkpoint%occupations).or..not.allocated(checkpoint%density))then
      code=ibset(code,3);return
    endif
    if(.not.ieee_is_finite(checkpoint%density_residual).or.&
       checkpoint%density_residual<0d0.or.checkpoint%density_residual>checkpoint%density_tolerance)&
      code=ibset(code,4)
    if(.not.ieee_is_finite(checkpoint%unmixed_density_residual).or.&
       checkpoint%unmixed_density_residual<0d0.or.&
       checkpoint%unmixed_density_residual>checkpoint%density_tolerance)code=ibset(code,5)
    if(.not.ieee_is_finite(checkpoint%charge_error).or.&
       abs(checkpoint%charge_error)>checkpoint%charge_tolerance)code=ibset(code,6)
    if(.not.ieee_is_finite(checkpoint%coefficient_residual).or.&
       checkpoint%coefficient_residual<0d0.or.&
       checkpoint%coefficient_residual>checkpoint%coefficient_tolerance)code=ibset(code,7)
    if(.not.ieee_is_finite(checkpoint%orthogonality_defect).or.&
       checkpoint%orthogonality_defect<0d0.or.&
       checkpoint%orthogonality_defect>checkpoint%orthogonality_tolerance)code=ibset(code,8)
    if(.not.ieee_is_finite(checkpoint%metric_condition).or.checkpoint%metric_condition<1d0.or.&
       checkpoint%metric_condition>checkpoint%condition_limit)code=ibset(code,9)
    if(.not.ieee_is_finite(checkpoint%symmetry_closure_residual).or.&
       checkpoint%symmetry_closure_residual<0d0.or.&
       checkpoint%symmetry_closure_residual>checkpoint%symmetry_tolerance)code=ibset(code,10)
    if(size(checkpoint%center_owner)<1.or.any(checkpoint%center_owner<0).or.&
       any(checkpoint%center_owner>=nproc).or.any(checkpoint%core_physical_ids<=0_int64).or.&
       any(checkpoint%density<0d0).or.any(checkpoint%occupations<0d0).or.&
       .not.finite_complex(checkpoint%overlap).or..not.finite_complex(checkpoint%hamiltonian0).or.&
       .not.finite_complex3(checkpoint%position).or..not.finite_complex3(checkpoint%velocity).or.&
       .not.finite_complex(checkpoint%coefficients))&
      code=ibset(code,11)
  end function

  logical function finite_complex(values)
    complex(8),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function

  logical function finite_complex3(values)
    complex(8),intent(in)::values(:,:,:)
    finite_complex3=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function

  subroutine versioned_shard_name(prefix,generation,fingerprint,transaction_id,rank,name)
    character(*),intent(in)::prefix
    integer,intent(in)::generation,rank
    integer(int64),intent(in)::fingerprint,transaction_id
    character(*),intent(out)::name
    write(name,'(a,".g",i0,"-",z16.16,".t",z16.16,".rank",i8.8)')trim(prefix),generation,&
      fingerprint,transaction_id,rank
  end subroutine

  subroutine choose_transaction_id(prefix,generation,fingerprint,transaction_id,reservation,ios)
    character(*),intent(in)::prefix
    integer,intent(in)::generation
    integer(int64),intent(in)::fingerprint
    integer(int64),intent(out)::transaction_id
    character(*),intent(out)::reservation
    integer,intent(out)::ios
    integer::rank,unit,values(8),attempt,close_ios
    logical::exists
    call system_clock(count=transaction_id)
    call date_and_time(values=values)
    transaction_id=ieor(transaction_id,ishft(int(values(8),int64),17))
    transaction_id=ieor(transaction_id,ishft(int(values(7),int64),33))
    transaction_id=ieor(transaction_id,ishft(int(values(6),int64),49))
    if(transaction_id==0_int64)transaction_id=1_int64
    rank=0
    ios=1
    do attempt=1,1024
      write(reservation,'(a,".transaction.",z16.16,".reservation")')trim(prefix),transaction_id
      open(newunit=unit,file=trim(reservation),status='new',action='write',iostat=ios)
      if(ios==0)then
        write(unit,'(z16.16)',iostat=ios)transaction_id
        close_ios=0;close(unit,iostat=close_ios)
        if(ios==0)ios=close_ios
        if(ios==0)exit
        call remove_file_if_present(reservation)
      else
        inquire(file=trim(reservation),exist=exists)
        if(.not.exists)return
      endif
      transaction_id=transaction_id+1_int64
      if(transaction_id==0_int64)transaction_id=1_int64
    enddo
  end subroutine

  subroutine abort_publication(rank,manifest_temporary,shard_temporary,shard,reservation,checkpoint,previous_id)
    integer,intent(in)::rank
    character(*),intent(in)::manifest_temporary,shard_temporary,shard,reservation
    type(s_dg_overlapping_wannier_checkpoint),intent(inout)::checkpoint
    integer(int64),intent(in)::previous_id
    call remove_file_if_present(shard_temporary)
    call remove_file_if_present(shard)
    if(rank==0)then
      call remove_file_if_present(manifest_temporary)
      call remove_file_if_present(reservation)
    endif
    checkpoint%publication_id=previous_id
  end subroutine

  subroutine remove_file_if_present(filename)
    character(*),intent(in)::filename
    integer::unit,ios
    open(newunit=unit,file=trim(filename),status='old',iostat=ios)
    if(ios==0)close(unit,status='delete')
  end subroutine
#endif
end module
