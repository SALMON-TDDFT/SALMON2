#include "config.h"
module rt_dg_hybrid_checkpoint
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_hybrid_checkpoint_v4,only:s_rt_dg_hybrid_v4_shard,write_rt_dg_hybrid_checkpoint_v4
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  character(16),parameter::occupied_magic="SALMON_DG_OCC02 "
  integer,parameter::occupied_version=2
  integer,parameter,public::rt_dg_hybrid_occupied_checkpoint_version=2
  public::write_rt_dg_hybrid_occupied_checkpoint,read_rt_dg_hybrid_occupied_checkpoint,&
    collective_rt_dg_hybrid_publication_precondition,&
    collective_rt_dg_hybrid_publication_mapping_precondition,&
    publish_rt_dg_hybrid_checkpoint_v4
  interface
    function c_rename(old_path,new_path) bind(C,name="rename") result(status)
      import::c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
      integer(c_int)::status
    end function c_rename
  end interface
contains
  subroutine collective_rt_dg_hybrid_publication_precondition(comm,local_valid,local_n,local_nocc,ok,message)
    integer,intent(in)::comm,local_n,local_nocc
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_bad,global_bad,local_signature(2),minimum_signature(2),maximum_signature(2)
    local_bad=merge(0,1,local_valid);local_signature=[local_n,local_nocc]
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v4 publication validity reduction failed';return;endif
    call MPI_Allreduce(local_signature,minimum_signature,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v4 publication minimum reduction failed';return;endif
    call MPI_Allreduce(local_signature,maximum_signature,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='terminal divided v4 publication maximum reduction failed';return;endif
    ok=global_bad==0.and.all(minimum_signature==maximum_signature)
    if(ok)then;message='';else;message='terminal divided v4 publication collective precondition failed';endif
#else
    ok=.false.;message='terminal divided v4 publication precondition requires MPI'
#endif
  end subroutine collective_rt_dg_hybrid_publication_precondition

  subroutine collective_rt_dg_hybrid_publication_mapping_precondition(comm,global_count,row_ids,row_owner,&
      occupied_row_ids,local_valid,ok,message)
    integer,intent(in)::comm,global_count,row_owner(:)
    integer(int64),intent(in)::row_ids(:),occupied_row_ids(:)
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,i,ierr,local_bad,global_bad,allocation_status
    integer,allocatable::row_counts(:),owner_minimum(:),owner_maximum(:)
    ok=.false.;message='terminal divided v4 publication stage=pre-gather global row mapping failed'
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(local_valid)then;local_bad=0;else;local_bad=1;endif
    if(global_count<1.or.size(row_owner)/=global_count.or.size(occupied_row_ids)/=size(row_ids))local_bad=1
    if(local_bad==0)then
      if(any(row_owner<0).or.any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or.&
        any(occupied_row_ids<1_int64).or.any(occupied_row_ids>int(global_count,int64)))local_bad=1
    endif
    if(local_bad==0.and.size(row_ids)>0)then
      if(any(row_owner(int(row_ids))/=rank).or.any(occupied_row_ids/=row_ids))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)return
    allocate(row_counts(global_count),owner_minimum(global_count),&
      owner_maximum(global_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)return
    row_counts=0
    do i=1,size(row_ids);row_counts(int(row_ids(i)))=row_counts(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(row_owner,owner_minimum,global_count,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(row_owner,owner_maximum,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,row_counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    ok=global_bad==0.and.all(row_counts==1).and.all(owner_minimum==owner_maximum)
    if(ok)message=''
#else
    ok=.false.;message='terminal divided v4 publication mapping precondition requires MPI'
#endif
  end subroutine collective_rt_dg_hybrid_publication_mapping_precondition

  subroutine publish_rt_dg_hybrid_checkpoint_v4(comm,path,global_count,noccupied,row_ids,row_owner,&
      occupied_row_ids,payload,local_valid,ok,message)
    integer,intent(in)::comm,global_count,noccupied,row_owner(:)
    character(*),intent(in)::path
    integer(int64),intent(in)::row_ids(:),occupied_row_ids(:)
    type(s_rt_dg_hybrid_v4_shard),intent(in)::payload
    logical,intent(in)::local_valid
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(512)::detail

    call collective_rt_dg_hybrid_publication_precondition(comm,local_valid,global_count,noccupied,ok,detail)
    if(.not.ok)then;message='distributed-v4 endpoint precondition failed: '//trim(detail);return;endif
    call collective_rt_dg_hybrid_publication_mapping_precondition(comm,global_count,row_ids,row_owner,&
      occupied_row_ids,local_valid,ok,detail)
    if(.not.ok)then;message='distributed-v4 endpoint row mapping failed: '//trim(detail);return;endif
    call write_rt_dg_hybrid_checkpoint_v4(comm,path,payload,ok,detail)
    if(.not.ok)then;message='distributed-v4 endpoint publication failed: '//trim(detail);return;endif
    message=''
  end subroutine publish_rt_dg_hybrid_checkpoint_v4

  subroutine write_rt_dg_hybrid_occupied_checkpoint(comm,path,global_count,row_ids,coefficients,occupations,eigenvalues,&
      catalog_fingerprint,basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,scf_receipts,&
      maximum_scf_residual,payload_fingerprint,ok,message)
    integer,intent(in)::comm,global_count
    character(*),intent(in)::path
    integer(int64),intent(in)::row_ids(:),catalog_fingerprint,basis_fingerprint,provenance_fingerprints(6),&
      operator_fingerprint,state_fingerprint
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),eigenvalues(:),scf_receipts(5),maximum_scf_residual
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,noccupied,nlocal,i,j,row,root,position,unit,io_status,close_status,bad,global_bad,status
    integer,allocatable::counts(:),owners(:),positions(:)
    complex(real64),allocatable::row_values(:)
    integer(int64)::minimum_i,maximum_i,bits
    real(real64)::minimum_r,maximum_r
    character(16)::probe
    character(:),allocatable::temporary_path
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    nlocal=size(row_ids);noccupied=size(coefficients,2);bad=0
    if(global_count<1.or.noccupied<1.or.size(coefficients,1)/=nlocal.or.size(occupations)/=noccupied.or.&
      size(eigenvalues)/=noccupied)bad=1
    call agree_int(global_count);call agree_int(noccupied)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint dimensions';return;endif
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or..not.finite_matrix(coefficients).or.&
      any(.not.ieee_is_finite(occupations)).or.any(.not.ieee_is_finite(eigenvalues)).or.&
      any(.not.ieee_is_finite(scf_receipts)).or.any(occupations<0d0))bad=1
    if(catalog_fingerprint==0_int64.or.basis_fingerprint==0_int64.or.any(provenance_fingerprints==0_int64).or.&
      operator_fingerprint==0_int64.or.state_fingerprint==0_int64)bad=1
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0.or.any(scf_receipts<0d0).or.&
      any(scf_receipts>maximum_scf_residual))bad=1
    call agree_i64(catalog_fingerprint);call agree_i64(basis_fingerprint);call agree_i64(operator_fingerprint);call agree_i64(state_fingerprint)
    do i=1,6;call agree_i64(provenance_fingerprints(i));enddo
    call agree_real(maximum_scf_residual)
    do i=1,noccupied;call agree_real(occupations(i));call agree_real(eigenvalues(i));enddo
    do i=1,5;call agree_real(scf_receipts(i));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint write contract';return;endif
    allocate(counts(global_count),owners(global_count),positions(global_count),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='cannot allocate occupied checkpoint workspace';return;endif
    counts=0;owners=-1;positions=0
    do i=1,nlocal;row=int(row_ids(i));counts(row)=counts(row)+1;owners(row)=rank;positions(row)=i;enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,positions,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(counts==1));call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='occupied checkpoint rows are not exactly once';return;endif
    temporary_path=trim(path)//'.occupied.tmp.'//trim(int64_string(state_fingerprint))
    if(rank==0)then
      open(newunit=unit,file=temporary_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status);opened=io_status==0
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    payload_fingerprint=catalog_fingerprint
    call occupied_hash(basis_fingerprint);do i=1,6;call occupied_hash(provenance_fingerprints(i));enddo
    call occupied_hash(operator_fingerprint);call occupied_hash(state_fingerprint)
    call occupied_hash(int(global_count,int64));call occupied_hash(int(noccupied,int64))
    do i=1,noccupied;call occupied_hash(transfer(occupations(i),bits));call occupied_hash(transfer(eigenvalues(i),bits));enddo
    do i=1,5;call occupied_hash(transfer(scf_receipts(i),bits));enddo
    if(rank==0)write(unit,iostat=io_status)occupied_magic,occupied_version,global_count,noccupied,catalog_fingerprint,&
      basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,occupations,eigenvalues,scf_receipts
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do row=1,global_count
      root=owners(row);position=positions(row);if(rank==root)row_values=coefficients(position,:)
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call occupied_hash(int(row,int64));do j=1,noccupied;call occupied_hash(transfer(real(row_values(j)),bits));call occupied_hash(transfer(aimag(row_values(j)),bits));enddo
      if(rank==0)write(unit,iostat=io_status)row_values
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)then
      write(unit,iostat=io_status)payload_fingerprint;close_status=0;close(unit,iostat=close_status);opened=.false.;if(io_status==0)io_status=close_status
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call clean;ok=.true.;return
900 message='occupied checkpoint MPI stream failed';if(rank==0.and.opened)close(unit);call clean;return
910 message='occupied checkpoint publication failed';if(rank==0.and.opened)close(unit);call clean;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_int(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_i64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,minimum_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minimum_i/=maximum_i)bad=1
    end subroutine
    subroutine agree_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,minimum_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.transfer(minimum_r,bits)/=transfer(maximum_r,bits))bad=1
    end subroutine
    subroutine occupied_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine clean
      if(allocated(counts))deallocate(counts);if(allocated(owners))deallocate(owners)
      if(allocated(positions))deallocate(positions);if(allocated(row_values))deallocate(row_values)
    end subroutine
#endif
  end subroutine write_rt_dg_hybrid_occupied_checkpoint

  subroutine read_rt_dg_hybrid_occupied_checkpoint(comm,path,expected_catalog,expected_basis,expected_operator,&
      expected_state,expected_provenance,expected_occupations,maximum_scf_residual,global_count,row_ids,coefficients,&
      occupations,eigenvalues,scf_receipts,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::expected_catalog,expected_basis,expected_operator,expected_state,expected_provenance(6)
    real(real64),intent(in)::expected_occupations(:),maximum_scf_residual
    integer,intent(out)::global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),allocatable,intent(out)::occupations(:),eigenvalues(:)
    real(real64),intent(out)::scf_receipts(5)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,bad,global_bad,noccupied,version,nlocal,row,position,j,status
    integer(int64)::catalog,basis,provenance(6),operator_receipt,state,stored_fingerprint,bits,lo_i,hi_i,file_size
    real(real64)::lo_r,hi_r
    complex(real64),allocatable::row_values(:)
    character(16)::magic,probe
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    bad=0
    call agree_expected(expected_catalog);call agree_expected(expected_basis);call agree_expected(expected_operator);call agree_expected(expected_state)
    do j=1,6;call agree_expected(expected_provenance(j));enddo
    call agree_expected_count(size(expected_occupations))
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    call agree_expected_real(maximum_scf_residual)
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0)bad=1
    do j=1,size(expected_occupations);call agree_expected_real(expected_occupations(j));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status);opened=io_status==0
      if(io_status==0)inquire(unit=unit,size=file_size,iostat=io_status)
      if(io_status==0)read(unit,iostat=io_status)magic,version,global_count,noccupied,catalog,basis,provenance,operator_receipt,state
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(magic,len(magic),MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(version,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(global_count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(noccupied,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(catalog,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(basis,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(operator_receipt,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(state,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(provenance,6,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=0
    if(magic/=occupied_magic.or.version/=occupied_version.or.global_count<1.or.noccupied<1)bad=1
    if(catalog/=expected_catalog.or.basis/=expected_basis.or.any(provenance/=expected_provenance).or.&
      operator_receipt/=expected_operator.or.state/=expected_state)bad=1
    if(global_count>huge(0)-nproc.or.noccupied>huge(0)/max(1,global_count))bad=1
    if(file_size<1_int64.or.int(noccupied,int64)>file_size/16_int64.or.&
      int(global_count,int64)>file_size/max(16_int64,16_int64*int(noccupied,int64)))bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='stale or incompatible occupied checkpoint';goto 920;endif
    nlocal=(global_count+nproc-1-rank)/nproc
    allocate(row_ids(nlocal),coefficients(nlocal,noccupied),occupations(noccupied),eigenvalues(noccupied),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate occupied checkpoint output';goto 920;endif
    if(rank==0)read(unit,iostat=io_status)occupations,eigenvalues,scf_receipts
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(occupations,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(eigenvalues,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(scf_receipts,5,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(ieee_is_finite(occupations)).and.all(occupations>=0d0).and.&
      all(ieee_is_finite(eigenvalues)).and.all(ieee_is_finite(scf_receipts)).and.all(scf_receipts>=0d0).and.&
      all(scf_receipts<=maximum_scf_residual).and.size(expected_occupations)==noccupied)
    if(bad==0)then;if(any(occupations/=expected_occupations))bad=1;endif
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    payload_fingerprint=catalog;call read_hash(basis);do j=1,6;call read_hash(provenance(j));enddo
    call read_hash(operator_receipt);call read_hash(state)
    call read_hash(int(global_count,int64));call read_hash(int(noccupied,int64))
    do j=1,noccupied;call read_hash(transfer(occupations(j),bits));call read_hash(transfer(eigenvalues(j),bits));enddo
    do j=1,5;call read_hash(transfer(scf_receipts(j),bits));enddo
    position=0
    do row=1,global_count
      if(rank==0)read(unit,iostat=io_status)row_values
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      if(.not.finite_vector(row_values))goto 930
      call read_hash(int(row,int64));do j=1,noccupied;call read_hash(transfer(real(row_values(j)),bits));call read_hash(transfer(aimag(row_values(j)),bits));enddo
      if(mod(row-1,nproc)==rank)then;position=position+1;row_ids(position)=row;coefficients(position,:)=row_values;endif
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)read(unit,iostat=io_status)stored_fingerprint
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(stored_fingerprint,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    if(rank==0)then;close(unit,iostat=io_status);opened=.false.;endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0.or.stored_fingerprint/=payload_fingerprint)goto 930
    if(allocated(row_values))deallocate(row_values);ok=.true.;return
900 message='occupied checkpoint MPI read failed';goto 920
910 message='occupied checkpoint file read failed';goto 920
930 message='corrupt occupied checkpoint payload'
920 if(rank==0.and.opened)close(unit);call cleanup_output;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0
#endif
  contains
#ifdef USE_MPI
    subroutine agree_expected(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,lo_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo_i/=hi_i.or.value==0_int64)bad=1
    end subroutine
    subroutine agree_expected_count(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_expected_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,lo_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.transfer(lo_r,bits)/=transfer(hi_r,bits))bad=1
    end subroutine
    subroutine read_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine cleanup_output
      if(allocated(row_ids))deallocate(row_ids);if(allocated(coefficients))deallocate(coefficients)
      if(allocated(occupations))deallocate(occupations);if(allocated(eigenvalues))deallocate(eigenvalues)
      if(allocated(row_values))deallocate(row_values);global_count=0;payload_fingerprint=0_int64
    end subroutine
#endif
  end subroutine read_rt_dg_hybrid_occupied_checkpoint

  subroutine validate_path(path,probe,comm,ierr)
    character(*),intent(in)::path;character(16),intent(out)::probe;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    integer::i,status
    integer(int64)::local_hash,minimum_hash,maximum_hash
    local_hash=int(len_trim(path),int64)
    do i=1,len_trim(path);local_hash=ieor(ishftc(local_hash,7),int(iachar(path(i:i)),int64));enddo
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.minimum_hash/=maximum_hash)ierr=MPI_ERR_OTHER
    probe=''
#else
    probe='';ierr=1
#endif
  end subroutine validate_path
  subroutine sync_io(io_status,comm,ierr)
    integer,intent(inout)::io_status;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
#else
    ierr=1
#endif
  end subroutine sync_io
  subroutine atomic_rename(old_path,new_path,status)
    character(*),intent(in)::old_path,new_path;integer,intent(out)::status
    character(c_char),allocatable::old_c(:),new_c(:)
    integer::i
    allocate(old_c(len_trim(old_path)+1),new_c(len_trim(new_path)+1))
    do i=1,len_trim(old_path);old_c(i)=old_path(i:i);enddo;old_c(size(old_c))=c_null_char
    do i=1,len_trim(new_path);new_c(i)=new_path(i:i);enddo;new_c(size(new_c))=c_null_char
    status=int(c_rename(old_c,new_c))
  end subroutine atomic_rename
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix

  pure integer(int64) function mix_hash(seed,value)
    integer(int64),intent(in)::seed,value
    mix_hash=ieor(ishftc(seed,9),value)
  end function mix_hash

  function int64_string(value) result(text)
    integer(int64),intent(in)::value;character(32)::text
    write(text,'(i0)')value
  end function int64_string
end module rt_dg_hybrid_checkpoint
