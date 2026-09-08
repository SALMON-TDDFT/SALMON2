#include "config.h"
#ifdef DG_GENERALIZED_PIPELINE_CONTRACT_SYNTAX
module dg_hybrid_projected_fragment_pipeline_red_contract
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  type,public::s_dg_hybrid_projection_factorization_receipt
    logical::valid=.false.
    integer::basis_generation=0,metric_rank=0,metric_factorization_count=0,projected_tile_count=0
    integer(int64)::wannier_fingerprint=0_int64,metric_fingerprint=0_int64
  end type s_dg_hybrid_projection_factorization_receipt
  interface
    subroutine build_dg_hybrid_projected_fragment_basis(comm,global_point_count,fragment_count,&
        fragment_id,core_ids,weights,core_wannier,core_coordinates,core_windows,buffer_ids,&
        buffer_wannier,buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,tile_width,&
        tolerance,wannier_fingerprint,basis,workspace_peak_bytes,fingerprint,ok,message,&
        basis_generation,projection_receipt)
      import::int64,real64,s_dg_hybrid_basis_catalog,s_dg_hybrid_fragment_basis,&
        s_dg_hybrid_projection_factorization_receipt
      integer,intent(in)::comm,global_point_count,fragment_count,fragment_id,tile_width,wannier_owner(:)
      integer(int64),intent(in)::core_ids(:),buffer_ids(:),wannier_fingerprint
      real(real64),intent(in)::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),&
        buffer_windows(:,:),g_vectors(:,:),tolerance
      complex(real64),intent(in)::core_wannier(:,:),buffer_wannier(:,:)
      type(s_dg_hybrid_basis_catalog),intent(in)::catalog
      type(s_dg_hybrid_fragment_basis),intent(out)::basis
      integer(int64),intent(out)::workspace_peak_bytes,fingerprint
      logical,intent(out)::ok
      character(*),intent(out)::message
      integer,intent(in),optional::basis_generation
      type(s_dg_hybrid_projection_factorization_receipt),intent(out),optional::projection_receipt
    end subroutine build_dg_hybrid_projected_fragment_basis
  end interface
end module dg_hybrid_projected_fragment_pipeline_red_contract
#endif

program test_dg_hybrid_projected_fragment_pipeline_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
#ifdef DG_GENERALIZED_PIPELINE_CONTRACT_SYNTAX
  use dg_hybrid_projected_fragment_pipeline_red_contract,only:build_dg_hybrid_projected_fragment_basis,&
    s_dg_hybrid_projection_factorization_receipt
#else
  use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_fragment_basis,&
    s_dg_hybrid_projection_factorization_receipt,build_dg_hybrid_projected_local_fragment_basis
#endif
  implicit none
  real(real64),parameter::global_point_weights(4)=[0.5_real64,1.25_real64,2.0_real64,0.75_real64]
  integer::comm,rank,nproc,ierr,fragment_id,nlocal,i,j,mismatched_generation
  integer(int64),allocatable::core_ids(:),buffer_ids(:)
  real(real64),allocatable::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),buffer_windows(:,:)
  real(real64),allocatable::g_vectors(:,:)
  complex(real64),allocatable::core_wannier(:,:),buffer_wannier(:,:),local_full(:,:),global_full(:,:)
  integer::wannier_owner(2)
  integer,allocatable::support_slots(:)
  type(s_dg_hybrid_basis_catalog)::catalog
  type(s_dg_hybrid_fragment_basis)::basis
  type(s_dg_hybrid_projection_factorization_receipt)::projection_receipt
  integer(int64)::workspace,fingerprint
  logical::ok,values_ok,owner_contract_ok,basis_metadata_ok,generation_failure_ok,&
    late_tile_failure_ok,global_metadata_ok,global_generation_failure_ok,global_late_tile_failure_ok
  complex(real64)::cross_overlap
  character(256)::message,generation_message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  fragment_id=0;if(rank<2)fragment_id=rank+1
  nlocal=0;if(rank<2)nlocal=2
  allocate(core_ids(nlocal),weights(nlocal),core_coordinates(3,nlocal),core_windows(2,nlocal),core_wannier(2,nlocal))
  if(rank==0)core_ids=[1_int64,2_int64]
  if(rank==1)core_ids=[3_int64,4_int64]
  core_coordinates=0d0;core_windows=0d0;core_wannier=(0d0,0d0)
  do i=1,nlocal
    weights(i)=global_point_weights(int(core_ids(i)))
    core_coordinates(1,i)=real(core_ids(i)-1_int64,real64)
    if(core_ids(i)<=2)core_windows(1,i)=1d0
    if(core_ids(i)>=3)core_windows(2,i)=1d0
    if(core_ids(i)==1)then
      core_wannier(1,i)=1d0/sqrt(global_point_weights(1))
      core_wannier(2,i)=1d0/sqrt(2d0*global_point_weights(1))
    endif
    if(core_ids(i)==3)core_wannier(2,i)=cmplx(0d0,1d0,real64)/sqrt(2d0*global_point_weights(3))
  enddo
  if(fragment_id>0)then
    allocate(buffer_ids(4),source=[1_int64,2_int64,3_int64,4_int64])
    allocate(buffer_coordinates(3,4),buffer_windows(2,4),buffer_wannier(2,4))
    buffer_coordinates=0d0;buffer_windows=0d0;buffer_wannier=(0d0,0d0)
    do i=1,4;buffer_coordinates(1,i)=real(i-1,real64);enddo
    buffer_windows(1,1:2)=1d0;buffer_windows(2,3:4)=1d0
    buffer_wannier(1,1)=1d0/sqrt(global_point_weights(1))
    buffer_wannier(2,1)=1d0/sqrt(2d0*global_point_weights(1))
    buffer_wannier(2,3)=cmplx(0d0,1d0,real64)/sqrt(2d0*global_point_weights(3))
  else
    allocate(buffer_ids(0),buffer_coordinates(3,0),buffer_windows(2,0),buffer_wannier(2,0))
  endif
  allocate(g_vectors(3,1));g_vectors=0d0;allocate(catalog%packets(2));catalog%valid=.true.
  do i=1,2
    catalog%packets(i)%fragment_id=i;catalog%packets(i)%star_id=1;catalog%packets(i)%owner_rank=i-1
    allocate(catalog%packets(i)%g_indices(1),source=[1])
  enddo
  catalog%packet_fingerprint=31_int64;catalog%catalog_fingerprint=37_int64;wannier_owner=[1,2]
  call build_dg_hybrid_projected_fragment_basis(comm,4,2,fragment_id,core_ids,weights,core_wannier,&
    core_coordinates,core_windows,buffer_ids,buffer_wannier,buffer_coordinates,buffer_windows,catalog,&
    g_vectors,wannier_owner,1,1d-12,41_int64,basis,workspace,fingerprint,ok,message,&
    basis_generation=17,projection_receipt=projection_receipt)
  call require(ok,'projected fragment pipeline failed: '//trim(message))
  owner_contract_ok=allocated(basis%global_ids).and.allocated(basis%sector).and.&
    allocated(basis%buffer_point_ids).and.allocated(basis%buffer_values)
  if(owner_contract_ok)then
    if(fragment_id==1)then
      owner_contract_ok=size(basis%global_ids)==2.and.size(basis%sector)==2.and.&
        size(basis%buffer_point_ids)==size(buffer_ids).and.&
        size(basis%buffer_values,1)==size(buffer_ids).and.size(basis%buffer_values,2)==2
      if(owner_contract_ok)owner_contract_ok=all(basis%global_ids==[1_int64,3_int64]).and.&
        all(basis%sector==[1,2]).and.all(basis%buffer_point_ids==buffer_ids).and.&
        bitwise_complex_vector_equal(basis%buffer_values(:,1),buffer_wannier(1,:))
    elseif(fragment_id==2)then
      owner_contract_ok=size(basis%global_ids)==2.and.size(basis%sector)==2.and.&
        size(basis%buffer_point_ids)==size(buffer_ids).and.&
        size(basis%buffer_values,1)==size(buffer_ids).and.size(basis%buffer_values,2)==2
      if(owner_contract_ok)owner_contract_ok=all(basis%global_ids==[2_int64,4_int64]).and.&
        all(basis%sector==[1,2]).and.all(basis%buffer_point_ids==buffer_ids).and.&
        bitwise_complex_vector_equal(basis%buffer_values(:,1),buffer_wannier(2,:))
    else
      owner_contract_ok=size(basis%global_ids)==0.and.size(basis%sector)==0.and.&
        size(basis%buffer_point_ids)==0.and.size(basis%buffer_values)==0
    endif
  endif
  call require(owner_contract_ok,&
    'pipeline changed owner IDs/sectors/buffer order/WF values or published data on an idle rank')
  if(fragment_id>0)then
    basis_metadata_ok=basis%fragment_id==fragment_id.and.basis%generation==17.and.&
      basis%provenance_fingerprint/=0_int64
  else
    basis_metadata_ok=basis%fragment_id==0.and.basis%generation==0.and.&
      basis%provenance_fingerprint==0_int64.and.allocated(basis%global_ids).and.&
      allocated(basis%sector).and.allocated(basis%buffer_point_ids).and.allocated(basis%buffer_values)
    if(basis_metadata_ok)basis_metadata_ok=size(basis%global_ids)==0.and.size(basis%sector)==0.and.&
      size(basis%buffer_point_ids)==0.and.size(basis%buffer_values)==0
  endif
  allocate(local_full(4,4),global_full(4,4));local_full=(0d0,0d0)
  do j=1,size(basis%global_ids);do i=1,size(buffer_ids)
    local_full(int(buffer_ids(i)),int(basis%global_ids(j)))=basis%buffer_values(i,j)
  enddo;enddo
  call MPI_Allreduce(local_full,global_full,16,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  cross_overlap=sum(global_point_weights*conjg(global_full(:,1))*global_full(:,2))
  values_ok=abs(cross_overlap)>0.1d0.and.&
    abs(sum(global_point_weights*conjg(global_full(:,1))*global_full(:,3)))<1d-12.and.&
    abs(sum(global_point_weights*conjg(global_full(:,2))*global_full(:,3)))<1d-12.and.&
    abs(sum(global_point_weights*conjg(global_full(:,1))*global_full(:,4)))<1d-12.and.&
    abs(sum(global_point_weights*conjg(global_full(:,2))*global_full(:,4)))<1d-12
  call require(values_ok,'pipeline did not use generalized projection for a nonidentity fragment-union Gram')
  call require(workspace<=512_int64,'projected fragment pipeline workspace is not tile bounded')
  call require(fingerprint/=0_int64,'projected fragment pipeline fingerprint is empty')
  call require((fragment_id==0.or.basis%generation==17).and.projection_receipt%valid.and.&
    projection_receipt%basis_generation==17,&
    'generalized projection factorization receipt has wrong basis generation')
  call require(projection_receipt%metric_rank==2.and.projection_receipt%metric_factorization_count==1.and.&
    projection_receipt%projected_tile_count==2,&
    'generalized WF Gram was not factorized exactly once and reused across both PW tiles')
  call require(projection_receipt%wannier_fingerprint==41_int64.and.&
    projection_receipt%metric_fingerprint/=0_int64,&
    'generalized metric factorization receipt lacks observable provenance')
  if(rank==0)write(*,'(a,i0,a,i0)')'PROJECTED_FRAGMENT ranks=',nproc,' fingerprint=',fingerprint

#ifndef DG_GENERALIZED_PIPELINE_CONTRACT_SYNTAX
  call test_variable_fragment_columns
  if(nproc==2)then
    call build_dg_hybrid_projected_local_fragment_basis(comm,4,2,fragment_id,core_ids,weights,&
      core_coordinates,core_windows,buffer_ids,buffer_wannier(fragment_id:fragment_id,:),&
      buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,1,1d-12,41_int64,&
      basis,workspace,fingerprint,ok,message,basis_generation=17,projection_receipt=projection_receipt)
    call require(ok,'local fragment input pipeline failed: '//trim(message))
    values_ok=basis%generation==17.and.projection_receipt%valid
    do j=1,size(basis%global_ids)
      values_ok=values_ok.and.maxval(abs(basis%buffer_values(:,j)-global_full(:,int(basis%global_ids(j)))))<1d-12
    enddo
    call require(values_ok,'local fragment input differs from full-union generalized projection')
    if(rank==0)then
      support_slots=[2,1]
    else
      support_slots=[3,1,4]
    endif
    call build_dg_hybrid_projected_local_fragment_basis(comm,4,2,fragment_id,core_ids,weights,&
      core_coordinates,core_windows,buffer_ids(support_slots),&
      buffer_wannier(fragment_id:fragment_id,support_slots),buffer_coordinates(:,support_slots),&
      buffer_windows(:,support_slots),catalog,g_vectors,wannier_owner,2,1d-12,41_int64,&
      basis,workspace,fingerprint,ok,message,basis_generation=17)
    call require(ok,'unequal reordered local WF support failed: '//trim(message))
    values_ok=all(basis%buffer_point_ids==buffer_ids(support_slots))
    do j=1,size(basis%global_ids)
      values_ok=values_ok.and.maxval(abs(basis%buffer_values(:,j)-&
        global_full(support_slots,int(basis%global_ids(j)))))<1d-12
    enddo
    call require(values_ok,'support exchange changed physical point ordering or zero extension')
    if(rank==0)buffer_ids(2)=buffer_ids(1)
    call build_dg_hybrid_projected_local_fragment_basis(comm,4,2,fragment_id,core_ids,weights,&
      core_coordinates,core_windows,buffer_ids,buffer_wannier(fragment_id:fragment_id,:),&
      buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,1,1d-12,41_int64,&
      basis,workspace,fingerprint,ok,message)
    call require(.not.ok.and..not.allocated(basis%global_ids),'duplicate source physical IDs accepted')
    buffer_ids(2)=2_int64
  endif
#endif

  mismatched_generation=17;if(rank==0)mismatched_generation=18
  call build_dg_hybrid_projected_fragment_basis(comm,4,2,fragment_id,core_ids,weights,core_wannier,&
    core_coordinates,core_windows,buffer_ids,buffer_wannier,buffer_coordinates,buffer_windows,catalog,&
    g_vectors,wannier_owner,1,1d-12,41_int64,basis,workspace,fingerprint,ok,message,&
    basis_generation=mismatched_generation,projection_receipt=projection_receipt)
  generation_failure_ok=.not.ok.and.index(lowercase(message),'generation')>0.and.&
    fragment_basis_unpublished(basis).and.projection_receipt_unpublished(projection_receipt)
  generation_message=message

  ! Packet one remains valid and is materialized/projected first.  Only the
  ! reciprocal index consumed by tile two is invalid, exercising rollback
  ! after a successful tile rather than initial-input rejection.
  catalog%packets(2)%g_indices(1)=size(g_vectors,2)+1
  call build_dg_hybrid_projected_fragment_basis(comm,4,2,fragment_id,core_ids,weights,core_wannier,&
    core_coordinates,core_windows,buffer_ids,buffer_wannier,buffer_coordinates,buffer_windows,catalog,&
    g_vectors,wannier_owner,1,1d-12,41_int64,basis,workspace,fingerprint,ok,message,&
    basis_generation=17,projection_receipt=projection_receipt)
  late_tile_failure_ok=.not.ok.and.index(lowercase(message),'reciprocal')>0.and.&
    fragment_basis_unpublished(basis).and.projection_receipt_unpublished(projection_receipt).and.&
    workspace==0_int64.and.fingerprint==0_int64
  call MPI_Allreduce(basis_metadata_ok,global_metadata_ok,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
  call MPI_Allreduce(generation_failure_ok,global_generation_failure_ok,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
  call MPI_Allreduce(late_tile_failure_ok,global_late_tile_failure_ok,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
  if(rank==0.and..not.global_metadata_ok)write(0,'(a)')&
    'EXPECTED RED: active/idle fragment basis metadata publication contract is not enforced'
  if(rank==0.and..not.global_generation_failure_ok)write(0,'(a,a)')&
    'EXPECTED RED: rank-dependent basis generation was not rejected without publication; diagnostic=',&
    trim(generation_message)
  if(rank==0.and..not.global_late_tile_failure_ok)write(0,'(a,a,a,i0,a,i0)')&
    'EXPECTED RED: tile-two failure leaked scalar outputs; diagnostic=',trim(message),&
    ' workspace=',workspace,' fingerprint=',fingerprint
  call require(global_metadata_ok.and.global_generation_failure_ok.and.global_late_tile_failure_ok,&
    'projected fragment metadata/generation/two-phase quality contracts are not implemented')
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid projected fragment pipeline on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
#ifndef DG_GENERALIZED_PIPELINE_CONTRACT_SYNTAX
  subroutine test_variable_fragment_columns
    type(s_dg_hybrid_basis_catalog)::pw_catalog
    type(s_dg_hybrid_fragment_basis)::reference,actual
    integer::nf,ng,nw,f,a,b,k,first,last,neighbor
    integer,allocatable::owners(:)
    integer(int64)::cids(2),bids(3)
    complex(real64),allocatable::wf(:,:)
    real(real64),allocatable::windows(:,:)
    real(real64)::cc(3,2),bc(3,3),gv(3,1)
    ! Every MPI rank owns exactly one fragment, in reversed rank order.
    ! Alternating one/two WF columns force multiple width-one tiles.
    nf=nproc;ng=2*nf;nw=sum([(1+modulo(a,2),a=1,nf)])
    allocate(owners(nw),wf(nw,ng),windows(nf,ng),pw_catalog%packets(nf))
    wf=(0d0,0d0);windows=0d0;k=0;f=nf-rank;first=0;last=0
    do a=1,nf
      neighbor=2*modulo(a,nf)+1
      if(a==f)first=k+1
      do b=1,1+modulo(a,2)
        k=k+1;owners(k)=a
        wf(k,2*a-2+b)=cmplx(1d0,0.1d0*b,real64)
        wf(k,neighbor)=cmplx(0.1d0*b,-0.05d0*a,real64)
      enddo
      if(a==f)last=k
      windows(a,2*a-1:2*a)=1d0
      pw_catalog%packets(a)%fragment_id=a
      pw_catalog%packets(a)%star_id=1;pw_catalog%packets(a)%owner_rank=nf-a
      allocate(pw_catalog%packets(a)%g_indices(1),source=[1])
    enddo
    pw_catalog%valid=.true.;pw_catalog%packet_fingerprint=31_int64;pw_catalog%catalog_fingerprint=37_int64
    cids=int([2*f,2*f-1],int64);bids=int([2*f,2*modulo(f,nf)+1,2*f-1],int64)
    cc=0d0;bc=0d0;gv=0d0
    cc(1,:)=real(cids-1_int64,real64);bc(1,:)=real(bids-1_int64,real64)
    call build_dg_hybrid_projected_fragment_basis(comm,ng,nf,f,cids,[1d0,1d0],wf(:,cids),&
      cc,windows(:,cids),bids,wf(:,bids),bc,windows(:,bids),pw_catalog,gv,owners,1,1d-12,41_int64,&
      reference,workspace,fingerprint,ok,message,basis_generation=5)
    call require(ok,'variable-column reference failed: '//trim(message))
    call build_dg_hybrid_projected_local_fragment_basis(comm,ng,nf,f,cids,[1d0,1d0],cc,windows(:,cids),&
      bids,wf(first:last,bids),bc,windows(:,bids),pw_catalog,gv,owners,1,1d-12,41_int64,&
      actual,workspace,fingerprint,ok,message,basis_generation=5)
    call require(ok,'variable-column local projection failed: '//trim(message))
    call require(all(actual%global_ids==reference%global_ids).and.&
      all(actual%buffer_point_ids==reference%buffer_point_ids).and.actual%generation==5.and.&
      maxval(abs(actual%buffer_values-reference%buffer_values))<1d-12,&
      'multiple local WF tiles changed column IDs, physical ordering or projected values')
  end subroutine test_variable_fragment_columns
#endif

  logical function bitwise_complex_vector_equal(left,right)result(equal)
    complex(real64),intent(in)::left(:),right(:)
    integer::k
    integer(int64)::left_bits,right_bits
    equal=size(left)==size(right);if(.not.equal)return
    do k=1,size(left)
      left_bits=transfer(real(left(k),real64),left_bits);right_bits=transfer(real(right(k),real64),right_bits)
      if(left_bits/=right_bits)then;equal=.false.;return;endif
      left_bits=transfer(aimag(left(k)),left_bits);right_bits=transfer(aimag(right(k)),right_bits)
      if(left_bits/=right_bits)then;equal=.false.;return;endif
    enddo
  end function bitwise_complex_vector_equal

  logical function fragment_basis_unpublished(candidate)result(unpublished)
    type(s_dg_hybrid_fragment_basis),intent(in)::candidate
    unpublished=candidate%fragment_id==0.and.candidate%generation==0.and.&
      candidate%provenance_fingerprint==0_int64
    if(allocated(candidate%global_ids))unpublished=unpublished.and.size(candidate%global_ids)==0
    if(allocated(candidate%sector))unpublished=unpublished.and.size(candidate%sector)==0
    if(allocated(candidate%buffer_point_ids))unpublished=unpublished.and.size(candidate%buffer_point_ids)==0
    if(allocated(candidate%buffer_values))unpublished=unpublished.and.size(candidate%buffer_values)==0
  end function fragment_basis_unpublished

  logical function projection_receipt_unpublished(receipt)result(unpublished)
    type(s_dg_hybrid_projection_factorization_receipt),intent(in)::receipt
    unpublished=.not.receipt%valid.and.receipt%basis_generation==0.and.receipt%metric_rank==0.and.&
      receipt%metric_factorization_count==0.and.receipt%projected_tile_count==0.and.&
      receipt%wannier_fingerprint==0_int64.and.receipt%metric_fingerprint==0_int64
  end function projection_receipt_unpublished

  pure function lowercase(text)result(lowered)
    character(*),intent(in)::text
    character(len(text))::lowered
    integer::k,code
    lowered=text
    do k=1,len(text)
      code=iachar(text(k:k));if(code>=iachar('A').and.code<=iachar('Z'))lowered(k:k)=achar(code+32)
    enddo
  end function lowercase

  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_projected_fragment_pipeline_mpi
