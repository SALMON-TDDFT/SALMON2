#include "config.h"
program test_rt_dg_hybrid_checkpoint_v5_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_checkpoint_v5,only:s_rt_dg_hybrid_v5_shard,&
    read_rt_dg_hybrid_checkpoint_v5,write_rt_dg_hybrid_checkpoint_v5,&
    checked_rt_dg_hybrid_extent_product
  use rt_dg_hybrid_checkpoint_v5,only:publish_rt_dg_hybrid_checkpoint_v5,&
    s_rt_dg_hybrid_v5_publication_authorization
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,initialize_rt_dg_hybrid_from_checkpoint,&
    fingerprint_rt_dg_hybrid_scope
  use dg_hybrid_continuation_controller,only:s_dg_hybrid_candidate_acceptance,&
    initialize_dg_hybrid_candidate_acceptance,record_dg_hybrid_complete_lcfo_solve,&
    record_dg_hybrid_occupation_policy,record_dg_hybrid_unconditional_gates,&
    record_dg_hybrid_spectral_certification,record_dg_hybrid_certified_rt_basis,&
    authorize_dg_hybrid_v5_publication
  implicit none
  type(s_rt_dg_hybrid_v5_shard)::written,loaded
  type(s_rt_dg_hybrid_state)::state
  type(s_rt_dg_hybrid_v5_publication_authorization)::authorization
  type(s_dg_hybrid_candidate_acceptance)::candidate
  integer::comm,rank,nproc,ierr,i,j,environment_status
  integer(int64)::extent_product
  integer,allocatable::row_owner(:)
  logical::ok
  character(256)::message
  character(512)::prefix,failure_prefix
  character(32)::test_mode
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call checked_rt_dg_hybrid_extent_product(7_int64,9_int64,extent_product,ok)
  call require(ok.and.extent_product==63_int64,'checked extent product changed a finite product')
  call checked_rt_dg_hybrid_extent_product(huge(0_int64),2_int64,extent_product,ok)
  call require(.not.ok.and.extent_product==0_int64,'checked extent product failed to reject int64 overflow')
  write(prefix,'(a,i0)')'/tmp/salmon-hybrid-v5-checkpoint-',nproc
  written%global_count=50*nproc;written%global_grid_count=100*nproc;written%nocc=7
  written%certified_rank=written%global_count;written%fragment_id=rank+1
  written%basis_fingerprint=7717_int64;written%operator_fingerprint=9919_int64
  written%operator_structure_fingerprint=1217_int64;written%scope_fingerprint=1811_int64
  written%payload_fingerprint=2027_int64
  written%system_fingerprint=[3037_int64,3038_int64,3039_int64,3040_int64]
  written%pseudopotential_fingerprint=4049_int64
  written%pseudopotential_digest=[4049_int64,4050_int64,4051_int64,4052_int64]
  allocate(row_owner(written%global_count))
  do i=1,written%global_count;row_owner(i)=(i-1)/50;enddo
  allocate(written%row_ids(50),written%metric_offsets(51),written%metric_columns(50),&
    written%metric_values(50),written%operator_offsets(51),written%operator_columns(50),&
    written%operator_values(50),written%kinetic_values(50),written%nonlocal_values(50),&
    written%local_values(50),written%sipg_values(50),written%position_values(3,50),&
    written%grid_ids(100),written%basis_point_offsets(101),written%basis_support_ids(100),&
    written%basis_support_values(100),written%grid_weights(100),written%density(100),&
    written%initial_occupied_amplitudes(50,7),written%occupations(7),written%eigenvalues(7),&
    written%scope_selectors(8),written%xc_types(1),written%acceptance_receipts(8),&
    written%pseudopotential_receipt(6),written%energy_receipt(7))
  do i=1,50
    written%row_ids(i)=int(rank*50+i,int64)
    written%metric_offsets(i)=i;written%operator_offsets(i)=i
    written%metric_columns(i)=rank*50+i;written%operator_columns(i)=rank*50+i
    written%metric_values(i)=(1d0,0d0)
    written%operator_values(i)=cmplx(merge(-0.5d0+0.1d0*(rank*50+i),&
      0.25d0+0.01d0*(rank*50+i),rank*50+i<=7),0d0,real64)
    written%kinetic_values(i)=written%operator_values(i)
    written%nonlocal_values(i)=(0d0,0d0);written%local_values(i)=(0d0,0d0)
    written%sipg_values(i)=(0d0,0d0);written%position_values(:,i)=(0d0,0d0)
    do j=1,7
      written%initial_occupied_amplitudes(i,j)=merge((1d0,0d0),(0d0,0d0),rank*50+i==j)
    enddo
  enddo
  written%metric_offsets(51)=51;written%operator_offsets(51)=51
  do i=1,100
    written%grid_ids(i)=int(rank*100+i,int64);written%grid_weights(i)=0.5d0
    written%basis_point_offsets(i)=i
    written%basis_support_ids(i)=rank*50+mod(i-1,50)+1
    written%basis_support_values(i)=cmplx(0.1d0+0.001d0*i,-0.02d0,real64)
    if(written%basis_support_ids(i)<=7)then
      written%density(i)=abs(written%basis_support_values(i))**2
    else
      written%density(i)=0d0
    endif
  enddo
  written%basis_point_offsets(101)=101
  do j=1,7;written%occupations(j)=1d0;written%eigenvalues(j)=-0.5d0+0.1d0*j;enddo
  written%scope_selectors=[1,1,1,0,0,0,0,0];written%xc_types=[1]
  written%scope_fingerprint=fingerprint_rt_dg_hybrid_scope(written%scope_selectors,written%xc_types)
  written%acceptance_receipts=[(0.01d0*i,i=1,8)]
  written%pseudopotential_receipt=[(1d0*i,i=1,6)]
  written%energy_receipt=[(2d0*i,i=1,7)]
  authorization%checkpoint_version=5;authorization%published_rank=written%global_count
  authorization%basis_fingerprint=written%basis_fingerprint
  authorization%operator_fingerprint=written%operator_fingerprint
  call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
  call require(.not.ok.and.index(message,'authorization')>0,&
    'common v5 endpoint accepted a payload without formal publication authorization')
  call initialize_dg_hybrid_candidate_acceptance(comm,written%global_count,-1d0,candidate,ok,message)
  call require(ok,'full-rank controller initialization failed: '//trim(message))
  call record_dg_hybrid_complete_lcfo_solve(comm,candidate,written%global_count,1,3001_int64,ok,message)
  call require(ok,'full-rank complete LCFO receipt failed: '//trim(message))
  call record_dg_hybrid_occupation_policy(comm,candidate,written%nocc,.true.,3002_int64,ok,message)
  call require(ok,'full-rank occupation receipt failed: '//trim(message))
  call record_dg_hybrid_unconditional_gates(comm,candidate,.true.,.true.,1d-13,2d-13,1d-10,ok,message)
  call require(ok,'full-rank physical gates failed: '//trim(message))
  call record_dg_hybrid_spectral_certification(comm,candidate,written%nocc,written%global_count,&
    written%global_count,.false.,.true.,.true.,3003_int64,ok,message)
  call require(ok.and.candidate%certified_rank==written%global_count,&
    'energy_window=-1 full-rank certification was rejected: '//trim(message))
  call record_dg_hybrid_certified_rt_basis(comm,candidate,written%global_count,written%basis_fingerprint,&
    written%operator_fingerprint,ok,message)
  call require(ok,'full-rank certified basis receipt failed: '//trim(message))
  call authorize_dg_hybrid_v5_publication(comm,candidate,5,written%global_count,.true.,ok,message)
  call require(ok.and.candidate%publication_authorized,'full-rank v5 publication authorization failed: '//trim(message))
  authorization%valid=candidate%publication_authorized
  authorization%checkpoint_version=candidate%checkpoint_version
  authorization%published_rank=candidate%published_rt_rank
  authorization%basis_fingerprint=candidate%basis_fingerprint
  authorization%operator_fingerprint=candidate%operator_fingerprint
  write(failure_prefix,'(a,i0,a)')'/tmp/salmon-hybrid-v5-open-failure-',nproc,'/missing/checkpoint'
  call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(failure_prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
  call require(.not.ok.and.index(message,'cannot atomically publish distributed-v5 rank shard')>0,&
    'v5 writer OPEN failure was not collectively rejected')
  call get_environment_variable('SALMON_TEST_V5_MANIFEST_OPEN_FAILURE',test_mode,status=environment_status)
  if(environment_status==0.and.trim(test_mode)=='1')then
    write(failure_prefix,'(a,i0)')'/tmp/salmon-hybrid-v5-manifest-open-failure-',nproc
    call write_rt_dg_hybrid_checkpoint_v5(comm,trim(failure_prefix),written,ok,message)
    call require(.not.ok.and.index(message,'cannot atomically publish distributed-v5 manifest')>0,&
      'v5 manifest OPEN failure was not collectively rejected')
    if(rank==0)write(*,'(a,i0)')'PASS v5 collective manifest OPEN failure ranks=',nproc
    call MPI_Finalize(ierr);stop
  endif
  call get_environment_variable('SALMON_TEST_V5_MANIFEST_FAILURE',test_mode,status=environment_status)
  if(environment_status==0.and.trim(test_mode)=='1')then
    write(failure_prefix,'(a,i0)')'/tmp/salmon-hybrid-v5-manifest-failure-',nproc
    call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(failure_prefix),written%global_count,written%nocc,&
      written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
    call require(.not.ok.and.index(message,'cannot atomically publish distributed-v5 manifest')>0,&
      'v5 manifest rename failure was not collectively rejected')
    if(rank==0)write(*,'(a,i0)')'PASS v5 collective manifest failure ranks=',nproc
    call MPI_Finalize(ierr);stop
  endif
  written%certified_rank=written%global_count+1
  call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
  call require(.not.ok.and.index(message,'invalid distributed-v5 rank shard payload')>0,&
    'certified rank greater than construction rank was accepted')
  written%certified_rank=written%nocc-1
  call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
  call require(.not.ok.and.index(message,'invalid distributed-v5 rank shard payload')>0,&
    'certified rank below occupied rank was accepted')
  written%certified_rank=written%global_count
  if(nproc>1)then
    written%scope_fingerprint=written%scope_fingerprint+rank
    call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
      written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
    call require(.not.ok.and.index(message,'rank-inconsistent')>0,&
      'rank-inconsistent common metadata reached shard publication')
    written%scope_fingerprint=fingerprint_rt_dg_hybrid_scope(written%scope_selectors,written%xc_types)
    written%system_fingerprint(1)=3037_int64+rank
    call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
      written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
    call require(.not.ok.and.index(message,'rank-inconsistent')>0,&
      'rank-inconsistent system identity reached shard publication')
    written%system_fingerprint=[3037_int64,3038_int64,3039_int64,3040_int64]
    written%pseudopotential_fingerprint=4049_int64+rank
    call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
      written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
    call require(.not.ok.and.index(message,'rank-inconsistent')>0,&
      'rank-inconsistent canonical PP identity reached shard publication')
    written%pseudopotential_fingerprint=4049_int64
    written%pseudopotential_digest(1)=4049_int64+rank
    call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
      written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
    call require(.not.ok.and.index(message,'rank-inconsistent')>0,&
      'rank-inconsistent authoritative PP digest reached shard publication')
    written%pseudopotential_digest=[4049_int64,4050_int64,4051_int64,4052_int64]
  endif
  call publish_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,authorization,.true.,ok,message)
  call require(ok,'v5 shard publication failed: '//trim(message))
  call read_rt_dg_hybrid_checkpoint_v5(comm,trim(prefix),loaded,ok,message)
  call require(ok,'v5 shard reload failed: '//trim(message))
  call require(loaded%global_count==written%global_count.and.loaded%nocc==written%nocc,&
    'v5 manifest dimensions changed')
  call require(loaded%fragment_id==written%fragment_id,'v5 rank-fragment mapping changed')
  call require(all(loaded%system_fingerprint==written%system_fingerprint).and.&
    all(loaded%pseudopotential_digest==written%pseudopotential_digest).and.&
    loaded%pseudopotential_fingerprint==written%pseudopotential_fingerprint,&
    'v5 physical-system identity changed')
  call require(all(loaded%row_ids==written%row_ids),'v5 owned rows changed')
  call require(all(loaded%metric_offsets==written%metric_offsets).and.&
    all(loaded%metric_columns==written%metric_columns).and.&
    maxval(abs(loaded%metric_values-written%metric_values))==0d0,'v5 metric CSR changed')
  call require(all(loaded%operator_offsets==written%operator_offsets).and.&
    all(loaded%operator_columns==written%operator_columns).and.&
    maxval(abs(loaded%operator_values-written%operator_values))==0d0,'v5 operator CSR changed')
  call require(maxval(abs(loaded%position_values-written%position_values))==0d0.and.&
    maxval(abs(loaded%local_values-written%local_values))==0d0,'v5 component CSR changed')
  call require(all(loaded%basis_point_offsets==written%basis_point_offsets).and.&
    all(loaded%basis_support_ids==written%basis_support_ids).and.&
    maxval(abs(loaded%basis_support_values-written%basis_support_values))==0d0,'v5 point CSR changed')
  call require(maxval(abs(loaded%initial_occupied_amplitudes-written%initial_occupied_amplitudes))==0d0,&
    'v5 distributed occupied coefficients changed')
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(prefix),'tddft_response',.true.,1,&
    .false.,.false.,.false.,.false.,.false.,[1],written%system_fingerprint,&
    written%pseudopotential_fingerprint,&
    written%pseudopotential_digest,&
    [1d-10,1d-10,1d-10,1d-10],state,ok,message)
  call require(ok.and.state%valid.and.state%initial_invariants_valid,&
    'common v5 endpoint did not initialize distributed Hybrid RT: '//trim(message))
  call require(state%startup_metric_defect<=1d-12.and.state%startup_orbital_residual<=1d-12,&
    'common v5 endpoint changed metric orthonormality or stationarity')
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(prefix),'tddft_response',.true.,1,&
    .false.,.false.,.false.,.false.,.false.,[1],written%system_fingerprint+[1_int64,0_int64,0_int64,0_int64],&
    written%pseudopotential_fingerprint,&
    written%pseudopotential_digest,&
    [1d-10,1d-10,1d-10,1d-10],state,ok,message)
  call require(.not.ok.and.index(message,'system identity mismatch')>0,&
    'common v5 endpoint accepted a different current physical system')
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(prefix),'tddft_response',.true.,1,&
    .false.,.false.,.false.,.false.,.false.,[1],written%system_fingerprint,&
    written%pseudopotential_fingerprint+1_int64,written%pseudopotential_digest,&
    [1d-10,1d-10,1d-10,1d-10],state,ok,message)
  call require(.not.ok.and.index(message,'system identity mismatch')>0,&
    'common v5 endpoint accepted a different current canonical pseudopotential')
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(prefix),'tddft_response',.true.,1,&
    .false.,.false.,.false.,.false.,.false.,[1],written%system_fingerprint,&
    written%pseudopotential_fingerprint,&
    written%pseudopotential_digest+[1_int64,0_int64,0_int64,0_int64],&
    [1d-10,1d-10,1d-10,1d-10],state,ok,message)
  call require(.not.ok.and.index(message,'system identity mismatch')>0,&
    'common endpoint accepted a different authoritative PP digest')
  if(rank==0)write(*,'(a,i0)')'PASS distributed-v5 shard manifest ranks=',nproc
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    integer::bad,global_bad,status
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.global_bad/=0)then
      if(.not.condition)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,status)
    endif
  end subroutine require
end program test_rt_dg_hybrid_checkpoint_v5_mpi
