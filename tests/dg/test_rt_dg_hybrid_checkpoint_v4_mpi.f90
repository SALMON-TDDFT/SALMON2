#include "config.h"
program test_rt_dg_hybrid_checkpoint_v4_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_checkpoint_v4,only:s_rt_dg_hybrid_v4_shard,&
    read_rt_dg_hybrid_checkpoint_v4
  use rt_dg_hybrid_checkpoint,only:publish_rt_dg_hybrid_checkpoint_v4
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,initialize_rt_dg_hybrid_from_checkpoint,&
    fingerprint_rt_dg_hybrid_scope
  implicit none
  type(s_rt_dg_hybrid_v4_shard)::written,loaded
  type(s_rt_dg_hybrid_state)::state
  integer::comm,rank,nproc,ierr,i,j
  integer,allocatable::row_owner(:)
  logical::ok
  character(256)::message
  character(512)::prefix
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  write(prefix,'(a,i0)')'/tmp/salmon-hybrid-v4-checkpoint-',nproc
  written%global_count=50*nproc;written%global_grid_count=100*nproc;written%nocc=7
  written%certified_rank=written%global_count-1;written%fragment_id=rank+1
  written%basis_fingerprint=7717_int64;written%operator_fingerprint=9919_int64
  written%operator_structure_fingerprint=1217_int64;written%scope_fingerprint=1811_int64
  written%payload_fingerprint=2027_int64
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
  call publish_rt_dg_hybrid_checkpoint_v4(comm,trim(prefix),written%global_count,written%nocc,&
    written%row_ids,row_owner,written%row_ids,written,.true.,ok,message)
  call require(ok,'v4 shard publication failed: '//trim(message))
  call read_rt_dg_hybrid_checkpoint_v4(comm,trim(prefix),loaded,ok,message)
  call require(ok,'v4 shard reload failed: '//trim(message))
  call require(loaded%global_count==written%global_count.and.loaded%nocc==written%nocc,&
    'v4 manifest dimensions changed')
  call require(loaded%fragment_id==written%fragment_id,'v4 rank-fragment mapping changed')
  call require(all(loaded%row_ids==written%row_ids),'v4 owned rows changed')
  call require(all(loaded%metric_offsets==written%metric_offsets).and.&
    all(loaded%metric_columns==written%metric_columns).and.&
    maxval(abs(loaded%metric_values-written%metric_values))==0d0,'v4 metric CSR changed')
  call require(all(loaded%operator_offsets==written%operator_offsets).and.&
    all(loaded%operator_columns==written%operator_columns).and.&
    maxval(abs(loaded%operator_values-written%operator_values))==0d0,'v4 operator CSR changed')
  call require(maxval(abs(loaded%position_values-written%position_values))==0d0.and.&
    maxval(abs(loaded%local_values-written%local_values))==0d0,'v4 component CSR changed')
  call require(all(loaded%basis_point_offsets==written%basis_point_offsets).and.&
    all(loaded%basis_support_ids==written%basis_support_ids).and.&
    maxval(abs(loaded%basis_support_values-written%basis_support_values))==0d0,'v4 point CSR changed')
  call require(maxval(abs(loaded%initial_occupied_amplitudes-written%initial_occupied_amplitudes))==0d0,&
    'v4 distributed occupied coefficients changed')
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(prefix),'tddft_response',.true.,1,&
    .false.,.false.,.false.,.false.,.false.,[1],[1d-10,1d-10,1d-10,1d-10],state,ok,message)
  call require(ok.and.state%valid.and.state%initial_invariants_valid,&
    'common v4 endpoint did not initialize distributed Hybrid RT: '//trim(message))
  call require(state%startup_metric_defect<=1d-12.and.state%startup_orbital_residual<=1d-12,&
    'common v4 endpoint changed metric orthonormality or stationarity')
  if(rank==0)write(*,'(a,i0)')'PASS distributed-v4 shard manifest ranks=',nproc
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
end program test_rt_dg_hybrid_checkpoint_v4_mpi
