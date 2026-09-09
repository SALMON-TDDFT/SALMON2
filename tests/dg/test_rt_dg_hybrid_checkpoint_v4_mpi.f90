#include "config.h"
program test_rt_dg_hybrid_checkpoint_v4_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_checkpoint_v4,only:s_rt_dg_hybrid_v4_shard,&
    write_rt_dg_hybrid_checkpoint_v4,read_rt_dg_hybrid_checkpoint_v4
  implicit none
  type(s_rt_dg_hybrid_v4_shard)::written,loaded
  integer::comm,rank,nproc,ierr,i,j
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
  allocate(written%row_ids(50),written%metric_offsets(51),written%metric_columns(150),&
    written%metric_values(150),written%operator_offsets(51),written%operator_columns(150),&
    written%operator_values(150),written%kinetic_values(150),written%nonlocal_values(150),&
    written%local_values(150),written%sipg_values(150),written%position_values(3,150),&
    written%grid_ids(100),written%basis_point_offsets(101),written%basis_support_ids(300),&
    written%basis_support_values(300),written%grid_weights(100),written%density(100),&
    written%initial_occupied_amplitudes(50,7),written%occupations(7),written%eigenvalues(7),&
    written%scope_selectors(8),written%xc_types(1),written%acceptance_receipts(8),&
    written%pseudopotential_receipt(6),written%energy_receipt(7))
  do i=1,50
    written%row_ids(i)=int(rank*50+i,int64)
    written%metric_offsets(i)=3*(i-1)+1;written%operator_offsets(i)=3*(i-1)+1
    do j=1,3
      written%metric_columns(3*(i-1)+j)=mod(rank*50+i+j-2,written%global_count)+1
      written%operator_columns(3*(i-1)+j)=written%metric_columns(3*(i-1)+j)
      written%metric_values(3*(i-1)+j)=cmplx(0.5d0/(i+j),0.125d0*j,real64)
      written%operator_values(3*(i-1)+j)=cmplx(-0.25d0/(i+j),0.0625d0*i,real64)
      written%kinetic_values(3*(i-1)+j)=0.25d0*written%operator_values(3*(i-1)+j)
      written%nonlocal_values(3*(i-1)+j)=0.125d0*written%operator_values(3*(i-1)+j)
      written%local_values(3*(i-1)+j)=0.375d0*written%operator_values(3*(i-1)+j)
      written%sipg_values(3*(i-1)+j)=0.25d0*written%operator_values(3*(i-1)+j)
      written%position_values(:,3*(i-1)+j)=[cmplx(0.1d0*i,0d0,real64),&
        cmplx(0.2d0*j,0d0,real64),cmplx(0.3d0*(i+j),0d0,real64)]
    enddo
    do j=1,7
      written%initial_occupied_amplitudes(i,j)=cmplx(0.01d0*i,0.02d0*j,real64)
    enddo
  enddo
  written%metric_offsets(51)=151;written%operator_offsets(51)=151
  do i=1,100
    written%grid_ids(i)=int(rank*100+i,int64);written%grid_weights(i)=0.5d0
    written%density(i)=0.125d0*i
    written%basis_point_offsets(i)=3*(i-1)+1
    do j=1,3
      written%basis_support_ids(3*(i-1)+j)=mod(rank*50+i+j-2,written%global_count)+1
      written%basis_support_values(3*(i-1)+j)=cmplx(0.03d0*i,-0.04d0*j,real64)
    enddo
  enddo
  written%basis_point_offsets(101)=301
  do j=1,7;written%occupations(j)=0.25d0*j;written%eigenvalues(j)=-0.5d0+0.1d0*j;enddo
  written%scope_selectors=[1,1,1,0,0,0,0,0];written%xc_types=[1]
  written%acceptance_receipts=[(0.01d0*i,i=1,8)]
  written%pseudopotential_receipt=[(1d0*i,i=1,6)]
  written%energy_receipt=[(2d0*i,i=1,7)]
  call write_rt_dg_hybrid_checkpoint_v4(comm,trim(prefix),written,ok,message)
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
