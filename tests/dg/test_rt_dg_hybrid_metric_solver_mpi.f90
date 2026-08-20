#include "config.h"
program test_rt_dg_hybrid_metric_solver_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_get_halting_mode,ieee_set_halting_mode,ieee_invalid,ieee_divide_by_zero
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use rt_dg_hybrid_metric_solver,only:solve_rt_dg_hybrid_metric
  implicit none
  integer,parameter::n=5,nrhs=2
  integer::comm,rank,nproc,ierr,nowned,row,i,j,k,pos,iterations
  integer(int64),allocatable::ids(:)
  complex(real64),allocatable::rhs(:,:),solution(:,:),dense_rhs(:,:),reference(:,:)
  complex(real64)::dense(n,n),work_matrix(n,n)
  integer::pivots(n),info
  type(s_dg_hybrid_sparse_metric)::metric
  integer(int64)::workspace,fingerprint,reference_fingerprint
  real(real64)::residual,defect
  logical::ok
  logical::halt_invalid,halt_zero
  character(256)::message
  external::zgesv
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  dense=(0d0,0d0)
  do i=1,n;dense(i,i)=2d0+0.1d0*i;enddo
  do i=1,n-1
    dense(i,i+1)=cmplx(-0.25d0,0.04d0,real64)
    dense(i+1,i)=conjg(dense(i,i+1))
  enddo
  call distribute_metric(dense,metric)
  nowned=size(metric%owned_row_ids);allocate(rhs(nowned,nrhs),dense_rhs(n,nrhs),reference(n,nrhs))
  dense_rhs(:,1)=[(1d0,0.2d0),(-0.3d0,0.4d0),(0.5d0,-0.1d0),(0.7d0,0.3d0),(-0.2d0,0.1d0)]
  dense_rhs(:,2)=[(0.2d0,-0.5d0),(0.4d0,0.1d0),(-0.6d0,0.2d0),(0.1d0,-0.3d0),(0.8d0,0d0)]
  do i=1,nowned;rhs(i,:)=dense_rhs(int(metric%owned_row_ids(i)),:);enddo
  work_matrix=dense;reference=dense_rhs
  call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
  call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
  call zgesv(n,nrhs,work_matrix,n,pivots,reference,n,info)
  call ieee_set_halting_mode(ieee_invalid,halt_invalid);call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero)
  call require(info==0,'dense metric oracle failed')
  call solve_rt_dg_hybrid_metric(comm,metric,rhs,1d-11,100,solution,iterations,residual,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  defect=0d0
  do i=1,nowned
    defect=max(defect,maxval(abs(solution(i,:)-reference(int(metric%owned_row_ids(i)),:))))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<2d-10.and.residual<1d-11,'matrix-free metric solve differs from dense LAPACK')
  call require(iterations>0.and.iterations<=100.and.workspace>0_int64,'solver receipts are invalid')

  call solve_rt_dg_hybrid_metric(comm,metric,rhs,1d-14,1,solution,iterations,residual,workspace,fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(solution),'iteration-cap failure did not cleanly reject')
  metric%values=0d0
  call solve_rt_dg_hybrid_metric(comm,metric,rhs,1d-11,100,solution,iterations,residual,workspace,fingerprint,ok,message)
  call require(.not.ok,'singular metric was accepted by solver')
  call distribute_metric(dense,metric)
  if(nproc>1)then
    call solve_rt_dg_hybrid_metric(comm,metric,rhs,merge(1d-11,2d-11,rank==0),100,solution,iterations,residual,&
      workspace,fingerprint,ok,message)
    call require(.not.ok,'rank-disagreeing solver tolerance was accepted')
  endif
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_METRIC_SOLVER ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid metric solver on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine distribute_metric(matrix,result)
    complex(real64),intent(in)::matrix(:,:)
    type(s_dg_hybrid_sparse_metric),intent(out)::result
    integer::nnz,column
    nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
    if(allocated(ids))deallocate(ids)
    allocate(ids(nowned));pos=0
    do row=n,1,-1;if(mod(row-1,nproc)==rank)then;pos=pos+1;ids(pos)=row;endif;enddo
    nnz=0
    do pos=1,nowned;row=int(ids(pos));do column=1,n;if(abs(matrix(row,column))>0d0)nnz=nnz+1;enddo;enddo
    result%valid=.true.;result%global_count=n;result%numerical_rank=n;result%condition_estimate=2d0
    result%maximum_value=maxval(abs(matrix));result%max_row_nnz=3;result%fingerprint=9191_int64
    allocate(result%active_rows(n),result%packet_ids(n),result%owned_row_ids(nowned),&
      result%row_offsets(nowned+1),result%column_ids(nnz),result%values(nnz))
    result%active_rows=.true.;result%packet_ids=1;result%owned_row_ids=ids;result%row_offsets(1)=1;k=0
    do pos=1,nowned
      row=int(ids(pos))
      do column=1,n
        if(abs(matrix(row,column))==0d0)cycle
        k=k+1;result%column_ids(k)=column;result%values(k)=matrix(row,column)
      enddo
      result%row_offsets(pos+1)=k+1
    enddo
  end subroutine distribute_metric
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_rt_dg_hybrid_metric_solver_mpi
