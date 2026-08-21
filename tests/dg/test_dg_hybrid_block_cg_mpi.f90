#include "config.h"
program test_dg_hybrid_block_cg_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_block_cg,only:solve_dg_hybrid_block_cg
  implicit none
  integer,parameter::n=6,m=2
  integer::comm,rank,nproc,ierr,nowned,row,i,j,position,iterations
  integer(int64),allocatable::row_ids(:)
  complex(real64),allocatable::initial(:,:),coefficients(:,:)
  complex(real64)::h(n,n),s(n,n),global_input(n,2*m)
  real(real64)::eigenvalues(m),residual
  integer(int64)::workspace,fingerprint,reference_fingerprint
  logical::ok
  character(64)::stop_reason
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  h=(0d0,0d0);s=(0d0,0d0)
  do i=1,n;h(i,i)=0.2d0*i;s(i,i)=1d0;enddo
  do i=1,n-1;h(i,i+1)=cmplx(-0.025d0,0.01d0,real64);h(i+1,i)=conjg(h(i,i+1));enddo
  nowned=count([(mod(row-1,nproc)==rank,row=1,n)]);allocate(row_ids(nowned),initial(nowned,m));position=0
  do row=n,1,-1
    if(mod(row-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=row
    initial(position,1)=cmplx(1d0/(row+1),0.03d0*row,real64)
    initial(position,2)=cmplx(0.2d0*row,-1d0/(row+2),real64)
  enddo
  call solve_dg_hybrid_block_cg(comm,n,row_ids,initial,apply_h,apply_s,1d-2,1d-11,24,&
    coefficients,eigenvalues,iterations,residual,stop_reason,workspace,fingerprint,ok,message)
  if(rank==0.and..not.ok)write(*,'(*(g0))')'loose block CG iterations=',iterations,&
    ' reason=',trim(stop_reason),' residual=',residual,' message=',trim(message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(iterations>0.and.iterations<24.and.residual<1d-3,'adaptive block CG did not stop at its loose outer target')
  initial=coefficients
  call solve_dg_hybrid_block_cg(comm,n,row_ids,initial,apply_h,apply_s,1d-8,1d-11,24,&
    coefficients,eigenvalues,iterations,residual,stop_reason,workspace,fingerprint,ok,message)
  if(rank==0.and..not.ok)write(*,'(*(g0))')'tight block CG iterations=',iterations,&
    ' reason=',trim(stop_reason),' residual=',residual,' message=',trim(message)
  call require(ok.and.residual<2d-9,'tightened block CG did not converge')
  initial=coefficients
  call solve_dg_hybrid_block_cg(comm,n,row_ids,initial,apply_h,apply_s,1d-8,1d-11,24,&
    coefficients,eigenvalues,iterations,residual,stop_reason,workspace,fingerprint,ok,message)
  call require(ok.and.iterations<=2,'converged warm-started block CG was over-iterated')
  do i=1,nowned;initial(i,1)=cmplx(1d0/(int(row_ids(i))+1),0.03d0*int(row_ids(i)),real64);enddo
  call solve_dg_hybrid_block_cg(comm,n,row_ids,initial,apply_h,apply_s,1d-10,1d-12,1,&
    coefficients,eigenvalues,iterations,residual,stop_reason,workspace,fingerprint,ok,message)
  call require(.not.ok.and.trim(stop_reason)=='iteration_cap','block CG iteration cap was not reported')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_BLOCK_CG ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid block CG on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine apply_h(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    call apply_dense(h,input,output);callback_ok=.true.
  end subroutine apply_h
  subroutine apply_s(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    call apply_dense(s,input,output);callback_ok=.true.
  end subroutine apply_s
  subroutine apply_dense(matrix,input,output)
    complex(real64),intent(in)::matrix(:,:),input(:,:);complex(real64),intent(out)::output(:,:)
    integer::column
    global_input(:,1:size(input,2))=(0d0,0d0)
    do i=1,nowned;global_input(int(row_ids(i)),1:size(input,2))=input(i,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_input,n*size(input,2),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,nowned;output(i,:)=matmul(matrix(int(row_ids(i)),:),global_input(:,1:size(input,2)));enddo
  end subroutine apply_dense
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_block_cg_mpi
