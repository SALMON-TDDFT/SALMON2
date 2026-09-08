#include "config.h"
program test_rt_dg_hybrid_occupied_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_checkpoint,only:write_rt_dg_hybrid_occupied_checkpoint,read_rt_dg_hybrid_occupied_checkpoint
  implicit none
  integer,parameter::n=5,m=2
  integer::comm,rank,nproc,ierr,row,position,nlocal,global_count,coefficient_bad
  integer(int64),allocatable::row_ids(:),read_ids(:)
  complex(real64),allocatable::coefficients(:,:),read_coefficients(:,:)
  real(real64)::occupations(m),eigenvalues(m),receipts(5),read_receipts(5)
  real(real64),allocatable::read_occupations(:),read_eigenvalues(:)
  integer(int64)::fingerprint,provenance(6)
  logical::ok
  character(256)::mode,path,message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD;call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_command_argument(1,mode);call get_command_argument(2,path)
  if(index(trim(mode),'write')==1)then
    nlocal=count([(mod(row-1,nproc)==rank,row=1,n)]);allocate(row_ids(nlocal),coefficients(nlocal,m));position=0
    do row=n,1,-1
      if(mod(row-1,nproc)/=rank)cycle
      position=position+1;row_ids(position)=row
      coefficients(position,1)=cmplx(0.1d0*row,-0.03d0*row,real64)
      coefficients(position,2)=cmplx(-0.04d0*row,0.07d0/(row+1),real64)
    enddo
    occupations=[2d0,1d0];eigenvalues=[-0.7d0,0.2d0];receipts=[1d-9,2d-10,3d-11,4d-12,5d-13]
    provenance=[501_int64,502_int64,503_int64,504_int64,505_int64,506_int64]
    if(trim(mode)=='write_bad_occupation')occupations(2)=-1d0
    if(trim(mode)=='write_stale_scf')receipts(1)=2d-8
    if(trim(mode)=='write_incomplete'.and.rank==0.and.nlocal>0)row_ids(1)=merge(2_int64,1_int64,row_ids(1)==1_int64)
    call write_rt_dg_hybrid_occupied_checkpoint(comm,trim(path),n,row_ids,coefficients,occupations,eigenvalues,&
      101_int64,202_int64,provenance,303_int64,404_int64,receipts,1d-8,fingerprint,ok,message)
    if(trim(mode)=='write')then;call require(ok,trim(message))
    else;call require(.not.ok,'invalid occupied state was checkpointed');endif
  else
    provenance=[501_int64,502_int64,503_int64,504_int64,505_int64,506_int64]
    if(trim(mode)=='stale_provenance')provenance(4)=999_int64
    occupations=[2d0,1d0];if(trim(mode)=='changed_occupation')occupations(2)=0.5d0
    call read_rt_dg_hybrid_occupied_checkpoint(comm,trim(path),101_int64,202_int64,303_int64,&
      merge(999_int64,404_int64,trim(mode)=='stale'),provenance,occupations,1d-8,global_count,read_ids,read_coefficients,read_occupations,&
      read_eigenvalues,read_receipts,fingerprint,ok,message)
    if(trim(mode)=='read')then
      call require(ok,trim(message));call require(global_count==n.and.size(read_coefficients,2)==m,'occupied checkpoint shape differs')
      call require(maxval(abs(read_occupations-[2d0,1d0]))<1d-15.and.&
        maxval(abs(read_eigenvalues-[-0.7d0,0.2d0]))<1d-15,'occupied metadata differs')
      coefficient_bad=0
      do position=1,size(read_ids)
        row=int(read_ids(position))
        if(abs(read_coefficients(position,1)-cmplx(0.1d0*row,-0.03d0*row,real64))>=1d-14)coefficient_bad=1
        if(abs(read_coefficients(position,2)-cmplx(-0.04d0*row,0.07d0/(row+1),real64))>=1d-14)coefficient_bad=1
      enddo
      call require(coefficient_bad==0,'occupied coefficients differ')
      call require(maxval(abs(read_receipts-[1d-9,2d-10,3d-11,4d-12,5d-13]))<1d-20,'SCF receipts differ')
      if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_OCCUPIED_CHECKPOINT ranks=',nproc,' fingerprint=',fingerprint
      if(rank==0)write(*,'(a,i0,a)')'PASS occupied checkpoint on ',nproc,' ranks'
    else
      call require(.not.ok.and..not.allocated(read_coefficients),'stale or old checkpoint was accepted')
    endif
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine
end program
