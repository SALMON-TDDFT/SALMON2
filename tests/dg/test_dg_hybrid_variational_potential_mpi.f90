#include "config.h"
program test_dg_hybrid_variational_potential_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_variational_potential,only:assemble_dg_hybrid_total_density,&
    combine_dg_hybrid_fragment_local_potential
  implicit none
  integer,parameter::nglobal=7
  integer::comm,rank,nproc,ierr,p,nlocal,nrequest,i
  integer(int64),allocatable::core_ids(:),request_ids(:)
  real(real64),allocatable::core_density(:),total_density(:),hartree(:),xc(:),ionic(:),combined(:)
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(core_ids(nlocal),core_density(nlocal));i=0
  do p=1,nglobal
    if(mod(p-1,nproc)/=rank)cycle
    i=i+1;core_ids(i)=int(p,int64);core_density(i)=0.2d0+0.03d0*p
  enddo
  call assemble_dg_hybrid_total_density(comm,nglobal,core_ids,core_density,total_density,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(total_density-[(0.2d0+0.03d0*p,p=1,nglobal)]))<1d-15,&
    'fragment core density was not assembled on the total grid')

  nrequest=3;allocate(request_ids(nrequest),hartree(nrequest),xc(nrequest),ionic(nrequest),combined(nrequest))
  request_ids=[int(modulo(2*rank,nglobal)+1,int64),int(modulo(2*rank+1,nglobal)+1,int64),&
    int(modulo(2*rank+2,nglobal)+1,int64)]
  do i=1,nrequest
    p=int(request_ids(i));hartree(i)=2d0*total_density(p)+0.01d0*p
    xc(i)=-0.3d0*sqrt(total_density(p));ionic(i)=-0.5d0+0.02d0*p
  enddo
  call combine_dg_hybrid_fragment_local_potential(comm,nglobal,request_ids,hartree,xc,ionic,.true.,&
    combined,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(combined-(hartree+xc+ionic)))<1d-15,&
    'fragment Hartree/XC/ionic local potential combination is incorrect')
  call combine_dg_hybrid_fragment_local_potential(comm,nglobal,request_ids,hartree,xc,ionic,.false.,&
    combined,ok,message)
  call require(.not.ok.and.index(message,'halo')>0,'incomplete semilocal XC halo was accepted')
  if(rank==0)write(*,'(a,i0,a)')'PASS variational potential distribution on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_variational_potential_mpi
