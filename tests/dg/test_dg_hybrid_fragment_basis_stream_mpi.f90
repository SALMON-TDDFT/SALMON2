#include "config.h"
program test_dg_hybrid_fragment_basis_stream_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_basis_stream,only:s_dg_hybrid_fragment_basis_stream,&
    initialize_dg_hybrid_fragment_basis_stream,append_dg_hybrid_projected_pw_tile,&
    finalize_dg_hybrid_fragment_basis_stream
  implicit none
  integer,parameter::nw=3,npw=4,npoint=5
  integer::comm,rank,nproc,ierr,fragment_id,i,j
  integer::wf_owner(nw),pw_owner(npw)
  integer(int64),allocatable::point_ids(:)
  complex(real64),allocatable::wannier_buffer(:,:),pw_tile(:,:)
  type(s_dg_hybrid_fragment_basis_stream)::stream
  type(s_dg_hybrid_fragment_basis)::basis
  integer(int64)::workspace,fingerprint
  logical::ok,values_ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  fragment_id=0;if(rank<2)fragment_id=rank+1
  wf_owner=[1,2,1];pw_owner=[2,1,2,1]
  if(fragment_id>0)then
    allocate(point_ids(npoint),wannier_buffer(nw,npoint))
    point_ids=int([(10*fragment_id+i,i=1,npoint)],int64)
    do j=1,nw;do i=1,npoint
      wannier_buffer(j,i)=cmplx(100*j+10*fragment_id+i,-j,real64)
    enddo;enddo
  else
    allocate(point_ids(0),wannier_buffer(nw,0))
  endif
  call initialize_dg_hybrid_fragment_basis_stream(comm,2,fragment_id,point_ids,wannier_buffer,&
    wf_owner,pw_owner,stream,basis,workspace,fingerprint,ok,message)
  call require(ok,'fragment stream initialization failed: '//trim(message))
  allocate(pw_tile(2,size(point_ids)))
  do j=1,2;do i=1,size(point_ids);pw_tile(j,i)=cmplx(1000*j+10*fragment_id+i,j,real64);enddo;enddo
  call append_dg_hybrid_projected_pw_tile(stream,1,pw_tile,basis,ok,message)
  call require(ok,'first projected PW tile append failed: '//trim(message))
  do j=1,2;do i=1,size(point_ids);pw_tile(j,i)=cmplx(1000*(j+2)+10*fragment_id+i,j+2,real64);enddo;enddo
  call append_dg_hybrid_projected_pw_tile(stream,3,pw_tile,basis,ok,message)
  call require(ok,'second projected PW tile append failed: '//trim(message))
  call finalize_dg_hybrid_fragment_basis_stream(comm,stream,basis,workspace,fingerprint,ok,message)
  call require(ok,'fragment stream finalize failed: '//trim(message))
  values_ok=.true.
  if(fragment_id>0)then
    values_ok=values_ok.and.all(basis%buffer_point_ids==point_ids)
    do j=1,size(basis%global_ids)
      if(basis%global_ids(j)<=nw)then
        values_ok=values_ok.and.wf_owner(int(basis%global_ids(j)))==fragment_id
      else
        values_ok=values_ok.and.pw_owner(int(basis%global_ids(j))-nw)==fragment_id
        do i=1,npoint
          values_ok=values_ok.and.abs(basis%buffer_values(i,j)-&
            cmplx(1000*(int(basis%global_ids(j))-nw)+10*fragment_id+i,&
            int(basis%global_ids(j))-nw,real64))<1d-12
        enddo
      endif
    enddo
  else
    values_ok=values_ok.and.size(basis%global_ids)==0
  endif
  call require(values_ok,'streamed fragment basis payload mismatch')
  call require(workspace<=int(2*npoint*16,8),'fragment stream retained more than one PW tile workspace')
  call require(fingerprint/=0_int64,'fragment stream fingerprint is empty')
  call append_dg_hybrid_projected_pw_tile(stream,3,pw_tile,basis,ok,message)
  call require(.not.ok,'duplicate projected PW tile append must fail')
  if(rank==0)write(*,'(a,i0,a,i0)')'FRAGMENT_BASIS_STREAM ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment basis stream on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_fragment_basis_stream_mpi
