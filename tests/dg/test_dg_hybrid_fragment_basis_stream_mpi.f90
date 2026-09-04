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
  call test_local_fragment_wannier
  if(rank==0)write(*,'(a,i0,a,i0)')'FRAGMENT_BASIS_STREAM ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment basis stream on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine test_local_fragment_wannier
    integer,allocatable::owners(:),pwowners(:)
    complex(real64),allocatable::wf(:,:),tile(:,:)
    integer(int64)::ids(npoint)
    integer::f,a,b,k,total,local_wf,first,generation
    ! Actual production layout: one rank per fragment, variable local WF count.
    f=nproc-rank;local_wf=f;total=nproc*(nproc+1)/2;first=f*(f-1)/2
    allocate(owners(total),pwowners(nproc),wf(local_wf,npoint),tile(nproc,npoint))
    k=0
    do a=1,nproc
      do b=1,a;k=k+1;owners(k)=a;enddo
      pwowners(a)=a
    enddo
    ids=int([(10*f+a,a=1,npoint)],int64)
    do a=1,npoint
      do b=1,local_wf;wf(b,a)=cmplx(100*b+10*f+a,-b,real64);enddo
      do b=1,nproc;tile(b,a)=cmplx(1000*b+10*f+a,b,real64);enddo
    enddo
    generation=7
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc,f,ids,wf,owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=.true.,basis_generation=generation)
    call require(ok,'local WF stream initialization failed: '//trim(message))
    call require(basis%generation==7.and.size(basis%global_ids)==local_wf+1,&
      'local WF stream changed fragment rank or generation')
    call append_dg_hybrid_projected_pw_tile(stream,1,tile,basis,ok,message)
    call require(ok,'local WF projected-PW append failed: '//trim(message))
    call finalize_dg_hybrid_fragment_basis_stream(comm,stream,basis,workspace,fingerprint,ok,message)
    call require(ok,'local WF stream finalization failed: '//trim(message))
    values_ok=all(basis%buffer_point_ids==ids)
    do b=1,local_wf
      values_ok=values_ok.and.basis%global_ids(b)==first+b.and.all(basis%buffer_values(:,b)==wf(b,:))
    enddo
    values_ok=values_ok.and.basis%global_ids(local_wf+1)==total+f.and.&
      all(basis%buffer_values(:,local_wf+1)==tile(f,:))
    call require(values_ok,'local WF/PW stream changed column order, values or global IDs')
    if(rank==0)generation=8
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc,f,ids,wf,owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=.true.,basis_generation=generation)
    call require(.not.ok.and..not.allocated(basis%global_ids),'rank-disagreeing generation accepted')
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc,f,ids,wf,owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=rank/=0,basis_generation=7)
    call require(.not.ok.and..not.allocated(basis%global_ids),'rank-disagreeing local WF mode accepted')
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc,f,ids,wf(:local_wf-1,:),owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=.true.)
    call require(.not.ok.and..not.allocated(basis%global_ids),'missing local WF column accepted')
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc,f,ids,wf,owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=.true.,basis_generation=0)
    call require(.not.ok.and..not.allocated(basis%global_ids),'nonpositive basis generation accepted')
    call initialize_dg_hybrid_fragment_basis_stream(comm,nproc+1,f,ids,wf,owners,pwowners,&
      stream,basis,workspace,fingerprint,ok,message,local_wannier_only=.true.)
    call require(.not.ok.and..not.allocated(basis%global_ids),'unequal fragment and rank counts accepted')
  end subroutine test_local_fragment_wannier

  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_fragment_basis_stream_mpi
