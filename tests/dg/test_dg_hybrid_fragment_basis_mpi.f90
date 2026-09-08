#include "config.h"
program test_dg_hybrid_fragment_basis_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis,build_dg_hybrid_fragment_basis
  implicit none
  integer::comm,rank,nproc,ierr,fragment_id,p,nwf,npw,iw,ip
  integer(int64),allocatable::wf_ids(:),pw_ids(:)
  complex(real64),allocatable::wf_values(:,:),pw_values(:,:)
  type(s_dg_hybrid_fragment_basis)::basis
  integer(int64)::fingerprints(2),combined
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  do fragment_id=1,2
    nwf=0;npw=0
    do p=1,4
      if(mod(p-1,nproc)/=rank)cycle
      if(p<=2)then;nwf=nwf+1;else;npw=npw+1;endif
    enddo
    allocate(wf_ids(nwf),pw_ids(npw),wf_values(2,nwf),pw_values(2,npw));iw=0;ip=0
    do p=1,4
      if(mod(p-1,nproc)/=rank)cycle
      if(p<=2)then
        iw=iw+1;wf_ids(iw)=int(4*(fragment_id-1)+p,int64)
        wf_values(:,iw)=[cmplx(real(p,real64),real(fragment_id,real64),real64),&
          cmplx(-real(p,real64),0.5d0*fragment_id,real64)]
      else
        ip=ip+1;pw_ids(ip)=int(4*(fragment_id-1)+p,int64)
        pw_values(:,ip)=[cmplx(real(p,real64),-real(fragment_id,real64),real64),&
          cmplx(0.25d0*p,real(fragment_id,real64),real64)]
      endif
    enddo
    call build_dg_hybrid_fragment_basis(comm,fragment_id,wf_ids,wf_values,pw_ids,pw_values,&
      2,3,basis,ok,message)
    call require(ok,trim(message))
    call require(basis%fragment_id==fragment_id.and.basis%generation==1,'fragment metadata mismatch')
    call require(size(basis%global_ids)==nwf+npw,'persistent catalog is not rank-local')
    call require(all(basis%sector(:nwf)==1).and.all(basis%sector(nwf+1:)==2),'WF/PW ordering mismatch')
    call require(size(basis%buffer_values,1)==2.and.size(basis%buffer_values,2)==nwf+npw,&
      'buffer storage shape mismatch')
    if(nwf>0)call require(all(basis%buffer_values(:,:nwf)==wf_values),'WF values changed')
    if(npw>0)call require(all(basis%buffer_values(:,nwf+1:)==pw_values),'PW values changed')
    fingerprints(fragment_id)=basis%provenance_fingerprint
    deallocate(wf_ids,pw_ids,wf_values,pw_values)
  enddo

  ! Duplicate ownership must fail collectively.
  allocate(wf_ids(1),pw_ids(1),wf_values(1,1),pw_values(1,1))
  wf_ids=91_int64;pw_ids=91_int64;wf_values=(1d0,0d0);pw_values=(2d0,0d0)
  call build_dg_hybrid_fragment_basis(comm,3,wf_ids,wf_values,pw_ids,pw_values,1,2,basis,ok,message)
  call require(.not.ok,'duplicate global IDs were accepted')
  pw_ids=92_int64
  call build_dg_hybrid_fragment_basis(comm,3,wf_ids,wf_values,pw_ids,pw_values,3,2,basis,ok,message)
  call require(.not.ok,'projector support larger than buffer was accepted')
  deallocate(wf_ids,pw_ids,wf_values,pw_values)

  combined=ieor(fingerprints(1),ishftc(fingerprints(2),17))
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_FRAGMENT_BASIS ranks=',nproc,' fingerprint=',combined
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment basis on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_fragment_basis_mpi
