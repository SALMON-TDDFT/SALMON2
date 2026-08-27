#include "config.h"
program test_dg_hybrid_sipg_operator_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use dg_hybrid_sipg_operator,only:s_dg_hybrid_sipg_face_operator,assemble_dg_hybrid_sipg_face,&
    scale_dg_hybrid_sipg_faces
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_coupling_envelope,&
    build_dg_hybrid_operator_envelope,set_dg_hybrid_operator_edge
  implicit none
  integer::comm,rank,nproc,ierr,i,j
  complex(real64)::value_minus(2),value_plus(1),derivative_minus(2),derivative_plus(1)
  complex(real64)::jump(3),average_derivative(3),reference(3,3)
  real(real64)::h,weight,eta,lambda
  type(s_dg_hybrid_sipg_face_operator)::face,scaled
  type(s_dg_hybrid_coupling_envelope)::envelope
  integer::basis_edges(2,3),nonlocal_edges(2,1),face_edges(2,1),position
  integer(8)::structure_fingerprint,schedule_fingerprint
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  value_minus=[cmplx(1d0,0.2d0,real64),cmplx(-0.3d0,0.4d0,real64)]
  value_plus=[cmplx(0.7d0,-0.1d0,real64)]
  derivative_minus=[cmplx(0.5d0,-0.2d0,real64),cmplx(-0.1d0,0.3d0,real64)]
  derivative_plus=[cmplx(-0.4d0,0.25d0,real64)]
  h=0.8d0;weight=1.25d0;eta=6d0
  call assemble_dg_hybrid_sipg_face(comm,17,0,value_minus,derivative_minus,value_plus,&
    derivative_plus,h,weight,eta,face,ok,message)
  call require(ok,trim(message))
  jump=[value_minus,-value_plus]
  average_derivative=0.5d0*[derivative_minus,derivative_plus]
  do j=1,3;do i=1,3
    reference(i,j)=0.5d0*weight*(-conjg(jump(i))*average_derivative(j)-&
      conjg(average_derivative(i))*jump(j)+(eta/h)*conjg(jump(i))*jump(j))
  enddo;enddo
  call require(maxval(abs(face%total-reference))<1d-13,'complete SIPG face matrix mismatch')
  call require(maxval(abs(face%total-conjg(transpose(face%total))))<1d-13,'SIPG face matrix is not Hermitian')
  call require(maxval(abs(face%total(1:2,3)))>1d-12.and.maxval(abs(face%total(3,1:2)))>1d-12,&
    'cross-fragment SIPG blocks are missing')
  call require(maxval(abs(face%physical_penalty-0.5d0*face%raw_penalty))<1d-13,&
    'SIPG penalty kinetic factor was not applied exactly once')
  lambda=0.375d0
  call scale_dg_hybrid_sipg_faces([face],lambda,[lambda],scaled,ok,message)
  call require(ok.and.maxval(abs(scaled%total-lambda*face%total))<1d-13,&
    'uniform lambda did not scale every SIPG block')
  call scale_dg_hybrid_sipg_faces([face],lambda,[0.5d0],scaled,ok,message)
  call require(.not.ok,'face-local lambda was accepted')

  basis_edges=reshape([1,1,2,2,3,3],[2,3])
  nonlocal_edges(:,1)=[2,3];face_edges(:,1)=[1,3]
  call build_dg_hybrid_operator_envelope(3,basis_edges,nonlocal_edges,face_edges,4,envelope,ok,message)
  call require(ok,trim(message))
  structure_fingerprint=envelope%structure_fingerprint
  schedule_fingerprint=envelope%schedule_fingerprint
  position=envelope%find_edge(1,3)
  call require(position>0.and.envelope%find_edge(3,1)>0,'complete reciprocal face edge is absent')
  call require(abs(envelope%component_values(2,position))==0d0,'allowed zero component edge was not explicit')
  call set_dg_hybrid_operator_edge(envelope,2,1,3,cmplx(0.2d0,-0.1d0,real64),ok,message)
  call require(ok.and.abs(envelope%component_values(2,position)-cmplx(0.2d0,-0.1d0,real64))<1d-14,&
    'allowed operator edge could not change from zero')
  call require(envelope%structure_fingerprint==structure_fingerprint.and.&
    envelope%schedule_fingerprint==schedule_fingerprint,'operator value update rebuilt the graph or schedule')
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid SIPG operator on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_sipg_operator_mpi
