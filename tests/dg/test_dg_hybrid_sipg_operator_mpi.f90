#include "config.h"
program test_dg_hybrid_sipg_operator_mpi
  use mpi, only: MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_Finalize, MPI_Init, MPI_INTEGER, MPI_MAX, &
    MPI_SUCCESS, MPI_SUM
  use,intrinsic::iso_fortran_env,only:real64
  use dg_hybrid_sipg_operator,only:s_dg_hybrid_sipg_face_operator,assemble_dg_hybrid_sipg_face,&
    scale_dg_hybrid_sipg_faces
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_coupling_envelope,&
    build_dg_hybrid_operator_envelope,set_dg_hybrid_operator_hermitian_edge
  use dg_nodal_sipg,only:s_dg_nodal_sipg_action,evaluate_dg_nodal_sipg_face
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange
  implicit none
  integer::icomm,id_rank,nproc,ierr,i,j,local_count,global_count
  complex(real64)::value_minus(2),value_plus(1),derivative_minus(2),derivative_plus(1)
  complex(real64)::jump(3),average_derivative(3),reference(3,3)
  real(real64)::h,weight,eta,lambda
  type(s_dg_hybrid_sipg_face_operator)::face,face2,scaled
  type(s_dg_nodal_sipg_action)::nodal_action
  integer::info
  type(s_dg_hybrid_coupling_envelope)::envelope
  type(s_rt_dg_sparse_exchange)::exchange_plan
  integer::basis_edges(2,3),nonlocal_edges(2,2),face_edges(2,2),position
  integer,allocatable::local_basis(:,:),local_nonlocal(:,:),local_face(:,:)
  integer(8),allocatable::owned_rows(:)
  integer(8)::structure_fingerprint
  logical::ok,edge_ok
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  value_minus=[cmplx(1d0,0.2d0,real64),cmplx(-0.3d0,0.4d0,real64)]
  value_plus=[cmplx(0.7d0,-0.1d0,real64)]
  derivative_minus=[cmplx(0.5d0,-0.2d0,real64),cmplx(-0.1d0,0.3d0,real64)]
  derivative_plus=[cmplx(-0.4d0,0.25d0,real64)]
  h=0.8d0;weight=1.25d0;eta=6d0
  call assemble_dg_hybrid_sipg_face(icomm,17,0,[0,0,0],[1,2],[3],value_minus,derivative_minus,value_plus,&
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
  do j=1,3
    call evaluate_dg_nodal_sipg_face(value_minus(min(j,2))*merge(1d0,0d0,j<=2),&
      value_plus(1)*merge(1d0,0d0,j==3),derivative_minus(min(j,2))*merge(1d0,0d0,j<=2),&
      derivative_plus(1)*merge(1d0,0d0,j==3),h,weight,eta,nodal_action,info)
    call require(info==0,'nodal SIPG reference rejected a projected basis trace')
    do i=1,3
      reference(i,j)=0.5d0*(conjg(merge(value_minus(min(i,2)),(0d0,0d0),i<=2))*nodal_action%total_value(1)+&
        conjg(merge(value_plus(1),(0d0,0d0),i==3))*nodal_action%total_value(2)+&
        conjg(merge(derivative_minus(min(i,2)),(0d0,0d0),i<=2))*nodal_action%total_normal(1)+&
        conjg(merge(derivative_plus(1),(0d0,0d0),i==3))*nodal_action%total_normal(2))
    enddo
  enddo
  call require(maxval(abs(face%total-reference))<1d-13,'projected SIPG disagrees with nodal SIPG action')
  if(nproc>1)then
    call assemble_dg_hybrid_sipg_face(icomm,17,id_rank,[0,0,0],[1,2],[3],value_minus,derivative_minus,value_plus,&
      derivative_plus,h,weight,eta,face2,ok,message)
    call require(.not.ok,'rank-disagreeing SIPG owner was accepted')
  endif
  call assemble_dg_hybrid_sipg_face(icomm,18,0,[1,0,0],[4],[5],value_plus,derivative_plus,value_plus,&
    derivative_plus,h,weight,eta,face2,ok,message)
  call require(ok,trim(message))
  call require(all(face2%periodic_shift==[1,0,0]),'physical periodic face shift was not preserved')
  lambda=0.375d0
  call scale_dg_hybrid_sipg_faces([face,face2],lambda,[lambda,lambda],scaled,ok,message)
  call require(ok.and.all(scaled%global_basis_ids==[1,2,3,4,5]),&
    'global basis IDs from distinct faces were not united')
  call require(maxval(abs(scaled%total(1:3,1:3)-lambda*face%total))<1d-13,&
    'uniform lambda did not scatter the first SIPG face')
  call scale_dg_hybrid_sipg_faces([face],lambda,[0.5d0],scaled,ok,message)
  call require(.not.ok,'face-local lambda was accepted')

  basis_edges=reshape([1,1,2,2,3,3],[2,3])
  nonlocal_edges=reshape([2,3,3,2],[2,2]);face_edges=reshape([1,3,3,1],[2,2])
  owned_rows=pack([1_8,2_8,3_8],[(mod(i-1,nproc)==id_rank,i=1,3)])
  call localize(basis_edges,owned_rows,local_basis);call localize(nonlocal_edges,owned_rows,local_nonlocal)
  call localize(face_edges,owned_rows,local_face)
  call build_dg_hybrid_operator_envelope(icomm,3,owned_rows,local_basis,local_nonlocal,local_face,4,envelope,ok,message)
  call require(ok,trim(message))
  structure_fingerprint=envelope%structure_fingerprint
  local_count=merge(1,0,envelope%find_edge(2,3)>0)
  call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,icomm,ierr)
  call require(global_count==1,'cross-fragment nonlocal projector edge was not assembled exactly once')
  call build_rt_dg_sparse_exchange(icomm,3,structure_fingerprint,owned_rows,envelope%column_ids,&
    exchange_plan,ok,message)
  call require(ok.and.exchange_plan%valid,'real sparse halo schedule was not constructed: '//trim(message))
  position=envelope%find_edge(1,3)
  call require((position>0).eqv.(mod(0,nproc)==id_rank),'row-owned face edge distribution is incorrect')
  call set_dg_hybrid_operator_hermitian_edge(icomm,envelope,2,1,3,cmplx(0.2d0,-0.1d0,real64),ok,message)
  call require(ok,trim(message))
  position=envelope%find_edge(1,3)
  edge_ok=.true.
  if(position>0)edge_ok=abs(envelope%component_values(2,position)-cmplx(0.2d0,-0.1d0,real64))<1d-14
  call require(edge_ok,'allowed operator edge could not change from zero')
  position=envelope%find_edge(3,1)
  edge_ok=.true.
  if(position>0)edge_ok=abs(envelope%component_values(2,position)-cmplx(0.2d0,0.1d0,real64))<1d-14
  call require(edge_ok,'reverse Hermitian edge was not updated')
  call require(envelope%structure_fingerprint==structure_fingerprint,'operator value update rebuilt the graph')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS hybrid SIPG operator on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine localize(global_edges,rows,local_edges)
    integer,intent(in)::global_edges(:,:)
    integer(8),intent(in)::rows(:)
    integer,allocatable,intent(out)::local_edges(:,:)
    logical,allocatable::keep(:)
    integer::q
    allocate(keep(size(global_edges,2)));keep=.false.
    do q=1,size(global_edges,2);keep(q)=any(rows==int(global_edges(1,q),8));enddo
    allocate(local_edges(2,count(keep)));local_edges=reshape(pack(global_edges,spread(keep,1,2)),[2,count(keep)])
  end subroutine localize
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(id_rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_sipg_operator_mpi
