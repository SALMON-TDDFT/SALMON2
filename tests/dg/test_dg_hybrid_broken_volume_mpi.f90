#include "config.h"
program test_dg_hybrid_broken_volume_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows,&
    assemble_dg_hybrid_local_potential_rows,assemble_dg_hybrid_exact_nonlocal_rows,&
    collect_dg_hybrid_exact_projector_overlaps
  implicit none
  integer,parameter::nbasis=5,npoint=5
  integer::comm,rank,nproc,ierr,i,j,p,nlocal,nrow,ip,ir
  integer::basis_fragment(nbasis),point_fragment(npoint)
  integer,allocatable::local_point_fragment(:)
  integer(int64),allocatable::point_ids(:),row_ids(:)
  real(real64),allocatable::weights(:),potential(:)
  complex(real64),allocatable::values(:,:),gradients(:,:,:),kinetic(:,:),local_rows(:,:)
  complex(real64),allocatable::local_only_rows(:,:)
  complex(real64),allocatable::nonlocal_rows(:,:)
  complex(real64),allocatable::partial_projector_overlap(:,:),collected_projector_overlap(:,:)
  integer(int64),allocatable::collected_projector_ids(:)
  integer,allocatable::collected_projector_owner(:)
  real(real64),allocatable::collected_projector_strength(:)
  complex(real64)::expected_t(nbasis,nbasis),expected_v(nbasis,nbasis),term
  complex(real64)::projector_overlap(nbasis,2),expected_nl(nbasis,nbasis)
  integer(int64)::projector_ids(2)
  integer::projector_owner(2)
  real(real64)::projector_strength(2),nonlocal_diagnostics(2)
  real(real64)::diagnostics(4)
  logical::ok,checks_ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  basis_fragment=[1,1,2,2,2];point_fragment=[1,1,2,2,2]
  nlocal=count([(mod(p-1,nproc)==rank,p=1,npoint)])
  nrow=count([(mod(i-1,nproc)==rank,i=1,nbasis)])
  allocate(point_ids(nlocal),local_point_fragment(nlocal),weights(nlocal),potential(nlocal),values(nbasis,nlocal),&
    gradients(3,nbasis,nlocal),row_ids(nrow))
  ip=0
  do p=1,npoint
    if(mod(p-1,nproc)/=rank)cycle
    ip=ip+1;point_ids(ip)=int(p,int64);local_point_fragment(ip)=point_fragment(p);weights(ip)=0.2d0+0.1d0*p
    potential(ip)=-0.4d0+0.07d0*p
    do i=1,nbasis
      values(i,ip)=cmplx(0.1d0*i+0.03d0*p,-0.02d0*i*p,real64)
      gradients(1,i,ip)=cmplx(0.04d0*i*p,0.01d0*i,real64)
      gradients(2,i,ip)=cmplx(-0.03d0*i,0.02d0*p,real64)
      gradients(3,i,ip)=cmplx(0.01d0*(i+p),-0.015d0*i,real64)
      if(basis_fragment(i)/=point_fragment(p))then
        values(i,ip)=cmplx(9d0+i,2d0+p,real64)
        gradients(:,i,ip)=cmplx(7d0+i,3d0+p,real64)
      endif
    enddo
  enddo
  ir=0
  do i=1,nbasis
    if(mod(i-1,nproc)/=rank)cycle
    ir=ir+1;row_ids(ir)=int(i,int64)
  enddo

  call assemble_dg_hybrid_broken_volume_rows(comm,nbasis,row_ids,basis_fragment,point_ids,&
    local_point_fragment,weights,values,gradients,potential,kinetic,local_rows,diagnostics,ok,message)
  call require(ok,trim(message))
  expected_t=(0d0,0d0);expected_v=(0d0,0d0)
  do p=1,npoint;do i=1,nbasis;do j=1,nbasis
    if(basis_fragment(i)/=point_fragment(p).or.basis_fragment(j)/=point_fragment(p))cycle
    term=sum(conjg([cmplx(0.04d0*i*p,0.01d0*i,real64),cmplx(-0.03d0*i,0.02d0*p,real64),&
      cmplx(0.01d0*(i+p),-0.015d0*i,real64)])*&
      [cmplx(0.04d0*j*p,0.01d0*j,real64),cmplx(-0.03d0*j,0.02d0*p,real64),&
      cmplx(0.01d0*(j+p),-0.015d0*j,real64)])
    expected_t(i,j)=expected_t(i,j)+0.5d0*(0.2d0+0.1d0*p)*term
    expected_v(i,j)=expected_v(i,j)+(0.2d0+0.1d0*p)*(-0.4d0+0.07d0*p)*&
      conjg(cmplx(0.1d0*i+0.03d0*p,-0.02d0*i*p,real64))*&
      cmplx(0.1d0*j+0.03d0*p,-0.02d0*j*p,real64)
  enddo;enddo;enddo
  checks_ok=.true.
  do ir=1,nrow
    i=int(row_ids(ir))
    checks_ok=checks_ok.and.maxval(abs(kinetic(ir,:)-expected_t(i,:)))<2d-13
    checks_ok=checks_ok.and.maxval(abs(local_rows(ir,:)-expected_v(i,:)))<2d-13
    do j=1,nbasis
      if(basis_fragment(i)/=basis_fragment(j))&
        checks_ok=checks_ok.and.abs(kinetic(ir,j))+abs(local_rows(ir,j))<1d-14
    enddo
  enddo
  call require(checks_ok,'broken-volume rows differ from the analytic fragment reference')
  call require(all(diagnostics>=0d0),'invalid broken-volume diagnostics')
  call assemble_dg_hybrid_local_potential_rows(comm,nbasis,row_ids,basis_fragment,point_ids,&
    local_point_fragment,weights,values,potential,local_only_rows,diagnostics(1:2),ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(local_only_rows-local_rows))<2d-13,&
    'local-potential-only rows differ from the combined assembler')
  allocate(partial_projector_overlap(nbasis,2));partial_projector_overlap=(0d0,0d0)
  do ip=1,nlocal
    p=int(point_ids(ip))
    do i=1,nbasis
      partial_projector_overlap(i,1)=partial_projector_overlap(i,1)+weights(ip)*&
        conjg(cmplx(0.2d0+0.01d0*p,-0.03d0*p,real64))*&
        cmplx(0.1d0*i+0.03d0*p,-0.02d0*i*p,real64)
      partial_projector_overlap(i,2)=partial_projector_overlap(i,2)+weights(ip)*&
        conjg(cmplx(-0.1d0+0.02d0*p,0.04d0*p,real64))*&
        cmplx(0.1d0*i+0.03d0*p,-0.02d0*i*p,real64)
    enddo
  enddo
  call collect_dg_hybrid_exact_projector_overlaps(comm,nbasis,[5,6],[1,2],[0.7d0,-0.25d0],&
    partial_projector_overlap,collected_projector_ids,collected_projector_owner,&
    collected_projector_strength,collected_projector_overlap,ok,message)
  call require(ok,trim(message))
  call require(size(collected_projector_ids)==2.and.&
    all(collected_projector_owner==[0,mod(1,nproc)]),&
    'cross-fragment projectors do not have one canonical owner')
  projector_overlap=collected_projector_overlap
  projector_ids=[1_int64,2_int64]
  projector_owner=collected_projector_owner
  projector_strength=collected_projector_strength
  call assemble_dg_hybrid_exact_nonlocal_rows(comm,nbasis,row_ids,projector_ids,projector_owner,&
    projector_strength,projector_overlap,nonlocal_rows,nonlocal_diagnostics,ok,message)
  call require(ok,trim(message))
  expected_nl=(0d0,0d0)
  do p=1,2;do i=1,nbasis;do j=1,nbasis
    expected_nl(i,j)=expected_nl(i,j)+projector_strength(p)*&
      conjg(projector_overlap(i,p))*projector_overlap(j,p)
  enddo;enddo;enddo
  checks_ok=.true.
  do ir=1,nrow
    checks_ok=checks_ok.and.maxval(abs(nonlocal_rows(ir,:)-expected_nl(int(row_ids(ir)),:)))<2d-13
  enddo
  checks_ok=checks_ok.and.abs(expected_nl(1,3))>1d-8
  call require(checks_ok,'crossing nonlocal projector was not accumulated exactly once')
  call require(all(nonlocal_diagnostics>=0d0),'invalid nonlocal diagnostics')
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid broken volume on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_broken_volume_mpi
