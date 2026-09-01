#include "config.h"
program test_dg_hybrid_localization_first_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_quiet_nan,ieee_value
  use dg_hybrid_localization_first,only:s_dg_hybrid_localization_receipt,&
    prepare_dg_hybrid_localization_first_seed,build_dg_hybrid_localization_receipt
  implicit none
  integer,parameter::npoint=8,nraw=3
  integer::comm,rank,nproc,ierr,nlocal,p,i,j,k,global_point,prepared_rank
  integer(int64),allocatable::point_ids(:)
  real(real64)::pi,x,closure_residual,nan_value
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::raw_seed(:,:),basis(:,:),global_basis(:,:),&
    localized_basis(:,:),global_localized_basis(:,:)
  complex(real64)::gram(nraw,nraw),identity(nraw,nraw)
  complex(real64)::projector(npoint,npoint),localized_projector(npoint,npoint),&
    expected_projector(npoint,npoint)
  complex(real64)::reflected(npoint),projected(npoint),phi_i,phi_j
  complex(real64)::transform(nraw,nraw),rank_loss_transform(nraw,nraw-1)
  real(real64)::centers(3,nraw),spreads(nraw)
  real(real64)::rank_loss_centers(3,nraw-1),rank_loss_spreads(nraw-1)
  integer(int64)::seed_fingerprint,reference_seed_fingerprint
  integer(int64)::reference_transform_fingerprint
  logical::ok
  character(256)::message
  type(s_dg_hybrid_localization_receipt)::receipt

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call require(mod(npoint,nproc)==0,'test requires an even block distribution')
  nlocal=npoint/nproc;pi=acos(-1d0)
  allocate(point_ids(nlocal),weights(nlocal),raw_seed(nraw,nlocal))
  weights=1d0/real(npoint,real64)
  do p=1,nlocal
    global_point=rank+1+(nlocal-p)*nproc
    point_ids(p)=int(global_point,int64)
    x=2d0*pi*real(global_point-1,real64)/real(npoint,real64)
    raw_seed(1,p)=cmplx(1d0,0d0,real64)
    raw_seed(2,p)=raw_seed(1,p)+2d0*exp(cmplx(0d0,x,real64))
    raw_seed(3,p)=cmplx(0.3d0,-0.2d0,real64)*exp(cmplx(0d0,x,real64))+&
      cmplx(1d0,1d0,real64)*exp(cmplx(0d0,2d0*x,real64))
  enddo

  call prepare_dg_hybrid_localization_first_seed(comm,point_ids,raw_seed,weights,1d-12,&
    basis,prepared_rank,seed_fingerprint,ok,message)
  call require(ok,trim(message))
  call require(prepared_rank==nraw.and.all(shape(basis)==[nraw,nlocal]),&
    'localization-first preparation changed the raw occupied+s+p rank')
  reference_seed_fingerprint=seed_fingerprint

  gram=matmul(basis*spread(weights,1,nraw),conjg(transpose(basis)))
  call MPI_Allreduce(MPI_IN_PLACE,gram,size(gram),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  identity=cmplx(0d0,0d0,real64)
  do i=1,nraw;identity(i,i)=cmplx(1d0,0d0,real64);enddo
  call require(maxval(abs(gram-identity))<1d-11,&
    'localization-first seed is not metric orthonormal')

  allocate(global_basis(nraw,npoint));global_basis=cmplx(0d0,0d0,real64)
  do p=1,nlocal;global_basis(:,int(point_ids(p)))=basis(:,p);enddo
  call MPI_Allreduce(MPI_IN_PLACE,global_basis,size(global_basis),MPI_DOUBLE_COMPLEX,&
    MPI_SUM,comm,ierr)
  call require(ierr==MPI_SUCCESS,'prepared seed reconstruction failed')
  projector=cmplx(0d0,0d0,real64)
  expected_projector=cmplx(0d0,0d0,real64)
  do i=1,npoint
    do j=1,npoint
      do k=1,nraw
        projector(i,j)=projector(i,j)+global_basis(k,i)*conjg(global_basis(k,j))/real(npoint,real64)
      enddo
      do k=0,nraw-1
        phi_i=exp(cmplx(0d0,2d0*pi*real(k*(i-1),real64)/real(npoint,real64),real64))
        phi_j=exp(cmplx(0d0,2d0*pi*real(k*(j-1),real64)/real(npoint,real64),real64))
        expected_projector(i,j)=expected_projector(i,j)+&
          phi_i*conjg(phi_j)/real(npoint,real64)
      enddo
    enddo
  enddo
  call require(maxval(abs(projector-expected_projector))<1d-10,&
    'metric preparation changed the weighted raw-span projector')

  do i=1,npoint
    x=2d0*pi*real(i-1,real64)/real(npoint,real64)
    reflected(i)=exp(cmplx(0d0,-x,real64))
  enddo
  projected=matmul(expected_projector,reflected)
  closure_residual=sqrt(sum(abs(reflected-projected)**2)/real(npoint,real64))
  call require(closure_residual>0.9d0,&
    'raw occupied+s+p fixture is accidentally closed under reflection')

  call set_unitary_transform(transform)
  allocate(localized_basis(nraw,nlocal),global_localized_basis(nraw,npoint))
  localized_basis=matmul(transpose(transform),basis)
  global_localized_basis=cmplx(0d0,0d0,real64)
  do p=1,nlocal;global_localized_basis(:,int(point_ids(p)))=localized_basis(:,p);enddo
  call MPI_Allreduce(MPI_IN_PLACE,global_localized_basis,size(global_localized_basis),&
    MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  localized_projector=cmplx(0d0,0d0,real64)
  do i=1,npoint;do j=1,npoint;do k=1,nraw
    localized_projector(i,j)=localized_projector(i,j)+&
      global_localized_basis(k,i)*conjg(global_localized_basis(k,j))/real(npoint,real64)
  enddo;enddo;enddo
  call require(maxval(abs(localized_projector-projector))<1d-10,&
    'retaining every localized column changed the construction-space projector')
  centers=reshape([0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,&
    0.7d0,0.8d0,0.9d0],[3,nraw])
  spreads=[1d200,2d200,3d200]
  call build_dg_hybrid_localization_receipt(comm,nraw,nraw,reference_seed_fingerprint,&
    transform,centers,spreads,.true.,17,1d-12,receipt,ok,message)
  call require(ok,trim(message))
  call require(receipt%valid.and..not.receipt%symmetry_constrained.and.receipt%converged.and.&
    receipt%raw_rank==nraw.and.receipt%retained_rank==nraw.and.receipt%iterations==17,&
    'localization-first receipt lost its fixed-rank unconstrained contract')
  call require(receipt%spread_min==1d200.and.receipt%spread_max==3d200.and.&
    abs(receipt%spread_mean/2d200-1d0)<1d-14.and.&
    abs(receipt%spread_total/6d200-1d0)<1d-14,&
    'finite large Wannier spreads were clipped or rejected')
  call require(receipt%transform_unitarity_defect<1d-12.and.&
    receipt%seed_fingerprint==reference_seed_fingerprint.and.&
    receipt%transform_fingerprint/=0_int64,'localization-first fingerprints are invalid')
  reference_transform_fingerprint=receipt%transform_fingerprint

  nan_value=ieee_value(0d0,ieee_quiet_nan)
  if(rank==0)centers(1,1)=nan_value
  call build_dg_hybrid_localization_receipt(comm,nraw,nraw,reference_seed_fingerprint,&
    transform,centers,spreads,.true.,17,1d-12,receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,'nonfinite Wannier center was accepted')
  if(rank==0)centers(1,1)=0.1d0

  if(rank==0)spreads(1)=nan_value
  call build_dg_hybrid_localization_receipt(comm,nraw,nraw,reference_seed_fingerprint,&
    transform,centers,spreads,.true.,17,1d-12,receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,'nonfinite Wannier spread was accepted')
  if(rank==0)spreads=[1d200,2d200,3d200]

  if(rank==0)transform(1,1)=2d0*transform(1,1)
  call build_dg_hybrid_localization_receipt(comm,nraw,nraw,reference_seed_fingerprint,&
    transform,centers,spreads,.true.,17,1d-12,receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,'nonunitary Wannier transform was accepted')
  call set_unitary_transform(transform)

  call build_dg_hybrid_localization_receipt(comm,nraw,nraw,reference_seed_fingerprint,&
    transform,centers,spreads,rank/=0,17,1d-12,receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,'nonconverged Wannier result was accepted')

  rank_loss_transform=transform(:,1:nraw-1)
  rank_loss_centers=centers(:,1:nraw-1);rank_loss_spreads=spreads(1:nraw-1)
  call build_dg_hybrid_localization_receipt(comm,nraw,nraw-1,reference_seed_fingerprint,&
    rank_loss_transform,rank_loss_centers,rank_loss_spreads,.true.,17,1d-12,&
    receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,'rank-losing Wannier transform was accepted')
  call build_dg_hybrid_localization_receipt(comm,nraw,merge(nraw-1,nraw,rank==0),&
    reference_seed_fingerprint,transform,centers,spreads,.true.,17,1d-12,&
    receipt,ok,message)
  call require(.not.ok.and..not.receipt%valid,&
    'rank-local retained-rank loss was not rejected collectively')

  raw_seed(3,:)=raw_seed(2,:)
  call prepare_dg_hybrid_localization_first_seed(comm,point_ids,raw_seed,weights,1d-12,&
    basis,prepared_rank,seed_fingerprint,ok,message)
  call require(.not.ok,'rank-deficient raw localization seed was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0,a,i0,a,es25.16e3)')'LOCALIZATION_FIRST ranks=',nproc,&
      ' seed=',reference_seed_fingerprint,' transform=',reference_transform_fingerprint,&
      ' spread_total=',6d200
    write(*,'(a,i0,a)')'PASS localization-first contract on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine set_unitary_transform(value)
    complex(real64),intent(out)::value(nraw,nraw)
    value=cmplx(0d0,0d0,real64)
    value(1,1)=1d0/sqrt(2d0);value(1,2)=cmplx(0d0,1d0,real64)/sqrt(2d0)
    value(2,1)=cmplx(0d0,1d0,real64)/sqrt(2d0);value(2,2)=1d0/sqrt(2d0)
    value(3,3)=cmplx(1d0,0d0,real64)
  end subroutine set_unitary_transform

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_localization_first_mpi
