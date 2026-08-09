program test_dg_overlapping_wannier_w90_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_w90,only:estimate_dg_w90_coordinator_bytes,&
    validate_dg_w90_result,setup_dg_w90_gamma_library,run_dg_w90_gamma_library,&
    assemble_dg_w90_gamma_matrices,apply_dg_w90_gamma_transform
  implicit none
  integer::ierr,rank,nproc,b,m,n,p,nlocal
  complex(8)::transform(2,2)
  real(8)::centers(3,2),spreads(2),spread(3)
  integer(8)::bytes
  logical::ok,matrix_matches
  character(256)::message
  complex(8),allocatable::local_values(:,:),local_anchors(:,:)
  complex(8),allocatable::assembled_m(:,:,:),assembled_a(:,:)
  real(8),allocatable::local_weights(:),local_fractional(:,:)
  integer::test_nncell(3,2),global_point
  integer(8)::matrix_peak,matrix_estimate
  complex(8)::local_m_reference(2,2,2),local_a_reference(2,2),m_reference(2,2,2),a_reference(2,2),phase
  real(8)::angle
  complex(8),allocatable::gauge_values(:,:),gauge_gradients(:,:,:)
  complex(8)::gauge_transform(2,2)
  real(8)::gauge_centers(3,2)
  integer(8),allocatable::gauge_ids(:)
#ifdef USE_WANNIER90
  integer::nntot
  integer,allocatable::nncell(:,:)
  complex(8),allocatable::m_matrix(:,:,:),a_matrix(:,:),library_transform(:,:)
  real(8),allocatable::library_centers(:,:),library_spreads(:)
  real(8)::lattice(3,3),reciprocal(3,3),atoms_cart(3,1),library_spread(3),eigenvalues(1)
  character(2)::atom_symbols(1)
#endif
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0
  centers=reshape([0.1d0,0.2d0,0.3d0,0.6d0,0.2d0,0.3d0],[3,2])
  spreads=[0.4d0,0.5d0];spread=[0.9d0,0.2d0,0.7d0]
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(ok,trim(message))
  call estimate_dg_w90_coordinator_bytes(384,384,12,1,bytes,ok,message)
  call require(ok.and.bytes>0_8,'finite Si64 Wannier90 byte estimate')
  call estimate_dg_w90_coordinator_bytes(huge(0),huge(0),12,1,bytes,ok,message)
  call require(.not.ok,'Wannier90 byte estimate rejects integer overflow')
  nlocal=8/nproc
  allocate(local_values(2,nlocal),local_anchors(2,nlocal),local_weights(nlocal),local_fractional(3,nlocal))
  do p=1,nlocal
    global_point=rank*nlocal+p
    local_values(:,p)=[cmplx(global_point,0d0,8),cmplx((-1d0)**global_point,0d0,8)]
    local_anchors(:,p)=[cmplx(1d0,0d0,8),cmplx(global_point,0d0,8)]
    local_fractional(:,p)=[real(global_point-1,8)/8d0,0.25d0*mod(global_point-1,4),0d0]
  enddo
  local_weights=0.125d0;test_nncell=reshape([1,0,0,0,1,0],[3,2])
  local_m_reference=(0d0,0d0);local_a_reference=(0d0,0d0)
  do b=1,2;do p=1,nlocal
    angle=-2d0*acos(-1d0)*dot_product(real(test_nncell(:,b),8),local_fractional(:,p))
    phase=cmplx(cos(angle),sin(angle),8)
    do n=1,2;do m=1,2
      local_m_reference(m,n,b)=local_m_reference(m,n,b)+&
        local_weights(p)*conjg(local_values(m,p))*phase*local_values(n,p)
    enddo;enddo
  enddo;enddo
  do p=1,nlocal;do n=1,2;do m=1,2
    local_a_reference(m,n)=local_a_reference(m,n)+&
      local_weights(p)*conjg(local_values(m,p))*local_anchors(n,p)
  enddo;enddo;enddo
  call MPI_Reduce(local_m_reference,m_reference,size(m_reference),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
    MPI_COMM_WORLD,ierr)
  call MPI_Reduce(local_a_reference,a_reference,size(a_reference),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
    MPI_COMM_WORLD,ierr)
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,huge(0_8),assembled_m,assembled_a,matrix_estimate,&
    matrix_peak,ok,message)
  call require(ok.and.matrix_estimate>0_8.and.matrix_peak<=matrix_estimate,trim(message))
  if(rank==0)then
    call require(all(shape(assembled_m)==[2,2,2]).and.all(shape(assembled_a)==[2,2]),&
      'coordinator owns complete Wannier90 M/A matrices')
    matrix_matches=maxval(abs(assembled_m-m_reference))<1d-12.and.&
      maxval(abs(assembled_a-a_reference))<1d-12
    write(*,'(a,1x,es24.16)')'W90_MATRIX_FINGERPRINT',sum(abs(assembled_m))+sum(abs(assembled_a))
  else
    call require(size(assembled_m)==0.and.size(assembled_a)==0,&
      'noncoordinator does not own complete Wannier90 M/A matrices')
    matrix_matches=.true.
  endif
  call require(matrix_matches,'distributed Wannier90 M/A match dense references')
  p=merge(nlocal,max(0,nlocal-1),rank==0)
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values(:,1:p),local_anchors(:,1:p),&
    local_weights(1:p),local_fractional(:,1:p),test_nncell,huge(0_8),assembled_m,assembled_a,&
    matrix_estimate,matrix_peak,ok,message)
  call require(ok,'Wannier90 matrix assembly accepts uneven local point ownership')
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,1_8,assembled_m,assembled_a,matrix_estimate,matrix_peak,&
    ok,message)
  call require(.not.ok,'Wannier90 matrix assembly rejects a low coordinator byte limit')
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,merge(1_8,huge(0_8),rank==0),assembled_m,assembled_a,&
    matrix_estimate,matrix_peak,ok,message)
  call require(.not.ok,'Wannier90 matrix assembly rejects rank-inconsistent byte limits')
  allocate(gauge_values,source=local_values)
  allocate(gauge_gradients(3,2,nlocal));gauge_gradients=(0d0,0d0)
  allocate(gauge_ids(nlocal));gauge_ids=[(int(rank*nlocal+p,8),p=1,nlocal)]
  gauge_transform=reshape([cmplx(0d0,0d0,8),cmplx(-1d0,0d0,8),&
    cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)],[2,2])
  gauge_centers=reshape([0.2d0,0d0,0d0,0.2d0,0d0,0d0],[3,2])
  call apply_dg_w90_gamma_transform(MPI_COMM_WORLD,gauge_ids,gauge_values,gauge_gradients,&
    gauge_transform,gauge_centers,1d-12,ok,message)
  call require(ok.and.maxval(abs(gauge_values(1,:)-local_values(1,:)))<1d-12.and.&
    maxval(abs(gauge_values(2,:)+local_values(2,:)))<1d-12.and.&
    abs(gauge_transform(1,1)-1d0)<1d-12.and.abs(gauge_transform(2,2)+1d0)<1d-12.and.&
    maxval(abs(gauge_transform-reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
      cmplx(0d0,0d0,8),cmplx(-1d0,0d0,8)],[2,2])))<1d-12,trim(message))
  transform(1,1)=2d0
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'nonunitary Wannier90 transform rejection')
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=cmplx(1d0,1d-4,8)
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'complex Gamma Wannier90 gauge rejection')
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0;spread(3)=0.9d0
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'increased Wannier90 gauge-dependent spread rejection')
  spread(3)=0.7d0;centers(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'nonfinite Wannier90 center rejection')
#ifdef USE_WANNIER90
  lattice=0d0;reciprocal=0d0
  lattice(1,1)=10d0;lattice(2,2)=10d0;lattice(3,3)=10d0
  reciprocal(1,1)=2d0*acos(-1d0)/10d0
  reciprocal(2,2)=reciprocal(1,1);reciprocal(3,3)=reciprocal(1,1)
  atoms_cart=0d0;atom_symbols(1)='H ';eigenvalues=0d0
  call setup_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,1,1,nntot,nncell,ok,message)
  call require(ok.and.nntot>0,trim(message))
  if(rank==0)then
    allocate(m_matrix(1,1,nntot),a_matrix(1,1));m_matrix=(1d0,0d0);a_matrix=(1d0,0d0)
  else
    allocate(m_matrix(0,0,0),a_matrix(0,0))
  endif
  call run_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d6,1d-10,library_transform,&
    library_centers,library_spreads,library_spread,ok,message)
  call require(ok,trim(message))
  call require(abs(abs(library_transform(1,1))-1d0)<1d-10,&
    'one-band Wannier90 library returns a unitary Gamma transform')
  if(rank==0)m_matrix(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call run_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d6,1d-10,library_transform,&
    library_centers,library_spreads,library_spread,ok,message)
  call require(.not.ok,'nonfinite Wannier90 M matrix is rejected before library entry')
#endif
  if(rank==0)write(*,'(a)')'PASS Wannier90 MLWF adapter validation'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_bad/=0)error stop label
  end subroutine
end program
