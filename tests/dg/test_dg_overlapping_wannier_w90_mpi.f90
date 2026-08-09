program test_dg_overlapping_wannier_w90_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_w90,only:estimate_dg_w90_coordinator_bytes,&
    validate_dg_w90_result,setup_dg_w90_gamma_library,run_dg_w90_gamma_library
  implicit none
  integer::ierr,rank
  complex(8)::transform(2,2)
  real(8)::centers(3,2),spreads(2),spread(3)
  integer(8)::bytes
  logical::ok
  character(256)::message
#ifdef USE_WANNIER90
  integer::nntot
  integer,allocatable::nncell(:,:)
  complex(8),allocatable::m_matrix(:,:,:),a_matrix(:,:),library_transform(:,:)
  real(8),allocatable::library_centers(:,:),library_spreads(:)
  real(8)::lattice(3,3),reciprocal(3,3),atoms_cart(3,1),library_spread(3),eigenvalues(1)
  character(2)::atom_symbols(1)
#endif
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0
  centers=reshape([0.1d0,0.2d0,0.3d0,0.6d0,0.2d0,0.3d0],[3,2])
  spreads=[0.4d0,0.5d0];spread=[0.9d0,0.2d0,0.7d0]
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(ok,trim(message))
  call estimate_dg_w90_coordinator_bytes(384,384,12,1,bytes,ok,message)
  call require(ok.and.bytes>0_8,'finite Si64 Wannier90 byte estimate')
  call estimate_dg_w90_coordinator_bytes(huge(0),huge(0),12,1,bytes,ok,message)
  call require(.not.ok,'Wannier90 byte estimate rejects integer overflow')
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
  allocate(m_matrix(1,1,nntot),a_matrix(1,1));m_matrix=(1d0,0d0);a_matrix=(1d0,0d0)
  call run_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d6,1d-10,library_transform,&
    library_centers,library_spreads,library_spread,ok,message)
  call require(ok,trim(message))
  call require(abs(abs(library_transform(1,1))-1d0)<1d-10,&
    'one-band Wannier90 library returns a unitary Gamma transform')
  m_matrix(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
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
