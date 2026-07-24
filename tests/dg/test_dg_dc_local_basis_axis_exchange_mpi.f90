program test_dg_dc_local_basis_axis_exchange_mpi
  use mpi
  use dg_dc_local_basis_ground_state
  implicit none
  complex(8), allocatable :: minus_value(:,:),minus_normal(:,:),plus_value(:,:),plus_normal(:,:)
  complex(8), allocatable :: neighbor_minus_value(:,:),neighbor_minus_normal(:,:)
  complex(8), allocatable :: neighbor_plus_value(:,:),neighbor_plus_normal(:,:)
  integer :: ierr,rank,nproc,minus_rank,plus_rank
  logical :: ok,global_ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=4) error stop 'test requires four ranks'
  minus_rank=modulo(rank-1,nproc)
  plus_rank=modulo(rank+1,nproc)
  allocate(minus_value(3,rank+1),minus_normal(3,rank+1),plus_value(3,rank+1),plus_normal(3,rank+1))
  allocate(neighbor_minus_value(3,minus_rank+1),neighbor_minus_normal(3,minus_rank+1))
  allocate(neighbor_plus_value(3,plus_rank+1),neighbor_plus_normal(3,plus_rank+1))
  minus_value=cmplx(100d0*rank+1d0,0.1d0*rank,8)
  plus_value=cmplx(100d0*rank+2d0,0.1d0*rank,8)
  minus_normal=-0.5d0*minus_value
  plus_normal=-0.5d0*plus_value
  call exchange_dg_dc_local_basis_axis_traces(minus_value,minus_normal,plus_value,plus_normal, &
    minus_rank,plus_rank,901,MPI_COMM_WORLD,neighbor_minus_value,neighbor_minus_normal, &
    neighbor_plus_value,neighbor_plus_normal,ok,message)
  ok=ok .and. maxval(abs(neighbor_minus_value-cmplx(100d0*minus_rank+2d0,0.1d0*minus_rank,8)))<1d-14
  ok=ok .and. maxval(abs(neighbor_plus_value-cmplx(100d0*plus_rank+1d0,0.1d0*plus_rank,8)))<1d-14
  ok=ok .and. maxval(abs(neighbor_minus_normal+0.5d0*neighbor_minus_value))<1d-14
  ok=ok .and. maxval(abs(neighbor_plus_normal+0.5d0*neighbor_plus_value))<1d-14
  call MPI_Allreduce(ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
  if(.not.global_ok) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  if(rank==0) print '(a)','PASS DG-DC local-basis four-rank axis exchange'
  call MPI_Finalize(ierr)
end program test_dg_dc_local_basis_axis_exchange_mpi
