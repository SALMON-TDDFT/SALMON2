program test_dg_hybrid_terminal_operator_guard_mpi
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use mpi
  use dg_hybrid_terminal_refinement, only: s_dg_hybrid_terminal_operator_guard, &
    initialize_dg_hybrid_terminal_operator_guard, validate_dg_hybrid_terminal_operator_guard
  implicit none

  integer :: ierr,rank,nproc,mutation
  integer,parameter :: nrow=2,ncol=3
  integer :: generations(nrow+1),owners(nrow+1)
  integer(int64) :: fixed_payload_fingerprint
  complex(real64) :: metric(nrow,ncol),kinetic(nrow,ncol),nonlocal(nrow,ncol),sipg(nrow,ncol)
  real(real64) :: seed_density(4)
  type(s_dg_hybrid_terminal_operator_guard) :: guard
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call reset_payload
  call initialize_dg_hybrid_terminal_operator_guard(MPI_COMM_WORLD,metric,kinetic,nonlocal,sipg,&
    generations,owners,fixed_payload_fingerprint,seed_density,guard,ok,message)
  call require(ok,'valid immutable payload was rejected: '//trim(message))
  call validate_dg_hybrid_terminal_operator_guard(MPI_COMM_WORLD,metric,kinetic,nonlocal,sipg,&
    generations,owners,fixed_payload_fingerprint,seed_density,guard,ok,message)
  call require(ok,'unchanged immutable payload was rejected: '//trim(message))

  do mutation=1,8
    call reset_payload
    if(rank==0)then
      select case(mutation)
      case(1);metric(1,1)=cmplx(nearest(real(metric(1,1)),1.0_real64),aimag(metric(1,1)),real64)
      case(2);kinetic(1,1)=cmplx(nearest(real(kinetic(1,1)),1.0_real64),aimag(kinetic(1,1)),real64)
      case(3);nonlocal(1,1)=cmplx(real(nonlocal(1,1)),nearest(aimag(nonlocal(1,1)),1.0_real64),real64)
      case(4);sipg(1,1)=cmplx(nearest(real(sipg(1,1)),1.0_real64),aimag(sipg(1,1)),real64)
      case(5);generations(1)=generations(1)+1
      case(6);owners(1)=owners(1)+1
      case(7);fixed_payload_fingerprint=ieor(fixed_payload_fingerprint,1_int64)
      case(8);seed_density(1)=nearest(seed_density(1),1.0_real64)
      end select
    endif
    call validate_dg_hybrid_terminal_operator_guard(MPI_COMM_WORLD,metric,kinetic,nonlocal,sipg,&
      generations,owners,fixed_payload_fingerprint,seed_density,guard,ok,message)
    call require(.not.ok,'immutable component mutation was accepted')
  enddo

  if(rank==0)write(*,'(a,i0,a)')'PASS terminal operator guard on ',nproc,' ranks'
  call MPI_Finalize(ierr)

contains

  subroutine reset_payload
    integer :: i,j
    do j=1,ncol
      do i=1,nrow
        metric(i,j)=cmplx(0.1_real64*i+0.01_real64*j,0.001_real64*(i+j),real64)
        kinetic(i,j)=cmplx(0.2_real64*i+0.02_real64*j,0.002_real64*(i+j),real64)
        nonlocal(i,j)=cmplx(0.3_real64*i+0.03_real64*j,0.003_real64*(i+j),real64)
        sipg(i,j)=cmplx(0.4_real64*i+0.04_real64*j,0.004_real64*(i+j),real64)
      enddo
    enddo
    generations=[11+rank,12+rank,13+rank]
    owners=[rank,modulo(rank+1,max(1,nproc)),modulo(rank+2,max(1,nproc))]
    fixed_payload_fingerprint=7001_int64
    seed_density=[0.11_real64,0.22_real64,0.33_real64,0.44_real64]+real(rank,real64)
  end subroutine reset_payload

  subroutine require(condition,why)
    logical,intent(in) :: condition
    character(*),intent(in) :: why
    if(.not.condition)then
      write(0,'(a,i0,2a)')'rank ',rank,': ',trim(why)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_terminal_operator_guard_mpi
