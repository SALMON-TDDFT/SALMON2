program test_dg_dc_local_basis_sipg_mpi
  use mpi
  use dg_dc_local_basis_ground_state
  implicit none
  complex(8) :: local_value(2,1),local_normal(2,1),neighbor_value(2,1),neighbor_normal(2,1)
  complex(8) :: hll(1,1),hln(1,1),hnl(1,1),hnn(1,1),zero11(1,1)
  complex(8) :: volume_left(1,1),volume_right(1,1),hamiltonian(2,2),dc_hamiltonian(2,2)
  complex(8) :: unequal_local_value(2,1),unequal_local_normal(2,1)
  complex(8) :: unequal_neighbor_value(2,2),unequal_neighbor_normal(2,2)
  complex(8) :: unequal_ll(1,1),unequal_ln(1,2),unequal_nl(2,1),unequal_nn(2,2)
  complex(8) :: unequal_volume_left(1,1),unequal_volume_right(2,2),unequal_h(3,3),unequal_dc(3,3)
  integer :: ierr,rank,nproc
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'

  local_value(:,1)=[cmplx(1d0+rank,0.1d0*rank,8),cmplx(0.5d0,0.2d0,8)]
  local_normal(:,1)=[cmplx(0.2d0,-0.1d0,8),cmplx(-0.3d0,0.05d0,8)]
  call MPI_Sendrecv(local_value,2,MPI_DOUBLE_COMPLEX,1-rank,31,neighbor_value,2,MPI_DOUBLE_COMPLEX, &
    1-rank,31,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Sendrecv(local_normal,2,MPI_DOUBLE_COMPLEX,1-rank,32,neighbor_normal,2,MPI_DOUBLE_COMPLEX, &
    1-rank,32,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  call assemble_dg_dc_local_basis_sipg_pair(local_value,local_normal,neighbor_value,neighbor_normal, &
    0.5d0,0.25d0,4d0,.true.,rank==0,.false.,[1,0,0],[1,0,0],2-rank,2-rank, &
    hll,hln,hnl,hnn,ok,message)
  call require(ok,'physical SIPG pair assembly')
  call require(rank/=0 .or. maxval(abs(hll-conjg(transpose(hll))))<1d-12, &
    'left diagonal block Hermitian')
  call require(rank/=0 .or. maxval(abs(hnn-conjg(transpose(hnn))))<1d-12, &
    'right diagonal block Hermitian')
  call require(rank/=0 .or. maxval(abs(hnl-conjg(transpose(hln))))<1d-12, &
    'neighbor blocks adjoint')
  call require(rank/=0 .or. maxval(abs([hll(1,1),hln(1,1),hnl(1,1),hnn(1,1)]))>1d-8, &
    'physical interface is nonzero')
  call require(rank==0 .or. maxval(abs([hll(1,1),hln(1,1),hnl(1,1),hnn(1,1)]))==0d0, &
    'noncanonical rank does not double count')
  volume_left(1,1)=2d0
  volume_right(1,1)=3d0
  call compose_dg_dc_local_basis_pair_hamiltonian(volume_left,volume_right,hll,hln,hnl,hnn, &
    0d0,dc_hamiltonian,ok,message)
  call require(ok .and. maxval(abs(dc_hamiltonian-reshape([cmplx(2d0,0d0,8),(0d0,0d0), &
    (0d0,0d0),cmplx(3d0,0d0,8)],[2,2])))<1d-12,'lambda zero reproduces DC volume blocks')
  call compose_dg_dc_local_basis_pair_hamiltonian(volume_left,volume_right,hll,hln,hnl,hnn, &
    1d0,hamiltonian,ok,message)
  call require(ok,'lambda one Hamiltonian composition')
  call require(maxval(abs(hamiltonian-conjg(transpose(hamiltonian))))<1d-12, &
    'assembled Hamiltonian Hermitian')
  call require(rank/=0 .or. maxval(abs(hamiltonian-dc_hamiltonian))>1d-8, &
    'lambda one includes SIPG')

  zero11=(9d0,1d0)
  call assemble_dg_dc_local_basis_sipg_pair(local_value,local_normal,neighbor_value,neighbor_normal, &
    0.5d0,0.25d0,4d0,.false.,.true.,.false.,[0,0,0],[0,0,0],2-rank,2-rank, &
    hll,hln,hnl,hnn,ok,message)
  call require(ok .and. maxval(abs([hll(1,1),hln(1,1),hnl(1,1),hnn(1,1)]))==0d0, &
    'auxiliary nonphysical boundary omitted')
  call assemble_dg_dc_local_basis_sipg_pair(local_value,local_normal,neighbor_value,neighbor_normal, &
    0.5d0,0.25d0,4d0,.true.,.true.,.true.,[1,0,0],[1,0,0],2-rank,2-rank, &
    hll,hln,hnl,hnn,ok,message)
  call require(ok .and. maxval(abs([hll(1,1),hln(1,1),hnl(1,1),hnn(1,1)]))==0d0, &
    'auxiliary periodic wrap omitted')

  call assemble_dg_dc_local_basis_sipg_pair(local_value,local_normal,neighbor_value,neighbor_normal, &
    0.5d0,0.25d0,4d0,.true.,.true.,.false.,[0,1,0],[1,0,0],2-rank,2-rank, &
    hll,hln,hnl,hnn,ok,message)
  call require(.not.ok,'wrong periodic neighbor image rejected')

  call assemble_dg_dc_local_basis_sipg_pair(local_value,local_normal,neighbor_value,neighbor_normal, &
    tiny(1d0),1d0,huge(1d0),.true.,.true.,.false.,[0,0,0],[0,0,0],2-rank,2-rank, &
    hll,hln,hnl,hnn,ok,message)
  call require(.not.ok,'nonfinite accumulated SIPG blocks rejected')

  unequal_local_value(:,1)=[cmplx(1d0,0.2d0,8),cmplx(0.4d0,-0.1d0,8)]
  unequal_local_normal(:,1)=[cmplx(0.3d0,-0.2d0,8),cmplx(-0.1d0,0.05d0,8)]
  unequal_neighbor_value(:,1)=[cmplx(0.7d0,0.1d0,8),cmplx(-0.2d0,0.3d0,8)]
  unequal_neighbor_value(:,2)=[cmplx(-0.4d0,0.2d0,8),cmplx(0.9d0,-0.15d0,8)]
  unequal_neighbor_normal(:,1)=[cmplx(-0.2d0,0.1d0,8),cmplx(0.25d0,0.05d0,8)]
  unequal_neighbor_normal(:,2)=[cmplx(0.15d0,-0.05d0,8),cmplx(-0.35d0,0.2d0,8)]
  call assemble_dg_dc_local_basis_sipg_pair(unequal_local_value,unequal_local_normal, &
    unequal_neighbor_value,unequal_neighbor_normal,0.5d0,0.25d0,4d0,.true.,.true.,.false., &
    [0,0,0],[0,0,0],2-rank,2-rank,unequal_ll,unequal_ln,unequal_nl,unequal_nn,ok,message)
  call require(ok,'unequal local-basis SIPG assembly')
  call require(maxval(abs(unequal_ll-conjg(transpose(unequal_ll))))<1d-12, &
    'unequal left diagonal Hermitian')
  call require(maxval(abs(unequal_nn-conjg(transpose(unequal_nn))))<1d-12, &
    'unequal right diagonal Hermitian')
  call require(maxval(abs(unequal_nl-conjg(transpose(unequal_ln))))<1d-12, &
    'unequal rectangular blocks adjoint')
  unequal_volume_left(1,1)=2d0
  unequal_volume_right=reshape([cmplx(3d0,0d0,8),(0d0,0d0),(0d0,0d0),cmplx(4d0,0d0,8)],[2,2])
  call compose_dg_dc_local_basis_pair_hamiltonian(unequal_volume_left,unequal_volume_right, &
    unequal_ll,unequal_ln,unequal_nl,unequal_nn,0d0,unequal_dc,ok,message)
  call require(ok .and. maxval(abs(unequal_dc(1,2:3)))==0d0 .and. &
    maxval(abs(unequal_dc(2:3,1)))==0d0,'unequal lambda zero block diagonal')
  call compose_dg_dc_local_basis_pair_hamiltonian(unequal_volume_left,unequal_volume_right, &
    unequal_ll,unequal_ln,unequal_nl,unequal_nn,1d0,unequal_h,ok,message)
  call require(ok .and. maxval(abs(unequal_h-conjg(transpose(unequal_h))))<1d-12, &
    'unequal lambda one Hamiltonian Hermitian')

  if(rank==0) print '(a)','PASS DG-DC local-basis SIPG pair assembly'
  call MPI_Finalize(ierr)

contains

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: global
    integer :: mpi_error
    call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if(.not.global) then
      if(rank==0) write(*,'(a,1x,a)') 'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    end if
  end subroutine require
end program test_dg_dc_local_basis_sipg_mpi
