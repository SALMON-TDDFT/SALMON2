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
  complex(8), allocatable :: exchanged_minus_value(:,:),exchanged_minus_normal(:,:)
  complex(8), allocatable :: exchanged_plus_value(:,:),exchanged_plus_normal(:,:)
  complex(8), allocatable :: neighbor_minus_value(:,:),neighbor_minus_normal(:,:)
  complex(8), allocatable :: neighbor_plus_value(:,:),neighbor_plus_normal(:,:)
  complex(8) :: trace_basis(3,2,2,1),packed_value(4,1),packed_normal(4,1)
  type(s_dg_dc_local_basis_layout) :: row_layout
  complex(8), allocatable :: interface_row(:,:),global_interface(:,:),row_local_block(:,:),row_neighbor_block(:,:)
  integer :: row_counts(2),row_displacements(2)
  integer :: origins8(3,8),sizes8(3,8),face_neighbors(3,2),face_shifts(3,2)
  integer :: fragment,axis,side_index,ix8,iy8,iz8,expected_neighbor,expected_shift
  complex(8), allocatable :: production_basis(:,:,:,:),production_rows(:,:),production_global(:,:)
  integer :: production_origins(3,2),production_sizes(3,2)
  integer :: ierr,rank,nproc,i
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'

  do fragment=1,8
    ix8=modulo(fragment-1,2)
    iy8=modulo((fragment-1)/2,2)
    iz8=(fragment-1)/4
    origins8(:,fragment)=[ix8,iy8,iz8]
    sizes8(:,fragment)=1
  end do
  do fragment=1,8
    call build_dg_dc_six_face_neighbors(fragment,origins8,sizes8,[2,2,2],face_neighbors,face_shifts,ok,message)
    call require(ok,'six-face topology construction')
    do axis=1,3
      do side_index=1,2
        ix8=origins8(1,fragment); iy8=origins8(2,fragment); iz8=origins8(3,fragment)
        select case(axis)
        case(1)
          ix8=1-ix8
        case(2)
          iy8=1-iy8
        case(3)
          iz8=1-iz8
        end select
        expected_neighbor=1+ix8+2*iy8+4*iz8
        expected_shift=0
        if(side_index==1 .and. origins8(axis,fragment)==0) expected_shift=-1
        if(side_index==2 .and. origins8(axis,fragment)==1) expected_shift=1
        call require(face_neighbors(axis,side_index)==expected_neighbor .and. &
          face_shifts(axis,side_index)==expected_shift,'signed periodic face neighbor')
        do i=1,3
          if(i/=axis) call require(origins8(i,face_neighbors(axis,side_index))==origins8(i,fragment), &
            'edge and corner fragments excluded from SIPG topology')
        end do
      end do
    end do
  end do

  do i=1,3
    trace_basis(i,:,:,1)=cmplx(real(i-1,8),0d0,8)
  end do
  call pack_dg_dc_local_basis_face_trace(trace_basis,1,-1,[1d0,-1d0],packed_value,packed_normal,ok,message)
  call require(ok .and. maxval(abs(packed_value))<1d-14 .and. &
    maxval(abs(packed_normal+1d0))<1d-14,'minus face value and outward normal derivative')
  call pack_dg_dc_local_basis_face_trace(trace_basis,1,1,[1d0,-1d0],packed_value,packed_normal,ok,message)
  call require(ok .and. maxval(abs(packed_value-2d0))<1d-14 .and. &
    maxval(abs(packed_normal-1d0))<1d-14,'plus face value and outward normal derivative')

  allocate(exchanged_minus_value(2,rank+1),exchanged_minus_normal(2,rank+1))
  allocate(exchanged_plus_value(2,rank+1),exchanged_plus_normal(2,rank+1))
  allocate(neighbor_minus_value(2,2-rank),neighbor_minus_normal(2,2-rank))
  allocate(neighbor_plus_value(2,2-rank),neighbor_plus_normal(2,2-rank))
  exchanged_minus_value=cmplx(real(10*rank+reshape([(i,i=1,2*(rank+1))],[2,rank+1]),8), &
    0.1d0*real(rank,8),8)
  exchanged_plus_value=exchanged_minus_value+100d0
  exchanged_minus_normal=-0.25d0*exchanged_minus_value
  exchanged_plus_normal=-0.25d0*exchanged_plus_value
  call exchange_dg_dc_local_basis_axis_traces(exchanged_minus_value,exchanged_minus_normal, &
    exchanged_plus_value,exchanged_plus_normal,1-rank,1-rank,77,MPI_COMM_WORLD, &
    neighbor_minus_value,neighbor_minus_normal,neighbor_plus_value,neighbor_plus_normal,ok,message)
  call require(ok,'unequal local-basis axis trace exchange')
  call require(size(neighbor_minus_value,2)==2-rank .and. &
    maxval(abs(neighbor_minus_normal+0.25d0*neighbor_minus_value))<1d-14 .and. &
    maxval(abs(neighbor_plus_normal+0.25d0*neighbor_plus_value))<1d-14, &
    'unequal neighbor trace shape and normal retained')
  deallocate(exchanged_minus_value,exchanged_minus_normal,exchanged_plus_value,exchanged_plus_normal, &
    neighbor_minus_value,neighbor_minus_normal,neighbor_plus_value,neighbor_plus_normal)

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

  call initialize_dg_dc_local_basis_layout(row_layout,rank+1,rank+1,2,71_8,81_8,MPI_COMM_WORLD,ok,message)
  call require(ok,'unequal interface row layout')
  allocate(production_basis(3,2,2,rank+1),production_rows(rank+1,3),production_global(3,3))
  do i=1,rank+1
    production_basis(:,:,:,i)=cmplx(real(rank+1,8)+0.1d0*real(i,8),0.05d0*real(i,8),8)
    production_basis(2,:,:,i)=production_basis(2,:,:,i)+0.2d0*i
    production_basis(3,:,:,i)=production_basis(3,:,:,i)+0.5d0*i
  end do
  production_origins=reshape([0,0,0,3,0,0],[3,2])
  production_sizes=reshape([3,2,2,3,2,2],[3,2])
  call assemble_dg_dc_local_basis_interface_rows(row_layout,production_basis,production_origins, &
    production_sizes,[6,2,2],[0.4d0,0.5d0,0.5d0],8d0,MPI_COMM_WORLD,production_rows,ok,message)
  call require(ok,'production six-face interface row assembly')
  row_counts=[1,2]; row_displacements=[0,1]
  do i=1,3
    call MPI_Allgatherv(production_rows(:,i),size(production_rows,1),MPI_DOUBLE_COMPLEX,production_global(:,i), &
      row_counts,row_displacements,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  end do
  call require(maxval(abs(production_global-conjg(transpose(production_global))))<1d-11, &
    'production interface rows are globally Hermitian')
  deallocate(production_basis,production_rows,production_global)
  allocate(interface_row(rank+1,3),global_interface(3,3))
  if(rank==0) then
    allocate(row_local_block(1,1),row_neighbor_block(1,2))
    row_local_block=unequal_ll
    row_neighbor_block=unequal_ln
  else
    allocate(row_local_block(2,2),row_neighbor_block(2,1))
    row_local_block=unequal_nn
    row_neighbor_block=unequal_nl
  end if
  interface_row=(0d0,0d0)
  call accumulate_dg_dc_local_basis_interface_rows(row_layout,2-rank,row_local_block,row_neighbor_block, &
    interface_row,ok,message)
  call require(ok,'unequal interface blocks inserted into owned rows')
  row_counts=[1,2]
  row_displacements=[0,1]
  do i=1,3
    call MPI_Allgatherv(interface_row(:,i),size(interface_row,1),MPI_DOUBLE_COMPLEX,global_interface(:,i), &
      row_counts,row_displacements,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  end do
  call require(maxval(abs(global_interface-conjg(transpose(global_interface))))<1d-12, &
    'distributed unequal interface rows form a Hermitian matrix')
  deallocate(interface_row,global_interface,row_local_block,row_neighbor_block)

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
