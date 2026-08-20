#include "config.h"
program test_rt_dg_hybrid_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_checkpoint,only:write_rt_dg_hybrid_checkpoint,read_rt_dg_hybrid_checkpoint
  implicit none
  integer,parameter::n=4
  integer::comm,rank,nproc,ierr,i,mode_length
  character(256)::mode,path,message
  type(s_dg_hybrid_sparse_metric)::metric
  type(s_dg_hybrid_sparse_operators)::operators
  complex(real64),allocatable::coefficients(:)
  integer(int64)::payload_fingerprint,expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
    expected_complement,expected_position,expected_operator
  real(real64)::observable,metric_observable,energy_observable,position_observable(3)
  logical::ok
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_command_argument(1,mode,length=mode_length);call get_command_argument(1,mode)
  call get_command_argument(2,path)
  if(trim(mode)=='write'.or.trim(mode)=='write_incomplete')then
    call construct_state(metric,operators,coefficients)
    if(trim(mode)=='write_incomplete')metric%packet_ids(n)=0
    call write_rt_dg_hybrid_checkpoint(comm,trim(path),6001_int64,metric,operators,coefficients,7001_int64,&
      payload_fingerprint,ok,message)
    if(trim(mode)=='write')then
      call require(ok,trim(message))
    else
      call require(.not.ok,'incomplete packet checkpoint was accepted')
    endif
  else
    expected_catalog=6001_int64;expected_state=7001_int64;expected_selection=101_int64
    expected_window=102_int64;expected_packet=103_int64;expected_complement=104_int64
    expected_position=105_int64;expected_operator=8181_int64
    select case(trim(mode))
    case('read_stale');expected_catalog=6002_int64
    case('read_stale_selection');expected_selection=999_int64
    case('read_stale_window');expected_window=999_int64
    case('read_stale_packet');expected_packet=999_int64
    case('read_stale_state');expected_state=999_int64
    case('read_stale_complement');expected_complement=999_int64
    case('read_stale_position');expected_position=999_int64
    case('read_stale_operator');expected_operator=999_int64
    case('read_rank_stale');
      if(rank==0)expected_operator=999_int64
    end select
    call read_rt_dg_hybrid_checkpoint(comm,trim(path),expected_catalog,expected_state,expected_selection,expected_window,&
      expected_packet,expected_complement,9191_int64,expected_position,expected_operator,metric,operators,coefficients,&
      payload_fingerprint,ok,message)
    if(trim(mode)=='read')then
      call require(ok,trim(message));observable=0d0
      do i=1,size(coefficients);observable=observable+abs(coefficients(i))**2;enddo
      call MPI_Allreduce(MPI_IN_PLACE,observable,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      call require(abs(observable-1.95d0)<1d-13,'restart coefficient observable differs')
      call require(metric%fingerprint==9191_int64.and.operators%fingerprint==8181_int64,&
        'restart provenance differs')
      call require(operators%window_fingerprint==102_int64.and.operators%packet_fingerprint==103_int64.and.&
        operators%complement_fingerprint==104_int64,'restart operator provenance was dropped')
      call restored_observables(metric,operators,coefficients,metric_observable,energy_observable,position_observable)
      call require(abs(metric_observable-1.9524d0)<1d-13.and.abs(energy_observable-0.8d0)<1d-13.and.&
        maxval(abs(position_observable-[0.4d0,-0.2d0,0.12d0]))<1d-13,&
        'restart S/H/Z observables differ')
    else
      call require(.not.ok,'stale or corrupt hybrid checkpoint was accepted')
      call require(.not.allocated(coefficients).and..not.allocated(metric%owned_row_ids).and.&
        .not.allocated(metric%values).and..not.allocated(operators%owned_row_ids).and.&
        .not.allocated(operators%hamiltonian_values),'rejected checkpoint retained output storage')
    endif
  endif
  if(rank==0.and.trim(mode)=='read')then
    write(*,'(a,i0,a,i0)')'HYBRID_CHECKPOINT ranks=',nproc,' fingerprint=',payload_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid checkpoint on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine construct_state(distributed_metric,distributed_operators,owned_coefficients)
    type(s_dg_hybrid_sparse_metric),intent(out)::distributed_metric
    type(s_dg_hybrid_sparse_operators),intent(out)::distributed_operators
    complex(real64),allocatable,intent(out)::owned_coefficients(:)
    complex(real64)::s(n,n),h(n,n),z(3,n,n),global_coefficients(n)
    integer::row,column,position,edge,nowned
    s=reshape([(1d0,0d0),(0.1d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0.1d0,0d0),(1.1d0,0d0),(0.05d0,0d0),(0d0,0d0),&
      (0d0,0d0),(0.05d0,0d0),(0.9d0,0d0),(0.08d0,0d0),&
      (0d0,0d0),(0d0,0d0),(0.08d0,0d0),(1.2d0,0d0)],[n,n])
    h=(0d0,0d0);z=(0d0,0d0)
    do row=1,n
      h(row,row)=0.2d0*row;z(1,row,row)=0.1d0*row;z(2,row,row)=-0.05d0*row;z(3,row,row)=0.03d0*row
    enddo
    global_coefficients=[(1d0,0.2d0),(-0.4d0,0.1d0),(0.5d0,-0.3d0),(0.6d0,0.2d0)]
    nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
    allocate(distributed_metric%owned_row_ids(nowned),distributed_metric%row_offsets(nowned+1),&
      distributed_metric%column_ids(nowned*n),distributed_metric%values(nowned*n),&
      distributed_metric%active_rows(n),distributed_metric%packet_ids(n),owned_coefficients(nowned))
    allocate(distributed_operators%owned_row_ids(nowned),distributed_operators%row_offsets(nowned+1),&
      distributed_operators%column_ids(nowned*n),distributed_operators%metric_values(nowned*n),&
      distributed_operators%hamiltonian_values(nowned*n),distributed_operators%position_values(3,nowned*n))
    position=0;edge=0;distributed_metric%row_offsets(1)=1
    do row=n,1,-1
      if(mod(row-1,nproc)/=rank)cycle
      position=position+1;distributed_metric%owned_row_ids(position)=row;owned_coefficients(position)=global_coefficients(row)
      do column=1,n
        edge=edge+1;distributed_metric%column_ids(edge)=column;distributed_metric%values(edge)=s(row,column)
        distributed_operators%column_ids(edge)=column;distributed_operators%metric_values(edge)=s(row,column)
        distributed_operators%hamiltonian_values(edge)=h(row,column);distributed_operators%position_values(:,edge)=z(:,row,column)
      enddo
      distributed_metric%row_offsets(position+1)=edge+1
    enddo
    distributed_metric%valid=.true.;distributed_metric%global_count=n;distributed_metric%numerical_rank=n
    distributed_metric%max_row_nnz=n;distributed_metric%maximum_value=1.2d0;distributed_metric%condition_estimate=2d0
    distributed_metric%fingerprint=9191_int64;distributed_metric%active_rows=.true.;distributed_metric%packet_ids=[1,1,2,2]
    distributed_operators%valid=.true.;distributed_operators%global_count=n
    distributed_operators%owned_row_ids=distributed_metric%owned_row_ids
    distributed_operators%row_offsets=distributed_metric%row_offsets
    distributed_operators%selection_fingerprint=101_int64;distributed_operators%window_fingerprint=102_int64
    distributed_operators%packet_fingerprint=103_int64;distributed_operators%complement_fingerprint=104_int64
    distributed_operators%metric_fingerprint=9191_int64;distributed_operators%position_convention_fingerprint=105_int64
    distributed_operators%fingerprint=8181_int64
  end subroutine construct_state
  subroutine restored_observables(distributed_metric,distributed_operators,owned_coefficients,&
      metric_value,hamiltonian_value,position_value)
    type(s_dg_hybrid_sparse_metric),intent(in)::distributed_metric
    type(s_dg_hybrid_sparse_operators),intent(in)::distributed_operators
    complex(real64),intent(in)::owned_coefficients(:)
    real(real64),intent(out)::metric_value,hamiltonian_value,position_value(3)
    complex(real64)::global_coefficients(n),local_metric,local_hamiltonian,local_position(3),applied
    integer::local_row,edge,component
    global_coefficients=(0d0,0d0)
    do local_row=1,size(owned_coefficients)
      global_coefficients(int(distributed_metric%owned_row_ids(local_row)))=owned_coefficients(local_row)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    local_metric=(0d0,0d0);local_hamiltonian=(0d0,0d0);local_position=(0d0,0d0)
    do local_row=1,size(owned_coefficients)
      applied=(0d0,0d0)
      do edge=distributed_metric%row_offsets(local_row),distributed_metric%row_offsets(local_row+1)-1
        applied=applied+distributed_metric%values(edge)*global_coefficients(distributed_metric%column_ids(edge))
      enddo
      local_metric=local_metric+conjg(owned_coefficients(local_row))*applied
      applied=(0d0,0d0)
      do edge=distributed_operators%row_offsets(local_row),distributed_operators%row_offsets(local_row+1)-1
        applied=applied+distributed_operators%hamiltonian_values(edge)*&
          global_coefficients(distributed_operators%column_ids(edge))
      enddo
      local_hamiltonian=local_hamiltonian+conjg(owned_coefficients(local_row))*applied
      do component=1,3
        applied=(0d0,0d0)
        do edge=distributed_operators%row_offsets(local_row),distributed_operators%row_offsets(local_row+1)-1
          applied=applied+distributed_operators%position_values(component,edge)*&
            global_coefficients(distributed_operators%column_ids(edge))
        enddo
        local_position(component)=local_position(component)+conjg(owned_coefficients(local_row))*applied
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_metric,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,local_hamiltonian,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,local_position,3,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    metric_value=real(local_metric);hamiltonian_value=real(local_hamiltonian);position_value=real(local_position)
  end subroutine restored_observables
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_rt_dg_hybrid_checkpoint_mpi
