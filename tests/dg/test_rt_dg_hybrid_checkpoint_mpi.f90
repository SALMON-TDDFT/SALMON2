#include "config.h"
program test_rt_dg_hybrid_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_checkpoint,only:write_rt_dg_hybrid_checkpoint,read_rt_dg_hybrid_checkpoint,&
    s_rt_dg_hybrid_ground_state_payload,write_rt_dg_hybrid_ground_state_checkpoint,&
    read_rt_dg_hybrid_ground_state_checkpoint,fingerprint_rt_dg_hybrid_component
  implicit none
  integer,parameter::n=4
  integer::comm,rank,nproc,ierr,i,mode_length
  character(256)::mode,path,message
  type(s_dg_hybrid_sparse_metric)::metric
  type(s_dg_hybrid_sparse_operators)::operators
  type(s_rt_dg_hybrid_ground_state_payload)::complete_payload,restored_payload
  complex(real64),allocatable::coefficients(:)
  integer(int64)::payload_fingerprint,expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
    expected_complement,expected_position,expected_operator
  real(real64)::observable,metric_observable,energy_observable,position_observable(3)
  logical::ok
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_command_argument(1,mode,length=mode_length);call get_command_argument(1,mode)
  call get_command_argument(2,path)
  if(trim(mode)=='write_complete'.or.trim(mode)=='write_bad_complete'.or.trim(mode)=='write_incomplete_complete'.or.&
      trim(mode)=='write_interrupted_complete')then
    call construct_complete_payload(complete_payload)
    if(trim(mode)=='write_bad_complete')complete_payload%hamiltonian_rows(1,1)=&
      complete_payload%hamiltonian_rows(1,1)+(1d0,0d0)
    if(trim(mode)=='write_incomplete_complete')deallocate(complete_payload%face_values)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),complete_payload,payload_fingerprint,ok,message,&
      interrupt_after_write=trim(mode)=='write_interrupted_complete')
    if(trim(mode)=='write_complete')then
      call require(ok,trim(message))
    else
      call require(.not.ok,'inconsistent or incomplete complete payload was published')
    endif
  else if(trim(mode)=='read_complete'.or.trim(mode)=='read_complete_corrupt')then
    call read_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),restored_payload,payload_fingerprint,ok,message)
    if(trim(mode)=='read_complete')then
      call require(ok,trim(message))
      call require(restored_payload%final_refresh_complete.and.restored_payload%analysis_complete,&
        'complete checkpoint lost final acceptance receipts')
      call require(all(restored_payload%hamiltonian_rows==restored_payload%kinetic_rows+&
        restored_payload%nonlocal_rows+restored_payload%local_rows+restored_payload%sipg_rows),&
        'complete checkpoint changed Hamiltonian component identity')
      call require(restored_payload%requested_ids(1)==1.and.restored_payload%effective_ids(4)==4.and.&
        restored_payload%added_ids(1)==4.and.restored_payload%closure_parent(1)==1,&
        'complete checkpoint lost requested/effective closure provenance')
      call require(size(restored_payload%basis_values,2)==size(restored_payload%grid_ids).and.&
        size(restored_payload%face_values,2)==size(restored_payload%face_ids).and.&
        size(restored_payload%nonlocal_values,2)==size(restored_payload%nonlocal_ids),&
        'complete checkpoint lost basis, face, or nonlocal payload')
      call require(size(restored_payload%metric_column_ids)/=size(restored_payload%operator_column_ids),&
        'complete checkpoint collapsed independent metric and operator graphs')
    else
      call require(.not.ok,'corrupt complete DG ground-state checkpoint was accepted')
    endif
  else if(trim(mode)=='write_legacy')then
    call write_legacy_checkpoint(trim(path))
  else if(trim(mode)=='write'.or.trim(mode)=='write_incomplete')then
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
      call require(size(operators%metric_values)==size(operators%column_ids).and.&
        maxval(abs(operators%metric_values-[((1d0,0d0),i=1,size(operators%metric_values))]))<0.21d0,&
        'authoritative metric was not projected onto the restored operator graph')
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
  if(rank==0.and.trim(mode)=='read_complete')write(*,'(a,i0,a)')'PASS complete hybrid checkpoint on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine construct_complete_payload(payload)
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::payload
    integer::row,point,local_row,local_point,nrow,npoint
    logical::fingerprint_ok
    payload%valid=.true.;payload%final_refresh_complete=.true.;payload%analysis_complete=.true.
    payload%identity_only=.false.;payload%global_count=n;payload%noccupied=2
    payload%operation_count=2;payload%nonidentity_operation_count=1
    payload%catalog_fingerprint=1101_int64;payload%state_fingerprint=1102_int64
    payload%metric_fingerprint=1103_int64;payload%operator_structure_fingerprint=1104_int64
    payload%operator_value_fingerprint=1105_int64;payload%kinetic_fingerprint=1106_int64
    payload%nonlocal_fingerprint=1107_int64;payload%local_fingerprint=1108_int64
    payload%sipg_fingerprint=1109_int64;payload%basis_fingerprint=1110_int64
    payload%face_fingerprint=1111_int64;payload%dc_seed_fingerprint=1112_int64
    payload%continuation_fingerprint=1113_int64;payload%scope_fingerprint=1114_int64
    payload%pseudopotential_fingerprint=1115_int64;payload%energy_fingerprint=1116_int64
    payload%analysis_fingerprint=1117_int64;payload%selection_fingerprint=1118_int64
    nrow=count([(mod(row-1,nproc)==rank,row=1,n)])
    npoint=count([(mod(point-1,nproc)==rank,point=1,n)])
    allocate(payload%row_ids(nrow),payload%metric_rows(nrow,n),payload%kinetic_rows(nrow,n),&
      payload%nonlocal_rows(nrow,n),payload%local_rows(nrow,n),payload%sipg_rows(nrow,n),&
      payload%hamiltonian_rows(nrow,n),payload%coefficients(nrow,2))
    allocate(payload%metric_row_offsets(nrow+1),payload%metric_column_ids(nrow*n),&
      payload%operator_row_offsets(nrow+1),payload%operator_column_ids(nrow))
    payload%metric_row_offsets=[(1+(row-1)*n,row=1,nrow+1)]
    do row=1,nrow;payload%metric_column_ids((row-1)*n+1:row*n)=[1,2,3,4];enddo
    payload%operator_row_offsets=[(row,row=1,nrow+1)]
    payload%metric_rows=(0d0,0d0);payload%kinetic_rows=(0d0,0d0);payload%nonlocal_rows=(0d0,0d0)
    payload%local_rows=(0d0,0d0);payload%sipg_rows=(0d0,0d0);payload%coefficients=(0d0,0d0)
    local_row=0
    do row=1,n
      if(mod(row-1,nproc)/=rank)cycle
      local_row=local_row+1;payload%row_ids(local_row)=row
      payload%operator_column_ids(local_row)=row
      payload%metric_rows(local_row,row)=cmplx(1d0+0.1d0*row,0d0,real64)
      payload%kinetic_rows(local_row,row)=cmplx(0.2d0*row,0d0,real64)
      payload%nonlocal_rows(local_row,row)=cmplx(-0.03d0*row,0d0,real64)
      payload%local_rows(local_row,row)=cmplx(0.07d0*row,0d0,real64)
      payload%sipg_rows(local_row,row)=cmplx(0.01d0*row,0d0,real64)
      payload%coefficients(local_row,1)=cmplx(0.1d0*row,0.02d0*row,real64)
      payload%coefficients(local_row,2)=cmplx(-0.03d0*row,0.04d0*row,real64)
    enddo
    payload%hamiltonian_rows=payload%kinetic_rows+payload%nonlocal_rows+payload%local_rows+payload%sipg_rows
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%kinetic_rows,&
      payload%kinetic_fingerprint,fingerprint_ok);call require(fingerprint_ok,'kinetic fixture fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%nonlocal_rows,&
      payload%nonlocal_fingerprint,fingerprint_ok);call require(fingerprint_ok,'nonlocal fixture fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%local_rows,&
      payload%local_fingerprint,fingerprint_ok);call require(fingerprint_ok,'local fixture fingerprint failed')
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%sipg_rows,&
      payload%sipg_fingerprint,fingerprint_ok);call require(fingerprint_ok,'SIPG fixture fingerprint failed')
    allocate(payload%grid_ids(npoint),payload%grid_weights(npoint),payload%partition_ids(npoint),&
      payload%basis_values(n,npoint),payload%density(npoint))
    local_point=0
    do point=1,n
      if(mod(point-1,nproc)/=rank)cycle
      local_point=local_point+1;payload%grid_ids(local_point)=100+point;payload%grid_weights(local_point)=0.25d0
      payload%partition_ids(local_point)=1+mod(point,2);payload%density(local_point)=0.5d0+0.01d0*point
      do row=1,n;payload%basis_values(row,local_point)=cmplx(0.01d0*row*point,-0.02d0*row,real64);enddo
    enddo
    allocate(payload%face_ids(1),payload%face_point_ids(1),payload%face_metadata(4,1),payload%face_offsets(2),&
      payload%face_value_offsets(2),payload%face_basis_ids(2),&
      payload%face_normals(3,1),payload%face_weights(1),payload%face_values(4,1),payload%interface_observables(3,1))
    payload%face_ids=200+rank;payload%face_point_ids=300+rank;payload%face_metadata(:,1)=[1,2,0,rank]
    payload%face_offsets=[1,2];payload%face_value_offsets=[1,5];payload%face_basis_ids=[1,2]
    payload%face_normals(:,1)=[1d0,0d0,0d0];payload%face_weights=0.5d0
    payload%face_values(:,1)=[(0.1d0,0.01d0),(0.2d0,0.02d0),(0.3d0,0.03d0),(0.4d0,0.04d0)]
    payload%interface_observables(:,1)=[(0.5d0,0d0),(0.6d0,0d0),(0.7d0,0d0)]
    allocate(payload%nonlocal_ids(1),payload%nonlocal_owner(1),payload%nonlocal_values(2,1))
    payload%nonlocal_ids=400+rank;payload%nonlocal_owner=rank
    payload%nonlocal_values(:,1)=[(0.11d0,0.02d0),(0.12d0,0.03d0)]
    allocate(payload%requested_ids(3),payload%effective_ids(4),payload%added_ids(1),payload%closure_parent(1),&
      payload%closure_reason(1),payload%closure_action(1),payload%scope_selectors(6),payload%xc_types(1))
    payload%requested_ids=[1,2,3];payload%effective_ids=[1,2,3,4];payload%added_ids=4
    payload%closure_parent=1;payload%closure_reason=2;payload%closure_action=2
    payload%scope_selectors=[1,1,0,0,0,0];payload%xc_types=4
    allocate(payload%occupations(2),payload%eigenvalues(2),payload%continuation_receipt(8),&
      payload%pseudopotential_receipt(3),payload%energy_receipt(4))
    payload%occupations=[2d0,2d0];payload%eigenvalues=[-0.5d0,-0.2d0]
    payload%continuation_receipt=[1d0,1d-9,2d-9,3d-9,4d-9,0d0,1d0,1d0]
    payload%pseudopotential_receipt=[1d0,2d0,3d0];payload%energy_receipt=[4d0,5d0,6d0,7d0]
  end subroutine construct_complete_payload

  subroutine write_legacy_checkpoint(checkpoint_path)
    character(*),intent(in)::checkpoint_path
    character(16),parameter::magic='SALMON_DG_HYB01 '
    integer,parameter::version=1
    integer::unit,io_status,row,column,component
    integer::metric_degrees(n),operator_degrees(n),metric_columns(n),operator_columns(1)
    integer::packet_ids(n)
    integer(int64)::fingerprint,bits
    logical::active_rows(n)
    complex(real64)::s(n,n),h,z(3),c(n),operator_metric(1)
    s=reshape([(1d0,0d0),(0.1d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0.1d0,0d0),(1.1d0,0d0),(0.05d0,0d0),(0d0,0d0),&
      (0d0,0d0),(0.05d0,0d0),(0.9d0,0d0),(0.08d0,0d0),&
      (0d0,0d0),(0d0,0d0),(0.08d0,0d0),(1.2d0,0d0)],[n,n])
    c=[(1d0,0.2d0),(-0.4d0,0.1d0),(0.5d0,-0.3d0),(0.6d0,0.2d0)]
    active_rows=.true.;packet_ids=[1,1,2,2];metric_degrees=n;operator_degrees=1
    metric_columns=[(column,column=1,n)]
    fingerprint=6001_int64
    call legacy_hash_int(fingerprint,7001_int64);call legacy_hash_int(fingerprint,9191_int64)
    call legacy_hash_int(fingerprint,8181_int64);call legacy_hash_int(fingerprint,int(n,int64))
    call legacy_hash_int(fingerprint,int(n,int64))
    bits=transfer(2d0,bits);call legacy_hash_int(fingerprint,bits)
    bits=transfer(1.2d0,bits);call legacy_hash_int(fingerprint,bits)
    call legacy_hash_int(fingerprint,101_int64);call legacy_hash_int(fingerprint,102_int64)
    call legacy_hash_int(fingerprint,103_int64);call legacy_hash_int(fingerprint,104_int64)
    call legacy_hash_int(fingerprint,9191_int64);call legacy_hash_int(fingerprint,105_int64)
    do row=1,n
      call legacy_hash_int(fingerprint,1_int64);call legacy_hash_int(fingerprint,int(packet_ids(row),int64))
      call legacy_hash_int(fingerprint,int(n,int64));call legacy_hash_int(fingerprint,1_int64)
    enddo
    do row=1,n
      h=cmplx(0.2d0*row,0d0,real64)
      z=[cmplx(0.1d0*row,0d0,real64),cmplx(-0.05d0*row,0d0,real64),cmplx(0.03d0*row,0d0,real64)]
      operator_columns(1)=row;operator_metric(1)=s(row,row)
      call legacy_hash_int(fingerprint,int(row,int64));call legacy_hash_int(fingerprint,int(n,int64))
      call legacy_hash_int(fingerprint,1_int64)
      do column=1,n
        call legacy_hash_int(fingerprint,int(column,int64));call legacy_hash_complex(fingerprint,s(row,column))
      enddo
      call legacy_hash_int(fingerprint,int(row,int64));call legacy_hash_complex(fingerprint,operator_metric(1))
      call legacy_hash_complex(fingerprint,h)
      do component=1,3;call legacy_hash_complex(fingerprint,z(component));enddo
      call legacy_hash_complex(fingerprint,c(row))
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    if(rank==0)then
      open(newunit=unit,file=checkpoint_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status)
      if(io_status==0)write(unit,iostat=io_status)magic,version,n,n,6001_int64,7001_int64,9191_int64,2d0,1.2d0,&
        101_int64,102_int64,103_int64,104_int64,9191_int64,105_int64,8181_int64,&
        active_rows,packet_ids,metric_degrees,operator_degrees
      if(io_status==0)then
        do row=1,n
          h=cmplx(0.2d0*row,0d0,real64)
          z=[cmplx(0.1d0*row,0d0,real64),cmplx(-0.05d0*row,0d0,real64),cmplx(0.03d0*row,0d0,real64)]
          operator_columns(1)=row;operator_metric(1)=s(row,row)
          write(unit,iostat=io_status)metric_columns,s(row,:),operator_columns,operator_metric,h,z,c(row)
          if(io_status/=0)exit
        enddo
      endif
      if(io_status==0)write(unit,iostat=io_status)fingerprint
      close(unit)
    else
      io_status=0
    endif
    call require(io_status==0,'cannot write legacy checkpoint fixture')
  end subroutine write_legacy_checkpoint

  subroutine legacy_hash_int(fingerprint,value)
    integer(int64),intent(inout)::fingerprint
    integer(int64),intent(in)::value
    fingerprint=ieor(ishftc(fingerprint,9),value)
  end subroutine legacy_hash_int

  subroutine legacy_hash_complex(fingerprint,value)
    integer(int64),intent(inout)::fingerprint
    complex(real64),intent(in)::value
    integer(int64)::bits
    bits=transfer(real(value,real64),bits);call legacy_hash_int(fingerprint,bits)
    bits=transfer(aimag(value),bits);call legacy_hash_int(fingerprint,bits)
  end subroutine legacy_hash_complex

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
      distributed_operators%column_ids(nowned),&
      distributed_operators%hamiltonian_values(nowned),distributed_operators%position_values(3,nowned))
    position=0;edge=0;distributed_metric%row_offsets(1)=1
    do row=n,1,-1
      if(mod(row-1,nproc)/=rank)cycle
      position=position+1;distributed_metric%owned_row_ids(position)=row;owned_coefficients(position)=global_coefficients(row)
      do column=1,n
        edge=edge+1;distributed_metric%column_ids(edge)=column;distributed_metric%values(edge)=s(row,column)
      enddo
      distributed_metric%row_offsets(position+1)=edge+1
      distributed_operators%row_offsets(position)=position
      distributed_operators%column_ids(position)=row
      distributed_operators%hamiltonian_values(position)=h(row,row)
      distributed_operators%position_values(:,position)=z(:,row,row)
    enddo
    distributed_operators%row_offsets(nowned+1)=nowned+1
    distributed_metric%valid=.true.;distributed_metric%global_count=n;distributed_metric%numerical_rank=n
    distributed_metric%max_row_nnz=n;distributed_metric%maximum_value=1.2d0;distributed_metric%condition_estimate=2d0
    distributed_metric%fingerprint=9191_int64;distributed_metric%active_rows=.true.;distributed_metric%packet_ids=[1,1,2,2]
    distributed_operators%valid=.true.;distributed_operators%global_count=n
    distributed_operators%owned_row_ids=distributed_metric%owned_row_ids
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
