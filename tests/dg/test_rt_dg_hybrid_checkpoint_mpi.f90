#include "config.h"
program test_rt_dg_hybrid_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_checkpoint,only:write_rt_dg_hybrid_checkpoint,read_rt_dg_hybrid_checkpoint,&
    s_rt_dg_hybrid_ground_state_payload,write_rt_dg_hybrid_ground_state_checkpoint,&
    read_rt_dg_hybrid_ground_state_checkpoint,read_rt_dg_hybrid_ground_state_checkpoint_coalesced,&
    authenticate_rt_dg_hybrid_ground_state_payload,fingerprint_rt_dg_hybrid_component,&
    collective_rt_dg_hybrid_publication_precondition,&
    rt_dg_hybrid_checkpoint_version,rt_dg_hybrid_occupied_checkpoint_version,&
    rt_dg_hybrid_ground_state_checkpoint_version,rt_dg_hybrid_energy_window_explicit,&
    rt_dg_hybrid_energy_window_legacy_dynamic,rt_dg_hybrid_vector_canonical_momentum
  implicit none
  integer,parameter::n=4,nrt=3,nocc=2,nface_global=2
  integer,parameter::face_point_degree(nface_global)=[4,6],&
    face_weight_degree(nface_global)=[2,3],face_basis_degree(nface_global)=[3,4],&
    face_value_degree(nface_global)=[12,24],face_observable_degree(nface_global)=[4,9]
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
  call exercise_publication_precondition
  call require(rt_dg_hybrid_checkpoint_version==2.and.rt_dg_hybrid_occupied_checkpoint_version==2.and.&
    rt_dg_hybrid_ground_state_checkpoint_version==3.and.rt_dg_hybrid_vector_canonical_momentum==1,&
    'checkpoint family versions or canonical-momentum vector slot changed unexpectedly')
  call get_command_argument(1,mode,length=mode_length);call get_command_argument(1,mode)
  call get_command_argument(2,path)
  if(trim(mode)=='write_complete_legacy_dynamic')then
    call require(nproc==1,'legacy-dynamic full-rank fixture requires one MPI rank')
    call construct_complete_payload(complete_payload)
    call promote_to_legacy_dynamic_full_rank(complete_payload)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),complete_payload,payload_fingerprint,ok,message)
    call require(ok,trim(message))
    if(rank==0)write(*,'(a,i0)')'HYBRID_GS_LEGACY_DYNAMIC phase=write fingerprint=',payload_fingerprint
  else if(trim(mode)=='read_complete_legacy_dynamic')then
    call require(nproc==1,'legacy-dynamic full-rank fixture requires one MPI rank')
    call read_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),restored_payload,payload_fingerprint,ok,message)
    call require(ok,trim(message))
    call verify_legacy_dynamic_full_rank(restored_payload,payload_fingerprint)
    if(rank==0)write(*,'(a,i0)')'HYBRID_GS_LEGACY_DYNAMIC phase=read fingerprint=',payload_fingerprint
  else if(trim(mode)=='write_complete'.or.trim(mode)=='write_complete_between_levels'.or.&
      trim(mode)=='write_complete_signed_spread'.or.trim(mode)=='write_complete_distinct_selection'.or.&
      trim(mode)=='write_bad_complete'.or.trim(mode)=='write_incomplete_complete'.or.&
      trim(mode)=='write_bad_grid_complete'.or.trim(mode)=='write_out_of_range_grid_complete'.or.&
      trim(mode)=='write_missing_position_convention_complete'.or.&
      trim(mode)=='write_bad_explicit_compat_complete'.or.&
      trim(mode)=='write_bad_explicit_proof_complete'.or.trim(mode)=='write_bad_proof_below_requested_complete'.or.&
      trim(mode)=='write_bad_extension_complete'.or.trim(mode)=='write_bad_boundary_complete'.or.&
      trim(mode)=='write_bad_window_mode_complete'.or.trim(mode)=='write_bad_face_id_complete'.or.&
      trim(mode)=='write_bad_nonlocal_id_complete'.or.trim(mode)=='write_duplicate_face_id_complete'.or.&
      trim(mode)=='write_duplicate_nonlocal_id_complete'.or.trim(mode)=='write_bad_electron_defect_complete'.or.&
      trim(mode)=='write_bad_omitted_tail_complete'.or.&
      trim(mode)=='write_bad_face_point_offset_complete'.or.&
      trim(mode)=='write_bad_face_weight_offset_complete'.or.&
      trim(mode)=='write_bad_face_basis_offset_complete'.or.&
      trim(mode)=='write_bad_face_value_offset_complete'.or.&
      trim(mode)=='write_bad_face_observable_offset_complete'.or.&
      trim(mode)=='write_bad_face_point_tail_complete'.or.&
      trim(mode)=='write_bad_face_weight_tail_complete'.or.&
      trim(mode)=='write_bad_face_basis_tail_complete'.or.&
      trim(mode)=='write_bad_face_value_tail_complete'.or.&
      trim(mode)=='write_bad_face_observable_tail_complete'.or.&
      trim(mode)=='write_interrupted_complete')then
    call construct_complete_payload(complete_payload)
    if(trim(mode)=='write_complete_between_levels')then
      complete_payload%energy_window%requested_rank=nrt
      complete_payload%energy_window%extension_states=0
      complete_payload%energy_window%window_size=0.7d0
      complete_payload%energy_window%requested_cutoff=&
        complete_payload%energy_window%e_homo+complete_payload%energy_window%window_size
      complete_payload%energy_window%extension_energy=0d0
    endif
    if(trim(mode)=='write_complete_signed_spread')then
      complete_payload%certified_basis%spreads_after=&
        complete_payload%certified_basis%spreads_before+0.1d0
      complete_payload%certified_basis%spread_after_total=&
        sum(complete_payload%certified_basis%spreads_after)
      complete_payload%certified_basis%spread_improvement=&
        complete_payload%certified_basis%spread_before_total-&
        complete_payload%certified_basis%spread_after_total
    endif
    if(trim(mode)=='write_complete_distinct_selection')then
      complete_payload%requested_ids=[11,12]
      complete_payload%effective_ids=[11,12,13]
      complete_payload%added_ids=[13]
      complete_payload%closure_parent=[11]
      complete_payload%closure_reason=[1]
      complete_payload%closure_action=[2]
    endif
    if(trim(mode)=='write_bad_complete'.and.size(complete_payload%hamiltonian_rows,1)>0)&
      complete_payload%hamiltonian_rows(1,1)=complete_payload%hamiltonian_rows(1,1)+(1d0,0d0)
    if(trim(mode)=='write_incomplete_complete')deallocate(complete_payload%face_values)
    if(trim(mode)=='write_bad_grid_complete'.and.size(complete_payload%grid_ids)>0)&
      complete_payload%grid_ids(1)=2_int64
    if(trim(mode)=='write_out_of_range_grid_complete'.and.size(complete_payload%grid_ids)>0)&
      complete_payload%grid_ids(1)=n+1_int64
    if(trim(mode)=='write_missing_position_convention_complete')then
      complete_payload%position_convention_fingerprint=0_int64
    endif
    if(trim(mode)=='write_bad_explicit_compat_complete')&
      complete_payload%energy_window%compatibility_dynamic_rank=.true.
    if(trim(mode)=='write_bad_explicit_proof_complete')then
      complete_payload%energy_window%proof_state_present=.false.
      complete_payload%energy_window%proof_status=0
    endif
    if(trim(mode)=='write_bad_proof_below_requested_complete')then
      complete_payload%energy_window%requested_rank=nrt
      complete_payload%energy_window%extension_states=0
      complete_payload%energy_window%window_size=1.2d0
      complete_payload%energy_window%requested_cutoff=&
        complete_payload%energy_window%e_homo+complete_payload%energy_window%window_size
      complete_payload%energy_window%extension_energy=0d0
    endif
    if(trim(mode)=='write_bad_extension_complete')then
      complete_payload%energy_window%extension_states=0
      complete_payload%energy_window%extension_energy=0d0
    endif
    if(trim(mode)=='write_bad_boundary_complete')complete_payload%energy_window%boundary_cluster_rank=2
    if(trim(mode)=='write_bad_window_mode_complete')&
      complete_payload%energy_window%mode=complete_payload%energy_window%mode+1
    if(trim(mode)=='write_bad_face_id_complete'.and.size(complete_payload%face_ids)>0)&
      complete_payload%face_ids(1)=0_int64
    if(trim(mode)=='write_bad_nonlocal_id_complete'.and.size(complete_payload%nonlocal_ids)>0)&
      complete_payload%nonlocal_ids(1)=0_int64
    if(trim(mode)=='write_duplicate_face_id_complete'.and.size(complete_payload%face_ids)>0)&
      complete_payload%face_ids(1)=9000000001_int64
    if(trim(mode)=='write_duplicate_nonlocal_id_complete'.and.size(complete_payload%nonlocal_ids)>0)&
      complete_payload%nonlocal_ids(1)=9000000001_int64
    if(trim(mode)=='write_bad_electron_defect_complete')&
      complete_payload%electron_count%defect=0.5d0*complete_payload%electron_count%tolerance
    if(trim(mode)=='write_bad_omitted_tail_complete')&
      complete_payload%electron_count%omitted_tail=2d0*complete_payload%electron_count%tolerance
    if(trim(mode)=='write_bad_face_point_offset_complete')&
      complete_payload%face_offsets(size(complete_payload%face_offsets))=&
        size(complete_payload%face_point_ids)+2
    if(trim(mode)=='write_bad_face_weight_offset_complete')&
      complete_payload%face_weight_offsets(size(complete_payload%face_weight_offsets))=&
        size(complete_payload%face_weights)+2
    if(trim(mode)=='write_bad_face_basis_offset_complete')&
      complete_payload%face_basis_offsets(size(complete_payload%face_basis_offsets))=&
        size(complete_payload%face_basis_ids)+2
    if(trim(mode)=='write_bad_face_value_offset_complete')&
      complete_payload%face_value_offsets(size(complete_payload%face_value_offsets))=&
        size(complete_payload%face_values,2)+2
    if(trim(mode)=='write_bad_face_observable_offset_complete')&
      complete_payload%face_observable_offsets(size(complete_payload%face_observable_offsets))=&
        size(complete_payload%interface_observables,2)+2
    if(trim(mode)=='write_bad_face_point_tail_complete')&
      complete_payload%face_point_ids=[complete_payload%face_point_ids,999999_int64]
    if(trim(mode)=='write_bad_face_weight_tail_complete')&
      complete_payload%face_weights=[complete_payload%face_weights,9d0]
    if(trim(mode)=='write_bad_face_basis_tail_complete')&
      complete_payload%face_basis_ids=[complete_payload%face_basis_ids,n]
    if(trim(mode)=='write_bad_face_value_tail_complete')&
      call append_complex_column(complete_payload%face_values,[(9d0,1d0)])
    if(trim(mode)=='write_bad_face_observable_tail_complete')&
      call append_complex_column(complete_payload%interface_observables,&
        [(9d0,1d0),(8d0,2d0),(7d0,3d0)])
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),complete_payload,payload_fingerprint,ok,message,&
      interrupt_after_write=trim(mode)=='write_interrupted_complete')
    if(trim(mode)=='write_complete'.or.trim(mode)=='write_complete_between_levels'.or.&
        trim(mode)=='write_complete_signed_spread'.or.trim(mode)=='write_complete_distinct_selection')then
      call require(ok,trim(message))
    else
      call require(.not.ok,'inconsistent or incomplete complete payload was published')
    endif
  else if(trim(mode)=='read_complete'.or.trim(mode)=='read_complete_coalesced'.or.&
      trim(mode)=='read_complete_coalesced_corrupt'.or.trim(mode)=='read_complete_auth'.or.&
      trim(mode)=='read_complete_corrupt')then
    if(trim(mode)=='read_complete_coalesced'.or.trim(mode)=='read_complete_coalesced_corrupt')then
      call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,trim(path),restored_payload,&
        payload_fingerprint,ok,message)
    else
      call read_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),restored_payload,payload_fingerprint,ok,message)
    endif
    if(trim(mode)=='read_complete'.or.trim(mode)=='read_complete_coalesced'.or.trim(mode)=='read_complete_auth')then
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
        size(restored_payload%nonlocal_values,2)==size(restored_payload%nonlocal_ids),&
        'complete checkpoint lost basis, face, or nonlocal payload')
      i=merge(1,0,size(restored_payload%metric_column_ids)/=size(restored_payload%operator_column_ids))
      call MPI_Allreduce(MPI_IN_PLACE,i,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      call require(i==1,'complete checkpoint collapsed independent metric and operator graphs')
      call verify_v3_payload(restored_payload,payload_fingerprint)
      if(trim(mode)=='read_complete_auth')call verify_named_tamper_rejection(restored_payload,payload_fingerprint)
    else if(trim(mode)=='read_complete_corrupt'.or.trim(mode)=='read_complete_coalesced_corrupt')then
      call require(.not.ok,'corrupt complete DG ground-state checkpoint was accepted')
      call require(.not.allocated(restored_payload%row_ids).and.&
        .not.allocated(restored_payload%certified_basis%c_cert).and.&
        .not.allocated(restored_payload%certified_basis%u_rt),&
        'rejected complete checkpoint retained partial output storage')
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
  if(rank==0.and.(trim(mode)=='write_complete'.or.trim(mode)=='read_complete'.or.&
      trim(mode)=='read_complete_coalesced'.or.trim(mode)=='read_complete_auth'))then
    write(*,'(a,i0,a,i0)')'HYBRID_GS_CHECKPOINT ranks=',nproc,' fingerprint=',payload_fingerprint
  endif
  if(rank==0.and.(trim(mode)=='read_complete'.or.trim(mode)=='read_complete_coalesced'.or.&
      trim(mode)=='read_complete_auth'))write(*,'(a,i0,a)')'PASS complete hybrid checkpoint on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine exercise_publication_precondition
    logical::collective_ok
    character(256)::collective_message
    call collective_rt_dg_hybrid_publication_precondition(comm,rank/=0,4,2,collective_ok,collective_message)
    call require(.not.collective_ok.and.index(collective_message,'precondition')>0,&
      'one-rank malformed v3 publication precondition was not collectively rejected')
    call collective_rt_dg_hybrid_publication_precondition(comm,.true.,4,2,collective_ok,collective_message)
    call require(collective_ok,'valid v3 publication precondition was rejected')
    call collective_rt_dg_hybrid_publication_precondition(comm,.true.,merge(5,4,nproc>1.and.rank==0),2,&
      collective_ok,collective_message)
    call require(nproc==1.or..not.collective_ok,'rank-disagreeing v3 publication extent was accepted')
  end subroutine exercise_publication_precondition
  subroutine construct_complete_payload(payload)
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::payload
    integer::row,point,local_row,local_point,nrow,npoint,face,nface,local_face,projector,nprojector,local_projector,&
      entry,component,face_point_count,face_weight_count,face_basis_count,face_value_count,face_observable_count,&
      face_point_position,face_weight_position,face_basis_position,face_value_position,face_observable_position
    logical::fingerprint_ok
    payload%valid=.true.;payload%final_refresh_complete=.true.;payload%analysis_complete=.true.
    payload%identity_only=.false.;payload%global_count=n;payload%noccupied=nocc
    payload%global_grid_count=n;payload%position_convention_fingerprint=1199_int64
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
    allocate(payload%position_rows(3,nrow,n),payload%symmetry_representation(n,n,2))
    payload%position_rows=(0d0,0d0);payload%symmetry_representation=(0d0,0d0)
    do i=1,n
      payload%symmetry_representation(i,i,1)=(1d0,0d0)
      payload%symmetry_representation(i,i,2)=merge((1d0,0d0),(-1d0,0d0),mod(i,2)==1)
    enddo
    allocate(payload%row_ids(nrow),payload%metric_rows(nrow,n),payload%kinetic_rows(nrow,n),&
      payload%nonlocal_rows(nrow,n),payload%local_rows(nrow,n),payload%sipg_rows(nrow,n),&
      payload%hamiltonian_rows(nrow,n),payload%coefficients(nrow,nocc))
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
      local_point=local_point+1;payload%grid_ids(local_point)=point;payload%grid_weights(local_point)=0.25d0
      payload%partition_ids(local_point)=1+mod(point,2);payload%density(local_point)=0.5d0+0.01d0*point
      do row=1,n;payload%basis_values(row,local_point)=cmplx(0.01d0*row*point,-0.02d0*row,real64);enddo
    enddo
    nface=count([(mod(face-1,nproc)==rank,face=1,nface_global)])
    face_point_count=0;face_weight_count=0;face_basis_count=0;face_value_count=0;face_observable_count=0
    do face=1,nface_global
      if(mod(face-1,nproc)/=rank)cycle
      face_point_count=face_point_count+face_point_degree(face)
      face_weight_count=face_weight_count+face_weight_degree(face)
      face_basis_count=face_basis_count+face_basis_degree(face)
      face_value_count=face_value_count+face_value_degree(face)
      face_observable_count=face_observable_count+face_observable_degree(face)
    enddo
    allocate(payload%face_ids(nface),payload%face_point_ids(face_point_count),payload%face_metadata(8,nface),&
      payload%face_offsets(nface+1),payload%face_weight_offsets(nface+1),payload%face_basis_offsets(nface+1),&
      payload%face_value_offsets(nface+1),payload%face_observable_offsets(nface+1),&
      payload%face_basis_ids(face_basis_count),payload%face_normals(3,nface),&
      payload%face_weights(face_weight_count),payload%face_values(1,face_value_count),&
      payload%interface_observables(3,face_observable_count))
    payload%face_offsets(1)=1;payload%face_weight_offsets(1)=1;payload%face_basis_offsets(1)=1
    payload%face_value_offsets(1)=1;payload%face_observable_offsets(1)=1
    local_face=0;face_point_position=0;face_weight_position=0;face_basis_position=0
    face_value_position=0;face_observable_position=0
    do face=nface_global,1,-1
      if(mod(face-1,nproc)/=rank)cycle
      local_face=local_face+1;payload%face_ids(local_face)=200+face
      select case(face)
      case(1)
        payload%face_metadata(:,local_face)=[1,2,0,0,0,2,1,2]
        payload%face_normals(:,local_face)=[1d0,0d0,0d0]
      case(2)
        payload%face_metadata(:,local_face)=[2,1,0,0,0,3,2,2]
        payload%face_normals(:,local_face)=[0d0,1d0,0d0]
      end select
      do entry=1,face_point_degree(face)
        face_point_position=face_point_position+1
        payload%face_point_ids(face_point_position)=int(10000*face+entry,int64)
      enddo
      payload%face_offsets(local_face+1)=face_point_position+1
      do entry=1,face_weight_degree(face)
        face_weight_position=face_weight_position+1
        payload%face_weights(face_weight_position)=0.1d0*face+0.01d0*entry
      enddo
      payload%face_weight_offsets(local_face+1)=face_weight_position+1
      do entry=1,face_basis_degree(face)
        face_basis_position=face_basis_position+1
        payload%face_basis_ids(face_basis_position)=1+mod(face+entry-2,n)
      enddo
      payload%face_basis_offsets(local_face+1)=face_basis_position+1
      do entry=1,face_value_degree(face)
        face_value_position=face_value_position+1
        payload%face_values(1,face_value_position)=cmplx(0.001d0*(100*face+entry),-0.0001d0*entry,real64)
      enddo
      payload%face_value_offsets(local_face+1)=face_value_position+1
      do entry=1,face_observable_degree(face)
        face_observable_position=face_observable_position+1
        do component=1,3
          payload%interface_observables(component,face_observable_position)=&
            cmplx(0.2d0*face+0.01d0*entry+0.001d0*component,-0.002d0*component,real64)
        enddo
      enddo
      payload%face_observable_offsets(local_face+1)=face_observable_position+1
    enddo
    nprojector=count([(mod(projector-1,nproc)==rank,projector=1,2)])
    allocate(payload%nonlocal_ids(nprojector),payload%nonlocal_owner(nprojector),&
      payload%nonlocal_values(2,nprojector));local_projector=0
    do projector=1,2
      if(mod(projector-1,nproc)/=rank)cycle
      local_projector=local_projector+1;payload%nonlocal_ids(local_projector)=400+projector
      payload%nonlocal_owner(local_projector)=20+projector
      payload%nonlocal_values(:,local_projector)=[cmplx(0.11d0*projector,0.02d0,real64),&
        cmplx(0.12d0*projector,0.03d0,real64)]
    enddo
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
    call construct_v3_payload(payload)
  end subroutine construct_complete_payload

  subroutine construct_v3_payload(payload)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::payload
    integer::row,local_row,nrtrow,point
    real(real64)::a
    complex(real64)::diagonal(nrt,nrt),hamiltonian(nrt,nrt),full_u_rt(nrt,nrt)

    payload%construction_catalog%valid=.true.
    payload%construction_catalog%global_count=n
    allocate(payload%construction_catalog%ids(n),payload%construction_catalog%generations(n),&
      payload%construction_catalog%ordering(n),payload%construction_catalog%ownership(n))
    payload%construction_catalog%ids=[101_int64,103_int64,107_int64,109_int64]
    payload%construction_catalog%generations=[2,4,6,8]
    payload%construction_catalog%ordering=[3,1,4,2]
    payload%construction_catalog%ownership=[11,11,12,12]
    payload%construction_catalog%ids_fingerprint=2101_int64
    payload%construction_catalog%generation_fingerprint=2102_int64
    payload%construction_catalog%ordering_fingerprint=2103_int64
    payload%construction_catalog%ownership_fingerprint=2104_int64
    payload%construction_catalog%provenance_fingerprint=2105_int64
    payload%construction_catalog%catalog_fingerprint=payload%catalog_fingerprint

    payload%certified_basis%valid=.true.
    payload%certified_basis%localization_converged=.true.
    payload%certified_basis%localization_symmetry_constrained=.false.
    payload%certified_basis%construction_count=n
    payload%certified_basis%certified_count=nrt
    payload%certified_basis%occupied_count=nocc
    payload%certified_basis%localization_iterations=7
    nrtrow=count([(mod(row-1,nproc)==rank,row=1,nrt)])
    allocate(payload%certified_basis%construction_row_ids(size(payload%row_ids)),&
      payload%certified_basis%transformation_row_ids(nrtrow),&
      payload%certified_basis%c_cert(size(payload%row_ids),nrt),&
      payload%certified_basis%b_rt(size(payload%row_ids),nrt),&
      payload%certified_basis%u_rt(nrtrow,nrt),&
      payload%certified_basis%certified_eigenvalues(nrt),payload%certified_basis%occupations(nocc),&
      payload%certified_basis%initial_occupied_amplitudes(nrt,nocc),&
      payload%certified_basis%centers(3,nrt),payload%certified_basis%spreads_before(nrt),&
      payload%certified_basis%spreads_after(nrt))
    a=1d0/sqrt(2d0);full_u_rt=(0d0,0d0)
    full_u_rt(:,1)=[cmplx(a,0d0,real64),cmplx(a,0d0,real64),(0d0,0d0)]
    full_u_rt(:,2)=[cmplx(0d0,a,real64),cmplx(0d0,-a,real64),(0d0,0d0)]
    full_u_rt(3,3)=(1d0,0d0)
    local_row=0
    do row=1,nrt
      if(mod(row-1,nproc)/=rank)cycle
      local_row=local_row+1
      payload%certified_basis%transformation_row_ids(local_row)=row
      payload%certified_basis%u_rt(local_row,:)=full_u_rt(row,:)
    enddo
    payload%certified_basis%c_cert=(0d0,0d0);payload%certified_basis%b_rt=(0d0,0d0)
    do local_row=1,size(payload%row_ids)
      row=int(payload%row_ids(local_row));payload%certified_basis%construction_row_ids(local_row)=row
      if(row<=nrt)payload%certified_basis%c_cert(local_row,row)=(1d0,0d0)
      payload%certified_basis%b_rt(local_row,:)=matmul(payload%certified_basis%c_cert(local_row,:),&
        full_u_rt)
    enddo
    payload%certified_basis%certified_eigenvalues=[-0.6d0,-0.2d0,0.4d0]
    payload%certified_basis%occupations=[2d0,2d0]
    payload%certified_basis%initial_occupied_amplitudes=&
      conjg(transpose(full_u_rt(1:nocc,:)))
    payload%certified_basis%centers=reshape([(0.1d0*row,row=1,3*nrt)],[3,nrt])
    payload%certified_basis%spreads_before=[1.2d0,1.1d0,1d0]
    payload%certified_basis%spreads_after=[0.9d0,0.8d0,0.7d0]
    payload%certified_basis%spread_before_total=3.3d0
    payload%certified_basis%spread_after_total=2.4d0
    payload%certified_basis%spread_improvement=0.9d0
    payload%certified_basis%transform_unitarity_defect=1d-14
    payload%certified_basis%certified_metric_defect=2d-14
    payload%certified_basis%rt_metric_defect=3d-14
    payload%certified_basis%embedding_defect=4d-14
    payload%certified_basis%projector_invariance_defect=5d-14
    payload%certified_basis%target_symmetry_defect_before=6d-14
    payload%certified_basis%target_symmetry_defect_after=6d-14
    payload%certified_basis%energy_symmetry_defect_before=7d-14
    payload%certified_basis%energy_symmetry_defect_after=7d-14
    payload%certified_basis%symmetry_defect_invariance=8d-14
    payload%certified_basis%scalar_covariance_defect=9d-14
    payload%certified_basis%vector_covariance_defect=1d-13
    payload%certified_basis%tensor_covariance_defect=1.1d-13
    payload%certified_basis%c_cert_fingerprint=2201_int64
    payload%certified_basis%u_rt_fingerprint=2202_int64
    payload%certified_basis%b_rt_fingerprint=2203_int64
    payload%certified_basis%initial_state_fingerprint=2204_int64
    payload%certified_basis%transformation_fingerprint=2205_int64
    payload%certified_basis%operator_fingerprint=2206_int64
    payload%certified_basis%fingerprint=2207_int64

    payload%electron_count%valid=.true.
    payload%electron_count%expected_count=4d0;payload%electron_count%actual_count=4d0
    payload%electron_count%tolerance=1d-11;payload%electron_count%defect=0d0
    payload%electron_count%omitted_tail=0d0;payload%electron_count%chemical_potential=-0.1d0
    payload%electron_count%fingerprint=2301_int64

    payload%rt_space%valid=.true.;payload%rt_space%rank=nrt
    payload%rt_space%operation_count=2;payload%rt_space%scalar_count=1
    payload%rt_space%vector_count=1;payload%rt_space%tensor_count=1
    allocate(payload%rt_space%row_ids(nrtrow),payload%rt_space%row_owner_keys(nrt),&
      payload%rt_space%grid_owner_keys(size(payload%grid_ids)),payload%rt_space%metric_rows(nrtrow,nrt),&
      payload%rt_space%kinetic_rows(nrtrow,nrt),payload%rt_space%nonlocal_rows(nrtrow,nrt),&
      payload%rt_space%local_rows(nrtrow,nrt),payload%rt_space%sipg_rows(nrtrow,nrt),&
      payload%rt_space%hamiltonian_rows(nrtrow,nrt),payload%rt_space%representation(nrt,nrt,2),&
      payload%rt_space%cartesian_rotations(3,3,2),payload%rt_space%scalar_operator_rows(nrtrow,nrt,1),&
      payload%rt_space%vector_operator_rows(nrtrow,nrt,3,1),&
      payload%rt_space%tensor_operator_rows(nrtrow,nrt,3,3,1),&
      payload%rt_space%basis_values(nrt,size(payload%grid_ids)),payload%rt_space%density(size(payload%grid_ids)))
    payload%rt_space%row_owner_keys=[31,32,33]
    do point=1,size(payload%grid_ids)
      payload%rt_space%grid_owner_keys(point)=40+mod(int(payload%grid_ids(point)),2)
      payload%rt_space%density(point)=payload%density(point)
      do row=1,nrt
        payload%rt_space%basis_values(row,point)=cmplx(0.03d0*row*real(payload%grid_ids(point)),&
          -0.01d0*row,real64)
      enddo
    enddo
    diagonal=(0d0,0d0);do row=1,nrt;diagonal(row,row)=payload%certified_basis%certified_eigenvalues(row);enddo
    hamiltonian=matmul(conjg(transpose(full_u_rt)),matmul(diagonal,full_u_rt))
    payload%rt_space%representation=(0d0,0d0);payload%rt_space%cartesian_rotations=0d0
    do row=1,nrt
      payload%rt_space%representation(row,row,1)=(1d0,0d0)
      payload%rt_space%representation(row,row,2)=merge((-1d0,0d0),(1d0,0d0),row<3)
    enddo
    do row=1,3
      payload%rt_space%cartesian_rotations(row,row,1)=1d0
      payload%rt_space%cartesian_rotations(row,row,2)=merge(-1d0,1d0,row<3)
    enddo
    payload%rt_space%metric_rows=(0d0,0d0);payload%rt_space%scalar_operator_rows=(0d0,0d0)
    payload%rt_space%vector_operator_rows=(0d0,0d0);payload%rt_space%tensor_operator_rows=(0d0,0d0)
    local_row=0
    do row=1,nrt
      if(mod(row-1,nproc)/=rank)cycle
      local_row=local_row+1;payload%rt_space%row_ids(local_row)=row
      payload%rt_space%metric_rows(local_row,row)=(1d0,0d0)
      payload%rt_space%hamiltonian_rows(local_row,:)=hamiltonian(row,:)
      payload%rt_space%kinetic_rows(local_row,:)=0.5d0*hamiltonian(row,:)
      payload%rt_space%nonlocal_rows(local_row,:)=0.25d0*hamiltonian(row,:)
      payload%rt_space%local_rows(local_row,:)=0.125d0*hamiltonian(row,:)
      payload%rt_space%sipg_rows(local_row,:)=0.125d0*hamiltonian(row,:)
      payload%rt_space%scalar_operator_rows(local_row,:,1)=hamiltonian(row,:)
      payload%rt_space%vector_operator_rows(local_row,:,1,1)=cmplx(0.01d0*row,0d0,real64)
      payload%rt_space%vector_operator_rows(local_row,:,2,1)=cmplx(0.02d0*row,0d0,real64)
      payload%rt_space%vector_operator_rows(local_row,:,3,1)=cmplx(0.03d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(local_row,:,1,1,1)=cmplx(0.04d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(local_row,:,2,2,1)=cmplx(0.05d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(local_row,:,3,3,1)=cmplx(0.06d0*row,0d0,real64)
    enddo
    payload%rt_space%metric_fingerprint=2401_int64;payload%rt_space%kinetic_fingerprint=2402_int64
    payload%rt_space%nonlocal_fingerprint=2403_int64;payload%rt_space%local_fingerprint=2404_int64
    payload%rt_space%sipg_fingerprint=2405_int64;payload%rt_space%hamiltonian_fingerprint=2406_int64
    payload%rt_space%basis_fingerprint=2407_int64;payload%rt_space%density_fingerprint=2408_int64
    payload%rt_space%ownership_fingerprint=2409_int64;payload%rt_space%scalar_fingerprint=2410_int64
    payload%rt_space%vector_fingerprint=2411_int64;payload%rt_space%tensor_fingerprint=2412_int64
    payload%rt_space%representation_fingerprint=2413_int64;payload%rt_space%fingerprint=2414_int64

    payload%energy_window%valid=.true.;payload%energy_window%compatibility_dynamic_rank=.false.
    payload%energy_window%proof_state_present=.true.;payload%energy_window%mode=rt_dg_hybrid_energy_window_explicit
    payload%energy_window%construction_rank=n;payload%energy_window%solved_rank=n
    payload%energy_window%occupied_rank=nocc;payload%energy_window%requested_rank=2
    payload%energy_window%certified_rank=nrt;payload%energy_window%extension_states=1
    payload%energy_window%boundary_cluster_rank=nrt;payload%energy_window%proof_status=1
    payload%energy_window%window_size=0.35d0;payload%energy_window%e_homo=-0.2d0
    payload%energy_window%requested_cutoff=payload%energy_window%e_homo+payload%energy_window%window_size
    payload%energy_window%certified_cutoff=payload%certified_basis%certified_eigenvalues(nrt)
    payload%energy_window%extension_energy=max(0d0,payload%energy_window%certified_cutoff-&
      payload%energy_window%requested_cutoff)
    payload%energy_window%proof_energy=0.9d0
    payload%energy_window%fingerprint=2501_int64

    payload%symmetry_receipt%valid=.true.;payload%symmetry_receipt%worst_operation=2
    payload%symmetry_receipt%occupied_subspace_defect=1d-13
    payload%symmetry_receipt%occupied_projector_defect=2d-13
    payload%symmetry_receipt%target_subspace_defect=3d-13
    payload%symmetry_receipt%target_energy_defect=4d-13
    payload%symmetry_receipt%density_defect=5d-13
    payload%symmetry_receipt%scalar_covariance_defect=6d-13
    payload%symmetry_receipt%vector_covariance_defect=7d-13
    payload%symmetry_receipt%tensor_covariance_defect=8d-13
    payload%symmetry_receipt%final_basis_defect=9d-13
    payload%symmetry_receipt%worst_operation_defect=1d-12
    payload%symmetry_receipt%maximum_physical_defect=1d-12
    payload%symmetry_receipt%fingerprint=2601_int64

    payload%handoff_receipts%valid=.true.
    payload%handoff_receipts%position_fingerprint=payload%position_convention_fingerprint
    payload%handoff_receipts%nonlocal_fingerprint=payload%nonlocal_fingerprint
    payload%handoff_receipts%face_fingerprint=payload%face_fingerprint
    payload%handoff_receipts%pseudopotential_fingerprint=payload%pseudopotential_fingerprint
    payload%handoff_receipts%transformation_fingerprint=payload%certified_basis%transformation_fingerprint
    payload%handoff_receipts%fingerprint=2701_int64
  end subroutine construct_v3_payload

  subroutine promote_to_legacy_dynamic_full_rank(payload)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::payload
    integer::row,point,npoint
    real(real64)::eigenvalue

    npoint=size(payload%grid_ids)
    payload%certified_basis%certified_count=n
    deallocate(payload%certified_basis%construction_row_ids,&
      payload%certified_basis%transformation_row_ids,payload%certified_basis%c_cert,&
      payload%certified_basis%u_rt,payload%certified_basis%b_rt,&
      payload%certified_basis%initial_occupied_amplitudes,&
      payload%certified_basis%certified_eigenvalues,payload%certified_basis%occupations,&
      payload%certified_basis%centers,payload%certified_basis%spreads_before,&
      payload%certified_basis%spreads_after)
    allocate(payload%certified_basis%construction_row_ids(n),&
      payload%certified_basis%transformation_row_ids(n),&
      payload%certified_basis%c_cert(n,n),payload%certified_basis%u_rt(n,n),&
      payload%certified_basis%b_rt(n,n),&
      payload%certified_basis%initial_occupied_amplitudes(n,nocc),&
      payload%certified_basis%certified_eigenvalues(n),&
      payload%certified_basis%occupations(nocc),payload%certified_basis%centers(3,n),&
      payload%certified_basis%spreads_before(n),payload%certified_basis%spreads_after(n))
    payload%certified_basis%construction_row_ids=[(int(row,int64),row=1,n)]
    payload%certified_basis%transformation_row_ids=[(int(row,int64),row=1,n)]
    payload%certified_basis%c_cert=(0d0,0d0);payload%certified_basis%u_rt=(0d0,0d0)
    payload%certified_basis%b_rt=(0d0,0d0)
    do row=1,n
      payload%certified_basis%c_cert(row,row)=(1d0,0d0)
      payload%certified_basis%u_rt(row,row)=(1d0,0d0)
      payload%certified_basis%b_rt(row,row)=(1d0,0d0)
    enddo
    payload%certified_basis%initial_occupied_amplitudes=(0d0,0d0)
    do row=1,nocc
      payload%certified_basis%initial_occupied_amplitudes(row,row)=(1d0,0d0)
    enddo
    payload%certified_basis%certified_eigenvalues=[-0.6d0,-0.2d0,0.4d0,0.9d0]
    payload%certified_basis%occupations=[2d0,2d0]
    payload%certified_basis%centers=reshape([(0.1d0*row,row=1,3*n)],[3,n])
    payload%certified_basis%spreads_before=[1.2d0,1.1d0,1d0,0.9d0]
    payload%certified_basis%spreads_after=[0.9d0,0.8d0,0.7d0,0.6d0]
    payload%certified_basis%spread_before_total=4.2d0
    payload%certified_basis%spread_after_total=3d0
    payload%certified_basis%spread_improvement=1.2d0
    payload%certified_basis%c_cert_fingerprint=3201_int64
    payload%certified_basis%u_rt_fingerprint=3202_int64
    payload%certified_basis%b_rt_fingerprint=3203_int64
    payload%certified_basis%initial_state_fingerprint=3204_int64
    payload%certified_basis%transformation_fingerprint=3205_int64
    payload%certified_basis%operator_fingerprint=3206_int64
    payload%certified_basis%fingerprint=3207_int64

    payload%rt_space%rank=n
    deallocate(payload%rt_space%row_ids,payload%rt_space%row_owner_keys,&
      payload%rt_space%metric_rows,payload%rt_space%kinetic_rows,&
      payload%rt_space%nonlocal_rows,payload%rt_space%local_rows,&
      payload%rt_space%sipg_rows,payload%rt_space%hamiltonian_rows,&
      payload%rt_space%representation,payload%rt_space%scalar_operator_rows,&
      payload%rt_space%vector_operator_rows,payload%rt_space%tensor_operator_rows,&
      payload%rt_space%basis_values,payload%rt_space%density)
    allocate(payload%rt_space%row_ids(n),payload%rt_space%row_owner_keys(n),&
      payload%rt_space%metric_rows(n,n),payload%rt_space%kinetic_rows(n,n),&
      payload%rt_space%nonlocal_rows(n,n),payload%rt_space%local_rows(n,n),&
      payload%rt_space%sipg_rows(n,n),payload%rt_space%hamiltonian_rows(n,n),&
      payload%rt_space%representation(n,n,payload%operation_count),&
      payload%rt_space%scalar_operator_rows(n,n,payload%rt_space%scalar_count),&
      payload%rt_space%vector_operator_rows(n,n,3,payload%rt_space%vector_count),&
      payload%rt_space%tensor_operator_rows(n,n,3,3,payload%rt_space%tensor_count),&
      payload%rt_space%basis_values(n,npoint),payload%rt_space%density(npoint))
    payload%rt_space%row_ids=[(int(row,int64),row=1,n)]
    payload%rt_space%row_owner_keys=[31,32,33,34]
    payload%rt_space%metric_rows=(0d0,0d0);payload%rt_space%kinetic_rows=(0d0,0d0)
    payload%rt_space%nonlocal_rows=(0d0,0d0);payload%rt_space%local_rows=(0d0,0d0)
    payload%rt_space%sipg_rows=(0d0,0d0);payload%rt_space%hamiltonian_rows=(0d0,0d0)
    payload%rt_space%representation=(0d0,0d0)
    payload%rt_space%scalar_operator_rows=(0d0,0d0)
    payload%rt_space%vector_operator_rows=(0d0,0d0)
    payload%rt_space%tensor_operator_rows=(0d0,0d0)
    do row=1,n
      eigenvalue=payload%certified_basis%certified_eigenvalues(row)
      payload%rt_space%metric_rows(row,row)=(1d0,0d0)
      payload%rt_space%kinetic_rows(row,row)=cmplx(0.5d0*eigenvalue,0d0,real64)
      payload%rt_space%nonlocal_rows(row,row)=cmplx(0.25d0*eigenvalue,0d0,real64)
      payload%rt_space%local_rows(row,row)=cmplx(0.125d0*eigenvalue,0d0,real64)
      payload%rt_space%sipg_rows(row,row)=cmplx(0.125d0*eigenvalue,0d0,real64)
      payload%rt_space%hamiltonian_rows(row,row)=cmplx(eigenvalue,0d0,real64)
      payload%rt_space%representation(row,row,1)=(1d0,0d0)
      payload%rt_space%representation(row,row,2)=merge((-1d0,0d0),(1d0,0d0),row<3)
      payload%rt_space%scalar_operator_rows(row,:,1)=payload%rt_space%hamiltonian_rows(row,:)
      payload%rt_space%vector_operator_rows(row,:,1,1)=cmplx(0.01d0*row,0d0,real64)
      payload%rt_space%vector_operator_rows(row,:,2,1)=cmplx(0.02d0*row,0d0,real64)
      payload%rt_space%vector_operator_rows(row,:,3,1)=cmplx(0.03d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(row,:,1,1,1)=cmplx(0.04d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(row,:,2,2,1)=cmplx(0.05d0*row,0d0,real64)
      payload%rt_space%tensor_operator_rows(row,:,3,3,1)=cmplx(0.06d0*row,0d0,real64)
    enddo
    do point=1,npoint
      payload%rt_space%density(point)=payload%density(point)
      do row=1,n
        payload%rt_space%basis_values(row,point)=cmplx(&
          0.03d0*row*real(payload%grid_ids(point),real64),-0.01d0*row,real64)
      enddo
    enddo
    payload%rt_space%metric_fingerprint=3401_int64
    payload%rt_space%kinetic_fingerprint=3402_int64
    payload%rt_space%nonlocal_fingerprint=3403_int64
    payload%rt_space%local_fingerprint=3404_int64
    payload%rt_space%sipg_fingerprint=3405_int64
    payload%rt_space%hamiltonian_fingerprint=3406_int64
    payload%rt_space%basis_fingerprint=3407_int64
    payload%rt_space%density_fingerprint=3408_int64
    payload%rt_space%ownership_fingerprint=3409_int64
    payload%rt_space%scalar_fingerprint=3410_int64
    payload%rt_space%vector_fingerprint=3411_int64
    payload%rt_space%tensor_fingerprint=3412_int64
    payload%rt_space%representation_fingerprint=3413_int64
    payload%rt_space%fingerprint=3414_int64

    payload%energy_window%mode=rt_dg_hybrid_energy_window_legacy_dynamic
    payload%energy_window%compatibility_dynamic_rank=.true.
    payload%energy_window%window_size=-1d0
    payload%energy_window%construction_rank=n
    payload%energy_window%solved_rank=n
    payload%energy_window%occupied_rank=nocc
    payload%energy_window%requested_rank=nocc
    payload%energy_window%certified_rank=n
    payload%energy_window%extension_states=n-nocc
    payload%energy_window%boundary_cluster_rank=n
    payload%energy_window%e_homo=payload%certified_basis%certified_eigenvalues(nocc)
    payload%energy_window%requested_cutoff=&
      payload%certified_basis%certified_eigenvalues(payload%energy_window%requested_rank)
    payload%energy_window%certified_cutoff=payload%certified_basis%certified_eigenvalues(n)
    payload%energy_window%extension_energy=max(0d0,payload%energy_window%certified_cutoff-&
      payload%energy_window%requested_cutoff)
    payload%energy_window%proof_state_present=.false.
    payload%energy_window%proof_status=0
    payload%energy_window%proof_energy=0d0
    payload%energy_window%fingerprint=3501_int64
    payload%handoff_receipts%transformation_fingerprint=&
      payload%certified_basis%transformation_fingerprint
  end subroutine promote_to_legacy_dynamic_full_rank

  subroutine verify_legacy_dynamic_full_rank(payload,fingerprint)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::fingerprint
    logical::authenticated
    character(256)::authentication_message
    call require(payload%certified_basis%certified_count==n.and.payload%rt_space%rank==n,&
      'legacy-dynamic v3 checkpoint did not retain the full certified rank')
    call require(size(payload%certified_basis%c_cert,2)==n.and.&
      size(payload%certified_basis%u_rt,1)==n.and.size(payload%certified_basis%u_rt,2)==n.and.&
      size(payload%certified_basis%b_rt,2)==n.and.size(payload%rt_space%metric_rows,2)==n,&
      'legacy-dynamic v3 checkpoint changed full-rank named payload dimensions')
    call require(all(payload%certified_basis%construction_row_ids==[(int(i,int64),i=1,n)]).and.&
      all(payload%certified_basis%transformation_row_ids==[(int(i,int64),i=1,n)]).and.&
      all(payload%rt_space%row_ids==[(int(i,int64),i=1,n)]),&
      'legacy-dynamic v3 checkpoint changed row-owned full-rank catalogs')
    call require(payload%energy_window%mode==rt_dg_hybrid_energy_window_legacy_dynamic.and.&
      payload%energy_window%compatibility_dynamic_rank.and.payload%energy_window%window_size==-1d0.and.&
      payload%energy_window%construction_rank==n.and.payload%energy_window%solved_rank==n.and.&
      payload%energy_window%certified_rank==n.and.payload%energy_window%boundary_cluster_rank==n,&
      'legacy-dynamic v3 checkpoint changed the full-rank energy-window receipt')
    call require(.not.payload%energy_window%proof_state_present.and.&
      payload%energy_window%proof_status==0.and.payload%energy_window%proof_energy==0d0,&
      'legacy-dynamic full-rank receipt invented a proof state')
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,payload,fingerprint,authenticated,&
      authentication_message)
    call require(authenticated,trim(authentication_message))
  end subroutine verify_legacy_dynamic_full_rank

  subroutine verify_v3_payload(payload,fingerprint)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::fingerprint
    integer::local_counts(n),global_counts(n),local_rt(nrt),global_rt(nrt),&
      local_u(nrt),global_u(nrt),p,row
    complex(real64)::local_u_rt(nrt,nrt),full_u_rt(nrt,nrt)
    logical::authenticated,local_ok
    character(256)::authentication_message
    local_counts=0;local_rt=0;local_u=0;local_u_rt=(0d0,0d0);local_ok=.true.
    do row=1,size(payload%certified_basis%transformation_row_ids)
      local_u(int(payload%certified_basis%transformation_row_ids(row)))=&
        local_u(int(payload%certified_basis%transformation_row_ids(row)))+1
      local_u_rt(int(payload%certified_basis%transformation_row_ids(row)),:)=payload%certified_basis%u_rt(row,:)
      local_ok=local_ok.and.mod(int(payload%certified_basis%transformation_row_ids(row))-1,nproc)==rank
    enddo
    call MPI_Allreduce(local_u,global_u,nrt,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_u_rt,full_u_rt,nrt*nrt,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do row=1,size(payload%certified_basis%construction_row_ids)
      local_counts(int(payload%certified_basis%construction_row_ids(row)))=&
      local_counts(int(payload%certified_basis%construction_row_ids(row)))+1
      local_ok=local_ok.and.maxval(abs(payload%certified_basis%b_rt(row,:)-&
        matmul(payload%certified_basis%c_cert(row,:),full_u_rt)))<1d-13
      local_ok=local_ok.and.mod(int(payload%certified_basis%construction_row_ids(row))-1,nproc)==rank
    enddo
    do row=1,size(payload%grid_ids)
      local_ok=local_ok.and.mod(int(payload%grid_ids(row))-1,nproc)==rank
    enddo
    do row=1,size(payload%face_ids)
      local_ok=local_ok.and.mod(int(payload%face_ids(row))-1,nproc)==rank
    enddo
    do row=1,size(payload%nonlocal_ids)
      local_ok=local_ok.and.mod(int(payload%nonlocal_ids(row))-1,nproc)==rank
    enddo
    do row=1,size(payload%rt_space%row_ids)
      local_rt(int(payload%rt_space%row_ids(row)))=local_rt(int(payload%rt_space%row_ids(row)))+1
      local_ok=local_ok.and.mod(int(payload%rt_space%row_ids(row))-1,nproc)==rank
    enddo
    call MPI_Allreduce(local_counts,global_counts,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_rt,global_rt,nrt,MPI_INTEGER,MPI_SUM,comm,ierr)
    call require(local_ok,'v3 checkpoint changed the embedding or canonical row layout')
    call require(all(global_counts==1).and.all(global_rt==1).and.all(global_u==1),&
      'v3 checkpoint lost construction, transformation, or RT row ownership')
    call require(payload%construction_catalog%valid.and.payload%construction_catalog%global_count==n.and.&
      payload%certified_basis%valid.and.payload%certified_basis%construction_count==n.and.&
      payload%certified_basis%certified_count==nrt.and.payload%certified_basis%occupied_count==nocc,&
      'v3 checkpoint changed construction/certified dimensions')
    call require(payload%rt_space%valid.and.payload%rt_space%rank==nrt.and.&
      payload%energy_window%valid.and.payload%energy_window%certified_rank==nrt.and.&
      payload%energy_window%proof_state_present.and.&
      payload%energy_window%proof_energy>payload%energy_window%certified_cutoff,&
      'v3 checkpoint lost certified-window proof')
    call require(maxval(abs(payload%certified_basis%initial_occupied_amplitudes-&
      conjg(transpose(full_u_rt(1:nocc,:)))))<1d-13,&
      'v3 checkpoint changed initial occupied amplitudes')
    call require(abs(sum(payload%certified_basis%occupations)-payload%electron_count%actual_count)<1d-13.and.&
      payload%electron_count%defect<=payload%electron_count%tolerance,&
      'v3 checkpoint changed the electron-count receipt')
    call require(size(payload%rt_space%scalar_operator_rows,3)==1.and.&
      size(payload%rt_space%vector_operator_rows,3)==3.and.size(payload%rt_space%vector_operator_rows,4)==1.and.&
      size(payload%rt_space%tensor_operator_rows,3)==3.and.size(payload%rt_space%tensor_operator_rows,4)==3.and.&
      size(payload%rt_space%tensor_operator_rows,5)==1,&
      'v3 checkpoint changed scalar/vector/tensor operator dimensions')
    call require(all(payload%rt_space%hamiltonian_rows==payload%rt_space%kinetic_rows+&
      payload%rt_space%nonlocal_rows+payload%rt_space%local_rows+payload%rt_space%sipg_rows),&
      'v3 checkpoint changed final RT Hamiltonian components')
    call require(payload%handoff_receipts%position_fingerprint==payload%position_convention_fingerprint.and.&
      payload%handoff_receipts%nonlocal_fingerprint==payload%nonlocal_fingerprint.and.&
      payload%handoff_receipts%face_fingerprint==payload%face_fingerprint.and.&
      payload%handoff_receipts%pseudopotential_fingerprint==payload%pseudopotential_fingerprint,&
      'v3 checkpoint lost named handoff provenance')
    call verify_face_payload(payload)
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,payload,fingerprint,authenticated,&
      authentication_message)
    call require(authenticated,trim(authentication_message))
    p=0
  end subroutine verify_v3_payload

  subroutine verify_face_payload(payload)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer::local_counts(nface_global),global_counts(nface_global),local_face,face,entry,component,position
    integer::expected_metadata(8)
    real(real64)::expected_normal(3),expected_weight
    complex(real64)::expected_value
    logical::local_ok
    local_counts=0;local_ok=allocated(payload%face_ids).and.allocated(payload%face_point_ids).and.&
      allocated(payload%face_metadata).and.allocated(payload%face_offsets).and.&
      allocated(payload%face_weight_offsets).and.allocated(payload%face_basis_offsets).and.&
      allocated(payload%face_value_offsets).and.allocated(payload%face_observable_offsets).and.&
      allocated(payload%face_basis_ids).and.allocated(payload%face_normals).and.&
      allocated(payload%face_weights).and.allocated(payload%face_values).and.&
      allocated(payload%interface_observables)
    if(local_ok)then
      local_ok=size(payload%face_metadata,1)==8.and.size(payload%face_metadata,2)==size(payload%face_ids).and.&
        size(payload%face_normals,1)==3.and.size(payload%face_normals,2)==size(payload%face_ids).and.&
        size(payload%face_offsets)==size(payload%face_ids)+1.and.&
        size(payload%face_weight_offsets)==size(payload%face_ids)+1.and.&
        size(payload%face_basis_offsets)==size(payload%face_ids)+1.and.&
        size(payload%face_value_offsets)==size(payload%face_ids)+1.and.&
        size(payload%face_observable_offsets)==size(payload%face_ids)+1.and.&
        size(payload%face_values,1)==1.and.size(payload%interface_observables,1)==3
    endif
    if(local_ok)then
      local_ok=payload%face_offsets(1)==1.and.payload%face_weight_offsets(1)==1.and.&
        payload%face_basis_offsets(1)==1.and.payload%face_value_offsets(1)==1.and.&
        payload%face_observable_offsets(1)==1.and.&
        payload%face_offsets(size(payload%face_offsets))==size(payload%face_point_ids)+1.and.&
        payload%face_weight_offsets(size(payload%face_weight_offsets))==size(payload%face_weights)+1.and.&
        payload%face_basis_offsets(size(payload%face_basis_offsets))==size(payload%face_basis_ids)+1.and.&
        payload%face_value_offsets(size(payload%face_value_offsets))==size(payload%face_values,2)+1.and.&
        payload%face_observable_offsets(size(payload%face_observable_offsets))==&
          size(payload%interface_observables,2)+1
    endif
    if(local_ok)then
      do local_face=1,size(payload%face_ids)
        face=int(payload%face_ids(local_face)-200_int64)
        if(face<1.or.face>nface_global)then;local_ok=.false.;cycle;endif
        local_counts(face)=local_counts(face)+1
        local_ok=local_ok.and.payload%face_offsets(local_face+1)-payload%face_offsets(local_face)==&
          face_point_degree(face)
        local_ok=local_ok.and.payload%face_weight_offsets(local_face+1)-payload%face_weight_offsets(local_face)==&
          face_weight_degree(face)
        local_ok=local_ok.and.payload%face_basis_offsets(local_face+1)-payload%face_basis_offsets(local_face)==&
          face_basis_degree(face)
        local_ok=local_ok.and.payload%face_value_offsets(local_face+1)-payload%face_value_offsets(local_face)==&
          face_value_degree(face)
        local_ok=local_ok.and.payload%face_observable_offsets(local_face+1)-&
          payload%face_observable_offsets(local_face)==face_observable_degree(face)
        if(face==1)then
          expected_metadata=[1,2,0,0,0,2,1,2];expected_normal=[1d0,0d0,0d0]
        else
          expected_metadata=[2,1,0,0,0,3,2,2];expected_normal=[0d0,1d0,0d0]
        endif
        local_ok=local_ok.and.all(payload%face_metadata(:,local_face)==expected_metadata).and.&
          maxval(abs(payload%face_normals(:,local_face)-expected_normal))<1d-15
        do entry=1,face_point_degree(face)
          position=payload%face_offsets(local_face)+entry-1
          local_ok=local_ok.and.payload%face_point_ids(position)==int(10000*face+entry,int64)
        enddo
        do entry=1,face_weight_degree(face)
          position=payload%face_weight_offsets(local_face)+entry-1
          expected_weight=0.1d0*face+0.01d0*entry
          local_ok=local_ok.and.abs(payload%face_weights(position)-expected_weight)<1d-15
        enddo
        do entry=1,face_basis_degree(face)
          position=payload%face_basis_offsets(local_face)+entry-1
          local_ok=local_ok.and.payload%face_basis_ids(position)==1+mod(face+entry-2,n)
        enddo
        do entry=1,face_value_degree(face)
          position=payload%face_value_offsets(local_face)+entry-1
          expected_value=cmplx(0.001d0*(100*face+entry),-0.0001d0*entry,real64)
          local_ok=local_ok.and.abs(payload%face_values(1,position)-expected_value)<1d-15
        enddo
        do entry=1,face_observable_degree(face)
          position=payload%face_observable_offsets(local_face)+entry-1
          do component=1,3
            expected_value=cmplx(0.2d0*face+0.01d0*entry+0.001d0*component,&
              -0.002d0*component,real64)
            local_ok=local_ok.and.abs(payload%interface_observables(component,position)-expected_value)<1d-15
          enddo
        enddo
      enddo
    endif
    call MPI_Allreduce(local_counts,global_counts,nface_global,MPI_INTEGER,MPI_SUM,comm,ierr)
    call require(local_ok.and.all(global_counts==1),&
      'v3 checkpoint changed variable-degree face offsets or payload content')
  end subroutine verify_face_payload

  subroutine verify_named_tamper_rejection(payload,fingerprint)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::fingerprint
    type(s_rt_dg_hybrid_ground_state_payload)::tampered
    integer::field,p,row
    logical::authenticated
    character(256)::authentication_message,label
    complex(real64)::local_u_rt(nrt,nrt),full_u_rt(nrt,nrt),swap_row(nrt)
    integer,parameter::named_field_count=138
    local_u_rt=(0d0,0d0)
    do row=1,size(payload%certified_basis%transformation_row_ids)
      local_u_rt(int(payload%certified_basis%transformation_row_ids(row)),:)=payload%certified_basis%u_rt(row,:)
    enddo
    call MPI_Allreduce(local_u_rt,full_u_rt,nrt*nrt,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do field=1,named_field_count
      tampered=payload
      select case(field)
      case(1);tampered%construction_catalog%valid=.false.
      case(2);tampered%construction_catalog%global_count=n+1
      case(3);tampered%construction_catalog%ids(1)=tampered%construction_catalog%ids(1)+1
      case(4);tampered%construction_catalog%generations(1)=tampered%construction_catalog%generations(1)+1
      case(5);tampered%construction_catalog%ordering(1)=tampered%construction_catalog%ordering(1)+1
      case(6);tampered%construction_catalog%ownership(1)=tampered%construction_catalog%ownership(1)+1
      case(7);tampered%construction_catalog%ids_fingerprint=tampered%construction_catalog%ids_fingerprint+1
      case(8);tampered%construction_catalog%generation_fingerprint=tampered%construction_catalog%generation_fingerprint+1
      case(9);tampered%construction_catalog%ordering_fingerprint=tampered%construction_catalog%ordering_fingerprint+1
      case(10);tampered%construction_catalog%ownership_fingerprint=tampered%construction_catalog%ownership_fingerprint+1
      case(11);tampered%construction_catalog%provenance_fingerprint=tampered%construction_catalog%provenance_fingerprint+1
      case(12);tampered%construction_catalog%catalog_fingerprint=tampered%construction_catalog%catalog_fingerprint+1
      case(13);tampered%certified_basis%valid=.false.
      case(14);tampered%certified_basis%localization_converged=.false.
      case(15);tampered%certified_basis%localization_symmetry_constrained=.true.
      case(16);tampered%certified_basis%construction_count=n+1
      case(17);tampered%certified_basis%certified_count=nrt+1
      case(18);tampered%certified_basis%occupied_count=nocc+1
      case(19);tampered%certified_basis%localization_iterations=tampered%certified_basis%localization_iterations+1
      case(20)
        p=owned_position(tampered%certified_basis%construction_row_ids,1_int64)
        if(p>0)tampered%certified_basis%construction_row_ids(p)=2_int64
      case(21)
        p=owned_position(tampered%certified_basis%construction_row_ids,1_int64)
        if(p>0)tampered%certified_basis%c_cert(p,1)=tampered%certified_basis%c_cert(p,1)+(0.01d0,0d0)
      case(22)
        p=owned_position(tampered%certified_basis%transformation_row_ids,1_int64)
        if(p>0)tampered%certified_basis%u_rt(p,1)=tampered%certified_basis%u_rt(p,1)+(0.01d0,0d0)
      case(23)
        p=owned_position(tampered%certified_basis%construction_row_ids,1_int64)
        if(p>0)tampered%certified_basis%b_rt(p,1)=tampered%certified_basis%b_rt(p,1)+(0.01d0,0d0)
      case(24);tampered%certified_basis%certified_eigenvalues(1)=tampered%certified_basis%certified_eigenvalues(1)+0.01d0
      case(25);tampered%certified_basis%occupations(1)=tampered%certified_basis%occupations(1)-0.01d0
      case(26);tampered%certified_basis%initial_occupied_amplitudes(1,1)=&
        tampered%certified_basis%initial_occupied_amplitudes(1,1)+(0.01d0,0d0)
      case(27);tampered%certified_basis%centers(1,1)=tampered%certified_basis%centers(1,1)+0.01d0
      case(28);tampered%certified_basis%spreads_before(1)=tampered%certified_basis%spreads_before(1)+0.01d0
      case(29);tampered%certified_basis%spreads_after(1)=tampered%certified_basis%spreads_after(1)+0.01d0
      case(30);tampered%certified_basis%spread_before_total=tampered%certified_basis%spread_before_total+0.01d0
      case(31);tampered%certified_basis%spread_after_total=tampered%certified_basis%spread_after_total+0.01d0
      case(32);tampered%certified_basis%spread_improvement=tampered%certified_basis%spread_improvement+0.01d0
      case(33);tampered%certified_basis%transform_unitarity_defect=tampered%certified_basis%transform_unitarity_defect+1d-14
      case(34);tampered%certified_basis%certified_metric_defect=tampered%certified_basis%certified_metric_defect+1d-14
      case(35);tampered%certified_basis%rt_metric_defect=tampered%certified_basis%rt_metric_defect+1d-14
      case(36);tampered%certified_basis%embedding_defect=tampered%certified_basis%embedding_defect+1d-14
      case(37);tampered%certified_basis%projector_invariance_defect=&
        tampered%certified_basis%projector_invariance_defect+1d-14
      case(38);tampered%certified_basis%target_symmetry_defect_before=&
        tampered%certified_basis%target_symmetry_defect_before+1d-14
      case(39);tampered%certified_basis%target_symmetry_defect_after=&
        tampered%certified_basis%target_symmetry_defect_after+1d-14
      case(40);tampered%certified_basis%energy_symmetry_defect_before=&
        tampered%certified_basis%energy_symmetry_defect_before+1d-14
      case(41);tampered%certified_basis%energy_symmetry_defect_after=&
        tampered%certified_basis%energy_symmetry_defect_after+1d-14
      case(42);tampered%certified_basis%symmetry_defect_invariance=&
        tampered%certified_basis%symmetry_defect_invariance+1d-14
      case(43);tampered%certified_basis%scalar_covariance_defect=&
        tampered%certified_basis%scalar_covariance_defect+1d-14
      case(44);tampered%certified_basis%vector_covariance_defect=&
        tampered%certified_basis%vector_covariance_defect+1d-14
      case(45);tampered%certified_basis%tensor_covariance_defect=&
        tampered%certified_basis%tensor_covariance_defect+1d-14
      case(46);tampered%certified_basis%c_cert_fingerprint=tampered%certified_basis%c_cert_fingerprint+1
      case(47);tampered%certified_basis%u_rt_fingerprint=tampered%certified_basis%u_rt_fingerprint+1
      case(48);tampered%certified_basis%b_rt_fingerprint=tampered%certified_basis%b_rt_fingerprint+1
      case(49);tampered%certified_basis%initial_state_fingerprint=tampered%certified_basis%initial_state_fingerprint+1
      case(50);tampered%certified_basis%transformation_fingerprint=&
        tampered%certified_basis%transformation_fingerprint+1
      case(51);tampered%certified_basis%operator_fingerprint=tampered%certified_basis%operator_fingerprint+1
      case(52);tampered%certified_basis%fingerprint=tampered%certified_basis%fingerprint+1
      case(53);tampered%electron_count%valid=.false.
      case(54);tampered%electron_count%expected_count=tampered%electron_count%expected_count+0.01d0
      case(55);tampered%electron_count%actual_count=tampered%electron_count%actual_count+0.01d0
      case(56);tampered%electron_count%tolerance=tampered%electron_count%tolerance*2d0
      case(57);tampered%electron_count%defect=tampered%electron_count%defect+1d-13
      case(58);tampered%electron_count%omitted_tail=tampered%electron_count%omitted_tail+1d-13
      case(59);tampered%electron_count%chemical_potential=tampered%electron_count%chemical_potential+0.01d0
      case(60);tampered%electron_count%fingerprint=tampered%electron_count%fingerprint+1
      case(61);tampered%rt_space%valid=.false.
      case(62);tampered%rt_space%rank=nrt+1
      case(63);tampered%rt_space%operation_count=tampered%rt_space%operation_count+1
      case(64);tampered%rt_space%scalar_count=tampered%rt_space%scalar_count+1
      case(65);tampered%rt_space%vector_count=tampered%rt_space%vector_count+1
      case(66);tampered%rt_space%tensor_count=tampered%rt_space%tensor_count+1
      case(67)
        p=owned_position(tampered%rt_space%row_ids,1_int64);if(p>0)tampered%rt_space%row_ids(p)=2_int64
      case(68);tampered%rt_space%row_owner_keys(1)=tampered%rt_space%row_owner_keys(1)+1
      case(69)
        p=owned_position(tampered%grid_ids,1_int64)
        if(p>0)tampered%rt_space%grid_owner_keys(p)=tampered%rt_space%grid_owner_keys(p)+1
      case(70)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%metric_rows(p,1)=tampered%rt_space%metric_rows(p,1)+(0.01d0,0d0)
      case(71)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%kinetic_rows(p,1)=tampered%rt_space%kinetic_rows(p,1)+(0.01d0,0d0)
      case(72)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%nonlocal_rows(p,1)=tampered%rt_space%nonlocal_rows(p,1)+(0.01d0,0d0)
      case(73)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%local_rows(p,1)=tampered%rt_space%local_rows(p,1)+(0.01d0,0d0)
      case(74)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%sipg_rows(p,1)=tampered%rt_space%sipg_rows(p,1)+(0.01d0,0d0)
      case(75)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%hamiltonian_rows(p,1)=tampered%rt_space%hamiltonian_rows(p,1)+(0.01d0,0d0)
      case(76);tampered%rt_space%representation(1,1,1)=tampered%rt_space%representation(1,1,1)+(0.01d0,0d0)
      case(77);tampered%rt_space%cartesian_rotations(1,1,1)=tampered%rt_space%cartesian_rotations(1,1,1)+0.01d0
      case(78)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%scalar_operator_rows(p,1,1)=&
          tampered%rt_space%scalar_operator_rows(p,1,1)+(0.01d0,0d0)
      case(79)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%vector_operator_rows(p,1,1,1)=&
          tampered%rt_space%vector_operator_rows(p,1,1,1)+(0.01d0,0d0)
      case(80)
        p=owned_position(tampered%rt_space%row_ids,1_int64)
        if(p>0)tampered%rt_space%tensor_operator_rows(p,1,1,1,1)=&
          tampered%rt_space%tensor_operator_rows(p,1,1,1,1)+(0.01d0,0d0)
      case(81)
        p=owned_position(tampered%grid_ids,1_int64)
        if(p>0)tampered%rt_space%basis_values(1,p)=tampered%rt_space%basis_values(1,p)+(0.01d0,0d0)
      case(82)
        p=owned_position(tampered%grid_ids,1_int64)
        if(p>0)tampered%rt_space%density(p)=tampered%rt_space%density(p)+0.01d0
      case(83);tampered%rt_space%metric_fingerprint=tampered%rt_space%metric_fingerprint+1
      case(84);tampered%rt_space%kinetic_fingerprint=tampered%rt_space%kinetic_fingerprint+1
      case(85);tampered%rt_space%nonlocal_fingerprint=tampered%rt_space%nonlocal_fingerprint+1
      case(86);tampered%rt_space%local_fingerprint=tampered%rt_space%local_fingerprint+1
      case(87);tampered%rt_space%sipg_fingerprint=tampered%rt_space%sipg_fingerprint+1
      case(88);tampered%rt_space%hamiltonian_fingerprint=tampered%rt_space%hamiltonian_fingerprint+1
      case(89);tampered%rt_space%basis_fingerprint=tampered%rt_space%basis_fingerprint+1
      case(90);tampered%rt_space%density_fingerprint=tampered%rt_space%density_fingerprint+1
      case(91);tampered%rt_space%ownership_fingerprint=tampered%rt_space%ownership_fingerprint+1
      case(92);tampered%rt_space%scalar_fingerprint=tampered%rt_space%scalar_fingerprint+1
      case(93);tampered%rt_space%vector_fingerprint=tampered%rt_space%vector_fingerprint+1
      case(94);tampered%rt_space%tensor_fingerprint=tampered%rt_space%tensor_fingerprint+1
      case(95);tampered%rt_space%representation_fingerprint=tampered%rt_space%representation_fingerprint+1
      case(96);tampered%rt_space%fingerprint=tampered%rt_space%fingerprint+1
      case(97);tampered%energy_window%valid=.false.
      case(98);tampered%energy_window%compatibility_dynamic_rank=.true.
      case(99);tampered%energy_window%proof_state_present=.false.
      case(100);tampered%energy_window%mode=tampered%energy_window%mode+1
      case(101);tampered%energy_window%construction_rank=tampered%energy_window%construction_rank+1
      case(102);tampered%energy_window%solved_rank=tampered%energy_window%solved_rank+1
      case(103);tampered%energy_window%occupied_rank=tampered%energy_window%occupied_rank+1
      case(104);tampered%energy_window%requested_rank=tampered%energy_window%requested_rank+1
      case(105);tampered%energy_window%certified_rank=tampered%energy_window%certified_rank+1
      case(106);tampered%energy_window%extension_states=tampered%energy_window%extension_states+1
      case(107);tampered%energy_window%boundary_cluster_rank=tampered%energy_window%boundary_cluster_rank+1
      case(108);tampered%energy_window%proof_status=tampered%energy_window%proof_status+1
      case(109);tampered%energy_window%window_size=tampered%energy_window%window_size+0.01d0
      case(110);tampered%energy_window%e_homo=tampered%energy_window%e_homo+0.01d0
      case(111);tampered%energy_window%requested_cutoff=tampered%energy_window%requested_cutoff+0.01d0
      case(112);tampered%energy_window%certified_cutoff=tampered%energy_window%certified_cutoff+0.01d0
      case(113);tampered%energy_window%extension_energy=tampered%energy_window%extension_energy+0.01d0
      case(114);tampered%energy_window%proof_energy=tampered%energy_window%proof_energy+0.01d0
      case(115);tampered%energy_window%fingerprint=tampered%energy_window%fingerprint+1
      case(116);tampered%symmetry_receipt%valid=.false.
      case(117);tampered%symmetry_receipt%worst_operation=tampered%symmetry_receipt%worst_operation+1
      case(118);tampered%symmetry_receipt%occupied_subspace_defect=&
        tampered%symmetry_receipt%occupied_subspace_defect+1d-14
      case(119);tampered%symmetry_receipt%occupied_projector_defect=&
        tampered%symmetry_receipt%occupied_projector_defect+1d-14
      case(120);tampered%symmetry_receipt%target_subspace_defect=&
        tampered%symmetry_receipt%target_subspace_defect+1d-14
      case(121);tampered%symmetry_receipt%target_energy_defect=tampered%symmetry_receipt%target_energy_defect+1d-14
      case(122);tampered%symmetry_receipt%density_defect=tampered%symmetry_receipt%density_defect+1d-14
      case(123);tampered%symmetry_receipt%scalar_covariance_defect=&
        tampered%symmetry_receipt%scalar_covariance_defect+1d-14
      case(124);tampered%symmetry_receipt%vector_covariance_defect=&
        tampered%symmetry_receipt%vector_covariance_defect+1d-14
      case(125);tampered%symmetry_receipt%tensor_covariance_defect=&
        tampered%symmetry_receipt%tensor_covariance_defect+1d-14
      case(126);tampered%symmetry_receipt%final_basis_defect=tampered%symmetry_receipt%final_basis_defect+1d-14
      case(127);tampered%symmetry_receipt%worst_operation_defect=&
        tampered%symmetry_receipt%worst_operation_defect+1d-14
      case(128);tampered%symmetry_receipt%maximum_physical_defect=&
        tampered%symmetry_receipt%maximum_physical_defect+1d-14
      case(129);tampered%symmetry_receipt%fingerprint=tampered%symmetry_receipt%fingerprint+1
      case(130);tampered%handoff_receipts%valid=.false.
      case(131);tampered%handoff_receipts%position_fingerprint=tampered%handoff_receipts%position_fingerprint+1
      case(132);tampered%handoff_receipts%nonlocal_fingerprint=tampered%handoff_receipts%nonlocal_fingerprint+1
      case(133);tampered%handoff_receipts%face_fingerprint=tampered%handoff_receipts%face_fingerprint+1
      case(134);tampered%handoff_receipts%pseudopotential_fingerprint=&
        tampered%handoff_receipts%pseudopotential_fingerprint+1
      case(135);tampered%handoff_receipts%transformation_fingerprint=&
        tampered%handoff_receipts%transformation_fingerprint+1
      case(136);tampered%handoff_receipts%fingerprint=tampered%handoff_receipts%fingerprint+1
      case(137)
        if(allocated(tampered%certified_basis%c_cert))deallocate(tampered%certified_basis%c_cert)
      case(138)
        p=owned_position(tampered%certified_basis%transformation_row_ids,1_int64)
        if(p>0)tampered%certified_basis%transformation_row_ids(p)=2_int64
      case default
        error stop 'named tamper table incomplete'
      end select
      call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,&
        authentication_message)
      write(label,'(a,i0)')'named v3 tamper accepted for field ',field
      call require(.not.authenticated,trim(label))
    enddo
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_offsets(2)=tampered%face_offsets(2)+1
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted a changed face-point offset')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_weight_offsets(2)=tampered%face_weight_offsets(2)+1
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted a changed face-weight offset')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_basis_offsets(2)=tampered%face_basis_offsets(2)+1
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted a changed face-basis offset')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_value_offsets(2)=tampered%face_value_offsets(2)+1
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted a changed face-value offset')
    tampered=payload
    if(size(tampered%face_ids)>0)&
      tampered%face_observable_offsets(2)=tampered%face_observable_offsets(2)+1
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted a changed face-observable offset')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_point_ids=[tampered%face_point_ids,999999_int64]
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted an unowned face-point tail')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_weights=[tampered%face_weights,9d0]
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted an unowned face-weight tail')
    tampered=payload
    if(size(tampered%face_ids)>0)tampered%face_basis_ids=[tampered%face_basis_ids,n]
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted an unowned face-basis tail')
    tampered=payload
    if(size(tampered%face_ids)>0)call append_complex_column(tampered%face_values,[(9d0,1d0)])
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted an unowned face-value tail')
    tampered=payload
    if(size(tampered%face_ids)>0)call append_complex_column(tampered%interface_observables,&
      [(9d0,1d0),(8d0,2d0),(7d0,3d0)])
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'authentication accepted an unowned face-observable tail')
    tampered=payload
    p=owned_position(tampered%certified_basis%construction_row_ids,1_int64)
    if(p>0)then
      tampered%certified_basis%c_cert(p,1)=tampered%certified_basis%c_cert(p,1)+(0.01d0,0d0)
      tampered%certified_basis%b_rt(p,:)=matmul(tampered%certified_basis%c_cert(p,:),full_u_rt)
    endif
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'digest omitted a relation-preserving C_cert/B_rt change')
    tampered=payload
    p=owned_position(tampered%rt_space%row_ids,1_int64)
    if(p>0)then
      tampered%rt_space%kinetic_rows(p,1)=tampered%rt_space%kinetic_rows(p,1)+(0.01d0,0d0)
      tampered%rt_space%hamiltonian_rows(p,:)=tampered%rt_space%kinetic_rows(p,:)+&
        tampered%rt_space%nonlocal_rows(p,:)+tampered%rt_space%local_rows(p,:)+tampered%rt_space%sipg_rows(p,:)
    endif
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'digest omitted a relation-preserving kinetic/Hamiltonian change')
    if(nproc==1)then
      tampered=payload
      swap_row=tampered%certified_basis%c_cert(1,:)
      tampered%certified_basis%c_cert(1,:)=tampered%certified_basis%c_cert(2,:)
      tampered%certified_basis%c_cert(2,:)=swap_row
      swap_row=tampered%certified_basis%b_rt(1,:)
      tampered%certified_basis%b_rt(1,:)=tampered%certified_basis%b_rt(2,:)
      tampered%certified_basis%b_rt(2,:)=swap_row
      call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
      call require(.not.authenticated,'digest accepted a relation-preserving two-row swap')
    endif
    tampered=payload
    tampered%certified_basis%embedding_defect=-1d-14
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'negative certified-basis defect was accepted')
    tampered=payload
    if(size(tampered%metric_row_offsets)>1)&
      tampered%metric_row_offsets(2)=size(tampered%metric_column_ids)+2
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'out-of-range CSR row offset was accepted')
    tampered=payload
    p=owned_position(tampered%rt_space%row_ids,1_int64)
    if(p>0)then
      tampered%rt_space%kinetic_rows(p,1)=cmplx(huge(0d0)/2d0,0d0,real64)
      tampered%rt_space%nonlocal_rows(p,1)=cmplx(huge(0d0)/2d0,0d0,real64)
      tampered%rt_space%hamiltonian_rows(p,1)=cmplx(huge(0d0),0d0,real64)
    endif
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,tampered,fingerprint,authenticated,authentication_message)
    call require(.not.authenticated,'overflowing Hamiltonian component sum was evaluated')
  end subroutine verify_named_tamper_rejection

  subroutine append_complex_column(values,column)
    complex(real64),allocatable,intent(inout)::values(:,:)
    complex(real64),intent(in)::column(:)
    complex(real64),allocatable::next_values(:,:)
    integer::old_columns
    old_columns=size(values,2)
    if(size(column)/=size(values,1))error stop 'invalid complex-column test fixture'
    allocate(next_values(size(values,1),old_columns+1))
    if(old_columns>0)next_values(:,:old_columns)=values
    next_values(:,old_columns+1)=column
    call move_alloc(next_values,values)
  end subroutine append_complex_column

  integer function owned_position(ids,id) result(position)
    integer(int64),intent(in)::ids(:),id
    integer::j
    position=0
    do j=1,size(ids)
      if(ids(j)==id)then;position=j;return;endif
    enddo
  end function owned_position

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
