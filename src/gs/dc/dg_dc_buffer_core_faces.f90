#include "config.h"
module dg_dc_buffer_core_faces
#ifdef USE_MPI
  use mpi
#endif
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use dg_buffer_window_projector,only:s_dg_buffer_projector_diagnostics,&
    build_dg_buffer_window_projector
  implicit none
  private
  integer(8),parameter :: maximum_quarter_mpi_count=536870911_8

  type,public :: s_dg_dc_buffer_core_face
    integer :: axis=0,side=0,local_fragment=0,neighbor_fragment=0,neighbor_rank=-1
    integer :: overlap_depth=0,point_count=0,generation=0
    integer(8),allocatable :: physical_grid_ids(:)
    integer,allocatable :: local_buffer_indices(:,:),neighbor_core_indices(:,:)
    real(8),allocatable :: neighbor_core_values(:,:),local_core_values(:,:)
    real(8),allocatable :: neighbor_core_normals(:,:),local_core_normals(:,:)
    real(8),allocatable :: local_buffer_values(:,:)
    real(8),allocatable :: coefficients(:,:),reconstructed_values(:,:),reconstructed_normals(:,:)
    type(s_dg_buffer_projector_diagnostics) :: diagnostics
    integer(8) :: operator_fingerprint=0_8
    real(8) :: forward_reverse_projection_mismatch=huge(1d0)
    logical :: physical_map_identity=.false.
  end type

  public :: build_dg_dc_buffer_core_faces,exchange_dg_dc_buffer_core_state_window
  public :: assemble_dg_dc_local_buffer_state_window
  public :: validate_dg_dc_buffer_core_faces,project_dg_dc_buffer_core_face
contains
  subroutine build_dg_dc_buffer_core_faces(local_fragment,fragment_origins,fragment_sizes,global_size,&
      buffer_width,generation,communicator,faces,ok,message,local_core_origin,local_core_size)
    integer,intent(in) :: local_fragment,fragment_origins(:,:),fragment_sizes(:,:),global_size(3)
    integer,intent(in) :: buffer_width,generation,communicator
    type(s_dg_dc_buffer_core_face),allocatable,intent(out) :: faces(:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer,optional,intent(in) :: local_core_origin(3),local_core_size(3)
    integer :: axis,side,iface,allocation_status,local_bad,global_bad,owned_origin(3),owned_size(3)
#ifdef USE_MPI
    integer :: ierr
#endif

    ok=.false.;message=''
    local_bad=0
    if(size(fragment_origins,1)/=3)local_bad=1
    if(any(shape(fragment_origins)/=shape(fragment_sizes)))local_bad=1
    if(local_fragment<1.or.local_fragment>size(fragment_origins,2))local_bad=1
    if(any(global_size<=0).or.any(fragment_sizes<=0).or.buffer_width<=0.or.generation<=0)local_bad=1
    if(any(fragment_origins<0))local_bad=1
    if(local_bad==0)then
      if(any(int(fragment_origins,8)+int(fragment_sizes,8)+int(buffer_width,8)>&
        int(huge(buffer_width),8)))local_bad=1
    endif
    owned_origin=0;owned_size=0
    if(local_bad==0)then
      owned_origin=fragment_origins(:,local_fragment)
      owned_size=fragment_sizes(:,local_fragment)
      if(present(local_core_origin))owned_origin=local_core_origin
      if(present(local_core_size))owned_size=local_core_size
      if(any(owned_size<=0))local_bad=1
    endif
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='DG buffer/core faces: invalid collective mapping context';return
    endif
    call validate_collective_mapping_metadata(fragment_origins,fragment_sizes,global_size,&
      buffer_width,generation,communicator,local_bad)
    if(local_bad/=0)then
      message='DG buffer/core faces: inconsistent collective mapping metadata';return
    endif
#else
    if(local_bad/=0.or.communicator<0.or.size(fragment_origins,2)/=1)then
      message='DG buffer/core faces: invalid serial mapping context';return
    endif
#endif
    call validate_partition(fragment_origins,fragment_sizes,global_size,local_bad)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then
      message='DG buffer/core faces: invalid, overlapping, or incomplete physical partition';return
    endif
    allocate(faces(6),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then
      message='DG buffer/core faces: face allocation failed collectively';return
    endif
    do axis=1,3
      do side=-1,1,2
        iface=face_slot(axis,side)
        call build_one_face(local_fragment,axis,side,fragment_origins,fragment_sizes,global_size,&
          buffer_width,generation,owned_origin,owned_size,faces(iface),local_bad)
        call collective_max(local_bad,communicator,global_bad)
        if(global_bad/=0)then
          message='DG buffer/core faces: physical face construction failed collectively';return
        endif
      enddo
    enddo
    ok=.true.
  end subroutine

#ifdef USE_MPI
  subroutine validate_collective_mapping_metadata(origins,sizes,global_size,buffer_width,generation,&
      communicator,bad)
    integer,intent(in) :: origins(:,:),sizes(:,:),global_size(3),buffer_width,generation,communicator
    integer,intent(out) :: bad
    integer :: ierr,nfragment_min,nfragment_max,scalar_min,scalar_max,allocation_status,global_bad
    integer,allocatable :: minimum(:),maximum(:),payload(:)
    bad=0
    call MPI_Allreduce(size(origins,2),nfragment_min,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
    call MPI_Allreduce(size(origins,2),nfragment_max,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.nfragment_min/=nfragment_max)then;bad=1;return;endif
    call MPI_Allreduce(buffer_width,scalar_min,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Allreduce(buffer_width,scalar_max,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.scalar_min/=scalar_max)bad=1
    call MPI_Allreduce(generation,scalar_min,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Allreduce(generation,scalar_max,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.scalar_min/=scalar_max)bad=1
    call MPI_Allreduce(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.bad/=0)return
    allocate(payload(3+6*nfragment_min),minimum(3+6*nfragment_min),&
      maximum(3+6*nfragment_min),stat=allocation_status)
    bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;bad=1;return;endif
    payload(1:3)=global_size
    payload(4:3+3*nfragment_min)=reshape(origins,[3*nfragment_min])
    payload(4+3*nfragment_min:3+6*nfragment_min)=reshape(sizes,[3*nfragment_min])
    call MPI_Allreduce(payload,minimum,size(payload),MPI_INTEGER,MPI_MIN,communicator,ierr)
    if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Allreduce(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.bad/=0)return
    call MPI_Allreduce(payload,maximum,size(payload),MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum))bad=1
  end subroutine
#endif

  subroutine build_one_face(local_fragment,axis,side,origins,sizes,global_size,buffer_width,generation,&
      owned_origin,owned_size,face,bad)
    integer,intent(in) :: local_fragment,axis,side,origins(:,:),sizes(:,:),global_size(3)
    integer,intent(in) :: buffer_width,generation
    integer,intent(in) :: owned_origin(3),owned_size(3)
    type(s_dg_dc_buffer_core_face),intent(out) :: face
    integer,intent(out) :: bad
    integer :: point_count,allocation_status,depth,t1,t2,tangent(2),point,position(3),wrapped(3)
    integer :: neighbor,found_neighbor
    integer(8) :: point_count_8
    tangent=pack([1,2,3],[1,2,3]/=axis)
    point_count_8=0_8
    if((side<0.and.owned_origin(axis)==origins(axis,local_fragment)).or.&
       (side>0.and.owned_origin(axis)+owned_size(axis)==&
        origins(axis,local_fragment)+sizes(axis,local_fragment)))then
      point_count_8=checked_grid_point_count([buffer_width,owned_size(tangent(1)),owned_size(tangent(2))])
    endif
    if(point_count_8<0_8.or.point_count_8>int(huge(point_count),8))then;bad=1;return;endif
    point_count=int(point_count_8)
    face%axis=axis;face%side=side;face%local_fragment=local_fragment
    face%overlap_depth=buffer_width;face%point_count=point_count;face%generation=generation
    allocate(face%physical_grid_ids(point_count),face%local_buffer_indices(3,point_count),&
      face%neighbor_core_indices(3,point_count),stat=allocation_status)
    bad=merge(0,1,allocation_status==0)
    if(bad/=0)return
    position=origins(:,local_fragment)
    position(axis)=merge(origins(axis,local_fragment)-1,&
      origins(axis,local_fragment)+sizes(axis,local_fragment),side<0)
    wrapped=modulo(position,global_size)
    found_neighbor=find_owner(wrapped,origins,sizes)
    if(found_neighbor<=0)then;bad=1;return;endif
    point=0
    if(point_count>0)then
    do depth=1,buffer_width
      do t2=0,owned_size(tangent(2))-1
        do t1=0,owned_size(tangent(1))-1
          point=point+1
          position=owned_origin
          position(tangent(1))=owned_origin(tangent(1))+t1
          position(tangent(2))=owned_origin(tangent(2))+t2
          if(side<0)then
            position(axis)=origins(axis,local_fragment)-depth
            face%local_buffer_indices(axis,point)=buffer_width+1-depth
          else
            position(axis)=origins(axis,local_fragment)+sizes(axis,local_fragment)-1+depth
            face%local_buffer_indices(axis,point)=buffer_width+sizes(axis,local_fragment)+depth
          endif
          face%local_buffer_indices(tangent(1),point)=buffer_width+&
            position(tangent(1))-origins(tangent(1),local_fragment)+1
          face%local_buffer_indices(tangent(2),point)=buffer_width+&
            position(tangent(2))-origins(tangent(2),local_fragment)+1
          wrapped=modulo(position,global_size)
          neighbor=find_owner(wrapped,origins,sizes)
          if(neighbor<=0)then;bad=1;return;endif
          if(neighbor/=found_neighbor)then;bad=1;return;endif
          face%physical_grid_ids(point)=physical_id(wrapped,global_size)
          face%neighbor_core_indices(:,point)=wrapped-origins(:,neighbor)+1
        enddo
      enddo
    enddo
    endif
    face%neighbor_fragment=found_neighbor
    face%neighbor_rank=-1
  end subroutine

  subroutine exchange_dg_dc_buffer_core_state_window(faces,fragment_origins,fragment_sizes,global_size,&
      local_core_origin,local_core_size,core_physical_grid_ids,core_values,core_gradients,&
      owned_state_ids,configured_states,communicator,ok,message)
    type(s_dg_dc_buffer_core_face),intent(inout) :: faces(:)
    integer,intent(in) :: fragment_origins(:,:),fragment_sizes(:,:)
    integer,intent(in) :: global_size(3),local_core_origin(3),local_core_size(3)
    integer(8),intent(in) :: core_physical_grid_ids(:,:,:)
    real(8),intent(in) :: core_values(:,:,:,:),core_gradients(:,:,:,:,:)
    integer,intent(in) :: owned_state_ids(:)
    integer,intent(in) :: configured_states,communicator
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: local_bad,global_bad,iface,point,state,allocation_status
    integer :: request_count,request,rank_count,rank_index,x,y,z
    integer(8) :: request_count_8,point_total_8,state_total_8
    type(s_dg_dc_buffer_core_face),allocatable :: candidate_windows(:)
    integer,allocatable :: rank_fragments(:),rank_core_origins(:,:),rank_core_sizes(:,:)
    integer,allocatable :: state_counts(:),state_displacements(:),request_owners(:)
    integer,allocatable :: request_faces(:),request_points(:),request_states(:)
    integer,allocatable :: owners_for_states(:)
#ifdef USE_MPI
    integer :: rank,ierr,min_states,max_states,position,owner
    integer(8) :: total_send_8,total_receive_8
    integer,allocatable :: send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:)
    integer,allocatable :: send_counts_2(:),receive_counts_2(:),send_displacements_2(:),receive_displacements_2(:)
    integer,allocatable :: send_counts_4(:),receive_counts_4(:),send_displacements_4(:),receive_displacements_4(:)
    integer,allocatable :: cursor(:),directory_state_ids(:)
    integer(8),allocatable :: request_send(:),request_receive(:)
    real(8),allocatable :: response_send(:),response_receive(:)
#else
    integer,allocatable :: directory_state_ids(:)
    real(8) :: local_record(4)
#endif
    ok=.false.;message=''
    local_bad=0
    if(size(faces)/=6.or.configured_states<=0.or.any(global_size<=0).or.any(local_core_size<=0))local_bad=1
    if(size(core_values,4)/=size(owned_state_ids))local_bad=1
    if(any(local_core_size/=[size(core_values,1),size(core_values,2),size(core_values,3)]))local_bad=1
    if(size(core_gradients,4)/=size(owned_state_ids).or.size(core_gradients,5)/=3)local_bad=1
    if(any([size(core_gradients,1),size(core_gradients,2),size(core_gradients,3)]/=&
      [size(core_values,1),size(core_values,2),size(core_values,3)]))local_bad=1
    if(any(shape(core_physical_grid_ids)/=&
      [size(core_values,1),size(core_values,2),size(core_values,3)]))local_bad=1
    if(any(owned_state_ids<1).or.any(owned_state_ids>configured_states))local_bad=1
    if(.not.all(ieee_is_finite(core_values)).or..not.all(ieee_is_finite(core_gradients)))local_bad=1
    if(local_bad==0)then
      do state=1,size(owned_state_ids)
        if(count(owned_state_ids==owned_state_ids(state))/=1)local_bad=1
      enddo
      do z=1,size(core_physical_grid_ids,3);do y=1,size(core_physical_grid_ids,2)
        do x=1,size(core_physical_grid_ids,1)
          if(core_physical_grid_ids(x,y,z)<1_8)local_bad=1
        enddo
      enddo;enddo
    endif
#ifdef USE_MPI
    call MPI_Allreduce(configured_states,min_states,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(configured_states,max_states,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS.or.min_states/=max_states)local_bad=1
#endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then
      message='DG buffer/core faces: inconsistent configured-state dimensions before communication';return
    endif
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(communicator,rank_count,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
#else
    rank_count=1
#endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: communicator query failed collectively';return;endif
    allocate(rank_fragments(rank_count),rank_core_origins(3,rank_count),rank_core_sizes(3,rank_count),&
      state_counts(rank_count),state_displacements(rank_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: directory allocation failed collectively';return;endif
#ifdef USE_MPI
    call MPI_Allgather(faces(1)%local_fragment,1,MPI_INTEGER,rank_fragments,1,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(local_core_origin,3,MPI_INTEGER,rank_core_origins,3,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(local_core_size,3,MPI_INTEGER,rank_core_sizes,3,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(owned_state_ids),1,MPI_INTEGER,state_counts,1,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
#else
    rank_fragments=faces(1)%local_fragment
    rank_core_origins(:,1)=local_core_origin;rank_core_sizes(:,1)=local_core_size
    state_counts=size(owned_state_ids)
#endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: ownership descriptor exchange failed';return;endif
    state_total_8=sum(int(state_counts,8))
    local_bad=merge(0,1,state_total_8>=0_8.and.state_total_8<=int(huge(request_count),8))
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: state directory exceeds integer capacity';return;endif
    state_displacements(1)=0
    do rank_index=2,rank_count
      state_displacements(rank_index)=state_displacements(rank_index-1)+state_counts(rank_index-1)
    enddo
    allocate(directory_state_ids(sum(state_counts)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: directory payload allocation failed collectively';return;endif
#ifdef USE_MPI
    call MPI_Allgatherv(owned_state_ids,size(owned_state_ids),MPI_INTEGER,directory_state_ids,&
      state_counts,state_displacements,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
#else
    directory_state_ids=owned_state_ids
#endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: ownership directory exchange failed';return;endif
    call validate_distributed_core_ownership(fragment_origins,fragment_sizes,rank_fragments,&
      rank_core_origins,rank_core_sizes,state_counts,state_displacements,directory_state_ids,&
      configured_states,local_bad)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then
      message='DG buffer/core faces: incomplete or overlapping distributed core ownership';return
    endif
    point_total_8=sum([(int(faces(iface)%point_count,8),iface=1,6)])
    if(point_total_8>0_8.and.int(configured_states,8)>huge(request_count_8)/point_total_8)then
      request_count_8=-1_8
    else
      request_count_8=int(configured_states,8)*point_total_8
    endif
    local_bad=merge(0,1,request_count_8>=0_8.and.request_count_8<=maximum_quarter_mpi_count)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: request count exceeds MPI integer capacity';return;endif
    request_count=int(request_count_8)
    allocate(request_owners(request_count),request_faces(request_count),request_points(request_count),&
      request_states(request_count),owners_for_states(configured_states),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    if(local_bad==0)then
      request=0
      do iface=1,6;do point=1,faces(iface)%point_count
        call index_record_owners(faces(iface)%neighbor_fragment,faces(iface)%physical_grid_ids(point),&
          global_size,rank_fragments,rank_core_origins,rank_core_sizes,state_counts,state_displacements,&
          directory_state_ids,owners_for_states,local_bad)
        do state=1,configured_states
          request=request+1
          request_owners(request)=owners_for_states(state)
          request_faces(request)=iface;request_points(request)=point;request_states(request)=state
        enddo
      enddo;enddo
    endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then
      message='DG buffer/core faces: missing or duplicate physical-grid/state owner';return
    endif
    allocate(candidate_windows(6),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    do iface=1,6
      if(local_bad/=0)exit
      allocate(candidate_windows(iface)%neighbor_core_values(faces(iface)%point_count,configured_states),&
        candidate_windows(iface)%neighbor_core_normals(faces(iface)%point_count,configured_states),&
        stat=allocation_status)
      if(allocation_status/=0)local_bad=1
      if(allocation_status==0)then
        candidate_windows(iface)%neighbor_core_values=0d0
        candidate_windows(iface)%neighbor_core_normals=0d0
      endif
    enddo
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: state-window allocation failed collectively';return;endif
#ifdef USE_MPI
    allocate(send_counts(rank_count),receive_counts(rank_count),send_displacements(rank_count),&
      receive_displacements(rank_count),cursor(rank_count),stat=allocation_status)
    if(allocation_status/=0)then;local_bad=1;else;send_counts=0;endif
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: request count allocation failed';return;endif
    do request=1,request_count
      send_counts(request_owners(request)+1)=send_counts(request_owners(request)+1)+1
    enddo
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,receive_counts,1,MPI_INTEGER,communicator,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: request count exchange failed';return;endif
    total_send_8=sum(int(send_counts,8));total_receive_8=sum(int(receive_counts,8))
    local_bad=merge(0,1,total_send_8<=maximum_quarter_mpi_count.and.&
      total_receive_8<=maximum_quarter_mpi_count)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: response count exceeds MPI integer capacity';return;endif
    send_displacements(1)=0;receive_displacements(1)=0
    do rank_index=2,rank_count
      send_displacements(rank_index)=send_displacements(rank_index-1)+send_counts(rank_index-1)
      receive_displacements(rank_index)=receive_displacements(rank_index-1)+receive_counts(rank_index-1)
    enddo
    allocate(request_send(2*sum(send_counts)),request_receive(2*sum(receive_counts)),&
      response_send(4*sum(receive_counts)),response_receive(4*sum(send_counts)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: payload allocation failed collectively';return;endif
      cursor=send_displacements
      do request=1,request_count
        owner=request_owners(request)+1;position=cursor(owner)+1;cursor(owner)=position
        request_send(2*position-1)=faces(request_faces(request))%physical_grid_ids(request_points(request))
        request_send(2*position)=int(request_states(request),8)
      enddo
      allocate(send_counts_2(rank_count),receive_counts_2(rank_count),send_displacements_2(rank_count),&
        receive_displacements_2(rank_count),send_counts_4(rank_count),receive_counts_4(rank_count),&
        send_displacements_4(rank_count),receive_displacements_4(rank_count),stat=allocation_status)
      local_bad=merge(0,1,allocation_status==0)
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: count allocation failed collectively';return;endif
      send_counts_2=2*send_counts;receive_counts_2=2*receive_counts
      send_displacements_2=2*send_displacements;receive_displacements_2=2*receive_displacements
      send_counts_4=4*receive_counts;receive_counts_4=4*send_counts
      send_displacements_4=4*receive_displacements;receive_displacements_4=4*send_displacements
      call MPI_Alltoallv(request_send,send_counts_2,send_displacements_2,MPI_INTEGER8,&
        request_receive,receive_counts_2,receive_displacements_2,MPI_INTEGER8,communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: request payload exchange failed';return;endif
      do request=1,sum(receive_counts)
        call lookup_local_record(request_receive(2*request-1),int(request_receive(2*request)),&
          global_size,local_core_origin,core_physical_grid_ids,owned_state_ids,core_values,core_gradients,&
          response_send(4*request-3:4*request),local_bad)
        if(local_bad/=0)exit
      enddo
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: local record lookup failed collectively';return;endif
      call MPI_Alltoallv(response_send,send_counts_4,send_displacements_4,MPI_DOUBLE_PRECISION,&
        response_receive,receive_counts_4,receive_displacements_4,MPI_DOUBLE_PRECISION,communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: response payload exchange failed';return;endif
      cursor=send_displacements
      do request=1,request_count
        owner=request_owners(request)+1;position=cursor(owner)+1;cursor(owner)=position
        iface=request_faces(request);point=request_points(request);state=request_states(request)
        candidate_windows(iface)%neighbor_core_values(point,state)=response_receive(4*position-3)
        candidate_windows(iface)%neighbor_core_normals(point,state)=-dble(faces(iface)%side)*&
          response_receive(4*position-3+faces(iface)%axis)
      enddo
#else
    do request=1,request_count
      iface=request_faces(request);point=request_points(request);state=request_states(request)
      call lookup_local_record(faces(iface)%physical_grid_ids(point),state,global_size,local_core_origin,&
        core_physical_grid_ids,owned_state_ids,core_values,core_gradients,local_record,local_bad)
      candidate_windows(iface)%neighbor_core_values(point,state)=local_record(1)
      candidate_windows(iface)%neighbor_core_normals(point,state)=-dble(faces(iface)%side)*&
        local_record(1+faces(iface)%axis)
    enddo
#endif
    call collective_max(local_bad,communicator,global_bad)
    ok=global_bad==0
    if(ok)then
      do iface=1,6
        call move_alloc(candidate_windows(iface)%neighbor_core_values,faces(iface)%neighbor_core_values)
        call move_alloc(candidate_windows(iface)%neighbor_core_normals,faces(iface)%neighbor_core_normals)
      enddo
      message=''
    else
      message='DG buffer/core faces: complete state-window exchange failed collectively'
    endif
  end subroutine

  subroutine assemble_dg_dc_local_buffer_state_window(faces,owned_face_values,owned_state_ids,&
      configured_states,communicator,ok,message)
    type(s_dg_dc_buffer_core_face),intent(inout) :: faces(:)
    real(8),intent(in) :: owned_face_values(:,:,:)
    integer,intent(in) :: owned_state_ids(:),configured_states,communicator
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: iface,state,allocation_status,local_bad,global_bad
    integer :: buffer_mpi_count
    integer(8) :: buffer_count_8
    type(s_dg_dc_buffer_core_face),allocatable :: candidate_buffers(:)
#ifdef USE_MPI
    integer :: rank,ierr,color,face_communicator,face_rank_count,member,rank_count,candidate,safe_gather_count
    integer(8) :: face_total_8
    integer,allocatable :: state_coverage(:),face_counts(:),face_displacements(:),face_fragments(:)
    integer(8),allocatable :: gathered_ids(:),all_face_ids(:)
#endif
    ok=.false.;message='';local_bad=0
    if(size(faces)/=6.or.configured_states<=0)local_bad=1
    if(size(owned_face_values,2)/=size(owned_state_ids).or.size(owned_face_values,3)/=6)local_bad=1
    if(any(owned_state_ids<1).or.any(owned_state_ids>configured_states))local_bad=1
    if(.not.all(ieee_is_finite(owned_face_values)))local_bad=1
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: invalid owned local-buffer window';return;endif
    allocate(candidate_buffers(6),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: buffer candidate allocation failed';return;endif
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(communicator,rank_count,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call collective_max(local_bad,communicator,global_bad)
    if(global_bad/=0)then;message='DG buffer/core faces: communicator query failed';return;endif
#endif
    do iface=1,6
#ifdef USE_MPI
      color=MPI_UNDEFINED
      allocate(face_counts(rank_count),face_displacements(rank_count),face_fragments(rank_count),&
        stat=allocation_status)
      local_bad=merge(0,1,allocation_status==0)
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: exact face-group allocation failed';return;endif
      call MPI_Allgather(faces(iface)%point_count,1,MPI_INTEGER,face_counts,1,MPI_INTEGER,communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call MPI_Allgather(faces(iface)%local_fragment,1,MPI_INTEGER,face_fragments,1,MPI_INTEGER,communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      face_total_8=sum(int(face_counts,8))
      if(face_total_8<0_8.or.face_total_8>int(huge(iface),8))local_bad=1
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: face metadata exceeds MPI capacity';return;endif
      face_displacements(1)=0
      do member=2,rank_count
        face_displacements(member)=face_displacements(member-1)+face_counts(member-1)
      enddo
      allocate(all_face_ids(sum(face_counts)),stat=allocation_status)
      if(allocation_status/=0)local_bad=1
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: exact face-group payload allocation failed';return;endif
      call MPI_Allgatherv(faces(iface)%physical_grid_ids,faces(iface)%point_count,MPI_INTEGER8,&
        all_face_ids,face_counts,face_displacements,MPI_INTEGER8,communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: exact face metadata exchange failed';return;endif
      if(faces(iface)%point_count>0)then
        color=rank+1
        do candidate=1,rank_count
          if(face_fragments(candidate)/=faces(iface)%local_fragment.or.&
            face_counts(candidate)/=faces(iface)%point_count)cycle
          if(all(all_face_ids(face_displacements(candidate)+1:&
            face_displacements(candidate)+face_counts(candidate))==faces(iface)%physical_grid_ids))then
            color=candidate;exit
          endif
        enddo
      endif
      call MPI_Comm_split(communicator,color,rank,face_communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      deallocate(face_counts,face_displacements,face_fragments,all_face_ids)
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: face-group split failed collectively';return;endif
      if(faces(iface)%point_count>0)then
      call MPI_Comm_size(face_communicator,face_rank_count,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      endif
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then;message='DG buffer/core faces: face communicator query failed';return;endif
      if(faces(iface)%point_count>0)then
      face_total_8=int(faces(iface)%point_count,8)*int(face_rank_count,8)
      safe_gather_count=0
      if(face_total_8>int(huge(safe_gather_count),8))local_bad=1
      if(local_bad==0)safe_gather_count=int(face_total_8)
      if(allocated(gathered_ids))deallocate(gathered_ids)
      if(allocated(state_coverage))deallocate(state_coverage)
      allocate(gathered_ids(safe_gather_count),&
        state_coverage(configured_states),stat=allocation_status)
      local_bad=merge(0,1,allocation_status==0)
      call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,face_communicator,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      if(local_bad==0)then
        call MPI_Allgather(faces(iface)%physical_grid_ids,faces(iface)%point_count,MPI_INTEGER8,&
          gathered_ids,faces(iface)%point_count,MPI_INTEGER8,face_communicator,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        do member=0,face_rank_count-1
          if(any(gathered_ids(member*faces(iface)%point_count+1:&
            (member+1)*faces(iface)%point_count)/=faces(iface)%physical_grid_ids))local_bad=1
        enddo
        state_coverage=0
        do state=1,size(owned_state_ids)
          state_coverage(owned_state_ids(state))=state_coverage(owned_state_ids(state))+1
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,state_coverage,configured_states,MPI_INTEGER,MPI_SUM,&
          face_communicator,ierr)
        if(ierr/=MPI_SUCCESS.or.any(state_coverage/=1))local_bad=1
      endif
#else
      if(faces(iface)%point_count==0)cycle
      if(size(owned_state_ids)/=configured_states)local_bad=1
#endif
      if(size(owned_face_values,1)<faces(iface)%point_count)local_bad=1
      buffer_count_8=int(faces(iface)%point_count,8)*int(configured_states,8)
      buffer_mpi_count=0
      if(buffer_count_8>int(huge(buffer_mpi_count),8))local_bad=1
      if(local_bad==0)buffer_mpi_count=int(buffer_count_8)
      allocate(candidate_buffers(iface)%local_buffer_values(faces(iface)%point_count,configured_states),&
        stat=allocation_status)
      if(allocation_status/=0)local_bad=1
#ifdef USE_MPI
      call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,face_communicator,ierr)
#endif
      if(local_bad==0)then
        candidate_buffers(iface)%local_buffer_values=0d0
        do state=1,size(owned_state_ids)
          candidate_buffers(iface)%local_buffer_values(:,owned_state_ids(state))=&
            owned_face_values(1:faces(iface)%point_count,state,iface)
        enddo
#ifdef USE_MPI
        call MPI_Allreduce(MPI_IN_PLACE,candidate_buffers(iface)%local_buffer_values,&
          buffer_mpi_count,MPI_DOUBLE_PRECISION,MPI_SUM,face_communicator,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
#endif
      endif
#ifdef USE_MPI
      call MPI_Comm_free(face_communicator,ierr)
      if(allocated(gathered_ids))deallocate(gathered_ids)
      if(allocated(state_coverage))deallocate(state_coverage)
      endif
      call collective_max(local_bad,communicator,global_bad)
      if(global_bad/=0)then
        message='DG buffer/core faces: face-group assembly failed collectively';return
      endif
#endif
    enddo
    call collective_max(local_bad,communicator,global_bad)
    ok=global_bad==0
    if(ok)then
      do iface=1,6
        if(allocated(candidate_buffers(iface)%local_buffer_values))&
          call move_alloc(candidate_buffers(iface)%local_buffer_values,faces(iface)%local_buffer_values)
      enddo
      message=''
    else
      message='DG buffer/core faces: local-buffer window assembly failed collectively'
    endif
  end subroutine

  subroutine validate_distributed_core_ownership(fragment_origins,fragment_sizes,rank_fragments,&
      rank_origins,rank_sizes,state_counts,state_displacements,state_ids,configured_states,bad)
    integer,intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),rank_fragments(:)
    integer,intent(in) :: rank_origins(:,:),rank_sizes(:,:),state_counts(:),state_displacements(:)
    integer,intent(in) :: state_ids(:),configured_states
    integer,intent(out) :: bad
    integer :: fragment,representative,other,rank_index,state,allocation_status
    integer,allocatable :: coverage(:)
    integer(8) :: covered_volume,fragment_volume,box_volume
    logical :: same_box
    bad=0
    if(size(fragment_origins,1)/=3.or.any(shape(fragment_origins)/=shape(fragment_sizes)))then
      bad=1;return
    endif
    allocate(coverage(configured_states),stat=allocation_status)
    if(allocation_status/=0)then;bad=1;return;endif
    do fragment=1,size(fragment_origins,2)
      covered_volume=0_8
      fragment_volume=checked_grid_point_count(fragment_sizes(:,fragment))
      if(fragment_volume<1_8)bad=1
      do representative=1,size(rank_fragments)
        if(rank_fragments(representative)/=fragment)cycle
        same_box=.false.
        do other=1,representative-1
          if(rank_fragments(other)==fragment.and.all(rank_origins(:,other)==rank_origins(:,representative)).and.&
            all(rank_sizes(:,other)==rank_sizes(:,representative)))same_box=.true.
        enddo
        if(same_box)cycle
        if(any(rank_sizes(:,representative)<=0).or..not.box_within_box(rank_origins(:,representative),&
          rank_sizes(:,representative),fragment_origins(:,fragment),fragment_sizes(:,fragment)))bad=1
        coverage=0
        do rank_index=1,size(rank_fragments)
          if(rank_fragments(rank_index)/=fragment)cycle
          if(.not.all(rank_origins(:,rank_index)==rank_origins(:,representative)).or.&
            .not.all(rank_sizes(:,rank_index)==rank_sizes(:,representative)))cycle
          do state=1,state_counts(rank_index)
            if(state_ids(state_displacements(rank_index)+state)<1.or.&
              state_ids(state_displacements(rank_index)+state)>configured_states)then
              bad=1
            else
              coverage(state_ids(state_displacements(rank_index)+state))=&
                coverage(state_ids(state_displacements(rank_index)+state))+1
            endif
          enddo
        enddo
        if(any(coverage/=1))bad=1
        do other=representative+1,size(rank_fragments)
          if(rank_fragments(other)/=fragment)cycle
          if(all(rank_origins(:,other)==rank_origins(:,representative)).and.&
            all(rank_sizes(:,other)==rank_sizes(:,representative)))cycle
          if(boxes_overlap(rank_origins(:,representative),rank_sizes(:,representative),&
            rank_origins(:,other),rank_sizes(:,other)))bad=1
        enddo
        box_volume=checked_grid_point_count(rank_sizes(:,representative))
        if(box_volume<1_8.or.covered_volume>huge(covered_volume)-max(box_volume,0_8))then
          bad=1
        else
          covered_volume=covered_volume+box_volume
        endif
      enddo
      if(covered_volume/=fragment_volume)bad=1
    enddo
  end subroutine

  subroutine index_record_owners(fragment,identifier,global_size,rank_fragments,rank_origins,&
      rank_sizes,state_counts,state_displacements,state_ids,owners,bad)
    integer,intent(in) :: fragment,global_size(3),rank_fragments(:),rank_origins(:,:),rank_sizes(:,:)
    integer,intent(in) :: state_counts(:),state_displacements(:),state_ids(:)
    integer(8),intent(in) :: identifier
    integer,intent(out) :: owners(:)
    integer,intent(inout) :: bad
    integer :: candidate,grid_position(3),state_index,state
    call position_from_physical_id(identifier,global_size,grid_position)
    owners=-1
    do candidate=1,size(rank_fragments)
      if(rank_fragments(candidate)/=fragment)cycle
      if(.not.point_in_box(grid_position,rank_origins(:,candidate),rank_sizes(:,candidate)))cycle
      do state_index=1,state_counts(candidate)
        state=state_ids(state_displacements(candidate)+state_index)
        if(state<1.or.state>size(owners))then;bad=1;cycle;endif
        if(owners(state)>=0)then;bad=1;else;owners(state)=candidate-1;endif
      enddo
    enddo
    if(any(owners<0))bad=1
  end subroutine

  subroutine lookup_local_record(identifier,state,global_size,local_origin,identifiers,state_ids,&
      values,gradients,record,bad)
    integer(8),intent(in) :: identifier
    integer,intent(in) :: state,global_size(3),local_origin(3),state_ids(:)
    integer(8),intent(in) :: identifiers(:,:,:)
    real(8),intent(in) :: values(:,:,:,:),gradients(:,:,:,:,:)
    real(8),intent(out) :: record(4)
    integer,intent(inout) :: bad
    integer :: column,index(3),grid_position(3),state_index
    record=0d0
    column=0
    do state_index=1,size(state_ids)
      if(state_ids(state_index)==state)column=state_index
    enddo
    if(column==0)then;bad=1;return;endif
    call position_from_physical_id(identifier,global_size,grid_position)
    index=grid_position-local_origin+1
    if(any(index<1).or.index(1)>size(identifiers,1).or.index(2)>size(identifiers,2).or.&
      index(3)>size(identifiers,3))then;bad=1;return;endif
    if(identifiers(index(1),index(2),index(3))/=identifier)then;bad=1;return;endif
    record(1)=values(index(1),index(2),index(3),column)
    record(2:4)=gradients(index(1),index(2),index(3),column,:)
  end subroutine

  subroutine position_from_physical_id(identifier,global_size,position)
    integer(8),intent(in) :: identifier
    integer,intent(in) :: global_size(3)
    integer,intent(out) :: position(3)
    integer(8) :: zero_based
    zero_based=identifier-1_8
    position(1)=int(modulo(zero_based,int(global_size(1),8)))
    zero_based=zero_based/int(global_size(1),8)
    position(2)=int(modulo(zero_based,int(global_size(2),8)))
    position(3)=int(zero_based/int(global_size(2),8))
  end subroutine

  subroutine validate_dg_dc_buffer_core_faces(faces,fragment_origins,fragment_sizes,global_size,&
      communicator,ok,message)
    type(s_dg_dc_buffer_core_face),intent(inout) :: faces(:)
    integer,intent(in) :: fragment_origins(:,:),fragment_sizes(:,:)
    integer,intent(in) :: global_size(3),communicator
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: iface,i,local_bad,global_bad,seen(6),position(3),expected_core(3)
    integer :: unwrapped(3),depth,tangent(2),local_fragment,neighbor_fragment,plane
    integer(8) :: global_point_count
    local_bad=merge(0,1,size(faces)==6.and.all(global_size>0))
    if(size(fragment_origins,1)/=3.or.any(shape(fragment_origins)/=shape(fragment_sizes)).or.&
      any(fragment_sizes<=0))local_bad=1
    global_point_count=checked_grid_point_count(global_size)
    if(global_point_count<1_8)local_bad=1
    seen=0
    if(local_bad==0)then
      do iface=1,6
        if(faces(iface)%axis<1.or.faces(iface)%axis>3.or.abs(faces(iface)%side)/=1.or.&
           faces(iface)%point_count<0.or.faces(iface)%overlap_depth<=0)local_bad=1
        if(local_bad/=0)exit
        seen(face_slot(faces(iface)%axis,faces(iface)%side))=&
          seen(face_slot(faces(iface)%axis,faces(iface)%side))+1
        if(.not.allocated(faces(iface)%physical_grid_ids))then
          local_bad=1;exit
        endif
        if(.not.allocated(faces(iface)%local_buffer_indices))then
          local_bad=1;exit
        endif
        if(.not.allocated(faces(iface)%neighbor_core_indices))then
          local_bad=1;exit
        endif
        if(size(faces(iface)%physical_grid_ids)/=faces(iface)%point_count.or.&
           any(shape(faces(iface)%local_buffer_indices)/=[3,faces(iface)%point_count]).or.&
           any(shape(faces(iface)%neighbor_core_indices)/=[3,faces(iface)%point_count]))then
          local_bad=1;exit
        endif
        if(any(faces(iface)%local_buffer_indices<1).or.any(faces(iface)%neighbor_core_indices<1))then
          local_bad=1;exit
        endif
        local_fragment=faces(iface)%local_fragment
        neighbor_fragment=faces(iface)%neighbor_fragment
        if(local_fragment<1.or.local_fragment>size(fragment_origins,2).or.&
          neighbor_fragment<1.or.neighbor_fragment>size(fragment_origins,2))then
          local_bad=1;exit
        endif
        tangent=pack([1,2,3],[1,2,3]/=faces(iface)%axis)
        if(faces(iface)%point_count>0)then
          if(modulo(faces(iface)%point_count,faces(iface)%overlap_depth)/=0)then
            local_bad=1;exit
          endif
          plane=faces(iface)%point_count/faces(iface)%overlap_depth
        else
          plane=1
        endif
        do i=1,faces(iface)%point_count
          if(faces(iface)%physical_grid_ids(i)<1_8.or.&
             faces(iface)%physical_grid_ids(i)>global_point_count)then
            local_bad=1;cycle
          endif
          call position_from_physical_id(faces(iface)%physical_grid_ids(i),global_size,position)
          unwrapped(tangent(1))=fragment_origins(tangent(1),local_fragment)+&
            faces(iface)%local_buffer_indices(tangent(1),i)-faces(iface)%overlap_depth-1
          unwrapped(tangent(2))=fragment_origins(tangent(2),local_fragment)+&
            faces(iface)%local_buffer_indices(tangent(2),i)-faces(iface)%overlap_depth-1
          if(faces(iface)%side<0)then
            unwrapped(faces(iface)%axis)=fragment_origins(faces(iface)%axis,local_fragment)-&
              (faces(iface)%overlap_depth+1-faces(iface)%local_buffer_indices(faces(iface)%axis,i))
          else
            unwrapped(faces(iface)%axis)=fragment_origins(faces(iface)%axis,local_fragment)+&
              fragment_sizes(faces(iface)%axis,local_fragment)-1+&
              (faces(iface)%local_buffer_indices(faces(iface)%axis,i)-&
              faces(iface)%overlap_depth-fragment_sizes(faces(iface)%axis,local_fragment))
          endif
          if(any(modulo(unwrapped,global_size)/=position))local_bad=1
          expected_core=position-fragment_origins(:,neighbor_fragment)+1
          if(any(expected_core/=faces(iface)%neighbor_core_indices(:,i)))local_bad=1
          depth=(i-1)/plane+1
          if(faces(iface)%side<0)then
            if(faces(iface)%local_buffer_indices(faces(iface)%axis,i)/=&
              faces(iface)%overlap_depth+1-depth)local_bad=1
          else
            if(faces(iface)%local_buffer_indices(faces(iface)%axis,i)/=&
              faces(iface)%overlap_depth+fragment_sizes(faces(iface)%axis,local_fragment)+depth)local_bad=1
          endif
        enddo
        if(contains_duplicate_ids(faces(iface)%physical_grid_ids))local_bad=1
      enddo
      if(any(seen/=1))local_bad=1
    endif
    call collective_max(local_bad,communicator,global_bad)
    ok=global_bad==0
    do iface=1,size(faces)
      faces(iface)%physical_map_identity=ok
    enddo
    if(ok)then;message='';else;message='DG buffer/core faces: invalid or duplicate physical grid IDs';endif
  end subroutine

  logical function contains_duplicate_ids(identifiers)result(duplicate)
    integer(8),intent(in) :: identifiers(:)
    integer(8),allocatable :: ordered(:)
    integer :: allocation_status,index
    duplicate=.true.
    allocate(ordered(size(identifiers)),stat=allocation_status)
    if(allocation_status/=0)return
    ordered=identifiers
    if(size(ordered)>1)call quicksort_int8(ordered,1,size(ordered))
    duplicate=.false.
    do index=2,size(ordered)
      if(ordered(index)==ordered(index-1))then
        duplicate=.true.;return
      endif
    enddo
  end function

  recursive subroutine quicksort_int8(values,left,right)
    integer(8),intent(inout) :: values(:)
    integer,intent(in) :: left,right
    integer :: lower,upper
    integer(8) :: pivot,temporary
    lower=left;upper=right;pivot=values(left+(right-left)/2)
    do
      do while(values(lower)<pivot);lower=lower+1;enddo
      do while(values(upper)>pivot);upper=upper-1;enddo
      if(lower<=upper)then
        temporary=values(lower);values(lower)=values(upper);values(upper)=temporary
        lower=lower+1;upper=upper-1
      endif
      if(lower>upper)exit
    enddo
    if(left<upper)call quicksort_int8(values,left,upper)
    if(lower<right)call quicksort_int8(values,lower,right)
  end subroutine

  subroutine project_dg_dc_buffer_core_face(face,buffer_values,weights,rank_tolerance,&
      minimum_retained_rank,maximum_projection_residual,maximum_escape_norm,&
      maximum_projection_mismatch,generation,operator_fingerprint,ok,message)
    type(s_dg_dc_buffer_core_face),intent(inout) :: face
    real(8),intent(in) :: buffer_values(:,:),weights(:),rank_tolerance
    integer,intent(in) :: minimum_retained_rank,generation
    real(8),intent(in) :: maximum_projection_residual,maximum_escape_norm,maximum_projection_mismatch
    integer(8),intent(in) :: operator_fingerprint
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: allocation_status
    real(8),allocatable :: reverse_coefficients(:,:),reverse_reconstructed(:,:)
    real(8),allocatable :: candidate_coefficients(:,:),candidate_values(:,:),candidate_normals(:,:)
    type(s_dg_buffer_projector_diagnostics) :: reverse_diagnostics,candidate_diagnostics
    logical :: reverse_ok
    character(256) :: reverse_message
    real(8) :: candidate_mismatch
    ok=.false.;message=''
    if(.not.allocated(face%neighbor_core_values).or..not.allocated(face%neighbor_core_normals))then
      message='DG buffer/core faces: missing exchanged neighbor state window';return
    endif
    if(generation/=face%generation.or..not.face%physical_map_identity.or.&
       generation<=0.or.operator_fingerprint==0_8)then
      message='DG buffer/core faces: stale generation or invalid projected trace context';return
    endif
    if(minimum_retained_rank<1.or.maximum_projection_residual<0d0.or.maximum_escape_norm<0d0.or.&
      maximum_projection_mismatch<0d0)then
      message='DG buffer/core faces: invalid projector acceptance policy';return
    endif
    if(size(buffer_values,2)/=size(face%neighbor_core_values,2).or.&
       size(buffer_values,1)/=face%point_count.or.size(weights)/=face%point_count)then
      message='DG buffer/core faces: invalid local buffer state-window shape';return
    endif
    allocate(candidate_coefficients(size(face%neighbor_core_values,2),size(buffer_values,2)),&
      candidate_values(size(buffer_values,1),size(buffer_values,2)),&
      candidate_normals(size(buffer_values,1),size(buffer_values,2)),&
      reverse_coefficients(size(buffer_values,2),size(face%neighbor_core_values,2)),&
      reverse_reconstructed(size(face%neighbor_core_values,1),size(face%neighbor_core_values,2)),&
      stat=allocation_status)
    if(allocation_status/=0)then;message='DG buffer/core faces: projector allocation failed';return;endif
    call build_dg_buffer_window_projector(face%neighbor_core_values,buffer_values,weights,rank_tolerance,&
      candidate_coefficients,candidate_values,candidate_diagnostics,ok,message)
    if(.not.ok)return
    call build_dg_buffer_window_projector(buffer_values,face%neighbor_core_values,weights,rank_tolerance,&
      reverse_coefficients,reverse_reconstructed,reverse_diagnostics,reverse_ok,reverse_message)
    if(.not.reverse_ok)then
      ok=.false.;message='DG buffer/core faces: reverse projector failed: '//trim(reverse_message);return
    endif
    candidate_mismatch=abs(candidate_diagnostics%projection_residual-&
      reverse_diagnostics%projection_residual)
    if(candidate_diagnostics%retained_rank<minimum_retained_rank.or.&
      reverse_diagnostics%retained_rank<minimum_retained_rank.or.&
      max(candidate_diagnostics%projection_residual,reverse_diagnostics%projection_residual)>&
      maximum_projection_residual.or.&
      max(candidate_diagnostics%escape_norm,reverse_diagnostics%escape_norm)>maximum_escape_norm.or.&
      candidate_mismatch>maximum_projection_mismatch)then
      ok=.false.;message='DG buffer/core faces: projector diagnostics rejected by acceptance policy';return
    endif
    candidate_normals=matmul(face%neighbor_core_normals,candidate_coefficients)
    if(.not.all(ieee_is_finite(candidate_normals)))then
      ok=.false.;message='DG buffer/core faces: nonfinite reconstructed normal trace';return
    endif
    call move_alloc(candidate_coefficients,face%coefficients)
    call move_alloc(candidate_values,face%reconstructed_values)
    call move_alloc(candidate_normals,face%reconstructed_normals)
    face%diagnostics=candidate_diagnostics
    face%forward_reverse_projection_mismatch=candidate_mismatch
    face%operator_fingerprint=operator_fingerprint
  end subroutine

  subroutine validate_partition(origins,sizes,global_size,bad)
    integer,intent(in) :: origins(:,:),sizes(:,:),global_size(3)
    integer,intent(out) :: bad
    integer :: first,second
    integer(8) :: covered_volume,global_volume
    bad=0
    covered_volume=0_8
    global_volume=checked_grid_point_count(global_size)
    if(global_volume<1_8)bad=1
    do first=1,size(origins,2)
      if(any(origins(:,first)<0).or..not.box_within_box(origins(:,first),sizes(:,first),&
        [0,0,0],global_size))bad=1
      if(checked_grid_point_count(sizes(:,first))<1_8.or.&
        covered_volume>huge(covered_volume)-checked_grid_point_count(sizes(:,first)))then
        bad=1
      else
        covered_volume=covered_volume+checked_grid_point_count(sizes(:,first))
      endif
      do second=first+1,size(origins,2)
        if(boxes_overlap(origins(:,first),sizes(:,first),origins(:,second),sizes(:,second)))bad=1
      enddo
    enddo
    if(covered_volume/=global_volume)bad=1
  end subroutine

  integer function find_owner(position,origins,sizes)result(owner)
    integer,intent(in) :: position(3),origins(:,:),sizes(:,:)
    integer :: fragment,count
    owner=0;count=0
    do fragment=1,size(origins,2)
      if(point_in_box(position,origins(:,fragment),sizes(:,fragment)))then
        owner=fragment;count=count+1
      endif
    enddo
    if(count/=1)owner=0
  end function

  integer(8) function physical_id(position,global_size)result(identifier)
    integer,intent(in) :: position(3),global_size(3)
    identifier=int(position(1),8)+int(global_size(1),8)*&
      (int(position(2),8)+int(global_size(2),8)*int(position(3),8))+1_8
  end function

  integer(8) function checked_grid_point_count(global_size)result(count)
    integer,intent(in) :: global_size(3)
    integer(8) :: partial,limit
    count=-1_8
    if(any(global_size<=0))return
    limit=huge(count)
    if(int(global_size(1),8)>limit/int(global_size(2),8))return
    partial=int(global_size(1),8)*int(global_size(2),8)
    if(partial>limit/int(global_size(3),8))return
    count=partial*int(global_size(3),8)
  end function

  logical function point_in_box(position,origin,extent)result(inside)
    integer,intent(in) :: position(3),origin(3),extent(3)
    inside=all(int(position,8)>=int(origin,8)).and.&
      all(int(position,8)<int(origin,8)+int(extent,8))
  end function

  logical function box_within_box(inner_origin,inner_extent,outer_origin,outer_extent)result(inside)
    integer,intent(in) :: inner_origin(3),inner_extent(3),outer_origin(3),outer_extent(3)
    inside=all(int(inner_origin,8)>=int(outer_origin,8)).and.&
      all(int(inner_origin,8)+int(inner_extent,8)<=int(outer_origin,8)+int(outer_extent,8))
  end function

  logical function boxes_overlap(first_origin,first_extent,second_origin,second_extent)result(overlap)
    integer,intent(in) :: first_origin(3),first_extent(3),second_origin(3),second_extent(3)
    overlap=all(int(first_origin,8)<int(second_origin,8)+int(second_extent,8)).and.&
      all(int(second_origin,8)<int(first_origin,8)+int(first_extent,8))
  end function

  integer function face_slot(axis,side)result(slot)
    integer,intent(in) :: axis,side
    slot=2*(axis-1)+merge(1,2,side<0)
  end function

  subroutine collective_max(local_value,communicator,global_value)
    integer,intent(in) :: local_value,communicator
    integer,intent(out) :: global_value
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS)global_value=1
#else
    global_value=local_value
    if(communicator<0)global_value=1
#endif
  end subroutine
end module
