program test_dg_dc_buffer_core_faces_mpi
  use mpi
  use dg_dc_buffer_core_faces,only:s_dg_dc_buffer_core_face,&
    build_dg_dc_buffer_core_faces,exchange_dg_dc_buffer_core_state_window,&
    validate_dg_dc_buffer_core_faces,assemble_dg_dc_local_buffer_state_window,&
    project_dg_dc_buffer_core_face
  implicit none
  type(s_dg_dc_buffer_core_face),allocatable :: faces(:)
  integer,allocatable :: origins(:,:),sizes(:,:)
  real(8),allocatable :: core_values(:,:,:,:),core_gradients(:,:,:,:,:)
  integer(8),allocatable :: core_physical_grid_ids(:,:,:)
  real(8),allocatable :: buffer_values(:,:),buffer_normals(:,:),weights(:)
  real(8),allocatable :: owned_buffer_face_values(:,:,:)
  integer,allocatable :: owned_state_ids(:)
  integer :: ierr,rank,nrank,fragment,buffer_width,configured_states,nfragment,lane,spatial_part,state_part
  integer :: test_generation
  integer :: global_size(3),fragment_grid(3),coords(3),local_size(3),local_origin(3)
  integer :: axis,side,iface,point,state,expected_fragment
  integer(8) :: saved_id,saved_id2
  integer :: saved_core1(3),saved_core2(3),face_plane
  logical :: ok,face_condition
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank==1)then
    fragment_grid=[1,1,1]
    global_size=[4,3,2]
    nfragment=1
    fragment=1
    spatial_part=0
    state_part=0
  else if(nrank==2)then
    fragment_grid=[2,1,1]
    global_size=[8,3,2]
    nfragment=2
    fragment=rank+1
    spatial_part=0
    state_part=0
  else if(nrank==8)then
    fragment_grid=[2,1,1]
    global_size=[8,4,2]
    nfragment=2
    fragment=modulo(rank,nfragment)+1
    lane=rank/nfragment
    spatial_part=modulo(lane,2)
    state_part=lane/2
  else
    error stop 'fixture requires one, two, or eight MPI ranks'
  endif
  allocate(origins(3,nfragment),sizes(3,nfragment))
  call build_uniform_partition(fragment_grid,global_size,origins,sizes)
  configured_states=3
  local_size=sizes(:,fragment)
  local_origin=origins(:,fragment)
  if(nrank==8)then
    local_size(2)=local_size(2)/2
    local_origin(2)=local_origin(2)+spatial_part*local_size(2)
    if(state_part==0)then
      allocate(owned_state_ids(2),source=[1,2])
    else
      allocate(owned_state_ids(1),source=[3])
    endif
  else
    allocate(owned_state_ids(configured_states))
    owned_state_ids=[(state,state=1,configured_states)]
  endif
  allocate(core_values(local_size(1),local_size(2),local_size(3),size(owned_state_ids)),source=0d0)
  allocate(core_gradients(local_size(1),local_size(2),local_size(3),size(owned_state_ids),3))
  allocate(core_physical_grid_ids(local_size(1),local_size(2),local_size(3)))
  call fill_core_ids(local_origin,local_size,global_size,core_physical_grid_ids)
  call fill_core_values(fragment,local_origin,local_size,global_size,core_values,owned_state_ids)
  core_gradients(:,:,:,:,1)=2d0*core_values
  core_gradients(:,:,:,:,2)=3d0*core_values
  core_gradients(:,:,:,:,3)=5d0*core_values

  do buffer_width=1,2
    call build_dg_dc_buffer_core_faces(fragment,origins,sizes,global_size,buffer_width,17, &
      MPI_COMM_WORLD,faces,ok,message,local_origin,local_size)
    call require_collective(ok,'canonical six-face construction: '//trim(message))
    call require_collective(size(faces)==6,'exactly six signed faces')
    do iface=1,6
      call require_collective(faces(iface)%overlap_depth==buffer_width,'buffer overlap depth')
      call require_collective(faces(iface)%generation==17,'mapping generation')
      call require_collective(faces(iface)%point_count==expected_face_points(faces(iface)%axis,&
        faces(iface)%side,buffer_width),'tangential point count')
    enddo
    call validate_dg_dc_buffer_core_faces(faces,origins,sizes,global_size,MPI_COMM_WORLD,ok,message)
    call require_collective(ok,'physical grid-map identity')
    call exchange_dg_dc_buffer_core_state_window(faces,origins,sizes,global_size,local_origin,local_size,&
      core_physical_grid_ids,core_values,core_gradients,owned_state_ids,configured_states, &
      MPI_COMM_WORLD,ok,message)
    call require_collective(ok,'complete configured-state exchange: '//trim(message))
    do iface=1,6
      do state=1,configured_states
        do point=1,faces(iface)%point_count
          if(abs(faces(iface)%neighbor_core_values(point,state)-&
            value_from_id(faces(iface)%physical_grid_ids(point),state))>1d-12) &
            call fail('physical buffer/core value mismatch')
          if(abs(faces(iface)%neighbor_core_normals(point,state)+dble(faces(iface)%side)*&
            gradient_factor(faces(iface)%axis)*value_from_id(faces(iface)%physical_grid_ids(point),state))&
            >1d-12)call fail('axis- and sign-dependent normal mismatch')
        enddo
      enddo
      do point=1,faces(iface)%point_count
        call require_face_point_identity(faces(iface),point,buffer_width)
      enddo
    enddo
    allocate(owned_buffer_face_values(maxval([(faces(iface)%point_count,iface=1,6)]),&
      size(owned_state_ids),6),source=0d0)
    do iface=1,6;do state=1,size(owned_state_ids);do point=1,faces(iface)%point_count
      owned_buffer_face_values(point,state,iface)=&
        value_from_id(faces(iface)%physical_grid_ids(point),owned_state_ids(state))
    enddo;enddo;enddo
    call assemble_dg_dc_local_buffer_state_window(faces,owned_buffer_face_values,owned_state_ids,&
      configured_states,MPI_COMM_WORLD,ok,message)
    call require_collective(ok,'disjoint-state local-buffer window assembly: '//trim(message))
    do iface=1,6
      ok=.true.;face_condition=.true.
      if(faces(iface)%point_count>0)then
        allocate(buffer_values(faces(iface)%point_count,configured_states),&
          buffer_normals(faces(iface)%point_count,configured_states),weights(faces(iface)%point_count))
        buffer_values=faces(iface)%local_buffer_values
        buffer_normals=faces(iface)%neighbor_core_normals
        weights=1d0
        call project_dg_dc_buffer_core_face(faces(iface),buffer_values,weights,1d-12,1,&
          1d-9,1d-9,1d-9,17,&
          991_8,ok,message)
      endif
      call require_collective(ok,'all-face value/normal projection: '//trim(message))
      if(faces(iface)%point_count>0)then
        face_condition=maxval(abs(faces(iface)%reconstructed_values-buffer_values))<1d-10.and.&
          maxval(abs(faces(iface)%reconstructed_normals-buffer_normals))<1d-10.and.&
          faces(iface)%diagnostics%configured_states==configured_states.and.&
          faces(iface)%diagnostics%retained_rank>0.and.faces(iface)%operator_fingerprint==991_8.and.&
          faces(iface)%physical_map_identity.and.&
          faces(iface)%forward_reverse_projection_mismatch<1d-10
        deallocate(buffer_values,buffer_normals,weights)
      endif
      call require_collective(face_condition,'all-face projected traces and diagnostics')
    enddo
    allocate(buffer_values(faces(1)%point_count,configured_states),weights(faces(1)%point_count))
    buffer_values=faces(1)%neighbor_core_values;weights=1d0
    buffer_values(1,1)=buffer_values(1,1)+0.1d0
    call project_dg_dc_buffer_core_face(faces(1),buffer_values,weights,1d-12,1,&
      1d-14,1d-14,1d-14,17,991_8,ok,message)
    call require_collective(.not.ok,'asymmetric projection diagnostics rejection')
    buffer_values=faces(1)%neighbor_core_values
    call project_dg_dc_buffer_core_face(faces(1),buffer_values,weights,1d-12,1,&
      1d-9,1d-9,1d-9,18,991_8,ok,message)
    call require_collective(.not.ok,'stale projector generation rejection')
    deallocate(buffer_values,weights)
    saved_id=faces(1)%physical_grid_ids(1)
    faces(1)%physical_grid_ids(1)=faces(1)%physical_grid_ids(2)
    call validate_dg_dc_buffer_core_faces(faces,origins,sizes,global_size,MPI_COMM_WORLD,ok,message)
    call require_collective(.not.ok,'duplicate physical grid ID rejection')
    faces(1)%physical_grid_ids(1)=int(product(global_size),8)+1_8
    call validate_dg_dc_buffer_core_faces(faces,origins,sizes,global_size,MPI_COMM_WORLD,ok,message)
    call require_collective(.not.ok,'out-of-range physical grid ID rejection')
    if(buffer_width==2)then
      face_plane=faces(1)%point_count/buffer_width
      saved_id2=faces(1)%physical_grid_ids(face_plane+1)
      saved_core1=faces(1)%neighbor_core_indices(:,1)
      saved_core2=faces(1)%neighbor_core_indices(:,face_plane+1)
      faces(1)%physical_grid_ids(1)=saved_id2
      faces(1)%neighbor_core_indices(:,1)=saved_core2
      faces(1)%physical_grid_ids(face_plane+1)=saved_id
      faces(1)%neighbor_core_indices(:,face_plane+1)=saved_core1
      call validate_dg_dc_buffer_core_faces(faces,origins,sizes,global_size,MPI_COMM_WORLD,ok,message)
      call require_collective(.not.ok,'normal-depth physical image mismatch rejection')
      faces(1)%physical_grid_ids(face_plane+1)=saved_id2
      faces(1)%neighbor_core_indices(:,face_plane+1)=saved_core2
      faces(1)%neighbor_core_indices(:,1)=saved_core1
    endif
    faces(1)%physical_grid_ids(1)=saved_id
    faces(1)%local_buffer_indices(1,1)=faces(1)%local_buffer_indices(1,1)+1
    call validate_dg_dc_buffer_core_faces(faces,origins,sizes,global_size,MPI_COMM_WORLD,ok,message)
    call require_collective(.not.ok,'buffer-index/physical-map mismatch rejection')
    faces(1)%local_buffer_indices(1,1)=faces(1)%local_buffer_indices(1,1)-1
    deallocate(owned_buffer_face_values)
    deallocate(faces)
  enddo

  call check_fragment_relabeling(origins,sizes,global_size,core_values,configured_states)

  call build_dg_dc_buffer_core_faces(fragment,origins,sizes,global_size,1,29, &
    MPI_COMM_WORLD,faces,ok,message,local_origin,local_size)
  if(rank==0)configured_states=4
  call exchange_dg_dc_buffer_core_state_window(faces,origins,sizes,global_size,local_origin,local_size,&
    core_physical_grid_ids,core_values,core_gradients,owned_state_ids,configured_states, &
    MPI_COMM_WORLD,ok,message)
  call require_collective(.not.ok,'rank disagreement rejected before communication')
  configured_states=3
  test_generation=41
  if(rank==0)test_generation=42
  call build_dg_dc_buffer_core_faces(fragment,origins,sizes,global_size,1,test_generation,&
    MPI_COMM_WORLD,faces,ok,message,local_origin,local_size)
  if(nrank>1)call require_collective(.not.ok,'mapping generation disagreement rejected before communication')

  if(rank==0)print '(a)','PASS buffer-to-neighbor-core physical face mapping'
  call MPI_Finalize(ierr)
contains
  real(8) function gradient_factor(face_axis)result(factor)
    integer,intent(in) :: face_axis
    select case(face_axis)
    case(1);factor=2d0
    case(2);factor=3d0
    case default;factor=5d0
    end select
  end function

  subroutine require_face_point_identity(face,face_point,depth)
    type(s_dg_dc_buffer_core_face),intent(in) :: face
    integer,intent(in) :: face_point,depth
    integer :: position(3),expected_core(3),expected_buffer_axis,depth_index,plane
    call position_from_id(face%physical_grid_ids(face_point),global_size,position)
    expected_core=position-origins(:,face%neighbor_fragment)+1
    if(any(face%neighbor_core_indices(:,face_point)/=expected_core))&
      call fail('neighbor core index/physical ID identity')
    plane=face%point_count/depth
    depth_index=(face_point-1)/plane+1
    expected_buffer_axis=merge(depth+1-depth_index,&
      depth+sizes(face%axis,fragment)+depth_index,face%side<0)
    if(face%local_buffer_indices(face%axis,face_point)/=expected_buffer_axis)&
      call fail('signed local buffer storage index')
  end subroutine

  subroutine position_from_id(identifier,global,position)
    integer(8),intent(in) :: identifier
    integer,intent(in) :: global(3)
    integer,intent(out) :: position(3)
    integer(8) :: zero_based
    zero_based=identifier-1_8
    position(1)=int(modulo(zero_based,int(global(1),8)))
    zero_based=zero_based/int(global(1),8)
    position(2)=int(modulo(zero_based,int(global(2),8)))
    position(3)=int(zero_based/int(global(2),8))
  end subroutine

  subroutine fill_core_ids(origin,extent,global,identifiers)
    integer,intent(in) :: origin(3),extent(3),global(3)
    integer(8),intent(out) :: identifiers(:,:,:)
    integer :: x,y,z
    do z=1,extent(3);do y=1,extent(2);do x=1,extent(1)
      identifiers(x,y,z)=grid_id(origin+[x-1,y-1,z-1],global)
    enddo;enddo;enddo
  end subroutine

  integer function expected_face_points(face_axis,face_side,depth)result(count)
    integer,intent(in) :: face_axis,face_side,depth
    integer :: tangential(2)
    tangential=pack([1,2,3],[1,2,3]/=face_axis)
    count=0
    if((face_side<0.and.local_origin(face_axis)==origins(face_axis,fragment)).or.&
       (face_side>0.and.local_origin(face_axis)+local_size(face_axis)==&
        origins(face_axis,fragment)+sizes(face_axis,fragment)))&
      count=depth*local_size(tangential(1))*local_size(tangential(2))
  end function

  subroutine build_uniform_partition(grid,global,partition_origins,partition_sizes)
    integer,intent(in) :: grid(3),global(3)
    integer,intent(out) :: partition_origins(:,:),partition_sizes(:,:)
    integer :: f,x,y,z
    f=0
    do x=0,grid(1)-1;do y=0,grid(2)-1;do z=0,grid(3)-1
      f=f+1
      partition_sizes(:,f)=global/grid
      partition_origins(:,f)=[x,y,z]*partition_sizes(:,f)
    enddo;enddo;enddo
  end subroutine

  subroutine fill_core_values(owner,origin,extent,global,values,state_ids)
    integer,intent(in) :: owner,origin(3),extent(3),global(3)
    integer,intent(in) :: state_ids(:)
    real(8),intent(out) :: values(:,:,:,:)
    integer :: x,y,z,s
    do s=1,size(values,4);do z=1,extent(3);do y=1,extent(2);do x=1,extent(1)
      values(x,y,z,s)=value_from_id(grid_id(origin+[x-1,y-1,z-1],global),state_ids(s))
    enddo;enddo;enddo;enddo
  end subroutine

  integer(8) function grid_id(position,global)result(identifier)
    integer,intent(in) :: position(3),global(3)
    integer :: wrapped(3)
    wrapped=modulo(position,global)
    identifier=int(wrapped(1)+global(1)*(wrapped(2)+global(2)*wrapped(3))+1,8)
  end function

  real(8) function value_from_id(identifier,state)result(value)
    integer(8),intent(in) :: identifier
    integer,intent(in) :: state
    select case(state)
    case(1);value=1d0
    case(2);value=dble(identifier)/50d0
    case(3);value=(dble(identifier)/50d0)**2
    case default;value=(dble(identifier)/50d0)**(state-1)
    end select
  end function

  subroutine check_fragment_relabeling(base_origins,base_sizes,global,values,nstate)
    integer,intent(in) :: base_origins(:,:),base_sizes(:,:),global(3),nstate
    real(8),intent(in) :: values(:,:,:,:)
    integer,allocatable :: relabeled_origins(:,:),relabeled_sizes(:,:),permutation(:)
    type(s_dg_dc_buffer_core_face),allocatable :: relabeled_faces(:)
    integer :: i,new_fragment
    allocate(relabeled_origins(3,nfragment),relabeled_sizes(3,nfragment),permutation(nfragment))
    do i=1,nfragment
      permutation(i)=nfragment-i+1
      relabeled_origins(:,i)=base_origins(:,permutation(i))
      relabeled_sizes(:,i)=base_sizes(:,permutation(i))
    enddo
    new_fragment=nfragment-fragment+1
    call build_dg_dc_buffer_core_faces(new_fragment,relabeled_origins,relabeled_sizes,global,1,23, &
      MPI_COMM_WORLD,relabeled_faces,ok,message)
    call require_collective(ok,'fragment relabeling')
    do i=1,6
      call require_collective(all(relabeled_faces(i)%physical_grid_ids==&
        faces_for_original(i,base_origins,base_sizes,global)),'relabeling physical identity')
    enddo
  end subroutine

  function faces_for_original(face_index,base_origins,base_sizes,global)result(ids)
    integer,intent(in) :: face_index,base_origins(:,:),base_sizes(:,:),global(3)
    integer(8),allocatable :: ids(:)
    type(s_dg_dc_buffer_core_face),allocatable :: original_faces(:)
    call build_dg_dc_buffer_core_faces(fragment,base_origins,base_sizes,global,1,23, &
      MPI_COMM_WORLD,original_faces,ok,message)
    ids=original_faces(face_index)%physical_grid_ids
  end function

  subroutine require_collective(condition,label)
    logical,intent(in) :: condition
    character(*),intent(in) :: label
    logical :: all_condition
    call MPI_Allreduce(condition,all_condition,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(.not.all_condition)call fail(label)
  end subroutine

  subroutine fail(label)
    character(*),intent(in) :: label
    if(rank==0)write(*,'(a,1x,a)')'FAIL',trim(label)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program
