program test_dg_dc_empty_state_surrogate_mpi
  use mpi
  use dg_dc_buffer_core_faces,only:s_dg_dc_buffer_core_face,build_dg_dc_buffer_core_faces,&
    exchange_dg_dc_buffer_core_state_window,assemble_dg_dc_local_buffer_state_window
  implicit none
  type(s_dg_dc_buffer_core_face),allocatable::faces(:)
  integer::ierr,rank,iface,point,state
  integer::origins(3,1),sizes(3,1),global_size(3),local_origin(3),local_size(3)
  integer,allocatable::state_ids(:)
  integer(8),allocatable::physical_ids(:,:,:)
  real(8),allocatable::values(:,:,:,:),gradients(:,:,:,:,:),owned_faces(:,:,:)
  logical::ok
  character(256)::message
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  origins(:,1)=[0,0,0];sizes(:,1)=[2,2,1];global_size=sizes(:,1)
  local_origin=origins(:,1);local_size=sizes(:,1)
  if(rank==0)then
    allocate(state_ids(2),source=[1,2])
  else
    allocate(state_ids(0))
  endif
  allocate(physical_ids(2,2,1),values(2,2,1,size(state_ids)),&
    gradients(2,2,1,size(state_ids),3))
  physical_ids(:,:,1)=reshape([1_8,2_8,3_8,4_8],[2,2])
  do state=1,size(state_ids)
    values(:,:,:,state)=dble(state_ids(state))
    gradients(:,:,:,state,:)=dble(state_ids(state))
  enddo
  call build_dg_dc_buffer_core_faces(1,origins,sizes,global_size,1,29,MPI_COMM_WORLD,&
    faces,ok,message,local_origin,local_size)
  if(.not.ok)error stop trim(message)
  call exchange_dg_dc_buffer_core_state_window(faces,origins,sizes,global_size,local_origin,local_size,&
    physical_ids,values,gradients,state_ids,2,MPI_COMM_WORLD,ok,message)
  if(.not.ok)error stop trim(message)
  allocate(owned_faces(maxval([(faces(iface)%point_count,iface=1,6)]),size(state_ids),6),source=0d0)
  do iface=1,6;do state=1,size(state_ids);do point=1,faces(iface)%point_count
    owned_faces(point,state,iface)=dble(state_ids(state))
  enddo;enddo;enddo
  call assemble_dg_dc_local_buffer_state_window(faces,owned_faces,state_ids,2,MPI_COMM_WORLD,ok,message)
  if(.not.ok)error stop trim(message)
  do iface=1,6
    if(faces(iface)%point_count>0.and.size(faces(iface)%local_buffer_values,2)/=2)&
      error stop 'empty-state surrogate did not receive the complete state window'
  enddo
  if(rank==0)print '(a)','PASS empty-state duplicate-core surrogate'
  call MPI_Finalize(ierr)
end program test_dg_dc_empty_state_surrogate_mpi
