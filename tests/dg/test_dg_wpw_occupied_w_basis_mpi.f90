program test_dg_wpw_occupied_w_basis_mpi
  use mpi
  use dg_wpw_occupied_w_basis,only:s_dg_wpw_occupied_w_basis,&
    initialize_dg_wpw_occupied_w_basis_collective
  implicit none
  type(s_dg_wpw_occupied_w_basis)::basis
  integer::rank,nrank,ierr,info,local_keys(5,1),grid_ids(2)
  integer(8)::saved_fingerprint
  complex(8)::core_values(2,1),buffer_values(2,1),buffer_gradients(3,2,1)

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  if(rank==0)then;local_keys(:,1)=[3,4,0,1,0]
  else;local_keys(:,1)=[1,2,0,0,0];endif
  grid_ids=[10*rank+1,10*rank+2];core_values=cmplx(rank+1d0,0d0,8)
  buffer_values=cmplx(rank+0.5d0,0d0,8);buffer_gradients=(0.1d0,0d0)
  call initialize_dg_wpw_occupied_w_basis_collective(basis,MPI_COMM_WORLD,rank+1,local_keys,&
    grid_ids,core_values,[1,1,1],[2,1,1],buffer_values,buffer_gradients,2d0,4,2,info)
  if(info/=0.or..not.basis%valid.or.basis%global_count/=2.or.basis%local_count/=1)&
    call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  if(rank==0.and.any(basis%owned_ids/=[2]))call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  if(rank==1.and.any(basis%owned_ids/=[1]))call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  saved_fingerprint=basis%fingerprint
  call initialize_dg_wpw_occupied_w_basis_collective(basis,MPI_COMM_WORLD,rank+1,local_keys,&
    grid_ids,core_values,[1,1,1],[2,1,1],buffer_values,buffer_gradients,2d0,5,3,info)
  if(info==0.or.basis%epoch/=4.or.basis%fingerprint/=saved_fingerprint)&
    call MPI_Abort(MPI_COMM_WORLD,5,ierr)
  if(rank==0)print '(a)','PASS collective occupied-W stable ID descriptor'
  call MPI_Finalize(ierr)
end program test_dg_wpw_occupied_w_basis_mpi
