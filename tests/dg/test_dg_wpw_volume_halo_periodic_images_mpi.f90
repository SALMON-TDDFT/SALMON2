program test_dg_wpw_volume_halo_periodic_images_mpi
  use mpi
  use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_send,s_dg_wpw_volume_halo_state,&
    exchange_dg_wpw_volume_halo_schedule,read_dg_wpw_volume_halo
  implicit none
  type(s_dg_wpw_volume_halo_send)::send(2)
  type(s_dg_wpw_volume_halo_state),allocatable::halo(:)
  integer::ierr,rank,nrank,peer,i,info
  complex(8)::value,gradient(3)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(ierr/=MPI_SUCCESS.or.nrank/=2)error stop 1
  peer=1-rank
  send(1)%peer=peer;send(1)%image_shift=[-1,0,0];send(1)%box_lo=[1,1,1];send(1)%box_hi=[1,1,1]
  send(2)%peer=peer;send(2)%image_shift=[ 1,0,0];send(2)%box_lo=[4,1,1];send(2)%box_hi=[4,1,1]
  do i=1,2
    allocate(send(i)%w_ids(1),send(i)%values(1,1),send(i)%gradients(3,1,1))
    send(i)%w_ids=[rank+1];send(i)%values(1,1)=cmplx(100*rank+i,0d0,8);send(i)%gradients=(0d0,0d0)
  enddo
  call exchange_dg_wpw_volume_halo_schedule(MPI_COMM_WORLD,21,send,halo,info)
  if(info/=0.or.size(halo)/=2.or.any(.not.halo%valid))error stop 2
  ! Local minus image receives the peer's plus-image payload, and vice versa.
  call read_dg_wpw_volume_halo(halo(1),peer+1,[4,1,1],21,value,gradient,info)
  if(info/=0.or.abs(value-cmplx(100*peer+2,0d0,8))>1d-13)error stop 3
  call read_dg_wpw_volume_halo(halo(2),peer+1,[1,1,1],21,value,gradient,info)
  if(info/=0.or.abs(value-cmplx(100*peer+1,0d0,8))>1d-13)error stop 4
  if(rank==0)print '(a)','PASS distinct two-fragment periodic volume-halo images'
  call MPI_Finalize(ierr)
end program test_dg_wpw_volume_halo_periodic_images_mpi
