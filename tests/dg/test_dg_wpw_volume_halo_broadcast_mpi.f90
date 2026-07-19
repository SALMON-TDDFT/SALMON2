program test_dg_wpw_volume_halo_broadcast_mpi
  use mpi
  use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_state,&
    broadcast_dg_wpw_volume_halos,read_dg_wpw_volume_halo,dg_wpw_volume_payload_checksum
  implicit none
  type(s_dg_wpw_volume_halo_state),allocatable::halos(:)
  complex(8)::value,gradient(3)
  integer::ierr,rank,nrank,info
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2)error stop 1
  if(rank==0)then
    allocate(halos(1));halos(1)%epoch=4;halos(1)%source_rank=7;halos(1)%image_shift=[-1,0,0]
    halos(1)%box_lo=[2,1,1];halos(1)%box_hi=[2,1,1];halos(1)%valid=.true.
    allocate(halos(1)%w_ids(1),halos(1)%values(1,1),halos(1)%gradients(3,1,1))
    halos(1)%w_ids=[8];halos(1)%values(1,1)=(3d0,0.5d0);halos(1)%gradients(:,1,1)=(2d0,0d0)
    halos(1)%checksum=dg_wpw_volume_payload_checksum(halos(1)%w_ids,halos(1)%box_lo,halos(1)%box_hi,&
      halos(1)%values,halos(1)%gradients)
  else
    allocate(halos(0))
  endif
  call broadcast_dg_wpw_volume_halos(MPI_COMM_WORLD,0,halos,info)
  if(info/=0.or.size(halos)/=1.or..not.halos(1)%valid.or.halos(1)%source_rank/=7.or.&
    any(halos(1)%image_shift/=[-1,0,0]))error stop 2
  call read_dg_wpw_volume_halo(halos(1),8,[2,1,1],4,value,gradient,info)
  if(info/=0.or.abs(value-(3d0,0.5d0))>1d-13.or.any(abs(gradient-(2d0,0d0))>1d-13))error stop 3
  if(rank==0)halos(1)%checksum=halos(1)%checksum+1
  call broadcast_dg_wpw_volume_halos(MPI_COMM_WORLD,0,halos,info)
  if(info==0.or.any([(halos(ierr)%valid,ierr=1,size(halos))]))error stop 4
  if(rank==0)print '(a)','PASS fragment-rank WPW volume-halo broadcast'
  call MPI_Finalize(ierr)
end program
