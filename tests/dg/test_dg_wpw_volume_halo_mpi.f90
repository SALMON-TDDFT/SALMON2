program test_dg_wpw_volume_halo_mpi
  use mpi
  use rt_dg_wpw_volume_halo_provider, only: s_dg_wpw_volume_halo_state, &
    exchange_dg_wpw_volume_halo, read_dg_wpw_volume_halo
  implicit none
  type(s_dg_wpw_volume_halo_state) :: halo
  integer :: ierr,rank,peer,info,send_w_ids(1),grid(3)
  integer :: box_lo(3),box_hi(3)
  complex(8) :: send_values(1,2),send_gradients(3,1,2),value,gradient(3)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,info,ierr)
  if(ierr/=MPI_SUCCESS.or.info/=2)error stop 1
  peer=1-rank
  send_w_ids=[rank+1]
  box_lo=[1,1,1];box_hi=[2,1,1]
  send_values(1,1)=cmplx(10*rank+1,rank,8)
  send_values(1,2)=cmplx(10*rank+2,-rank,8)
  send_gradients=(0d0,0d0)
  send_gradients(1,1,1)=cmplx(rank+0.25d0,0.5d0,8)
  send_gradients(1,1,2)=cmplx(rank+0.5d0,-0.25d0,8)

  call exchange_dg_wpw_volume_halo(MPI_COMM_WORLD,7,peer,send_w_ids,box_lo,box_hi, &
    send_values,send_gradients,halo,info)
  if(info/=0.or..not.halo%valid.or.halo%epoch/=7)error stop 2
  grid=[2,1,1]
  call read_dg_wpw_volume_halo(halo,peer+1,grid,7,value,gradient,info)
  if(info/=0)error stop 3
  if(abs(value-cmplx(10*peer+2,-peer,8))>1d-13)error stop 4
  if(maxval(abs(gradient-[cmplx(peer+0.5d0,-0.25d0,8),(0d0,0d0),(0d0,0d0)]))>1d-13)error stop 5
  call read_dg_wpw_volume_halo(halo,peer+1,grid,6,value,gradient,info)
  if(info==0)error stop 6
  if(rank==0)print '(a)','PASS two-rank WPW volume-support halo'
  call MPI_Finalize(ierr)
end program test_dg_wpw_volume_halo_mpi
