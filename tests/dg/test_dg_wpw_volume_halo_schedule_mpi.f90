program test_dg_wpw_volume_halo_schedule_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use rt_dg_wpw_volume_halo_provider, only: s_dg_wpw_volume_halo_send, &
    s_dg_wpw_volume_halo_state, exchange_dg_wpw_volume_halo_schedule, read_dg_wpw_volume_halo
  implicit none
  type(s_dg_wpw_volume_halo_send),allocatable::send(:)
  type(s_dg_wpw_volume_halo_state),allocatable::halo(:)
  integer::ierr,rank,nrank,info,i,peer,grid(3)
  complex(8)::value,gradient(3)

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(ierr/=MPI_SUCCESS.or.nrank/=4)error stop 1
  allocate(send(2))
  ! Every rank lists its clockwise neighbor first.  A sequential blocking
  ! exchange in this order deadlocks; the schedule API must post all receives.
  send(1)%peer=modulo(rank+1,nrank)
  send(2)%peer=modulo(rank-1,nrank)
  do i=1,2
    peer=send(i)%peer
    send(i)%box_lo=[1,1,1];send(i)%box_hi=[1,1,1]
    allocate(send(i)%w_ids(1),send(i)%values(1,1),send(i)%gradients(3,1,1))
    send(i)%w_ids=[rank+1]
    send(i)%values(1,1)=cmplx(100*rank+peer,rank,8)
    send(i)%gradients(:,1,1)=[cmplx(rank,peer,8),(0d0,0d0),(0d0,0d0)]
  enddo
  call exchange_dg_wpw_volume_halo_schedule(MPI_COMM_WORLD,11,send,halo,info)
  if(info/=0.or.size(halo)/=2.or.any(.not.halo%valid))error stop 2
  grid=[1,1,1]
  do i=1,2
    peer=send(i)%peer
    call read_dg_wpw_volume_halo(halo(i),peer+1,grid,11,value,gradient,info)
    if(info/=0)error stop 3
    if(abs(value-cmplx(100*peer+rank,peer,8))>1d-13)error stop 4
  enddo
  if(rank==0)print '(a)','PASS multi-neighbor nonblocking WPW volume halo schedule'

  if(rank==0)send(1)%values(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call exchange_dg_wpw_volume_halo_schedule(MPI_COMM_WORLD,12,send,halo,info)
  if(info==0.or.size(halo)/=0)error stop 5
  call MPI_Finalize(ierr)
end program test_dg_wpw_volume_halo_schedule_mpi
