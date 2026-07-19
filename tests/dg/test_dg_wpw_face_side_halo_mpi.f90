program test_dg_wpw_face_side_halo_mpi
  use mpi
  use rt_dg_wpw_face_side_halo,only:s_dg_wpw_face_side_send,s_dg_wpw_face_side_state,&
    exchange_dg_wpw_face_side_schedule
  implicit none
  type(s_dg_wpw_face_side_send)::send(2)
  type(s_dg_wpw_face_side_state),allocatable::remote(:)
  integer::ierr,rank,nrank,info,i,side,peer
  complex(8)::expected
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2)error stop 1
  peer=1-rank
  do i=1,2
    side=merge(-1,1,i==1)
    send(i)%peer=peer;send(i)%epoch=5;send(i)%local_fragment=rank+1
    send(i)%neighbor_fragment=peer+1;send(i)%axis=1;send(i)%side_from_local=side
    allocate(send(i)%w_ids(1),send(i)%p_ids(1));send(i)%w_ids=[1];send(i)%p_ids=[1]
    allocate(send(i)%grid(3,2),send(i)%w(1,2),send(i)%grad_w(3,1,2))
    allocate(send(i)%p(1,2),send(i)%grad_p(3,1,2))
    send(i)%grid=reshape([rank*2+merge(1,2,side<0),1,1,&
      rank*2+merge(1,2,side<0),2,1],[3,2])
    send(i)%w=cmplx(100*rank+10*i,rank,8);send(i)%grad_w=0
    send(i)%p=cmplx(1000*rank+10*i,0d0,8);send(i)%grad_p=0
  enddo
  call exchange_dg_wpw_face_side_schedule(MPI_COMM_WORLD,5,send,remote,info)
  if(info/=0.or.size(remote)/=2.or.any(.not.remote%valid))error stop 2
  do i=1,2
    if(remote(i)%side_from_remote/=-send(i)%side_from_local)error stop 3
    if(any(remote(i)%w_ids/=[1]).or.any(remote(i)%p_ids/=[1]))error stop 3
    expected=cmplx(100*peer+10*(3-i),peer,8)
    if(any(abs(remote(i)%w-expected)>1d-13))error stop 4
  enddo
  if(rank==0)print '(a)','PASS two-fragment periodic plus/minus side traces remain distinct'
  do i=1,2;send(i)%epoch=6;enddo
  if(rank==1)send(2)%grid(2,1)=9
  call exchange_dg_wpw_face_side_schedule(MPI_COMM_WORLD,6,send,remote,info)
  if(info==0)error stop 5
  call MPI_Finalize(ierr)
end program
