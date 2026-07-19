program test_dg_wpw_face_side_reduce_mpi
  use mpi
  use rt_dg_wpw_face_side_halo,only:s_dg_wpw_face_side_send,reduce_dg_wpw_face_side_parts
  implicit none
  type(s_dg_wpw_face_side_send)::side
  integer::ierr,rank,info,coverage(2),grid(3,2)
  integer::w_ids(1),p_ids(1)
  complex(8)::w(1,2),gw(3,1,2),p(1,2),gp(3,1,2)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  grid=reshape([2,1,1,2,2,1],[3,2]);w_ids=[4];p_ids=[7]
  coverage=0;coverage(rank+1)=1;w=0;gw=0;p=0;gp=0
  w(1,rank+1)=cmplx(rank+1,0d0,8);p(1,rank+1)=cmplx(10+rank,0d0,8)
  gw(1,1,rank+1)=cmplx(20+rank,0d0,8);gp(1,1,rank+1)=cmplx(30+rank,0d0,8)
  call reduce_dg_wpw_face_side_parts(MPI_COMM_WORLD,0,1,3,1,2,1,1,w_ids,p_ids,grid,coverage,&
    w,gw,p,gp,side,info)
  if(info/=0)error stop 1
  if(rank==0)then
    if(side%epoch/=3.or.side%peer/=1.or.any(side%w_ids/=w_ids).or.any(side%p_ids/=p_ids))error stop 2
    if(abs(side%w(1,1)-(1d0,0d0))>1d-13.or.abs(side%w(1,2)-(2d0,0d0))>1d-13)error stop 3
  endif
  coverage=0;if(rank==0)coverage(1)=1
  call reduce_dg_wpw_face_side_parts(MPI_COMM_WORLD,0,1,4,1,2,1,1,w_ids,p_ids,grid,coverage,&
    w,gw,p,gp,side,info)
  if(info==0)error stop 4
  if(rank==0.and.side%epoch/=3)error stop 5
  if(rank==0)print '(a)','PASS rank-local face side publishes only after complete coverage'
  call MPI_Finalize(ierr)
end program
