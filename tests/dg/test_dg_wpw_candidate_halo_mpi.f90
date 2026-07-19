program test_dg_wpw_candidate_halo_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use rt_dg_wpw_candidate_halo,only:s_dg_wpw_staged_candidate,s_dg_wpw_owned_candidates,&
    route_dg_wpw_candidate_halo,wpw_candidate_kind_wp,wpw_candidate_kind_pp
  implicit none
  type(s_dg_wpw_staged_candidate),allocatable::staged(:)
  type(s_dg_wpw_owned_candidates)::owned
  integer,parameter::nlarge=100000
  integer::ierr,rank,nrank,info,i,owner,p_row,source_rank
  integer::support(4)
  complex(8)::expected

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(ierr/=MPI_SUCCESS.or.nrank/=4)error stop 1
  support=[1,2,3,4]
  allocate(staged(6))
  do i=1,3
    owner=modulo(rank+i,nrank)
    p_row=2*owner+1
    staged(2*i-1)%kind=wpw_candidate_kind_wp;staged(2*i-1)%image_id=i
    staged(2*i-1)%wp_w=rank+1;staged(2*i-1)%wp_p=p_row
    staged(2*i-1)%wp_h=cmplx(10*(rank+1)+i,1d0,8);staged(2*i-1)%wp_s=cmplx(rank+1,i,8)
    staged(2*i)%kind=wpw_candidate_kind_pp;staged(2*i)%image_id=i
    staged(2*i)%pp_r=p_row;staged(2*i)%pp_c=2*rank+2
    staged(2*i)%pp_h=cmplx(100*(rank+1)+i,-1d0,8);staged(2*i)%pp_s=cmplx(i,rank+1,8)
  enddo
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,7,4,2,support,staged,owned,info,16)
  if(info/=0.or.owned%epoch/=7.or.size(owned%wp_w)/=3.or.size(owned%pp_r)/=3)error stop 2
  if(any(owned%wp_p/=2*rank+1).or.any(owned%pp_r/=2*rank+1))error stop 3
  do i=1,3
    owner=modulo(rank-i,nrank)
    expected=cmplx(10*(owner+1)+i,1d0,8)
    if(.not.any(abs(owned%wp_h-expected)<1d-13))error stop 4
  enddo
  if(rank==0)print '(a)','PASS bounded WPW candidate halo routes P rows to fragment owners'

  if(rank==0)staged(1)%wp_h=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,8,4,2,support,staged,owned,info,16)
  if(info==0)error stop 5
  if(rank==0)staged(1)%wp_h=cmplx(11d0,1d0,8)

  staged(3)=staged(1)
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,9,4,2,support,staged,owned,info,16)
  if(info==0)error stop 6

  deallocate(staged);allocate(staged(4))
  owner=modulo(rank+1,nrank);p_row=2*owner+1
  do i=1,2
    staged(2*i-1)%kind=wpw_candidate_kind_wp;staged(2*i-1)%image_id=i
    staged(2*i-1)%wp_w=rank+1;staged(2*i-1)%wp_p=p_row
    staged(2*i-1)%wp_h=cmplx(dble(i),0d0,8);staged(2*i-1)%wp_s=cmplx(0d0,dble(i),8)
    staged(2*i)%kind=wpw_candidate_kind_pp;staged(2*i)%image_id=i
    staged(2*i)%pp_r=p_row;staged(2*i)%pp_c=2*rank+2
    staged(2*i)%pp_h=cmplx(10d0*dble(i),0d0,8);staged(2*i)%pp_s=cmplx(0d0,10d0*dble(i),8)
  enddo
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,10,4,2,support,staged,owned,info,16)
  if(info/=0.or.size(owned%wp_w)/=1.or.size(owned%pp_r)/=1)error stop 7
  if(abs(owned%wp_h(1)-(3d0,0d0))>1d-13.or.abs(owned%pp_h(1)-(30d0,0d0))>1d-13)error stop 8

  deallocate(staged)
  if(rank==0)then
    allocate(staged(4))
    do i=1,2
      staged(2*i-1)%kind=wpw_candidate_kind_wp;staged(2*i-1)%image_id=i
      staged(2*i-1)%wp_w=i;staged(2*i-1)%wp_p=3
      staged(2*i-1)%wp_h=cmplx(dble(i),0d0,8);staged(2*i-1)%wp_s=cmplx(dble(i),0d0,8)
      staged(2*i)%kind=wpw_candidate_kind_pp;staged(2*i)%image_id=i
      staged(2*i)%pp_r=3;staged(2*i)%pp_c=i
      staged(2*i)%pp_h=cmplx(dble(i),0d0,8);staged(2*i)%pp_s=cmplx(dble(i),0d0,8)
    enddo
  else
    allocate(staged(0))
  endif
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,11,4,2,support,staged,owned,info,merge(4,1,rank==0))
  if(info/=0)error stop 9
  if(rank==1.and.(size(owned%wp_w)/=2.or.size(owned%pp_r)/=2))error stop 10

  deallocate(staged);allocate(staged(nlarge))
  owner=modulo(rank+1,nrank);p_row=2*owner+1
  do i=1,nlarge
    staged(i)%kind=wpw_candidate_kind_wp;staged(i)%image_id=i
    staged(i)%wp_w=rank*nlarge+nlarge-i+1;staged(i)%wp_p=p_row
    staged(i)%wp_h=cmplx(dble(i),0d0,8);staged(i)%wp_s=(0d0,0d0)
  enddo
  call route_dg_wpw_candidate_halo(MPI_COMM_WORLD,12,4,2,support,staged,owned,info,nlarge)
  if(info/=0.or.size(owned%wp_w)/=nlarge.or.size(owned%pp_r)/=0)error stop 11
  source_rank=modulo(rank-1,nrank)
  if(owned%wp_w(1)/=source_rank*nlarge+1.or.&
     owned%wp_w(nlarge)/=(source_rank+1)*nlarge)error stop 12
  if(rank==0)print '(a)','PASS large reverse-order candidate route is subquadratic'
  call MPI_Finalize(ierr)
end program
