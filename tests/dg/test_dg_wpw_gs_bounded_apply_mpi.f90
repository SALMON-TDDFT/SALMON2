program test_dg_wpw_gs_bounded_apply_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_bounded_operator,only:s_dg_wpw_bounded_operator,initialize_dg_wpw_bounded_operator,&
    apply_h_dg_wpw_bounded,apply_s_dg_wpw_bounded,global_gram_dg_wpw_bounded,&
    fetch_dg_wpw_support_coefficients,release_dg_wpw_bounded_operator
  use dg_wpw_owner_exchange,only:peer_sets_equal,exchange_sizes_fit,collective_exchange_sizes_fit
  implicit none
  type(s_dg_wpw_bounded_operator)::op
  type(s_dg_wpw_bounded_operator)::bad_op
  integer::ierr,rank,info,i
  integer,allocatable::owned_w(:),owned_p(:),required_w(:),required_p(:),peers(:)
  integer,allocatable::ww_r(:),ww_c(:),wp_w(:),wp_p(:),pp_r(:),pp_c(:)
  integer,allocatable::missing_w(:)
  complex(8),allocatable::ww_h(:),ww_s(:),wp_h(:),wp_s(:),pp_h(:),pp_s(:)
  complex(8),allocatable::xw(:,:),xp(:,:),yw(:,:),yp(:,:),sw(:,:),sp(:,:)
  complex(8),allocatable::rw(:,:),rp(:,:)
  complex(8)::x_global(5),y_oracle(5),h(5,5),gram(1,1),gram_expected,local_gram
  complex(8),allocatable::bad_gram(:,:),bad_gram_counts(:,:)
  complex(8),allocatable::bad_xw(:,:),bad_xp(:,:),bad_yw(:,:),bad_yp(:,:)

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(.not.peer_sets_equal([1,3,5],[5,1,3]))error stop 17
  if(peer_sets_equal([1,3,5],[5,1,4]))error stop 18
  if(.not.exchange_sizes_fit([2,3],[4,1],2))error stop 19
  if(exchange_sizes_fit([huge(1)],[1],2))error stop 20
  if(rank==0)then
    owned_w=[1,2];owned_p=[1];required_w=[1,2,3];required_p=[1,2]
    ww_r=[1,2,2,3];ww_c=[1,2,3,1];ww_h=[(2d0,0d0),(3d0,0d0),(-0.1d0,0.05d0),(0.07d0,0d0)]
    ww_s=[(1d0,0d0),(1d0,0d0),(0d0,0d0),(0d0,0d0)]
    wp_w=[1,3];wp_p=[1,1];wp_h=[(0.4d0,-0.2d0),(0.3d0,0.1d0)]
    wp_s=[(0d0,0d0),(0d0,0d0)]
    pp_r=[1,1];pp_c=[1,2];pp_h=[(4d0,0d0),(0.15d0,-0.05d0)]
    pp_s=[(1d0,0d0),(0d0,0d0)]
  else
    owned_w=[3];owned_p=[2];required_w=[2,3];required_p=[1,2]
    ww_r=[3,3];ww_c=[2,3];ww_h=[conjg((-0.1d0,0.05d0)),(5d0,0d0)]
    ww_s=[(0d0,0d0),(1d0,0d0)]
    wp_w=[2,3];wp_p=[2,2];wp_h=[(-0.25d0,0.1d0),(0.35d0,-0.15d0)]
    wp_s=[(0d0,0d0),(0d0,0d0)]
    pp_r=[2,2];pp_c=[1,2];pp_h=[conjg((0.15d0,-0.05d0)),(6d0,0d0)]
    pp_s=[(0d0,0d0),(1d0,0d0)]
  endif
  peers=[1-rank]
  call initialize_dg_wpw_bounded_operator(op,MPI_COMM_WORLD,7,87123_8,'fragment_root_v1',&
    'orthonormal_ww','windowed_kg_sipg_v1',peers,owned_w,owned_p,required_w,required_p,&
    ww_r,ww_c,ww_h,ww_s,wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0)error stop 1
  call collective_exchange_sizes_fit(op%w_schedule,huge(1)/2+1,info)
  if(info==0)error stop 21
  allocate(xw(size(owned_w),1),xp(size(owned_p),1),yw(size(owned_w),1),yp(size(owned_p),1),&
    sw(size(owned_w),1),sp(size(owned_p),1))
  do i=1,size(owned_w);xw(i,1)=cmplx(0.1d0*owned_w(i),-0.03d0*owned_w(i),8);enddo
  do i=1,size(owned_p);xp(i,1)=cmplx(0.2d0*owned_p(i),0.04d0*owned_p(i),8);enddo
  allocate(rw(size(required_w),1),rp(size(required_p),1))
  call fetch_dg_wpw_support_coefficients(op,xw,xp,rw,rp,info);if(info/=0)error stop 22
  do i=1,size(required_w)
    if(abs(rw(i,1)-cmplx(0.1d0*required_w(i),-0.03d0*required_w(i),8))>1d-13)error stop 23
  enddo
  do i=1,size(required_p)
    if(abs(rp(i,1)-cmplx(0.2d0*required_p(i),0.04d0*required_p(i),8))>1d-13)error stop 24
  enddo
  call apply_h_dg_wpw_bounded(op,7,87123_8,xw,xp,yw,yp,info);if(info/=0)error stop 2
  call apply_s_dg_wpw_bounded(op,7,87123_8,xw,xp,sw,sp,info);if(info/=0)error stop 3
  if(maxval(abs(sw-xw))>1d-13.or.maxval(abs(sp-xp))>1d-13)error stop 4

  x_global(1:3)=[(cmplx(0.1d0*i,-0.03d0*i,8),i=1,3)]
  x_global(4:5)=[(cmplx(0.2d0*i,0.04d0*i,8),i=1,2)]
  h=(0d0,0d0)
  h(1,1)=2;h(2,2)=3;h(3,3)=5
  h(2,3)=(-0.1d0,0.05d0);h(3,2)=conjg(h(2,3))
  h(3,1)=h(3,1)+(0.07d0,0d0)
  h(1,4)=(0.4d0,-0.2d0);h(3,4)=(0.3d0,0.1d0)
  h(2,5)=(-0.25d0,0.1d0);h(3,5)=(0.35d0,-0.15d0)
  h(4,1)=conjg(h(1,4));h(4,3)=conjg(h(3,4));h(5,2)=conjg(h(2,5));h(5,3)=conjg(h(3,5))
  h(4,4)=4;h(4,5)=(0.15d0,-0.05d0);h(5,4)=conjg(h(4,5));h(5,5)=6
  y_oracle=matmul(h,x_global)
  do i=1,size(owned_w);if(abs(yw(i,1)-y_oracle(owned_w(i)))>1d-12)error stop 5;enddo
  do i=1,size(owned_p);if(abs(yp(i,1)-y_oracle(3+owned_p(i)))>1d-12)error stop 6;enddo
  call global_gram_dg_wpw_bounded(op,xw,xp,yw,yp,gram,info);if(info/=0)error stop 7
  local_gram=sum(conjg(xw(:,1))*yw(:,1))+sum(conjg(xp(:,1))*yp(:,1))
  call MPI_Allreduce(local_gram,gram_expected,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(abs(gram(1,1)-gram_expected)>1d-13)error stop 8
  if(rank==0)then;allocate(bad_gram(2,1));else;allocate(bad_gram(1,1));endif
  call global_gram_dg_wpw_bounded(op,xw,xp,yw,yp,bad_gram,info)
  if(info==0.or.any(abs(bad_gram)>1d-13))error stop 9
  i=merge(2,1,rank==0)
  allocate(bad_xw(size(owned_w),i),bad_xp(size(owned_p),i),bad_yw(size(owned_w),i),bad_yp(size(owned_p),i))
  bad_xw=1;bad_xp=1;bad_yw=(9d0,0d0);bad_yp=(9d0,0d0)
  call apply_h_dg_wpw_bounded(op,7,87123_8,bad_xw,bad_xp,bad_yw,bad_yp,info)
  if(info==0.or.any(abs(bad_yw)>1d-13).or.any(abs(bad_yp)>1d-13))error stop 10
  allocate(bad_gram_counts(i,i))
  call global_gram_dg_wpw_bounded(op,bad_xw,bad_xp,bad_yw,bad_yp,bad_gram_counts,info)
  if(info==0.or.any(abs(bad_gram_counts)>1d-13))error stop 11
  yw=(9d0,0d0);yp=(9d0,0d0)
  call apply_h_dg_wpw_bounded(op,merge(6,7,rank==0),87123_8,xw,xp,yw,yp,info)
  if(info==0.or.any(abs(yw)>1d-13).or.any(abs(yp)>1d-13))error stop 12
  yw=(9d0,0d0);yp=(9d0,0d0)
  call apply_h_dg_wpw_bounded(op,7,merge(87122_8,87123_8,rank==0),xw,xp,yw,yp,info)
  if(info==0.or.any(abs(yw)>1d-13).or.any(abs(yp)>1d-13))error stop 13
  if(rank==0)xw(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  yw=(9d0,0d0);yp=(9d0,0d0)
  call apply_h_dg_wpw_bounded(op,7,87123_8,xw,xp,yw,yp,info)
  if(info==0.or.any(abs(yw)>1d-13).or.any(abs(yp)>1d-13))error stop 14
  if(rank==0)xw(1,1)=cmplx(0.1d0,-0.03d0,8)
  if(rank==0)then;missing_w=[required_w,99];else;missing_w=required_w;endif
  call initialize_dg_wpw_bounded_operator(bad_op,MPI_COMM_WORLD,8,87124_8,'fragment_root_v1',&
    'orthonormal_ww','windowed_kg_sipg_v1',peers,owned_w,owned_p,missing_w,required_p,&
    ww_r,ww_c,ww_h,ww_s,wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info==0.or.bad_op%valid)error stop 15
  if(rank==0)then;peers=[1];else;deallocate(peers);allocate(peers(0));endif
  call initialize_dg_wpw_bounded_operator(bad_op,MPI_COMM_WORLD,8,87124_8,'fragment_root_v1',&
    'orthonormal_ww','windowed_kg_sipg_v1',peers,owned_w,owned_p,required_w,required_p,&
    ww_r,ww_c,ww_h,ww_s,wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info==0.or.bad_op%valid)error stop 16
  call release_dg_wpw_bounded_operator(op)
  call release_dg_wpw_bounded_operator(bad_op)
  if(rank==0)print '(a)','PASS bounded GS-neutral WW/WP/PP apply matches dense oracle'
  call MPI_Finalize(ierr)
end program
