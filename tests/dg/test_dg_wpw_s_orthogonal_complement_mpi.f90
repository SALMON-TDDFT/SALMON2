program test_dg_wpw_s_orthogonal_complement_mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use mpi
  use dg_wpw_bounded_operator,only:s_dg_wpw_bounded_operator,initialize_dg_wpw_bounded_operator,&
    release_dg_wpw_bounded_operator,set_dg_wpw_interface_lambda
  use dg_wpw_s_orthogonal_complement,only:s_dg_wpw_s_orthogonal_complement,&
    initialize_dg_wpw_s_orthogonal_complement,apply_h_dg_wpw_s_orthogonal_complement,&
    apply_s_dg_wpw_s_orthogonal_complement,map_dg_wpw_complement_to_original,&
    map_dg_wpw_original_to_complement,validate_dg_wpw_s_orthogonal_complement,&
    release_dg_wpw_s_orthogonal_complement
  implicit none
  type(s_dg_wpw_bounded_operator)::op,bad_op,rank_op,deflated_op,dynamic_op
  type(s_dg_wpw_s_orthogonal_complement)::transform,failed,rank_transform,deflated_transform,dynamic_transform
  integer::rank,ierr,info,owned,irow,ww_r(2),ww_c(2),wp_w(2),wp_p(2),pp_r(2),pp_c(2),ipiv(2)
  integer::peers(1),owned_w(1),owned_p(1),required(2)
  complex(8)::sww(2,2),swp(2,2),spp(2,2),hww(2,2),hwp(2,2),hpp(2,2),s(4,4),h(4,4),st(4,4),ht(4,4)
  complex(8)::ww_s(2),wp_s(2),pp_s(2),ww_h(2),wp_h(2),pp_h(2),a(2,2),rhs(2,2),t(4,4),spp_rank(2,2)
  complex(8)::x(4,2),original(4,2),roundtrip(4,2),yw(1,2),yp(1,2),local(4,2),actual(4,2)
  complex(8),allocatable::xw(:,:),xp(:,:),mw(:,:),mp(:,:)
  integer,allocatable::owned_w3(:),ww_r3(:),ww_c3(:)
  complex(8),allocatable::ww_s3(:),ww_h3(:)
  integer::required_w3(3),wp_w3(3),wp_p3(3),pp_r3(2),pp_c3(2)
  complex(8)::wp_s3(3),wp_h3(3),pp_s3(2),pp_h3(2),expected_a3(3,2)
  real(8)::scale

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(rank==0)then;peers=[1];else;peers=[0];endif
  owned=rank+1;owned_w=[owned];owned_p=[owned];required=[1,2]
  sww=reshape([(2.0d0,0d0),(0.25d0,-0.10d0),(0.25d0,0.10d0),(1.5d0,0d0)],[2,2])
  swp=reshape([(0.30d0,0.05d0),(-0.12d0,0.08d0),(0.15d0,-0.04d0),(0.22d0,0.03d0)],[2,2])
  spp=reshape([(1.30d0,0d0),(0.07d0,-0.02d0),(0.07d0,0.02d0),(1.10d0,0d0)],[2,2])
  hww=reshape([(3.0d0,0d0),(0.2d0,0.1d0),(0.2d0,-0.1d0),(2.4d0,0d0)],[2,2])
  hwp=reshape([(0.8d0,0.1d0),(-0.3d0,0.2d0),(0.25d0,-0.15d0),(0.6d0,0.05d0)],[2,2])
  hpp=reshape([(1.8d0,0d0),(0.12d0,0.03d0),(0.12d0,-0.03d0),(2.1d0,0d0)],[2,2])
  ww_r=owned;ww_c=[1,2];ww_s=sww(owned,:);ww_h=hww(owned,:)
  wp_w=[1,2];wp_p=owned;wp_s=swp(:,owned);wp_h=hwp(:,owned)
  pp_r=owned;pp_c=[1,2];pp_s=spp(owned,:);pp_h=hpp(owned,:)
  call initialize_dg_wpw_bounded_operator(op,MPI_COMM_WORLD,7,12345_8,'canonical_owner',&
    'real_space_h1','bounded_hs',peers,owned_w,owned_p,required,required,ww_r,ww_c,ww_h,ww_s,&
    wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0)error stop 1
  call initialize_dg_wpw_s_orthogonal_complement(op,1d-12,transform,info)
  if(info/=0.or..not.transform%valid)error stop 2
  if(transform%solve_residual>1d-11.or.transform%cross_metric_defect>1d-11)error stop 3
  rhs=swp;a=sww;call zgesv(2,2,a,2,ipiv,rhs,2,info);if(info/=0)error stop 4;a=rhs
  t=0;t(1,1)=1;t(2,2)=1;t(3,3)=1;t(4,4)=1;t(1:2,3:4)=-a
  s=0;s(1:2,1:2)=sww;s(1:2,3:4)=swp;s(3:4,1:2)=conjg(transpose(swp));s(3:4,3:4)=spp
  h=0;h(1:2,1:2)=hww;h(1:2,3:4)=hwp;h(3:4,1:2)=conjg(transpose(hwp));h(3:4,3:4)=hpp
  scale=max(1d0,maxval(abs(s)))
  st=matmul(conjg(transpose(t)),matmul(s,t));ht=matmul(conjg(transpose(t)),matmul(h,t))
  if(maxval(abs(st(1:2,3:4)))/scale>1d-11)error stop 5
  x=reshape([(0.2d0,0.1d0),(-0.3d0,0.4d0),(0.7d0,-0.2d0),(0.1d0,0.6d0),&
    (-0.5d0,0.2d0),(0.9d0,-0.1d0),(0.4d0,0.3d0),(-0.2d0,-0.4d0)],[4,2])
  allocate(xw(1,2),xp(1,2),mw(1,2),mp(1,2))
  xw(1,:)=x(owned,:);xp(1,:)=x(2+owned,:)
  call map_dg_wpw_complement_to_original(op,transform,xw,xp,mw,mp,info);if(info/=0)error stop 6
  local=0;local(owned,:)=mw(1,:);local(2+owned,:)=mp(1,:)
  call MPI_Allreduce(local,original,8,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(original-matmul(t,x)))>1d-11)error stop 7
  call map_dg_wpw_original_to_complement(op,transform,mw,mp,xw,xp,info);if(info/=0)error stop 8
  local=0;local(owned,:)=xw(1,:);local(2+owned,:)=xp(1,:)
  call MPI_Allreduce(local,roundtrip,8,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(roundtrip-x))>1d-11)error stop 9
  xw(1,:)=x(owned,:);xp(1,:)=x(2+owned,:)
  call apply_s_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,yw,yp,info);if(info/=0)error stop 10
  local=0;local(owned,:)=yw(1,:);local(2+owned,:)=yp(1,:)
  call MPI_Allreduce(local,actual,8,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(actual-matmul(conjg(transpose(t)),matmul(s,matmul(t,x)))))>1d-11)error stop 11
  call apply_h_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,yw,yp,info);if(info/=0)error stop 12
  local=0;local(owned,:)=yw(1,:);local(2+owned,:)=yp(1,:)
  call MPI_Allreduce(local,actual,8,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(actual-matmul(conjg(transpose(t)),matmul(h,matmul(t,x)))))>1d-11)error stop 13
  if(maxval(abs(ht(1:2,3:4)))<1d-3)error stop 14
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)error stop 15
  op%wp_h_dense=op%wp_h_dense+(0.125d0,0d0)
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)error stop 151
  op%wp_h_dense=op%wp_h_dense-(0.125d0,0d0)
  call set_dg_wpw_interface_lambda(op,0.5d0,info);if(info/=0)error stop 152
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)error stop 153
  if(rank==0)op%wp_s_dense(1,1)=op%wp_s_dense(1,1)+(1d-4,0d0)
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info==0)error stop 16
  if(rank==0)op%wp_s_dense(1,1)=op%wp_s_dense(1,1)-(1d-4,0d0)
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)error stop 17
  op%pp_s_dense=op%pp_s_dense+(1d-4,0d0)
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info==0)error stop 171
  op%pp_s_dense=op%pp_s_dense-(1d-4,0d0)
  call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)error stop 172
  bad_op=op
  if(rank==0)bad_op%ww_s_dense(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call initialize_dg_wpw_s_orthogonal_complement(bad_op,1d-12,failed,info);if(info==0)error stop 18
  bad_op=op
  if(rank==0)bad_op%ww_s_dense(1,2)=bad_op%ww_s_dense(1,2)+(0.1d0,0.05d0)
  call initialize_dg_wpw_s_orthogonal_complement(bad_op,1d-12,failed,info);if(info==0)error stop 181
  bad_op=op;bad_op%ww_s_dense=1
  call initialize_dg_wpw_s_orthogonal_complement(bad_op,1d-12,failed,info);if(info==0)error stop 19
  bad_op=op;bad_op%ww_s_dense=0
  bad_op%ww_s_dense(owned,owned)=merge(1d0,1d-14,rank==0)
  call initialize_dg_wpw_s_orthogonal_complement(bad_op,1d-12,failed,info);if(info==0)error stop 190
  spp_rank=matmul(conjg(transpose(swp)),a)
  spp_rank=spp_rank+reshape([(1d0,0d0),(0.4d0,-0.2d0),(0.4d0,0.2d0),(0.2d0,0d0)],[2,2])
  pp_s=spp_rank(owned,:)
  call initialize_dg_wpw_bounded_operator(rank_op,MPI_COMM_WORLD,8,12346_8,'canonical_owner',&
    'real_space_h1','bounded_hs',peers,owned_w,owned_p,required,required,ww_r,ww_c,ww_h,ww_s,&
    wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0)error stop 191
  call initialize_dg_wpw_s_orthogonal_complement(rank_op,1d-12,rank_transform,info)
  if(info/=0.or.rank_transform%numerical_p_rank/=1)error stop 192
  call release_dg_wpw_s_orthogonal_complement(rank_transform)
  call release_dg_wpw_bounded_operator(rank_op)
  swp(:,2)=0;rhs=swp;a=sww;call zgesv(2,2,a,2,ipiv,rhs,2,info);if(info/=0)error stop 193;a=rhs
  spp_rank=matmul(conjg(transpose(swp)),a);spp_rank(1,1)=spp_rank(1,1)+1;spp_rank(2,2)=spp_rank(2,2)+1
  wp_s=swp(:,owned);pp_s=spp_rank(owned,:)
  call initialize_dg_wpw_bounded_operator(deflated_op,MPI_COMM_WORLD,9,12347_8,'canonical_owner',&
    'real_space_h1','bounded_hs',peers,owned_w,owned_p,required,required,ww_r,ww_c,ww_h,ww_s,&
    wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0)error stop 194
  call initialize_dg_wpw_s_orthogonal_complement(deflated_op,1d-12,deflated_transform,info)
  if(info/=0.or.maxval(abs(deflated_transform%a_owned_w_global_p(:,2)))>1d-12)error stop 195
  call release_dg_wpw_s_orthogonal_complement(deflated_transform)
  call release_dg_wpw_bounded_operator(deflated_op)
  required_w3=[1,2,3];wp_w3=required_w3;wp_p3=rank+1;pp_r3=rank+1;pp_c3=[1,2]
  if(rank==0)then
    allocate(owned_w3(2),ww_r3(6),ww_c3(6),ww_s3(6),ww_h3(6));owned_w3=[1,2]
    ww_r3=[1,1,1,2,2,2];ww_c3=[1,2,3,1,2,3]
    ww_s3=[(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),(2d0,0d0),(0.5d0,0d0)]
  else
    allocate(owned_w3(1),ww_r3(3),ww_c3(3),ww_s3(3),ww_h3(3));owned_w3=[3]
    ww_r3=[3,3,3];ww_c3=[1,2,3];ww_s3=[(0d0,0d0),(0.5d0,0d0),(4d0,0d0)]
  endif
  ww_h3=ww_s3
  if(rank==0)then
    wp_s3=[(1d0,0d0),(0d0,0d0),(0d0,0d0)];pp_s3=[(2d0,0d0),(0d0,0d0)]
  else
    wp_s3=[(0d0,0d0),(1d0,0d0),(1d0,0d0)];pp_s3=[(0d0,0d0),(1.6451612903225806d0,0d0)]
  endif
  wp_h3=wp_s3;pp_h3=pp_s3
  call initialize_dg_wpw_bounded_operator(dynamic_op,MPI_COMM_WORLD,11,12349_8,'canonical_owner',&
    'real_space_h1','bounded_hs',peers,owned_w3,owned_p,required_w3,required,ww_r3,ww_c3,ww_h3,ww_s3,&
    wp_w3,wp_p3,wp_h3,wp_s3,pp_r3,pp_c3,pp_h3,pp_s3,info)
  if(info/=0)error stop 199
  call initialize_dg_wpw_s_orthogonal_complement(dynamic_op,1d-12,dynamic_transform,info)
  if(info/=0)error stop 1991
  expected_a3=0;expected_a3(1,1)=1;expected_a3(2,2)=3.5d0/7.75d0;expected_a3(3,2)=1.5d0/7.75d0
  do irow=1,size(owned_w3)
    if(maxval(abs(dynamic_transform%a_owned_w_global_p(irow,:)-expected_a3(owned_w3(irow),:)))>1d-11)&
      error stop 1992
  enddo
  call release_dg_wpw_s_orthogonal_complement(dynamic_transform)
  call release_dg_wpw_bounded_operator(dynamic_op)
  swp(:,2)=swp(:,1);rhs=swp;a=sww;call zgesv(2,2,a,2,ipiv,rhs,2,info);if(info/=0)error stop 196;a=rhs
  spp_rank=matmul(conjg(transpose(swp)),a);spp_rank(1,1)=spp_rank(1,1)+1;spp_rank(2,2)=spp_rank(2,2)+1
  wp_s=swp(:,owned);pp_s=spp_rank(owned,:)
  call initialize_dg_wpw_bounded_operator(deflated_op,MPI_COMM_WORLD,10,12348_8,'canonical_owner',&
    'real_space_h1','bounded_hs',peers,owned_w,owned_p,required,required,ww_r,ww_c,ww_h,ww_s,&
    wp_w,wp_p,wp_h,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0)error stop 197
  call initialize_dg_wpw_s_orthogonal_complement(deflated_op,1d-12,deflated_transform,info)
  if(info/=0.or.maxval(abs(deflated_transform%a_owned_w_global_p(:,1)-&
    deflated_transform%a_owned_w_global_p(:,2)))>1d-12)error stop 198
  call release_dg_wpw_s_orthogonal_complement(deflated_transform)
  call release_dg_wpw_bounded_operator(deflated_op)
  deallocate(xw,xp,mw,mp);allocate(xw(1,merge(1,2,rank==0)),xp(1,merge(1,2,rank==0)),&
    mw(1,merge(1,2,rank==0)),mp(1,merge(1,2,rank==0)))
  xw=0;xp=0
  call map_dg_wpw_complement_to_original(op,transform,xw,xp,mw,mp,info);if(info==0)error stop 20
  deallocate(mw);allocate(mw(2,size(xw,2)))
  call apply_s_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,mw,mp,info)
  if(info==0)error stop 201
  call release_dg_wpw_s_orthogonal_complement(transform);if(transform%valid)error stop 21
  deallocate(xw,xp,mw,mp);allocate(xw(1,1),xp(1,1),mw(1,1),mp(1,1));xw=0;xp=0
  call apply_s_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,mw,mp,info)
  if(info==0)error stop 22
  call release_dg_wpw_bounded_operator(op)
  if(rank==0)write(*,'(a)')'PASS distributed S-orthogonal PW complement'
  call MPI_Finalize(ierr)
end program
