program test_dg_wpw_production_operator_mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use mpi
  use dg_wpw_production_context, only: s_dg_wpw_production_context,&
    s_dg_wpw_production_context_snapshot,initialize_dg_wpw_fragment_root_context,&
    snapshot_dg_wpw_production_context,validate_dg_wpw_production_context_snapshot,&
    release_dg_wpw_production_context_snapshot
  use dg_wpw_lcfo_operator_adapter,only:s_dg_wpw_lcfo_ww_components,import_dg_wpw_lcfo_ww_components,&
    publish_dg_wpw_lcfo_ww_components
  use dg_wpw_bounded_operator,only:apply_h_dg_wpw_bounded,apply_s_dg_wpw_bounded,&
    set_dg_wpw_interface_lambda,release_dg_wpw_bounded_operator,&
    s_dg_wpw_bounded_operator_snapshot,snapshot_dg_wpw_bounded_operator,&
    validate_dg_wpw_bounded_operator_snapshot,release_dg_wpw_bounded_operator_snapshot,&
    reduce_dg_wpw_metric_rhs_partials,s_dg_wpw_fragment_block_preconditioner,&
    initialize_dg_wpw_fragment_block_preconditioner,apply_dg_wpw_fragment_block_preconditioner,&
    release_dg_wpw_fragment_block_preconditioner
  use rt_dg_wpw_trace_halo_provider, only: s_dg_wpw_trace_halo_state, prepare_dg_wpw_trace_halo
  use rt_dg_wpw_face_trace_provider, only: s_wpw_face_trace_provider
  use rt_dg_wpw_production_builder, only: build_dg_wpw_rank_local_quadrature, &
    scan_and_route_dg_wpw_canonical_faces, replace_dg_wpw_potential_volume, &
    build_dg_wpw_production_operator, bind_dg_wpw_hs_callbacks,install_dg_wpw_projector_nonlocal
  implicit none
  type(s_dg_wpw_production_context) :: context
  type(s_dg_wpw_bounded_operator_snapshot)::frozen_operator
  type(s_dg_wpw_fragment_block_preconditioner)::block_preconditioner
  type(s_dg_wpw_production_context_snapshot)::frozen_context
  type(s_dg_wpw_lcfo_ww_components)::ww_components
  type(s_wpw_face_trace_provider) :: provider
  type(s_dg_wpw_trace_halo_state),target :: halo
  integer :: ierr, rank, info, owned_p
  integer(8)::old_fingerprint
  real(8)::old_ww_potential
  integer,allocatable::fkminus(:),fkplus(:),faxis(:),fside(:)
  integer,allocatable::fminus_lo(:,:),fminus_hi(:,:),fplus_lo(:,:),fplus_hi(:,:)
  integer :: wp_w(2), wp_p(2), pp_r(2), pp_c(2)
  complex(8) :: wp_h_volume(2), wp_s(2), pp_h(2), pp_s(2)
  complex(8) :: saved_face(2),saved_total(2),new_wp_volume(2),new_pp_volume(2)
  complex(8) :: x(4), yh_local(4), ys_local(4), yh(4), ys(4)
  complex(8) :: h_oracle(4,4), s_oracle(4,4), wp(2,2), swp(2,2), pp(2,2), spp(2,2)
  complex(8) :: h0_oracle(4,4),interface_oracle(4,4)
  integer,allocatable::saved_owned_w(:),saved_required_w(:),saved_owned_p(:),saved_required_p(:)
  complex(8) :: face_delta_local(2,2),face_delta(2,2)
  complex(8)::ww_nl_local(1,1),ww_nl_cross(1),wp_nl(2),pp_nl(2)
  complex(8)::saved_ww_nl_local(1,1),saved_ww_nl_cross(1),saved_wp_nl(2),saved_pp_nl(2)
  complex(8)::saved_wp_volume(2),saved_wp_total(2)
  complex(8)::saved_h0_cache,saved_interface_cache,saved_trial_h
  complex(8)::saved_context_value
  integer::saved_transport_id
  complex(8)::bxw(1,1),bxp(1,1),byw(1,1),byp(1,1),bsw(1,1),bsp(1,1),bounded_local(4),bounded(4)
  complex(8)::block_rhs_w(1,2),block_rhs_p(1,2),block_z_w(1,2),block_z_p(1,2),&
    block_expected(2,2),block_solution(2,2),block_matrix(2,2),saved_block_z_w(1,2),&
    saved_block_z_p(1,2),saved_metric_value
  complex(8)::rhs_partial_w(2,2),rhs_partial_p(2,2),rhs_owned_w(1,2),rhs_owned_p(1,2)
  complex(8) :: wm(2,1),wplus(2,1),gwm(3,2,1),gwplus(3,2,1)
  complex(8) :: pm(2,1),pplus(2,1),gpm(3,2,1),gpplus(3,2,1)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(ierr/=MPI_SUCCESS) error stop 10
  call initialize_dg_wpw_fragment_root_context(context,MPI_COMM_WORLD,2,1,2,rank+1,[1,2],[rank+1],[1,2],info)
  if(info/=0) error stop 11
  owned_p=rank+1
  wp=reshape([(0.8d0,0.1d0),(0.2d0,0.3d0),(0.1d0,-0.2d0),(0.6d0,0.1d0)],[2,2])
  swp=0.05d0*wp
  pp=reshape([(2d0,0d0),(0d0,0.1d0),(0d0,-0.1d0),(3d0,0d0)],[2,2])
  spp=reshape([(1d0,0d0),(0d0,0.01d0),(0d0,-0.01d0),(1.1d0,0d0)],[2,2])
  wp_w=[1,2]; wp_p=owned_p; pp_r=owned_p; pp_c=[1,2]
  wp_h_volume=wp(:,owned_p)
  wp_s=swp(:,owned_p); pp_h=pp(owned_p,:); pp_s=spp(owned_p,:)
  call build_dg_wpw_rank_local_quadrature(context,wp_w,wp_p,wp_h_volume,wp_s,pp_r,pp_c,pp_h,pp_s,info)
  if(info/=0) error stop 12
  if(maxval(abs(context%wp_h_volume-wp_h_volume))>1d-14.or.&
    maxval(abs(context%wp_h_face))>1d-14.or.maxval(abs(context%wp_h-wp_h_volume))>1d-14)error stop 121
  if(maxval(abs(context%pp_h_volume-pp_h))>1d-14.or.maxval(abs(context%pp_h-pp_h))>1d-14)error stop 122
  wm(:,1)=[(1d0,0d0),(0.5d0,0.1d0)];wplus(:,1)=wm(:,1)+[(0.2d0,0d0),(-0.1d0,0.05d0)]
  gwm=(0d0,0d0);gwplus=(0d0,0d0);gwm(1,:,1)=[(0.1d0,0d0),(0.2d0,0d0)];gwplus(1,:,1)=-gwm(1,:,1)
  pm(:,1)=[(0.8d0,0.2d0),(0.4d0,-0.1d0)];pplus=pm;gpm=(0d0,0d0);gpplus=(0d0,0d0)
  gpm(1,:,1)=[(0.05d0,0d0),(-0.02d0,0d0)];gpplus=gpm
  call prepare_dg_wpw_trace_halo(context,halo,provider,1,1,2,1,1,reshape([1,1,1],[3,1]),&
    wm(1:1,:),wplus(1:1,:),gwm(:,1:1,:),gwplus(:,1:1,:),pm,pplus,gpm,gpplus,info)
  if(info/=0) error stop 13
  if(rank==0)then
    allocate(fkminus(1),fkplus(1),faxis(1),fside(1),fminus_lo(3,1),fminus_hi(3,1),&
      fplus_lo(3,1),fplus_hi(3,1))
    fkminus=1;fkplus=2;faxis=1;fside=1
    fminus_lo=1;fminus_hi=1;fplus_lo=1;fplus_hi=1
  else
    allocate(fkminus(0),fkplus(0),faxis(0),fside(0),fminus_lo(3,0),fminus_hi(3,0),&
      fplus_lo(3,0),fplus_hi(3,0))
  endif
  call scan_and_route_dg_wpw_canonical_faces(context,1,[1,2],provider,fkminus,fkplus,faxis,fside,&
    fminus_lo,fminus_hi,fplus_lo,fplus_hi,[1d0,1d0,1d0],wp_w(1:1),[1,2],8,info)
  if(info/=0) error stop 14
  if(maxval(abs(context%wp_h-context%wp_h_volume-context%wp_h_face))>1d-14)error stop 143
  if(maxval(abs(context%wp_h_volume-wp_h_volume))>1d-14)error stop 144
  face_delta_local=(0d0,0d0);face_delta_local(:,owned_p)=context%wp_h-wp_h_volume
  call MPI_Allreduce(face_delta_local,face_delta,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  wp=wp+face_delta
  call import_dg_wpw_lcfo_ww_components(rank+1,[1,1],reshape([1,2],[1,2]),reshape([1d0],[1,1]),&
    reshape([2d0],[1,1]),reshape([3d0],[1,1]),reshape([4d0],[1,1]),[integer::],[integer::],&
    [integer::],[integer::],[integer::],reshape([integer::],[3,0]),[real(8)::],'orthonormal_ww',ww_components,info)
  if(info/=0)error stop 141
  call publish_dg_wpw_lcfo_ww_components(context,ww_components,info)
  if(info/=0)error stop 142
  ww_nl_local=(0.5d0,0d0);ww_nl_cross=(0.25d0,0d0);wp_nl=0;pp_nl=0
  call install_dg_wpw_projector_nonlocal(context,ww_nl_local,[rank+1],[2-rank],ww_nl_cross,&
    wp_nl,pp_nl,info)
  if(info/=0)error stop 149
  saved_ww_nl_local=context%ww_projector_nonlocal
  saved_ww_nl_cross=context%ww_projector_cross_value
  saved_wp_nl=context%wp_h_nonlocal;saved_pp_nl=context%pp_h_nonlocal
  if(rank==0)ww_nl_cross(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call install_dg_wpw_projector_nonlocal(context,ww_nl_local,[rank+1],[2-rank],ww_nl_cross,&
    wp_nl,pp_nl,info)
  if(info==0.or.any(context%ww_projector_nonlocal/=saved_ww_nl_local).or.&
    any(context%ww_projector_cross_value/=saved_ww_nl_cross).or.&
    any(context%wp_h_nonlocal/=saved_wp_nl).or.any(context%pp_h_nonlocal/=saved_pp_nl))error stop 1491
  ww_nl_cross=(0.25d0,0d0)
  saved_wp_volume=context%wp_h_volume;saved_wp_total=context%wp_h
  context%wp_h_volume=cmplx(huge(0d0),0d0,8);wp_nl=cmplx(huge(0d0),0d0,8)
  call install_dg_wpw_projector_nonlocal(context,ww_nl_local,[rank+1],[2-rank],ww_nl_cross,&
    wp_nl,pp_nl,info)
  if(info==0.or.any(context%wp_h/=saved_wp_total).or.&
    any(context%ww_projector_nonlocal/=saved_ww_nl_local))error stop 1492
  context%wp_h_volume=saved_wp_volume;context%wp_h=saved_wp_total;wp_nl=0
  call build_dg_wpw_production_operator(context,info)
  if(info/=0) error stop 15
  if(.not.context%bounded_operator%valid)error stop 151
  call initialize_dg_wpw_fragment_block_preconditioner(context%bounded_operator,1d-8,&
    block_preconditioner,info)
  if(info/=0.or..not.block_preconditioner%valid.or.block_preconditioner%dimension/=2.or.&
    block_preconditioner%condition_number<1d0)error stop 1501
  block_rhs_w=reshape([(1d0,0.25d0),(-0.5d0,0.2d0)],[1,2])
  block_rhs_p=reshape([(0.3d0,-0.1d0),(0.8d0,0.4d0)],[1,2])
  call apply_dg_wpw_fragment_block_preconditioner(context%bounded_operator,block_preconditioner,&
    block_rhs_w,block_rhs_p,block_z_w,block_z_p,info)
  if(info/=0)error stop 1502
  block_matrix=0
  block_matrix(1,1)=context%bounded_operator%ww_s_dense(rank+1,rank+1)
  block_matrix(1,2)=context%bounded_operator%wp_s_dense(rank+1,1)
  block_matrix(2,1)=conjg(block_matrix(1,2))
  block_matrix(2,2)=context%bounded_operator%pp_s_dense(1,rank+1)
  block_expected(1,:)=block_rhs_w(1,:)
  block_expected(2,:)=block_rhs_p(1,:)
  block_solution(1,:)=block_z_w(1,:);block_solution(2,:)=block_z_p(1,:)
  if(maxval(abs(matmul(block_matrix,block_solution)-block_expected))>1d-11)error stop 1503
  saved_block_z_w=block_z_w;saved_block_z_p=block_z_p
  saved_metric_value=context%bounded_operator%ww_s_dense(rank+1,rank+1)
  if(rank==0)context%bounded_operator%ww_s_dense(1,1)=cmplx(-1d0,0d0,8)
  call initialize_dg_wpw_fragment_block_preconditioner(context%bounded_operator,1d-8,&
    block_preconditioner,info)
  if(info==0.or..not.block_preconditioner%valid)error stop 1504
  context%bounded_operator%ww_s_dense(rank+1,rank+1)=saved_metric_value
  call apply_dg_wpw_fragment_block_preconditioner(context%bounded_operator,block_preconditioner,&
    block_rhs_w,block_rhs_p,block_z_w,block_z_p,info)
  if(info/=0.or.maxval(abs(block_z_w-saved_block_z_w))>1d-13.or.&
    maxval(abs(block_z_p-saved_block_z_p))>1d-13)error stop 1505
  rhs_partial_w=reshape([cmplx(1+rank,0d0,8),cmplx(10*(1+rank),0d0,8),&
    cmplx(2*(1+rank),0d0,8),cmplx(20*(1+rank),0d0,8)],[2,2])
  rhs_partial_p=cmplx(3d0,0d0,8)*rhs_partial_w
  call reduce_dg_wpw_metric_rhs_partials(context%bounded_operator,rhs_partial_w,rhs_partial_p,&
    rhs_owned_w,rhs_owned_p,info)
  if(info/=0)error stop 1511
  if(maxval(abs(rhs_owned_w(1,:)-merge([cmplx(3d0,0d0,8),cmplx(6d0,0d0,8)],&
    [cmplx(30d0,0d0,8),cmplx(60d0,0d0,8)],rank==0)))>1d-12)error stop 1512
  if(maxval(abs(rhs_owned_p-3d0*rhs_owned_w))>1d-12)error stop 1513
  call snapshot_dg_wpw_bounded_operator(context%bounded_operator,frozen_operator,info)
  if(info/=0)error stop 152
  saved_h0_cache=context%bounded_operator%ww_h0_dense(1,1)
  saved_interface_cache=context%bounded_operator%wp_interface_dense(1,1)
  saved_transport_id=context%bounded_operator%owned_w_ids(1)
  if(rank==0)context%bounded_operator%ww_h0_dense(1,1)=&
    context%bounded_operator%ww_h0_dense(1,1)+(0.125d0,0d0)
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info==0)error stop 153
  if(rank==0)context%bounded_operator%ww_h0_dense(1,1)=saved_h0_cache
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info/=0)error stop 154
  if(rank==0)then
    deallocate(context%bounded_operator%ww_h0_dense)
    allocate(context%bounded_operator%ww_h0_dense(0,0))
  endif
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info==0)error stop 1541
  if(rank==0)then
    deallocate(context%bounded_operator%ww_h0_dense)
    allocate(context%bounded_operator%ww_h0_dense,source=frozen_operator%ww_h0_dense)
  endif
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info/=0)error stop 1542
  if(rank==1)context%bounded_operator%wp_interface_dense(1,1)=&
    context%bounded_operator%wp_interface_dense(1,1)+(0.25d0,0d0)
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info==0)error stop 155
  if(rank==1)context%bounded_operator%wp_interface_dense(1,1)=saved_interface_cache
  if(rank==0)context%bounded_operator%owned_w_ids(1)=context%bounded_operator%owned_w_ids(1)+100
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info==0)error stop 156
  if(rank==0)context%bounded_operator%owned_w_ids(1)=saved_transport_id
  call validate_dg_wpw_bounded_operator_snapshot(context%bounded_operator,frozen_operator,info)
  if(info/=0)error stop 157
  call bind_dg_wpw_hs_callbacks(context,info)
  if(info/=0) error stop 16
  call snapshot_dg_wpw_production_context(context,frozen_context,info)
  if(info/=0)error stop 158
  saved_context_value=context%ww_projector_nonlocal(1,1)
  if(rank==0)context%ww_projector_nonlocal(1,1)=saved_context_value+(0.375d0,0d0)
  call validate_dg_wpw_production_context_snapshot(context,frozen_context,info)
  if(info==0)error stop 159
  if(rank==0)context%ww_projector_nonlocal(1,1)=saved_context_value
  saved_context_value=context%ww_projector_cross_value(1)
  if(rank==1)context%ww_projector_cross_value(1)=saved_context_value+(0.5d0,0d0)
  call validate_dg_wpw_production_context_snapshot(context,frozen_context,info)
  if(info==0)error stop 160
  if(rank==1)context%ww_projector_cross_value(1)=saved_context_value
  saved_context_value=context%wp_h_nonlocal(1)
  if(rank==0)context%wp_h_nonlocal(1)=saved_context_value+(0.625d0,0d0)
  call validate_dg_wpw_production_context_snapshot(context,frozen_context,info)
  if(info==0)error stop 161
  if(rank==0)context%wp_h_nonlocal(1)=saved_context_value
  if(rank==1)context%callbacks_bound=.false.
  call validate_dg_wpw_production_context_snapshot(context,frozen_context,info)
  if(info==0)error stop 162
  if(rank==1)context%callbacks_bound=.true.
  call validate_dg_wpw_production_context_snapshot(context,frozen_context,info)
  if(info/=0)error stop 163
  x=[(1d0,0.2d0),(-0.3d0,0.4d0),(0.7d0,-0.1d0),(-0.2d0,0.5d0)]
  call context%apply_h(x,yh_local,info); if(info/=0) error stop 17
  call context%apply_s(x,ys_local,info); if(info/=0) error stop 18
  call MPI_Allreduce(yh_local,yh,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(ys_local,ys,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  h_oracle=(0d0,0d0); s_oracle=(0d0,0d0)
  h_oracle(1,1)=(7.5d0,0d0);h_oracle(2,2)=(7.5d0,0d0)
  h_oracle(1,2)=(0.25d0,0d0);h_oracle(2,1)=(0.25d0,0d0)
  s_oracle(1,1)=(1d0,0d0);s_oracle(2,2)=(1d0,0d0)
  h_oracle(1:2,3:4)=wp; h_oracle(3:4,1:2)=transpose(conjg(wp)); h_oracle(3:4,3:4)=pp
  s_oracle(1:2,3:4)=swp; s_oracle(3:4,1:2)=transpose(conjg(swp)); s_oracle(3:4,3:4)=spp
  if(maxval(abs(yh-matmul(h_oracle,x)))>1d-12) error stop 19
  if(maxval(abs(ys-matmul(s_oracle,x)))>1d-12) error stop 20
  bxw(1,1)=x(rank+1);bxp(1,1)=x(2+owned_p)
  call apply_h_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
    context%bounded_operator%layout_fingerprint,bxw,bxp,byw,byp,info)
  if(info/=0)error stop 21
  bounded_local=0;bounded_local(rank+1)=byw(1,1);bounded_local(2+owned_p)=byp(1,1)
  call MPI_Allreduce(bounded_local,bounded,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(bounded-matmul(h_oracle,x)))>1d-12)error stop 22
  h0_oracle=h_oracle;h0_oracle(1,1)=(3.5d0,0d0);h0_oracle(2,2)=(3.5d0,0d0)
  h0_oracle(1:2,3:4)=h0_oracle(1:2,3:4)-face_delta
  h0_oracle(3:4,1:2)=transpose(conjg(h0_oracle(1:2,3:4)))
  interface_oracle=h_oracle-h0_oracle
  saved_owned_w=context%bounded_operator%owned_w_ids
  saved_required_w=context%bounded_operator%required_w_ids
  saved_owned_p=context%bounded_operator%owned_p_ids
  saved_required_p=context%bounded_operator%required_p_ids
  call set_dg_wpw_interface_lambda(context%bounded_operator,0d0,info)
  if(info/=0)error stop 221
  call apply_h_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
    context%bounded_operator%layout_fingerprint,bxw,bxp,byw,byp,info)
  if(info/=0)error stop 222
  bounded_local=0;bounded_local(rank+1)=byw(1,1);bounded_local(2+owned_p)=byp(1,1)
  call MPI_Allreduce(bounded_local,bounded,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(bounded-matmul(h0_oracle,x)))>1d-12)error stop 223
  call set_dg_wpw_interface_lambda(context%bounded_operator,0.5d0,info)
  if(info/=0)error stop 224
  call apply_h_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
    context%bounded_operator%layout_fingerprint,bxw,bxp,byw,byp,info)
  if(info/=0)error stop 225
  bounded_local=0;bounded_local(rank+1)=byw(1,1);bounded_local(2+owned_p)=byp(1,1)
  call MPI_Allreduce(bounded_local,bounded,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(bounded-matmul(h0_oracle+0.5d0*interface_oracle,x)))>1d-12)error stop 226
  if(any(context%bounded_operator%owned_w_ids/=saved_owned_w).or.&
     any(context%bounded_operator%required_w_ids/=saved_required_w).or.&
     any(context%bounded_operator%owned_p_ids/=saved_owned_p).or.&
     any(context%bounded_operator%required_p_ids/=saved_required_p))error stop 227
  saved_h0_cache=context%bounded_operator%ww_h0_dense(1,1)
  saved_interface_cache=context%bounded_operator%ww_interface_dense(1,1)
  saved_trial_h=context%bounded_operator%ww_h_dense(1,1)
  saved_transport_id=context%bounded_operator%owned_w_ids(1)
  call set_dg_wpw_interface_lambda(context%bounded_operator,1.25d0,info)
  if(info==0.or.context%bounded_operator%interface_lambda/=0.5d0.or.&
    context%bounded_operator%ww_h0_dense(1,1)/=saved_h0_cache.or.&
    context%bounded_operator%ww_interface_dense(1,1)/=saved_interface_cache.or.&
    context%bounded_operator%ww_h_dense(1,1)/=saved_trial_h.or.&
    context%bounded_operator%owned_w_ids(1)/=saved_transport_id)error stop 2271
  call set_dg_wpw_interface_lambda(context%bounded_operator,1d0,info)
  if(info/=0)error stop 228
  call apply_s_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
    context%bounded_operator%layout_fingerprint,bxw,bxp,bsw,bsp,info)
  if(info/=0)error stop 23
  bounded_local=0;bounded_local(rank+1)=bsw(1,1);bounded_local(2+owned_p)=bsp(1,1)
  call MPI_Allreduce(bounded_local,bounded,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(bounded-matmul(s_oracle,x)))>1d-12)error stop 24
  old_fingerprint=context%bounded_operator%layout_fingerprint;old_ww_potential=context%ww_potential(1,1)
  call import_dg_wpw_lcfo_ww_components(rank+1,[1,1],reshape([1,2],[1,2]),reshape([1d0],[1,1]),&
    reshape([5d0],[1,1]),reshape([3d0],[1,1]),reshape([4d0],[1,1]),[201],[99],[rank+1],&
    [1],[1],reshape([0,0,0],[3,1]),[0.1d0],'orthonormal_ww',ww_components,info)
  if(info/=0)error stop 241
  call publish_dg_wpw_lcfo_ww_components(context,ww_components,info)
  if(info/=0.or..not.context%bounded_operator%valid.or..not.context%callbacks_bound)error stop 242
  call build_dg_wpw_production_operator(context,info)
  if(info==0.or..not.context%bounded_operator%valid.or..not.context%callbacks_bound.or.&
    context%bounded_operator%layout_fingerprint/=old_fingerprint.or.&
    context%ww_potential(1,1)/=old_ww_potential)error stop 243
  call import_dg_wpw_lcfo_ww_components(rank+1,[1,1],reshape([1,2],[1,2]),reshape([1d0],[1,1]),&
    reshape([5d0],[1,1]),reshape([3d0],[1,1]),reshape([4d0],[1,1]),[integer::],[integer::],&
    [integer::],[integer::],[integer::],reshape([integer::],[3,0]),[real(8)::],'orthonormal_ww',ww_components,info)
  if(info/=0)error stop 244
  call publish_dg_wpw_lcfo_ww_components(context,ww_components,info)
  if(info/=0.or..not.context%bounded_operator%valid)error stop 245
  call build_dg_wpw_production_operator(context,info);if(info/=0.or.context%ww_potential(1,1)/=5d0)error stop 246
  old_fingerprint=context%bounded_operator%layout_fingerprint
  saved_face=context%wp_h_face;new_wp_volume=context%wp_h_volume;new_pp_volume=context%pp_h_volume
  new_wp_volume(1)=new_wp_volume(1)+(0.01d0,0d0)
  call replace_dg_wpw_potential_volume(context,new_wp_volume,new_pp_volume,info)
  if(info/=0.or.maxval(abs(context%wp_h_face-saved_face))>1d-14.or.&
    maxval(abs(context%wp_h-context%wp_h_volume-context%wp_h_face))>1d-14)error stop 247
  saved_total=context%wp_h
  call replace_dg_wpw_potential_volume(context,new_wp_volume(1:1),new_pp_volume,info)
  if(info==0.or.maxval(abs(context%wp_h-saved_total))>1d-14)error stop 248
  call build_dg_wpw_production_operator(context,info)
  if(info/=0.or.context%bounded_operator%layout_fingerprint==old_fingerprint)error stop 25
  if(rank==0) print '(a)','PASS two-rank rank-local production operator matches dense oracle'
  call release_dg_wpw_bounded_operator_snapshot(frozen_operator)
  call release_dg_wpw_production_context_snapshot(frozen_context)
  call release_dg_wpw_bounded_operator(context%bounded_operator)
  call MPI_Finalize(ierr)
end program
