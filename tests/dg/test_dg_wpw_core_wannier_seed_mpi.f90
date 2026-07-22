program test_dg_wpw_core_wannier_seed_mpi
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_positive_inf
  use mpi
  use lcfo_wannier_sawf_seed,only:canonicalize_sawf_wannier_center,&
    transform_sawf_wannier_occupation,build_sawf_projected_wannier,&
    build_sawf_projected_wannier_from_overlap,&
    apply_sawf_projected_wannier_transform,&
    apply_sawf_projected_wannier_gradient_transform,&
    canonicalize_sawf_bond_identity,build_sawf_wannier_density,&
    qualify_sawf_wannier_density_projection,diagnose_sawf_discrete_wannier_spread,&
    assemble_sawf_diagonal_periodic_links,normalize_sawf_projected_wannier_columns
  use dg_wpw_wannier_tail_halo,only:validate_sawf_wannier_tail_schedule,&
    exchange_sawf_wannier_tail_values,exchange_sawf_discovered_wannier_tails
  use dg_wpw_wannier_tail_halo,only:sawf_tail_records_unique
  use dg_wpw_wannier_tail_halo,only:locate_sawf_wannier_tail_core
  use dg_wpw_wannier_tail_halo,only:locate_sawf_wannier_tail_rank
  use dg_wpw_wannier_tail_halo,only:qualify_sawf_wannier_buffer_tail
  use dg_wpw_wannier_tail_halo,only:is_sawf_outer_buffer_shell
  implicit none
  integer::ierr,rank,info,owner_count,image(3),destination_fragment,destination_local(3),destination_rank
  integer::fragment_start(3,2),fragment_shape(3,2)
  integer::rank_fragment(2),rank_orbital_lane(2),rank_grid_lo(3,2),rank_grid_hi(3,2)
  integer::bond_atoms(2),bond_image(3)
  integer,allocatable::send_source(:,:),send_destination(:),send_point(:),recv_source(:,:),recv_point(:)
  integer,allocatable::send_image(:,:),discovered_image(:,:)
  integer,allocatable::many_source(:,:),many_destination(:),many_point(:),many_image(:,:)
  real(8),parameter::tol=1d-10
  real(8)::lower(3),upper(3),wrapped(3)
  complex(8)::c(2,2),t(2,2),f_source(2,2),f_q(2,2),q(2,2),rho_source(2,2),rho_q(2,2)
  complex(8)::occupied(4,3),occupied_rotated(4,3),trial(4,2),wannier(4,2),wannier_rotated(4,2)
  complex(8)::projection_overlap(3,2),wannier_from_overlap(4,2),polar_transform(2,2)
  complex(8)::sample_values(2,3),transformed_samples(2,2)
  complex(8)::occupied_gradients(3,3,3),gradient_partial(3,3,2),gradient_total(3,3,2)
  complex(8)::expected_gradient(3,2),derivative(3,3)
  integer::ix
  complex(8)::rotation(3,3),rank_bad_trial(4,2),rank_bad_wannier(4,2)
  complex(8),allocatable::send_value(:),recv_value(:)
  integer,allocatable::discovered_source(:,:),discovered_point(:)
  complex(8),allocatable::discovered_value(:)
  complex(8)::source_on_core,w_value,p_value,overlap_local(2),overlap_full(2),&
    overlap_missing_tail(2),overlap_missing_w(2),overlap_missing_p(2)
  complex(8)::point_values(2,2),density_occupation(2,2),bad_occupation(2,2)
  real(8)::density(2),charge
  real(8)::projection_residual,normalization_residual,projection_charge_error
  real(8)::outer_shell_ratio,omitted_tail_ratio
  real(8)::condition,condition_rotated,condition_bad,theta
  complex(8)::spread_link(2,6)
  real(8)::spread_norm(2),spread_bvec(3,6),spread_weight(6),spread_center(3,2),spread_omega(2)
  real(8)::spread_center_expected(3),spread_sigma,spread_expected
  logical::spread_center_valid(2)
  complex(8)::link_values(4,2),link_local(2,6),link_global(2,6),link_expected(2,6)
  real(8)::link_coordinates(3,4),link_norm_local(2),link_norm_global(2),link_norm_expected(2)
  real(8)::column_norm(2)
  complex(8)::normalization_polar(2,2),normalization_core(3,2),normalization_p(4,2),&
    normalization_polar_before(2,2),normalization_core_before(3,2),normalization_p_before(4,2)
  integer::ip,ib
  logical::owned

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(ierr/=MPI_SUCCESS)error stop 10
  if(rank>1)error stop 11
  normalization_polar=reshape([(cmplx(dble(info),-0.1d0*dble(info),8),info=1,4)],[2,2])
  normalization_core=reshape([(cmplx(0.2d0*dble(info),0.05d0*dble(info),8),info=1,6)],[3,2])
  normalization_p=reshape([(cmplx(-0.1d0*dble(info),0.03d0*dble(info),8),info=1,8)],[4,2])
  normalization_polar_before=normalization_polar;normalization_core_before=normalization_core
  normalization_p_before=normalization_p;column_norm=[4d0,9d0]
  call normalize_sawf_projected_wannier_columns(column_norm,normalization_polar,&
    normalization_core,normalization_p,info)
  if(info/=0.or.maxval(abs(normalization_polar(:,1)-0.5d0*normalization_polar_before(:,1)))>tol.or.&
    maxval(abs(normalization_polar(:,2)-normalization_polar_before(:,2)/3d0))>tol.or.&
    maxval(abs(normalization_core(:,1)-0.5d0*normalization_core_before(:,1)))>tol.or.&
    maxval(abs(normalization_core(:,2)-normalization_core_before(:,2)/3d0))>tol.or.&
    maxval(abs(normalization_p(:,1)-0.5d0*normalization_p_before(:,1)))>tol.or.&
    maxval(abs(normalization_p(:,2)-normalization_p_before(:,2)/3d0))>tol)error stop 123
  normalization_polar_before=normalization_polar;normalization_core_before=normalization_core
  normalization_p_before=normalization_p;column_norm=1d0
  call normalize_sawf_projected_wannier_columns(column_norm,normalization_polar,&
    normalization_core,normalization_p,info)
  if(info/=0.or.maxval(abs(normalization_polar-normalization_polar_before))>0d0.or.&
    maxval(abs(normalization_core-normalization_core_before))>0d0.or.&
    maxval(abs(normalization_p-normalization_p_before))>0d0)error stop 124
  column_norm=[1d0,0d0]
  call normalize_sawf_projected_wannier_columns(column_norm,normalization_polar,&
    normalization_core,normalization_p,info)
  if(info==0.or.maxval(abs(normalization_polar-normalization_polar_before))>0d0.or.&
    maxval(abs(normalization_core-normalization_core_before))>0d0.or.&
    maxval(abs(normalization_p-normalization_p_before))>0d0)error stop 125
  column_norm=[1d0,ieee_value(0d0,ieee_positive_inf)]
  call normalize_sawf_projected_wannier_columns(column_norm,normalization_polar,&
    normalization_core,normalization_p,info)
  if(info==0.or.maxval(abs(normalization_polar-normalization_polar_before))>0d0)error stop 126
  column_norm=1d0
  call normalize_sawf_projected_wannier_columns(column_norm,normalization_polar(:,1:1),&
    normalization_core,normalization_p,info)
  if(info==0.or.maxval(abs(normalization_core-normalization_core_before))>0d0)error stop 127
  link_values=reshape([(cmplx(0.1d0*dble(info),-0.03d0*dble(info),8),info=1,8)],[4,2])
  link_coordinates=reshape([0d0,0d0,0d0, 0.4d0,0.2d0,0.1d0, 0.9d0,0.3d0,0.7d0, &
    1.4d0,0.8d0,1.1d0],[3,4])
  spread_bvec=0d0
  spread_bvec(1,1:2)=[1d0,-1d0]
  spread_bvec(2,3:4)=[1d0,-1d0]
  spread_bvec(3,5:6)=[1d0,-1d0]
  spread_weight=0.5d0
  link_expected=(0d0,0d0);link_norm_expected=0d0
  do ib=1,6;do ip=1,4
    link_expected(:,ib)=link_expected(:,ib)+abs(link_values(ip,:))**2*&
      exp(cmplx(0d0,-dot_product(spread_bvec(:,ib),link_coordinates(:,ip)),8))*0.25d0
    link_norm_expected=link_norm_expected+merge(abs(link_values(ip,:))**2*0.25d0,0d0,ib==1)
  enddo;enddo
  link_local=(0d0,0d0);link_norm_local=0d0
  if(rank==0)call assemble_sawf_diagonal_periodic_links(link_values,link_coordinates,spread_bvec,&
    0.25d0,link_norm_local,link_local,info)
  if(rank==0.and.info/=0)error stop 121
  call MPI_Allreduce(link_local,link_global,size(link_global),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(link_norm_local,link_norm_global,size(link_norm_global),MPI_DOUBLE_PRECISION,&
    MPI_SUM,MPI_COMM_WORLD,ierr)
  if(ierr/=MPI_SUCCESS.or.maxval(abs(link_global-link_expected))>1d-14.or.&
    maxval(abs(link_norm_global-link_norm_expected))>1d-14)error stop 122
  spread_center_expected=[0.2d0,-0.15d0,0.1d0]
  spread_sigma=0.3d0
  do ix=1,6
    spread_link(1,ix)=exp(cmplx(-0.5d0*spread_sigma**2*sum(spread_bvec(:,ix)**2),&
      -dot_product(spread_bvec(:,ix),spread_center_expected),8))
  enddo
  spread_norm=[1d0,4d0]
  spread_link(2,:)=4d0*spread_link(1,:)
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  spread_expected=3d0*(1d0-exp(-spread_sigma**2))
  if(info/=0.or..not.all(spread_center_valid).or.&
    maxval(abs(spread_center-spread(spread_center_expected,2,2)))>1d-12.or.&
    maxval(abs(spread_omega-spread_expected))>1d-12)error stop 110
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info,require_unit_norm=.true.)
  if(info==0)error stop 114
  spread_norm=1d0;spread_link(2,:)=spread_link(1,:)
  spread_center_expected=[acos(-1d0)-0.1d0,0d0,0d0]
  do ix=1,6
    spread_link(:,ix)=exp(cmplx(0d0,-dot_product(spread_bvec(:,ix),spread_center_expected),8))
  enddo
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info,require_unit_norm=.true.)
  if(info/=0.or..not.all(spread_center_valid).or.&
    maxval(abs(spread_center-spread(spread_center_expected,2,2)))>1d-12.or.&
    maxval(abs(spread_omega))>1d-12)error stop 115
  spread_link=cmplx(0d0,0d0,8)
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info/=0.or.any(spread_center_valid).or.any(spread_omega<huge(1d0)))error stop 116
  spread_link=cmplx(1d0,0d0,8)
  spread_link(1,1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info==0)error stop 117
  spread_link=cmplx(1d0,0d0,8);spread_weight(1)=0.4d0
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info==0)error stop 118
  spread_weight=0.5d0
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec(:,1:5),spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info==0)error stop 119
  spread_link(1,:)=2d0;spread_link(2,:)=2d0
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info==0)error stop 120
  spread_center_expected=[0.2d0,-0.15d0,0.1d0]
  do ix=1,6
    spread_link(:,ix)=exp(cmplx(-0.5d0*spread_sigma**2*sum(spread_bvec(:,ix)**2),&
      -dot_product(spread_bvec(:,ix),spread_center_expected),8))
  enddo
  spread_link(1,2)=spread_link(1,1)
  call diagnose_sawf_discrete_wannier_spread(spread_link,spread_norm,spread_bvec,spread_weight,&
    spread_center,spread_omega,spread_center_valid,info)
  if(info==0)error stop 113
  allocate(many_source(5,100000),many_destination(100000),many_point(100000),many_image(3,100000))
  many_source=0;many_source(1,:)=[(info,info=1,100000)]
  many_destination=1-rank;many_point=[(info,info=1,100000)];many_image=0
  if(.not.sawf_tail_records_unique(many_source,many_destination,many_point,many_image))error stop 111
  many_image(:,100000)=many_image(:,99999);many_source(:,100000)=many_source(:,99999)
  many_point(100000)=many_point(99999)
  if(sawf_tail_records_unique(many_source,many_destination,many_point,many_image))error stop 112
  deallocate(many_source,many_destination,many_point,many_image)
  lower=[0.5d0*dble(rank),0d0,0d0]
  upper=[0.5d0*dble(rank+1),1d0,1d0]
  fragment_start=reshape([1,1,1,5,1,1],[3,2]);fragment_shape=4
  call locate_sawf_wannier_tail_core([8,2,3],[8,4,4],fragment_start,fragment_shape,&
    destination_fragment,destination_local,info)
  if(info/=0.or.destination_fragment/=2.or.any(destination_local/=[4,2,3]))error stop 12
  call locate_sawf_wannier_tail_core([9,2,3],[8,4,4],fragment_start,fragment_shape,&
    destination_fragment,destination_local,info)
  if(info/=0.or.destination_fragment/=1.or.any(destination_local/=[1,2,3]))error stop 13
  rank_fragment=[1,2];rank_orbital_lane=0;rank_grid_lo=1;rank_grid_hi=4
  call locate_sawf_wannier_tail_rank(2,[4,2,3],rank_fragment,rank_orbital_lane,&
    rank_grid_lo,rank_grid_hi,destination_rank,info)
  if(info/=0.or.destination_rank/=1)error stop 14
  call qualify_sawf_wannier_buffer_tail(2d0,1d-14,2d-14,1d-14,1d-10,&
    outer_shell_ratio,omitted_tail_ratio,info)
  if(info/=0.or.abs(outer_shell_ratio-1d-14)>1d-20.or.&
    abs(omitted_tail_ratio-5d-15)>1d-20)error stop 15
  call qualify_sawf_wannier_buffer_tail(2d0,0d0,4d-8,0d0,1d-10,&
    outer_shell_ratio,omitted_tail_ratio,info)
  if(info==0.or.outer_shell_ratio<=1d-10)error stop 16
  if(.not.is_sawf_outer_buffer_shell([1,12,24],[24,24,24],[24,24,24]))error stop 17
  if(.not.is_sawf_outer_buffer_shell([1,12,18],[18,18,18],[24,24,24]))error stop 18
  if(.not.is_sawf_outer_buffer_shell([32,12,12],[32,32,32],[24,24,24]))error stop 181

  call canonicalize_sawf_wannier_center([0.5d0,0.25d0,0.75d0],[1d0,1d0,1d0],&
    lower,upper,tol,wrapped,image,owned,info)
  if(info/=0.or.any(abs(wrapped-[0.5d0,0.25d0,0.75d0])>tol).or.any(image/=0))error stop 20
  call MPI_Allreduce(merge(1,0,owned),owner_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(ierr/=MPI_SUCCESS.or.owner_count/=1.or.(rank==0.and.owned).or.(rank==1.and..not.owned))error stop 21

  call canonicalize_sawf_wannier_center([1d0,0.25d0,0.75d0],[1d0,1d0,1d0],&
    lower,upper,tol,wrapped,image,owned,info)
  if(info/=0.or.any(abs(wrapped-[0d0,0.25d0,0.75d0])>tol).or.any(image/=[1,0,0]))error stop 22
  call MPI_Allreduce(merge(1,0,owned),owner_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(ierr/=MPI_SUCCESS.or.owner_count/=1.or.(rank==0.and..not.owned).or.(rank==1.and.owned))error stop 23

  call canonicalize_sawf_wannier_center([2.25d0,0.25d0,0.75d0],[1d0,1d0,1d0],&
    lower,upper,tol,wrapped,image,owned,info)
  if(info/=0.or.any(abs(wrapped-[0.25d0,0.25d0,0.75d0])>tol).or.any(image/=[2,0,0]))error stop 24

  call canonicalize_sawf_wannier_center([0.5d0-2d0*tol,0.25d0,0.75d0],[1d0,1d0,1d0],&
    lower,upper,tol,wrapped,image,owned,info)
  if(info/=0.or.(rank==0.and..not.owned).or.(rank==1.and.owned))error stop 25
  call canonicalize_sawf_wannier_center([0.5d0+2d0*tol,0.25d0,0.75d0],[1d0,1d0,1d0],&
    lower,upper,tol,wrapped,image,owned,info)
  if(info/=0.or.(rank==0.and.owned).or.(rank==1.and..not.owned))error stop 26

  c=reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0.5d0,0d0,8),&
    cmplx(sqrt(0.75d0),0d0,8)],[2,2])
  t=reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(-0.5d0/sqrt(0.75d0),0d0,8),cmplx(1d0/sqrt(0.75d0),0d0,8)],[2,2])
  f_source=0;f_source(1,1)=2;f_source(2,2)=2
  call transform_sawf_wannier_occupation(t,f_source,f_q,info)
  if(info/=0.or.maxval(abs(f_q-conjg(transpose(f_q))))>1d-12)error stop 30
  q=matmul(c,t)
  rho_source=matmul(c,matmul(f_source,conjg(transpose(c))))
  rho_q=matmul(q,matmul(f_q,conjg(transpose(q))))
  if(maxval(abs(rho_source-rho_q))>1d-12)error stop 31
  if(maxval(abs(rho_source-2d0*matmul(q,conjg(transpose(q)))))<1d-3)error stop 32

  occupied=0;occupied(1,1)=1;occupied(2,2)=1;occupied(3,3)=1
  trial=0;trial(:,1)=[cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0.2d0,0d0,8),cmplx(0.4d0,0d0,8)]
  trial(:,2)=[cmplx(0d0,0d0,8),cmplx(1d0,0d0,8),cmplx(0.1d0,0d0,8),cmplx(-0.3d0,0d0,8)]
  call build_sawf_projected_wannier(occupied,trial,1d0,1d-12,wannier,condition,info)
  if(info/=0.or.condition<1d0.or.maxval(abs(matmul(conjg(transpose(wannier)),wannier)-&
    reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),&
      cmplx(1d0,0d0,8)],[2,2])))>1d-12)error stop 40
  projection_overlap=matmul(conjg(transpose(occupied)),trial)
  call build_sawf_projected_wannier_from_overlap(occupied,projection_overlap,1d0,1d-12,&
    wannier_from_overlap,polar_transform,condition_bad,info)
  if(info/=0.or.abs(condition-condition_bad)>1d-12.or.&
    maxval(abs(wannier-wannier_from_overlap))>1d-12.or.&
    maxval(abs(matmul(matmul(occupied,projection_overlap),polar_transform)-wannier))>1d-12)error stop 401
  sample_values=reshape([(cmplx(dble(info),-0.25d0*dble(info),8),info=1,6)],[2,3])
  call apply_sawf_projected_wannier_transform(sample_values,projection_overlap,polar_transform,&
    transformed_samples,info)
  if(info/=0.or.maxval(abs(transformed_samples-&
    matmul(matmul(sample_values,projection_overlap),polar_transform)))>1d-12)error stop 402
  transformed_samples=cmplx(7d0,-3d0,8)
  call apply_sawf_projected_wannier_transform(sample_values(:,1:2),projection_overlap,&
    polar_transform,transformed_samples,info)
  if(info==0.or.maxval(abs(transformed_samples))>0d0)error stop 403
  occupied_gradients=0
  derivative=0;derivative(:,1)=[cmplx(1d0,0d0,8),cmplx(2d0,0d0,8),cmplx(3d0,0d0,8)]
  derivative(:,2)=[cmplx(-2d0,0d0,8),cmplx(0.5d0,0d0,8),cmplx(1d0,0d0,8)]
  do ix=1,3
    occupied_gradients(ix,:,rank+1)=spread(derivative(ix,rank+1),1,3)
  enddo
  call apply_sawf_projected_wannier_gradient_transform(occupied_gradients,projection_overlap,&
    polar_transform,gradient_partial,info)
  if(info/=0)error stop 404
  call MPI_Allreduce(gradient_partial,gradient_total,size(gradient_total),MPI_DOUBLE_COMPLEX,&
    MPI_SUM,MPI_COMM_WORLD,ierr)
  do ix=1,3
    expected_gradient(ix,:)=matmul(derivative(ix,:),matmul(projection_overlap,polar_transform))
    if(maxval(abs(gradient_total(ix,:,:)-spread(expected_gradient(ix,:),1,3)))>1d-12)error stop 405
  enddo
  theta=0.37d0;rotation=0;rotation(1,1)=cos(theta);rotation(2,1)=sin(theta)
  rotation(1,2)=-sin(theta);rotation(2,2)=cos(theta);rotation(3,3)=1
  occupied_rotated=matmul(occupied,rotation)
  call build_sawf_projected_wannier(occupied_rotated,trial,1d0,1d-12,&
    wannier_rotated,condition_rotated,info)
  if(info/=0.or.abs(condition-condition_rotated)>1d-12.or.&
    maxval(abs(wannier-wannier_rotated))>1d-12)error stop 41
  rank_bad_trial(:,1)=trial(:,1);rank_bad_trial(:,2)=trial(:,1)
  call build_sawf_projected_wannier(occupied,rank_bad_trial,1d0,1d-12,&
    rank_bad_wannier,condition_bad,info)
  if(info==0)error stop 42

  call canonicalize_sawf_bond_identity(9,3,[1,-2,0],bond_atoms,bond_image,info)
  if(info/=0.or.any(bond_atoms/=[3,9]).or.any(bond_image/=[-1,2,0]))error stop 50
  call canonicalize_sawf_bond_identity(3,9,[-1,2,0],bond_atoms,bond_image,info)
  if(info/=0.or.any(bond_atoms/=[3,9]).or.any(bond_image/=[-1,2,0]))error stop 51
  call canonicalize_sawf_bond_identity(3,3,[0,0,0],bond_atoms,bond_image,info)
  if(info==0)error stop 52

  allocate(send_source(5,1),send_destination(1),send_point(1),recv_source(5,1),recv_point(1))
  if(rank==0)then
    send_source(:,1)=[1,2,1,0,0];send_destination=1;send_point=7
    recv_source(:,1)=[3,4,-1,0,0];recv_point=9
  else
    send_source(:,1)=[3,4,-1,0,0];send_destination=0;send_point=9
    recv_source(:,1)=[1,2,1,0,0];recv_point=7
  endif
  call validate_sawf_wannier_tail_schedule(MPI_COMM_WORLD,send_source,send_destination,send_point,&
    recv_source,recv_point,info)
  if(info/=0)error stop 60
  allocate(send_value(1),recv_value(1))
  send_value(1)=merge(cmplx(2d0,3d0,8),cmplx(-1d0,0.5d0,8),rank==0)
  call exchange_sawf_wannier_tail_values(MPI_COMM_WORLD,send_source,send_destination,send_point,&
    send_value,recv_source,recv_point,recv_value,info)
  if(info/=0.or.abs(recv_value(1)-merge(cmplx(-1d0,0.5d0,8),&
    cmplx(2d0,3d0,8),rank==0))>1d-14)error stop 601

  ! A basis row on either core needs the neighboring source tail.  Both W and
  ! P blocks must therefore be accumulated after the keyed halo exchange.
  source_on_core=send_value(1)+recv_value(1)
  w_value=cmplx(1d0+rank,0.25d0,8);p_value=cmplx(-0.5d0,1d0+rank,8)
  overlap_local=[conjg(w_value)*source_on_core,conjg(p_value)*source_on_core]
  call MPI_Allreduce(overlap_local,overlap_full,2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  overlap_local=[conjg(w_value)*send_value(1),conjg(p_value)*send_value(1)]
  call MPI_Allreduce(overlap_local,overlap_missing_tail,2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  overlap_local=[cmplx(0d0,0d0,8),conjg(p_value)*source_on_core]
  call MPI_Allreduce(overlap_local,overlap_missing_w,2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  overlap_local=[conjg(w_value)*source_on_core,cmplx(0d0,0d0,8)]
  call MPI_Allreduce(overlap_local,overlap_missing_p,2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(ierr/=MPI_SUCCESS.or.any(abs(overlap_full-overlap_missing_tail)<1d-12).or.&
    abs(overlap_missing_w(1))>1d-14.or.abs(overlap_missing_p(2))>1d-14.or.&
    abs(overlap_full(1))<1d-12.or.abs(overlap_full(2))<1d-12)error stop 602
  send_value(1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call exchange_sawf_wannier_tail_values(MPI_COMM_WORLD,send_source,send_destination,send_point,&
    send_value,recv_source,recv_point,recv_value,info)
  if(info==0)error stop 603
  send_value(1)=merge(cmplx(2d0,3d0,8),cmplx(-1d0,0.5d0,8),rank==0)
  call exchange_sawf_discovered_wannier_tails(MPI_COMM_WORLD,send_source,send_destination,&
    send_point,reshape([0,0,0],[3,1]),send_value,discovered_source,discovered_point,&
    discovered_image,discovered_value,info)
  if(info/=0.or.size(discovered_point)/=1.or.discovered_point(1)/=recv_point(1).or.&
    any(discovered_source(:,1)/=recv_source(:,1)).or.abs(discovered_value(1)-&
      merge(cmplx(-1d0,0.5d0,8),cmplx(2d0,3d0,8),rank==0))>1d-14)&
    error stop 604
  if(any(discovered_image(:,1)/=0))error stop 605

  if(rank==1)recv_source(3,1)=2
  call validate_sawf_wannier_tail_schedule(MPI_COMM_WORLD,send_source,send_destination,send_point,&
    recv_source,recv_point,info)
  if(info==0)error stop 61
  if(rank==1)recv_source(3,1)=1

  deallocate(send_source,send_destination,send_point)
  if(rank==0)then
    allocate(send_source(5,2),send_destination(2),send_point(2))
    send_source(:,1)=[1,2,1,0,0];send_source(:,2)=send_source(:,1)
    send_destination=1;send_point=7
  else
    allocate(send_source(5,1),send_destination(1),send_point(1))
    send_source(:,1)=[3,4,-1,0,0];send_destination=0;send_point=9
  endif
  call validate_sawf_wannier_tail_schedule(MPI_COMM_WORLD,send_source,send_destination,send_point,&
    recv_source,recv_point,info)
  if(info==0)error stop 62

  deallocate(send_source,send_destination,send_point,send_value)
  if(rank==0)then
    allocate(send_source(5,2),send_destination(2),send_point(2),send_image(3,2),send_value(2))
    send_source(:,1)=[1,2,1,0,0];send_source(:,2)=send_source(:,1)
    send_destination=1;send_point=7
    send_image(:,1)=[0,0,0];send_image(:,2)=[1,0,0]
    send_value=[cmplx(2d0,0d0,8),cmplx(3d0,0d0,8)]
  else
    allocate(send_source(5,0),send_destination(0),send_point(0),send_image(3,0),send_value(0))
  endif
  call exchange_sawf_discovered_wannier_tails(MPI_COMM_WORLD,send_source,send_destination,&
    send_point,send_image,send_value,discovered_source,discovered_point,discovered_image,&
    discovered_value,info)
  if(info/=0)error stop 63
  if(rank==1)then
    if(size(discovered_point)/=2.or.any(discovered_point/=7).or.&
      abs(sum(discovered_value)-cmplx(5d0,0d0,8))>1d-14)error stop 64
  else if(size(discovered_point)/=0)then
    error stop 65
  endif
  if(rank==0)send_image(:,2)=send_image(:,1)
  call exchange_sawf_discovered_wannier_tails(MPI_COMM_WORLD,send_source,send_destination,&
    send_point,send_image,send_value,discovered_source,discovered_point,discovered_image,&
    discovered_value,info)
  if(info==0)error stop 66

  point_values=reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0.5d0,0d0,8),cmplx(0d0,1d0,8)],[2,2])
  density_occupation=reshape([cmplx(2d0,0d0,8),cmplx(0.25d0,0d0,8),&
    cmplx(0.25d0,0d0,8),cmplx(1d0,0d0,8)],[2,2])
  call build_sawf_wannier_density(point_values,density_occupation,[0.4d0,0.6d0],density,charge,info)
  if(info/=0.or.maxval(abs(density-[2.5d0,1d0]))>1d-14.or.abs(charge-1.6d0)>1d-14)error stop 70
  bad_occupation=density_occupation;bad_occupation(1,2)=bad_occupation(1,2)+(0d0,0.1d0)
  call build_sawf_wannier_density(point_values,bad_occupation,[0.4d0,0.6d0],density,charge,info)
  if(info==0)error stop 71
  call qualify_sawf_wannier_density_projection([2.5d0,1d0],[2.5d0,1d0],[2.5d0,1d0],&
    [0.4d0,0.6d0],1.6d0,1d0,1d-12,projection_residual,normalization_residual,&
    projection_charge_error,info)
  if(info/=0.or.maxval(abs([projection_residual,normalization_residual,&
    projection_charge_error]))>1d-14)error stop 72
  call qualify_sawf_wannier_density_projection([2.5d0,1d0],[2d0,1d0],[2d0,1d0],&
    [0.4d0,0.6d0],1.6d0,0.99d0,1d-3,projection_residual,normalization_residual,&
    projection_charge_error,info)
  if(info==0.or.projection_residual<=1d-3.or.projection_charge_error<=1d-3)error stop 73

  if(rank==0)write(*,'(a)')'PASS DG WPW core Wannier seed MPI oracle'
  call MPI_Finalize(ierr)
end program
