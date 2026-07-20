program test_dg_wpw_core_wannier_seed_mpi
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_positive_inf
  use mpi
  use lcfo_wannier_sawf_seed,only:canonicalize_sawf_wannier_center,&
    transform_sawf_wannier_occupation,build_sawf_projected_wannier,&
    build_sawf_projected_wannier_from_overlap,&
    apply_sawf_projected_wannier_transform,&
    canonicalize_sawf_bond_identity,build_sawf_wannier_density,&
    qualify_sawf_wannier_density_projection
  use dg_wpw_wannier_tail_halo,only:validate_sawf_wannier_tail_schedule,&
    exchange_sawf_wannier_tail_values,exchange_sawf_discovered_wannier_tails
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
  real(8),parameter::tol=1d-10
  real(8)::lower(3),upper(3),wrapped(3)
  complex(8)::c(2,2),t(2,2),f_source(2,2),f_q(2,2),q(2,2),rho_source(2,2),rho_q(2,2)
  complex(8)::occupied(4,3),occupied_rotated(4,3),trial(4,2),wannier(4,2),wannier_rotated(4,2)
  complex(8)::projection_overlap(3,2),wannier_from_overlap(4,2),polar_transform(2,2)
  complex(8)::sample_values(2,3),transformed_samples(2,2)
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
  logical::owned

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(ierr/=MPI_SUCCESS)error stop 10
  if(rank>1)error stop 11
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
  if(is_sawf_outer_buffer_shell([1,12,24],[24,24,24],[24,24,24]))error stop 17
  if(.not.is_sawf_outer_buffer_shell([1,12,18],[18,18,18],[24,24,24]))error stop 18

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
    send_point,send_value,discovered_source,discovered_point,discovered_value,info)
  if(info/=0.or.size(discovered_point)/=1.or.discovered_point(1)/=recv_point(1).or.&
    any(discovered_source(:,1)/=recv_source(:,1)).or.abs(discovered_value(1)-&
      merge(cmplx(-1d0,0.5d0,8),cmplx(2d0,3d0,8),rank==0))>1d-14)&
    error stop 604

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
