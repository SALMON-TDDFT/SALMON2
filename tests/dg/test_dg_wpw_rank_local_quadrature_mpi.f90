program test_dg_wpw_rank_local_quadrature_mpi
  use mpi
  use dg_wpw_rank_local_quadrature,only:build_dg_wpw_rank_local_quadrature,&
    s_dg_wpw_volume_accumulator,initialize_dg_wpw_volume_accumulator,add_dg_wpw_core_point
  implicit none
  complex(8)::w(1,1),gw(3,1,1),po(1,1),gpo(3,1,1),ps(2,1),gps(3,2,1)
  complex(8)::wp_h(1,1),wp_s(1,1),pp_h(1,2),pp_s(1,2)
  complex(8)::ref_wp_h(1,1),ref_wp_s(1,1),ref_pp_h(1,2),ref_pp_s(1,2)
  real(8)::potential(1),weight(1)
  integer::ierr,rank,nrank,info
  type(s_dg_wpw_volume_accumulator)::accumulator

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2)error stop 1
  w(1,1)=merge((1d0,0.2d0),(0.5d0,-0.1d0),rank==0);gw=(0d0,0d0)
  gw(1,1,1)=merge((0.2d0,0.1d0),(-0.1d0,0.05d0),rank==0)
  po(1,1)=merge((0.7d0,-0.2d0),(0.4d0,0.3d0),rank==0);gpo=(0d0,0d0)
  gpo(1,1,1)=merge((0.3d0,0d0),(0.1d0,-0.2d0),rank==0)
  ps(1,1)=po(1,1);ps(2,1)=merge((-0.2d0,0.4d0),(0.6d0,0.1d0),rank==0);gps=(0d0,0d0)
  gps(1,1,1)=gpo(1,1,1);gps(2,2,1)=merge((0.25d0,0d0),(-0.15d0,0.1d0),rank==0)
  potential(1)=merge(-0.4d0,0.7d0,rank==0);weight(1)=merge(0.25d0,0.5d0,rank==0)
  call build_dg_wpw_rank_local_quadrature(MPI_COMM_WORLD,0,w,gw,po,gpo,ps,gps,potential,weight,&
    wp_h,wp_s,pp_h,pp_s,info)
  if(info/=0)error stop 2
  if(rank==0)then
    call dense_reference(ref_wp_h,ref_wp_s,ref_pp_h,ref_pp_s)
    if(maxval(abs(wp_h-ref_wp_h))>1d-13.or.maxval(abs(wp_s-ref_wp_s))>1d-13)error stop 3
    if(maxval(abs(pp_h-ref_pp_h))>1d-13.or.maxval(abs(pp_s-ref_pp_s))>1d-13)error stop 4
  else
    if(any(abs(wp_h)>0d0).or.any(abs(pp_h)>0d0))error stop 5
  endif
  call initialize_dg_wpw_volume_accumulator(accumulator,1,1,2,info,point_capacity=1)
  call add_dg_wpw_core_point(accumulator,w(:,1),gw(:,:,1),po(:,1),gpo(:,:,1),ps(:,1),gps(:,:,1),&
    potential(1),weight(1),info)
  call build_dg_wpw_rank_local_quadrature(MPI_COMM_WORLD,0,accumulator,wp_h,wp_s,pp_h,pp_s,info)
  if(info/=0)error stop 7
  if(rank==0)then
    call dense_reference(ref_wp_h,ref_wp_s,ref_pp_h,ref_pp_s)
    if(maxval(abs(wp_h-ref_wp_h))>1d-13.or.maxval(abs(pp_h-ref_pp_h))>1d-13)error stop 8
  endif
  if(rank==1)weight(1)=-1d0
  call build_dg_wpw_rank_local_quadrature(MPI_COMM_WORLD,0,w,gw,po,gpo,ps,gps,potential,weight,&
    wp_h,wp_s,pp_h,pp_s,info)
  if(info==0.or.any(abs(wp_h)>0d0).or.any(abs(wp_s)>0d0).or.any(abs(pp_h)>0d0).or.any(abs(pp_s)>0d0))&
    error stop 6
  if(rank==0)print '(a)','PASS fragment-rank WPW quadrature reduction is transactional'
  call MPI_Finalize(ierr)
contains
  subroutine dense_reference(wh,ws,ph,psout)
    complex(8),intent(out)::wh(1,1),ws(1,1),ph(1,2),psout(1,2)
    complex(8)::ww(2),gww(2),ppo(2),gppo(2),pps(2,2),gpps(2,2)
    real(8)::v(2),wt(2);integer::i,j
    ww=[(1d0,0.2d0),(0.5d0,-0.1d0)];gww=[(0.2d0,0.1d0),(-0.1d0,0.05d0)]
    ppo=[(0.7d0,-0.2d0),(0.4d0,0.3d0)];gppo=[(0.3d0,0d0),(0.1d0,-0.2d0)]
    pps(1,:)=ppo;pps(2,:)=[(-0.2d0,0.4d0),(0.6d0,0.1d0)]
    gpps(1,:)=gppo;gpps(2,:)=[(0.25d0,0d0),(-0.15d0,0.1d0)]
    v=[-0.4d0,0.7d0];wt=[0.25d0,0.5d0];wh=0;ws=0;ph=0;psout=0
    do i=1,2
      ws(1,1)=ws(1,1)+wt(i)*conjg(ww(i))*ppo(i)
      wh(1,1)=wh(1,1)+wt(i)*(0.5d0*conjg(gww(i))*gppo(i)+v(i)*conjg(ww(i))*ppo(i))
      do j=1,2
        psout(1,j)=psout(1,j)+wt(i)*conjg(ppo(i))*pps(j,i)
        ph(1,j)=ph(1,j)+wt(i)*(0.5d0*merge(conjg(gppo(i))*gpps(j,i),(0d0,0d0),j==1)+&
          v(i)*conjg(ppo(i))*pps(j,i))
      enddo
    enddo
  end subroutine
end program
