program test_dg_wpw_rank_local_quadrature
  use dg_wpw_rank_local_quadrature, only: accumulate_dg_wpw_core_volume,s_dg_wpw_volume_accumulator,&
    initialize_dg_wpw_volume_accumulator,add_dg_wpw_core_point,finalize_dg_wpw_volume_accumulator
  implicit none
  integer,parameter::nw=1,npo=1,nps=2,np=2
  complex(8)::w(nw,np),gw(3,nw,np),po(npo,np),gpo(3,npo,np),ps(nps,np),gps(3,nps,np)
  complex(8)::wp_h(nw,npo),wp_s(nw,npo),pp_h(npo,nps),pp_s(npo,nps)
  complex(8)::wp_h_ref(nw,npo),wp_s_ref(nw,npo),pp_h_ref(npo,nps),pp_s_ref(npo,nps)
  real(8)::potential(np),weight(np)
  integer::info,ipoint,j
  type(s_dg_wpw_volume_accumulator)::accumulator

  w(1,:)=[(1d0,0.2d0),(0.5d0,-0.1d0)];gw=(0d0,0d0)
  gw(1,1,:)=[(0.2d0,0.1d0),(-0.1d0,0.05d0)]
  po(1,:)=[(0.7d0,-0.2d0),(0.4d0,0.3d0)];gpo=(0d0,0d0)
  gpo(1,1,:)=[(0.3d0,0d0),(0.1d0,-0.2d0)]
  ps(1,:)=po(1,:);ps(2,:)=[(-0.2d0,0.4d0),(0.6d0,0.1d0)];gps=(0d0,0d0)
  gps(1,1,:)=gpo(1,1,:);gps(2,2,:)=[(0.25d0,0d0),(-0.15d0,0.1d0)]
  potential=[-0.4d0,0.7d0];weight=[0.25d0,0.5d0]

  call accumulate_dg_wpw_core_volume(w,gw,po,gpo,ps,gps,potential,weight,wp_h,wp_s,pp_h,pp_s,info)
  if(info/=0)error stop 1
  wp_h_ref=(0d0,0d0);wp_s_ref=(0d0,0d0);pp_h_ref=(0d0,0d0);pp_s_ref=(0d0,0d0)
  do ipoint=1,np
    wp_s_ref(1,1)=wp_s_ref(1,1)+weight(ipoint)*conjg(w(1,ipoint))*po(1,ipoint)
    wp_h_ref(1,1)=wp_h_ref(1,1)+weight(ipoint)*(0.5d0*sum(conjg(gw(:,1,ipoint))*gpo(:,1,ipoint))+&
      potential(ipoint)*conjg(w(1,ipoint))*po(1,ipoint))
    do j=1,nps
      pp_s_ref(1,j)=pp_s_ref(1,j)+weight(ipoint)*conjg(po(1,ipoint))*ps(j,ipoint)
      pp_h_ref(1,j)=pp_h_ref(1,j)+weight(ipoint)*(0.5d0*sum(conjg(gpo(:,1,ipoint))*gps(:,j,ipoint))+&
        potential(ipoint)*conjg(po(1,ipoint))*ps(j,ipoint))
    enddo
  enddo
  if(maxval(abs(wp_h-wp_h_ref))>1d-13.or.maxval(abs(wp_s-wp_s_ref))>1d-13)error stop 2
  if(maxval(abs(pp_h-pp_h_ref))>1d-13.or.maxval(abs(pp_s-pp_s_ref))>1d-13)error stop 3
  call initialize_dg_wpw_volume_accumulator(accumulator,nw,npo,nps,info,point_capacity=np)
  if(info/=0)error stop 7
  do ipoint=1,np
    call add_dg_wpw_core_point(accumulator,w(:,ipoint),gw(:,:,ipoint),po(:,ipoint),gpo(:,:,ipoint),&
      ps(:,ipoint),gps(:,:,ipoint),potential(ipoint),weight(ipoint),info,grid_id=100+ipoint,&
      density=merge(0.6d0,0.9d0,ipoint==1))
    if(info/=0)error stop 8
  enddo
  call finalize_dg_wpw_volume_accumulator(accumulator,wp_h,wp_s,pp_h,pp_s,info)
  if(info/=0.or.maxval(abs(wp_h-wp_h_ref))>1d-13.or.maxval(abs(pp_h-pp_h_ref))>1d-13)error stop 9
  if(any(accumulator%grid_ids/=[101,102]).or.maxval(abs(accumulator%w_points-w))>1d-14.or.&
    maxval(abs(accumulator%p_points-ps))>1d-14.or.maxval(abs(accumulator%grad_w_points-gw))>1d-14.or.&
    maxval(abs(accumulator%grad_p_points-gps))>1d-14.or.any(accumulator%weights/=weight).or.&
    any(accumulator%densities/=[0.6d0,0.9d0]))error stop 12
  call accumulate_dg_wpw_core_volume(w(:,1:1),gw(:,:,1:1),po,gpo,ps,gps,potential,weight,&
    wp_h,wp_s,pp_h,pp_s,info)
  if(info==0)error stop 4
  call initialize_dg_wpw_volume_accumulator(accumulator,nw,npo,nps,info,point_capacity=1)
  call add_dg_wpw_core_point(accumulator,w(:,1),gw(:,:,1),po(:,1),gpo(:,:,1),ps(:,1),gps(:,:,1),&
    potential(1),-1d0,info)
  if(info==0)error stop 10
  call finalize_dg_wpw_volume_accumulator(accumulator,wp_h,wp_s,pp_h,pp_s,info)
  if(info==0.or.any(abs(wp_h)>0d0).or.any(abs(pp_h)>0d0))error stop 11
  print '(a)','PASS uniquely-owned core WPW volume quadrature'
end program test_dg_wpw_rank_local_quadrature
