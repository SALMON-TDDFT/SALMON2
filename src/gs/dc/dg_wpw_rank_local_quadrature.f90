module dg_wpw_rank_local_quadrature
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_wpw_quadrature_assembler,only:assemble_wpw_volume_point
  use rt_dg_wpw_weak_form_evaluator,only:wpw_volume_weak_pair
  implicit none
  private
  type,public::s_dg_wpw_volume_accumulator
    logical::valid=.false.,failed=.false.
    integer::npoint=0,point_capacity=0
    complex(8),allocatable::wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:)
    integer,allocatable::grid_ids(:)
    real(8),allocatable::potentials(:),weights(:),densities(:)
    complex(8),allocatable::w_points(:,:),p_points(:,:)
    complex(8),allocatable::grad_w_points(:,:,:),grad_p_points(:,:,:)
  end type
  public::accumulate_dg_wpw_core_volume
  public::initialize_dg_wpw_volume_accumulator,add_dg_wpw_core_point
  public::finalize_dg_wpw_volume_accumulator
  public::build_dg_wpw_rank_local_quadrature
  interface build_dg_wpw_rank_local_quadrature
    module procedure build_dg_wpw_rank_local_quadrature_batch
    module procedure build_dg_wpw_rank_local_quadrature_accumulator
  end interface
contains
  subroutine build_dg_wpw_rank_local_quadrature_batch(comm,root,w,grad_w,p_owned,grad_p_owned,p_support,&
      grad_p_support,potential,weight,wp_h,wp_s,pp_h,pp_s,info)
    use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Reduce,MPI_INTEGER,MPI_MAX,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_SUCCESS
    integer,intent(in)::comm,root
    complex(8),intent(in)::w(:,:),grad_w(:,:,:),p_owned(:,:),grad_p_owned(:,:,:)
    complex(8),intent(in)::p_support(:,:),grad_p_support(:,:,:)
    real(8),intent(in)::potential(:),weight(:)
    complex(8),intent(out)::wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:)
    integer,intent(out)::info
    complex(8),allocatable::local_wp_h(:,:),local_wp_s(:,:),local_pp_h(:,:),local_pp_s(:,:)
    integer::rank,nrank,ierr,local_info,local_bad,global_bad

    wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);info=1
    call MPI_Comm_rank(comm,rank,ierr);local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank)local_bad=1
    allocate(local_wp_h,source=wp_h);allocate(local_wp_s,source=wp_s)
    allocate(local_pp_h,source=pp_h);allocate(local_pp_s,source=pp_s)
    call accumulate_dg_wpw_core_volume(w,grad_w,p_owned,grad_p_owned,p_support,grad_p_support,&
      potential,weight,local_wp_h,local_wp_s,local_pp_h,local_pp_s,local_info)
    if(local_info/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    local_bad=0
    call MPI_Reduce(local_wp_h,wp_h,size(wp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(local_wp_s,wp_s,size(wp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(local_pp_h,pp_h,size(pp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(local_pp_s,pp_s,size(pp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);return
    endif
    info=0
  end subroutine build_dg_wpw_rank_local_quadrature_batch

  subroutine build_dg_wpw_rank_local_quadrature_accumulator(comm,root,accumulator,wp_h,wp_s,pp_h,pp_s,info)
    use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Reduce,MPI_INTEGER,MPI_MAX,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_SUCCESS
    integer,intent(in)::comm,root
    type(s_dg_wpw_volume_accumulator),intent(in)::accumulator
    complex(8),intent(out)::wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,global_npoint

    wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);info=1
    call MPI_Comm_rank(comm,rank,ierr);local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank.or..not.accumulator%valid.or.accumulator%failed.or.accumulator%npoint<0)&
      local_bad=1
    if(accumulator%valid)then
      if(any(shape(wp_h)/=shape(accumulator%wp_h)).or.any(shape(wp_s)/=shape(accumulator%wp_s)).or.&
         any(shape(pp_h)/=shape(accumulator%pp_h)).or.any(shape(pp_s)/=shape(accumulator%pp_s)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(accumulator%npoint,global_npoint,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_npoint<=0)return
    local_bad=0
    call MPI_Reduce(accumulator%wp_h,wp_h,size(wp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(accumulator%wp_s,wp_s,size(wp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(accumulator%pp_h,pp_h,size(pp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Reduce(accumulator%pp_s,pp_s,size(pp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);return
    endif
    info=0
  end subroutine build_dg_wpw_rank_local_quadrature_accumulator

  subroutine initialize_dg_wpw_volume_accumulator(accumulator,nw,npo,nps,info,point_capacity)
    type(s_dg_wpw_volume_accumulator),intent(out)::accumulator
    integer,intent(in)::nw,npo,nps
    integer,intent(out)::info
    integer,intent(in),optional::point_capacity
    integer::capacity
    info=1
    capacity=0;if(present(point_capacity))capacity=point_capacity
    if(nw<=0.or.npo<=0.or.nps<=0.or.capacity<=0)return
    allocate(accumulator%wp_h(nw,npo),accumulator%wp_s(nw,npo),&
      accumulator%pp_h(npo,nps),accumulator%pp_s(npo,nps))
    allocate(accumulator%grid_ids(capacity),accumulator%potentials(capacity),accumulator%weights(capacity),&
      accumulator%densities(capacity))
    allocate(accumulator%w_points(nw,capacity),accumulator%p_points(nps,capacity))
    allocate(accumulator%grad_w_points(3,nw,capacity),accumulator%grad_p_points(3,nps,capacity))
    accumulator%wp_h=(0d0,0d0);accumulator%wp_s=(0d0,0d0)
    accumulator%pp_h=(0d0,0d0);accumulator%pp_s=(0d0,0d0)
    accumulator%point_capacity=capacity
    accumulator%valid=.true.;accumulator%failed=.false.;accumulator%npoint=0;info=0
  end subroutine

  subroutine add_dg_wpw_core_point(accumulator,w,grad_w,p_owned,grad_p_owned,p_support,&
      grad_p_support,potential,weight,info,grid_id,density)
    type(s_dg_wpw_volume_accumulator),intent(inout)::accumulator
    complex(8),intent(in)::w(:),grad_w(:,:),p_owned(:),grad_p_owned(:,:),p_support(:),grad_p_support(:,:)
    real(8),intent(in)::potential,weight
    integer,intent(out)::info
    integer,intent(in),optional::grid_id
    real(8),intent(in),optional::density
    complex(8)::overlap,hvalue
    integer::iw,ip,jp,pair_info,next_point,point_id
    real(8)::point_density
    info=1
    point_density=0d0;if(present(density))point_density=density
    if(.not.accumulator%valid.or.accumulator%failed)return
    if(.not.ieee_is_finite(point_density))then;accumulator%failed=.true.;return;endif
    if(any(shape(grad_w)/=[3,size(w)]).or.any(shape(grad_p_owned)/=[3,size(p_owned)]).or.&
       any(shape(grad_p_support)/=[3,size(p_support)]).or.&
       any(shape(accumulator%wp_h)/=[size(w),size(p_owned)]).or.&
       any(shape(accumulator%wp_s)/=shape(accumulator%wp_h)).or.&
       any(shape(accumulator%pp_h)/=[size(p_owned),size(p_support)]).or.&
       any(shape(accumulator%pp_s)/=shape(accumulator%pp_h)))then
      accumulator%failed=.true.;return
    endif
    ! Validate the complete point before publishing any partial contribution.
    do ip=1,size(p_owned)
      do iw=1,size(w)
        call wpw_volume_weak_pair(w(iw),grad_w(:,iw),p_owned(ip),grad_p_owned(:,ip),&
          potential,weight,overlap,hvalue,pair_info)
        if(pair_info/=0)then;accumulator%failed=.true.;return;endif
      enddo
      do jp=1,size(p_support)
        call wpw_volume_weak_pair(p_owned(ip),grad_p_owned(:,ip),p_support(jp),grad_p_support(:,jp),&
          potential,weight,overlap,hvalue,pair_info)
        if(pair_info/=0)then;accumulator%failed=.true.;return;endif
      enddo
    enddo
    next_point=accumulator%npoint+1;point_id=0;if(present(grid_id))point_id=grid_id
    if(next_point>accumulator%point_capacity)then;accumulator%failed=.true.;return;endif
    accumulator%grid_ids(next_point)=point_id;accumulator%potentials(next_point)=potential
    accumulator%weights(next_point)=weight;accumulator%densities(next_point)=point_density
    accumulator%w_points(:,next_point)=w
    accumulator%p_points(:,next_point)=p_support;accumulator%grad_w_points(:,:,next_point)=grad_w
    accumulator%grad_p_points(:,:,next_point)=grad_p_support
    do ip=1,size(p_owned)
      do iw=1,size(w)
        call wpw_volume_weak_pair(w(iw),grad_w(:,iw),p_owned(ip),grad_p_owned(:,ip),&
          potential,weight,overlap,hvalue,pair_info)
        accumulator%wp_s(iw,ip)=accumulator%wp_s(iw,ip)+overlap
        accumulator%wp_h(iw,ip)=accumulator%wp_h(iw,ip)+hvalue
      enddo
      do jp=1,size(p_support)
        call wpw_volume_weak_pair(p_owned(ip),grad_p_owned(:,ip),p_support(jp),grad_p_support(:,jp),&
          potential,weight,overlap,hvalue,pair_info)
        accumulator%pp_s(ip,jp)=accumulator%pp_s(ip,jp)+overlap
        accumulator%pp_h(ip,jp)=accumulator%pp_h(ip,jp)+hvalue
      enddo
    enddo
    accumulator%npoint=next_point
    info=0
  end subroutine

  subroutine finalize_dg_wpw_volume_accumulator(accumulator,wp_h,wp_s,pp_h,pp_s,info)
    type(s_dg_wpw_volume_accumulator),intent(in)::accumulator
    complex(8),intent(out)::wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:)
    integer,intent(out)::info
    wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);info=1
    if(.not.accumulator%valid.or.accumulator%failed.or.accumulator%npoint<=0)return
    if(any(shape(wp_h)/=shape(accumulator%wp_h)).or.any(shape(wp_s)/=shape(accumulator%wp_s)).or.&
       any(shape(pp_h)/=shape(accumulator%pp_h)).or.any(shape(pp_s)/=shape(accumulator%pp_s)))return
    wp_h=accumulator%wp_h;wp_s=accumulator%wp_s;pp_h=accumulator%pp_h;pp_s=accumulator%pp_s;info=0
  end subroutine

  subroutine accumulate_dg_wpw_core_volume(w,grad_w,p_owned,grad_p_owned,p_support,grad_p_support,&
      potential,weight,wp_h,wp_s,pp_h,pp_s,info)
    complex(8),intent(in)::w(:,:),grad_w(:,:,:),p_owned(:,:),grad_p_owned(:,:,:)
    complex(8),intent(in)::p_support(:,:),grad_p_support(:,:,:)
    real(8),intent(in)::potential(:),weight(:)
    complex(8),intent(out)::wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:)
    integer,intent(out)::info
    integer::ipoint,npoint,point_info
    complex(8),allocatable::point_wp_h(:,:),point_wp_s(:,:),point_pp_h(:,:),point_pp_s(:,:)
    complex(8),allocatable::sum_wp_h(:,:),sum_wp_s(:,:),sum_pp_h(:,:),sum_pp_s(:,:)

    wp_h=(0d0,0d0);wp_s=(0d0,0d0);pp_h=(0d0,0d0);pp_s=(0d0,0d0);info=1
    npoint=size(potential)
    if(npoint<=0.or.size(weight)/=npoint.or.size(w,2)/=npoint.or.size(p_owned,2)/=npoint.or.&
       size(p_support,2)/=npoint.or.any(shape(grad_w)/=[3,size(w,1),npoint]).or.&
       any(shape(grad_p_owned)/=[3,size(p_owned,1),npoint]).or.&
       any(shape(grad_p_support)/=[3,size(p_support,1),npoint]).or.&
       any(shape(wp_h)/=[size(w,1),size(p_owned,1)]).or.any(shape(wp_s)/=shape(wp_h)).or.&
       any(shape(pp_h)/=[size(p_owned,1),size(p_support,1)]).or.any(shape(pp_s)/=shape(pp_h)))return
    allocate(point_wp_h,source=wp_h);allocate(point_wp_s,source=wp_s)
    allocate(point_pp_h,source=pp_h);allocate(point_pp_s,source=pp_s)
    allocate(sum_wp_h,source=wp_h);allocate(sum_wp_s,source=wp_s)
    allocate(sum_pp_h,source=pp_h);allocate(sum_pp_s,source=pp_s)
    do ipoint=1,npoint
      call assemble_wpw_volume_point(w(:,ipoint),grad_w(:,:,ipoint),p_owned(:,ipoint),&
        grad_p_owned(:,:,ipoint),p_support(:,ipoint),grad_p_support(:,:,ipoint),potential(ipoint),&
        weight(ipoint),point_wp_h,point_wp_s,point_pp_h,point_pp_s,point_info)
      if(point_info/=0)then;info=2;return;endif
      sum_wp_h=sum_wp_h+point_wp_h;sum_wp_s=sum_wp_s+point_wp_s
      sum_pp_h=sum_pp_h+point_pp_h;sum_pp_s=sum_pp_s+point_pp_s
    enddo
    wp_h=sum_wp_h;wp_s=sum_wp_s;pp_h=sum_pp_h;pp_s=sum_pp_s;info=0
  end subroutine
end module dg_wpw_rank_local_quadrature
