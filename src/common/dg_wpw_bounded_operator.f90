module dg_wpw_bounded_operator
  use mpi,only:MPI_Allreduce,MPI_INTEGER,MPI_MAX,MPI_MIN,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_SUCCESS
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_wpw_owner_exchange,only:s_dg_wpw_owner_schedule,initialize_dg_wpw_owner_schedule,&
    fetch_rows_from_owners,reduce_w_partial_to_owners,release_dg_wpw_owner_schedule
  implicit none
  private
  type,public::s_dg_wpw_bounded_operator
    integer::operator_epoch=-1
    integer(8)::layout_fingerprint=0
    real(8)::interface_lambda=1d0
    character(32)::ownership_kind='',metric_convention='',operator_convention=''
    logical::valid=.false.
    integer,allocatable::owned_w_ids(:),owned_p_ids(:),required_w_ids(:),required_p_ids(:)
    integer,allocatable::ww_r(:),ww_c(:),wp_w(:),wp_p(:),pp_r(:),pp_c(:)
    integer,allocatable::ww_ri(:),ww_ci(:),wp_wi(:),wp_pi(:),pp_ri(:),pp_ci(:)
    complex(8),allocatable::ww_h(:),ww_s(:),wp_h(:),wp_s(:),pp_h(:),pp_s(:)
    complex(8),allocatable::ww_h_dense(:,:),ww_s_dense(:,:),wp_h_dense(:,:),wp_s_dense(:,:),&
      pp_h_dense(:,:),pp_s_dense(:,:)
    complex(8),allocatable::ww_h0_dense(:,:),ww_interface_dense(:,:),&
      wp_h0_dense(:,:),wp_interface_dense(:,:),pp_h0_dense(:,:),pp_interface_dense(:,:)
    type(s_dg_wpw_owner_schedule)::w_schedule,p_schedule
  end type
  public::initialize_dg_wpw_bounded_operator,apply_h_dg_wpw_bounded,apply_s_dg_wpw_bounded
  public::global_gram_dg_wpw_bounded
  public::fetch_dg_wpw_support_coefficients
  public::release_dg_wpw_bounded_operator
  public::set_dg_wpw_interface_lambda
contains
  subroutine initialize_dg_wpw_bounded_operator(op,comm,epoch,fingerprint,ownership,metric,convention,&
      peers,owned_w,owned_p,required_w,required_p,ww_r,ww_c,ww_h,ww_s,wp_w,wp_p,wp_h,wp_s,&
      pp_r,pp_c,pp_h,pp_s,info,ww_h0,ww_interface,wp_h0,wp_interface,pp_h0,pp_interface)
    type(s_dg_wpw_bounded_operator),intent(inout)::op
    integer,intent(in)::comm,epoch,peers(:),owned_w(:),owned_p(:),required_w(:),required_p(:)
    integer(8),intent(in)::fingerprint
    character(*),intent(in)::ownership,metric,convention
    integer,intent(in)::ww_r(:),ww_c(:),wp_w(:),wp_p(:),pp_r(:),pp_c(:)
    complex(8),intent(in)::ww_h(:),ww_s(:),wp_h(:),wp_s(:),pp_h(:),pp_s(:)
    integer,intent(out)::info
    complex(8),intent(in),optional::ww_h0(:),ww_interface(:),wp_h0(:),wp_interface(:),&
      pp_h0(:),pp_interface(:)
    type(s_dg_wpw_bounded_operator)::candidate
    integer::i,local_bad,global_bad,ierr,astat
    info=1;local_bad=0
    if(epoch<=0.or.fingerprint==0.or.len_trim(ownership)==0.or.len_trim(metric)==0.or.&
       len_trim(convention)==0.or.size(ww_r)/=size(ww_c).or.size(ww_r)/=size(ww_h).or.&
       size(ww_r)/=size(ww_s).or.size(wp_w)/=size(wp_p).or.size(wp_w)/=size(wp_h).or.&
       size(wp_w)/=size(wp_s).or.size(pp_r)/=size(pp_c).or.size(pp_r)/=size(pp_h).or.&
       size(pp_r)/=size(pp_s).or..not.finite1(ww_h).or..not.finite1(ww_s).or.&
       .not.finite1(wp_h).or..not.finite1(wp_s).or..not.finite1(pp_h).or..not.finite1(pp_s).or.&
       .not.strictly_increasing(owned_w).or..not.strictly_increasing(owned_p).or.&
       .not.strictly_increasing(required_w).or..not.strictly_increasing(required_p))local_bad=1
    call initialize_dg_wpw_owner_schedule(candidate%w_schedule,comm,peers,owned_w,required_w,info)
    if(info/=0)local_bad=1
    call initialize_dg_wpw_owner_schedule(candidate%p_schedule,comm,peers,owned_p,required_p,info)
    if(info/=0)local_bad=1
    allocate(candidate%ww_ri(size(ww_r)),candidate%ww_ci(size(ww_c)),&
      candidate%wp_wi(size(wp_w)),candidate%wp_pi(size(wp_p)),&
      candidate%pp_ri(size(pp_r)),candidate%pp_ci(size(pp_c)))
    astat=0
    allocate(candidate%ww_h_dense(size(required_w),size(required_w)),&
      candidate%ww_s_dense(size(required_w),size(required_w)),&
      candidate%wp_h_dense(size(required_w),size(owned_p)),&
      candidate%wp_s_dense(size(required_w),size(owned_p)),&
      candidate%pp_h_dense(size(owned_p),size(required_p)),&
      candidate%pp_s_dense(size(owned_p),size(required_p)),&
      candidate%ww_h0_dense(size(required_w),size(required_w)),&
      candidate%ww_interface_dense(size(required_w),size(required_w)),&
      candidate%wp_h0_dense(size(required_w),size(owned_p)),&
      candidate%wp_interface_dense(size(required_w),size(owned_p)),&
      candidate%pp_h0_dense(size(owned_p),size(required_p)),&
      candidate%pp_interface_dense(size(owned_p),size(required_p)),stat=astat)
    if(astat/=0)local_bad=1
    candidate%ww_ri=0;candidate%ww_ci=0;candidate%wp_wi=0;candidate%wp_pi=0
    candidate%pp_ri=0;candidate%pp_ci=0
    if(astat==0)then
      candidate%ww_h_dense=0;candidate%ww_s_dense=0
      candidate%wp_h_dense=0;candidate%wp_s_dense=0
      candidate%pp_h_dense=0;candidate%pp_s_dense=0
      candidate%ww_h0_dense=0;candidate%ww_interface_dense=0
      candidate%wp_h0_dense=0;candidate%wp_interface_dense=0
      candidate%pp_h0_dense=0;candidate%pp_interface_dense=0
    endif
    if(present(ww_h0).neqv.present(ww_interface).or.present(wp_h0).neqv.present(wp_interface).or.&
       present(pp_h0).neqv.present(pp_interface).or.present(ww_h0).neqv.present(wp_h0).or.&
       present(ww_h0).neqv.present(pp_h0))local_bad=1
    if(present(ww_h0))then
      if(size(ww_h0)/=size(ww_h).or.size(ww_interface)/=size(ww_h).or.&
         size(wp_h0)/=size(wp_h).or.size(wp_interface)/=size(wp_h).or.&
         size(pp_h0)/=size(pp_h).or.size(pp_interface)/=size(pp_h).or.&
         .not.finite1(ww_h0).or..not.finite1(ww_interface).or..not.finite1(wp_h0).or.&
         .not.finite1(wp_interface).or..not.finite1(pp_h0).or..not.finite1(pp_interface))local_bad=1
    endif
    if(local_bad==0)then
      do i=1,size(ww_r)
        candidate%ww_ri(i)=find_id_sorted(required_w,ww_r(i))
        candidate%ww_ci(i)=find_id_sorted(required_w,ww_c(i))
        if(candidate%ww_ri(i)<=0.or.candidate%ww_ci(i)<=0)local_bad=1
      enddo
      do i=1,size(wp_w)
        candidate%wp_wi(i)=find_id_sorted(required_w,wp_w(i))
        candidate%wp_pi(i)=find_id_sorted(owned_p,wp_p(i))
        if(candidate%wp_wi(i)<=0.or.candidate%wp_pi(i)<=0)local_bad=1
      enddo
      do i=1,size(pp_r)
        candidate%pp_ri(i)=find_id_sorted(owned_p,pp_r(i))
        candidate%pp_ci(i)=find_id_sorted(required_p,pp_c(i))
        if(candidate%pp_ri(i)<=0.or.candidate%pp_ci(i)<=0)local_bad=1
      enddo
      if(local_bad==0)then
        do i=1,size(ww_r)
          candidate%ww_h_dense(candidate%ww_ri(i),candidate%ww_ci(i))=&
            candidate%ww_h_dense(candidate%ww_ri(i),candidate%ww_ci(i))+ww_h(i)
          candidate%ww_s_dense(candidate%ww_ri(i),candidate%ww_ci(i))=&
            candidate%ww_s_dense(candidate%ww_ri(i),candidate%ww_ci(i))+ww_s(i)
          if(present(ww_h0))then
            candidate%ww_h0_dense(candidate%ww_ri(i),candidate%ww_ci(i))=&
              candidate%ww_h0_dense(candidate%ww_ri(i),candidate%ww_ci(i))+ww_h0(i)
            candidate%ww_interface_dense(candidate%ww_ri(i),candidate%ww_ci(i))=&
              candidate%ww_interface_dense(candidate%ww_ri(i),candidate%ww_ci(i))+ww_interface(i)
          else
            candidate%ww_h0_dense(candidate%ww_ri(i),candidate%ww_ci(i))=&
              candidate%ww_h0_dense(candidate%ww_ri(i),candidate%ww_ci(i))+ww_h(i)
          endif
        enddo
        do i=1,size(wp_w)
          candidate%wp_h_dense(candidate%wp_wi(i),candidate%wp_pi(i))=&
            candidate%wp_h_dense(candidate%wp_wi(i),candidate%wp_pi(i))+wp_h(i)
          candidate%wp_s_dense(candidate%wp_wi(i),candidate%wp_pi(i))=&
            candidate%wp_s_dense(candidate%wp_wi(i),candidate%wp_pi(i))+wp_s(i)
          if(present(wp_h0))then
            candidate%wp_h0_dense(candidate%wp_wi(i),candidate%wp_pi(i))=&
              candidate%wp_h0_dense(candidate%wp_wi(i),candidate%wp_pi(i))+wp_h0(i)
            candidate%wp_interface_dense(candidate%wp_wi(i),candidate%wp_pi(i))=&
              candidate%wp_interface_dense(candidate%wp_wi(i),candidate%wp_pi(i))+wp_interface(i)
          else
            candidate%wp_h0_dense(candidate%wp_wi(i),candidate%wp_pi(i))=&
              candidate%wp_h0_dense(candidate%wp_wi(i),candidate%wp_pi(i))+wp_h(i)
          endif
        enddo
        do i=1,size(pp_r)
          candidate%pp_h_dense(candidate%pp_ri(i),candidate%pp_ci(i))=&
            candidate%pp_h_dense(candidate%pp_ri(i),candidate%pp_ci(i))+pp_h(i)
          candidate%pp_s_dense(candidate%pp_ri(i),candidate%pp_ci(i))=&
            candidate%pp_s_dense(candidate%pp_ri(i),candidate%pp_ci(i))+pp_s(i)
          if(present(pp_h0))then
            candidate%pp_h0_dense(candidate%pp_ri(i),candidate%pp_ci(i))=&
              candidate%pp_h0_dense(candidate%pp_ri(i),candidate%pp_ci(i))+pp_h0(i)
            candidate%pp_interface_dense(candidate%pp_ri(i),candidate%pp_ci(i))=&
              candidate%pp_interface_dense(candidate%pp_ri(i),candidate%pp_ci(i))+pp_interface(i)
          else
            candidate%pp_h0_dense(candidate%pp_ri(i),candidate%pp_ci(i))=&
              candidate%pp_h0_dense(candidate%pp_ri(i),candidate%pp_ci(i))+pp_h(i)
          endif
        enddo
        if(maxval(abs(candidate%ww_h_dense-candidate%ww_h0_dense-candidate%ww_interface_dense))>&
             1d-12*max(1d0,maxval(abs(candidate%ww_h_dense))).or.&
           maxval(abs(candidate%wp_h_dense-candidate%wp_h0_dense-candidate%wp_interface_dense))>&
             1d-12*max(1d0,maxval(abs(candidate%wp_h_dense))).or.&
           maxval(abs(candidate%pp_h_dense-candidate%pp_h0_dense-candidate%pp_interface_dense))>&
             1d-12*max(1d0,maxval(abs(candidate%pp_h_dense))))local_bad=1
        if(.not.finite2(candidate%ww_h_dense).or..not.finite2(candidate%ww_s_dense).or.&
           .not.finite2(candidate%wp_h_dense).or..not.finite2(candidate%wp_s_dense).or.&
           .not.finite2(candidate%pp_h_dense).or..not.finite2(candidate%pp_s_dense))local_bad=1
      endif
    endif
    info=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%w_schedule)
      call release_dg_wpw_owner_schedule(candidate%p_schedule)
      return
    endif
    candidate%operator_epoch=epoch;candidate%layout_fingerprint=fingerprint
    candidate%ownership_kind=ownership;candidate%metric_convention=metric
    candidate%operator_convention=convention
    candidate%owned_w_ids=owned_w;candidate%owned_p_ids=owned_p
    candidate%required_w_ids=required_w;candidate%required_p_ids=required_p
    candidate%ww_r=ww_r;candidate%ww_c=ww_c;candidate%ww_h=ww_h;candidate%ww_s=ww_s
    candidate%wp_w=wp_w;candidate%wp_p=wp_p;candidate%wp_h=wp_h;candidate%wp_s=wp_s
    candidate%pp_r=pp_r;candidate%pp_c=pp_c;candidate%pp_h=pp_h;candidate%pp_s=pp_s
    candidate%valid=.true.
    call release_dg_wpw_bounded_operator(op)
    op=candidate;info=0
  end subroutine

  subroutine set_dg_wpw_interface_lambda(op,lambda,info)
    type(s_dg_wpw_bounded_operator),intent(inout)::op
    real(8),intent(in)::lambda
    integer,intent(out)::info
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,op%valid.and.ieee_is_finite(lambda).and.lambda>=0d0.and.lambda<=1d0.and.&
      allocated(op%ww_h0_dense).and.allocated(op%ww_interface_dense).and.&
      allocated(op%wp_h0_dense).and.allocated(op%wp_interface_dense).and.&
      allocated(op%pp_h0_dense).and.allocated(op%pp_interface_dense))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    op%ww_h_dense=op%ww_h0_dense+lambda*op%ww_interface_dense
    op%wp_h_dense=op%wp_h0_dense+lambda*op%wp_interface_dense
    op%pp_h_dense=op%pp_h0_dense+lambda*op%pp_interface_dense
    op%interface_lambda=lambda
    info=0
  end subroutine set_dg_wpw_interface_lambda

  subroutine release_dg_wpw_bounded_operator(op)
    type(s_dg_wpw_bounded_operator),intent(inout)::op
    call release_dg_wpw_owner_schedule(op%w_schedule)
    call release_dg_wpw_owner_schedule(op%p_schedule)
    op=s_dg_wpw_bounded_operator()
  end subroutine

  subroutine apply_h_dg_wpw_bounded(op,expected_epoch,expected_fingerprint,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer,intent(in)::expected_epoch;integer(8),intent(in)::expected_fingerprint
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::info
    call apply_blocks(op,expected_epoch,expected_fingerprint,op%ww_h_dense,op%wp_h_dense,&
      op%pp_h_dense,xw,xp,yw,yp,info)
  end subroutine
  subroutine fetch_dg_wpw_support_coefficients(op,xw,xp,rw,rp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    complex(8),intent(in)::xw(:,:),xp(:,:)
    complex(8),intent(out)::rw(:,:),rp(:,:)
    integer,intent(out)::info
    integer::local_bad,global_bad,ierr
    rw=0;rp=0;info=1
    local_bad=0
    if(.not.op%valid.or.size(xw,1)/=size(op%owned_w_ids).or.size(xp,1)/=size(op%owned_p_ids).or.&
       size(xp,2)/=size(xw,2).or.any(shape(rw)/=[size(op%required_w_ids),size(xw,2)]).or.&
       any(shape(rp)/=[size(op%required_p_ids),size(xw,2)]).or..not.finite2(xw).or.&
       .not.finite2(xp))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call fetch_rows_from_owners(op%w_schedule,xw,rw,info);if(info/=0)return
    call fetch_rows_from_owners(op%p_schedule,xp,rp,info)
    if(info/=0)then;rw=0;rp=0;return;endif
    local_bad=merge(0,1,finite2(rw).and.finite2(rp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;rw=0;rp=0;info=1;return;endif
    info=0
  end subroutine
  subroutine apply_s_dg_wpw_bounded(op,expected_epoch,expected_fingerprint,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer,intent(in)::expected_epoch;integer(8),intent(in)::expected_fingerprint
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::info
    call apply_blocks(op,expected_epoch,expected_fingerprint,op%ww_s_dense,op%wp_s_dense,&
      op%pp_s_dense,xw,xp,yw,yp,info)
  end subroutine

  subroutine apply_blocks(op,expected_epoch,expected_fingerprint,ww_dense,wp_dense,pp_dense,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer,intent(in)::expected_epoch;integer(8),intent(in)::expected_fingerprint
    complex(8),intent(in)::ww_dense(:,:),wp_dense(:,:),pp_dense(:,:),xw(:,:),xp(:,:)
    complex(8),intent(out)::yw(:,:),yp(:,:);integer,intent(out)::info
    complex(8),allocatable::rw(:,:),rp(:,:),wpartial(:,:),wowned(:,:),candidate_yw(:,:),candidate_yp(:,:)
    integer::nvec,nvec_min,nvec_max,local_bad,global_bad,ierr
    yw=0;yp=0;info=1;local_bad=0;nvec=size(xw,2)
    if(.not.op%valid.or.expected_epoch/=op%operator_epoch.or.&
       expected_fingerprint/=op%layout_fingerprint.or.any(shape(xw)/=[size(op%owned_w_ids),nvec]).or.&
       any(shape(xp)/=[size(op%owned_p_ids),nvec]).or.any(shape(yw)/=shape(xw)).or.&
       any(shape(yp)/=shape(xp)).or.size(xp,2)/=nvec.or.&
       any(shape(ww_dense)/=[size(op%required_w_ids),size(op%required_w_ids)]).or.&
       any(shape(wp_dense)/=[size(op%required_w_ids),size(op%owned_p_ids)]).or.&
       any(shape(pp_dense)/=[size(op%owned_p_ids),size(op%required_p_ids)]).or.&
       .not.finite2(xw).or..not.finite2(xp))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(nvec,nvec_min,1,MPI_INTEGER,MPI_MIN,op%w_schedule%comm,ierr)
    call MPI_Allreduce(nvec,nvec_max,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nvec_min/=nvec_max)return
    allocate(rw(size(op%required_w_ids),nvec),rp(size(op%required_p_ids),nvec))
    allocate(wpartial(size(op%required_w_ids),nvec),wowned(size(op%owned_w_ids),nvec));wpartial=0
    allocate(candidate_yw(size(yw,1),size(yw,2)),candidate_yp(size(yp,1),size(yp,2)))
    candidate_yw=0;candidate_yp=0
    call fetch_rows_from_owners(op%w_schedule,xw,rw,info);if(info/=0)return
    call fetch_rows_from_owners(op%p_schedule,xp,rp,info);if(info/=0)return
    wpartial=matmul(ww_dense,rw)+matmul(wp_dense,xp)
    candidate_yp=matmul(conjg(transpose(wp_dense)),rw)+matmul(pp_dense,rp)
    call reduce_w_partial_to_owners(op%w_schedule,wpartial,wowned,info);if(info/=0)return
    candidate_yw=candidate_yw+wowned
    local_bad=merge(0,1,finite2(candidate_yw).and.finite2(candidate_yp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    yw=candidate_yw;yp=candidate_yp;info=0
  end subroutine

  subroutine global_gram_dg_wpw_bounded(op,xw,xp,yw,yp,g,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    complex(8),intent(in)::xw(:,:),xp(:,:),yw(:,:),yp(:,:);complex(8),intent(out)::g(:,:)
    integer,intent(out)::info
    complex(8)::local(size(xw,2),size(yw,2))
    integer::ierr,local_bad,global_bad,nx_min,nx_max,ny_min,ny_max
    info=1;g=0;local_bad=0
    if(.not.op%valid.or.size(xw,1)/=size(op%owned_w_ids).or.size(xp,1)/=size(op%owned_p_ids).or.&
       size(yw,1)/=size(op%owned_w_ids).or.size(yp,1)/=size(op%owned_p_ids).or.&
       size(xp,2)/=size(xw,2).or.size(yp,2)/=size(yw,2).or.&
       any(shape(g)/=[size(xw,2),size(yw,2)]).or..not.finite2(xw).or..not.finite2(xp).or.&
       .not.finite2(yw).or..not.finite2(yp))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(size(xw,2),nx_min,1,MPI_INTEGER,MPI_MIN,op%w_schedule%comm,ierr)
    call MPI_Allreduce(size(xw,2),nx_max,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    call MPI_Allreduce(size(yw,2),ny_min,1,MPI_INTEGER,MPI_MIN,op%w_schedule%comm,ierr)
    call MPI_Allreduce(size(yw,2),ny_max,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nx_min/=nx_max.or.ny_min/=ny_max)return
    local=matmul(conjg(transpose(xw)),yw)+matmul(conjg(transpose(xp)),yp)
    call MPI_Allreduce(local,g,size(g),MPI_DOUBLE_COMPLEX,MPI_SUM,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS)info=0
  end subroutine
  integer function find_id_sorted(ids,target)result(pos)
    integer,intent(in)::ids(:),target
    integer::lo,hi,mid
    pos=0;lo=1;hi=size(ids)
    do while(lo<=hi)
      mid=lo+(hi-lo)/2
      if(ids(mid)<target)then;lo=mid+1
      else if(ids(mid)>target)then;hi=mid-1
      else;pos=mid;return;endif
    enddo
  end function
  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:)
    integer::i
    ok=.true.
    do i=1,size(ids)
      if(ids(i)<=0)then;ok=.false.;return;endif
      if(i>1.and.ids(i)<=ids(i-1))then;ok=.false.;return;endif
    enddo
  end function
  logical function finite1(x)result(ok)
    complex(8),intent(in)::x(:);ok=all(ieee_is_finite(real(x,8))).and.all(ieee_is_finite(aimag(x)))
  end function
  logical function finite2(x)result(ok)
    complex(8),intent(in)::x(:,:);ok=all(ieee_is_finite(real(x,8))).and.all(ieee_is_finite(aimag(x)))
  end function
end module
