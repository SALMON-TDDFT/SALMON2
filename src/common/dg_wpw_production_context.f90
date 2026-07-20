module dg_wpw_production_context
  use,intrinsic::iso_fortran_env,only:int64
  use dg_wpw_bounded_operator,only:s_dg_wpw_bounded_operator,&
    s_dg_wpw_fragment_block_preconditioner,&
    s_dg_wpw_bounded_operator_snapshot,snapshot_dg_wpw_bounded_operator,&
    validate_dg_wpw_bounded_operator_snapshot,release_dg_wpw_bounded_operator_snapshot
  use mpi, only: MPI_COMM_NULL, MPI_Comm_rank, MPI_Comm_size, MPI_SUCCESS,&
    MPI_Allreduce,MPI_INTEGER,MPI_MAX
  implicit none
  private
  type,public::s_dg_wpw_production_context
    integer::comm=MPI_COMM_NULL,rank_id=-1,nrank=0,n_fragments=0,n_g=0,n_w=0
    integer::owned_fragment_id=0
    character(16)::ownership_kind=''
    integer::halo_epoch=0,scan_epoch=-1,operator_epoch=-1
    logical::quadrature_valid=.false.,face_valid=.false.,operator_valid=.false.,callbacks_bound=.false.
    logical::ww_valid=.false.
    character(32)::ww_metric_convention=''
    integer(int64)::ww_provenance_fingerprint=0_int64
    real(8),allocatable::ww_kinetic(:,:),ww_potential(:,:),ww_nonlocal(:,:),ww_face_self(:,:)
    logical::ww_projector_nonlocal_valid=.false.
    complex(8),allocatable::ww_projector_nonlocal(:,:),ww_projector_cross_value(:)
    integer,allocatable::ww_projector_cross_row_id(:),ww_projector_cross_col_id(:)
    integer,allocatable::ww_cross_face_id(:),ww_cross_row_id(:),ww_cross_col_id(:),ww_cross_axis(:),&
      ww_cross_side(:),ww_cross_image(:,:)
    real(8),allocatable::ww_cross_value(:)
    logical::pending_ww_valid=.false.
    character(32)::pending_ww_metric_convention=''
    integer(int64)::pending_ww_provenance_fingerprint=0_int64
    integer,allocatable::pending_ww_owned_w_ids(:)
    real(8),allocatable::pending_ww_kinetic(:,:),pending_ww_potential(:,:),pending_ww_nonlocal(:,:),&
      pending_ww_face_self(:,:)
    integer,allocatable::pending_ww_cross_face_id(:),pending_ww_cross_row_id(:),pending_ww_cross_col_id(:),&
      pending_ww_cross_axis(:),pending_ww_cross_side(:),pending_ww_cross_image(:,:)
    real(8),allocatable::pending_ww_cross_value(:)
    integer,allocatable::owned_column_ids(:),support_column_ids(:),owned_w_ids(:),support_w_ids(:)
    integer,allocatable::support_fragment_ids(:)
    type(s_dg_wpw_bounded_operator)::bounded_operator
    type(s_dg_wpw_fragment_block_preconditioner)::metric_block_preconditioner
    integer,allocatable::wp_w(:),wp_p(:),wp_origin(:),pp_r(:),pp_c(:),pp_origin(:)
    complex(8),allocatable::wp_h(:),wp_h_volume(:),wp_h_nonlocal(:),wp_h_face(:),wp_s(:)
    complex(8),allocatable::pp_h(:),pp_h_volume(:),pp_h_nonlocal(:),pp_s(:)
    integer,allocatable::h_wp_w(:),h_wp_p(:),h_pp_r(:),h_pp_c(:)
    integer,allocatable::s_wp_w(:),s_wp_p(:),s_pp_r(:),s_pp_c(:)
    complex(8),allocatable::h_wp(:),h_pp(:),s_wp(:),s_pp(:)
  contains
    procedure::apply_h=>apply_h_context
    procedure::apply_s=>apply_s_context
  end type
  type,public::s_dg_wpw_production_context_snapshot
    integer::halo_epoch=-1,scan_epoch=-1,operator_epoch=-1
    integer(int64)::ww_provenance_fingerprint=0_int64
    logical::valid=.false.,callbacks_bound=.false.,operator_valid=.false.,&
      ww_projector_nonlocal_valid=.false.
    character(16)::ownership_kind=''
    character(32)::ww_metric_convention=''
    complex(8),allocatable::wp_h_volume(:),wp_h_nonlocal(:),wp_h_face(:),&
      pp_h_volume(:),pp_h_nonlocal(:),ww_projector_nonlocal(:,:),&
      ww_projector_cross_value(:)
    integer,allocatable::ww_projector_cross_row_id(:),ww_projector_cross_col_id(:)
    type(s_dg_wpw_bounded_operator_snapshot)::bounded_operator
  end type
  public::initialize_dg_wpw_production_context,initialize_dg_wpw_fragment_root_context
  public::consume_dg_wpw_bounded_subspace
  public::snapshot_dg_wpw_production_context,validate_dg_wpw_production_context_snapshot
  public::release_dg_wpw_production_context_snapshot
contains
  subroutine snapshot_dg_wpw_production_context(ctx,snapshot,info)
    type(s_dg_wpw_production_context),intent(in)::ctx
    type(s_dg_wpw_production_context_snapshot),intent(inout)::snapshot
    integer,intent(out)::info
    integer::astat,local_bad,global_bad,ierr
    call release_dg_wpw_production_context_snapshot(snapshot)
    local_bad=merge(0,1,ctx%comm/=MPI_COMM_NULL.and.ctx%operator_valid.and.&
      ctx%callbacks_bound.and.ctx%ww_projector_nonlocal_valid)
    if(local_bad==0)then
      allocate(snapshot%wp_h_volume,source=ctx%wp_h_volume,stat=astat);if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%wp_h_nonlocal,source=ctx%wp_h_nonlocal,stat=astat);if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%wp_h_face,source=ctx%wp_h_face,stat=astat);if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%pp_h_volume,source=ctx%pp_h_volume,stat=astat);if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%pp_h_nonlocal,source=ctx%pp_h_nonlocal,stat=astat);if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%ww_projector_nonlocal,source=ctx%ww_projector_nonlocal,stat=astat)
      if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%ww_projector_cross_value,source=ctx%ww_projector_cross_value,stat=astat)
      if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%ww_projector_cross_row_id,source=ctx%ww_projector_cross_row_id,stat=astat)
      if(astat/=0)local_bad=1
    endif
    if(local_bad==0)then
      allocate(snapshot%ww_projector_cross_col_id,source=ctx%ww_projector_cross_col_id,stat=astat)
      if(astat/=0)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_production_context_snapshot(snapshot);info=1;return
    endif
    call snapshot_dg_wpw_bounded_operator(ctx%bounded_operator,snapshot%bounded_operator,info)
    if(info/=0)then;call release_dg_wpw_production_context_snapshot(snapshot);return;endif
    snapshot%halo_epoch=ctx%halo_epoch;snapshot%scan_epoch=ctx%scan_epoch
    snapshot%operator_epoch=ctx%operator_epoch
    snapshot%ww_provenance_fingerprint=ctx%ww_provenance_fingerprint
    snapshot%callbacks_bound=ctx%callbacks_bound;snapshot%operator_valid=ctx%operator_valid
    snapshot%ww_projector_nonlocal_valid=ctx%ww_projector_nonlocal_valid
    snapshot%ownership_kind=ctx%ownership_kind
    snapshot%ww_metric_convention=ctx%ww_metric_convention
    snapshot%valid=.true.;info=0
  end subroutine snapshot_dg_wpw_production_context

  subroutine validate_dg_wpw_production_context_snapshot(ctx,snapshot,info,allow_interface_lambda_change)
    type(s_dg_wpw_production_context),intent(in)::ctx
    type(s_dg_wpw_production_context_snapshot),intent(in)::snapshot
    integer,intent(out)::info
    logical,intent(in),optional::allow_interface_lambda_change
    integer::local_bad,global_bad,ierr,operator_info
    logical::allow_lambda
    allow_lambda=.false.;if(present(allow_interface_lambda_change))allow_lambda=allow_interface_lambda_change
    local_bad=merge(0,1,snapshot%valid.and.ctx%comm/=MPI_COMM_NULL.and.&
      snapshot%halo_epoch==ctx%halo_epoch.and.snapshot%scan_epoch==ctx%scan_epoch.and.&
      snapshot%operator_epoch==ctx%operator_epoch.and.&
      snapshot%ww_provenance_fingerprint==ctx%ww_provenance_fingerprint.and.&
      (snapshot%callbacks_bound.eqv.ctx%callbacks_bound).and.&
      (snapshot%operator_valid.eqv.ctx%operator_valid).and.&
      (snapshot%ww_projector_nonlocal_valid.eqv.ctx%ww_projector_nonlocal_valid).and.&
      snapshot%ownership_kind==ctx%ownership_kind.and.&
      snapshot%ww_metric_convention==ctx%ww_metric_convention)
    if(local_bad==0.and..not.same_context_z1(snapshot%wp_h_volume,ctx%wp_h_volume))local_bad=1
    if(local_bad==0.and..not.same_context_z1(snapshot%wp_h_nonlocal,ctx%wp_h_nonlocal))local_bad=1
    if(local_bad==0.and..not.same_context_z1(snapshot%wp_h_face,ctx%wp_h_face))local_bad=1
    if(local_bad==0.and..not.same_context_z1(snapshot%pp_h_volume,ctx%pp_h_volume))local_bad=1
    if(local_bad==0.and..not.same_context_z1(snapshot%pp_h_nonlocal,ctx%pp_h_nonlocal))local_bad=1
    if(local_bad==0.and..not.same_context_z2(snapshot%ww_projector_nonlocal,&
      ctx%ww_projector_nonlocal))local_bad=1
    if(local_bad==0.and..not.same_context_z1(snapshot%ww_projector_cross_value,&
      ctx%ww_projector_cross_value))local_bad=1
    if(local_bad==0.and..not.same_context_i1(snapshot%ww_projector_cross_row_id,&
      ctx%ww_projector_cross_row_id))local_bad=1
    if(local_bad==0.and..not.same_context_i1(snapshot%ww_projector_cross_col_id,&
      ctx%ww_projector_cross_col_id))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    call validate_dg_wpw_bounded_operator_snapshot(ctx%bounded_operator,&
      snapshot%bounded_operator,operator_info,allow_lambda)
    info=merge(0,1,operator_info==0)
  end subroutine validate_dg_wpw_production_context_snapshot

  logical function same_context_z1(left,right)result(ok)
    complex(8),allocatable,intent(in)::left(:),right(:)
    ok=allocated(left).and.allocated(right);if(.not.ok)return
    ok=size(left)==size(right);if(.not.ok)return
    ok=all(left==right)
  end function same_context_z1

  logical function same_context_z2(left,right)result(ok)
    complex(8),allocatable,intent(in)::left(:,:),right(:,:)
    ok=allocated(left).and.allocated(right);if(.not.ok)return
    ok=all(shape(left)==shape(right));if(.not.ok)return
    ok=all(left==right)
  end function same_context_z2

  logical function same_context_i1(left,right)result(ok)
    integer,allocatable,intent(in)::left(:),right(:)
    ok=allocated(left).and.allocated(right);if(.not.ok)return
    ok=size(left)==size(right);if(.not.ok)return
    ok=all(left==right)
  end function same_context_i1

  subroutine release_dg_wpw_production_context_snapshot(snapshot)
    type(s_dg_wpw_production_context_snapshot),intent(inout)::snapshot
    if(allocated(snapshot%wp_h_volume))deallocate(snapshot%wp_h_volume)
    if(allocated(snapshot%wp_h_nonlocal))deallocate(snapshot%wp_h_nonlocal)
    if(allocated(snapshot%wp_h_face))deallocate(snapshot%wp_h_face)
    if(allocated(snapshot%pp_h_volume))deallocate(snapshot%pp_h_volume)
    if(allocated(snapshot%pp_h_nonlocal))deallocate(snapshot%pp_h_nonlocal)
    if(allocated(snapshot%ww_projector_nonlocal))deallocate(snapshot%ww_projector_nonlocal)
    if(allocated(snapshot%ww_projector_cross_value))deallocate(snapshot%ww_projector_cross_value)
    if(allocated(snapshot%ww_projector_cross_row_id))deallocate(snapshot%ww_projector_cross_row_id)
    if(allocated(snapshot%ww_projector_cross_col_id))deallocate(snapshot%ww_projector_cross_col_id)
    call release_dg_wpw_bounded_operator_snapshot(snapshot%bounded_operator)
    snapshot%valid=.false.
  end subroutine release_dg_wpw_production_context_snapshot

  subroutine initialize_dg_wpw_production_context(ctx,comm,n_fragments,n_g,n_w,owned_w_ids,support_w_ids,info)
    type(s_dg_wpw_production_context),intent(out)::ctx
    integer,intent(in)::comm,n_fragments,n_g,n_w,owned_w_ids(:),support_w_ids(:)
    integer,intent(out)::info
    integer::ierr,first,last,i,ncol
    info=0
    if(comm==MPI_COMM_NULL.or.n_fragments<=0.or.n_g<=0.or.n_w<=0.or.&
       .not.strictly_increasing(owned_w_ids).or..not.strictly_increasing(support_w_ids).or.&
       any(owned_w_ids<1).or.any(owned_w_ids>n_w).or.any(support_w_ids<1).or.any(support_w_ids>n_w))then
      info=1;return
    endif
    do i=1,size(owned_w_ids)
      if(.not.any(support_w_ids==owned_w_ids(i)))then;info=1;return;endif
    enddo
    ctx%comm=comm;ctx%n_fragments=n_fragments;ctx%n_g=n_g;ctx%n_w=n_w
    ctx%ownership_kind='arithmetic'
    call MPI_Comm_rank(comm,ctx%rank_id,ierr);if(ierr/=MPI_SUCCESS)then;info=2;return;endif
    call MPI_Comm_size(comm,ctx%nrank,ierr);if(ierr/=MPI_SUCCESS.or.ctx%nrank<=0)then;info=2;return;endif
    ncol=n_fragments*n_g
    first=(ctx%rank_id*ncol+ctx%nrank-1)/ctx%nrank+1
    last=((ctx%rank_id+1)*ncol+ctx%nrank-1)/ctx%nrank
    allocate(ctx%owned_column_ids(max(0,last-first+1)),ctx%owned_w_ids(size(owned_w_ids)),&
      ctx%support_w_ids(size(support_w_ids)))
    ctx%owned_w_ids=owned_w_ids;ctx%support_w_ids=support_w_ids
    do i=1,size(ctx%owned_column_ids);ctx%owned_column_ids(i)=first+i-1;enddo
    ! Stable production column id is (K-1)*n_G+G_id.
  end subroutine

  subroutine initialize_dg_wpw_fragment_root_context(ctx,comm,n_fragments,n_g,n_w,owned_fragment_id, &
      support_fragment_ids,owned_w_ids,support_w_ids,info)
    type(s_dg_wpw_production_context),intent(out)::ctx
    integer,intent(in)::comm,n_fragments,n_g,n_w,owned_fragment_id,support_fragment_ids(:),&
      owned_w_ids(:),support_w_ids(:)
    integer,intent(out)::info
    integer::ierr,i,first
    integer(8)::ncol64,last64

    info=0
    if(comm==MPI_COMM_NULL.or.n_fragments<=0.or.n_g<=0.or.n_w<=0.or.&
       owned_fragment_id<1.or.owned_fragment_id>n_fragments.or.&
       .not.strictly_increasing(support_fragment_ids).or.any(support_fragment_ids<1).or.&
       any(support_fragment_ids>n_fragments).or..not.any(support_fragment_ids==owned_fragment_id).or.&
       .not.strictly_increasing(owned_w_ids).or..not.strictly_increasing(support_w_ids).or.&
       any(owned_w_ids<1).or.any(owned_w_ids>n_w).or.any(support_w_ids<1).or.any(support_w_ids>n_w))then
      info=1;return
    endif
    do i=1,size(owned_w_ids)
      if(.not.any(support_w_ids==owned_w_ids(i)))then;info=1;return;endif
    enddo
    ncol64=int(n_fragments,8)*int(n_g,8)
    last64=int(owned_fragment_id,8)*int(n_g,8)
    if(ncol64>int(huge(ctx%n_fragments),8).or.last64>int(huge(ctx%n_fragments),8))then
      info=1;return
    endif
    ctx%comm=comm;ctx%n_fragments=n_fragments;ctx%n_g=n_g;ctx%n_w=n_w
    ctx%owned_fragment_id=owned_fragment_id;ctx%ownership_kind='fragment_root'
    call MPI_Comm_rank(comm,ctx%rank_id,ierr);if(ierr/=MPI_SUCCESS)then;info=2;return;endif
    call MPI_Comm_size(comm,ctx%nrank,ierr);if(ierr/=MPI_SUCCESS.or.ctx%nrank/=n_fragments.or.&
       ctx%rank_id/=owned_fragment_id-1)then;info=2;return;endif
    allocate(ctx%owned_column_ids(n_g),ctx%support_column_ids(size(support_fragment_ids)*n_g),&
      ctx%owned_w_ids(size(owned_w_ids)),&
      ctx%support_w_ids(size(support_w_ids)),ctx%support_fragment_ids(size(support_fragment_ids)))
    ctx%owned_w_ids=owned_w_ids;ctx%support_w_ids=support_w_ids
    ctx%support_fragment_ids=support_fragment_ids
    first=(owned_fragment_id-1)*n_g+1
    do i=1,n_g;ctx%owned_column_ids(i)=first+i-1;enddo
    do i=1,size(support_fragment_ids)*n_g
      ctx%support_column_ids(i)=(support_fragment_ids((i-1)/n_g+1)-1)*n_g+modulo(i-1,n_g)+1
    enddo
  end subroutine initialize_dg_wpw_fragment_root_context
  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::j
    ok=.true.
    do j=2,size(ids);if(ids(j)<=ids(j-1))then;ok=.false.;return;endif;enddo
  end function
  subroutine consume_dg_wpw_bounded_subspace(ctx,nretain,info)
    type(s_dg_wpw_production_context),intent(in)::ctx
    integer,intent(in)::nretain
    integer,intent(out)::info
    info=merge(0,1,ctx%callbacks_bound.and.nretain>0.and.nretain<=ctx%n_w+ctx%n_fragments*ctx%n_g)
  end subroutine
  subroutine apply_h_context(ctx,x,y,info)
    class(s_dg_wpw_production_context),intent(in)::ctx;complex(8),intent(in)::x(:);complex(8),intent(out)::y(:);integer,intent(out)::info
    call apply_blocks(ctx,ctx%h_wp_w,ctx%h_wp_p,ctx%h_wp,ctx%h_pp_r,ctx%h_pp_c,ctx%h_pp,x,y,info)
    if(info==0.and.trim(ctx%ownership_kind)=='fragment_root')call add_ww_h(ctx,x,y,info)
  end subroutine
  subroutine apply_s_context(ctx,x,y,info)
    class(s_dg_wpw_production_context),intent(in)::ctx;complex(8),intent(in)::x(:);complex(8),intent(out)::y(:);integer,intent(out)::info
    call apply_blocks(ctx,ctx%s_wp_w,ctx%s_wp_p,ctx%s_wp,ctx%s_pp_r,ctx%s_pp_c,ctx%s_pp,x,y,info)
    if(info==0.and.trim(ctx%ownership_kind)=='fragment_root')call add_ww_s(ctx,x,y,info)
  end subroutine
  subroutine add_ww_h(ctx,x,y,info)
    class(s_dg_wpw_production_context),intent(in)::ctx
    complex(8),intent(in)::x(:);complex(8),intent(inout)::y(:);integer,intent(out)::info
    integer::i,j,row,col
    info=1
    if(.not.ctx%ww_valid.or.any(shape(ctx%ww_kinetic)/=[size(ctx%owned_w_ids),size(ctx%owned_w_ids)]))then
      y=0;return
    endif
    do j=1,size(ctx%owned_w_ids);col=ctx%owned_w_ids(j)
      do i=1,size(ctx%owned_w_ids);row=ctx%owned_w_ids(i)
        if(ctx%ww_projector_nonlocal_valid)then
          y(row)=y(row)+(cmplx(ctx%ww_kinetic(i,j)+ctx%ww_potential(i,j)+ctx%ww_face_self(i,j),0d0,8)+&
            ctx%ww_projector_nonlocal(i,j))*x(col)
        else
          y(row)=y(row)+cmplx(ctx%ww_kinetic(i,j)+ctx%ww_potential(i,j)+ctx%ww_nonlocal(i,j)+&
            ctx%ww_face_self(i,j),0d0,8)*x(col)
        endif
      enddo
    enddo
    do i=1,size(ctx%ww_cross_value)
      row=ctx%ww_cross_row_id(i);col=ctx%ww_cross_col_id(i)
      if(row<1.or.row>ctx%n_w.or.col<1.or.col>ctx%n_w)then;y=0;return;endif
      y(row)=y(row)+ctx%ww_cross_value(i)*x(col)
      y(col)=y(col)+ctx%ww_cross_value(i)*x(row)
    enddo
    if(ctx%ww_projector_nonlocal_valid)then
      do i=1,size(ctx%ww_projector_cross_value)
        row=ctx%ww_projector_cross_row_id(i);col=ctx%ww_projector_cross_col_id(i)
        if(row<1.or.row>ctx%n_w.or.col<1.or.col>ctx%n_w)then;y=0;return;endif
        y(row)=y(row)+ctx%ww_projector_cross_value(i)*x(col)
      enddo
    endif
    info=0
  end subroutine
  subroutine add_ww_s(ctx,x,y,info)
    class(s_dg_wpw_production_context),intent(in)::ctx
    complex(8),intent(in)::x(:);complex(8),intent(inout)::y(:);integer,intent(out)::info
    integer::i
    if(.not.ctx%ww_valid)then;y=0;info=1;return;endif
    do i=1,size(ctx%owned_w_ids);y(ctx%owned_w_ids(i))=y(ctx%owned_w_ids(i))+x(ctx%owned_w_ids(i));enddo
    info=0
  end subroutine
  subroutine apply_blocks(ctx,wr,pc,wv,pr,qc,pv,x,y,info)
    class(s_dg_wpw_production_context),intent(in)::ctx
    integer,intent(in)::wr(:),pc(:),pr(:),qc(:);complex(8),intent(in)::wv(:),pv(:),x(:);complex(8),intent(out)::y(:);integer,intent(out)::info
    integer::i,np
    y=(0d0,0d0);info=0;np=ctx%n_fragments*ctx%n_g
    if(.not.ctx%callbacks_bound.or.size(x)/=ctx%n_w+np.or.size(y)/=size(x))then;info=1;return;endif
    do i=1,size(wv)
      y(wr(i))=y(wr(i))+wv(i)*x(ctx%n_w+pc(i))
      y(ctx%n_w+pc(i))=y(ctx%n_w+pc(i))+conjg(wv(i))*x(wr(i))
    enddo
    do i=1,size(pv);y(ctx%n_w+pr(i))=y(ctx%n_w+pr(i))+pv(i)*x(ctx%n_w+qc(i));enddo
  end subroutine
end module
