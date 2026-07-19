module rt_dg_wpw_production_builder
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use mpi,only:MPI_Allreduce,MPI_INTEGER,MPI_MAX,MPI_SUCCESS
  use dg_wpw_production_context,only:s_dg_wpw_production_context
  use dg_wpw_bounded_operator,only:initialize_dg_wpw_bounded_operator
  use rt_dg_wpw_column_layout,only:s_dg_wpw_column_layout,initialize_wpw_column_layout,&
    initialize_wpw_fragment_root_layout
  use rt_dg_wpw_sparse_builder,only:s_dg_wpw_sparse_candidates,build_windowed_sparse_wpw_operators,&
    wpw_candidate_volume,wpw_candidate_volume_face
  use rt_dg_wpw_sparse_blocks,only:s_dg_wpw_sparse_blocks
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider
  use rt_dg_wpw_face_trace_scanner,only:assemble_wpw_canonical_face_grid
  use rt_dg_wpw_candidate_halo,only:s_dg_wpw_staged_candidate,s_dg_wpw_owned_candidates,&
    route_dg_wpw_candidate_halo,wpw_candidate_kind_wp,wpw_candidate_kind_pp
  implicit none
  private
  public::build_dg_wpw_rank_local_quadrature,scan_dg_wpw_canonical_faces
  public::route_dg_wpw_staged_candidates
  public::route_and_replace_dg_wpw_potential_volume
  public::scan_and_route_dg_wpw_canonical_face
  public::scan_and_route_dg_wpw_canonical_faces
  public::replace_dg_wpw_potential_volume
  public::install_dg_wpw_nonlocal_volume
  public::install_dg_wpw_projector_nonlocal
  public::build_dg_wpw_production_operator,bind_dg_wpw_hs_callbacks
contains
  subroutine replace_dg_wpw_potential_volume(ctx,wp_h_volume,pp_h_volume,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    complex(8),intent(in)::wp_h_volume(:),pp_h_volume(:)
    integer,intent(out)::info
    integer::local_bad,global_bad,ierr

    local_bad=0
    if(.not.ctx%quadrature_valid.or..not.allocated(ctx%wp_h_volume).or.&
       .not.allocated(ctx%wp_h_face).or..not.allocated(ctx%pp_h_volume))then
      local_bad=1
    else if(size(wp_h_volume)/=size(ctx%wp_h_volume).or.&
            size(pp_h_volume)/=size(ctx%pp_h_volume))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    ctx%wp_h_volume=wp_h_volume
    ctx%wp_h=ctx%wp_h_volume+ctx%wp_h_nonlocal+ctx%wp_h_face
    ctx%pp_h_volume=pp_h_volume
    ctx%pp_h=ctx%pp_h_volume+ctx%pp_h_nonlocal
    ctx%operator_valid=.false.;ctx%callbacks_bound=.false.
    info=0
  end subroutine replace_dg_wpw_potential_volume

  subroutine install_dg_wpw_nonlocal_volume(ctx,wp_h_nonlocal,pp_h_nonlocal,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    complex(8),intent(in)::wp_h_nonlocal(:),pp_h_nonlocal(:)
    integer,intent(out)::info
    integer::local_bad,global_bad,ierr

    local_bad=0
    if(.not.ctx%quadrature_valid.or..not.allocated(ctx%wp_h_nonlocal).or.&
       .not.allocated(ctx%pp_h_nonlocal))then
      local_bad=1
    else if(size(wp_h_nonlocal)/=size(ctx%wp_h_nonlocal).or.&
            size(pp_h_nonlocal)/=size(ctx%pp_h_nonlocal))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    ctx%wp_h_nonlocal=wp_h_nonlocal
    ctx%pp_h_nonlocal=pp_h_nonlocal
    ctx%wp_h=ctx%wp_h_volume+ctx%wp_h_nonlocal+ctx%wp_h_face
    ctx%pp_h=ctx%pp_h_volume+ctx%pp_h_nonlocal
    ctx%operator_valid=.false.;ctx%callbacks_bound=.false.
    info=0
  end subroutine install_dg_wpw_nonlocal_volume

  subroutine install_dg_wpw_projector_nonlocal(ctx,ww_local,ww_cross_row,ww_cross_col,ww_cross_value,&
      wp_nonlocal,pp_nonlocal,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    complex(8),intent(in)::ww_local(:,:),ww_cross_value(:),wp_nonlocal(:),pp_nonlocal(:)
    integer,intent(in)::ww_cross_row(:),ww_cross_col(:)
    integer,intent(out)::info
    integer::i,j,local_bad,global_bad,ierr,astat
    integer,allocatable::new_cross_row(:),new_cross_col(:)
    complex(8),allocatable::new_ww(:,:),new_cross(:),new_wp(:),new_pp(:),new_wp_total(:),new_pp_total(:)

    local_bad=0
    if(trim(ctx%ownership_kind)/='fragment_root'.or..not.ctx%quadrature_valid.or.&
       .not.(ctx%ww_valid.or.ctx%pending_ww_valid).or..not.allocated(ctx%owned_w_ids).or.&
       .not.allocated(ctx%support_w_ids).or..not.allocated(ctx%wp_h_nonlocal).or.&
       .not.allocated(ctx%pp_h_nonlocal).or..not.allocated(ctx%wp_h_volume).or.&
       .not.allocated(ctx%wp_h_face).or..not.allocated(ctx%pp_h_volume).or.&
       .not.allocated(ctx%wp_h).or..not.allocated(ctx%pp_h))then
      local_bad=1
    else
      if(any(shape(ww_local)/=[size(ctx%owned_w_ids),size(ctx%owned_w_ids)]).or.&
         size(ww_cross_row)/=size(ww_cross_col).or.size(ww_cross_row)/=size(ww_cross_value).or.&
         size(wp_nonlocal)/=size(ctx%wp_h_nonlocal).or.size(pp_nonlocal)/=size(ctx%pp_h_nonlocal).or.&
         size(ctx%wp_h_volume)/=size(wp_nonlocal).or.size(ctx%wp_h_face)/=size(wp_nonlocal).or.&
         size(ctx%pp_h_volume)/=size(pp_nonlocal))local_bad=1
      do i=1,size(ww_cross_row)
        if(.not.any(ctx%owned_w_ids==ww_cross_row(i)).or..not.any(ctx%support_w_ids==ww_cross_col(i)))local_bad=1
        do j=1,i-1
          if(ww_cross_row(j)==ww_cross_row(i).and.ww_cross_col(j)==ww_cross_col(i))local_bad=1
        enddo
      enddo
    endif
    if(.not.all_finite_complex_2d(ww_local).or..not.all_finite_complex(ww_cross_value).or.&
       .not.all_finite_complex(wp_nonlocal).or..not.all_finite_complex(pp_nonlocal))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    astat=0
    allocate(new_ww(size(ww_local,1),size(ww_local,2)),new_cross(size(ww_cross_value)),&
      new_wp(size(wp_nonlocal)),new_pp(size(pp_nonlocal)),new_wp_total(size(wp_nonlocal)),&
      new_pp_total(size(pp_nonlocal)),new_cross_row(size(ww_cross_row)),&
      new_cross_col(size(ww_cross_col)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    new_ww=ww_local;new_cross=ww_cross_value;new_cross_row=ww_cross_row;new_cross_col=ww_cross_col
    new_wp=wp_nonlocal;new_pp=pp_nonlocal
    new_wp_total=ctx%wp_h_volume+new_wp+ctx%wp_h_face
    new_pp_total=ctx%pp_h_volume+new_pp
    local_bad=merge(0,1,all_finite_complex(new_wp_total).and.all_finite_complex(new_pp_total))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    call move_alloc(new_ww,ctx%ww_projector_nonlocal)
    call move_alloc(new_cross_row,ctx%ww_projector_cross_row_id)
    call move_alloc(new_cross_col,ctx%ww_projector_cross_col_id)
    call move_alloc(new_cross,ctx%ww_projector_cross_value);ctx%ww_projector_nonlocal_valid=.true.
    call move_alloc(new_wp,ctx%wp_h_nonlocal);call move_alloc(new_pp,ctx%pp_h_nonlocal)
    call move_alloc(new_wp_total,ctx%wp_h);call move_alloc(new_pp_total,ctx%pp_h)
    ctx%operator_valid=.false.;ctx%callbacks_bound=.false.;info=0
  end subroutine install_dg_wpw_projector_nonlocal

  pure logical function all_finite_complex(values)
    complex(8),intent(in)::values(:)
    all_finite_complex=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function all_finite_complex

  pure logical function all_finite_complex_2d(values)
    complex(8),intent(in)::values(:,:)
    all_finite_complex_2d=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function all_finite_complex_2d

  subroutine scan_and_route_dg_wpw_canonical_faces(ctx,epoch,support_fragments,provider,k_minus,&
      k_plus,axes,sides,minus_lo,minus_hi,plus_lo,plus_hi,hgs,w_ids,p_ids,max_candidates_per_peer,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::epoch,support_fragments(:),k_minus(:),k_plus(:),axes(:),sides(:)
    integer,intent(in)::minus_lo(:,:),minus_hi(:,:),plus_lo(:,:),plus_hi(:,:)
    integer,intent(in)::w_ids(:),p_ids(:),max_candidates_per_peer
    real(8),intent(in)::hgs(3)
    type(s_wpw_face_trace_provider),intent(inout)::provider
    integer,intent(out)::info
    type(s_dg_wpw_staged_candidate),allocatable::staged(:)
    type(s_dg_wpw_owned_candidates)::owned
    complex(8),allocatable::face_h(:,:),candidate_face(:)
    integer,allocatable::candidate_origin(:),match(:)
    integer::nface,iface,iw,ip,i,entry,nmatch,local_bad,global_bad,ierr

    info=1;local_bad=0;nface=size(k_minus)
    if(.not.ctx%quadrature_valid.or.epoch/=ctx%halo_epoch.or.size(k_plus)/=nface.or.&
       size(axes)/=nface.or.size(sides)/=nface.or.any(shape(minus_lo)/=[3,nface]).or.&
       any(shape(minus_hi)/=[3,nface]).or.any(shape(plus_lo)/=[3,nface]).or.&
       any(shape(plus_hi)/=[3,nface]).or.size(w_ids)<=0.or.size(p_ids)<=0)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(staged(nface*size(w_ids)*size(p_ids)));i=0
    do iface=1,nface
      allocate(face_h(size(w_ids),size(p_ids)))
      call assemble_wpw_canonical_face_grid(provider,k_minus(iface),k_plus(iface),axes(iface),&
        sides(iface),minus_lo(:,iface),minus_hi(:,iface),plus_lo(:,iface),plus_hi(:,iface),&
        hgs,w_ids,p_ids,face_h,info)
      if(info/=0)then;local_bad=1;deallocate(face_h);cycle;endif
      do ip=1,size(p_ids);do iw=1,size(w_ids)
        i=i+1;staged(i)%kind=wpw_candidate_kind_wp;staged(i)%image_id=iface
        staged(i)%wp_w=w_ids(iw);staged(i)%wp_p=p_ids(ip)
        staged(i)%wp_h=face_h(iw,ip);staged(i)%wp_s=(0d0,0d0)
      enddo;enddo
      deallocate(face_h)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call route_dg_wpw_candidate_halo(ctx%comm,epoch,ctx%n_fragments,ctx%n_g,support_fragments,&
      staged,owned,info,max_candidates_per_peer)
    if(info/=0)return
    allocate(match(size(owned%wp_w)));match=0
    do i=1,size(owned%wp_w)
      nmatch=0
      do entry=1,size(ctx%wp_w)
        if(ctx%wp_w(entry)==owned%wp_w(i).and.ctx%wp_p(entry)==owned%wp_p(i))then
          match(i)=entry;nmatch=nmatch+1
        endif
      enddo
      if(nmatch/=1)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif
    allocate(candidate_face,source=ctx%wp_h_face);allocate(candidate_origin,source=ctx%wp_origin)
    do i=1,size(owned%wp_w)
      candidate_face(match(i))=candidate_face(match(i))+owned%wp_h(i)
      candidate_origin(match(i))=wpw_candidate_volume_face
    enddo
    ctx%wp_h_face=candidate_face;ctx%wp_h=ctx%wp_h_volume+ctx%wp_h_nonlocal+ctx%wp_h_face;ctx%wp_origin=candidate_origin
    ctx%scan_epoch=epoch;ctx%face_valid=.true.;ctx%operator_valid=.false.;ctx%callbacks_bound=.false.
    info=0
  end subroutine scan_and_route_dg_wpw_canonical_faces

  subroutine scan_and_route_dg_wpw_canonical_face(ctx,epoch,support_fragments,canonical_owner,&
      provider,k_minus,k_plus,axis,side_from_k_minus,minus_lo,minus_hi,plus_lo,plus_hi,hgs,&
      w_ids,p_ids,max_candidates_per_peer,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::epoch,support_fragments(:),k_minus,k_plus,axis,side_from_k_minus
    integer,intent(in)::minus_lo(3),minus_hi(3),plus_lo(3),plus_hi(3),w_ids(:),p_ids(:)
    integer,intent(in)::max_candidates_per_peer
    logical,intent(in)::canonical_owner
    real(8),intent(in)::hgs(3)
    type(s_wpw_face_trace_provider),intent(inout)::provider
    integer,intent(out)::info
    type(s_dg_wpw_staged_candidate),allocatable::staged(:)
    type(s_dg_wpw_owned_candidates)::owned
    complex(8),allocatable::face_h(:,:),candidate_face(:)
    integer,allocatable::candidate_origin(:),match(:)
    integer::iw,ip,i,entry,nmatch,local_bad,global_bad,ierr

    info=1;local_bad=0
    if(.not.ctx%quadrature_valid.or.epoch/=ctx%halo_epoch)then;local_bad=1;endif
    if(canonical_owner)then
      allocate(face_h(size(w_ids),size(p_ids)))
      call assemble_wpw_canonical_face_grid(provider,k_minus,k_plus,axis,side_from_k_minus,&
        minus_lo,minus_hi,plus_lo,plus_hi,hgs,w_ids,p_ids,face_h,info)
      if(info/=0)then;local_bad=1;allocate(staged(0))
      else
        allocate(staged(size(w_ids)*size(p_ids)));i=0
        do ip=1,size(p_ids);do iw=1,size(w_ids)
          i=i+1;staged(i)%kind=wpw_candidate_kind_wp
          staged(i)%image_id=1+2*(axis-1)+merge(1,0,side_from_k_minus>0)
          staged(i)%wp_w=w_ids(iw);staged(i)%wp_p=p_ids(ip)
          staged(i)%wp_h=face_h(iw,ip);staged(i)%wp_s=(0d0,0d0)
        enddo;enddo
      endif
    else
      allocate(staged(0))
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call route_dg_wpw_candidate_halo(ctx%comm,epoch,ctx%n_fragments,ctx%n_g,support_fragments,&
      staged,owned,info,max_candidates_per_peer)
    if(info/=0)return
    allocate(match(size(owned%wp_w)));match=0
    do i=1,size(owned%wp_w)
      nmatch=0
      do entry=1,size(ctx%wp_w)
        if(ctx%wp_w(entry)==owned%wp_w(i).and.ctx%wp_p(entry)==owned%wp_p(i))then
          match(i)=entry;nmatch=nmatch+1
        endif
      enddo
      if(nmatch/=1)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif
    allocate(candidate_face,source=ctx%wp_h_face);allocate(candidate_origin,source=ctx%wp_origin)
    do i=1,size(owned%wp_w)
      candidate_face(match(i))=candidate_face(match(i))+owned%wp_h(i)
      candidate_origin(match(i))=wpw_candidate_volume_face
    enddo
    ctx%wp_h_face=candidate_face;ctx%wp_h=ctx%wp_h_volume+ctx%wp_h_nonlocal+ctx%wp_h_face;ctx%wp_origin=candidate_origin
    ctx%scan_epoch=epoch;ctx%face_valid=.true.;ctx%operator_valid=.false.;ctx%callbacks_bound=.false.
    info=0
  end subroutine scan_and_route_dg_wpw_canonical_face

  subroutine route_dg_wpw_staged_candidates(ctx,epoch,support_fragments,wr,pc,wh,ws,pr,qc,ph,ps,&
      max_candidates_per_peer,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::epoch,support_fragments(:),wr(:),pc(:),pr(:),qc(:),max_candidates_per_peer
    complex(8),intent(in)::wh(:),ws(:),ph(:),ps(:)
    integer,intent(out)::info
    type(s_dg_wpw_staged_candidate),allocatable::staged(:)
    type(s_dg_wpw_owned_candidates)::owned
    integer::i,nwp
    info=1;nwp=size(wr)
    if(size(pc)/=nwp.or.size(wh)/=nwp.or.size(ws)/=nwp.or.size(qc)/=size(pr).or.&
       size(ph)/=size(pr).or.size(ps)/=size(pr))return
    allocate(staged(nwp+size(pr)))
    do i=1,nwp
      staged(i)%kind=wpw_candidate_kind_wp;staged(i)%image_id=i
      staged(i)%wp_w=wr(i);staged(i)%wp_p=pc(i);staged(i)%wp_h=wh(i);staged(i)%wp_s=ws(i)
    enddo
    do i=1,size(pr)
      staged(nwp+i)%kind=wpw_candidate_kind_pp;staged(nwp+i)%image_id=nwp+i
      staged(nwp+i)%pp_r=pr(i);staged(nwp+i)%pp_c=qc(i)
      staged(nwp+i)%pp_h=ph(i);staged(nwp+i)%pp_s=ps(i)
    enddo
    call route_dg_wpw_candidate_halo(ctx%comm,epoch,ctx%n_fragments,ctx%n_g,support_fragments,&
      staged,owned,info,max_candidates_per_peer)
    if(info/=0)return
    call build_dg_wpw_rank_local_quadrature(ctx,owned%wp_w,owned%wp_p,owned%wp_h,owned%wp_s,&
      owned%pp_r,owned%pp_c,owned%pp_h,owned%pp_s,info)
  end subroutine route_dg_wpw_staged_candidates

  subroutine route_and_replace_dg_wpw_potential_volume(ctx,epoch,support_fragments,wr,pc,wh,pr,qc,ph,&
      max_candidates_per_peer,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::epoch,support_fragments(:),wr(:),pc(:),pr(:),qc(:),max_candidates_per_peer
    complex(8),intent(in)::wh(:),ph(:)
    integer,intent(out)::info
    type(s_dg_wpw_staged_candidate),allocatable::staged(:)
    type(s_dg_wpw_owned_candidates)::owned
    complex(8),allocatable::wp_volume(:),pp_volume(:)
    integer::i,j,nmatch,local_bad,global_bad,ierr
    info=1;local_bad=0
    if(size(wr)/=size(pc).or.size(wr)/=size(wh).or.size(pr)/=size(qc).or.size(pr)/=size(ph))local_bad=1
    allocate(staged(size(wr)+size(pr)))
    do i=1,size(wr)
      staged(i)%kind=wpw_candidate_kind_wp;staged(i)%image_id=i
      staged(i)%wp_w=wr(i);staged(i)%wp_p=pc(i);staged(i)%wp_h=wh(i);staged(i)%wp_s=(0d0,0d0)
    enddo
    do i=1,size(pr)
      j=size(wr)+i;staged(j)%kind=wpw_candidate_kind_pp;staged(j)%image_id=j
      staged(j)%pp_r=pr(i);staged(j)%pp_c=qc(i);staged(j)%pp_h=ph(i);staged(j)%pp_s=(0d0,0d0)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call route_dg_wpw_candidate_halo(ctx%comm,epoch,ctx%n_fragments,ctx%n_g,support_fragments,staged,owned,info,&
      max_candidates_per_peer)
    if(info/=0)return
    allocate(wp_volume(size(ctx%wp_h_volume)),pp_volume(size(ctx%pp_h_volume)));wp_volume=0;pp_volume=0
    do i=1,size(owned%wp_w)
      nmatch=0
      do j=1,size(ctx%wp_w)
        if(ctx%wp_w(j)==owned%wp_w(i).and.ctx%wp_p(j)==owned%wp_p(i))then
          wp_volume(j)=owned%wp_h(i);nmatch=nmatch+1
        endif
      enddo
      if(nmatch/=1)local_bad=1
    enddo
    do i=1,size(owned%pp_r)
      nmatch=0
      do j=1,size(ctx%pp_r)
        if(ctx%pp_r(j)==owned%pp_r(i).and.ctx%pp_c(j)==owned%pp_c(i))then
          pp_volume(j)=owned%pp_h(i);nmatch=nmatch+1
        endif
      enddo
      if(nmatch/=1)local_bad=1
    enddo
    if(size(owned%wp_w)/=size(ctx%wp_w).or.size(owned%pp_r)/=size(ctx%pp_r))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,ctx%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif
    call replace_dg_wpw_potential_volume(ctx,wp_volume,pp_volume,info)
  end subroutine route_and_replace_dg_wpw_potential_volume

  subroutine build_dg_wpw_rank_local_quadrature(ctx,wr,pc,wh,ws,pr,qc,ph,ps,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::wr(:),pc(:),pr(:),qc(:);complex(8),intent(in)::wh(:),ws(:),ph(:),ps(:);integer,intent(out)::info
    integer::i
    integer,allocatable::new_wp_w(:),new_wp_p(:),new_wp_origin(:),new_pp_r(:),new_pp_c(:),new_pp_origin(:)
    complex(8),allocatable::new_wp_h(:),new_wp_h_volume(:),new_wp_h_nonlocal(:),new_wp_h_face(:),new_wp_s(:)
    complex(8),allocatable::new_pp_h(:),new_pp_h_volume(:),new_pp_h_nonlocal(:),new_pp_s(:)
    info=0
    if(size(wr)/=size(pc).or.size(wr)/=size(wh).or.size(wr)/=size(ws).or.&
       size(pr)/=size(qc).or.size(pr)/=size(ph).or.size(pr)/=size(ps))then;info=1;return;endif
    do i=1,size(pc);if(.not.any(ctx%owned_column_ids==pc(i)))then;info=2;return;endif;enddo
    do i=1,size(pr);if(.not.any(ctx%owned_column_ids==pr(i)))then;info=2;return;endif;enddo
    allocate(new_wp_w,source=wr);allocate(new_wp_p,source=pc);allocate(new_wp_origin(size(wr)))
    allocate(new_wp_h,source=wh);allocate(new_wp_h_volume,source=wh)
    allocate(new_wp_h_nonlocal(size(wh)),source=(0d0,0d0))
    allocate(new_wp_h_face(size(wh)),source=(0d0,0d0));allocate(new_wp_s,source=ws)
    allocate(new_pp_r,source=pr);allocate(new_pp_c,source=qc);allocate(new_pp_origin(size(pr)))
    allocate(new_pp_h,source=ph);allocate(new_pp_h_volume,source=ph);allocate(new_pp_s,source=ps)
    allocate(new_pp_h_nonlocal(size(ph)),source=(0d0,0d0))
    new_wp_origin=wpw_candidate_volume;new_pp_origin=wpw_candidate_volume
    if(allocated(ctx%wp_w))deallocate(ctx%wp_w,ctx%wp_p,ctx%wp_origin,ctx%wp_h,ctx%wp_h_volume,&
      ctx%wp_h_nonlocal,ctx%wp_h_face,ctx%wp_s)
    if(allocated(ctx%pp_r))deallocate(ctx%pp_r,ctx%pp_c,ctx%pp_origin,ctx%pp_h,ctx%pp_h_volume,&
      ctx%pp_h_nonlocal,ctx%pp_s)
    call move_alloc(new_wp_w,ctx%wp_w);call move_alloc(new_wp_p,ctx%wp_p)
    call move_alloc(new_wp_origin,ctx%wp_origin);call move_alloc(new_wp_h,ctx%wp_h)
    call move_alloc(new_wp_h_volume,ctx%wp_h_volume);call move_alloc(new_wp_h_face,ctx%wp_h_face)
    call move_alloc(new_wp_h_nonlocal,ctx%wp_h_nonlocal)
    call move_alloc(new_wp_s,ctx%wp_s)
    call move_alloc(new_pp_r,ctx%pp_r);call move_alloc(new_pp_c,ctx%pp_c)
    call move_alloc(new_pp_origin,ctx%pp_origin);call move_alloc(new_pp_h,ctx%pp_h)
    call move_alloc(new_pp_h_volume,ctx%pp_h_volume);call move_alloc(new_pp_s,ctx%pp_s)
    call move_alloc(new_pp_h_nonlocal,ctx%pp_h_nonlocal)
    ctx%quadrature_valid=.true.;ctx%face_valid=.false.;ctx%operator_valid=.false.;ctx%callbacks_bound=.false.
  end subroutine
  subroutine scan_dg_wpw_canonical_faces(ctx,epoch,provider,k_minus,k_plus,axis,side_from_k_minus,&
      minus_lo,minus_hi,plus_lo,plus_hi,hgs,w_ids,p_ids,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(in)::epoch,k_minus,k_plus,axis,side_from_k_minus
    integer,intent(in)::minus_lo(3),minus_hi(3),plus_lo(3),plus_hi(3),w_ids(:),p_ids(:)
    real(8),intent(in)::hgs(3)
    type(s_wpw_face_trace_provider),intent(inout)::provider
    integer,intent(out)::info
    integer::iw,ip,ientry,nmatch
    complex(8),allocatable::face_h(:,:),candidate_wp_h_face(:)
    integer,allocatable::candidate_origin(:),match_entry(:,:)
    info=0
    if(.not.ctx%quadrature_valid.or.epoch/=ctx%halo_epoch)then;info=1;return;endif
    allocate(face_h(size(w_ids),size(p_ids)),match_entry(size(w_ids),size(p_ids)))
    call assemble_wpw_canonical_face_grid(provider,k_minus,k_plus,axis,side_from_k_minus,&
      minus_lo,minus_hi,plus_lo,plus_hi,hgs,w_ids,p_ids,face_h,info)
    if(info/=0)return
    match_entry=0
    do ip=1,size(p_ids)
      do iw=1,size(w_ids)
        nmatch=0
        do ientry=1,size(ctx%wp_h)
          if(ctx%wp_w(ientry)==w_ids(iw).and.ctx%wp_p(ientry)==p_ids(ip))then
            match_entry(iw,ip)=ientry
            nmatch=nmatch+1
          endif
        enddo
        if(nmatch/=1)then;info=2;return;endif
      enddo
    enddo
    allocate(candidate_wp_h_face,source=ctx%wp_h_face);allocate(candidate_origin,source=ctx%wp_origin)
    do ip=1,size(p_ids);do iw=1,size(w_ids)
      ientry=match_entry(iw,ip);candidate_wp_h_face(ientry)=candidate_wp_h_face(ientry)+face_h(iw,ip)
      candidate_origin(ientry)=wpw_candidate_volume_face
    enddo;enddo
    ctx%wp_h_face=candidate_wp_h_face;ctx%wp_h=ctx%wp_h_volume+ctx%wp_h_nonlocal+ctx%wp_h_face;ctx%wp_origin=candidate_origin
    ctx%scan_epoch=epoch;ctx%face_valid=.true.
  end subroutine
  subroutine build_dg_wpw_production_operator(ctx,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    integer,intent(out)::info
    type(s_dg_wpw_column_layout)::layout
    type(s_dg_wpw_sparse_candidates)::c
    type(s_dg_wpw_sparse_blocks)::hb,sb
    real(8)::stage_clock
    call cpu_time(stage_clock)
    info=0
    if(.not.ctx%face_valid.or.ctx%scan_epoch/=ctx%halo_epoch.or.&
       (trim(ctx%ownership_kind)=='fragment_root'.and..not.(ctx%ww_valid.or.ctx%pending_ww_valid)))then
      info=1;return
    endif
    select case(trim(ctx%ownership_kind))
    case('fragment_root')
      call initialize_wpw_fragment_root_layout(layout,ctx%n_fragments,ctx%n_g,ctx%owned_fragment_id,&
        ctx%rank_id,ctx%nrank,info)
    case('arithmetic')
      call initialize_wpw_column_layout(layout,ctx%n_fragments,ctx%n_g,ctx%rank_id,ctx%nrank,info)
    case default
      info=2
    end select
    if(info/=0)return
    allocate(c%wp_w_row_ids(size(ctx%wp_w)),c%wp_pw_col_ids(size(ctx%wp_p)),c%wp_origin(size(ctx%wp_origin)))
    allocate(c%wp_h_values(size(ctx%wp_h),1),c%wp_s_values(size(ctx%wp_s),1))
    allocate(c%pp_pw_row_ids(size(ctx%pp_r)),c%pp_pw_col_ids(size(ctx%pp_c)),c%pp_origin(size(ctx%pp_origin)))
    allocate(c%pp_h_values(size(ctx%pp_h),1),c%pp_s_values(size(ctx%pp_s),1))
    c%wp_w_row_ids=ctx%wp_w;c%wp_pw_col_ids=ctx%wp_p;c%wp_origin=ctx%wp_origin;c%wp_h_values(:,1)=ctx%wp_h;c%wp_s_values(:,1)=ctx%wp_s
    c%pp_pw_row_ids=ctx%pp_r;c%pp_pw_col_ids=ctx%pp_c;c%pp_origin=ctx%pp_origin;c%pp_h_values(:,1)=ctx%pp_h;c%pp_s_values(:,1)=ctx%pp_s
    call trace_production_builder('candidate_copy',ctx,stage_clock,size(c%wp_w_row_ids),size(c%pp_pw_row_ids))
    call build_windowed_sparse_wpw_operators(layout,ctx%rank_id,ctx%nrank,c,hb,sb,info);if(info/=0)return
    call trace_production_builder('sparse_blocks',ctx,stage_clock,size(hb%wp_w_row_ids),size(hb%pp_pw_row_ids))
    if(trim(ctx%ownership_kind)=='fragment_root')then
      call build_bounded_operator(ctx,hb,sb,info);if(info/=0)return
      call trace_production_builder('bounded_operator',ctx,stage_clock,size(hb%wp_w_row_ids),size(hb%pp_pw_row_ids))
      if(ctx%pending_ww_valid)call commit_pending_ww(ctx)
    endif
    call copy_blocks(ctx,hb,sb);ctx%operator_epoch=ctx%halo_epoch
    call trace_production_builder('block_publication',ctx,stage_clock,size(hb%wp_w_row_ids),size(hb%pp_pw_row_ids))
    ctx%operator_valid=.true.;ctx%callbacks_bound=.false.
  end subroutine
  subroutine trace_production_builder(stage,ctx,stage_clock,nwp,npp)
    character(*),intent(in)::stage
    type(s_dg_wpw_production_context),intent(in)::ctx
    real(8),intent(inout)::stage_clock
    integer,intent(in)::nwp,npp
    real(8)::now
    call cpu_time(now)
    write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,f10.3)')'[DG-WPW-BUILDER-STAGE] stage=',trim(stage),&
      ' rank=',ctx%rank_id,' nwp=',nwp,' npp=',npp,' cpu_seconds=',now-stage_clock
    flush(6);stage_clock=now
  end subroutine
  subroutine build_bounded_operator(ctx,hb,sb,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    type(s_dg_wpw_sparse_blocks),intent(in)::hb,sb
    integer,intent(out)::info
    integer::i,j,k,nlocal,ncross,nprojector_cross,nww,npeer
    integer(8)::fingerprint
    integer,allocatable::peers(:),ww_r(:),ww_c(:)
    complex(8),allocatable::ww_h(:),ww_s(:),ww_h0(:),ww_interface(:),&
      wp_h0(:),wp_interface(:),pp_h0(:),pp_interface(:)
    real(8),allocatable::kinetic(:,:),potential(:,:),nonlocal(:,:),face_self(:,:),cross_value(:)
    integer,allocatable::cross_row(:),cross_col(:)
    character(32)::metric
    if(ctx%pending_ww_valid)then
      kinetic=ctx%pending_ww_kinetic;potential=ctx%pending_ww_potential
      nonlocal=ctx%pending_ww_nonlocal;face_self=ctx%pending_ww_face_self
      cross_value=ctx%pending_ww_cross_value;cross_row=ctx%pending_ww_cross_row_id
      cross_col=ctx%pending_ww_cross_col_id;metric=ctx%pending_ww_metric_convention
      fingerprint=ctx%pending_ww_provenance_fingerprint
    else
      kinetic=ctx%ww_kinetic;potential=ctx%ww_potential;nonlocal=ctx%ww_nonlocal;face_self=ctx%ww_face_self
      cross_value=ctx%ww_cross_value;cross_row=ctx%ww_cross_row_id;cross_col=ctx%ww_cross_col_id
      metric=ctx%ww_metric_convention;fingerprint=ctx%ww_provenance_fingerprint
    endif
    nlocal=size(ctx%owned_w_ids);ncross=size(cross_value);nprojector_cross=0
    if(ctx%ww_projector_nonlocal_valid)nprojector_cross=size(ctx%ww_projector_cross_value)
    nww=nlocal*nlocal+2*ncross+nprojector_cross
    npeer=count(ctx%support_fragment_ids/=ctx%owned_fragment_id);allocate(peers(npeer));k=0
    do i=1,size(ctx%support_fragment_ids)
      if(ctx%support_fragment_ids(i)==ctx%owned_fragment_id)cycle
      k=k+1;peers(k)=ctx%support_fragment_ids(i)-1
    enddo
    allocate(ww_r(nww),ww_c(nww),ww_h(nww),ww_s(nww),ww_h0(nww),ww_interface(nww));k=0
    do j=1,nlocal;do i=1,nlocal
      k=k+1;ww_r(k)=ctx%owned_w_ids(i);ww_c(k)=ctx%owned_w_ids(j)
      if(ctx%ww_projector_nonlocal_valid)then
        ww_h0(k)=cmplx(kinetic(i,j)+potential(i,j),0d0,8)+ctx%ww_projector_nonlocal(i,j)
      else
        ww_h0(k)=cmplx(kinetic(i,j)+potential(i,j)+nonlocal(i,j),0d0,8)
      endif
      ww_interface(k)=cmplx(face_self(i,j),0d0,8);ww_h(k)=ww_h0(k)+ww_interface(k)
      ww_s(k)=merge((1d0,0d0),(0d0,0d0),i==j)
    enddo;enddo
    do i=1,ncross
      k=k+1;ww_r(k)=cross_row(i);ww_c(k)=cross_col(i)
      ww_h0(k)=0;ww_interface(k)=cmplx(cross_value(i),0d0,8)
      ww_h(k)=ww_interface(k);ww_s(k)=0
      k=k+1;ww_r(k)=cross_col(i);ww_c(k)=cross_row(i)
      ww_h0(k)=0;ww_interface(k)=cmplx(cross_value(i),0d0,8)
      ww_h(k)=ww_interface(k);ww_s(k)=0
    enddo
    if(ctx%ww_projector_nonlocal_valid)then
      do i=1,nprojector_cross
        k=k+1;ww_r(k)=ctx%ww_projector_cross_row_id(i);ww_c(k)=ctx%ww_projector_cross_col_id(i)
        ww_h0(k)=ctx%ww_projector_cross_value(i);ww_interface(k)=0
        ww_h(k)=ww_h0(k);ww_s(k)=0
      enddo
    endif
    allocate(wp_h0(size(hb%wp_values,1)),wp_interface(size(hb%wp_values,1)),&
      pp_h0(size(hb%pp_values,1)),pp_interface(size(hb%pp_values,1)))
    wp_h0=ctx%wp_h_volume+ctx%wp_h_nonlocal;wp_interface=ctx%wp_h_face
    pp_h0=ctx%pp_h_volume+ctx%pp_h_nonlocal;pp_interface=0
    fingerprint=ieor(fingerprint,int(ctx%halo_epoch,8))
    do i=1,size(hb%wp_w_row_ids)
      call mix64(fingerprint,int(hb%wp_w_row_ids(i),8));call mix64(fingerprint,int(hb%wp_pw_col_ids(i),8))
      call mix_complex(fingerprint,hb%wp_values(i,1));call mix_complex(fingerprint,sb%wp_values(i,1))
    enddo
    do i=1,size(hb%pp_pw_row_ids)
      call mix64(fingerprint,int(hb%pp_pw_row_ids(i),8));call mix64(fingerprint,int(hb%pp_pw_col_ids(i),8))
      call mix_complex(fingerprint,hb%pp_values(i,1));call mix_complex(fingerprint,sb%pp_values(i,1))
    enddo
    do i=1,size(ctx%owned_w_ids);call mix64(fingerprint,int(ctx%owned_w_ids(i),8));enddo
    do i=1,size(ctx%owned_column_ids);call mix64(fingerprint,int(ctx%owned_column_ids(i),8));enddo
    do i=1,size(ctx%support_w_ids);call mix64(fingerprint,int(ctx%support_w_ids(i),8));enddo
    do i=1,size(ctx%support_column_ids);call mix64(fingerprint,int(ctx%support_column_ids(i),8));enddo
    if(fingerprint==0)fingerprint=1
    call initialize_dg_wpw_bounded_operator(ctx%bounded_operator,ctx%comm,ctx%halo_epoch,fingerprint,&
      'fragment_root_v1',metric,'windowed_kg_sipg_v1',peers,ctx%owned_w_ids,&
      ctx%owned_column_ids,ctx%support_w_ids,ctx%support_column_ids,ww_r,ww_c,ww_h,ww_s,&
      hb%wp_w_row_ids,hb%wp_pw_col_ids,hb%wp_values(:,1),sb%wp_values(:,1),&
      hb%pp_pw_row_ids,hb%pp_pw_col_ids,hb%pp_values(:,1),sb%pp_values(:,1),info,&
      ww_h0=ww_h0,ww_interface=ww_interface,wp_h0=wp_h0,wp_interface=wp_interface,&
      pp_h0=pp_h0,pp_interface=pp_interface)
  end subroutine
  subroutine commit_pending_ww(ctx)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    ctx%ww_kinetic=ctx%pending_ww_kinetic;ctx%ww_potential=ctx%pending_ww_potential
    ctx%ww_nonlocal=ctx%pending_ww_nonlocal;ctx%ww_face_self=ctx%pending_ww_face_self
    ctx%ww_cross_face_id=ctx%pending_ww_cross_face_id;ctx%ww_cross_row_id=ctx%pending_ww_cross_row_id
    ctx%ww_cross_col_id=ctx%pending_ww_cross_col_id;ctx%ww_cross_axis=ctx%pending_ww_cross_axis
    ctx%ww_cross_side=ctx%pending_ww_cross_side;ctx%ww_cross_image=ctx%pending_ww_cross_image
    ctx%ww_cross_value=ctx%pending_ww_cross_value;ctx%ww_metric_convention=ctx%pending_ww_metric_convention
    ctx%ww_provenance_fingerprint=ctx%pending_ww_provenance_fingerprint;ctx%ww_valid=.true.
    ctx%pending_ww_valid=.false.
  end subroutine
  subroutine mix_complex(hash,value)
    integer(8),intent(inout)::hash
    complex(8),intent(in)::value
    call mix64(hash,transfer(real(value,8),hash));call mix64(hash,transfer(aimag(value),hash))
  end subroutine
  subroutine mix64(hash,value)
    integer(8),intent(inout)::hash
    integer(8),intent(in)::value
    hash=ieor(hash,value);hash=ieor(hash,shiftl(hash,13));hash=ieor(hash,shiftr(hash,7))
  end subroutine
  subroutine copy_blocks(ctx,hb,sb)
    type(s_dg_wpw_production_context),intent(inout)::ctx;type(s_dg_wpw_sparse_blocks),intent(in)::hb,sb
    if(allocated(ctx%h_wp_w))deallocate(ctx%h_wp_w,ctx%h_wp_p,ctx%h_wp,ctx%h_pp_r,ctx%h_pp_c,ctx%h_pp)
    if(allocated(ctx%s_wp_w))deallocate(ctx%s_wp_w,ctx%s_wp_p,ctx%s_wp,ctx%s_pp_r,ctx%s_pp_c,ctx%s_pp)
    allocate(ctx%h_wp_w(size(hb%wp_w_row_ids)),ctx%h_wp_p(size(hb%wp_pw_col_ids)),ctx%h_wp(size(hb%wp_values,1)))
    allocate(ctx%h_pp_r(size(hb%pp_pw_row_ids)),ctx%h_pp_c(size(hb%pp_pw_col_ids)),ctx%h_pp(size(hb%pp_values,1)))
    allocate(ctx%s_wp_w(size(sb%wp_w_row_ids)),ctx%s_wp_p(size(sb%wp_pw_col_ids)),ctx%s_wp(size(sb%wp_values,1)))
    allocate(ctx%s_pp_r(size(sb%pp_pw_row_ids)),ctx%s_pp_c(size(sb%pp_pw_col_ids)),ctx%s_pp(size(sb%pp_values,1)))
    ctx%h_wp_w=hb%wp_w_row_ids;ctx%h_wp_p=hb%wp_pw_col_ids;ctx%h_wp=hb%wp_values(:,1)
    ctx%h_pp_r=hb%pp_pw_row_ids;ctx%h_pp_c=hb%pp_pw_col_ids;ctx%h_pp=hb%pp_values(:,1)
    ctx%s_wp_w=sb%wp_w_row_ids;ctx%s_wp_p=sb%wp_pw_col_ids;ctx%s_wp=sb%wp_values(:,1)
    ctx%s_pp_r=sb%pp_pw_row_ids;ctx%s_pp_c=sb%pp_pw_col_ids;ctx%s_pp=sb%pp_values(:,1)
  end subroutine
  subroutine bind_dg_wpw_hs_callbacks(ctx,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx;integer,intent(out)::info
    info=0
    if(.not.ctx%operator_valid.or.ctx%operator_epoch/=ctx%halo_epoch)then;info=1;return;endif
    ctx%callbacks_bound=.true.
  end subroutine
end module
