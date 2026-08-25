module dg_hybrid_projected_fragment_pipeline
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_windowed_pw_basis,only:materialize_dg_hybrid_windowed_pw_columns
  use dg_hybrid_wannier_complement,only:compute_dg_hybrid_wannier_projection_tile,&
    materialize_dg_hybrid_projected_pw_tile
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_basis_stream,only:s_dg_hybrid_fragment_basis_stream,&
    initialize_dg_hybrid_fragment_basis_stream,append_dg_hybrid_projected_pw_tile,&
    finalize_dg_hybrid_fragment_basis_stream
  implicit none
  private
  public::build_dg_hybrid_projected_fragment_basis
contains
  subroutine build_dg_hybrid_projected_fragment_basis(comm,global_point_count,fragment_count,&
      fragment_id,core_ids,weights,core_wannier,core_coordinates,core_windows,buffer_ids,&
      buffer_wannier,buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,tile_width,&
      tolerance,wannier_fingerprint,basis,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_id,tile_width,wannier_owner(:)
    integer(int64),intent(in)::core_ids(:),buffer_ids(:),wannier_fingerprint
    real(real64),intent(in)::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),&
      buffer_windows(:,:),g_vectors(:,:),tolerance
    complex(real64),intent(in)::core_wannier(:,:),buffer_wannier(:,:)
    type(s_dg_hybrid_basis_catalog),intent(in)::catalog
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_fragment_basis_stream)::stream
    integer::packet,g,npw,column,width,i
    integer,allocatable::pw_owner(:)
    real(real64),allocatable::normalized_buffer_windows(:,:)
    complex(real64),allocatable::core_raw(:,:),buffer_raw(:,:),coefficients(:,:),projected_buffer(:,:)
    integer(int64)::stage_workspace,stage_fingerprint,stream_workspace,stream_fingerprint
    real(real64)::scale
    logical::stage_ok
    character(256)::stage_message
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;npw=0
    if(.not.catalog%valid.or.tile_width<1.or.size(core_wannier,1)/=size(wannier_owner).or.&
        size(core_wannier,2)/=size(core_ids).or.size(buffer_wannier,1)/=size(wannier_owner).or.&
        size(buffer_wannier,2)/=size(buffer_ids).or.size(weights)/=size(core_ids).or.&
        any(shape(core_coordinates)/=[3,size(core_ids)]).or.any(shape(buffer_coordinates)/=[3,size(buffer_ids)]).or.&
        size(core_windows,2)/=size(core_ids).or.size(buffer_windows,2)/=size(buffer_ids).or.&
        size(core_windows,1)/=fragment_count.or.size(buffer_windows,1)/=fragment_count)then
      message='invalid projected fragment pipeline shape';return
    endif
    do packet=1,size(catalog%packets);npw=npw+size(catalog%packets(packet)%g_indices);enddo
    if(npw<1)then;message='empty projected fragment PW catalog';return;endif
    allocate(pw_owner(npw))
    allocate(normalized_buffer_windows,source=buffer_windows)
    column=0
    do packet=1,size(catalog%packets)
      if(catalog%packets(packet)%fragment_id<1.or.catalog%packets(packet)%fragment_id>fragment_count)then
        message='invalid projected fragment packet owner';return
      endif
      do g=1,size(catalog%packets(packet)%g_indices)
        column=column+1;pw_owner(column)=catalog%packets(packet)%fragment_id
      enddo
    enddo
    do i=1,size(buffer_ids)
      scale=sqrt(sum(normalized_buffer_windows(:,i)**2))
      if(.not.ieee_is_finite(scale).or.scale<=0d0)then
        message='fragment buffer lies outside all PW windows';return
      endif
      normalized_buffer_windows(:,i)=normalized_buffer_windows(:,i)/scale
    enddo
    call initialize_dg_hybrid_fragment_basis_stream(comm,fragment_count,fragment_id,buffer_ids,&
      buffer_wannier,wannier_owner,pw_owner,stream,basis,stream_workspace,stream_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do column=1,npw,tile_width
      width=min(tile_width,npw-column+1)
      call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,core_coordinates,core_windows,&
        column,width,core_raw,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      call compute_dg_hybrid_wannier_projection_tile(comm,global_point_count,core_ids,weights,core_wannier,&
        core_raw,wannier_fingerprint,catalog%packet_fingerprint,column,tolerance,coefficients,&
        stage_workspace,stage_fingerprint,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      workspace_peak_bytes=max(workspace_peak_bytes,stage_workspace)
      call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,buffer_coordinates,&
        normalized_buffer_windows,column,width,buffer_raw,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      call materialize_dg_hybrid_projected_pw_tile(buffer_wannier,buffer_raw,coefficients,&
        projected_buffer,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      call append_dg_hybrid_projected_pw_tile(stream,column,projected_buffer,basis,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      workspace_peak_bytes=max(workspace_peak_bytes,16_int64*int(size(core_raw)+size(buffer_raw)+&
        size(coefficients)+size(projected_buffer),int64))
      fingerprint=ieor(fingerprint,ishftc(stage_fingerprint,modulo(column,63)))
      deallocate(core_raw,buffer_raw,coefficients,projected_buffer)
    enddo
    call finalize_dg_hybrid_fragment_basis_stream(comm,stream,basis,stream_workspace,&
      stream_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    workspace_peak_bytes=max(workspace_peak_bytes,stream_workspace)
    fingerprint=ieor(fingerprint,stream_fingerprint)
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(pw_owner,normalized_buffer_windows);ok=.true.
  end subroutine build_dg_hybrid_projected_fragment_basis
end module dg_hybrid_projected_fragment_pipeline
