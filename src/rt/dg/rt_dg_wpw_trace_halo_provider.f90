module rt_dg_wpw_trace_halo_provider
  use dg_wpw_production_context,only:s_dg_wpw_production_context
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider,bind_wpw_face_trace_provider
  implicit none
  private
  type,public::s_dg_wpw_trace_halo_state
    integer::epoch=-1,k_minus=0,k_plus=0,axis=0,side=0
    logical::complete=.false.
    integer,allocatable::grid(:,:)
    complex(8),allocatable::wm(:,:),wp(:,:),gwm(:,:,:),gwp(:,:,:)
    complex(8),allocatable::pm(:,:),pp(:,:),gpm(:,:,:),gpp(:,:,:)
  end type
  type,public::s_dg_wpw_trace_halo_set
    integer::epoch=-1,nface=0
    type(s_dg_wpw_trace_halo_state),allocatable::faces(:)
  end type
  public::prepare_dg_wpw_trace_halo
  public::prepare_dg_wpw_trace_halo_face
  public::reduce_dg_wpw_face_trace_parts
contains
  subroutine reduce_dg_wpw_face_trace_parts(comm,root,ctx,set,provider,epoch,k_minus,k_plus,axis,&
      side,grid,coverage,wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
    use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Reduce,MPI_Bcast,MPI_INTEGER,&
      MPI_MAX,MPI_SUM,MPI_DOUBLE_COMPLEX,MPI_SUCCESS
    integer,intent(in)::comm,root,epoch,k_minus,k_plus,axis,side,grid(:,:),coverage(:)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    type(s_dg_wpw_trace_halo_set),target,intent(inout)::set
    type(s_wpw_face_trace_provider),intent(inout)::provider
    complex(8),intent(in)::wm(:,:),wp(:,:),gwm(:,:,:),gwp(:,:,:)
    complex(8),intent(in)::pm(:,:),pp(:,:),gpm(:,:,:),gpp(:,:,:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,npoint,root_info
    integer::local_shape(4),reference_shape(4)
    integer,allocatable::reference_grid(:,:),total_coverage(:)
    complex(8),allocatable::rwm(:,:),rwp(:,:),rgwm(:,:,:),rgwp(:,:,:)
    complex(8),allocatable::rpm(:,:),rpp(:,:),rgpm(:,:,:),rgpp(:,:,:)

    info=1;local_bad=0;npoint=size(grid,2)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank.or.epoch<=0.or.k_minus<=0.or.k_minus>=k_plus.or.&
       axis<1.or.axis>3.or.abs(side)/=1.or.size(grid,1)/=3.or.npoint<=0.or.&
       size(coverage)/=npoint.or.any(coverage<0).or.any(coverage>1).or.&
       size(wm,2)/=npoint.or.any(shape(wp)/=shape(wm)).or.&
       any(shape(gwm)/=[3,size(wm,1),npoint]).or.any(shape(gwp)/=shape(gwm)).or.&
       size(pm,2)/=npoint.or.any(shape(pp)/=shape(pm)).or.&
       any(shape(gpm)/=[3,size(pm,1),npoint]).or.any(shape(gpp)/=shape(gpm)).or.&
       .not.all(finite_complex_3(wm)).or..not.all(finite_complex_3(wp)).or.&
       .not.all(finite_complex_3(gwm)).or..not.all(finite_complex_3(gwp)).or.&
       .not.all(finite_complex_3(pm)).or..not.all(finite_complex_3(pp)).or.&
       .not.all(finite_complex_3(gpm)).or..not.all(finite_complex_3(gpp)))local_bad=1
    local_shape=[size(grid,1),npoint,size(wm,1),size(pm,1)];reference_shape=local_shape
    call MPI_Bcast(reference_shape,4,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(reference_shape/=local_shape))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(reference_grid(3,npoint));if(rank==root)reference_grid=grid
    call MPI_Bcast(reference_grid,size(reference_grid),MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(reference_grid/=grid))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return

    allocate(total_coverage(npoint));allocate(rwm,source=wm);allocate(rwp,source=wp)
    allocate(rgwm,source=gwm);allocate(rgwp,source=gwp)
    allocate(rpm,source=pm);allocate(rpp,source=pp)
    allocate(rgpm,source=gpm);allocate(rgpp,source=gpp)
    call MPI_Reduce(coverage,total_coverage,npoint,MPI_INTEGER,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(wm,rwm,size(wm),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(wp,rwp,size(wp),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(gwm,rgwm,size(gwm),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(gwp,rgwp,size(gwp),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(pm,rpm,size(pm),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(pp,rpp,size(pp),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(gpm,rgpm,size(gpm),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(gpp,rgpp,size(gpp),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    root_info=1
    if(rank==root.and.ierr==MPI_SUCCESS.and.all(total_coverage==1))then
      call prepare_dg_wpw_trace_halo_face(ctx,set,provider,epoch,k_minus,k_plus,axis,side,grid,&
        rwm,rwp,rgwm,rgwp,rpm,rpp,rgpm,rgpp,root_info)
    endif
    call MPI_Bcast(root_info,1,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.root_info/=0)return
    info=0
  end subroutine reduce_dg_wpw_face_trace_parts

  elemental logical function finite_complex_3(value)result(ok)
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function finite_complex_3

  subroutine prepare_dg_wpw_trace_halo_face(ctx,set,provider,epoch,k_minus,k_plus,axis,side,grid,&
      wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    type(s_dg_wpw_trace_halo_set),target,intent(inout)::set
    type(s_wpw_face_trace_provider),intent(inout)::provider
    integer,intent(in)::epoch,k_minus,k_plus,axis,side,grid(:,:)
    complex(8),intent(in)::wm(:,:),wp(:,:),gwm(:,:,:),gwp(:,:,:)
    complex(8),intent(in)::pm(:,:),pp(:,:),gpm(:,:,:),gpp(:,:,:)
    integer,intent(out)::info
    type(s_dg_wpw_trace_halo_state)::face
    type(s_dg_wpw_trace_halo_state),allocatable::work(:)
    class(*),pointer::bound_context
    integer::npoint,i

    info=1;npoint=size(grid,2)
    if(epoch<ctx%halo_epoch.or.k_minus<=0.or.k_minus>=k_plus.or.axis<1.or.axis>3.or.abs(side)/=1.or.&
       size(grid,1)/=3.or.npoint<=0.or.size(wm,2)/=npoint.or.any(shape(wp)/=shape(wm)).or.&
       any(shape(gwm)/=[3,size(wm,1),npoint]).or.any(shape(gwp)/=shape(gwm)).or.&
       size(pm,2)/=npoint.or.any(shape(pp)/=shape(pm)).or.&
       any(shape(gpm)/=[3,size(pm,1),npoint]).or.any(shape(gpp)/=shape(gpm)))return
    if(epoch>set%epoch)then
      if(allocated(set%faces))deallocate(set%faces)
      set%epoch=epoch;set%nface=0
    else if(epoch<set%epoch)then
      return
    endif
    do i=1,set%nface
      if(set%faces(i)%k_minus==k_minus.and.set%faces(i)%k_plus==k_plus.and.&
         set%faces(i)%axis==axis.and.set%faces(i)%side==side)return
    enddo
    face%epoch=epoch;face%k_minus=k_minus;face%k_plus=k_plus;face%axis=axis;face%side=side
    face%complete=.true.;face%grid=grid;face%wm=wm;face%wp=wp;face%gwm=gwm;face%gwp=gwp
    face%pm=pm;face%pp=pp;face%gpm=gpm;face%gpp=gpp
    allocate(work(set%nface+1))
    if(set%nface>0)work(1:set%nface)=set%faces
    work(set%nface+1)=face;call move_alloc(work,set%faces);set%nface=set%nface+1
    bound_context=>set
    call bind_wpw_face_trace_provider(provider,bound_context,read_prepared_trace_set,info)
    if(info/=0)return
    ctx%halo_epoch=epoch;info=0
  end subroutine

  subroutine prepare_dg_wpw_trace_halo(ctx,halo,provider,epoch,k_minus,k_plus,axis,side,grid,&
      wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    type(s_dg_wpw_trace_halo_state),target,intent(inout)::halo
    type(s_wpw_face_trace_provider),intent(inout)::provider
    integer,intent(in)::epoch,k_minus,k_plus,axis,side,grid(:,:)
    complex(8),intent(in)::wm(:,:),wp(:,:),gwm(:,:,:),gwp(:,:,:)
    complex(8),intent(in)::pm(:,:),pp(:,:),gpm(:,:,:),gpp(:,:,:)
    integer,intent(out)::info
    type(s_dg_wpw_trace_halo_state)::candidate
    class(*),pointer::bound_context
    integer::npoint
    info=0;npoint=size(grid,2)
    if(epoch<=ctx%halo_epoch.or.k_minus<=0.or.k_minus>=k_plus.or.axis<1.or.axis>3.or.abs(side)/=1.or.&
       size(grid,1)/=3.or.npoint<=0.or.size(wm,2)/=npoint.or.any(shape(wp)/=shape(wm)).or.&
       any(shape(gwm)/=[3,size(wm,1),npoint]).or.any(shape(gwp)/=shape(gwm)).or.&
       size(pm,2)/=npoint.or.any(shape(pp)/=shape(pm)).or.&
       any(shape(gpm)/=[3,size(pm,1),npoint]).or.any(shape(gpp)/=shape(gpm)))then
      info=1;return
    endif
    candidate%epoch=epoch;candidate%k_minus=k_minus;candidate%k_plus=k_plus
    candidate%axis=axis;candidate%side=side;candidate%complete=.true.
    candidate%grid=grid;candidate%wm=wm;candidate%wp=wp;candidate%gwm=gwm;candidate%gwp=gwp
    candidate%pm=pm;candidate%pp=pp;candidate%gpm=gpm;candidate%gpp=gpp
    halo=candidate
    bound_context=>halo
    call bind_wpw_face_trace_provider(provider,bound_context,read_prepared_trace,info)
    if(info/=0)return
    ctx%halo_epoch=epoch
  end subroutine
  subroutine read_prepared_trace(user_context,k_minus,k_plus,axis,side,grid,&
      wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
    class(*),pointer,intent(inout)::user_context
    integer,intent(in)::k_minus,k_plus,axis,side,grid(3)
    complex(8),intent(out)::wm(:),wp(:),gwm(:,:),gwp(:,:),pm(:),pp(:),gpm(:,:),gpp(:,:)
    integer,intent(out)::info
    integer::i
    info=1
    select type(halo=>user_context)
    type is(s_dg_wpw_trace_halo_state)
      if(.not.halo%complete.or.k_minus/=halo%k_minus.or.k_plus/=halo%k_plus.or.&
         axis/=halo%axis.or.side/=halo%side) return
      do i=1,size(halo%grid,2)
        if(all(grid==halo%grid(:,i)))then
          if(size(wm)/=size(halo%wm,1).or.size(pm)/=size(halo%pm,1))return
          wm=halo%wm(:,i);wp=halo%wp(:,i);gwm=halo%gwm(:,:,i);gwp=halo%gwp(:,:,i)
          pm=halo%pm(:,i);pp=halo%pp(:,i);gpm=halo%gpm(:,:,i);gpp=halo%gpp(:,:,i)
          info=0;return
        endif
      enddo
    end select
  end subroutine

  subroutine read_prepared_trace_set(user_context,k_minus,k_plus,axis,side,grid,&
      wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
    class(*),pointer,intent(inout)::user_context
    integer,intent(in)::k_minus,k_plus,axis,side,grid(3)
    complex(8),intent(out)::wm(:),wp(:),gwm(:,:),gwp(:,:),pm(:),pp(:),gpm(:,:),gpp(:,:)
    integer,intent(out)::info
    integer::iface,i
    info=1
    select type(set=>user_context)
    type is(s_dg_wpw_trace_halo_set)
      do iface=1,set%nface
        if(.not.set%faces(iface)%complete.or.k_minus/=set%faces(iface)%k_minus.or.&
           k_plus/=set%faces(iface)%k_plus.or.axis/=set%faces(iface)%axis.or.&
           side/=set%faces(iface)%side)cycle
        do i=1,size(set%faces(iface)%grid,2)
          if(all(grid==set%faces(iface)%grid(:,i)))then
            if(size(wm)/=size(set%faces(iface)%wm,1).or.size(pm)/=size(set%faces(iface)%pm,1))return
            wm=set%faces(iface)%wm(:,i);wp=set%faces(iface)%wp(:,i)
            gwm=set%faces(iface)%gwm(:,:,i);gwp=set%faces(iface)%gwp(:,:,i)
            pm=set%faces(iface)%pm(:,i);pp=set%faces(iface)%pp(:,i)
            gpm=set%faces(iface)%gpm(:,:,i);gpp=set%faces(iface)%gpp(:,:,i)
            info=0;return
          endif
        enddo
      enddo
    end select
  end subroutine read_prepared_trace_set
end module
