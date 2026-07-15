module rt_dg_wpw_grid_point_adapter
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_wpw_window, only: map_global_to_phi_box_coord_pw
  use rt_dg_wpw_point_evaluator, only: evaluate_windowed_kg_point, evaluate_wannier_point
  implicit none
  private
  public :: evaluate_local_wannier_grid_point
  public :: evaluate_windowed_kg_grid_point

contains

  subroutine evaluate_windowed_kg_grid_point(dg_frag,k_fragment,g_mode_id,grid,omega_cell,value,gradient,info)
    type(s_dg_fragment_rt),intent(in)::dg_frag
    integer,intent(in)::k_fragment,g_mode_id,grid(3)
    real(8),intent(in)::omega_cell
    complex(8),intent(out)::value,gradient(3)
    integer,intent(out)::info
    integer::ilocal,nlocal,ix,iy,iz,window_extent(3)
    real(8)::chi,grad_chi(3),position(3)
    value=(0d0,0d0); gradient=(0d0,0d0); info=0
    if(dg_frag%ifrag_start<1 .or. dg_frag%ifrag_end<dg_frag%ifrag_start .or. &
       dg_frag%ifrag_end>dg_frag%n_frag .or. k_fragment<1 .or. k_fragment>dg_frag%n_frag) then
      info=1
      return
    end if
    nlocal=dg_frag%ifrag_end-dg_frag%ifrag_start+1
    ilocal=k_fragment-dg_frag%ifrag_start+1
    if(.not.dg_frag%has_wpw_window .or. ilocal<1 .or. ilocal>nlocal .or. &
       .not.allocated(dg_frag%k_pw) .or. .not.allocated(dg_frag%wpw_window_box_lo) .or. &
       .not.allocated(dg_frag%wpw_window_box_hi) .or. .not.allocated(dg_frag%wpw_chi) .or. &
       .not.allocated(dg_frag%wpw_grad_chi) .or. any(dg_frag%hgs<=0d0)) then
      info=1
      return
    end if
    if(size(dg_frag%k_pw,1)/=3 .or. g_mode_id<1 .or. g_mode_id>size(dg_frag%k_pw,2)) then
      info=1
      return
    end if
    ! Only precomputed windows owned by this rank are admissible here.  A remote
    ! K trace must be supplied by the later support-halo/provider layer; never
    ! reconstruct it by scanning all fragments at every quadrature point.
    if(size(dg_frag%wpw_window_box_lo,1)/=3 .or. size(dg_frag%wpw_window_box_hi,1)/=3 .or. &
       size(dg_frag%wpw_window_box_lo,2)/=nlocal .or. size(dg_frag%wpw_window_box_hi,2)/=nlocal .or. &
       size(dg_frag%wpw_chi,4)/=nlocal .or. size(dg_frag%wpw_grad_chi,1)/=3 .or. &
       size(dg_frag%wpw_grad_chi,5)/=nlocal .or. &
       size(dg_frag%wpw_chi,1)/=size(dg_frag%wpw_grad_chi,2) .or. &
       size(dg_frag%wpw_chi,2)/=size(dg_frag%wpw_grad_chi,3) .or. &
       size(dg_frag%wpw_chi,3)/=size(dg_frag%wpw_grad_chi,4) .or. &
       any(dg_frag%wpw_window_box_hi(:,ilocal)<dg_frag%wpw_window_box_lo(:,ilocal)) .or. &
       any(grid<dg_frag%wpw_window_box_lo(:,ilocal)) .or. &
       any(grid>dg_frag%wpw_window_box_hi(:,ilocal))) then
      info=1
      return
    end if
    window_extent=dg_frag%wpw_window_box_hi(:,ilocal)-dg_frag%wpw_window_box_lo(:,ilocal)+1
    if(window_extent(1)>size(dg_frag%wpw_chi,1) .or. &
       window_extent(2)>size(dg_frag%wpw_chi,2) .or. &
       window_extent(3)>size(dg_frag%wpw_chi,3) .or. &
       window_extent(1)>size(dg_frag%wpw_grad_chi,2) .or. &
       window_extent(2)>size(dg_frag%wpw_grad_chi,3) .or. &
       window_extent(3)>size(dg_frag%wpw_grad_chi,4)) then
      info=1
      return
    end if
    ix=grid(1)-dg_frag%wpw_window_box_lo(1,ilocal)+1
    iy=grid(2)-dg_frag%wpw_window_box_lo(2,ilocal)+1
    iz=grid(3)-dg_frag%wpw_window_box_lo(3,ilocal)+1
    if(ix<1 .or. iy<1 .or. iz<1 .or. ix>size(dg_frag%wpw_chi,1) .or. &
       iy>size(dg_frag%wpw_chi,2) .or. iz>size(dg_frag%wpw_chi,3) .or. &
       ix>size(dg_frag%wpw_grad_chi,2) .or. iy>size(dg_frag%wpw_grad_chi,3) .or. &
       iz>size(dg_frag%wpw_grad_chi,4)) then
      info=1
      return
    end if
    chi=dg_frag%wpw_chi(ix,iy,iz,ilocal)
    grad_chi=dg_frag%wpw_grad_chi(:,ix,iy,iz,ilocal)
    position=dble(grid)*dg_frag%hgs
    call evaluate_windowed_kg_point(chi,grad_chi,dg_frag%k_pw(:,g_mode_id),position,omega_cell,value,gradient,info)
  end subroutine evaluate_windowed_kg_grid_point

  subroutine evaluate_local_wannier_grid_point(dg_frag,ifrag,ispin,w_row_ids,grid,values,gradients,info)
    type(s_dg_fragment_rt),intent(in)::dg_frag
    integer,intent(in)::ifrag,ispin,w_row_ids(:),grid(3)
    complex(8),intent(out)::values(:),gradients(:,:)
    integer,intent(out)::info
    integer::ilocal,iw,slot,ib,nbf,bx,by,bz,ix,iy,iz
    integer::phi_lo(3),phi_hi(3)
    integer,allocatable::slots(:)
    complex(8),allocatable::coef(:,:),basis_values(:),basis_gradients(:,:)

    values=(0d0,0d0); gradients=(0d0,0d0); info=0
    ilocal=ifrag-dg_frag%ifrag_start+1
    if(.not.dg_frag%has_global_wannier_local_basis .or. .not.dg_frag%gradient_basis_cache_valid .or. &
       ilocal<1 .or. ifrag>dg_frag%ifrag_end .or. ispin<1 .or. ispin>dg_frag%nspin .or. &
       size(values)/=size(w_row_ids) .or. any(shape(gradients)/=[3,size(w_row_ids)]) .or. &
       .not.strictly_increasing(w_row_ids)) then
      info=1
      return
    end if
    if(.not.allocated(dg_frag%n_basis) .or. .not.allocated(dg_frag%ixyz_frag) .or. &
       .not.allocated(dg_frag%nxyz_domain) .or. .not.allocated(dg_frag%global_wannier_local_ids) .or. &
       .not.allocated(dg_frag%global_wannier_local_nkeep) .or. &
       .not.allocated(dg_frag%global_wannier_local_coef) .or. &
       .not.allocated(dg_frag%gradient_basis_cache) .or. &
       (.not.allocated(dg_frag%phi_frag) .and. .not.allocated(dg_frag%phi_frag_c))) then
      info=1
      return
    end if
    if(.not.grid_context_extents_ready(dg_frag,ifrag,ispin,ilocal)) then
      info=1
      return
    end if
    allocate(slots(size(w_row_ids)))
    do iw=1,size(w_row_ids)
      slots(iw)=0
      do slot=1,min(dg_frag%global_wannier_local_nkeep(ilocal),size(dg_frag%global_wannier_local_ids,1))
        if(dg_frag%global_wannier_local_ids(slot,ilocal)==w_row_ids(iw)) then
          slots(iw)=slot
          exit
        end if
      end do
      if(slots(iw)==0) then; info=2; return; end if
    end do
    if(allocated(dg_frag%phi_frag_c)) then
      nbf=min(dg_frag%n_basis(ifrag,ispin),size(dg_frag%phi_frag_c,4),size(dg_frag%global_wannier_local_coef,1))
      phi_lo=[lbound(dg_frag%phi_frag_c,1),lbound(dg_frag%phi_frag_c,2),lbound(dg_frag%phi_frag_c,3)]
      phi_hi=[ubound(dg_frag%phi_frag_c,1),ubound(dg_frag%phi_frag_c,2),ubound(dg_frag%phi_frag_c,3)]
      bx=map_global_to_phi_box_coord_pw(grid(1),lbound(dg_frag%phi_frag_c,1),ubound(dg_frag%phi_frag_c,1),dg_frag%lgnum_total(1))
      by=map_global_to_phi_box_coord_pw(grid(2),lbound(dg_frag%phi_frag_c,2),ubound(dg_frag%phi_frag_c,2),dg_frag%lgnum_total(2))
      bz=map_global_to_phi_box_coord_pw(grid(3),lbound(dg_frag%phi_frag_c,3),ubound(dg_frag%phi_frag_c,3),dg_frag%lgnum_total(3))
    else
      nbf=min(dg_frag%n_basis(ifrag,ispin),size(dg_frag%phi_frag,4),size(dg_frag%global_wannier_local_coef,1))
      phi_lo=[lbound(dg_frag%phi_frag,1),lbound(dg_frag%phi_frag,2),lbound(dg_frag%phi_frag,3)]
      phi_hi=[ubound(dg_frag%phi_frag,1),ubound(dg_frag%phi_frag,2),ubound(dg_frag%phi_frag,3)]
      bx=map_global_to_phi_box_coord_pw(grid(1),lbound(dg_frag%phi_frag,1),ubound(dg_frag%phi_frag,1),dg_frag%lgnum_total(1))
      by=map_global_to_phi_box_coord_pw(grid(2),lbound(dg_frag%phi_frag,2),ubound(dg_frag%phi_frag,2),dg_frag%lgnum_total(2))
      bz=map_global_to_phi_box_coord_pw(grid(3),lbound(dg_frag%phi_frag,3),ubound(dg_frag%phi_frag,3),dg_frag%lgnum_total(3))
    end if
    ix=grid(1)-dg_frag%ixyz_frag(1,ifrag)+1
    iy=grid(2)-dg_frag%ixyz_frag(2,ifrag)+1
    iz=grid(3)-dg_frag%ixyz_frag(3,ifrag)+1
    if(nbf<=0 .or. bx<phi_lo(1) .or. bx>phi_hi(1) .or. by<phi_lo(2) .or. by>phi_hi(2) .or. &
       bz<phi_lo(3) .or. bz>phi_hi(3) .or. ix<1 .or. iy<1 .or. iz<1 .or. &
       ix>size(dg_frag%gradient_basis_cache,1) .or. iy>size(dg_frag%gradient_basis_cache,2) .or. &
       iz>size(dg_frag%gradient_basis_cache,3)) then; info=3; return; end if
    allocate(coef(nbf,size(w_row_ids)),basis_values(nbf),basis_gradients(3,nbf))
    do iw=1,size(w_row_ids)
      coef(:,iw)=dg_frag%global_wannier_local_coef(1:nbf,slots(iw),ispin,ilocal)
    end do
    do ib=1,nbf
      if(allocated(dg_frag%phi_frag_c)) then
        basis_values(ib)=dg_frag%phi_frag_c(bx,by,bz,ib,ilocal)
      else
        basis_values(ib)=cmplx(dg_frag%phi_frag(bx,by,bz,ib,ilocal),0d0,kind=8)
      end if
      basis_gradients(:,ib)=cmplx(dg_frag%gradient_basis_cache(ix,iy,iz,:,ib,ilocal),0d0,kind=8)
    end do
    call evaluate_wannier_point(coef,basis_values,basis_gradients,values,gradients,info)
  end subroutine evaluate_local_wannier_grid_point

  logical function strictly_increasing(ids) result(ok)
    integer,intent(in)::ids(:)
    integer::i
    ok=.true.
    do i=2,size(ids)
      if(ids(i)<=ids(i-1)) then; ok=.false.; return; end if
    end do
  end function strictly_increasing

  logical function grid_context_extents_ready(dg_frag,ifrag,ispin,ilocal) result(ok)
    type(s_dg_fragment_rt),intent(in)::dg_frag
    integer,intent(in)::ifrag,ispin,ilocal
    ok=.false.
    if(size(dg_frag%n_basis,1)<ifrag .or. size(dg_frag%n_basis,2)<ispin) return
    if(size(dg_frag%ixyz_frag,1)<3 .or. size(dg_frag%ixyz_frag,2)<ifrag) return
    if(size(dg_frag%nxyz_domain,1)<3 .or. size(dg_frag%nxyz_domain,2)<ifrag) return
    if(any(dg_frag%lgnum_total<=0)) return
    if(size(dg_frag%global_wannier_local_nkeep)<ilocal) return
    if(size(dg_frag%global_wannier_local_ids,2)<ilocal) return
    if(size(dg_frag%global_wannier_local_coef,2)<max(1,dg_frag%global_wannier_local_nkeep(ilocal))) return
    if(size(dg_frag%global_wannier_local_coef,3)<ispin .or. &
       size(dg_frag%global_wannier_local_coef,4)<ilocal) return
    if(size(dg_frag%gradient_basis_cache,4)<3 .or. &
       size(dg_frag%gradient_basis_cache,5)<dg_frag%n_basis(ifrag,ispin) .or. &
       size(dg_frag%gradient_basis_cache,6)<ilocal) return
    ok=.true.
    if(allocated(dg_frag%phi_frag_c)) then
      ok=size(dg_frag%phi_frag_c,4)>=dg_frag%n_basis(ifrag,ispin) .and. size(dg_frag%phi_frag_c,5)>=ilocal
    else
      ok=size(dg_frag%phi_frag,4)>=dg_frag%n_basis(ifrag,ispin) .and. size(dg_frag%phi_frag,5)>=ilocal
    end if
  end function grid_context_extents_ready
end module rt_dg_wpw_grid_point_adapter
