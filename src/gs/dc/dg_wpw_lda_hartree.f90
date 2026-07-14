module dg_wpw_lda_hartree
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: validate_core_ownership, integrate_core_lda_terms
  public :: hartree_energy_global, update_wpw_lda_hartree
contains
  subroutine validate_core_ownership(core_owner,npoint,nfragment,info,bad_point)
    integer,intent(in) :: npoint,nfragment
    logical,intent(in) :: core_owner(npoint,nfragment)
    integer,intent(out) :: info,bad_point
    integer :: ipoint
    info=0; bad_point=0
    if(npoint<1 .or. nfragment<1)then
      info=1; return
    endif
    do ipoint=1,npoint
      if(count(core_owner(ipoint,:))/=1)then
        info=2; bad_point=ipoint; return
      endif
    enddo
  end subroutine validate_core_ownership

  subroutine integrate_core_lda_terms(rho,vxc,exc_density,weight,core_owner,npoint,nfragment,e_xc,n_vxc,info)
    integer,intent(in) :: npoint,nfragment
    real(8),intent(in) :: rho(npoint,nfragment),vxc(npoint,nfragment),exc_density(npoint,nfragment),weight(npoint)
    logical,intent(in) :: core_owner(npoint,nfragment)
    real(8),intent(out) :: e_xc,n_vxc
    integer,intent(out) :: info
    integer :: ipoint,ifragment,bad_point
    e_xc=0d0; n_vxc=0d0; info=0
    call validate_core_ownership(core_owner,npoint,nfragment,info,bad_point)
    if(info/=0)return
    if(any(weight<0d0).or.any(.not.ieee_is_finite(weight)))then
      info=3; return
    endif
    do ifragment=1,nfragment
      do ipoint=1,npoint
        if(.not.core_owner(ipoint,ifragment))cycle
        if(.not.ieee_is_finite(rho(ipoint,ifragment)).or. &
           .not.ieee_is_finite(vxc(ipoint,ifragment)).or. &
           .not.ieee_is_finite(exc_density(ipoint,ifragment)))then
          info=4; return
        endif
        e_xc=e_xc+weight(ipoint)*exc_density(ipoint,ifragment)
        n_vxc=n_vxc+weight(ipoint)*rho(ipoint,ifragment)*vxc(ipoint,ifragment)
      enddo
    enddo
  end subroutine integrate_core_lda_terms

  subroutine hartree_energy_global(rho,vh,weight,npoint,comm,e_h,info)
    use communication, only: comm_summation
    integer,intent(in) :: npoint
    real(8),intent(in) :: rho(npoint),vh(npoint),weight(npoint)
    integer,intent(in) :: comm
    real(8),intent(out) :: e_h
    integer,intent(out) :: info
    real(8) :: e_local
    integer :: local_bad,global_bad
    e_h=0d0; e_local=0d0; info=0
    if(npoint<1 .or. any(weight<0d0) .or. any(.not.ieee_is_finite(weight)) .or. &
       any(.not.ieee_is_finite(rho)) .or. any(.not.ieee_is_finite(vh)))then
      info=1
    endif
    local_bad=merge(1,0,info/=0)
    call comm_summation(local_bad,global_bad,comm)
    if(global_bad/=0)then
      if(info==0)info=3
      return
    endif
    e_local=0.5d0*sum(weight*rho*vh)
    call comm_summation(e_local,e_h,comm)
    if(.not.ieee_is_finite(e_h))info=2
  end subroutine hartree_energy_global

  subroutine update_wpw_lda_hartree(lg,mg,system,info,stencil,xc_func,pp,ppn,spsi, &
      srg,srg_scalar,poisson,fg,rho,rho_s,vh,vxc,rho_core,core_owner,npoint,nfragment, &
      e_xc,e_xc_core,n_vxc_core,status,bad_point)
    use structures
    use hartree_sub, only: hartree
    use salmon_xc, only: exchange_correlation,eexc_tmp
    use salmon_global, only: yn_hse
    use communication, only: comm_summation
    type(s_rgrid),intent(in) :: lg,mg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_xc_functional),intent(in) :: xc_func
    type(s_pp_info),intent(in) :: pp
    type(s_pp_nlcc),intent(in) :: ppn
    type(s_orbital),intent(inout) :: spsi
    type(s_sendrecv_grid),intent(inout) :: srg,srg_scalar
    type(s_poisson),intent(inout) :: poisson
    type(s_reciprocal_grid),intent(inout) :: fg
    type(s_scalar),intent(inout) :: rho,rho_s(system%nspin),vh,vxc(system%nspin)
    integer,intent(in) :: npoint,nfragment
    logical,intent(in) :: core_owner(npoint,nfragment)
    real(8),intent(in) :: rho_core(npoint,nfragment,system%nspin)
    real(8),intent(out) :: e_xc,e_xc_core,n_vxc_core
    integer,intent(out) :: status,bad_point
    integer :: ispin,ipoint,ifragment,local_bad,global_bad
    real(8) :: e_xc_local,n_vxc_local
    real(8),allocatable :: rho_flat(:)
    ! ipoint is frozen to Fortran array-element order of rho_s(:)%f.  The
    ! caller supplies per-fragment values in exactly that order. Hartree remains
    ! one global Poisson solve; XC remains SALMON's configured LDA implementation.
    e_xc=0d0; e_xc_core=0d0; n_vxc_core=0d0; status=0; bad_point=0
    if(xc_func%use_gradient .or. xc_func%use_laplacian .or. &
       xc_func%use_kinetic_energy .or. xc_func%use_current .or. yn_hse=='y')then
      status=10
    endif
    if(status==0)call validate_core_ownership(core_owner,npoint,nfragment,status,bad_point)
    if(status==0 .and. npoint/=size(rho_s(1)%f))status=12
    if(status==0)then
      do ispin=1,system%nspin
        do ifragment=1,nfragment
          do ipoint=1,npoint
            if(core_owner(ipoint,ifragment) .and. .not.ieee_is_finite(rho_core(ipoint,ifragment,ispin)))status=11
          enddo
        enddo
      enddo
    endif
    local_bad=merge(1,0,status/=0)
    call comm_summation(local_bad,global_bad,info%icomm_r)
    if(global_bad/=0)then
      if(status==0)status=13
      return
    endif
    allocate(rho_flat(npoint))
    do ispin=1,system%nspin
      rho_flat=0d0
      do ifragment=1,nfragment
        do ipoint=1,npoint
          if(core_owner(ipoint,ifragment))rho_flat(ipoint)=rho_core(ipoint,ifragment,ispin)
        enddo
      enddo
      rho_s(ispin)%f=reshape(rho_flat,shape(rho_s(ispin)%f))
    enddo
    deallocate(rho_flat)
    rho%f=0d0
    do ispin=1,system%nspin
      rho%f=rho%f+rho_s(ispin)%f
    enddo
    call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,vh)
    call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi,stencil,vxc,e_xc)
    local_bad=0
    if(.not.allocated(eexc_tmp))then
      status=15; local_bad=1
    else if(any(shape(eexc_tmp)/=mg%num))then
      status=15; local_bad=1
    endif
    call comm_summation(local_bad,global_bad,info%icomm_r)
    if(global_bad/=0)then
      if(status==0)status=16
      return
    endif
    e_xc_local=sum(eexc_tmp)*system%hvol
    n_vxc_local=0d0
    do ispin=1,system%nspin
      n_vxc_local=n_vxc_local+sum(rho_s(ispin)%f*vxc(ispin)%f)*system%hvol
    enddo
    call comm_summation(e_xc_local,e_xc_core,info%icomm_r)
    call comm_summation(n_vxc_local,n_vxc_core,info%icomm_r)
    if(abs(e_xc_core-e_xc)>1d-11*max(1d0,abs(e_xc)))status=14
  end subroutine update_wpw_lda_hartree
end module dg_wpw_lda_hartree
