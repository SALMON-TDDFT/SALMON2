!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Adapter from replicated SALMON real-space data to one core-owned nodal
! fragment per MPI rank. This is the correctness-first stage-1 layout.
!=======================================================================
module rt_dg_nodal_salmon_adapter
  use structures, only: s_dft_system, s_parallel_info, s_orbital, s_scalar, s_stencil
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state, allocate_nodal_state, allocate_nodal_face_buffers
  use rt_dg_nodal_halo, only: initialize_nodal_face_topology
  implicit none
  private
  public :: initialize_nodal_from_full_orbital, initialize_nodal_from_dg_coefficients
  public :: build_nodal_local_potential
  public :: get_nodal_stencil_coefficients

contains

  subroutine initialize_nodal_from_dg_coefficients(dg_frag,system,info,halo_width,state)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    integer, intent(in) :: halo_width
    type(s_dg_nodal_state), intent(inout) :: state
    integer :: ifrag,i_local,nx,ny,nz,nstate,nocc_max,ispin,istate,ibasis,gid,lrow
    integer :: ix,iy,iz,gx,gy,gz,iface,ilocal,igrid,nbasis,fragment_coords(3),fragment_rank(dg_frag%n_frag)
    complex(8), allocatable :: phi_matrix(:,:),coef_matrix(:,:),psi_matrix(:,:)

    if (dg_frag%nproc_frag /= 1 .or. info%nporbital /= 1) &
      stop 'nodal DG coefficient seed requires one MPI rank per fragment'
    if (any(info%nprgrid(1:3) /= dg_frag%num_fragment(1:3))) &
      stop 'nodal DG coefficient seed requires parent domains to match fragments'
    if (.not. allocated(dg_frag%phi_frag) .or. .not. allocated(dg_frag%coef) .or. &
        .not. allocated(dg_frag%coef_global_to_local)) &
      stop 'nodal DG coefficient seed data are incomplete'
    ifrag=dg_frag%ifrag_group
    i_local=ifrag-dg_frag%ifrag_start+1
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag,5)) &
      stop 'nodal DG coefficient seed local fragment index is invalid'
    nocc_max=0
    do ispin=1,system%nspin
      nocc_max=max(nocc_max,count(system%rocc(:,1,ispin) > 1.0d-12))
    end do
    if (nocc_max < 1) stop 'nodal DG coefficient seed has no occupied states'
    nstate=min(nocc_max,size(dg_frag%coef,2))
    nx=dg_frag%nxyz_domain(1,ifrag); ny=dg_frag%nxyz_domain(2,ifrag); nz=dg_frag%nxyz_domain(3,ifrag)
    call allocate_nodal_state(state,ifrag,[nx,ny,nz],halo_width,nstate,system%nspin)
    do ilocal=1,dg_frag%n_frag
      fragment_coords(1)=(ilocal-1)/(dg_frag%num_fragment(2)*dg_frag%num_fragment(3))
      fragment_coords(2)=modulo((ilocal-1)/dg_frag%num_fragment(3),dg_frag%num_fragment(2))
      fragment_coords(3)=modulo(ilocal-1,dg_frag%num_fragment(3))
      fragment_rank(ilocal)=info%imap(fragment_coords(1),fragment_coords(2),fragment_coords(3), &
                                      lbound(info%imap,4),lbound(info%imap,5))
    end do
    fragment_coords(1)=(ifrag-1)/(dg_frag%num_fragment(2)*dg_frag%num_fragment(3))
    fragment_coords(2)=modulo((ifrag-1)/dg_frag%num_fragment(3),dg_frag%num_fragment(2))
    fragment_coords(3)=modulo(ifrag-1,dg_frag%num_fragment(3))
    call initialize_nodal_face_topology(state%faces,fragment_coords,dg_frag%num_fragment,fragment_rank,halo_width)
    do iface=1,size(state%faces)
      call allocate_nodal_face_buffers(state%faces(iface),state%core_size,halo_width,nstate,system%nspin)
    end do
    state%psi_core=(0.0d0,0.0d0)
    do ispin=1,system%nspin
      nbasis=dg_frag%n_basis(ifrag,ispin)
      allocate(phi_matrix(nx*ny*nz,nbasis),coef_matrix(nbasis,nstate),psi_matrix(nx*ny*nz,nstate))
      phi_matrix=(0.0d0,0.0d0); coef_matrix=(0.0d0,0.0d0)
      do ibasis=1,nbasis
        gid=dg_frag%index_basis(ibasis,ifrag,ispin)
        if (gid < 1 .or. gid > size(dg_frag%coef_global_to_local,1)) cycle
        lrow=dg_frag%coef_global_to_local(gid,ispin)
        if (lrow < 1 .or. lrow > size(dg_frag%coef,1)) cycle
        coef_matrix(ibasis,:)=dg_frag%coef(lrow,1:nstate,ispin)
        do iz=1,nz
          gz=dg_frag%ixyz_frag(3,ifrag)+iz-1
          do iy=1,ny
            gy=dg_frag%ixyz_frag(2,ifrag)+iy-1
            do ix=1,nx
              gx=dg_frag%ixyz_frag(1,ifrag)+ix-1
              if (gx < lbound(dg_frag%phi_frag,1) .or. gx > ubound(dg_frag%phi_frag,1) .or. &
                  gy < lbound(dg_frag%phi_frag,2) .or. gy > ubound(dg_frag%phi_frag,2) .or. &
                  gz < lbound(dg_frag%phi_frag,3) .or. gz > ubound(dg_frag%phi_frag,3)) &
                stop 'nodal DG coefficient seed point is outside phi_frag core box'
              igrid=ix+nx*((iy-1)+ny*(iz-1))
              phi_matrix(igrid,ibasis)=cmplx(dg_frag%phi_frag(gx,gy,gz,ibasis,i_local),0.0d0,8)
            end do
          end do
        end do
      end do
      psi_matrix=matmul(phi_matrix,coef_matrix)
      do istate=1,nstate
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
              igrid=ix+nx*((iy-1)+ny*(iz-1))
              state%psi_core(ix,iy,iz,istate,ispin)=psi_matrix(igrid,istate)
            end do
          end do
        end do
      end do
      deallocate(phi_matrix,coef_matrix,psi_matrix)
    end do
    state%enabled=.true.
  end subroutine initialize_nodal_from_dg_coefficients

  subroutine initialize_nodal_from_full_orbital(dg_frag, system, info, spsi, halo_width, state)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_orbital), intent(in) :: spsi
    integer, intent(in) :: halo_width
    type(s_dg_nodal_state), intent(inout) :: state
    integer :: ifrag, fragment_coords(3), fragment_rank(dg_frag%n_frag)
    integer :: nx, ny, nz, nstate, nstate_available, nocc_max, ispin, ilocal, io, ix, iy, iz, ixg, iyg, izg
    integer :: ik, im, iface
    logical :: use_complex, use_real

    if (dg_frag%nproc_frag /= 1 .or. info%nporbital /= 1) &
      stop 'nodal DG stage-1 requires one MPI rank per fragment and no orbital decomposition'
    if (any(info%nprgrid(1:3) /= dg_frag%num_fragment(1:3))) &
      stop 'nodal DG stage-1 requires parent real-space domains to match fragment cores'
    if (.not. allocated(info%imap)) stop 'nodal DG: parent rank-address map is unavailable'
    if (dg_frag%n_frag /= product(dg_frag%num_fragment)) &
      stop 'nodal DG: inconsistent fragment topology'
    ifrag = dg_frag%ifrag_group
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) stop 'nodal DG: invalid local fragment id'
    use_complex = allocated(spsi%zwf)
    use_real = allocated(spsi%rwf)
    if (.not. use_complex .and. .not. use_real) stop 'nodal DG: missing full real-space seed orbital'

    nx = dg_frag%nxyz_domain(1,ifrag)
    ny = dg_frag%nxyz_domain(2,ifrag)
    nz = dg_frag%nxyz_domain(3,ifrag)
    if (use_complex) then
      nstate_available = size(spsi%zwf,5)
      ik = lbound(spsi%zwf,6)
      im = lbound(spsi%zwf,7)
    else
      nstate_available = size(spsi%rwf,5)
      ik = lbound(spsi%rwf,6)
      im = lbound(spsi%rwf,7)
    end if
    nocc_max = 0
    if (allocated(system%rocc)) then
      do ispin=1,system%nspin
        nocc_max=max(nocc_max,count(system%rocc(:,1,ispin) > 1.0d-12))
      end do
    end if
    if (nocc_max < 1) stop 'nodal DG: no occupied seed orbitals were found'
    nstate=min(nstate_available,nocc_max)
    call allocate_nodal_state(state,ifrag,[nx,ny,nz],halo_width,nstate,system%nspin)

    fragment_coords(1) = (ifrag-1)/(dg_frag%num_fragment(2)*dg_frag%num_fragment(3))
    fragment_coords(2) = modulo((ifrag-1)/dg_frag%num_fragment(3),dg_frag%num_fragment(2))
    fragment_coords(3) = modulo(ifrag-1,dg_frag%num_fragment(3))
    do ilocal = 1, dg_frag%n_frag
      fragment_coords(1) = (ilocal-1)/(dg_frag%num_fragment(2)*dg_frag%num_fragment(3))
      fragment_coords(2) = modulo((ilocal-1)/dg_frag%num_fragment(3),dg_frag%num_fragment(2))
      fragment_coords(3) = modulo(ilocal-1,dg_frag%num_fragment(3))
      fragment_rank(ilocal) = info%imap(fragment_coords(1),fragment_coords(2),fragment_coords(3), &
                                        lbound(info%imap,4),lbound(info%imap,5))
    end do
    fragment_coords(1) = (ifrag-1)/(dg_frag%num_fragment(2)*dg_frag%num_fragment(3))
    fragment_coords(2) = modulo((ifrag-1)/dg_frag%num_fragment(3),dg_frag%num_fragment(2))
    fragment_coords(3) = modulo(ifrag-1,dg_frag%num_fragment(3))
    call initialize_nodal_face_topology(state%faces,fragment_coords,dg_frag%num_fragment, &
                                        fragment_rank,halo_width)
    do iface = 1, size(state%faces)
      call allocate_nodal_face_buffers(state%faces(iface),state%core_size,halo_width,nstate,system%nspin)
    end do

    do ispin = 1, system%nspin
      do ilocal = 1, nstate
        if (use_complex) then
          io = lbound(spsi%zwf,5) + ilocal - 1
        else
          io = lbound(spsi%rwf,5) + ilocal - 1
        end if
        do iz = 1, nz
          izg = modulo(dg_frag%ixyz_frag(3,ifrag)+iz-2,dg_frag%lgnum_total(3))+1
          do iy = 1, ny
            iyg = modulo(dg_frag%ixyz_frag(2,ifrag)+iy-2,dg_frag%lgnum_total(2))+1
            do ix = 1, nx
              ixg = modulo(dg_frag%ixyz_frag(1,ifrag)+ix-2,dg_frag%lgnum_total(1))+1
              if (use_complex) then
                if (ixg < lbound(spsi%zwf,1) .or. ixg > ubound(spsi%zwf,1) .or. &
                    iyg < lbound(spsi%zwf,2) .or. iyg > ubound(spsi%zwf,2) .or. &
                    izg < lbound(spsi%zwf,3) .or. izg > ubound(spsi%zwf,3)) &
                  stop 'nodal DG: seed orbital grid is not replicated'
                state%psi_core(ix,iy,iz,ilocal,ispin) = spsi%zwf(ixg,iyg,izg,ispin,io,ik,im)
              else
                if (ixg < lbound(spsi%rwf,1) .or. ixg > ubound(spsi%rwf,1) .or. &
                    iyg < lbound(spsi%rwf,2) .or. iyg > ubound(spsi%rwf,2) .or. &
                    izg < lbound(spsi%rwf,3) .or. izg > ubound(spsi%rwf,3)) &
                  stop 'nodal DG: seed orbital grid is not replicated'
                state%psi_core(ix,iy,iz,ilocal,ispin) = &
                  cmplx(spsi%rwf(ixg,iyg,izg,ispin,io,ik,im),0.0d0,8)
              end if
            end do
          end do
        end do
      end do
    end do
    state%enabled = .true.
  end subroutine initialize_nodal_from_full_orbital

  subroutine build_nodal_local_potential(dg_frag, state, Vh, Vxc, Vpsl, vlocal)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dg_nodal_state), intent(in) :: state
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), allocatable, intent(out) :: vlocal(:,:,:,:)
    integer :: ifrag, ix, iy, iz, ixg, iyg, izg, ispin

    ifrag = state%fragment
    allocate(vlocal(state%core_size(1),state%core_size(2),state%core_size(3),state%nspin))
    do ispin = 1, state%nspin
      do iz = 1, state%core_size(3)
        izg = modulo(dg_frag%ixyz_frag(3,ifrag)+iz-2,dg_frag%lgnum_total(3))+1
        do iy = 1, state%core_size(2)
          iyg = modulo(dg_frag%ixyz_frag(2,ifrag)+iy-2,dg_frag%lgnum_total(2))+1
          do ix = 1, state%core_size(1)
            ixg = modulo(dg_frag%ixyz_frag(1,ifrag)+ix-2,dg_frag%lgnum_total(1))+1
            if (ixg < lbound(Vh%f,1) .or. ixg > ubound(Vh%f,1) .or. &
                iyg < lbound(Vh%f,2) .or. iyg > ubound(Vh%f,2) .or. &
                izg < lbound(Vh%f,3) .or. izg > ubound(Vh%f,3)) &
              stop 'nodal DG: local potential grid is not replicated'
            vlocal(ix,iy,iz,ispin) = Vh%f(ixg,iyg,izg) + Vxc(ispin)%f(ixg,iyg,izg) + Vpsl%f(ixg,iyg,izg)
          end do
        end do
      end do
    end do
  end subroutine build_nodal_local_potential

  subroutine get_nodal_stencil_coefficients(stencil, kinetic_center, kinetic_offset, gradient_offset)
    type(s_stencil), intent(in) :: stencil
    real(8), intent(out) :: kinetic_center, kinetic_offset(4,3), gradient_offset(4,3)

    if (.not. stencil%if_orthogonal) stop 'nodal DG stage-1 requires an orthogonal real-space grid'
    kinetic_center = stencil%coef_lap0
    kinetic_offset = -0.5d0 * stencil%coef_lap
    gradient_offset = stencil%coef_nab
  end subroutine get_nodal_stencil_coefficients

end module rt_dg_nodal_salmon_adapter
