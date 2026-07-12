!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Velocity-gauge current for core-owned nodal orbitals.
!=======================================================================
module rt_dg_nodal_current
#ifdef USE_MPI
  use mpi
#endif
  use structures, only: s_dft_system,s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state,nodal_face_slot
  use rt_dg_nodal_mpi, only: exchange_nodal_face_halos
  use rt_dg_nodal_nonlocal, only: global_to_local_core_index
  implicit none
  private
  public :: calculate_nodal_velocity_current_mpi
contains

  subroutine calculate_nodal_velocity_current_mpi(state,dg_frag,system,ppg,ik,gradient_offset, &
                                                    communicator,current)
    type(s_dg_nodal_state), intent(inout) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik,communicator
    real(8), intent(in) :: gradient_offset(:,:)
    real(8), intent(out) :: current(3)
    complex(8), allocatable :: work(:,:,:,:,:),projector_local(:,:,:),projector_global(:,:,:)
    complex(8) :: derivative(3),u_r(3),psi0
    real(8) :: current_local(3),current_global(3),local_cart(3),weight,kavec(3),bt(3,3)
    integer :: ix,iy,iz,istate,ispin,axis,iorder,ilma,ia,j,lx,ly,lz,ifrag
#ifdef USE_MPI
    integer :: ierr
#endif
    if (ik < 1 .or. ik > system%nk) stop 'nodal DG: current k index out of bounds'
    if (size(gradient_offset,2) /= 3) stop 'nodal DG: invalid current derivative stencil'
    call exchange_nodal_face_halos(state,communicator)
    call build_halo_work(state,work)
    current_local=0.0d0
    bt=transpose(system%rmatrix_B)
    kavec=system%vec_k(:,ik)+system%vec_Ac
    do ispin=1,state%nspin
      do istate=1,state%nstate
        weight=system%rocc(istate,ik,ispin)*system%wtk(ik)
        do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
          psi0=conjg(state%psi_core(ix,iy,iz,istate,ispin))
          derivative=(0.0d0,0.0d0)
          do axis=1,3
            do iorder=1,size(gradient_offset,1)
              derivative(axis)=derivative(axis)+gradient_offset(iorder,axis)*psi0* &
                positive_grid_value(work,state%halo_width,ix,iy,iz,axis,iorder,istate,ispin)
            end do
          end do
          local_cart=matmul(bt,2.0d0*aimag(derivative))
          current_local=current_local+weight*(kavec*abs(state%psi_core(ix,iy,iz,istate,ispin))**2+local_cart)
        end do; end do; end do
      end do
    end do
    deallocate(work)

    if (ppg%Nlma > 0) then
      allocate(projector_local(ppg%Nlma,state%nstate,state%nspin), &
               projector_global(ppg%Nlma,state%nstate,state%nspin))
      projector_local=(0.0d0,0.0d0)
      ifrag=state%fragment
      do ilma=1,ppg%Nlma
        ia=ppg%ia_tbl(ilma)
        do j=1,ppg%mps(ia)
          lx=global_to_local_core_index(ppg%jxyz(1,j,ia),dg_frag%ixyz_frag(1,ifrag), &
                                        state%core_size(1),dg_frag%lgnum_total(1))
          ly=global_to_local_core_index(ppg%jxyz(2,j,ia),dg_frag%ixyz_frag(2,ifrag), &
                                        state%core_size(2),dg_frag%lgnum_total(2))
          lz=global_to_local_core_index(ppg%jxyz(3,j,ia),dg_frag%ixyz_frag(3,ifrag), &
                                        state%core_size(3),dg_frag%lgnum_total(3))
          if (lx == 0 .or. ly == 0 .or. lz == 0) cycle
          do ispin=1,state%nspin; do istate=1,state%nstate
            projector_local(ilma,istate,ispin)=projector_local(ilma,istate,ispin)+ &
              conjg(ppg%zekr_uV(j,ilma,ik))*state%psi_core(lx,ly,lz,istate,ispin)
          end do; end do
        end do
      end do
#ifdef USE_MPI
      call MPI_Allreduce(projector_local,projector_global,size(projector_local),MPI_DOUBLE_COMPLEX, &
                         MPI_SUM,communicator,ierr)
      if (ierr /= MPI_SUCCESS) stop 'nodal DG: current projector reduction failed'
#else
      projector_global=projector_local
#endif
      do ilma=1,ppg%Nlma
        ia=ppg%ia_tbl(ilma)
        do ispin=1,state%nspin; do istate=1,state%nstate
          u_r=(0.0d0,0.0d0)
          do j=1,ppg%mps(ia)
            lx=global_to_local_core_index(ppg%jxyz(1,j,ia),dg_frag%ixyz_frag(1,ifrag), &
                                          state%core_size(1),dg_frag%lgnum_total(1))
            ly=global_to_local_core_index(ppg%jxyz(2,j,ia),dg_frag%ixyz_frag(2,ifrag), &
                                          state%core_size(2),dg_frag%lgnum_total(2))
            lz=global_to_local_core_index(ppg%jxyz(3,j,ia),dg_frag%ixyz_frag(3,ifrag), &
                                          state%core_size(3),dg_frag%lgnum_total(3))
            if (lx == 0 .or. ly == 0 .or. lz == 0) cycle
            u_r=u_r+conjg(ppg%zekr_uV(j,ilma,ik))*ppg%rxyz(:,j,ia)* &
              state%psi_core(lx,ly,lz,istate,ispin)
          end do
          weight=system%rocc(istate,ik,ispin)*system%wtk(ik)
          current_local=current_local+weight*2.0d0*aimag(conjg(u_r)* &
            (ppg%rinv_uvu(ilma)*projector_global(ilma,istate,ispin)))
        end do; end do
      end do
      deallocate(projector_local,projector_global)
    end if
#ifdef USE_MPI
    call MPI_Allreduce(current_local,current_global,3,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: current reduction failed'
#else
    current_global=current_local
#endif
    current=current_global/(dble(system%ngrid)*system%hvol)
  end subroutine calculate_nodal_velocity_current_mpi

  subroutine build_halo_work(state,work)
    type(s_dg_nodal_state), intent(in) :: state
    complex(8), allocatable, intent(out) :: work(:,:,:,:,:)
    integer :: axis,side,iface,layer,i1,i2,idx,nh,nx,ny,nz
    nh=state%halo_width; nx=state%core_size(1); ny=state%core_size(2); nz=state%core_size(3)
    allocate(work(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh,state%nstate,state%nspin))
    work=(0.0d0,0.0d0); work(1:nx,1:ny,1:nz,:,:)=state%psi_core
    do axis=1,3; do side=-1,1,2
      iface=nodal_face_slot(axis,side)
      do layer=1,nh
        idx=merge(1-layer,state%core_size(axis)+layer,side < 0)
        do i2=1,state%faces(iface)%tangential_size(2); do i1=1,state%faces(iface)%tangential_size(1)
          select case(axis)
          case(1); work(idx,i1,i2,:,:)=state%faces(iface)%recv_value(layer,i1,i2,:,:)
          case(2); work(i1,idx,i2,:,:)=state%faces(iface)%recv_value(layer,i1,i2,:,:)
          case(3); work(i1,i2,idx,:,:)=state%faces(iface)%recv_value(layer,i1,i2,:,:)
          end select
        end do; end do
      end do
    end do; end do
  end subroutine build_halo_work

  pure complex(8) function positive_grid_value(work,nh,ix,iy,iz,axis,offset,istate,ispin) result(value)
    integer, intent(in) :: nh,ix,iy,iz,axis,offset,istate,ispin
    complex(8), intent(in) :: work(1-nh:,1-nh:,1-nh:,:,:)
    integer :: point(3)
    point=[ix,iy,iz]; point(axis)=point(axis)+offset
    value=work(point(1),point(2),point(3),istate,ispin)
  end function positive_grid_value

end module rt_dg_nodal_current
