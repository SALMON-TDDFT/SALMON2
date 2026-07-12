!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Matrix-free local Hamiltonian action on core-owned nodal DG states.
!=======================================================================
module rt_dg_nodal_hamiltonian
  use rt_dg_nodal_types, only: s_dg_nodal_state, nodal_face_slot
  implicit none
  private
  public :: apply_nodal_local_hamiltonian

contains

  subroutine apply_nodal_local_hamiltonian(state, vlocal, kinetic_center, kinetic_offset, &
                                            gradient_offset, avec, hpsi)
    type(s_dg_nodal_state), intent(in) :: state
    real(8), intent(in) :: vlocal(:,:,:,:)
    real(8), intent(in) :: kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:)
    real(8), intent(in) :: gradient_offset(:,:)
    real(8), intent(in) :: avec(3)
    complex(8), intent(out) :: hpsi(:,:,:,:,:)
    complex(8), allocatable :: work(:,:,:,:,:)
    complex(8) :: value
    real(8) :: a2
    integer :: nx, ny, nz, nh, norder
    integer :: ix, iy, iz, istate, ispin, axis, iorder

    if (.not. state%initialized) stop 'nodal DG: Hamiltonian applied before state initialization'
    nx = state%core_size(1)
    ny = state%core_size(2)
    nz = state%core_size(3)
    nh = state%halo_width
    norder = size(kinetic_offset, 1)
    if (size(kinetic_offset,2) /= 3 .or. any(shape(gradient_offset) /= shape(kinetic_offset))) &
      stop 'nodal DG: inconsistent stencil coefficients'
    if (nh < norder) stop 'nodal DG: halo is narrower than the finite-difference stencil'
    if (any(shape(vlocal) /= [nx,ny,nz,state%nspin])) stop 'nodal DG: local-potential shape mismatch'
    if (any(shape(hpsi) /= [nx,ny,nz,state%nstate,state%nspin])) stop 'nodal DG: output shape mismatch'

    allocate(work(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh,state%nstate,state%nspin))
    work = (0.0d0, 0.0d0)
    work(1:nx,1:ny,1:nz,:,:) = state%psi_core
    call load_face_halos(state, work)

    a2 = 0.5d0 * dot_product(avec, avec)
    do ispin = 1, state%nspin
      do istate = 1, state%nstate
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              value = (kinetic_center + vlocal(ix,iy,iz,ispin) + a2) * work(ix,iy,iz,istate,ispin)
              do axis = 1, 3
                do iorder = 1, norder
                  if (kinetic_offset(iorder,axis) == 0.0d0 .and. &
                      gradient_offset(iorder,axis) == 0.0d0) cycle
                  value = value + kinetic_offset(iorder,axis) * &
                    (grid_value(work,nh,ix,iy,iz,axis,+iorder,istate,ispin) + &
                     grid_value(work,nh,ix,iy,iz,axis,-iorder,istate,ispin))
                  value = value - (0.0d0,1.0d0) * avec(axis) * gradient_offset(iorder,axis) * &
                    (grid_value(work,nh,ix,iy,iz,axis,+iorder,istate,ispin) - &
                     grid_value(work,nh,ix,iy,iz,axis,-iorder,istate,ispin))
                end do
              end do
              hpsi(ix,iy,iz,istate,ispin) = value
            end do
          end do
        end do
      end do
    end do
    deallocate(work)
  end subroutine apply_nodal_local_hamiltonian

  subroutine load_face_halos(state, work)
    type(s_dg_nodal_state), intent(in) :: state
    complex(8), intent(inout) :: work(1-state%halo_width:,1-state%halo_width:,1-state%halo_width:,:,:)
    integer :: axis, side, iface, layer, i1, i2, idx

    do axis = 1, 3
      do side = -1, 1, 2
        iface = nodal_face_slot(axis, side)
        if (.not. allocated(state%faces(iface)%recv_value)) stop 'nodal DG: face halo is not allocated'
        do layer = 1, state%halo_width
          idx = merge(1-layer, state%core_size(axis)+layer, side < 0)
          do i2 = 1, state%faces(iface)%tangential_size(2)
            do i1 = 1, state%faces(iface)%tangential_size(1)
              select case(axis)
              case(1)
                work(idx,i1,i2,:,:) = state%faces(iface)%recv_value(layer,i1,i2,:,:)
              case(2)
                work(i1,idx,i2,:,:) = state%faces(iface)%recv_value(layer,i1,i2,:,:)
              case(3)
                work(i1,i2,idx,:,:) = state%faces(iface)%recv_value(layer,i1,i2,:,:)
              end select
            end do
          end do
        end do
      end do
    end do
  end subroutine load_face_halos

  pure complex(8) function grid_value(work, nh, ix, iy, iz, axis, offset, istate, ispin) result(value)
    integer, intent(in) :: nh
    complex(8), intent(in) :: work(1-nh:,1-nh:,1-nh:,:,:)
    integer, intent(in) :: ix, iy, iz, axis, offset, istate, ispin
    integer :: point(3)

    point = [ix,iy,iz]
    point(axis) = point(axis) + offset
    value = work(point(1),point(2),point(3),istate,ispin)
  end function grid_value

end module rt_dg_nodal_hamiltonian
