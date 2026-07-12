!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Nonlocal pseudopotential action on core-owned nodal fragment states.
!=======================================================================
module rt_dg_nodal_nonlocal
#ifdef USE_MPI
  use mpi
#endif
  use structures, only: s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  implicit none
  private
  public :: apply_nodal_nonlocal_potential_mpi, global_to_local_core_index

contains

  subroutine apply_nodal_nonlocal_potential_mpi(state, dg_frag, ppg, ik, communicator, hpsi)
    type(s_dg_nodal_state), intent(in) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik, communicator
    complex(8), intent(inout) :: hpsi(:,:,:,:,:)
    complex(8), allocatable :: projector_local(:,:,:), projector_global(:,:,:)
    integer :: ilma, ia, j, ix, iy, iz, istate, ispin, ifrag
#ifdef USE_MPI
    integer :: ierr
#endif

    if (ppg%Nlma <= 0) return
    if (.not. allocated(ppg%zekr_uV) .or. .not. allocated(ppg%rinv_uvu)) &
      stop 'nodal DG: nonlocal projector data are missing'
    if (ik < lbound(ppg%zekr_uV,3) .or. ik > ubound(ppg%zekr_uV,3)) &
      stop 'nodal DG: nonlocal projector k index out of bounds'
    ifrag = state%fragment
    allocate(projector_local(ppg%Nlma,state%nstate,state%nspin))
    allocate(projector_global(ppg%Nlma,state%nstate,state%nspin))
    projector_local = (0.0d0,0.0d0)

    do ilma = 1, ppg%Nlma
      ia = ppg%ia_tbl(ilma)
      do j = 1, ppg%mps(ia)
        ix = global_to_local_core_index(ppg%jxyz(1,j,ia),dg_frag%ixyz_frag(1,ifrag), &
                                        state%core_size(1),dg_frag%lgnum_total(1))
        iy = global_to_local_core_index(ppg%jxyz(2,j,ia),dg_frag%ixyz_frag(2,ifrag), &
                                        state%core_size(2),dg_frag%lgnum_total(2))
        iz = global_to_local_core_index(ppg%jxyz(3,j,ia),dg_frag%ixyz_frag(3,ifrag), &
                                        state%core_size(3),dg_frag%lgnum_total(3))
        if (ix == 0 .or. iy == 0 .or. iz == 0) cycle
        do ispin = 1, state%nspin
          do istate = 1, state%nstate
            projector_local(ilma,istate,ispin) = projector_local(ilma,istate,ispin) + &
              conjg(ppg%zekr_uV(j,ilma,ik)) * state%psi_core(ix,iy,iz,istate,ispin)
          end do
        end do
      end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(projector_local,projector_global,size(projector_local),MPI_DOUBLE_COMPLEX, &
                       MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: nonlocal projector reduction failed'
#else
    projector_global = projector_local
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
    do ilma = 1, ppg%Nlma
      projector_global(ilma,:,:) = ppg%rinv_uvu(ilma) * projector_global(ilma,:,:)
      ia = ppg%ia_tbl(ilma)
      do j = 1, ppg%mps(ia)
        ix = global_to_local_core_index(ppg%jxyz(1,j,ia),dg_frag%ixyz_frag(1,ifrag), &
                                        state%core_size(1),dg_frag%lgnum_total(1))
        iy = global_to_local_core_index(ppg%jxyz(2,j,ia),dg_frag%ixyz_frag(2,ifrag), &
                                        state%core_size(2),dg_frag%lgnum_total(2))
        iz = global_to_local_core_index(ppg%jxyz(3,j,ia),dg_frag%ixyz_frag(3,ifrag), &
                                        state%core_size(3),dg_frag%lgnum_total(3))
        if (ix == 0 .or. iy == 0 .or. iz == 0) cycle
        do ispin = 1, state%nspin
          do istate = 1, state%nstate
            hpsi(ix,iy,iz,istate,ispin) = hpsi(ix,iy,iz,istate,ispin) + &
              ppg%zekr_uV(j,ilma,ik) * projector_global(ilma,istate,ispin)
          end do
        end do
      end do
    end do
    deallocate(projector_local,projector_global)
  end subroutine apply_nodal_nonlocal_potential_mpi

  pure integer function global_to_local_core_index(global_index,origin,ncore,ngrid) result(local_index)
    integer, intent(in) :: global_index, origin, ncore, ngrid
    integer :: offset, wrapped_origin, wrapped_global

    local_index = 0
    if (ngrid < 1 .or. ncore < 1) return
    wrapped_origin = modulo(origin-1,ngrid)+1
    wrapped_global = modulo(global_index-1,ngrid)+1
    offset = modulo(wrapped_global-wrapped_origin,ngrid)
    if (offset < ncore) local_index = offset+1
  end function global_to_local_core_index

end module rt_dg_nodal_nonlocal
