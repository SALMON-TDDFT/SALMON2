!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
module rt_dg_nodal_diagnostics
#ifdef USE_MPI
  use mpi
#endif
  use rt_dg_nodal_types, only: s_dg_nodal_state
  implicit none
  private
  public :: calculate_nodal_norm_diagnostics_mpi
contains

  subroutine calculate_nodal_norm_diagnostics_mpi(state,communicator,total_norm,max_orbital_drift)
    type(s_dg_nodal_state), intent(in) :: state
    integer, intent(in) :: communicator
    real(8), intent(out) :: total_norm,max_orbital_drift
    real(8), allocatable :: norm_local(:,:),norm_global(:,:)
    integer :: istate,ispin
#ifdef USE_MPI
    integer :: ierr
#endif
    allocate(norm_local(state%nstate,state%nspin),norm_global(state%nstate,state%nspin))
    do ispin=1,state%nspin
      do istate=1,state%nstate
        norm_local(istate,ispin)=sum(abs(state%psi_core(:,:,:,istate,ispin))**2)
      end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(norm_local,norm_global,size(norm_local),MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: norm diagnostic reduction failed'
#else
    norm_global=norm_local
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
    total_norm=sum(norm_global)
    max_orbital_drift=maxval(abs(norm_global-1.0d0))
    deallocate(norm_local,norm_global)
  end subroutine calculate_nodal_norm_diagnostics_mpi

end module rt_dg_nodal_diagnostics
