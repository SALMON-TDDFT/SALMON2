!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Reconstruct a replicated all-system density from core-owned nodal states.
!=======================================================================
module rt_dg_nodal_density
#ifdef USE_MPI
  use mpi
#endif
  use structures, only: s_dft_system,s_scalar
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use sym_rho_sub, only: sym_rho
  implicit none
  private
  public :: reconstruct_nodal_density_mpi
contains

  subroutine reconstruct_nodal_density_mpi(state,dg_frag,system,ik,communicator,rho_s,rho,electron_number)
    type(s_dg_nodal_state), intent(in) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    integer, intent(in) :: ik,communicator
    type(s_scalar), intent(inout) :: rho_s(system%nspin),rho
    real(8), intent(out) :: electron_number
    real(8), allocatable :: density_local(:,:,:)
    integer :: ix,iy,iz,ixg,iyg,izg,istate,ispin,ifrag
#ifdef USE_MPI
    integer :: ierr
#endif
    if (state%nspin /= system%nspin) stop 'nodal DG: density spin count mismatch'
    if (state%nstate > system%no) stop 'nodal DG: density state count exceeds SALMON states'
    ifrag=state%fragment
    do ispin=1,state%nspin
      allocate(density_local(lbound(rho_s(ispin)%f,1):ubound(rho_s(ispin)%f,1), &
                             lbound(rho_s(ispin)%f,2):ubound(rho_s(ispin)%f,2), &
                             lbound(rho_s(ispin)%f,3):ubound(rho_s(ispin)%f,3)))
      density_local=0.0d0
      do iz=1,state%core_size(3)
        izg=modulo(dg_frag%ixyz_frag(3,ifrag)+iz-2,dg_frag%lgnum_total(3))+1
        do iy=1,state%core_size(2)
          iyg=modulo(dg_frag%ixyz_frag(2,ifrag)+iy-2,dg_frag%lgnum_total(2))+1
          do ix=1,state%core_size(1)
            ixg=modulo(dg_frag%ixyz_frag(1,ifrag)+ix-2,dg_frag%lgnum_total(1))+1
            do istate=1,state%nstate
              density_local(ixg,iyg,izg)=density_local(ixg,iyg,izg)+ &
                system%rocc(istate,ik,ispin)*system%wtk(ik)* &
                abs(state%psi_core(ix,iy,iz,istate,ispin))**2/system%hvol
            end do
          end do
        end do
      end do
#ifdef USE_MPI
      call MPI_Allreduce(density_local,rho_s(ispin)%f,size(density_local),MPI_DOUBLE_PRECISION, &
                         MPI_SUM,communicator,ierr)
      if (ierr /= MPI_SUCCESS) stop 'nodal DG: density reduction failed'
#else
      rho_s(ispin)%f=density_local
#endif
      call sym_rho(rho_s(ispin)%f)
      deallocate(density_local)
    end do
    rho%f=0.0d0
    do ispin=1,state%nspin
      rho%f=rho%f+rho_s(ispin)%f
    end do
    electron_number=system%hvol*sum(rho%f)
  end subroutine reconstruct_nodal_density_mpi

end module rt_dg_nodal_density
