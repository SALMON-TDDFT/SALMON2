!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
module rt_dg_nodal_rayleigh_ritz
#ifdef USE_MPI
  use mpi
#endif
  use eigen_subdiag_sub, only: eigen_zheev
  use rt_dg_nodal_types, only: s_dg_nodal_state
  implicit none
  private
  public :: rayleigh_ritz_nodal_subspace_mpi
  logical, save :: trace_hermiticity = .true.
contains
  subroutine rayleigh_ritz_nodal_subspace_mpi(state,hpsi,eigenvalues,communicator)
    type(s_dg_nodal_state), intent(inout) :: state
    complex(8), intent(inout) :: hpsi(:,:,:,:,:)
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin)
    integer, intent(in) :: communicator
    complex(8), allocatable :: hsub_local(:,:),hsub_global(:,:),evec(:,:)
    complex(8), allocatable :: psi_matrix(:,:),hpsi_matrix(:,:)
    integer :: ispin,i,j,ix,iy,iz,igrid,ngrid,myrank
    real(8) :: antiherm_rel
#ifdef USE_MPI
    integer :: ierr
#endif
    myrank=0
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,myrank,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: Rayleigh-Ritz rank query failed'
#endif
    ngrid=product(state%core_size)
    allocate(hsub_local(state%nstate,state%nstate),hsub_global(state%nstate,state%nstate))
    allocate(evec(state%nstate,state%nstate),psi_matrix(ngrid,state%nstate),hpsi_matrix(ngrid,state%nstate))
    do ispin=1,state%nspin
      hsub_local=(0.0d0,0.0d0)
      do j=1,state%nstate
        do i=1,state%nstate
          hsub_local(i,j)=sum(conjg(state%psi_core(:,:,:,i,ispin))*hpsi(:,:,:,j,ispin))
        end do
      end do
#ifdef USE_MPI
      call MPI_Allreduce(hsub_local,hsub_global,size(hsub_local),MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
      if (ierr /= MPI_SUCCESS) stop 'nodal DG: Rayleigh-Ritz reduction failed'
#else
      hsub_global=hsub_local
      if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
      antiherm_rel=maxval(abs(hsub_global-conjg(transpose(hsub_global)))) / &
        max(maxval(abs(hsub_global)),tiny(1.0d0))
      if (trace_hermiticity .and. myrank == 0) then
        write(*,'(1x,a,i0,a,1pe13.5)') '[DG-NODAL-HERMITICITY] ispin=',ispin, &
          ' antiherm_rel=',antiherm_rel
        flush(6)
      end if
      hsub_global=0.5d0*(hsub_global+conjg(transpose(hsub_global)))
      call eigen_zheev(hsub_global,eigenvalues(:,ispin),evec)
      do j=1,state%nstate
        do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
          igrid=ix+state%core_size(1)*((iy-1)+state%core_size(2)*(iz-1))
          psi_matrix(igrid,j)=state%psi_core(ix,iy,iz,j,ispin)
          hpsi_matrix(igrid,j)=hpsi(ix,iy,iz,j,ispin)
        end do; end do; end do
      end do
      psi_matrix=matmul(psi_matrix,evec)
      hpsi_matrix=matmul(hpsi_matrix,evec)
      do j=1,state%nstate
        do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
          igrid=ix+state%core_size(1)*((iy-1)+state%core_size(2)*(iz-1))
          state%psi_core(ix,iy,iz,j,ispin)=psi_matrix(igrid,j)
          hpsi(ix,iy,iz,j,ispin)=hpsi_matrix(igrid,j)
        end do; end do; end do
      end do
    end do
    trace_hermiticity=.false.
    deallocate(hsub_local,hsub_global,evec,psi_matrix,hpsi_matrix)
  end subroutine rayleigh_ritz_nodal_subspace_mpi
end module rt_dg_nodal_rayleigh_ritz
