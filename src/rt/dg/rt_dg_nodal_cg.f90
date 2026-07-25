!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Conjugate-gradient eigensolver for a matrix-free nodal Hamiltonian.
!=======================================================================
module rt_dg_nodal_cg
#ifdef USE_MPI
  use mpi
#endif
  use rt_dg_nodal_types, only: s_dg_nodal_state, accept_nodal_dg_ground_state
  use dg_nodal_interfaces, only: nodal_hamiltonian_action => nodal_complete_h_action
  use dg_nodal_interfaces, only: nodal_subspace_rotation
  use rt_dg_nodal_ground_state_solver, only: orthonormalize_nodal_states_mpi
  implicit none
  private
  public :: solve_nodal_ground_state_cg_mpi

contains

  subroutine solve_nodal_ground_state_cg_mpi(state,apply_hamiltonian,max_iteration,tolerance, &
                                               communicator,eigenvalues,max_residual,niteration, &
                                               rotate_subspace,require_convergence)
    type(s_dg_nodal_state), intent(inout) :: state
    procedure(nodal_hamiltonian_action) :: apply_hamiltonian
    integer, intent(in) :: max_iteration, communicator
    real(8), intent(in) :: tolerance
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin),max_residual
    integer, intent(out) :: niteration
    procedure(nodal_subspace_rotation) :: rotate_subspace
    logical, intent(in), optional :: require_convergence
    complex(8), allocatable :: hpsi(:,:,:,:,:),residual(:,:,:,:,:),search_direction(:,:,:,:,:)
    complex(8), allocatable :: hsearch(:,:,:,:,:),psi_saved(:,:,:,:,:)
    real(8), allocatable :: residual_norm(:,:),previous_residual_norm(:,:)
    complex(8) :: h2(2,2),s2(2,2),work(9),phase
    real(8) :: eval2(2),rwork(9),beta,pnorm_local,pnorm_global
    complex(8) :: overlap_local,overlap_global
    integer :: iteration,istate,ispin,ierr,myrank
    logical :: must_converge

    if (max_iteration < 1 .or. tolerance <= 0.0d0) stop 'nodal DG: invalid CG eigensolver control'
    must_converge=.true.
    if(present(require_convergence)) must_converge=require_convergence
    state%ground_state_ready=.false.
    allocate(hpsi,mold=state%psi_core)
    allocate(residual,mold=state%psi_core)
    allocate(search_direction,mold=state%psi_core)
    allocate(hsearch,mold=state%psi_core)
    allocate(psi_saved,mold=state%psi_core)
    hpsi=(0.0d0,0.0d0); residual=(0.0d0,0.0d0); search_direction=(0.0d0,0.0d0)
    hsearch=(0.0d0,0.0d0); psi_saved=(0.0d0,0.0d0)
    allocate(residual_norm(state%nstate,state%nspin),previous_residual_norm(state%nstate,state%nspin))
    previous_residual_norm=0.0d0
    myrank=0
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,myrank,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: CG rank query failed'
#endif
    call orthonormalize_nodal_states_mpi(state,communicator)
    call apply_hamiltonian(state,hpsi)
    call rotate_subspace(state,hpsi,eigenvalues,communicator)

    do iteration=1,max_iteration
      call build_residuals(state,hpsi,eigenvalues,residual,residual_norm,max_residual,communicator)
      if (myrank == 0 .and. (iteration == 1 .or. mod(iteration,10) == 0 .or. max_residual <= tolerance)) then
        write(*,'(1x,a,i0,3(a,1pe13.5))') '[DG-NODAL-CG-ITER] iter=',iteration, &
          ' residual=',max_residual,' eigen_min=',minval(eigenvalues),' eigen_max=',maxval(eigenvalues)
        flush(6)
      end if
      if (max_residual <= tolerance) exit

      do ispin=1,state%nspin
        do istate=1,state%nstate
          if (previous_residual_norm(istate,ispin) > tiny(1.0d0)) then
            beta=residual_norm(istate,ispin)/previous_residual_norm(istate,ispin)
          else
            beta=0.0d0
          end if
          search_direction(:,:,:,istate,ispin)=-residual(:,:,:,istate,ispin)+ &
            beta*search_direction(:,:,:,istate,ispin)
        end do
      end do
      previous_residual_norm=residual_norm

      do ispin=1,state%nspin
        do istate=1,state%nstate
          pnorm_local=sum(abs(search_direction(:,:,:,istate,ispin))**2)
          call global_real_sum(pnorm_local,pnorm_global,communicator)
          if (pnorm_global > tiny(1.0d0)) search_direction(:,:,:,istate,ispin)= &
            search_direction(:,:,:,istate,ispin)/sqrt(pnorm_global)
        end do
      end do
      psi_saved=state%psi_core
      state%psi_core=search_direction
      call apply_hamiltonian(state,hsearch)
      state%psi_core=psi_saved

      do ispin=1,state%nspin
        do istate=1,state%nstate
          h2=(0.0d0,0.0d0); s2=(0.0d0,0.0d0)
          h2(1,1)=cmplx(eigenvalues(istate,ispin),0.0d0,8)
          overlap_local=sum(conjg(state%psi_core(:,:,:,istate,ispin))*hsearch(:,:,:,istate,ispin))
          call global_complex_sum(overlap_local,overlap_global,communicator)
          h2(1,2)=overlap_global; h2(2,1)=conjg(overlap_global)
          overlap_local=sum(conjg(search_direction(:,:,:,istate,ispin))*hsearch(:,:,:,istate,ispin))
          call global_complex_sum(overlap_local,overlap_global,communicator)
          h2(2,2)=cmplx(real(overlap_global,8),0.0d0,8)
          s2(1,1)=(1.0d0,0.0d0); s2(2,2)=(1.0d0,0.0d0)
          call zhegv(1,'V','U',2,h2,2,s2,2,eval2,work,9,rwork,ierr)
          if (ierr /= 0) stop 'nodal DG: CG two-vector diagonalization failed'
          if (abs(h2(1,1)) > tiny(1.0d0)) then
            phase=conjg(h2(1,1))/abs(h2(1,1))
          else
            phase=(1.0d0,0.0d0)
          end if
          h2(:,1)=h2(:,1)*phase
          state%psi_core(:,:,:,istate,ispin)=h2(1,1)*state%psi_core(:,:,:,istate,ispin)+ &
            h2(2,1)*search_direction(:,:,:,istate,ispin)
          hpsi(:,:,:,istate,ispin)=h2(1,1)*hpsi(:,:,:,istate,ispin)+ &
            h2(2,1)*hsearch(:,:,:,istate,ispin)
          eigenvalues(istate,ispin)=eval2(1)
        end do
      end do
    end do
    niteration=min(iteration,max_iteration)
    call orthonormalize_nodal_states_mpi(state,communicator)
    call apply_hamiltonian(state,hpsi)
    call rotate_subspace(state,hpsi,eigenvalues,communicator)
    call build_residuals(state,hpsi,eigenvalues,residual,residual_norm,max_residual,communicator)
    if(max_residual>tolerance) then
      if(must_converge) stop 'nodal DG: callback CG eigensolver did not converge'
    else
      call accept_nodal_dg_ground_state(state,max_residual,tolerance)
    end if
    deallocate(hpsi,residual,search_direction,hsearch,psi_saved,residual_norm,previous_residual_norm)
  end subroutine solve_nodal_ground_state_cg_mpi

  subroutine build_residuals(state,hpsi,eigenvalues,residual,residual_norm,max_residual,communicator)
    type(s_dg_nodal_state), intent(in) :: state
    complex(8), intent(in) :: hpsi(:,:,:,:,:)
    real(8), intent(inout) :: eigenvalues(state%nstate,state%nspin)
    complex(8), intent(out) :: residual(:,:,:,:,:)
    real(8), intent(out) :: residual_norm(state%nstate,state%nspin),max_residual
    integer, intent(in) :: communicator
    complex(8) :: local_value,global_value
    real(8) :: local_norm,global_norm
    integer :: istate,ispin
    max_residual=0.0d0
    do ispin=1,state%nspin; do istate=1,state%nstate
      local_value=sum(conjg(state%psi_core(:,:,:,istate,ispin))*hpsi(:,:,:,istate,ispin))
      call global_complex_sum(local_value,global_value,communicator)
      eigenvalues(istate,ispin)=real(global_value,8)
      residual(:,:,:,istate,ispin)=hpsi(:,:,:,istate,ispin)- &
        eigenvalues(istate,ispin)*state%psi_core(:,:,:,istate,ispin)
      local_norm=sum(abs(residual(:,:,:,istate,ispin))**2)
      call global_real_sum(local_norm,global_norm,communicator)
      residual_norm(istate,ispin)=global_norm
      max_residual=max(max_residual,sqrt(global_norm))
    end do; end do
  end subroutine build_residuals

  subroutine global_real_sum(local_value,global_value,communicator)
    real(8), intent(in) :: local_value
    real(8), intent(out) :: global_value
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: CG real reduction failed'
#else
    global_value=local_value
#endif
  end subroutine global_real_sum

  subroutine global_complex_sum(local_value,global_value,communicator)
    complex(8), intent(in) :: local_value
    complex(8), intent(out) :: global_value
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: CG complex reduction failed'
#else
    global_value=local_value
#endif
  end subroutine global_complex_sum

end module rt_dg_nodal_cg
