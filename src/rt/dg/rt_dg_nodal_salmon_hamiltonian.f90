!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Complete SALMON nodal Hamiltonian action used by GS and RT adapters.
!=======================================================================
module rt_dg_nodal_salmon_hamiltonian
  use structures, only: s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use rt_dg_nodal_mpi, only: exchange_nodal_face_halos
  use rt_dg_nodal_hamiltonian, only: apply_nodal_local_hamiltonian
  use rt_dg_nodal_nonlocal, only: apply_nodal_nonlocal_potential_mpi
  implicit none
  private
  public :: apply_nodal_salmon_hamiltonian_mpi
  logical, save :: trace_first_action = .true.

contains

  subroutine apply_nodal_salmon_hamiltonian_mpi(state,dg_frag,ppg,ik,vlocal,kinetic_center, &
                                                 kinetic_offset,gradient_offset,avec,communicator,hpsi)
    type(s_dg_nodal_state), intent(inout) :: state
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik, communicator
    real(8), intent(in) :: vlocal(:,:,:,:), kinetic_center
    real(8), intent(in) :: kinetic_offset(:,:), gradient_offset(:,:), avec(3)
    complex(8), intent(out) :: hpsi(:,:,:,:,:)
    real(8) :: t0,t1,t2,t3

    call cpu_time(t0)
    call exchange_nodal_face_halos(state,communicator)
    call cpu_time(t1)
    if (trace_first_action .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe13.5)') '[DG-NODAL-H-PROFILE] halo=',t1-t0
      flush(6)
    end if
    call apply_nodal_local_hamiltonian(state,vlocal,kinetic_center,kinetic_offset,gradient_offset,avec,hpsi)
    call cpu_time(t2)
    if (trace_first_action .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe13.5)') '[DG-NODAL-H-PROFILE] local=',t2-t1
      flush(6)
    end if
    call apply_nodal_nonlocal_potential_mpi(state,dg_frag,ppg,ik,communicator,hpsi)
    call cpu_time(t3)
    if (trace_first_action .and. dg_frag%id == 0) then
      write(*,'(1x,a,3(a,1pe13.5))') '[DG-NODAL-H-PROFILE]', &
        ' halo=',t1-t0,' local=',t2-t1,' nonlocal=',t3-t2
      flush(6)
    end if
    trace_first_action=.false.
  end subroutine apply_nodal_salmon_hamiltonian_mpi

end module rt_dg_nodal_salmon_hamiltonian
