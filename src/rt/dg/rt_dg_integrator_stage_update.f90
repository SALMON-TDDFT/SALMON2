  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_global, only: yn_dual_rho_vh_only
    use communication, only: comm_summation, comm_get_min
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use hartree_sub, only: build_hartree_density_from_rho
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: Ac_tot(3)
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_xc_functional),  intent(in)    :: xc_func
    type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_pp_nlcc),        intent(in)    :: ppn
    type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
    type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
    type(s_dft_energy),     intent(inout) :: energy
    integer :: n_frag, n_pw, ispin
    integer :: n_floor_local, n_floor_global
    real(8) :: t_stage0, t_stage1
    real(8) :: rho_floor_min_local, rho_floor_min_global
    real(8) :: rho_floor_min_buf(1)
    real(8), allocatable :: rho_buffer(:,:,:), Vh_buffer(:,:,:)
    logical :: use_rank_buffered_potential
    logical :: use_dual_rho_vh
    logical, parameter :: enable_stage_update_trace = .false.
    logical, parameter :: enable_stage_update_progress = .false.
    type(s_scalar) :: rho_h

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        density-hmat stage trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if

    call cpu_time(t_stage0)
    if (itt == 1) then
      write(*,'(1x,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_DENSITY_RECON_CALL rank=', dg_frag%id, ' itt=', itt, &
        ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt)
    if (itt == 1) then
      write(*,'(1x,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_DENSITY_RECON_CALL rank=', dg_frag%id, ' itt=', itt, &
        ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    call cpu_time(t_stage1)

    dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (system%nspin > 0) then
      dg_frag%rho_s_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_frag = 1, system%nspin
        dg_frag%rho_s_frag(:, :, :, n_frag) = rho_s(n_frag)%f(:, :, :)
      end do
    end if

    use_dual_rho_vh = (yn_dual_rho_vh_only == 'y')
    if (use_dual_rho_vh) then
      call allocate_scalar(mg, rho_h)
      if (itt == 1 .and. dg_frag%has_seed_rho_h) then
        rho_h%f(:, :, :) = dg_frag%rho_h_frag(:, :, :)
      else
        call build_hartree_density_from_rho(info, rho, rho_h)
        dg_frag%rho_h_frag(:, :, :) = rho_h%f(:, :, :)
      end if
    end if

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    if (use_rank_buffered_potential) then
      allocate(rho_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, rho, rho_buffer)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-density dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-hartree"
      flush(6)
    end if
    call cpu_time(t_stage0)
    if (use_dual_rho_vh) then
      call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho_h, Vh)
    else
      call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    end if
    call cpu_time(t_stage1)
    dg_frag%Vh_frag(:, :, :) = Vh%f(:, :, :)
    if (use_rank_buffered_potential) then
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vh, Vh_buffer)
    end if
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-hartree dt=", t_stage1 - t_stage0
      flush(6)
    end if
    if (use_dual_rho_vh) then
      call deallocate_scalar(rho_h)
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-xc"
      flush(6)
    end if
    ! Guard against tiny negative rho before XC; this prevents non-physical Vxc NaN under aggressive FP reassociation.
    ! Validation note: O3 and O3+no-unsafe-math runs matched key observables (J_para/J_total/Ne_raw) in this case.
    n_floor_local = 0
    rho_floor_min_local = huge(1.0d0)
    do ispin = 1, system%nspin
      n_floor_local = n_floor_local + count(rho_s(ispin)%f(:, :, :) < 0.0d0)
      if (any(rho_s(ispin)%f(:, :, :) < 0.0d0)) then
        rho_floor_min_local = min(rho_floor_min_local, minval(rho_s(ispin)%f(:, :, :), &
          mask=rho_s(ispin)%f(:, :, :) < 0.0d0))
      end if
      where (rho_s(ispin)%f(:, :, :) < 0.0d0) rho_s(ispin)%f(:, :, :) = 0.0d0
    end do
    call comm_summation(n_floor_local, n_floor_global, dg_frag%icomm)
    if (n_floor_global > 0) then
      rho_floor_min_buf(1) = rho_floor_min_local
      call comm_get_min(rho_floor_min_buf, rho_floor_min_buf, 1, dg_frag%icomm)
      rho_floor_min_global = rho_floor_min_buf(1)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a,1pe14.6)') '[RHO-FLOOR] itt=', itt, ' clamped_points=', n_floor_global, &
          ' min_rho_before=', rho_floor_min_global
        flush(6)
      end if
    end if
    call cpu_time(t_stage0)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    call cpu_time(t_stage1)
    if (system%nspin > 0) then
      dg_frag%Vxc_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_frag = 1, system%nspin
        dg_frag%Vxc_frag(:, :, :, n_frag) = Vxc(n_frag)%f(:, :, :)
      end do
    end if
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-xc dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if (dg_frag%use_plusu .and. PLUS_U_ON) then
      call calc_density_matrix_and_energy_plusU(rt%tpsi0, ppg, info, system, energy%E_U)
    else
      energy%E_U = 0.0d0
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=before-reconstruct"
      flush(6)
    end if
    call cpu_time(t_stage0)
    if (use_rank_buffered_potential) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, Vh_buffer)
    else
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    end if
    call cpu_time(t_stage1)
    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') "        density-hmat stage trace: stage=after-reconstruct dt=", t_stage1 - t_stage0
      flush(6)
    end if

    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      n_frag = dg_frag%n_mat_max
      n_pw = dg_frag%n_plane_waves

      if (.not. allocated(dg_frag%S_mat_frag_pw)) then
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      else if (size(dg_frag%S_mat_frag_pw,1) /= n_frag .or. size(dg_frag%S_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%S_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%S_mat_frag_pw)
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      end if

      if (.not. allocated(dg_frag%H_mat_frag_pw)) then
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%H_mat_frag_pw)
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      end if

      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, dg_frag%H_mat_frag_pw)
      call build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, dg_frag%S_mat_frag_pw, dg_frag%H_mat_frag_pw)
      if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        density-hmat stage trace: stage=after-mixed-refresh"
        flush(6)
      end if
    end if

    if ((enable_stage_update_trace .or. enable_stage_update_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density-hmat stage trace: stage=exit"
      flush(6)
    end if

    if (allocated(rho_buffer)) deallocate(rho_buffer)
    if (allocated(Vh_buffer)) deallocate(Vh_buffer)

  end subroutine update_density_hamiltonian_stage
