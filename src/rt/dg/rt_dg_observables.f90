  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt)
    use salmon_global, only: nelec, theory
    use structures
    use communication, only: comm_summation
    use rt_dg_fragment_ops, only: apply_momentum_blocks
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    
    integer :: io, ispin, idir, n, nocc
    integer :: ifrag, jfrag, ib, jb, i_idx, j_idx
    logical :: do_interface_check
    real(8), allocatable :: interface_flow(:,:), dndt_frag(:)
    real(8) :: pair_residual, max_pair_residual, charge_balance_residual
    real(8) :: current_tmp, energy_tmp, pw_weight_local
    real(8) :: Ac_tot(3), A_squared
    real(8) :: current_local(3), energy_local
    complex(8) :: minus_i
    complex(8), allocatable :: op_mat(:,:), tmp_mat(:,:)
    logical :: has_nonlocal, use_hmat_complex
    logical :: use_spatial_A
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    ! Calculate local observables (only for assigned fragments)
    ! MPI aggregation will sum across all ranks
    current_local = 0.0d0
    energy_local = 0.0d0
    pw_weight_local = 0.0d0

    n = dg_frag%n_mat_max
    ! Occupied orbitals per spin channel (closed-shell): nelec = 2 * nocc.
    ! Using nelec/nspin here over-occupies states when nspin=1.
    nocc = min(int(nelec / 2.0d0 + 1.0d-12), min(dg_frag%nstate_tot, n))
    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)

    do_interface_check = .false.
    if (do_interface_check) then
      allocate(interface_flow(dg_frag%n_frag, dg_frag%n_frag))
      interface_flow = 0.0d0
    end if
    if (n <= 0) then
      current_local = 0.0d0
      energy_local = 0.0d0
      goto 1000
    end if

    allocate(op_mat(n, n), tmp_mat(n, nocc))
    minus_i = cmplx(0.0d0, -1.0d0, kind=8)

    ! Current calculation via momentum operator matrix (velocity gauge)
    ! Following conventional RT implementation in density_matrix.f90:
    !   - momentum_mat stores <φ|∇|φ> (gradient operator)
    !   - Current: j = Im[<ψ|∇|ψ>] with factor 2 and normalization by ngrid
    !   - Sign: Testing -2.0 to match conventional RT direction
    do ispin = 1, dg_frag%nspin
      do idir = 1, 3
        ! momentum_mat = <φ|∇|φ>, need to apply -i via aimag() and include factor 2
        if (allocated(dg_frag%momentum_blocks)) then
          tmp_mat(:, :) = (0.0d0, 0.0d0)
          Ac_tot = 0.0d0
          Ac_tot(idir) = 1.0d0
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, dg_frag%coef(1:n, 1:nocc, ispin), tmp_mat)
        else if (allocated(dg_frag%momentum_mat_c)) then
          op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     dg_frag%coef(1:n, 1:nocc, ispin), n, (0.0d0, 0.0d0), tmp_mat, n)
        else
          op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     dg_frag%coef(1:n, 1:nocc, ispin), n, (0.0d0, 0.0d0), tmp_mat, n)
        end if
        
        current_tmp = 0.0d0
        do io = 1, nocc
          ! Factor -2.0: -1 for operator sign convention, 2 for Im[ψ*∇ψ] normalization
          current_tmp = current_tmp + sum(aimag(conjg(dg_frag%coef(1:n, io, ispin)) * tmp_mat(1:n, io)))
        end do
        current_local(idir) = current_local(idir) - 2.0d0 * current_tmp
      end do
    end do
    
    ! Get vector potential at current time for energy calculation
    Ac_tot = rt%Ac_tot(:, itt)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    
    call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
    has_nonlocal = dg_frag%has_nl_cache

    ! Calculate total energy: E = <ψ|H(t)|ψ>
    ! H(t) = H_0 - i*A(t)·∇ + A²(t)/2 + V_NL(A)
    do ispin = 1, dg_frag%nspin
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      if (use_hmat_complex) then
        op_mat(:, :) = dg_frag%H_mat_c(1:n, 1:n, ispin)
      else
        op_mat(:, :) = cmplx(dg_frag%H_mat(1:n, 1:n, ispin), 0.0d0, kind=8)
      end if
      if (has_nonlocal) then
        op_mat(:, :) = op_mat(:, :) + dg_frag%H_nl_cache(1:n, 1:n, ispin)
      end if
      if (use_spatial_A) then
        if (.not. allocated(Ap_mat)) allocate(Ap_mat(n, n), A2_mat(n, n))
        call build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
        op_mat(:, :) = op_mat(:, :) + cmplx(A2_mat(:, :), 0.0d0, kind=8)
        op_mat(:, :) = op_mat(:, :) + minus_i * cmplx(Ap_mat(:, :), 0.0d0, kind=8)
      else
        do io = 1, n
          op_mat(io, io) = op_mat(io, io) + 0.5d0 * A_squared
        end do
        if (allocated(dg_frag%momentum_blocks)) then
          tmp_mat(:, :) = (0.0d0, 0.0d0)
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, dg_frag%coef(1:n, 1:nocc, ispin), tmp_mat)
          tmp_mat(:, :) = minus_i * tmp_mat(:, :)
        else
          do idir = 1, 3
            if (allocated(dg_frag%momentum_mat_c)) then
              op_mat(:, :) = op_mat(:, :) + minus_i * Ac_tot(idir) * dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
            else
              op_mat(:, :) = op_mat(:, :) + minus_i * Ac_tot(idir) * dg_frag%momentum_mat(idir, 1:n, 1:n, ispin)
            end if
          end do
        end if
      end if

      if (do_interface_check) then
        do ifrag = 1, dg_frag%n_frag
          do jfrag = 1, dg_frag%n_frag
            if (jfrag == ifrag) cycle
            do io = 1, nocc
              do ib = 1, dg_frag%n_basis(ifrag, ispin)
                i_idx = dg_frag%index_basis(ib, ifrag, ispin)
                if (i_idx < 1 .or. i_idx > n) cycle
                do jb = 1, dg_frag%n_basis(jfrag, ispin)
                  j_idx = dg_frag%index_basis(jb, jfrag, ispin)
                  if (j_idx < 1 .or. j_idx > n) cycle
                  interface_flow(ifrag, jfrag) = interface_flow(ifrag, jfrag) + &
                    2.0d0 * aimag(conjg(dg_frag%coef(i_idx, io, ispin)) * op_mat(i_idx, j_idx) * &
                                  dg_frag%coef(j_idx, io, ispin))
                end do
              end do
            end do
          end do
        end do
      end if

      if (.not. allocated(dg_frag%momentum_blocks) .or. use_spatial_A) then
        call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                   dg_frag%coef(1:n, 1:nocc, ispin), n, (0.0d0, 0.0d0), tmp_mat, n)
      else
        call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                   dg_frag%coef(1:n, 1:nocc, ispin), n, (1.0d0, 0.0d0), tmp_mat, n)
      end if
      energy_tmp = 0.0d0
      do io = 1, nocc
        energy_tmp = energy_tmp + sum(real(conjg(dg_frag%coef(1:n, io, ispin)) * tmp_mat(1:n, io)))
      end do
      energy_local = energy_local + energy_tmp
    end do

    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
      do ispin = 1, dg_frag%nspin
        do io = 1, nocc
          pw_weight_local = pw_weight_local + sum(abs(dg_frag%coef_pw(:, io, ispin))**2)
        end do
      end do
    end if

    if (do_interface_check) then
      allocate(dndt_frag(dg_frag%n_frag))
      dndt_frag = 0.0d0
      max_pair_residual = 0.0d0
      do ifrag = 1, dg_frag%n_frag
        do jfrag = 1, dg_frag%n_frag
          if (jfrag == ifrag) cycle
          dndt_frag(ifrag) = dndt_frag(ifrag) - interface_flow(ifrag, jfrag)
        end do
      end do

      do ifrag = 1, dg_frag%n_frag - 1
        do jfrag = ifrag + 1, dg_frag%n_frag
          pair_residual = interface_flow(ifrag, jfrag) + interface_flow(jfrag, ifrag)
          max_pair_residual = max(max_pair_residual, abs(pair_residual))
        end do
      end do
      charge_balance_residual = abs(sum(dndt_frag))

      deallocate(dndt_frag, interface_flow)
    end if

    if (allocated(Ap_mat)) deallocate(Ap_mat)
    if (allocated(A2_mat)) deallocate(A2_mat)
    deallocate(op_mat, tmp_mat)


    ! Cache retained for reuse

  1000 continue
    
    ! MPI aggregation: sum local contributions from all ranks
    call comm_summation(current_local, dg_frag%current, 3, dg_frag%icomm)
    call comm_summation(energy_local, dg_frag%total_energy, dg_frag%icomm)
    call comm_summation(pw_weight_local, dg_frag%pw_weight_raw, dg_frag%icomm)

    ! coef/momentum/H are synchronized globally on all ranks. Each rank therefore
    ! evaluates the same full-system observable; use process-average to avoid
    ! fragment-count-dependent overcounting after MPI summation.
    dg_frag%current(:) = dg_frag%current(:) / real(max(1, dg_frag%isize), 8)
    dg_frag%total_energy = dg_frag%total_energy / real(max(1, dg_frag%isize), 8)
    dg_frag%pw_weight_raw = dg_frag%pw_weight_raw / real(max(1, dg_frag%isize), 8)
    
    ! Normalize by global grid count exactly as conventional calc_current().
    ! This avoids decomposition-dependent scaling from local/grid-view differences.
    dg_frag%current(:) = dg_frag%current(:) / dble(system%ngrid)
    
    ! Store in rt structure for output
    rt%curr(:, itt) = dg_frag%current(:)
    
  end subroutine calculate_observables
