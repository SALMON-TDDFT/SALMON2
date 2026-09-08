if(NOT DEFINED WANNIER90_SOURCE_DIR)
  message(FATAL_ERROR "WANNIER90_SOURCE_DIR is required")
endif()

set(source "${WANNIER90_SOURCE_DIR}/src/sitesym.F90")
file(READ "${source}" contents)

# The DMN stream contains every stabilizer operation, including identity.
# Keep Wannier90's exact one-pass Reynolds average.  Repeatedly applying this
# already-idempotent full-group projector only adds dense GEMMs and roundoff.
if(NOT contents MATCHES "grad\\(:, :, ik\\) = grad_total/ngk")
  message(FATAL_ERROR "Wannier90 one-pass full-group gradient projection is missing")
endif()

# Large symmetry-adapted retained spaces can converge monotonically but need
# more than the upstream fixed cap of 100 projection iterations.  Si64 reaches
# only O(1e-8) by iteration 100 and satisfies the requested 1e-10 tolerance
# when the same iteration is allowed to continue.
set(old_sitesym_iteration_cap "    integer, parameter :: niter = 100\n")
set(new_sitesym_iteration_cap "    integer, parameter :: niter = 1000\n")
string(REPLACE "${old_sitesym_iteration_cap}" "${new_sitesym_iteration_cap}"
       contents "${contents}")

# symmetrize_eps is an elementwise representation tolerance.  Summing the
# whole residual matrix makes the acceptance threshold grow as num_wann**2
# and rejects large, otherwise compatible retained spaces.
string(REPLACE "      diff = sum(abs(cmat2))\n"
               "      diff = maxval(abs(cmat2))\n"
               contents "${contents}")

if(contents MATCHES "integer, parameter :: niter = 100\\n" OR
   contents MATCHES "diff = sum\\(abs\\(cmat2\\)\\)")
  message(FATAL_ERROR "Failed to patch Wannier90 generator symmetry projection")
endif()

file(WRITE "${source}" "${contents}")

set(wannierise_source "${WANNIER90_SOURCE_DIR}/src/wannierise.F90")
file(READ "${wannierise_source}" wannierise_contents)

set(old_line_search_state [=[    real(kind=dp) :: doda0
]=])
set(new_line_search_state [=[    real(kind=dp) :: doda0, symmetry_backtracking_step, symmetry_backtracking_tolerance
    real(kind=dp) :: symmetry_backtracking_min_step
    logical :: symmetry_backtracking_accepted
]=])
string(REPLACE "${old_line_search_state}" "${new_line_search_state}"
       wannierise_contents "${wannierise_contents}")

set(old_investigation_state [=[    real(kind=dp) :: doda0, unprojected_doda0, projected_doda0, projected_direction_norm, projected_antihermitian_defect
]=])
string(REPLACE "${old_investigation_state}" "${new_line_search_state}"
       wannierise_contents "${wannierise_contents}")
string(REPLACE "      unprojected_doda0 = doda0\n" ""
       wannierise_contents "${wannierise_contents}")

set(old_gradient_projection [=[      if (lsitesymmetry) call sitesym_symmetrize_gradient(2, cdq) !RS:
]=])
set(new_gradient_projection [=[      if (lsitesymmetry) then
        ! internal_search_direction updates the distributed cdq_loc array.
        ! Synchronize it before projection and return the projected owned slice
        ! so that the line search uses the symmetry-preserving direction.
        call comms_gatherv(cdq_loc, num_wann*num_wann*counts(my_node_id), &
                           cdq, num_wann*num_wann*counts, num_wann*num_wann*displs)
        call comms_bcast(cdq(1, 1, 1), num_wann*num_wann*num_kpts)
        call sitesym_symmetrize_gradient(2, cdq)
        cdq_loc(:, :, 1:counts(my_node_id)) = cdq(:, :, 1 + displs(my_node_id): &
          displs(my_node_id) + counts(my_node_id))
      endif
]=])
string(REPLACE "${old_gradient_projection}" "${new_gradient_projection}"
       wannierise_contents "${wannierise_contents}")

set(old_investigation_projection [=[        cdq_loc(:, :, 1:counts(my_node_id)) = cdq(:, :, 1 + displs(my_node_id): &
          displs(my_node_id) + counts(my_node_id))
        ! The projection changes the direction after internal_search_direction
        ! has evaluated dOmega/dalpha.  Recompute the slope from the direction
        ! that the trial step will actually apply.
        projected_doda0 = -real(sum(conjg(cdodq_loc(:, :, 1:counts(my_node_id)))* &
          cdq_loc(:, :, 1:counts(my_node_id))), dp)
        call comms_allreduce(projected_doda0, 1, 'SUM')
        doda0 = projected_doda0/(4.0_dp*wbtot)
        projected_direction_norm = sum(abs(cdq_loc(:, :, 1:counts(my_node_id)))**2)
        call comms_allreduce(projected_direction_norm, 1, 'SUM')
        projected_antihermitian_defect = 0.0_dp
        do nkp_loc = 1, counts(my_node_id)
          do n = 1, num_wann
            do i = 1, num_wann
              projected_antihermitian_defect = projected_antihermitian_defect + &
                abs(cdq_loc(i, n, nkp_loc) + conjg(cdq_loc(n, i, nkp_loc)))**2
            enddo
          enddo
        enddo
        call comms_allreduce(projected_antihermitian_defect, 1, 'SUM')
        projected_antihermitian_defect = sqrt(projected_antihermitian_defect/ &
          max(projected_direction_norm, tiny(1.0_dp)))
        if (lprint .and. iprint > 2 .and. on_root) then
          write (stdout, *) ' LINE --> Pre-symmetry slope             :', unprojected_doda0*lenconfac**2
          write (stdout, *) ' LINE --> Symmetry-projected slope       :', doda0*lenconfac**2
          write (stdout, *) ' LINE --> ||projected direction||^2      :', &
            projected_direction_norm*lenconfac**2
          write (stdout, *) ' LINE --> Projected anti-Hermitian defect:', &
            projected_antihermitian_defect
          write (stdout, *) ' LINE --> Line-search wbtot              :', wbtot
        endif
]=])
set(new_investigation_projection [=[        cdq_loc(:, :, 1:counts(my_node_id)) = cdq(:, :, 1 + displs(my_node_id): &
          displs(my_node_id) + counts(my_node_id))
]=])
string(REPLACE "${old_investigation_projection}" "${new_investigation_projection}"
       wannierise_contents "${wannierise_contents}")

set(old_line_search_initialization [=[    gcnorm1 = 0.0_dp; gcnorm0 = 0.0_dp
]=])
set(new_line_search_initialization [=[    gcnorm1 = 0.0_dp; gcnorm0 = 0.0_dp
    symmetry_backtracking_step = trial_step
]=])
string(REPLACE "${old_line_search_initialization}" "${new_line_search_initialization}"
       wannierise_contents "${wannierise_contents}")

set(old_trial_line_search [=[        ! take trial step
        cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(trial_step/(4.0_dp*wbtot))

        ! store original U and M before rotating
        u0_loc = u_matrix_loc

        if (optimisation <= 0) then
!             write(page_unit)   m_matrix
          write (page_unit) m_matrix_loc
          rewind (page_unit)
        else
          m0_loc = m_matrix_loc
        endif

        ! update U and M
        call internal_new_u_and_m()

        ! calculate spread at trial step
        call wann_omega(csheet, sheet, rave, r2ave, rave2, trial_spread)

        ! Calculate optimal step (alphamin)
        call internal_optimal_step()
]=])
set(new_trial_line_search [=[        ! Store the original state once and reuse it for every retry.
        u0_loc = u_matrix_loc
        if (optimisation <= 0) then
          write (page_unit) m_matrix_loc
          rewind (page_unit)
        else
          m0_loc = m_matrix_loc
        endif

        if (lsitesymmetry) then
          symmetry_backtracking_min_step = sqrt(epsilon(1.0_dp))*max(1.0_dp, trial_step)
          symmetry_backtracking_tolerance = 16.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, abs(wann_spread%om_tot))
          symmetry_backtracking_accepted = .false.
          do
            u_matrix_loc = u0_loc
            if (optimisation <= 0) then
              read (page_unit) m_matrix_loc
              rewind (page_unit)
            else
              m_matrix_loc = m0_loc
            endif
            cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)* &
              (symmetry_backtracking_step/(4.0_dp*wbtot))
            call internal_new_u_and_m()
            call wann_omega(csheet, sheet, rave, r2ave, rave2, trial_spread)
            if (trial_spread%om_tot <= wann_spread%om_tot + symmetry_backtracking_tolerance) then
              symmetry_backtracking_accepted = .true.
              exit
            endif
            if (symmetry_backtracking_step <= symmetry_backtracking_min_step) exit
            symmetry_backtracking_step = 0.5_dp*symmetry_backtracking_step
          enddo
          if (.not. symmetry_backtracking_accepted) then
            u_matrix_loc = u0_loc
            if (optimisation <= 0) then
              read (page_unit) m_matrix_loc
              rewind (page_unit)
            else
              m_matrix_loc = m0_loc
            endif
            symmetry_backtracking_step = 0.0_dp
            ncg = 0
            gcfac = 0.0_dp
            call wann_omega(csheet, sheet, rave, r2ave, rave2, trial_spread)
          endif
          alphamin = symmetry_backtracking_step
          falphamin = trial_spread%om_tot
          lquad = .false.
        else
          cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(trial_step/(4.0_dp*wbtot))
          call internal_new_u_and_m()
          call wann_omega(csheet, sheet, rave, r2ave, rave2, trial_spread)
          call internal_optimal_step()
        endif
]=])
string(REPLACE "${old_trial_line_search}" "${new_trial_line_search}"
       wannierise_contents "${wannierise_contents}")

if(NOT wannierise_contents MATCHES "Accepted symmetry step")
  string(REPLACE
    "          lquad = .false.\n        else\n          cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(trial_step/(4.0_dp*wbtot))\n"
    "          lquad = .false.\n          if (lprint .and. iprint > 2 .and. on_root) &\n            write (stdout, *) ' LINE --> Accepted symmetry step         :', symmetry_backtracking_step\n        else\n          cdq_loc(:, :, :) = cdqkeep_loc(:, :, :)*(trial_step/(4.0_dp*wbtot))\n"
    wannierise_contents "${wannierise_contents}")
endif()

# Accelerate's complex-return BLAS ABI is not compatible with gfortran's
# external ZDOTC declaration on Apple Silicon.  These are simple local
# Frobenius inner products, so express them directly and keep GEMM on BLAS.
set(old_preconditioned_norm [=[        gcnorm1 = real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_precond_loc, 1, cdodq_loc, 1), dp)
]=])
set(new_preconditioned_norm [=[        gcnorm1 = real(sum(conjg(cdodq_precond_loc(:, :, 1:counts(my_node_id)))* &
                           cdodq_loc(:, :, 1:counts(my_node_id))), dp)
]=])
string(REPLACE "${old_preconditioned_norm}" "${new_preconditioned_norm}"
       wannierise_contents "${wannierise_contents}")

set(old_gradient_norm [=[        gcnorm1 = real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_loc, 1, cdodq_loc, 1), dp)
]=])
set(new_gradient_norm [=[        gcnorm1 = sum(abs(cdodq_loc(:, :, 1:counts(my_node_id)))**2)
]=])
string(REPLACE "${old_gradient_norm}" "${new_gradient_norm}"
       wannierise_contents "${wannierise_contents}")

set(old_directional_derivative [=[      doda0 = -real(zdotc(counts(my_node_id)*num_wann*num_wann, cdodq_loc, 1, cdq_loc, 1), dp)
]=])
set(new_directional_derivative [=[      doda0 = -real(sum(conjg(cdodq_loc(:, :, 1:counts(my_node_id)))* &
                         cdq_loc(:, :, 1:counts(my_node_id))), dp)
]=])
string(REPLACE "${old_directional_derivative}" "${new_directional_derivative}"
       wannierise_contents "${wannierise_contents}")

set(old_optimal_step [=[      if (abs(eqa/(fac*wann_spread%om_tot)) .gt. epsilon(1.0_dp)) then
        lquad = .true.
]=])
set(new_optimal_step [=[      if (abs(fac*wann_spread%om_tot) .gt. tiny(1.0_dp)) then
        lquad = abs(eqa/(fac*wann_spread%om_tot)) .gt. epsilon(1.0_dp)
      else
        lquad = .false.
      endif
      if (lquad) then
]=])
string(REPLACE "${old_optimal_step}" "${new_optimal_step}"
       wannierise_contents "${wannierise_contents}")

if(wannierise_contents MATCHES
   "if \(lsitesymmetry\) call sitesym_symmetrize_gradient\(2, cdq\)")
  message(FATAL_ERROR "Failed to patch Wannier90 distributed gradient synchronization")
endif()
if(wannierise_contents MATCHES
   "abs\(eqa/\(fac\*wann_spread%om_tot\)\)")
  message(FATAL_ERROR "Failed to patch Wannier90 zero-spread line search")
endif()
file(WRITE "${wannierise_source}" "${wannierise_contents}")

set(overlap_source "${WANNIER90_SOURCE_DIR}/src/overlap.F90")
file(READ "${overlap_source}" overlap_contents)

set(old_gamma_import "      m_matrix_orig, u_matrix_opt, cp_pp, use_bloch_phases, gamma_only, & ![ysl]\n")
set(new_gamma_import "      m_matrix_orig, u_matrix_opt, cp_pp, use_bloch_phases, gamma_only, lsitesymmetry, & ![ysl]\n")
string(REPLACE "${old_gamma_import}" "${new_gamma_import}"
       overlap_contents "${overlap_contents}")

set(old_gamma_branch "      if (.not. gamma_only) then\n")
set(new_gamma_branch "      if (.not. gamma_only .or. lsitesymmetry) then\n")
string(REPLACE "${old_gamma_branch}" "${new_gamma_branch}"
       overlap_contents "${overlap_contents}")

if(overlap_contents MATCHES "if \\(\\.not\\. gamma_only\\) then" OR
   NOT overlap_contents MATCHES "gamma_only, lsitesymmetry")
  message(FATAL_ERROR "Failed to patch Wannier90 symmetry-adapted Gamma initialization")
endif()

file(WRITE "${overlap_source}" "${overlap_contents}")

set(library_source "${WANNIER90_SOURCE_DIR}/src/wannier_lib.F90")
file(READ "${library_source}" library_contents)

set(old_library_kmesh_use "  use w90_kmesh\n")
set(new_library_kmesh_use "  use w90_kmesh\n  use w90_sitesym, only: sitesym_read\n")
if(NOT library_contents MATCHES "use w90_sitesym, only: sitesym_read")
  string(REPLACE "${old_library_kmesh_use}" "${new_library_kmesh_use}"
         library_contents "${library_contents}")
endif()

set(old_library_kmesh_call [=[  call kmesh_get()

  time2 = io_time()
]=])
set(new_library_kmesh_call [=[  call kmesh_get()
  if (lsitesymmetry) call sitesym_read()

  time2 = io_time()
]=])
string(REPLACE "${old_library_kmesh_call}" "${new_library_kmesh_call}"
       library_contents "${library_contents}")

set(old_library_projection [=[    if (gamma_only) then
      call overlap_project_gamma()
    else
      call overlap_project()
    endif
]=])
set(new_library_projection [=[    if (gamma_only .and. .not. lsitesymmetry) then
      call overlap_project_gamma()
    else
      call overlap_project()
    endif
]=])
string(REPLACE "${old_library_projection}" "${new_library_projection}"
       library_contents "${library_contents}")

set(old_library_optimizer [=[  if (gamma_only) then
    call wann_main_gamma()
  else
    call wann_main()
  endif
]=])
set(new_library_optimizer [=[  if (gamma_only .and. .not. lsitesymmetry) then
    call wann_main_gamma()
  else
    call wann_main()
  endif
]=])
string(REPLACE "${old_library_optimizer}" "${new_library_optimizer}"
       library_contents "${library_contents}")

if(NOT library_contents MATCHES
   "if \\(gamma_only \\.and\\. \\.not\\. lsitesymmetry\\) then" OR
   NOT library_contents MATCHES "if \\(lsitesymmetry\\) call sitesym_read")
  message(FATAL_ERROR "Failed to patch Wannier90 library Gamma symmetry dispatch")
endif()

file(WRITE "${library_source}" "${library_contents}")

set(program_source "${WANNIER90_SOURCE_DIR}/src/wannier_prog.F90")
file(READ "${program_source}" program_contents)

set(old_program_optimizer "  if (.not. gamma_only) then\n")
set(new_program_optimizer "  if (.not. gamma_only .or. lsitesymmetry) then\n")
string(REPLACE "${old_program_optimizer}" "${new_program_optimizer}"
       program_contents "${program_contents}")
if(NOT program_contents MATCHES
   "if \\(\\.not\\. gamma_only \\.or\\. lsitesymmetry\\) then")
  message(FATAL_ERROR "Failed to patch Wannier90 executable Gamma symmetry dispatch")
endif()

file(WRITE "${program_source}" "${program_contents}")
