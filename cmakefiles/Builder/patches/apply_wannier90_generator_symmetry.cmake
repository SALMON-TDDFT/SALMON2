if(NOT DEFINED WANNIER90_SOURCE_DIR)
  message(FATAL_ERROR "WANNIER90_SOURCE_DIR is required")
endif()

set(source "${WANNIER90_SOURCE_DIR}/src/sitesym.F90")
file(READ "${source}" contents)

set(old_use "    use w90_parameters, only: num_wann, num_kpts\n")
set(new_use "    use w90_parameters, only: num_wann, num_kpts, symmetrize_eps\n")
string(REPLACE "${old_use}" "${new_use}" contents "${contents}")

set(old_decl "    integer :: ik, ir, isym, irk, ngk\n\n    complex(kind=dp) :: grad_total(num_wann, num_wann)\n")
set(new_decl "    integer :: ik, ir, isym, irk, ngk, iter\n    integer, parameter :: generator_projection_iterations = 100\n    real(kind=dp) :: generator_projection_diff\n\n    complex(kind=dp) :: grad_total(num_wann, num_wann)\n    complex(kind=dp) :: grad_previous(num_wann, num_wann)\n")
string(REPLACE "${old_decl}" "${new_decl}" contents "${contents}")

set(old_loop [=[    do ir = 1, nkptirr
      ik = ir2ik(ir)
      ngk = count(kptsym(:, ir) .eq. ik)
      if (ngk .eq. 1) cycle
      grad_total = grad(:, :, ik)
      do isym = 2, nsymmetry
        if (kptsym(isym, ir) .ne. ik) cycle
        !
        ! calculate cmat1 = D^{+}(R,k) G(Rk) D(R,k)
        !
        ! step 1: cmat2 =  G(Rk) D(R,k)
        call utility_zgemm(cmat2, grad(:, :, ik), 'N', &
                           d_matrix_wann(:, :, isym, ir), 'N', num_wann)
        ! step 2: cmat1 = D^{+}(R,k) * cmat2
        call utility_zgemm(cmat1, d_matrix_wann(:, :, isym, ir), 'C', &
                           cmat2, 'N', num_wann)
        grad_total = grad_total + cmat1
      enddo
      grad(:, :, ik) = grad_total/ngk
    enddo
]=])
set(new_loop [=[    do ir = 1, nkptirr
      ik = ir2ik(ir)
      ngk = count(kptsym(:, ir) .eq. ik)
      if (ngk .eq. 1) cycle
      do iter = 1, generator_projection_iterations
        grad_previous = grad(:, :, ik)
        grad_total = grad_previous
        do isym = 2, nsymmetry
          if (kptsym(isym, ir) .ne. ik) cycle
          ! Repeated averaging over identity plus a generating set converges
          ! to the same common fixed space as the full finite group Reynolds
          ! projector, without retaining every group element.
          call utility_zgemm(cmat2, grad_previous, 'N', &
                             d_matrix_wann(:, :, isym, ir), 'N', num_wann)
          call utility_zgemm(cmat1, d_matrix_wann(:, :, isym, ir), 'C', &
                             cmat2, 'N', num_wann)
          grad_total = grad_total + cmat1
        enddo
        grad_total = grad_total/ngk
        generator_projection_diff = sum(abs(grad_total - grad_previous))
        grad(:, :, ik) = grad_total
        if (generator_projection_diff .lt. symmetrize_eps) exit
      enddo
      if (iter .gt. generator_projection_iterations) then
        write (stdout, "(a,2e20.10)") &
          'generator symmetry gradient projection did not converge: diff,eps=', &
          generator_projection_diff, symmetrize_eps
        call io_error('sitesym_symmetrize_gradient: generator projection not converged')
      endif
    enddo
]=])
string(REPLACE "${old_loop}" "${new_loop}" contents "${contents}")

# Large symmetry-adapted retained spaces can converge monotonically but need
# more than the upstream fixed cap of 100 projection iterations.  Si64 reaches
# only O(1e-8) by iteration 100 and satisfies the requested 1e-10 tolerance
# when the same iteration is allowed to continue.
set(old_sitesym_iteration_cap "    integer, parameter :: niter = 100\n")
set(new_sitesym_iteration_cap "    integer, parameter :: niter = 1000\n")
string(REPLACE "${old_sitesym_iteration_cap}" "${new_sitesym_iteration_cap}"
       contents "${contents}")

if(contents MATCHES "integer :: ik, ir, isym, irk, ngk\\n" OR
   contents MATCHES "grad\\(:, :, ik\\) = grad_total/ngk" OR
   contents MATCHES "integer, parameter :: niter = 100\\n")
  message(FATAL_ERROR "Failed to patch Wannier90 generator symmetry projection")
endif()

file(WRITE "${source}" "${contents}")

set(wannierise_source "${WANNIER90_SOURCE_DIR}/src/wannierise.F90")
file(READ "${wannierise_source}" wannierise_contents)

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
string(REPLACE "${old_library_kmesh_use}" "${new_library_kmesh_use}"
       library_contents "${library_contents}")

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
