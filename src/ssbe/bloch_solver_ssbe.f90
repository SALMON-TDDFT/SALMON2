module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use communication, only: comm_get_groupinfo, comm_summation, &
                             comm_isend, comm_irecv, comm_wait_all
    use gs_info_ssbe
    use util_ssbe, only: split_range
    use sbe_lg_mode_ssbe, only: uses_prod_dk, uses_xfull_links, uses_fd_gicov, uses_integral_gicov
    implicit none



    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction
        ! VG completion (yn_sbe_vnl_exact='y'): replace the first-order A.rvnl
        ! nonlocal correction by the all-order DeltaV = V_nl(k+A)-V_nl(k)
        ! (propagation commutator) and W_i(k+A) (current readout), interpolated
        ! (centered 4-point Lagrange) from the gs%vnl_V/vnl_W kappa-stencil.
        ! Mutually exclusive with flag_vnl_correction (checker-enforced).
        logical :: flag_vnl_exact = .false.
        ! Window f-sum-rule deficiency tensor D_ij of the truncated band basis
        ! (yn_sbe_gs_current_subtract='y'; see build_fsum_deficiency_tensor and
        ! the calc_current_bloch header for the derivation).  Built ONCE per
        ! solver at init; the velocity-gauge readout then applies
        ! J(t) -= matmul(fsum_D, A(t)) -- per-step cost is one 3x3 matvec.
        real(8) :: fsum_D(1:3, 1:3) = 0d0
        logical :: fsum_D_built = .false.
        complex(8), allocatable :: qnm(:, :, :)
        complex(8), allocatable :: grad_qnm(:, :, :, :)
        complex(8), allocatable :: qnm_new(:, :, :)
        complex(8), allocatable :: dqnm_stock(:, :, :, :)
        complex(8), allocatable :: dnm_i(:, :, :, :)
        real(8), allocatable    :: abs_dnm(:, :, :)
        complex(8), allocatable :: exp_iphi(:, :, :)
        ! T2 Delta-omega dephasing gate weight (design: gw/gw_design/plans/
        ! 2026-07-08-t2-gate-shape.md), precomputed ONCE at init from the
        ! STATIC gs%delta_omega (t2_gate_weight, gs_info_ssbe.f90) -- gicov_rhs
        ! then just Hadamard-multiplies it into the phenomenological -rho/t_2
        ! term every step.  Diagonal forced to 0 (no population dephasing).
        real(8), allocatable    :: t2_gate_w(:, :, :)
    end type

    ! --- lightweight gicov perf timers (accumulate in gicov_rhs; print via sbe_print_timers) ---
    real(8),    save :: sbe_tmr_cov  = 0d0    ! covariant_grad_block wall [s], summed over RHS calls
    real(8),    save :: sbe_tmr_comm = 0d0    ! rho halo-exchange (was: full allreduce) wall [s]
    integer(8), save :: sbe_tmr_nrhs = 0      ! number of gicov_rhs calls

    ! --- TEST hook (test_gicov_hex.f90 orthogonal-equivalence check): force the
    ! GENERAL (non-orthogonal, reduced-axis c_i = (E.a_i)/2pi) field-term path
    ! in gicov_rhs even when b_matrix is strictly diagonal, so the two paths can
    ! be compared on the same fixture.  NEVER set in production: the default
    ! .false. keeps every diagonal-b campaign on the legacy branch bit-for-bit.
    ! PRIVATE (codex review, cd9ddc49): mutable module state that gicov_rhs
    ! honors unconditionally -- if it were a public variable, production code
    ! could accidentally force the non-bit-for-bit general branch even on
    ! diagonal cells.  Flip it only via the public set_gicov_force_general_field
    ! setter below (used by test_gicov_hex.f90); there is no other writer.
    logical, save, private :: gicov_force_general_field = .false.

    ! Single-density current evaluation (the legacy calc_current_bloch body):
    ! PRIVATE (codex review): every caller must go through the public
    ! calc_current_bloch wrapper so the yn_sbe_gs_current_subtract readout
    ! cannot be silently bypassed by a future direct call to the core.
    private :: calc_current_bloch_core

    ! --- gicov halo-exchange plan (built once on first gicov_rhs call) ---
    ! Replaces the full-nk rho allreduce: each rank ships only the +-m_max-shell
    ! neighbour-k densities the covariant stencil actually reads (plan derived
    ! from the stencil's own chain-walk => coverage by construction; superset
    ! over all 3 axes => static under any axis_active/field direction).
    ! Deterministic ascending-k packing on both sides (see build_halo_lists).
    ! Single-solver module state (same assumption as the timers above).
    type t_gh_buf
      complex(8), allocatable :: b(:,:,:)     ! (nb, nb, cnt) per-partner message
    end type
    logical,    save :: gh_built = .false.
    integer,    save :: gh_nb = 0, gh_nk = 0
    integer,    save :: gh_nsrc = 0, gh_ndst = 0
    integer,    allocatable, save :: gh_src_rank(:), gh_src_cnt(:), gh_src_ofs(:), gh_src_k(:)
    integer,    allocatable, save :: gh_dst_rank(:), gh_dst_cnt(:), gh_dst_ofs(:), gh_dst_k(:)
    type(t_gh_buf), allocatable, save :: gh_sb(:), gh_rb(:)
    integer,    allocatable, save :: gh_reqs(:), gh_reqr(:)
    complex(8), allocatable, save :: gh_rho(:,:,:)    ! persistent rho_full (nb,nb,nk); non-needed entries stay 0 from build
    complex(8), allocatable, save :: gh_Dq(:,:,:,:)   ! persistent Dq (nb,nb,3,nk); only local slice defined per call

    ! --- integral (covariant-Houston) transport cache (sbe_lg_degen='gicov_int') ---
    ! Built ONCE by build_gicov_integral_cache from the whole-pulse trajectory:
    ! the driven reduced axis, the pulse span j_max, and -- for every LOCAL k and
    ! every integer shift n in [-j_max,j_max] -- the bounded transported band-
    ! energy matrix H~_n and the three transported velocity components v~_{n,i}.
    ! The single-step Wilson links (gs%u_transport) are replicated full-nk on
    ! every rank, so the transport chain W(kappa,kappa+n) is assembled ENTIRELY
    ! ON-RANK (no halo): gi_iknb_p/gi_iknb_m are the +/- reduced-axis neighbour
    ! index tables.  Cache is complex(8) over the LOCAL k-slice only.
    !
    ! gi_jmax is the CACHED span and gi_jmax_phys the PHYSICAL pulse span; they
    ! differ by the degree-p interpolation halo (gicov_int_jmax_cache), because
    ! the moving evaluation point x = kappa - a(t) sits BETWEEN nodes and its
    ! degree-p stencil reaches past the bracketing pair.  Everything that indexes
    ! the cache (allocation, stencil span check, memory estimate) uses gi_jmax;
    ! gi_jmax_phys is kept only to report the pulse span it was derived from.
    logical,    save :: gi_built = .false.
    integer,    save :: gi_nb = 0, gi_axis = 0, gi_jmax = 0, gi_jmax_phys = 0
    integer,    save :: gi_ikmin = 0, gi_ikmax = 0
    complex(8), allocatable, save :: gi_Ht(:,:,:,:)      ! (nb,nb,-jmax:jmax, ik_min:ik_max)
    ! shift index BEFORE the Cartesian component, so that gi_vt(:,:,:,i,ik) -- the
    ! whole cached shift range of one component of one k-block, i.e. exactly what
    ! the interpolation stencil consumes -- is CONTIGUOUS.  With the component
    ! index inner (the natural-looking order) that section is strided, and passing
    ! it to the kernel would build a copy-in/copy-out temporary on every k of
    ! every step, inside the OpenMP loop.
    complex(8), allocatable, save :: gi_vt(:,:,:,:,:)    ! (nb,nb,-jmax:jmax,3, ik_min:ik_max)
    ! co-moving GROUND-STATE occupation reference F~0_n = W diag(f0(k+n)) W^dag.
    ! _sbe_occ.data is an EXCITATION (rho - f0), so the integral path needs f0 in
    ! the SAME instantaneous frame as rho~.  It must be transported and
    ! interpolated exactly like H~: at a node H~ and F~0 share an eigenbasis and
    ! the occupation could be recovered from eigenvalue order, but BETWEEN nodes
    ! the interpolated matrices no longer commute, so the reference has to be
    ! carried explicitly.  Also the only correct metal reference: gs%occup is the
    ! real fractional, k-varying Fermi-Dirac GS occupation, which a rigid
    ! "1 below nvb" block would misreport as a spurious excitation.
    complex(8), allocatable, save :: gi_Ft(:,:,:,:)      ! (nb,nb,-jmax:jmax, ik_min:ik_max)
    integer,    allocatable, save :: gi_iknb_p(:), gi_iknb_m(:)   ! (nk) +/- driven-axis neighbour
    ! per-thread eigensolver work bank (heap allocated ONCE, outside the OMP
    ! k-loop -- frtpx discipline: no per-thread automatic arrays, no heap in loop)
    real(8),    allocatable, save :: gi_w_eps(:,:)       ! (nb, 0:nth-1)
    real(8),    allocatable, save :: gi_w_rw(:,:)        ! (3*nb, 0:nth-1)
    complex(8), allocatable, save :: gi_w_P(:,:,:), gi_w_R(:,:,:)   ! (nb,nb,0:nth-1)
    complex(8), allocatable, save :: gi_w_cw(:,:)        ! (lcwork, 0:nth-1)
    complex(8), allocatable, save :: gi_w_H(:,:,:)       ! (nb,nb,0:nth-1) interpolated H~
    complex(8), allocatable, save :: gi_w_rout(:,:,:)    ! (nb,nb,0:nth-1) step output
    complex(8), allocatable, save :: gi_w_v(:,:,:,:)     ! (nb,nb,3,0:nth-1) interpolated v~
    integer,    allocatable, save :: gi_w_blk(:,:)       ! (nb, 0:nth-1) instantaneous degenerate-block ids
    integer,    save :: gi_lcwork = 0

    ! --- persistent RK4 stage vectors of the FD-gicov propagator ---
    ! Allocated once on the LOCAL k-slice instead of four full-nk allocate +
    ! deallocate per step (the malloc churn dominated the fixed per-step cost at
    ! production nk).  Values and arithmetic are untouched.
    logical,    save :: gk_built = .false.
    integer,    save :: gk_nb = 0, gk_lo = 0, gk_hi = -1
    complex(8), allocatable, save :: gk_rho0(:,:,:), gk_rnew(:,:,:)
    complex(8), allocatable, save :: gk_k1(:,:,:), gk_k2(:,:,:), gk_k3(:,:,:), gk_k4(:,:,:)

    ! --- k-roughness diagnostic (SBE_KROUGH), see sbe_krough_diag ---
    complex(8), allocatable, save :: kr_rho0(:,:,:)      ! (nb,nb,ik_min:ik_max) rho~(t=0)
    logical,    save :: kr_built = .false.

    ! --- driven-axis transport halo (build_gicov_integral_cache; freed after) ---
    integer, save :: gih_nsrc = 0, gih_ndst = 0, gih_nh = 0
    integer, allocatable, save :: gih_src_rank(:), gih_src_cnt(:), gih_src_ofs(:), gih_src_k(:)
    integer, allocatable, save :: gih_dst_rank(:), gih_dst_cnt(:), gih_dst_ofs(:), gih_dst_k(:)
    integer, allocatable, save :: gih_map(:)             ! (nk) remote k -> compact velocity-halo slot
    complex(8), allocatable, save :: gih_p(:, :, :, :)   ! (nb,nb,3,nh) p_tm at the remote nodes
    complex(8), allocatable, save :: gih_r(:, :, :, :)   ! (nb,nb,3,nh) rvnl at the remote nodes

    ! --- one-layer forward halo plan of the k-roughness probe (sbe_krough_diag) ---
    logical, save :: krp_built = .false.
    integer, save :: krp_nsrc = 0, krp_ndst = 0, krp_nx = 0
    integer, allocatable, save :: krp_src_rank(:), krp_src_cnt(:), krp_src_ofs(:), krp_src_k(:)
    integer, allocatable, save :: krp_dst_rank(:), krp_dst_cnt(:), krp_dst_ofs(:), krp_dst_k(:)
    integer, allocatable, save :: krp_map(:)             ! (nk) k -> probe slot
    integer, allocatable, save :: krp_nbr(:, :)          ! (3, nk) forward +e_axis neighbour

contains

  !-------------------------------------------------------------------
  ! Public setter for the gicov_force_general_field TEST hook (declared
  ! private above): production code has no way to write the flag directly,
  ! only test_gicov_hex.f90 flips it, and only through this setter.
  !-------------------------------------------------------------------
  subroutine set_gicov_force_general_field(flag)
    implicit none
    logical, intent(in) :: flag
    gicov_force_general_field = flag
  end subroutine set_gicov_force_general_field

  !-------------------------------------------------------------------
  ! Print the accumulated gicov RHS timers (covariant kernel vs rho-gather
  ! comm), reduced max/avg over icomm, once on rank 0. system_clock only,
  ! no barriers were added in the hot path. Used for the scaling study.
  !-------------------------------------------------------------------
  subroutine sbe_print_timers(icomm)
    use communication, only: comm_get_max
    implicit none
    integer, intent(in) :: icomm
    integer :: irank, nproc
    real(8) :: cin(1), cout(1), din(1), dout(1), cov_sum, comm_sum, rn
    call comm_get_groupinfo(icomm, irank, nproc)
    if (sbe_tmr_nrhs == 0) return
    rn = 1d0 / dble(sbe_tmr_nrhs)
    cin(1) = sbe_tmr_cov;  din(1) = sbe_tmr_comm
    call comm_get_max(cin, cout, 1, icomm)
    call comm_get_max(din, dout, 1, icomm)
    call comm_summation(sbe_tmr_cov,  cov_sum,  icomm)
    call comm_summation(sbe_tmr_comm, comm_sum, icomm)
    if (irank == 0) then
      write(*, '(a)') "=== gicov RHS timers (covariant_grad_block vs rho halo exchange) ==="
      write(*, '(a,i0,a,i0)') "  n_rhs=", sbe_tmr_nrhs, "  nproc=", nproc
      write(*, '(a,es12.4,a,es12.4)') "  covariant [s/call]: max_rank ", cout(1)*rn, "  avg_rank ", cov_sum/dble(nproc)*rn
      write(*, '(a,es12.4,a,es12.4)') "  comm      [s/call]: max_rank ", dout(1)*rn, "  avg_rank ", comm_sum/dble(nproc)*rn
    end if
  end subroutine sbe_print_timers

  !-------------------------------------------------------------------
  ! One-time construction of the gicov rho halo-exchange plan.
  ! Every rank derives every rank's needed-set from replicated tables
  ! (split_range + the stencil chain-walk), so send/recv lists agree with
  ! no metadata exchange. Fail-closed: coverage of my own needed-set by
  ! (local + recv) is verified here; a gap is a hard error, not silent.
  ! At nproc=1 all lists are empty => the exchange issues NO point-to-point
  ! calls (required: the no-MPI dummy build aborts on isend/irecv to a real
  ! rank, and the standalone tests run exactly that build).
  !-------------------------------------------------------------------
  subroutine build_gicov_halo(sbe, gs, icomm)
    use salmon_global, only: num_kgrid
    use degenerate_block_ssbe, only: covariant_halo_needed, build_halo_lists
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    integer, intent(in) :: icomm
    integer :: irank, nproc, o, p, ik, nk, nb, jj, i, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)
    logical, allocatable :: needed_all(:,:), covered(:)
    type(t_gh_buf), allocatable :: usb(:), urb(:)
    integer, allocatable :: ureqs(:), ureqr(:)

    nk = sbe%nk;  nb = sbe%nb
    call comm_get_groupinfo(icomm, irank, nproc)
    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
    call split_range(1, nk, nproc, itbl_min, itbl_max)
    if (itbl_min(irank) /= sbe%ik_min .or. itbl_max(irank) /= sbe%ik_max) then
      write(*,*) "ERROR build_gicov_halo: split_range mismatch with solver ik range"
      error stop
    end if

    allocate(needed_all(nk, 0:nproc-1))
    do o = 0, nproc - 1
      call covariant_halo_needed(nk, gs%nbvec, gs%bvec, num_kgrid, &
                                 itbl_min(o), itbl_max(o), needed_all(:, o))
    end do
    call build_halo_lists(nk, nproc, irank, itbl_min, itbl_max, needed_all, &
                          gh_nsrc, gh_src_rank, gh_src_cnt, gh_src_ofs, gh_src_k, &
                          gh_ndst, gh_dst_rank, gh_dst_cnt, gh_dst_ofs, gh_dst_k)

    ! fail-closed coverage check: every k I need is local or in a recv list, AND
    ! it has a storage slot (gs%kmap) -- the two must agree, since gs%kmap was
    ! built from this very same covariant_halo_needed set.
    allocate(covered(nk));  covered = .false.
    covered(sbe%ik_min:sbe%ik_max) = .true.
    do p = 1, gh_nsrc
      do ik = gh_src_ofs(p) + 1, gh_src_ofs(p) + gh_src_cnt(p)
        covered(gh_src_k(ik)) = .true.
      end do
    end do
    do ik = 1, nk
      if (needed_all(ik, irank) .and. .not. covered(ik)) then
        write(*,*) "ERROR build_gicov_halo: uncovered needed k =", ik
        error stop
      end if
      if (needed_all(ik, irank) .and. gs%kmap(ik) == 0) then
        write(*,*) "ERROR build_gicov_halo: needed k has no storage slot (kmap=0), k =", ik
        error stop
      end if
    end do
    deallocate(covered)

    allocate(gh_rb(gh_nsrc), gh_sb(gh_ndst))
    do p = 1, gh_nsrc
      allocate(gh_rb(p)%b(nb, nb, gh_src_cnt(p)))
    end do
    do p = 1, gh_ndst
      allocate(gh_sb(p)%b(nb, nb, gh_dst_cnt(p)))
    end do
    allocate(gh_reqr(max(1, gh_nsrc)), gh_reqs(max(1, gh_ndst)))
    ! gh_rho lives on the SAME slots as the static fields (gs%kmap): owned k
    ! first, then the halo k this plan ships.  Replicated (kmap = identity,
    ! nks = nk) reproduces the old full-nk buffer exactly.
    allocate(gh_rho(nb, nb, max(1, gs%nks)), gh_Dq(nb, nb, 3, sbe%ik_min:sbe%ik_max))
    gh_rho = (0d0, 0d0)   ! one-time zero: never-needed entries stay 0 forever
    gh_nb = nb;  gh_nk = nk

    ! ---- ONE-TIME halo exchange of the STATIC Wilson links.  gs%u_transport is
    ! built per-k (build_block_transport is k-local), so under klocal each rank
    ! holds only its own k; the covariant stencil, however, walks U out to the
    ! +-m_max shells.  Ship those once here -- the links never change in time,
    ! so this is a setup cost, not a per-step one (contrast gh_rho, which is
    ! re-exchanged every RHS call because the density does change).
    ! Skipped when the links are replicated: every k is already local.
    if (gs%klocal) then
      allocate(urb(gh_nsrc), usb(gh_ndst))
      do p = 1, gh_nsrc
        allocate(urb(p)%b(nb, nb, 3 * gh_src_cnt(p)))
      end do
      do p = 1, gh_ndst
        allocate(usb(p)%b(nb, nb, 3 * gh_dst_cnt(p)))
      end do
      allocate(ureqr(max(1, gh_nsrc)), ureqs(max(1, gh_ndst)))
      do p = 1, gh_nsrc
        ureqr(p) = comm_irecv(urb(p)%b, gh_src_rank(p), 1, icomm)
      end do
      do p = 1, gh_ndst
        do jj = 1, gh_dst_cnt(p)
          do i = 1, 3
            usb(p)%b(:, :, 3 * (jj - 1) + i) = &
              & gs%u_transport(:, :, i, gs%kmap(gh_dst_k(gh_dst_ofs(p) + jj)))
          end do
        end do
        ureqs(p) = comm_isend(usb(p)%b, gh_dst_rank(p), 1, icomm)
      end do
      if (gh_nsrc > 0) then
        call comm_wait_all(ureqr(1:gh_nsrc))
        do p = 1, gh_nsrc
          do jj = 1, gh_src_cnt(p)
            do i = 1, 3
              gs%u_transport(:, :, i, gs%kmap(gh_src_k(gh_src_ofs(p) + jj))) = &
                & urb(p)%b(:, :, 3 * (jj - 1) + i)
            end do
          end do
        end do
      end if
      if (gh_ndst > 0) call comm_wait_all(ureqs(1:gh_ndst))
      do p = 1, gh_nsrc
        deallocate(urb(p)%b)
      end do
      do p = 1, gh_ndst
        deallocate(usb(p)%b)
      end do
      deallocate(urb, usb, ureqr, ureqs)
      ierr = 0
    end if

    gh_built = .true.
  end subroutine build_gicov_halo



subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication
    use salmon_global, only: yn_sbe_gs_current_subtract, yn_vnl_correction, yn_sbe_vnl_exact, &
                              sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_t2_gate_width, &
                              sbe_lg_degen_floor, sbe_lg_degen
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, jb, nk_proc, irank, nproc, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    ! Fail-closed: with rank-local static fields the gs k-slice IS the solver's
    ! k-slice.  Both come from split_range on the same nk, but they are split on
    ! the communicator each was HANDED -- a caller that builds gs on one
    ! communicator and the solver on another (multiscale does exactly that, which
    ! is why klocal is off for maxwell_sbe) would otherwise read k it does not own.
    if (gs%klocal) then
        if (gs%ik_lo /= sbe%ik_min .or. gs%ik_hi /= sbe%ik_max) then
            write(*, '(a)') "ERROR(init_sbe_bloch_solver): the rank-local gs k-slice does not " // &
                & "match the solver k-slice (init_sbe_gs_info and init_sbe_bloch_solver were " // &
                & "given different communicators)."
            error stop
        end if
    end if

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))

    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do

    ! T2 Delta-omega dephasing gate weight: precompute ONCE from the STATIC
    ! gs%delta_omega (design: gw/gw_design/plans/2026-07-08-t2-gate-shape.md).
    ! gicov_rhs Hadamard-multiplies this into the -rho/t_2 dephasing term every
    ! step instead of recomputing it. Diagonal forced to 0 (no population
    ! dephasing) regardless of what t2_gate_weight would return for delta_omega=0.
    ! GUARDED to sbe_lg_degen=='gicov': only gicov_rhs's 'gauss' branch reads
    ! this array, and gauss is a gicov-only construction. The 'step' branch
    ! (default, and every non-gicov length-gauge path) never touches
    ! sbe%t2_gate_w, so allocating the nb^2 * nk_local array there would waste
    ! memory for no consumer.
    if (uses_fd_gicov(sbe_lg_degen)) then
        allocate(sbe%t2_gate_w(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
        do ik = sbe%ik_min, sbe%ik_max
            do jb = 1, sbe%nb
                do ib = 1, sbe%nb
                    if (ib == jb) then
                        sbe%t2_gate_w(ib, jb, ik) = 0d0
                    else
                        sbe%t2_gate_w(ib, jb, ik) = t2_gate_weight(gs%delta_omega(ib, jb, ik), &
                            & sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_t2_gate_width, sbe_lg_degen_floor)
                    end if
                end do
            end do
        end do
    end if

    sbe%flag_vnl_correction = .false.

    ! VG completion: wire the exact-nonlocal flag from the SAME global the
    ! checker validated, and fail closed on the two structural preconditions:
    ! the gs-side stencil must have been read (read_vnl_kappa_data), and its
    ! k-slice must coincide with this solver's split (both use split_range on
    ! the same icomm, so a mismatch is a wiring bug, not a configuration).
    sbe%flag_vnl_exact = (yn_sbe_vnl_exact == 'y')
    if (sbe%flag_vnl_exact) then
        if (.not. gs%vnl_exact) then
            write(*, '(a)') "ERROR(init_sbe_bloch_solver): yn_sbe_vnl_exact='y' but the " // &
                & "kappa-stencil was not loaded into gs (read_vnl_kappa_data wiring bug)."
            error stop
        end if
        if (gs%vnl_ik_min /= sbe%ik_min .or. gs%vnl_ik_max /= sbe%ik_max) then
            write(*, '(a)') "ERROR(init_sbe_bloch_solver): vnl_kappa k-slice differs from " // &
                & "the solver k-slice (split_range mismatch)."
            error stop
        end if
    end if

    ! Window f-sum-rule deficiency tensor (yn_sbe_gs_current_subtract='y'):
    ! built once here, applied per step by calc_current_bloch.  Reset FIRST
    ! unconditionally (codex review): a REUSED solver object must not carry a
    ! stale flag-on tensor from an earlier init into a flag-off -> flag-on
    ! toggle, where it would pass the fail-closed readout check.  The vnl choice
    ! is taken from the SAME global every production caller assigns to
    ! flag_vnl_correction right after this init (realtime_ssbe / multiscale_ssbe
    ! both set it to yn_vnl_correction=='y'); a caller that diverges from that
    ! convention must rebuild via build_fsum_deficiency_tensor explicitly.
    sbe%fsum_D(:, :) = 0d0
    sbe%fsum_D_built = .false.
    if (yn_sbe_gs_current_subtract == 'y') then
        ! pi convention: p + rvnl whenever the nonlocal-velocity physics is ON
        ! (first-order rvnl OR the all-order exact mode -- at A=0, where the
        ! linear-response tensor lives, both use the same velocity p + W_0;
        ! in exact mode gs%rvnl_tm_matrix has been overwritten by the
        ! binary-precision W_0, so the SAME array serves both).  In exact mode
        ! the builder additionally adds the nonlocal sum-rule term T (see its
        ! header), reduced over icomm because gs%vnl_W is k-sliced.
        call build_fsum_deficiency_tensor(sbe, gs, &
            & yn_vnl_correction == 'y' .or. yn_sbe_vnl_exact == 'y', icomm)
    end if
end subroutine


!===================================================================
! VG completion: locate A(t) on the 1D kappa-stencil and build the centered
! 4-point Lagrange weights.
!   a  = e_dir . A   (the trajectory is validated to be collinear with e_dir)
!   x  = a/h, base interval [q0, q0+1), q0 = floor(x), nodes q0-1 .. q0+2,
!   wl = cubic Lagrange weights at t = x - q0 in [0,1).
! Exact at the nodes (t=0 => wl = (0,1,0,0)), so DeltaV(A=0) == 0 identically.
! Fail-closed: a transverse A component or a stencil overrun is an error stop
! (both are excluded up front by sbe_vnl_validate_trajectory; the checks here
! are the last line of defense for future callers with dynamic A).
!===================================================================
subroutine sbe_vnl_interp_weights(gs, Ac, wl, q0)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in)  :: Ac(1:3)
    real(8), intent(out) :: wl(1:4)
    integer, intent(out) :: q0
    real(8) :: a, t, trans(1:3)

    a = dot_product(gs%vnl_dir(1:3), Ac(1:3))
    trans(1:3) = Ac(1:3) - a * gs%vnl_dir(1:3)
    if (norm2(trans) > 1d-8 * max(1d0, abs(a))) then
        write(*, '(a,es12.4)') "ERROR(sbe_vnl_interp_weights): A(t) has a transverse " // &
            & "component w.r.t. the kappa-stencil axis: |A_perp| = ", norm2(trans)
        error stop
    end if
    q0 = floor(a / gs%vnl_h)
    if (q0 - 1 < -gs%vnl_ns .or. q0 + 2 > gs%vnl_ns) then
        write(*, '(a,es12.4,a)') "ERROR(sbe_vnl_interp_weights): |e.A| = ", abs(a), &
            & " exceeds the kappa-stencil range (increase sbe_vnl_kappa_amax and re-export)."
        error stop
    end if
    t = a / gs%vnl_h - dble(q0)
    wl(1) = -t * (t - 1d0) * (t - 2d0) / 6d0
    wl(2) =  (t * t - 1d0) * (t - 2d0) / 2d0
    wl(3) = -t * (t + 1d0) * (t - 2d0) / 2d0
    wl(4) =  t * (t * t - 1d0) / 6d0
end subroutine sbe_vnl_interp_weights


!===================================================================
! VG completion: pre-validate the WHOLE precomputed A(t) trajectory against
! the kappa-stencil BEFORE propagation starts (realtime_ssbe builds
! Ac_ext_t(-1..nt+1) up front), so an under-sized stencil or a mis-aligned
! polarization axis fails at t=0 with a clear message instead of mid-run.
!===================================================================
subroutine sbe_vnl_validate_trajectory(gs, Ac_ext_t, n_lo, n_hi, irank)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: n_lo, n_hi, irank
    real(8), intent(in) :: Ac_ext_t(1:3, n_lo:n_hi)
    integer :: it, q0
    real(8) :: a, amax_traj, trans(1:3), tmax

    if (.not. gs%vnl_exact) then
        write(*, '(a)') "ERROR(sbe_vnl_validate_trajectory): kappa-stencil not loaded."
        error stop
    end if
    amax_traj = 0d0
    tmax = 0d0
    do it = n_lo, n_hi
        a = dot_product(gs%vnl_dir(1:3), Ac_ext_t(1:3, it))
        trans(1:3) = Ac_ext_t(1:3, it) - a * gs%vnl_dir(1:3)
        tmax = max(tmax, norm2(trans) / max(1d0, abs(a)))
        if (norm2(trans) > 1d-8 * max(1d0, abs(a))) then
            if (irank == 0) write(*, '(a,i0,a,es12.4)') &
                & "ERROR(sbe_vnl_validate_trajectory): A(t) at step ", it, &
                & " has a transverse component w.r.t. the stencil axis: |A_perp| = ", norm2(trans)
            error stop
        end if
        q0 = floor(a / gs%vnl_h)
        if (q0 - 1 < -gs%vnl_ns .or. q0 + 2 > gs%vnl_ns) then
            if (irank == 0) then
                write(*, '(a,i0,a,es12.4)') "ERROR(sbe_vnl_validate_trajectory): |e.A| at step ", &
                    & it, " = ", abs(a)
                write(*, '(a,es12.4,a)') "  exceeds the usable stencil range ", &
                    & (gs%vnl_ns - 2) * gs%vnl_h, &
                    & " (increase sbe_vnl_kappa_amax and re-export the stencil)."
            end if
            error stop
        end if
        amax_traj = max(amax_traj, abs(a))
    end do
    if (irank == 0) then
        write(*, '(a,es12.4,a,es12.4,a,es10.2)') "# vnl_kappa trajectory: max|e.A| = ", &
            & amax_traj, " of usable range ", (gs%vnl_ns - 2) * gs%vnl_h, &
            & ", max transverse ratio = ", tmax
    end if
end subroutine sbe_vnl_validate_trajectory


!===================================================================
! Velocity-gauge current readout.  Default (yn_sbe_gs_current_subtract='n'):
! the legacy evaluation J = Tr[rho(t) v(k+A)] (calc_current_bloch_core,
! byte-identical code path).  With 'y', the WINDOW F-SUM-RULE DEFICIENCY
! current is subtracted:
!
!   J(t) -= D * A(t)          (3x3 tensor D = sbe%fsum_D, built once at init)
!
! WHY (v2; replaces the v1 frozen-GS subtraction of commit cda60ee4):
! Truncating the basis at nstate_sbe breaks the f-sum rule: the propagated rho
! can only build the paramagnetic counter-term to the diamagnetic A(t)*Ne
! readout term from occupied->WINDOW transitions, so a pseudo-linear current
! D*A(t) survives, where D is exactly the part of the f-sum rule the window
! CANNOT capture (a filled-band insulator in a complete basis carries no
! A-linear current).  The v1 remedy subtracted the whole frozen-GS current
! Tr[rho0 v(k+A)] = A*Ne/V (at norder=0 on a TR-symmetric grid the
! paramagnetic trace Tr[rho0 p] vanishes identically), i.e. it removed the
! FULL diamagnetic term -- including the part that the retained window
! response legitimately cancels dynamically.  Measured on Si@nb=32 VG
! (JID49460686): J_v1 - J_legacy = -(Ne/V)*A = -0.0296284*A to 6 digits,
! and peak |Jz| moved AWAY from the LG oracle (2.0517e-4 -> 9.4142e-4).
! v2 subtracts only the response the window is MISSING (see
! build_fsum_deficiency_tensor for the derivation and normalization):
! complete window => D -> 0 (nothing subtracted); empty conduction window =>
! D -> (Ne/V)*delta_ij (the v1 limit, then correct).  D*A vanishes for A=0,
! so post-pulse observables are unchanged, and D is A-independent, so the
! nonlinear response is untouched.  In EXACT mode (yn_sbe_vnl_exact='y') D
! additionally carries the nonlocal sum-rule term T of the A-dependent
! readout velocity W(k+A) -- see build_fsum_deficiency_tensor; without it
! the complete-window limit of D is -T/V instead of 0.
!
! Scope: derived for the norder_correction=0 readout (the production default).
! At norder_correction>=1 the core's own A-linear velocity-operator correction
! (which at rho ~ rho0 supplies a bare-p, 1d-3-floored variant of the
! captured-window response -S*A) overlaps with the response D is built
! against; the checker WARNs on that combination.
!===================================================================
subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    use salmon_global, only: yn_sbe_gs_current_subtract
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm

    call calc_current_bloch_core(sbe, gs, Ac, jmat, icomm)

    if (yn_sbe_gs_current_subtract == 'y') then
        ! Fail-closed: the deficiency tensor is built by init_sbe_bloch_solver
        ! (or explicitly by build_fsum_deficiency_tensor); a solver that
        ! reaches the readout without it is a wiring bug, not a soft default.
        if (.not. sbe%fsum_D_built) then
            write(*, '(a)') "ERROR(calc_current_bloch): yn_sbe_gs_current_subtract='y' " // &
                & "but the f-sum deficiency tensor has not been built for this solver " // &
                & "(init_sbe_bloch_solver / build_fsum_deficiency_tensor)."
            error stop
        end if
        jmat(1:3) = jmat(1:3) - matmul(sbe%fsum_D(1:3, 1:3), Ac(1:3))
    end if

    return
end subroutine calc_current_bloch


!===================================================================
! Build the window f-sum-rule deficiency tensor D (sbe%fsum_D) from the GS
! data -- called once per solver (setup time), zero per-step cost.
!
! MODEL: the VG-SBE propagates  i drho/dt = [H, rho],  H_k(t) = diag(eps_k)
! + A(t).pi_k  with pi = p_tm (+ rvnl_tm iff flag_vnl), and reads out
!   J(t) = ( <Tr[rho pi]>_k + A(t) <Tr[rho]>_k ) / V,
! where <X>_k = sum_k w_k X_k / sum_k w_k (calc_current_bloch_core, norder=0).
!
! DERIVATION: first-order response of rho to dH = A.pi (adiabatic limit of
! the VG Kubo kernel; equal to static perturbation theory of the shifted
! Fermi sea):  rho1_ab = (f_b - f_a) (A.pi)_ab / (eps_b - eps_a),  a /= b,
! with f = gs%occup.  The A-linear part of the readout is then
!
!   J_lin,i = D_ij A_j,
!   D_ij = ( Ne delta_ij
!            - < sum_{v: f_v>0} f_v sum_{m: f_m=0}
!                2 Re[ pi^i_vm pi^j_mv ] / (eps_m - eps_v) >_k ) / V
!
! NORMALIZATION / SPIN CONVENTION: f = gs%occup as set by init_sbe_gs_info
! (spinless: 2 per filled band; spinor/noncollinear SO-SBE: 1 per filled
! band), and
! Ne = <sum_b f_b>_k = Tr[rho0] -- the same trace the diamagnetic readout
! term A*calc_trace() carries (the SBE propagation conserves the trace).
! With the occupation factor f_v kept EXPLICIT under the sum, the
! Thomas-Reiche-Kuhn sum rule per band,
!   sum_{m/=v} 2 Re[ p^i_vm p^j_mv ] / (eps_m - eps_v) = delta_ij
! (one electron; spin enters only through f_v), gives in a COMPLETE basis
!   < sum_v f_v delta_ij >_k = Ne delta_ij   =>   D -> 0 identically --
! BOTH occupation conventions are served by the same expression.
!
! EXACT MODE (yn_sbe_vnl_exact='y') -- nonlocal sum-rule term T: with a
! NONLOCAL potential the TRK sum rule generalizes to the effective-mass form
!   sum_m 2 Re[ pi^i_vm pi^j_mv ] / (eps_m - eps_v)
!     = delta_ij + <v| d2 V_nl(kappa)/dkappa_i dkappa_j |v> - d2 eps_v/dk_i dk_j
! (kappa-derivatives at frozen |u_vk>; the band-mass term integrates to zero
! over the BZ for filled bands), and the exact readout velocity p + W(k+A)
! carries the matching A-linear GROUND-STATE term Tr[rho0 dW_i/dkappa_j] A_j
! that the legacy frozen-velocity readout does not have.  The window
! deficiency is therefore
!   D_ij = ( Ne delta_ij + T_ij - S_ij ) / V,
!   T_ij = < sum_v f_v Re <v| dW_i/dkappa_j |v> >_k,
! and OMITTING T over-subtracts by -T/V even in a COMPLETE basis (measured
! on Si 8^3 @ nb=32, z-stencil: T_zz/V = -1.82e-3 a.u. = 46% of peak |J|).
! The 1D kappa-stencil provides the derivative along its axis only, i.e. the
! e_dir COLUMN of T -- exactly the part ever contracted with A(t), because
! exact mode pre-validates A(t) // vnl_dir (sbe_vnl_validate_trajectory).
! dW/dkappa is evaluated by the O(h^4) five-point first-derivative stencil at
! the s=0 node (the reader enforces ns >= 4, so s = +-1, +-2 always exist).
! gs%vnl_W is stored k-sliced (unlike the replicated gs%p/rvnl), so the T
! part is a sliced partial sum reduced over icomm; legacy modes are
! communication-free and bit-identical to the pre-T behavior.
! Occupied<->occupied pairs are excluded: they carry (f_v - f_m) = 0 in the
! response (ssbe occupations are exactly 0 or focc by construction), and in the
! TRK rearrangement they cancel pairwise (numerator symmetric, denominator
! antisymmetric under v<->m) -- restricting m to the unoccupied window also
! keeps degenerate occupied pairs (eps_m - eps_v -> 0) out of the sum by
! construction.  Occupied x unoccupied near-degeneracies (Fermi-crossing,
! metals) are guarded by the theta_off floor: the adiabatic response of such
! a pair is not perturbative, so it is dropped rather than amplified.
!
! On a time-reversal-symmetric k-grid D is a symmetric real tensor
! (pi(-k) = -conjg(pi(k)) makes each pair's contribution real and i<->j
! symmetric after the k-sum); the full 3x3 is computed without assuming it.
!
! The gs%* arrays are replicated on every rank (comm_bcast in
! init_sbe_gs_info), so the FULL k-sum of the S/Ne part is evaluated
! redundantly on each rank: bit-identical, no communication.  The exact-mode
! T part is the one exception (k-sliced gs%vnl_W, comm_summation reduction --
! see the EXACT MODE paragraph above); legacy modes never reach it.
!===================================================================
subroutine build_fsum_deficiency_tensor(sbe, gs, flag_vnl, icomm)
    use degenerate_block_ssbe, only: theta_off
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    logical, intent(in) :: flag_vnl
    integer, intent(in) :: icomm
    integer :: ik, iv, im, i, j
    real(8) :: s(1:3, 1:3), fv, de, w, wsum, ne_w
    real(8) :: t_col(1:3), t_col_l(1:3)
    complex(8) :: pvm(1:3), pmv(1:3)
    real(8), parameter :: occ_eps = 1d-12

    ! Fail-closed (metal): the occ_eps binary threshold below ("occupied<->
    ! occupied: no net response") assumes occupations are exactly 0 or focc by
    ! construction (see this routine's own header derivation).  For a metal
    ! with genuinely fractional (finite-temperature) occupation this silently
    ! DROPS any valence<->conduction pair where BOTH bands are fractional
    ! (nonzero response weight f_v - f_m gets skipped rather than summed).
    ! Not yet generalized to (f_v-f_m)-weighted summation (follow-up) -- fail
    ! loud rather than silently under-count until then.
    if (gs%is_metal) then
        write(*, '(a)') "ERROR(build_fsum_deficiency_tensor): metal-like occupation detected; " // &
            & "yn_sbe_gs_current_subtract='y' is not yet generalized to fractional/k-varying " // &
            & "occupation (see this routine's occ_eps threshold) -- set " // &
            & "yn_sbe_gs_current_subtract='n' for a metal."
        error stop
    end if

    s(:, :) = 0d0
    ne_w = 0d0
    do ik = 1, gs%nk
        w = gs%kweight(ik)
        do iv = 1, sbe%nb
            fv = gs%occup(iv, ik)
            ne_w = ne_w + w * fv
            if (fv <= occ_eps) cycle
            do im = 1, sbe%nb
                if (gs%occup(im, ik) > occ_eps) cycle   ! occupied<->occupied: no net response
                de = gs%eigen(im, ik) - gs%eigen(iv, ik)
                if (abs(de) < theta_off) cycle          ! Fermi-crossing near-degeneracy guard
                pvm(1:3) = gs%p_tm_matrix(iv, im, 1:3, ik)
                pmv(1:3) = gs%p_tm_matrix(im, iv, 1:3, ik)
                if (flag_vnl) then
                    ! same velocity operator as the propagation + readout
                    pvm(1:3) = pvm(1:3) + gs%rvnl_tm_matrix(iv, im, 1:3, ik)
                    pmv(1:3) = pmv(1:3) + gs%rvnl_tm_matrix(im, iv, 1:3, ik)
                end if
                do j = 1, 3
                    do i = 1, 3
                        s(i, j) = s(i, j) + w * fv * 2d0 * dble(pvm(i) * pmv(j)) / de
                    end do
                end do
            end do
        end do
    end do

    wsum = sum(gs%kweight(:))
    sbe%fsum_D(:, :) = -s(:, :) / wsum
    do i = 1, 3
        sbe%fsum_D(i, i) = sbe%fsum_D(i, i) + ne_w / wsum
    end do
    sbe%fsum_D(:, :) = sbe%fsum_D(:, :) / gs%volume

    ! Exact mode: add the nonlocal sum-rule term T (e_dir column; see the
    ! header).  D(i, :) += (T_i / V) * e_j with
    !   T_i = < sum_v f_v Re[ (dW_i/da)_vv at s=0 ] >_k,   a = e_dir.kappa,
    ! from the O(h^4) five-point derivative of the exact-mode W stencil.
    ! Sliced k-sum (gs%vnl_W holds this rank's slice only) + allreduce.
    if (sbe%flag_vnl_exact) then
        t_col_l(1:3) = 0d0
        do ik = gs%vnl_ik_min, gs%vnl_ik_max
            w = gs%kweight(ik)
            do iv = 1, sbe%nb
                fv = gs%occup(iv, ik)
                if (fv <= occ_eps) cycle
                do i = 1, 3
                    t_col_l(i) = t_col_l(i) + w * fv * dble( &
                        &  ( 8d0 * (gs%vnl_W(iv, iv, i,  1, ik) - gs%vnl_W(iv, iv, i, -1, ik))  &
                        &        - (gs%vnl_W(iv, iv, i,  2, ik) - gs%vnl_W(iv, iv, i, -2, ik)) ) &
                        & / (12d0 * gs%vnl_h) )
                end do
            end do
        end do
        call comm_summation(t_col_l, t_col, 3, icomm)
        do i = 1, 3
            sbe%fsum_D(i, 1:3) = sbe%fsum_D(i, 1:3) &
                & + (t_col(i) / (wsum * gs%volume)) * gs%vnl_dir(1:3)
        end do
    end if

    sbe%fsum_D_built = .true.
end subroutine build_fsum_deficiency_tensor


! Single-density current evaluation (the pre-subtraction legacy body, moved
! verbatim from calc_current_bloch; see the wrapper above).
subroutine calc_current_bloch_core(sbe, gs, Ac, jmat, icomm)
    use salmon_global, only: norder_correction
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, kb, nb, ierr
    complex(8) :: tmp1(1:3), tmp(1:3)
    complex(8) :: sum1(1:3)
    complex(8) :: pnn(1:3),pni(1:3),pin(1:3)
    complex(8) :: pnj(1:3),pjn(1:3),pnk(1:3),pkn(1:3)
    complex(8) :: pij(1:3),pji(1:3),pik(1:3),pki(1:3),pjk(1:3),pkj(1:3)
    complex(8) :: pnn_Ac,pni_Ac,pin_Ac
    complex(8) :: pnj_Ac,pjn_Ac,pnk_Ac,pkn_Ac
    complex(8) :: pij_Ac,pji_Ac,pik_Ac,pki_Ac,pjk_Ac,pkj_Ac
    ! VG completion: interpolated all-order nonlocal velocity W_i(k+A)
    logical :: do_exact
    real(8) :: wl(1:4)
    integer :: q0

    do_exact = sbe%flag_vnl_exact
    if (do_exact) call sbe_vnl_interp_weights(gs, Ac, wl, q0)

    tmp1(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                        & gs%p_tm_matrix(ib, jb, idir, ik) &
                    & )
                    if (sbe%flag_vnl_correction) then
                        tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                            & gs%rvnl_tm_matrix(ib, jb, idir, ik) &
                        & )
                    endif
                    if (do_exact) then
                        ! v^nl_i(k+A): all 3 Cartesian components are exact at
                        ! every stencil node (W_i = dV/dkappa_i evaluated at
                        ! kappa_s in the export), interpolated along the axis.
                        tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                            & wl(1) * gs%vnl_W(ib, jb, idir, q0-1, ik) &
                            & + wl(2) * gs%vnl_W(ib, jb, idir, q0,   ik) &
                            & + wl(3) * gs%vnl_W(ib, jb, idir, q0+1, ik) &
                            & + wl(4) * gs%vnl_W(ib, jb, idir, q0+2, ik) &
                        & )
                    endif
                end do
            end do
        end do
    end do
    !$omp end parallel do

    if(norder_correction>=1)then
        ! Fail-closed (metal): this term's occupied-band loop is a HARD
        ! 1..gs%nvb cutoff and never reads gs%occup inside the loop (unlike
        ! the norder=0 term above) -- for a metal the true occupied-band set
        ! varies by k (Fermi surface), so the same nvb cutoff at every k
        ! silently mixes in wrong bands. Not yet generalized to an
        ! occupation-weighted sum over the full window (follow-up); the
        ! production default is norder_correction=0, so this costs nothing
        ! for the metal-safe path.
        if (gs%is_metal) then
            write(*, '(a)') "ERROR(calc_current_bloch_core): metal-like occupation detected; " // &
                & "norder_correction>=1 is not yet generalized to a k-varying occupied-band " // &
                & "window (hard 1..nvb cutoff below) -- set norder_correction=0 for a metal."
            error stop
        end if
        do ik = sbe%ik_min, sbe%ik_max
            do idir = 1, 3
                ! occupied bands 1..gs%nvb (= gs%ne/2 spinless, = gs%ne
                ! spinor); the occupation itself is carried by rho(nb,nb) and
                ! the 2.d0 below is the c.c.-pair factor, NOT spin degeneracy.
                do nb = 1, gs%nvb
                    do ib = 1, sbe%nb
                        pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
                        pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
                        if(nb /= ib) then
                            if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                                tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                    2.d0 * sbe%rho(nb, nb, ik) * &
                                    dble(pni_Ac*pin(idir)) / gs%delta_omega(ib, nb, ik) 
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end if

    if(norder_correction>=2)then
        do ik = sbe%ik_min, sbe%ik_max
            do idir = 1, 3
                do nb = 1, sbe%nb
                    pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
                    pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                             gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                             gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
                    do ib = 1, sbe%nb
                        pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
                        pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
                        pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
                        do jb = 1, sbe%nb
                            pjn(idir) = gs%p_tm_matrix(jb, nb, idir, ik)
                            pij(idir) = gs%p_tm_matrix(ib, jb, idir, ik)
                            pij_Ac = gs%p_tm_matrix(ib, jb, 1, ik) * Ac(1) + &
                                     gs%p_tm_matrix(ib, jb, 2, ik) * Ac(2) + &
                                     gs%p_tm_matrix(ib, jb, 3, ik) * Ac(3)
                            pjn_Ac = gs%p_tm_matrix(jb, nb, 1, ik) * Ac(1) + &
                                     gs%p_tm_matrix(jb, nb, 2, ik) * Ac(2) + &
                                     gs%p_tm_matrix(jb, nb, 3, ik) * Ac(3)
                            if(nb /= ib .and. nb /= jb) then
                                if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                                   abs(gs%delta_omega(jb, nb, ik))> 1.d-3)then
                                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * &
                                        sbe%rho(nb, nb, ik) / &
                                        gs%delta_omega(ib, nb, ik) / &
                                        gs%delta_omega(jb, nb, ik) * &
                                        (pjn(idir) * pij_Ac * pni_Ac &
                                       + pjn_Ac * (pij(idir) * pni_Ac  + pni(idir) * pij_Ac))
                                end if
                            end if
                        end do
                        if(nb /= ib) then
                            if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                                tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                    sbe%rho(nb, nb, ik) * &
                                    (2.d0 * pnn_Ac * dble(pni(idir) * pin_Ac) + &
                                    pnn(idir)*abs(pni_Ac)**2) / &
                                    gs%delta_omega(ib, nb, ik)**2
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end if

    if(norder_correction>=3)then
      do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
          do nb = 1, sbe%nb
            pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
            pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                     gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                     gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
            do ib = 1, sbe%nb
              pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
              pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
              pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
              pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
              do jb = 1, sbe%nb
                pnj(idir) = gs%p_tm_matrix(nb, jb, idir, ik)
                pjn(idir) = gs%p_tm_matrix(jb, nb, idir, ik)
                pij(idir) = gs%p_tm_matrix(ib, jb, idir, ik)
                pji(idir) = gs%p_tm_matrix(jb, ib, idir, ik)
                pnj_Ac = gs%p_tm_matrix(nb, jb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(nb, jb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(nb, jb, 3, ik) * Ac(3)
                pjn_Ac = gs%p_tm_matrix(jb, nb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(jb, nb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(jb, nb, 3, ik) * Ac(3)
                pij_Ac = gs%p_tm_matrix(ib, jb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(ib, jb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(ib, jb, 3, ik) * Ac(3)
                pji_Ac = gs%p_tm_matrix(jb, ib, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(jb, ib, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(jb, ib, 3, ik) * Ac(3)
                do kb = 1, sbe%nb
                  pnk(idir) = gs%p_tm_matrix(nb, kb, idir, ik)
                  pkn(idir) = gs%p_tm_matrix(kb, nb, idir, ik)
                  pik(idir) = gs%p_tm_matrix(ib, kb, idir, ik)
                  pki(idir) = gs%p_tm_matrix(kb, ib, idir, ik)
                  pjk(idir) = gs%p_tm_matrix(jb, kb, idir, ik)
                  pkj(idir) = gs%p_tm_matrix(kb, jb, idir, ik)
                  pnk_Ac = gs%p_tm_matrix(nb, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(nb, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(nb, kb, 3, ik) * Ac(3)
                  pkn_Ac = gs%p_tm_matrix(kb, nb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, nb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, nb, 3, ik) * Ac(3)
                  pik_Ac = gs%p_tm_matrix(ib, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(ib, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(ib, kb, 3, ik) * Ac(3)
                  pki_Ac = gs%p_tm_matrix(kb, ib, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, ib, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, ib, 3, ik) * Ac(3)
                  pjk_Ac = gs%p_tm_matrix(jb, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(jb, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(jb, kb, 3, ik) * Ac(3)
                  pkj_Ac = gs%p_tm_matrix(kb, jb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, jb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, jb, 3, ik) * Ac(3)
                  if(nb /= ib .and. nb /= jb .and. nb /= kb) then
                    if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                       abs(gs%delta_omega(jb, nb, ik))> 1.d-3 .and. &
                       abs(gs%delta_omega(kb, nb, ik))> 1.d-3)then
                      tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                   sbe%rho(nb, nb, ik) / &
                                   gs%delta_omega(ib, nb, ik) / &
                                   gs%delta_omega(jb, nb, ik) / &
                                   gs%delta_omega(kb, nb, ik) * &
                                   dble(pji_Ac * ( pnk_Ac * ( pin(idir) * pkj_Ac + &
                                                              pkj(idir) * pin_Ac ) + &
                                                   pnk(idir) * pin_Ac * pkj_Ac ) + &
                                        pji(idir) * pin_Ac * pkj_Ac * pnk_Ac)
                    end if
                  end if
                end do
                if(nb /= ib .and. nb /= jb) then
                  if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                     abs(gs%delta_omega(jb, nb, ik))> 1.d-3)then
                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * &
                                 sbe%rho(nb, nb, ik) * &
                                 (gs%delta_omega(ib, nb, ik) + gs%delta_omega(jb, nb, ik)) / &
                                 2.d0 / &
                                 gs%delta_omega(ib, nb, ik)**2 / &
                                 gs%delta_omega(jb, nb, ik)**2 * &
                                 (pji(idir) * pnn_Ac * pin_Ac * pnj_Ac + &
                                  pjn(idir) * pni_Ac * (pin_Ac * pnj_Ac + &
                                                        pnn_Ac * pij_Ac) + &
                                  pji_Ac * (pnn_Ac * (pin(idir) * pnj_Ac + &
                                                      pnj(idir) * pin_Ac) + &
                                            pnn(idir) * pin_Ac * pnj_Ac) + &
                                  pjn_Ac * (pni(idir) * ( pin_Ac * pnj_Ac + &
                                                          pnn_Ac * pij_Ac ) + &
                                            pni_Ac * (pij(idir) * pnn_Ac + &
                                                      pin(idir) * pnj_Ac + &
                                                      pnj(idir) * pin_Ac + &
                                                      pnn(idir) * pij_Ac)))
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
      do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
          do nb = 1, sbe%nb
            pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
            pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                     gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                     gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
            sum1 = 0.d0
            do ib = 1, sbe%nb
              pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
              pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
              pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
              pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
              if(nb /= ib) then
                if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                  sum1(idir) = sum1(idir) + &
                          (pni(idir) * pnn_Ac * pin_Ac + &
                           pni_Ac * (pin(idir) * pnn_Ac + &
                                    2.d0 * pnn(idir) * pin_Ac)) / &
                           gs%delta_omega(ib, nb, ik)**3
                end if
              end if
            end do
            tmp1(idir) = tmp1(idir) - gs%kweight(ik) * sbe%rho(nb, nb, ik) * pnn_Ac * sum1(idir)
          end do
        end do
      end do
    end if

    call comm_summation(tmp1, tmp, 3, icomm)

    jmat(:) = (real(tmp(1:3)) / sum(gs%kweight(:)) &
        & + Ac * calc_trace(sbe, gs, sbe%nb, icomm)) / gs%volume

    return
end subroutine calc_current_bloch_core


subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    use salmon_global, only: yn_sbe_gw_collision, sbe_deph_mode
    use sbe_collision_gw, only: add_collision_vg
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: nb, nk, ik

    complex(8) :: hrho1_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho2_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho3_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho4_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3)
    ! VG completion: all-order nonlocal difference DeltaV = V(k+A) - V(k),
    ! interpolated once per (ik, step) -- the SAME A serves all 4 Taylor
    ! stages (existing convention: Ac is constant within the step).
    complex(8) :: dv_k(1:sbe%nb, 1:sbe%nb)
    logical :: do_exact
    real(8) :: wl(1:4)
    integer :: q0

    nb = sbe%nb
    nk = sbe%nk

    do_exact = sbe%flag_vnl_exact
    if (do_exact) call sbe_vnl_interp_weights(gs, Ac, wl, q0)

    !$omp parallel do default(shared) private(ik, p_rvnl_k, dv_k, hrho1_k, hrho2_k, hrho3_k, hrho4_k)
    do ik = sbe%ik_min, sbe%ik_max
        p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) = gs%p_tm_matrix(1:sbe%nb, 1:sbe%nb, 1:3, ik)
        if (sbe%flag_vnl_correction) then
            p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) =  p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) &
                & + gs%rvnl_tm_matrix(1:sbe%nb, 1:sbe%nb, 1:3, ik)
        end if
        if (do_exact) then
            ! Hermitian by construction (real weights x Hermitized stencil);
            ! exactly zero at A=0 (node-exact Lagrange weights).
            dv_k(:, :) = wl(1) * gs%vnl_V(:, :, q0-1, ik) &
                     & + wl(2) * gs%vnl_V(:, :, q0,   ik) &
                     & + wl(3) * gs%vnl_V(:, :, q0+1, ik) &
                     & + wl(4) * gs%vnl_V(:, :, q0+2, ik) &
                     & - gs%vnl_V(:, :, 0, ik)
        end if

        call calc_hrho_bloch_k(ik, sbe%rho(:, :, ik), p_rvnl_k, dv_k, hrho1_k)
        call calc_hrho_bloch_k(ik, hrho1_k, p_rvnl_k, dv_k, hrho2_k)
        call calc_hrho_bloch_k(ik, hrho2_k, p_rvnl_k, dv_k, hrho3_k)
        call calc_hrho_bloch_k(ik, hrho3_k, p_rvnl_k, dv_k, hrho4_k)

        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho1_k * (- zi * dt)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho2_k * (- zi * dt) ** 2 * (1d0 / 2d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho3_k * (- zi * dt) ** 3 * (1d0 / 6d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho4_k * (- zi * dt) ** 4 * (1d0 / 24d0)
    end do
    ! GW collision term (Phase 2), velocity gauge: explicit forward-Euler step
    ! applied to rho after the Taylor propagation.  OFF path untouched.
    if (yn_sbe_gw_collision == 'y') then
        call add_collision_vg(sbe%rho, gs%gamma_gw, gs%f0_ref, dt, &
          & sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
    end if
    return

contains


    !Calculate [H, rho] commutation:
    ! dv_k = the all-order nonlocal difference DeltaV(k,A) (VG completion);
    ! passed as an ARGUMENT (not host-associated) because it is an OpenMP
    ! PRIVATE copy of the enclosing loop -- host association would read the
    ! shared original.  Only referenced when do_exact (shared read-only).
    subroutine calc_hrho_bloch_k(ik, rho_k, p_k, dv_k, hrho_k)
        implicit none
        integer, intent(in) :: ik
        complex(8), intent(in) :: rho_k(nb, nb)
        complex(8), intent(in) :: p_k(nb, nb, 1:3)
        complex(8), intent(in) :: dv_k(nb, nb)
        complex(8), intent(out) :: hrho_k(nb, nb)
        integer :: idir
        !hrho = hrho + Ac(t) * (p * rho - rho * p)
        hrho_k(1:nb, 1:nb) = gs%delta_omega(1:nb, 1:nb, ik) * rho_k(1:nb, 1:nb)
        do idir=1, 3 !1:x, 2:y, 3:z
            ! hrho(1:nb, 1:nb, ik) = hrho(1:nb, 1:nb, ik) + Ac(idir) * (&
            ! & + matmul(gs%p_mod_matrix(1:nb, 1:nb, idir, ik), rho(1:nb, 1:nb, ik)) &
            ! & - matmul(rho(1:nb, 1:nb, ik), gs%p_mod_matrix(1:nb, 1:nb, idir, ik)) &
            ! & )

            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(+Ac(idir), 0d0), &
                p_k(:, :, idir),nb, &
                rho_k(:, :), nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :),nb)

            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(-Ac(idir), 0d0), &
                rho_k(:, :), nb, &
                p_k(:, :, idir),nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :), nb)

        end do !idir

        ! VG completion: hrho += [DeltaV, rho] (trace-conserving commutator;
        ! DeltaV is Hermitian, so the propagation stays unitary-consistent)
        if (do_exact) then
            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(+1d0, 0d0), &
                dv_k(:, :), nb, &
                rho_k(:, :), nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :), nb)
            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(-1d0, 0d0), &
                rho_k(:, :), nb, &
                dv_k(:, :), nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :), nb)
        end if
        return
    end subroutine calc_hrho_bloch_k
end subroutine

function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication
    use salmon_global
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    select case(trim(gauge_sbe))
    case ("velocity_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    case ("length_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    end select
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)

    return
end function calc_trace


!===================================================================
! Band-resolved excited population (yn_sbe_out_occ='y').
!   nex_b(b) = sum_k w_k ( Re rho_bb(k,t) - f0_b(k) ) / sum_k w_k
! with f0 = gs%occup (the ground-state occupation; equal to gs%f0_ref at init,
! but gs%occup is always allocated whereas gs%f0_ref exists only under the GW
! collision term).  Gauge-aware diagonal source MIRRORING calc_trace: rho for
! velocity_gauge, qnm_new for length_gauge (in LG sbe%rho is frozen at its
! init value -- only qnm/qnm_new are propagated -- but their diagonals coincide
! by construction; see prepare_qnm).  MPI: each rank sums its own k-slice, then
! comm_summation reduces the full nb-vector (collective; every rank must call).
! Normalized by sum(kweight) exactly as calc_trace, so summing the conduction
! columns reproduces the nelec column of the *_sbe_nex.data file.
!===================================================================
subroutine calc_band_population(sbe, gs, nex_b, icomm)
    use communication, only: comm_summation
    use salmon_global, only: gauge_sbe
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(out) :: nex_b(1:sbe%nb)
    integer, intent(in) :: icomm
    integer :: ik, ib
    real(8) :: wsum
    real(8) :: tmp_l(1:sbe%nb)

    tmp_l(:) = 0d0
    select case(trim(gauge_sbe))
    case ("length_gauge")
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, sbe%nb
                tmp_l(ib) = tmp_l(ib) + gs%kweight(ik) * &
                    & (real(sbe%qnm_new(ib, ib, ik)) - gs%occup(ib, ik))
            end do
        end do
    case default   ! velocity_gauge
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, sbe%nb
                tmp_l(ib) = tmp_l(ib) + gs%kweight(ik) * &
                    & (real(sbe%rho(ib, ib, ik)) - gs%occup(ib, ik))
            end do
        end do
    end select
    call comm_summation(tmp_l, nex_b, sbe%nb, icomm)
    wsum = sum(gs%kweight(:))
    nex_b(:) = nex_b(:) / wsum

    return
end subroutine calc_band_population


!===================================================================
! Energy grid for the nonequilibrium distribution (yn_sbe_out_occ='y').
! Built ONCE from gs%eigen (replicated on every rank, so the grid is identical
! everywhere without communication).  nbin bin centers span the window band
! energies padded by +-6*sigma; sigma = broadening (a few bin widths).
!===================================================================
subroutine sbe_edist_grid(gs, nb, nbin, e_lo, de, sigma)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb, nbin
    real(8), intent(out) :: e_lo, de, sigma
    real(8) :: emin, emax, espan, e_hi

    emin = minval(gs%eigen(1:nb, 1:gs%nk))
    emax = maxval(gs%eigen(1:nb, 1:gs%nk))
    espan = max(emax - emin, 1d-6)
    ! Fixed relative padding, then tie the broadening to the ACTUAL bin width
    ! (sigma = 2*de) so the Gaussian is always well-resolved on the grid,
    ! independent of nbin -- avoids the sigma<<de undersampling that would
    ! collapse the deposited weight for a small nbin.
    e_lo = emin - 0.05d0 * espan
    e_hi = emax + 0.05d0 * espan
    de = (e_hi - e_lo) / dble(nbin - 1)
    sigma = 2d0 * de

    return
end subroutine sbe_edist_grid


!===================================================================
! Energy-resolved nonequilibrium occupation distribution (yn_sbe_out_occ='y').
!   dist(j) = sum_{b,k} w_k Re rho_bb(k,t) g(e_j - eps_b(k)) / sum_k w_k
! e_j = e_lo + (j-1)*de, g = unit-normalized Gaussian of width sigma.  Each
! (b,k) contribution is deposited only on the +-5*sigma bin window around its
! eigenvalue (cost ~ nb*nk_slice*O(sigma/de), not nb*nk_slice*nbin).  Same
! gauge-aware diagonal source as calc_band_population.  MPI: local k-slice
! partial histogram, comm_summation over nbin (collective; all ranks call).
! Normalized by sum(kweight) so the energy integral gives electrons/cell.
!===================================================================
subroutine calc_energy_distribution(sbe, gs, e_lo, de, nbin, sigma, dist, icomm)
    use communication, only: comm_summation
    use salmon_global, only: gauge_sbe
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: e_lo, de, sigma
    integer, intent(in) :: nbin
    real(8), intent(out) :: dist(1:nbin)
    integer, intent(in) :: icomm
    integer :: ik, ib, j, jc, jlo, jhi, nsig
    real(8) :: occ, eps, ec, arg, gnorm, inv2s2, wsum
    logical :: is_lg
    real(8), allocatable :: dist_l(:)

    allocate(dist_l(1:nbin))
    dist_l(:) = 0d0
    is_lg = (trim(gauge_sbe) == "length_gauge")
    inv2s2 = 1d0 / (2d0 * sigma * sigma)
    gnorm = 1d0 / (sqrt(2d0 * pi) * sigma)
    nsig = int(5d0 * sigma / de) + 1   ! +-5 sigma deposition window (bins)

    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            if (is_lg) then
                occ = real(sbe%qnm_new(ib, ib, ik))
            else
                occ = real(sbe%rho(ib, ib, ik))
            end if
            eps = gs%eigen(ib, ik)
            jc = nint((eps - e_lo) / de) + 1
            jlo = max(1, jc - nsig)
            jhi = min(nbin, jc + nsig)
            do j = jlo, jhi
                ec = e_lo + dble(j - 1) * de
                arg = ec - eps
                dist_l(j) = dist_l(j) + gs%kweight(ik) * occ * gnorm * exp(-arg * arg * inv2s2)
            end do
        end do
    end do
    call comm_summation(dist_l, dist, nbin, icomm)
    wsum = sum(gs%kweight(:))
    dist(:) = dist(:) / wsum
    deallocate(dist_l)

    return
end subroutine calc_energy_distribution


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    ! real(8) :: kvec(1:3)
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                !& + dot_product(kvec(:), Ac(:))
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return
end function calc_energy

subroutine prepare_qnm(sbe, gs, icomm)
  use salmon_global, only: epdir_re1, am_s, sbe_lg_diag, sbe_lg_degen
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  integer, intent(in) :: icomm
  complex(8) :: dnm(sbe%nb, sbe%nb, sbe%nk)
  integer :: nb,nk
  integer :: ik,ib,jb,ii,jj

  nk = sbe%nk
  nb = sbe%nb

  allocate(sbe%qnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  if (.not. uses_xfull_links(sbe_lg_degen)) then
    allocate(sbe%grad_qnm(sbe%nb, sbe%nb, 1:3, sbe%nk))
  end if
  allocate(sbe%qnm_new(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%dqnm_stock(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max, am_s))
  allocate(sbe%dnm_i(sbe%nb, sbe%nb, 1:3, sbe%ik_min:sbe%ik_max))
  allocate(sbe%abs_dnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%exp_iphi(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
    sbe%qnm(ib, jb, ik) = 0.d0
    sbe%qnm_new(ib, jb, ik) = 0.d0
    dnm(ib, jb, ik)=0.d0
    sbe%abs_dnm(ib, jb, ik)=0.d0
    sbe%exp_iphi(ib, jb, ik)=0.d0
  end do
  end do
  end do

  do ii=1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
      sbe%dqnm_stock(ib, jb, ik, ii) = 0.d0
  end do
  end do
  end do
  end do

  if (.not. uses_xfull_links(sbe_lg_degen)) then
    !$omp parallel do default(shared) private(ik, ib, jb, jj)
    do ik = sbe%ik_min, sbe%ik_max
      do jj=1,3
        do ib=1,nb
          do jb=1,nb
            dnm(ib, jb, ik) = dnm(ib, jb, ik) + epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
          end do
        end do
      end do
    end do

    !$omp parallel do default(shared) private(ik, ib, jb, jj) collapse(4)
    do ik = sbe%ik_min, sbe%ik_max
      do jj=1,3
        do ib=1,nb
          do jb=1,nb
            sbe%dnm_i(ib, jb, jj, ik) = epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
          end do
        end do
      end do
    end do
  end if

  sbe%abs_dnm=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        sbe%abs_dnm(ib, jb, ik) = abs(dnm(ib, jb, ik))
      end do
    end do
  end do

  ! LG-SBE diagnostic knockout (2): suppress near-degenerate dipole couplings.
  ! gs%delta_omega is (nb,nb,1:nk) on every rank; sbe%abs_dnm is the local
  ! k-slice (nb,nb,ik_min:ik_max), so slice delta_omega to make the WHERE
  ! conform (identical to the full array when nproc==1).
  if (sbe_lg_diag == 2 .or. sbe_lg_diag == 3) then
    where (abs(gs%delta_omega(:, :, sbe%ik_min:sbe%ik_max)) < 1.d-3) sbe%abs_dnm = 0.d0
  end if

  ! gicov X-full (R-1 representation): ALL off-diagonal coherences are carried
  ! in the PHYSICAL basis (qnm == rho).  Force exp_iphi=1 and abs_dnm=0 on
  ! every off-diagonal pair BEFORE the phase-fill below, so (i) the
  ! phase-fill's abs_dnm>=1.d-13 gate skips them and leaves exp_iphi=1, giving
  ! qnm_new == rho for every off-diagonal pair, and (ii) the dipole source
  ! terms (which multiply abs_dnm) vanish everywhere -- the interband coupling
  ! is supplied instead by the full-band gauge-covariant gradient (gicov_rhs),
  ! not a dipole source.  X-full has a SINGLE full-band block (block_id === 1),
  ! so "same-block" and "ib/=jb" now coincide; the loop body no longer
  ! references block_id.  gi/gifix/off are untouched (this whole block is
  ! guarded to the X-full modes gicov / gicov_int).
  if (uses_xfull_links(sbe_lg_degen)) then
    !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
    do ik = sbe%ik_min, sbe%ik_max
      do ib=1,nb
        do jb=1,nb
          if (ib /= jb) then
            sbe%abs_dnm(ib, jb, ik) = 0.d0
            sbe%exp_iphi(ib, jb, ik) = (1.d0, 0.d0)
          end if
        end do
      end do
    end do
  end if

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        if(sbe%abs_dnm(ib, jb, ik) >= 1.d-13) then
          sbe%exp_iphi(ib, jb, ik) = dnm(ib, jb, ik)/sbe%abs_dnm(ib, jb, ik)
        end if
        if(ib==jb)then
          sbe%qnm_new(ib, jb, ik) = sbe%rho(ib, jb, ik)
        else
          sbe%qnm_new(ib, jb, ik) = conjg(sbe%exp_iphi(ib, jb, ik)) *  &
                                 &  sbe%rho(ib, jb, ik)
        end if
      end do
    end do
  end do

end subroutine prepare_qnm

!===================================================================
! R-1 representation bridge (orientation-explicit).  qnm is the propagated
! variable; the covariant gradient needs PHYSICAL rho.
!   rho(i,j) =              qnm(i,j)                (i==j)
!            = exp_iphi(i,j) * qnm(i,j)             (i/=j)
!   qnm(i,j) =              rho(i,j)                (i==j)
!            = conjg(exp_iphi(i,j)) * rho(i,j)      (i/=j)
! In gicov (X-full), prepare_qnm sets exp_iphi=1 on EVERY off-diagonal pair
! (ib/=jb), not just same-block ones, so qnm == rho for all off-diagonal
! coherences -- there is no out-of-block dipole-derived phase left; the
! interband coupling is carried entirely by the full-band gauge-covariant
! gradient (gicov_rhs), not by a dipole-phase-carrying qnm/rho distinction.
!===================================================================
pure function rho_ij_from_q(sbe, ib, jb, ik) result(r)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  integer, intent(in) :: ib, jb, ik
  complex(8) :: r
  if (ib == jb) then
    r = sbe%qnm(ib, jb, ik)
  else
    r = sbe%exp_iphi(ib, jb, ik) * sbe%qnm(ib, jb, ik)
  end if
end function rho_ij_from_q

pure function q_ij_from_rho(sbe, rho_ij, ib, jb, ik) result(q)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  complex(8), intent(in) :: rho_ij
  integer, intent(in) :: ib, jb, ik
  complex(8) :: q
  if (ib == jb) then
    q = rho_ij
  else
    q = conjg(sbe%exp_iphi(ib, jb, ik)) * rho_ij
  end if
end function q_ij_from_rho

!===================================================================
! gicov RHS operator (Phase 3, X-full): instantaneous COHERENT drho/dt of the
! gauge-covariant length-gauge propagator, as a callable routine.  NO
! integrator / AB4 / dqnm_stock here (that is Task 5b) -- this returns the
! instantaneous physical drho only, so a RHS bug can be told apart from an
! integrator-stability issue.
!
!   drho = + sum_i c_i * D_cov[rho]_i          covariant transport (WHOLE field term) (1)
!          - i * delta_omega .* rho           band-energy coherent term       (2a; off-diag)
!          - rho / t_2                         legacy interband dephasing      (2b; off-diag)
!
! (1) covariant_grad_block returns D_cov rho = d_k rho - i[xi,rho] on the FULL
!     nb x nb density; it REPLACES the legacy grad_qnm term.  gs%u_transport is
!     fed AS-IS (orientation U, Task 2 -- NOT conjugate-transposed).  The field
!     weights c_i and axis spacings dk are lattice-dependent (DUAL PATH, see
!     the inline comment at the field-term block below): on a strictly
!     diagonal b_matrix, c_i = Efield(i) and the Cartesian k-spacing
!     dk_a = b_matrix(a,a)/num_kgrid(a) reproduce the legacy
!     grad_k_array_nb2d_dcomplex spacing (nabt = bNmat(:,4)/(b(a,a)/N_a))
!     bit-for-bit, so with U==I this reduces EXACTLY to the legacy gradient;
!     on a non-diagonal (hexagonal, ...) b_matrix, c_i = (E . a_i)/(2 pi) and
!     dk_i = 1/num_kgrid(i) project E . grad_k onto the reduced axes.  The
!     "+E*grad" field prefactor/sign matches dt_evolve_bloch_lg (:663,:685).
!     X-full (block_id === 1, the full nb x nb overlap polar-factored as ONE
!     block): xi is now the FULL-BAND gauge-covariant connection, so its
!     off-diagonal entries (xi_inter) already ARE the physical interband
!     dipole -- the separate analytic dipole commutator that Phase-3-fixed-block
!     gicov used to add here would double-count it, so it is DELETED (proven
!     equivalent, not merely "unneeded", by src/ssbe/test/test_gicov_xfull.f90's
!     Test D: a REAL sigma_y rotation gives zero diagonal Berry connection, so
!     xi_full = xi_inter, and covariant_grad_block(full) is shown to equal
!     covariant_grad_block(bare) - i[xi_inter,rho] to stencil-matched tol).
!     gs%d_matrix is NO LONGER READ by this routine.
! (2) energy -i*delta_omega*rho (= -i[H0,rho]) and the 1/t_2 dephasing are kept
!     exactly as the legacy LG RHS (:680,:690); when GW dephasing is active the
!     t_2 term is suppressed (5b applies the GW collision as a separate substep).
!
! rho is reconstructed from sbe%qnm via rho_ij_from_q and gathered across ranks
! (comm_summation) because the covariant stencil is non-local in k.  Returns
! PHYSICAL drho (full nb x nb x nk).  5b maps drho -> dqnm elementwise via
! q_ij_from_rho if it steps qnm (diag: dqnm=drho; off-diag: dqnm=conjg(exp_iphi)*drho).
!
! NOTE (deviation from the brief signature): icomm is added as a trailing
! argument.  The covariant gradient must see neighbouring-k densities that may
! live on other ranks, so a gather is mandatory -- exactly as the legacy
! dt_evolve_bloch_lg gathers qnm via comm_summation before grad_k.  drho is
! returned for the full k range (nb,nb,nk); a caller slices its local k-range.
!===================================================================
subroutine gicov_rhs(sbe, gs, Efield, drho, icomm)
  use salmon_global, only: num_kgrid, t_2, yn_sbe_gw_collision, sbe_deph_mode, &
                            sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_lg_degen_floor
  use degenerate_block_ssbe, only: covariant_grad_block
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: Efield(1:3)
  ! only the LOCAL k-slice of drho was ever produced or read (the RK4 stage
  ! vectors are local-only); declaring the full nk was a per-step nb^2*nk
  ! complex allocation on every rank for nothing.
  complex(8), intent(out) :: drho(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max)
  integer, intent(in) :: icomm
  integer :: nb, nk, ik, ib, jb, axis
  real(8) :: dk(1:3), cvec(1:3)
  logical :: deph_by_gw, use_general
  logical :: axis_active(3)
  complex(8) :: gterm
  integer :: p, jj
  integer(8) :: tc1, tc2, trate

  nb = sbe%nb
  nk = sbe%nk

  ! ---- (0) physical rho on the local k-slice + halo exchange.
  ! Replaces the full-nk allreduce: only the +-m_max-shell neighbour-k rho the
  ! covariant stencil reads is shipped (point-to-point, plan built once).
  ! gh_rho persists across calls; entries neither local nor halo stay 0 from
  ! the build and are never read.
  if (.not. gh_built) call build_gicov_halo(sbe, gs, icomm)
  if (gh_nb /= sbe%nb .or. gh_nk /= sbe%nk) then
    write(*,*) "ERROR gicov_rhs: halo plan built for a different solver size"
    error stop
  end if

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
      do jb = 1, nb
        gh_rho(ib, jb, gs%kmap(ik)) = rho_ij_from_q(sbe, ib, jb, ik)
      end do
    end do
  end do

  sbe_tmr_nrhs = sbe_tmr_nrhs + 1
  call system_clock(tc1, trate)
  do p = 1, gh_nsrc
    gh_reqr(p) = comm_irecv(gh_rb(p)%b, gh_src_rank(p), 0, icomm)
  end do
  do p = 1, gh_ndst
    do jj = 1, gh_dst_cnt(p)
      gh_sb(p)%b(:, :, jj) = gh_rho(:, :, gs%kmap(gh_dst_k(gh_dst_ofs(p) + jj)))
    end do
    gh_reqs(p) = comm_isend(gh_sb(p)%b, gh_dst_rank(p), 0, icomm)
  end do
  if (gh_nsrc > 0) then
    call comm_wait_all(gh_reqr(1:gh_nsrc))
    do p = 1, gh_nsrc
      do jj = 1, gh_src_cnt(p)
        gh_rho(:, :, gs%kmap(gh_src_k(gh_src_ofs(p) + jj))) = gh_rb(p)%b(:, :, jj)
      end do
    end do
  end if
  if (gh_ndst > 0) call comm_wait_all(gh_reqs(1:gh_ndst))
  call system_clock(tc2);  sbe_tmr_comm = sbe_tmr_comm + dble(tc2 - tc1) / dble(trate)

  ! ---- (1) covariant transport (physical), WHOLE field term: + sum_i c_i D_cov[rho]_i
  ! Field-term projection, DUAL PATH (non-orthogonal lattice support; design:
  ! gw/gw_design/specs/2026-07-05-gicov-nonorthogonal-bmatrix.md):
  !
  !  (legacy) b_matrix STRICTLY diagonal -- all six off-diagonal entries
  !    exactly 0.d0 (every current orthorhombic campaign): the ORIGINAL
  !    formula verbatim,
  !      dk(axis)   = b_matrix(axis,axis)/num_kgrid(axis)   (Cartesian spacing)
  !      cvec(axis) = Efield(axis)
  !    cvec(axis)*Dq below is then the SAME multiplication Efield(axis)*Dq
  !    used to be => bit-for-bit identical to the pre-dual-path code.
  !
  !  (general) b_matrix non-diagonal (hexagonal, ...): reduced-axis
  !    projection.  k = sum_i s_i b_i (s_i = reduced coordinate) and
  !    a_i . b_j = 2 pi delta_ij give grad_k s_i = a_i/(2 pi), hence
  !      E . grad_k rho = sum_i c_i drho/ds_i,   c_i = (E . a_i)/(2 pi),
  !    with the stencil differentiating in the REDUCED coordinate:
  !      dk(i) = 1/num_kgrid(i)   (= Delta s_i).
  !    a_i = gs%a_matrix(:, i) (columns; set unconditionally by
  !    calc_lattice_info at gs init).  axis_active is judged on c_i, NOT on
  !    Efield: for a hexagonal cell, E || x already activates BOTH in-plane
  !    reduced axes.  In the orthogonal limit this path is ALGEBRAICALLY
  !    identical to the legacy one (c_i Dq_red = (E_i/b_ii)(b_ii Dq_cart))
  !    but differs in ULPs (1/b vs a/2pi) -- hence the strict dual path.
  !    Orthogonal equivalence (~1e-12) + a 60-degree hexagonal analytic
  !    directional-derivative check: src/ssbe/test/test_gicov_hex.f90.
  !
  ! Skipping field-inactive axes (cvec(axis)==0 exactly): c*D_cov for them is 0,
  ! so their covariant gradient is not computed (linear pol on a diagonal
  ! lattice => 1 of 3 axes, circular in-plane => 2). Bit-for-bit identical to
  ! computing all axes.  covariant_grad_block itself is lattice-agnostic (a
  ! reduced-grid stencil scaled by dk); halo/needed are index-based (unchanged).
  ! Exact `/= 0.d0` comparison, INTENTIONAL (not an "epsilon would be safer"
  ! bug): for every orthorhombic cell, calc_lattice_info (gs_info_ssbe.f90)
  ! writes the six off-diagonal b_matrix entries as literal 0.d0, so this is
  ! an exact-vs-exact test that reliably selects the legacy branch -- and
  ! test_gicov_hex.f90's Part O pins that selection to be bit-for-bit with
  ! the pre-dual-path code (1e-12 gate) for exactly that reason: a tolerance
  ! here would let the legacy branch silently start taking the (algebraically
  ! equivalent, ULP-different) general path and break that guarantee.  A
  ! nominally-orthorhombic input whose off-diagonal b_matrix entries are only
  ! FP-noise-close to zero (not exactly 0.d0) will therefore take the GENERAL
  ! branch instead; that is by design (see the DUAL PATH comment above: both
  ! branches are mathematically equivalent in the diagonal limit, differing
  ! only in ULPs), not a defect to special-case with an epsilon.
  use_general = gicov_force_general_field .or. &
    & (gs%b_matrix(1, 2) /= 0.d0) .or. (gs%b_matrix(1, 3) /= 0.d0) .or. &
    & (gs%b_matrix(2, 1) /= 0.d0) .or. (gs%b_matrix(2, 3) /= 0.d0) .or. &
    & (gs%b_matrix(3, 1) /= 0.d0) .or. (gs%b_matrix(3, 2) /= 0.d0)
  if (.not. use_general) then
    do axis = 1, 3
      dk(axis) = gs%b_matrix(axis, axis) / dble(num_kgrid(axis))
      cvec(axis) = Efield(axis)
      axis_active(axis) = (Efield(axis) /= 0.d0)
    end do
  else
    do axis = 1, 3
      dk(axis) = 1.d0 / dble(num_kgrid(axis))
      cvec(axis) = (Efield(1) * gs%a_matrix(1, axis) &
        &         + Efield(2) * gs%a_matrix(2, axis) &
        &         + Efield(3) * gs%a_matrix(3, axis)) / (2.d0 * pi)
      axis_active(axis) = (cvec(axis) /= 0.d0)
    end do
  end if
  call system_clock(tc1)
  call covariant_grad_block(nb, nk, gs%nbvec, gs%bvec, num_kgrid, &
                            gs%u_transport, gh_rho, dk, gh_Dq, sbe%ik_min, sbe%ik_max, axis_active, &
                            gs%nks, gs%kmap)
  call system_clock(tc2);  sbe_tmr_cov = sbe_tmr_cov + dble(tc2 - tc1) / dble(trate)

  deph_by_gw = (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw')

  ! only the local k-slice of drho is produced (the caller uses only ik_min:ik_max);
  ! Dq was likewise computed only on [ik_min:ik_max] by covariant_grad_block.
  ! Zero only that slice -- drho outside it is undefined and never read
  ! (dt_evolve_bloch_lg_gicov's k1..k4 consumers are local-only), and the
  ! full-nk memset was a measurable per-call fixed cost.
  drho(:, :, sbe%ik_min:sbe%ik_max) = (0.d0, 0.d0)
  !$omp parallel do default(shared) private(ik, ib, jb, gterm)
  do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
      do jb = 1, nb
        ! (1) covariant transport (intraband + interband, full-band).
        ! cvec == Efield verbatim on the legacy (diagonal-b) path, so this sum
        ! is bit-for-bit the old Efield(1..3)*Dq expression there; on the
        ! general path cvec(i) = (E . a_i)/(2 pi) are the reduced-axis weights.
        gterm = cvec(1) * gh_Dq(ib, jb, 1, ik) + &
                cvec(2) * gh_Dq(ib, jb, 2, ik) + &
                cvec(3) * gh_Dq(ib, jb, 3, ik)
        drho(ib, jb, ik) = gterm
        ! (2) coherent band energy + phenomenological dephasing (off-diagonal only)
        if (ib /= jb) then
          drho(ib, jb, ik) = drho(ib, jb, ik) &
            & - zi * gs%delta_omega(ib, jb, ik) * gh_rho(ib, jb, gs%kmap(ik))
          ! Covariant T2: dephase only ENERGY-DISTINCT pairs. Inside a
          ! (near-)degenerate manifold the scalar -rho/t_2 relaxation is NOT
          ! U(N)-gauge-covariant -- it splits diagonal/off-diagonal, a
          ! partition U(N) rotations mix -- so it is skipped there;
          ! intra-manifold relaxation belongs to the covariant (GW collision)
          ! channel, not to a phenomenological scalar T2.
          !
          ! Shape (design: gw/gw_design/plans/2026-07-08-t2-gate-shape.md):
          ! 'step' reuses the ORIGINAL hard-skip arithmetic VERBATIM (same
          ! gh_rho/t_2 division order => byte-identical with
          ! sbe_t2_gate_theta defaulting to the old theta_off literal 2d-3;
          ! do NOT route it through sbe%t2_gate_w, which would change the
          ! rounding). 'gauss' Hadamard-multiplies the precomputed weight
          ! (built once in init_sbe_bloch_solver from the STATIC
          ! gs%delta_omega) into the same term.
          !
          ! Exact-degeneracy floor clamp (|Delta-omega| <= sbe_lg_degen_floor
          ! => weight 0) is spec-mandated for BOTH shapes. gauss gets it inside
          ! t2_gate_weight (precomputed). step needs it EXPLICITLY here: the
          ! bare `abs(dw) > theta` test alone does NOT protect exact-degeneracy
          ! pairs when theta < floor (e.g. the checker-valid theta=0d0, where
          ! 0 < |dw| <= floor would otherwise be dephased). The `.and.` keeps
          ! the DEFAULT path bit-identical: with theta=2d-3 > floor=1d-9,
          ! `|dw| > theta` already implies `|dw| > floor`, so the added clause
          ! is always .true. there and never changes the arithmetic.
          if (.not. deph_by_gw) then
            if (trim(sbe_t2_gate_shape) == 'step') then
              if (abs(gs%delta_omega(ib, jb, ik)) > sbe_t2_gate_theta .and. &
                & abs(gs%delta_omega(ib, jb, ik)) > sbe_lg_degen_floor) then
                drho(ib, jb, ik) = drho(ib, jb, ik) - gh_rho(ib, jb, gs%kmap(ik)) / t_2
              end if
            else   ! 'gauss'
              drho(ib, jb, ik) = drho(ib, jb, ik) &
                & - sbe%t2_gate_w(ib, jb, ik) * gh_rho(ib, jb, gs%kmap(ik)) / t_2
            end if
          end if
        end if
      end do
    end do
  end do
  !$omp end parallel do
end subroutine gicov_rhs

subroutine dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
  use salmon_global, only: am_s, num_kgrid, t_2, epdir_re1, &
    & yn_sbe_gw_collision, sbe_deph_mode, sbe_lg_diag, sbe_lg_degen
  use sbe_collision_gw, only: add_collision_diag, add_collision_offdiag
  use common_ssbe, only: grad_k_array_nb2d_dcomplex
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: E(1:3)
  real(8), intent(in) :: bj_am(8,8)
  real(8), intent(in) :: dt
  integer, intent(in) :: icomm
  integer :: nb,nk
  integer :: ik,ib,jb
  integer :: ii,lb
  complex(8) :: shift_vector(3)
  real(8) :: abs_E, proj_E
  complex(8) :: qnm_tmp1(sbe%nb, sbe%nb, sbe%nk)
  complex(8) :: qnm_tmp(sbe%nb, sbe%nb, sbe%nk)

  nb = sbe%nb
  nk = sbe%nk

  ! gicov (R-3, Task 5b): self-starting RK4 on the PHYSICAL density via gicov_rhs
  ! (imaginary-axis stability 2*sqrt(2)~2.83), NO dqnm_stock / Adams-Moulton history,
  ! NO bare Euler; the GW collision is a separate VG-style substep.  gi/gifix/off/VG
  ! fall through to the AB4 path below (byte-unchanged).
  if (uses_fd_gicov(sbe_lg_degen)) then
    call dt_evolve_bloch_lg_gicov(sbe, gs, E, dt, icomm)
    return
  end if
  if (uses_integral_gicov(sbe_lg_degen)) then
    ! integral (covariant-Houston) transport propagation is driven from
    ! realtime_ssbe via dt_evolve_bloch_lg_integral (it needs the precomputed
    ! mesh shift a(t) at the step midpoint, not just E); reaching the generic
    ! length-gauge dispatch here is a wiring bug.
    write(*, '(a)') "ERROR(dt_evolve_bloch_lg): sbe_lg_degen='gicov_int' must be " // &
        & "propagated via dt_evolve_bloch_lg_integral (realtime_ssbe wiring)."
    error stop 1
  end if

  shift_vector(:) = 0.d0 ! shift_vector is set to be 0

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm(ib, jb, ik) = sbe%qnm_new(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s-1
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, ii) = sbe%dqnm_stock(ib, jb, ik, ii+1)
  end do
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, am_s) = 0.d0
  end do
  end do
  end do

  qnm_tmp1=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    qnm_tmp1(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  call comm_summation(qnm_tmp1, qnm_tmp, nb*nb*nk, icomm)
  call grad_k_array_nb2d_dcomplex(nb,nk,gs%b_matrix,qnm_tmp,sbe%grad_qnm)

  ! LG-SBE diagnostic knockout (1): suppress the intraband k-gradient term
  if (sbe_lg_diag == 1 .or. sbe_lg_diag == 3) sbe%grad_qnm = 0.d0

  abs_E = sqrt(E(1)**2+E(2)**2+E(3)**2)
  proj_E = E(1)*epdir_re1(1) + E(2)*epdir_re1(2) + E(3)*epdir_re1(3) 

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    if(ib == jb) then
    ! qnn (diagonal part)
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, ib, ik, am_s) = E(1) * sbe%grad_qnm(ib, ib, 1, ik) + &
                                           E(2) * sbe%grad_qnm(ib, ib, 2, ik) + &
                                           E(3) * sbe%grad_qnm(ib, ib, 3, ik)
        do lb = 1,nb
          if(ib /= lb)then
            sbe%dqnm_stock(ib, ib, ik, am_s) = &
              & sbe%dqnm_stock(ib, ib, ik, am_s) + &
              & 2.d0 * proj_E * sbe%abs_dnm(ib, lb, ik) * &
              &        aimag(sbe%qnm(lb, ib, ik))
          end if
        end do
      end if
    else
    ! qnm (off-diagonal part)
      if (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw') then
        sbe%dqnm_stock(ib, jb, ik, am_s) = 0.d0      ! t_2 replaced by GW (added below)
      else
        sbe%dqnm_stock(ib, jb, ik, am_s) = - sbe%qnm(ib, jb, ik)/t_2
      end if
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & E(1) * sbe%grad_qnm(ib, jb, 1, ik) + &
          & E(2) * sbe%grad_qnm(ib, jb, 2, ik) + &
          & E(3) * sbe%grad_qnm(ib, jb, 3, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * gs%delta_omega(ib, jb, ik) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * (E(1) * shift_vector(1) + &
                  E(2) * shift_vector(2) + &
                  E(3) * shift_vector(3)) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & zi * proj_E * sbe%abs_dnm(ib, jb, ik) * &
          &     (sbe%qnm(ib, ib, ik) - sbe%qnm(jb, jb, ik))
        do lb = 1,nb
          if(lb /= ib .and. lb /= jb) then
            sbe%dqnm_stock(ib, jb, ik, am_s) = &
              & sbe%dqnm_stock(ib, jb, ik, am_s) - &
              & zi * proj_E * &
              &     (sbe%abs_dnm(ib, lb, ik) * &
              &      sbe%qnm(lb, jb, ik) - &
              &      sbe%abs_dnm(lb, jb, ik) * &
              &      sbe%qnm(ib, lb, ik)) * &
              &      sbe%exp_iphi(ib, lb, ik) * &
              &      sbe%exp_iphi(lb, jb, ik) * &
              &      conjg(sbe%exp_iphi(ib, jb, ik))
          end if
        end do
      end if
    end if
  end do
  end do
  end do

  ! GW collision term (Phase 2): diagonal QP-lifetime loss + off-diagonal
  ! dephasing into the same Adams-Moulton RHS slot.  OFF path untouched.
  if (yn_sbe_gw_collision == 'y') then
    call add_collision_diag(sbe%dqnm_stock(:, :, :, am_s), sbe%qnm, &
      & gs%gamma_gw, gs%f0_ref, sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max)
    call add_collision_offdiag(sbe%dqnm_stock(:, :, :, am_s), sbe%qnm, &
      & gs%gamma_gw, sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
  end if

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm_new(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
      sbe%qnm_new(ib, jb, ik) = sbe%qnm_new(ib, jb, ik) + &
       &  bj_am(am_s+1-ii,am_s) * sbe%dqnm_stock(ib, jb, ik, ii) * dt
  end do
  end do
  end do
  end do

end subroutine dt_evolve_bloch_lg

!===================================================================
! gicov integrator (Phase 3, Task 5b).  Steps the PHYSICAL density rho one dt with
! a self-starting classical RK4 built on the (Task 5a machine-exact) gicov_rhs:
!
!   rho0 = rho_ij_from_q(qnm)                        (physical rho at step start)
!   k1 = f(rho0);  k2 = f(rho0 + dt/2 k1)
!   k3 = f(rho0 + dt/2 k2);  k4 = f(rho0 + dt k3)     f == gicov_rhs
!   rho_new = rho0 + dt/6 (k1 + 2 k2 + 2 k3 + k4)
!
! RK4's imaginary-axis stability limit 2*sqrt(2)~2.83 (vs the legacy AB4's 0.43) is
! the R-3 choice; RK4 is SELF-STARTING so there is NO dqnm_stock / Adams-Moulton
! history and NO bare Euler.  gicov_rhs reads sbe%qnm, reconstructs physical rho
! via rho_ij_from_q, and gathers neighbouring-k slices internally (comm_summation),
! so each stage writes its intermediate physical rho back into sbe%qnm through the
! R-1 bridge q_ij_from_rho before its gicov_rhs call.  Because prepare_qnm made
! exp_iphi a time-constant UNIT phase (same-block exp_iphi=1), the qnm<->rho bridge
! is lossless, so stepping rho through qnm is exact.
!
! The GW collision is applied as a SEPARATE substep on the physical rho AFTER the
! coherent RK4 (mirrors dt_evolve_bloch's post-Taylor add_collision_vg); it is a
! no-op when yn_sbe_gw_collision/='y' (e.g. the U3 gate: t_2=1e30, GW off).
!===================================================================
subroutine dt_evolve_bloch_lg_gicov(sbe, gs, E, dt, icomm)
  use salmon_global, only: yn_sbe_gw_collision, sbe_deph_mode
  use sbe_collision_gw, only: add_collision_vg
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: E(1:3)
  real(8), intent(in) :: dt
  integer, intent(in) :: icomm
  integer :: nb, nk, ik, ib, jb

  nb = sbe%nb
  nk = sbe%nk

  ! RK4 stage vectors: PERSISTENT (module-save), allocated once on the LOCAL
  ! k-slice.  They used to be allocate/deallocate'd on EVERY step at the FULL
  ! nk -- four nb^2*nk complex arrays of malloc churn per step (~361 MB/step at
  ! the graphene k99 production size) for buffers that were only ever written
  ! and read on [ik_min:ik_max].  Nothing about the arithmetic changes.
  if (.not. gk_built .or. gk_nb /= nb .or. gk_lo /= sbe%ik_min .or. gk_hi /= sbe%ik_max) then
    if (allocated(gk_rho0)) deallocate(gk_rho0, gk_rnew, gk_k1, gk_k2, gk_k3, gk_k4)
    allocate(gk_rho0(nb, nb, sbe%ik_min:sbe%ik_max))
    allocate(gk_rnew(nb, nb, sbe%ik_min:sbe%ik_max))
    allocate(gk_k1(nb, nb, sbe%ik_min:sbe%ik_max))
    allocate(gk_k2(nb, nb, sbe%ik_min:sbe%ik_max))
    allocate(gk_k3(nb, nb, sbe%ik_min:sbe%ik_max))
    allocate(gk_k4(nb, nb, sbe%ik_min:sbe%ik_max))
    gk_nb = nb;  gk_lo = sbe%ik_min;  gk_hi = sbe%ik_max;  gk_built = .true.
  end if

  ! (0) advance the propagated state (previous step's qnm_new) into qnm; then
  !     reconstruct the physical gk_rho0 at step start on the local k-slice.
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    sbe%qnm(ib, jb, ik) = sbe%qnm_new(ib, jb, ik)
  end do
  end do
  end do
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    gk_rho0(ib, jb, ik) = rho_ij_from_q(sbe, ib, jb, ik)
  end do
  end do
  end do

  ! (1) RK4 stage 1 (sbe%qnm already holds gk_rho0's representation)
  call gicov_rhs(sbe, gs, E, gk_k1, icomm)
  ! (2) stage 2 at gk_rho0 + (dt/2) gk_k1
  call load_stage_qnm(sbe, gk_rho0, gk_k1, 0.5d0 * dt)
  call gicov_rhs(sbe, gs, E, gk_k2, icomm)
  ! (3) stage 3 at gk_rho0 + (dt/2) gk_k2
  call load_stage_qnm(sbe, gk_rho0, gk_k2, 0.5d0 * dt)
  call gicov_rhs(sbe, gs, E, gk_k3, icomm)
  ! (4) stage 4 at gk_rho0 + dt gk_k3
  call load_stage_qnm(sbe, gk_rho0, gk_k3, dt)
  call gicov_rhs(sbe, gs, E, gk_k4, icomm)

  ! (5) combine on the local slice
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    gk_rnew(ib, jb, ik) = gk_rho0(ib, jb, ik) + (dt / 6.d0) * &
      & (gk_k1(ib, jb, ik) + 2.d0 * gk_k2(ib, jb, ik) + 2.d0 * gk_k3(ib, jb, ik) + gk_k4(ib, jb, ik))
  end do
  end do
  end do

  ! (6) GW collision as a SEPARATE substep on the physical rho (VG-style).
  if (yn_sbe_gw_collision == 'y') then
    call add_collision_vg(gk_rnew, gs%gamma_gw, gs%f0_ref, dt, &
      & nb, nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
  end if

  ! (7) write stepped physical rho into qnm_new; restore qnm to the step-start
  !     representation (AB4 leaves qnm holding the step-start state).
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    sbe%qnm_new(ib, jb, ik) = q_ij_from_rho(sbe, gk_rnew(ib, jb, ik), ib, jb, ik)
    sbe%qnm(ib, jb, ik)     = q_ij_from_rho(sbe, gk_rho0(ib, jb, ik),    ib, jb, ik)
  end do
  end do
  end do

contains

  ! Write the RK stage density rho0 + a*kX into sbe%qnm (local slice) through the
  ! R-1 bridge so the next gicov_rhs sees the intermediate physical rho.
  subroutine load_stage_qnm(sbe, rho0, kX, a)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    complex(8), intent(in) :: rho0(nb, nb, sbe%ik_min:sbe%ik_max)
    complex(8), intent(in) :: kX(nb, nb, sbe%ik_min:sbe%ik_max)
    real(8), intent(in) :: a
    integer :: ik, ib, jb
    !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
    do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
    do jb = 1, nb
      sbe%qnm(ib, jb, ik) = q_ij_from_rho(sbe, rho0(ib, jb, ik) + a * kX(ib, jb, ik), ib, jb, ik)
    end do
    end do
    end do
  end subroutine load_stage_qnm

end subroutine dt_evolve_bloch_lg_gicov

!===================================================================
! Integral (covariant-Houston) transport propagation (sbe_lg_degen='gicov_int').
!
! Driven by realtime_ssbe, which analyses the WHOLE precomputed field
! trajectory: it confirms single-axis linear polarization (gicov_int_axis_single
! on q_i(t) = num_kgrid(i)*(a(t).a_i)/2pi), fixes the driven reduced axis and
! the pulse span j_max = ceil(max_t|q_axis(t)|), builds the transport cache
! once, then per step propagates with the moving-frame shift q_axis(t) evaluated
! at the STEP MIDPOINT.  The co-moving density rho~ lives in place in sbe%rho
! (== the physical rho on the X-full path: prepare_qnm set exp_iphi=1, so
! qnm==rho; the initial rho~(0)=rho(0)=gs occupation since a(0)=0 -> W=I).
!
! The single-step Wilson links and the velocities are needed at k out to
! +-j_cache along the driven axis.  They are NOT replicated (that full-nk
! replication is exactly what caps the k-mesh), so the chain is marched outward
! with a ROLLING one-step halo shift -- see the comment at the assembly loop.
!===================================================================
subroutine build_gicov_integral_cache(sbe, gs, axis, jmax, icomm)
  use salmon_global, only: num_kgrid
  use degenerate_block_ssbe, only: build_ik_neighbor, polar_unitary
  use gicov_integral_ssbe, only: gicov_int_transport_op, gicov_int_cache_bytes, &
                                 gicov_int_jmax_cache
  !$ use omp_lib
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  ! inout: the driven-axis halo grows gs%kmap / gs%u_transport by the remote
  ! nodes the chain walk visits (the link array itself is the halo store)
  type(s_sbe_gs_info), intent(inout) :: gs
  integer, intent(in) :: axis, jmax, icomm
  integer :: nb, nk, ik, n, i, istep, kk, krem, ierr, nth, jc
  integer :: bplus(3, 1), bminus(3, 1)
  integer, allocatable :: nbtab(:, :)
  complex(8), allocatable :: Wc(:, :), acc(:, :), Ohost(:, :), Yt(:, :)
  real(8) :: sigma_min
  integer(8) :: nbytes

  nb = sbe%nb;  nk = sbe%nk
  ! jmax is the PHYSICAL pulse span; the cache is padded by the interpolation
  ! halo so the degree-p stencil of an off-node x never reaches past its end.
  jc = gicov_int_jmax_cache(jmax)
  gi_nb = nb;  gi_axis = axis;  gi_jmax = jc;  gi_jmax_phys = jmax
  gi_ikmin = sbe%ik_min;  gi_ikmax = sbe%ik_max

  ! +/- driven-axis neighbour tables on the reduced uniform grid
  bplus = 0;   bplus(axis, 1) = 1
  bminus = 0;  bminus(axis, 1) = -1
  allocate(nbtab(1, nk))
  if (.not. allocated(gi_iknb_p)) allocate(gi_iknb_p(nk))
  if (.not. allocated(gi_iknb_m)) allocate(gi_iknb_m(nk))
  call build_ik_neighbor(num_kgrid, bplus,  1, nk, nbtab);  gi_iknb_p(:) = nbtab(1, :)
  call build_ik_neighbor(num_kgrid, bminus, 1, nk, nbtab);  gi_iknb_m(:) = nbtab(1, :)
  deallocate(nbtab)

  ! transport cache over the LOCAL k-slice only (never full nk); int64 sizing.
  ! Sized on the PADDED span jc -- the halo shifts are really allocated, so the
  ! estimate must count them (the lint memory gate applies the same padding).
  nbytes = gicov_int_cache_bytes(jc, nb, gi_ikmax - gi_ikmin + 1)
  allocate(gi_Ht(nb, nb, -jc:jc, gi_ikmin:gi_ikmax), stat=ierr)
  if (ierr /= 0) then
    write(*, '(a,i0,a)') "ERROR build_gicov_integral_cache: gi_Ht alloc failed (", &
      & int(nbytes / 2_8**20, 4), " MiB/rank requested)"
    error stop 1
  end if
  allocate(gi_vt(nb, nb, -jc:jc, 3, gi_ikmin:gi_ikmax), stat=ierr)
  if (ierr /= 0) then
    write(*, '(a)') "ERROR build_gicov_integral_cache: gi_vt alloc failed"
    error stop 1
  end if
  allocate(gi_Ft(nb, nb, -jc:jc, gi_ikmin:gi_ikmax), stat=ierr)
  if (ierr /= 0) then
    write(*, '(a)') "ERROR build_gicov_integral_cache: gi_Ft alloc failed"
    error stop 1
  end if

  ! ---- transport-chain assembly over the driven-axis halo --------------
  ! The chain W(kappa, kappa+n) needs the single-step links AND the velocity at
  ! REMOTE k, out to +-j_cache along the driven axis.  Those fields are no longer
  ! replicated full-nk (that replication is the ~248*nb^2*nk B/rank term this
  ! change removes), so the exact set the walk will touch is fetched ONCE here.
  !
  ! The walk BELOW is the original one, statement for statement: same operand
  ! order, same matmul expressions, same polar_unitary input.  Only the k-index
  ! is redirected (gs%kmap for the links; the p/rvnl halo for a non-owned node).
  ! That is deliberate and load-bearing: re-associating the chain into a rolling
  ! recursion is algebraically identical but NOT bit-identical -- gfortran emits a
  ! different matmul for a conj-transposed local-array operand than for the
  ! derived-type-component operand, and the resulting 1e-15 drift in W propagated
  ! into every cached H~ and v~ (measured; the forward walk, which has no
  ! conj-transpose, stayed exact).  Bit-identity is the requirement, so the
  ! arithmetic does not move.
  call gi_build_axis_halo(gs, nk, axis, jc, gi_ikmin, gi_ikmax, gi_iknb_p, gi_iknb_m, &
                        & sbe%flag_vnl_correction, icomm)

  allocate(Wc(nb, nb), acc(nb, nb), Ohost(nb, nb), Yt(nb, nb))
  do ik = gi_ikmin, gi_ikmax
    do n = -jc, jc
      ! W(kappa, kappa+n) by walking |n| single-step links
      acc = (0d0, 0d0)
      do i = 1, nb
        acc(i, i) = (1d0, 0d0)
      end do
      kk = ik
      if (n > 0) then
        do istep = 1, n
          acc = matmul(acc, gs%u_transport(:, :, axis, gs%kmap(kk)))
          kk = gi_iknb_p(kk)
        end do
      else if (n < 0) then
        do istep = 1, -n
          kk = gi_iknb_m(kk)
          acc = matmul(acc, conjg(transpose(gs%u_transport(:, :, axis, gs%kmap(kk)))))
        end do
      end if
      krem = kk
      call polar_unitary(acc, nb, Wc, sigma_min, ierr)
      if (ierr /= 0) Wc = acc                 ! transient window leak: keep raw product
      ! transported band-energy matrix H~_n = Wc diag(eigen(remote)) Wc^dag
      Ohost = (0d0, 0d0)
      do i = 1, nb
        Ohost(i, i) = cmplx(gs%eigen(i, krem), 0d0, 8)
      end do
      call gicov_int_transport_op(Wc, Ohost, nb, Yt)
      gi_Ht(:, :, n, ik) = Yt
      ! transported GS occupation reference F~0_n = Wc diag(f0(remote)) Wc^dag.
      ! gs%occup is the REAL per-(band,k) ground-state occupation (fractional and
      ! k-varying for a metal), so this is the metal-correct f0(x) as well.
      Ohost = (0d0, 0d0)
      do i = 1, nb
        Ohost(i, i) = cmplx(gs%occup(i, krem), 0d0, 8)
      end do
      call gicov_int_transport_op(Wc, Ohost, nb, Yt)
      gi_Ft(:, :, n, ik) = Yt
      ! transported velocity v~_{n,i} = Wc v_i(remote) Wc^dag, i=1..3
      do i = 1, 3
        if (krem >= gi_ikmin .and. krem <= gi_ikmax) then
          if (sbe%flag_vnl_correction) then
            Ohost = gs%p_tm_matrix(:, :, i, krem) + gs%rvnl_tm_matrix(:, :, i, krem)
          else
            Ohost = gs%p_tm_matrix(:, :, i, krem)
          end if
        else
          ! non-owned node: the SAME two matrices, fetched from the halo (an
          ! elementwise sum -- no summation order to perturb)
          if (sbe%flag_vnl_correction) then
            Ohost = gih_p(:, :, i, gih_map(krem)) + gih_r(:, :, i, gih_map(krem))
          else
            Ohost = gih_p(:, :, i, gih_map(krem))
          end if
        end if
        call gicov_int_transport_op(Wc, Ohost, nb, Yt)
        gi_vt(:, :, n, i, ik) = Yt
      end do
    end do
  end do
  deallocate(Wc, acc, Ohost, Yt)
  call gi_free_axis_halo()

  ! per-thread eigensolver + interp work bank: allocated ONCE here, OUTSIDE any
  ! OMP loop (frtpx: no heap inside the loop, no per-thread automatic arrays)
  nth = 1
  !$ nth = omp_get_max_threads()
  gi_lcwork = max(1, 2 * nb - 1)
  allocate(gi_w_eps(nb, 0:nth-1), gi_w_rw(max(1, 3*nb-2), 0:nth-1))
  allocate(gi_w_P(nb, nb, 0:nth-1), gi_w_R(nb, nb, 0:nth-1))
  allocate(gi_w_cw(gi_lcwork, 0:nth-1))
  allocate(gi_w_H(nb, nb, 0:nth-1), gi_w_rout(nb, nb, 0:nth-1), gi_w_v(nb, nb, 3, 0:nth-1))
  allocate(gi_w_blk(nb, 0:nth-1))



  gi_built = .true.
end subroutine build_gicov_integral_cache


!===================================================================
! Driven-axis transport halo for build_gicov_integral_cache.
!
! The chain walk reads, for every OWNED k, the single-step links and the
! velocity at kappa + n*e_axis for n in [-j_cache, j_cache].  Under klocal those
! nodes may belong to another rank, so fetch exactly that set, once, before the
! walk:
!   * the LINKS are appended as extra SLOTS of gs%u_transport (gs%kmap grows to
!     cover them), so the walk keeps indexing the very same derived-type
!     component it always did -- only the slot lookup is new.  That is what
!     makes the chain bit-identical (see the note at the walk).
!   * the VELOCITIES (p_tm, and rvnl when the nonlocal correction is on) go into
!     side buffers with their own compact map; they are only ever summed
!     elementwise, so no expression shape matters there.
!
! Size: the set is the axis shell of the owned slice, i.e. ~2*j_cache extra k per
! boundary row when the driven axis is the fastest one (the production case) --
! kilobytes.  It is freed as soon as the cache is built; from then on gicov_int
! reads NOTHING but its own local cache.
!===================================================================
subroutine gi_build_axis_halo(gs, nk, axis, jc, ik_lo, ik_hi, knb_p, knb_m, want_rvnl, icomm)
  use degenerate_block_ssbe, only: build_halo_lists
  implicit none
  type(s_sbe_gs_info), intent(inout) :: gs
  integer, intent(in) :: nk, axis, jc, ik_lo, ik_hi, icomm
  integer, intent(in) :: knb_p(nk), knb_m(nk)
  logical, intent(in) :: want_rvnl
  integer :: irank, nproc, nb, o, ik, kk, n, p, jj, i, nslot0, nh
  integer, allocatable :: itbl_min(:), itbl_max(:)
  logical, allocatable :: needed_all(:, :)
  complex(8), allocatable :: utmp(:, :, :, :)
  type(t_gh_buf), allocatable :: sb(:), rb(:)
  integer, allocatable :: reqs(:), reqr(:)
  integer :: ncomp

  nb = gs%nb
  call comm_get_groupinfo(icomm, irank, nproc)
  allocate(gih_map(nk));  gih_map(:) = 0
  gih_nh = 0
  if (.not. gs%klocal) return          ! replicated: every node is already here

  ! ---- needed-set of EVERY rank (deterministic, no metadata exchange) ----
  allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
  call split_range(1, nk, nproc, itbl_min, itbl_max)
  allocate(needed_all(nk, 0:nproc-1))
  needed_all = .false.
  do o = 0, nproc - 1
    if (itbl_max(o) < itbl_min(o)) cycle
    do ik = itbl_min(o), itbl_max(o)
      needed_all(ik, o) = .true.
      kk = ik
      do n = 1, jc
        kk = knb_p(kk);  needed_all(kk, o) = .true.
      end do
      kk = ik
      do n = 1, jc
        kk = knb_m(kk);  needed_all(kk, o) = .true.
      end do
    end do
  end do
  call build_halo_lists(nk, nproc, irank, itbl_min, itbl_max, needed_all, &
                        gih_nsrc, gih_src_rank, gih_src_cnt, gih_src_ofs, gih_src_k, &
                        gih_ndst, gih_dst_rank, gih_dst_cnt, gih_dst_ofs, gih_dst_k)
  deallocate(itbl_min, itbl_max, needed_all)

  ! ---- grow gs%kmap / gs%u_transport by the received nodes ----
  nslot0 = gs%nks
  do p = 1, gih_nsrc
    do jj = 1, gih_src_cnt(p)
      ik = gih_src_k(gih_src_ofs(p) + jj)
      if (gs%kmap(ik) == 0) then
        gs%nks = gs%nks + 1
        gs%kmap(ik) = gs%nks
      end if
      if (gih_map(ik) == 0) then
        gih_nh = gih_nh + 1
        gih_map(ik) = gih_nh
      end if
    end do
  end do
  if (gs%nks > nslot0) then
    allocate(utmp(nb, nb, 3, gs%nks))
    utmp(:, :, :, 1:nslot0) = gs%u_transport(:, :, :, 1:nslot0)
    utmp(:, :, :, nslot0+1:gs%nks) = (0d0, 0d0)
    call move_alloc(utmp, gs%u_transport)
  end if
  nh = max(1, gih_nh)
  allocate(gih_p(nb, nb, 3, nh))
  allocate(gih_r(nb, nb, 3, nh))
  gih_p = (0d0, 0d0);  gih_r = (0d0, 0d0)

  ! ---- ONE exchange carrying the links AND the velocity components ----
  ncomp = 6
  if (want_rvnl) ncomp = 9
  allocate(rb(gih_nsrc), sb(gih_ndst))
  do p = 1, gih_nsrc
    allocate(rb(p)%b(nb, nb, ncomp * gih_src_cnt(p)))
  end do
  do p = 1, gih_ndst
    allocate(sb(p)%b(nb, nb, ncomp * gih_dst_cnt(p)))
  end do
  allocate(reqr(max(1, gih_nsrc)), reqs(max(1, gih_ndst)))
  do p = 1, gih_nsrc
    reqr(p) = comm_irecv(rb(p)%b, gih_src_rank(p), 4, icomm)
  end do
  do p = 1, gih_ndst
    do jj = 1, gih_dst_cnt(p)
      ik = gih_dst_k(gih_dst_ofs(p) + jj)
      do i = 1, 3
        sb(p)%b(:, :, ncomp*(jj-1) + i)     = gs%u_transport(:, :, i, gs%kmap(ik))
        sb(p)%b(:, :, ncomp*(jj-1) + 3 + i) = gs%p_tm_matrix(:, :, i, ik)
        if (want_rvnl) sb(p)%b(:, :, ncomp*(jj-1) + 6 + i) = gs%rvnl_tm_matrix(:, :, i, ik)
      end do
    end do
    reqs(p) = comm_isend(sb(p)%b, gih_dst_rank(p), 4, icomm)
  end do
  if (gih_nsrc > 0) then
    call comm_wait_all(reqr(1:gih_nsrc))
    do p = 1, gih_nsrc
      do jj = 1, gih_src_cnt(p)
        ik = gih_src_k(gih_src_ofs(p) + jj)
        do i = 1, 3
          gs%u_transport(:, :, i, gs%kmap(ik)) = rb(p)%b(:, :, ncomp*(jj-1) + i)
          gih_p(:, :, i, gih_map(ik))          = rb(p)%b(:, :, ncomp*(jj-1) + 3 + i)
          if (want_rvnl) gih_r(:, :, i, gih_map(ik)) = rb(p)%b(:, :, ncomp*(jj-1) + 6 + i)
        end do
      end do
    end do
  end if
  if (gih_ndst > 0) call comm_wait_all(reqs(1:gih_ndst))
  do p = 1, gih_nsrc
    deallocate(rb(p)%b)
  end do
  do p = 1, gih_ndst
    deallocate(sb(p)%b)
  end do
  deallocate(rb, sb, reqr, reqs)

  ! fail-closed: every node the walk will visit must now resolve
  do ik = ik_lo, ik_hi
    kk = ik
    do n = 1, jc
      kk = knb_p(kk)
      if (gs%kmap(kk) == 0) then
        write(*,*) "ERROR gi_build_axis_halo: +axis node has no link slot, k =", kk
        error stop 1
      end if
    end do
    kk = ik
    do n = 1, jc
      kk = knb_m(kk)
      if (gs%kmap(kk) == 0) then
        write(*,*) "ERROR gi_build_axis_halo: -axis node has no link slot, k =", kk
        error stop 1
      end if
    end do
  end do

  if (irank == 0) write(*, '(a,i0,a)') &
    & "# gicov_int: driven-axis transport halo = ", gih_nh, " remote k / rank (freed after the cache build)"
end subroutine gi_build_axis_halo


subroutine gi_free_axis_halo()
  implicit none
  if (allocated(gih_p))        deallocate(gih_p)
  if (allocated(gih_r))        deallocate(gih_r)
  if (allocated(gih_map))      deallocate(gih_map)
  if (allocated(gih_src_rank)) deallocate(gih_src_rank, gih_src_cnt, gih_src_ofs, gih_src_k)
  if (allocated(gih_dst_rank)) deallocate(gih_dst_rank, gih_dst_cnt, gih_dst_ofs, gih_dst_k)
  gih_nsrc = 0;  gih_ndst = 0;  gih_nh = 0
end subroutine gi_free_axis_halo


!===================================================================
! k-ROUGHNESS OF THE DENSITY, measured with the Wilson shift (SBE_KROUGH).
!
! What it answers: the covariant noise floor is set by how much of the density
! sits at the SHORT-wavelength end of the k-mesh -- the part the interpolation
! (or the finite-difference stencil) cannot represent.  A raw FFT of rho(k)
! cannot measure that, because rho(k) is expressed in a k-dependent band gauge:
! its Fourier content is gauge NOISE as much as physics.  The Wilson shift is
! the gauge-covariant replacement for the plain shift:
!
!     (S_a X)(k) = U_a(k) X(k + e_a) U_a(k)^dagger ,
!     L_a        = 2I - S_a - S_a^dagger            (covariant discrete -d^2/dk^2)
!     G_a        = <X, L_a X> / (4 <X, X>)          in [0, 1]
!
! G_a = 0 means X is perfectly smooth along axis a (in the parallel-transported
! frame); G_a = 1 means it alternates at the Nyquist wavelength.  Both S_a and
! the inner product are invariant under a k-local gauge rotation
! (X -> W^H X W, U -> W(k)^H U W(k+e)), so G_a is a PHYSICAL number.
!
! Using <X, L X> = 2(<X,X> - Re<X, S X>) the whole thing costs ONE forward
! Wilson shift, i.e. a single one-layer halo exchange -- no FFT, no eigensolve.
!
! Reported for two probes, because they answer different questions:
!   X = drho = rho~ - rho~(0)   the EXCITED density -- what actually radiates;
!   X = [H~, rho~]              the local generator of the dynamics -- the part
!                               the propagator differentiates/interpolates, so
!                               its roughness is what the floor tracks.
!
! R = the fraction of |X|^2 living in the upper half of the covariant spectrum;
! G alone bounds it (Markov/Chebyshev on the spectral measure of L in [0,4]):
!     max(0, (4G-1)/3)  <=  R  <=  min(1, 4G).
! Both bounds are printed, so a floor claim can be read off without inverting a
! spectral density that was never computed.
!===================================================================
subroutine sbe_krough_diag(sbe, gs, q_now, gval, rlo, rhi, icomm)
  use salmon_global, only: num_kgrid, sbe_lg_degen
  use degenerate_block_ssbe, only: build_ik_neighbor, find_bvec, build_halo_lists
  use gicov_integral_ssbe, only: gicov_int_stencil, gicov_int_interp_p, gicov_int_nsten
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: q_now
  real(8), intent(out) :: gval(3, 2), rlo(3, 2), rhi(3, 2)   ! (axis, probe): 1 = drho, 2 = [H,rho]
  integer, intent(in) :: icomm
  integer :: nb, nk, ik, ib, jb, axis, iv, iprobe, ierr
  integer :: iv_axis(3)
  integer :: nodes(gicov_int_nsten), nsten
  real(8) :: wts(gicov_int_nsten)
  real(8) :: acc_l(7), acc_g(7), xx, sx
  complex(8), allocatable :: X(:, :, :), Hk(:, :), tmp1(:, :), tmp2(:, :)
  complex(8) :: rr

  nb = sbe%nb;  nk = sbe%nk
  gval = 0d0;  rlo = 0d0;  rhi = 0d0

  call kr_build_plan(sbe, gs, icomm)
  allocate(X(nb, nb, max(1, krp_nx)), Hk(nb, nb), tmp1(nb, nb), tmp2(nb, nb))

  iv_axis(1) = find_bvec(gs%bvec, gs%nbvec, 1, 0, 0)
  iv_axis(2) = find_bvec(gs%bvec, gs%nbvec, 0, 1, 0)
  iv_axis(3) = find_bvec(gs%bvec, gs%nbvec, 0, 0, 1)

  ! rho~(t=0): captured by sbe_krough_init BEFORE the first step.  Capturing it
  ! lazily on first USE would make the excitation probe read zero at its own
  ! first sample (the reference would be the state at that sample), which is a
  ! silently plausible and completely wrong curve.
  if (.not. kr_built) then
    write(*, '(a)') "ERROR(sbe_krough_diag): rho~(0) reference not captured " // &
      & "(sbe_krough_init must run before the propagation)."
    error stop 1
  end if

  if (uses_integral_gicov(sbe_lg_degen)) then
    call gicov_int_stencil(q_now, 1d0, gi_jmax, nodes, wts, nsten, ierr)
    if (ierr /= 0) then
      write(*, '(a)') "ERROR(sbe_krough_diag): mesh shift leaves the cached span."
      error stop 1
    end if
  end if

  do iprobe = 1, 2
    ! ---- build the probe field X on the owned k, in the SLOT layout ----
    do ik = sbe%ik_min, sbe%ik_max
      call kr_rho_at(sbe, nb, ik, tmp1)
      if (iprobe == 1) then
        X(:, :, krp_map(ik)) = tmp1(:, :) - kr_rho0(:, :, ik)
      else
        if (uses_integral_gicov(sbe_lg_degen)) then
          ! instantaneous co-moving Hamiltonian at x = kappa - a(t)
          call gicov_int_interp_p(gi_Ht(:, :, :, ik), nb, -gi_jmax, gi_jmax, &
                                & nodes, wts, nsten, Hk)
          X(:, :, krp_map(ik)) = matmul(Hk, tmp1) - matmul(tmp1, Hk)
        else
          ! FD gicov: H0 = diag(eigen), so [H0, rho]_ij = (eps_i - eps_j) rho_ij
          do jb = 1, nb
            do ib = 1, nb
              X(ib, jb, krp_map(ik)) = gs%delta_omega(ib, jb, ik) * tmp1(ib, jb)
            end do
          end do
        end if
      end if
    end do

    call kr_exchange(nb, X, icomm)

    ! ---- <X,X> and Re <X, S_a X> for the three axes ----
    acc_l(:) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      rr = (0d0, 0d0)
      do jb = 1, nb
        do ib = 1, nb
          rr = rr + conjg(X(ib, jb, krp_map(ik))) * X(ib, jb, krp_map(ik))
        end do
      end do
      acc_l(1) = acc_l(1) + gs%kweight(ik) * dble(rr)
      do axis = 1, 3
        iv = iv_axis(axis)
        if (iv == 0) cycle
        if (num_kgrid(axis) < 2) cycle          ! singleton axis: no shift exists
        ! S_a X (k) = U_a(k) X(k+e_a) U_a(k)^H
        call zgemm('N', 'N', nb, nb, nb, (1d0, 0d0), &
                 & gs%u_transport(:, :, axis, gs%kmap(ik)), nb, &
                 & X(:, :, krp_map(krp_nbr(axis, ik))), nb, (0d0, 0d0), tmp1, nb)
        call zgemm('N', 'C', nb, nb, nb, (1d0, 0d0), tmp1, nb, &
                 & gs%u_transport(:, :, axis, gs%kmap(ik)), nb, (0d0, 0d0), tmp2, nb)
        rr = (0d0, 0d0)
        do jb = 1, nb
          do ib = 1, nb
            rr = rr + conjg(X(ib, jb, krp_map(ik))) * tmp2(ib, jb)
          end do
        end do
        acc_l(1 + axis) = acc_l(1 + axis) + gs%kweight(ik) * dble(rr)
      end do
    end do
    call comm_summation(acc_l, acc_g, 7, icomm)

    xx = acc_g(1)
    do axis = 1, 3
      if (xx <= 0d0) cycle
      ! an axis with no +shift column, or a singleton axis, carries S_a = I, so
      ! <X, L_a X> = 0 exactly: report G = 0 rather than the 0.5 an unaccumulated
      ! sx would fake.
      if (iv_axis(axis) == 0 .or. num_kgrid(axis) < 2) then
        gval(axis, iprobe) = 0d0;  rlo(axis, iprobe) = 0d0;  rhi(axis, iprobe) = 0d0
        cycle
      end if
      sx = acc_g(1 + axis)
      gval(axis, iprobe) = (xx - sx) / (2d0 * xx)
      rlo(axis, iprobe) = max(0d0, (4d0 * gval(axis, iprobe) - 1d0) / 3d0)
      rhi(axis, iprobe) = min(1d0, 4d0 * gval(axis, iprobe))
    end do
  end do

  deallocate(X, Hk, tmp1, tmp2)
end subroutine sbe_krough_diag


! Capture the t=0 co-moving density, the reference of the excitation probe.
! Called once, before the first step, when SBE_KROUGH is on.
subroutine sbe_krough_init(sbe)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  integer :: ik
  if (kr_built) return
  allocate(kr_rho0(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  do ik = sbe%ik_min, sbe%ik_max
    call kr_rho_at(sbe, sbe%nb, ik, kr_rho0(:, :, ik))
  end do
  kr_built = .true.
end subroutine sbe_krough_init


! Physical rho~ of the CURRENT propagated state on one owned k.  X-full carries
! qnm == rho (prepare_qnm forces exp_iphi = 1), and the integral propagator
! mirrors its result into rho / qnm / qnm_new alike, so qnm_new is the state in
! BOTH X-full modes.
subroutine kr_rho_at(sbe, nb, ik, rk)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  integer, intent(in) :: nb, ik
  complex(8), intent(out) :: rk(nb, nb)
  integer :: ib, jb
  do jb = 1, nb
    do ib = 1, nb
      rk(ib, jb) = rho_ij_from_q_new(sbe, ib, jb, ik)
    end do
  end do
end subroutine kr_rho_at

pure function rho_ij_from_q_new(sbe, ib, jb, ik) result(r)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  integer, intent(in) :: ib, jb, ik
  complex(8) :: r
  if (ib == jb) then
    r = sbe%qnm_new(ib, jb, ik)
  else
    r = sbe%exp_iphi(ib, jb, ik) * sbe%qnm_new(ib, jb, ik)
  end if
end function rho_ij_from_q_new


! One-layer FORWARD (+e_a, all three axes) halo plan for the roughness probe.
! Independent of the propagation halos: gicov_int has no resident rho halo at
! all, and the FD plan carries m_max shells this only needs one of.
subroutine kr_build_plan(sbe, gs, icomm)
  use salmon_global, only: num_kgrid
  use degenerate_block_ssbe, only: build_ik_neighbor, find_bvec, build_halo_lists
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  integer, intent(in) :: icomm
  integer :: irank, nproc, nk, o, ik, axis, iv, p, jj, iv_axis(3)
  integer, allocatable :: itbl_min(:), itbl_max(:), ikn(:, :)
  logical, allocatable :: needed_all(:, :)

  if (krp_built) return
  nk = sbe%nk
  call comm_get_groupinfo(icomm, irank, nproc)
  allocate(ikn(gs%nbvec, nk))
  call build_ik_neighbor(num_kgrid, gs%bvec, gs%nbvec, nk, ikn)
  iv_axis(1) = find_bvec(gs%bvec, gs%nbvec, 1, 0, 0)
  iv_axis(2) = find_bvec(gs%bvec, gs%nbvec, 0, 1, 0)
  iv_axis(3) = find_bvec(gs%bvec, gs%nbvec, 0, 0, 1)

  allocate(krp_nbr(3, nk))
  do ik = 1, nk
    do axis = 1, 3
      iv = iv_axis(axis)
      if (iv == 0 .or. num_kgrid(axis) < 2) then
        krp_nbr(axis, ik) = ik            ! no shift on this axis: S_a = identity
      else
        krp_nbr(axis, ik) = ikn(iv, ik)
      end if
    end do
  end do

  allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
  call split_range(1, nk, nproc, itbl_min, itbl_max)
  allocate(needed_all(nk, 0:nproc-1))
  needed_all = .false.
  do o = 0, nproc - 1
    if (itbl_max(o) < itbl_min(o)) cycle
    do ik = itbl_min(o), itbl_max(o)
      needed_all(ik, o) = .true.
      do axis = 1, 3
        needed_all(krp_nbr(axis, ik), o) = .true.
      end do
    end do
  end do
  call build_halo_lists(nk, nproc, irank, itbl_min, itbl_max, needed_all, &
                        krp_nsrc, krp_src_rank, krp_src_cnt, krp_src_ofs, krp_src_k, &
                        krp_ndst, krp_dst_rank, krp_dst_cnt, krp_dst_ofs, krp_dst_k)

  ! slot map: owned k first (ascending), then the received halo k
  allocate(krp_map(nk))
  krp_map(:) = 0
  krp_nx = 0
  do ik = sbe%ik_min, sbe%ik_max
    krp_nx = krp_nx + 1
    krp_map(ik) = krp_nx
  end do
  do p = 1, krp_nsrc
    do jj = 1, krp_src_cnt(p)
      ik = krp_src_k(krp_src_ofs(p) + jj)
      if (krp_map(ik) == 0) then
        krp_nx = krp_nx + 1
        krp_map(ik) = krp_nx
      end if
    end do
  end do
  do ik = sbe%ik_min, sbe%ik_max
    do axis = 1, 3
      if (krp_map(krp_nbr(axis, ik)) == 0) then
        write(*,*) "ERROR kr_build_plan: shifted k not covered, k =", ik, " axis =", axis
        error stop 1
      end if
    end do
  end do

  deallocate(itbl_min, itbl_max, needed_all, ikn)
  krp_built = .true.
end subroutine kr_build_plan


subroutine kr_exchange(nb, X, icomm)
  implicit none
  integer, intent(in) :: nb, icomm
  complex(8), intent(inout) :: X(nb, nb, *)
  integer :: p, jj, ik
  type(t_gh_buf), allocatable :: sb(:), rb(:)
  integer, allocatable :: reqs(:), reqr(:)

  allocate(rb(krp_nsrc), sb(krp_ndst))
  do p = 1, krp_nsrc
    allocate(rb(p)%b(nb, nb, krp_src_cnt(p)))
  end do
  do p = 1, krp_ndst
    allocate(sb(p)%b(nb, nb, krp_dst_cnt(p)))
  end do
  allocate(reqr(max(1, krp_nsrc)), reqs(max(1, krp_ndst)))
  do p = 1, krp_nsrc
    reqr(p) = comm_irecv(rb(p)%b, krp_src_rank(p), 3, icomm)
  end do
  do p = 1, krp_ndst
    do jj = 1, krp_dst_cnt(p)
      ik = krp_dst_k(krp_dst_ofs(p) + jj)
      sb(p)%b(:, :, jj) = X(:, :, krp_map(ik))
    end do
    reqs(p) = comm_isend(sb(p)%b, krp_dst_rank(p), 3, icomm)
  end do
  if (krp_nsrc > 0) then
    call comm_wait_all(reqr(1:krp_nsrc))
    do p = 1, krp_nsrc
      do jj = 1, krp_src_cnt(p)
        ik = krp_src_k(krp_src_ofs(p) + jj)
        X(:, :, krp_map(ik)) = rb(p)%b(:, :, jj)
      end do
    end do
  end if
  if (krp_ndst > 0) call comm_wait_all(reqs(1:krp_ndst))
  do p = 1, krp_nsrc
    deallocate(rb(p)%b)
  end do
  do p = 1, krp_ndst
    deallocate(sb(p)%b)
  end do
  deallocate(rb, sb, reqr, reqs)
end subroutine kr_exchange


! Output interval M from the environment (SBE_KROUGH=M); 0 = off (the default,
! so no run that does not ask for it changes by a single byte).
integer function sbe_krough_every() result(m)
  implicit none
  character(16) :: val
  integer :: ln, ios
  m = 0
  call get_environment_variable('SBE_KROUGH', val, ln)
  if (ln <= 0) return
  read(val, *, iostat=ios) m
  if (ios /= 0 .or. m < 0) m = 0
end function sbe_krough_every


! Midpoint-frozen exact-exponential step of the whole local k-slice.  q_mid =
! q_axis(t + dt/2) (reduced-mesh index units); the kernel evaluates x=kappa-a by
! degree-p Lagrange interpolation of the bounded cached H~ (node-exact when x
! lands on a shift) and steps each k-block exactly with the moving-gap-gated T2.
!
! The stencil is resolved ONCE per step (it depends only on q_mid, not on k) and
! reused across the whole k-loop.
subroutine dt_evolve_bloch_lg_integral(sbe, gs, q_mid, dt, icomm)
  use salmon_global, only: t_2, sbe_t2_gate_shape, sbe_t2_gate_theta, &
                            sbe_t2_gate_width, sbe_lg_degen_floor
  use gicov_integral_ssbe, only: gicov_int_stencil, gicov_int_interp_p, &
                                 gicov_int_step_k, gicov_int_nsten
  !$ use omp_lib
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: q_mid, dt
  integer, intent(in) :: icomm
  integer :: ik, nb, tid, ierr, nbad_l, nbad
  integer :: nodes(gicov_int_nsten), nsten
  real(8) :: wts(gicov_int_nsten), gamma

  if (.not. gi_built) then
    write(*, '(a)') "ERROR(dt_evolve_bloch_lg_integral): transport cache not built."
    error stop 1
  end if
  nb = sbe%nb
  gamma = 1d0 / t_2
  call gicov_int_stencil(q_mid, 1d0, gi_jmax, nodes, wts, nsten, ierr)
  if (ierr /= 0) then
    write(*, '(a,es16.8,a,i0)') "ERROR(dt_evolve_bloch_lg_integral): mesh shift q=", q_mid, &
      & " leaves the cached transport span +/-", gi_jmax
    error stop 1
  end if

  nbad_l = 0
  !$omp parallel do default(shared) private(ik, tid, ierr) reduction(+: nbad_l)
  do ik = gi_ikmin, gi_ikmax
    tid = 0
    !$ tid = omp_get_thread_num()
    call gicov_int_interp_p(gi_Ht(:, :, :, ik), nb, -gi_jmax, gi_jmax, &
                          & nodes, wts, nsten, gi_w_H(:, :, tid))
    call gicov_int_step_k(gi_w_H(:, :, tid), sbe%rho(:, :, ik), nb, dt, gamma, &
                        & sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_t2_gate_width, &
                        & sbe_lg_degen_floor, gi_w_eps(:, tid), gi_w_P(:, :, tid), &
                        & gi_w_R(:, :, tid), gi_w_cw(:, tid), gi_lcwork, gi_w_rw(:, tid), &
                        & gi_w_blk(:, tid), gi_w_rout(:, :, tid), ierr)
    if (ierr /= 0) then
      nbad_l = nbad_l + 1
    else
      sbe%rho(:, :, ik) = gi_w_rout(:, :, tid)
      ! X-full carries qnm == rho (exp_iphi=1), so mirror the co-moving density
      ! into BOTH qnm and qnm_new -- the length-gauge calc_trace reads qnm_new
      ! and the _sbe_diag herm/trace monitor reads qnm (the legacy AB4 path
      ! refreshes qnm itself at :1746; the integral path must do it here, else
      ! the monitor would report the t=0 density forever).
      sbe%qnm_new(:, :, ik) = gi_w_rout(:, :, tid)
      sbe%qnm(:, :, ik)     = gi_w_rout(:, :, tid)
    end if
  end do
  !$omp end parallel do

  ! COLLECTIVE fail: an instantaneous eigensolve that broke (LAPACK info /= 0 or
  ! non-finite eigenvalues) must abort the whole communicator, not be absorbed.
  ! Both former fallbacks were trace-preserving, so the Ne monitor could not see
  ! them -- silence here is precisely the failure mode this guard exists for.
  call comm_summation(nbad_l, nbad, icomm)
  if (nbad > 0) then
    write(*, '(a,i0,a)') "ERROR(dt_evolve_bloch_lg_integral): instantaneous eigensolve " // &
      & "failed on ", nbad, " k-block(s) (LAPACK info /= 0 or non-finite eigenvalues). " // &
      & "The co-moving Hamiltonian is corrupted; aborting rather than propagating it."
    error stop 1
  end if
end subroutine dt_evolve_bloch_lg_integral

! Integral-form current J_i = (1/Omega sum w) sum_kappa w Tr[rho~ v~_i], with
! v~_i the cached TRANSPORTED velocity interpolated to x = kappa - a(t).  No
! D.A subtraction (that is a velocity-gauge concept; the length gauge is Tr(v rho)).
subroutine calc_current_bloch_lg_integral(sbe, gs, q_now, jmat, icomm)
  use gicov_integral_ssbe, only: gicov_int_stencil, gicov_int_interp_p, &
                                 gicov_int_current_k, gicov_int_nsten
  !$ use omp_lib
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: q_now
  real(8), intent(out) :: jmat(3)
  integer, intent(in) :: icomm
  integer :: ik, nb, tid, i, ierr
  integer :: nodes(gicov_int_nsten), nsten
  real(8) :: wts(gicov_int_nsten), jk(3), tmp1(3), tmp(3)

  nb = sbe%nb
  call gicov_int_stencil(q_now, 1d0, gi_jmax, nodes, wts, nsten, ierr)
  if (ierr /= 0) then
    write(*, '(a)') "ERROR(calc_current_bloch_lg_integral): mesh shift leaves the cached span."
    error stop 1
  end if
  tmp1(:) = 0d0
  !$omp parallel do default(shared) private(ik, tid, i, jk) reduction(+:tmp1)
  do ik = gi_ikmin, gi_ikmax
    tid = 0
    !$ tid = omp_get_thread_num()
    do i = 1, 3
      call gicov_int_interp_p(gi_vt(:, :, :, i, ik), nb, -gi_jmax, gi_jmax, &
                            & nodes, wts, nsten, gi_w_v(:, :, i, tid))
    end do
    call gicov_int_current_k(sbe%rho(:, :, ik), gi_w_v(:, :, :, tid), nb, jk)
    tmp1(1:3) = tmp1(1:3) + gs%kweight(ik) * jk(1:3)
  end do
  !$omp end parallel do
  call comm_summation(tmp1, tmp, 3, icomm)
  jmat(1:3) = tmp(1:3) / sum(gs%kweight(:)) / gs%volume
end subroutine calc_current_bloch_lg_integral

! Band-resolved EXCITATION at x = kappa - a(t) (yn_sbe_out_occ='y'), the
! integral-transport counterpart of calc_band_population:
!
!   dn_a(x) = ( P^dag [ rho~ - F~0(x) ] P )_aa ,   H~(x) = P diag(eps) P^dag
!
! summed over the local k-slice and reduced.  TWO things are essential here and
! both were wrong before:
!
!  (1) the projector must be the INSTANTANEOUS eigenbasis of the interpolated
!      H~ (diag(rho~) in the frozen band basis is not transport-invariant), and
!  (2) the reference must be the co-moving GS occupation F~0(x) -- the same
!      gs%occup, transported by the same W and interpolated to the same x.
!      Reporting the ABSOLUTE population instead (the previous behaviour) makes
!      an unexcited state read out as fully occupied, i.e. it is not an
!      excitation at all; and a rigid "1 below nvb" reference would manufacture
!      a spurious excitation for a metal, whose f0 is fractional and k-varying.
!
! Degenerate members are reported as a block sum (gicov_int_occupation_k), since
! the individual eigen-columns inside a degenerate manifold are arbitrary.
subroutine calc_band_population_integral(sbe, gs, q_now, nex_b, icomm)
  use salmon_global, only: sbe_lg_degen_floor
  use gicov_integral_ssbe, only: gicov_int_stencil, gicov_int_interp_p, &
                                 gicov_int_occupation_k, gicov_int_nsten
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: q_now
  real(8), intent(out) :: nex_b(1:sbe%nb)
  integer, intent(in) :: icomm
  integer :: ik, nb, ib, ierr, nbad_l, nbad
  integer :: nodes(gicov_int_nsten), nsten
  real(8) :: wts(gicov_int_nsten)
  real(8), allocatable :: acc(:), dn(:)
  complex(8), allocatable :: F0(:, :), D(:, :)

  nb = sbe%nb
  call gicov_int_stencil(q_now, 1d0, gi_jmax, nodes, wts, nsten, ierr)
  if (ierr /= 0) then
    write(*, '(a)') "ERROR(calc_band_population_integral): mesh shift leaves the cached span."
    error stop 1
  end if
  ! serial over k (periodic diagnostic, not the propagation hot path): reuses
  ! the thread-0 work slot; heap work allocated OUTSIDE any OMP region.
  allocate(acc(nb), dn(nb), F0(nb, nb), D(nb, nb));  acc(:) = 0d0
  nbad_l = 0
  do ik = gi_ikmin, gi_ikmax
    call gicov_int_interp_p(gi_Ht(:, :, :, ik), nb, -gi_jmax, gi_jmax, &
                          & nodes, wts, nsten, gi_w_H(:, :, 0))
    ! co-moving GS reference at the SAME x, interpolated the SAME way as H~
    call gicov_int_interp_p(gi_Ft(:, :, :, ik), nb, -gi_jmax, gi_jmax, &
                          & nodes, wts, nsten, F0)
    D(:, :) = sbe%rho(:, :, ik) - F0(:, :)      ! Hermitian: the excitation
    call gicov_int_occupation_k(gi_w_H(:, :, 0), D, nb, sbe_lg_degen_floor, &
                              & gi_w_eps(:, 0), gi_w_P(:, :, 0), gi_w_R(:, :, 0), &
                              & gi_w_cw(:, 0), gi_lcwork, gi_w_rw(:, 0), dn, &
                              & gi_w_blk(:, 0), ierr)
    if (ierr /= 0) then
      nbad_l = nbad_l + 1
      cycle
    end if
    do ib = 1, nb
      acc(ib) = acc(ib) + gs%kweight(ik) * dn(ib)
    end do
  end do
  ! COLLECTIVE fail (same contract as the propagator): a broken instantaneous
  ! eigenproblem must abort, never fall back to the transport-non-invariant
  ! diag(rho~) -- that fallback is trace-correct and therefore undetectable.
  call comm_summation(nbad_l, nbad, icomm)
  if (nbad > 0) then
    write(*, '(a,i0,a)') "ERROR(calc_band_population_integral): instantaneous eigensolve " // &
      & "failed on ", nbad, " k-block(s) (LAPACK info /= 0 or non-finite eigenvalues)."
    error stop 1
  end if
  call comm_summation(acc, nex_b, nb, icomm)
  nex_b(1:nb) = nex_b(1:nb) / sum(gs%kweight(:))
  deallocate(acc, dn, F0, D)
end subroutine calc_band_population_integral


! Transport-invariant VALENCE trace at x = kappa - a(t), for _sbe_nex.data.
!
!   N_v(x) = sum_kappa w_kappa sum_{a=1..nvb} ( P^dag rho~ P )_aa / sum w
!
! The frozen-basis partial trace it replaces -- calc_trace(sbe, gs, gs%nvb),
! i.e. the sum of the leading nvb DIAGONAL entries of rho~ -- is NOT transport
! invariant.  A Wilson rotation preserves the FULL trace of rho~ but not the
! trace of its leading principal block, so a pure transport (no field-driven
! excitation whatsoever) already moves weight across the nvb cut and is reported
! as excitation.  Projecting onto the INSTANTANEOUS eigenstates and summing the
! lowest nvb of them is the invariant statement; the full trace (n_electrons)
! stays valid as a plain Tr rho~ and keeps using calc_trace.
!
! The eigenvalues come back ascending from zheev, so "the lowest nvb
! instantaneous states" is the leading nvb of the projected diagonal.  Values
! are block-symmetrized (gicov_int_occupation_k), so a degenerate manifold
! straddling the nvb cut contributes its invariant proportional share rather
! than an arbitrary basis-dependent split.
function calc_valence_trace_integral(sbe, gs, q_now, nvb, icomm) result(trv)
  use salmon_global, only: sbe_lg_degen_floor
  use gicov_integral_ssbe, only: gicov_int_stencil, gicov_int_interp_p, &
                                 gicov_int_occupation_k, gicov_int_nsten
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: q_now
  integer, intent(in) :: nvb, icomm
  real(8) :: trv
  integer :: ik, nb, ib, ierr, nbad_l, nbad
  integer :: nodes(gicov_int_nsten), nsten
  real(8) :: wts(gicov_int_nsten), acc, tmp
  real(8), allocatable :: nocc(:)

  nb = sbe%nb
  call gicov_int_stencil(q_now, 1d0, gi_jmax, nodes, wts, nsten, ierr)
  if (ierr /= 0) then
    write(*, '(a)') "ERROR(calc_valence_trace_integral): mesh shift leaves the cached span."
    error stop 1
  end if
  allocate(nocc(nb))
  acc = 0d0
  nbad_l = 0
  do ik = gi_ikmin, gi_ikmax
    call gicov_int_interp_p(gi_Ht(:, :, :, ik), nb, -gi_jmax, gi_jmax, &
                          & nodes, wts, nsten, gi_w_H(:, :, 0))
    call gicov_int_occupation_k(gi_w_H(:, :, 0), sbe%rho(:, :, ik), nb, sbe_lg_degen_floor, &
                              & gi_w_eps(:, 0), gi_w_P(:, :, 0), gi_w_R(:, :, 0), &
                              & gi_w_cw(:, 0), gi_lcwork, gi_w_rw(:, 0), nocc, &
                              & gi_w_blk(:, 0), ierr)
    if (ierr /= 0) then
      nbad_l = nbad_l + 1
      cycle
    end if
    do ib = 1, min(nvb, nb)
      acc = acc + gs%kweight(ik) * nocc(ib)
    end do
  end do
  call comm_summation(nbad_l, nbad, icomm)
  if (nbad > 0) then
    write(*, '(a,i0,a)') "ERROR(calc_valence_trace_integral): instantaneous eigensolve " // &
      & "failed on ", nbad, " k-block(s) (LAPACK info /= 0 or non-finite eigenvalues)."
    error stop 1
  end if
  call comm_summation(acc, tmp, icomm)
  trv = tmp / sum(gs%kweight(:))
  deallocate(nocc)
end function calc_valence_trace_integral

subroutine calc_current_bloch_lg(sbe, gs, jmat, icomm)
    use salmon_global, only: sbe_lg_degen
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(out) :: jmat(3)
    integer, intent(in) :: icomm
    integer :: ik, ib, jb, jj
    integer :: nk, nb
    complex(8) :: tmp1(3),tmp(3)
    complex(8) :: rho_ji

    nk = sbe%nk
    nb = sbe%nb

    tmp1(:) = 0.d0

    ! Pb4 (GI current): when sbe_lg_degen=='gi'/'gifix'/'gicov' the Pb3 xi substitution
    ! (gi/gifix) or the same-block exp_iphi=1 / covariant-transport representation (gicov)
    ! makes the dnm_i / exp_iphi phase machinery inconsistent inside degenerate blocks, so
    ! the delta_omega*|dnm|*aimag(qnm) current below is unreliable there.  Instead reconstruct
    ! the physical density matrix rho from the propagated qnm (rho = exp_iphi*qnm off the
    ! diagonal, rho = qnm on the diagonal; inverse of prepare_qnm's qnm=conjg(exp_iphi)*rho,
    ! == rho_ij_from_q) and contract it with the velocity v = p_tm_matrix (+ rvnl_tm_matrix
    ! when flag_vnl_correction), using the exact index/sign convention of the VG current
    ! calc_current_bloch (paramagnetic part; the length gauge carries no Ac*N term).  This
    ! reconstruction is Tr(v rho), gauge-invariant under a consistent GS gauge, so it is the
    ! minimal correct gicov current (subsumes Task 6's in-block current reconstruction).
    ! gi/gifix/gicov reconstruct Tr(v rho) at the frozen k.  gicov_int does NOT
    ! use this branch: its current is Tr[rho~ v~] with the TRANSPORTED velocity
    ! v~(kappa-a) from the integral-form cache (calc_current_bloch_lg_integral),
    ! so it is excluded here (uses_prod_dk minus the integral mode).
    if (uses_prod_dk(sbe_lg_degen) .and. .not. uses_integral_gicov(sbe_lg_degen)) then
        !$omp parallel do default(shared) private(ik,jj,ib,jb,rho_ji) reduction(+:tmp1)
        do ik = sbe%ik_min, sbe%ik_max
        do jj = 1,3
        do ib = 1,nb
        do jb = 1,nb
            if (ib == jb) then
                rho_ji = sbe%qnm_new(ib, ib, ik)
            else
                rho_ji = sbe%exp_iphi(jb, ib, ik) * sbe%qnm_new(jb, ib, ik)
            end if
            tmp1(jj) = tmp1(jj) + gs%kweight(ik) * rho_ji * gs%p_tm_matrix(ib, jb, jj, ik)
            if (sbe%flag_vnl_correction) then
                tmp1(jj) = tmp1(jj) + gs%kweight(ik) * rho_ji * gs%rvnl_tm_matrix(ib, jb, jj, ik)
            end if
        end do
        end do
        end do
        end do
        !$omp end parallel do

        call comm_summation(tmp1, tmp, 3, icomm)
        jmat(:) = dble(tmp(:)) / sum(gs%kweight(:)) / gs%volume
        return
    end if

    !$omp parallel do default(shared) private(ik,ib,jb,jj) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
    do jj = 1,3
    do ib = 1,nb
    do jb = 1,nb
        if(ib == jb) then
            tmp1(jj) = tmp1(jj) + gs%grad_k_eigen(ib, jj, ik) *  &
                               & dble(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
        else
            tmp1(jj) = tmp1(jj) + gs%delta_omega(ib, jb, ik) * &
                               & abs(sbe%dnm_i(ib, jb, jj, ik)) * &
                               & aimag(sbe%qnm_new(jb, ib, ik)) * gs%kweight(ik)
        end if
    end do
    end do
    end do
    end do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = -dble(tmp(:)) / sum(gs%kweight(:)) / gs%volume

end subroutine calc_current_bloch_lg

subroutine adams_moulton_coefs(bj_am)
  implicit none
  real(8) :: bj_am(8,8)

  bj_am(1,1) = 1.d0

  !bj_am(1,2) = 3.d0/2.d0
  !bj_am(2,2) = -1.d0/2.d0

  !bj_am(1,3) = 23.d0/12.d0
  !bj_am(2,3) = -4.d0/3.d0
  !bj_am(3,3) = 5.d0/12.d0

  bj_am(1,4) = 55.d0/24.d0
  bj_am(2,4) = -59.d0/24.d0
  bj_am(3,4) = 37.d0/24.d0
  bj_am(4,4) = -3.d0/8.d0

  !bj_am(1,5) = 1901.d0/720.d0
  !bj_am(2,5) = -1387.d0/360.d0
  !bj_am(3,5) = 109.d0/30.d0
  !bj_am(4,5) = -637.d0/360.d0
  !bj_am(5,5) = 251.d0/720.d0

  !bj_am(1,6) = 4277.d0/1440.d0
  !bj_am(2,6) = -2641.d0/480.d0
  !bj_am(3,6) = 4991.d0/720.d0
  !bj_am(4,6) = -3649.d0/720.d0
  !bj_am(5,6) = 959.d0/480.d0
  !bj_am(6,6) = -95.d0/288.d0

  !bj_am(1,7) = 198721.d0/60480.d0
  !bj_am(2,7) = -18637.d0/2520.d0
  !bj_am(3,7) = 235183.d0/20160.d0
  !bj_am(4,7) = -10754.d0/945.d0
  !bj_am(5,7) = 135713.d0/20160.d0
  !bj_am(6,7) = -5603.d0/2520.d0
  !bj_am(7,7) = 19087.d0/60480.d0

  bj_am(1,8) = 16083.d0/4480.d0
  bj_am(2,8) = -1152169.d0/120960.d0
  bj_am(3,8) = 242653.d0/13440.d0
  bj_am(4,8) = -296053.d0/13440.d0
  bj_am(5,8) = 2102243.d0/120960.d0
  bj_am(6,8) = -115747.d0/13440.d0
  bj_am(7,8) = 32863.d0/13440.d0
  bj_am(8,8) = -5257.d0/17280.d0

end subroutine adams_moulton_coefs


end module

