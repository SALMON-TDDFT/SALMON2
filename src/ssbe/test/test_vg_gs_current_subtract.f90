! src/ssbe/test/test_vg_gs_current_subtract.f90
! Velocity-gauge window f-sum-deficiency current subtraction
! (yn_sbe_gs_current_subtract, v2): analytic-fixture gate.
!
! WHY THIS FILE EXISTS
!   Truncating the VG-SBE basis at nstate_sbe breaks the f-sum rule, so the
!   readout J = Tr[rho(t) v(k+A)] carries a spurious current proportional to
!   A(t) (a filled-band insulator in a COMPLETE basis carries none; observed as
!   the non-converging N_b-ladder of peak |J| in the gw campaign).  v1 (commit
!   cda60ee4) subtracted the whole frozen-GS current Tr[rho0 v(k+A)] = A*Ne/V,
!   which removes the FULL diamagnetic term -- including the part the retained
!   window response legitimately cancels -- and overcorrected on real Si@nb=32
!   (peak |Jz| 2.05e-4 -> 9.41e-4, away from the LG oracle; JID49460686).
!   v2 subtracts only what the window is MISSING: the f-sum deficiency tensor
!
!     D_ij = ( Ne d_ij - < sum_{v occ} f_v sum_{m unocc}
!              2 Re[ pi^i_vm pi^j_mv ] / (eps_m - eps_v) >_k ) / V,
!     J(t) -= D * A(t)
!
!   (f = gs%occup, spinless convention occup=2; Ne = <sum f>_k = Tr[rho0];
!   pi = p_tm (+ rvnl iff vnl flag); <>_k = kweight average; see
!   build_fsum_deficiency_tensor).  This test pins the defining properties on
!   small analytic fixtures plus one independent evaluation of the legacy
!   formula:
!
!     0. OFF path unchanged: flag='n' reproduces an INDEPENDENTLY coded
!        evaluation of the legacy norder=0 formula
!        J = (Re sum_k w tr[rho p] / sum w + A * tr[rho]_w) / V   (refactor
!        witness: the wrapper edit did not move the legacy math).
!     a. COMPLETE basis => D == 0 to machine precision: fixture saturating the
!        Thomas-Reiche-Kuhn sum per occupied band exactly by construction
!        (each Cartesian direction coupled to its own partner band with
!        2|p|^2/(eps_m - eps_v) = 1).
!     b. TRUNCATED window => D_zz > 0 with the hand value (2 - 2*gamma^2)/V,
!        and the on/off readout difference equals -D*A(t) pointwise for every
!        A in a sweep x norder 0..3 x vnl on/off, on a generic excited rho
!        (populations + coherences); D built through the PUBLIC init path.
!     c. A == 0 => the correction vanishes IDENTICALLY (J_on == J_off
!        bitwise): post-pulse observables are untouched by the flag.
!     d. Analytic 2-band fixture: all 9 components of D match the hand
!        formula D_ij = (2 d_ij - 4 Re[pi^i_12 pi^j_21]/Delta)/V, for both
!        vnl off (pi = p) and vnl on (pi = p + rvnl); D symmetric (TR grid).
!     e. Guards: degenerate occupied<->occupied pairs are excluded (no NaN,
!        no contribution) and occupied x unoccupied pairs closer than the
!        theta_off floor are dropped, not amplified.
!
!   k-grids are exact time-reversal pairs (p(k2) = -conjg(p(k1)), same
!   eigenvalues), as for any real material's full-BZ grid.
!
! BUILD (already-built ninja tree at build_local/; single-process communication
! dummy).  Compile from a CLEAN dir (a stale .mod file in the repo root shadows
! the fresh build_local one -- use -I build_local ONLY).  Link the SAME objects
! the salmon executable built, minus main.f90.o, via an @objs.txt response
! file:
!
!   find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o' > objs.txt
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<clean_dir> \
!     -c <repo>/src/ssbe/test/test_vg_gs_current_subtract.f90 -o <clean_dir>/test_vg_gs_current_subtract.o
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     @objs.txt <clean_dir>/test_vg_gs_current_subtract.o -o <clean_dir>/test_vg_gs_current_subtract \
!     -framework Accelerate -lm -ldl
!   <clean_dir>/test_vg_gs_current_subtract
!
program test_vg_gs_current_subtract
  use gs_info_ssbe,      only: s_sbe_gs_info, sbe_gs_set_replicated_kmap
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                               calc_current_bloch, build_fsum_deficiency_tensor
  use degenerate_block_ssbe, only: theta_off
  use salmon_global,     only: gauge_sbe, norder_correction, &
                               yn_sbe_gs_current_subtract, yn_vnl_correction
  implicit none

  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: vol = 150.0d0
  ! A-sweep: zero, small/large along z, oblique, negative -- the subtraction
  ! must hold pointwise in A, not just for one field.
  integer, parameter :: nA = 5
  real(8), parameter :: Aset(3, nA) = reshape( (/ &
    &  0.00d0, 0.00d0, 0.00d0, &
    &  0.00d0, 0.00d0, 0.05d0, &
    &  0.00d0, 0.00d0, 0.80d0, &
    &  0.30d0, -0.20d0, 0.55d0, &
    & -0.15d0, 0.10d0, -0.40d0 /), (/3, nA/) )
  ! truncation factor of the z-channel in the windowed fixture (check b)
  real(8), parameter :: gamma_z = 0.6d0
  integer :: nfail

  nfail = 0
  call set_globals()

  call test_off_path_analytic(nfail)          ! check 0 (refactor witness)
  call test_complete_basis_d_zero(nfail)      ! check a
  call test_truncated_window(nfail)           ! check b
  call test_a0_identity(nfail)                ! check c
  call test_two_band_analytic(nfail)          ! check d
  call test_guards(nfail)                     ! check e

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_vg_gs_current_subtract checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    gauge_sbe = 'velocity_gauge'      ! calc_current_bloch is the VG readout
    norder_correction = 0             ! per-check loops override
    yn_sbe_gs_current_subtract = 'n'  ! per-check code overrides
    yn_vnl_correction = 'n'           ! init-path D build: vnl off default
  end subroutine set_globals

  !======================= assert helper (test_gicov_* style) =================
  subroutine check_true(cond, label, nfail)
    implicit none
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (cond) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a)') "FAIL  " // label
      nfail = nfail + 1
    end if
  end subroutine check_true

  !======================= gs fixture allocator ===============================
  subroutine alloc_gs(gs, nb, nk, ne)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    integer, intent(in) :: nb, nk, ne
    gs%nk = nk
    gs%nb = nb
    call sbe_gs_set_replicated_kmap(gs, nk)   ! replicated k layout (kmap = identity)
    gs%ne = ne
    gs%volume = vol
    allocate(gs%kweight(nk), gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_tm_matrix(nb, nb, 3, nk), gs%rvnl_tm_matrix(nb, nb, 3, nk))
    gs%kweight(:) = 1.0d0
    gs%occup(:, :) = 0d0
    gs%occup(1:ne/2, :) = 2d0
    gs%p_tm_matrix(:, :, :, :) = (0d0, 0d0)
    gs%rvnl_tm_matrix(:, :, :, :) = (0d0, 0d0)
  end subroutine alloc_gs

  subroutine finish_gs(gs)
    ! delta_omega from eigen; k2 = exact TR partner of k1 for every operator
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    integer :: ik, ib, jb, idir
    do ik = 1, gs%nk
      do jb = 1, gs%nb
        do ib = 1, gs%nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do
    end do
    if (gs%nk == 2) then
      do idir = 1, 3
        gs%p_tm_matrix(:, :, idir, 2)    = -conjg(gs%p_tm_matrix(:, :, idir, 1))
        gs%rvnl_tm_matrix(:, :, idir, 2) = -conjg(gs%rvnl_tm_matrix(:, :, idir, 1))
      end do
    end if
  end subroutine finish_gs

  subroutine free_gs(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    deallocate(gs%kweight, gs%eigen, gs%occup, gs%delta_omega, &
             & gs%p_tm_matrix, gs%rvnl_tm_matrix)
  end subroutine free_gs

  !======================= TRK-saturating fixture (checks a/b) ================
  ! nb=4, nk=2 (TR pair), band 1 occupied (occup=2, ne=2).  Each Cartesian
  ! direction d couples band 1 to its OWN partner band 1+d with
  ! |p^d_{1,1+d}|^2 = (eps_{1+d} - eps_1)/2, i.e. 2|p|^2/Delta = 1 exactly:
  ! the one-electron TRK sum is saturated per direction, off-diagonal D
  ! channels vanish (no shared partners), so D == 0 to machine precision.
  ! scale_z multiplies the z-channel element: sum captured = scale_z^2 of the
  ! rule => D_zz = (2 - 2*scale_z^2)/V (hand value used by check b).
  ! Nonzero band-diagonal p entries are included (they must not enter D).
  subroutine build_gs_trk(gs, scale_z)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: scale_z
    integer, parameter :: nb = 4, nk = 2, ne = 2
    real(8) :: delta_d, c, re, im
    complex(8) :: phase, h(nb, nb)
    integer :: d, ib, jb

    call alloc_gs(gs, nb, nk, ne)
    gs%eigen(1, :) = -0.30d0
    gs%eigen(2, :) =  0.20d0
    gs%eigen(3, :) =  0.45d0
    gs%eigen(4, :) =  0.70d0

    do d = 1, 3
      delta_d = gs%eigen(1 + d, 1) - gs%eigen(1, 1)
      c = sqrt(delta_d / 2d0)
      if (d == 3) c = c * scale_z
      phase = exp(zi_ * (0.4d0 * dble(d)))    ! nontrivial complex phase
      gs%p_tm_matrix(1, 1 + d, d, 1) = c * phase
      gs%p_tm_matrix(1 + d, 1, d, 1) = c * conjg(phase)
      ! band-diagonal (intraband) elements: irrelevant for D by construction
      do ib = 1, nb
        gs%p_tm_matrix(ib, ib, d, 1) = dcmplx(0.03d0 * dble(ib * d), 0d0)
      end do
    end do

    ! Hermitian rvnl (nonlocal-velocity correction): NOT part of the TRK
    ! saturation above -- checks a/b build D with flag_vnl = .false., so this
    ! only feeds the vnl branch of the READOUT core (check b's vnl loop).
    do d = 1, 3
      do jb = 1, nb
        do ib = 1, nb
          re = 0.017d0 * dble(ib + jb) - 0.009d0 * dble(d) &
             & + 0.004d0 * dble(ib * jb)
          im = 0.006d0 * dble(ib - jb) * dble(4 - d)
          h(ib, jb) = dcmplx(re, im)
        end do
      end do
      gs%rvnl_tm_matrix(:, :, d, 1) = 0.5d0 * (h + transpose(conjg(h)))
    end do
    call finish_gs(gs)
  end subroutine build_gs_trk

  !======================= solver fixtures ====================================
  ! Ground-state solver through the PUBLIC path: init_sbe_bloch_solver sets
  ! rho = diag(gs%occup) and, when yn_sbe_gs_current_subtract=='y' at init
  ! time, builds the deficiency tensor sbe%fsum_D.
  subroutine make_solver_gs(sbe, gs, vnl)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    logical, intent(in) :: vnl
    integer :: icomm
    icomm = 0
    call init_sbe_bloch_solver(sbe, gs, gs%nb, icomm)
    sbe%flag_vnl_correction = vnl
  end subroutine make_solver_gs

  ! Generic excited-state density: Hermitian, populations perturbed away from
  ! occup (trace NOT preserved on purpose) plus dense coherences, different
  ! at the two k (no TR symmetry imposed on rho: none is required).
  subroutine make_solver_excited(sbe, gs, vnl)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    logical, intent(in) :: vnl
    integer :: ik, ib, jb
    complex(8) :: c
    call make_solver_gs(sbe, gs, vnl)
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, gs%nb
        sbe%rho(ib, ib, ik) = gs%occup(ib, ik) &
          & + 0.17d0 * dble(ik) / dble(ib) - 0.11d0 * dble(ib - 2)
      end do
      do jb = 1, gs%nb
        do ib = jb + 1, gs%nb
          c = dcmplx(0.13d0 * dble(ib) - 0.05d0 * dble(jb) + 0.02d0 * dble(ik), &
                   & 0.09d0 * dble(ib - jb) - 0.03d0 * dble(ik))
          sbe%rho(ib, jb, ik) = c
          sbe%rho(jb, ib, ik) = conjg(c)
        end do
      end do
    end do
  end subroutine make_solver_excited

  !======================= independent legacy formula (norder=0) ==============
  ! Deliberately re-coded from the physics (NOT from calc_current_bloch_core):
  ! J = ( Re sum_k w tr[rho p] / sum_k w  +  A * (sum_k w tr[rho] / sum_k w) ) / V
  subroutine current_reference_norder0(sbe, gs, Ac, jref)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(3)
    real(8), intent(out) :: jref(3)
    integer :: ik, ib, jb, idir
    complex(8) :: s(3), pmat
    real(8) :: tr
    s(:) = (0d0, 0d0)
    tr = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, gs%nb
        tr = tr + dble(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        do jb = 1, gs%nb
          do idir = 1, 3
            pmat = gs%p_tm_matrix(ib, jb, idir, ik)
            if (sbe%flag_vnl_correction) pmat = pmat + gs%rvnl_tm_matrix(ib, jb, idir, ik)
            s(idir) = s(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * pmat
          end do
        end do
      end do
    end do
    jref(:) = (dble(s(:)) / sum(gs%kweight) + Ac(:) * (tr / sum(gs%kweight))) / vol
  end subroutine current_reference_norder0

  !======================= check 0: OFF path analytic =========================
  subroutine test_off_path_analytic(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: jm(3), jref(3), errmax
    integer :: ia, icomm, ivnl

    icomm = 0
    call build_gs_trk(gs, gamma_z)
    errmax = 0d0
    norder_correction = 0
    yn_sbe_gs_current_subtract = 'n'
    do ivnl = 0, 1
      call make_solver_excited(sbe, gs, ivnl == 1)
      do ia = 1, nA
        call calc_current_bloch(sbe, gs, Aset(:, ia), jm, icomm)
        call current_reference_norder0(sbe, gs, Aset(:, ia), jref)
        errmax = max(errmax, maxval(abs(jm - jref)))
      end do
      deallocate(sbe%rho)
    end do
    call free_gs(gs)
    call check_true(errmax < 1d-13, "off-path == independent legacy formula (norder=0)", nfail)
  end subroutine test_off_path_analytic

  !======================= check a: complete basis => D == 0 ==================
  subroutine test_complete_basis_d_zero(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: dmax

    call build_gs_trk(gs, 1.0d0)              ! scale_z = 1: TRK fully saturated
    call make_solver_gs(sbe, gs, .false.)
    call build_fsum_deficiency_tensor(sbe, gs, .false., 0)
    dmax = maxval(abs(sbe%fsum_D))
    deallocate(sbe%rho)
    call free_gs(gs)
    ! exact cancellation Ne - <sum f_v * 1> up to sqrt/square roundoff (~ulp)
    call check_true(dmax < 5d-15, "complete basis (TRK saturated): D == 0 to machine precision", nfail)
  end subroutine test_complete_basis_d_zero

  !======================= check b: truncated window ==========================
  subroutine test_truncated_window(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: j_on(3), j_off(3), dref(3, 3), errd, errmax, submax
    integer :: ia, no, icomm, ivnl, i

    icomm = 0
    call build_gs_trk(gs, gamma_z)

    ! D through the PUBLIC init path (flag set at init time)
    yn_sbe_gs_current_subtract = 'y'
    call make_solver_gs(sbe, gs, .false.)
    yn_sbe_gs_current_subtract = 'n'
    call check_true(sbe%fsum_D_built, "init path builds D when the flag is on at init", nfail)

    ! hand values: x/y channels saturated => 0; z channel captures gamma^2 of
    ! the rule => D_zz = (2 - 2*gamma^2)/V; all off-diagonals 0.
    dref(:, :) = 0d0
    dref(3, 3) = (2d0 - 2d0 * gamma_z**2) / vol
    errd = maxval(abs(sbe%fsum_D - dref))
    call check_true(errd < 5d-15, "truncated window: D == hand value ((2-2g^2)/V on zz)", nfail)
    call check_true(sbe%fsum_D(3, 3) > 1d-3, "truncated window: D_zz > 0 (deficiency present)", nfail)
    deallocate(sbe%rho)

    ! on/off difference == -D*A pointwise: generic excited rho, A-sweep x
    ! norder 0..3 x vnl on/off.  D is built at init (flag on), then the flag
    ! is toggled around the readout calls only.
    errmax = 0d0
    submax = 0d0
    do ivnl = 0, 1
      yn_sbe_gs_current_subtract = 'y'
      call make_solver_excited(sbe, gs, ivnl == 1)   ! init under 'y': D built
      do no = 0, 3
        norder_correction = no
        do ia = 1, nA
          yn_sbe_gs_current_subtract = 'n'
          call calc_current_bloch(sbe, gs, Aset(:, ia), j_off, icomm)
          yn_sbe_gs_current_subtract = 'y'
          call calc_current_bloch(sbe, gs, Aset(:, ia), j_on, icomm)
          errmax = max(errmax, maxval(abs((j_on - j_off) &
            & + matmul(sbe%fsum_D, Aset(:, ia)))))
          ! non-vacuous: for A /= 0 the subtracted term must be well above
          ! roundoff, or the identity above is the trivial 0 == 0.
          if (ia > 1) submax = max(submax, maxval(abs(matmul(sbe%fsum_D, Aset(:, ia)))))
        end do
      end do
      deallocate(sbe%rho)
    end do
    norder_correction = 0
    yn_sbe_gs_current_subtract = 'n'
    call free_gs(gs)
    call check_true(errmax < 1d-15, "flag on/off difference == -D*A (all A/norder/vnl)", nfail)
    call check_true(submax > 1d-4, "non-vacuous: |D*A| well above roundoff for A/=0", nfail)
  end subroutine test_truncated_window

  !======================= check c: A=0 => correction == 0 bitwise ============
  subroutine test_a0_identity(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: j_on(3), j_off(3), a0(3), errmax
    integer :: no, icomm, ivnl

    icomm = 0
    a0(:) = 0d0
    call build_gs_trk(gs, gamma_z)
    errmax = 0d0
    do ivnl = 0, 1
      yn_sbe_gs_current_subtract = 'y'
      call make_solver_excited(sbe, gs, ivnl == 1)   ! init under 'y': D built
      do no = 0, 3
        norder_correction = no
        yn_sbe_gs_current_subtract = 'n'
        call calc_current_bloch(sbe, gs, a0, j_off, icomm)
        yn_sbe_gs_current_subtract = 'y'
        call calc_current_bloch(sbe, gs, a0, j_on, icomm)
        errmax = max(errmax, maxval(abs(j_on - j_off)))
      end do
      deallocate(sbe%rho)
    end do
    norder_correction = 0
    yn_sbe_gs_current_subtract = 'n'
    call free_gs(gs)
    ! D*0 == 0 exactly and x - 0 == x bitwise: the flag is an exact no-op at
    ! A=0 (stronger than the v1 roundoff bound).
    call check_true(errmax == 0d0, "A=0: correction is identically zero (bitwise)", nfail)
  end subroutine test_a0_identity

  !======================= check d: analytic 2-band ===========================
  subroutine test_two_band_analytic(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    complex(8) :: p12(3), r12(3), q12(3)
    real(8) :: delta, dref(3, 3), errd(0:1), symerr
    integer :: i, j, ivnl

    call alloc_gs(gs, 2, 2, 2)                 ! nb=2, nk=2 (TR pair), band 1 occ
    gs%eigen(1, :) = -0.20d0
    gs%eigen(2, :) =  0.30d0
    delta = gs%eigen(2, 1) - gs%eigen(1, 1)

    p12 = (/ dcmplx(0.31d0,  0.12d0), dcmplx(-0.07d0, 0.22d0), dcmplx(0.40d0, -0.05d0) /)
    r12 = (/ dcmplx(0.04d0, -0.02d0), dcmplx( 0.06d0, 0.01d0), dcmplx(-0.03d0, 0.05d0) /)
    do i = 1, 3
      gs%p_tm_matrix(1, 2, i, 1) = p12(i)
      gs%p_tm_matrix(2, 1, i, 1) = conjg(p12(i))
      gs%rvnl_tm_matrix(1, 2, i, 1) = r12(i)
      gs%rvnl_tm_matrix(2, 1, i, 1) = conjg(r12(i))
    end do
    call finish_gs(gs)

    call make_solver_gs(sbe, gs, .false.)
    do ivnl = 0, 1
      call build_fsum_deficiency_tensor(sbe, gs, ivnl == 1, 0)
      ! hand formula: D_ij = (2 d_ij - f * 2 Re[pi^i_12 pi^j_21] / Delta)/V,
      ! f = 2 (both TR k carry the same real product, so the k-average drops).
      q12 = p12
      if (ivnl == 1) q12 = p12 + r12
      do j = 1, 3
        do i = 1, 3
          dref(i, j) = -2d0 * 2d0 * dble(q12(i) * conjg(q12(j))) / delta / vol
        end do
        dref(j, j) = dref(j, j) + 2d0 / vol
      end do
      errd(ivnl) = maxval(abs(sbe%fsum_D - dref))
      symerr = maxval(abs(sbe%fsum_D - transpose(sbe%fsum_D)))
      call check_true(symerr < 1d-16, "2-band: D symmetric on the TR grid", nfail)
      ! production-style PUBLIC init wiring (codex review): re-init the SAME
      ! solver with yn_vnl_correction matching ivnl -- init must pick the vnl
      ! choice up from the global and reproduce the explicit build exactly.
      deallocate(sbe%rho)
      yn_sbe_gs_current_subtract = 'y'
      if (ivnl == 1) yn_vnl_correction = 'y'
      call make_solver_gs(sbe, gs, ivnl == 1)
      yn_vnl_correction = 'n'
      yn_sbe_gs_current_subtract = 'n'
      errd(ivnl) = max(errd(ivnl), maxval(abs(sbe%fsum_D - dref)))
    end do
    ! stale-state regression (codex review): a flag-off re-init of a solver
    ! that HAD a flag-on tensor must clear it (no stale D past the readout's
    ! fail-closed gate after an off->on toggle).
    deallocate(sbe%rho)
    call make_solver_gs(sbe, gs, .false.)     ! yn flag is 'n' here
    call check_true(.not. sbe%fsum_D_built .and. maxval(abs(sbe%fsum_D)) == 0d0, &
                  & "flag-off re-init clears a previously built D (no stale state)", nfail)
    deallocate(sbe%rho)
    call free_gs(gs)
    call check_true(errd(0) < 1d-15, "2-band: D == hand formula (vnl off, pi = p; explicit + init path)", nfail)
    call check_true(errd(1) < 1d-15, "2-band: D == hand formula (vnl on, pi = p + rvnl; explicit + init path)", nfail)
  end subroutine test_two_band_analytic

  !======================= check e: degeneracy guards =========================
  subroutine test_guards(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: d_with(3, 3), d_without(3, 3)
    integer :: i

    ! (e1) degenerate occupied<->occupied pair with a large p element between
    ! the two occupied bands: excluded by construction (no NaN, no change).
    call alloc_gs(gs, 3, 2, 4)                 ! bands 1,2 occupied (degenerate)
    gs%eigen(1, :) = -0.20d0
    gs%eigen(2, :) = -0.20d0                   ! exact degeneracy inside occ manifold
    gs%eigen(3, :) =  0.40d0
    do i = 1, 3
      gs%p_tm_matrix(1, 3, i, 1) = dcmplx(0.2d0 * i, 0.1d0)
      gs%p_tm_matrix(3, 1, i, 1) = conjg(gs%p_tm_matrix(1, 3, i, 1))
    end do
    call finish_gs(gs)
    call make_solver_gs(sbe, gs, .false.)
    call build_fsum_deficiency_tensor(sbe, gs, .false., 0)
    d_without = sbe%fsum_D
    ! now add a HUGE matrix element between the degenerate occupied bands
    do i = 1, 3
      gs%p_tm_matrix(1, 2, i, 1) = dcmplx(1d6, -2d6)
      gs%p_tm_matrix(2, 1, i, 1) = conjg(gs%p_tm_matrix(1, 2, i, 1))
    end do
    call finish_gs(gs)
    call build_fsum_deficiency_tensor(sbe, gs, .false., 0)
    d_with = sbe%fsum_D
    call check_true(all(d_with == d_with) .and. &                 ! no NaN
                  & maxval(abs(d_with - d_without)) == 0d0, &
                  & "occ<->occ pair (degenerate, huge p): excluded exactly", nfail)
    deallocate(sbe%rho)
    call free_gs(gs)

    ! (e2) occupied x unoccupied pair below the theta_off floor: dropped.
    call alloc_gs(gs, 3, 2, 2)                 ! band 1 occupied
    gs%eigen(1, :) = -0.20d0
    gs%eigen(2, :) = -0.20d0 + 0.5d0 * theta_off   ! Fermi-crossing near-degeneracy
    gs%eigen(3, :) =  0.40d0
    do i = 1, 3
      gs%p_tm_matrix(1, 3, i, 1) = dcmplx(0.2d0 * i, 0.1d0)
      gs%p_tm_matrix(3, 1, i, 1) = conjg(gs%p_tm_matrix(1, 3, i, 1))
    end do
    call finish_gs(gs)
    call make_solver_gs(sbe, gs, .false.)
    call build_fsum_deficiency_tensor(sbe, gs, .false., 0)
    d_without = sbe%fsum_D
    ! huge element on the sub-floor occ x unocc pair: must be dropped, not
    ! amplified (1/de would be ~1e3 * 1e12 without the guard)
    do i = 1, 3
      gs%p_tm_matrix(1, 2, i, 1) = dcmplx(1d6, -2d6)
      gs%p_tm_matrix(2, 1, i, 1) = conjg(gs%p_tm_matrix(1, 2, i, 1))
    end do
    call finish_gs(gs)
    call build_fsum_deficiency_tensor(sbe, gs, .false., 0)
    d_with = sbe%fsum_D
    call check_true(maxval(abs(d_with - d_without)) == 0d0, &
                  & "occ x unocc pair below theta_off floor: dropped exactly", nfail)
    deallocate(sbe%rho)
    call free_gs(gs)
  end subroutine test_guards

end program test_vg_gs_current_subtract
