! src/ssbe/test/test_vnl_kappa.f90
! VG completion (yn_sbe_vnl_exact): all-order nonlocal V_nl(k+A) via
! kappa-stencil interpolation -- analytic-fixture gate for the RUNTIME
! production path (spec: gw_design/specs/2026-07-07-vg-vnl-kappa-interpolation.md).
!
! WHY THIS FILE EXISTS
!   The VG-SBE propagation approximates <u|e^{-i(k+A)r} V_nl e^{i(k+A)r}|u>
!   by the first-order A.rvnl commutator term; the missing higher orders are
!   the 12.3% LG-VG current residual (GC-SBE.qmd).  The exact mode replaces
!   that term by DeltaV = V(k+A)-V(k) in the propagation and by W_i(k+A) in
!   the readout, both interpolated (centered 4-point Lagrange) from a
!   precomputed kappa-stencil.  This test drives the PRODUCTION solver
!   (dt_evolve_bloch / calc_current_bloch / init_sbe_bloch_solver) on an
!   analytic band-limited fixture V(a) = Hvol f(a)^H R f(a) built from the
!   same separable structure as the exporter, so every claimed identity is
!   pinned on the code path production runs:
!
!     1. Lagrange weights: node exactness (t=0 => (0,1,0,0)) + partition of
!        unity at off-node points.
!     2. A=0 BITWISE identity: exact-mode propagation with A=0 equals the
!        legacy path bit for bit (DeltaV(0) == 0 by node-exact weights).
!     3. Conservation: finite-A exact propagation conserves Tr rho_k and
!        Hermiticity to roundoff over many steps (commutator structure).
!     4. Readout interpolation convergence: J from calc_current_bloch at an
!        off-node A converges O(h^4) to the analytic-W reference (ns-scan),
!        and is exact (roundoff) at a stencil node.
!     5. Propagation accuracy: production interpolated propagation matches an
!        INDEPENDENT dense Taylor-4 reference built from the analytic V(a);
!        residual shrinks ~2^4 per ns doubling.
!     6. Small-A legacy consistency, order-separated: propagated rho differs
!        O(A^2) (the legacy A.rvnl term is the first-order expansion of
!        DeltaV) while the READOUT differs O(A) (frozen p+W(0) vs p+W(A);
!        the O(A) coefficient is the nonlocal sum-rule term T of check 7).
!     7. fsum_D wiring: with yn_sbe_vnl_exact='y' the PUBLIC init path builds
!        the deficiency tensor with pi = p + rvnl (rvnl := W_0 convention)
!        PLUS the nonlocal sum-rule term T_i = <sum_v f_v (dW_i/da)_vv>_k / V
!        on the stencil-axis column (checked against an independent analytic
!        off-node reference; legacy modes carry no T -- regression-guarded).
!
! BUILD (already-built ninja tree at build_local/; single-process
! communication dummy).  Compile from a CLEAN dir; link the SAME objects the
! salmon executable built, minus main.f90.o, via an @objs.txt response file:
!
!   find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o' > objs.txt
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<clean_dir> \
!     -c <repo>/src/ssbe/test/test_vnl_kappa.f90 -o <clean_dir>/test_vnl_kappa.o
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     @objs.txt <clean_dir>/test_vnl_kappa.o -o <clean_dir>/test_vnl_kappa \
!     -framework Accelerate -lm -ldl
!   <clean_dir>/test_vnl_kappa
!
program test_vnl_kappa
  use gs_info_ssbe,      only: s_sbe_gs_info
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                               dt_evolve_bloch, calc_current_bloch, &
                               sbe_vnl_interp_weights, build_fsum_deficiency_tensor
  use salmon_global,     only: gauge_sbe, norder_correction, yn_vnl_correction, &
                               yn_sbe_vnl_exact, yn_sbe_gs_current_subtract, &
                               yn_sbe_gw_collision
  implicit none

  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: pi_ = 3.14159265358979323846d0
  ! fixture geometry
  integer, parameter :: nb = 3, nk = 2, ne = 2
  real(8), parameter :: vol = 120.0d0
  real(8), parameter :: hvol_fix = 0.7d0     ! "Hvol" of the analytic model
  ! analytic projector model: nchan channels x npts support points
  integer, parameter :: nchan = 2, npts = 5
  real(8) :: rpts(3, npts)                   ! support coordinates (x,y,z)
  real(8) :: rchan(nchan)                    ! KB weights R(alpha) = +-hvol
  complex(8) :: cofs(npts, nchan, nb, nk)    ! model coefficients
  ! stencil axis = z
  real(8), parameter :: edir(3) = (/ 0d0, 0d0, 1d0 /)
  real(8), parameter :: amax_fix = 0.5d0
  integer :: nfail

  nfail = 0
  call set_globals()
  call build_model()

  call test_weights(nfail)                 ! check 1
  call test_a0_bitwise(nfail)              ! check 2
  call test_conservation(nfail)            ! check 3
  call test_readout_convergence(nfail)     ! check 4
  call test_propagation_reference(nfail)   ! check 5
  call test_small_a_legacy(nfail)          ! check 6
  call test_fsum_wiring(nfail)             ! check 7

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_vnl_kappa checks passed."
  end if

contains

  subroutine set_globals()
    implicit none
    gauge_sbe = 'velocity_gauge'
    norder_correction = 0
    yn_vnl_correction = 'n'
    yn_sbe_vnl_exact = 'n'            ! per-check code overrides
    yn_sbe_gs_current_subtract = 'n'
    yn_sbe_gw_collision = 'n'
  end subroutine set_globals

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

  !======================= analytic separable model ===========================
  ! fc_{alpha,m}(k; a) = sum_j cofs(j,alpha,m,ik) e^{i a rz_j}    (band-limited)
  ! g^{(i)}            = sum_j i r_{i,j} cofs e^{i a rz_j}        (= dfc/da for i=z)
  ! V(a)   = hvol_fix * fc^H R fc          (Hermitian by construction)
  ! W_i(a) = hvol_fix * (g_i^H R fc + fc^H R g_i)
  ! Deterministic pseudo-random coefficients (no RNG state dependence).
  subroutine build_model()
    implicit none
    integer :: j, a, m, ik
    real(8) :: s1, s2
    rpts(1, :) = (/ -1.7d0,  0.4d0,  1.1d0, -0.6d0,  2.0d0 /)
    rpts(2, :) = (/  0.9d0, -1.3d0,  0.2d0,  1.8d0, -0.5d0 /)
    rpts(3, :) = (/ -2.1d0, -0.8d0,  0.3d0,  1.2d0,  2.4d0 /)   ! parallel (z)
    rchan(1) = +hvol_fix
    rchan(2) = -hvol_fix
    do ik = 1, nk
      do m = 1, nb
        do a = 1, nchan
          do j = 1, npts
            s1 = sin(dble(3*j + 5*a + 7*m + 11*ik)) * 0.6d0
            s2 = cos(dble(2*j + 3*a + 5*m + 13*ik)) * 0.6d0
            cofs(j, a, m, ik) = dcmplx(s1, s2)
          end do
        end do
      end do
    end do
  end subroutine build_model

  subroutine model_fc(ik, a, fc)
    implicit none
    integer, intent(in) :: ik
    real(8), intent(in) :: a
    complex(8), intent(out) :: fc(nchan, nb)
    integer :: j, al, m
    complex(8) :: ph
    fc = (0d0, 0d0)
    do m = 1, nb
      do al = 1, nchan
        do j = 1, npts
          ph = exp(zi_ * a * rpts(3, j))
          fc(al, m) = fc(al, m) + cofs(j, al, m, ik) * ph
        end do
      end do
    end do
  end subroutine model_fc

  subroutine model_gc(ik, a, i, gc)
    implicit none
    integer, intent(in) :: ik, i
    real(8), intent(in) :: a
    complex(8), intent(out) :: gc(nchan, nb)
    integer :: j, al, m
    complex(8) :: ph
    gc = (0d0, 0d0)
    do m = 1, nb
      do al = 1, nchan
        do j = 1, npts
          ph = exp(zi_ * a * rpts(3, j))
          gc(al, m) = gc(al, m) + zi_ * rpts(i, j) * cofs(j, al, m, ik) * ph
        end do
      end do
    end do
  end subroutine model_gc

  subroutine model_V(ik, a, v)
    implicit none
    integer, intent(in) :: ik
    real(8), intent(in) :: a
    complex(8), intent(out) :: v(nb, nb)
    complex(8) :: fc(nchan, nb)
    integer :: n, m, al
    call model_fc(ik, a, fc)
    do m = 1, nb
      do n = 1, nb
        v(n, m) = (0d0, 0d0)
        do al = 1, nchan
          v(n, m) = v(n, m) + conjg(fc(al, n)) * rchan(al) * fc(al, m)
        end do
        v(n, m) = hvol_fix * v(n, m)
      end do
    end do
  end subroutine model_V

  subroutine model_W(ik, a, i, w)
    implicit none
    integer, intent(in) :: ik, i
    real(8), intent(in) :: a
    complex(8), intent(out) :: w(nb, nb)
    complex(8) :: fc(nchan, nb), gc(nchan, nb)
    integer :: n, m, al
    call model_fc(ik, a, fc)
    call model_gc(ik, a, i, gc)
    do m = 1, nb
      do n = 1, nb
        w(n, m) = (0d0, 0d0)
        do al = 1, nchan
          w(n, m) = w(n, m) + conjg(gc(al, n)) * rchan(al) * fc(al, m) &
                          & + conjg(fc(al, n)) * rchan(al) * gc(al, m)
        end do
        w(n, m) = hvol_fix * w(n, m)
      end do
    end do
  end subroutine model_W

  !======================= gs fixture =========================================
  subroutine build_gs(gs, ns)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    integer, intent(in) :: ns
    integer :: ik, ib, jb, s, i
    complex(8) :: v(nb, nb), w(nb, nb)
    real(8) :: h

    gs%nk = nk;  gs%nb = nb;  gs%ne = ne
    gs%focc = 2d0;  gs%nvb = ne / 2
    gs%volume = vol
    if (.not. allocated(gs%kweight)) then
      allocate(gs%kweight(nk), gs%eigen(nb, nk), gs%occup(nb, nk))
      allocate(gs%delta_omega(nb, nb, nk))
      allocate(gs%p_tm_matrix(nb, nb, 3, nk), gs%rvnl_tm_matrix(nb, nb, 3, nk))
    end if
    gs%kweight(:) = 1.0d0
    gs%eigen(1, :) = -0.25d0
    gs%eigen(2, :) =  0.15d0
    gs%eigen(3, :) =  0.55d0
    gs%occup(:, :) = 0d0
    gs%occup(1:ne/2, :) = 2d0
    do ik = 1, nk
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do
    end do
    ! p: Hermitian, deterministic
    do ik = 1, nk
      do i = 1, 3
        do jb = 1, nb
          do ib = 1, jb
            gs%p_tm_matrix(ib, jb, i, ik) = dcmplx( &
              & 0.11d0 * sin(dble(ib + 2*jb + 3*i + ik)), &
              & 0.07d0 * cos(dble(2*ib + jb + i + 2*ik)))
            if (ib == jb) gs%p_tm_matrix(ib, jb, i, ik) = &
              & dcmplx(dble(gs%p_tm_matrix(ib, jb, i, ik)), 0d0)
            gs%p_tm_matrix(jb, ib, i, ik) = conjg(gs%p_tm_matrix(ib, jb, i, ik))
          end do
        end do
      end do
    end do
    ! rvnl := analytic W_i(0)  (the exact-mode reader-overwrite convention)
    do ik = 1, nk
      do i = 1, 3
        call model_W(ik, 0d0, i, w)
        gs%rvnl_tm_matrix(:, :, i, ik) = w(:, :)
      end do
    end do

    ! kappa-stencil tables at the nodes
    h = amax_fix / dble(ns)
    gs%vnl_exact  = .true.
    gs%vnl_ns     = ns
    gs%vnl_h      = h
    gs%vnl_dir(1:3) = edir(1:3)
    gs%vnl_ik_min = 1
    gs%vnl_ik_max = nk
    if (allocated(gs%vnl_V)) deallocate(gs%vnl_V)
    if (allocated(gs%vnl_W)) deallocate(gs%vnl_W)
    allocate(gs%vnl_V(nb, nb, -ns:ns, 1:nk))
    allocate(gs%vnl_W(nb, nb, 3, -ns:ns, 1:nk))
    do ik = 1, nk
      do s = -ns, ns
        call model_V(ik, dble(s) * h, v)
        gs%vnl_V(:, :, s, ik) = v(:, :)
        do i = 1, 3
          call model_W(ik, dble(s) * h, i, w)
          gs%vnl_W(:, :, i, s, ik) = w(:, :)
        end do
      end do
    end do
  end subroutine build_gs

  ! independent norder=0 VG readout (the legacy formula, coded here, with a
  ! caller-supplied velocity matrix vmat = p + vnl_velocity(a))
  subroutine reference_current(gs, rho, vmat, Ac, jm)
    implicit none
    type(s_sbe_gs_info), intent(in) :: gs
    complex(8), intent(in) :: rho(nb, nb, nk), vmat(nb, nb, 3, nk)
    real(8), intent(in) :: Ac(3)
    real(8), intent(out) :: jm(3)
    complex(8) :: acc(3)
    real(8) :: tr
    integer :: ik, i, n, m
    acc = (0d0, 0d0);  tr = 0d0
    do ik = 1, nk
      do i = 1, 3
        do n = 1, nb
          do m = 1, nb
            acc(i) = acc(i) + gs%kweight(ik) * rho(m, n, ik) * vmat(n, m, i, ik)
          end do
        end do
      end do
      do n = 1, nb
        tr = tr + gs%kweight(ik) * dble(rho(n, n, ik))
      end do
    end do
    jm(1:3) = (dble(acc(1:3)) / sum(gs%kweight) + Ac(1:3) * tr / sum(gs%kweight)) / gs%volume
  end subroutine reference_current

  !======================= check 1: Lagrange weights ==========================
  subroutine test_weights(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    real(8) :: wl(4), a0(3)
    integer :: q0
    call build_gs(gs, 8)
    ! node a = 0
    a0 = 0d0
    call sbe_vnl_interp_weights(gs, a0, wl, q0)
    call check_true(q0 == 0 .and. abs(wl(2) - 1d0) < 1d-15 .and. &
      & abs(wl(1)) + abs(wl(3)) + abs(wl(4)) < 1d-15, &
      & "weights: exact at node a=0 -> (0,1,0,0)", nfail)
    ! node a = +2h
    a0 = 0d0;  a0(3) = 2d0 * gs%vnl_h
    call sbe_vnl_interp_weights(gs, a0, wl, q0)
    call check_true(q0 == 2 .and. abs(wl(2) - 1d0) < 1d-14 .and. &
      & abs(wl(1)) + abs(wl(3)) + abs(wl(4)) < 1d-14, &
      & "weights: exact at node a=+2h", nfail)
    ! off-node, negative a: partition of unity
    a0 = 0d0;  a0(3) = -0.37d0 * amax_fix
    call sbe_vnl_interp_weights(gs, a0, wl, q0)
    call check_true(abs(sum(wl) - 1d0) < 1d-14, &
      & "weights: partition of unity off-node (a<0)", nfail)
  end subroutine test_weights

  !======================= check 2: A=0 bitwise identity ======================
  subroutine test_a0_bitwise(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe_e, sbe_l
    real(8) :: a0(3)
    integer :: icomm, it, ik, n, m
    logical :: same
    call build_gs(gs, 8)
    icomm = 0
    yn_sbe_vnl_exact = 'y'
    call init_sbe_bloch_solver(sbe_e, gs, nb, icomm)
    yn_sbe_vnl_exact = 'n'
    call init_sbe_bloch_solver(sbe_l, gs, nb, icomm)
    sbe_e%flag_vnl_correction = .false.
    sbe_l%flag_vnl_correction = .false.
    ! identical generic excited state (Hermitian, trace = ne)
    call seed_rho(sbe_e)
    call seed_rho(sbe_l)
    a0 = 0d0
    do it = 1, 10
      call dt_evolve_bloch(sbe_e, gs, a0, 0.01d0)
      call dt_evolve_bloch(sbe_l, gs, a0, 0.01d0)
    end do
    same = .true.
    do ik = 1, nk
      do n = 1, nb
        do m = 1, nb
          if (sbe_e%rho(m, n, ik) /= sbe_l%rho(m, n, ik)) same = .false.
        end do
      end do
    end do
    call check_true(same, "A=0: exact-mode propagation BITWISE equal to legacy", nfail)
  end subroutine test_a0_bitwise

  subroutine seed_rho(sbe)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    integer :: ik, n, m
    do ik = 1, nk
      sbe%rho(:, :, ik) = (0d0, 0d0)
      sbe%rho(1, 1, ik) = dcmplx(1.7d0, 0d0)
      sbe%rho(2, 2, ik) = dcmplx(0.25d0, 0d0)
      sbe%rho(3, 3, ik) = dcmplx(0.05d0, 0d0)
      do n = 1, nb
        do m = n+1, nb
          sbe%rho(n, m, ik) = dcmplx(0.03d0 * dble(n), -0.02d0 * dble(m + ik))
          sbe%rho(m, n, ik) = conjg(sbe%rho(n, m, ik))
        end do
      end do
    end do
  end subroutine seed_rho

  !======================= check 3: conservation ==============================
  subroutine test_conservation(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: ac(3), tr0(nk), tr, herm
    integer :: icomm, it, ik, n, m
    call build_gs(gs, 16)
    icomm = 0
    yn_sbe_vnl_exact = 'y'
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    yn_sbe_vnl_exact = 'n'
    sbe%flag_vnl_correction = .false.
    call seed_rho(sbe)
    do ik = 1, nk
      tr0(ik) = 0d0
      do n = 1, nb
        tr0(ik) = tr0(ik) + dble(sbe%rho(n, n, ik))
      end do
    end do
    do it = 1, 200
      ac = 0d0
      ac(3) = 0.6d0 * amax_fix * sin(0.15d0 * dble(it))
      call dt_evolve_bloch(sbe, gs, ac, 0.005d0)
    end do
    tr = 0d0;  herm = 0d0
    do ik = 1, nk
      do n = 1, nb
        tr = max(tr, abs(dble(sum_diag(sbe%rho(:, :, ik))) - tr0(ik)))
        do m = 1, nb
          herm = max(herm, abs(sbe%rho(n, m, ik) - conjg(sbe%rho(m, n, ik))))
        end do
      end do
    end do
    call check_true(tr < 1d-10, "conservation: per-k trace conserved over 200 exact steps", nfail)
    call check_true(herm < 1d-10, "conservation: Hermiticity preserved over 200 exact steps", nfail)
  end subroutine test_conservation

  function sum_diag(a) result(t)
    implicit none
    complex(8), intent(in) :: a(nb, nb)
    complex(8) :: t
    integer :: n
    t = (0d0, 0d0)
    do n = 1, nb
      t = t + a(n, n)
    end do
  end function sum_diag

  !======================= check 4: readout convergence =======================
  subroutine test_readout_convergence(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: ac(3), jm(3), jref(3), err(4), rat
    complex(8) :: vmat(nb, nb, 3, nk), w(nb, nb)
    complex(8) :: rho_fix(nb, nb, nk)
    integer :: icomm, ik, i, ins, nss(4)
    icomm = 0
    nss = (/ 4, 8, 16, 32 /)
    do ins = 1, 4
      call build_gs(gs, nss(ins))
      call free_solver(sbe)
      yn_sbe_vnl_exact = 'y'
      call init_sbe_bloch_solver(sbe, gs, nb, icomm)
      yn_sbe_vnl_exact = 'n'
      sbe%flag_vnl_correction = .false.
      call seed_rho(sbe)
      rho_fix(:, :, :) = sbe%rho(:, :, 1:nk)
      ! off-node for every ns in the scan: a = 0.2345*amax
      ac = 0d0;  ac(3) = 0.2345d0 * amax_fix
      call calc_current_bloch(sbe, gs, ac, jm, icomm)
      ! reference with the ANALYTIC W(a) in an independent readout
      do ik = 1, nk
        do i = 1, 3
          call model_W(ik, ac(3), i, w)
          vmat(:, :, i, ik) = gs%p_tm_matrix(:, :, i, ik) + w(:, :)
        end do
      end do
      call reference_current(gs, rho_fix, vmat, ac, jref)
      err(ins) = maxval(abs(jm - jref))
    end do
    write(*, '(a,4es12.4)') "  readout interpolation errors (ns=4,8,16,32): ", err
    ! O(h^4) magnitude: (B*h)^4-scale with B = max|r_z| = 2.4, so ns=32
    ! (h=amax/32) sits at ~1e-8; the exact 16x-per-doubling is modulated by
    ! the Lagrange pi(t) factor (t moves with ns), hence the >8 gates.
    rat = err(2) / max(err(3), 1d-300)
    call check_true(err(4) < err(1) .and. err(4) < 1d-7, &
      & "readout: interpolation error decreases and is O((Bh)^4)-small at ns=32", nfail)
    call check_true(rat > 8d0 .and. err(3) / max(err(4), 1d-300) > 8d0, &
      & "readout: O(h^4) convergence (error ratio > 8 per ns doubling)", nfail)

    ! node exactness through the production readout: a = +3h (ns=32 fixture)
    ac = 0d0;  ac(3) = 3d0 * gs%vnl_h
    call calc_current_bloch(sbe, gs, ac, jm, icomm)
    do ik = 1, nk
      do i = 1, 3
        call model_W(ik, ac(3), i, w)
        vmat(:, :, i, ik) = gs%p_tm_matrix(:, :, i, ik) + w(:, :)
      end do
    end do
    call reference_current(gs, rho_fix, vmat, ac, jref)
    call check_true(maxval(abs(jm - jref)) < 1d-13, &
      & "readout: exact (roundoff) at a stencil node", nfail)
  end subroutine test_readout_convergence

  !======================= check 5: propagation vs dense reference ============
  subroutine test_propagation_reference(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    complex(8) :: rho_ref(nb, nb, nk), hmat(nb, nb), v0(nb, nb), va(nb, nb)
    complex(8) :: c1(nb, nb), c2(nb, nb), c3(nb, nb), c4(nb, nb)
    real(8) :: ac(3), dt, errs(2)
    integer :: icomm, it, ik, i, ins, nss(2), n
    icomm = 0
    dt = 0.005d0
    nss = (/ 8, 16 /)
    do ins = 1, 2
      call build_gs(gs, nss(ins))
      call free_solver(sbe)
      yn_sbe_vnl_exact = 'y'
      call init_sbe_bloch_solver(sbe, gs, nb, icomm)
      yn_sbe_vnl_exact = 'n'
      sbe%flag_vnl_correction = .false.
      call seed_rho(sbe)
      ! independent reference: same Taylor-4, dense H from the ANALYTIC V(a)
      rho_ref(:, :, :) = sbe%rho(:, :, 1:nk)
      do it = 1, 100
        ac = 0d0
        ac(3) = 0.55d0 * amax_fix * sin(0.2d0 * dble(it))
        call dt_evolve_bloch(sbe, gs, ac, dt)
        do ik = 1, nk
          call model_V(ik, 0d0, v0)
          call model_V(ik, ac(3), va)
          hmat = va - v0
          do n = 1, nb
            hmat(n, n) = hmat(n, n) + gs%eigen(n, ik)
          end do
          do i = 1, 3
            hmat = hmat + ac(i) * gs%p_tm_matrix(:, :, i, ik)
          end do
          c1 = matmul(hmat, rho_ref(:, :, ik)) - matmul(rho_ref(:, :, ik), hmat)
          c2 = matmul(hmat, c1) - matmul(c1, hmat)
          c3 = matmul(hmat, c2) - matmul(c2, hmat)
          c4 = matmul(hmat, c3) - matmul(c3, hmat)
          rho_ref(:, :, ik) = rho_ref(:, :, ik) &
            & + c1 * (-zi_ * dt) + c2 * (-zi_ * dt)**2 / 2d0 &
            & + c3 * (-zi_ * dt)**3 / 6d0 + c4 * (-zi_ * dt)**4 / 24d0
        end do
      end do
      errs(ins) = 0d0
      do ik = 1, nk
        errs(ins) = max(errs(ins), maxval(abs(sbe%rho(:, :, ik) - rho_ref(:, :, ik))))
      end do
    end do
    write(*, '(a,2es12.4)') "  propagation vs dense-analytic reference (ns=8,16): ", errs
    call check_true(errs(2) < 5d-6, &
      & "propagation: matches independent dense reference at ns=16", nfail)
    call check_true(errs(1) / max(errs(2), 1d-300) > 8d0, &
      & "propagation: residual shrinks O(h^4) with stencil refinement", nfail)
  end subroutine test_propagation_reference

  !======================= check 6: small-A legacy consistency ================
  ! The legacy first-order mode differs from the exact mode in TWO cleanly
  ! separable orders (this is the physics of the 12.3% gate, decomposed):
  !   (a) PROPAGATION: A.rvnl == A.W(0) is the first-order expansion of
  !       DeltaV(A), so the propagated rho differs at O(A^2);
  !   (b) READOUT: legacy norder=0 uses the FROZEN velocity p + W(0), while
  !       the exact readout uses p + W(A) -- an O(A) gap (during the pulse;
  !       post-pulse A=0 both readouts coincide and only (a) survives).
  subroutine test_small_a_legacy(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe_e, sbe_l
    real(8) :: ac(3), jm_e(3), jm_l(3), rmax(2), dread(2), eps, rat_rho, rat_read
    integer :: icomm, it, ia, ik
    icomm = 0
    call build_gs(gs, 16)
    do ia = 1, 2
      eps = 0.02d0 / dble(2**(ia-1))     ! 0.02, 0.01
      yn_sbe_vnl_exact = 'y'
      call init_sbe_bloch_solver(sbe_e, gs, nb, icomm)
      yn_sbe_vnl_exact = 'n'
      call init_sbe_bloch_solver(sbe_l, gs, nb, icomm)
      sbe_e%flag_vnl_correction = .false.
      sbe_l%flag_vnl_correction = .true.   ! legacy first-order A.rvnl (rvnl = W_0)
      call seed_rho(sbe_e)
      call seed_rho(sbe_l)
      ! (b) readout gap at FIXED rho (same seeded state, no propagation).
      ! Smaller eps than the propagation part: at eps=0.02 the O(A^2) tail of
      ! W(A)-W(0) still contributes ~30% and blurs the O(A) ratio.
      ac = 0d0;  ac(3) = 0.2d0 * eps
      call calc_current_bloch(sbe_e, gs, ac, jm_e, icomm)
      call calc_current_bloch(sbe_l, gs, ac, jm_l, icomm)
      dread(ia) = maxval(abs(jm_e - jm_l))
      ! (a) propagation gap: max |rho_e - rho_l| over the run
      rmax(ia) = 0d0
      do it = 1, 100
        ac = 0d0
        ac(3) = eps * sin(0.2d0 * dble(it))
        call dt_evolve_bloch(sbe_e, gs, ac, 0.005d0)
        call dt_evolve_bloch(sbe_l, gs, ac, 0.005d0)
        do ik = 1, nk
          rmax(ia) = max(rmax(ia), maxval(abs(sbe_e%rho(:, :, ik) - sbe_l%rho(:, :, ik))))
        end do
      end do
      call free_solver(sbe_e)
      call free_solver(sbe_l)
    end do
    rat_rho  = rmax(1) / max(rmax(2), 1d-300)
    rat_read = dread(1) / max(dread(2), 1d-300)
    write(*, '(a,2es12.4,a,f7.3)') "  propagation |rho_e - rho_l| at eps, eps/2: ", rmax, &
      & "  ratio: ", rat_rho
    write(*, '(a,2es12.4,a,f7.3)') "  readout gap at fixed rho at eps, eps/2:    ", dread, &
      & "  ratio: ", rat_read
    call check_true(rat_rho > 3.2d0 .and. rat_rho < 5.0d0, &
      & "small-A: propagation difference scales O(A^2) (legacy = 1st-order DeltaV)", nfail)
    call check_true(rat_read > 1.7d0 .and. rat_read < 2.4d0, &
      & "small-A: readout gap scales O(A) (frozen W(0) vs exact W(A))", nfail)
  end subroutine test_small_a_legacy

  subroutine free_solver(sbe)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    if (allocated(sbe%rho)) deallocate(sbe%rho)
  end subroutine free_solver

  !======================= check 7: fsum_D exact-mode wiring ==================
  ! (a) exact mode builds D through the public init with pi = p + W_0 AND the
  !     nonlocal sum-rule term T (e_dir column): D_exact = D_legacy + (T/V) e^T,
  !     T_i = < sum_v f_v Re (dW_i/da)_vv |_{a=0} >_k, checked against an
  !     INDEPENDENT reference (central difference of the ANALYTIC model_W at
  !     a = +-delta, delta << h -- not the production 5-point node stencil).
  ! (b) legacy (non-exact) D is unchanged by the T machinery (regression).
  subroutine test_fsum_wiring(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe_a, sbe_b
    integer :: icomm, ik, iv, i, j
    real(8) :: tref(3), dexp(3, 3), delta, wsum, derr
    complex(8) :: wp(nb, nb), wm(nb, nb)
    icomm = 0
    ! ns=64 (h = amax/64): the production T uses the O(h^4) five-point node
    ! stencil, so the gap to the off-node analytic reference is ~(B*h)^4/30
    ! ~ 7e-8 here -- comfortably inside the 1e-6 assertion; at the coarse
    ! ns=8 fixture it would be ~3e-4 and the tolerance would mask real bugs.
    call build_gs(gs, 64)
    ! PUBLIC init path with exact mode + subtraction: D must be built with
    ! pi = p + rvnl (rvnl carries W_0 by the reader-overwrite convention)
    yn_sbe_vnl_exact = 'y'
    yn_vnl_correction = 'n'
    yn_sbe_gs_current_subtract = 'y'
    call init_sbe_bloch_solver(sbe_a, gs, nb, icomm)
    yn_sbe_vnl_exact = 'n'
    yn_sbe_gs_current_subtract = 'n'
    ! explicit reference: flag_vnl = .true., non-exact solver (no T term)
    call init_sbe_bloch_solver(sbe_b, gs, nb, icomm)
    call build_fsum_deficiency_tensor(sbe_b, gs, .true., icomm)
    ! independent T reference from the analytic model (off-node FD, O(delta^2))
    delta = 1d-5
    tref(:) = 0d0
    wsum = sum(gs%kweight(1:nk))
    do ik = 1, nk
      do i = 1, 3
        call model_W(ik, +delta, i, wp)
        call model_W(ik, -delta, i, wm)
        do iv = 1, gs%nvb
          tref(i) = tref(i) + gs%kweight(ik) * gs%occup(iv, ik) &
            & * dble(wp(iv, iv) - wm(iv, iv)) / (2d0 * delta)
        end do
      end do
    end do
    tref(:) = tref(:) / (wsum * gs%volume)
    do j = 1, 3
      do i = 1, 3
        dexp(i, j) = sbe_b%fsum_D(i, j) + tref(i) * edir(j)
      end do
    end do
    derr = maxval(abs(sbe_a%fsum_D - dexp)) / max(maxval(abs(dexp)), 1d-300)
    write(*, '(a,3es12.4)') "  fsum_D exact-mode T column (analytic ref): ", tref
    write(*, '(a,es12.4)')  "  fsum_D exact vs legacy+T rel residual:     ", derr
    call check_true(sbe_a%fsum_D_built .and. derr < 1d-6, &
      & "fsum_D: exact mode adds the nonlocal sum-rule term T to pi = p + W_0", nfail)
    call check_true(maxval(abs(tref)) > 1d-12 .and. &
      & maxval(abs(sbe_a%fsum_D - sbe_b%fsum_D)) > 1d-12, &
      & "fsum_D: fixture T is nonzero (the check cannot pass vacuously)", nfail)
  end subroutine test_fsum_wiring

end program test_vnl_kappa
