!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!===================================================================
! Integral (covariant-Houston) transport propagator kernels.
!
! Purpose: replace the finite-difference gauge-covariant k-derivative of the
! length-gauge SBE by its EXACT path-ordered (Wilson) integral form, for a
! single reciprocal-mesh axis of linear polarization.  In the co-moving frame
! the density obeys the k-LOCAL, purely unitary law
!
!     d rho~/dt = -i [ H~(kappa,t), rho~ ],
!     H~(kappa,t) = W^dagger H0(kappa - a(t)) W,     x = kappa - a(t)
!
! with W the window-projected Wilson transport (a bounded unitary applied
! MULTIPLICATIVELY, never as a finite difference), a(t) the mesh shift built
! from the external field (a = -[A_ext - A_ext(0)], da/dt = E; see the field
! convention note in realtime_ssbe).  The singular Berry connection near a band
! touching (norm ~ 1/(2r)) enters only through the bounded, unitary W, so the
! mesh-density growth of the finite-difference stencil floor is removed by
! construction.
!
! This module holds ONLY the pure / LAPACK-only numerical kernels, with NO
! SALMON communication / filesystem / global-variable dependency (the same
! standalone-testable discipline as degenerate_block_ssbe.f90): it compiles and
! unit-tests on its own against a 2-band analytic model
!   gfortran sbe_lg_mode_ssbe.f90 degenerate_block_ssbe.f90 \
!            gicov_integral_ssbe.f90 \
!            test/test_gicov_integral.f90 -o t -framework Accelerate
! The MPI k-halo gather that assembles the per-axis single-step link chain
! across ranks, and the SALMON wiring (a(t) cache, dispatch, current), live in
! bloch_solver_ssbe.f90 / gs_info_ssbe.f90 / realtime_ssbe.f90 and consume
! these kernels.
!===================================================================
module gicov_integral_ssbe
  use degenerate_block_ssbe, only: polar_unitary
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: gicov_int_jmax           ! j_max = ceil(a_max / dk) mesh shifts spanned by the pulse
  public :: gicov_int_cache_bytes    ! int64 transport-cache footprint per rank (H~_j + 3 v~_j)
  public :: gicov_int_floor_shift    ! (a,dk) -> floor node shift n_int + fractional weight (floor(), not int())
  public :: gicov_int_bracket        ! (a,dk,jmax) -> the two cached nodes bracketing x, span-checked
  public :: gicov_int_axis_single    ! runtime all-trajectory single-axis / linear-polarization guard
  public :: gicov_int_gate_weight    ! Delta-omega T2 gate weight (mirror of gs_info t2_gate_weight)
  public :: gicov_int_degen_blocks   ! connected closure of |eps_a-eps_b|<=floor over the sorted spectrum
  public :: gicov_int_chain          ! ordered product of single-step Wilson links -> W(kappa,kappa+n), polar
  public :: gicov_int_transport_op   ! Y = W O W^dagger  (transport a band-frame operator into the kappa frame)
  public :: gicov_int_interp         ! linear interpolation of a bounded transported operator between nodes
  public :: gicov_int_step_k         ! midpoint-frozen EXACT-exponential step of one k-block + gated T2
  public :: gicov_int_current_k      ! Re Tr[rho~ v~_i] contribution of one k-block (i=1..3)
  public :: gicov_int_occupation_k   ! instantaneous-eigenbasis populations n_a(x) = (P^dag rho~ P)_aa

contains

  !-------------------------------------------------------------------
  ! Number of integer mesh shifts the pulse spans on the driven axis:
  !   j_max = ceil( a_max / dk ),   a_max = max_t |a(t).a_axis|,
  ! dk = reduced-mesh spacing on that axis.  The transport cache then covers
  ! j in [-j_max, +j_max] (2*j_max+1 shifts).  graphene k99 production:
  ! a_max=0.102424, dk=0.00788381 -> j_max = 13 (NOT ~7; a whole-pulse bound,
  ! not a per-step one).
  !-------------------------------------------------------------------
  pure integer function gicov_int_jmax(a_max, dk) result(jmax)
    implicit none
    real(8), intent(in) :: a_max, dk
    jmax = ceiling(abs(a_max) / dk)
    if (jmax < 1) jmax = 1
  end function gicov_int_jmax

  !-------------------------------------------------------------------
  ! Per-rank transport-cache footprint in BYTES, computed in integer(8) so it
  ! never wraps for production (nb up to ~192, nk_local up to ~O(10^4)).  The
  ! cache keeps FOUR complex(8) nb x nb matrices per shift -- H~_j and the
  ! three velocity components v~_{j,1..3} -- over 2*j_max+1 shifts and the
  ! LOCAL k-slice only (ik_min:ik_max = nk_local), never the full nk:
  !     bytes = 16 * 4 * (2*j_max+1) * nb^2 * nk_local
  !           = 64 * (2*j_max+1) * nb^2 * nk_local.
  ! (The factor is 64 = 4 matrices x 16 bytes, NOT 32; and the k extent is
  ! nk_local, NOT nk -- a full-nk cache would be ~9-12 GiB/rank at production
  ! scale and must be rejected by the lint memory gate.)
  !-------------------------------------------------------------------
  pure integer(8) function gicov_int_cache_bytes(jmax, nb, nk_local) result(nbytes)
    implicit none
    integer, intent(in) :: jmax, nb, nk_local
    nbytes = 64_8 * int(2 * jmax + 1, 8) * int(nb, 8) * int(nb, 8) * int(nk_local, 8)
  end function gicov_int_cache_bytes

  !-------------------------------------------------------------------
  ! Bracketing integer node shift and fractional weight for the off-mesh
  ! evaluation point x = kappa - a on the driven axis.  The physical offset
  ! from kappa is -a, i.e. s = -a/dk mesh steps, so
  !   n_int = FLOOR(s),   frac = s - n_int  in [0,1),
  ! and the bounded transported operator is interpolated as
  !   Y(x) = (1-frac) Y[n_int] + frac Y[n_int+1].
  ! FLOOR (round toward -infinity) is mandatory: INT()/truncation rounds toward
  ! zero and would jump the bracket by one for a < 0, biasing the shift.
  !-------------------------------------------------------------------
  pure subroutine gicov_int_floor_shift(a, dk, n_int, frac)
    implicit none
    real(8), intent(in)  :: a, dk
    integer, intent(out) :: n_int
    real(8), intent(out) :: frac
    real(8) :: s
    s = -a / dk
    n_int = floor(s)
    frac = s - real(n_int, 8)
  end subroutine gicov_int_floor_shift

  !-------------------------------------------------------------------
  ! The two CACHED nodes bracketing the evaluation point x = kappa - a, with the
  ! interpolation weight, span-checked against the cache extent [-jmax, jmax].
  !
  ! Wraps gicov_int_floor_shift and fixes its endpoint: when x lands exactly on
  ! a node (frac = 0) the upper bracket must be that node ITSELF, not n_lo+1.
  ! Otherwise the legal endpoint a = +j_max*dk (i.e. q = -j_max, the far turning
  ! point of the pulse the cache was sized for) asks for node j_max+1 -- one past
  ! the cache -- and the run aborts (or, unguarded, reads out of bounds), even
  ! though that node carries weight frac = 0 and is never actually needed.
  !
  ! ierr /= 0 iff the bracket leaves the cached span (a genuine sizing error:
  ! the pulse moved further than j_max was built for) -- reported, never clamped.
  !-------------------------------------------------------------------
  pure subroutine gicov_int_bracket(a, dk, jmax, n_lo, n_hi, frac, ierr)
    implicit none
    real(8), intent(in)  :: a, dk
    integer, intent(in)  :: jmax
    integer, intent(out) :: n_lo, n_hi, ierr
    real(8), intent(out) :: frac
    call gicov_int_floor_shift(a, dk, n_lo, frac)
    if (frac == 0d0) then
      n_hi = n_lo            ! exactly on a node: collapse, never reach past it
    else
      n_hi = n_lo + 1
    end if
    if (n_lo < -jmax .or. n_hi > jmax) then
      ierr = 1
    else
      ierr = 0
    end if
  end subroutine gicov_int_bracket

  !-------------------------------------------------------------------
  ! Runtime linear-polarization / single-axis guard, evaluated on the WHOLE
  ! precomputed trajectory (NOT a deck vector check): given the reduced-mesh
  ! displacement q(i,t) = N_i * (a(t).a_i) / 2pi for the three reciprocal axes
  ! i=1..3 at every time sample t=1..nt, require exactly ONE axis to move.
  ! qmax(i) = max_t |q(i,t)|; the moved set = { i : qmax(i) > tol }.  Accept
  ! (axis = that i) iff the moved set has cardinality 1.  Reject (axis = 0):
  !   - 0 moved axes  (no field -> nothing to transport; caller may treat as
  !     trivial, but a genuine pulse always moves one axis),
  !   - >= 2 moved axes ([110] non-axis linear polarization, or elliptic /
  !     circular where every component is time-nonzero).
  ! tol is a scale-aware absolute tolerance the caller sets from max qmax.
  !-------------------------------------------------------------------
  pure subroutine gicov_int_axis_single(q, naxis, nt, tol, ok, axis)
    implicit none
    integer, intent(in)  :: naxis, nt
    real(8), intent(in)  :: q(naxis, nt)
    real(8), intent(in)  :: tol
    logical, intent(out) :: ok
    integer, intent(out) :: axis
    integer :: i, it, nmoved, last
    real(8) :: qm
    nmoved = 0
    last = 0
    do i = 1, naxis
      qm = 0d0
      do it = 1, nt
        qm = max(qm, abs(q(i, it)))
      end do
      if (qm > tol) then
        nmoved = nmoved + 1
        last = i
      end if
    end do
    ok = (nmoved == 1)
    if (ok) then
      axis = last
    else
      axis = 0
    end if
  end subroutine gicov_int_axis_single

  !-------------------------------------------------------------------
  ! Delta-omega T2 gate weight, evaluated on the INSTANTANEOUS (co-moving)
  ! band gap Delta-omega = eps_a - eps_b.  Bit-identical semantics to
  ! gs_info_ssbe's t2_gate_weight (which is only reachable through the SALMON
  ! dependency tree, hence replicated here to keep this module standalone):
  !   |Delta-omega| <= floor      -> 0   (exact-degeneracy clamp: NEVER dephase
  !                                        within a Kramers / touching block ->
  !                                        gauge covariance g(0)=0)
  !   shape 'gauss'               -> 1 - exp(-(Delta-omega/width)^2)
  !   shape 'step' (default)      -> 1 if |Delta-omega| > theta else 0.
  ! The RATE is NEVER interpolated: it is recomputed here from the moving gap.
  !-------------------------------------------------------------------
  pure real(8) function gicov_int_gate_weight(delta_omega, shape, theta, width, floor) result(w)
    implicit none
    real(8),      intent(in) :: delta_omega, theta, width, floor
    character(*), intent(in) :: shape
    real(8) :: adw
    adw = abs(delta_omega)
    if (adw <= floor) then
      w = 0d0
    else if (trim(shape) == 'gauss') then
      w = 1d0 - exp( -(delta_omega / width)**2 )
    else
      w = merge(1d0, 0d0, adw > theta)
    end if
  end function gicov_int_gate_weight

  !-------------------------------------------------------------------
  ! Instantaneous degenerate blocks = the CONNECTED COMPONENTS of the relation
  ! |eps_a - eps_b| <= floor.  blk(a) is the (1-based, ascending) block index of
  ! eigenvalue a.
  !
  ! Why the connected closure and not the pairwise test: "within floor of" is
  ! NOT transitive, so a near-degenerate CHAIN eps = (0, 0.75*floor, 1.5*floor)
  ! has both adjacent pairs inside the floor while the END pair (0, 1.5*floor)
  ! sits outside it.  Gating pair-by-pair therefore dephases the end pair of a
  ! manifold whose members are physically indistinguishable -- and the answer
  ! then depends on which basis the eigensolver happened to pick INSIDE that
  ! manifold, which is exactly the gauge dependence the floor exists to kill.
  ! The dephasing rate must be a function of the BLOCK, not of the pair.
  !
  ! zheev returns eps ASCENDING, so the connected components are precisely the
  ! runs delimited by an adjacent gap greater than floor -- one linear pass.
  ! (Ascending order is what makes the single pass exact: any two members of a
  ! run are connected THROUGH the run, and any two members of different runs are
  ! separated by at least one gap > floor at every path between them.)
  !-------------------------------------------------------------------
  pure subroutine gicov_int_degen_blocks(eps, nb, floor, blk)
    implicit none
    integer, intent(in)  :: nb
    real(8), intent(in)  :: eps(nb), floor
    integer, intent(out) :: blk(nb)
    integer :: a
    if (nb < 1) return
    blk(1) = 1
    do a = 2, nb
      if (eps(a) - eps(a - 1) > floor) then
        blk(a) = blk(a - 1) + 1     ! adjacent gap breaks the chain: new block
      else
        blk(a) = blk(a - 1)         ! still connected to the running block
      end if
    end do
  end subroutine gicov_int_degen_blocks

  !-------------------------------------------------------------------
  ! Ordered product of single-step Wilson links into the transport
  !   W(kappa, kappa+n) = V(kappa)^dagger V(kappa+n)   (telescoped),
  ! with a final polar re-projection (round-off / window-leak cleanup, the same
  ! polar_unitary the single-step links themselves were built with).
  !
  ! link(:, :, s) is the FORWARD single-step link of the node s mesh steps
  ! ahead of kappa toward kappa+n, i.e. link(:,:,s) = polar(V(kappa+s-1)^dagger
  ! V(kappa+s)) for s = 1..|n|.  For n > 0 the product runs low-to-high; for
  ! n < 0 the transport toward kappa-|n| is the adjoint chain built from the
  ! links on the kappa-|n| .. kappa side (the caller supplies THOSE links in
  ! link(:,:,1..|n|) already oriented as the backward-side forward links, and
  ! sets back=.true.).  n = 0 returns the identity.
  !-------------------------------------------------------------------
  subroutine gicov_int_chain(link, nb, nabs, back, W)
    implicit none
    integer,    intent(in)  :: nb, nabs
    complex(8), intent(in)  :: link(nb, nb, nabs)
    logical,    intent(in)  :: back
    complex(8), intent(out) :: W(nb, nb)
    complex(8) :: acc(nb, nb), tmp(nb, nb)
    real(8)    :: sigma_min
    integer    :: s, i, ierr
    acc = (0d0, 0d0)
    do i = 1, nb
      acc(i, i) = (1d0, 0d0)
    end do
    if (nabs >= 1) then
      if (.not. back) then
        do s = 1, nabs
          tmp = matmul(acc, link(:, :, s))
          acc = tmp
        end do
      else
        ! backward transport: adjoint links, applied from the far node inward
        do s = nabs, 1, -1
          tmp = matmul(acc, conjg(transpose(link(:, :, s))))
          acc = tmp
        end do
      end if
    end if
    call polar_unitary(acc, nb, W, sigma_min, ierr)
    if (ierr /= 0) then
      ! near-singular telescoped overlap: fall back to the raw (still unitary
      ! to round-off) product so a transient window leak does not abort the run
      W = acc
    end if
  end subroutine gicov_int_chain

  !-------------------------------------------------------------------
  ! Transport a band-frame operator O at the remote node kappa+n INTO the
  ! kappa frame:  Y = W O W^dagger,  W = W(kappa, kappa+n).  Orientation fixed
  ! by matching the exact analytic co-moving operator off-axis (the on-axis row
  ! is blind to the orientation, so the sandwich side is validated off-axis):
  ! W O W^dagger, NOT W^dagger O W.
  !-------------------------------------------------------------------
  pure subroutine gicov_int_transport_op(W, O, nb, Y)
    implicit none
    integer,    intent(in)  :: nb
    complex(8), intent(in)  :: W(nb, nb), O(nb, nb)
    complex(8), intent(out) :: Y(nb, nb)
    Y = matmul(matmul(W, O), conjg(transpose(W)))
  end subroutine gicov_int_transport_op

  !-------------------------------------------------------------------
  ! Linear interpolation of a BOUNDED transported operator between the two
  ! mesh nodes bracketing x = kappa - a:  Y = (1-frac) Ylo + frac Yhi.  The
  ! interpolated object is bounded (||Y|| = ||O||), so the error is O(dk^2) on
  ! a smooth field -- never the 1/dk amplification of interpolating the
  ! singular connection or the finite-difference rate.
  !-------------------------------------------------------------------
  pure subroutine gicov_int_interp(Ylo, Yhi, frac, nb, Y)
    implicit none
    integer,    intent(in)  :: nb
    complex(8), intent(in)  :: Ylo(nb, nb), Yhi(nb, nb)
    real(8),    intent(in)  :: frac
    complex(8), intent(out) :: Y(nb, nb)
    Y = (1d0 - frac) * Ylo + frac * Yhi
  end subroutine gicov_int_interp

  !-------------------------------------------------------------------
  ! Midpoint-frozen EXACT-exponential update of one k-block over a step dt,
  ! plus the Delta-omega-gated phenomenological T2 dephasing, all in the
  ! instantaneous eigenbasis of the co-moving Hamiltonian H~ evaluated at the
  ! step midpoint (H = H~(kappa, t+dt/2), supplied by the caller):
  !
  !   H = P diag(eps) P^dagger        (Hermitian eigensolve, zheev)
  !   R = P^dagger rho~ P
  !   R_ab <- exp(-i (eps_a - eps_b) dt)
  !            * exp( -gamma * g(eps_a - eps_b) * dt ) * R_ab
  !   rho~_out = P R P^dagger
  !
  ! gamma = 1/T2, g = gicov_int_gate_weight on the MOVING gap (never
  ! interpolated).  Diagonal a=b: gap 0 <= floor -> g=0 and phase 0, so
  ! populations are left EXACTLY invariant (trace- and Ne-exact by
  ! construction, independent of mesh); an exactly degenerate off-diagonal pair
  ! (|eps_a-eps_b| <= floor) is likewise never dephased (gauge covariance).
  ! The whole update is unitary-similarity + a Hermitian-preserving Hadamard
  ! decay, so Hermiticity is machine-exact.
  !
  ! Work arrays (eps, P, R, cwork, rwork) are passed in so the routine performs
  ! NO heap allocation and NO per-thread automatic array: the caller allocates
  ! one work bank per OpenMP thread OUTSIDE the k-loop and passes the thread's
  ! slice, per the frtpx heap-outside-OMP discipline.  lcwork must be >= 2*nb-1.
  !
  ! ierr /= 0 means the instantaneous eigenproblem BROKE and rho_out is NOT to
  ! be trusted: the caller MUST abort the whole communicator (see
  ! gicov_int_eig_status).  It is deliberately NOT absorbed here -- the former
  ! "leave the block unchanged" fallback preserved the trace exactly, so a
  ! corrupted k-block slipped silently past the Ne/trace monitor.
  !-------------------------------------------------------------------
  subroutine gicov_int_step_k(H, rho, nb, dt, gamma, shape, theta, width, floor, &
                            & eps, P, R, cwork, lcwork, rwork, blk, rho_out, ierr)
    implicit none
    integer,      intent(in)    :: nb, lcwork
    complex(8),   intent(in)    :: H(nb, nb), rho(nb, nb)
    real(8),      intent(in)    :: dt, gamma, theta, width, floor
    character(*), intent(in)    :: shape
    real(8),      intent(inout) :: eps(nb), rwork(*)
    complex(8),   intent(inout) :: P(nb, nb), R(nb, nb), cwork(*)
    integer,      intent(inout) :: blk(nb)
    complex(8),   intent(out)   :: rho_out(nb, nb)
    integer,      intent(out)   :: ierr
    integer    :: a, b, info
    real(8)    :: dw, g, decay
    complex(8) :: fac
    P = H
    call zheev('V', 'U', nb, P, nb, eps, cwork, lcwork, rwork, info)
    ierr = gicov_int_eig_status(info, eps, nb)
    if (ierr /= 0) then
      ! REPORT, never absorb.  rho_out is defined only so the dummy is not left
      ! undefined; the caller aborts on ierr /= 0 and must not read it.
      rho_out = rho
      return
    end if
    ! connected near-degenerate manifolds of the INSTANTANEOUS spectrum: the
    ! dephasing rate is a function of the BLOCK PAIR, never of the raw gap, so
    ! that no coherence inside a physically indistinguishable manifold decays.
    call gicov_int_degen_blocks(eps, nb, floor, blk)
    R = matmul(matmul(conjg(transpose(P)), rho), P)
    do b = 1, nb
      do a = 1, nb
        dw = eps(a) - eps(b)
        if (blk(a) == blk(b)) then
          g = 0d0          ! same connected manifold (incl. a==b): NEVER dephase
        else
          g = gicov_int_gate_weight(dw, shape, theta, width, floor)
        end if
        decay = exp( -gamma * g * dt )
        fac = exp( cmplx(0d0, -dw * dt, 8) ) * decay
        R(a, b) = fac * R(a, b)
      end do
    end do
    rho_out = matmul(matmul(P, R), conjg(transpose(P)))
  end subroutine gicov_int_step_k

  !-------------------------------------------------------------------
  ! Current contribution of one k-block: j_i = Re Tr[ rho~ v~_i ], i=1..3, with
  ! v~_i the cached transported velocity interpolated to x = kappa - a (same
  ! frame/label as H~).  The caller multiplies by the k-weight and divides by
  ! the cell volume and the weight sum.  No D.A subtraction is applied (that is
  ! a velocity-gauge concept; the length-gauge current is Tr(v rho)).
  !-------------------------------------------------------------------
  pure subroutine gicov_int_current_k(rho, v, nb, jk)
    implicit none
    integer,    intent(in)  :: nb
    complex(8), intent(in)  :: rho(nb, nb), v(nb, nb, 3)
    real(8),    intent(out) :: jk(3)
    integer    :: a, b, i
    complex(8) :: tr
    do i = 1, 3
      tr = (0d0, 0d0)
      do b = 1, nb
        do a = 1, nb
          tr = tr + rho(a, b) * v(b, a, i)
        end do
      end do
      jk(i) = real(tr, 8)
    end do
  end subroutine gicov_int_current_k

  !-------------------------------------------------------------------
  ! Instantaneous-eigenbasis populations of one k-block:
  !   n_a(x) = (P^dagger rho~ P)_aa,   H~(x) = P diag(eps) P^dagger.
  ! The occupation of a co-moving state is its projection onto the INSTANTANEOUS
  ! eigenvector of the moving Hamiltonian, NOT diag(rho~) in the frozen band
  ! basis (which is basis-dependent under transport).  Degenerate members share
  ! a block; the caller sums n_a over each degenerate block for a basis-
  ! invariant readout.  Work arrays passed in (no heap / no automatic array).
  !
  ! Within a degenerate block the eigensolver's choice of columns is ARBITRARY,
  ! so the individual n_a are not observable -- only the block total is.  The
  ! readout is therefore symmetrized over each connected block (blk, the same
  ! closure the T2 gate uses): every member reports blockTotal/blockSize.  That
  ! keeps the column basis-invariant AND preserves the sum rule (summing all
  ! bands still gives the total electron count), unlike reporting the raw
  ! diagonal of the instantaneous eigenbasis.
  !
  ! ierr /= 0 => the eigenproblem broke and nocc is meaningless; the caller MUST
  ! abort.  The former fallback -- reading diag(rho~) in the frozen band basis --
  ! is FORBIDDEN: it is exactly the transport-non-invariant quantity this routine
  ! exists to replace, and it silently returns a plausible, trace-correct number.
  !-------------------------------------------------------------------
  subroutine gicov_int_occupation_k(H, rho, nb, floor, eps, P, R, cwork, lcwork, &
                                  & rwork, nocc, blk, ierr)
    implicit none
    integer,    intent(in)    :: nb, lcwork
    complex(8), intent(in)    :: H(nb, nb), rho(nb, nb)
    real(8),    intent(in)    :: floor
    real(8),    intent(inout) :: eps(nb), rwork(*)
    complex(8), intent(inout) :: P(nb, nb), R(nb, nb), cwork(*)
    real(8),    intent(out)   :: nocc(nb)
    integer,    intent(inout) :: blk(nb)
    integer,    intent(out)   :: ierr
    integer :: a, b, info, m
    real(8) :: s
    P = H
    call zheev('V', 'U', nb, P, nb, eps, cwork, lcwork, rwork, info)
    ierr = gicov_int_eig_status(info, eps, nb)
    if (ierr /= 0) then
      nocc(1:nb) = 0d0        ! NOT diag(rho~): the caller aborts on ierr /= 0
      return
    end if
    R = matmul(matmul(conjg(transpose(P)), rho), P)
    do a = 1, nb
      nocc(a) = real(R(a, a), 8)
    end do
    ! symmetrize over each connected degenerate block (basis-invariant readout)
    call gicov_int_degen_blocks(eps, nb, floor, blk)
    a = 1
    do while (a <= nb)
      b = a
      do while (b < nb)
        if (blk(b + 1) /= blk(a)) exit
        b = b + 1
      end do
      m = b - a + 1
      if (m > 1) then
        s = sum(nocc(a:b)) / real(m, 8)
        nocc(a:b) = s
      end if
      a = b + 1
    end do
  end subroutine gicov_int_occupation_k

  !-------------------------------------------------------------------
  ! Status of one instantaneous Hermitian eigensolve.  0 = usable.
  !
  ! Checking LAPACK's info alone is NOT sufficient: zheev returns info=0 for a
  ! NaN-poisoned Hermitian matrix (verified on Accelerate), so a corrupted H~
  ! would yield NaN eigenvalues, a NaN density, and -- because the unitary-
  ! similarity update is trace-preserving in form -- no complaint from the
  ! Ne/trace monitor.  Reject non-finite eigenvalues explicitly.
  !   info /= 0        -> ierr = info   (LAPACK-reported failure)
  !   non-finite eps   -> ierr = -999   (silent corruption)
  !-------------------------------------------------------------------
  pure integer function gicov_int_eig_status(info, eps, nb) result(ierr)
    implicit none
    integer, intent(in) :: info, nb
    real(8), intent(in) :: eps(nb)
    integer :: a
    ierr = info
    if (ierr /= 0) return
    do a = 1, nb
      if (.not. ieee_is_finite(eps(a))) then
        ierr = -999
        return
      end if
    end do
  end function gicov_int_eig_status

end module gicov_integral_ssbe
