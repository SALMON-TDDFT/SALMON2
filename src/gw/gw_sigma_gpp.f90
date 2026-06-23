!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Dynamic one-shot screened-Coulomb self-energy via the generalized
! plasmon-pole (GPP) model -- the proof-of-concept rung (ladder rung 3).
!
! This rung restores the FREQUENCY dependence that the static COHSEX rung threw
! away.  The full GW correlation self-energy is the convolution
!   Sigma_c(E) = (i/2pi) integral dw  G(E-w) [ W(w) - v ],
! where the screened interaction W(w) carries the dynamics of the electronic
! response.  Evaluating that convolution exactly needs W on the whole frequency
! axis; the generalized plasmon-pole model instead REPRESENTS the inverse
! dielectric matrix by a single effective plasmon mode per (G,G') channel,
!
!   eps^{-1}_{GG'}(q,w) = delta_{GG'}
!         + Omega^2_{GG'}(q) / ( w^2 - ( wt_{GG'}(q) - i eta )^2 ),
!
! whose strength Omega^2 and pole position wt are FIXED by two exact
! constraints: the static limit eps^{-1}(q,0) (the COHSEX input, from
! calc_epsinv) and the f-sum rule on the density-density response.  With one
! pole the frequency convolution is done analytically and yields a closed form
! for Sigma_c(E) with explicit, well-defined E dependence -- hence a nontrivial
! renormalization Z and a one-shot G0W0 quasiparticle shift that, unlike static
! COHSEX, does NOT over-open the gap.
!
! Physics -- the GPP parameters (atomic units)
! --------------------------------------------
! Total-density Fourier components (the f-sum-rule input):
!   rho(G'') = (1/Omega) integral_cell n(r) e^{-i G''.r} d3r
!            = (1/Ngrid) sum_{grid r} n(r) e^{-i G''.r},
! built by a direct discrete Fourier sum of the TOTAL ground-state charge
! density n(r) = rho%f (summed over spin) on the same real-space slab the
! orbitals live on; hvol = Omega/Ngrid folds the cell integral into the average
! (the same (1/Omega) integral -> (1/Ngrid) sum bookkeeping as calc_mtxel and
! calc_vxc_expect).  The G''=0 component is the average density rho(0) =
! N_elec/Omega, and the plasma frequency squared is wp^2 = 4 pi rho(0).
!
! Sum-rule oscillator strength (the standard generalized-plasmon-pole form):
!   Omega^2_{GG'}(q) = wp^2 * [ (q+G).(q+G') / |q+G|^2 ] * [ rho(G-G') / rho(0) ].
! The bracket [(q+G).(q+G')/|q+G|^2] is the longitudinal f-sum-rule projector
! (it reduces to 1 for G=G'); rho(G-G')/rho(0) carries the local-field
! (off-diagonal) structure.
!
! Pole position from the static limit.  Setting w=0 in the GPP ansatz,
!   eps^{-1}_{GG'}(q,0) = delta_{GG'} - Omega^2_{GG'}(q) / wt^2_{GG'}(q),
! so the effective mode frequency is fixed WITHOUT any empty-state sum:
!   wt^2_{GG'}(q) = Omega^2_{GG'}(q) / ( delta_{GG'} - eps^{-1}_{GG'}(q,0) ).
! Unphysical poles.  A genuine plasmon needs wt^2 > 0, i.e. Omega^2 and
! ( delta - eps^{-1} ) must share sign and the denominator must not vanish.
! Where  ( delta - eps^{-1} ) <= 0  (or the resulting wt^2 <= 0) the channel has
! no physical pole; following the standard plasmon-pole handling that (G,G')
! contribution is ZEROED and COUNTED, and the skipped fraction is reported (it sits in the
! high-|G-G'| / weakly screened corner and is small).
!
! Physics -- the correlation self-energy (single-pole screened-Coulomb form)
! ----------------------------------------------------------------
! Doing the single-pole frequency convolution analytically gives, for state
! (n,k) at probe energy E, a sum over ALL intermediate bands n' (occupied AND
! unoccupied) at the folded point k-q:
!
!   Sigma_c(nk;E) = (1/(N_k Omega)) sum_{n'} sum_q sum_{G,G'}
!        M*_{nn',G}(k,q) M_{nn',G'}(k,q) v(q+G')
!        ( Omega^2_{GG'}(q) / (2 wt_{GG'}(q)) ) *
!        [ (1 - f_{n',k-q}) / ( E - eps_{n',k-q} - wt_{GG'}(q) + i eta )
!        +       f_{n',k-q}  / ( E - eps_{n',k-q} + wt_{GG'}(q) - i eta ) ].
!
! Occupied/unoccupied pole signs.  The two terms are the two halves of the
! contour integral: the (1-f) term is the EMPTY-state branch (the pole sits at
! E = eps + wt, above the state), the f term is the OCCUPIED branch (pole at
! E = eps - wt, below).  For a gap state probed at E = eps^KS this combination is
! Coulomb-hole dominated and gives Sigma_c < 0, with |Sigma_c| smaller than the
! static COH because the dynamics partially cancels it -- which is exactly why
! one-shot GPP G0W0 lands BELOW static COHSEX.
!
! M-element bridge -- identical to calc_sigma_x
! --------------------------------------------
! M_{nn',G}(k,q) = < n,k | e^{i(q+G).r} | n',k-q >, with n' now over the FULL
! band set at k-q (Sigma_x kept only occupied v; here the band sum is complete,
! occ + unocc).  We reuse the exact same machinery:
!   find_kpq(system, ik, -q, ikm, g0vec)              => k - q = k_{ikm} + G0,
!   imap(ig) = gindex_of( gvec(:,ig) - g0vec )         (the +(G-G0) convention),
!   mblk(ig', bra@ikm, ket@ik) = calc_mtxel(... ik, ikm ...),
!   M*_{nn',G(ig)} = mblk( imap(ig), n', n ),   M_{nn',G'(jg)} = conjg(that).
! Entries whose (G-G0) leaves the cutoff are unavailable and contribute zero
! (boundary truncation; harmless where v(q+G') is largest).
!
! Quasiparticle solve -- the linearized one-shot G0W0 with Z
! ----------------------------------------------------------
! Sigma_c is frequency dependent, so the renormalization is NOT unity:
!   Z(nk) = ( 1 - d ReSigma_c/dE |_{E=eps^KS} )^{-1},   Z in [0.7,0.9] expected,
! computed by a central finite difference Sigma_c(eps^KS +- dE).  The QP energy
! (combined with the bare exchange Sigma_x from calc_sigma_x in gw_main) is
!   eps^QP(nk) = eps^KS(nk) + Z(nk) ( Sigma_x(nk) + ReSigma_c(nk;eps^KS)
!                                     - <V_xc>(nk) ).
! Only the REAL part of Sigma_c enters the QP energy; the i eta (a small
! broadening, eta_ev eV in a.u.) keeps the on-shell denominators finite and is
! discarded by taking Re.
!
! Prefactor 1/(N_k Omega) -- identical reasoning to calc_sigma_x / cohsex
! ----------------------------------------------------------------------
! The BZ integral (1/(2 pi)^3) int_BZ dq -> (1/(N_k Omega)) sum_q, and v(q+G')
! from build_vcoul carries no 1/Omega, so the cell volume Omega = system%det_a
! enters the prefactor here.  Omega^2_{GG'} (the oscillator strength) is itself
! O(wp^2) with no extra cell factor (rho(G-G')/rho(0) and the f-sum projector are
! both dimensionless), so it does not disturb the 1/(N_k Omega) bookkeeping.
!
! Parallel layout
! ---------------
! Launched with nproc_k=1 (all k local); calc_epsinv and calc_mtxel are
! collective (comm_summation over icomm_rko) and return identical global results
! on every rank, and the density Fourier sum is reduced over icomm_rko, so the
! whole Sigma_c accumulation is a deterministic local reduction afterwards.
!
! No proper nouns appear in this file.
!
module gw_sigma_gpp_sub
  implicit none
  private

  public :: calc_sigma_gpp
  public :: calc_sigma_gpp_qcache
  public :: calc_sigma_gpp_sym

  real(8), parameter :: occ_thr     = 1.0d-6   ! occupied if rocc > occ_thr
  real(8), parameter :: g_match_tol = 1.0d-6   ! |G|-match tolerance (bohr^-1)^2
  integer, parameter :: nsub_head   = 10       ! mini-BZ sub-sampling for v head
  integer, parameter :: ngmax_l     = 200000   ! capacity of the larger G-set
  real(8), parameter :: gset_factor = 4.0d0    ! larger cutoff = 4 * eps cutoff
  real(8), parameter :: eta_ev      = 0.1d0    ! pole broadening eta (eV)
  real(8), parameter :: dE_ev       = 0.1d0    ! finite-diff step dE for Z (eV)
  real(8), parameter :: ha2ev       = 27.21138505d0  ! eV per Hartree
  real(8), parameter :: wt2_tiny    = 1.0d-12  ! |delta-eps^{-1}| floor (unphysical)

contains

  ! --------------------------------------------------------------------------
  ! gindex_of
  !
  ! List index jg with gv(:,jg) == gtarget (within tolerance), 0 if gtarget is
  ! outside the list (cutoff sphere).  Local copy (kept private so the module is
  ! self-contained); identical role to the helper in gw_epsilon / gw_sigma_x /
  ! gw_sigma_cohsex.
  ! --------------------------------------------------------------------------
  integer function gindex_of(ngl, gv, gtarget) result(jg)
    implicit none
    integer, intent(in) :: ngl
    real(8), intent(in) :: gv(3,ngl)
    real(8), intent(in) :: gtarget(3)
    integer :: i
    real(8) :: d1, d2, d3, d2sum
    jg = 0
    do i = 1, ngl
      d1 = gv(1,i) - gtarget(1)
      d2 = gv(2,i) - gtarget(2)
      d3 = gv(3,i) - gtarget(3)
      d2sum = d1*d1 + d2*d2 + d3*d3
      if (d2sum < g_match_tol) then
        jg = i
        return
      end if
    end do
  end function gindex_of

  ! --------------------------------------------------------------------------
  ! build_density_ft
  !
  ! Total-density Fourier components rho(G'') = (1/Ngrid) sum_grid n(r)
  ! e^{-i G''.r} on a G''-set gvl(:,1..ngl), by a direct discrete Fourier sum of
  ! the TOTAL ground-state charge density n(r) = rho%f over the local mg slab,
  ! reduced over icomm_rko.  rho%f and the orbitals share the SAME global grid
  ! index space (mg%is..mg%ie) and lg%coordinate gives the Cartesian r, exactly
  ! as in calc_mtxel / calc_vxc_expect.  (1/Ngrid) = hvol/Omega turns the cell
  ! integral into the (1/Omega) average that defines rho(G'').  rhoft(1) (the
  ! G''=0 entry; build_gvectors puts G=0 first) equals the average density
  ! N_elec/Omega and pins wp^2 = 4 pi rhoft(1).
  ! --------------------------------------------------------------------------
  subroutine build_density_ft(system, info, mg, lg, rho, gvl, ngl, rhoft, local_only)
    use structures,    only: s_dft_system, s_parallel_info, s_rgrid, s_scalar
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_scalar),        intent(in)  :: rho           ! total density rho%f
    integer,               intent(in)  :: ngl
    real(8),               intent(in)  :: gvl(3,ngl)
    complex(8),            intent(out) :: rhoft(ngl)
    ! With nproc_rgrid=1 each rank already holds the whole grid, so the local sum
    ! IS the full density-FT; local_only skips the assembling collective (which,
    ! over a k-replicated communicator, would multiply-count the identical sum).
    logical, optional,     intent(in)  :: local_only

    complex(8), allocatable :: rloc(:)
    complex(8) :: zacc
    real(8) :: hvol, omega, rx, ry, rz, phase, nr, inv_ngrid
    integer :: ig, ix, iy, iz
    logical :: ll

    ll = .false.
    if (present(local_only)) ll = local_only

    omega = abs(system%det_a)
    hvol  = system%hvol
    inv_ngrid = hvol / omega          ! = 1/Ngrid : cell integral -> (1/Omega) avg

    allocate(rloc(ngl))
    rloc(:) = (0.0d0, 0.0d0)

    ! Direct Fourier sum of the local-slab total density onto each requested G''.
    ! rho%f is real; the kernel e^{-i G''.r} makes the running sum complex.
!$omp parallel do private(ig,iz,iy,ix,zacc,rx,ry,rz,nr,phase)
    do ig = 1, ngl
      zacc = (0.0d0, 0.0d0)
      do iz = mg%is(3), mg%ie(3)
        rz = lg%coordinate(iz,3)
      do iy = mg%is(2), mg%ie(2)
        ry = lg%coordinate(iy,2)
      do ix = mg%is(1), mg%ie(1)
        rx = lg%coordinate(ix,1)
        nr = rho%f(ix,iy,iz)
        phase = -( gvl(1,ig)*rx + gvl(2,ig)*ry + gvl(3,ig)*rz )
        zacc = zacc + nr * cmplx(cos(phase), sin(phase), 8)
      end do
      end do
      end do
      rloc(ig) = zacc * inv_ngrid
    end do
!$omp end parallel do

    ! Assemble across the r-space (and any k/orbital) partitions (same
    ! communicator calc_mtxel / calc_vxc_expect use).
    if (ll) then
      rhoft(:) = rloc(:)
    else
      call comm_summation(rloc, rhoft, ngl, info%icomm_rko)
    end if

    deallocate(rloc)

  end subroutine build_density_ft

  ! --------------------------------------------------------------------------
  ! sigma_c_at_energy
  !
  ! Inner kernel: given the precomputed, q-resolved GPP weight wgpp(ig,jg) =
  ! v(q+G'(jg)) * Omega^2_{GG'}/(2 wt_{GG'}), the per-channel mode frequency
  ! wt(ig,jg), a per-channel physical-pole mask phys(ig,jg), the M-element gather
  ! msig(ig, n') = M*_{nn',G(ig)} for a FIXED output band n (and the band
  ! energies/occupations at k-q), accumulate the (complex) correlation
  ! self-energy contribution at probe energy E into ssig:
  !
  !   ssig += sum_{n'} sum_{G,G'} msig(G,n') wgpp(G,G') conjg(msig(G',n'))
  !           [ (1-f)/(E - e_{n'} - wt + i eta) + f/(E - e_{n'} + wt - i eta) ].
  !
  ! Note v(q+G') and the M*_G M_{G'} contraction live in wgpp/msig; the band-
  ! and energy-dependent pole bracket is applied per (n',G,G').  The 1/(N_k
  ! Omega) prefactor is applied by the caller.  Channels with phys=.false. are
  ! skipped (their wgpp is already zero, but the mask guards the bracket too).
  ! --------------------------------------------------------------------------
  subroutine sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, &
                               eprobe, eta, ssig)
    implicit none
    integer,    intent(in)  :: ng, no
    complex(8), intent(in)  :: wgpp(ng,ng)
    real(8),    intent(in)  :: wt(ng,ng)
    logical,    intent(in)  :: phys(ng,ng)
    complex(8), intent(in)  :: msig(ng, no)      ! M*_{nn',G(ig)} for fixed n
    real(8),    intent(in)  :: eband(no)         ! eps_{n',k-q}
    real(8),    intent(in)  :: focc(no)          ! occupation in [0,1]
    real(8),    intent(in)  :: eprobe            ! probe energy E
    real(8),    intent(in)  :: eta               ! broadening (a.u.)
    complex(8), intent(out) :: ssig

    integer :: ig, jg, inp
    complex(8) :: zeta, mji, racc, brk, ww
    real(8) :: en, f1, fo, de

    zeta = cmplx(0.0d0, eta, 8)
    ssig = (0.0d0, 0.0d0)

    do inp = 1, no
      en = eband(inp)
      fo = focc(inp)
      if (fo < 0.0d0) fo = 0.0d0
      if (fo > 1.0d0) fo = 1.0d0
      f1 = 1.0d0 - fo
      de = eprobe - en
      ! Contract the double-G GPP weight with M*_G M_{G'} for this n', applying
      ! the per-channel plasmon-pole denominators.  brk depends on (G,G') only
      ! through wt(ig,jg); we fold it into the inner sum.
      racc = (0.0d0, 0.0d0)
      do jg = 1, ng
        mji = conjg(msig(jg,inp))            ! M_{nn',G'(jg)}
        if (mji == (0.0d0,0.0d0)) cycle
        do ig = 1, ng
          if (.not. phys(ig,jg)) cycle       ! unphysical pole: zeroed + counted
          ww  = wgpp(ig,jg)
          if (ww == (0.0d0,0.0d0)) cycle
          ! plasmon-pole bracket for channel (ig,jg):
          !   (1-f)/(de - wt + i eta) + f/(de + wt - i eta)
          brk = f1 / ( cmplx(de - wt(ig,jg), 0.0d0, 8) + zeta ) &
              + fo / ( cmplx(de + wt(ig,jg), 0.0d0, 8) - zeta )
          racc = racc + msig(ig,inp) * ww * mji * brk
        end do
      end do
      ssig = ssig + racc
    end do

  end subroutine sigma_c_at_energy

  ! --------------------------------------------------------------------------
  ! calc_sigma_gpp
  !
  ! Diagonal dynamic correlation self-energy Re Sigma_c(nk; eps^KS) and the
  ! renormalization Z(nk) for a band window [nb_lo,nb_hi] and all k-points, spin
  ! ispin (a.u.).  The bare exchange Sigma_x comes from calc_sigma_x and is
  ! combined with these in gw_main.  Internally: build rho(G'') on a larger
  ! G-set, then per q call calc_epsinv for eps^{-1}(q,0), form Omega^2 and wt,
  ! gather the M block, and accumulate Sigma_c at E = eps^KS and E +- dE (for Z).
  !
  ! q-set: q = vec_k(:,iq) - vec_k(:,ik), iq=1..nk (the BZ-mesh differences,
  ! including q=0 with its mini-BZ head), matching calc_sigma_x / calc_sigma_cohsex.
  !
  ! Arguments:
  !   system   [in]  -- vec_k, rocc, det_a (Omega), nk, no, nspin, hvol, lattice.
  !   info     [in]  -- parallel ownership / communicator (to calc_mtxel/epsinv).
  !   mg, lg   [in]  -- real-space grid descriptors.
  !   spsi     [in]  -- cell-periodic orbitals.
  !   esp(no,nk,nspin)[in] -- KS eigenvalues (a.u.); pass energy%esp.
  !   rho      [in]  -- total ground-state charge density rho%f (s_scalar).
  !   gvec(3,ng)[in] -- Cartesian G vectors (bohr^-1), the eps G-set.
  !   gg(ng)   [in]  -- |G|^2 (bohr^-2), the eps G-set.
  !   ng       [in]  -- number of eps G vectors (eps^{-1} matrix dimension).
  !   ispin    [in]  -- spin channel (1 for the non-magnetic case).
  !   nb_lo,nb_hi[in]-- output band window (states n).
  !   sigc_re(nb_lo:nb_hi, nk) [out] -- Re Sigma_c(nk; eps^KS) (a.u.).
  !   zfac   (nb_lo:nb_hi, nk) [out] -- renormalization Z(nk).
  !   skip_frac [out,opt] -- fraction of (G,G') channels with no physical pole
  !                          (delta-eps^{-1}<=0 or wt^2<=0), averaged over the
  !                          representable q-points.
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_gpp(system, info, mg, lg, spsi, esp, rho, gvec, gg, ng, &
                            ispin, nb_lo, nb_hi, sigc_re, zfac, skip_frac, local_only)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, s_scalar
    use gw_mtxel_sub,   only: calc_mtxel
    use gw_coulomb_sub, only: build_vcoul, build_gvectors
    use gw_epsilon_sub, only: find_kpq, calc_epsinv
    use communication,  only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_orbital),       intent(in)  :: spsi
    real(8),               intent(in)  :: esp(system%no, system%nk, system%nspin)
    type(s_scalar),        intent(in)  :: rho
    integer,               intent(in)  :: ng
    real(8),               intent(in)  :: gvec(3,ng)
    real(8),               intent(in)  :: gg(ng)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nb_lo, nb_hi
    real(8),               intent(out) :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out) :: zfac   (nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: skip_frac
    ! Node-parallel mode: replicated orbitals + the output-state k-loop split over
    ! info%icomm_k, with one reduction of the per-k results at the end.  Absent =
    ! the original collective path (all k per rank).
    logical,    optional,  intent(in)  :: local_only

    complex(8), allocatable :: epsinv(:,:)       ! eps^{-1}_{GG'}(q,0)
    complex(8), allocatable :: mblk(:,:,:)       ! < u_{.,ikm}|e^{-iG.r}|u_{.,ik} >
    complex(8), allocatable :: msig(:,:)         ! M*_{nn',G(ig)} for fixed n
    complex(8), allocatable :: wgpp(:,:)         ! v(q+G') Omega^2/(2 wt)  per (G,G')
    complex(8), allocatable :: rhoft(:)          ! rho(G'') on the larger set
    real(8),    allocatable :: vcoul(:)
    real(8),    allocatable :: wt(:,:)           ! mode frequency wt_{GG'} (>0)
    logical,    allocatable :: phys(:,:)         ! physical-pole mask per (G,G')
    integer,    allocatable :: imap(:)           ! eps-set G -> mtxel column
    integer,    allocatable :: idiff(:,:)        ! larger-set index of G(ig)-G(jg)
    real(8),    allocatable :: gvecl(:,:), ggl(:)
    real(8),    allocatable :: eband(:), focc(:)

    real(8),    allocatable :: e0(:,:)           ! eps^KS over the window (probe)
    complex(8), allocatable :: sc0(:), scp(:), scm(:)  ! Sigma_c at E0, E0+dE, E0-dE
    real(8),    allocatable :: sigc_g(:,:), zfac_g(:,:)  ! reduction buffers

    integer :: no, nk, nq, ik, iq, ikm, inp, in, ig, jg, kg, ik_lo, ik_hi
    integer :: ngl, nchan_phys_tot, nchan_unphys_tot
    integer :: ncnt(2), ncnt_g(2)
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, ecut_l, gmax2, fourpi, pi
    real(8) :: rho0, wp2, qgi(3), qg2i, qgj(3), proj, om2, denom, wt2
    real(8) :: eta_au, de_au, e0nk, dsig
    complex(8) :: rhod, s0, sp, sm
    logical   :: q_ok
    logical   :: ll

    ll = .false.
    if (present(local_only)) ll = local_only

    no    = system%no
    nk    = system%nk
    nq    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)  ! 1/(N_k*Omega): v(q+G)=4pi/|q+G|^2 carries
                                        ! no 1/Omega, so the cell volume enters here
    pi     = acos(-1.0d0)
    fourpi = 4.0d0 * pi
    eta_au = eta_ev / ha2ev
    de_au  = dE_ev  / ha2ev

    sigc_re(:,:) = 0.0d0
    ! In the parallel mode zfac is summed over icomm_k (each k filled by its sole
    ! owner), so it must start at zero; the serial path keeps the 1.0 default.
    if (ll) then
      zfac(:,:) = 0.0d0
    else
      zfac(:,:) = 1.0d0
    end if

    ! ---- larger G-set for rho(G-G') (the off-diagonal density structure) ----
    ! G and G' each range over the eps G-list, so the argument G-G' can reach
    ! |G-G'|^2 up to ~4 (eps cutoff) -- outside the eps list.  Build rho(G'') on
    ! a SEPARATE, LARGER set (cutoff = gset_factor * max(gg)); a G-G' that still
    ! falls outside it marks the channel unphysical (skipped + counted) since its
    ! oscillator strength is unavailable.
    gmax2 = 0.0d0
    do ig = 1, ng
      if (gg(ig) > gmax2) gmax2 = gg(ig)
    end do
    ecut_l = gset_factor * gmax2

    allocate(gvecl(3, ngmax_l), ggl(ngmax_l))
    call build_gvectors(system%primitive_b, ecut_l, ngmax_l, ngl, gvecl, ggl)

    allocate(rhoft(ngl))
    call build_density_ft(system, info, mg, lg, rho, gvecl, ngl, rhoft, local_only=ll)

    ! Average density and plasma frequency.  rhoft(1) is the G''=0 component
    ! (build_gvectors lists G=0 first) = N_elec/Omega; it must be real & positive.
    rho0 = dble(rhoft(1))
    wp2  = fourpi * rho0

    allocate(epsinv(ng,ng), vcoul(ng), imap(ng))
    allocate(wgpp(ng,ng), wt(ng,ng), phys(ng,ng), idiff(ng,ng))
    allocate(mblk(ng, no, no), msig(ng, no))
    allocate(eband(no), focc(no))
    allocate(e0(nb_lo:nb_hi, nk))
    allocate(sc0(nb_lo:nb_hi), scp(nb_lo:nb_hi), scm(nb_lo:nb_hi))

    ! ---- larger-set lookup of G(ig)-G(jg) (the rho argument).  rho(G''=G-G')
    ! is rhoft( gindexL(G(ig)-G(jg)) ); precompute once (q-independent: depends
    ! only on the fixed eps G-list).  0 marks a difference outside the larger set.
    do jg = 1, ng
      do ig = 1, ng
        gtarget(1) = gvec(1,ig) - gvec(1,jg)
        gtarget(2) = gvec(2,ig) - gvec(2,jg)
        gtarget(3) = gvec(3,ig) - gvec(3,jg)
        idiff(ig,jg) = gindex_of(ngl, gvecl, gtarget)
      end do
    end do

    nchan_phys_tot   = 0
    nchan_unphys_tot = 0

    ! Cache the probe (KS) energies over the band window for every k.
    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! Output-state k-points owned by this rank (its share when parallel).
    if (ll) then
      ik_lo = info%ik_s
      ik_hi = info%ik_e
    else
      ik_lo = 1
      ik_hi = nk
    end if

    ! ---- output-state loop -------------------------------------------------
    do ik = ik_lo, ik_hi

      sc0(:) = (0.0d0, 0.0d0)
      scp(:) = (0.0d0, 0.0d0)
      scm(:) = (0.0d0, 0.0d0)

      ! ---- BZ q-sum: q = vec_k(:,iq) - vec_k(:,ik) ------------------------
      do iq = 1, nq
        qvec(1:3)  =  system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        mqvec(1:3) = -qvec(1:3)

        ! eps^{-1}_{GG'}(q,0) on the eps G-set (collective unless local_only).
        call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, iq, &
                         qvec, ispin, epsinv, ok=q_ok, local_only=ll)
        if (.not. q_ok) cycle              ! q not representable on this mesh

        ! v(q+G'): bare Coulomb kernel; q->0 head averaged over the mini-BZ.
        call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

        ! ---- GPP parameters per (G,G') for this q ------------------------
        ! Omega^2_{GG'} = wp^2 [ (q+G).(q+G')/|q+G|^2 ] [ rho(G-G')/rho0 ];
        ! wt^2_{GG'} = Omega^2 / ( delta - eps^{-1} ); the channel weight that
        ! enters Sigma_c is v(q+G') Omega^2/(2 wt).  Unphysical poles (denom<=0
        ! or wt^2<=0 or rho(G-G') unavailable) are masked off and counted.
        do ig = 1, ng
          qgi(1) = qvec(1) + gvec(1,ig)
          qgi(2) = qvec(2) + gvec(2,ig)
          qgi(3) = qvec(3) + gvec(3,ig)
          qg2i   = qgi(1)*qgi(1) + qgi(2)*qgi(2) + qgi(3)*qgi(3)
          do jg = 1, ng
            phys(ig,jg) = .false.
            wt(ig,jg)   = 0.0d0
            wgpp(ig,jg) = (0.0d0, 0.0d0)

            if (qg2i < 1.0d-12) cycle      ! |q+G|->0 head: no longitudinal pole
            kg = idiff(ig,jg)
            if (kg == 0) cycle             ! rho(G-G') outside larger set: skip

            qgj(1) = qvec(1) + gvec(1,jg)
            qgj(2) = qvec(2) + gvec(2,jg)
            qgj(3) = qvec(3) + gvec(3,jg)
            proj   = ( qgi(1)*qgj(1) + qgi(2)*qgj(2) + qgi(3)*qgj(3) ) / qg2i

            rhod = rhoft(kg)               ! rho(G-G'); real part carries strength
            om2  = wp2 * proj * ( dble(rhod) / rho0 )

            ! denom = delta_{GG'} - eps^{-1}_{GG'}(q,0)
            denom = -dble(epsinv(ig,jg))
            if (ig == jg) denom = denom + 1.0d0

            ! Physical pole requires wt^2 = Omega^2/denom > 0 with denom away
            ! from zero -> Omega^2 and denom must share sign.
            if (abs(denom) < wt2_tiny) cycle
            wt2 = om2 / denom
            if (wt2 <= 0.0d0) cycle

            phys(ig,jg) = .true.
            wt(ig,jg)   = sqrt(wt2)
            ! channel weight: v(q+G') Omega^2/(2 wt).  v on the COLUMN G' (jg),
            ! consistent with calc_sigma_cohsex (wv(ig,jg)=epsinv(ig,jg) v(jg)).
            wgpp(ig,jg) = cmplx( vcoul(jg) * om2 / (2.0d0 * wt(ig,jg)), 0.0d0, 8 )
          end do
        end do

        ! ---- M block: ket index of k-q (= k+(-q) folded) and umklapp G0 ----
        call find_kpq(system, ik, mqvec, ikm, g0vec)
        if (ikm == 0) cycle                ! no mesh partner (skip this q)

        ! tally the physical/unphysical channels for the skip fraction (only for
        ! q-points that actually enter the Sigma_c sum, i.e. ikm /= 0).
        do jg = 1, ng
          do ig = 1, ng
            if (phys(ig,jg)) then
              nchan_phys_tot = nchan_phys_tot + 1
            else
              nchan_unphys_tot = nchan_unphys_tot + 1
            end if
          end do
        end do

        ! Per-(k,q) eps-set index map: imap(ig) -> calc_mtxel column giving
        ! M*_{nn',G(ig)} (output n at ik, intermediate n' at ikm).  Same +(G-G0)
        ! convention as calc_sigma_x.
        do ig = 1, ng
          gtarget(1) = gvec(1,ig) - g0vec(1)
          gtarget(2) = gvec(2,ig) - g0vec(2)
          gtarget(3) = gvec(3,ig) - g0vec(3)
          imap(ig) = gindex_of(ng, gvec, gtarget)   ! 0 if (G-G0) outside cutoff
        end do

        ! Full M block at (ik,ikm): mblk(ig', bra@ikm, ket@ik) on the eps set.
        call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                        no, no, mblk, local_only=ll)

        ! Intermediate-band energies and occupations at k-q (ALL bands n').
        do inp = 1, no
          eband(inp) = esp(inp,ikm,ispin)
          focc(inp)  = system%rocc(inp,ikm,ispin)
        end do
        ! Normalize occupations to [0,1] (rocc folds spin multiplicity; the GPP
        ! occupied/empty split needs a fractional occupation).  For nspin=1 the
        ! closed-shell rocc is ~2 per occupied state -> divide by the maximum.
        call normalize_occ(no, focc)

        ! ---- band-sum Sigma_c at E0 and E0 +- dE for every output band ----
        do in = nb_lo, nb_hi

          ! gather M*_{nn',G(ig)} = mblk(imap(ig), n', n) for fixed n.
          do inp = 1, no
            do ig = 1, ng
              jg = imap(ig)
              if (jg == 0) then
                msig(ig,inp) = (0.0d0, 0.0d0)   ! (G-G0) outside cutoff -> truncate
              else
                msig(ig,inp) = mblk(jg, inp, in)
              end if
            end do
          end do

          e0nk = e0(in,ik)
          call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, &
                                 e0nk,         eta_au, s0)
          call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, &
                                 e0nk + de_au, eta_au, sp)
          call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, &
                                 e0nk - de_au, eta_au, sm)
          sc0(in) = sc0(in) + s0
          scp(in) = scp(in) + sp
          scm(in) = scm(in) + sm

        end do  ! in

      end do  ! iq

      ! Apply the shared 1/(N_k Omega) factor, take the real part, and form Z by
      ! a central finite difference of Re Sigma_c in E.
      do in = nb_lo, nb_hi
        sigc_re(in,ik) = rnk * dble(sc0(in))
        dsig = rnk * ( dble(scp(in)) - dble(scm(in)) ) / (2.0d0 * de_au)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do

    end do  ! ik

    ! ---- assemble per-k results across the output-state partition ----------
    ! Each rank filled only its own k columns (disjoint), so the sum is exact.
    ! The pole-skip counters are likewise partial and are reduced for the global
    ! fraction below.
    if (ll) then
      allocate(sigc_g(nb_lo:nb_hi, nk), zfac_g(nb_lo:nb_hi, nk))
      call comm_summation(sigc_re, sigc_g, size(sigc_re), info%icomm_k)
      call comm_summation(zfac,    zfac_g, size(zfac),    info%icomm_k)
      sigc_re(:,:) = sigc_g(:,:)
      zfac(:,:)    = zfac_g(:,:)
      deallocate(sigc_g, zfac_g)
      ncnt(1) = nchan_phys_tot
      ncnt(2) = nchan_unphys_tot
      call comm_summation(ncnt, ncnt_g, 2, info%icomm_k)
      nchan_phys_tot   = ncnt_g(1)
      nchan_unphys_tot = ncnt_g(2)
    end if

    ! Channel skip fraction (physical pole missing): averaged over the
    ! representable q-points and all (G,G').
    if (present(skip_frac)) then
      if (nchan_phys_tot + nchan_unphys_tot > 0) then
        skip_frac = dble(nchan_unphys_tot) &
                  / dble(nchan_phys_tot + nchan_unphys_tot)
      else
        skip_frac = 0.0d0
      end if
    end if

    deallocate(epsinv, vcoul, imap, wgpp, wt, phys, idiff)
    deallocate(mblk, msig, rhoft, gvecl, ggl)
    deallocate(eband, focc, e0, sc0, scp, scm)

  end subroutine calc_sigma_gpp

  ! --------------------------------------------------------------------------
  ! calc_sigma_gpp_qcache
  !
  ! Same physics as calc_sigma_gpp, reordered so eps^{-1}(q), v(q) and the GPP
  ! weights {wgpp,wt,phys} -- all functions of the momentum transfer q ALONE --
  ! are built ONCE per distinct q rather than once per (ik,iq) pair.  The
  ! distinct q (the difference set vec_k(iq)-vec_k(ik)) are distributed over
  ! info%icomm_k; the inner per-pair body is byte-for-byte that of
  ! calc_sigma_gpp, so the answer matches it up to reduction order.  Selected by
  ! yn_gw_qcache='y'.  Requires replicated orbitals (spsi spans all k); the inner
  ! kernels run in local mode.  Cost O(nk^2) vs the default O(nk^3).
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_gpp_qcache(system, info, mg, lg, spsi, esp, rho, gvec, gg, &
                                   ng, ispin, nb_lo, nb_hi, sigc_re, zfac, skip_frac, &
                                   do_remainder, nb_sigma, nb_eps)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, s_scalar
    use gw_mtxel_sub,   only: calc_mtxel, build_state_density_ft
    use gw_coulomb_sub, only: build_vcoul, build_gvectors
    use gw_epsilon_sub, only: find_kpq, calc_epsinv
    use communication,  only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_orbital),       intent(in)  :: spsi
    real(8),               intent(in)  :: esp(system%no, system%nk, system%nspin)
    type(s_scalar),        intent(in)  :: rho
    integer,               intent(in)  :: ng
    real(8),               intent(in)  :: gvec(3,ng)
    real(8),               intent(in)  :: gg(ng)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nb_lo, nb_hi
    real(8),               intent(out) :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out) :: zfac   (nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: skip_frac
    logical,    optional,  intent(in)  :: do_remainder
    ! band caps: nb_sigma = Sigma_c intermediate sum, nb_eps = chi0/eps unocc sum.
    integer,    optional,  intent(in)  :: nb_sigma, nb_eps

    complex(8), allocatable :: epsinv(:,:), mblk(:,:,:), msig(:,:), wgpp(:,:), rhoft(:)
    real(8),    allocatable :: vcoul(:), wt(:,:), gvecl(:,:), ggl(:), eband(:), focc(:)
    real(8),    allocatable :: e0(:,:)
    logical,    allocatable :: phys(:,:)
    integer,    allocatable :: imap(:), idiff(:,:), qid(:,:)
    real(8),    allocatable :: qrep(:,:)
    complex(8), allocatable :: sc_all(:,:), scp_all(:,:), scm_all(:,:), scg(:,:)
    ! static-remainder (Coulomb-hole band-truncation correction)
    complex(8), allocatable :: rem_all(:,:), rho_nft(:,:,:), wc0(:,:)
    complex(8) :: pgg, rterm, zt
    real(8)    :: rem_re, rem_im
    integer    :: kk
    logical    :: lrem

    integer :: no, nk, ik, iq, ikm, inp, in, ig, jg, kg, nsig, neps
    integer :: ngl, nchan_phys_tot, nchan_unphys_tot, ncnt(2), ncnt_g(2)
    integer :: nqd, iqd, qd_lo, qd_hi, nper
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, ecut_l, gmax2, fourpi, pi
    real(8) :: rho0, wp2, qgi(3), qg2i, qgj(3), proj, om2, denom, wt2
    real(8) :: eta_au, de_au, e0nk, dsig, qtol
    complex(8) :: rhod, s0, sp, sm
    logical   :: q_ok

    no    = system%no
    nsig  = no;  if (present(nb_sigma)) nsig = min(nb_sigma, no)
    neps  = no;  if (present(nb_eps))   neps = min(nb_eps, no)
    nk    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)
    pi     = acos(-1.0d0)
    fourpi = 4.0d0 * pi
    eta_au = eta_ev / ha2ev
    de_au  = dE_ev  / ha2ev
    qtol   = 1.0d-6                       ! |q1-q2| (bohr^-1) below which q's match

    sigc_re(:,:) = 0.0d0
    zfac(:,:)    = 1.0d0

    lrem = .false.
    if (present(do_remainder)) lrem = do_remainder

    ! ---- larger G-set + rho(G'') (q-independent) -- same as calc_sigma_gpp --
    gmax2 = 0.0d0
    do ig = 1, ng
      if (gg(ig) > gmax2) gmax2 = gg(ig)
    end do
    ecut_l = gset_factor * gmax2
    allocate(gvecl(3, ngmax_l), ggl(ngmax_l))
    call build_gvectors(system%primitive_b, ecut_l, ngmax_l, ngl, gvecl, ggl)
    allocate(rhoft(ngl))
    call build_density_ft(system, info, mg, lg, rho, gvecl, ngl, rhoft, local_only=.true.)
    rho0 = dble(rhoft(1))
    wp2  = fourpi * rho0

    allocate(epsinv(ng,ng), vcoul(ng), imap(ng))
    allocate(wgpp(ng,ng), wt(ng,ng), phys(ng,ng), idiff(ng,ng))
    ! intermediate (bra) dim capped at nsig; output (ket) dim capped at the QP
    ! window upper bound nb_hi (output bands used are only [nb_lo,nb_hi]).
    allocate(mblk(ng, nsig, nb_hi), msig(ng, nsig))
    allocate(eband(nsig), focc(nsig), e0(nb_lo:nb_hi, nk))
    if (lrem) allocate(wc0(ng,ng))

    do jg = 1, ng
      do ig = 1, ng
        gtarget(1) = gvec(1,ig) - gvec(1,jg)
        gtarget(2) = gvec(2,ig) - gvec(2,jg)
        gtarget(3) = gvec(3,ig) - gvec(3,jg)
        idiff(ig,jg) = gindex_of(ngl, gvecl, gtarget)
      end do
    end do

    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! ---- static remainder: per-state density rho_nk(G) on the larger G-set ---
    ! Completeness limit of Sum_n' M*_nn'(G) M_nn'(G') (all bands) = rho_nk(G-G').
    ! Replicated orbitals -> every rank holds all (in,ik); cheap, no comm.
    if (lrem) then
      allocate(rho_nft(ngl, nb_lo:nb_hi, nk), rem_all(nb_lo:nb_hi, nk))
      rem_all(:,:) = (0.0d0, 0.0d0)
      do ik = 1, nk
        do in = nb_lo, nb_hi
          call build_state_density_ft(system, info, mg, lg, spsi, ispin, in, ik, &
                                       gvecl, ngl, rho_nft(:,in,ik), local_only=.true.)
        end do
      end do
    end if

    ! ---- (1) group all (ik,iq) pairs by distinct q = vec_k(iq)-vec_k(ik) ----
    ! Exact for a uniform mesh: equal q's give identical qvec to machine
    ! precision, distinct q's differ by >= |b|/N >> qtol.
    allocate(qid(nk,nk), qrep(3, nk*nk))
    nqd = 0
    do ik = 1, nk
      do iq = 1, nk
        qvec(1:3) = system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        qid(ik,iq) = 0
        do iqd = 1, nqd
          if ( abs(qvec(1)-qrep(1,iqd)) + abs(qvec(2)-qrep(2,iqd)) &
             + abs(qvec(3)-qrep(3,iqd)) < qtol ) then
            qid(ik,iq) = iqd
            exit
          end if
        end do
        if (qid(ik,iq) == 0) then
          nqd = nqd + 1
          qid(ik,iq) = nqd
          qrep(1:3,nqd) = qvec(1:3)
        end if
      end do
    end do

    ! ---- (2) block-distribute the distinct q over icomm_k ------------------
    nper  = (nqd + info%isize_k - 1) / info%isize_k
    qd_lo = info%id_k * nper + 1
    qd_hi = min((info%id_k + 1) * nper, nqd)

    nchan_phys_tot   = 0
    nchan_unphys_tot = 0
    allocate(sc_all (nb_lo:nb_hi,nk), scp_all(nb_lo:nb_hi,nk), scm_all(nb_lo:nb_hi,nk))
    sc_all  = (0.0d0,0.0d0)
    scp_all = (0.0d0,0.0d0)
    scm_all = (0.0d0,0.0d0)

    ! ---- (3) per distinct q: eps^{-1}/v/W once, then its (ik,iq) pairs ------
    do iqd = qd_lo, qd_hi
      qvec(1:3)  =  qrep(1:3,iqd)
      mqvec(1:3) = -qvec(1:3)

      call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, iqd, &
                       qvec, ispin, epsinv, ok=q_ok, local_only=.true., nb_eps=neps)
      if (.not. q_ok) cycle
      call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

      ! full static screened-Coulomb correlation W^c_GG'(q,0) = (eps^{-1}-I) v,
      ! for the static remainder (no GPP pole mask; head finite via nsub_head).
      if (lrem) then
        do jg = 1, ng
          do ig = 1, ng
            wc0(ig,jg) = epsinv(ig,jg) * cmplx(vcoul(jg), 0.0d0, 8)
            if (ig == jg) wc0(ig,jg) = wc0(ig,jg) - cmplx(vcoul(jg), 0.0d0, 8)
          end do
        end do
      end if

      do ig = 1, ng
        qgi(1) = qvec(1) + gvec(1,ig)
        qgi(2) = qvec(2) + gvec(2,ig)
        qgi(3) = qvec(3) + gvec(3,ig)
        qg2i   = qgi(1)*qgi(1) + qgi(2)*qgi(2) + qgi(3)*qgi(3)
        do jg = 1, ng
          phys(ig,jg) = .false.
          wt(ig,jg)   = 0.0d0
          wgpp(ig,jg) = (0.0d0, 0.0d0)
          if (qg2i < 1.0d-12) cycle
          kg = idiff(ig,jg)
          if (kg == 0) cycle
          qgj(1) = qvec(1) + gvec(1,jg)
          qgj(2) = qvec(2) + gvec(2,jg)
          qgj(3) = qvec(3) + gvec(3,jg)
          proj   = ( qgi(1)*qgj(1) + qgi(2)*qgj(2) + qgi(3)*qgj(3) ) / qg2i
          rhod = rhoft(kg)
          om2  = wp2 * proj * ( dble(rhod) / rho0 )
          denom = -dble(epsinv(ig,jg))
          if (ig == jg) denom = denom + 1.0d0
          if (abs(denom) < wt2_tiny) cycle
          wt2 = om2 / denom
          if (wt2 <= 0.0d0) cycle
          phys(ig,jg) = .true.
          wt(ig,jg)   = sqrt(wt2)
          wgpp(ig,jg) = cmplx( vcoul(jg) * om2 / (2.0d0 * wt(ig,jg)), 0.0d0, 8 )
        end do
      end do

      do ik = 1, nk
        do iq = 1, nk
          if (qid(ik,iq) /= iqd) cycle

          call find_kpq(system, ik, mqvec, ikm, g0vec)
          if (ikm == 0) cycle

          do jg = 1, ng
            do ig = 1, ng
              if (phys(ig,jg)) then
                nchan_phys_tot = nchan_phys_tot + 1
              else
                nchan_unphys_tot = nchan_unphys_tot + 1
              end if
            end do
          end do

          do ig = 1, ng
            gtarget(1) = gvec(1,ig) - g0vec(1)
            gtarget(2) = gvec(2,ig) - g0vec(2)
            gtarget(3) = gvec(3,ig) - g0vec(3)
            imap(ig) = gindex_of(ng, gvec, gtarget)
          end do

          ! bra (intermediate) = nsig; ket (output) = nb_hi (QP window upper bound)
          call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                          nsig, nb_hi, mblk, local_only=.true.)

          do inp = 1, nsig
            eband(inp) = esp(inp,ikm,ispin)
            focc(inp)  = system%rocc(inp,ikm,ispin)
          end do
          call normalize_occ(nsig, focc)

          do in = nb_lo, nb_hi
            do inp = 1, nsig
              do ig = 1, ng
                jg = imap(ig)
                if (jg == 0) then
                  msig(ig,inp) = (0.0d0, 0.0d0)
                else
                  msig(ig,inp) = mblk(jg, inp, in)
                end if
              end do
            end do
            e0nk = e0(in,ik)
            call sigma_c_at_energy(ng, nsig, wgpp, wt, phys, msig, eband, focc, &
                                   e0nk,         eta_au, s0)
            call sigma_c_at_energy(ng, nsig, wgpp, wt, phys, msig, eband, focc, &
                                   e0nk + de_au, eta_au, sp)
            call sigma_c_at_energy(ng, nsig, wgpp, wt, phys, msig, eband, focc, &
                                   e0nk - de_au, eta_au, sm)
            sc_all (in,ik) = sc_all (in,ik) + s0
            scp_all(in,ik) = scp_all(in,ik) + sp
            scm_all(in,ik) = scm_all(in,ik) + sm

            ! static remainder: Sum_GG' W^c_GG'(q,0) [ rho_nk(G_jg-G_ig) - P(ig,jg) ]
            ! with P(ig,jg) = Sum_n' M_nn'(G_ig) conjg(M_nn'(G_jg)).  Same contraction
            ! order as wgpp; rho term indexed by idiff(jg,ig) (umklapp cancels).
            if (lrem) then
              rem_re = 0.0d0; rem_im = 0.0d0
!$omp parallel do collapse(2) default(shared) &
!$omp   private(jg,ig,kk,rterm,pgg,inp,zt) reduction(+:rem_re,rem_im)
              do jg = 1, ng
                do ig = 1, ng
                  if (wc0(ig,jg) == (0.0d0,0.0d0)) cycle
                  kk = idiff(jg,ig)
                  if (kk == 0) then
                    rterm = (0.0d0, 0.0d0)
                  else
                    rterm = rho_nft(kk,in,ik)
                  end if
                  pgg = (0.0d0, 0.0d0)
                  do inp = 1, nsig
                    pgg = pgg + msig(ig,inp) * conjg(msig(jg,inp))
                  end do
                  zt = wc0(ig,jg) * (rterm - pgg)
                  rem_re = rem_re + dble(zt)
                  rem_im = rem_im + aimag(zt)
                end do
              end do
!$omp end parallel do
              rem_all(in,ik) = rem_all(in,ik) + cmplx(rem_re, rem_im, 8)
            end if
          end do

        end do
      end do
    end do

    ! ---- (4) assemble the per-(band,k) sums across the q distribution -------
    allocate(scg(nb_lo:nb_hi,nk))
    call comm_summation(sc_all,  scg, size(sc_all),  info%icomm_k);  sc_all  = scg
    call comm_summation(scp_all, scg, size(scp_all), info%icomm_k);  scp_all = scg
    call comm_summation(scm_all, scg, size(scm_all), info%icomm_k);  scm_all = scg
    if (lrem) then
      call comm_summation(rem_all, scg, size(rem_all), info%icomm_k);  rem_all = scg
    end if
    deallocate(scg)
    ncnt(1) = nchan_phys_tot
    ncnt(2) = nchan_unphys_tot
    call comm_summation(ncnt, ncnt_g, 2, info%icomm_k)
    nchan_phys_tot   = ncnt_g(1)
    nchan_unphys_tot = ncnt_g(2)

    ! ---- (5) Re Sigma_c and Z for every output k ---------------------------
    do ik = 1, nk
      do in = nb_lo, nb_hi
        sigc_re(in,ik) = rnk * dble(sc_all(in,ik))
        ! static remainder (energy-independent -> does not enter Z): + 1/4 rnk Re(rem)
        if (lrem) sigc_re(in,ik) = sigc_re(in,ik) + 0.25d0 * rnk * dble(rem_all(in,ik))
        dsig = rnk * ( dble(scp_all(in,ik)) - dble(scm_all(in,ik)) ) / (2.0d0 * de_au)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do
    end do

    if (present(skip_frac)) then
      if (nchan_phys_tot + nchan_unphys_tot > 0) then
        skip_frac = dble(nchan_unphys_tot) / dble(nchan_phys_tot + nchan_unphys_tot)
      else
        skip_frac = 0.0d0
      end if
    end if

    deallocate(epsinv, vcoul, imap, wgpp, wt, phys, idiff)
    deallocate(mblk, msig, rhoft, gvecl, ggl, eband, focc, e0)
    deallocate(qid, qrep, sc_all, scp_all, scm_all)
    if (lrem) deallocate(wc0, rho_nft, rem_all)

  end subroutine calc_sigma_gpp_qcache

  ! --------------------------------------------------------------------------
  ! calc_sigma_gpp_sym
  !
  ! Point-group symmetry-reduced version of calc_sigma_gpp_qcache.  Requires the
  ! orbitals to have been made exactly symmetric first (gw_symmetrize_orbitals),
  ! so the dielectric matrix obeys eps^{-1}(q)=rot[eps^{-1}(q_irr)] to machine
  ! precision.  Two reductions, each ~ n_sym:
  !   q-IBZ : eps^{-1} is solved only for the irreducible momentum transfers; for
  !           a star member q = R q_irr it is the G-index rotation of the rep
  !           (eps^{-1}(q)_{a,b} = eps^{-1}(q_irr)_{iperm(a),iperm(b)}).
  !   out-IBZ: Re Sigma_c and Z are formed only for irreducible output k; the rest
  !           are recovered from Sigma(Rk)=Sigma(k_irr) (a band invariant).
  ! The per-pair body is identical to calc_sigma_gpp_qcache.  Validate against the
  ! full-BZ symmetric result (qcache + symmetrised orbitals).
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_gpp_sym(system, info, mg, lg, spsi, esp, rho, gvec, gg, &
                                ng, ispin, nb_lo, nb_hi, sigc_re, zfac, skip_frac)
    use structures,      only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, s_scalar
    use gw_mtxel_sub,    only: calc_mtxel
    use gw_coulomb_sub,  only: build_vcoul, build_gvectors
    use gw_epsilon_sub,  only: find_kpq, calc_epsinv
    use communication,   only: comm_summation
    use gw_symmetry_sub, only: gw_sym_init_ops, build_g_perm, build_ibz_map
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg, lg
    type(s_orbital),       intent(in)  :: spsi
    real(8),               intent(in)  :: esp(system%no, system%nk, system%nspin)
    type(s_scalar),        intent(in)  :: rho
    integer,               intent(in)  :: ng, ispin, nb_lo, nb_hi
    real(8),               intent(in)  :: gvec(3,ng), gg(ng)
    real(8),               intent(out) :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out) :: zfac   (nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: skip_frac

    complex(8), allocatable :: epsinv(:,:), epsi_irr(:,:), mblk(:,:,:), msig(:,:)
    complex(8), allocatable :: wgpp(:,:), rhoft(:), sc_all(:,:), scp_all(:,:), scm_all(:,:), scg(:,:)
    real(8),    allocatable :: vcoul(:), wt(:,:), gvecl(:,:), ggl(:), eband(:), focc(:), e0(:,:)
    logical,    allocatable :: phys(:,:)
    integer,    allocatable :: imap(:), idiff(:,:), qid(:,:), gperm(:,:), iperm(:,:)
    real(8),    allocatable :: qrep(:,:)
    integer,    allocatable :: qrep_i(:), qop(:), krep(:), kop(:), ibzq(:)

    integer :: no, nk, ik, iq, ikm, inp, in, ig, jg, kg, a, b
    integer :: ngl, nchan_phys_tot, nchan_unphys_tot, ncnt(2), ncnt_g(2)
    integer :: nqd, iqd, nsym, isym, nqirr, nkirr, irep, irlo, irhi, nper
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, ecut_l, gmax2, fourpi, pi
    real(8) :: rho0, wp2, qgi(3), qg2i, qgj(3), proj, om2, denom, wt2
    real(8) :: eta_au, de_au, e0nk, dsig, qtol
    complex(8) :: rhod, s0, sp, sm
    logical   :: q_ok

    no    = system%no
    nk    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)
    pi     = acos(-1.0d0)
    fourpi = 4.0d0 * pi
    eta_au = eta_ev / ha2ev
    de_au  = dE_ev  / ha2ev
    qtol   = 1.0d-6

    sigc_re(:,:) = 0.0d0
    zfac(:,:)    = 1.0d0

    ! larger G-set + rho(G'') (q-independent) -- same as the qcache path
    gmax2 = 0.0d0
    do ig = 1, ng
      if (gg(ig) > gmax2) gmax2 = gg(ig)
    end do
    ecut_l = gset_factor * gmax2
    allocate(gvecl(3, ngmax_l), ggl(ngmax_l))
    call build_gvectors(system%primitive_b, ecut_l, ngmax_l, ngl, gvecl, ggl)
    allocate(rhoft(ngl))
    call build_density_ft(system, info, mg, lg, rho, gvecl, ngl, rhoft, local_only=.true.)
    rho0 = dble(rhoft(1))
    wp2  = fourpi * rho0

    allocate(epsinv(ng,ng), epsi_irr(ng,ng), vcoul(ng), imap(ng))
    allocate(wgpp(ng,ng), wt(ng,ng), phys(ng,ng), idiff(ng,ng))
    allocate(mblk(ng, no, no), msig(ng, no))
    allocate(eband(no), focc(no), e0(nb_lo:nb_hi, nk))

    do jg = 1, ng
      do ig = 1, ng
        gtarget(1) = gvec(1,ig) - gvec(1,jg)
        gtarget(2) = gvec(2,ig) - gvec(2,jg)
        gtarget(3) = gvec(3,ig) - gvec(3,jg)
        idiff(ig,jg) = gindex_of(ngl, gvecl, gtarget)
      end do
    end do
    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! ---- group (ik,iq) pairs by distinct q (same as qcache) -----------------
    allocate(qid(nk,nk), qrep(3, nk*nk))
    nqd = 0
    do ik = 1, nk
      do iq = 1, nk
        qvec(1:3) = system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        qid(ik,iq) = 0
        do iqd = 1, nqd
          if ( abs(qvec(1)-qrep(1,iqd)) + abs(qvec(2)-qrep(2,iqd)) &
             + abs(qvec(3)-qrep(3,iqd)) < qtol ) then
            qid(ik,iq) = iqd; exit
          end if
        end do
        if (qid(ik,iq) == 0) then
          nqd = nqd + 1; qid(ik,iq) = nqd; qrep(1:3,nqd) = qvec(1:3)
        end if
      end do
    end do

    ! ---- symmetry: ops, G-rotation perm (+inverse), q-IBZ and output-k IBZ ---
    call gw_sym_init_ops(system, nsym)
    call build_g_perm(system%primitive_b, gvec, ng, gperm)   ! gperm(ng,nsym)
    allocate(iperm(ng,nsym))
    do isym = 1, nsym
      do ig = 1, ng
        iperm(gperm(ig,isym), isym) = ig
      end do
    end do
    allocate(qrep_i(nqd), qop(nqd))
    call build_ibz_map(system%primitive_b, qrep, nqd, qrep_i, qop, nqirr)
    allocate(krep(nk), kop(nk))
    call build_ibz_map(system%primitive_b, system%vec_k, nk, krep, kop, nkirr)

    ! list of irreducible-q representatives, then block-distribute over icomm_k
    allocate(ibzq(nqd))
    nqirr = 0
    do iqd = 1, nqd
      if (qop(iqd) == 0) then; nqirr = nqirr + 1; ibzq(nqirr) = iqd; end if
    end do
    nper = (nqirr + info%isize_k - 1) / info%isize_k
    irlo = info%id_k * nper + 1
    irhi = min((info%id_k + 1) * nper, nqirr)

    nchan_phys_tot = 0; nchan_unphys_tot = 0
    allocate(sc_all(nb_lo:nb_hi,nk), scp_all(nb_lo:nb_hi,nk), scm_all(nb_lo:nb_hi,nk))
    sc_all = (0d0,0d0); scp_all = (0d0,0d0); scm_all = (0d0,0d0)

    ! ---- per irreducible q: eps^{-1} ONCE, rotate for the star ---------------
    do irep = irlo, irhi
      ! eps^{-1} at the representative q
      call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, ibzq(irep), &
                       qrep(1:3,ibzq(irep)), ispin, epsi_irr, ok=q_ok, local_only=.true.)
      if (.not. q_ok) cycle

      do iqd = 1, nqd
        if (qrep_i(iqd) /= ibzq(irep)) cycle      ! only this rep's star

        ! eps^{-1} at q = R q_irr : identity for the rep, G-rotation otherwise
        if (qop(iqd) == 0) then
          epsinv(:,:) = epsi_irr(:,:)
        else
          do b = 1, ng
            do a = 1, ng
              epsinv(a,b) = epsi_irr(iperm(a,qop(iqd)), iperm(b,qop(iqd)))
            end do
          end do
        end if

        qvec(1:3)  =  qrep(1:3,iqd)
        mqvec(1:3) = -qvec(1:3)
        call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

        ! GPP weights {wgpp,wt,phys} (identical to the qcache path)
        do ig = 1, ng
          qgi(1) = qvec(1) + gvec(1,ig)
          qgi(2) = qvec(2) + gvec(2,ig)
          qgi(3) = qvec(3) + gvec(3,ig)
          qg2i   = qgi(1)*qgi(1) + qgi(2)*qgi(2) + qgi(3)*qgi(3)
          do jg = 1, ng
            phys(ig,jg) = .false.; wt(ig,jg) = 0.0d0; wgpp(ig,jg) = (0.0d0,0.0d0)
            if (qg2i < 1.0d-12) cycle
            kg = idiff(ig,jg)
            if (kg == 0) cycle
            qgj(1) = qvec(1) + gvec(1,jg)
            qgj(2) = qvec(2) + gvec(2,jg)
            qgj(3) = qvec(3) + gvec(3,jg)
            proj   = ( qgi(1)*qgj(1) + qgi(2)*qgj(2) + qgi(3)*qgj(3) ) / qg2i
            rhod = rhoft(kg)
            om2  = wp2 * proj * ( dble(rhod) / rho0 )
            denom = -dble(epsinv(ig,jg))
            if (ig == jg) denom = denom + 1.0d0
            if (abs(denom) < wt2_tiny) cycle
            wt2 = om2 / denom
            if (wt2 <= 0.0d0) cycle
            phys(ig,jg) = .true.
            wt(ig,jg)   = sqrt(wt2)
            wgpp(ig,jg) = cmplx( vcoul(jg) * om2 / (2.0d0 * wt(ig,jg)), 0.0d0, 8 )
          end do
        end do

        ! the (ik,iq) pairs of this q, but only irreducible output k
        do ik = 1, nk
          if (kop(ik) /= 0) cycle                 ! output-IBZ: rep output k only
          do iq = 1, nk
            if (qid(ik,iq) /= iqd) cycle

            call find_kpq(system, ik, mqvec, ikm, g0vec)
            if (ikm == 0) cycle

            do jg = 1, ng
              do ig = 1, ng
                if (phys(ig,jg)) then; nchan_phys_tot = nchan_phys_tot + 1
                else; nchan_unphys_tot = nchan_unphys_tot + 1; end if
              end do
            end do

            do ig = 1, ng
              gtarget(1) = gvec(1,ig) - g0vec(1)
              gtarget(2) = gvec(2,ig) - g0vec(2)
              gtarget(3) = gvec(3,ig) - g0vec(3)
              imap(ig) = gindex_of(ng, gvec, gtarget)
            end do

            call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                            no, no, mblk, local_only=.true.)

            do inp = 1, no
              eband(inp) = esp(inp,ikm,ispin)
              focc(inp)  = system%rocc(inp,ikm,ispin)
            end do
            call normalize_occ(no, focc)

            do in = nb_lo, nb_hi
              do inp = 1, no
                do ig = 1, ng
                  jg = imap(ig)
                  if (jg == 0) then; msig(ig,inp) = (0.0d0,0.0d0)
                  else; msig(ig,inp) = mblk(jg, inp, in); end if
                end do
              end do
              e0nk = e0(in,ik)
              call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, e0nk,         eta_au, s0)
              call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, e0nk + de_au, eta_au, sp)
              call sigma_c_at_energy(ng, no, wgpp, wt, phys, msig, eband, focc, e0nk - de_au, eta_au, sm)
              sc_all (in,ik) = sc_all (in,ik) + s0
              scp_all(in,ik) = scp_all(in,ik) + sp
              scm_all(in,ik) = scm_all(in,ik) + sm
            end do
          end do
        end do
      end do
    end do

    ! ---- assemble over the q distribution -----------------------------------
    allocate(scg(nb_lo:nb_hi,nk))
    call comm_summation(sc_all,  scg, size(sc_all),  info%icomm_k);  sc_all  = scg
    call comm_summation(scp_all, scg, size(scp_all), info%icomm_k);  scp_all = scg
    call comm_summation(scm_all, scg, size(scm_all), info%icomm_k);  scm_all = scg
    deallocate(scg)
    ncnt(1) = nchan_phys_tot; ncnt(2) = nchan_unphys_tot
    call comm_summation(ncnt, ncnt_g, 2, info%icomm_k)
    nchan_phys_tot = ncnt_g(1); nchan_unphys_tot = ncnt_g(2)

    ! ---- Re Sigma_c and Z for irreducible output k, then symmetry recovery --
    do ik = 1, nk
      if (kop(ik) /= 0) cycle
      do in = nb_lo, nb_hi
        sigc_re(in,ik) = rnk * dble(sc_all(in,ik))
        dsig = rnk * ( dble(scp_all(in,ik)) - dble(scm_all(in,ik)) ) / (2.0d0 * de_au)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do
    end do
    do ik = 1, nk
      if (kop(ik) == 0) cycle                     ! star member: copy from its rep
      sigc_re(:,ik) = sigc_re(:,krep(ik))
      zfac(:,ik)    = zfac(:,krep(ik))
    end do

    if (present(skip_frac)) then
      if (nchan_phys_tot + nchan_unphys_tot > 0) then
        skip_frac = dble(nchan_unphys_tot) / dble(nchan_phys_tot + nchan_unphys_tot)
      else
        skip_frac = 0.0d0
      end if
    end if

    deallocate(epsinv, epsi_irr, vcoul, imap, wgpp, wt, phys, idiff)
    deallocate(mblk, msig, rhoft, gvecl, ggl, eband, focc, e0)
    deallocate(qid, qrep, sc_all, scp_all, scm_all)
    deallocate(gperm, iperm, qrep_i, qop, krep, kop, ibzq)

  end subroutine calc_sigma_gpp_sym

  ! --------------------------------------------------------------------------
  ! normalize_occ
  !
  ! Map system%rocc (which folds the spin multiplicity: ~2 per occupied state in
  ! the closed-shell nspin=1 case, ~1 for nspin=2) onto a fractional occupation
  ! f in [0,1] for the GPP occupied/empty pole split.  Scale by the maximum
  ! occupation present (the closed-shell weight); fully occupied -> 1, empty ->
  ! 0.  If all occupations are ~0 (no electrons, pathological) leave them.
  ! --------------------------------------------------------------------------
  subroutine normalize_occ(no, focc)
    implicit none
    integer, intent(in)    :: no
    real(8), intent(inout) :: focc(no)
    integer :: i
    real(8) :: fmax
    fmax = 0.0d0
    do i = 1, no
      if (focc(i) > fmax) fmax = focc(i)
    end do
    if (fmax > occ_thr) then
      do i = 1, no
        focc(i) = focc(i) / fmax
      end do
    end if
  end subroutine normalize_occ

end module gw_sigma_gpp_sub
