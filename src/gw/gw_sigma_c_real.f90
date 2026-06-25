!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Real-axis full-frequency correlation self-energy Sigma_c (spec-b1-4).
!
!   Sigma_c(nk;E) = (1/(N_k Omega)) sum_{n'} sum_q sum_{G,G'}
!        M*_{nn',G}(k,q) M_{nn',G'}(k,q)
!        INT_0^inf dw' B_{GG'}(q,w')
!        [ (1 - f_{n',k-q}) / ( E - eps_{n',k-q} - w' + i eta )
!        +       f_{n',k-q}  / ( E - eps_{n',k-q} + w' - i eta ) ]
!
! with the spectral function of the SCREENED-minus-bare interaction
!   W^c_{GG'}(q,w) = ( eps^{-1}_{GG'}(q,w) - delta_{GG'} ) v(q+G'),
!   B_{GG'}(q,w)   = -(1/pi) Im W^c_{GG'}(q,w)       (w >= 0).
!
! This GENERALISES the generalized-plasmon-pole Sigma_c (gw_sigma_gpp.f90):
! GPP replaces B(w') by a single delta at the plasmon frequency wt with residue
! v Omega^2/(2 wt); here B(w') is the real-frequency continuum taken straight
! from the full-frequency eps^{-1}(w) engine (calc_chi0_freq -> calc_w_freq).
! Removing the single-pole approximation is exactly what lets Im Sigma (the
! scattering / impact-ionisation rate, spec-b2) be read off later; for the gap
! it must reproduce the GPP value to ~0.1-0.2 eV (gate (i)), the residual being
! the true single-pole-vs-full-frequency physics difference.
!
! Normalisation reuses calc_sigma_gpp exactly: the BZ q-sum carries 1/(N_k Omega)
! (v(q+G') from build_vcoul has no cell volume), and W^c inherits the F1(1/N_k)
! + F2(factor-2) chi0 prefactor through eps^{-1}(w), the same normalisation that
! puts the velocity-head eps_inf at the experimental ~11.9 -- so no extra factor
! enters here.  Only Re Sigma_c enters the QP energy; i eta keeps the on-shell
! denominators finite and is discarded by taking Re.
!
! M-element bridge identical to calc_sigma_gpp/calc_sigma_x:
!   find_kpq(system, ik, -q, ikm, g0vec)  => k - q = k_{ikm} + G0,
!   imap(ig) = gindex_of( gvec(:,ig) - g0vec ),
!   M_{nn',G(ig)} = mblk( imap(ig), n', n ).
!
! Launched with nproc_k=1 (all k local); calc_mtxel / calc_chi0_freq / calc_w_freq
! return identical global results on every rank, so the accumulation is a
! deterministic local reduction.  The optional output-k split over info%icomm_k
! mirrors calc_sigma_gpp (local_only).
!
! No proper nouns appear in this file.
!
module gw_sigma_c_real_sub
  implicit none
  private

  public :: calc_sigma_c_real
  public :: calc_sigma_c_real_qcache
  public :: calc_sigma_c_real_sym

  real(8), parameter :: ha2ev   = 27.21138505d0  ! eV per Hartree
  real(8), parameter :: dE_ev   = 0.1d0          ! finite-diff step dE for Z (eV)
  integer, parameter :: ngmax_l     = 200000     ! capacity of the larger G-set (remainder)
  real(8), parameter :: gset_factor = 4.0d0      ! larger cutoff = 4 * eps cutoff

contains

  ! --------------------------------------------------------------------------
  ! build_spectral_weight
  !
  ! M-contracted spectral weight  swt(n',iw) = sum_{G,G'} M_{nn',G} bspec_{GG'}(iw)
  ! M*_{nn',G'}  for a fixed output band n (msig).  This is the O(ng^2) cost of
  ! the correlation self-energy; doing it ONCE here (then reusing it for the three
  ! probe energies, and stripping the omega'-grid complex divides out of the G,G'
  ! double loop) is the key optimisation over evaluating the full integral per
  ! (probe,G,G').  OMP over the (disjoint) omega' slices.
  ! --------------------------------------------------------------------------
  subroutine build_spectral_weight(ng, no, nomega, bspec, msig, swt)
    implicit none
    integer,    intent(in)  :: ng, no, nomega
    real(8),    intent(in)  :: bspec(ng, ng, nomega)
    complex(8), intent(in)  :: msig(ng, no)
    complex(8), intent(out) :: swt(no, nomega)
    integer    :: ig, jg, inp, iw
    complex(8) :: mji, acc, bm
!$omp parallel do default(shared) private(iw, inp, jg, ig, mji, acc, bm) schedule(dynamic)
    do iw = 1, nomega
      do inp = 1, no
        acc = (0.0d0, 0.0d0)
        do jg = 1, ng
          mji = conjg(msig(jg,inp))
          if (mji == (0.0d0,0.0d0)) cycle
          bm = (0.0d0, 0.0d0)
          do ig = 1, ng
            bm = bm + msig(ig,inp) * bspec(ig,jg,iw)
          end do
          acc = acc + bm * mji
        end do
        swt(inp,iw) = acc
      end do
    end do
!$omp end parallel do
  end subroutine build_spectral_weight

  ! --------------------------------------------------------------------------
  ! sigma_c_probe
  !
  ! Sigma_c at one probe energy E from the precomputed swt:
  !   Sigma_c(E) = sum_{n'} sum_{iw} domg * swt(n',iw)
  !                [ (1-f)/(E-eps_n'-w'+i eta) + f/(E-eps_n'+w'-i eta) ].
  ! The complex divides (the dominant cost) appear only here -- once per (n',iw),
  ! NOT per (G,G') -- so this is cheap (no ng^2 factor).
  ! --------------------------------------------------------------------------
  subroutine sigma_c_probe(no, nomega, omg, domg, swt, eband, focc, eprobe, eta, ssig)
    implicit none
    integer,    intent(in)  :: no, nomega
    real(8),    intent(in)  :: omg(nomega), domg
    complex(8), intent(in)  :: swt(no, nomega)
    real(8),    intent(in)  :: eband(no), focc(no)
    real(8),    intent(in)  :: eprobe, eta
    complex(8), intent(out) :: ssig
    integer    :: inp, iw
    complex(8) :: zeta, acc
    real(8)    :: f1, fo, de
    zeta = cmplx(0.0d0, eta, 8)
    ssig = (0.0d0, 0.0d0)
    do inp = 1, no
      fo = focc(inp)
      if (fo < 0.0d0) fo = 0.0d0
      if (fo > 1.0d0) fo = 1.0d0
      f1 = 1.0d0 - fo
      de = eprobe - eband(inp)
      acc = (0.0d0, 0.0d0)
      do iw = 1, nomega
        acc = acc + swt(inp,iw) * ( f1 / ( cmplx(de - omg(iw), 0.0d0, 8) + zeta ) &
                                  + fo / ( cmplx(de + omg(iw), 0.0d0, 8) - zeta ) )
      end do
      ssig = ssig + acc * domg
    end do
  end subroutine sigma_c_probe

  ! --------------------------------------------------------------------------
  ! calc_sigma_c_real
  !
  ! Full-frequency real-axis Re Sigma_c(eps^KS) over the band window for every
  ! output k, plus the renormalisation Z from a central finite difference.
  ! Mirrors calc_sigma_gpp's output-k / q double loop and M assembly, but builds
  ! the screened interaction from the full-frequency eps^{-1}(w) engine instead
  ! of the plasmon-pole model.  No density Fourier transform / Omega^2 / wt
  ! machinery is needed here (that is the GPP single-pole construction).
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_c_real(system, info, mg, lg, spsi, esp, gvec, gg, ng, &
                               ispin, nb_lo, nb_hi, nomega, omega_grid, eta_au, &
                               sigc_re, zfac, local_only)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,   only: calc_mtxel
    use gw_epsilon_sub, only: find_kpq, calc_chi0_freq, calc_w_freq
    use communication,  only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)    :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_rgrid),         intent(in)    :: lg
    type(s_orbital),       intent(in)    :: spsi
    real(8),               intent(in)    :: esp(system%no, system%nk, system%nspin)
    integer,               intent(in)    :: ng
    real(8),               intent(in)    :: gvec(3,ng)
    real(8),               intent(in)    :: gg(ng)
    integer,               intent(in)    :: ispin, nb_lo, nb_hi, nomega
    real(8),               intent(in)    :: omega_grid(nomega)
    real(8),               intent(in)    :: eta_au
    real(8),               intent(out)   :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out)   :: zfac   (nb_lo:nb_hi, system%nk)
    logical,    optional,  intent(in)    :: local_only

    complex(8), allocatable :: mblk(:,:,:)            ! < u_{.,ikm}|e^{-iG.r}|u_{.,ik} >
    complex(8), allocatable :: msig(:,:)              ! M_{nn',G(ig)} for fixed n
    complex(8), allocatable :: chi0_w(:,:,:)          ! chi0_{GG'}(q,w)
    complex(8), allocatable :: epsinv_w(:,:,:)        ! eps^{-1}_{GG'}(q,w)
    real(8),    allocatable :: bspec(:,:,:)           ! -Im W^c / pi
    complex(8), allocatable :: swt(:,:)               ! M-contracted spectral weight (no,nomega)
    real(8),    allocatable :: vcoul(:)
    integer,    allocatable :: imap(:)
    real(8),    allocatable :: eband(:), focc(:)
    real(8),    allocatable :: e0(:,:)
    real(8),    allocatable :: sigc_g(:,:), zfac_g(:,:)

    integer :: no, nk, nq, ik, iq, ikm, inp, in, ig, jg, iw, ik_lo, ik_hi
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, pi, domg, de_au, e0nk, dsig
    complex(8) :: wc, s0, sp, sm
    logical    :: q_ok, ll

    ll = .false.
    if (present(local_only)) ll = local_only

    no    = system%no
    nk    = system%nk
    nq    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)
    pi    = acos(-1.0d0)
    de_au = dE_ev / ha2ev
    ! uniform omega' grid spacing (calc_chi0_freq builds a uniform 0..omega_max)
    if (nomega > 1) then
      domg = omega_grid(2) - omega_grid(1)
    else
      domg = 1.0d0
    end if

    ! zfac accumulates the q-summed derivative dReSigma_c/dE (seed 0); converted
    ! to Z = 1/(1 - dReSigma_c/dE) after the loop (and the icomm_k reduction).
    sigc_re(:,:) = 0.0d0
    zfac(:,:)    = 0.0d0

    allocate(mblk(ng, no, no), msig(ng, no), imap(ng))
    allocate(chi0_w(ng,ng,nomega), epsinv_w(ng,ng,nomega), vcoul(ng))
    allocate(bspec(ng,ng,nomega), swt(no,nomega))
    allocate(eband(no), focc(no))
    allocate(e0(nb_lo:nb_hi, nk))

    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! output-k range: split over icomm_k in the parallel mode, else all k.
    if (ll) then
      ik_lo = info%ik_s
      ik_hi = info%ik_e
    else
      ik_lo = 1
      ik_hi = nk
    end if

    do ik = ik_lo, ik_hi
      do iq = 1, nq
        qvec(1:3)  = system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        mqvec(1:3) = -qvec(1:3)

        ! k - q = k_{ikm} + G0, and the M G-bridge (same convention as GPP)
        call find_kpq(system, ik, mqvec, ikm, g0vec)
        if (ikm == 0) cycle
        do ig = 1, ng
          gtarget(1) = gvec(1,ig) - g0vec(1)
          gtarget(2) = gvec(2,ig) - g0vec(2)
          gtarget(3) = gvec(3,ig) - g0vec(3)
          imap(ig) = gindex_local(ng, gvec, gtarget)
        end do

        ! M block for this (ik, ikm)
        call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                        no, no, mblk, local_only=ll)

        ! full-frequency screened interaction at this q:
        ! chi0(w) -> eps^{-1}(w) -> W^c = (eps^{-1}-I) v -> spectral B = -Im W^c/pi
        call calc_chi0_freq(system, info, mg, lg, spsi, esp, gvec, ng, iq, qvec, &
                            ispin, nomega, omega_grid, eta_au, chi0_w, q_ok, &
                            local_only=ll)
        if (.not. q_ok) cycle
        call calc_w_freq(system, gvec, gg, ng, qvec, nomega, chi0_w, epsinv_w, &
                         vcoul, q_ok)
        if (.not. q_ok) cycle

        do iw = 1, nomega
          do jg = 1, ng
            do ig = 1, ng
              wc = epsinv_w(ig,jg,iw) * vcoul(jg)
              if (ig == jg) wc = wc - vcoul(jg)        ! W^c = (eps^{-1}-I) v
              bspec(ig,jg,iw) = -aimag(wc) / pi
            end do
          end do
        end do

        ! per output band n: assemble M_{nn',G(ig)}, accumulate Re Sigma_c at the
        ! three probes (E0, E0+-dE) for the Z finite difference.
        do inp = 1, no
          eband(inp) = esp(inp,ikm,ispin)
          focc(inp)  = system%rocc(inp,ikm,ispin)
        end do
        call normalize_occ(no, focc)

        do in = nb_lo, nb_hi
          do inp = 1, no
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
          ! O(ng^2) M-contraction ONCE (OMP), then the three probes are cheap.
          call build_spectral_weight(ng, no, nomega, bspec, msig, swt)
          call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, &
                             e0nk,         eta_au, s0)
          call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, &
                             e0nk + de_au, eta_au, sp)
          call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, &
                             e0nk - de_au, eta_au, sm)
          sigc_re(in,ik) = sigc_re(in,ik) + rnk * dble(s0)
          ! d ReSigma_c/dE accumulates the q-sum derivative; Z formed after the loop
          zfac(in,ik) = zfac(in,ik) + rnk * dble(sp - sm) / (2.0d0 * de_au)
        end do
      end do
    end do

    ! reduce the q-distributed partial sums (parallel mode) then form Z.
    if (ll) then
      allocate(sigc_g(nb_lo:nb_hi,nk), zfac_g(nb_lo:nb_hi,nk))
      call comm_summation(sigc_re, sigc_g, size(sigc_re), info%icomm_k); sigc_re = sigc_g
      call comm_summation(zfac,    zfac_g, size(zfac),    info%icomm_k); zfac    = zfac_g
      deallocate(sigc_g, zfac_g)
    end if
    ! finite-difference Z = 1 / (1 - dReSigma_c/dE).
    do ik = 1, nk
      do in = nb_lo, nb_hi
        dsig = zfac(in,ik)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do
    end do

    deallocate(mblk, msig, imap, chi0_w, epsinv_w, vcoul, bspec, swt, eband, focc, e0)

  end subroutine calc_sigma_c_real

  ! --------------------------------------------------------------------------
  ! calc_sigma_c_real_qcache
  !
  ! Production variant: the screened interaction W(q,w) depends ONLY on q, so
  ! group the (ik,iq) pairs by distinct q, block-distribute the distinct q over
  ! info%icomm_k, and build chi0(w) -> W(w) -> spectral B ONCE per distinct q
  ! instead of once per (ik,iq).  This removes the O(nk*nq) redundant full-
  ! frequency dielectric inversions (the dominant cost) and is node-scalable.
  ! Mirrors calc_sigma_gpp_qcache; the per-pair body is the same build_spectral_
  ! weight + sigma_c_probe as the base path, so the result is identical.
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_c_real_qcache(system, info, mg, lg, spsi, esp, gvec, gg, ng, &
                                      ispin, nb_lo, nb_hi, nomega, omega_grid, eta_au, &
                                      sigc_re, zfac, nw_scan, w_scan, k_scan, sigc_scan, &
                                      nb_sigma, nb_eps, do_remainder)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,   only: calc_mtxel
    use gw_epsilon_sub, only: find_kpq, calc_chi0_freq, calc_w_freq
    use gw_coulomb_sub, only: build_gvectors, gw_loop_progress
    use gw_mtxel_sub,     only: build_state_density_ft
    use communication,  only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)    :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg, lg
    type(s_orbital),       intent(in)    :: spsi
    real(8),               intent(in)    :: esp(system%no, system%nk, system%nspin)
    integer,               intent(in)    :: ng
    real(8),               intent(in)    :: gvec(3,ng), gg(ng)
    integer,               intent(in)    :: ispin, nb_lo, nb_hi, nomega
    real(8),               intent(in)    :: omega_grid(nomega), eta_au
    real(8),               intent(out)   :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out)   :: zfac   (nb_lo:nb_hi, system%nk)
    ! Optional spectral scan (yn_out_gw_spectral): complex Sigma_c(in, w) at output
    ! k = k_scan over the probe grid w_scan (Im Sigma_c = scattering rate, the
    ! basis of A(k,w) -- spec-b2 / sp3).  All four present together.
    integer,    optional,  intent(in)    :: nw_scan, k_scan
    real(8),    optional,  intent(in)    :: w_scan(:)
    complex(8), optional,  intent(out)   :: sigc_scan(nb_lo:nb_hi, *)
    ! band caps + static remainder (CH band-truncation correction)
    integer,    optional,  intent(in)    :: nb_sigma, nb_eps
    logical,    optional,  intent(in)    :: do_remainder

    complex(8), allocatable :: mblk(:,:,:), msig(:,:), chi0_w(:,:,:), epsinv_w(:,:,:), swt(:,:)
    real(8),    allocatable :: bspec(:,:,:), vcoul(:)
    integer,    allocatable :: imap(:), qid(:,:)
    real(8),    allocatable :: qrep(:,:), eband(:), focc(:), e0(:,:)
    complex(8), allocatable :: sc_all(:,:), scp_all(:,:), scm_all(:,:), scg(:,:)
    complex(8), allocatable :: scan_all(:,:), scan_g(:,:)
    ! static remainder
    real(8),    allocatable :: gvecl(:,:), ggl(:)
    integer,    allocatable :: idiff(:,:)
    complex(8), allocatable :: rho_nft(:,:,:), wc0(:,:), rem_all(:,:), rem_g(:,:)
    complex(8) :: pgg, rterm, zt
    real(8)    :: rem_re, rem_im, ecut_l, gmax2
    integer    :: ngl, kk, iw0
    logical    :: lrem

    integer :: no, nk, ik, iq, iqd, ikm, inp, in, ig, jg, iw, nqd, nper, qd_lo, qd_hi
    integer :: iws, nws, nsig, neps
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, pi, domg, de_au, e0nk, dsig
    real(8), parameter :: qtol = 1.0d-6
    logical :: do_scan
    real(8) :: tpg_s, tpg_l
    complex(8) :: ssc
    complex(8) :: wc, s0, sp, sm
    logical    :: q_ok

    no    = system%no
    nsig  = no;  if (present(nb_sigma)) nsig = min(nb_sigma, no)
    neps  = no;  if (present(nb_eps))   neps = min(nb_eps, no)
    lrem  = .false.;  if (present(do_remainder)) lrem = do_remainder
    nk    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)
    pi    = acos(-1.0d0)
    de_au = dE_ev / ha2ev
    if (nomega > 1) then
      domg = omega_grid(2) - omega_grid(1)
    else
      domg = 1.0d0
    end if
    ! omega'=0 index for the static remainder W^c(q,0): find the grid point
    ! closest to 0 (don't assume the grid starts at 0).
    iw0 = 1
    do iw = 1, nomega
      if (abs(omega_grid(iw)) < abs(omega_grid(iw0))) iw0 = iw
    end do

    ! intermediate (bra) dim capped at nsig; output (ket) dim capped at nb_hi
    allocate(mblk(ng,nsig,nb_hi), msig(ng,nsig), imap(ng))
    allocate(chi0_w(ng,ng,nomega), epsinv_w(ng,ng,nomega), vcoul(ng))
    allocate(bspec(ng,ng,nomega), swt(nsig,nomega))
    allocate(eband(nsig), focc(nsig), e0(nb_lo:nb_hi,nk))
    allocate(sc_all(nb_lo:nb_hi,nk), scp_all(nb_lo:nb_hi,nk), scm_all(nb_lo:nb_hi,nk))
    sc_all = (0d0,0d0); scp_all = (0d0,0d0); scm_all = (0d0,0d0)

    do_scan = present(sigc_scan) .and. present(w_scan) .and. present(k_scan) .and. present(nw_scan)
    nws = 0
    if (do_scan) then
      nws = nw_scan
      allocate(scan_all(nb_lo:nb_hi, nws)); scan_all = (0d0,0d0)
    end if

    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! ---- static remainder setup: larger G-set, idiff(G-G'), per-state rho_nk ----
    if (lrem) then
      gmax2 = 0.0d0
      do ig = 1, ng
        if (gg(ig) > gmax2) gmax2 = gg(ig)
      end do
      ecut_l = gset_factor * gmax2
      allocate(gvecl(3, ngmax_l), ggl(ngmax_l))
      call build_gvectors(system%primitive_b, ecut_l, ngmax_l, ngl, gvecl, ggl)
      allocate(idiff(ng,ng), wc0(ng,ng))
      do jg = 1, ng
        do ig = 1, ng
          gtarget(1) = gvec(1,ig) - gvec(1,jg)
          gtarget(2) = gvec(2,ig) - gvec(2,jg)
          gtarget(3) = gvec(3,ig) - gvec(3,jg)
          idiff(ig,jg) = gindex_local(ngl, gvecl, gtarget)
        end do
      end do
      allocate(rho_nft(ngl, nb_lo:nb_hi, nk), rem_all(nb_lo:nb_hi, nk))
      rem_all(:,:) = (0.0d0, 0.0d0)
      do ik = 1, nk
        do in = nb_lo, nb_hi
          call build_state_density_ft(system, info, mg, lg, spsi, ispin, in, ik, &
                                       gvecl, ngl, rho_nft(:,in,ik), local_only=.true.)
        end do
      end do
    end if

    ! (1) group (ik,iq) by distinct q = vec_k(iq) - vec_k(ik)
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

    ! (2) block-distribute the distinct q over icomm_k
    nper  = (nqd + info%isize_k - 1) / info%isize_k
    qd_lo = info%id_k * nper + 1
    qd_hi = min((info%id_k + 1) * nper, nqd)

    ! (3) per distinct q: chi0(w) -> W(w) -> spectral B ONCE, then its (ik,iq) pairs
    tpg_s = -1.0d0;  tpg_l = -1.0d0
    call gw_loop_progress('Sigma_c(w) q-loop ', 0, qd_hi-qd_lo+1, tpg_s, tpg_l)
    do iqd = qd_lo, qd_hi
      qvec(1:3)  =  qrep(1:3,iqd)
      mqvec(1:3) = -qvec(1:3)
      call calc_chi0_freq(system, info, mg, lg, spsi, esp, gvec, ng, iqd, qvec, &
                          ispin, nomega, omega_grid, eta_au, chi0_w, q_ok, &
                          local_only=.true., nb_eps=neps)
      if (.not. q_ok) cycle
      call calc_w_freq(system, gvec, gg, ng, qvec, nomega, chi0_w, epsinv_w, vcoul, q_ok)
      if (.not. q_ok) cycle
      do iw = 1, nomega
        do jg = 1, ng
          do ig = 1, ng
            wc = epsinv_w(ig,jg,iw) * vcoul(jg)
            if (ig == jg) wc = wc - vcoul(jg)
            bspec(ig,jg,iw) = -aimag(wc) / pi
            ! static W^c(q,0) = (epsinv_w(0)-I) v  (same convention as GPP)
            if (lrem .and. iw == iw0) wc0(ig,jg) = wc
          end do
        end do
      end do

      do ik = 1, nk
        do iq = 1, nk
          if (qid(ik,iq) /= iqd) cycle
          call find_kpq(system, ik, mqvec, ikm, g0vec)
          if (ikm == 0) cycle
          do ig = 1, ng
            gtarget(1) = gvec(1,ig) - g0vec(1)
            gtarget(2) = gvec(2,ig) - g0vec(2)
            gtarget(3) = gvec(3,ig) - g0vec(3)
            imap(ig) = gindex_local(ng, gvec, gtarget)
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
            call build_spectral_weight(ng, nsig, nomega, bspec, msig, swt)
            e0nk = e0(in,ik)
            call sigma_c_probe(nsig, nomega, omega_grid, domg, swt, eband, focc, e0nk,         eta_au, s0)
            call sigma_c_probe(nsig, nomega, omega_grid, domg, swt, eband, focc, e0nk + de_au, eta_au, sp)
            call sigma_c_probe(nsig, nomega, omega_grid, domg, swt, eband, focc, e0nk - de_au, eta_au, sm)
            sc_all (in,ik) = sc_all (in,ik) + s0
            scp_all(in,ik) = scp_all(in,ik) + sp
            scm_all(in,ik) = scm_all(in,ik) + sm
            ! static remainder: Sum_GG' W^c(q,0)[ rho_nk(G_jg-G_ig) - Sum_n'^nsig M M* ]
            ! (energy-independent -> does not enter Z); contraction order as bspec.
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
            ! spectral scan: complex Sigma_c(in, k_scan; w_scan) -- reuses swt.
            if (do_scan .and. ik == k_scan) then
              do iws = 1, nws
                call sigma_c_probe(nsig, nomega, omega_grid, domg, swt, eband, focc, &
                                   w_scan(iws), eta_au, ssc)
                scan_all(in,iws) = scan_all(in,iws) + ssc
              end do
            end if
          end do
        end do
      end do
      call gw_loop_progress('Sigma_c(w) q-loop ', iqd-qd_lo+1, qd_hi-qd_lo+1, tpg_s, tpg_l)
    end do

    ! (4) assemble the per-(band,k) sums across the q distribution
    allocate(scg(nb_lo:nb_hi,nk))
    call comm_summation(sc_all,  scg, size(sc_all),  info%icomm_k); sc_all  = scg
    call comm_summation(scp_all, scg, size(scp_all), info%icomm_k); scp_all = scg
    call comm_summation(scm_all, scg, size(scm_all), info%icomm_k); scm_all = scg
    if (lrem) then
      allocate(rem_g(nb_lo:nb_hi,nk))
      call comm_summation(rem_all, rem_g, size(rem_all), info%icomm_k); rem_all = rem_g
      deallocate(rem_g)
    end if
    deallocate(scg)
    if (do_scan) then
      allocate(scan_g(nb_lo:nb_hi,nws))
      call comm_summation(scan_all, scan_g, size(scan_all), info%icomm_k)
      scan_all = scan_g; deallocate(scan_g)
      do iws = 1, nws
        do in = nb_lo, nb_hi
          sigc_scan(in,iws) = rnk * scan_all(in,iws)   ! complex Sigma_c(in,k_scan;w_scan)
        end do
      end do
      deallocate(scan_all)
    end if

    ! (5) Re Sigma_c and Z for every output k
    do ik = 1, nk
      do in = nb_lo, nb_hi
        sigc_re(in,ik) = rnk * dble(sc_all(in,ik))
        ! static remainder (energy-independent -> does not enter Z): + 1/4 rnk Re(rem)
        if (lrem) sigc_re(in,ik) = sigc_re(in,ik) + 0.25d0 * rnk * dble(rem_all(in,ik))
        dsig = rnk * ( dble(scp_all(in,ik)) - dble(scm_all(in,ik)) ) / (2.0d0 * de_au)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do
    end do

    deallocate(mblk, msig, imap, chi0_w, epsinv_w, vcoul, bspec, swt, eband, focc, e0)
    deallocate(qid, qrep, sc_all, scp_all, scm_all)
    if (lrem) deallocate(gvecl, ggl, idiff, wc0, rho_nft, rem_all)

  end subroutine calc_sigma_c_real_qcache

  ! --------------------------------------------------------------------------
  ! calc_sigma_c_real_sym
  !
  ! Point-group symmetry-reduced variant: compute the full-frequency eps^{-1}(q,w)
  ! ONLY for the irreducible q, reconstruct each star member by the (omega-
  ! independent) G-rotation eps^{-1}(Rq)_{GG'} = eps^{-1}(q)_{R^-1 G, R^-1 G'}, and
  ! accumulate Sigma_c only for the irreducible output k, recovering Sigma_c(Rk) =
  ! Sigma_c(k_irr) at the end.  Mirrors calc_sigma_gpp_sym; the dominant chi0(w)->
  ! W(w) inversions drop by the point-group order (~24) on top of the q-cache.
  ! Needs the symmetrised orbitals (gw_symmetrize_orbitals in gw_main) + sym.dat.
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_c_real_sym(system, info, mg, lg, spsi, esp, gvec, gg, ng, &
                                   ispin, nb_lo, nb_hi, nomega, omega_grid, eta_au, &
                                   sigc_re, zfac)
    use structures,      only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,    only: calc_mtxel
    use gw_coulomb_sub,  only: build_vcoul
    use gw_epsilon_sub,  only: find_kpq, calc_chi0_freq, calc_w_freq
    use gw_symmetry_sub, only: gw_sym_init_ops, build_g_perm, build_ibz_map
    use communication,   only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)    :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg, lg
    type(s_orbital),       intent(in)    :: spsi
    real(8),               intent(in)    :: esp(system%no, system%nk, system%nspin)
    integer,               intent(in)    :: ng
    real(8),               intent(in)    :: gvec(3,ng), gg(ng)
    integer,               intent(in)    :: ispin, nb_lo, nb_hi, nomega
    real(8),               intent(in)    :: omega_grid(nomega), eta_au
    real(8),               intent(out)   :: sigc_re(nb_lo:nb_hi, system%nk)
    real(8),               intent(out)   :: zfac   (nb_lo:nb_hi, system%nk)

    complex(8), allocatable :: mblk(:,:,:), msig(:,:), chi0_w(:,:,:), swt(:,:)
    complex(8), allocatable :: epsi_irr(:,:,:), epsinv(:,:,:)
    real(8),    allocatable :: bspec(:,:,:), vcoul(:)
    integer,    allocatable :: imap(:), qid(:,:), gperm(:,:), iperm(:,:)
    integer,    allocatable :: qrep_i(:), qop(:), krep(:), kop(:), ibzq(:)
    real(8),    allocatable :: qrep(:,:), eband(:), focc(:), e0(:,:)
    complex(8), allocatable :: sc_all(:,:), scp_all(:,:), scm_all(:,:), scg(:,:)

    integer :: no, nk, ik, iq, iqd, ikm, inp, in, ig, jg, iw, nqd, nper, irlo, irhi
    integer :: nsym, isym, irep, nqirr, nkirr, a, b, op
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, pi, domg, de_au, e0nk, dsig
    real(8), parameter :: qtol = 1.0d-6
    integer, parameter :: nsub_head = 10
    complex(8) :: wc, s0, sp, sm
    logical    :: q_ok

    no    = system%no
    nk    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)
    pi    = acos(-1.0d0)
    de_au = dE_ev / ha2ev
    if (nomega > 1) then
      domg = omega_grid(2) - omega_grid(1)
    else
      domg = 1.0d0
    end if

    allocate(mblk(ng,no,no), msig(ng,no), imap(ng))
    allocate(chi0_w(ng,ng,nomega), epsi_irr(ng,ng,nomega), epsinv(ng,ng,nomega), vcoul(ng))
    allocate(bspec(ng,ng,nomega), swt(no,nomega))
    allocate(eband(no), focc(no), e0(nb_lo:nb_hi,nk))
    allocate(sc_all(nb_lo:nb_hi,nk), scp_all(nb_lo:nb_hi,nk), scm_all(nb_lo:nb_hi,nk))
    sc_all = (0d0,0d0); scp_all = (0d0,0d0); scm_all = (0d0,0d0)

    do ik = 1, nk
      do in = nb_lo, nb_hi
        e0(in,ik) = esp(in,ik,ispin)
      end do
    end do

    ! group (ik,iq) by distinct q
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

    ! symmetry ops, G-rotation perm (+inverse), q-IBZ and output-k IBZ
    call gw_sym_init_ops(system, nsym)
    allocate(gperm(ng,nsym), iperm(ng,nsym))
    call build_g_perm(system%primitive_b, gvec, ng, gperm)
    do isym = 1, nsym
      do ig = 1, ng
        iperm(gperm(ig,isym), isym) = ig
      end do
    end do
    allocate(qrep_i(nqd), qop(nqd), krep(nk), kop(nk))
    call build_ibz_map(system%primitive_b, qrep, nqd, qrep_i, qop, nqirr)
    call build_ibz_map(system%primitive_b, system%vec_k, nk, krep, kop, nkirr)
    allocate(ibzq(nqd))
    nqirr = 0
    do iqd = 1, nqd
      if (qop(iqd) == 0) then; nqirr = nqirr + 1; ibzq(nqirr) = iqd; end if
    end do
    nper = (nqirr + info%isize_k - 1) / info%isize_k
    irlo = info%id_k * nper + 1
    irhi = min((info%id_k + 1) * nper, nqirr)

    ! per irreducible q: chi0(w) -> eps^{-1}(w) ONCE; rotate for the star
    do irep = irlo, irhi
      qvec(1:3) = qrep(1:3, ibzq(irep))
      call calc_chi0_freq(system, info, mg, lg, spsi, esp, gvec, ng, ibzq(irep), qvec, &
                          ispin, nomega, omega_grid, eta_au, chi0_w, q_ok, local_only=.true.)
      if (.not. q_ok) cycle
      call calc_w_freq(system, gvec, gg, ng, qvec, nomega, chi0_w, epsi_irr, vcoul, q_ok)
      if (.not. q_ok) cycle

      do iqd = 1, nqd
        if (qrep_i(iqd) /= ibzq(irep)) cycle      ! only this rep's star

        op = qop(iqd)
        if (op == 0) then
          epsinv(:,:,:) = epsi_irr(:,:,:)
        else
          do iw = 1, nomega
            do b = 1, ng
              do a = 1, ng
                epsinv(a,b,iw) = epsi_irr(iperm(a,op), iperm(b,op), iw)
              end do
            end do
          end do
        end if

        qvec(1:3)  =  qrep(1:3,iqd)
        mqvec(1:3) = -qvec(1:3)
        call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)
        do iw = 1, nomega
          do jg = 1, ng
            do ig = 1, ng
              wc = epsinv(ig,jg,iw) * vcoul(jg)
              if (ig == jg) wc = wc - vcoul(jg)
              bspec(ig,jg,iw) = -aimag(wc) / pi
            end do
          end do
        end do

        ! pairs of this q with IRREDUCIBLE output k only
        do ik = 1, nk
          if (kop(ik) /= 0) cycle
          do iq = 1, nk
            if (qid(ik,iq) /= iqd) cycle
            call find_kpq(system, ik, mqvec, ikm, g0vec)
            if (ikm == 0) cycle
            do ig = 1, ng
              gtarget(1) = gvec(1,ig) - g0vec(1)
              gtarget(2) = gvec(2,ig) - g0vec(2)
              gtarget(3) = gvec(3,ig) - g0vec(3)
              imap(ig) = gindex_local(ng, gvec, gtarget)
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
                  if (jg == 0) then
                    msig(ig,inp) = (0.0d0, 0.0d0)
                  else
                    msig(ig,inp) = mblk(jg, inp, in)
                  end if
                end do
              end do
              call build_spectral_weight(ng, no, nomega, bspec, msig, swt)
              e0nk = e0(in,ik)
              call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, e0nk,         eta_au, s0)
              call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, e0nk + de_au, eta_au, sp)
              call sigma_c_probe(no, nomega, omega_grid, domg, swt, eband, focc, e0nk - de_au, eta_au, sm)
              sc_all (in,ik) = sc_all (in,ik) + s0
              scp_all(in,ik) = scp_all(in,ik) + sp
              scm_all(in,ik) = scm_all(in,ik) + sm
            end do
          end do
        end do
      end do
    end do

    ! assemble over the q distribution
    allocate(scg(nb_lo:nb_hi,nk))
    call comm_summation(sc_all,  scg, size(sc_all),  info%icomm_k); sc_all  = scg
    call comm_summation(scp_all, scg, size(scp_all), info%icomm_k); scp_all = scg
    call comm_summation(scm_all, scg, size(scm_all), info%icomm_k); scm_all = scg
    deallocate(scg)

    ! Re Sigma_c and Z for irreducible output k, then symmetry recovery
    do ik = 1, nk
      if (kop(ik) /= 0) cycle
      do in = nb_lo, nb_hi
        sigc_re(in,ik) = rnk * dble(sc_all(in,ik))
        dsig = rnk * ( dble(scp_all(in,ik)) - dble(scm_all(in,ik)) ) / (2.0d0 * de_au)
        zfac(in,ik) = 1.0d0 / (1.0d0 - dsig)
      end do
    end do
    do ik = 1, nk
      if (kop(ik) == 0) cycle
      sigc_re(:,ik) = sigc_re(:,krep(ik))
      zfac(:,ik)    = zfac(:,krep(ik))
    end do

    deallocate(mblk, msig, imap, chi0_w, epsi_irr, epsinv, vcoul, bspec, swt, eband, focc, e0)
    deallocate(qid, qrep, sc_all, scp_all, scm_all)
    deallocate(gperm, iperm, qrep_i, qop, krep, kop, ibzq)

  end subroutine calc_sigma_c_real_sym

  ! --------------------------------------------------------------------------
  ! gindex_local : list index jg with gvec(:,jg) == gtarget (tolerance), 0 if
  ! outside the cutoff sphere.  Private copy (same role as the helpers in
  ! gw_epsilon / gw_sigma_gpp); keeps the module self-contained.
  ! --------------------------------------------------------------------------
  integer function gindex_local(ng, gv, gtarget) result(jg)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in) :: gv(3,ng), gtarget(3)
    integer :: i
    real(8), parameter :: tol = 1.0d-6
    jg = 0
    do i = 1, ng
      if ( abs(gv(1,i)-gtarget(1)) + abs(gv(2,i)-gtarget(2)) &
         + abs(gv(3,i)-gtarget(3)) < tol ) then
        jg = i
        return
      end if
    end do
  end function gindex_local

  ! occupation -> [0,1] per spin (rocc is 2 per doubly-filled state for nspin=1).
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
    if (fmax > 1.0d0) then
      do i = 1, no
        focc(i) = focc(i) / fmax
      end do
    end if
  end subroutine normalize_occ

end module gw_sigma_c_real_sub
