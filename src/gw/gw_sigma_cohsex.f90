!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Static Coulomb-hole + screened-exchange (COHSEX) self-energy for the GW
! module -- the first SCREENED rung (ladder rung 2).
!
! This rung promotes the bare-exchange Sigma_x of rung 1 to a screened
! interaction W = eps^{-1} v and adds the static Coulomb-hole term.  It consumes
! the inverse dielectric matrix eps^{-1}_{GG'}(q,0) of the static-RPA rung
! (calc_epsinv) and reuses the bare Coulomb kernel v(q+G)=4 pi/|q+G|^2
! (build_vcoul, which carries NO 1/Omega -- the cell volume enters the
! self-energy prefactor instead, exactly as in calc_sigma_x).
!
! Physics (atomic units)
! ----------------------
! Static screened interaction:  W_{GG'}(q) = eps^{-1}_{GG'}(q,0) * v(q+G').
!
! Screened exchange (SEX) -- like Sigma_x but bare v -> W:
!
!   <nk|Sigma_SEX|nk> = - (1/(N_k Omega)) sum_q sum_v^occ sum_{G,G'}
!         M*_{nv,G}(k,q) * eps^{-1}_{GG'}(q) * v(q+G') * M_{nv,G'}(k,q),
!
! with  M_{nv,G}(k,q) = < n,k | e^{i(q+G).r} | v,k-q >.  Setting eps^{-1}=delta
! recovers Sigma_x exactly (the double-G sum collapses to sum_G v(q+G)|M_G|^2),
! which fixes both the overall sign (Sigma_SEX<0) and the 1/(N_k Omega) factor.
! Screening reduces |Sigma_SEX| below |Sigma_x|.
!
! Coulomb hole (COH), static, via completeness (no slow empty-band sum):
!
!   <nk|Sigma_COH|nk> = (1/2)(1/(N_k Omega)) sum_q sum_{G,G'}
!         [ eps^{-1}_{GG'}(q) - delta_{GG'} ] * v(q+G') * rho_n(k; G-G'),
!
! where the state density form factor is
!
!   rho_n(k; G'') = < u_{nk} | e^{i G''.r} | u_{nk} >
!                 = (Fourier component of |u_{nk}(r)|^2 at G'').
!
! Derivation of the COH closure.  The frequency integral of the (W - v) part of
! GW, taken at the static (plasmon-pole -> static) limit, leaves a sum over ALL
! intermediate bands n':
!   Sigma_COH = (1/2) sum_{n'} sum_{G,G'} M*_{nn',G} [eps^{-1}-delta]_{GG'} v M_{nn',G'}
!   (times the mesh+cell prefactor).  The completeness relation of the full band
! set at fixed (k,q=0) gives
!   sum_{n'} M*_{nn',G}(k,0) M_{nn',G'}(k,0)
!     = sum_{n'} <u_nk|e^{-iG.r}|u_n'k><u_n'k|e^{iG'.r}|u_nk>
!     = <u_nk| e^{-iG.r} ( sum_{n'} |u_n'k><u_n'k| ) e^{iG'.r} |u_nk>
!     = <u_nk| e^{-i(G-G').r} |u_nk>
!     = rho_n(k; G-G'),
! using sum_{n'}|u_n'k><u_n'k| = 1 (closure over the complete band set).  This
! removes the explicit band sum: COH needs only the diagonal density form factor
! of the single state n, evaluated at the differences G-G'.  (Note the static
! COH is q-independent in its band structure -- the only q-dependence is through
! eps^{-1}(q) and v(q+G'); the form factor is taken at q=0, ik=ikq=k.)
!
! Larger G-set for rho_n
! ----------------------
! G and G' each range over the eps G-list (|G|^2 <= eps_cutoff), so the argument
! G-G' can reach |G-G'|^2 up to ~4 eps_cutoff -- well outside the eps list.  We
! therefore build rho_n on a SEPARATE, LARGER G-set (cutoff ~ 4 eps_cutoff: the
! eps cutoff is max(gg) here, scaled by gset_factor=4) and look G-G' up there.
! If a needed G-G' still falls outside the larger set, that (G,G') pair is
! SKIPPED and COUNTED (the skipped fraction is returned and should be small --
! it lives in the high-|G-G'| corner where v(q+G') is small).
!
! Computing rho_n from calc_mtxel.  calc_mtxel at q=0 (ik=ikq=k) returns
!   mtxL(ig, n, n) = < u_{nk} | e^{-i G_L(ig).r} | u_{nk} > = rho_n(k; -G_L(ig)),
! on the LARGER G-set G_L.  Hence
!   rho_n(k; G'') = mtxL( gindexL(-G''), n, n ),
! and for the COH argument G''=G-G' we look up gindexL(-(G-G')) = gindexL(G'-G).
! The list G_L is closed under negation, so equivalently rho_n(k;G'') and its
! negation partner are both present whenever G'' is inside the larger sphere.
!
! Total and QP
! ------------
!   Sigma_COHSEX = Sigma_SEX + Sigma_COH.  (Do NOT add bare Sigma_x on top: the
!   SEX term already IS the screened exchange.)  Linearised QP, Z=1 (static, no
!   frequency dependence):  eps^QP = eps^KS + (Sigma_COHSEX - <V_xc>).
!
! Sign / magnitude sanity
! -----------------------
! Sigma_SEX < 0 (attractive, |.| < |Sigma_x| because of screening); Sigma_COH<0
! (the Coulomb hole lowers every level).  For a semiconductor the COHSEX gap
! sits ABOVE the KS gap (COHSEX is known to OVER-open relative to one-shot GW).
!
! Prefactor 1/(N_k Omega) -- identical reasoning to calc_sigma_x
! -------------------------------------------------------------
! The BZ integral (1/(2 pi)^3) int_BZ dq -> (1/(N_k Omega)) sum_q, and v(q+G)
! from build_vcoul carries no 1/Omega, so the cell volume Omega = system%det_a
! enters the prefactor explicitly here (this was the fixed Sigma_x bug).  Both
! SEX and the SEX-equivalent part of COH share the same 1/(N_k Omega); COH
! additionally carries the standard 1/2.
!
! Parallel layout
! ---------------
! Launched with nproc_k=1 (all k local); calc_epsinv and calc_mtxel are
! collective (comm_summation over icomm_rko) and return identical global results
! on every rank, so the COHSEX accumulation is a deterministic local reduction.
!
! No proper nouns appear in this file.
!
module gw_sigma_cohsex_sub
  implicit none
  private

  public :: calc_sigma_cohsex

  real(8), parameter :: occ_thr     = 1.0d-6   ! occupied if rocc > occ_thr
  real(8), parameter :: g_match_tol = 1.0d-6   ! |G|-match tolerance (bohr^-1)^2
  integer, parameter :: nsub_head   = 10       ! mini-BZ sub-sampling for v head
  integer, parameter :: ngmax_l     = 200000   ! capacity of the larger G-set
  real(8), parameter :: gset_factor = 4.0d0    ! larger cutoff = 4 * eps cutoff

contains

  ! --------------------------------------------------------------------------
  ! gindex_of
  !
  ! List index jg with gv(:,jg) == gtarget (within tolerance), 0 if gtarget is
  ! outside the list (cutoff sphere).  Local copy (kept private so the module is
  ! self-contained); identical role to the helper in gw_epsilon / gw_sigma_x.
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
  ! calc_sigma_cohsex
  !
  ! Diagonal static COHSEX self-energy expectation <n k|Sigma_COHSEX|n k> =
  ! Sigma_SEX + Sigma_COH for a band window [nb_lo,nb_hi] and all k-points, spin
  ! ispin.  Returns real (a.u.).  Internally calls calc_epsinv per q for
  ! eps^{-1}(q) and calc_mtxel for the M blocks and the rho_n form factors.
  !
  !   sigc(n,k) = - (1/(N_k Omega)) sum_q sum_v^occ sum_{G,G'}
  !                   M*_{nv,G} eps^{-1}_{GG'}(q) v(q+G') M_{nv,G'}
  !             + (1/2)(1/(N_k Omega)) sum_q sum_{G,G'}
  !                   [eps^{-1}_{GG'}(q)-delta] v(q+G') rho_n(k; G-G').
  !
  ! q-set: q = vec_k(:,iq) - vec_k(:,ik), iq=1..nk (the BZ-mesh differences,
  ! including q=0 with its mini-BZ head), matching calc_sigma_x / [gw][t4].
  !
  ! Arguments:
  !   system   [in]  -- vec_k, rocc, det_a (Omega), nk, no, nspin, hvol, lattice.
  !   info     [in]  -- parallel ownership / communicator (to calc_mtxel/epsinv).
  !   mg, lg   [in]  -- real-space grid descriptors.
  !   spsi     [in]  -- cell-periodic orbitals.
  !   esp(no,nk,nspin)[in] -- KS eigenvalues (a.u.); pass energy%esp (to epsinv).
  !   gvec(3,ng)[in] -- Cartesian G vectors (bohr^-1), the eps G-set.
  !   gg(ng)   [in]  -- |G|^2 (bohr^-2), the eps G-set.
  !   ng       [in]  -- number of eps G vectors (eps^{-1} matrix dimension).
  !   ispin    [in]  -- spin channel (1 for the non-magnetic case).
  !   nb_lo,nb_hi[in]-- output band window (states n).
  !   sigc(nb_lo:nb_hi, system%nk) [out] -- Sigma_COHSEX expectation (a.u.).
  !   sex_out, coh_out (nb_lo:nb_hi, nk) [out,opt] -- the SEX and COH parts.
  !   skip_frac [out,opt] -- fraction of (G,G') pairs skipped in COH (G-G'
  !                          outside the larger set), averaged over q.
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_cohsex(system, info, mg, lg, spsi, esp, gvec, gg, ng, &
                               ispin, nb_lo, nb_hi, sigc, sex_out, coh_out, skip_frac)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,   only: calc_mtxel
    use gw_coulomb_sub, only: build_vcoul, build_gvectors
    use gw_epsilon_sub, only: find_kpq, calc_epsinv
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_orbital),       intent(in)  :: spsi
    real(8),               intent(in)  :: esp(system%no, system%nk, system%nspin)
    integer,               intent(in)  :: ng
    real(8),               intent(in)  :: gvec(3,ng)
    real(8),               intent(in)  :: gg(ng)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nb_lo, nb_hi
    real(8),               intent(out) :: sigc(nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: sex_out(nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: coh_out(nb_lo:nb_hi, system%nk)
    real(8),    optional,  intent(out) :: skip_frac

    complex(8), allocatable :: epsinv(:,:)       ! eps^{-1}_{GG'}(q)
    complex(8), allocatable :: mblk(:,:,:)       ! < u_{.,ikm}|e^{-iG.r}|u_{.,ik} >
    complex(8), allocatable :: msex(:,:)         ! M*_{nv,G(ig)} per output band n
    complex(8), allocatable :: wv(:,:)           ! eps^{-1}_{GG'} v(q+G') (W-like)
    complex(8), allocatable :: mtxl(:,:,:)       ! rho form factors on larger set
    real(8),    allocatable :: vcoul(:)
    integer,    allocatable :: imap(:)           ! eps-set G->mtxel column (SEX)
    integer,    allocatable :: idiff(:,:)        ! larger-set index of G(jg)-G(ig)
    real(8),    allocatable :: gvecl(:,:), ggl(:)

    complex(8), allocatable :: sexk(:), cohk(:)  ! per-k accumulators (band window)

    integer :: no, nk, nq, ik, iq, ikm, iv, in, ig, jg, kg
    integer :: ngl, npair_tot, npair_skip
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, ecut_l, gmax2
    complex(8) :: ssex, scoh, zw, mg_ig, mg_jg
    logical   :: q_ok

    no    = system%no
    nk    = system%nk
    nq    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)  ! 1/(N_k*Omega): v(q+G)=4pi/|q+G|^2 carries
                                        ! no 1/Omega, so the cell volume enters here

    sigc(:,:) = 0.0d0
    if (present(sex_out)) sex_out(:,:) = 0.0d0
    if (present(coh_out)) coh_out(:,:) = 0.0d0

    ! ---- larger G-set for the rho_n(k;G-G') form factors -------------------
    ! cutoff ~ 4 * (eps cutoff).  The eps cutoff in bohr^-2 is max(gg(1:ng));
    ! build_gvectors takes the cutoff as |G|^2 (Ry == bohr^-2 here), so pass
    ! gset_factor * max(gg).
    gmax2  = 0.0d0
    do ig = 1, ng
      if (gg(ig) > gmax2) gmax2 = gg(ig)
    end do
    ecut_l = gset_factor * gmax2

    allocate(gvecl(3, ngmax_l), ggl(ngmax_l))
    call build_gvectors(system%primitive_b, ecut_l, ngmax_l, ngl, gvecl, ggl)

    allocate(epsinv(ng,ng), vcoul(ng), imap(ng), wv(ng,ng))
    allocate(mblk(ng, no, no))
    allocate(msex(ng, nb_lo:nb_hi))
    ! mtxl carries calc_mtxel's full 1..nb_hi band block (its output is indexed
    ! from 1); only the diagonal columns in_in are read in the COH sum, but the
    ! allocation/shape must match what calc_mtxel writes.
    allocate(mtxl(ngl, nb_hi, nb_hi))
    allocate(idiff(ng,ng))
    allocate(sexk(nb_lo:nb_hi), cohk(nb_lo:nb_hi))

    ! ---- larger-set lookup of G(jg)-G(ig) (the COH argument G-G', here the
    ! "G" of e^{iG''.r} is G-G' with rows ig, cols jg).  rho_n(k;G(ig)-G(jg))
    ! lives at column gindexL(-(G(ig)-G(jg))) = gindexL(G(jg)-G(ig)).  Precompute
    ! the larger-set column index once (q-independent: it depends only on the
    ! eps G-list, which is fixed).  0 marks a difference outside the larger set.
    npair_skip = 0
    npair_tot  = ng * ng
    do jg = 1, ng
      do ig = 1, ng
        gtarget(1) = gvec(1,jg) - gvec(1,ig)
        gtarget(2) = gvec(2,jg) - gvec(2,ig)
        gtarget(3) = gvec(3,jg) - gvec(3,ig)
        kg = gindex_of(ngl, gvecl, gtarget)
        idiff(ig,jg) = kg
        if (kg == 0) npair_skip = npair_skip + 1
      end do
    end do

    ! ---- output-state loop -------------------------------------------------
    do ik = 1, nk

      sexk(:) = (0.0d0, 0.0d0)
      cohk(:) = (0.0d0, 0.0d0)

      ! rho_n(k;G'') form factors for the COH state-density: diagonal of
      ! calc_mtxel at q=0 (ikq=ik) on the LARGER set.  mtxl(:, n, n) is used.
      call calc_mtxel(system, info, mg, lg, spsi, gvecl, ngl, ik, ik, ispin, &
                      nb_hi, nb_hi, mtxl)

      ! ---- BZ q-sum: q = vec_k(:,iq) - vec_k(:,ik) ------------------------
      do iq = 1, nq
        qvec(1:3)  =  system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        mqvec(1:3) = -qvec(1:3)

        ! eps^{-1}_{GG'}(q) on the eps G-set (collective inside).
        call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, iq, &
                         qvec, ispin, epsinv, ok=q_ok)
        if (.not. q_ok) cycle              ! q not representable on this mesh

        ! v(q+G'): bare Coulomb kernel; q->0 head averaged over the mini-BZ.
        call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

        ! W-like matrix W_{GG'} = eps^{-1}_{GG'} * v(q+G') (v on the COLUMN G').
        ! Shared by SEX (contracted with M on both sides) and COH (its G-G' form
        ! factor weight).  For COH the COULOMB-HOLE bracket is [eps^{-1}-delta].
        do jg = 1, ng
          do ig = 1, ng
            wv(ig,jg) = epsinv(ig,jg) * vcoul(jg)
          end do
        end do

        ! ================= COH part (state-density closure) ================
        ! Sigma_COH += (1/2) sum_{G,G'} [eps^{-1}_{GG'}-delta] v(q+G')
        !                                 rho_n(k; G-G').
        ! [eps^{-1}-delta]*v = wv with the diagonal v(q+G) subtracted.
        do in = nb_lo, nb_hi
          scoh = (0.0d0, 0.0d0)
          do jg = 1, ng
            do ig = 1, ng
              kg = idiff(ig,jg)
              if (kg == 0) cycle           ! G-G' outside larger set: skip+count
              zw = wv(ig,jg)
              if (ig == jg) zw = zw - vcoul(jg)   ! subtract the delta_{GG'} v
              ! rho_n(k; G(ig)-G(jg)) = mtxl(idiff(ig,jg), n, n)
              scoh = scoh + zw * mtxl(kg, in, in)
            end do
          end do
          cohk(in) = cohk(in) + 0.5d0 * scoh
        end do

        ! ================= SEX part (screened exchange) ====================
        ! ket index of k-q (= k + (-q) folded) and umklapp G0.
        call find_kpq(system, ik, mqvec, ikm, g0vec)
        if (ikm == 0) cycle                ! no mesh partner (skip SEX for this q)

        ! Per-(k,q) eps-set index map: imap(ig) -> calc_mtxel column giving
        ! M*_{nv,G(ig)} (output n at ik, occupied v at ikm).  Same +(G-G0)
        ! convention as calc_sigma_x.
        do ig = 1, ng
          gtarget(1) = gvec(1,ig) - g0vec(1)
          gtarget(2) = gvec(2,ig) - g0vec(2)
          gtarget(3) = gvec(3,ig) - g0vec(3)
          imap(ig) = gindex_of(ng, gvec, gtarget)   ! 0 if (G-G0) outside cutoff
        end do

        ! Full M block at (ik,ikm): mblk(ig', bra@ikm, ket@ik) on the eps set.
        call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                        no, no, mblk)

        do iv = 1, no
          if (system%rocc(iv,ikm,ispin) <= occ_thr) cycle   ! v occupied at k-q

          ! Gather M*_{nv,G(ig)} = mblk(imap(ig), v, n) for every output band n.
          ! Entries whose (G-G0) leaves the cutoff are unavailable -> 0.
          do in = nb_lo, nb_hi
            do ig = 1, ng
              jg = imap(ig)
              if (jg == 0) then
                msex(ig,in) = (0.0d0, 0.0d0)
              else
                msex(ig,in) = mblk(jg, iv, in)            ! = M*_{nv,G(ig)}
              end if
            end do
          end do

          ! Sigma_SEX += - sum_{G,G'} M*_G eps^{-1}_{GG'} v(q+G') M_{G'}.
          ! With A_ig = msex(ig,in) = M*_{nv,G(ig)}, M_{nv,G'(jg)} = conjg(A_jg),
          ! and wv(ig,jg) = eps^{-1}_{ig,jg} v(q+G'(jg)) :
          !   sum_{ig,jg} A_ig wv(ig,jg) conjg(A_jg).
          do in = nb_lo, nb_hi
            ssex = (0.0d0, 0.0d0)
            do jg = 1, ng
              mg_jg = conjg(msex(jg,in))                  ! M_{nv,G'(jg)}
              if (mg_jg == (0.0d0,0.0d0)) cycle
              do ig = 1, ng
                mg_ig = msex(ig,in)                       ! M*_{nv,G(ig)}
                ssex = ssex + mg_ig * wv(ig,jg) * mg_jg
              end do
            end do
            sexk(in) = sexk(in) - ssex
          end do

        end do  ! iv

      end do  ! iq

      ! Apply the shared 1/(N_k Omega) factor and store (real part; the diagonal
      ! self-energy is real up to rounding).
      do in = nb_lo, nb_hi
        if (present(sex_out)) sex_out(in,ik) = rnk * dble(sexk(in))
        if (present(coh_out)) coh_out(in,ik) = rnk * dble(cohk(in))
        sigc(in,ik) = rnk * ( dble(sexk(in)) + dble(cohk(in)) )
      end do

    end do  ! ik

    ! Skipped-(G,G') fraction (q- and k-independent: depends only on the eps
    ! G-list geometry vs the larger sphere).
    if (present(skip_frac)) then
      if (npair_tot > 0) then
        skip_frac = dble(npair_skip) / dble(npair_tot)
      else
        skip_frac = 0.0d0
      end if
    end if

    deallocate(epsinv, vcoul, imap, wv, mblk, msex, mtxl, idiff)
    deallocate(gvecl, ggl, sexk, cohk)

  end subroutine calc_sigma_cohsex

end module gw_sigma_cohsex_sub
