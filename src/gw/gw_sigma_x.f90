!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Bare exchange (Fock) self-energy for the GW module -- ladder rung 1.
!
! Physics
! -------
! The bare-exchange (Fock) self-energy is the static, unscreened limit of the
! GW self-energy (W -> v, the bare Coulomb interaction).  Its diagonal
! expectation value for state (n,k), spin ispin, on a regular k-mesh is
!
!   < n k | Sigma_x | n k > = - (1/N_k) sum_q^BZ sum_v^occ sum_G
!                                 v(q+G) * | M_{nv}(k,q,G) |^2,
!
! with
!   M_{nv}(k,q,G) = < n, k | e^{ i (q+G).r } | v, k-q >
!                 = < u_{n,k} | e^{ i G.r } | u_{v,k-q} >   (the e^{i q.r} cancels
!                   against the Bloch phases of the bra/ket; same as in chi0).
!
! Here N_k = number of k-points (system%nk), v(q+G) = 4 pi / |q+G|^2 (the bare
! Coulomb kernel from build_vcoul; its q->0, G=0 head is the mini-BZ average),
! and the q-sum runs over the BZ k-mesh.  Only OCCUPIED bands v contribute, so
! Sigma_x is independent of the number of empty bands carried in the run.
!
! Sign / magnitude: |M|^2 >= 0 and v(q+G) > 0, so Sigma_x < 0 (attractive
! exchange).  For valence states of a semiconductor it is O(-10 .. -25 eV).
! The exchange correction Sigma_x - V_xc widens the gap toward Hartree-Fock.
!
! Normalisation:  1/N_k, NOT 1/Omega
! ----------------------------------
! The mesh prefactor here is 1/N_k with NO explicit cell-volume factor.  The
! reasoning: the Brillouin-zone integral that defines Fock exchange is
!   (1/(2 pi)^3) integral_BZ dq ...  ->  (1/(N_k Omega)) sum_q ... ,
! because each mesh point owns a BZ volume (2 pi)^3 / (N_k Omega).  The 1/Omega
! that this produces is cancelled exactly by the Omega carried in v(q+G):
! v(q+G) = 4 pi / |q+G|^2 is the bare Coulomb kernel per unit cell charge, and
! the matrix element M is a CELL-AVERAGED (dimensionless) quantity
! [M_{nn}(k,0,0) = 1, pinned by the calc_mtxel q=0 check].  Collecting factors,
! only 1/N_k survives.  (Contrast chi0 in gw_epsilon, which keeps an explicit
! 1/Omega: there M appears once with conjg(M)*M but the response is a density-
! density correlation that carries an extra 1/Omega; the Fock self-energy does
! not.)  This is verified independently by the n_empty invariance and the
! O(-10 eV) magnitude in the [gw][t5] sanity block.
!
! k-q folding (umklapp) and the calc_mtxel bridge
! -----------------------------------------------
! The ket lives at k-q.  We obtain its stored mesh index ikm and umklapp G0 from
!   find_kpq(system, ik, -qvec, ikm, g0vec)   =>   k - q = k_{ikm} + G0,
! i.e. find_kpq is called with -q (it folds k + (argument)).  Then
! u_{v,k-q}(r) = e^{-i G0.r} u_{v,ikm}(r).
!
! calc_mtxel(..., ik, ikm, ...) returns the e^{-iG'.r} kernel with the BRA at
! ikm and the KET at ik:
!   mblk(ig', nb, nk2) = < u_{nb, ikm} | e^{-i G'(ig').r} | u_{nk2, ik} >.
! Put nb = v (occupied, at k-q) and nk2 = n (output state, at k).  Because |M|^2
! is taken, we may work with M*_{nv} :
!   M*_{nv,G} = < u_{v,k-q} | e^{-i G.r} | u_{n,k} >
!             = < u_{v,ikm} | e^{ i G0.r} e^{-i G.r} | u_{n,k} >
!             = < u_{v,ikm} | e^{-i (G - G0).r} | u_{n,k} >
!             = mblk( gindex_of(G - G0), v, n ).
! Hence the per-q index map is
!   imap(ig) = gindex_of( gvec(:,ig) - g0vec ),
! and  | M_{nv}(k,q,G(ig)) |^2 = | mblk(imap(ig), v, n) |^2.
!
! Note the SIGN difference from chi0's map.  In chi0 the bra sits at the FOLDED
! point (k+q) and the kernel target is -(G + G0); here the bra (v) sits at the
! folded point k-q reached via -q, and the e^{-iG.r} convention plus the G0
! sign flip give +(G - G0) instead.  For G0 = 0 both reduce to the identity map
! (imap(ig) = ig).  Entries whose (G - G0) leaves the cutoff sphere are
! unavailable and contribute zero (boundary truncation; harmless deep inside
! the sphere where v(q+G) is largest).
!
! Parallel layout
! ---------------
! Launched with nproc_k = 1 (all k local), so every (ik, ikm) pair is available
! and calc_mtxel (collective over icomm_rko) returns the full block on every
! rank.  The Fock accumulation is then a deterministic local reduction; no extra
! communication is needed beyond what calc_mtxel already performs.
!
! No proper nouns appear in this file.
!
module gw_sigma_x_sub
  implicit none
  private

  public :: calc_sigma_x

  real(8), parameter :: occ_thr   = 1.0d-6   ! occupied if rocc > occ_thr
  real(8), parameter :: g_match_tol = 1.0d-6 ! |G|-match tolerance (bohr^-1)^2
  integer, parameter :: nsub_head = 10       ! mini-BZ sub-sampling for v head

contains

  ! --------------------------------------------------------------------------
  ! gindex_of
  !
  ! List index jg with gvec(:,jg) == gtarget (within tolerance), 0 if gtarget
  ! lies outside the cutoff sphere.  (Local copy; identical role to the helper
  ! in gw_epsilon, kept private here so this module is self-contained.)
  ! --------------------------------------------------------------------------
  integer function gindex_of(ng, gvec, gtarget) result(jg)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in) :: gvec(3,ng)
    real(8), intent(in) :: gtarget(3)
    integer :: i
    real(8) :: d1, d2, d3, d2sum
    jg = 0
    do i = 1, ng
      d1 = gvec(1,i) - gtarget(1)
      d2 = gvec(2,i) - gtarget(2)
      d3 = gvec(3,i) - gtarget(3)
      d2sum = d1*d1 + d2*d2 + d3*d3
      if (d2sum < g_match_tol) then
        jg = i
        return
      end if
    end do
  end function gindex_of

  ! --------------------------------------------------------------------------
  ! calc_sigma_x
  !
  ! Diagonal bare-exchange self-energy expectation < n k | Sigma_x | n k > for a
  ! band window [nb_lo, nb_hi] and all k-points, spin ispin.  Returns real,
  ! negative values (a.u.).
  !
  !   sigx(n,k) = - (1/N_k) sum_q sum_v^occ sum_G v(q+G) |M_{nv}(k,q,G)|^2.
  !
  ! Loop structure: for each k (output states n) and each momentum transfer q
  ! (taken as the mesh-point differences q = vec_k(:,iq) - vec_k(:,ik), so the
  ! q-mesh is the full BZ mesh and includes q = 0 with its mini-BZ head), find
  ! the folded ket index ikm = k - q via find_kpq(ik, -q), fetch the M block
  ! once with calc_mtxel, then accumulate over occupied v and all G.
  !
  ! Arguments:
  !   system  [in]  -- vec_k, rocc, det_a (Omega), nk, no, nspin, hvol, lattice.
  !   info    [in]  -- parallel ownership / communicator (passed to calc_mtxel).
  !   mg, lg  [in]  -- real-space grid descriptors (passed to calc_mtxel).
  !   spsi    [in]  -- cell-periodic orbitals.
  !   gvec(3,ng)[in]-- Cartesian G vectors (bohr^-1), from build_gvectors.
  !   gg(ng)  [in]  -- |G|^2 (bohr^-2), from build_gvectors.
  !   ng      [in]  -- number of G vectors.
  !   ispin   [in]  -- spin channel (1 for the non-magnetic case).
  !   nb_lo,nb_hi[in]-- output band window (states n).
  !   sigx(nb_lo:nb_hi, system%nk) [out] -- exchange expectation (a.u.).
  ! --------------------------------------------------------------------------
  subroutine calc_sigma_x(system, info, mg, lg, spsi, gvec, gg, ng, &
                          ispin, nb_lo, nb_hi, sigx)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,   only: calc_mtxel
    use gw_coulomb_sub, only: build_vcoul
    use gw_epsilon_sub, only: find_kpq
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_orbital),       intent(in)  :: spsi
    integer,               intent(in)  :: ng
    real(8),               intent(in)  :: gvec(3,ng)
    real(8),               intent(in)  :: gg(ng)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nb_lo, nb_hi
    real(8),               intent(out) :: sigx(nb_lo:nb_hi, system%nk)

    complex(8), allocatable :: mblk(:,:,:)   ! < u_{.,ikm}|e^{-iG.r}|u_{.,ik} >
    real(8),    allocatable :: vcoul(:)
    integer,    allocatable :: imap(:)

    integer :: no, nk, nq, ik, iq, ikm, iv, in, ig, jg
    real(8) :: qvec(3), mqvec(3), g0vec(3), gtarget(3)
    real(8) :: omega, rnk, acc, m2

    no    = system%no
    nk    = system%nk
    nq    = system%nk
    omega = abs(system%det_a)
    rnk   = 1.0d0 / (dble(nk) * omega)  ! Fock mesh+cell factor 1/(N_k*Omega):
                                        ! v(q+G)=4pi/|q+G|^2 (build_vcoul) carries
                                        ! no 1/Omega, so the cell volume enters here

    sigx(:,:) = 0.0d0

    allocate(mblk(ng, no, no), vcoul(ng), imap(ng))

    ! ---- output-state loop (n at k) ---------------------------------------
    do ik = 1, nk

      ! ---- BZ q-sum: q = vec_k(:,iq) - vec_k(:,ik) ------------------------
      do iq = 1, nq
        qvec(1:3)  =  system%vec_k(1:3,iq) - system%vec_k(1:3,ik)
        mqvec(1:3) = -qvec(1:3)

        ! ket index of k - q (= k + (-q) folded) and umklapp G0.
        call find_kpq(system, ik, mqvec, ikm, g0vec)
        if (ikm == 0) cycle            ! no mesh partner (symmetry-reduced mesh)

        ! v(q+G): bare Coulomb kernel; q->0 head averaged over the mini-BZ.
        call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

        ! Per-(k,q) index map: imap(ig) -> calc_mtxel column giving M*_{nv,G(ig)}
        ! for the OUTPUT state n (ket at ik) and occupied v (bra at ikm).
        ! M*_{nv,G} = < u_{v,ikm}| e^{-i(G - G0).r} | u_{n,ik} >, so the e^{-iG'.r}
        ! convention needs G' = G - G0 (note the + sign on gvec, - on g0vec).
        do ig = 1, ng
          gtarget(1) = gvec(1,ig) - g0vec(1)
          gtarget(2) = gvec(2,ig) - g0vec(2)
          gtarget(3) = gvec(3,ig) - g0vec(3)
          imap(ig) = gindex_of(ng, gvec, gtarget)   ! 0 if (G-G0) outside cutoff
        end do

        ! Full M block at (ik, ikm): mblk(ig', nb=bra@ikm, nk2=ket@ik).
        ! calc_mtxel is collective -> entered on all ranks.
        call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikm, ispin, &
                        no, no, mblk)

        ! ---- occupied-v sum and G sum -------------------------------------
        do iv = 1, no
          if (system%rocc(iv,ikm,ispin) <= occ_thr) cycle   ! v occupied at k-q

          do in = nb_lo, nb_hi
            acc = 0.0d0
            do ig = 1, ng
              jg = imap(ig)
              if (jg == 0) cycle                ! (G-G0) outside cutoff: truncate
              ! |M_{nv}(k,q,G(ig))|^2 = |mblk(imap(ig), v, n)|^2
              m2  = dble( mblk(jg,iv,in) * conjg(mblk(jg,iv,in)) )
              acc = acc + vcoul(ig) * m2
            end do
            ! - (1/N_k) accumulation (q and v summed across iterations)
            sigx(in,ik) = sigx(in,ik) - rnk * acc
          end do

        end do  ! iv

      end do  ! iq
    end do  ! ik

    deallocate(mblk, vcoul, imap)

  end subroutine calc_sigma_x

end module gw_sigma_x_sub
