!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Static (omega=0) RPA dielectric matrix for the GW module.
!
! Physics
! -------
! The static random-phase-approximation polarizability (independent-particle
! response) in a plane-wave basis, at omega=0, is:
!
!   chi0_{G G'}(q,0) = (1/Omega) * sum_k sum_{v,c}
!         ( f_{v k} - f_{c,k+q} ) * conjg(M_{cv,G}(k,q)) * M_{cv,G'}(k,q)
!         / ( eps_{v k} - eps_{c,k+q} )
!
! where v runs over occupied bands, c over unoccupied bands, f the occupation
! (system%rocc, so the spin multiplicity is already folded in), eps the
! Kohn-Sham eigenvalue (energy%esp), Omega the cell volume (system%det_a), and
!
!   M_{cv,G}(k,q) = < u_{c,k+q} | e^{i G.r} | u_{v,k} >.
!
! Sign / prefactor sanity (occupied v, unoccupied c):
!   numerator   (f_v - f_c) > 0          (f_v ~ 2, f_c ~ 0)
!   denominator (eps_v - eps_c) < 0      (occupied lies below unoccupied)
!   diagonal weight conjg(M)*M = |M|^2 >= 0
! so chi0_{GG} < 0 (the response is "negative-definite-ish"), and therefore
!   eps_{GG}(q) = 1 - v(q+G)*chi0_{GG}(q) = 1 - (positive)*(negative) > 1,
! i.e. the macroscopic dielectric constant exceeds 1, as it must.
!
! Convention bridge to calc_mtxel
! -------------------------------
! calc_mtxel returns < u_{nb,ikq} | e^{-i G.r} | u_{nk2,ik} >, i.e. the
! e^{-iG.r} kernel.  With nb=c at ikq, nk2=v at ik this equals M_{cv,-G(ig)}.
! Because the G-vector list from build_gvectors is closed under G->-G, we form
! the permutation iGm(ig) = (index of -gvec(:,ig)) and then
!   M_{cv,G(ig)} = calc_mtxel-result( iGm(ig) ).
!
! k+q folding (umklapp)
! ---------------------
! For a momentum transfer q the response needs band c at k+q.  On the stored
! mesh, k+q equals some mesh point k_{ikq} plus a reciprocal-lattice vector G0
! (the umklapp vector):  k + q = k_{ikq} + G0.  Then
!   < u_{c,k+q} | e^{i G.r} | u_{v,k} > = < u_{c,k_{ikq}} | e^{i (G+G0).r} | u_{v,k} >,
! so the requested G is rigidly shifted by G0.  In the calc_mtxel (e^{-iG.r})
! convention the entry that yields M_{cv,G(ig)} for folded q is the gvec index
! of -(G(ig)+G0).  We precompute that mapping per q (imap); entries whose
! -(G+G0) falls outside the cutoff sphere are unavailable and set to zero
! (boundary truncation).  For G0=0 the map reduces to iGm and is always exact.
!
! find_kpq locates ikq and G0 directly from the Cartesian k-vectors
! system%vec_k (bohr^-1, filled by init_kvector as vec_k = primitive_b . k_frac):
! it searches the stored mesh for the point whose offset from k+q is a
! reciprocal-lattice vector.  This is robust to the storage order of vec_k.  On
! a full (unreduced) Gamma-centred mesh -- the layout the GW controller uses
! (use_symmetry off) -- every k+q has such a partner.  If the mesh has been
! symmetry-reduced so that no partner exists, find_kpq returns ikq=0 and the
! caller flags/skips that q.
!
! Dielectric matrix and inverse
! -----------------------------
!   eps_{G G'}(q) = delta_{G G'} - v(q+G) * chi0_{G G'}(q)
! with v(q+G) (the bare Coulomb kernel from build_vcoul) multiplying the ROW
! index G.  The inverse is obtained with LAPACK zgetrf + zgetri (the same calls
! used by matrix_inverse_complex in src/math/salmon_math.f90).
!
! Parallel note
! -------------
! The GW job is launched with a single MPI rank (nproc_k=1, all k local), so
! every rank holds every orbital and calc_mtxel works for any (ik,ikq).
! calc_mtxel is collective (comm_summation), so it is entered on all ranks; the
! resulting global M, and hence chi0, eps and its inverse, are identical on
! every rank (the build is deterministic), so no extra broadcast is needed.
!
! No proper nouns appear in this file.
!
module gw_epsilon_sub
  implicit none
  private

  public :: build_gminus
  public :: find_kpq
  public :: calc_epsinv
  public :: calc_chi0_freq
  public :: calc_w_freq

  real(8), parameter :: g_match_tol = 1.0d-6   ! |G|-match tolerance (bohr^-1)^2
  real(8), parameter :: occ_thr     = 1.0d-6   ! occupied if rocc > occ_thr
  real(8), parameter :: denom_tiny  = 1.0d-8   ! skip near-degenerate transitions

contains

  ! --------------------------------------------------------------------------
  ! build_gminus
  !
  ! For each G(ig) find the list index of -G and return it in iGm(ig).  The
  ! gvec list from build_gvectors is closed under negation (the search loops a
  ! symmetric integer box and keeps |G|^2 <= cutoff), so a partner must exist
  ! for every ig; a missing partner is a fatal inconsistency and is reported.
  !
  ! Arguments:
  !   ng         [in]  -- number of G vectors.
  !   gvec(3,ng) [in]  -- Cartesian G vectors (bohr^-1).
  !   iGm(ng)    [out] -- index of -gvec(:,ig); iGm(ig)=0 on (unexpected) miss.
  ! --------------------------------------------------------------------------
  subroutine build_gminus(ng, gvec, iGm)
    implicit none
    integer, intent(in)  :: ng
    real(8), intent(in)  :: gvec(3,ng)
    integer, intent(out) :: iGm(ng)

    integer :: ig, jg
    real(8) :: d1, d2, d3, d2sum
    logical :: found

    iGm(:) = 0
    do ig = 1, ng
      found = .false.
      do jg = 1, ng
        d1 = gvec(1,jg) + gvec(1,ig)
        d2 = gvec(2,jg) + gvec(2,ig)
        d3 = gvec(3,jg) + gvec(3,ig)
        d2sum = d1*d1 + d2*d2 + d3*d3
        if (d2sum < g_match_tol) then
          iGm(ig) = jg
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        write(*,*) "[gw][t4] FATAL: build_gminus found no -G partner for ig=", ig
      end if
    end do

  end subroutine build_gminus

  ! --------------------------------------------------------------------------
  ! find_kpq
  !
  ! Locate the stored mesh index ikq and umklapp vector g0vec such that
  !   vec_k(:,ik) + qvec = vec_k(:,ikq) + g0vec,
  ! where g0vec is a reciprocal-lattice vector (an integer combination of the
  ! reciprocal primitive vectors system%primitive_b).  All vectors are
  ! Cartesian (bohr^-1).
  !
  ! Strategy: for each candidate ikq, form the residual
  !   res = vec_k(:,ik) + qvec - vec_k(:,ikq)
  ! and test whether res is (numerically) a reciprocal-lattice vector by
  ! projecting onto the dual (real primitive) basis: the fractional components
  ! n_j = res . a_j / (2 pi) must be (near) integers, because b_i . a_j =
  ! 2 pi delta_ij.  If so, g0vec = round(n) . b and we are done.
  !
  ! Arguments:
  !   system     [in]  -- provides vec_k, nk, primitive_a, primitive_b.
  !   ik         [in]  -- index of k.
  !   qvec(3)    [in]  -- momentum transfer (Cartesian, bohr^-1).
  !   ikq        [out] -- stored index of the folded k+q; 0 if none found.
  !   g0vec(3)   [out] -- umklapp vector (Cartesian, bohr^-1); 0 if ikq=0.
  ! --------------------------------------------------------------------------
  subroutine find_kpq(system, ik, qvec, ikq, g0vec)
    use structures, only: s_dft_system
    implicit none
    type(s_dft_system), intent(in)  :: system
    integer,            intent(in)  :: ik
    real(8),            intent(in)  :: qvec(3)
    integer,            intent(out) :: ikq
    real(8),            intent(out) :: g0vec(3)

    integer :: jk, j
    real(8) :: kpq(3), res(3), fn(3), pi, twopi
    integer :: n1, n2, n3
    real(8) :: r1, r2, r3

    pi    = acos(-1.0d0)
    twopi = 2.0d0 * pi

    ikq      = 0
    g0vec(:) = 0.0d0

    kpq(1) = system%vec_k(1,ik) + qvec(1)
    kpq(2) = system%vec_k(2,ik) + qvec(2)
    kpq(3) = system%vec_k(3,ik) + qvec(3)

    do jk = 1, system%nk
      res(1) = kpq(1) - system%vec_k(1,jk)
      res(2) = kpq(2) - system%vec_k(2,jk)
      res(3) = kpq(3) - system%vec_k(3,jk)

      ! fractional reciprocal components: n_j = res . a_j / (2 pi)
      ! (primitive_a(:,j) is the j-th real primitive vector; b_i . a_j = 2 pi delta_ij)
      do j = 1, 3
        fn(j) = ( res(1)*system%primitive_a(1,j) &
                + res(2)*system%primitive_a(2,j) &
                + res(3)*system%primitive_a(3,j) ) / twopi
      end do

      if ( abs(fn(1) - anint(fn(1))) < 1.0d-5 .and. &
           abs(fn(2) - anint(fn(2))) < 1.0d-5 .and. &
           abs(fn(3) - anint(fn(3))) < 1.0d-5 ) then
        n1 = nint(fn(1)); n2 = nint(fn(2)); n3 = nint(fn(3))
        ! reconstruct g0 = sum_j n_j * b_j (Cartesian)
        r1 = dble(n1)*system%primitive_b(1,1) + dble(n2)*system%primitive_b(1,2) &
           + dble(n3)*system%primitive_b(1,3)
        r2 = dble(n1)*system%primitive_b(2,1) + dble(n2)*system%primitive_b(2,2) &
           + dble(n3)*system%primitive_b(2,3)
        r3 = dble(n1)*system%primitive_b(3,1) + dble(n2)*system%primitive_b(3,2) &
           + dble(n3)*system%primitive_b(3,3)
        ikq      = jk
        g0vec(1) = r1
        g0vec(2) = r2
        g0vec(3) = r3
        return
      end if
    end do

  end subroutine find_kpq

  ! --------------------------------------------------------------------------
  ! gindex_of
  !
  ! Return the list index jg with gvec(:,jg) == gtarget (within tolerance), or
  ! 0 if gtarget lies outside the cutoff sphere / is not on the list.
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
  ! calc_epsinv
  !
  ! Assemble the static RPA chi0 for one momentum transfer q, form the
  ! dielectric matrix eps = 1 - v(q+G) chi0, and invert it with LAPACK.
  !
  ! The k-loop runs over all stored k-points; for each k the partner ikq and
  ! umklapp g0vec are obtained from find_kpq.  If a partner does not exist
  ! (symmetry-reduced mesh) the whole q is unrepresentable: epsinv is set to the
  ! identity and ok=.false. is returned so the caller can flag/skip it.
  !
  ! For each k the occupied (v) and unoccupied (c) band sets are determined from
  ! rocc, the full block of M_{cv,G} for that (k,ikq) is fetched once with
  ! calc_mtxel (e^{-iG.r} convention; remapped to M_{cv,G} via imap that also
  ! carries the umklapp shift), and the static polarizability sum is accumulated.
  !
  ! Arguments:
  !   system   [in]  -- vec_k, rocc, det_a (Omega), nk, no, nspin, hvol, lattice.
  !   info     [in]  -- parallel ownership / communicator (passed to calc_mtxel).
  !   mg, lg   [in]  -- real-space grid descriptors (passed to calc_mtxel).
  !   spsi     [in]  -- cell-periodic orbitals.
  !   esp(no,nk,nspin)[in] -- Kohn-Sham eigenvalues (a.u.); pass energy%esp.
  !   gvec(3,ng)[in] -- Cartesian G vectors (bohr^-1), from build_gvectors.
  !   gg(ng)   [in]  -- |G|^2 (bohr^-2), from build_gvectors.
  !   ng       [in]  -- number of G vectors (matrix dimension).
  !   iq       [in]  -- q index (for diagnostics only).
  !   qvec(3)  [in]  -- momentum transfer (Cartesian, bohr^-1).
  !   ispin    [in]  -- spin channel (1 for the non-magnetic case).
  !   epsinv(ng,ng)[out] -- inverse dielectric matrix eps^{-1}(q).
  !   eps_diag(ng) [out,opt] -- diagonal of eps (before inversion), for checks.
  !   residual [out,opt]     -- max| eps . eps^{-1} - I |, an inversion check.
  !   ok       [out,opt]     -- .true. if q was representable and inverted.
  ! --------------------------------------------------------------------------
  subroutine calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, iq, qvec, &
                         ispin, epsinv, eps_diag, residual, ok, local_only)
    use structures,    only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,  only: calc_mtxel
    use gw_coulomb_sub,only: build_vcoul
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
    integer,               intent(in)  :: iq
    real(8),               intent(in)  :: qvec(3)
    integer,               intent(in)  :: ispin
    complex(8),            intent(out) :: epsinv(ng,ng)
    real(8),    optional,  intent(out) :: eps_diag(ng)
    real(8),    optional,  intent(out) :: residual
    logical,    optional,  intent(out) :: ok
    ! Forwarded to calc_mtxel: when set, the orbitals are replicated and this
    ! rank builds the whole chi0/eps locally (no inner collective).
    logical,    optional,  intent(in)  :: local_only

    complex(8), allocatable :: chi0(:,:), eps(:,:), mtxel(:,:,:)
    complex(8), allocatable :: mcv(:,:)        ! M_{cv,G}(ig, ipair), remapped
    complex(8), allocatable :: prod(:,:)
    real(8),    allocatable :: vcoul(:)
    integer,    allocatable :: iGm(:), imap(:)
    integer,    allocatable :: ipiv(:)
    complex(8), allocatable :: zwork(:)

    integer :: no, nk, ik, ikq, iv, ic, ig, jg, ipair, npair, nfill
    integer :: lwork, linfo
    integer :: nsub_head
    real(8) :: g0vec(3), gtarget(3), omega, inv_omega
    real(8) :: fv, fc, dw, de, wgt
    complex(8) :: zfac
    logical :: q_ok
    logical :: ll

    ll = .false.
    if (present(local_only)) ll = local_only

    no    = system%no
    nk    = system%nk
    omega = abs(system%det_a)
    ! chi0 prefactor 2/(N_k*Omega).  Two factors of 2 enter the static RPA chi0
    ! summed once over occ->unocc: the spin 2 (here in dw = f_v-f_c = 2) and the
    ! static-response 2 (the antiresonant pole, lumped into the prefactor).  The
    ! reference convention is fact = 4/(N_k*Omega*nspin) summed occ->unocc once;
    ! with dw=2 this needs inv_omega = 2/(N_k*Omega).  (Confirmed by the f-sum
    ! rule: with 1/(N_k Omega) the absorption f-sum was ~0.48 of (pi/2)wp^2.)
    ! 1/N_k was also missing earlier (eps_M grew with mesh at fixed |q|).
    inv_omega = 2.0d0 / (dble(nk) * omega)
    nsub_head = 10            ! mini-BZ sub-sampling density for the v head

    allocate(chi0(ng,ng), eps(ng,ng))
    allocate(iGm(ng), imap(ng), vcoul(ng))
    chi0(:,:) = (0.0d0, 0.0d0)

    ! G -> -G permutation (closed list; asserted in build_gminus).  Used directly
    ! for the unfolded (G0=0) case and as the closure assertion: if any partner is
    ! missing the gvec list is not closed under negation, which would corrupt M.
    call build_gminus(ng, gvec, iGm)
    if (any(iGm(:) == 0)) then
      write(*,*) "[gw][t4] FATAL: gvec list not closed under G->-G (iq=", iq, ")"
    end if

    ! v(q+G), q-dependent (head averaged over the mini-BZ for q->0)
    call build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub_head, vcoul)

    q_ok = .true.

    ! --- k-loop: accumulate chi0 -------------------------------------------
    do ik = 1, nk

      call find_kpq(system, ik, qvec, ikq, g0vec)
      if (ikq == 0) then
        q_ok = .false.
        exit                                   ! q not representable on this mesh
      end if

      ! Per-(k,q) index map: imap(ig) = calc_mtxel column that gives M_{cv,G(ig)}.
      ! M_{cv,G}(folded) = < u_{c,ikq} | e^{i (G+G0).r} | u_{v,k} >, and calc_mtxel
      ! stores the e^{-iG'.r} kernel, so we need the gvec index of -(G(ig)+G0).
      do ig = 1, ng
        gtarget(1) = -( gvec(1,ig) + g0vec(1) )
        gtarget(2) = -( gvec(2,ig) + g0vec(2) )
        gtarget(3) = -( gvec(3,ig) + g0vec(3) )
        jg = gindex_of(ng, gvec, gtarget)
        imap(ig) = jg                          ! 0 if -(G+G0) is outside cutoff
      end do

      ! Fetch the full M block for this (ik,ikq): mtxel(ig, c, v).
      ! calc_mtxel is collective -> entered on all ranks unconditionally.
      allocate(mtxel(ng, no, no))
      call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikq, ispin, &
                      no, no, mtxel, local_only=ll)

      ! Remapped, umklapp-shifted M for every occ->unocc pair at this k.
      ! Count pairs first.
      npair = 0
      do iv = 1, no
        if (system%rocc(iv,ik,ispin) <= occ_thr) cycle
        do ic = 1, no
          if (system%rocc(ic,ikq,ispin) > occ_thr) cycle
          npair = npair + 1
        end do
      end do

      if (npair > 0) then
        ! npair upper-bounds the allocation; nfill counts the pairs actually
        ! built (near-degenerate transitions are skipped), and only the first
        ! nfill columns of mcv/prod are valid.
        allocate(mcv(ng, npair), prod(ng, npair))
        ipair = 0
        do iv = 1, no
          fv = system%rocc(iv,ik,ispin)
          if (fv <= occ_thr) cycle
          do ic = 1, no
            fc = system%rocc(ic,ikq,ispin)
            if (fc > occ_thr) cycle
            de = esp(iv,ik,ispin) - esp(ic,ikq,ispin)
            if (abs(de) < denom_tiny) cycle    ! skip ill-defined transition

            ipair = ipair + 1
            dw  = fv - fc
            wgt = inv_omega * dw / de           ! (1/Omega)(f_v-f_c)/(e_v-e_c)

            ! Build M_{cv,G(ig)} for this pair, applying the G->-G+umklapp map.
            ! mtxel(jg, ic, iv) = M_{cv,-G(jg)}; with jg=imap(ig) this is M_{cv,G(ig)}.
            do ig = 1, ng
              jg = imap(ig)
              if (jg == 0) then
                mcv(ig,ipair) = (0.0d0, 0.0d0)  ! G+G0 outside cutoff -> truncate
              else
                mcv(ig,ipair) = mtxel(jg, ic, iv)
              end if
            end do

            ! Pre-scale a working copy by the (real) transition weight so that
            ! chi0 += conjg(mcv) (outer) (wgt*mcv).
            do ig = 1, ng
              prod(ig,ipair) = wgt * mcv(ig,ipair)
            end do
          end do
        end do
        nfill = ipair                          ! pairs actually built

        ! Rank-update: chi0_{G G'} += sum_pair conjg(mcv(G,p)) * prod(G',p).
        ! (prod already carries the (1/Omega)(f_v-f_c)/(e_v-e_c) weight.)
        do ipair = 1, nfill
          do jg = 1, ng
            zfac = prod(jg,ipair)
            if (zfac == (0.0d0,0.0d0)) cycle
            do ig = 1, ng
              chi0(ig,jg) = chi0(ig,jg) + conjg(mcv(ig,ipair)) * zfac
            end do
          end do
        end do

        deallocate(mcv, prod)
      end if

      deallocate(mtxel)

    end do  ! ik

    ! --- not representable: return identity and flag -----------------------
    if (.not. q_ok) then
      epsinv(:,:) = (0.0d0, 0.0d0)
      do ig = 1, ng
        epsinv(ig,ig) = (1.0d0, 0.0d0)
      end do
      if (present(eps_diag)) eps_diag(:) = 1.0d0
      if (present(residual)) residual = 0.0d0
      if (present(ok))       ok = .false.
      deallocate(chi0, eps, iGm, imap, vcoul)
      return
    end if

    ! --- eps_{G G'} = delta_{G G'} - v(q+G) * chi0_{G G'} ------------------
    ! v multiplies the ROW index G.
    do jg = 1, ng
      do ig = 1, ng
        eps(ig,jg) = - vcoul(ig) * chi0(ig,jg)
      end do
      eps(jg,jg) = eps(jg,jg) + (1.0d0, 0.0d0)
    end do

    if (present(eps_diag)) then
      do ig = 1, ng
        eps_diag(ig) = dble(eps(ig,ig))
      end do
    end if

    ! --- LAPACK inversion: epsinv = eps^{-1} (zgetrf + zgetri) -------------
    ! Mirror src/math/salmon_math.f90: factorize then invert in place.  Invert a
    ! copy so the residual check below can use the original eps.
    epsinv(:,:) = eps(:,:)
    lwork = ng * max(ng, 64)
    allocate(ipiv(ng), zwork(lwork))
    call zgetrf(ng, ng, epsinv, ng, ipiv, linfo)
    if (linfo /= 0) then
      write(*,*) "[gw][t4] WARNING: zgetrf failed for iq=", iq, " info=", linfo
    end if
    call zgetri(ng, epsinv, ng, ipiv, zwork, lwork, linfo)
    if (linfo /= 0) then
      write(*,*) "[gw][t4] WARNING: zgetri failed for iq=", iq, " info=", linfo
    end if
    deallocate(ipiv, zwork)

    ! --- inversion residual: max| eps . epsinv - I | (optional) -----------
    if (present(residual)) then
      residual = 0.0d0
      do jg = 1, ng
        do ig = 1, ng
          zfac = (0.0d0, 0.0d0)
          do ipair = 1, ng
            zfac = zfac + eps(ig,ipair) * epsinv(ipair,jg)
          end do
          if (ig == jg) zfac = zfac - (1.0d0, 0.0d0)
          residual = max(residual, abs(zfac))
        end do
      end do
    end if

    if (present(ok)) ok = .true.

    deallocate(chi0, eps, iGm, imap, vcoul)

  end subroutine calc_epsinv


  ! ----------------------------------------------------------------------
  ! Real-frequency RPA polarizability chi0_{GG'}(q,omega) on a real-omega
  ! grid (spec-b1).  Same occ->unocc assembly as calc_epsinv, but the single
  ! static denominator (f_v-f_c)/(e_v-e_c) is generalised to the omega-
  ! dependent two-pole form, normalised so its omega=0, eta->0 limit
  ! reproduces calc_epsinv's static weight EXACTLY (regression gate iii):
  !
  !   wgt(w) = (1/Omega)(f_v-f_c) * (1/2)[ 1/(w - Delta + i eta)
  !                                      - 1/(w + Delta - i eta) ]
  !   Delta = e_c - e_v > 0 ;  at w=0,eta->0 : (1/2)(-1/Delta - 1/Delta)
  !                                          = -1/Delta = 1/(e_v-e_c).
  !
  ! NOTE on the prefactor: the static weight here sums occ->unocc ONCE with a
  ! single denominator, i.e. the 1/2 normalisation of the symmetric two-pole
  ! form.  The absolute scale of Im eps (vs RT-TDDFT, spec gate ii) is what
  ! arbitrates whether this convention needs the historical x2 (the long-open
  ! chi0 prefactor / eps_inf question).  Here we only guarantee the w=0 limit
  ! matches the existing static path so the QP/GPP gap is reproduced.
  !
  ! This returns the BARE polarizability only (no vcoul, no inversion); the
  ! dielectric inversion eps = 1 - v chi0 -> eps^{-1}, W is calc_w_freq.
  ! ----------------------------------------------------------------------
  subroutine calc_chi0_freq(system, info, mg, lg, spsi, esp, gvec, ng, iq, qvec, &
                            ispin, nomega, omega_grid, eta, chi0_w, ok, local_only, run_sanity)
    use structures,    only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_mtxel_sub,  only: calc_mtxel
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_rgrid),         intent(in)  :: lg
    type(s_orbital),       intent(in)  :: spsi
    real(8),               intent(in)  :: esp(system%no, system%nk, system%nspin)
    integer,               intent(in)  :: ng
    real(8),               intent(in)  :: gvec(3,ng)
    integer,               intent(in)  :: iq
    real(8),               intent(in)  :: qvec(3)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nomega
    real(8),               intent(in)  :: omega_grid(nomega)   ! a.u.
    real(8),               intent(in)  :: eta                  ! a.u.
    complex(8),            intent(out) :: chi0_w(ng,ng,nomega)
    logical,    optional,  intent(out) :: ok
    logical,    optional,  intent(in)  :: local_only
    logical,    optional,  intent(in)  :: run_sanity

    complex(8), allocatable :: mtxel(:,:,:)
    complex(8), allocatable :: mcv(:,:)
    real(8),    allocatable :: dwp(:), delp(:)        ! per-pair (f_v-f_c) and Delta
    integer,    allocatable :: iGm(:), imap(:)
    integer :: no, nk, ik, ikq, iv, ic, ig, jg, ipair, npair, nfill, iw, iw0
    real(8) :: g0vec(3), gtarget(3), omega, inv_omega
    real(8) :: fv, fc, de, smax
    complex(8) :: zi, zw, zden, wgt_w, zfac
    logical :: q_ok, ll, dosan

    zi   = (0.0d0, 1.0d0)
    ll   = .false.;  if (present(local_only)) ll = local_only
    dosan= .false.;  if (present(run_sanity)) dosan = run_sanity
    no   = system%no;  nk = system%nk
    omega= abs(system%det_a);  inv_omega = 2.0d0/(dble(nk)*omega)  ! 2/(N_k Omega): spin(in dw)+static-response 2

    allocate(iGm(ng), imap(ng))
    chi0_w(:,:,:) = (0.0d0, 0.0d0)
    call build_gminus(ng, gvec, iGm)
    if (any(iGm(:)==0)) write(*,*) "[gw][chi0w] FATAL: gvec not closed under -G (iq=", iq, ")"
    q_ok = .true.

    iw0 = 0
    do iw = 1, nomega
      if (abs(omega_grid(iw)) < 1.0d-12) iw0 = iw     ! locate the omega=0 point
    end do
    smax = 0.0d0

    do ik = 1, nk
      call find_kpq(system, ik, qvec, ikq, g0vec)
      if (ikq == 0) then
        q_ok = .false.; exit
      end if
      do ig = 1, ng
        gtarget(1) = -( gvec(1,ig) + g0vec(1) )
        gtarget(2) = -( gvec(2,ig) + g0vec(2) )
        gtarget(3) = -( gvec(3,ig) + g0vec(3) )
        imap(ig) = gindex_of(ng, gvec, gtarget)
      end do

      allocate(mtxel(ng, no, no))
      call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ikq, ispin, &
                      no, no, mtxel, local_only=ll)

      npair = 0
      do iv = 1, no
        if (system%rocc(iv,ik,ispin) <= occ_thr) cycle
        do ic = 1, no
          if (system%rocc(ic,ikq,ispin) > occ_thr) cycle
          npair = npair + 1
        end do
      end do

      if (npair > 0) then
        allocate(mcv(ng,npair), dwp(npair), delp(npair))
        ! Pass 1: build remapped M and the (omega-independent) per-pair data.
        ipair = 0
        do iv = 1, no
          fv = system%rocc(iv,ik,ispin)
          if (fv <= occ_thr) cycle
          do ic = 1, no
            fc = system%rocc(ic,ikq,ispin)
            if (fc > occ_thr) cycle
            de = esp(iv,ik,ispin) - esp(ic,ikq,ispin)
            if (abs(de) < denom_tiny) cycle
            ipair = ipair + 1
            dwp(ipair)  = fv - fc
            delp(ipair) = -de                       ! Delta = e_c - e_v > 0
            do ig = 1, ng
              jg = imap(ig)
              if (jg == 0) then
                mcv(ig,ipair) = (0.0d0, 0.0d0)
              else
                mcv(ig,ipair) = mtxel(jg, ic, iv)
              end if
            end do
          end do
        end do
        nfill = ipair

        ! Pass 2: per-omega weighted rank-1 accumulation.  OMP over omega -- each
        ! omega writes a disjoint chi0_w(:,:,iw) slice, so there is no race; this
        ! is the work-dominant loop (nomega * nfill * ng^2 per k).
!$omp parallel do default(shared) schedule(dynamic) &
!$omp   private(iw, zw, ipair, zden, wgt_w, jg, ig, zfac) reduction(max:smax)
        do iw = 1, nomega
          zw = cmplx(omega_grid(iw), 0.0d0, 8)
          do ipair = 1, nfill
            zden  = 0.5d0*( 1.0d0/(zw - delp(ipair) + zi*eta) &
                         -  1.0d0/(zw + delp(ipair) - zi*eta) )
            wgt_w = inv_omega * dwp(ipair) * zden
            if (dosan .and. iw == iw0) &
              smax = max(smax, abs(wgt_w - cmplx(-inv_omega*dwp(ipair)/delp(ipair), 0.0d0, 8)))
            do jg = 1, ng
              zfac = wgt_w * mcv(jg,ipair)
              if (zfac == (0.0d0,0.0d0)) cycle
              do ig = 1, ng
                chi0_w(ig,jg,iw) = chi0_w(ig,jg,iw) + conjg(mcv(ig,ipair)) * zfac
              end do
            end do
          end do
        end do
!$omp end parallel do

        deallocate(mcv, dwp, delp)
      end if

      deallocate(mtxel)
    end do  ! ik

    if (.not. q_ok) chi0_w(:,:,:) = (0.0d0, 0.0d0)
    if (present(ok)) ok = q_ok
    if (dosan .and. iw0 > 0) &
      write(*,'(A,ES12.4)') "  [gw][chi0w] max|wgt_w(0)-wgt_static| = ", smax
    if (dosan .and. iw0 == 0) &
      write(*,*) "  [gw][chi0w] WARN: omega_grid has no 0 point; static gate skipped"
    deallocate(iGm, imap)
  end subroutine calc_chi0_freq


  ! ----------------------------------------------------------------------
  ! Real-frequency dielectric inversion (spec-b1).  For each omega:
  !   eps_{GG'}(q,w) = delta_{GG'} - v(q+G) chi0_{GG'}(q,w)   (v on the ROW
  !                    index G, the same convention as calc_epsinv)
  ! invert with LAPACK (zgetrf+zgetri) -> epsinv_w(:,:,iw).  vcoul is returned
  ! so callers can form W = epsinv . v (the absorption head needs only
  ! eps_M = 1/epsinv_w(0,0); the correlation self-energy needs W).  chi0_w is
  ! supplied by calc_chi0_freq.  At w=0,eta->0 the result equals calc_epsinv's
  ! static epsinv to machine precision (regression gate iii at the W level).
  ! ----------------------------------------------------------------------
  subroutine calc_w_freq(system, gvec, gg, ng, qvec, nomega, chi0_w, epsinv_w, vcoul, ok)
    use structures,    only: s_dft_system
    use gw_coulomb_sub,only: build_vcoul
    implicit none
    type(s_dft_system), intent(in)  :: system
    integer,            intent(in)  :: ng
    real(8),            intent(in)  :: gvec(3,ng)
    real(8),            intent(in)  :: gg(ng)
    real(8),            intent(in)  :: qvec(3)
    integer,            intent(in)  :: nomega
    complex(8),         intent(in)  :: chi0_w(ng,ng,nomega)
    complex(8),         intent(out) :: epsinv_w(ng,ng,nomega)
    real(8),            intent(out) :: vcoul(ng)
    logical, optional,  intent(out) :: ok

    complex(8), allocatable :: epsm(:,:), zwork(:)
    integer,    allocatable :: ipiv(:)
    integer :: ig, jg, iw, lwork, linfo, nsub_head
    real(8) :: omega

    omega = abs(system%det_a)
    nsub_head = 10
    call build_vcoul(ng, gvec, gg, qvec, omega, system%nk, nsub_head, vcoul)

    allocate(epsm(ng,ng), ipiv(ng))
    lwork = ng * max(ng, 64)
    allocate(zwork(lwork))

    do iw = 1, nomega
      do jg = 1, ng
        do ig = 1, ng
          epsm(ig,jg) = - vcoul(ig) * chi0_w(ig,jg,iw)
        end do
        epsm(jg,jg) = epsm(jg,jg) + (1.0d0, 0.0d0)
      end do
      call zgetrf(ng, ng, epsm, ng, ipiv, linfo)
      if (linfo /= 0) write(*,*) "[gw][wfreq] zgetrf info=", linfo, " iw=", iw
      call zgetri(ng, epsm, ng, ipiv, zwork, lwork, linfo)
      if (linfo /= 0) write(*,*) "[gw][wfreq] zgetri info=", linfo, " iw=", iw
      epsinv_w(:,:,iw) = epsm(:,:)
    end do

    deallocate(epsm, ipiv, zwork)
    if (present(ok)) ok = .true.
  end subroutine calc_w_freq

end module gw_epsilon_sub
