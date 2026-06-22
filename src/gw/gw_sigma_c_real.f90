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

  real(8), parameter :: ha2ev   = 27.21138505d0  ! eV per Hartree
  real(8), parameter :: dE_ev   = 0.1d0          ! finite-diff step dE for Z (eV)

contains

  ! --------------------------------------------------------------------------
  ! sigma_c_real_at_energy
  !
  ! Spectral-integral Sigma_c at one probe energy E for one output band n.
  ! bspec(ig,jg,iw) = -(1/pi) Im W^c_{G(ig) G(jg)}(q, omg(iw)) is the precomputed
  ! spectral function for this q (band-independent); domg is the (uniform) grid
  ! spacing.  Mirrors sigma_c_at_energy of gw_sigma_gpp with the single-pole
  ! bracket replaced by the omega'-grid integral.
  ! --------------------------------------------------------------------------
  subroutine sigma_c_real_at_energy(ng, no, nomega, omg, domg, bspec, &
                                    msig, eband, focc, eprobe, eta, ssig)
    implicit none
    integer,    intent(in)  :: ng, no, nomega
    real(8),    intent(in)  :: omg(nomega)        ! omega' grid (a.u., >= 0)
    real(8),    intent(in)  :: domg               ! grid spacing (a.u.)
    real(8),    intent(in)  :: bspec(ng, ng, nomega)
    complex(8), intent(in)  :: msig(ng, no)       ! M_{nn',G(ig)} for fixed n
    real(8),    intent(in)  :: eband(no)          ! eps_{n',k-q}
    real(8),    intent(in)  :: focc(no)           ! occupation in [0,1]
    real(8),    intent(in)  :: eprobe             ! probe energy E
    real(8),    intent(in)  :: eta                ! broadening (a.u.)
    complex(8), intent(out) :: ssig

    integer    :: ig, jg, inp, iw
    complex(8) :: zeta, mji, racc, brk
    real(8)    :: en, f1, fo, de, b

    zeta = cmplx(0.0d0, eta, 8)
    ssig = (0.0d0, 0.0d0)

    do inp = 1, no
      en = eband(inp)
      fo = focc(inp)
      if (fo < 0.0d0) fo = 0.0d0
      if (fo > 1.0d0) fo = 1.0d0
      f1 = 1.0d0 - fo
      de = eprobe - en
      racc = (0.0d0, 0.0d0)
      do jg = 1, ng
        mji = conjg(msig(jg,inp))               ! M_{nn',G'(jg)}
        if (mji == (0.0d0,0.0d0)) cycle
        do ig = 1, ng
          if (msig(ig,inp) == (0.0d0,0.0d0)) cycle
          ! spectral integral of the (ig,jg) channel:
          !   INT dw' B(w') [ (1-f)/(de - w' + i eta) + f/(de + w' - i eta) ]
          brk = (0.0d0, 0.0d0)
          do iw = 1, nomega
            b = bspec(ig,jg,iw)
            if (b == 0.0d0) cycle
            brk = brk + b * ( f1 / ( cmplx(de - omg(iw), 0.0d0, 8) + zeta ) &
                            + fo / ( cmplx(de + omg(iw), 0.0d0, 8) - zeta ) )
          end do
          racc = racc + msig(ig,inp) * (brk * domg) * mji
        end do
      end do
      ssig = ssig + racc
    end do

  end subroutine sigma_c_real_at_energy

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
    allocate(bspec(ng,ng,nomega))
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
          call sigma_c_real_at_energy(ng, no, nomega, omega_grid, domg, bspec, &
                                      msig, eband, focc, e0nk,         eta_au, s0)
          call sigma_c_real_at_energy(ng, no, nomega, omega_grid, domg, bspec, &
                                      msig, eband, focc, e0nk + de_au, eta_au, sp)
          call sigma_c_real_at_energy(ng, no, nomega, omega_grid, domg, bspec, &
                                      msig, eband, focc, e0nk - de_au, eta_au, sm)
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

    deallocate(mblk, msig, imap, chi0_w, epsinv_w, vcoul, bspec, eband, focc, e0)

  end subroutine calc_sigma_c_real

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
