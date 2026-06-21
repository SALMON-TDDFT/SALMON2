!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Plane-wave matrix elements for the GW module.
!
! Physics
! -------
! Bloch states are written  psi_{n,k}(r) = e^{i k.r} u_{n,k}(r),  where u is
! the cell-periodic part.  The plane-wave matrix element used by GW is
!
!   M_{n n'}(k,q,G) = < n, k+q | e^{i (q+G).r} | n', k >
!                   = (1/Omega) * integral_cell  u*_{n,k+q}(r) u_{n',k}(r) e^{i G.r} d3r
!
! i.e. the Fourier component (at reciprocal vector G) of the cell-periodic
! transition density  rho(r) = u*_{n,k+q}(r) u_{n',k}(r).
!
! Orbital convention (verified in the source)
! -------------------------------------------
! The real-space orbitals are stored in  spsi%zwf(ix,iy,iz,is,io,ik,im).  The
! kinetic operator (see the periodic stencil in the Hamiltonian) acts as
! (-i grad + k)^2 / 2 on zwf: it adds the constant 0.5*|k|^2 to the diagonal
! and the cross term k.(-i grad) through the k-dependent first-derivative
! coefficients, while the bare Laplacian acts on zwf.  The non-local
! projectors are stored pre-multiplied by e^{i k.r} (the "ekr" projectors) so
! that they pair with zwf.  Both facts show that zwf holds the cell-periodic
! part u_{n,k}, NOT the full Bloch function.  Therefore the transition density
! is exactly  rho(r) = conjg(zwf_{n,k+q}) * zwf_{n',k}  and the formula above
! applies directly with no extra e^{-i q.r} phase.
!
! Normalisation (pinned by the q=0 sanity check)
! ----------------------------------------------
! The orbitals are normalised so that  sum_grid |zwf|^2 * hvol = 1, with
! hvol = Omega / Ngrid (Omega = cell volume = system%det_a, Ngrid = number of
! real-space grid points).  Hence  sum_grid |u|^2 = Ngrid / Omega.
!
! The reused forward FFT (the periodic Hartree solver) returns
!   zrhoG_ele(G) = (1/Ngrid) * sum_grid rho(r) e^{-i G.r}.
! For the diagonal q=0, G=0 element, rho = |u|^2, so
!   zrhoG_ele(0) = (1/Ngrid) * (Ngrid/Omega) = 1/Omega.
! To obtain M_{nn}(k,0,0) = 1 we must therefore multiply by Omega:
!   M(ig) = zrhoG_ele(gmap(ig)) * Omega        (Omega = system%det_a).
! This is the physically correct factor: it turns the discretised
! (1/Ngrid) sum into the cell average that defines M.
!
! Sign of G: the reused FFT uses the e^{-i G.r} kernel and the reciprocal grid
! fg%vec_G is built with the matching convention, so reading zrhoG_ele at the
! grid point whose fg%vec_G equals gvec(:,ig) returns the Fourier coefficient
! that pairs with the stored G.  For G=0 (the validation) the sign is
! irrelevant; for G/=0 the convention is internally consistent provided Task 4
! uses the same gvec list for both M and the screened interaction v(q+G).
!
! Parallel layout
! ---------------
! The reciprocal-space density zrhoG_ele and fg%vec_G are distributed: each
! rank holds only the FFT grid points of its local mg slab.  build_gmap maps
! the global gvec list onto LOCAL grid points (returns iown=.false. for G not
! owned by this rank).  calc_mtxel fills the locally owned components and a
! comm_summation over icomm_rko assembles the full M across the r-space, k and
! orbital partitions.
!
! No proper nouns appear in this file.
!
module gw_mtxel_sub
  implicit none
  private

  public :: build_gmap
  public :: calc_mtxel

contains

  ! --------------------------------------------------------------------------
  ! build_gmap
  !
  ! One-time map from the global gvec list (Task 2) to local FFT grid points.
  !
  ! For each ig the routine searches the LOCAL part of fg%vec_G (the part this
  ! rank owns) for the grid point whose Cartesian G matches -gvec(:,ig).  The
  ! minus sign is deliberate: the reused forward FFT returns
  !   zrhoG_ele(vec_G) = (1/Ngrid) sum_r rho(r) e^{-i vec_G . r},
  ! whereas the GW matrix element is defined with the +i G.r kernel
  !   M(k,q,G) = (1/Omega) integral rho(r) e^{+i G.r} d3r.
  ! Reading zrhoG_ele at the grid point where vec_G = -G therefore yields the
  ! coefficient of e^{+i G.r}, i.e. M with the convention of the brief.  For
  ! G=0 the sign is immaterial (the validation), and because the gvec list is
  ! closed under negation the map is well defined for every ig.
  !
  ! The match is exact up to round-off because fg%vec_G is built from the same
  ! reciprocal primitive vectors as gvec; the tolerance is a small fraction of
  ! the shortest reciprocal primitive vector.
  !
  ! Arguments:
  !   fg          [in]  -- reciprocal grid (fg%vec_G distributed over local mg).
  !   mg          [in]  -- local real/recip grid descriptor (index ranges).
  !   ng          [in]  -- number of G vectors in the gvec list.
  !   gvec(3,ng)  [in]  -- Cartesian G vectors (bohr^-1), from build_gvectors.
  !   gmap(3,ng)  [out] -- FFT grid index (kx,ky,kz) of each ig if owned here,
  !                         else (0,0,0).
  !   iown(ng)    [out] -- .true. if this rank owns the grid point for ig.
  !
  ! An unmatched G that is NOT inside the local slab is normal (another rank
  ! owns it): iown(ig)=.false.  A G that should exist on the global FFT grid
  ! but is found nowhere by ANY rank indicates ng is not a subset of the FFT
  ! grid (cutoff too large for the chosen mesh); that global condition is not
  ! checked here but flagged by the caller's q=0 sanity (M_nn would not be 1).
  ! --------------------------------------------------------------------------
  subroutine build_gmap(fg, mg, ng, gvec, gmap, iown)
    use structures, only: s_reciprocal_grid, s_rgrid
    implicit none
    type(s_reciprocal_grid), intent(in)  :: fg
    type(s_rgrid),           intent(in)  :: mg
    integer,                 intent(in)  :: ng
    real(8),                 intent(in)  :: gvec(3,ng)
    integer,                 intent(out) :: gmap(3,ng)
    logical,                 intent(out) :: iown(ng)

    integer :: ig, kx, ky, kz
    real(8) :: gx, gy, gz, d2, tol2, gmin
    integer :: kxlo, kxhi, kylo, kyhi, kzlo, kzhi

    ! Local index window where BOTH fg%vec_G and zrhoG_ele are valid: the local
    ! mg slab in all three directions (zrhoG_ele is allocated over mg).
    kxlo = mg%is(1); kxhi = mg%ie(1)
    kylo = mg%is(2); kyhi = mg%ie(2)
    kzlo = mg%is(3); kzhi = mg%ie(3)

    ! Match tolerance: a small fraction of the shortest grid spacing in G.
    ! Use the smallest non-zero |G| present on this slab as the length scale;
    ! fall back to 1e-3 bohr^-1 if the slab is degenerate.
    gmin = huge(1.0d0)
    do kz = kzlo, kzhi
    do ky = kylo, kyhi
    do kx = kxlo, kxhi
      d2 = fg%vec_G(1,kx,ky,kz)**2 + fg%vec_G(2,kx,ky,kz)**2 + fg%vec_G(3,kx,ky,kz)**2
      if (d2 > 1.0d-12 .and. d2 < gmin) gmin = d2
    end do
    end do
    end do
    if (gmin >= huge(1.0d0)) then
      tol2 = 1.0d-6
    else
      tol2 = (1.0d-4 * sqrt(gmin))**2   ! (1e-4 * |G|_min)^2
    end if

    gmap = 0
    iown = .false.

    do ig = 1, ng
      ! Match -G so that mtxel(ig) carries the +i G.r convention (see header).
      gx = -gvec(1,ig); gy = -gvec(2,ig); gz = -gvec(3,ig)
      do kz = kzlo, kzhi
      do ky = kylo, kyhi
      do kx = kxlo, kxhi
        d2 = (fg%vec_G(1,kx,ky,kz)-gx)**2 &
           + (fg%vec_G(2,kx,ky,kz)-gy)**2 &
           + (fg%vec_G(3,kx,ky,kz)-gz)**2
        if (d2 <= tol2) then
          gmap(1,ig) = kx
          gmap(2,ig) = ky
          gmap(3,ig) = kz
          iown(ig)   = .true.
          go to 100
        end if
      end do
      end do
      end do
100   continue
    end do

  end subroutine build_gmap

  ! --------------------------------------------------------------------------
  ! calc_mtxel
  !
  ! Plane-wave matrix elements M_{n n'}(k,q,G) for a block of band pairs at a
  ! fixed (ik, ikq), where ikq is the stored index of k+q and ik the index of
  ! k.  Returns
  !   mtxel(ig, nb, nk') = M_{ (nb at ikq) , (nk' at ik) }(G_ig)
  ! for nb = 1..nb_bra (bra band index n) and nk' = 1..nb_ket (ket band n').
  !
  ! For each band pair the routine forms the local-slab transition density
  !   rho(r) = conjg(zwf_{nb,ikq}(r)) * zwf_{nk',ik}(r),
  ! Fourier transforms it via the periodic Hartree forward FFT (real and
  ! imaginary parts separately, by linearity), multiplies by Omega = det_a to
  ! pin M_nn=1, picks the locally owned G components, and a single
  ! comm_summation over icomm_rko assembles the global mtxel.
  !
  ! Spin: the spin index is fixed by the argument ispin (1 for the
  ! non-magnetic test); kept explicit so a spin-polarised caller can loop.
  !
  ! Distribution note (q=0 validation / Task 4):
  ! The product needs BOTH orbitals (at ikq and ik) on the same r-space slab.
  ! This is satisfied when each k-point's bands are locally available on the
  ! ranks owning that k (e.g. nporbital=1, the small test layout).  Bands not
  ! in [io_s,io_e] or k not in [ik_s,ik_e] contribute zero on this rank and are
  ! supplied by the owning rank through the final comm_summation.  A general
  ! cross-rank orbital gather (when k+q lives on a different orbital/k rank than
  ! k) is an orthogonal parallelisation issue deferred to Task 4; the q=0
  ! validation uses ikq=ik so no gather is needed.
  ! --------------------------------------------------------------------------
  subroutine calc_mtxel(system, info, mg, lg, fg, poisson, spsi, &
                        gmap, iown, ng, ik, ikq, ispin, nb_bra, nb_ket, mtxel)
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, &
                          s_reciprocal_grid, s_poisson, s_orbital, s_scalar, &
                          allocate_scalar
    use communication, only: comm_summation
    use poisson_periodic, only: poisson_ft, poisson_ffte
#ifdef USE_FFTW
    use poisson_periodic, only: poisson_fftw
#endif
    use salmon_global, only: yn_ffte, yn_fftw
    implicit none
    type(s_dft_system),     intent(in) :: system
    type(s_parallel_info),  intent(in) :: info
    type(s_rgrid),          intent(in) :: mg
    type(s_rgrid),          intent(in) :: lg
    type(s_reciprocal_grid),intent(in) :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_orbital),        intent(in) :: spsi
    integer,                intent(in) :: gmap(3,ng)
    logical,                intent(in) :: iown(ng)
    integer,                intent(in) :: ng
    integer,                intent(in) :: ik, ikq, ispin
    integer,                intent(in) :: nb_bra, nb_ket
    complex(8),             intent(out):: mtxel(ng, nb_bra, nb_ket)

    type(s_scalar) :: rho_re, rho_im, vh_scr
    complex(8), allocatable :: mloc(:,:,:)
    complex(8) :: zd
    real(8) :: omega
    integer :: nb, nk2, ig, ix, iy, iz
    integer :: im
    logical :: have_k, have_kq
    logical :: bra_owned, ket_owned

    omega = system%det_a
    im    = info%im_s   ! single Maxwell replica for the GW use case

    ! Scratch real fields over the local mg slab (match poisson_* expectations).
    call allocate_scalar(mg, rho_re)
    call allocate_scalar(mg, rho_im)
    call allocate_scalar(mg, vh_scr)

    allocate(mloc(ng, nb_bra, nb_ket))
    mloc = (0.0d0, 0.0d0)

    ! Does this rank hold the k and k+q orbital blocks?
    have_k  = (ik  >= info%ik_s .and. ik  <= info%ik_e)
    have_kq = (ikq >= info%ik_s .and. ikq <= info%ik_e)

    if (have_k .and. have_kq) then
      do nk2 = 1, nb_ket
        ket_owned = (nk2 >= info%io_s .and. nk2 <= info%io_e)
        if (.not. ket_owned) cycle
        do nb = 1, nb_bra
          bra_owned = (nb >= info%io_s .and. nb <= info%io_e)
          if (.not. bra_owned) cycle

          ! Transition density on the local slab, split into Re and Im parts.
!$omp parallel do collapse(2) private(iz,iy,ix,zd)
          do iz = mg%is(3), mg%ie(3)
          do iy = mg%is(2), mg%ie(2)
          do ix = mg%is(1), mg%ie(1)
            zd = conjg(spsi%zwf(ix,iy,iz,ispin,nb ,ikq,im)) &
               *       spsi%zwf(ix,iy,iz,ispin,nk2,ik ,im)
            rho_re%f(ix,iy,iz) = dble(zd)
            rho_im%f(ix,iy,iz) = aimag(zd)
          end do
          end do
          end do

          ! Forward FFT of the real part -> zrhoG_ele; capture owned G.
          call forward_fft(rho_re, vh_scr)
          do ig = 1, ng
            if (iown(ig)) &
              mloc(ig,nb,nk2) = mloc(ig,nb,nk2) &
                + poisson%zrhoG_ele(gmap(1,ig),gmap(2,ig),gmap(3,ig)) * omega
          end do

          ! Forward FFT of the imaginary part -> add i * result.
          call forward_fft(rho_im, vh_scr)
          do ig = 1, ng
            if (iown(ig)) &
              mloc(ig,nb,nk2) = mloc(ig,nb,nk2) &
                + (0.0d0,1.0d0) &
                * poisson%zrhoG_ele(gmap(1,ig),gmap(2,ig),gmap(3,ig)) * omega
          end do

        end do
      end do
    end if

    ! Assemble the full M over the r-space, k and orbital partitions.
    call comm_summation(mloc, mtxel, ng*nb_bra*nb_ket, info%icomm_rko)

    deallocate(mloc)
    deallocate(rho_re%f, rho_im%f, vh_scr%f)

  contains

    ! Dispatch the periodic forward FFT for one real field, leaving the
    ! G-space result in poisson%zrhoG_ele.  Mirrors the backend selection of
    ! the Hartree driver; the returned Vh is unused scratch.
    subroutine forward_fft(rho_in, vh_out)
      type(s_scalar), intent(in)    :: rho_in
      type(s_scalar), intent(inout) :: vh_out
#ifdef USE_FFTW
      if (yn_fftw == 'y') then
        call poisson_fftw(lg, mg, info, fg, rho_in, vh_out, poisson)
        return
      end if
#endif
      if (yn_ffte == 'y') then
        call poisson_ffte(lg, mg, info, fg, rho_in, vh_out, poisson)
      else
        call poisson_ft(lg, mg, info, fg, rho_in, vh_out, poisson)
      end if
    end subroutine forward_fft

  end subroutine calc_mtxel

end module gw_mtxel_sub
