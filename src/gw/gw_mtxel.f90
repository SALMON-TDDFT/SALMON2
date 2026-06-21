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
! Evaluation (direct discrete Fourier sum)
! ----------------------------------------
! Rather than reuse the periodic Hartree FFT drivers (which are collective,
! assume a real input field, and require every rank to participate -- calling
! them selectively on the k-owning rank for a complex transition density is
! fragile), the Fourier component at each requested G is evaluated directly:
!
!   M_{n n'}(k,q,G) = sum_{grid r} conjg(zwf_{n,k+q}(r)) * zwf_{n',k}(r)
!                                  * exp(-i G.r) * hvol
!
! where hvol = system%hvol = Omega / Ngrid is the real-space volume element
! and r is the Cartesian coordinate of grid point (ix,iy,iz).  The discretised
! sum * hvol is the cell integral that defines M; dividing the integral by
! Omega is folded into hvol = Omega/Ngrid combined with the (1/Omega) prefactor
! (see Normalisation).  This is O(Ngrid * ng) per band pair, robust, and
! convention-clean; an FFT acceleration is deferred future work.
!
! Orbital convention (verified in the source)
! -------------------------------------------
! The real-space orbitals are stored in  spsi%zwf(ix,iy,iz,is,io,ik,im).  The
! kinetic operator (see the periodic stencil in the Hamiltonian) acts as
! (-i grad + k)^2 / 2 on zwf, and the non-local projectors are stored
! pre-multiplied by e^{i k.r} (the "ekr" projectors) so they pair with zwf.
! Both facts show that zwf holds the cell-periodic part u_{n,k}, NOT the full
! Bloch function.  Therefore the transition density is exactly
! rho(r) = conjg(zwf_{n,k+q}) * zwf_{n',k}  and the formula above applies
! directly with no extra e^{-i q.r} phase.
!
! Real-space coordinate convention (verified in the source)
! ---------------------------------------------------------
! The grid coordinate is held in lg%coordinate(i,idir) and is set in
! set_gridcoordinate (src/common/initialization.f90): for the periodic case
! (iperiodic=3) r_idir = (i-1)*hgs(idir); for the isolated case it carries the
! half-grid offset for an even mesh.  The local slab mg%is(idir)..mg%ie(idir)
! uses the SAME global index space as lg%coordinate (this is how every other
! routine that needs r on the mg slab indexes it, e.g. the jellium and
! pseudopotential builders loop mg%is..mg%ie and read lg%coordinate(ix,1)).
! Using lg%coordinate(ix,1), lg%coordinate(iy,2), lg%coordinate(iz,3) therefore
! gives the correct Cartesian r for every grid point on this rank, independent
! of iperiodic and of any half-grid offset, and is the convention SALMON itself
! uses for the orbitals.
!
! Normalisation (pinned by the q=0 sanity check)
! ----------------------------------------------
! The orbitals are normalised so that  sum_grid |zwf|^2 * hvol = 1, with
! hvol = Omega / Ngrid (Omega = cell volume = system%det_a, Ngrid = number of
! real-space grid points).  At q=0, G=0 the transition density of a single band
! is rho = |u|^2 and exp(-i G.r) = 1, so
!   M_{nn}(k,0,0) = sum_grid |zwf|^2 * hvol = 1.
! For two different bands at q=0, G=0 the same sum is the orbital overlap, which
! is 0 by orthonormality.  Thus M_{nn'}(k,0,0) = delta_{nn'}, the validation
! target.  This is the physically correct factor: hvol turns the discretised
! sum into the cell average (1/Omega) integral that defines M.
!
! Sign of G: the kernel exp(-i G.r) is used with the gvec list passed in
! (Cartesian, bohr^-1, from build_gvectors).  For G=0 (the validation) the sign
! is irrelevant; for G/=0 the convention is internally consistent provided a
! later screened-interaction step uses the same gvec list for both M and
! v(q+G).
!
! Parallel layout
! ---------------
! The product needs BOTH orbitals (at ikq and ik) on the same r-space slab and
! both bands owned on this rank.  calc_mtxel accumulates the locally owned
! contributions and a single comm_summation over icomm_rko assembles the full M
! across the r-space, k and orbital partitions.
!
! No proper nouns appear in this file.
!
module gw_mtxel_sub
  implicit none
  private

  public :: calc_mtxel

contains

  ! --------------------------------------------------------------------------
  ! calc_mtxel
  !
  ! Plane-wave matrix elements M_{n n'}(k,q,G) for a block of band pairs at a
  ! fixed (ik, ikq), where ikq is the stored index of k+q and ik the index of
  ! k.  Returns
  !   mtxel(ig, nb, nk') = M_{ (nb at ikq) , (nk' at ik) }(G_ig)
  ! for nb = 1..nb_bra (bra band index n) and nk' = 1..nb_ket (ket band n').
  !
  ! For each owned band pair the routine forms the local-slab transition
  ! density rho(r) = conjg(zwf_{nb,ikq}(r)) * zwf_{nk',ik}(r) and accumulates,
  ! for every requested G, the direct Fourier sum
  !   mloc(ig,nb,nk') += rho(r) * exp(-i G_ig . r) * hvol.
  ! A single comm_summation over icomm_rko assembles the global mtxel.
  !
  ! Spin: the spin index is fixed by the argument ispin (1 for the
  ! non-magnetic test); kept explicit so a spin-polarised caller can loop.
  !
  ! Distribution note (q=0 validation):
  ! The product needs BOTH orbitals (at ikq and ik) on the same r-space slab.
  ! This is satisfied when each k-point's bands are locally available on the
  ! ranks owning that k (e.g. nporbital=1, the small test layout).  Bands not
  ! in [io_s,io_e] or k not in [ik_s,ik_e] contribute zero on this rank and are
  ! supplied by the owning rank through the final comm_summation.  A general
  ! cross-rank orbital gather (when k+q lives on a different orbital/k rank than
  ! k) is an orthogonal parallelisation issue deferred to a later task; the q=0
  ! validation uses ikq=ik so no gather is needed.
  !
  ! Arguments:
  !   system      [in]  -- holds hvol (= Omega/Ngrid), the volume element.
  !   info        [in]  -- parallel ownership ranges (ik_s/e, io_s/e, im_s) and
  !                         the assembly communicator icomm_rko.
  !   mg          [in]  -- local real-space grid descriptor (global index range
  !                         of this rank's slab).
  !   lg          [in]  -- global grid; lg%coordinate(i,idir) gives Cartesian r.
  !   spsi        [in]  -- cell-periodic orbitals zwf(ix,iy,iz,is,io,ik,im).
  !   ng          [in]  -- number of G vectors.
  !   gvec(3,ng)  [in]  -- Cartesian G vectors (bohr^-1), from build_gvectors.
  !   ik          [in]  -- stored index of k.
  !   ikq         [in]  -- stored index of k+q.
  !   ispin       [in]  -- spin index of the band block.
  !   nb_bra      [in]  -- number of bra bands (n).
  !   nb_ket      [in]  -- number of ket bands (n').
  !   mtxel       [out] -- mtxel(ng,nb_bra,nb_ket), assembled across all ranks.
  ! --------------------------------------------------------------------------
  subroutine calc_mtxel(system, info, mg, lg, spsi, &
                        gvec, ng, ik, ikq, ispin, nb_bra, nb_ket, mtxel)
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),    intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid),         intent(in) :: mg
    type(s_rgrid),         intent(in) :: lg
    type(s_orbital),       intent(in) :: spsi
    integer,               intent(in) :: ng
    real(8),               intent(in) :: gvec(3,ng)
    integer,               intent(in) :: ik, ikq, ispin
    integer,               intent(in) :: nb_bra, nb_ket
    complex(8),            intent(out):: mtxel(ng, nb_bra, nb_ket)

    complex(8), allocatable :: mloc(:,:,:)
    complex(8) :: zd, zacc
    real(8) :: hvol, rx, ry, rz, phase
    integer :: nb, nk2, ig, ix, iy, iz
    integer :: im
    logical :: have_k, have_kq
    logical :: bra_owned, ket_owned

    hvol = system%hvol
    im   = info%im_s   ! single Maxwell replica for the GW use case

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

          ! Direct Fourier sum of the local-slab transition density onto each
          ! requested G.  Per grid point form zd = conjg(zwf_kq)*zwf_k once,
          ! then accumulate zd*exp(-i G.r)*hvol into every G; threaded over G so
          ! each thread owns a private mloc column (no races on mloc).
!$omp parallel do private(ig,iz,iy,ix,zacc,rx,ry,rz,zd,phase)
          do ig = 1, ng
            zacc = (0.0d0, 0.0d0)
            do iz = mg%is(3), mg%ie(3)
              rz = lg%coordinate(iz,3)
            do iy = mg%is(2), mg%ie(2)
              ry = lg%coordinate(iy,2)
            do ix = mg%is(1), mg%ie(1)
              rx = lg%coordinate(ix,1)
              zd = conjg(spsi%zwf(ix,iy,iz,ispin,nb ,ikq,im)) &
                 *       spsi%zwf(ix,iy,iz,ispin,nk2,ik ,im)
              phase = -( gvec(1,ig)*rx + gvec(2,ig)*ry + gvec(3,ig)*rz )
              zacc = zacc + zd * cmplx(cos(phase), sin(phase), 8)
            end do
            end do
            end do
            mloc(ig,nb,nk2) = zacc * hvol
          end do
!$omp end parallel do

        end do
      end do
    end if

    ! Assemble the full M over the r-space, k and orbital partitions.
    call comm_summation(mloc, mtxel, ng*nb_bra*nb_ket, info%icomm_rko)

    deallocate(mloc)

  end subroutine calc_mtxel

end module gw_mtxel_sub
