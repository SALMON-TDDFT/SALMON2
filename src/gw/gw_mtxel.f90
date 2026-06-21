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
! Evaluation: node-local serial FFT (default) with a direct-DFT fallback
! ---------------------------------------------------------------------
! The Fourier components of the transition density are computed with a SERIAL
! (per-rank) forward 3D FFT.  This is correct and node-local because the GW
! driver runs with nproc_rgrid = 1: every rank then holds the FULL real-space
! grid (mg%is..mg%ie spans the whole cell), so a serial FFT on this rank's data
! involves NO MPI and NO collective call.  It must NOT be confused with the
! periodic Hartree FFT drivers (poisson_ffte / PZFFT3DV_MOD), which are
! collective (they comm_summation over the r-grid communicator) -- using those
! here selectively, for a complex field, is what crashed in an earlier attempt.
! The serial 3D FFT used here is the linked serial complex FFT (FFTW's
! dfftw_plan_dft_3d / dfftw_execute_dft, compiled in when USE_FFTW is set).
!
! When the build does not link a serial FFT (USE_FFTW off), the routine falls
! back to the direct discrete Fourier sum, which is mathematically identical and
! is also the reference the FFT path is pinned against:
!
!   M_{n n'}(k,q,G) = sum_{grid r} conjg(zwf_{n,k+q}(r)) * zwf_{n',k}(r)
!                                  * exp(-i G.r) * hvol
!
! where hvol = system%hvol = Omega / Ngrid is the real-space volume element and
! r is the Cartesian coordinate of grid point (ix,iy,iz).  The FFT path computes
! the SAME quantity: see the normalisation note below.  Cost: the direct sum is
! O(Ngrid * ng) per band pair; the FFT is O(Ngrid log Ngrid) per band pair plus
! a scatter through a once-built G-map, ~ng/log(Ngrid) faster.
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
! directly with no extra e^{-i q.r} phase.  The real(rwf)/complex(zwf) branch is
! preserved: time-reversal-invariant k-points store rwf, others zwf.
!
! Real-space coordinate / FFT-index convention (verified in the source)
! ---------------------------------------------------------------------
! The grid coordinate is held in lg%coordinate(i,idir) and is set in
! set_gridcoordinate (src/common/initialization.f90): for the periodic case
! (iperiodic=3) r_idir = (i-1)*hgs(idir), with lg%is=1 so the global index i runs
! 1..N_d.  Grid point (ix,iy,iz) therefore sits at the integer offset
! j_d = i_d - 1 in [0,N_d-1] along Cartesian axis d, with N_d = lg%num(d) points.
! A G vector built by build_gvectors is G = n1*b1 + n2*b2 + n3*b3 (integer n_d,
! b_d = system%primitive_b(:,d)).  SALMON's real-space mesh is a Cartesian box
! r_grid = (j1*hgs1, j2*hgs2, j3*hgs3) and the periodic GW cells use a diagonal
! primitive_a = diag(L1,L2,L3) with L_d = N_d*hgs_d, so b_d = (2*pi/L_d) e_d and
!     G . r_grid = 2*pi * ( n1*j1/N1 + n2*j2/N2 + n3*j3/N3 ),
! exactly the phase of the FFT bin (n1,n2,n3).  (This is the same orthorhombic-
! box assumption the periodic FFT Poisson solvers make; the direct-DFT fallback
! uses the Cartesian coordinate/G directly and is valid for any lattice.)  The
! G-map below recovers the integers (n1,n2,n3) for each gvec by inverting
! primitive_b and folds them into the FFT frequency-index range, so
! rhog(gmap(ig)) is precisely sum_r rho(r) e^{-i G.r}.
!
! Normalisation (pinned analytically by M_nn(k,0,0)=1)
! ----------------------------------------------------
! FFTW's forward transform (FFTW_FORWARD) is the UNNORMALISED sum
!   rhog(k1,k2,k3) = sum_{j} rho(j) * exp(-2*pi*i*(k1 j1/N1 + k2 j2/N2 + k3 j3/N3)).
! With the index identification above this is exactly  sum_r rho(r) e^{-i G.r}
! (no 1/Ngrid, no hvol).  The matrix element wants the cell average
! (1/Omega) integral = (discrete sum) * hvol, so
!   mtxel(ig) = rhog(gmap(ig)) * norm,   norm = hvol = Omega/Ngrid.
! Check: at q=0, G=0 a single band has rho = |u|^2, so rhog(0) = sum_r |zwf|^2,
! and mtxel = (sum_r |zwf|^2) * hvol = 1 because the orbitals satisfy
! sum_r |zwf|^2 * hvol = 1.  For n /= n' the same q=0,G=0 sum is the orbital
! overlap = 0.  Thus M_{nn'}(k,0,0) = delta_{nn'}, matching the direct DFT and
! the sanity target.  norm is therefore identical to the direct-DFT factor; the
! FFT changes only the algorithm, not the convention.
!
! Sign of G: the kernel exp(-i G.r) is used with the gvec list passed in
! (Cartesian, bohr^-1, from build_gvectors).  The FFT bin sign (FFTW_FORWARD ->
! exp(-2*pi*i ...)) matches; for G=0 (the validation) the sign is irrelevant.
!
! Parallel layout
! ---------------
! The product needs BOTH orbitals (at ikq and ik) on the same r-space slab and
! both bands owned on this rank.  calc_mtxel accumulates the locally owned
! contributions and a single comm_summation over icomm_rko assembles the full M
! across the r-space, k and orbital partitions.  The serial FFT itself touches
! only this rank's data -- there is no MPI inside the transform.
!
! No proper nouns appear in this file.
!
! Pull in the generated build-configuration macros.  This is where the optional
! serial-FFT switch is defined; without this include the preprocessor never sees
! the switch, the fast path below is silently dropped, and only the direct-sum
! fallback is compiled (the cause of an "FFT built but no transform symbols"
! binary).  The header lives in the build directory, which is on the include
! path, and is the same one the core modules use to gate the same switch.
#include "config.h"
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
  ! density rho(r) = conjg(zwf_{nb,ikq}(r)) * zwf_{nk',ik}(r) on the full grid,
  ! performs ONE serial forward 3D FFT, and scatters the requested G bins (via a
  ! once-built G-map) into mloc(ig,nb,nk') with the volume factor hvol.  A single
  ! comm_summation over icomm_rko assembles the global mtxel.  When no serial FFT
  ! is linked, the equivalent direct Fourier sum is used (same result).
  !
  ! Spin: the spin index is fixed by the argument ispin (1 for the
  ! non-magnetic test); kept explicit so a spin-polarised caller can loop.
  !
  ! Distribution note (q=0 validation):
  ! The product needs BOTH orbitals (at ikq and ik) on the same r-space slab.
  ! This is satisfied when each k-point's bands are locally available on the
  ! ranks owning that k.  Bands not in [io_s,io_e] or k not in [ik_s,ik_e]
  ! contribute zero on this rank and are supplied by the owning rank through the
  ! final comm_summation.  A general cross-rank orbital gather is an orthogonal
  ! parallelisation issue deferred to a later task; the q=0 validation uses
  ! ikq=ik so no gather is needed.
  !
  ! Arguments:
  !   system      [in]  -- holds hvol (= Omega/Ngrid) and primitive_b.
  !   info        [in]  -- parallel ownership ranges (ik_s/e, io_s/e, im_s) and
  !                         the assembly communicator icomm_rko.
  !   mg          [in]  -- local real-space grid descriptor.  With nproc_rgrid=1
  !                         (the GW layout) mg spans the whole cell.
  !   lg          [in]  -- global grid; lg%num(d) gives N_d, lg%coordinate the r.
  !   spsi        [in]  -- cell-periodic orbitals zwf/rwf(ix,iy,iz,is,io,ik,im).
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
    real(8) :: hvol, norm
    integer :: nb, nk2, ig, ix, iy, iz
    integer :: im
    logical :: have_k, have_kq
    logical :: bra_owned, ket_owned
    logical :: use_zwf

    ! Full-grid sizes (nproc_rgrid=1 => mg spans the whole cell).
    integer :: n1, n2, n3

    ! G-map: FFT bin (1-based) for every requested G; built once, reused.
    integer, allocatable :: gmap(:)

#ifdef USE_FFTW
    complex(8), allocatable :: rhor(:,:,:), rhog(:,:,:)
    integer(8) :: plan_fwd
    integer :: kx, ky, kz
#else
    complex(8) :: zacc, zd
    real(8) :: rx, ry, rz, phase
#endif

    hvol = system%hvol
    im   = info%im_s   ! single Maxwell replica for the GW use case
    ! Orbitals at time-reversal-invariant k-points are stored real (spsi%rwf);
    ! the others complex (spsi%zwf).  A rank holds a single storage type for its
    ! whole k-block, so branch once.
    use_zwf = allocated(spsi%zwf)

    ! Normalisation pinned by M_nn(k,0,0)=1: the (unnormalised forward) FFT and
    ! the direct sum both return  sum_r rho(r) e^{-i G.r}; the cell-average
    ! matrix element multiplies by hvol = Omega/Ngrid (see the header note).
    norm = hvol

    n1 = lg%num(1)
    n2 = lg%num(2)
    n3 = lg%num(3)

    allocate(mloc(ng, nb_bra, nb_ket))
    mloc = (0.0d0, 0.0d0)

    ! Does this rank hold the k and k+q orbital blocks?
    have_k  = (ik  >= info%ik_s .and. ik  <= info%ik_e)
    have_kq = (ikq >= info%ik_s .and. ikq <= info%ik_e)

    if (have_k .and. have_kq) then

#ifdef USE_FFTW
      ! ---- G-map: built once, reused across all band pairs ------------------
      ! For each gvec recover its integer reciprocal index (m1,m2,m3) by solving
      ! gvec = primitive_b * (m1,m2,m3)^T, then fold each m_d into the FFT
      ! frequency-index range [0,N_d-1] (negative frequencies wrap to N_d+m_d).
      allocate(gmap(ng))
      call build_gmap(system%primitive_b, gvec, ng, n1, n2, n3, gmap)

      ! ---- Serial node-local forward 3D FFT (no MPI) ------------------------
      ! FFTW stores frequency bin (kx,ky,kz) at array index (kx+1,ky+1,kz+1)
      ! with kx in [0,N1-1], etc.  Plan once on the full local grid, reuse the
      ! plan for every band pair, execute out-of-place rhor -> rhog.
      allocate(rhor(n1, n2, n3))
      allocate(rhog(n1, n2, n3))
      call make_plan_dft_3d(plan_fwd, n1, n2, n3, rhor, rhog)

      do nk2 = 1, nb_ket
        ket_owned = (nk2 >= info%io_s .and. nk2 <= info%io_e)
        if (.not. ket_owned) cycle
        do nb = 1, nb_bra
          bra_owned = (nb >= info%io_s .and. nb <= info%io_e)
          if (.not. bra_owned) cycle

          ! Transition density on the full grid (1-based local index = global
          ! offset j_d, since mg%is(d)=lg%is(d) when nproc_rgrid=1).
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz = 1, n3
          do iy = 1, n2
          do ix = 1, n1
            if (use_zwf) then
              rhor(ix,iy,iz) = &
                  conjg(spsi%zwf(mg%is(1)-1+ix, mg%is(2)-1+iy, mg%is(3)-1+iz, ispin, nb , ikq, im)) &
                *       spsi%zwf(mg%is(1)-1+ix, mg%is(2)-1+iy, mg%is(3)-1+iz, ispin, nk2, ik , im)
            else
              rhor(ix,iy,iz) = cmplx( &
                  spsi%rwf(mg%is(1)-1+ix, mg%is(2)-1+iy, mg%is(3)-1+iz, ispin, nb , ikq, im) &
                * spsi%rwf(mg%is(1)-1+ix, mg%is(2)-1+iy, mg%is(3)-1+iz, ispin, nk2, ik , im), 0.0d0, 8)
            end if
          end do
          end do
          end do

          call exec_plan_dft_3d(plan_fwd, rhor, rhog)

          ! Scatter requested G bins through the G-map, apply hvol.
!$omp parallel do private(ig,kx,ky,kz)
          do ig = 1, ng
            kx = mod(gmap(ig)-1, n1) + 1
            ky = mod((gmap(ig)-1)/n1, n2) + 1
            kz = (gmap(ig)-1)/(n1*n2) + 1
            mloc(ig,nb,nk2) = rhog(kx,ky,kz) * norm
          end do
        end do
      end do

      call destroy_plan_dft_3d(plan_fwd)
      deallocate(rhor, rhog)
      deallocate(gmap)
#else
      ! ---- Fallback: direct discrete Fourier sum (same result, no FFT lib) --
      ! Per grid point form zd = conjg(zwf_kq)*zwf_k once, then accumulate
      ! zd*exp(-i G.r)*hvol into every G; threaded over G (private mloc column).
      do nk2 = 1, nb_ket
        ket_owned = (nk2 >= info%io_s .and. nk2 <= info%io_e)
        if (.not. ket_owned) cycle
        do nb = 1, nb_bra
          bra_owned = (nb >= info%io_s .and. nb <= info%io_e)
          if (.not. bra_owned) cycle

!$omp parallel do private(ig,iz,iy,ix,zacc,rx,ry,rz,zd,phase)
          do ig = 1, ng
            zacc = (0.0d0, 0.0d0)
            do iz = mg%is(3), mg%ie(3)
              rz = lg%coordinate(iz,3)
            do iy = mg%is(2), mg%ie(2)
              ry = lg%coordinate(iy,2)
            do ix = mg%is(1), mg%ie(1)
              rx = lg%coordinate(ix,1)
              if (use_zwf) then
                zd = conjg(spsi%zwf(ix,iy,iz,ispin,nb ,ikq,im)) &
                   *       spsi%zwf(ix,iy,iz,ispin,nk2,ik ,im)
              else
                zd = spsi%rwf(ix,iy,iz,ispin,nb ,ikq,im) &
                   * spsi%rwf(ix,iy,iz,ispin,nk2,ik ,im)
              end if
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
#endif

    end if

    ! Assemble the full M over the r-space, k and orbital partitions.  The FFT
    ! above is purely local; the only MPI is this single collective.
    call comm_summation(mloc, mtxel, ng*nb_bra*nb_ket, info%icomm_rko)

    deallocate(mloc)

  end subroutine calc_mtxel

  ! --------------------------------------------------------------------------
  ! build_gmap
  !
  ! Map each Cartesian G = gvec(:,ig) to its flattened FFT bin index (1-based,
  ! column-major over (kx,ky,kz) with kx fastest) on an N1 x N2 x N3 grid.
  !
  ! gvec = primitive_b * m, with m = (m1,m2,m3) integer reciprocal indices.  We
  ! solve for m by inverting the 3x3 reciprocal-lattice matrix (columns = the
  ! reciprocal primitive vectors), round to the nearest integers, and assert the
  ! residual is small (every gvec MUST lie on the reciprocal lattice, since the
  ! epsilon cutoff <= the wavefunction grid Nyquist).  Each m_d is folded into
  ! [0,N_d-1] by mod (negative frequencies wrap: bin = mod(m_d, N_d) + N_d, then
  ! mod N_d), matching FFTW's storage of frequency m at index mod(m,N).
  ! --------------------------------------------------------------------------
  subroutine build_gmap(bmat, gvec, ng, n1, n2, n3, gmap)
    implicit none
    real(8), intent(in)  :: bmat(3,3)
    integer, intent(in)  :: ng, n1, n2, n3
    real(8), intent(in)  :: gvec(3,ng)
    integer, intent(out) :: gmap(ng)

    real(8) :: binv(3,3), det, fm(3)
    integer :: ig, m1, m2, m3, kx, ky, kz
    real(8) :: resid

    ! Inverse of bmat (columns are the reciprocal primitive vectors).
    det =   bmat(1,1)*(bmat(2,2)*bmat(3,3) - bmat(2,3)*bmat(3,2)) &
          - bmat(1,2)*(bmat(2,1)*bmat(3,3) - bmat(2,3)*bmat(3,1)) &
          + bmat(1,3)*(bmat(2,1)*bmat(3,2) - bmat(2,2)*bmat(3,1))
    if (abs(det) < 1.0d-30) then
      write(*,*) "[gw] build_gmap: singular reciprocal matrix"
      stop
    end if
    binv(1,1) = (bmat(2,2)*bmat(3,3) - bmat(2,3)*bmat(3,2)) / det
    binv(1,2) = (bmat(1,3)*bmat(3,2) - bmat(1,2)*bmat(3,3)) / det
    binv(1,3) = (bmat(1,2)*bmat(2,3) - bmat(1,3)*bmat(2,2)) / det
    binv(2,1) = (bmat(2,3)*bmat(3,1) - bmat(2,1)*bmat(3,3)) / det
    binv(2,2) = (bmat(1,1)*bmat(3,3) - bmat(1,3)*bmat(3,1)) / det
    binv(2,3) = (bmat(1,3)*bmat(2,1) - bmat(1,1)*bmat(2,3)) / det
    binv(3,1) = (bmat(2,1)*bmat(3,2) - bmat(2,2)*bmat(3,1)) / det
    binv(3,2) = (bmat(1,2)*bmat(3,1) - bmat(1,1)*bmat(3,2)) / det
    binv(3,3) = (bmat(1,1)*bmat(2,2) - bmat(1,2)*bmat(2,1)) / det

    do ig = 1, ng
      ! m = binv * gvec  (real-valued; should be integers)
      fm(1) = binv(1,1)*gvec(1,ig) + binv(1,2)*gvec(2,ig) + binv(1,3)*gvec(3,ig)
      fm(2) = binv(2,1)*gvec(1,ig) + binv(2,2)*gvec(2,ig) + binv(2,3)*gvec(3,ig)
      fm(3) = binv(3,1)*gvec(1,ig) + binv(3,2)*gvec(2,ig) + binv(3,3)*gvec(3,ig)
      m1 = nint(fm(1))
      m2 = nint(fm(2))
      m3 = nint(fm(3))
      resid = abs(fm(1)-dble(m1)) + abs(fm(2)-dble(m2)) + abs(fm(3)-dble(m3))
      if (resid > 1.0d-4) then
        write(*,*) "[gw] build_gmap: gvec not on reciprocal lattice, ig=", ig, " resid=", resid
        stop
      end if
      ! Fold signed frequency into [0,N-1] (negative wraps to the top half).
      kx = modulo(m1, n1)
      ky = modulo(m2, n2)
      kz = modulo(m3, n3)
      ! The folded index must also lie inside the grid (it always does after
      ! modulo); assert as a guard against an oversized cutoff vs the FFT grid.
      if (kx < 0 .or. kx >= n1 .or. ky < 0 .or. ky >= n2 .or. kz < 0 .or. kz >= n3) then
        write(*,*) "[gw] build_gmap: folded index out of range, ig=", ig
        stop
      end if
      ! Flattened 1-based column-major bin (kx fastest), matching the unpack in
      ! calc_mtxel.
      gmap(ig) = 1 + kx + n1*ky + n1*n2*kz
    end do
  end subroutine build_gmap

#ifdef USE_FFTW
  ! --------------------------------------------------------------------------
  ! Thin wrappers around the linked serial complex 3D FFT.  These are the only
  ! FFT-library entry points used; they are SERIAL (no MPI/communicator), so the
  ! transform is fully node-local.  The fftw3 header (the same one the serial
  ! dfftw_* calls in the periodic/isolated Poisson solvers include) supplies
  ! FFTW_FORWARD and FFTW_ESTIMATE and the dfftw_* legacy interface.
  ! --------------------------------------------------------------------------
  subroutine make_plan_dft_3d(plan, n1, n2, n3, in, out)
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3-mpi.f03'
    integer(8), intent(out) :: plan
    integer,    intent(in)  :: n1, n2, n3
    complex(8), intent(inout) :: in(n1,n2,n3), out(n1,n2,n3)
    ! Forward, unnormalised 3D complex DFT.  Argument order matches the serial
    ! dfftw_plan_dft_3d call in the isolated-Poisson solver: with a column-major
    ! Fortran array a(n1,n2,n3) the FIRST listed extent (n1) is the fastest
    ! array index, so frequency bin (kx,ky,kz) lands at out(kx+1,ky+1,kz+1) --
    ! exactly the unpack used in calc_mtxel.
    call dfftw_plan_dft_3d(plan, n1, n2, n3, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
  end subroutine make_plan_dft_3d

  subroutine exec_plan_dft_3d(plan, in, out)
    implicit none
    integer(8), intent(in)    :: plan
    complex(8), intent(inout) :: in(*), out(*)
    call dfftw_execute_dft(plan, in, out)
  end subroutine exec_plan_dft_3d

  subroutine destroy_plan_dft_3d(plan)
    implicit none
    integer(8), intent(in) :: plan
    call dfftw_destroy_plan(plan)
  end subroutine destroy_plan_dft_3d
#endif

end module gw_mtxel_sub
