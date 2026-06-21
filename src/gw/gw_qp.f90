!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Quasiparticle solve for the GW module -- V_xc expectation and the linearised
! QP energy update.
!
! V_xc expectation
! ----------------
! The exchange-correlation potential matrix element for state (n,k), spin
! ispin, is the diagonal expectation of the (multiplicative, real-space) xc
! potential built by initialization2_dft and stored in Vxc(ispin)%f(ix,iy,iz)
! (atomic units):
!
!   < n k | V_xc | n k > = integral_cell |u_{n,k}(r)|^2 V_xc(r) d3r
!                        = sum_{grid r} |zwf_{n,k}(r)|^2 V_xc(r) * hvol,
!
! with hvol = system%hvol = Omega / Ngrid the real-space volume element.  Since
! the orbitals are normalised as sum_grid |zwf|^2 hvol = 1, this is the cell
! average of V_xc weighted by the orbital density -- exactly the term that the
! KS Hamiltonian contributes and that GW replaces by Sigma.  The sum runs over
! the local mg slab and the locally owned (k,band) block; a comm_summation over
! icomm_rko assembles it across the r-space, k and orbital partitions (the same
! communicator calc_mtxel uses).  Vxc(ispin)%f and the orbitals
! share the SAME global grid index space (mg%is..mg%ie), the same convention
! used by every routine that combines a scalar field with the orbitals.
!
! QP solve (linearised, static)
! -----------------------------
! The linearised quasiparticle equation, replacing V_xc by the self-energy at
! the KS energy, is
!
!   eps^QP(n,k) = eps^KS(n,k) + Z * ( < Sigma(n,k) > - < V_xc(n,k) > ),
!
! with Sigma = Sigma_x (+ Sigma_c).  For the bare-exchange rung there is no
! frequency dependence, so the renormalisation Z = 1 exactly (no dSigma/domega
! term), and Sigma_c = 0.  The correction Sigma_x - V_xc is the exchange piece
! of the GW correction and widens the gap toward Hartree-Fock.
!
! No proper nouns appear in this file.
!
module gw_qp_sub
  implicit none
  private

  public :: calc_vxc_expect
  public :: solve_qp

contains

  ! --------------------------------------------------------------------------
  ! calc_vxc_expect
  !
  ! Diagonal V_xc expectation < n k | V_xc | n k > for a band window
  ! [nb_lo, nb_hi] and all k-points, spin ispin (a.u.).
  !
  ! Arguments:
  !   system  [in]  -- nk, no, nspin, hvol.
  !   info    [in]  -- parallel ownership (ik_s/e, io_s/e, im_s) + icomm_r.
  !   mg      [in]  -- local real-space slab descriptor.
  !   spsi    [in]  -- cell-periodic orbitals zwf / rwf.
  !   Vxc(:)  [in]  -- xc potential scalar fields, Vxc(ispin)%f(ix,iy,iz) (a.u.).
  !   ispin   [in]  -- spin channel.
  !   nb_lo,nb_hi[in]-- band window (states n).
  !   vxc_nk(nb_lo:nb_hi, system%nk) [out] -- V_xc expectation (a.u.).
  ! --------------------------------------------------------------------------
  subroutine calc_vxc_expect(system, info, mg, spsi, Vxc, ispin, nb_lo, nb_hi, vxc_nk)
    use structures,    only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, s_scalar
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),    intent(in)  :: system
    type(s_parallel_info), intent(in)  :: info
    type(s_rgrid),         intent(in)  :: mg
    type(s_orbital),       intent(in)  :: spsi
    type(s_scalar),        intent(in)  :: Vxc(system%nspin)
    integer,               intent(in)  :: ispin
    integer,               intent(in)  :: nb_lo, nb_hi
    real(8),               intent(out) :: vxc_nk(nb_lo:nb_hi, system%nk)

    real(8), allocatable :: vloc(:,:)     ! local-slab partial sums (n,k)
    real(8) :: hvol, acc, dens
    integer :: nk, ik, in, im, ix, iy, iz
    logical :: have_k, band_owned, use_zwf

    nk   = system%nk
    hvol = system%hvol
    im   = info%im_s
    use_zwf = allocated(spsi%zwf)

    allocate(vloc(nb_lo:nb_hi, nk))
    vloc(:,:) = 0.0d0

    do ik = 1, nk
      have_k = (ik >= info%ik_s .and. ik <= info%ik_e)
      if (.not. have_k) cycle
      do in = nb_lo, nb_hi
        band_owned = (in >= info%io_s .and. in <= info%io_e)
        if (.not. band_owned) cycle

        acc = 0.0d0
        do iz = mg%is(3), mg%ie(3)
        do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          if (use_zwf) then
            dens = dble( spsi%zwf(ix,iy,iz,ispin,in,ik,im) &
                       * conjg(spsi%zwf(ix,iy,iz,ispin,in,ik,im)) )
          else
            dens = spsi%rwf(ix,iy,iz,ispin,in,ik,im) ** 2
          end if
          acc = acc + dens * Vxc(ispin)%f(ix,iy,iz)
        end do
        end do
        end do
        vloc(in,ik) = acc * hvol
      end do
    end do

    ! Assemble across the r-space / k / orbital partitions.  With nproc_k=1 the
    ! grid slab is whole and each k/band is local, but the comm_summation keeps
    ! it correct under any partitioning (zero contributions elsewhere).
    call comm_summation(vloc, vxc_nk, (nb_hi-nb_lo+1)*nk, info%icomm_rko)

    deallocate(vloc)

  end subroutine calc_vxc_expect

  ! --------------------------------------------------------------------------
  ! solve_qp
  !
  ! Elementwise linearised quasiparticle update (one band/k/spin entry):
  !   eqp = eks + zfac * ( sigx + sigc - vxc ).
  ! For the bare-exchange rung pass zfac = 1 and sigc = 0.
  ! --------------------------------------------------------------------------
  elemental subroutine solve_qp(eks, sigx, sigc, vxc, zfac, eqp)
    implicit none
    real(8), intent(in)  :: eks, sigx, sigc, vxc, zfac
    real(8), intent(out) :: eqp
    eqp = eks + zfac * ( sigx + sigc - vxc )
  end subroutine solve_qp

end module gw_qp_sub
