!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Bare Coulomb interaction in reciprocal space for the GW module.
!
! In Hartree atomic units the bare Coulomb matrix element in a plane-wave
! basis is:
!
!   v(q+G) = 4*pi / |q+G|^2
!
! The G-vector basis is defined by all reciprocal lattice vectors whose
! kinetic energy lies below the cutoff:
!
!   |G|^2  <=  epsilon_cutoff    (epsilon_cutoff in Rydberg; |G| in bohr^-1;
!                                  kinetic energy  = |G|^2/2 Hartree = |G|^2 Rydberg)
!
! At q+G = 0 (the "head" element: G=0 and q->0) the Coulomb kernel diverges.
! It is regularised by averaging 4*pi/|q'|^2 over the mini Brillouin zone
! (mini-BZ) surrounding the singular point, whose volume is
!
!   v_mBZ = (2*pi)^3 / (Omega * Nk)
!
! where Omega is the unit-cell volume and Nk is the number of k-points.
! The average is evaluated numerically by sub-sampling a cubic region of
! side length L = v_mBZ^(1/3) on a uniform grid, excluding the exact origin.
! The sampling density (number of sub-intervals per axis, `nsub`) is a
! caller-supplied argument so that convergence can be verified.
!
! Conventions from structures.f90 (s_dft_system):
!   system%primitive_b(1:3, j)   -- j-th reciprocal primitive vector (Cartesian, bohr^-1)
!   system%det_a                 -- unit-cell volume (bohr^3)
!   system%nk                    -- number of k-points
!
! No proper nouns appear in this file.
!
module gw_coulomb_sub
  implicit none
  private

  public :: build_gvectors
  public :: build_vcoul
  public :: gw_loop_progress

contains

  ! --------------------------------------------------------------------------
  ! Throttled progress for the (otherwise silent) q-loops: root rank only,
  ! prints at the first iteration, every ~10% of the loop, and at least every
  ! 300 s, with a wall-time ETA.  A handful of lines total -- not stdout spam.
  ! Usage: set t_start=-1; call once with idone=0 before the loop to latch the
  ! start time, then call with idone=1..itot at the bottom of each iteration.
  subroutine gw_loop_progress(label, idone, itot, t_start, t_last)
    use parallelization, only: nproc_id_global
    use communication,   only: comm_is_root
    implicit none
    character(*), intent(in)    :: label
    integer,      intent(in)    :: idone, itot
    real(8),      intent(inout) :: t_start, t_last
    integer(8) :: c, cr
    real(8)    :: now, eta
    if (.not. comm_is_root(nproc_id_global)) return
    call system_clock(c, cr)
    now = dble(c) / dble(cr)
    if (t_start < 0.0d0) then        ! latch loop start (idone=0 pre-loop call)
      t_start = now;  t_last = now;  return
    end if
    if ( idone == 1 .or. idone == itot .or. now - t_last >= 300.0d0 .or. &
         (itot >= 10 .and. mod(idone, max(1, itot/10)) == 0) ) then
      eta = (now - t_start) / dble(max(idone,1)) * dble(max(itot - idone, 0))
      write(*,'(A,A,I5,A,I5,A,F7.1,A,F7.1,A)') &
        '   [gw] ', trim(label), idone, ' /', itot, &
        '   elapsed ', (now - t_start)/60.0d0, ' min   ETA ~', eta/60.0d0, ' min'
      flush(6)
      t_last = now
    end if
  end subroutine gw_loop_progress

  ! --------------------------------------------------------------------------
  ! build_gvectors
  !
  ! Build the ordered list of reciprocal lattice vectors G whose kinetic energy
  ! |G|^2 does not exceed ecut_ry (given in Rydberg).
  !
  ! Arguments:
  !   bvec(3,3)  [in]  -- reciprocal primitive vectors in Cartesian (bohr^-1):
  !                        bvec(:,j) = j-th vector.  Pass system%primitive_b.
  !   ecut_ry    [in]  -- plane-wave cutoff in Rydberg (= epsilon_cutoff).
  !   ngmax      [in]  -- maximum number of G vectors that gvec/gg can hold.
  !   ng         [out] -- actual number of G vectors found.
  !   gvec(3,ngmax) [out] -- Cartesian G vectors (bohr^-1), G=0 first, sorted
  !                           by |G|^2 ascending.
  !   gg(ngmax)  [out] -- |G|^2 for each entry (bohr^-2).
  !
  ! Notes:
  !   The integer search range along each reciprocal axis is determined by the
  !   minimum |b_j| and the cutoff: |n_j| <= ceil( sqrt(ecut_ry) / |b_j| ).
  !   Sorting is done with a simple insertion sort (ng is typically small:
  !   a few hundred to a few thousand for typical cutoffs).
  ! --------------------------------------------------------------------------
  subroutine build_gvectors(bvec, ecut_ry, ngmax, ng, gvec, gg)
    implicit none
    real(8), intent(in)  :: bvec(3,3)
    real(8), intent(in)  :: ecut_ry
    integer, intent(in)  :: ngmax
    integer, intent(out) :: ng
    real(8), intent(out) :: gvec(3,ngmax)
    real(8), intent(out) :: gg(ngmax)

    integer :: n1, n2, n3, nmax1, nmax2, nmax3
    real(8) :: gx, gy, gz, g2, blen1, blen2, blen3, sqcut
    integer :: ig, jg
    real(8) :: gtmp(3), g2tmp
    real(8) :: pi

    pi   = acos(-1.0d0)
    sqcut = sqrt(ecut_ry)

    ! Search range: |n_j| <= ceil( sqrt(ecut) / |b_j| )
    blen1 = sqrt(bvec(1,1)**2 + bvec(2,1)**2 + bvec(3,1)**2)
    blen2 = sqrt(bvec(1,2)**2 + bvec(2,2)**2 + bvec(3,2)**2)
    blen3 = sqrt(bvec(1,3)**2 + bvec(2,3)**2 + bvec(3,3)**2)

    nmax1 = ceiling(sqcut / blen1) + 1
    nmax2 = ceiling(sqcut / blen2) + 1
    nmax3 = ceiling(sqcut / blen3) + 1

    ng = 0
    do n1 = -nmax1, nmax1
    do n2 = -nmax2, nmax2
    do n3 = -nmax3, nmax3
      gx = dble(n1)*bvec(1,1) + dble(n2)*bvec(1,2) + dble(n3)*bvec(1,3)
      gy = dble(n1)*bvec(2,1) + dble(n2)*bvec(2,2) + dble(n3)*bvec(2,3)
      gz = dble(n1)*bvec(3,1) + dble(n2)*bvec(3,2) + dble(n3)*bvec(3,3)
      g2 = gx*gx + gy*gy + gz*gz

      if (g2 > ecut_ry) cycle

      if (ng >= ngmax) then
        write(*,*) "[gw][t2] WARNING: build_gvectors ngmax exceeded; increase ngmax."
        return
      end if

      ng = ng + 1
      gvec(1,ng) = gx
      gvec(2,ng) = gy
      gvec(3,ng) = gz
      gg(ng)     = g2
    end do
    end do
    end do

    ! Insertion sort by |G|^2 ascending (G=0 at index 1 naturally first)
    do ig = 2, ng
      gtmp(1:3) = gvec(1:3, ig)
      g2tmp     = gg(ig)
      jg = ig - 1
      do while (jg >= 1 .and. gg(jg) > g2tmp)
        gvec(1:3, jg+1) = gvec(1:3, jg)
        gg(jg+1)        = gg(jg)
        jg = jg - 1
      end do
      gvec(1:3, jg+1) = gtmp(1:3)
      gg(jg+1)        = g2tmp
    end do

  end subroutine build_gvectors

  ! --------------------------------------------------------------------------
  ! build_vcoul
  !
  ! Compute the bare Coulomb matrix elements v(q+G) = 4*pi/|q+G|^2 for a
  ! given crystal momentum qvec (Cartesian, bohr^-1).
  !
  ! The singular head element (|q+G|^2 < sing_thr) is replaced by the
  ! mini-BZ-averaged value computed by numerical sub-sampling.
  !
  ! Arguments:
  !   ng         [in]  -- number of G vectors.
  !   gvec(3,ng) [in]  -- Cartesian G vectors (bohr^-1).
  !   gg(ng)     [in]  -- |G|^2 for each entry (bohr^-2).
  !   qvec(3)    [in]  -- crystal momentum in Cartesian (bohr^-1).
  !   omega      [in]  -- unit-cell volume (bohr^3); pass system%det_a.
  !   nk         [in]  -- number of k-points; pass system%nk.
  !   nsub       [in]  -- number of sub-intervals per axis for mini-BZ average.
  !                        Use nsub=10 for production; nsub=5 for a quick check.
  !   vcoul(ng)  [out] -- Coulomb matrix elements (Hartree * bohr^3).
  !
  ! Notes:
  !   sing_thr = 1.0e-8 bohr^-2; any |q+G|^2 below this is treated as
  !   singular.  For non-zero q with G=0 the normal formula applies.
  ! --------------------------------------------------------------------------
  subroutine build_vcoul(ng, gvec, gg, qvec, omega, nk, nsub, vcoul)
    implicit none
    integer, intent(in)  :: ng
    real(8), intent(in)  :: gvec(3,ng)
    real(8), intent(in)  :: gg(ng)
    real(8), intent(in)  :: qvec(3)
    real(8), intent(in)  :: omega
    integer, intent(in)  :: nk
    integer, intent(in)  :: nsub
    real(8), intent(out) :: vcoul(ng)

    integer :: ig, i1, i2, i3, npts
    real(8) :: qg2, pi, fourpi, head_val, L, dq, qx, qy, qz, r2
    real(8), parameter :: sing_thr = 1.0d-8

    pi     = acos(-1.0d0)
    fourpi = 4.0d0 * pi

    ! mini-BZ side length: v_mBZ = (2*pi)^3 / (Omega * Nk), L = v_mBZ^(1/3)
    L  = (2.0d0*pi)**3 / (omega * dble(nk))
    L  = L**(1.0d0/3.0d0)
    dq = L / dble(nsub)   ! sub-interval length

    ! --- non-singular entries ---
    do ig = 1, ng
      qg2 = (qvec(1)+gvec(1,ig))**2 + (qvec(2)+gvec(2,ig))**2 + (qvec(3)+gvec(3,ig))**2
      if (qg2 < sing_thr) then
        vcoul(ig) = 0.0d0   ! placeholder; will be overwritten below
      else
        vcoul(ig) = fourpi / qg2
      end if
    end do

    ! --- mini-BZ average for the singular head ---
    ! Check whether any entry is singular; typically only G=0 at q->0.
    do ig = 1, ng
      qg2 = (qvec(1)+gvec(1,ig))**2 + (qvec(2)+gvec(2,ig))**2 + (qvec(3)+gvec(3,ig))**2
      if (qg2 >= sing_thr) cycle

      ! Numerical average of 4*pi/|q'|^2 over a cube [-L/2, L/2]^3,
      ! excluding the exact origin, using nsub^3 uniform sub-cells.
      ! q' at the centre of each sub-cell: q'_i = (i - 0.5*(nsub+1)) * dq
      head_val = 0.0d0
      npts     = 0
      do i1 = 1, nsub
      do i2 = 1, nsub
      do i3 = 1, nsub
        qx = (dble(i1) - 0.5d0*(dble(nsub)+1.0d0)) * dq
        qy = (dble(i2) - 0.5d0*(dble(nsub)+1.0d0)) * dq
        qz = (dble(i3) - 0.5d0*(dble(nsub)+1.0d0)) * dq
        r2 = qx*qx + qy*qy + qz*qz
        if (r2 < 1.0d-30) cycle   ! exact centre of the cube (nsub odd)
        head_val = head_val + fourpi / r2
        npts = npts + 1
      end do
      end do
      end do

      if (npts > 0) then
        vcoul(ig) = head_val / dble(npts)
      else
        ! Fallback: nsub=1 corner-centred estimate
        vcoul(ig) = fourpi / (3.0d0 * (0.5d0*dq)**2)
      end if

    end do

  end subroutine build_vcoul

end module gw_coulomb_sub
