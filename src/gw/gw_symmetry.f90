!
! Point-group symmetry helpers for the GW Brillouin-zone sum.
!
! The dielectric matrix obeys, for a crystal point-group operation R,
!   eps_{G,G'}(R q) = eps_{R^{-1}G, R^{-1}G'}(q),  hence
!   eps_{G,G'}(q)   = eps_{R G, R G'}(R q).
! So eps^{-1} at a general q is a G-index permutation of eps^{-1} at the
! irreducible representative R q -- no orbital rotation, an integer permutation
! of the G-list (a point-group op preserves |G|, so it permutes the cutoff
! sphere).  The operations are symmorphic for the test cell (origin on an atom),
! so there is NO non-symmorphic phase.
!
! This module provides the G-index rotation and a self-test that checks the
! rotation against a direct calc_epsinv, isolating the (subtle) convention
! before it is used in the self-energy.  No proper nouns.
!
module gw_symmetry_sub
  implicit none
  private
  public :: gw_sym_selftest

contains

  ! Inverse of a 3x3 matrix (columns = reciprocal vectors), analytic cofactor.
  subroutine inv3(a, ainv)
    implicit none
    real(8), intent(in)  :: a(3,3)
    real(8), intent(out) :: ainv(3,3)
    real(8) :: det
    det =   a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
          - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
          + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    ainv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
    ainv(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/det
    ainv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))/det
    ainv(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/det
    ainv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
    ainv(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/det
    ainv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
    ainv(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/det
    ainv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
  end subroutine inv3

  ! --------------------------------------------------------------------------
  ! gw_sym_selftest
  !
  ! Reads the symmetry operations (init_sym_sub, with use_symmetry toggled on so
  ! the ops are built without disturbing the already-set full k-mesh), then for
  ! one non-trivial op R checks
  !   eps^{-1}_{G,G'}(q_test)  ==  eps^{-1}_{R G, R G'}(R q_test)
  ! by computing both sides directly (calc_epsinv) and comparing.  Root prints
  ! the operation count and the max discrepancy.  Run on a single rank.
  ! --------------------------------------------------------------------------
  subroutine gw_sym_selftest(system, info, mg, lg, spsi, esp, ecut)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital
    use gw_coulomb_sub, only: build_gvectors
    use gw_epsilon_sub, only: calc_epsinv
    use sym_sub,        only: SymMatB, init_sym_sub, use_symmetry, read_sw_symmetry
    use parallelization,only: nproc_id_global
    use communication,  only: comm_is_root
    implicit none
    type(s_dft_system),    intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid),         intent(in) :: mg, lg
    type(s_orbital),       intent(in) :: spsi
    real(8),               intent(in) :: esp(system%no, system%nk, system%nspin)
    real(8),               intent(in) :: ecut

    integer, parameter :: ngmax = 200000
    real(8), allocatable :: gvec(:,:), gg(:)
    complex(8), allocatable :: eps_test(:,:), eps_irr(:,:)
    integer, allocatable :: mind(:,:), gperm(:)
    real(8) :: binv(3,3), mr(3), gtgt(3), qtest(3), qirr(3), qf(3)
    integer :: ng, nsym, isym, ig, jg, kg, m1, m2, m3
    real(8) :: resid, errmax
    logical :: ok
    logical :: rootp

    rootp = comm_is_root(nproc_id_global)

    ! Build the operations (full k-mesh already set; this only fills SymMatB).
    ! read_sw_symmetry sets BOTH use_symmetry and the per-direction filter flags
    ! that init_sym_sub uses to keep operations (without it only the identity
    ! survives the filter).  Called here, after init_dft, so the full k-mesh that
    ! init_dft already built is undisturbed.
    call read_sw_symmetry('yyy')
    call init_sym_sub(system%primitive_a, system%primitive_b)
    nsym = size(SymMatB, 3)
    use_symmetry = .false.   ! ops persist in SymMatB; keep the rest of GW unaffected

    ! eps G-set.
    allocate(gvec(3,ngmax), gg(ngmax))
    call build_gvectors(system%primitive_b, ecut, ngmax, ng, gvec, gg)

    ! Integer reciprocal index m of each G (G = primitive_b . m).
    call inv3(system%primitive_b, binv)
    allocate(mind(3,ng))
    do ig = 1, ng
      mind(1,ig) = nint(binv(1,1)*gvec(1,ig)+binv(1,2)*gvec(2,ig)+binv(1,3)*gvec(3,ig))
      mind(2,ig) = nint(binv(2,1)*gvec(1,ig)+binv(2,2)*gvec(2,ig)+binv(2,3)*gvec(3,ig))
      mind(3,ig) = nint(binv(3,1)*gvec(1,ig)+binv(3,2)*gvec(2,ig)+binv(3,3)*gvec(3,ig))
    end do

    ! Use the first non-trivial operation.
    isym = min(2, nsym)

    ! G-permutation: gperm(ig) = index of R.G_ig (m' = SymMatB . m, integer).
    allocate(gperm(ng))
    gperm(:) = 0
    do ig = 1, ng
      mr(1) = SymMatB(1,1,isym)*mind(1,ig)+SymMatB(1,2,isym)*mind(2,ig)+SymMatB(1,3,isym)*mind(3,ig)
      mr(2) = SymMatB(2,1,isym)*mind(1,ig)+SymMatB(2,2,isym)*mind(2,ig)+SymMatB(2,3,isym)*mind(3,ig)
      mr(3) = SymMatB(3,1,isym)*mind(1,ig)+SymMatB(3,2,isym)*mind(2,ig)+SymMatB(3,3,isym)*mind(3,ig)
      m1 = nint(mr(1)); m2 = nint(mr(2)); m3 = nint(mr(3))
      do jg = 1, ng
        if (mind(1,jg)==m1 .and. mind(2,jg)==m2 .and. mind(3,jg)==m3) then
          gperm(ig) = jg; exit
        end if
      end do
      if (gperm(ig) == 0) then
        if (rootp) write(*,*) "[gw][sym] FATAL: R.G not in G-set, ig=", ig
        deallocate(gvec,gg,mind,gperm); return
      end if
    end do

    ! q_test (a non-zero mesh q) and q_irr = R q_test.
    qtest(1:3) = system%vec_k(1:3,2) - system%vec_k(1:3,1)
    qf(1) = binv(1,1)*qtest(1)+binv(1,2)*qtest(2)+binv(1,3)*qtest(3)
    qf(2) = binv(2,1)*qtest(1)+binv(2,2)*qtest(2)+binv(2,3)*qtest(3)
    qf(3) = binv(3,1)*qtest(1)+binv(3,2)*qtest(2)+binv(3,3)*qtest(3)
    ! mr = SymMatB . qf  (rotated fractional), then back to Cartesian.
    mr(1) = SymMatB(1,1,isym)*qf(1)+SymMatB(1,2,isym)*qf(2)+SymMatB(1,3,isym)*qf(3)
    mr(2) = SymMatB(2,1,isym)*qf(1)+SymMatB(2,2,isym)*qf(2)+SymMatB(2,3,isym)*qf(3)
    mr(3) = SymMatB(3,1,isym)*qf(1)+SymMatB(3,2,isym)*qf(2)+SymMatB(3,3,isym)*qf(3)
    qirr(1) = system%primitive_b(1,1)*mr(1)+system%primitive_b(1,2)*mr(2)+system%primitive_b(1,3)*mr(3)
    qirr(2) = system%primitive_b(2,1)*mr(1)+system%primitive_b(2,2)*mr(2)+system%primitive_b(2,3)*mr(3)
    qirr(3) = system%primitive_b(3,1)*mr(1)+system%primitive_b(3,2)*mr(2)+system%primitive_b(3,3)*mr(3)

    allocate(eps_test(ng,ng), eps_irr(ng,ng))
    call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, 1, qtest, 1, eps_test, ok=ok)
    call calc_epsinv(system, info, mg, lg, spsi, esp, gvec, gg, ng, 1, qirr,  1, eps_irr,  ok=ok)

    ! eps^{-1}(q_test)_{G,G'} should equal eps^{-1}(q_irr)_{R G, R G'}.
    errmax = 0.0d0
    do jg = 1, ng
      do ig = 1, ng
        resid = abs( eps_test(ig,jg) - eps_irr(gperm(ig),gperm(jg)) )
        if (resid > errmax) errmax = resid
      end do
    end do

    if (rootp) then
      write(*,*)
      write(*,*) "[gw][sym] operations read (nsym) =", nsym
      write(*,'(A,I3,A,I6)') "  [gw][sym] test op isym=", isym, "   eps G-set ng=", ng
      write(*,'(A,3F10.5)') "  [gw][sym] q_test  =", qtest
      write(*,'(A,3F10.5)') "  [gw][sym] R.q_test=", qirr
      write(*,'(A,ES14.6)') "  [gw][sym] max| eps^-1(q) - rot[eps^-1(Rq)] | =", errmax
      write(*,'(A,L2)')     "  [gw][sym] G-rotation convention OK (<1e-8)? ", (errmax < 1.0d-8)
      write(*,*)
    end if

    deallocate(gvec, gg, mind, gperm, eps_test, eps_irr)
  end subroutine gw_sym_selftest

end module gw_symmetry_sub
