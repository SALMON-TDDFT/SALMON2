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
  public :: gw_grid_perm_selftest
  public :: gw_symmetrize_orbitals
  public :: gw_sym_init_ops
  public :: build_g_perm
  public :: build_ibz_map

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

  ! --------------------------------------------------------------------------
  ! gw_grid_perm_selftest
  !
  ! Build, for each point-group operation R, the real-space grid-index map that
  ! realises u_{Rk}(r) = u_k(R^{-1} r):  perm_R(i) = grid index of R^{-1} r_i,
  ! with r_i the fractional coordinate (i-1)/N of grid point i (origin at i=1,
  ! periodic).  For the cubic grid + signed-permutation (Td) operations this is
  ! an exact index permutation (sign flips fold via modulo).  Then check it is a
  ! BIJECTION (catches a wrong fold/reflection) and that the identity operation
  ! gives the identity map, and that the group closes under composition.  This
  ! isolates the trickiest part of the orbital rotation before it touches data.
  ! Requires a cubic grid (N1=N2=N3) so axis swaps map the grid to itself.
  ! --------------------------------------------------------------------------
  subroutine gw_grid_perm_selftest(system, lg)
    use structures,     only: s_dft_system, s_rgrid
    use sym_sub,        only: SymMatA, init_sym_sub, use_symmetry, read_sw_symmetry
    use parallelization,only: nproc_id_global
    use communication,  only: comm_is_root
    implicit none
    type(s_dft_system), intent(in) :: system
    type(s_rgrid),      intent(in) :: lg

    integer, allocatable :: perm(:,:), hit(:), comp(:)
    real(8) :: sa(3,3), sainv(3,3), f(3), fp(3)
    integer :: n1, n2, n3, ntot, nsym, isym, jsym, ksym
    integer :: ix, iy, iz, jx, jy, jz, i, nbad, nident, nclose, nfound
    logical :: rootp, ok

    rootp = comm_is_root(nproc_id_global)

    call read_sw_symmetry('yyy')
    call init_sym_sub(system%primitive_a, system%primitive_b)
    use_symmetry = .false.
    nsym = size(SymMatA, 3)

    n1 = lg%num(1); n2 = lg%num(2); n3 = lg%num(3)
    ntot = n1*n2*n3
    if (n1 /= n2 .or. n2 /= n3) then
      if (rootp) write(*,*) "[gw][sym2] non-cubic grid; axis-swap ops need N1=N2=N3"
    end if

    allocate(perm(ntot,nsym), hit(ntot), comp(ntot))

    do isym = 1, nsym
      sa = SymMatA(1:3,1:3,isym)
      call inv3(sa, sainv)
      do iz = 1, n3
      do iy = 1, n2
      do ix = 1, n1
        i = ix + (iy-1)*n1 + (iz-1)*n1*n2
        f(1) = dble(ix-1)/dble(n1)
        f(2) = dble(iy-1)/dble(n2)
        f(3) = dble(iz-1)/dble(n3)
        fp(1) = sainv(1,1)*f(1)+sainv(1,2)*f(2)+sainv(1,3)*f(3)
        fp(2) = sainv(2,1)*f(1)+sainv(2,2)*f(2)+sainv(2,3)*f(3)
        fp(3) = sainv(3,1)*f(1)+sainv(3,2)*f(2)+sainv(3,3)*f(3)
        jx = modulo(nint(fp(1)*dble(n1)), n1) + 1
        jy = modulo(nint(fp(2)*dble(n2)), n2) + 1
        jz = modulo(nint(fp(3)*dble(n3)), n3) + 1
        perm(i,isym) = jx + (jy-1)*n1 + (jz-1)*n1*n2
      end do
      end do
      end do
    end do

    ! (1) each map a bijection?
    nbad = 0
    do isym = 1, nsym
      hit(:) = 0
      do i = 1, ntot
        hit(perm(i,isym)) = hit(perm(i,isym)) + 1
      end do
      if (any(hit(:) /= 1)) nbad = nbad + 1
    end do

    ! (2) the identity operation -> identity map?
    nident = 0
    do isym = 1, nsym
      ok = .true.
      do i = 1, ntot
        if (abs(SymMatA(1,1,isym)-1d0)+abs(SymMatA(2,2,isym)-1d0)+abs(SymMatA(3,3,isym)-1d0) < 1d-8 &
            .and. abs(SymMatA(1,2,isym))+abs(SymMatA(1,3,isym))+abs(SymMatA(2,3,isym)) < 1d-8) then
          if (perm(i,isym) /= i) ok = .false.
        end if
      end do
      if (.not. ok) nident = nident + 1
    end do

    ! (3) closure: perm_isym ( perm_jsym (i) ) must equal some perm_ksym for all i.
    nclose = 0
    do isym = 1, nsym
      do jsym = 1, nsym
        do i = 1, ntot
          comp(i) = perm(perm(i,jsym), isym)   ! apply jsym then isym
        end do
        nfound = 0
        do ksym = 1, nsym
          if (all(comp(:) == perm(:,ksym))) then; nfound = 1; exit; end if
        end do
        if (nfound == 0) nclose = nclose + 1
      end do
    end do

    if (rootp) then
      write(*,*)
      write(*,'(A,I3,A,3I4)') " [gw][sym2] grid-perm: nsym=", nsym, "  grid=", n1,n2,n3
      write(*,'(A,I4,A)')     " [gw][sym2] non-bijective maps     =", nbad,   " (must be 0)"
      write(*,'(A,I4,A)')     " [gw][sym2] identity-op mismatches  =", nident, " (must be 0)"
      write(*,'(A,I4,A)')     " [gw][sym2] non-closed compositions =", nclose, " (must be 0)"
      write(*,'(A,L2)')       " [gw][sym2] grid permutation group OK? ", (nbad==0 .and. nident==0 .and. nclose==0)
      write(*,*)
    end if

    deallocate(perm, hit, comp)
  end subroutine gw_grid_perm_selftest

  ! --------------------------------------------------------------------------
  ! gw_symmetrize_orbitals
  !
  ! Make the full-mesh orbitals EXACTLY symmetric by overwriting every star
  ! member with the grid-rotated copy of its star representative:
  !   u_k(r) = u_{k_irr}(R^{-1} r)   for k = R k_irr,
  ! realised on the grid via the validated index permutation perm_R, and the
  ! star energies/occupations set equal to the representative's (these are
  ! point-group invariants).  After this the dielectric matrix is symmetric to
  ! machine precision, so the G-rotation reduction is EXACT.  This is the
  ! orbital-rotation that replaces the (broken on the full mesh) density
  ! symmetrisation.  Single rank (full mesh on one rank).
  ! --------------------------------------------------------------------------
  subroutine gw_symmetrize_orbitals(system, info, lg, spsi, energy)
    use structures,     only: s_dft_system, s_parallel_info, s_rgrid, s_orbital, s_dft_energy
    use sym_sub,        only: SymMatA, SymMatB, init_sym_sub, use_symmetry, read_sw_symmetry
    use parallelization,only: nproc_id_global
    use communication,  only: comm_is_root
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: lg
    type(s_orbital),       intent(inout) :: spsi
    type(s_dft_energy),    intent(inout) :: energy

    integer, allocatable :: perm(:,:), rep(:), kop(:), gx(:), gy(:), gz(:)
    logical, allocatable :: assigned(:)
    real(8) :: sa(3,3), sainv(3,3), binv(3,3), f(3), fp(3), kf(3), rf(3), d
    integer :: n1, n2, n3, ntot, nsym, isym, ix, iy, iz, i, j
    integer :: nk, ik, kk, k, r, is, io, im, nstar, nrep
    integer :: jx, jy, jz
    logical :: rootp, use_zwf
    real(8), allocatable :: rbuf(:)
    complex(8), allocatable :: zbuf(:)

    rootp = comm_is_root(nproc_id_global)
    call read_sw_symmetry('yyy')
    call init_sym_sub(system%primitive_a, system%primitive_b)
    use_symmetry = .false.
    nsym = size(SymMatB, 3)

    n1 = lg%num(1); n2 = lg%num(2); n3 = lg%num(3); ntot = n1*n2*n3
    nk = system%nk
    use_zwf = allocated(spsi%zwf)

    ! ---- grid-index permutation perm_R(i) = grid point of R^{-1} r_i ----------
    allocate(perm(ntot,nsym), gx(ntot), gy(ntot), gz(ntot))
    do iz = 1, n3
    do iy = 1, n2
    do ix = 1, n1
      i = ix + (iy-1)*n1 + (iz-1)*n1*n2
      gx(i) = ix; gy(i) = iy; gz(i) = iz
    end do
    end do
    end do
    do isym = 1, nsym
      sa = SymMatA(1:3,1:3,isym)
      call inv3(sa, sainv)
      do i = 1, ntot
        f(1) = dble(gx(i)-1)/dble(n1)
        f(2) = dble(gy(i)-1)/dble(n2)
        f(3) = dble(gz(i)-1)/dble(n3)
        fp(1) = sainv(1,1)*f(1)+sainv(1,2)*f(2)+sainv(1,3)*f(3)
        fp(2) = sainv(2,1)*f(1)+sainv(2,2)*f(2)+sainv(2,3)*f(3)
        fp(3) = sainv(3,1)*f(1)+sainv(3,2)*f(2)+sainv(3,3)*f(3)
        jx = modulo(nint(fp(1)*dble(n1)), n1) + 1
        jy = modulo(nint(fp(2)*dble(n2)), n2) + 1
        jz = modulo(nint(fp(3)*dble(n3)), n3) + 1
        perm(i,isym) = jx + (jy-1)*n1 + (jz-1)*n1*n2
      end do
    end do

    ! ---- star map: each k -> (rep, op) with k = R_op . rep ; rep -> kop=0 -----
    call inv3(system%primitive_b, binv)
    allocate(rep(nk), kop(nk), assigned(nk))
    assigned(:) = .false.
    rep(:) = 0; kop(:) = 0
    do ik = 1, nk
      if (assigned(ik)) cycle
      rep(ik) = ik; kop(ik) = 0; assigned(ik) = .true.
      kf(1) = binv(1,1)*system%vec_k(1,ik)+binv(1,2)*system%vec_k(2,ik)+binv(1,3)*system%vec_k(3,ik)
      kf(2) = binv(2,1)*system%vec_k(1,ik)+binv(2,2)*system%vec_k(2,ik)+binv(2,3)*system%vec_k(3,ik)
      kf(3) = binv(3,1)*system%vec_k(1,ik)+binv(3,2)*system%vec_k(2,ik)+binv(3,3)*system%vec_k(3,ik)
      do isym = 1, nsym
        rf(1) = SymMatB(1,1,isym)*kf(1)+SymMatB(1,2,isym)*kf(2)+SymMatB(1,3,isym)*kf(3)
        rf(2) = SymMatB(2,1,isym)*kf(1)+SymMatB(2,2,isym)*kf(2)+SymMatB(2,3,isym)*kf(3)
        rf(3) = SymMatB(3,1,isym)*kf(1)+SymMatB(3,2,isym)*kf(2)+SymMatB(3,3,isym)*kf(3)
        ! find the mesh point with fractional coord rf (mod 1)
        do kk = 1, nk
          f(1) = binv(1,1)*system%vec_k(1,kk)+binv(1,2)*system%vec_k(2,kk)+binv(1,3)*system%vec_k(3,kk)
          f(2) = binv(2,1)*system%vec_k(1,kk)+binv(2,2)*system%vec_k(2,kk)+binv(2,3)*system%vec_k(3,kk)
          f(3) = binv(3,1)*system%vec_k(1,kk)+binv(3,2)*system%vec_k(2,kk)+binv(3,3)*system%vec_k(3,kk)
          d = abs(modulo(rf(1)-f(1)+0.5d0,1d0)-0.5d0) &
            + abs(modulo(rf(2)-f(2)+0.5d0,1d0)-0.5d0) &
            + abs(modulo(rf(3)-f(3)+0.5d0,1d0)-0.5d0)
          if (d < 1d-6) then
            if (.not. assigned(kk)) then
              rep(kk) = ik; kop(kk) = isym; assigned(kk) = .true.
            end if
            exit
          end if
        end do
      end do
    end do
    nrep  = count(kop(:) == 0)
    nstar = nk - nrep

    ! ---- overwrite each non-representative k by the rotated representative ----
    allocate(zbuf(ntot), rbuf(ntot))
    im = info%im_s
    do k = 1, nk
      if (kop(k) == 0) cycle           ! representative: leave as is
      r = rep(k); isym = kop(k)
      do is = 1, system%nspin
        do io = 1, system%no
          if (use_zwf) then
            do i = 1, ntot
              zbuf(i) = spsi%zwf(gx(perm(i,isym)), gy(perm(i,isym)), gz(perm(i,isym)), is, io, r, im)
            end do
            do i = 1, ntot
              spsi%zwf(gx(i), gy(i), gz(i), is, io, k, im) = zbuf(i)
            end do
          else
            do i = 1, ntot
              rbuf(i) = spsi%rwf(gx(perm(i,isym)), gy(perm(i,isym)), gz(perm(i,isym)), is, io, r, im)
            end do
            do i = 1, ntot
              spsi%rwf(gx(i), gy(i), gz(i), is, io, k, im) = rbuf(i)
            end do
          end if
        end do
        ! energies and occupations are point-group invariants
        energy%esp(:,k,is)   = energy%esp(:,r,is)
        system%rocc(:,k,is)  = system%rocc(:,r,is)
      end do
    end do

    if (rootp) then
      write(*,*)
      write(*,'(A,I4,A,I4,A,I4)') " [gw][sym3] orbital symmetrise: nk=", nk, &
        "  representatives=", nrep, "  star members rotated=", nstar
      write(*,*)
    end if

    deallocate(perm, gx, gy, gz, rep, kop, assigned, zbuf, rbuf)
  end subroutine gw_symmetrize_orbitals

  ! --------------------------------------------------------------------------
  ! gw_sym_init_ops : populate SymMatA/SymMatB from sym.dat without disturbing
  ! the already-built full k-mesh (read_sw_symmetry sets the per-direction
  ! filter; init_sym_sub reads the operations).  Returns the operation count.
  ! Safe to call repeatedly (init_sym_sub guards on its own flag).
  ! --------------------------------------------------------------------------
  subroutine gw_sym_init_ops(system, nsym)
    use structures, only: s_dft_system
    use sym_sub,    only: SymMatB, init_sym_sub, use_symmetry, read_sw_symmetry
    implicit none
    type(s_dft_system), intent(in)  :: system
    integer,            intent(out) :: nsym
    call read_sw_symmetry('yyy')
    call init_sym_sub(system%primitive_a, system%primitive_b)
    use_symmetry = .false.
    nsym = size(SymMatB, 3)
  end subroutine gw_sym_init_ops

  ! --------------------------------------------------------------------------
  ! build_g_perm : gperm(ig,isym) = index in the G-list of R_isym . G_ig, i.e.
  ! the integer reciprocal index m' = SymMatB(isym) . m_ig.  Point-group ops
  ! preserve |G| so they permute the cutoff sphere.  Used for the dielectric
  ! G-rotation eps^-1_{G,G'}(q) = eps^-1_{R G, R G'}(q_irr) (validated to 3e-15).
  ! Requires SymMatB populated (gw_sym_init_ops).
  ! --------------------------------------------------------------------------
  subroutine build_g_perm(bmat, gvec, ng, gperm)
    use sym_sub, only: SymMatB
    implicit none
    real(8), intent(in)  :: bmat(3,3)
    integer, intent(in)  :: ng
    real(8), intent(in)  :: gvec(3,ng)
    integer, allocatable, intent(out) :: gperm(:,:)
    real(8) :: binv(3,3), mr(3)
    integer, allocatable :: mind(:,:)
    integer :: nsym, isym, ig, jg, m1, m2, m3
    nsym = size(SymMatB,3)
    call inv3(bmat, binv)
    allocate(mind(3,ng), gperm(ng,nsym))
    do ig = 1, ng
      mind(1,ig) = nint(binv(1,1)*gvec(1,ig)+binv(1,2)*gvec(2,ig)+binv(1,3)*gvec(3,ig))
      mind(2,ig) = nint(binv(2,1)*gvec(1,ig)+binv(2,2)*gvec(2,ig)+binv(2,3)*gvec(3,ig))
      mind(3,ig) = nint(binv(3,1)*gvec(1,ig)+binv(3,2)*gvec(2,ig)+binv(3,3)*gvec(3,ig))
    end do
    gperm(:,:) = 0
    do isym = 1, nsym
      do ig = 1, ng
        mr(1) = SymMatB(1,1,isym)*mind(1,ig)+SymMatB(1,2,isym)*mind(2,ig)+SymMatB(1,3,isym)*mind(3,ig)
        mr(2) = SymMatB(2,1,isym)*mind(1,ig)+SymMatB(2,2,isym)*mind(2,ig)+SymMatB(2,3,isym)*mind(3,ig)
        mr(3) = SymMatB(3,1,isym)*mind(1,ig)+SymMatB(3,2,isym)*mind(2,ig)+SymMatB(3,3,isym)*mind(3,ig)
        m1 = nint(mr(1)); m2 = nint(mr(2)); m3 = nint(mr(3))
        do jg = 1, ng
          if (mind(1,jg)==m1 .and. mind(2,jg)==m2 .and. mind(3,jg)==m3) then
            gperm(ig,isym) = jg; exit
          end if
        end do
      end do
    end do
    deallocate(mind)
  end subroutine build_g_perm

  ! --------------------------------------------------------------------------
  ! build_ibz_map : reduce a set of momentum vectors vec(3,nv) (Cartesian) to
  ! the irreducible set.  rep(i) = the representative index, kop(i) = the op with
  ! vec(:,i) = R_kop . vec(:,rep(i))  (mod reciprocal lattice).  kop=0 marks a
  ! representative.  Requires SymMatB populated.
  ! --------------------------------------------------------------------------
  subroutine build_ibz_map(bmat, vec, nv, rep, kop, nrep)
    use sym_sub, only: SymMatB
    implicit none
    real(8), intent(in)  :: bmat(3,3)
    integer, intent(in)  :: nv
    real(8), intent(in)  :: vec(3,nv)
    integer, intent(out) :: rep(nv), kop(nv)
    integer, intent(out) :: nrep
    real(8) :: binv(3,3), vf(3), rf(3), f(3), d
    integer :: nsym, isym, i, j
    logical, allocatable :: assigned(:)
    nsym = size(SymMatB,3)
    call inv3(bmat, binv)
    allocate(assigned(nv)); assigned(:) = .false.
    rep(:) = 0; kop(:) = 0; nrep = 0
    do i = 1, nv
      if (assigned(i)) cycle
      rep(i) = i; kop(i) = 0; assigned(i) = .true.; nrep = nrep + 1
      vf(1) = binv(1,1)*vec(1,i)+binv(1,2)*vec(2,i)+binv(1,3)*vec(3,i)
      vf(2) = binv(2,1)*vec(1,i)+binv(2,2)*vec(2,i)+binv(2,3)*vec(3,i)
      vf(3) = binv(3,1)*vec(1,i)+binv(3,2)*vec(2,i)+binv(3,3)*vec(3,i)
      do isym = 1, nsym
        rf(1) = SymMatB(1,1,isym)*vf(1)+SymMatB(1,2,isym)*vf(2)+SymMatB(1,3,isym)*vf(3)
        rf(2) = SymMatB(2,1,isym)*vf(1)+SymMatB(2,2,isym)*vf(2)+SymMatB(2,3,isym)*vf(3)
        rf(3) = SymMatB(3,1,isym)*vf(1)+SymMatB(3,2,isym)*vf(2)+SymMatB(3,3,isym)*vf(3)
        do j = 1, nv
          f(1) = binv(1,1)*vec(1,j)+binv(1,2)*vec(2,j)+binv(1,3)*vec(3,j)
          f(2) = binv(2,1)*vec(1,j)+binv(2,2)*vec(2,j)+binv(2,3)*vec(3,j)
          f(3) = binv(3,1)*vec(1,j)+binv(3,2)*vec(2,j)+binv(3,3)*vec(3,j)
          ! EXACT match (no mod-reciprocal-lattice fold): on the MP-shifted mesh R
          ! maps every point to another mesh point exactly, so grouping by the
          ! exact vector avoids umklapp ambiguity (q and q+G both present), which
          ! the pure G-rotation could not undo.
          d = abs(rf(1)-f(1)) + abs(rf(2)-f(2)) + abs(rf(3)-f(3))
          if (d < 1d-6) then
            if (.not. assigned(j)) then
              rep(j) = i; kop(j) = isym; assigned(j) = .true.
            end if
            exit
          end if
        end do
      end do
    end do
    deallocate(assigned)
  end subroutine build_ibz_map

end module gw_symmetry_sub
