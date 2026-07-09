! src/ssbe/test/test_gicov_halo.f90
! Unit test (pure logic, NO MPI) for the gicov rho halo-exchange PLAN:
!   covariant_halo_needed  -- marks every k the covariant stencil reads
!   build_halo_lists       -- per-rank send/recv k-lists from all needed-sets
! (degenerate_block_ssbe). The MPI exchange itself cannot run here (the
! standalone build is the no-MPI dummy); its correctness rests on these
! invariants + the deterministic ascending-k packing + the Fugaku nproc_k>1
! smoke (trace_re invariant). What IS provable locally:
!
!   C (coverage):    every k READ by an INDEPENDENT re-walk of the stencil
!                    chains for rank r's slice is in local_r UNION recv_r.
!                    The re-walk is written here in the test (not calling the
!                    library's walk), so a bug in covariant_halo_needed's
!                    mirroring of the stencil is caught, not reproduced.
!   S (symmetry):    recv list of r from o == send list of o to r,
!                    elementwise (same count, same ascending k) -- this is
!                    what lets sender/receiver agree with no metadata.
!   O (ownership):   every recv k lies in the source rank's slice; every
!                    send k lies in MY slice; no self-partners.
!   E (edges):       nproc=1 -> no partners at all (required: the no-MPI
!                    dummy build aborts on isend/irecv to a real rank);
!                    nproc > nk -> empty slices handled (empty lists).
!   V (volume):      max recv k per rank printed (the halo's point: << nk).
!
! Grids tested: (8,8,8) nk=512 with nproc = 1, 4, 32, 100, 600 (600 > nk
! forces empty slices), and anisotropic (8,8,16) nk=1024 with nproc = 32
! (mixed m_max = 3,3,4 exercises per-axis shell caps).
!
! BUILD (standalone; no LAPACK call is reached but the module links it):
!   gfortran -O2 util_ssbe.f90 degenerate_block_ssbe.f90 \
!            test/test_gicov_halo.f90 -o /tmp/t_halo -framework Accelerate && /tmp/t_halo
!
program test_gicov_halo
  use util_ssbe,             only: split_range
  use degenerate_block_ssbe, only: covariant_halo_needed, build_halo_lists, &
                                    build_ik_neighbor
  implicit none

  type t_plan
    integer :: nsrc, ndst
    integer, allocatable :: src_rank(:), src_cnt(:), src_ofs(:), src_k(:)
    integer, allocatable :: dst_rank(:), dst_cnt(:), dst_ofs(:), dst_k(:)
  end type

  integer :: nfail

  nfail = 0
  call run_case((/8, 8, 8/),  1,   nfail)
  call run_case((/8, 8, 8/),  4,   nfail)
  call run_case((/8, 8, 8/),  32,  nfail)
  call run_case((/8, 8, 8/),  100, nfail)
  call run_case((/8, 8, 8/),  600, nfail)   ! nproc > nk: empty slices
  call run_case((/8, 8, 16/), 32,  nfail)   ! anisotropic: m_max = (3,3,4)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_halo checks passed."
  end if

contains

  subroutine check_true(cond, label, nfail)
    implicit none
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (cond) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a)') "FAIL  " // label
      nfail = nfail + 1
    end if
  end subroutine check_true

  subroutine run_case(num_kgrid, nproc, nfail)
    implicit none
    integer, intent(in) :: num_kgrid(3), nproc
    integer, intent(inout) :: nfail
    integer, parameter :: nbvec = 3
    integer :: bvec(3, nbvec)
    integer :: nk, r, o, p, p2, q, ik, m, axis, iv, kfwd, kbwd
    integer :: m_max(3), maxrecv, i1, i2
    integer, allocatable :: itbl_min(:), itbl_max(:), ik_neighbor(:,:), bwd(:)
    logical, allocatable :: needed_all(:,:), inrecv(:)
    logical :: okC, okS, okO, found
    character(64) :: tag
    type(t_plan), allocatable :: plan(:)

    nk = num_kgrid(1) * num_kgrid(2) * num_kgrid(3)
    bvec = 0
    bvec(1, 1) = 1;  bvec(2, 2) = 1;  bvec(3, 3) = 1
    write(tag, '(a,i0,a,i0,a,i0,a,i0)') "grid ", num_kgrid(1), "x", num_kgrid(2), &
      & "x", num_kgrid(3), " np=", nproc

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
    call split_range(1, nk, nproc, itbl_min, itbl_max)

    allocate(needed_all(nk, 0:nproc-1))
    do o = 0, nproc - 1
      call covariant_halo_needed(nk, nbvec, bvec, num_kgrid, &
                                 itbl_min(o), itbl_max(o), needed_all(:, o))
    end do

    allocate(plan(0:nproc-1))
    do r = 0, nproc - 1
      call build_halo_lists(nk, nproc, r, itbl_min, itbl_max, needed_all, &
        plan(r)%nsrc, plan(r)%src_rank, plan(r)%src_cnt, plan(r)%src_ofs, plan(r)%src_k, &
        plan(r)%ndst, plan(r)%dst_rank, plan(r)%dst_cnt, plan(r)%dst_ofs, plan(r)%dst_k)
    end do

    ! ---- C: coverage via an INDEPENDENT re-walk of the stencil chains ----
    allocate(ik_neighbor(nbvec, nk), bwd(nk), inrecv(nk))
    call build_ik_neighbor(num_kgrid, bvec, nbvec, nk, ik_neighbor)
    do axis = 1, 3
      m_max(axis) = min(4, (num_kgrid(axis) - 1) / 2)
    end do
    okC = .true.
    maxrecv = 0
    do r = 0, nproc - 1
      inrecv = .false.
      if (itbl_min(r) <= itbl_max(r)) inrecv(itbl_min(r):itbl_max(r)) = .true.  ! local
      do p = 1, plan(r)%nsrc
        do q = plan(r)%src_ofs(p) + 1, plan(r)%src_ofs(p) + plan(r)%src_cnt(p)
          inrecv(plan(r)%src_k(q)) = .true.                                     ! halo
        end do
      end do
      maxrecv = max(maxrecv, size(plan(r)%src_k))
      do axis = 1, 3
        iv = axis                     ! identity bvec: the +axis column IS column axis
        if (m_max(axis) == 0) cycle
        do ik = 1, nk
          bwd(ik_neighbor(iv, ik)) = ik
        end do
        do ik = itbl_min(r), itbl_max(r)
          kfwd = ik;  kbwd = ik
          do m = 1, m_max(axis)
            kfwd = ik_neighbor(iv, kfwd)
            kbwd = bwd(kbwd)
            if (.not. inrecv(kfwd)) okC = .false.
            if (.not. inrecv(kbwd)) okC = .false.
          end do
        end do
      end do
    end do
    call check_true(okC, "C coverage: independent re-walk within local+recv [" // trim(tag) // "]", nfail)

    ! ---- S: recv(r <- o) == send(o -> r), elementwise ----
    okS = .true.
    do r = 0, nproc - 1
      ! every recv pairing must have an exactly matching send pairing
      do p = 1, plan(r)%nsrc
        o = plan(r)%src_rank(p)
        found = .false.
        do p2 = 1, plan(o)%ndst
          if (plan(o)%dst_rank(p2) == r) then
            found = .true.
            if (plan(o)%dst_cnt(p2) /= plan(r)%src_cnt(p)) then
              okS = .false.
            else
              do q = 1, plan(r)%src_cnt(p)
                i1 = plan(r)%src_k(plan(r)%src_ofs(p) + q)
                i2 = plan(o)%dst_k(plan(o)%dst_ofs(p2) + q)
                if (i1 /= i2) okS = .false.
              end do
            end if
          end if
        end do
        if (.not. found) okS = .false.
      end do
      ! and every send pairing must have a matching recv pairing (no orphans)
      do p = 1, plan(r)%ndst
        o = plan(r)%dst_rank(p)
        found = .false.
        do p2 = 1, plan(o)%nsrc
          if (plan(o)%src_rank(p2) == r) found = .true.
        end do
        if (.not. found) okS = .false.
      end do
    end do
    call check_true(okS, "S symmetry: recv(r<-o) == send(o->r) elementwise  [" // trim(tag) // "]", nfail)

    ! ---- O: ownership + no self ----
    okO = .true.
    do r = 0, nproc - 1
      do p = 1, plan(r)%nsrc
        o = plan(r)%src_rank(p)
        if (o == r) okO = .false.
        do q = plan(r)%src_ofs(p) + 1, plan(r)%src_ofs(p) + plan(r)%src_cnt(p)
          if (plan(r)%src_k(q) < itbl_min(o) .or. plan(r)%src_k(q) > itbl_max(o)) okO = .false.
        end do
      end do
      do p = 1, plan(r)%ndst
        if (plan(r)%dst_rank(p) == r) okO = .false.
        do q = plan(r)%dst_ofs(p) + 1, plan(r)%dst_ofs(p) + plan(r)%dst_cnt(p)
          if (plan(r)%dst_k(q) < itbl_min(r) .or. plan(r)%dst_k(q) > itbl_max(r)) okO = .false.
        end do
      end do
    end do
    call check_true(okO, "O ownership: recv k owned by src, send k owned by me [" // trim(tag) // "]", nfail)

    ! ---- E: nproc=1 -> no partners ----
    if (nproc == 1) then
      call check_true(plan(0)%nsrc == 0 .and. plan(0)%ndst == 0, &
        "E edge: nproc=1 has zero partners (dummy-build safe) [" // trim(tag) // "]", nfail)
    end if

    ! ---- V: volume for the record ----
    if (nproc > 1) then
      write(*, '(a,i0,a,i0,a,f6.1,a)') "      volume: max recv k/rank = ", maxrecv, &
        & " of nk=", nk, "  (", 100d0 * dble(maxrecv) / dble(nk), "% of the old full gather)"
    end if

    deallocate(itbl_min, itbl_max, needed_all, ik_neighbor, bwd, inrecv, plan)
  end subroutine run_case

end program test_gicov_halo
