program test_overlap_closed_blocks
  use degenerate_block_ssbe
  implicit none
  integer, parameter :: nb = 5, nk = 4, nbvec = 7
  integer :: bvec(3, nbvec), num_kgrid(3)
  complex(8) :: prod_dk(nb, nb, nbvec, nk)
  real(8) :: eigen(nb, nk)
  integer :: fixed_block_id(nb), ik, iv, n, ivp1
  logical :: isolated_ok
  integer :: nfail
  nfail = 0

  ! --- k-independent stencil: index 1=(0,0,0), 2=(+1,0,0), 3=(-1,0,0),
  !     4=(0,+1,0),5=(0,-1,0),6=(0,0,+1),7=(0,0,-1) ---
  bvec(:,1)=[0,0,0]; bvec(:,2)=[1,0,0]; bvec(:,3)=[-1,0,0]
  bvec(:,4)=[0,1,0]; bvec(:,5)=[0,-1,0]; bvec(:,6)=[0,0,1]; bvec(:,7)=[0,0,-1]
  num_kgrid = [4,1,1]
  ivp1 = 2   ! the +i1 axis column

  ! --- energies: {2,3} degenerate; 1 separated below; 4,5 above ---
  do ik = 1, nk
    eigen(1, ik) = -0.02d0
    eigen(2, ik) =  0.00d0
    eigen(3, ik) =  0.00d0
    eigen(4, ik) =  0.05d0
    eigen(5, ik) =  0.50d0
  end do

  ! --- prod_dk = identity everywhere ---
  prod_dk = (0d0, 0d0)
  do ik = 1, nk
    do iv = 1, nbvec
      do n = 1, nb
        prod_dk(n, n, iv, ik) = (1d0, 0d0)
      end do
    end do
  end do

  ! --- scramble: at ik=2, +i1 link, rotate bands (1,2) by pi/2 (full swap) ---
  !     rotation R(1,2)=[[cos, sin],[-sin, cos]], theta=pi/2 -> [[0,1],[-1,0]]
  prod_dk(1, 1, ivp1, 2) = (0d0, 0d0)
  prod_dk(1, 2, ivp1, 2) = (1d0, 0d0)
  prod_dk(2, 1, ivp1, 2) = (-1d0, 0d0)
  prod_dk(2, 2, ivp1, 2) = (0d0, 0d0)
  ! band 3 untouched (self-overlap already 1); bands 4,5 identity

  ! === T-closed-1a: default wclose merges {1,2,3} into one block ===
  call build_fixed_blocks_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, &
                                 fixed_block_id, isolated_ok)
  if (fixed_block_id(1) /= fixed_block_id(2)) then
    write(*,*) 'FAIL T1a: band 1 not merged with block'; nfail = nfail + 1
  end if
  if (fixed_block_id(2) /= fixed_block_id(3)) then
    write(*,*) 'FAIL T1a: energy block {2,3} broken'; nfail = nfail + 1
  end if
  if (fixed_block_id(4) == fixed_block_id(3)) then
    write(*,*) 'FAIL T1a: band 4 wrongly merged'; nfail = nfail + 1
  end if
  if (fixed_block_id(5) == fixed_block_id(4)) then
    write(*,*) 'FAIL T1a: band 5 wrongly merged'; nfail = nfail + 1
  end if

  ! === T-closed-1b (negative control): wclose above the injected mixing (=1)
  !     -> no overlap union; block stays the energy partition {2,3} ===
  call build_fixed_blocks_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, &
                                 fixed_block_id, isolated_ok, wclose_in=1.5d0)
  if (fixed_block_id(1) == fixed_block_id(2)) then
    write(*,*) 'FAIL T1b: band 1 merged despite high wclose (vacuous closure)'; nfail = nfail + 1
  end if
  if (fixed_block_id(2) /= fixed_block_id(3)) then
    write(*,*) 'FAIL T1b: energy block {2,3} broken'; nfail = nfail + 1
  end if

  ! === T-closed-3: build_blocks_fixed_closed caches same_block == block_id ===
  block
    integer :: block_id2(nb, nk)
    call build_blocks_fixed_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, block_id2)
    ! per-k broadcast of the k-independent partition
    if (block_id2(1, 1) /= block_id2(2, 1)) then
      write(*,*) 'FAIL T3: closed partition not reflected in block_id'; nfail = nfail + 1
    end if
    ! same_block cache must agree with the returned block_id
    if (.not. same_block(1, 2, 3)) then
      write(*,*) 'FAIL T3: same_block(1,2) false but block_id merged'; nfail = nfail + 1
    end if
    if (same_block(2, 4, 3)) then
      write(*,*) 'FAIL T3: same_block(2,4) true but they are different blocks'; nfail = nfail + 1
    end if
    if (block_id2(3, 2) /= block_id2(3, 4)) then
      write(*,*) 'FAIL T3: partition not k-independent'; nfail = nfail + 1
    end if
  end block

  ! === T-closed-1c: multi-hop closure REQUIRES >=2 fixed-point sweeps ===
  ! Band 1 -> band 2 strong (|.|^2=0.25) merges {1,2,3} in sweep 1. Band 4 gets
  ! only sub-threshold leaks from bands 1 and 2 individually (0.0625 each < 0.1),
  ! but the COMBINED {1,2,3} column weight w(4)=0.125>0.1 fires ONLY in sweep 2.
  ! A single-pass (non-iterating) closure would leave band 4 separate -> FAIL.
  block
    complex(8) :: pd(nb, nb, nbvec, nk)
    real(8) :: eg(nb, nk)
    integer :: bid(nb), nn, kk, vv
    logical :: iso
    do kk = 1, nk
      eg(1, kk) = -0.02d0; eg(2, kk) = 0.00d0; eg(3, kk) = 0.00d0
      eg(4, kk) =  0.05d0; eg(5, kk) = 0.50d0
    end do
    pd = (0d0, 0d0)
    do kk = 1, nk
      do vv = 1, nbvec
        do nn = 1, nb
          pd(nn, nn, vv, kk) = (1d0, 0d0)
        end do
      end do
    end do
    pd(1, 2, ivp1, 2) = (0.5d0,  0d0)   ! |.|^2 = 0.25   -> sweep-1 merge {1,2,3}
    pd(1, 4, ivp1, 2) = (0.25d0, 0d0)   ! |.|^2 = 0.0625 (< wclose alone)
    pd(2, 4, ivp1, 2) = (0.25d0, 0d0)   ! |.|^2 = 0.0625 (< wclose alone); sum 0.125 > wclose
    call build_fixed_blocks_closed(nb, nk, nbvec, bvec, pd, eg, num_kgrid, bid, iso)
    if (bid(1) /= bid(2) .or. bid(2) /= bid(3)) then
      write(*,*) 'FAIL T1c: sweep-1 merge {1,2,3} did not happen'; nfail = nfail + 1
    end if
    if (bid(4) /= bid(1)) then
      write(*,*) 'FAIL T1c: band 4 not merged (single-pass bug? needs 2nd sweep)'; nfail = nfail + 1
    end if
    if (bid(5) == bid(1)) then
      write(*,*) 'FAIL T1c: band 5 wrongly merged'; nfail = nfail + 1
    end if
  end block

  if (nfail == 0) then
    write(*,*) 'PASS test_overlap_closed_blocks (T-closed-1)'
  else
    write(*,*) 'FAILURES:', nfail; error stop 1
  end if
end program
