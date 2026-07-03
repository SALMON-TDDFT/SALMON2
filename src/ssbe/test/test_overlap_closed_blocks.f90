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

  if (nfail == 0) then
    write(*,*) 'PASS test_overlap_closed_blocks (T-closed-1)'
  else
    write(*,*) 'FAILURES:', nfail; error stop 1
  end if
end program
