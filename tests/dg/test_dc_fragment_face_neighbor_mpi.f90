program test_dc_fragment_face_neighbor_mpi
  use mpi_f08
  use dc_fragment_geometry, only: find_fragment_face_neighbor
  implicit none
  integer :: ierr, rank, current, neighbor, matches
  integer :: origins2(3,2), extents2(3,2), total2(3)
  integer :: origins8(3,8), extents8(3,8), total8(3)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (rank > 1) error stop 'fixture requires exactly two ranks'

  total2 = [6, 1, 1]
  origins2(:,1) = [0, 0, 0]
  origins2(:,2) = [2, 0, 0]
  extents2(:,1) = [2, 1, 1]
  extents2(:,2) = [4, 1, 1]
  current = rank + 1
  call find_fragment_face_neighbor(origins2, extents2, total2, current, [1,0,0], neighbor, matches)
  call require(matches == 1 .and. neighbor == 3-current, 'unequal +/wrap neighbor')
  call find_fragment_face_neighbor(origins2, extents2, total2, current, [-1,0,0], neighbor, matches)
  call require(matches == 1 .and. neighbor == 3-current, 'unequal -/wrap neighbor')

  if (rank == 0) then
    extents2(1,:) = [3,3]
    origins2(1,:) = [0,3]
    call find_fragment_face_neighbor(origins2, extents2, total2, 1, [1,0,0], neighbor, matches)
    call require(matches == 1 .and. neighbor == 2, 'equal-width + neighbor')
    call find_fragment_face_neighbor(origins2, extents2, total2, 1, [-1,0,0], neighbor, matches)
    call require(matches == 1 .and. neighbor == 2, 'equal-width wrap neighbor')

    call make_3d_fixture(origins8, extents8)
    total8 = [6,8,10]
    call find_fragment_face_neighbor(origins8, extents8, total8, 1, [1,1,1], neighbor, matches)
    call require(matches == 1 .and. neighbor == 8, '3D positive corner')
    call find_fragment_face_neighbor(origins8, extents8, total8, 1, [-1,1,0], neighbor, matches)
    call require(matches == 1 .and. neighbor == 7, '3D mixed wrap edge')
    call find_fragment_face_neighbor(origins8, extents8, total8, 8, [1,1,1], neighbor, matches)
    call require(matches == 1 .and. neighbor == 1, '3D full wrap corner')

    origins2(1,:) = [0,3]
    extents2(1,:) = [2,3]
    call find_fragment_face_neighbor(origins2, extents2, total2, 1, [1,0,0], neighbor, matches)
    call require(matches == 0, 'gap rejection')

    block
      integer :: origins3(3,3), extents3(3,3)
      origins3(:,1) = [0,0,0]
      origins3(:,2) = [2,0,0]
      origins3(:,3) = [2,0,0]
      extents3(:,1) = [2,1,1]
      extents3(:,2) = [4,1,1]
      extents3(:,3) = [4,1,1]
      call find_fragment_face_neighbor(origins3, extents3, total2, 1, [1,0,0], neighbor, matches)
      call require(matches == 2, 'ambiguous overlap rejection')
    end block
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (rank == 0) print *, 'PASS production fragment face topology on 2 MPI ranks'
  call MPI_Finalize(ierr)

contains

  subroutine make_3d_fixture(origins, extents)
    integer, intent(out) :: origins(3,8), extents(3,8)
    integer :: ix, iy, iz, f
    f = 0
    do ix = 1,2
      do iy = 1,2
        do iz = 1,2
          f = f + 1
          origins(:,f) = [merge(0,2,ix==1), merge(0,3,iy==1), merge(0,4,iz==1)]
          extents(:,f) = [merge(2,4,ix==1), merge(3,5,iy==1), merge(4,6,iz==1)]
        end do
      end do
    end do
  end subroutine make_3d_fixture

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a,i0,2a)') 'rank ', rank, ': FAIL ', trim(label)
      error stop 1
    end if
  end subroutine require

end program test_dc_fragment_face_neighbor_mpi
