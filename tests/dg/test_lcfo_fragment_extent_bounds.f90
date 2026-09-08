program test_lcfo_fragment_extent_bounds
  implicit none
  call exercise([2,3,2],[3,2,3])
  call exercise([4,2,3],[3,3,2])
  call exercise([3,3,3],[3,3,3])
  write(*,'(a)') 'PASS smaller, larger, and nominal LCFO fragment extents'
contains
  subroutine exercise(actual,nominal)
    integer, intent(in) :: actual(3),nominal(3)
    real(8), allocatable :: f_basis(:,:,:),orbital(:,:,:),hf(:,:,:)
    real(8) :: expected,block
    integer :: ix,iy,iz

    allocate(f_basis(actual(1),actual(2),actual(3)))
    allocate(orbital(max(actual(1),nominal(1)),max(actual(2),nominal(2)),max(actual(3),nominal(3))))
    allocate(hf(max(actual(1),nominal(1)),max(actual(2),nominal(2)),max(actual(3),nominal(3))))
    f_basis=1d0
    orbital=-1d0
    hf=2d0
    do iz=1,actual(3)
      do iy=1,actual(2)
        do ix=1,actual(1)
          orbital(ix,iy,iz)=f_basis(ix,iy,iz)
        end do
      end do
    end do
    block=sum(f_basis(1:actual(1),1:actual(2),1:actual(3))* &
      & hf(1:actual(1),1:actual(2),1:actual(3)))
    expected=2d0*real(product(actual),8)
    if(abs(block-expected)>1d-12) error stop 'fragment Hamiltonian extent truncated'
    if(any(orbital(1:actual(1),1:actual(2),1:actual(3))/=1d0)) &
      error stop 'fragment basis copy lost points'
    deallocate(hf,orbital,f_basis)
  end subroutine exercise
end program test_lcfo_fragment_extent_bounds
