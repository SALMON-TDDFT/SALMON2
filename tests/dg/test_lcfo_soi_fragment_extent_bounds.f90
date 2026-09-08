program test_lcfo_soi_fragment_extent_bounds
  implicit none
  call exercise([2,3,2],[3,2,3])
  call exercise([4,2,3],[3,3,2])
  call exercise([3,3,3],[3,3,3])
  write(*,'(a)') 'PASS smaller, larger, and nominal LCFO-SOI complex spinor extents'
contains
  subroutine exercise(actual,nominal)
    integer, intent(in) :: actual(3),nominal(3)
    integer, parameter :: nspin=2
    complex(8), allocatable :: f_basis(:,:,:,:),orbital(:,:,:,:),hf(:,:,:,:)
    complex(8) :: block,expected
    integer :: ix,iy,iz,ispin

    allocate(f_basis(actual(1),actual(2),actual(3),nspin))
    allocate(orbital(max(actual(1),nominal(1)),max(actual(2),nominal(2)),&
      max(actual(3),nominal(3)),nspin))
    allocate(hf(max(actual(1),nominal(1)),max(actual(2),nominal(2)),&
      max(actual(3),nominal(3)),nspin))
    do ispin=1,nspin
      f_basis(:,:,:,ispin)=cmplx(real(ispin,8),0.25d0*real(ispin,8),8)
      hf(:,:,:,ispin)=cmplx(-7d0,3d0,8)
      hf(1:actual(1),1:actual(2),1:actual(3),ispin)=2d0*f_basis(:,:,:,ispin)
    end do
    orbital=cmplx(-99d0,77d0,8)
    do ispin=1,nspin
      do iz=1,actual(3);do iy=1,actual(2);do ix=1,actual(1)
        orbital(ix,iy,iz,ispin)=f_basis(ix,iy,iz,ispin)
      end do;end do;end do
    end do
    block=(0d0,0d0)
    do ispin=1,nspin
      block=block+sum(conjg(f_basis(:,:,:,ispin))*&
        hf(1:actual(1),1:actual(2),1:actual(3),ispin))
    end do
    expected=cmplx(2.125d0*sum([(real(ispin*ispin,8),ispin=1,nspin)])*&
      real(product(actual),8),0d0,8)
    if(abs(block-expected)>1d-11)error stop 'SOI Hamiltonian extent truncated'
    do ispin=1,nspin
      if(any(orbital(1:actual(1),1:actual(2),1:actual(3),ispin)/=f_basis(:,:,:,ispin))) &
        error stop 'SOI spinor basis copy lost points'
    end do
    deallocate(hf,orbital,f_basis)
  end subroutine exercise
end program test_lcfo_soi_fragment_extent_bounds
