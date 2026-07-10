module lcfo_wannier_sawf_basis
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: close_sawf_candidate_basis

contains

  subroutine close_sawf_candidate_basis(candidate,npoint,ncandidate,hvol, &
      rank_tolerance,max_basis,basis,nbasis,singular_values,ok,message)
    integer, intent(in) :: npoint,ncandidate,max_basis
    real(8), intent(in) :: candidate(npoint,ncandidate),hvol,rank_tolerance
    real(8), allocatable, intent(out) :: basis(:,:),singular_values(:)
    integer, intent(out) :: nbasis
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), allocatable :: gram(:,:),eigenvalues(:),work(:)
    real(8), allocatable :: basis_tmp(:,:),singular_tmp(:)
    real(8) :: alpha,norm2
    integer :: i,j,index,info,lwork,allocation_status
    external :: dsyev

    ok=.false.; message=''; nbasis=0
    if(allocated(basis)) deallocate(basis)
    if(allocated(singular_values)) deallocate(singular_values)
    if(npoint <= 0 .or. ncandidate <= 0 .or. max_basis <= 0 .or. &
        .not.ieee_is_finite(hvol) .or. hvol <= 0.0d0 .or. &
        .not.ieee_is_finite(rank_tolerance) .or. rank_tolerance <= 0.0d0) then
      message='SAWF candidate-basis dimensions or tolerances are invalid'
      return
    end if
    if(.not.all(ieee_is_finite(candidate))) then
      message='SAWF candidate basis contains a non-finite value'
      return
    end if

    allocate(gram(ncandidate,ncandidate),eigenvalues(ncandidate),work(1), &
      stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF candidate Gram allocation failed'
      return
    end if
    gram=hvol*matmul(transpose(candidate),candidate)
    gram=0.5d0*(gram+transpose(gram))

    lwork=-1
    call dsyev('V','U',ncandidate,gram,ncandidate,eigenvalues,work,lwork,info)
    if(info /= 0 .or. .not.ieee_is_finite(work(1))) then
      message='SAWF candidate Gram eigensolver workspace query failed'
      deallocate(gram,eigenvalues,work)
      return
    end if
    lwork=max(1,int(work(1)))
    deallocate(work)
    allocate(work(lwork),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF candidate Gram workspace allocation failed'
      deallocate(gram,eigenvalues)
      return
    end if
    call dsyev('V','U',ncandidate,gram,ncandidate,eigenvalues,work,lwork,info)
    deallocate(work)
    if(info /= 0 .or. .not.all(ieee_is_finite(eigenvalues)) .or. &
        .not.all(ieee_is_finite(gram))) then
      message='SAWF candidate Gram eigensolver failed'
      deallocate(gram,eigenvalues)
      return
    end if

    nbasis=count(eigenvalues > rank_tolerance*rank_tolerance)
    if(nbasis <= 0) then
      message='SAWF candidate basis has zero numerical rank'
      nbasis=0
      deallocate(gram,eigenvalues)
      return
    end if
    if(nbasis > max_basis) then
      write(message,'(a,i0,a,i0)') 'SAWF symmetry-closed basis exceeds capacity: rank=', &
        nbasis,' capacity=',max_basis
      nbasis=0
      deallocate(gram,eigenvalues)
      return
    end if

    allocate(basis_tmp(npoint,nbasis),singular_tmp(nbasis),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF closed-basis output allocation failed'
      nbasis=0
      deallocate(gram,eigenvalues)
      return
    end if
    do i=1,nbasis
      index=ncandidate-i+1
      singular_tmp(i)=sqrt(max(0.0d0,eigenvalues(index)))
      basis_tmp(:,i)=matmul(candidate,gram(:,index))/singular_tmp(i)
      do j=1,i-1
        alpha=hvol*dot_product(basis_tmp(:,j),basis_tmp(:,i))
        basis_tmp(:,i)=basis_tmp(:,i)-alpha*basis_tmp(:,j)
      end do
      do j=1,i-1
        alpha=hvol*dot_product(basis_tmp(:,j),basis_tmp(:,i))
        basis_tmp(:,i)=basis_tmp(:,i)-alpha*basis_tmp(:,j)
      end do
      norm2=hvol*dot_product(basis_tmp(:,i),basis_tmp(:,i))
      if(.not.ieee_is_finite(norm2) .or. norm2 <= rank_tolerance*rank_tolerance) then
        message='SAWF closed-basis reorthogonalization lost a retained direction'
        nbasis=0
        deallocate(gram,eigenvalues,basis_tmp,singular_tmp)
        return
      end if
      basis_tmp(:,i)=basis_tmp(:,i)/sqrt(norm2)
    end do
    deallocate(gram,eigenvalues)
    call move_alloc(basis_tmp,basis)
    call move_alloc(singular_tmp,singular_values)
    ok=.true.
  end subroutine close_sawf_candidate_basis

end module lcfo_wannier_sawf_basis
