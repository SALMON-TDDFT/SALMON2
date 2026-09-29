! Exact nonzero storage of already masked, spatially local orbital columns.
module exx_sparse_orbitals
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_sparse_orbitals,sparse_pack,sparse_clear,sparse_valid,sparse_column,sparse_dot,sparse_norms
  type s_sparse_orbitals
    integer :: ng=0,no=0
    integer,allocatable :: offset(:),row(:)
    complex(8),allocatable :: value(:)
  end type
contains
  subroutine sparse_clear(a)
    implicit none
    type(s_sparse_orbitals),intent(inout) :: a
    if(allocated(a%offset))deallocate(a%offset)
    if(allocated(a%row))deallocate(a%row)
    if(allocated(a%value))deallocate(a%value)
    a%ng=0;a%no=0
  end subroutine
  subroutine sparse_pack(a,dense)
    implicit none
    type(s_sparse_orbitals),intent(inout) :: a
    complex(8),intent(in) :: dense(:,:)
    integer :: i,j,k
    call sparse_clear(a)
    a%ng=size(dense,1);a%no=size(dense,2)
    allocate(a%offset(a%no+1),a%row(count(dense/=(0d0,0d0))),a%value(count(dense/=(0d0,0d0))))
    k=1
    do j=1,a%no
      a%offset(j)=k
      do i=1,a%ng
        if(dense(i,j)==(0d0,0d0))cycle
        a%row(k)=i;a%value(k)=dense(i,j);k=k+1
      enddo
    enddo
    a%offset(a%no+1)=k
  end subroutine
  logical function sparse_valid(a,ng,no) result(valid)
    implicit none
    type(s_sparse_orbitals),intent(in) :: a
    integer,intent(in) :: ng,no
    integer :: j,k,previous
    valid=.false.
    if(a%ng/=ng.or.a%no/=no.or.ng<0.or.no<0)return
    if(.not.allocated(a%offset).or..not.allocated(a%row).or..not.allocated(a%value))return
    if(size(a%offset)/=no+1.or.size(a%row)/=size(a%value))return
    if(a%offset(1)/=1.or.a%offset(no+1)/=size(a%value)+1)return
    do j=1,no
      if(a%offset(j)<1.or.a%offset(j+1)<a%offset(j).or.a%offset(j+1)>size(a%value)+1)return
      previous=0
      do k=a%offset(j),a%offset(j+1)-1
        if(a%row(k)<=previous.or.a%row(k)>ng)return
        if(.not.ieee_is_finite(real(a%value(k))).or..not.ieee_is_finite(aimag(a%value(k))))return
        previous=a%row(k)
      enddo
    enddo
    valid=.true.
  end function
  subroutine sparse_column(a,j,column)
    implicit none
    type(s_sparse_orbitals),intent(in) :: a
    integer,intent(in) :: j
    complex(8),intent(out) :: column(:)
    integer :: k
    column=0d0
    do k=a%offset(j),a%offset(j+1)-1
      column(a%row(k))=a%value(k)
    enddo
  end subroutine
  complex(8) function sparse_dot(a,j,column) result(value)
    implicit none
    type(s_sparse_orbitals),intent(in) :: a
    integer,intent(in) :: j
    complex(8),intent(in) :: column(:)
    integer :: k
    value=0d0
    do k=a%offset(j),a%offset(j+1)-1
      value=value+conjg(a%value(k))*column(a%row(k))
    enddo
  end function
  subroutine sparse_norms(a,scale,norms)
    implicit none
    type(s_sparse_orbitals),intent(in) :: a
    real(8),intent(in) :: scale
    real(8),intent(out) :: norms(:)
    integer :: j,k
    norms=0d0
    if(scale<=0d0)return
    do j=1,a%no
      do k=a%offset(j),a%offset(j+1)-1
        norms(j)=norms(j)+abs(a%value(k)/scale)**2
      enddo
    enddo
  end subroutine
end module
