program test_sparse_orbitals
  use exx_sparse_orbitals
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(s_sparse_orbitals) :: a,b
  complex(8) :: dense(8,3),column(8),probe(8)
  real(8) :: norms(3)
  integer :: j,k
  dense=0d0
  dense(2,1)=(1d0,2d0);dense(8,1)=(1d-200,-1d-200)
  dense(4,3)=(-3d0,.5d0)
  do k=1,8
    probe(k)=cmplx(k,1-k,8)
  enddo
  call sparse_pack(a,dense)
  if(.not.sparse_valid(a,8,3))error stop 'valid storage rejected'
  if(size(a%value)/=3)error stop 'nonzero or tiny coefficient lost'
  do j=1,3
    call sparse_column(a,j,column)
    if(any(column/=dense(:,j)))error stop 'roundtrip mismatch'
    if(abs(sparse_dot(a,j,probe)-dot_product(dense(:,j),probe))>1d-14)error stop 'dot mismatch'
  enddo
  call sparse_norms(a,4d0,norms)
  if(maxval(abs(norms-sum(abs(dense/4d0)**2,dim=1)))>1d-14)error stop 'norm mismatch'
  b=a
  call sparse_clear(a)
  call sparse_column(b,1,column)
  if(any(column/=dense(:,1)))error stop 'copy aliases cleared storage'
  if(sparse_valid(a,8,3))error stop 'cleared storage accepted'
  b%row(1)=0
  if(sparse_valid(b,8,3))error stop 'invalid row accepted'
  b%row(1)=2;b%offset(2)=0
  if(sparse_valid(b,8,3))error stop 'invalid offset accepted'
  call sparse_pack(b,dense)
  b%value(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  if(sparse_valid(b,8,3))error stop 'nonfinite value accepted'
  call sparse_pack(a,dense(:,:0))
  if(.not.sparse_valid(a,8,0))error stop 'empty orbital rank rejected'
  call sparse_pack(a,dense(:0,:))
  if(.not.sparse_valid(a,0,3))error stop 'empty spatial rank rejected'
  call sparse_clear(a);call sparse_clear(b)
  print *, 'PASS exact support storage, tiny values, empty columns, deep copy and invalid inputs'
end program
