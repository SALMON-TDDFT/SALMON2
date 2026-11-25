! Diagnostic only: complete polynomials expose radial terms missing from solid harmonics.
! No change to molecular Hamiltonian or production field representation.
subroutine sp_polynomial_probe(n, nv, degree, xyz, val, query, predicted, rank, status) bind(C)
 use iso_c_binding
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer(c_int), value :: n,nv,degree
 real(c_double), intent(in) :: xyz(3,n),val(n),query(3,nv)
 real(c_double), intent(out) :: predicted(nv)
 integer(c_int), intent(out) :: rank,status
 real(c_double), allocatable :: a(:,:),b(:,:),s(:),work(:)
 real(c_double) :: item,worksize(1)
 integer :: p,i,j,k,l,m,info,lwork,ast
 rank=0;status=1
 if(n<1.or.nv<1.or.degree<0.or.degree>8) return
 p=(degree+1)*(degree+2)*(degree+3)/6
 if(n<p) return
 allocate(a(n,p),b(max(n,p),1),s(min(n,p)),stat=ast)
 if(ast/=0) return
 do i=1,n
  if(.not.ieee_is_finite(val(i))) return
  do j=1,3
   if(.not.ieee_is_finite(xyz(j,i))) return
  enddo
 enddo
 do i=1,nv
  do j=1,3
   if(.not.ieee_is_finite(query(j,i))) return
  enddo
 enddo
 call matrix()
 call dgelss(n,p,1,a,n,b,max(n,p),s,1d-12,rank,worksize,-1,info)
 if(info/=0) return
 lwork=ceiling(worksize(1));allocate(work(lwork),stat=ast)
 if(ast/=0) return
 call matrix()
 call dgelss(n,p,1,a,n,b,max(n,p),s,1d-12,rank,work,lwork,info)
 if(info/=0.or.rank/=p) return
 do i=1,nv
  item=0d0;m=0
  do j=0,degree
   do k=0,degree-j
    do l=0,degree-j-k
     m=m+1
     item=item+b(m,1)*query(1,i)**j*query(2,i)**k*query(3,i)**l
    enddo
   enddo
  enddo
  if(.not.ieee_is_finite(item)) return
  predicted(i)=item
 enddo
 status=0
contains
 subroutine matrix()
  m=0
  do j=0,degree
   do k=0,degree-j
    do l=0,degree-j-k
     m=m+1
     do i=1,n
      a(i,m)=xyz(1,i)**j*xyz(2,i)**k*xyz(3,i)**l
     enddo
    enddo
   enddo
  enddo
  b=0d0;b(1:n,1)=val
 end subroutine
end subroutine
