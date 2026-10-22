program test
use sym_sub
use sym_vector_sub
implicit none
integer::n,i,j,k
real(8)::v(3)
read(*,*) n
allocate(SymMatA(3,4,n),SymMatB(3,4,n))
do j=1,3
do i=1,3
read(*,*) Amat(i,j)
enddo
enddo
do j=1,3
do i=1,3
read(*,*) Ainv(i,j)
enddo
enddo
do k=1,n
do j=1,4
do i=1,3
read(*,*) SymMatA(i,j,k)
enddo
enddo
enddo
do k=1,n
do j=1,4
do i=1,3
read(*,*) SymMatB(i,j,k)
enddo
enddo
enddo
do i=1,3
read(*,*) v(i)
enddo
call sym_vector_xyz(v)
write(*,'(3es25.16)') v
end program
