program test
use sym_sub
use sym_spin_density_sub
implicit none
integer::n,i,j,k
real(8)::v(3),re,im
integer::g(3),l,h
complex(8),allocatable::dm(:,:,:,:,:)
read(*,*) n
do i=1,3
read(*,*) g(i)
enddo
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
allocate(dm(g(1),g(2),g(3),2,2))
do h=1,2
do l=1,2
do k=1,g(3)
do j=1,g(2)
do i=1,g(1)
read(*,*) re
read(*,*) im
dm(i,j,k,l,h)=cmplx(re,im,8)
enddo
enddo
enddo
enddo
enddo
call sym_spin_density_grid(dm)
do h=1,2
do l=1,2
do k=1,g(3)
do j=1,g(2)
do i=1,g(1)
write(*,'(2es25.16)') real(dm(i,j,k,l,h)),aimag(dm(i,j,k,l,h))
enddo
enddo
enddo
enddo
enddo
end program
