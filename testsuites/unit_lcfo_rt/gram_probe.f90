program gram_probe
 use mpi
 use, intrinsic :: ieee_arithmetic, only:ieee_value,ieee_quiet_nan
 use lcfo_gram, only:lcfo_gram_error
 implicit none
 complex(8),allocatable :: c(:,:),full(:,:),local(:,:),total(:,:)
 integer :: rank,np,ierr,n,m,lo,j,k,trial
 real(8) :: got,expected,pi
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 pi=acos(-1d0)
 do trial=1,4
  n=5;if(trial==2)n=129
  m=2*n+1;allocate(full(m,n),local(n,n),total(n,n))
  do j=1,n;do k=1,m
   full(k,j)=exp(cmplx(0d0,2*pi*(k-1)*(j-1)/m,8))/sqrt(dble(m))
  enddo;enddo
  if(trial>=3)full(:,n)=full(:,n)+(.1d0,.2d0)*full(:,1)
  ! Empty rank0, remaining uneven rows.
  lo=rank*m/(np-1)
  if(rank==0)then
   allocate(c(0,n))
  else
   lo=(rank-1)*m/(np-1)
   c=full(lo+1:rank*m/(np-1),:)
  endif
  local=matmul(conjg(transpose(c)),c)
  call MPI_Allreduce(local,total,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  do j=1,n;total(j,j)=total(j,j)-1d0;enddo
  expected=maxval(abs(total))
  if(trial==4.and.rank==np-1)c(1,1)=cmplx(ieee_value(1d0,ieee_quiet_nan),0d0,8)
  call lcfo_gram_error(c,MPI_COMM_WORLD,got)
  if(trial==4)then
   if(got/=huge(1d0))error stop 'nonfinite Gram not rejected'
  else
   if(abs(got-expected)>1d-12)error stop 'packed Gram differs from dense'
  endif
  if(rank==0)write(*,*)'packed Gram trial/dimension/error:',trial,n,got
  deallocate(c,full,local,total)
 enddo
 call MPI_Finalize(ierr)
end program
