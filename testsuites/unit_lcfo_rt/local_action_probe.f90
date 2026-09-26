program local_action_probe
 use mpi
 use lcfo_ace_local,only:lcfo_ace_local_action
 use hse_ace,only:hse_ace_state,hse_ace_apply
 implicit none
 type(hse_ace_state) :: ace
 integer :: rank,nproc,ierr,n,lo,hi,i,j,k,nf
 complex(8),allocatable :: b(:,:),psi(:,:),h(:,:),old(:,:),expected(:,:),c(:,:),global(:,:),w(:,:,:)
 real(8) :: dv,pi,error,allerror
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
 if(nproc/=2)error stop 'run with 2 ranks'
 n=rank+2;lo=1+2*rank;hi=lo+n-1;dv=.5d0;pi=acos(-1d0)
 allocate(b(8,n),psi(8,2),h(8,2),old(8,2),expected(8,2),c(5,2),global(5,2),w(5,2,1))
 do j=1,n;do i=1,8
 b(i,j)=exp(cmplx(0d0,2*pi*(i-1)*(j-1)/8,8))/sqrt(8*dv)
 enddo;enddo
 do j=1,2;do i=1,8
 psi(i,j)=cmplx(sin(real(i+2*j+rank,8)),cos(real(2*i-j+rank,8)),8)
 old(i,j)=cmplx(cos(real(i*j+rank,8)),sin(real(i-j+2*rank,8)),8)
 enddo;enddo
 c=0;c(lo:hi,:)=matmul(conjg(transpose(b)),psi)*dv
 call MPI_Allreduce(c,global,10,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 do k=1,3
 nf=2*k;allocate(ace%factors(5,nf,1));ace%dv=.7d0
 do j=1,nf;do i=1,5
 ace%factors(i,j,1)=cmplx(sin(real(2*i+j,8)),cos(real(i+3*j,8)),8)/3
 enddo;enddo
 call hse_ace_apply(ace,reshape(global,[5,2,1]),w,ierr)
 if(ierr/=0)error stop 'dense ACE reference'
 expected=old+matmul(b,w(lo:hi,:,1))
 expected=matmul(b,matmul(conjg(transpose(b)),expected))*dv
 h=old
 call lcfo_ace_local_action(b,psi,h,ace%factors(lo:hi,:,1),dv,ace%dv,MPI_COMM_WORLD)
 error=maxval(abs(h-expected))
 call MPI_Allreduce(error,allerror,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(allerror>2d-12)error stop 'distributed/fused ACE differs from dense projection'
 if(rank==0)write(*,*)'unequal partitions, ACE rank, max error:',nf,allerror
 deallocate(ace%factors)
 enddo
 call MPI_Finalize(ierr)
end program
