program link_probe
 use mpi
 use iso_fortran_env,only:int64
 use lcfo_mlwf_links
 implicit none
 complex(8),allocatable :: grid(:,:),shifted(:,:),raw(:,:,:,:),local(:,:),dense(:,:)
 real(8),allocatable :: position(:,:)
 real(8) :: lengths(3),dv,delta
 integer :: rank,np,ierr,m,n,i,j,a,width
 integer(int64) :: scratch
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 n=129;m=3+rank;if(rank==0)m=0
 allocate(grid(m,n),position(3,m),shifted(m,n),local(n,n),dense(n,n))
 lengths=[11d0,13d0,17d0];dv=.037d0
 do i=1,m
  position(:,i)=[.3d0*i+.2d0*rank,.7d0*i,.9d0*rank+.4d0*i]
  do j=1,n
   grid(i,j)=cmplx(sin(.13d0*(i+j+rank)),cos(.19d0*(i-j+2*rank)),8)
  enddo
 enddo
 do width=1,64,21
  call lcfo_initial_links(grid,position,lengths,dv,MPI_COMM_WORLD,raw,width,scratch)
  if(rank==0)then
   if(any(shape(raw)/=[n,n,6,1]))error stop 'root shape'
  else
   if(size(raw)/=0)error stop 'replicated raw storage'
  endif
  if(scratch/=int((m+2*n)*min(n,width)+m,int64))error stop 'unbounded scratch'
  do a=1,3
   delta=2*acos(-1d0)/lengths(a)
   do j=1,n
    shifted(:,j)=grid(:,j)*exp(cmplx(0d0,-delta*position(a,:),8))
   enddo
   local=matmul(conjg(transpose(grid)),shifted)*dv
   call MPI_Allreduce(local,dense,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
   if(rank==0)then
    if(maxval(abs(raw(:,:,a,1)-dense))>1d-12)error stop 'dense link mismatch'
    if(maxval(abs(raw(:,:,a+3,1)-conjg(transpose(dense))))>1d-12)error stop 'adjoint mismatch'
   endif
  enddo
 enddo
 if(rank==0)print *,'Root-only tiled initial links passed'
 call MPI_Finalize(ierr)
end program
