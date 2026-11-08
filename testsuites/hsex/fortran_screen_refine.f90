program test_refine
  use mpi
  use hsex_q_response, only: hsex_q_solve
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer::rank,np,ierr,n,g,nc,unit,team,i,j,c,bad,total_bad
  integer::counts(1),status(1)
  complex(8),allocatable::h(:,:),b(:,:),adj(:,:),expected(:,:,:),result(:,:,:,:),u(:,:),psi(:,:,:),phase(:,:,:)
  real(8),allocatable::d(:),e(:),etas(:)
  real(8)::level(1,1),occ(1,1),residual(1),ratio,maximum,global_max
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  open(newunit=unit,file='oracle.bin',access='stream',form='unformatted',status='old')
  read(unit)n,g,nc
  allocate(h(n,n),b(n,g),adj(n,n),d(n),e(n-1),etas(nc),expected(g,g,nc))
  allocate(result(g,g,nc,1),u(n,0),psi(n,1,1),phase(n,g,1))
  read(unit)h,b,adj,d,e,level,etas,expected
  close(unit)
  psi=1d0;phase(:,:,1)=b;counts=g;occ=2d0
  call hsex_q_solve(n,1,1,g,0,nc,12,1,adj,d,e,u,psi,phase,counts,level,occ,etas, &
    result,residual,status,team,1,h)
  bad=0;maximum=0d0
  if(any(status/=0).or.team/=12)bad=1
  do c=1,nc
    do j=1,g
      do i=1,g
        ratio=abs(result(i,j,c,1)-expected(i,j,c))/(2d-12+2d-12*abs(expected(i,j,c)))
        if(.not.ieee_is_finite(ratio))then
          bad=1
        else
          maximum=max(maximum,ratio)
          if(ratio>1d0)bad=1
        endif
      enddo
    enddo
  enddo
  call MPI_Allreduce(bad,total_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(maximum,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(total_bad/=0)then
    print *, 'REFINEMENT_ORACLE_FAILURE',rank,status,team,maximum
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif
  if(rank==0)print *, 'REFINEMENT_ORACLE_PASS mpi/omp/maximum_tolerance_ratio',np,team,global_max
  call MPI_Finalize(ierr)
end program
