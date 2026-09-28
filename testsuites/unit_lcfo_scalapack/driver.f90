program probe
  use mpi
  use lcfo_scalapack
  implicit none
  type(lcfo_dense_state) :: s
  complex(8),allocatable :: a(:,:),v(:,:),work(:),rows(:,:)
  real(8),allocatable :: e(:),rw(:)
  real(8) :: he,oe,re
  integer :: rank,np,err,n,i,j,status,comm,first,count,case_id
  integer,parameter :: sizes(6)=[1,2,3,7,13,19]
  call MPI_Init(err)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
  call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call MPI_Comm_split(MPI_COMM_WORLD,0,np-1-rank,comm,err)
  call MPI_Comm_rank(comm,rank,err)
  do case_id=1,size(sizes)
    n=sizes(case_id)
    allocate(a(n,n),v(n,n),work(2*n),e(n),rw(max(1,3*n-2)),rows(n,n))
    do j=1,n;do i=1,n
      a(i,j)=cmplx(cos(.17d0*(i+j)),sin(.13d0*(i-j)),8)
      if(i==j)a(i,j)=a(i,j)+i
    enddo;enddo
    call lcfo_dense_init(s,n,comm,status)
    if(status/=0)error stop 'init'
    call lcfo_dense_add(s,1,1,a)
    call lcfo_dense_solve(s,n,he,oe,re,status)
    if(status/=0.or.max(he,oe,re)>1d-10)then
      print *, 'FAILED status/diagnostics/values/vectors',status,he,oe,re,s%values,s%vectors
      error stop 'solve diagnostics'
    endif
    v=a
    call zheev('V','U',n,v,n,e,work,size(work),rw,status)
    if(status/=0.or.maxval(abs(e-s%values))>1d-10)error stop 'eigenvalues'
    call lcfo_dense_rows(s,1,n,n,0,rows)
    if(rank==0)then
      do i=1,n
        if(sqrt(sum(abs(matmul(a,rows(:,i))-e(i)*rows(:,i))**2))>1d-9)error stop 'exported residual'
      enddo
      print *, 'PASS distributed LCFO',np,n,he,oe,re
    endif
    first=max(1,n/2);count=n-first+1
    call lcfo_dense_rows(s,first,count,n,np-1,rows(:count,:))
    if(rank==np-1)then
      do j=1,n
        ! Eigenvector phases differ; compare entry magnitudes for the nondegenerate test spectrum.
        if(maxval(abs(abs(rows(:count,j))-abs(v(first:n,j))))>1d-9)error stop 'partial row extraction'
      enddo
    endif
    ! Malformed Hermitian input must fail collectively before diagonalization.
    s%h=0d0
    if(s%myrow==0.and.s%mycol==0)s%h(1,1)=(0d0,1d0)
    call lcfo_dense_solve(s,n,he,oe,re,status)
    if(status==0)error stop 'nonhermitian matrix accepted'
    call lcfo_dense_free(s)
    deallocate(a,v,work,e,rw,rows)
  enddo
  call MPI_Comm_free(comm,err)
  call MPI_Finalize(err)
end program
