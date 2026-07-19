program test_dg_wpw_full_rank_seed
  use mpi
  use dg_wpw_matrix_free_scf,only:initialize_dg_wpw_deterministic_seed
  implicit none
  integer,parameter::nw=24,np=16,nretain=160
  integer::rank,nrank,i,j,ierr,info,lwork
  integer::w_ids(nw),p_ids(np)
  complex(8)::qw(nw,nretain),qp(np,nretain)
  complex(8)::local_gram(nretain,nretain),gram(nretain,nretain)
  complex(8),allocatable::work(:)
  real(8)::eigenvalues(nretain),rwork(3*nretain-2)
  complex(8)::work_query(1)
  external zheev

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=8)error stop 'full-rank seed fixture requires 8 MPI ranks'

  do i=1,nw
    w_ids(i)=rank*nw+i
  enddo
  do i=1,np
    p_ids(i)=100000+rank*np+i
  enddo
  call initialize_dg_wpw_deterministic_seed(w_ids,p_ids,qw,qp,info)
  if(info/=0)error stop 'deterministic seed initialization failed'

  local_gram=matmul(conjg(transpose(qw)),qw)+matmul(conjg(transpose(qp)),qp)
  call MPI_Allreduce(local_gram,gram,size(gram),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'seed Gram reduction failed'
  if(rank==0)then
    call zheev('N','U',nretain,gram,nretain,eigenvalues,work_query,-1,rwork,info)
    if(info/=0)error stop 'seed Gram workspace query failed'
    lwork=max(1,int(real(work_query(1),8)));allocate(work(lwork))
    call zheev('N','U',nretain,gram,nretain,eigenvalues,work,lwork,rwork,info)
    if(info/=0)error stop 'seed Gram diagonalization failed'
    if(eigenvalues(1)<=1d-10*eigenvalues(nretain))then
      write(*,'(a,2es16.8)')'seed metric extrema: ',eigenvalues(1),eigenvalues(nretain)
      error stop 'deterministic production seed is rank deficient at the production cutoff'
    endif
    print '(a)','PASS Si64-sized deterministic WPW seed has full retained metric rank'
  endif
  call MPI_Finalize(ierr)
end program
