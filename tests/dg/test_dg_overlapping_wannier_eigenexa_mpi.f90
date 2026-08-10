program test_dg_overlapping_wannier_eigenexa_mpi
  use mpi
  use structures,only:s_parallel_info
  use eigen_libs_mod
  use dg_overlapping_wannier_solver,only:solve_dg_overlapping_wannier_generalized_eigenexa
  implicit none
  type(s_parallel_info)::info
  integer::comm,rank,nproc,ierr,i,p,nlocal
  integer(8),allocatable::row_ids(:)
  complex(8),allocatable::hrows(:,:),srows(:,:),coeff(:,:)
  real(8)::metric_diagonal(4),target_eigenvalues(4),residual,orthogonality,condition,gamma_defect
  real(8),allocatable::eigenvalues(:)
  integer(8)::workspace,signature
  logical::ok
  character(256)::message
  character(32)::case_name
  real(8)::solve_tolerance

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
  call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(4,info%nrow_local,info%ncol_local)
  info%flag_eigenexa_init=.true.
  case_name='normal';if(command_argument_count()>=1)call get_command_argument(1,case_name)
  solve_tolerance=1d-10
  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(row_ids(nlocal),hrows(nlocal,4),srows(nlocal,4),coeff(4,3),eigenvalues(3))
  metric_diagonal=[2d0,1.5d0,1.2d0,0.9d0];target_eigenvalues=[0.2d0,0.6d0,1.1d0,2d0]
  hrows=(0d0,0d0);srows=(0d0,0d0);p=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;row_ids(p)=i;srows(p,i)=metric_diagonal(i)
    hrows(p,i)=metric_diagonal(i)*target_eigenvalues(i)
  enddo
  select case(trim(case_name))
  case('degenerate')
    target_eigenvalues=[0d0,1d0,2d0,2d0]
    solve_tolerance=1d-6
    do p=1,nlocal;hrows(p,:)=0d0;hrows(p,int(row_ids(p)))=&
      metric_diagonal(int(row_ids(p)))*target_eigenvalues(int(row_ids(p)));enddo
  case('nonreal')
    if(nlocal>0.and.row_ids(1)==1_8)hrows(1,1)=hrows(1,1)+cmplx(0d0,1d-4,8)
  case('illmetric')
    do p=1,nlocal;if(row_ids(p)==4_8)srows(p,4)=1d-16;enddo
  case('residual')
    solve_tolerance=1d-30
    do p=1,nlocal
      if(row_ids(p)==1_8)hrows(p,2)=0.123456789d0
      if(row_ids(p)==2_8)hrows(p,1)=0.123456789d0
    enddo
  end select
  call solve_dg_overlapping_wannier_generalized_eigenexa(info,comm,row_ids,hrows,srows,3,&
    solve_tolerance,1d-12,1d-12,coeff,eigenvalues,residual,orthogonality,condition,gamma_defect,&
    workspace,ok,message)
  if(trim(case_name)/='normal')then
    call require(.not.ok,'negative generalized EigenExa case must reject')
    if(rank==0)write(*,'(3a,i0)')'REJECT ',trim(case_name),' ranks=',nproc
    call eigen_free();call MPI_Finalize(ierr);stop
  endif
  call require(ok,trim(message));call require(maxval(abs(eigenvalues-target_eigenvalues(1:3)))<1d-10,&
    'known generalized EigenExa spectrum')
  call require(residual<1d-10.and.orthogonality<1d-10,'generalized EigenExa quality receipts')
  call require(workspace>0_8.and.gamma_defect==0d0,'measured workspace and Gamma-real receipts')
  signature=nint(sum(eigenvalues*[1d0,3d0,7d0])*1d12,8)
  if(rank==0)write(*,'(a,i0,a,i0)')'EIGENEXA ranks=',nproc,' signature=',signature
  call eigen_free();call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure,error
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,error)
    if(global_failure/=0)error stop label
  end subroutine
end program
