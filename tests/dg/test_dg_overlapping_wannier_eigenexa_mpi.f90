program test_dg_overlapping_wannier_eigenexa_mpi
  use mpi
  use structures,only:s_parallel_info
  use eigen_libs_mod
  use dg_overlapping_wannier_solver,only:solve_dg_overlapping_wannier_generalized_eigenexa
  use dg_overlapping_wannier_construction,only:build_dg_group_averaged_occupied_candidates_eigenexa,&
    build_dg_cocycle_averaged_occupied_candidates_eigenexa,measure_dg_rank_fixed_symmetry_residuals
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
  case_name='normal';if(command_argument_count()>=1)call get_command_argument(1,case_name)
  if(index(trim(case_name),'average')==1)then;call run_average_case();call MPI_Finalize(ierr);stop;endif
  if(trim(case_name)=='cocycle')then;call run_cocycle_case();call MPI_Finalize(ierr);stop;endif
  call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
  call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(4,info%nrow_local,info%ncol_local)
  info%flag_eigenexa_init=.true.
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
  subroutine run_average_case()
    complex(8),allocatable::average_occupied(:,:),average_candidates(:,:)
    real(8),allocatable::average_spectrum(:)
    real(8),allocatable::average_weights(:)
    integer(8),allocatable::average_maps(:,:)
    integer::average_product(2,2),ii,global_point,global_count,average_rank,average_requested
    real(8)::average_trace,average_closure,average_gamma,average_selected_edge,&
      average_rejected_edge,average_cluster_gap
    integer(8)::average_workspace,average_signature
    logical::average_ok
    character(256)::average_message
    allocate(average_occupied(1,2),average_weights(2),average_maps(2,2))
    average_occupied=(0d0,0d0);average_weights=1d0;global_count=2*nproc
    do ii=1,2
      global_point=2*rank+ii
      average_maps(ii,1)=global_point
      average_maps(ii,2)=modulo(global_point-1+global_count/2,global_count)+1
      if(global_point==1)average_occupied(1,ii)=1d0
    enddo
    if(trim(case_name)=='average_unique')average_occupied=1d0/sqrt(real(global_count,8))
    if(trim(case_name)=='average_nonorthogonal')average_occupied=2d0/sqrt(real(global_count,8))
    average_product=reshape([1,2,2,1],[2,2])
    average_requested=merge(2,1,trim(case_name)=='average')
    call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
    call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(2,info%nrow_local,info%ncol_local)
    info%flag_eigenexa_init=.true.
    call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,average_occupied,&
      average_weights,average_maps,average_product,1,average_requested,1d-12,average_candidates,average_spectrum,&
      average_rank,average_trace,average_closure,average_gamma,average_workspace,average_ok,average_message,&
      average_selected_edge,average_rejected_edge,average_cluster_gap)
    if(trim(case_name)=='average_split')then
      call require(.not.average_ok,'group-average selection must reject a split degenerate block')
      if(rank==0)write(*,'(a,i0)')'REJECT average_split ranks=',nproc
      call eigen_free();return
    endif
    if(trim(case_name)=='average_nonorthogonal')then
      call require(.not.average_ok,'group-average input occupied space must be metric orthonormal')
      if(rank==0)write(*,'(a,i0)')'REJECT average_nonorthogonal ranks=',nproc
      call eigen_free();return
    endif
    if(trim(case_name)=='average_unique')then
      call require(average_ok.and.average_rank==1.and.abs(average_spectrum(1)-1d0)<1d-12.and.&
        abs(average_trace-1d0)<1d-12.and.average_closure<1d-12.and.&
        abs(average_selected_edge-1d0)<1d-12.and.abs(average_rejected_edge)<1d-12.and.&
        abs(average_cluster_gap-1d0)<1d-12,&
        'unique invariant group-average rank is accepted')
      average_signature=nint(average_spectrum(1)*1d12,8)
      if(rank==0)write(*,'(a,i0,a,i0)')'AVERAGE_UNIQUE ranks=',nproc,' signature=',average_signature
      call eigen_free();return
    endif
    call require(average_ok,trim(average_message))
    call require(average_rank==2.and.maxval(abs(average_spectrum-[0.5d0,0.5d0]))<1d-12,&
      'distributed group-average spectrum')
    call require(abs(average_trace-1d0)<1d-12.and.average_closure<1d-12.and.&
      average_gamma==0d0.and.average_workspace>0_8,'distributed group-average receipts')
    average_signature=nint(sum(average_spectrum*[1d0,3d0])*1d12,8)
    if(rank==0)write(*,'(a,i0,a,i0)')'AVERAGE ranks=',nproc,' signature=',average_signature
    if(rank==0)average_maps(:,2)=1_8
    call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,average_occupied,&
      average_weights,average_maps,average_product,1,2,1d-12,average_candidates,average_spectrum,&
      average_rank,average_trace,average_closure,average_gamma,average_workspace,average_ok,average_message)
    call require(.not.average_ok,'distributed group-average rejects a non-group point action')
    call eigen_free()
  end subroutine

  subroutine run_cocycle_case()
    complex(8),allocatable::occupied(:,:),candidates(:,:)
    real(8),allocatable::weights(:),candidate_spectrum(:)
    integer(8),allocatable::translation_maps(:,:),representative_maps(:,:)
    integer(8),allocatable::full_maps(:,:)
    real(8),allocatable::full_total(:),full_boundary(:),full_interior(:)
    logical,allocatable::no_boundary(:)
    integer::point_product(2,2),cocycle(2,2),ii,global_point,global_count,candidate_rank
    real(8)::projector_trace,closure,gamma
    integer(8)::workspace,cocycle_signature
    logical::cocycle_ok
    character(256)::cocycle_message

    global_count=4*nproc
    allocate(occupied(1,4),weights(4),translation_maps(4,2),representative_maps(4,2),&
      full_maps(4,4),full_total(4),full_boundary(4),full_interior(4),no_boundary(4))
    occupied=(0d0,0d0);weights=1d0
    do ii=1,4
      global_point=4*rank+ii
      translation_maps(ii,1)=global_point
      translation_maps(ii,2)=modulo(global_point-1+global_count/2,global_count)+1
      representative_maps(ii,1)=global_point
      representative_maps(ii,2)=modulo(global_point-1+global_count/4,global_count)+1
      do p=1,4
        full_maps(ii,p)=modulo(global_point-1+(p-1)*global_count/4,global_count)+1
      enddo
      if(global_point==1.or.global_point==1+global_count/2)occupied(1,ii)=1d0/sqrt(2d0)
    enddo
    point_product=reshape([1,2,2,1],[2,2])
    cocycle=reshape([1,1,1,2],[2,2])
    call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
    call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(2,info%nrow_local,info%ncol_local)
    info%flag_eigenexa_init=.true.
    call build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      translation_maps,representative_maps,point_product,cocycle,1,2,1d-12,candidates,&
      candidate_spectrum,candidate_rank,projector_trace,closure,gamma,workspace,cocycle_ok,cocycle_message)
    call require(cocycle_ok,trim(cocycle_message))
    call require(candidate_rank==2.and.maxval(abs(candidate_spectrum-[0.5d0,0.5d0]))<1d-12,&
      'cocycle representative average matches the explicit affine orbit spectrum')
    call require(abs(projector_trace-1d0)<1d-12.and.closure<1d-12.and.gamma==0d0.and.workspace>0_8,&
      'cocycle representative average receipts')
    no_boundary=.false.
    call measure_dg_rank_fixed_symmetry_residuals(comm,candidates,weights,full_maps,no_boundary,&
      total_residual=full_total,boundary_residual=full_boundary,interior_residual=full_interior,&
      ok=cocycle_ok,message=cocycle_message)
    call require(cocycle_ok.and.maxval(full_total)<1d-12,&
      'cocycle representative average equals the explicit full-affine projector')
    cocycle_signature=nint(sum(candidate_spectrum*[1d0,3d0])*1d12,8)
    if(rank==0)write(*,'(a,i0,a,i0)')'COCYCLE ranks=',nproc,' signature=',cocycle_signature
    cocycle(2,2)=1
    call build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      translation_maps,representative_maps,point_product,cocycle,1,2,1d-12,candidates,&
      candidate_spectrum,candidate_rank,projector_trace,closure,gamma,workspace,cocycle_ok,cocycle_message)
    call require(.not.cocycle_ok,'corrupt representative cocycle must reject')
    call eigen_free()
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure,error
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,error)
    if(global_failure/=0)error stop label
  end subroutine
end program
