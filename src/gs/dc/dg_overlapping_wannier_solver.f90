#include "config.h"
module dg_overlapping_wannier_solver
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
#ifdef USE_EIGENEXA
  use structures,only:s_parallel_info
  use eigen_eigenexa,only:eigen_pdsyevd_ex_distributed_blocks
  use eigen_libs_mod,only:eigen_owner_node,eigen_translate_g2l,eigen_translate_l2g,&
    eigen_loop_start,eigen_loop_end
#endif
  implicit none
  private
  public::solve_dg_overlapping_wannier_coefficients,select_dg_symmetry_complete_occupied_states
#ifdef USE_EIGENEXA
  public::solve_dg_overlapping_wannier_generalized_eigenexa
#endif
contains
  subroutine select_dg_symmetry_complete_occupied_states(eigenvalues,generator_representation,target_count,&
      degeneracy_tolerance,gamma_real_tolerance,coefficients,canonical_coefficients,selected_count,&
      gamma_real_defect,ok,message)
    real(8),intent(in)::eigenvalues(:),degeneracy_tolerance,gamma_real_tolerance
    integer,intent(in)::target_count
    complex(8),intent(in)::generator_representation(:,:,:),coefficients(:,:)
    complex(8),intent(out)::canonical_coefficients(:,:)
    integer,intent(out)::selected_count
    real(8),intent(out)::gamma_real_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,g,pivot,n
    real(8)::scale,imaginary_defect,spectral_scale
    complex(8)::phase
    complex(8),allocatable::identity(:,:),energy_commutator(:,:)

    ok=.false.;message='';selected_count=0;gamma_real_defect=huge(1d0)
    canonical_coefficients=(0d0,0d0);n=size(eigenvalues)
    if(n<1.or.size(coefficients,1)/=n.or.size(coefficients,2)/=n.or.&
        any(shape(canonical_coefficients)/=shape(coefficients)).or.&
        size(generator_representation,1)/=n.or.size(generator_representation,2)/=n.or.&
        size(generator_representation,3)<1.or.target_count<1.or.&
        target_count>size(eigenvalues).or.degeneracy_tolerance<=0d0.or.&
        gamma_real_tolerance<=0d0.or..not.all(ieee_is_finite(eigenvalues)).or.&
        .not.finite_matrix(coefficients).or..not.finite_matrix(&
        reshape(generator_representation,[n,n*size(generator_representation,3)])))then
      message='invalid symmetry-complete occupied selection contract';return
    endif
    allocate(identity(n,n),energy_commutator(n,n));identity=(0d0,0d0)
    do i=1,n;identity(i,i)=1d0;enddo
    spectral_scale=max(1d0,maxval(abs(eigenvalues)))
    do g=1,size(generator_representation,3)
      if(maxval(abs(matmul(conjg(transpose(generator_representation(:,:,g))),&
          generator_representation(:,:,g))-identity))>degeneracy_tolerance)then
        message='occupied symmetry generator is not unitary';return
      endif
      do j=1,n;do i=1,n
        energy_commutator(i,j)=(eigenvalues(i)-eigenvalues(j))*generator_representation(i,j,g)
      enddo;enddo
      if(maxval(abs(energy_commutator))>degeneracy_tolerance*spectral_scale)then
        message='occupied symmetry generator does not preserve generalized eigenspaces';return
      endif
      if(target_count<n)then
        if(maxval(abs(generator_representation(1:target_count,target_count+1:n,g)))>&
            degeneracy_tolerance.or.maxval(abs(generator_representation(target_count+1:n,&
            1:target_count,g)))>degeneracy_tolerance)then
          message='occupied boundary splits a generator-connected symmetry block';return
        endif
      endif
    enddo
    do i=2,size(eigenvalues)
      if(eigenvalues(i)<eigenvalues(i-1)-degeneracy_tolerance)then
        message='generalized eigenvalues are not ordered';return
      endif
    enddo
    if(target_count<size(eigenvalues))then
      if(abs(eigenvalues(target_count+1)-eigenvalues(target_count))<=degeneracy_tolerance)then
        message='occupied boundary splits a degenerate cluster';return
      endif
    endif
    canonical_coefficients=coefficients;gamma_real_defect=0d0
    do j=1,target_count
      pivot=maxloc(abs(coefficients(:,j)),dim=1)
      scale=maxval(abs(coefficients(:,j)))
      if(scale<=0d0)then;message='zero occupied generalized eigenvector';return;endif
      phase=conjg(coefficients(pivot,j))/abs(coefficients(pivot,j))
      imaginary_defect=maxval(abs(aimag(phase*coefficients(:,j))))/scale
      gamma_real_defect=max(gamma_real_defect,imaginary_defect)
      if(imaginary_defect>gamma_real_tolerance)then
        message='non-real Gamma occupied coefficient';return
      endif
      canonical_coefficients(:,j)=phase*coefficients(:,j)
      canonical_coefficients(:,j)=cmplx(real(canonical_coefficients(:,j),8),0d0,8)
    enddo
    selected_count=target_count;ok=.true.
  end subroutine

#if defined(USE_MPI) && defined(USE_EIGENEXA)
  subroutine solve_dg_overlapping_wannier_generalized_eigenexa(info,comm,row_ids,hrows,srows,nstate,&
      tolerance,metric_tolerance,gamma_real_tolerance,coefficients,eigenvalues,maximum_residual,&
      orthogonality_defect,metric_condition,gamma_real_defect,workspace_peak_bytes,ok,message)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::comm,nstate
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    real(8),intent(in)::tolerance,metric_tolerance,gamma_real_tolerance
    complex(8),intent(out)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:),maximum_residual,orthogonality_defect,metric_condition,&
      gamma_real_defect
    integer(8),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,parameter::column_tile=32
    real(8),allocatable::cyclic_h(:,:),cyclic_vectors(:,:),all_values(:),&
      local_metric_vectors(:,:),eigenvector_row(:),eigenvector_column(:),&
      h_times_columns(:,:),local_transformed(:,:),global_transformed(:,:),transformed_rows(:,:)
    real(8)::matrix_scale,imaginary_scale,local_scale,global_scale,local_trace,global_trace,&
      pivot,minimum_pivot,maximum_pivot
    integer,allocatable::ownership(:)
    integer::n,nlocal,i,j,k,grow,gcol,lrow,lcol,first,count,ierr,rank,pivot_owner,pivot_local
    logical::eigen_ok
    character(256)::detail

    ok=.false.;message='';workspace_peak_bytes=0_8;maximum_residual=huge(1d0)
    orthogonality_defect=huge(1d0);metric_condition=huge(1d0);gamma_real_defect=huge(1d0)
    coefficients=(0d0,0d0);eigenvalues=0d0;n=size(hrows,2);nlocal=size(row_ids)
    if(.not.info%flag_eigenexa_init.or.n<1.or.nstate<1.or.nstate>n.or.&
        size(hrows,1)/=nlocal.or.any(shape(srows)/=shape(hrows)).or.&
        any(shape(coefficients)/=[n,nstate]).or.size(eigenvalues)/=nstate.or.&
        any(row_ids<1_8).or.any(row_ids>int(n,8)).or.tolerance<=0d0.or.&
        metric_tolerance<=0d0.or.gamma_real_tolerance<=0d0.or.&
        .not.finite_matrix(hrows).or..not.finite_matrix(srows))then
      message='invalid distributed EigenExa generalized solve contract';return
    endif
    matrix_scale=max(1d0,maxval(abs(srows)))
    imaginary_scale=max(maxval(abs(aimag(hrows))),maxval(abs(aimag(srows))))
    call MPI_Allreduce(MPI_IN_PLACE,matrix_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,imaginary_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    imaginary_scale=imaginary_scale/matrix_scale
    if(imaginary_scale>gamma_real_tolerance)then
      message='non-real Gamma stitched generalized pencil';return
    endif
    call distributed_hermiticity(comm,row_ids,hrows,tolerance,i)
    call distributed_hermiticity(comm,row_ids,srows,metric_tolerance,j)
    if(i/=0.or.j/=0)then;message='non-Hermitian distributed EigenExa generalized pencil';return;endif
    allocate(ownership(n));ownership=0
    do i=1,nlocal;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      message='EigenExa coefficient rows do not form a unique global partition';return
    endif
    call MPI_Comm_rank(comm,rank,ierr)
    allocate(cyclic_h(info%nrow_local,info%ncol_local),&
      cyclic_vectors(info%nrow_local,info%ncol_local),all_values(n))
    cyclic_h=0d0
    allocate(eigenvector_row(n),eigenvector_column(n))
    allocate(local_metric_vectors(nlocal,n));local_metric_vectors=0d0
    minimum_pivot=huge(1d0);maximum_pivot=0d0
    do k=1,n
      pivot_owner=-1;pivot_local=0
      do i=1,nlocal
        if(int(row_ids(i))==k)then;pivot_owner=rank;pivot_local=i;endif
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,pivot_owner,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      eigenvector_row=0d0;pivot=-huge(1d0)
      if(rank==pivot_owner)then
        if(k>1)eigenvector_row(1:k-1)=local_metric_vectors(pivot_local,1:k-1)
        pivot=real(srows(pivot_local,k),8)-sum(eigenvector_row(1:k-1)**2)
        if(pivot>0d0)then
          pivot=sqrt(pivot);eigenvector_row(k)=pivot;local_metric_vectors(pivot_local,k)=pivot
        endif
      endif
      call MPI_Bcast(pivot,1,MPI_DOUBLE_PRECISION,pivot_owner,comm,ierr)
      call MPI_Bcast(eigenvector_row,n,MPI_DOUBLE_PRECISION,pivot_owner,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.pivot<=sqrt(metric_tolerance*matrix_scale).or..not.ieee_is_finite(pivot))then
        message='distributed overlap Cholesky rank or conditioning gate failed';return
      endif
      minimum_pivot=min(minimum_pivot,pivot);maximum_pivot=max(maximum_pivot,pivot)
      do i=1,nlocal
        if(row_ids(i)<=int(k,8))cycle
        local_metric_vectors(i,k)=(real(srows(i,k),8)-&
          sum(local_metric_vectors(i,1:k-1)*eigenvector_row(1:k-1)))/pivot
      enddo
    enddo
    metric_condition=(maximum_pivot/minimum_pivot)**2
    allocate(h_times_columns(nlocal,n),transformed_rows(nlocal,n));h_times_columns=0d0;transformed_rows=0d0
    ! Y=L^{-1}H.
    do first=1,n,column_tile
      count=min(column_tile,n-first+1)
      allocate(global_transformed(n,count));global_transformed=0d0
      do i=1,nlocal;global_transformed(int(row_ids(i)),:)=real(hrows(i,first:first+count-1),8);enddo
      call MPI_Allreduce(MPI_IN_PLACE,global_transformed,n*count,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      call solve_lower_tile(global_transformed,count,ok,detail)
      if(.not.ok)then;message=trim(detail);return;endif
      do i=1,nlocal;h_times_columns(i,first:first+count-1)=global_transformed(int(row_ids(i)),:);enddo
      deallocate(global_transformed)
    enddo
    ! B=Y L^{-T}; solve L B^T=Y^T.
    do first=1,n,column_tile
      count=min(column_tile,n-first+1)
      allocate(global_transformed(n,count));global_transformed=0d0
      do j=1,count
        grow=first+j-1;eigenvector_row=0d0
        do i=1,nlocal
          if(int(row_ids(i))==grow)eigenvector_row=h_times_columns(i,:)
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,eigenvector_row,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        global_transformed(:,j)=eigenvector_row
      enddo
      call solve_lower_tile(global_transformed,count,ok,detail)
      if(.not.ok)then;message=trim(detail);return;endif
      do i=1,nlocal;transformed_rows(i,first:first+count-1)=global_transformed(int(row_ids(i)),:);enddo
      deallocate(global_transformed)
    enddo
    cyclic_h=0d0
    do grow=1,n
      eigenvector_row=0d0
      do i=1,nlocal
        if(int(row_ids(i))==grow)eigenvector_row=transformed_rows(i,:)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,eigenvector_row,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(eigen_owner_node(grow,info%nprow,info%myrow)/=info%myrow)cycle
      lrow=eigen_translate_g2l(grow,info%nprow,info%myrow)
      do gcol=1,n
        if(eigen_owner_node(gcol,info%npcol,info%mycol)/=info%mycol)cycle
        lcol=eigen_translate_g2l(gcol,info%npcol,info%mycol)
        cyclic_h(lrow,lcol)=eigenvector_row(gcol)
      enddo
    enddo
    local_scale=maxval(abs(cyclic_h))
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(global_scale<=0d0)then;message='EigenExa transformed Hamiltonian is zero';return;endif
    local_trace=0d0
    do grow=1,n
      if(eigen_owner_node(grow,info%nprow,info%myrow)/=info%myrow.or.&
          eigen_owner_node(grow,info%npcol,info%mycol)/=info%mycol)cycle
      lrow=eigen_translate_g2l(grow,info%nprow,info%myrow)
      lcol=eigen_translate_g2l(grow,info%npcol,info%mycol)
      local_trace=local_trace+abs(cyclic_h(lrow,lcol))
    enddo
    call MPI_Allreduce(local_trace,global_trace,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(global_trace<=0d0)then;message='EigenExa transformed Hamiltonian active diagonal is zero';return;endif
    ok=.false.
    call eigen_pdsyevd_ex_distributed_blocks(info,n,cyclic_h,all_values,cyclic_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;message='EigenExa transformed-H diagonalization: '//trim(detail);return;endif
    if(nstate<n)then
      if(abs(all_values(nstate+1)-all_values(nstate))<=tolerance*max(1d0,maxval(abs(all_values))))then
        message='occupied boundary splits a symmetry-degenerate generalized eigenspace';return
      endif
    endif
    allocate(global_transformed(n,nstate))
    do j=1,nstate
      call gather_cyclic_column(j,cyclic_vectors,eigenvector_column,ok,detail)
      if(.not.ok)then;message=trim(detail);return;endif
      global_transformed(:,j)=eigenvector_column
    enddo
    coefficients=(0d0,0d0)
    do grow=n,1,-1
      allocate(local_transformed(1,nstate));local_transformed=0d0
      do i=1,nlocal
        if(row_ids(i)<=int(grow,8))cycle
        local_transformed(1,:)=local_transformed(1,:)+&
          local_metric_vectors(i,grow)*real(coefficients(int(row_ids(i)),:),8)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,local_transformed,nstate,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      pivot=0d0
      do i=1,nlocal
        if(int(row_ids(i))==grow)pivot=local_metric_vectors(i,grow)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,pivot,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      coefficients(grow,:)=cmplx((global_transformed(grow,:)-local_transformed(1,:))/pivot,0d0,8)
      deallocate(local_transformed)
    enddo
    eigenvalues=all_values(1:nstate);gamma_real_defect=0d0
    call canonicalize_real_columns(coefficients)
    call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
      maximum_residual,orthogonality_defect)
    if(maximum_residual>10d0*tolerance.or.orthogonality_defect>10d0*tolerance)then
      ok=.false.
      message='EigenExa generalized residual or S-orthonormality gate failed';return
    endif
    workspace_peak_bytes=8_8*int(size(cyclic_h)+size(cyclic_vectors)+size(all_values)+&
      size(local_metric_vectors)+size(h_times_columns)+size(transformed_rows)+n*(2+column_tile),8)+&
      16_8*int(size(coefficients),8)
    ok=.true.;message=''
  contains
    subroutine solve_lower_tile(rhs,tile_count,tile_ok,tile_message)
      real(8),intent(inout)::rhs(:,:)
      integer,intent(in)::tile_count
      logical,intent(out)::tile_ok
      character(*),intent(out)::tile_message
      integer::row,ii,error
      real(8)::diagonal
      tile_ok=.false.;tile_message=''
      if(size(rhs,1)/=n.or.size(rhs,2)/=tile_count)then
        tile_message='invalid distributed lower-triangular tile';return
      endif
      do row=1,n
        eigenvector_row=0d0;diagonal=0d0
        do ii=1,nlocal
          if(int(row_ids(ii))/=row)cycle
          eigenvector_row=local_metric_vectors(ii,:);diagonal=local_metric_vectors(ii,row)
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,eigenvector_row,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,error)
        call MPI_Allreduce(MPI_IN_PLACE,diagonal,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,error)
        if(error/=MPI_SUCCESS.or.diagonal<=0d0)then
          tile_message='distributed lower-triangular solve failed';return
        endif
        if(row>1)rhs(row,:)=(rhs(row,:)-matmul(eigenvector_row(1:row-1),rhs(1:row-1,:)))/diagonal
        if(row==1)rhs(row,:)=rhs(row,:)/diagonal
      enddo
      tile_ok=.true.
    end subroutine

    subroutine gather_cyclic_column(column,cyclic,values,column_ok,column_message)
      integer,intent(in)::column
      real(8),intent(in)::cyclic(:,:)
      real(8),intent(out)::values(:)
      logical,intent(out)::column_ok
      character(*),intent(out)::column_message
      integer::row,global_row,global_column,lr,lc,error,row_start,row_end,column_start,column_end
      values=0d0
      row_start=eigen_loop_start(1,info%nprow,info%myrow)
      row_end=eigen_loop_end(n,info%nprow,info%myrow)
      column_start=eigen_loop_start(1,info%npcol,info%mycol)
      column_end=eigen_loop_end(n,info%npcol,info%mycol)
      do lc=column_start,column_end
        global_column=eigen_translate_l2g(lc,info%npcol,info%mycol)
        if(global_column/=column)cycle
        do lr=row_start,row_end
          global_row=eigen_translate_l2g(lr,info%nprow,info%myrow)
          values(global_row)=cyclic(lr,lc)
        enddo
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,values,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,error)
      column_ok=error==MPI_SUCCESS
      if(column_ok)then;column_message='';else;column_message='cyclic column gather failed';endif
    end subroutine
    subroutine canonicalize_real_columns(matrix)
      complex(8),intent(inout)::matrix(:,:)
      integer::column,p
      do column=1,size(matrix,2)
        p=maxloc(abs(matrix(:,column)),dim=1)
        if(real(matrix(p,column),8)<0d0)matrix(:,column)=-matrix(:,column)
      enddo
    end subroutine
  end subroutine
#endif

  subroutine solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,nstate,max_iterations,&
      tolerance,metric_tolerance,coefficients,eigenvalues,maximum_residual,orthogonality_defect,&
      metric_condition,ok,message,initial_coefficients)
    integer,intent(in)::comm,nstate,max_iterations
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    real(8),intent(in)::tolerance,metric_tolerance
    complex(8),intent(out)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:),maximum_residual,orthogonality_defect,metric_condition
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),intent(in),optional::initial_coefficients(:,:)
#ifdef USE_MPI
    complex(8),allocatable::qlocal(:,:),vector(:),hvector(:),svector(:),residual(:),direction(:),&
      hdirection(:),sdirection(:),gram(:,:),block_basis(:,:),block_coefficients(:,:),&
      block_h(:,:),block_s(:,:),block_residual(:,:),candidate_block(:,:),previous_direction(:,:),old_q(:,:)
    complex(8),allocatable::reference_initial(:,:)
    complex(8),allocatable::full_coefficients(:,:)
    complex(8)::value,reduced_h(2,2),reduced_s(2,2),ritz_vector(2)
    real(8),allocatable::block_eigenvalues(:)
    real(8),allocatable::full_eigenvalues(:)
    real(8)::ritz_value,pivot,residual_norm,block_maximum_residual,warm_defect,global_warm_defect
    integer,allocatable::ownership(:)
    integer::n,nlocal,i,j,band,iteration,ierr,local_bad,global_bad,block_size,maximum_block_size,&
      warm_local,warm_min,warm_max
    integer::scalar_local(3),scalar_min(3),scalar_max(3)
    real(8)::real_local(2),real_min(2),real_max(2)

    ok=.false.;message='';coefficients=(0d0,0d0);eigenvalues=0d0
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);metric_condition=huge(1d0)
    n=size(hrows,2);nlocal=size(row_ids)
    scalar_local=[n,nstate,max_iterations]
    call MPI_Allreduce(scalar_local,scalar_min,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(scalar_local,scalar_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    real_local=[tolerance,metric_tolerance]
    call MPI_Allreduce(real_local,real_min,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(real_local,real_max,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    warm_local=merge(1,0,present(initial_coefficients))
    call MPI_Allreduce(warm_local,warm_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(warm_local,warm_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=0
    if(any(scalar_min/=scalar_max).or.any(real_min/=real_max).or.warm_min/=warm_max)local_bad=1
    if(n<1.or.nstate<1.or.nstate>n.or.max_iterations<1)local_bad=1
    if(tolerance<=0d0.or.metric_tolerance<=0d0.or.metric_tolerance>=1d0)local_bad=1
    if(.not.all(ieee_is_finite(real_local)))local_bad=1
    if(size(hrows,1)/=nlocal.or.any(shape(srows)/=shape(hrows)))local_bad=1
    if(any(shape(coefficients)/=[n,nstate]).or.size(eigenvalues)/=nstate)local_bad=1
    if(any(row_ids<1_8).or.any(row_ids>int(n,8)))local_bad=1
    if(.not.finite_matrix(hrows).or..not.finite_matrix(srows))local_bad=1
    if(present(initial_coefficients))then
      if(any(shape(initial_coefficients)/=[n,nstate]).or..not.finite_matrix(initial_coefficients))local_bad=1
    endif
    if(int(n,int64)>0_int64.and.int(n,int64)>int(huge(1),int64)/int(n,int64))local_bad=1
    if(int(nstate,int64)>0_int64.and.int(nstate,int64)>int(huge(1),int64)/int(nstate,int64))local_bad=1
    if(int(nstate,int64)>int(huge(1),int64)/3_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed overlapping-Wannier coefficient solve contract';return;endif
    if(warm_max==1)then
      allocate(reference_initial,source=initial_coefficients)
      call MPI_Bcast(reference_initial,n*nstate,MPI_DOUBLE_COMPLEX,0,comm,ierr)
      warm_defect=maxval(abs(initial_coefficients-reference_initial))
      call MPI_Allreduce(warm_defect,global_warm_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(global_warm_defect>metric_tolerance)then
        message='rank-inconsistent warm-start coefficient block';return
      endif
    endif

    allocate(ownership(n));ownership=0
    do i=1,nlocal;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(ownership/=1))then;message='coefficient rows do not form a unique global partition';return;endif
    call distributed_hermiticity(comm,row_ids,hrows,tolerance,local_bad)
    call distributed_hermiticity(comm,row_ids,srows,metric_tolerance,global_bad)
    if(local_bad/=0.or.global_bad/=0)then;message='non-Hermitian distributed H or S';return;endif

    allocate(qlocal(nlocal,n),vector(n),hvector(nlocal),svector(nlocal),residual(nlocal),&
      direction(nlocal),hdirection(nlocal),sdirection(nlocal),gram(nstate,nstate))
    qlocal=(0d0,0d0)
    do j=1,n
      do i=1,nlocal
        if(row_ids(i)==int(j,8))qlocal(i,j)=1d0
      enddo
      call s_orthogonalize(comm,row_ids,srows,qlocal(:,j),qlocal(:,1:j-1),pivot,ok)
      if(.not.ok)then;message='overlap metric is indefinite or numerically singular';return;endif
    enddo
    call estimate_metric_condition(comm,row_ids,srows,qlocal,metric_condition,ok)
    if(.not.ok.or.metric_condition*metric_tolerance>=1d0)then
      ok=.false.;message='overlap metric condition gate failed';return
    endif
    if(.not.present(initial_coefficients))then
      allocate(full_coefficients(n,n),full_eigenvalues(n))
      do j=1,n;call gather_vector(comm,row_ids,qlocal(:,j),full_coefficients(:,j));end do
      call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal,full_coefficients,&
        full_eigenvalues,ok)
      if(.not.ok)then;message='cold full-space Rayleigh-Ritz failed';return;end if
      coefficients=full_coefficients(:,1:nstate);eigenvalues=full_eigenvalues(1:nstate)
      call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
        maximum_residual,orthogonality_defect)
      if(maximum_residual>10d0*tolerance.or.orthogonality_defect>10d0*tolerance)then
        ok=.false.;message='cold full-space generalized residual exceeds numerical quality gate';return
      end if
      ok=.true.;message='';return
    end if
    if(present(initial_coefficients))then
      do j=1,nstate
        do i=1,nlocal
          qlocal(i,j)=initial_coefficients(int(row_ids(i)),j)
        enddo
        call s_orthogonalize(comm,row_ids,srows,qlocal(:,j),qlocal(:,1:j-1),pivot,ok)
        if(.not.ok)then;message='warm-start coefficient block lost metric rank';return;endif
      enddo
    endif

    maximum_block_size=min(n,3*nstate)
    allocate(block_basis(nlocal,maximum_block_size),block_coefficients(n,maximum_block_size),&
      block_eigenvalues(maximum_block_size),block_h(nlocal,nstate),block_s(nlocal,nstate),&
      block_residual(nlocal,nstate),candidate_block(nlocal,2*nstate),&
      previous_direction(nlocal,nstate),old_q(nlocal,nstate))
    previous_direction=(0d0,0d0)
    do iteration=1,max_iterations
      do band=1,nstate
        call gather_vector(comm,row_ids,qlocal(:,band),coefficients(:,band))
      enddo
      call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal(:,1:nstate),coefficients,eigenvalues,ok)
      if(.not.ok)then;message='coefficient block Rayleigh-Ritz failed';return;endif
      block_h=matmul(hrows,coefficients);block_s=matmul(srows,coefficients)
      block_maximum_residual=0d0
      do band=1,nstate
        block_residual(:,band)=block_h(:,band)-eigenvalues(band)*block_s(:,band)
        call relative_residual(comm,block_residual(:,band),block_h(:,band),block_s(:,band),&
          eigenvalues(band),residual_norm)
        if(.not.ieee_is_finite(residual_norm))then
          ok=.false.;message='nonfinite coefficient block residual';return
        endif
        block_maximum_residual=max(block_maximum_residual,residual_norm)
      enddo
      if(block_maximum_residual<=tolerance)exit
      block_basis=(0d0,0d0);block_basis(:,1:nstate)=qlocal(:,1:nstate);block_size=nstate
      old_q=qlocal(:,1:nstate)
      candidate_block(:,1:nstate)=-block_residual
      candidate_block(:,nstate+1:2*nstate)=previous_direction
      call append_s_orthonormal_block(comm,row_ids,srows,&
        candidate_block(:,1:merge(2*nstate,nstate,iteration>1)),block_basis,block_size,ok)
      if(.not.ok)then;message='coefficient block residual orthogonalization failed';return;endif
      if(block_size==nstate)then
        ok=.false.;message='coefficient block residual space lost metric rank';return
      endif
      do j=1,block_size
        call gather_vector(comm,row_ids,block_basis(:,j),block_coefficients(:,j))
      enddo
      call block_rayleigh_ritz(comm,row_ids,hrows,srows,block_basis(:,1:block_size),&
        block_coefficients(:,1:block_size),block_eigenvalues(1:block_size),ok)
      if(.not.ok)then;message='coefficient expanded block Rayleigh-Ritz failed';return;endif
      previous_direction=old_q
      qlocal(:,1:nstate)=block_basis(:,1:nstate)
      if(block_size==n)then
        eigenvalues=block_eigenvalues(1:nstate)
        coefficients=block_coefficients(:,1:nstate)
        call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
          maximum_residual,orthogonality_defect)
        if(maximum_residual<=10d0*tolerance.and.orthogonality_defect<=10d0*tolerance)exit
        ok=.false.;message='full-space Ritz residual exceeds numerical quality gate';return
      endif
    enddo
    if(iteration>max_iterations)then
      ok=.false.
      write(message,'(a,i0,a,es12.4,a,es12.4)')&
        'coefficient block iteration did not converge: iterations=',max_iterations,&
        ' residual=',block_maximum_residual,' metric_condition=',metric_condition
      return
    endif
    do band=1,nstate
      call gather_vector(comm,row_ids,qlocal(:,band),coefficients(:,band))
    enddo
    call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal(:,1:nstate),coefficients,eigenvalues,ok)
    if(.not.ok)then;message='final block Rayleigh-Ritz failed';return;endif
    call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
      maximum_residual,orthogonality_defect)
    if(maximum_residual>10d0*tolerance.or.orthogonality_defect>10d0*tolerance)then
      ok=.false.;message='coefficient residual or S-orthonormality gate failed';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier coefficient solve requires MPI'
#endif
  end subroutine

#ifdef USE_MPI
  subroutine gather_vector(comm,row_ids,local,global)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::local(:)
    complex(8),intent(out)::global(:)
    complex(8),allocatable::staged(:)
    integer::i,ierr
    allocate(staged(size(global)));staged=(0d0,0d0)
    do i=1,size(local);staged(int(row_ids(i)))=local(i);enddo
    call MPI_Allreduce(staged,global,size(global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine inner(comm,left,right,value)
    integer,intent(in)::comm
    complex(8),intent(in)::left(:),right(:)
    complex(8),intent(out)::value
    complex(8)::local
    integer::ierr
    local=sum(conjg(left)*right)
    call MPI_Allreduce(local,value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine s_orthogonalize(comm,row_ids,srows,vector,previous,pivot,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    real(8),intent(out)::pivot
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
    call inner(comm,vector,svector,value);pivot=real(value,8)
    ok=ieee_is_finite(pivot).and.pivot>0d0.and.abs(aimag(value))<=1d-10*max(1d0,abs(value))
    if(ok)vector=vector/sqrt(pivot)
  end subroutine

  subroutine s_project_out(comm,row_ids,srows,vector,previous,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    ok=finite_matrix(reshape(vector,[size(vector),1]))
  end subroutine

  subroutine distributed_hermiticity(comm,row_ids,rows,tolerance,bad)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::rows(:,:)
    real(8),intent(in)::tolerance
    integer,intent(out)::bad
    integer,parameter::row_batch_size=32
    complex(8),allocatable::local_batch(:,:),global_batch(:,:)
    real(8)::local_scale,scale,local_defect,defect
    integer::i,j,ierr,n,first,last,nbatch
    n=size(rows,2)
    allocate(local_batch(min(row_batch_size,n),n),global_batch(min(row_batch_size,n),n))
    local_scale=1d0
    if(size(rows)>0)local_scale=max(local_scale,maxval(abs(rows)))
    local_defect=0d0
    do first=1,n,row_batch_size
      last=min(n,first+row_batch_size-1);nbatch=last-first+1
      local_batch=(0d0,0d0)
      do i=1,size(row_ids)
        if(row_ids(i)>=int(first,8).and.row_ids(i)<=int(last,8))&
          local_batch(int(row_ids(i))-first+1,:)=rows(i,:)
      enddo
      call MPI_Allreduce(local_batch,global_batch,size(local_batch),&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      do i=1,size(row_ids)
        do j=first,last
          local_defect=max(local_defect,&
            abs(rows(i,j)-conjg(global_batch(j-first+1,int(row_ids(i))))))
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_scale,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_defect,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    bad=merge(1,0,defect>tolerance*scale)
  end subroutine

  subroutine relative_residual(comm,residual,hvector,svector,eigenvalue,value)
    integer,intent(in)::comm
    complex(8),intent(in)::residual(:),hvector(:),svector(:)
    real(8),intent(in)::eigenvalue
    real(8),intent(out)::value
    complex(8)::z
    real(8)::nr,nh,ns
    call inner(comm,residual,residual,z);nr=sqrt(max(0d0,real(z,8)))
    call inner(comm,hvector,hvector,z);nh=sqrt(max(0d0,real(z,8)))
    call inner(comm,svector,svector,z);ns=sqrt(max(0d0,real(z,8)))
    value=nr/max(tiny(1d0),nh+abs(eigenvalue)*ns)
  end subroutine

  subroutine estimate_metric_condition(comm,row_ids,srows,q,condition,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),q(:,:)
    real(8),intent(out)::condition
    logical,intent(out)::ok
    complex(8),allocatable::gathered_column(:),root_q(:,:)
    integer(8),allocatable::gathered_ids(:)
    complex(8)::gram_value
    real(8)::squares_s,squares_inverse,local_squares,local_scale,scale_s,&
      norm_s,norm_inverse
    real(8)::factor_s
    integer,allocatable::counts(:),displacements(:)
    integer::ierr,n,i,j,k,rank,nproc,nlocal
    n=size(q,2)
    nlocal=size(row_ids)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    local_scale=0d0
    if(size(srows)>0)local_scale=maxval(abs(srows))
    call MPI_Allreduce(local_scale,scale_s,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(scale_s<=0d0.or..not.ieee_is_finite(scale_s))then;ok=.false.;return;endif
    local_squares=sum((abs(srows)/scale_s)**2)
    call MPI_Allreduce(local_squares,squares_s,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    allocate(counts(nproc),displacements(nproc))
    call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm,ierr)
    if(rank==0)then
      displacements(1)=0
      do i=2,nproc;displacements(i)=displacements(i-1)+counts(i-1);enddo
      allocate(gathered_ids(n),gathered_column(n),root_q(n,n))
    else
      allocate(gathered_ids(1),gathered_column(1),root_q(1,1))
    endif
    call MPI_Gatherv(row_ids,nlocal,MPI_INTEGER8,gathered_ids,counts,displacements,&
      MPI_INTEGER8,0,comm,ierr)
    do j=1,n
      call MPI_Gatherv(q(:,j),nlocal,MPI_DOUBLE_COMPLEX,gathered_column,counts,displacements,&
        MPI_DOUBLE_COMPLEX,0,comm,ierr)
      if(rank==0)then
        do k=1,n;root_q(int(gathered_ids(k)),j)=gathered_column(k);enddo
      endif
    enddo
    squares_inverse=0d0
    if(rank==0)then
      do j=1,n;do i=1,n
        gram_value=dot_product(root_q(:,i),root_q(:,j))
        squares_inverse=squares_inverse+abs(gram_value)**2
      enddo;enddo
    endif
    call MPI_Bcast(squares_inverse,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    factor_s=sqrt(squares_s)
    if(factor_s<=0d0)then;ok=.false.;return;endif
    if(squares_inverse<=0d0.or..not.ieee_is_finite(squares_inverse))then;ok=.false.;return;endif
    if(factor_s>1d0)then
      if(scale_s>huge(1d0)/factor_s)then
        condition=huge(1d0);ok=.false.;return
      endif
    endif
    norm_s=scale_s*factor_s;norm_inverse=sqrt(squares_inverse)
    if(norm_inverse>1d0)then
      if(norm_s>huge(1d0)/norm_inverse)then
        condition=huge(1d0);ok=.false.;return
      else
        condition=norm_s*norm_inverse
      endif
    else
      condition=norm_s*norm_inverse
    endif
    ok=ieee_is_finite(condition).and.condition>=1d0
  end subroutine

  subroutine block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal,coefficients,eigenvalues,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    complex(8),intent(inout)::qlocal(:,:)
    complex(8),intent(inout)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:)
    logical,intent(out)::ok
    complex(8),allocatable::hq(:,:),projected(:,:),local_projected(:,:),rotated_q(:,:),work(:)
    real(8),allocatable::rwork(:)
    integer::p,ierr,n,lapack_info,lwork
    external::zheev
    n=size(eigenvalues)
    hq=matmul(hrows,coefficients)
    local_projected=matmul(conjg(transpose(qlocal)),hq)
    allocate(projected(size(eigenvalues),size(eigenvalues)))
    call MPI_Allreduce(local_projected,projected,size(eigenvalues)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    projected=0.5d0*(projected+conjg(transpose(projected)))
    if(.not.finite_matrix(projected))then;ok=.false.;return;endif
    allocate(rwork(max(1,3*n-2)))
    lwork=max(1,2*n-1)
    allocate(work(lwork))
    call zheev('V','U',n,projected,n,eigenvalues,work,lwork,rwork,lapack_info)
    if(lapack_info/=0.or..not.finite_matrix(projected))then;ok=.false.;return;endif
    allocate(rotated_q(size(qlocal,1),n))
    rotated_q=matmul(qlocal,projected)
    qlocal=rotated_q
    do p=1,size(eigenvalues);call gather_vector(comm,row_ids,qlocal(:,p),coefficients(:,p));enddo
    ok=all(ieee_is_finite(eigenvalues)).and.finite_matrix(coefficients)
  end subroutine

  subroutine append_s_orthonormal_block(comm,row_ids,srows,candidates,basis,basis_size,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),candidates(:,:)
    complex(8),intent(inout)::basis(:,:)
    integer,intent(inout)::basis_size
    logical,intent(out)::ok
    complex(8),allocatable::work_block(:,:),global_block(:,:),sblock(:,:),cross(:,:),local_cross(:,:),&
      gram(:,:),local_gram(:,:),work(:),rotated(:,:)
    real(8),allocatable::spectrum(:),rwork(:)
    integer::ncandidate,n,keep,i,j,ierr,info,lwork,projection_pass
    external::zheev
    ncandidate=size(candidates,2);n=size(srows,2)
    allocate(work_block(size(candidates,1),ncandidate),global_block(n,ncandidate),&
      sblock(size(candidates,1),ncandidate))
    work_block=candidates
    if(basis_size>0)then
      allocate(local_cross(basis_size,ncandidate),cross(basis_size,ncandidate))
      do projection_pass=1,2
        call gather_block(comm,row_ids,work_block,global_block)
        sblock=matmul(srows,global_block)
        local_cross=matmul(conjg(transpose(basis(:,1:basis_size))),sblock)
        call MPI_Allreduce(local_cross,cross,basis_size*ncandidate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        work_block=work_block-matmul(basis(:,1:basis_size),cross)
      enddo
    endif
    call gather_block(comm,row_ids,work_block,global_block)
    sblock=matmul(srows,global_block)
    allocate(local_gram(ncandidate,ncandidate),gram(ncandidate,ncandidate),spectrum(ncandidate))
    local_gram=matmul(conjg(transpose(work_block)),sblock)
    call MPI_Allreduce(local_gram,gram,ncandidate*ncandidate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    gram=0.5d0*(gram+conjg(transpose(gram)))
    allocate(rwork(max(1,3*ncandidate-2)))
    lwork=max(1,2*ncandidate-1);allocate(work(lwork))
    call zheev('V','U',ncandidate,gram,ncandidate,spectrum,work,lwork,rwork,info)
    if(info/=0.or..not.all(ieee_is_finite(spectrum)))then;ok=.false.;return;endif
    keep=count(spectrum>max(tiny(1d0),1d-12*max(0d0,maxval(spectrum))))
    keep=min(keep,size(basis,2)-basis_size)
    if(keep<1)then;ok=.true.;return;endif
    allocate(rotated(size(work_block,1),keep))
    do j=1,keep
      i=ncandidate-keep+j
      rotated(:,j)=matmul(work_block,gram(:,i))/sqrt(spectrum(i))
    enddo
    basis(:,basis_size+1:basis_size+keep)=rotated
    basis_size=basis_size+keep
    ok=finite_matrix(rotated)
  end subroutine

  subroutine gather_block(comm,row_ids,local,global)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::local(:,:)
    complex(8),intent(out)::global(:,:)
    complex(8),allocatable::staged(:,:)
    integer::i,ierr
    allocate(staged(size(global,1),size(global,2)));staged=(0d0,0d0)
    do i=1,size(local,1);staged(int(row_ids(i)),:)=local(i,:);enddo
    call MPI_Allreduce(staged,global,size(global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine coefficient_diagnostics(comm,row_ids,hrows,srows,c,e,residual,orthogonality)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:),c(:,:)
    real(8),intent(in)::e(:)
    real(8),intent(out)::residual,orthogonality
    complex(8),allocatable::rlocal(:,:),gram(:,:),local_gram(:,:),hv(:),sv(:)
    real(8)::band_residual
    integer::j,ierr
    allocate(rlocal(size(row_ids),size(e)),gram(size(e),size(e)),local_gram(size(e),size(e)))
    allocate(hv(size(row_ids)),sv(size(row_ids)));residual=0d0
    do j=1,size(e)
      hv=matmul(hrows,c(:,j));sv=matmul(srows,c(:,j));rlocal(:,j)=hv-e(j)*sv
      call relative_residual(comm,rlocal(:,j),hv,sv,e(j),band_residual)
      residual=max(residual,band_residual)
    enddo
    local_gram=matmul(conjg(transpose(c(int(row_ids),:))),matmul(srows,c))
    call MPI_Allreduce(local_gram,gram,size(e)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do j=1,size(e);gram(j,j)=gram(j,j)-1d0;enddo
    orthogonality=maxval(abs(gram))
    if(.not.all(ieee_is_finite([residual,orthogonality])))then
      residual=huge(1d0);orthogonality=huge(1d0)
    endif
  end subroutine

  subroutine lowest_generalized_2x2(h,s,eigenvalue,eigenvector,ok)
    complex(8),intent(in)::h(2,2),s(2,2)
    real(8),intent(out)::eigenvalue
    complex(8),intent(out)::eigenvector(2)
    logical,intent(out)::ok
    complex(8)::l(2,2),linv(2,2),a(2,2),u(2),phase,hs(2,2),temporary(2,2)
    real(8)::l11,l22,s22_root,disc_scaled,norm,scale,delta,t,c,absb_scaled,&
      difference_scaled,eigenvalue_scaled,xdisc,ydisc,metric_ratio,phase_threshold
    real(8)::hamiltonian_scale,total_scale
    complex(8)::quotient
    integer::i,j,k
    ok=.false.;eigenvalue=0d0;eigenvector=(0d0,0d0)
    if(.not.finite_matrix(h).or..not.finite_matrix(s))return
    if(real(s(1,1),8)<=0d0.or.real(s(2,2),8)<=0d0)return
    l11=sqrt(real(s(1,1),8))
    l=(0d0,0d0);l(1,1)=l11
    call safe_complex_divide(s(2,1),l11,l(2,1),ok);if(.not.ok)return
    s22_root=sqrt(real(s(2,2),8))
    if(safe_complex_abs(l(2,1))>=s22_root)return
    metric_ratio=safe_complex_abs(l(2,1))/s22_root
    l22=s22_root*sqrt((1d0-metric_ratio)*(1d0+metric_ratio));l(2,2)=l22
    if(.not.all(ieee_is_finite([l11,l22])).or.min(l11,l22)<=0d0)return
    if(l11<1d0/huge(1d0).or.l22<1d0/huge(1d0))return
    linv=(0d0,0d0);linv(1,1)=1d0/l11;linv(2,2)=1d0/l22
    call safe_complex_divide(l(2,1),l11,quotient,ok);if(.not.ok)return
    call safe_complex_divide(-quotient,l22,linv(2,1),ok);if(.not.ok)return
    hamiltonian_scale=max(tiny(1d0),max(maxval(abs(real(h,8))),maxval(abs(aimag(h)))))
    hs=h/hamiltonian_scale
    temporary=(0d0,0d0);a=(0d0,0d0)
    do j=1,2;do i=1,2;do k=1,2
      temporary(i,j)=temporary(i,j)+linv(i,k)*hs(k,j)
    enddo;enddo;enddo
    do j=1,2;do i=1,2;do k=1,2
      a(i,j)=a(i,j)+temporary(i,k)*conjg(linv(j,k))
    enddo;enddo;enddo
    a=0.5d0*(a+conjg(transpose(a)))
    if(.not.finite_matrix(a))return
    scale=max(tiny(1d0),max(maxval(abs(real(a,8))),maxval(abs(aimag(a)))))
    delta=real(a(1,1),8)/scale-real(a(2,2),8)/scale
    xdisc=abs(delta)
    absb_scaled=safe_real_hypot(real(a(1,2),8)/scale,aimag(a(1,2))/scale)
    ydisc=2d0*absb_scaled
    disc_scaled=sqrt(xdisc*xdisc+ydisc*ydisc)
    eigenvalue_scaled=0.5d0*(real(a(1,1),8)/scale+real(a(2,2),8)/scale-disc_scaled)
    if(scale>1d0)then
      if(hamiltonian_scale>huge(1d0)/scale)return
    endif
    total_scale=scale*hamiltonian_scale
    if(abs(eigenvalue_scaled)>1d0)then
      if(total_scale>huge(1d0)/abs(eigenvalue_scaled))return
    endif
    eigenvalue=total_scale*eigenvalue_scaled
    phase_threshold=epsilon(1d0)*max(1d0,1d0/scale)
    if(absb_scaled>phase_threshold)then
      phase=(a(1,2)/scale)/absb_scaled
      if(real(a(1,1),8)<=real(a(2,2),8))then
        difference_scaled=real(a(2,2),8)/scale-real(a(1,1),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[cmplx(c,0d0,8),-t*c*conjg(phase)]
      else
        difference_scaled=real(a(1,1),8)/scale-real(a(2,2),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[-t*c*phase,cmplx(c,0d0,8)]
      endif
    else
      u=(0d0,0d0);u(merge(1,2,real(a(1,1),8)<=real(a(2,2),8)))=1d0
    endif
    scale=max(safe_complex_abs(u(1)),safe_complex_abs(u(2)))
    if(scale<=0d0.or..not.ieee_is_finite(scale))return
    norm=scale*safe_real_hypot(safe_complex_abs(u(1))/scale,safe_complex_abs(u(2))/scale)
    if(norm<=0d0.or..not.ieee_is_finite(norm))return
    u=u/norm;eigenvector=matmul(conjg(transpose(linv)),u)
    ok=all(ieee_is_finite(real(eigenvector,8))).and.all(ieee_is_finite(aimag(eigenvector)))
  end subroutine
#endif

  real(8) function safe_complex_abs(value)
    complex(8),intent(in)::value
    safe_complex_abs=safe_real_hypot(real(value,8),aimag(value))
  end function

  real(8) function safe_real_hypot(x,y)
    real(8),intent(in)::x,y
    real(8)::larger,ratio,factor
    larger=max(abs(x),abs(y))
    if(larger==0d0)then
      safe_real_hypot=0d0
      return
    endif
    ratio=min(abs(x),abs(y))/larger
    factor=sqrt(1d0+ratio*ratio)
    if(larger>huge(1d0)/factor)then
      safe_real_hypot=huge(1d0)
    else
      safe_real_hypot=larger*factor
    endif
  end function

  subroutine safe_complex_divide(value,divisor,result,ok)
    complex(8),intent(in)::value
    real(8),intent(in)::divisor
    complex(8),intent(out)::result
    logical,intent(out)::ok
    real(8)::magnitude
    result=(0d0,0d0);ok=.false.
    if(divisor<=0d0.or..not.ieee_is_finite(divisor))return
    magnitude=safe_complex_abs(value)
    if(divisor<1d0)then
      if(magnitude>huge(1d0)*divisor)return
    endif
    result=value/divisor
    ok=ieee_is_finite(real(result,8))
    if(ok)ok=ieee_is_finite(aimag(result))
  end subroutine

  logical function finite_matrix(matrix)
    complex(8),intent(in)::matrix(:,:)
    integer::i,j
    real(8)::real_part,imaginary_part
    finite_matrix=.true.
    do j=1,size(matrix,2)
      do i=1,size(matrix,1)
        real_part=real(matrix(i,j),8)
        imaginary_part=aimag(matrix(i,j))
        if(real_part/=real_part)then
          finite_matrix=.false.
          return
        endif
        if(imaginary_part/=imaginary_part)then
          finite_matrix=.false.
          return
        endif
        if(abs(real_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
        if(abs(imaginary_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
      enddo
    enddo
  end function
end module
