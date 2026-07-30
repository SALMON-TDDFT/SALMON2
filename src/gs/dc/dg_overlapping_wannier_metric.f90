#include "config.h"
module dg_overlapping_wannier_metric
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_metric,assemble_dg_overlapping_wannier_metric_rows
contains
  subroutine assemble_dg_overlapping_wannier_metric_rows(comm,nwann,row_ids,core_ids,weights,values,&
      pairs,expected_core_count,relative_threshold,metric_rows,retained_spectrum,minimum_eigenvalue,&
      condition_number,rejected_rank,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::row_ids(:),core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),relative_threshold
    complex(real64),intent(in)::values(:,:)
    logical,intent(in)::pairs(:,:)
    complex(real64),allocatable,intent(out)::metric_rows(:,:)
    real(real64),allocatable,intent(out)::retained_spectrum(:)
    real(real64),intent(out)::minimum_eigenvalue,condition_number
    integer,intent(out)::rejected_rank,ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,nwann_min,nwann_max,total_count,total_rows,&
      r,i,j,p,nrows,batch_first,batch_count,status,info,lwork,nretained
    integer,allocatable::core_counts(:),core_displs(:),row_counts(:),row_displs(:)
    integer(int64),allocatable::all_core_ids(:),all_row_ids(:),validation_ids(:)
    integer(int64)::expected_min,expected_max
    complex(real64),allocatable::partial(:,:),reduced(:,:),block(:,:),root_metric(:,:),work(:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::threshold_min,threshold_max,scale,defect,largest,cutoff
    logical::shape_ok,finite_values
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character(1),intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface

    ok=.false.;message='';minimum_eigenvalue=0d0;condition_number=huge(1d0)
    rejected_rank=0;ownership_count=0;local_bad=0;status=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(relative_threshold,threshold_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(relative_threshold,threshold_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(nwann_min/=nwann_max.or.expected_min/=expected_max.or.threshold_min/=threshold_max)local_bad=1
    shape_ok=size(weights)==size(core_ids).and.size(values,1)==nwann.and.&
      size(values,2)==size(core_ids).and.size(pairs,1)==nwann.and.size(pairs,2)==size(core_ids)
    finite_values=all(ieee_is_finite(weights))
    if(shape_ok)finite_values=finite_values.and.all(ieee_is_finite(real(values))).and.&
      all(ieee_is_finite(aimag(values)))
    if(nwann<=0.or.expected_core_count<=0_int64.or.relative_threshold<=0d0.or..not.shape_ok.or.&
        .not.finite_values.or.any(core_ids<=0_int64).or.any(weights<=0d0).or..not.all(pairs).or.&
        any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid or inconsistent row-owned metric contract';return
    endif

    allocate(core_counts(nproc),core_displs(nproc),row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then;message='metric ownership metadata collective failed';return;endif
    total_count=0;total_rows=0
    do r=1,nproc
      core_displs(r)=total_count;row_displs(r)=total_rows
      if(core_counts(r)<0.or.row_counts(r)<0.or.total_count>huge(total_count)-core_counts(r).or.&
          total_rows>huge(total_rows)-row_counts(r))local_bad=1
      if(local_bad==0)then
        total_count=total_count+core_counts(r);total_rows=total_rows+row_counts(r)
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS.or.int(total_count,int64)/=expected_core_count.or.&
        total_rows/=nwann)then;message='missing or extra row/core owner in metric assembly';return;endif
    allocate(all_core_ids(total_count),all_row_ids(total_rows),validation_ids(max(total_count,total_rows)))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_core_ids,core_counts,core_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then;message='metric ownership payload collective failed';return;endif
    validation_ids(1:total_count)=all_core_ids;call sort_ids(validation_ids(1:total_count))
    do i=1,total_count;if(validation_ids(i)/=int(i,int64))local_bad=1;enddo
    validation_ids(1:total_rows)=all_row_ids;call sort_ids(validation_ids(1:total_rows))
    do i=1,total_rows;if(validation_ids(i)/=int(i,int64))local_bad=1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='duplicate or missing row/core owner in metric assembly';return
    endif

    allocate(metric_rows(size(row_ids),nwann));metric_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial(batch_count,nwann),reduced(batch_count,nwann));partial=(0d0,0d0)
        do p=1,size(core_ids);do j=1,nwann;do i=1,batch_count
          partial(i,j)=partial(i,j)+weights(p)*conjg(values(&
            int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*values(j,p)
        enddo;enddo;enddo
        call MPI_Reduce(partial,reduced,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==r)metric_rows(batch_first:batch_first+batch_count-1,:)=reduced
        deallocate(partial,reduced)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then;message='metric row reduction failed';return;endif

    if(rank==0)allocate(root_metric(nwann,nwann))
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(block(batch_count,nwann));block=(0d0,0d0)
        if(rank==r)block=metric_rows(batch_first:batch_first+batch_count-1,:)
        call MPI_Bcast(block,batch_count*nwann,MPI_DOUBLE_COMPLEX,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==0)then
          do i=1,batch_count
            root_metric(int(all_row_ids(row_displs(r+1)+batch_first+i-1)),:)=block(i,:)
          enddo
        endif
        deallocate(block)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then;message='metric root diagnostic gather failed';return;endif
    allocate(eigenvalues(nwann));eigenvalues=0d0
    if(rank==0)then
      if(.not.all(ieee_is_finite(real(root_metric))).or.&
          .not.all(ieee_is_finite(aimag(root_metric))))then
        status=6
      else
        scale=max(1d0,maxval(abs(root_metric)))
        defect=maxval(abs(root_metric-conjg(transpose(root_metric))))
        if(defect>relative_threshold*scale)then
          status=1
        else
          root_metric=0.5d0*(root_metric+conjg(transpose(root_metric)))
          allocate(rwork(max(1,3*nwann-2)),work(1));lwork=-1
          call zheev('N','U',nwann,root_metric,nwann,eigenvalues,work,lwork,rwork,info)
          if(info/=0)then
            status=2
          else
            lwork=max(1,int(real(work(1))));deallocate(work);allocate(work(lwork))
            call zheev('N','U',nwann,root_metric,nwann,eigenvalues,work,lwork,rwork,info)
            if(info/=0)status=2
          endif
        endif
        if(status==0)then
          largest=maxval(eigenvalues);cutoff=relative_threshold*largest
          if(largest<=0d0.or.minval(eigenvalues)<-cutoff)then
            status=3
          else
            nretained=count(eigenvalues>cutoff)
            if(nretained==0)then
              status=4
            else
              rejected_rank=nwann-nretained
              minimum_eigenvalue=eigenvalues(rejected_rank+1)
              condition_number=largest/minimum_eigenvalue
            endif
          endif
        endif
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(rejected_rank,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(minimum_eigenvalue,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(condition_number,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Bcast(eigenvalues,nwann,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)status=5
    if(status/=0)then
      select case(status)
      case(1);message='overlapping-Wannier metric Hermiticity defect exceeds tolerance'
      case(2);message='overlapping-Wannier metric eigensolve failed'
      case(3);message='overlapping-Wannier metric is not positive semidefinite'
      case(4);message='overlapping-Wannier metric has no positive rank'
      case(6);message='overlapping-Wannier assembled metric is nonfinite'
      case default;message='metric spectral diagnostic broadcast failed'
      end select
      return
    endif
    nretained=nwann-rejected_rank;allocate(retained_spectrum(nretained))
    retained_spectrum=eigenvalues(rejected_rank+1:nwann)
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='row-owned overlapping-Wannier metric assembly requires MPI'
    minimum_eigenvalue=0d0;condition_number=huge(1d0);rejected_rank=0;ownership_count=0
#endif
  end subroutine

  subroutine assemble_dg_overlapping_wannier_metric(comm,nwann,core_ids,weights,values,pairs,&
      expected_core_count,relative_threshold,metric,retained_vectors,retained_spectrum,&
      minimum_eigenvalue,condition_number,rejected_rank,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),relative_threshold
    complex(real64),intent(in)::values(:,:)
    logical,intent(in)::pairs(:,:)
    complex(real64),allocatable,intent(out)::metric(:,:),retained_vectors(:,:)
    real(real64),allocatable,intent(out)::retained_spectrum(:)
    real(real64),intent(out)::minimum_eigenvalue,condition_number
    integer,intent(out)::rejected_rank,ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,i,j,p,total_count,info,lwork,nretained,matrix_count,&
      nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    complex(real64),allocatable::local_metric(:,:),eigenvectors(:,:),work(:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::cutoff,largest,threshold_min,threshold_max,hermiticity_defect,scale
    integer(int64)::matrix_count64,expected_min,expected_max
    logical::finite_values
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character(1),intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface

    ok=.false.;message='';minimum_eigenvalue=0d0;condition_number=huge(1d0)
    rejected_rank=0;ownership_count=0
    call MPI_Comm_size(comm,nproc,ierr)
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(relative_threshold,threshold_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(relative_threshold,threshold_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max.or.threshold_min/=threshold_max)then
      message='inconsistent metric assembly contract across ranks';return
    endif
    finite_values=.true.
    do p=1,size(values,2)
      do i=1,size(values,1)
        finite_values=finite_values.and.ieee_is_finite(real(values(i,p))).and.&
          ieee_is_finite(aimag(values(i,p)))
      enddo
    enddo
    local_bad=0
    if(nwann<=0.or.expected_core_count<=0_int64.or.relative_threshold<=0d0) local_bad=1
    if(size(weights)/=size(core_ids).or.size(values,1)/=nwann.or.&
        size(values,2)/=size(core_ids))local_bad=1
    if(size(pairs,1)/=nwann.or.size(pairs,2)/=size(core_ids))local_bad=1
    if(any(core_ids<=0_int64).or.any(weights<=0d0).or..not.all(ieee_is_finite(weights)))local_bad=1
    if(.not.finite_values.or..not.all(pairs))local_bad=1
    if(nwann<=0)then
      matrix_count64=0_int64
    else if(int(nwann,int64)>huge(1_int64)/int(nwann,int64))then
      local_bad=1;matrix_count64=0_int64
    else
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='invalid or missing unique-core owner pairs before metric payload collective';return
    endif

    allocate(counts(nproc),displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_count=0;displs(1)=0
    do i=1,nproc
      if(counts(i)<0.or.total_count>huge(total_count)-counts(i))then
        local_bad=1;exit
      endif
      if(i>1)displs(i)=total_count
      total_count=total_count+counts(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_core_count)then
      message='missing or extra unique-core quadrature owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_ids,counts,displs,MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then
        message='duplicate or missing unique-core quadrature owner';return
      endif
    enddo

    allocate(local_metric(nwann,nwann),metric(nwann,nwann));local_metric=(0d0,0d0)
    do p=1,size(core_ids)
      do j=1,nwann
        do i=1,nwann
          local_metric(i,j)=local_metric(i,j)+weights(p)*conjg(values(i,p))*values(j,p)
        enddo
      enddo
    enddo
    matrix_count=int(matrix_count64)
    call MPI_Allreduce(local_metric,metric,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale=max(1d0,maxval(abs(metric)))
    hermiticity_defect=maxval(abs(metric-conjg(transpose(metric))))
    if(hermiticity_defect>relative_threshold*scale)then
      message='overlapping-Wannier metric Hermiticity defect exceeds tolerance';return
    endif
    metric=0.5d0*(metric+conjg(transpose(metric)))
    allocate(eigenvectors,source=metric)
    allocate(eigenvalues(nwann),rwork(max(1,3*nwann-2)),work(1))
    lwork=-1
    call zheev('V','U',nwann,eigenvectors,nwann,eigenvalues,work,lwork,rwork,info)
    if(info/=0)then;message='LAPACK workspace query failed';return;endif
    lwork=max(1,int(real(work(1))))
    deallocate(work);allocate(work(lwork))
    call zheev('V','U',nwann,eigenvectors,nwann,eigenvalues,work,lwork,rwork,info)
    if(info/=0)then;message='overlapping-Wannier metric eigensolve failed';return;endif
    largest=maxval(eigenvalues);cutoff=relative_threshold*largest
    if(largest<=0d0.or.minval(eigenvalues)<-cutoff)then
      message='overlapping-Wannier metric is not positive semidefinite';return
    endif
    nretained=count(eigenvalues>cutoff)
    if(nretained==0)then;message='overlapping-Wannier metric has no positive rank';return;endif
    rejected_rank=nwann-nretained
    minimum_eigenvalue=eigenvalues(rejected_rank+1)
    condition_number=largest/minimum_eigenvalue
    allocate(retained_spectrum(nretained),retained_vectors(nwann,nretained))
    retained_spectrum=eigenvalues(rejected_rank+1:nwann)
    retained_vectors=eigenvectors(:,rejected_rank+1:nwann)
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='overlapping-Wannier metric assembly requires MPI'
    minimum_eigenvalue=0d0;condition_number=huge(1d0);rejected_rank=0;ownership_count=0
#endif
  end subroutine

  subroutine sort_ids(ids)
    integer(int64),intent(inout)::ids(:)
    integer::i,j
    integer(int64)::key
    do i=2,size(ids)
      key=ids(i);j=i-1
      do while(j>=1)
        if(ids(j)<=key)exit
        ids(j+1)=ids(j);j=j-1
      enddo
      ids(j+1)=key
    enddo
  end subroutine
end module dg_overlapping_wannier_metric
