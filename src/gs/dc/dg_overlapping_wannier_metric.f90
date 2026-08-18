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
  public::assemble_dg_eigenexa_cyclic_metric_block
  public::assemble_dg_stitched_overlap_density_rows
contains
  subroutine assemble_dg_stitched_overlap_density_rows(comm,nbasis,row_ids,physical_ids,&
      partition_weight,basis_values,density_values,cell_volume,expected_physical_count,&
      expected_electrons,tolerance,electron_count_tolerance,srows,rhorows,electron_count,&
      s_hermiticity,rho_hermiticity,&
      minimum_cholesky_pivot,pivot_condition,peak_elements,ok,message)
    integer,intent(in)::comm,nbasis
    integer(int64),intent(in)::row_ids(:),physical_ids(:),expected_physical_count
    real(real64),intent(in)::partition_weight(:),density_values(:),cell_volume,expected_electrons,tolerance,&
      electron_count_tolerance
    complex(real64),intent(in)::basis_values(:,:)
    complex(real64),allocatable,intent(out)::srows(:,:),rhorows(:,:)
    real(real64),intent(out)::electron_count,s_hermiticity,rho_hermiticity
    real(real64),intent(out)::minimum_cholesky_pivot,pivot_condition
    integer(int64),intent(out)::peak_elements
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,r,total_rows,nrows,batch_first,batch_count,&
      i,j,k,p,local_row,pivot_owner,pivot_local,owner,total_send,total_recv,slot,nowned
    integer,allocatable::row_counts(:),row_displs(:),send_counts(:),recv_counts(:),&
      send_displs(:),recv_displs(:),send_cursor(:)
    integer(int64),allocatable::all_row_ids(:),sorted_row_ids(:),send_ids(:),recv_ids(:)
    complex(real64),allocatable::partial_s(:,:),reduced_s(:,:),partial_rho(:,:),reduced_rho(:,:),&
      block_s(:,:),block_rho(:,:)
    complex(real64),allocatable::cholesky_rows(:,:),pivot_row(:)
    real(real64),allocatable::send_weights(:),recv_weights(:),owned_coverage(:)
    real(real64)::local_electrons,local_moments(3),global_moments(3),expected_moments(3),&
      local_s_defect,local_rho_defect,scale_s,scale_rho,pivot_value,maximum_cholesky_pivot
    ok=.false.;message='';electron_count=huge(1d0);s_hermiticity=huge(1d0)
    rho_hermiticity=huge(1d0);peak_elements=0_int64;local_bad=0
    minimum_cholesky_pivot=0d0;pivot_condition=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS.or.nbasis<1.or.expected_physical_count<1_int64.or.cell_volume<=0d0.or.&
        tolerance<=0d0.or.electron_count_tolerance<=0d0.or.&
        .not.ieee_is_finite(electron_count_tolerance).or.&
        size(partition_weight)/=size(physical_ids).or.&
        size(density_values)/=size(physical_ids).or.size(basis_values,1)/=nbasis.or.&
        size(basis_values,2)/=size(physical_ids).or.any(physical_ids<1_int64).or.&
        any(physical_ids>expected_physical_count).or.any(partition_weight<0d0).or.&
        any(row_ids<1_int64).or.any(row_ids>int(nbasis,int64)))local_bad=1
    if(.not.all(ieee_is_finite(partition_weight)).or..not.all(ieee_is_finite(density_values)).or.&
        .not.all(ieee_is_finite(real(basis_values))).or..not.all(ieee_is_finite(aimag(basis_values))).or.&
        .not.ieee_is_finite(expected_electrons))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid stitched overlap-density contract';return
    endif
    allocate(row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    total_rows=0
    do r=1,nproc;row_displs(r)=total_rows;total_rows=total_rows+row_counts(r);enddo
    if(total_rows/=nbasis)local_bad=1
    allocate(all_row_ids(total_rows),sorted_row_ids(total_rows))
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    sorted_row_ids=all_row_ids;call sort_ids(sorted_row_ids)
    do i=1,total_rows;if(sorted_row_ids(i)/=int(i,int64))local_bad=1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='duplicate or missing stitched matrix row owner';return
    endif
    allocate(send_counts(nproc),recv_counts(nproc),send_displs(nproc),recv_displs(nproc),&
      send_cursor(nproc));send_counts=0
    do p=1,size(physical_ids)
      owner=int(modulo(physical_ids(p)-1_int64,int(nproc,int64)))+1
      send_counts(owner)=send_counts(owner)+1
    enddo
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr)
    total_send=0;total_recv=0
    do r=1,nproc
      send_displs(r)=total_send;recv_displs(r)=total_recv
      total_send=total_send+send_counts(r);total_recv=total_recv+recv_counts(r)
    enddo
    allocate(send_ids(total_send),recv_ids(total_recv),send_weights(total_send),recv_weights(total_recv))
    send_cursor=send_displs
    do p=1,size(physical_ids)
      owner=int(modulo(physical_ids(p)-1_int64,int(nproc,int64)))+1
      send_cursor(owner)=send_cursor(owner)+1
      send_ids(send_cursor(owner))=physical_ids(p);send_weights(send_cursor(owner))=partition_weight(p)
    enddo
    call MPI_Alltoallv(send_ids,send_counts,send_displs,MPI_INTEGER8,recv_ids,recv_counts,recv_displs,&
      MPI_INTEGER8,comm,ierr)
    call MPI_Alltoallv(send_weights,send_counts,send_displs,MPI_DOUBLE_PRECISION,recv_weights,&
      recv_counts,recv_displs,MPI_DOUBLE_PRECISION,comm,ierr)
    nowned=int((expected_physical_count-int(rank,int64)+int(nproc,int64)-1_int64)/int(nproc,int64))
    allocate(owned_coverage(nowned));owned_coverage=0d0
    do p=1,total_recv
      if(modulo(recv_ids(p)-1_int64,int(nproc,int64))/=int(rank,int64))then
        local_bad=1;cycle
      endif
      slot=int((recv_ids(p)-1_int64)/int(nproc,int64))+1
      if(slot<1.or.slot>nowned)then;local_bad=1;cycle;endif
      owned_coverage(slot)=owned_coverage(slot)+recv_weights(p)
    enddo
    if(any(abs(owned_coverage-1d0)>tolerance))local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='stitched partition has missing, excess, or nonunit physical-grid coverage';return
    endif
    peak_elements=max(peak_elements,int(7*nproc+2*nbasis+2*total_send+2*total_recv+nowned,int64))
    deallocate(send_counts,recv_counts,send_displs,recv_displs,send_cursor,send_ids,recv_ids,&
      send_weights,recv_weights,owned_coverage)
    local_moments=0d0
    do p=1,size(physical_ids)
      local_moments(1)=local_moments(1)+partition_weight(p)
      local_moments(2)=local_moments(2)+partition_weight(p)*real(physical_ids(p),real64)
      local_moments(3)=local_moments(3)+partition_weight(p)*real(physical_ids(p),real64)**2
    enddo
    call MPI_Allreduce(local_moments,global_moments,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    expected_moments=[real(expected_physical_count,real64),&
      0.5d0*real(expected_physical_count,real64)*real(expected_physical_count+1_int64,real64),&
      real(expected_physical_count,real64)*real(expected_physical_count+1_int64,real64)*&
      real(2_int64*expected_physical_count+1_int64,real64)/6d0]
    if(maxval(abs(global_moments-expected_moments))>tolerance*max(1d0,maxval(expected_moments)))then
      message='stitched partition has missing or excess physical-grid coverage';return
    endif
    local_electrons=cell_volume*sum(partition_weight*density_values)
    call MPI_Allreduce(local_electrons,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.abs(electron_count-expected_electrons)>&
        electron_count_tolerance*max(1d0,abs(expected_electrons)))then
      message='stitched density does not preserve electron count';return
    endif
    allocate(srows(size(row_ids),nbasis),rhorows(size(row_ids),nbasis));srows=0d0;rhorows=0d0
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial_s(batch_count,nbasis),reduced_s(batch_count,nbasis),&
          partial_rho(batch_count,nbasis),reduced_rho(batch_count,nbasis))
        partial_s=0d0;partial_rho=0d0
        do p=1,size(physical_ids);do j=1,nbasis;do i=1,batch_count
          local_row=int(all_row_ids(row_displs(r+1)+batch_first+i-1))
          partial_s(i,j)=partial_s(i,j)+cell_volume*partition_weight(p)*&
            conjg(basis_values(local_row,p))*basis_values(j,p)
          partial_rho(i,j)=partial_rho(i,j)+cell_volume*partition_weight(p)*density_values(p)*&
            conjg(basis_values(local_row,p))*basis_values(j,p)
        enddo;enddo;enddo
        call MPI_Reduce(partial_s,reduced_s,batch_count*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        call MPI_Reduce(partial_rho,reduced_rho,batch_count*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(rank==r)then
          srows(batch_first:batch_first+batch_count-1,:)=reduced_s
          rhorows(batch_first:batch_first+batch_count-1,:)=reduced_rho
        endif
        peak_elements=max(peak_elements,int(2*size(srows)+2*size(partial_s)+2*size(reduced_s)+&
          2*nproc+2*nbasis,int64))
        deallocate(partial_s,reduced_s,partial_rho,reduced_rho)
      enddo
    enddo
    local_s_defect=0d0;local_rho_defect=0d0;scale_s=1d0;scale_rho=1d0
    if(size(srows)>0)then;scale_s=max(1d0,maxval(abs(srows)));scale_rho=max(1d0,maxval(abs(rhorows)));endif
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(block_s(batch_count,nbasis),block_rho(batch_count,nbasis))
        if(rank==r)then
          block_s=srows(batch_first:batch_first+batch_count-1,:)
          block_rho=rhorows(batch_first:batch_first+batch_count-1,:)
        endif
        call MPI_Bcast(block_s,batch_count*nbasis,MPI_DOUBLE_COMPLEX,r,comm,ierr)
        call MPI_Bcast(block_rho,batch_count*nbasis,MPI_DOUBLE_COMPLEX,r,comm,ierr)
        do local_row=1,size(row_ids);do i=1,batch_count
          j=int(all_row_ids(row_displs(r+1)+batch_first+i-1))
          local_s_defect=max(local_s_defect,abs(srows(local_row,j)-conjg(block_s(i,int(row_ids(local_row))))))
          local_rho_defect=max(local_rho_defect,&
            abs(rhorows(local_row,j)-conjg(block_rho(i,int(row_ids(local_row))))))
        enddo;enddo
        deallocate(block_s,block_rho)
      enddo
    enddo
    call MPI_Allreduce(local_s_defect,s_hermiticity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_rho_defect,rho_hermiticity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,scale_s,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,scale_rho,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.s_hermiticity>tolerance*scale_s.or.rho_hermiticity>tolerance*scale_rho)then
      message='stitched overlap or density tile is not Hermitian';return
    endif
    allocate(cholesky_rows(size(row_ids),nbasis),pivot_row(nbasis));cholesky_rows=0d0
    minimum_cholesky_pivot=huge(1d0);maximum_cholesky_pivot=0d0
    do k=1,nbasis
      pivot_owner=-1;pivot_local=0
      do r=0,nproc-1
        do i=1,row_counts(r+1)
          if(all_row_ids(row_displs(r+1)+i)==int(k,int64))then
            pivot_owner=r
            if(rank==r)pivot_local=i
          endif
        enddo
      enddo
      pivot_row=0d0;pivot_value=-huge(1d0)
      if(rank==pivot_owner)then
        if(k>1)pivot_row(1:k-1)=cholesky_rows(pivot_local,1:k-1)
        pivot_value=real(srows(pivot_local,k)-sum(pivot_row(1:k-1)*conjg(pivot_row(1:k-1))),real64)
        if(pivot_value>0d0)then
          pivot_value=sqrt(pivot_value);pivot_row(k)=cmplx(pivot_value,0d0,real64)
          cholesky_rows(pivot_local,k)=pivot_row(k)
        endif
      endif
      call MPI_Bcast(pivot_value,1,MPI_DOUBLE_PRECISION,pivot_owner,comm,ierr)
      call MPI_Bcast(pivot_row,nbasis,MPI_DOUBLE_COMPLEX,pivot_owner,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.pivot_value<=sqrt(tolerance*scale_s).or..not.ieee_is_finite(pivot_value))then
        message='stitched overlap lost positive-definite rank';return
      endif
      minimum_cholesky_pivot=min(minimum_cholesky_pivot,pivot_value)
      maximum_cholesky_pivot=max(maximum_cholesky_pivot,pivot_value)
      do local_row=1,size(row_ids)
        if(row_ids(local_row)<=int(k,int64))cycle
        cholesky_rows(local_row,k)=(srows(local_row,k)-&
          sum(cholesky_rows(local_row,1:k-1)*conjg(pivot_row(1:k-1))))/pivot_value
      enddo
    enddo
    pivot_condition=maximum_cholesky_pivot/minimum_cholesky_pivot
    peak_elements=max(peak_elements,int(2*size(srows)+size(cholesky_rows)+size(pivot_row)+&
      2*nproc+2*nbasis,int64))
    call MPI_Allreduce(MPI_IN_PLACE,peak_elements,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='stitched workspace receipt reduction failed';return;endif
    ok=.true.
#else
    ok=.false.;message='stitched overlap-density assembly requires MPI'
    electron_count=huge(1d0);s_hermiticity=huge(1d0);rho_hermiticity=huge(1d0)
    minimum_cholesky_pivot=0d0;pivot_condition=huge(1d0);peak_elements=0_int64
#endif
  end subroutine assemble_dg_stitched_overlap_density_rows

  subroutine assemble_dg_eigenexa_cyclic_metric_block(comm,nprow,npcol,myrow,mycol,&
      local_row_capacity,local_col_capacity,values,weights,&
      local_metric,peak_elements,ok,message)
    integer,intent(in)::comm,nprow,npcol,myrow,mycol,local_row_capacity,local_col_capacity
    complex(real64),intent(in)::values(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),allocatable,intent(out)::local_metric(:,:)
    integer(int64),intent(out)::peak_elements
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,norbital,norbital_min,norbital_max,nlocal,nrowlocal,ncollocal,&
      i,j,ilocal,jlocal,r,batch_first,batch_count
    integer::local_bad,global_bad
    integer,allocatable::rows(:),cols(:),coordinate_owner(:,:)
    real(real64)::scale,local_scale,imaginary_max,local_imaginary_max,gamma_tolerance
    real(real64),allocatable::local_block(:,:),global_block(:,:)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    norbital=size(values,1);nlocal=size(values,2);ok=.false.;message='';peak_elements=0_int64
    local_scale=1d0;local_imaginary_max=0d0
    if(size(values)>0)then
      local_scale=max(1d0,maxval(abs(real(values))))
      local_imaginary_max=maxval(abs(aimag(values)))
    endif
    call MPI_Allreduce(local_scale,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_imaginary_max,imaginary_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(norbital,norbital_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(norbital,norbital_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    gamma_tolerance=1024d0*epsilon(1d0)*scale
    local_bad=0
    if(norbital<=0.or.norbital_min/=norbital_max.or.size(weights)/=nlocal.or.&
       nprow<=0.or.npcol<=0.or.local_row_capacity<=0.or.local_col_capacity<=0)then
      local_bad=1
    else
      if(nprow>huge(nprow)/npcol.or.nprow*npcol/=nproc)local_bad=1
    endif
    if(myrow<1.or.myrow>nprow.or.mycol<1.or.mycol>npcol.or.any(weights<0d0).or.&
       .not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(real(values))).or.&
       .not.all(ieee_is_finite(aimag(values))).or.imaginary_max>gamma_tolerance)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='invalid Gamma-real EigenExa cyclic metric contract';return
    endif
    allocate(rows(nproc),cols(nproc),coordinate_owner(nprow,npcol));coordinate_owner=-1
    call MPI_Allgather(myrow,1,MPI_INTEGER,rows,1,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(mycol,1,MPI_INTEGER,cols,1,MPI_INTEGER,comm,ierr)
    local_bad=0
    do r=1,nproc
      if(rows(r)<1.or.rows(r)>nprow.or.cols(r)<1.or.cols(r)>npcol)then
        local_bad=1
      else if(coordinate_owner(rows(r),cols(r))/=-1)then
        local_bad=1
      else
        coordinate_owner(rows(r),cols(r))=r-1
      endif
    enddo
    if(any(coordinate_owner<0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='duplicate or missing EigenExa cyclic process coordinate';return
    endif
    nrowlocal=count([(mod(i-1,nprow)==myrow-1,i=1,norbital)])
    ncollocal=count([(mod(i-1,npcol)==mycol-1,i=1,norbital)])
    if(local_row_capacity<nrowlocal.or.local_col_capacity<ncollocal)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then
      message='EigenExa local cyclic capacity is smaller than owned metric block';return
    endif
    allocate(local_metric(local_row_capacity,local_col_capacity));local_metric=0d0
    peak_elements=int(size(local_metric),int64)
    do batch_first=1,norbital,row_batch_size
      batch_count=min(row_batch_size,norbital-batch_first+1)
      allocate(local_block(batch_count,norbital),global_block(batch_count,norbital))
      do i=1,batch_count
        local_block(i,:)=real(matmul(conjg(values(batch_first+i-1,:))*weights,transpose(values)))
      enddo
      call MPI_Allreduce(local_block,global_block,batch_count*norbital,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      do i=1,batch_count
        if(mod(batch_first+i-2,nprow)/=myrow-1)cycle
        ilocal=(batch_first+i-2)/nprow+1
        do j=1,norbital
          if(mod(j-1,npcol)/=mycol-1)cycle
          jlocal=(j-1)/npcol+1
          local_metric(ilocal,jlocal)=global_block(i,j)
        enddo
      enddo
      peak_elements=max(peak_elements,int(size(local_metric)+size(local_block)+size(global_block),int64))
      deallocate(local_block,global_block)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or..not.all(ieee_is_finite(local_metric)))then
      message='EigenExa cyclic metric tiled reduction failed';return
    endif
    ok=.true.
#else
    ok=.false.;message='EigenExa cyclic metric assembly requires MPI';peak_elements=0_int64
#endif
  end subroutine assemble_dg_eigenexa_cyclic_metric_block

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
