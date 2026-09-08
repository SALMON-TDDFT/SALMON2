#include "config.h"
module dg_overlapping_wannier_operators
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_weak_operators,&
    assemble_dg_overlapping_wannier_weak_operator_rows,assemble_dg_stitched_weak_operator_rows
contains
  subroutine assemble_dg_stitched_weak_operator_rows(comm,nbasis,row_ids,physical_ids,&
      partition_weight,partition_gradient,basis_values,basis_gradients,local_potential,cell_volume,&
      kinetic_rows,potential_rows,weight_gradient_rows,kinetic_hermiticity,potential_hermiticity,&
      weight_gradient_trace,peak_elements,ok,message)
    integer,intent(in)::comm,nbasis
    integer(int64),intent(in)::row_ids(:),physical_ids(:)
    real(real64),intent(in)::partition_weight(:),partition_gradient(:,:),local_potential(:),cell_volume
    complex(real64),intent(in)::basis_values(:,:),basis_gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic_rows(:,:),potential_rows(:,:),weight_gradient_rows(:,:)
    real(real64),intent(out)::kinetic_hermiticity,potential_hermiticity
    real(real64),intent(out)::weight_gradient_trace
    integer(int64),intent(out)::peak_elements
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,total_rows,r,nrows,batch_first,batch_count,&
      i,j,p,row_index
    integer,allocatable::row_counts(:),row_displs(:)
    integer(int64),allocatable::all_row_ids(:),sorted_row_ids(:)
    complex(real64),allocatable::partial_t(:,:),partial_v(:,:),partial_w(:,:),reduced_t(:,:),&
      reduced_v(:,:),reduced_w(:,:),block_t(:,:),block_v(:,:)
    complex(real64)::weighted_value_i,weighted_value_j,weighted_gradient_i(3),weighted_gradient_j(3)
    real(real64)::sqrt_weight,local_t_defect,local_v_defect,t_scale,v_scale,local_weight_gradient_energy
    ok=.false.;message='';kinetic_hermiticity=huge(1d0);potential_hermiticity=huge(1d0)
    peak_elements=0_int64;local_bad=0
    weight_gradient_trace=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS.or.nbasis<1.or.cell_volume<=0d0.or.&
        size(partition_weight)/=size(physical_ids).or.size(local_potential)/=size(physical_ids).or.&
        any(shape(partition_gradient)/=[3,size(physical_ids)]).or.&
        any(shape(basis_values)/=[nbasis,size(physical_ids)]).or.&
        any(shape(basis_gradients)/=[3,nbasis,size(physical_ids)]).or.&
        any(row_ids<1_int64).or.any(row_ids>int(nbasis,int64)).or.any(physical_ids<1_int64).or.&
        any(partition_weight<0d0).or.any(partition_weight==0d0.and.&
        maxval(abs(partition_gradient),dim=1)>0d0))local_bad=1
    if(.not.all(ieee_is_finite(partition_weight)).or..not.all(ieee_is_finite(partition_gradient)).or.&
        .not.all(ieee_is_finite(local_potential)).or..not.all(ieee_is_finite(real(basis_values))).or.&
        .not.all(ieee_is_finite(aimag(basis_values))).or.&
        .not.all(ieee_is_finite(real(basis_gradients))).or.&
        .not.all(ieee_is_finite(aimag(basis_gradients))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='invalid stitched weak-operator contract';return;endif
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
      message='duplicate or missing stitched weak-operator row owner';return
    endif
    local_weight_gradient_energy=0d0
    do p=1,size(physical_ids)
      if(partition_weight(p)==0d0)cycle
      sqrt_weight=sqrt(partition_weight(p))
      do i=1,nbasis
        weighted_gradient_i=sqrt_weight*basis_gradients(:,i,p)+&
          0.5d0*partition_gradient(:,p)*basis_values(i,p)/sqrt_weight
        local_weight_gradient_energy=local_weight_gradient_energy+0.5d0*cell_volume*&
          (sum(abs(weighted_gradient_i)**2)-partition_weight(p)*sum(abs(basis_gradients(:,i,p))**2))
      enddo
    enddo
    call MPI_Allreduce(local_weight_gradient_energy,weight_gradient_trace,1,MPI_DOUBLE_PRECISION,&
      MPI_SUM,comm,ierr)
    allocate(kinetic_rows(size(row_ids),nbasis),potential_rows(size(row_ids),nbasis),&
      weight_gradient_rows(size(row_ids),nbasis))
    kinetic_rows=0d0;potential_rows=0d0;weight_gradient_rows=0d0
    peak_elements=int(size(kinetic_rows)+size(potential_rows)+size(weight_gradient_rows)+&
      2*nproc+2*nbasis+8,int64)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial_t(batch_count,nbasis),partial_v(batch_count,nbasis),&
          partial_w(batch_count,nbasis),reduced_t(batch_count,nbasis),&
          reduced_v(batch_count,nbasis),reduced_w(batch_count,nbasis))
        partial_t=0d0;partial_v=0d0;partial_w=0d0
        do p=1,size(physical_ids)
          if(partition_weight(p)==0d0)cycle
          sqrt_weight=sqrt(partition_weight(p))
          do j=1,nbasis
            weighted_value_j=sqrt_weight*basis_values(j,p)
            weighted_gradient_j=sqrt_weight*basis_gradients(:,j,p)+&
              0.5d0*partition_gradient(:,p)*basis_values(j,p)/sqrt_weight
            do i=1,batch_count
          row_index=int(all_row_ids(row_displs(r+1)+batch_first+i-1))
          weighted_value_i=sqrt_weight*basis_values(row_index,p)
          weighted_gradient_i=sqrt_weight*basis_gradients(:,row_index,p)+&
            0.5d0*partition_gradient(:,p)*basis_values(row_index,p)/sqrt_weight
          partial_t(i,j)=partial_t(i,j)+0.5d0*cell_volume*&
            sum(conjg(weighted_gradient_i)*weighted_gradient_j)
          partial_w(i,j)=partial_w(i,j)+0.5d0*cell_volume*&
            (sum(conjg(weighted_gradient_i)*weighted_gradient_j)-partition_weight(p)*&
            sum(conjg(basis_gradients(:,row_index,p))*basis_gradients(:,j,p)))
          partial_v(i,j)=partial_v(i,j)+cell_volume*local_potential(p)*&
            conjg(weighted_value_i)*weighted_value_j
            enddo
          enddo
        enddo
        call MPI_Reduce(partial_t,reduced_t,batch_count*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        call MPI_Reduce(partial_v,reduced_v,batch_count*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        call MPI_Reduce(partial_w,reduced_w,batch_count*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(rank==r)then
          kinetic_rows(batch_first:batch_first+batch_count-1,:)=reduced_t
          potential_rows(batch_first:batch_first+batch_count-1,:)=reduced_v
          weight_gradient_rows(batch_first:batch_first+batch_count-1,:)=reduced_w
        endif
        peak_elements=max(peak_elements,int(size(kinetic_rows)+size(potential_rows)+&
          size(weight_gradient_rows)+3*size(partial_t)+3*size(reduced_t)+&
          2*nproc+2*nbasis,int64))
        deallocate(partial_t,partial_v,partial_w,reduced_t,reduced_v,reduced_w)
      enddo
    enddo
    local_t_defect=0d0;local_v_defect=0d0;t_scale=1d0;v_scale=1d0
    if(size(kinetic_rows)>0)then
      t_scale=max(1d0,maxval(abs(kinetic_rows)));v_scale=max(1d0,maxval(abs(potential_rows)))
    endif
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(block_t(batch_count,nbasis),block_v(batch_count,nbasis))
        if(rank==r)then
          block_t=kinetic_rows(batch_first:batch_first+batch_count-1,:)
          block_v=potential_rows(batch_first:batch_first+batch_count-1,:)
        endif
        call MPI_Bcast(block_t,batch_count*nbasis,MPI_DOUBLE_COMPLEX,r,comm,ierr)
        call MPI_Bcast(block_v,batch_count*nbasis,MPI_DOUBLE_COMPLEX,r,comm,ierr)
        do j=1,size(row_ids);do i=1,batch_count
          row_index=int(all_row_ids(row_displs(r+1)+batch_first+i-1))
          local_t_defect=max(local_t_defect,&
            abs(kinetic_rows(j,row_index)-conjg(block_t(i,int(row_ids(j))))))
          local_v_defect=max(local_v_defect,&
            abs(potential_rows(j,row_index)-conjg(block_v(i,int(row_ids(j))))))
        enddo;enddo
        deallocate(block_t,block_v)
      enddo
    enddo
    call MPI_Allreduce(local_t_defect,kinetic_hermiticity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_v_defect,potential_hermiticity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,t_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,v_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,peak_elements,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.kinetic_hermiticity>1d-12*t_scale.or.&
        potential_hermiticity>1d-12*v_scale)then
      message='stitched weak operator is not Hermitian';return
    endif
    ok=.true.
#else
    ok=.false.;message='stitched weak operators require MPI';kinetic_hermiticity=huge(1d0)
    potential_hermiticity=huge(1d0);weight_gradient_trace=huge(1d0);peak_elements=0_int64
#endif
  end subroutine assemble_dg_stitched_weak_operator_rows

  subroutine assemble_dg_overlapping_wannier_weak_operator_rows(comm,nwann,row_ids,core_ids,weights,&
      values,gradients,local_potential,expected_core_count,kinetic_rows,potential_rows,&
      ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::row_ids(:),core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::values(:,:),gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic_rows(:,:),potential_rows(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,total_count,total_rows,p,i,j,r,nrows,&
      batch_first,batch_count,nwann_min,nwann_max
    integer,allocatable::core_counts(:),core_displs(:),row_counts(:),row_displs(:)
    integer(int64),allocatable::all_core_ids(:),all_row_ids(:),validation_ids(:)
    integer(int64)::expected_min,expected_max
    complex(real64),allocatable::partial_t(:,:),partial_v(:,:),reduced_t(:,:),reduced_v(:,:)
    logical::finite_payload,shapes_valid

    ok=.false.;message='';ownership_count=0;local_bad=0
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
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)local_bad=1
    shapes_valid=size(values,1)==nwann.and.size(values,2)==size(core_ids).and.&
      size(gradients,1)==3.and.size(gradients,2)==nwann.and.size(gradients,3)==size(core_ids)
    finite_payload=all(ieee_is_finite(weights)).and.all(ieee_is_finite(local_potential))
    if(nwann<=0.or.expected_core_count<=0_int64.or.size(weights)/=size(core_ids).or.&
        size(local_potential)/=size(core_ids).or..not.shapes_valid) local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)).or.any(core_ids<=0_int64).or.&
        any(weights<=0d0))local_bad=1
    if(shapes_valid)then
      finite_payload=finite_payload.and.all(ieee_is_finite(real(values))).and.&
        all(ieee_is_finite(aimag(values))).and.all(ieee_is_finite(real(gradients))).and.&
        all(ieee_is_finite(aimag(gradients)))
    endif
    if(.not.finite_payload)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid row-owned weak-operator payload';return;endif

    allocate(core_counts(nproc),core_displs(nproc),row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator ownership metadata collective failed';return;endif
    total_count=0;total_rows=0
    do r=1,nproc
      core_displs(r)=total_count;row_displs(r)=total_rows
      if(core_counts(r)<0.or.row_counts(r)<0.or.total_count>huge(total_count)-core_counts(r).or.&
          total_rows>huge(total_rows)-row_counts(r))local_bad=1
      if(local_bad==0)then
        total_count=total_count+core_counts(r);total_rows=total_rows+row_counts(r)
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_core_count.or.total_rows/=nwann)then
      message='missing or extra row/core owner in weak-operator assembly';return
    endif
    allocate(all_core_ids(total_count),all_row_ids(total_rows))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_core_ids,core_counts,core_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator ownership payload collective failed';return;endif
    allocate(validation_ids(max(total_count,total_rows)))
    validation_ids(1:total_count)=all_core_ids;call sort_ids(validation_ids(1:total_count))
    do i=1,total_count
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    validation_ids(1:total_rows)=all_row_ids;call sort_ids(validation_ids(1:total_rows))
    do i=1,total_rows
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='duplicate or missing row/core owner in weak-operator assembly';return;endif

    allocate(kinetic_rows(size(row_ids),nwann),potential_rows(size(row_ids),nwann))
    kinetic_rows=(0d0,0d0);potential_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial_t(batch_count,nwann),partial_v(batch_count,nwann),&
          reduced_t(batch_count,nwann),reduced_v(batch_count,nwann))
        partial_t=(0d0,0d0);partial_v=(0d0,0d0)
        do p=1,size(core_ids)
          do j=1,nwann;do i=1,batch_count
            partial_t(i,j)=partial_t(i,j)+0.5d0*weights(p)*sum(conjg(gradients(:,&
              int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*gradients(:,j,p))
            partial_v(i,j)=partial_v(i,j)+weights(p)*local_potential(p)*conjg(values(&
              int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*values(j,p)
          enddo;enddo
        enddo
        call MPI_Reduce(partial_t,reduced_t,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        call MPI_Reduce(partial_v,reduced_v,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==r)then
          kinetic_rows(batch_first:batch_first+batch_count-1,:)=reduced_t
          potential_rows(batch_first:batch_first+batch_count-1,:)=reduced_v
        endif
        deallocate(partial_t,partial_v,reduced_t,reduced_v)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator row reduction failed';return;endif
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='row-owned weak overlapping-Wannier operators require MPI'
    ownership_count=0
#endif
  end subroutine

  subroutine assemble_dg_overlapping_wannier_weak_operators(comm,nwann,core_ids,weights,values,&
      gradients,local_potential,expected_core_count,kinetic,potential,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::values(:,:),gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic(:,:),potential(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,matrix_count,nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64,expected_min,expected_max
    complex(real64),allocatable::local_kinetic(:,:),local_potential_matrix(:,:)
    logical::finite_payload,shapes_valid
    real(real64)::hermiticity_defect,scale
    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)then
      message='inconsistent weak-operator assembly contract across ranks';return
    endif
    finite_payload=all(ieee_is_finite(weights)).and.all(ieee_is_finite(local_potential))
    if(nwann<=0.or.expected_core_count<=0_int64)local_bad=1
    if(size(weights)/=size(core_ids).or.size(local_potential)/=size(core_ids))local_bad=1
    shapes_valid=size(values,1)==nwann.and.size(values,2)==size(core_ids).and.&
      size(gradients,1)==3.and.size(gradients,2)==nwann.and.size(gradients,3)==size(core_ids)
    if(.not.shapes_valid)local_bad=1
    if(shapes_valid)then
      do p=1,size(values,2);do i=1,size(values,1)
        finite_payload=finite_payload.and.ieee_is_finite(real(values(i,p))).and.&
          ieee_is_finite(aimag(values(i,p)))
        finite_payload=finite_payload.and.all(ieee_is_finite(real(gradients(:,i,p)))).and.&
          all(ieee_is_finite(aimag(gradients(:,i,p))))
      enddo;enddo
    endif
    if(any(core_ids<=0_int64).or.any(weights<=0d0).or..not.finite_payload)local_bad=1
    if(nwann>0.and.int(nwann,int64)<=huge(1_int64)/int(nwann,int64))then
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    else
      matrix_count64=0_int64;local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid weak-operator unique-core payload';return;endif
    matrix_count=int(matrix_count64)

    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displs(nproc))
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
      message='missing or extra weak-operator core owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_ids,counts,displs,MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then
        message='duplicate or missing weak-operator core owner';return
      endif
    enddo

    allocate(local_kinetic(nwann,nwann),local_potential_matrix(nwann,nwann))
    local_kinetic=(0d0,0d0);local_potential_matrix=(0d0,0d0)
    do p=1,size(core_ids)
      do j=1,nwann;do i=1,nwann
        local_kinetic(i,j)=local_kinetic(i,j)+0.5d0*weights(p)*&
          sum(conjg(gradients(:,i,p))*gradients(:,j,p))
        local_potential_matrix(i,j)=local_potential_matrix(i,j)+weights(p)*local_potential(p)*&
          conjg(values(i,p))*values(j,p)
      enddo;enddo
    enddo
    allocate(kinetic(nwann,nwann),potential(nwann,nwann))
    call MPI_Allreduce(local_kinetic,kinetic,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_potential_matrix,potential,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale=max(1d0,max(maxval(abs(kinetic)),maxval(abs(potential))))
    hermiticity_defect=max(maxval(abs(kinetic-conjg(transpose(kinetic)))),&
      maxval(abs(potential-conjg(transpose(potential)))))
    if(hermiticity_defect>1d-12*scale)then
      message='weak-operator Hermiticity defect exceeds tolerance';return
    endif
    kinetic=0.5d0*(kinetic+conjg(transpose(kinetic)))
    potential=0.5d0*(potential+conjg(transpose(potential)))
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='weak overlapping-Wannier operators require MPI';ownership_count=0
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
end module dg_overlapping_wannier_operators
