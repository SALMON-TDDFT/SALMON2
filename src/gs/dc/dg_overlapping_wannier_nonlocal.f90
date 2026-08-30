#include "config.h"
module dg_overlapping_wannier_nonlocal
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_nonlocal,assemble_dg_overlapping_wannier_nonlocal_rows
  public::collect_dg_overlapping_wannier_projector_overlaps
  public::apply_dg_overlapping_wannier_nonlocal_action
contains
  subroutine apply_dg_overlapping_wannier_nonlocal_action(comm,nglobal,ncore,core_positions,&
      projector_positions,projector_values,action_strength,overlap,action,ok,message)
    integer,intent(in)::comm,nglobal,ncore,core_positions(:),projector_positions(:)
    complex(real64),intent(in)::projector_values(:),overlap(:,:)
    real(real64),intent(in)::action_strength(:)
    complex(real64),allocatable,intent(out)::action(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,basis,ierr,local_bad,global_bad
    logical::shape_ok
    ok=.false.;message='';local_bad=0
    shape_ok=size(core_positions)==size(projector_positions).and.&
      size(core_positions)==size(projector_values).and.size(overlap,1)==nglobal.and.&
      size(overlap,2)==size(action_strength)
    if(nglobal<1.or.ncore<1.or..not.shape_ok)local_bad=1
    if(shape_ok)then
      if(any(core_positions<1).or.any(core_positions>ncore).or.any(projector_positions<1).or.&
          any(projector_positions>size(action_strength)).or.any(.not.ieee_is_finite(action_strength)).or.&
          .not.all(ieee_is_finite(real(projector_values))).or.&
          .not.all(ieee_is_finite(aimag(projector_values))).or.&
          .not.all(ieee_is_finite(real(overlap))).or..not.all(ieee_is_finite(aimag(overlap))))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid or missing nonlocal projector support';return;endif
    allocate(action(nglobal,ncore));action=(0d0,0d0)
    do i=1,size(core_positions)
      do basis=1,nglobal
        action(basis,core_positions(i))=action(basis,core_positions(i))+&
          action_strength(projector_positions(i))*projector_values(i)*overlap(basis,projector_positions(i))
      enddo
    enddo
    ok=.true.
#else
    ok=.false.;message='overlapping-Wannier nonlocal action requires MPI'
#endif
  end subroutine apply_dg_overlapping_wannier_nonlocal_action

  subroutine collect_dg_overlapping_wannier_projector_overlaps(comm,nwann,atom_ids,ordinals,&
      matrix_strength,action_strength,partial_overlap,projector_ids,owned_matrix_strength,owned_overlap,&
      expected_projector_count,ok,message,complete_atom_ids,complete_ordinals,complete_matrix_strength,&
      complete_action_strength,complete_overlap)
    integer,intent(in)::comm,nwann,atom_ids(:),ordinals(:)
    real(real64),intent(in)::matrix_strength(:),action_strength(:)
    complex(real64),intent(in)::partial_overlap(:,:)
    integer(int64),allocatable,intent(out)::projector_ids(:)
    real(real64),allocatable,intent(out)::owned_matrix_strength(:)
    complex(real64),allocatable,intent(out)::owned_overlap(:,:)
    integer,intent(out)::expected_projector_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable,optional,intent(out)::complete_atom_ids(:),complete_ordinals(:)
    real(real64),allocatable,optional,intent(out)::complete_matrix_strength(:),complete_action_strength(:)
    complex(real64),allocatable,optional,intent(out)::complete_overlap(:,:)
#ifdef USE_MPI
    integer::rank,nproc,ierr,r,p,q,total_records,nowned,unique_count,local_bad,global_bad
    integer,allocatable::counts(:),displacements(:),complex_counts(:),complex_displacements(:),&
      all_atom_ids(:),all_ordinals(:),unique_ids(:),owner_ranks(:)
    real(real64),allocatable::all_matrix_strength(:),all_action_strength(:)
    complex(real64),allocatable::all_overlap(:,:)
    logical::matched
    real(real64)::reference_action_strength
    ok=.false.;message='';expected_projector_count=0;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(nwann<1.or.size(atom_ids)/=size(ordinals).or.size(atom_ids)/=size(matrix_strength).or.&
        size(atom_ids)/=size(action_strength).or.&
        size(partial_overlap,1)/=nwann.or.&
        size(partial_overlap,2)/=size(atom_ids))local_bad=1
    if(local_bad==0)then
      if(any(atom_ids<1).or.any(ordinals<1).or.&
          any(.not.ieee_is_finite(matrix_strength)).or.any(.not.ieee_is_finite(action_strength)).or.&
          .not.all(ieee_is_finite(real(partial_overlap))).or.&
          .not.all(ieee_is_finite(aimag(partial_overlap))))local_bad=1
    end if
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid fragment projector overlap payload';return;end if
    allocate(counts(nproc),displacements(nproc),complex_counts(nproc),complex_displacements(nproc))
    call MPI_Allgather(size(atom_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_records=0
    do r=1,nproc
      displacements(r)=total_records;total_records=total_records+counts(r)
      complex_counts(r)=nwann*counts(r);complex_displacements(r)=nwann*displacements(r)
    end do
    allocate(all_atom_ids(total_records),all_ordinals(total_records),all_matrix_strength(total_records),&
      all_action_strength(total_records),&
      all_overlap(nwann,total_records),unique_ids(total_records),owner_ranks(total_records))
    call MPI_Allgatherv(atom_ids,size(atom_ids),MPI_INTEGER,all_atom_ids,counts,displacements,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgatherv(ordinals,size(ordinals),MPI_INTEGER,all_ordinals,counts,displacements,&
      MPI_INTEGER,comm,ierr)
    call MPI_Allgatherv(matrix_strength,size(matrix_strength),MPI_DOUBLE_PRECISION,all_matrix_strength,&
      counts,displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    call MPI_Allgatherv(action_strength,size(action_strength),MPI_DOUBLE_PRECISION,all_action_strength,&
      counts,displacements,&
      MPI_DOUBLE_PRECISION,comm,ierr)
    call MPI_Allgatherv(partial_overlap,size(partial_overlap),MPI_DOUBLE_COMPLEX,all_overlap,&
      complex_counts,complex_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    unique_count=0
    do p=1,total_records
      unique_ids(p)=0
      do q=1,p-1
        if(all_atom_ids(q)==all_atom_ids(p).and.all_ordinals(q)==all_ordinals(p))then
          unique_ids(p)=unique_ids(q);exit
        end if
      end do
      if(unique_ids(p)==0)then
        unique_count=unique_count+1;unique_ids(p)=unique_count
        owner_ranks(unique_count)=0
        do r=1,nproc
          if(p>displacements(r).and.p<=displacements(r)+counts(r))then
            owner_ranks(unique_count)=r-1;exit
          end if
        end do
      end if
    end do
    expected_projector_count=unique_count;nowned=count(owner_ranks(1:unique_count)==rank)
    if(present(complete_atom_ids))then
      allocate(complete_atom_ids(unique_count));complete_atom_ids=0
    endif
    if(present(complete_ordinals))then
      allocate(complete_ordinals(unique_count));complete_ordinals=0
    endif
    if(present(complete_matrix_strength))then
      allocate(complete_matrix_strength(unique_count));complete_matrix_strength=0d0
    endif
    if(present(complete_action_strength))then
      allocate(complete_action_strength(unique_count));complete_action_strength=0d0
    endif
    if(present(complete_overlap))then
      allocate(complete_overlap(nwann,unique_count));complete_overlap=(0d0,0d0)
    endif
    do p=1,total_records
      q=unique_ids(p)
      if(present(complete_atom_ids))complete_atom_ids(q)=all_atom_ids(p)
      if(present(complete_ordinals))complete_ordinals(q)=all_ordinals(p)
      if(present(complete_matrix_strength))complete_matrix_strength(q)=all_matrix_strength(p)
      if(present(complete_action_strength))complete_action_strength(q)=all_action_strength(p)
      if(present(complete_overlap))complete_overlap(:,q)=complete_overlap(:,q)+all_overlap(:,p)
    enddo
    allocate(projector_ids(nowned),owned_matrix_strength(nowned),owned_overlap(nwann,nowned))
    owned_overlap=(0d0,0d0);nowned=0
    do q=1,unique_count
      if(owner_ranks(q)/=rank)cycle
      nowned=nowned+1;projector_ids(nowned)=int(q,int64);matched=.false.
      do p=1,total_records
        if(unique_ids(p)/=q)cycle
        if(.not.matched)then
          owned_matrix_strength(nowned)=all_matrix_strength(p)
          reference_action_strength=all_action_strength(p);matched=.true.
        else if(abs(all_matrix_strength(p)-owned_matrix_strength(nowned))>&
            1024d0*epsilon(1d0)*max(1d0,abs(owned_matrix_strength(nowned))).or.&
            abs(all_action_strength(p)-reference_action_strength)>&
            1024d0*epsilon(1d0)*max(1d0,abs(reference_action_strength)))then
          local_bad=1
        end if
        owned_overlap(:,nowned)=owned_overlap(:,nowned)+all_overlap(:,p)
      end do
    end do
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='fragment copies disagree on nonlocal projector strength';return;end if
    ok=.true.;message=''
#else
    ok=.false.;message='fragment projector overlap collection requires MPI';expected_projector_count=0
#endif
  end subroutine collect_dg_overlapping_wannier_projector_overlaps

  subroutine assemble_dg_overlapping_wannier_nonlocal_rows(comm,nwann,row_ids,projector_ids,strength,&
      overlap,complete_tail_overlap,expected_projector_count,matrix_rows,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::row_ids(:),projector_ids(:),expected_projector_count
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::overlap(:,:)
    logical,intent(in)::complete_tail_overlap(:,:)
    complex(real64),allocatable,intent(out)::matrix_rows(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,total_projectors,total_rows,r,i,j,p,nrows,&
      batch_first,batch_count,nwann_min,nwann_max
    integer,allocatable::projector_counts(:),projector_displs(:),row_counts(:),row_displs(:)
    integer(int64),allocatable::all_projector_ids(:),all_row_ids(:),validation_ids(:)
    integer(int64)::expected_min,expected_max
    complex(real64),allocatable::partial(:,:),reduced(:,:)
    logical::shape_ok

    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_projector_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_projector_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)local_bad=1
    shape_ok=size(strength)==size(projector_ids).and.size(overlap,1)==nwann.and.&
      size(overlap,2)==size(projector_ids).and.size(complete_tail_overlap,1)==nwann.and.&
      size(complete_tail_overlap,2)==size(projector_ids)
    if(nwann<=0.or.expected_projector_count<0_int64.or..not.shape_ok)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)).or.any(projector_ids<=0_int64).or.&
        any(.not.ieee_is_finite(strength)))local_bad=1
    if(shape_ok)then
      if(.not.all(complete_tail_overlap).or..not.all(ieee_is_finite(real(overlap))).or.&
          .not.all(ieee_is_finite(aimag(overlap))))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid row-owned nonlocal payload';return;endif
    allocate(projector_counts(nproc),projector_displs(nproc),row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(projector_ids),1,MPI_INTEGER,projector_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal ownership metadata collective failed';return;endif
    total_projectors=0;total_rows=0
    do r=1,nproc
      projector_displs(r)=total_projectors;row_displs(r)=total_rows
      if(projector_counts(r)<0.or.row_counts(r)<0.or.&
          total_projectors>huge(total_projectors)-projector_counts(r).or.&
          total_rows>huge(total_rows)-row_counts(r))local_bad=1
      if(local_bad==0)then
        total_projectors=total_projectors+projector_counts(r);total_rows=total_rows+row_counts(r)
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_projectors,int64)/=expected_projector_count.or.total_rows/=nwann)then
      message='missing or extra row/projector owner in nonlocal assembly';return
    endif
    allocate(all_projector_ids(total_projectors),all_row_ids(total_rows))
    call MPI_Allgatherv(projector_ids,size(projector_ids),MPI_INTEGER8,all_projector_ids,&
      projector_counts,projector_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal ownership payload collective failed';return;endif
    allocate(validation_ids(max(total_projectors,total_rows)))
    validation_ids(1:total_projectors)=all_projector_ids
    call sort_ids(validation_ids(1:total_projectors))
    do i=1,total_projectors
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    validation_ids(1:total_rows)=all_row_ids
    call sort_ids(validation_ids(1:total_rows))
    do i=1,total_rows
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='duplicate or missing row/projector owner in nonlocal assembly';return;endif
    allocate(matrix_rows(size(row_ids),nwann));matrix_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial(batch_count,nwann),reduced(batch_count,nwann));partial=(0d0,0d0)
        do p=1,size(projector_ids);do j=1,nwann;do i=1,batch_count
          partial(i,j)=partial(i,j)+strength(p)*conjg(overlap(&
            int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*overlap(j,p)
        enddo;enddo;enddo
        call MPI_Reduce(partial,reduced,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==r)matrix_rows(batch_first:batch_first+batch_count-1,:)=reduced
        deallocate(partial,reduced)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal row reduction failed';return;endif
    ownership_count=total_projectors;ok=.true.
#else
    ok=.false.;message='row-owned overlapping-Wannier nonlocal assembly requires MPI'
    ownership_count=0
#endif
  end subroutine

  subroutine assemble_dg_overlapping_wannier_nonlocal(comm,nwann,projector_ids,strength,overlap,&
      complete_tail_overlap,expected_projector_count,matrix,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::projector_ids(:),expected_projector_count
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::overlap(:,:)
    logical,intent(in)::complete_tail_overlap(:,:)
    complex(real64),allocatable,intent(out)::matrix(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,matrix_count,nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64,expected_min,expected_max
    complex(real64),allocatable::local_matrix(:,:)
    logical::shape_ok,finite_overlap
    real(real64)::hermiticity_defect,scale
    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_projector_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_projector_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)then
      message='inconsistent nonlocal assembly contract across ranks';return
    endif
    shape_ok=size(strength)==size(projector_ids).and.size(overlap,1)==nwann.and.&
      size(overlap,2)==size(projector_ids).and.size(complete_tail_overlap,1)==nwann.and.&
      size(complete_tail_overlap,2)==size(projector_ids)
    if(nwann<=0.or.expected_projector_count<=0_int64.or..not.shape_ok)local_bad=1
    if(any(projector_ids<=0_int64).or.any(.not.ieee_is_finite(strength)))local_bad=1
    finite_overlap=.true.
    if(shape_ok)then
      do p=1,size(projector_ids);do i=1,nwann
        finite_overlap=finite_overlap.and.ieee_is_finite(real(overlap(i,p))).and.&
          ieee_is_finite(aimag(overlap(i,p)))
      enddo;enddo
      if(.not.all(complete_tail_overlap))local_bad=1
    endif
    if(.not.finite_overlap)local_bad=1
    if(nwann>0.and.int(nwann,int64)<=huge(1_int64)/int(nwann,int64))then
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    else
      matrix_count64=0_int64;local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid or incomplete tail-projector overlaps';return;endif
    matrix_count=int(matrix_count64)
    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displs(nproc))
    call MPI_Allgather(size(projector_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_count=0;displs(1)=0
    do i=1,nproc
      if(counts(i)<0.or.total_count>huge(total_count)-counts(i))then;local_bad=1;exit;endif
      if(i>1)displs(i)=total_count
      total_count=total_count+counts(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_projector_count)then
      message='missing or extra atom/projector owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(projector_ids,size(projector_ids),MPI_INTEGER8,all_ids,counts,displs,&
      MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then;message='duplicate or missing atom/projector owner';return;endif
    enddo
    allocate(local_matrix(nwann,nwann),matrix(nwann,nwann));local_matrix=(0d0,0d0)
    do p=1,size(projector_ids);do j=1,nwann;do i=1,nwann
      local_matrix(i,j)=local_matrix(i,j)+strength(p)*conjg(overlap(i,p))*overlap(j,p)
    enddo;enddo;enddo
    call MPI_Allreduce(local_matrix,matrix,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale=max(1d0,maxval(abs(matrix)))
    hermiticity_defect=maxval(abs(matrix-conjg(transpose(matrix))))
    if(hermiticity_defect>1d-12*scale)then
      message='nonlocal Hermiticity defect exceeds tolerance';return
    endif
    matrix=0.5d0*(matrix+conjg(transpose(matrix)))
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='overlapping-Wannier nonlocal assembly requires MPI';ownership_count=0
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
end module dg_overlapping_wannier_nonlocal
