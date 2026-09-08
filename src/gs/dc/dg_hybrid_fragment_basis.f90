module dg_hybrid_fragment_basis
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private
  type,public::s_dg_hybrid_fragment_basis
    integer::fragment_id=0,generation=0
    integer(int64),allocatable::global_ids(:)
    integer(int64),allocatable::buffer_point_ids(:)
    integer,allocatable::sector(:)
    complex(real64),allocatable::buffer_values(:,:)
    integer(int64)::provenance_fingerprint=0_int64
  end type s_dg_hybrid_fragment_basis
  public::build_dg_hybrid_fragment_basis
contains
  subroutine build_dg_hybrid_fragment_basis(comm,fragment_id,wf_ids,wf_values,pw_ids,pw_values,&
      projector_radius,buffer_radius,basis,ok,message)
    integer,intent(in)::comm,fragment_id,projector_radius,buffer_radius
    integer(int64),intent(in)::wf_ids(:),pw_ids(:)
    complex(real64),intent(in)::wf_values(:,:),pw_values(:,:)
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,ierr,nproc,nlocal,ntotal,metadata(4),metadata_min(4),metadata_max(4)
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::local_ids(:),all_ids(:)
    integer(int64)::local_fingerprint,value_bits,global_fingerprint

    ok=.false.;message='';basis%fragment_id=0;basis%generation=0
    metadata=[fragment_id,projector_radius,buffer_radius,size(wf_values,1)]
    call MPI_Allreduce(metadata,metadata_min,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(metadata,metadata_max,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(metadata_min/=metadata_max))then
      message='fragment basis metadata differs between ranks';return
    endif
    if(fragment_id<=0.or.projector_radius<0.or.buffer_radius<0)then
      message='fragment basis metadata must be nonnegative';return
    endif
    if(projector_radius>buffer_radius)then
      message='projector support radius exceeds fragment buffer';return
    endif
    if(size(wf_values,2)/=size(wf_ids).or.size(pw_values,2)/=size(pw_ids).or.&
        size(pw_values,1)/=size(wf_values,1))then
      message='fragment basis ID/value shape mismatch';return
    endif

    nlocal=size(wf_ids)+size(pw_ids)
    allocate(local_ids(nlocal));local_ids=[wf_ids,pw_ids]
    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displacements(nproc))
    call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis count exchange failed';return;endif
    displacements(1)=0
    do i=2,nproc;displacements(i)=displacements(i-1)+counts(i-1);enddo
    ntotal=sum(counts);allocate(all_ids(ntotal))
    call MPI_Allgatherv(local_ids,nlocal,MPI_INTEGER8,all_ids,counts,displacements,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis ownership exchange failed';return;endif
    do i=1,ntotal
      if(all_ids(i)<=0_int64)then;message='fragment basis global IDs must be positive';return;endif
      do j=i+1,ntotal
        if(all_ids(i)==all_ids(j))then;message='fragment basis global ID has duplicate ownership';return;endif
      enddo
    enddo

    allocate(basis%global_ids(nlocal),basis%sector(nlocal),&
      basis%buffer_values(size(wf_values,1),nlocal),basis%buffer_point_ids(size(wf_values,1)))
    basis%fragment_id=fragment_id;basis%generation=1;basis%global_ids=local_ids
    basis%sector=[(1,i=1,size(wf_ids)),(2,i=1,size(pw_ids))]
    basis%buffer_point_ids=[(int(i,int64),i=1,size(wf_values,1))]
    if(size(wf_ids)>0)basis%buffer_values(:,:size(wf_ids))=wf_values
    if(size(pw_ids)>0)basis%buffer_values(:,size(wf_ids)+1:)=pw_values

    local_fingerprint=0_int64
    do i=1,nlocal
      local_fingerprint=ieor(local_fingerprint,ishftc(basis%global_ids(i),&
        modulo(7*int(basis%global_ids(i)),63)))
      local_fingerprint=ieor(local_fingerprint,ishftc(int(basis%sector(i),int64),&
        modulo(13*int(basis%global_ids(i)),63)))
      do j=1,size(basis%buffer_values,1)
        value_bits=transfer(real(basis%buffer_values(j,i),real64),value_bits)
        local_fingerprint=ieor(local_fingerprint,ishftc(value_bits,&
          modulo(11*int(basis%global_ids(i))+3*j,63)))
        value_bits=transfer(aimag(basis%buffer_values(j,i)),value_bits)
        local_fingerprint=ieor(local_fingerprint,ishftc(value_bits,&
          modulo(17*int(basis%global_ids(i))+5*j,63)))
      enddo
    enddo
    call MPI_Allreduce(local_fingerprint,global_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis fingerprint reduction failed';return;endif
    basis%provenance_fingerprint=ieor(global_fingerprint,ishftc(int(fragment_id,int64),29))
    ok=.true.;message=''
  end subroutine build_dg_hybrid_fragment_basis
end module dg_hybrid_fragment_basis
