module dg_hybrid_fragment_basis_stream
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  private
  type,public::s_dg_hybrid_fragment_basis_stream
    logical::initialized=.false.,finalized=.false.
    integer::fragment_id=0,nwannier=0,npw=0,npoint=0
    integer(int64)::maximum_tile_elements=0_int64
    integer,allocatable::pw_owner(:),pw_slot(:)
    logical,allocatable::pw_filled(:)
  end type
  public::initialize_dg_hybrid_fragment_basis_stream,append_dg_hybrid_projected_pw_tile,&
    finalize_dg_hybrid_fragment_basis_stream
contains
  subroutine initialize_dg_hybrid_fragment_basis_stream(comm,fragment_count,fragment_id,point_ids,&
      wannier_buffer,wannier_owner,pw_owner,stream,basis,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,fragment_count,fragment_id,wannier_owner(:),pw_owner(:)
    integer(int64),intent(in)::point_ids(:)
    complex(real64),intent(in)::wannier_buffer(:,:)
    type(s_dg_hybrid_fragment_basis_stream),intent(out)::stream
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,slot,nlocal,local_bad,global_bad,ierr,minimum,maximum,allocation_status
    integer,allocatable::fragment_presence(:),ownership(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    call agree_integer(fragment_count,minimum,maximum,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='inconsistent stream fragment count';return;endif
    call agree_integer(size(wannier_owner),minimum,maximum,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='inconsistent stream Wannier count';return;endif
    call agree_integer(size(pw_owner),minimum,maximum,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='inconsistent stream PW count';return;endif
    local_bad=merge(0,1,fragment_count>0.and.fragment_id>=0.and.fragment_id<=fragment_count.and.&
      size(wannier_owner)>0.and.size(pw_owner)>0.and.size(wannier_buffer,1)==size(wannier_owner).and.&
      size(wannier_buffer,2)==size(point_ids).and.all(wannier_owner>=1).and.&
      all(wannier_owner<=fragment_count).and.all(pw_owner>=1).and.all(pw_owner<=fragment_count).and.&
      all(point_ids>0_int64).and.finite_complex(wannier_buffer))
    if(fragment_id==0.and.size(point_ids)/=0)local_bad=1
    do i=1,size(wannier_owner)
      call agree_integer(wannier_owner(i),minimum,maximum,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum/=maximum)local_bad=1
    enddo
    do i=1,size(pw_owner)
      call agree_integer(pw_owner(i),minimum,maximum,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum/=maximum)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid fragment basis stream input';return;endif
    allocate(fragment_presence(fragment_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);global_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='cannot allocate stream fragment-presence workspace';return
    endif
    fragment_presence=0
    if(fragment_id>0)fragment_presence(fragment_id)=1
    call MPI_Allreduce(MPI_IN_PLACE,fragment_presence,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(fragment_presence/=1))then
      message='duplicate or missing stream fragment owner';return
    endif
    nlocal=count(wannier_owner==fragment_id)+count(pw_owner==fragment_id)
    allocate(basis%global_ids(nlocal),basis%sector(nlocal),basis%buffer_point_ids(size(point_ids)),&
      basis%buffer_values(size(point_ids),nlocal),stream%pw_owner(size(pw_owner)),&
      stream%pw_slot(size(pw_owner)),stream%pw_filled(size(pw_owner)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);global_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call clear_fragment_basis_stream_outputs(stream,basis)
      message='cannot allocate fragment basis stream payload';return
    endif
    basis%fragment_id=fragment_id;basis%generation=1;basis%buffer_point_ids=point_ids
    basis%buffer_values=(0d0,0d0);stream%pw_owner=pw_owner;stream%pw_slot=0;stream%pw_filled=.false.
    slot=0
    do i=1,size(wannier_owner)
      if(wannier_owner(i)/=fragment_id)cycle
      slot=slot+1;basis%global_ids(slot)=int(i,int64);basis%sector(slot)=1
      basis%buffer_values(:,slot)=wannier_buffer(i,:)
    enddo
    do i=1,size(pw_owner)
      if(pw_owner(i)/=fragment_id)cycle
      slot=slot+1;basis%global_ids(slot)=int(size(wannier_owner)+i,int64);basis%sector(slot)=2
      stream%pw_slot(i)=slot
    enddo
    allocate(ownership(size(wannier_owner)+size(pw_owner)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);global_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call clear_fragment_basis_stream_outputs(stream,basis)
      message='cannot allocate stream basis-ownership workspace';return
    endif
    ownership=0
    do i=1,nlocal;ownership(int(basis%global_ids(i)))=ownership(int(basis%global_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,size(ownership),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then;message='stream basis IDs are not owned once';return;endif
    stream%initialized=.true.;stream%fragment_id=fragment_id;stream%nwannier=size(wannier_owner)
    stream%npw=size(pw_owner);stream%npoint=size(point_ids)
    fingerprint=int(z'1F83D9ABFB41BD6B',int64)
    do i=1,size(wannier_owner);fingerprint=ieor(ishftc(fingerprint,7),int(wannier_owner(i),int64));enddo
    do i=1,size(pw_owner);fingerprint=ieor(ishftc(fingerprint,7),int(pw_owner(i),int64));enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(fragment_presence,ownership);ok=.true.
  end subroutine initialize_dg_hybrid_fragment_basis_stream

  subroutine clear_fragment_basis_stream_outputs(stream,basis)
    type(s_dg_hybrid_fragment_basis_stream),intent(inout)::stream
    type(s_dg_hybrid_fragment_basis),intent(inout)::basis
    if(allocated(stream%pw_owner))deallocate(stream%pw_owner)
    if(allocated(stream%pw_slot))deallocate(stream%pw_slot)
    if(allocated(stream%pw_filled))deallocate(stream%pw_filled)
    if(allocated(basis%global_ids))deallocate(basis%global_ids)
    if(allocated(basis%buffer_point_ids))deallocate(basis%buffer_point_ids)
    if(allocated(basis%sector))deallocate(basis%sector)
    if(allocated(basis%buffer_values))deallocate(basis%buffer_values)
    stream%initialized=.false.;stream%finalized=.false.;stream%fragment_id=0
    stream%nwannier=0;stream%npw=0;stream%npoint=0;stream%maximum_tile_elements=0_int64
    basis%fragment_id=0;basis%generation=0;basis%provenance_fingerprint=0_int64
  end subroutine clear_fragment_basis_stream_outputs

  subroutine append_dg_hybrid_projected_pw_tile(stream,first_column,projected_pw_tile,basis,ok,message)
    type(s_dg_hybrid_fragment_basis_stream),intent(inout)::stream
    integer,intent(in)::first_column
    complex(real64),intent(in)::projected_pw_tile(:,:)
    type(s_dg_hybrid_fragment_basis),intent(inout)::basis
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::j,column,slot
    ok=.false.;message=''
    if(.not.stream%initialized.or.stream%finalized.or.first_column<1.or.size(projected_pw_tile,1)<1.or.&
        first_column+size(projected_pw_tile,1)-1>stream%npw.or.&
        size(projected_pw_tile,2)/=stream%npoint.or..not.finite_complex(projected_pw_tile))then
      message='invalid projected PW stream tile';return
    endif
    if(any(stream%pw_filled(first_column:first_column+size(projected_pw_tile,1)-1)))then
      message='duplicate projected PW stream tile';return
    endif
    do j=1,size(projected_pw_tile,1)
      column=first_column+j-1;stream%pw_filled(column)=.true.;slot=stream%pw_slot(column)
      if(slot>0)basis%buffer_values(:,slot)=projected_pw_tile(j,:)
    enddo
    stream%maximum_tile_elements=max(stream%maximum_tile_elements,int(size(projected_pw_tile),int64))
    ok=.true.
  end subroutine append_dg_hybrid_projected_pw_tile

  subroutine finalize_dg_hybrid_fragment_basis_stream(comm,stream,basis,workspace_peak_bytes,&
      fingerprint,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis_stream),intent(inout)::stream
    type(s_dg_hybrid_fragment_basis),intent(inout)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,ierr,local_bad,global_bad
    integer(int64)::local_hash,bits
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    local_bad=merge(0,1,stream%initialized.and..not.stream%finalized.and.all(stream%pw_filled).and.&
      allocated(basis%global_ids).and.allocated(basis%buffer_values).and.finite_complex(basis%buffer_values))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='incomplete fragment basis stream';return;endif
    local_hash=int(z'5BE0CD19137E2179',int64)
    do j=1,size(basis%global_ids)
      local_hash=ieor(ishftc(local_hash,7),basis%global_ids(j))
      do i=1,size(basis%buffer_values,1)
        bits=transfer(real(basis%buffer_values(i,j)),bits);local_hash=ieor(ishftc(local_hash,7),bits)
        bits=transfer(aimag(basis%buffer_values(i,j)),bits);local_hash=ieor(ishftc(local_hash,11),bits)
      enddo
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment stream fingerprint reduction failed';return;endif
    if(fingerprint==0_int64)fingerprint=1_int64
    basis%provenance_fingerprint=fingerprint;stream%finalized=.true.
    workspace_peak_bytes=16_int64*stream%maximum_tile_elements;ok=.true.
  end subroutine finalize_dg_hybrid_fragment_basis_stream

  subroutine agree_integer(value,minimum,maximum,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum,maximum,ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex
end module dg_hybrid_fragment_basis_stream
