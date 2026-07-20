module dg_wpw_occupied_w_basis
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_occupied_w_basis
    logical::valid=.false.
    integer::local_fragment=0,local_count=0,global_count=0,epoch=-1
    integer(int64)::fingerprint=0_int64
    real(8)::source_condition=huge(1d0)
    integer::buffer_lo(3)=0,buffer_hi(3)=-1
    integer,allocatable::owned_ids(:),stable_keys(:,:),core_grid_ids(:)
    complex(8),allocatable::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
  end type
  public::initialize_dg_wpw_occupied_w_basis
  public::initialize_dg_wpw_occupied_w_basis_collective
contains
  subroutine initialize_dg_wpw_occupied_w_basis_collective(basis,comm,local_fragment,local_keys,&
      core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,buffer_gradients,&
      source_condition,epoch,expected_global_count,info)
    use mpi,only:MPI_Comm_size,MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_INTEGER,MPI_MAX,&
      MPI_SUM,MPI_SUCCESS
    type(s_dg_wpw_occupied_w_basis),intent(inout)::basis
    integer,intent(in)::comm,local_fragment,local_keys(:,:),core_grid_ids(:),buffer_lo(3),buffer_hi(3)
    complex(8),intent(in)::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
    real(8),intent(in)::source_condition
    integer,intent(in)::epoch,expected_global_count
    integer,intent(out)::info
    type(s_dg_wpw_occupied_w_basis)::candidate
    integer::nrank,ierr,local_bad,global_bad,nlocal,nglobal,i
    integer,allocatable::counts(:),displacements(:),key_counts(:),key_displacements(:)
    integer,allocatable::flat_keys(:),global_keys(:,:),local_owners(:),global_owners(:)

    info=1;local_bad=0;nlocal=size(local_keys,2)
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(local_fragment<=0.or.nlocal<=0.or.size(local_keys,1)/=5.or.expected_global_count<=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(counts(nrank),displacements(nrank),key_counts(nrank),key_displacements(nrank))
    call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    displacements(1)=0
    do i=2,nrank;displacements(i)=displacements(i-1)+counts(i-1);enddo
    nglobal=sum(counts)
    key_counts=5*counts;key_displacements=5*displacements
    allocate(flat_keys(5*nglobal),global_keys(5,nglobal),local_owners(nlocal),global_owners(nglobal))
    call MPI_Allgatherv(local_keys,5*nlocal,MPI_INTEGER,flat_keys,key_counts,key_displacements,&
      MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    global_keys=reshape(flat_keys,[5,nglobal]);local_owners=local_fragment
    call MPI_Allgatherv(local_owners,nlocal,MPI_INTEGER,global_owners,counts,displacements,&
      MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,nglobal==expected_global_count)
    candidate=basis
    if(local_bad==0)call initialize_dg_wpw_occupied_w_basis(candidate,local_fragment,global_keys,&
      global_owners,local_keys,core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,&
      buffer_gradients,source_condition,epoch,local_bad)
    local_bad=merge(0,1,local_bad==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    basis=candidate;info=0
  end subroutine initialize_dg_wpw_occupied_w_basis_collective

  subroutine initialize_dg_wpw_occupied_w_basis(basis,local_fragment,global_keys,global_owners,&
      local_keys,core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,buffer_gradients,&
      source_condition,epoch,info)
    type(s_dg_wpw_occupied_w_basis),intent(inout)::basis
    integer,intent(in)::local_fragment,global_keys(:,:),global_owners(:),local_keys(:,:),core_grid_ids(:)
    integer,intent(in)::buffer_lo(3),buffer_hi(3),epoch
    complex(8),intent(in)::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
    real(8),intent(in)::source_condition
    integer,intent(out)::info
    type(s_dg_wpw_occupied_w_basis)::candidate
    integer,allocatable::order(:),expected_ids(:),expected_keys(:,:)
    integer::i,j,k,nlocal,nglobal,nbuffer,tmp

    info=1;nlocal=size(local_keys,2);nglobal=size(global_keys,2)
    if(local_fragment<=0.or.epoch<=0.or.nglobal<=0.or.nlocal<=0.or.&
       size(global_keys,1)/=5.or.size(local_keys,1)/=5.or.size(global_owners)/=nglobal.or.&
       any(global_owners<=0).or.count(global_owners==local_fragment)/=nlocal.or.&
       size(core_grid_ids)<=0.or.any(shape(core_values)/=[size(core_grid_ids),nlocal]).or.&
       any(buffer_hi<buffer_lo).or..not.ieee_is_finite(source_condition).or.source_condition<=0d0.or.&
       .not.all(finite_complex(core_values)).or..not.all(finite_complex(buffer_values)).or.&
       .not.all(finite_complex(buffer_gradients)))return
    nbuffer=product(buffer_hi-buffer_lo+1)
    if(any(shape(buffer_values)/=[nbuffer,nlocal]).or.&
       any(shape(buffer_gradients)/=[3,nbuffer,nlocal]).or..not.strictly_increasing(core_grid_ids))return
    allocate(order(nglobal));order=[(i,i=1,nglobal)]
    do i=1,nglobal-1
      k=i
      do j=i+1,nglobal
        if(key_less(global_keys(:,order(j)),global_keys(:,order(k))))k=j
      enddo
      if(k/=i)then;tmp=order(i);order(i)=order(k);order(k)=tmp;endif
    enddo
    do i=2,nglobal
      if(all(global_keys(:,order(i))==global_keys(:,order(i-1))))return
    enddo
    allocate(expected_ids(nlocal),expected_keys(5,nlocal));j=0
    do i=1,nglobal
      if(global_owners(order(i))/=local_fragment)cycle
      j=j+1;expected_ids(j)=i;expected_keys(:,j)=global_keys(:,order(i))
    enddo
    if(j/=nlocal.or.any(expected_keys/=local_keys))return
    candidate%local_fragment=local_fragment;candidate%local_count=nlocal
    candidate%global_count=nglobal;candidate%epoch=epoch;candidate%source_condition=source_condition
    candidate%buffer_lo=buffer_lo;candidate%buffer_hi=buffer_hi
    candidate%owned_ids=expected_ids;candidate%stable_keys=local_keys
    candidate%core_grid_ids=core_grid_ids;candidate%core_values=core_values
    candidate%buffer_values=buffer_values;candidate%buffer_gradients=buffer_gradients
    candidate%fingerprint=basis_fingerprint(candidate)
    if(candidate%fingerprint==0_int64)candidate%fingerprint=1_int64
    candidate%valid=.true.;basis=candidate;info=0
  end subroutine initialize_dg_wpw_occupied_w_basis

  logical function key_less(left,right)result(less)
    integer,intent(in)::left(:),right(:)
    integer::i
    less=.false.
    do i=1,size(left)
      if(left(i)<right(i))then;less=.true.;return
      elseif(left(i)>right(i))then;return
      endif
    enddo
  end function key_less

  logical function strictly_increasing(values)result(ok)
    integer,intent(in)::values(:)
    integer::i
    ok=size(values)>0
    do i=2,size(values)
      if(values(i)<=values(i-1))then;ok=.false.;return;endif
    enddo
  end function strictly_increasing

  elemental logical function finite_complex(value)result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function finite_complex

  integer(int64) function basis_fingerprint(basis)result(hash)
    type(s_dg_wpw_occupied_w_basis),intent(in)::basis
    integer::i,j,k
    hash=1469598103934665603_int64
    call mix(hash,int(basis%local_fragment,int64));call mix(hash,int(basis%epoch,int64))
    call mix(hash,transfer(basis%source_condition,hash))
    do j=1,size(basis%stable_keys,2);do i=1,5
      call mix(hash,int(basis%stable_keys(i,j),int64))
    enddo;enddo
    do j=1,size(basis%core_values,2);do i=1,size(basis%core_values,1)
      call mix(hash,transfer(real(basis%core_values(i,j),8),hash))
      call mix(hash,transfer(aimag(basis%core_values(i,j)),hash))
    enddo;enddo
    do k=1,size(basis%buffer_gradients,3);do j=1,size(basis%buffer_gradients,2);do i=1,3
      call mix(hash,transfer(real(basis%buffer_gradients(i,j,k),8),hash))
      call mix(hash,transfer(aimag(basis%buffer_gradients(i,j,k)),hash))
    enddo;enddo;enddo
  end function basis_fingerprint

  subroutine mix(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    hash=ieor(hash,value);hash=ieor(hash,shiftl(hash,13));hash=ieor(hash,shiftr(hash,7))
  end subroutine mix
end module dg_wpw_occupied_w_basis
