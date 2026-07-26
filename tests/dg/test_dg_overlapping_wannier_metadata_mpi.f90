#include "config.h"
program test_dg_overlapping_wannier_metadata_mpi
#ifdef USE_MPI
  use mpi
#endif
  use dg_overlapping_wannier_types, only: s_dg_overlapping_wannier_basis, &
    initialize_dg_wannier_tail,validate_dg_overlapping_wannier_basis, &
    checked_dg_wannier_extent_product,release_dg_overlapping_wannier_basis
  implicit none
  integer :: comm,rank,nproc,ierr,scenario,ownership_count
  integer(8) :: fingerprint,product_value
  logical :: ok
  character(256) :: message
  type(s_dg_overlapping_wannier_basis) :: basis

#ifdef USE_MPI
  call MPI_Init(ierr)
  comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr)
  call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif

  call checked_dg_wannier_extent_product([2_8,3_8,4_8],product_value,ok)
  call require(ok.and.product_value==24_8,'checked extent product')
  call checked_dg_wannier_extent_product([huge(1_8),2_8],product_value,ok)
  call require(.not.ok,'overflow rejection')

  do scenario=2,4,2
    call make_valid_basis(scenario,basis)
    call validate_dg_overlapping_wannier_basis(basis,comm,int(scenario,8),ok,message,&
      ownership_count,fingerprint)
    call require(ok,trim(message))
    call require(ownership_count==scenario,'unique physical core ownership count')
    call require(fingerprint/=0_8,'nonzero deterministic fingerprint')
    call require_same_int(ownership_count,'ownership count agreement')
    call require_same_int8(fingerprint,'fingerprint agreement')
    call release_dg_overlapping_wannier_basis(basis)
  enddo

  call make_valid_basis(2,basis)
  call add_duplicate_center(basis)
  call validate_dg_overlapping_wannier_basis(basis,comm,2_8,ok,message)
  call require(.not.ok,'duplicate center owner rejection')
  call release_dg_overlapping_wannier_basis(basis)

  call make_valid_basis(2,basis)
  call add_duplicate_core(basis)
  call validate_dg_overlapping_wannier_basis(basis,comm,2_8,ok,message)
  call require(.not.ok,'duplicate core owner rejection')
  call release_dg_overlapping_wannier_basis(basis)

  call make_valid_basis(2,basis)
  if(rank==mod(1,nproc))basis%owned_core_physical_ids=[integer(8)::]
  call validate_dg_overlapping_wannier_basis(basis,comm,2_8,ok,message)
  call require(.not.ok,'missing core owner rejection')
  call release_dg_overlapping_wannier_basis(basis)

  call make_valid_basis(2,basis)
  if(rank==0.and.size(basis%tail)>0)basis%tail(1)%generation=basis%generation+1
  call validate_dg_overlapping_wannier_basis(basis,comm,2_8,ok,message)
  call require(.not.ok,'stale tail rejection')
  call release_dg_overlapping_wannier_basis(basis)

  call make_valid_basis(2,basis)
  if(rank==0.and.size(basis%tail)>0)then
    deallocate(basis%tail(1)%gradient)
    allocate(basis%tail(1)%gradient(3,1))
    basis%tail(1)%gradient=0d0
  endif
  call validate_dg_overlapping_wannier_basis(basis,comm,2_8,ok,message)
  call require(.not.ok,'incomplete tail map rejection')
  call release_dg_overlapping_wannier_basis(basis)

  if(rank==0)write(*,'(a,i0,a,i0)')'PASS overlapping-Wannier metadata on ',nproc, &
    ' ranks fingerprint=',fingerprint
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif

contains

  subroutine make_valid_basis(nfragment,b)
    integer,intent(in)::nfragment
    type(s_dg_overlapping_wannier_basis),intent(out)::b
    integer::fragment,ntail,index,component
    integer(8)::ids(2)
    real(8)::values(2),gradients(3,2)
    b%candidate_rank=nfragment+2
    b%target_rank=nfragment
    b%retained_rank=nfragment
    b%generation=17
    b%geometry_fingerprint=9001_8+int(nfragment,8)
    b%basis_fingerprint=7001_8+int(nfragment,8)
    allocate(b%owned_core_physical_ids(count([(mod(fragment-1,nproc)==rank,fragment=1,nfragment)])))
    index=0
    do fragment=1,nfragment
      if(mod(fragment-1,nproc)==rank)then
        index=index+1;b%owned_core_physical_ids(index)=int(fragment,8)
      endif
    enddo
    ntail=count([(mod(fragment-1,nproc)==rank,fragment=1,nfragment)])
    allocate(b%tail(ntail));index=0
    do fragment=1,nfragment
      if(mod(fragment-1,nproc)/=rank)cycle
      index=index+1
      ids=[int(fragment,8),int(mod(fragment,nfragment)+1,8)]
      values=[1d0+0.1d0*fragment,0.25d0*fragment]
      gradients=reshape([(0.01d0*dble(fragment+component),component=1,6)],[3,2])
      call initialize_dg_wannier_tail(b%tail(index),fragment,rank,1,b%generation,&
        b%geometry_fingerprint,ids,values,gradients,ok,message)
      if(.not.ok)error stop trim(message)
    enddo
  end subroutine make_valid_basis

  subroutine add_duplicate_center(b)
    type(s_dg_overlapping_wannier_basis),intent(inout)::b
    type(s_dg_overlapping_wannier_basis)::copy
    integer::n
    if(rank/=0)return
    n=size(b%tail)
    if(n<1)return
    allocate(copy%tail(n+1))
    copy%tail(1:n)=b%tail
    copy%tail(n+1)=b%tail(1)
    call move_alloc(copy%tail,b%tail)
  end subroutine add_duplicate_center

  subroutine add_duplicate_core(b)
    type(s_dg_overlapping_wannier_basis),intent(inout)::b
    integer(8),allocatable::copy(:)
    integer::target,n
    target=merge(1,0,nproc>1)
    if(rank/=target)return
    n=size(b%owned_core_physical_ids)
    allocate(copy(n+1))
    if(n>0)copy(1:n)=b%owned_core_physical_ids
    copy(n+1)=1_8
    call move_alloc(copy,b%owned_core_physical_ids)
  end subroutine add_duplicate_core

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
#ifdef USE_MPI
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    global_failure=local_failure
#endif
    if(global_failure/=0)error stop trim(label)
  end subroutine require

  subroutine require_same_int(value,label)
    integer,intent(in)::value
    character(*),intent(in)::label
    integer::minimum,maximum
#ifdef USE_MPI
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    minimum=value;maximum=value
#endif
    call require(minimum==maximum,label)
  end subroutine require_same_int

  subroutine require_same_int8(value,label)
    integer(8),intent(in)::value
    character(*),intent(in)::label
    integer(8)::minimum,maximum
#ifdef USE_MPI
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
#else
    minimum=value;maximum=value
#endif
    call require(minimum==maximum,label)
  end subroutine require_same_int8
end program test_dg_overlapping_wannier_metadata_mpi
