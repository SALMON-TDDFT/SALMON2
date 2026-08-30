#include "config.h"
module dg_hybrid_variational_payload
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_fixed_payload
    logical::frozen=.false.
    integer::global_basis_count=0
    integer(int64),allocatable::row_ids(:)
    complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),interface_rows(:,:)
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64,interface_fingerprint=0_int64
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_fixed_payload
  type,public::s_dg_hybrid_variational_iterate
    integer::epoch=0
    real(real64)::lambda=0d0
    complex(real64),allocatable::local_rows(:,:),hamiltonian_rows(:,:)
    integer(int64)::fixed_fingerprint=0_int64
  end type s_dg_hybrid_variational_iterate
  type,public::s_dg_hybrid_accepted_variational_state
    logical::valid=.false.
    type(s_dg_hybrid_variational_iterate)::iterate
  end type s_dg_hybrid_accepted_variational_state
  public::freeze_dg_hybrid_variational_payload,compose_dg_hybrid_variational_hamiltonian
contains
  subroutine freeze_dg_hybrid_variational_payload(comm,global_basis_count,row_ids,metric_rows,kinetic_rows,&
      nonlocal_rows,interface_rows,basis_fingerprint,metric_fingerprint,interface_fingerprint,payload,ok,message)
    integer,intent(in)::comm,global_basis_count
    integer(int64),intent(in)::row_ids(:),basis_fingerprint,metric_fingerprint,interface_fingerprint
    complex(real64),intent(in)::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),interface_rows(:,:)
    type(s_dg_hybrid_fixed_payload),intent(out)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,local_bad,global_bad,total_rows
    integer,allocatable::ownership(:)
    ok=.false.;message='';local_bad=0
    if(global_basis_count<1.or.any(row_ids<1_int64).or.any(row_ids>int(global_basis_count,int64)).or.&
        any(shape(metric_rows)/=[size(row_ids),global_basis_count]).or.&
        any(shape(kinetic_rows)/=shape(metric_rows)).or.any(shape(nonlocal_rows)/=shape(metric_rows)).or.&
        any(shape(interface_rows)/=shape(metric_rows)))local_bad=ibset(local_bad,0)
    if(basis_fingerprint==0_int64.or.metric_fingerprint==0_int64.or.&
      interface_fingerprint==0_int64)local_bad=ibset(local_bad,1)
    if(.not.finite_matrix(metric_rows))local_bad=ibset(local_bad,2)
    if(.not.finite_matrix(kinetic_rows))local_bad=ibset(local_bad,3)
    if(.not.finite_matrix(nonlocal_rows))local_bad=ibset(local_bad,4)
    if(.not.finite_matrix(interface_rows))local_bad=ibset(local_bad,5)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_BOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='variational fixed payload validation reduction failed';return;endif
    if(btest(global_bad,0))then;message='invalid variational fixed payload extent';return;endif
    if(btest(global_bad,1))then;message='invalid variational fixed payload fingerprint';return;endif
    if(btest(global_bad,2))then;message='variational metric rows contain nonfinite values';return;endif
    if(btest(global_bad,3))then;message='variational kinetic rows contain nonfinite values';return;endif
    if(btest(global_bad,4))then;message='variational nonlocal rows contain nonfinite values';return;endif
    if(btest(global_bad,5))then;message='variational interface rows contain nonfinite values';return;endif
    call MPI_Allreduce(size(row_ids),total_rows,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    allocate(ownership(global_basis_count));ownership=0
    do i=1,size(row_ids);ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.total_rows/=global_basis_count.or.any(ownership/=1))then
      message='variational fixed rows are not owned exactly once';return
    endif
    payload%global_basis_count=global_basis_count
    allocate(payload%row_ids,source=row_ids);allocate(payload%metric_rows,source=metric_rows)
    allocate(payload%kinetic_rows,source=kinetic_rows);allocate(payload%nonlocal_rows,source=nonlocal_rows)
    allocate(payload%interface_rows,source=interface_rows)
    payload%basis_fingerprint=basis_fingerprint;payload%metric_fingerprint=metric_fingerprint
    payload%interface_fingerprint=interface_fingerprint
    call compute_payload_fingerprint(comm,payload,payload%fingerprint,ok)
    if(.not.ok.or.payload%fingerprint==0_int64)then;message='variational fixed fingerprint failed';return;endif
    payload%frozen=.true.;message=''
#else
    ok=.false.;message='MPI is required for variational payload freezing'
#endif
  end subroutine freeze_dg_hybrid_variational_payload

  subroutine compose_dg_hybrid_variational_hamiltonian(comm,payload,local_rows,lambda,epoch,iterate,ok,message)
    integer,intent(in)::comm,epoch
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    complex(real64),intent(in)::local_rows(:,:)
    real(real64),intent(in)::lambda
    type(s_dg_hybrid_variational_iterate),intent(out)::iterate
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_bad,global_bad
    integer(int64)::current_fingerprint
    ok=.false.;message='';local_bad=0
    if(.not.payload%frozen.or.epoch<1.or..not.ieee_is_finite(lambda).or.lambda<0d0.or.lambda>1d0.or.&
        any(shape(local_rows)/=shape(payload%metric_rows)).or..not.finite_matrix(local_rows))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid variational Hamiltonian composition';return;endif
    call compute_payload_fingerprint(comm,payload,current_fingerprint,ok)
    if(.not.ok.or.current_fingerprint/=payload%fingerprint)then
      ok=.false.;message='fixed variational payload fingerprint changed';return
    endif
    allocate(iterate%local_rows,source=local_rows)
    allocate(iterate%hamiltonian_rows(size(local_rows,1),size(local_rows,2)))
    iterate%hamiltonian_rows=payload%kinetic_rows+payload%nonlocal_rows+local_rows+lambda*payload%interface_rows
    iterate%lambda=lambda;iterate%epoch=epoch;iterate%fixed_fingerprint=payload%fingerprint
    ok=.true.;message=''
#else
    ok=.false.;message='MPI is required for variational Hamiltonian composition'
#endif
  end subroutine compose_dg_hybrid_variational_hamiltonian

#ifdef USE_MPI
  subroutine compute_payload_fingerprint(comm,payload,fingerprint,ok)
    integer,intent(in)::comm
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr
    integer(int64)::local_hash,bits
    local_hash=ieor(payload%basis_fingerprint,ishftc(payload%metric_fingerprint,11))
    local_hash=ieor(local_hash,ishftc(payload%interface_fingerprint,23))
    do i=1,size(payload%row_ids);do j=1,payload%global_basis_count
      call hash_complex(local_hash,payload%row_ids(i),j,payload%metric_rows(i,j))
      call hash_complex(local_hash,payload%row_ids(i),j+payload%global_basis_count,payload%kinetic_rows(i,j))
      call hash_complex(local_hash,payload%row_ids(i),j+2*payload%global_basis_count,payload%nonlocal_rows(i,j))
      call hash_complex(local_hash,payload%row_ids(i),j+3*payload%global_basis_count,payload%interface_rows(i,j))
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    bits=int(payload%global_basis_count,int64);fingerprint=ieor(fingerprint,ishftc(bits,37))
    if(fingerprint==0_int64)fingerprint=1543_int64
    ok=ierr==MPI_SUCCESS
  end subroutine compute_payload_fingerprint

  subroutine hash_complex(hash,row,column,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::row
    integer,intent(in)::column
    complex(real64),intent(in)::value
    integer(int64)::bits
    bits=transfer(real(value,real64),bits);hash=ieor(hash,ishftc(ieor(bits,row),mod(7*column,63)))
    bits=transfer(aimag(value),bits);hash=ieor(hash,ishftc(ieor(bits,ishftc(row,9)),mod(13*column,63)))
  end subroutine hash_complex
#endif

  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
end module dg_hybrid_variational_payload
