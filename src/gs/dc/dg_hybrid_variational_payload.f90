#include "config.h"
module dg_hybrid_variational_payload
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_fixed_payload
    logical::frozen=.false.
    integer::global_basis_count=0
    integer(int64),allocatable::row_ids(:)
    integer(int64),allocatable::row_fingerprints(:)
    complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),interface_rows(:,:)
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64,interface_fingerprint=0_int64
    integer(int64)::basis_directory_fingerprint=0_int64
    integer(int64)::metadata_fingerprint=0_int64
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
  public::verify_dg_hybrid_variational_payload_rows
  public::write_dg_hybrid_variational_payload_bundle,read_dg_hybrid_variational_payload_bundle
  integer,parameter::variational_payload_bundle_version=1
  interface
    integer(c_int) function c_rename(old_path,new_path) bind(C,name='rename')
      import::c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
    end function c_rename
  end interface
contains
  subroutine write_dg_hybrid_variational_payload_bundle(comm,prefix,global_basis_count,row_ids,&
      metric_rows,kinetic_rows,nonlocal_rows,interface_rows,basis_fingerprint,metric_fingerprint,&
      interface_fingerprint,ok,message)
    integer,intent(in)::comm,global_basis_count
    character(*),intent(in)::prefix
    integer(int64),intent(in)::row_ids(:),basis_fingerprint,metric_fingerprint,interface_fingerprint
    complex(real64),intent(in)::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),interface_rows(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,ios,local_bad,global_bad,dims(2,4)
    logical::exists
    character(1024)::shard,temporary,manifest,manifest_temporary
    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    write(shard,'(a,".rank",i6.6)')trim(prefix),rank
    temporary=trim(shard)//'.tmp';manifest=trim(prefix)//'.manifest'
    manifest_temporary=trim(manifest)//'.tmp'
    inquire(file=trim(shard),exist=exists);if(exists)local_bad=1
    if(rank==0)then;inquire(file=trim(manifest),exist=exists);if(exists)local_bad=1;endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='variational payload bundle already exists';return;endif
    dims(:,1)=shape(metric_rows);dims(:,2)=shape(kinetic_rows)
    dims(:,3)=shape(nonlocal_rows);dims(:,4)=shape(interface_rows)
    open(newunit=unit,file=trim(temporary),status='replace',access='stream',form='unformatted',&
      action='write',iostat=ios)
    if(ios==0)write(unit,iostat=ios)variational_payload_bundle_version,nproc,rank,global_basis_count,&
      size(row_ids),dims,basis_fingerprint,metric_fingerprint,interface_fingerprint
    if(ios==0)write(unit,iostat=ios)row_ids,metric_rows,kinetic_rows,nonlocal_rows,interface_rows
    if(ios==0)close(unit,iostat=ios)
    if(ios==0)call atomic_rename(temporary,shard,ios)
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='variational payload rank shard publication failed';return;endif
    call MPI_Barrier(comm,ierr)
    ios=0
    if(rank==0)then
      open(newunit=unit,file=trim(manifest_temporary),status='replace',access='stream',form='unformatted',&
        action='write',iostat=ios)
      if(ios==0)write(unit,iostat=ios)variational_payload_bundle_version,nproc
      if(ios==0)close(unit,iostat=ios)
      if(ios==0)call atomic_rename(manifest_temporary,manifest,ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.ios/=0)then;message='variational payload manifest publication failed';return;endif
    ok=.true.
#else
    ok=.false.;message='MPI is required for variational payload bundle writing'
#endif
  end subroutine write_dg_hybrid_variational_payload_bundle

  subroutine read_dg_hybrid_variational_payload_bundle(comm,prefix,global_basis_count,row_ids,&
      metric_rows,kinetic_rows,nonlocal_rows,interface_rows,basis_fingerprint,metric_fingerprint,&
      interface_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    integer,intent(out)::global_basis_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    complex(real64),allocatable,intent(out)::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      interface_rows(:,:)
    integer(int64),intent(out)::basis_fingerprint,metric_fingerprint,interface_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,ios,local_bad,global_bad,version,stored_nproc,stored_rank,nrows,dims(2,4)
    integer::minimum_integer,maximum_integer
    integer(int64)::minimum_fingerprint,maximum_fingerprint
    character(1024)::shard,manifest
    ok=.false.;message='';global_basis_count=0;basis_fingerprint=0_int64
    metric_fingerprint=0_int64;interface_fingerprint=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    manifest=trim(prefix)//'.manifest'
    open(newunit=unit,file=trim(manifest),status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios==0)read(unit,iostat=ios)version,stored_nproc
    if(ios==0)close(unit,iostat=ios)
    if(ios/=0.or.version/=variational_payload_bundle_version.or.stored_nproc/=nproc)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid variational payload bundle manifest';return;endif
    write(shard,'(a,".rank",i6.6)')trim(prefix),rank
    open(newunit=unit,file=trim(shard),status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios==0)read(unit,iostat=ios)version,stored_nproc,stored_rank,global_basis_count,nrows,dims,&
      basis_fingerprint,metric_fingerprint,interface_fingerprint
    if(ios/=0.or.version/=variational_payload_bundle_version.or.stored_nproc/=nproc.or.&
      stored_rank/=rank.or.global_basis_count<0.or.nrows<0.or.any(dims<0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(ios==0)close(unit);message='invalid variational payload rank shard header';return
    endif
    allocate(row_ids(nrows),metric_rows(dims(1,1),dims(2,1)),kinetic_rows(dims(1,2),dims(2,2)),&
      nonlocal_rows(dims(1,3),dims(2,3)),interface_rows(dims(1,4),dims(2,4)),stat=ios)
    if(ios==0)read(unit,iostat=ios)row_ids,metric_rows,kinetic_rows,nonlocal_rows,interface_rows
    if(ios==0)close(unit,iostat=ios)
    local_bad=merge(0,1,ios==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='truncated variational payload rank shard';return;endif
    call MPI_Allreduce(global_basis_count,minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(global_basis_count,maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(basis_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(basis_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer.or.minimum_fingerprint/=maximum_fingerprint)then
      message='variational payload shard metadata disagree across ranks';return
    endif
    ok=.true.
#else
    ok=.false.;message='MPI is required for variational payload bundle reading'
#endif
  end subroutine read_dg_hybrid_variational_payload_bundle

  subroutine atomic_rename(old_path,new_path,status)
    character(*),intent(in)::old_path,new_path
    integer,intent(out)::status
    character(c_char),allocatable::old_c(:),new_c(:)
    integer::i
    allocate(old_c(len_trim(old_path)+1),new_c(len_trim(new_path)+1))
    do i=1,len_trim(old_path);old_c(i)=old_path(i:i);enddo;old_c(size(old_c))=c_null_char
    do i=1,len_trim(new_path);new_c(i)=new_path(i:i);enddo;new_c(size(new_c))=c_null_char
    status=int(c_rename(old_c,new_c))
  end subroutine atomic_rename

  subroutine freeze_dg_hybrid_variational_payload(comm,global_basis_count,row_ids,metric_rows,kinetic_rows,&
      nonlocal_rows,interface_rows,basis_fingerprint,metric_fingerprint,interface_fingerprint,payload,ok,message,&
      basis_directory_fingerprint)
    integer,intent(in)::comm,global_basis_count
    integer(int64),intent(in)::row_ids(:),basis_fingerprint,metric_fingerprint,interface_fingerprint
    complex(real64),intent(in)::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),interface_rows(:,:)
    type(s_dg_hybrid_fixed_payload),intent(out)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(in),optional::basis_directory_fingerprint
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
    if(present(basis_directory_fingerprint))then
      if(basis_directory_fingerprint==0_int64)local_bad=ibset(local_bad,6)
    endif
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
    if(btest(global_bad,6))then;message='invalid variational basis-directory fingerprint';return;endif
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
    payload%basis_directory_fingerprint=0_int64
    if(present(basis_directory_fingerprint))payload%basis_directory_fingerprint=basis_directory_fingerprint
    payload%metadata_fingerprint=compute_payload_metadata_fingerprint(payload)
    allocate(payload%row_fingerprints(size(payload%row_ids)))
    do i=1,size(payload%row_ids)
      payload%row_fingerprints(i)=compute_payload_row_fingerprint(payload,i)
    enddo
    call compute_payload_fingerprint(comm,payload,payload%fingerprint,ok)
    if(.not.ok.or.payload%fingerprint==0_int64)then;message='variational fixed fingerprint failed';return;endif
    payload%frozen=.true.;message=''
#else
    ok=.false.;message='MPI is required for variational payload freezing'
#endif
  end subroutine freeze_dg_hybrid_variational_payload

  subroutine verify_dg_hybrid_variational_payload_rows(comm,payload,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,local_bad,global_bad,nrows
    logical::storage_ready
    integer(int64)::current_metadata_fingerprint
    ok=.false.;message='';local_bad=0;nrows=0
    if(.not.payload%frozen)local_bad=ibset(local_bad,0)
    if(payload%global_basis_count<1.or.payload%basis_fingerprint==0_int64.or.&
        payload%metric_fingerprint==0_int64.or.payload%interface_fingerprint==0_int64.or.&
        payload%metadata_fingerprint==0_int64.or.payload%fingerprint==0_int64)local_bad=ibset(local_bad,1)
    storage_ready=allocated(payload%row_ids).and.allocated(payload%row_fingerprints).and.&
      allocated(payload%metric_rows).and.allocated(payload%kinetic_rows).and.&
      allocated(payload%nonlocal_rows).and.allocated(payload%interface_rows)
    if(.not.storage_ready)then
      local_bad=ibset(local_bad,2)
    else
      nrows=size(payload%row_ids)
      if(size(payload%row_fingerprints)/=nrows.or.&
          any(shape(payload%metric_rows)/=[nrows,payload%global_basis_count]).or.&
          any(shape(payload%kinetic_rows)/=shape(payload%metric_rows)).or.&
          any(shape(payload%nonlocal_rows)/=shape(payload%metric_rows)).or.&
          any(shape(payload%interface_rows)/=shape(payload%metric_rows)))local_bad=ibset(local_bad,2)
    endif
    if(storage_ready.and..not.btest(local_bad,2))then
      if(any(payload%row_ids<1_int64).or.any(payload%row_ids>int(payload%global_basis_count,int64)))&
        local_bad=ibset(local_bad,1)
      if(.not.finite_matrix(payload%metric_rows).or..not.finite_matrix(payload%kinetic_rows).or.&
          .not.finite_matrix(payload%nonlocal_rows).or..not.finite_matrix(payload%interface_rows))&
        local_bad=ibset(local_bad,3)
      current_metadata_fingerprint=compute_payload_metadata_fingerprint(payload)
      if(current_metadata_fingerprint/=payload%metadata_fingerprint)local_bad=ibset(local_bad,4)
      do i=1,nrows
        if(compute_payload_row_fingerprint(payload,i)/=payload%row_fingerprints(i))&
          local_bad=ibset(local_bad,5)
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_BOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fixed variational payload row verification reduction failed';return;endif
    if(btest(global_bad,0))then;message='fixed variational payload is not frozen';return;endif
    if(btest(global_bad,1))then;message='invalid fixed variational payload metadata';return;endif
    if(btest(global_bad,2))then;message='invalid fixed variational payload row storage';return;endif
    if(btest(global_bad,3))then;message='fixed variational payload rows contain nonfinite values';return;endif
    if(btest(global_bad,4))then;message='fixed variational payload metadata fingerprint changed';return;endif
    if(btest(global_bad,5))then;message='fixed variational payload row fingerprint changed';return;endif
    ok=.true.
#else
    ok=.false.;message='MPI is required for variational payload row verification'
#endif
  end subroutine verify_dg_hybrid_variational_payload_rows

  subroutine compose_dg_hybrid_variational_hamiltonian(comm,payload,local_rows,lambda,epoch,iterate,ok,message)
    integer,intent(in)::comm,epoch
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    complex(real64),intent(in)::local_rows(:,:)
    real(real64),intent(in)::lambda
    type(s_dg_hybrid_variational_iterate),intent(out)::iterate
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_bad,global_bad,allocation_status
    integer(int64)::current_fingerprint
    complex(real64),allocatable::working_local_rows(:,:),working_hamiltonian_rows(:,:)
    ok=.false.;message='';local_bad=0
    call verify_dg_hybrid_variational_payload_rows(comm,payload,ok,message)
    if(.not.ok)return
    ok=.false.
    if(epoch<1.or..not.ieee_is_finite(lambda).or.lambda<0d0.or.lambda>1d0.or.&
        any(shape(local_rows)/=shape(payload%metric_rows)).or..not.finite_matrix(local_rows))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid variational Hamiltonian composition';return;endif
    call compute_payload_fingerprint(comm,payload,current_fingerprint,ok)
    if(.not.ok.or.current_fingerprint/=payload%fingerprint)then
      ok=.false.;message='fixed variational payload fingerprint changed';return
    endif
    ok=.false.
    allocation_status=0
    allocate(working_local_rows,source=local_rows,stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='variational local-row allocation failed collectively';return
    endif
    allocation_status=0
    allocate(working_hamiltonian_rows(size(local_rows,1),size(local_rows,2)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='variational Hamiltonian-row allocation failed collectively';return
    endif
    working_hamiltonian_rows=payload%kinetic_rows+payload%nonlocal_rows+local_rows+lambda*payload%interface_rows
    call move_alloc(working_local_rows,iterate%local_rows)
    call move_alloc(working_hamiltonian_rows,iterate%hamiltonian_rows)
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
    integer::i,ierr,shift
    integer(int64)::local_hash,row_hash
    fingerprint=0_int64;local_hash=0_int64
    do i=1,size(payload%row_ids)
      row_hash=compute_payload_row_fingerprint(payload,i)
      shift=int(modulo(payload%row_ids(i),63_int64))
      local_hash=ieor(local_hash,ishftc(row_hash,shift))
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr==MPI_SUCCESS)call mix_hash_word(fingerprint,compute_payload_metadata_fingerprint(payload),97_int64)
    if(fingerprint==0_int64)fingerprint=1543_int64
    ok=ierr==MPI_SUCCESS
  end subroutine compute_payload_fingerprint

  pure integer(int64) function compute_payload_metadata_fingerprint(payload) result(fingerprint)
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    fingerprint=int(z'243F6A8885A308D3',int64)
    call mix_hash_word(fingerprint,int(payload%global_basis_count,int64),1_int64)
    call mix_hash_word(fingerprint,payload%basis_fingerprint,2_int64)
    call mix_hash_word(fingerprint,payload%metric_fingerprint,3_int64)
    call mix_hash_word(fingerprint,payload%interface_fingerprint,4_int64)
    call mix_hash_word(fingerprint,payload%basis_directory_fingerprint,5_int64)
    if(fingerprint==0_int64)fingerprint=1543_int64
  end function compute_payload_metadata_fingerprint

  pure integer(int64) function compute_payload_row_fingerprint(payload,index) result(fingerprint)
    type(s_dg_hybrid_fixed_payload),intent(in)::payload
    integer,intent(in)::index
    integer::column
    fingerprint=compute_payload_metadata_fingerprint(payload)
    call mix_hash_word(fingerprint,payload%row_ids(index),11_int64)
    do column=1,payload%global_basis_count
      call hash_complex(fingerprint,column,1,payload%metric_rows(index,column))
      call hash_complex(fingerprint,column,2,payload%kinetic_rows(index,column))
      call hash_complex(fingerprint,column,3,payload%nonlocal_rows(index,column))
      call hash_complex(fingerprint,column,4,payload%interface_rows(index,column))
    enddo
    if(fingerprint==0_int64)fingerprint=1543_int64
  end function compute_payload_row_fingerprint

  pure subroutine hash_complex(hash,column,component,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::column,component
    complex(real64),intent(in)::value
    integer(int64)::bits,position
    position=int(column,int64)
    call mix_hash_word(hash,position,100_int64+int(component,int64))
    bits=transfer(real(value,real64),bits)
    call mix_hash_word(hash,bits,200_int64+int(component,int64))
    bits=transfer(aimag(value),bits)
    call mix_hash_word(hash,bits,300_int64+int(component,int64))
  end subroutine hash_complex

  pure subroutine mix_hash_word(hash,word,tag)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::word,tag
    integer::shift
    shift=1+int(modulo(ieor(word,ishftc(tag,7)),63_int64))
    hash=ieor(ishftc(hash,shift),word)
    hash=ieor(hash,ishftc(tag,modulo(shift+23,64)))
    hash=ieor(hash,int(z'13198A2E03707344',int64))
  end subroutine mix_hash_word
#endif

  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
end module dg_hybrid_variational_payload
