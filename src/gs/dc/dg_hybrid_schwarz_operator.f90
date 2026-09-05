module dg_hybrid_schwarz_operator
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  type,public::s_dg_hybrid_schwarz_schedule
    logical::valid=.false.
    integer::fragment_id=0,fragment_count=0,basis_generation=0,peer_count=0
    integer(int64)::directory_fingerprint=0_int64,face_fingerprint=0_int64
    integer(int64)::mapping_fingerprint=0_int64,fingerprint=0_int64
    integer,allocatable::peers(:),receive_offsets(:),send_offsets(:),send_local_slots(:)
    integer(int64),allocatable::receive_global_ids(:)
  end type

  public::build_dg_hybrid_schwarz_schedule
contains
  subroutine build_dg_hybrid_schwarz_schedule(comm,fragment_id,basis_generation,row_ids,basis_owner,&
      basis_fragment,basis_local_slot,basis_generations,interface_rows,directory_fingerprint,&
      face_fingerprint,mapping_fingerprint,schedule,ok,message,declared_neighbors)
    integer,intent(in)::comm,fragment_id,basis_generation
    integer(int64),intent(in)::row_ids(:)
    integer,intent(in)::basis_owner(:),basis_fragment(:),basis_local_slot(:),basis_generations(:)
    complex(real64),intent(in)::interface_rows(:,:)
    integer(int64),intent(in)::directory_fingerprint,face_fingerprint,mapping_fingerprint
    type(s_dg_hybrid_schwarz_schedule),intent(inout)::schedule
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::declared_neighbors(:)
    type(s_dg_hybrid_schwarz_schedule)::work
    integer::rank,nproc,ierr,global_count,p,q,g,peer,index,stat,receive_count,send_count
    integer::controls(2),minimum_controls(2),maximum_controls(2)
    integer(int64)::fingerprints(3),minimum_fingerprints(3),maximum_fingerprints(3)
    integer,allocatable::minimum_directory(:,:),maximum_directory(:,:),presence(:),requests(:,:),adjacency(:,:)
    integer,allocatable::peers(:),receive_offsets(:),send_offsets(:),send_slots(:)
    integer(int64),allocatable::receive_ids(:)
    logical::valid

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule size query failed';return;endif
    global_count=size(basis_owner)
    controls=[basis_generation,global_count]
    call MPI_Allreduce(controls,minimum_controls,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule control minimum failed';return;endif
    call MPI_Allreduce(controls,maximum_controls,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule control maximum failed';return;endif
    fingerprints=[directory_fingerprint,face_fingerprint,mapping_fingerprint]
    call MPI_Allreduce(fingerprints,minimum_fingerprints,3,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule fingerprint minimum failed';return;endif
    call MPI_Allreduce(fingerprints,maximum_fingerprints,3,MPI_INTEGER8,MPI_MAX,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(controls==minimum_controls).and.all(controls==maximum_controls).and.&
      all(fingerprints==minimum_fingerprints).and.all(fingerprints==maximum_fingerprints).and.&
      fragment_id==rank+1.and.basis_generation>0.and.global_count>0.and.all(fingerprints/=0_int64)
    call collective_gate(comm,valid,'invalid or rank-disagreeing Schwarz schedule controls',ok,message)
    if(.not.ok)return
    valid=size(basis_fragment)==global_count.and.size(basis_local_slot)==global_count.and.&
      size(basis_generations)==global_count.and.size(interface_rows,1)==size(row_ids).and.&
      size(interface_rows,2)==global_count.and.all(row_ids>=1_int64).and.&
      all(row_ids<=int(global_count,int64)).and.unique_int64(row_ids).and.finite_matrix(interface_rows)
    call collective_gate(comm,valid,'invalid Schwarz schedule directory or interface shape',ok,message)
    if(.not.ok)return
    allocate(minimum_directory(global_count,4),maximum_directory(global_count,4),presence(global_count),&
      requests(global_count,nproc),adjacency(nproc,nproc),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz schedule workspace allocation failed',ok,message)
    if(.not.ok)return
    minimum_directory(:,1)=basis_owner;minimum_directory(:,2)=basis_fragment
    minimum_directory(:,3)=basis_local_slot;minimum_directory(:,4)=basis_generations
    maximum_directory=minimum_directory
    call MPI_Allreduce(MPI_IN_PLACE,minimum_directory,4*global_count,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz directory minimum failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,maximum_directory,4*global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(minimum_directory==maximum_directory)
    if(valid)valid=all(basis_fragment>=1).and.all(basis_fragment<=nproc).and.&
      all(basis_owner==basis_fragment-1).and.all(basis_local_slot>=1).and.&
      all(basis_generations==basis_generation)
    do g=1,nproc
      if(valid)valid=all(pack(basis_local_slot,basis_fragment==g)==&
        [(p,p=1,count(basis_fragment==g))])
    enddo
    call collective_gate(comm,valid,'stale or inconsistent Schwarz basis directory',ok,message)
    if(.not.ok)return
    presence=0
    do p=1,size(row_ids);presence(int(row_ids(p)))=presence(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,presence,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(presence==1)
    if(valid)valid=all([(basis_fragment(int(row_ids(p)))==fragment_id,p=1,size(row_ids))])
    call collective_gate(comm,valid,'Schwarz operator rows are not owned exactly once',ok,message)
    if(.not.ok)return
    requests=0
    do p=1,size(row_ids)
      do q=1,global_count
        if(basis_fragment(q)/=fragment_id.and.abs(interface_rows(p,q))>0d0)requests(q,fragment_id)=1
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,requests,global_count*nproc,MPI_INTEGER,MPI_SUM,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(requests>=0).and.all(requests<=1)
    adjacency=0
    if(valid)then
      do g=1,nproc
        do q=1,global_count
          if(requests(q,g)==1)adjacency(g,basis_fragment(q))=1
        enddo
      enddo
      valid=all([(adjacency(g,g)==0,g=1,nproc)]).and.all(adjacency==transpose(adjacency))
    endif
    call collective_gate(comm,valid,'missing, duplicate, or nonreciprocal Schwarz face request',ok,message)
    if(.not.ok)return
    allocate(peers(count(adjacency(fragment_id,:)==1)),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz peer allocation failed',ok,message)
    if(.not.ok)return
    index=0
    do g=1,nproc
      if(adjacency(fragment_id,g)==1)then;index=index+1;peers(index)=g;endif
    enddo
    valid=.true.
    if(present(declared_neighbors))then
      valid=all(declared_neighbors>=1).and.all(declared_neighbors<=nproc).and.&
        all(declared_neighbors/=fragment_id).and.unique_integer(declared_neighbors)
      if(valid)valid=size(declared_neighbors)==size(peers)
      if(valid.and.size(peers)>0)valid=all(declared_neighbors==peers)
    endif
    call collective_gate(comm,valid,'declared Schwarz neighbors are missing, duplicate, or non-neighbor',ok,message)
    if(.not.ok)return
    receive_count=0;send_count=0
    do p=1,size(peers)
      peer=peers(p)
      receive_count=receive_count+count(requests(:,fragment_id)==1.and.basis_fragment==peer)
      send_count=send_count+count(requests(:,peer)==1.and.basis_fragment==fragment_id)
    enddo
    allocate(receive_offsets(size(peers)+1),send_offsets(size(peers)+1),receive_ids(receive_count),&
      send_slots(send_count),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz transfer schedule allocation failed',ok,message)
    if(.not.ok)return
    receive_offsets(1)=1;send_offsets(1)=1;receive_count=0;send_count=0
    do p=1,size(peers)
      peer=peers(p)
      do q=1,global_count
        if(requests(q,fragment_id)==1.and.basis_fragment(q)==peer)then
          receive_count=receive_count+1;receive_ids(receive_count)=int(q,int64)
        endif
        if(requests(q,peer)==1.and.basis_fragment(q)==fragment_id)then
          send_count=send_count+1;send_slots(send_count)=basis_local_slot(q)
        endif
      enddo
      receive_offsets(p+1)=receive_count+1;send_offsets(p+1)=send_count+1
    enddo
    work%valid=.true.;work%fragment_id=fragment_id;work%fragment_count=nproc
    work%basis_generation=basis_generation;work%peer_count=size(peers)
    work%directory_fingerprint=directory_fingerprint;work%face_fingerprint=face_fingerprint
    work%mapping_fingerprint=mapping_fingerprint
    call move_alloc(peers,work%peers);call move_alloc(receive_offsets,work%receive_offsets)
    call move_alloc(send_offsets,work%send_offsets);call move_alloc(receive_ids,work%receive_global_ids)
    call move_alloc(send_slots,work%send_local_slots)
    work%fingerprint=schedule_fingerprint(basis_generation,directory_fingerprint,face_fingerprint,&
      mapping_fingerprint,basis_owner,basis_fragment,basis_local_slot,requests)
    call collective_gate(comm,work%fingerprint/=0_int64,'Schwarz schedule fingerprint is zero',ok,message)
    if(.not.ok)return
    schedule=work;ok=.true.;message=''
  end subroutine build_dg_hybrid_schwarz_schedule

  integer(int64) function schedule_fingerprint(generation,directory_fp,face_fp,mapping_fp,owners,&
      fragments,slots,requests)result(hash)
    integer,intent(in)::generation,owners(:),fragments(:),slots(:),requests(:,:)
    integer(int64),intent(in)::directory_fp,face_fp,mapping_fp
    integer::p,g
    hash=int(z'A54FF53A5F1D36F1',int64)
    hash=mix_hash(hash,int(generation,int64));hash=mix_hash(hash,directory_fp)
    hash=mix_hash(hash,face_fp);hash=mix_hash(hash,mapping_fp)
    do p=1,size(owners)
      hash=mix_hash(hash,int(owners(p)+1,int64));hash=mix_hash(hash,int(fragments(p),int64))
      hash=mix_hash(hash,int(slots(p),int64))
      do g=1,size(requests,2);hash=mix_hash(hash,int(requests(p,g),int64));enddo
    enddo
    if(hash==0_int64)hash=1_int64
  end function schedule_fingerprint

  pure integer(int64) function mix_hash(hash,value)result(mixed)
    integer(int64),intent(in)::hash,value
    mixed=ieor(ishftc(hash,19),value);mixed=ieor(mixed,ishftc(value,37))
  end function mix_hash

  pure logical function unique_integer(values)result(unique)
    integer,intent(in)::values(:)
    integer::p
    unique=.true.
    do p=1,size(values)
      if(count(values==values(p))/=1)then;unique=.false.;return;endif
    enddo
  end function unique_integer

  pure logical function unique_int64(values)result(unique)
    integer(int64),intent(in)::values(:)
    integer::p
    unique=.true.
    do p=1,size(values)
      if(count(values==values(p))/=1)then;unique=.false.;return;endif
    enddo
  end function unique_int64

  pure logical function finite_matrix(values)result(finite)
    complex(real64),intent(in)::values(:,:)
    finite=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix

  subroutine collective_gate(comm,local_ok,detail,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    character(*),intent(in)::detail
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,failed,first_failed,ierr
    character(512)::shared
    ok=.false.;message='Schwarz schedule status rank query failed'
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    failed=huge(0);if(.not.local_ok)failed=rank
    call MPI_Allreduce(failed,first_failed,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule status reduction failed';return;endif
    if(first_failed==huge(0))then;ok=.true.;message='';return;endif
    shared='';if(rank==first_failed)shared=detail
    call MPI_Bcast(shared,len(shared),MPI_CHARACTER,first_failed,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz schedule diagnostic broadcast failed';return;endif
    message=trim(shared)
  end subroutine collective_gate
end module dg_hybrid_schwarz_operator
