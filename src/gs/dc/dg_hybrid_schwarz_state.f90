module dg_hybrid_schwarz_state
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  real(real64),parameter::boltzmann_hartree_per_kelvin=3.166811563d-6

  type,public::s_dg_hybrid_schwarz_state
    logical::valid=.false.
    integer::fragment_id=0,fragment_count=0,basis_generation=0
    integer::local_basis_count=0,trial_count=0,candidate_count=0,thermal_tail_count=0
    integer::coefficient_epoch=0
    integer(int64)::mapping_fingerprint=0_int64,candidate_fingerprint=0_int64,fingerprint=0_int64
    real(real64)::electron_target=0d0,temperature=0d0,wspin=0d0,chemical_potential=0d0
    real(real64)::electron_count=0d0,electron_defect=huge(1d0)
    integer(int64),allocatable::column_ids(:),source_candidate_ids(:)
    complex(real64),allocatable::coefficients(:,:)
    real(real64),allocatable::energies(:),occupations(:)
  end type

  public::initialize_dg_hybrid_schwarz_state,extend_dg_hybrid_schwarz_state,&
    validate_dg_hybrid_schwarz_dynamic_receipt
contains
  subroutine initialize_dg_hybrid_schwarz_state(comm,fragment_id,fragment_count,basis_generation,&
      electron_target,temperature,wspin,guard_count,occupation_tail_tolerance,degeneracy_tolerance,&
      mapping_fingerprint,candidate_fingerprint,candidate_ids,candidate_energies,candidate_vectors,&
      state,ok,message,local_publish_ok)
    integer,intent(in)::comm,fragment_id,fragment_count,basis_generation,guard_count
    real(real64),intent(in)::electron_target,temperature,wspin,occupation_tail_tolerance,degeneracy_tolerance
    integer(int64),intent(in)::mapping_fingerprint,candidate_fingerprint,candidate_ids(:)
    real(real64),intent(in)::candidate_energies(:)
    complex(real64),intent(in)::candidate_vectors(:,:)
    type(s_dg_hybrid_schwarz_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::local_publish_ok
    type(s_dg_hybrid_schwarz_state)::work
    integer::rank,nproc,ierr,ncandidate,ntrial,thermal_count,stat
    integer::integer_controls(4),minimum_controls(4),maximum_controls(4)
    integer(int64)::fingerprint_controls(2),minimum_fingerprints(2),maximum_fingerprints(2)
    real(real64)::real_controls(5),minimum_reals(5),maximum_reals(5),mu
    real(real64),allocatable::minimum_energies(:),maximum_energies(:),occupations(:)
    integer(int64),allocatable::common_ids(:)
    logical::valid,publish

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz state rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz state size query failed';return;endif
    ncandidate=size(candidate_ids)
    integer_controls=[fragment_count,basis_generation,guard_count,ncandidate]
    call MPI_Allreduce(integer_controls,minimum_controls,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz integer control minimum failed';return;endif
    call MPI_Allreduce(integer_controls,maximum_controls,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz integer control maximum failed';return;endif
    fingerprint_controls=[mapping_fingerprint,candidate_fingerprint]
    call MPI_Allreduce(fingerprint_controls,minimum_fingerprints,2,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz fingerprint minimum failed';return;endif
    call MPI_Allreduce(fingerprint_controls,maximum_fingerprints,2,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz fingerprint maximum failed';return;endif
    real_controls=[electron_target,temperature,wspin,occupation_tail_tolerance,degeneracy_tolerance]
    call MPI_Allreduce(real_controls,minimum_reals,5,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz real control minimum failed';return;endif
    call MPI_Allreduce(real_controls,maximum_reals,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz real control maximum failed';return;endif
    valid=all(integer_controls==minimum_controls).and.all(integer_controls==maximum_controls).and.&
      all(fingerprint_controls==minimum_fingerprints).and.all(fingerprint_controls==maximum_fingerprints).and.&
      all(real_controls==minimum_reals).and.all(real_controls==maximum_reals)
    valid=valid.and.fragment_count==nproc.and.fragment_id==rank+1.and.basis_generation>0.and.guard_count>=0.and.&
      mapping_fingerprint/=0_int64.and.candidate_fingerprint/=0_int64.and.ncandidate>0.and.&
      electron_target>0d0.and.temperature>=0d0.and.wspin>0d0.and.occupation_tail_tolerance>0d0.and.&
      occupation_tail_tolerance<1d0.and.degeneracy_tolerance>=0d0.and.all(ieee_is_finite(real_controls))
    call collective_gate(comm,valid,'invalid or rank-disagreeing Schwarz state controls',ok,message)
    if(.not.ok)return
    valid=size(candidate_energies)==ncandidate.and.size(candidate_vectors,1)>0.and.&
      size(candidate_vectors,2)==ncandidate.and.all(candidate_ids>0_int64).and.unique_ids(candidate_ids).and.&
      all(ieee_is_finite(candidate_energies)).and.finite_matrix(candidate_vectors)
    call collective_gate(comm,valid,'invalid Schwarz candidate catalog',ok,message)
    if(.not.ok)return
    allocate(minimum_energies(ncandidate),maximum_energies(ncandidate),occupations(ncandidate),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz thermal workspace allocation failed',ok,message)
    if(.not.ok)return
    call MPI_Allreduce(candidate_energies,minimum_energies,ncandidate,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz energy minimum failed';return;endif
    call MPI_Allreduce(candidate_energies,maximum_energies,ncandidate,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.all(minimum_energies==maximum_energies),&
      'Schwarz reference energies differ between ranks',ok,message)
    if(.not.ok)return
    call thermal_inventory(candidate_energies,electron_target,temperature,wspin,guard_count,&
      occupation_tail_tolerance,degeneracy_tolerance,ntrial,thermal_count,mu,occupations,valid)
    call collective_gate(comm,valid,'Schwarz candidate capacity does not resolve the 300 K occupation tail',ok,message)
    if(.not.ok)return
    call build_common_column_ids(comm,candidate_ids,ntrial,common_ids,ok,message)
    if(.not.ok)return
    work%valid=.true.;work%fragment_id=fragment_id;work%fragment_count=fragment_count
    work%basis_generation=basis_generation;work%local_basis_count=size(candidate_vectors,1)
    work%trial_count=ntrial;work%candidate_count=ncandidate;work%thermal_tail_count=thermal_count
    work%coefficient_epoch=0
    work%mapping_fingerprint=mapping_fingerprint;work%candidate_fingerprint=candidate_fingerprint
    work%electron_target=electron_target;work%temperature=temperature;work%wspin=wspin
    work%chemical_potential=mu
    allocate(work%column_ids(ntrial),work%source_candidate_ids(ncandidate),&
      work%coefficients(size(candidate_vectors,1),ntrial),work%energies(ntrial),&
      work%occupations(ntrial),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz state staging allocation failed',ok,message)
    if(.not.ok)return
    work%column_ids=common_ids;work%source_candidate_ids=candidate_ids
    work%coefficients=candidate_vectors(:,1:ntrial)
    work%energies=candidate_energies(:ntrial);work%occupations=occupations(:ntrial)
    work%electron_count=wspin*sum(work%occupations)
    work%electron_defect=abs(work%electron_count-electron_target)
    work%fingerprint=state_fingerprint(work)
    publish=.true.;if(present(local_publish_ok))publish=local_publish_ok
    valid=publish.and.work%fingerprint/=0_int64.and.finite_matrix(work%coefficients)
    call collective_gate(comm,valid,'Schwarz state publication rolled back',ok,message)
    if(.not.ok)return
    state=work;ok=.true.;message=''
  end subroutine initialize_dg_hybrid_schwarz_state

  subroutine extend_dg_hybrid_schwarz_state(comm,basis_generation,requested_count,mapping_fingerprint,&
      candidate_fingerprint,candidate_ids,candidate_vectors,state,ok,message,local_publish_ok)
    integer,intent(in)::comm,basis_generation,requested_count
    integer(int64),intent(in)::mapping_fingerprint,candidate_fingerprint,candidate_ids(:)
    complex(real64),intent(in)::candidate_vectors(:,:)
    type(s_dg_hybrid_schwarz_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::local_publish_ok
    type(s_dg_hybrid_schwarz_state)::work
    integer::rank,nproc,ierr,minimum_count,maximum_count,stat
    integer(int64),allocatable::common_ids(:)
    complex(real64),allocatable::extended_coefficients(:,:)
    logical::valid,publish

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz extension rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz extension size query failed';return;endif
    call MPI_Allreduce(requested_count,minimum_count,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz extension count minimum failed';return;endif
    call MPI_Allreduce(requested_count,maximum_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.minimum_count==maximum_count.and.state%valid.and.&
      state%fragment_count==nproc.and.state%fragment_id==rank+1.and.&
      state%basis_generation==basis_generation.and.state%mapping_fingerprint==mapping_fingerprint.and.&
      state%candidate_fingerprint==candidate_fingerprint.and.requested_count>state%trial_count.and.&
      requested_count<=state%candidate_count.and.size(candidate_ids)==state%candidate_count.and.&
      size(candidate_vectors,1)==state%local_basis_count.and.size(candidate_vectors,2)==state%candidate_count.and.&
      allocated(state%column_ids).and.allocated(state%source_candidate_ids).and.allocated(state%coefficients)
    if(valid)valid=all(candidate_ids==state%source_candidate_ids).and.unique_ids(candidate_ids).and.&
      finite_matrix(candidate_vectors).and.all(shape(state%coefficients)==[state%local_basis_count,state%trial_count])
    call collective_gate(comm,valid,'invalid or stale Schwarz extension context',ok,message)
    if(.not.ok)return
    call build_common_column_ids(comm,candidate_ids,requested_count,common_ids,ok,message)
    if(.not.ok)return
    valid=all(common_ids(1:state%trial_count)==state%column_ids)
    call collective_gate(comm,valid,'Schwarz extension changed accepted column bindings',ok,message)
    if(.not.ok)return
    allocate(extended_coefficients(state%local_basis_count,requested_count),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz extension staging allocation failed',ok,message)
    if(.not.ok)return
    extended_coefficients(:,1:state%trial_count)=state%coefficients
    extended_coefficients(:,state%trial_count+1:requested_count)=&
      candidate_vectors(:,state%trial_count+1:requested_count)
    work=state
    deallocate(work%column_ids,work%coefficients)
    allocate(work%column_ids(requested_count),work%coefficients(state%local_basis_count,requested_count),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz extension publication allocation failed',ok,message)
    if(.not.ok)return
    work%trial_count=requested_count;work%column_ids=common_ids;work%coefficients=extended_coefficients
    if(allocated(work%energies))deallocate(work%energies)
    if(allocated(work%occupations))deallocate(work%occupations)
    work%electron_count=0d0;work%electron_defect=huge(1d0)
    work%fingerprint=state_fingerprint(work)
    publish=.true.;if(present(local_publish_ok))publish=local_publish_ok
    valid=publish.and.work%fingerprint/=0_int64.and.finite_matrix(work%coefficients)
    call collective_gate(comm,valid,'Schwarz extension publication rolled back',ok,message)
    if(.not.ok)return
    state=work;ok=.true.;message=''
  end subroutine extend_dg_hybrid_schwarz_state

  subroutine validate_dg_hybrid_schwarz_dynamic_receipt(comm,state,receipt,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_schwarz_state),intent(in)::state
    integer(int64),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::local_receipt,minimum_receipt,maximum_receipt
    integer::ierr,j,rank,nproc
    logical::valid

    receipt=0_int64;ok=.false.;message='';rank=-1;nproc=-1
    call MPI_Comm_rank(comm,rank,ierr)
    valid=ierr==MPI_SUCCESS
    call MPI_Comm_size(comm,nproc,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
    valid=valid.and.state%valid.and.state%fragment_count==nproc.and.state%fragment_id==rank+1.and.&
      state%coefficient_epoch>=0.and.&
      state%trial_count>0.and.state%candidate_count>=state%trial_count.and.&
      state%thermal_tail_count>=0.and.state%thermal_tail_count<=state%trial_count.and.&
      state%fingerprint/=0_int64.and.allocated(state%energies).and.allocated(state%occupations)
    if(valid)valid=size(state%energies)==state%trial_count.and.size(state%occupations)==state%trial_count
    if(valid)valid=all(ieee_is_finite(state%energies)).and.all(ieee_is_finite(state%occupations)).and.&
      ieee_is_finite(state%chemical_potential).and.ieee_is_finite(state%electron_count).and.&
      ieee_is_finite(state%electron_defect).and.ieee_is_finite(state%temperature).and.&
      ieee_is_finite(state%wspin)
    call collective_gate(comm,valid,'invalid Schwarz dynamic state receipt',ok,message)
    if(.not.ok)return
    local_receipt=mix_hash(int(z'510E527FADE682D1',int64),state%fingerprint)
    local_receipt=mix_hash(local_receipt,int(state%fragment_count,int64))
    local_receipt=mix_hash(local_receipt,int(state%coefficient_epoch,int64))
    local_receipt=mix_hash(local_receipt,int(state%trial_count,int64))
    local_receipt=mix_hash(local_receipt,int(state%candidate_count,int64))
    local_receipt=mix_hash(local_receipt,int(state%thermal_tail_count,int64))
    local_receipt=mix_hash(local_receipt,transfer(state%chemical_potential,0_int64))
    local_receipt=mix_hash(local_receipt,transfer(state%electron_count,0_int64))
    local_receipt=mix_hash(local_receipt,transfer(state%electron_defect,0_int64))
    local_receipt=mix_hash(local_receipt,transfer(state%temperature,0_int64))
    local_receipt=mix_hash(local_receipt,transfer(state%wspin,0_int64))
    do j=1,state%trial_count
      local_receipt=mix_hash(local_receipt,transfer(state%energies(j),0_int64))
      local_receipt=mix_hash(local_receipt,transfer(state%occupations(j),0_int64))
    enddo
    if(local_receipt==0_int64)local_receipt=1_int64
    call MPI_Allreduce(local_receipt,minimum_receipt,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    valid=ierr==MPI_SUCCESS
    call MPI_Allreduce(local_receipt,maximum_receipt,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    valid=valid.and.ierr==MPI_SUCCESS.and.minimum_receipt==maximum_receipt.and.minimum_receipt/=0_int64
    call collective_gate(comm,valid,'rank-disagreeing Schwarz dynamic state receipt',ok,message)
    if(.not.ok)return
    receipt=minimum_receipt;message=''
  end subroutine validate_dg_hybrid_schwarz_dynamic_receipt

  subroutine thermal_inventory(energies,target,temperature,wspin,guard,tail_tolerance,degeneracy_tolerance,&
      count,thermal_count,mu,occupations,ok)
    real(real64),intent(in)::energies(:),target,temperature,wspin,tail_tolerance,degeneracy_tolerance
    integer,intent(in)::guard
    integer,intent(out)::count,thermal_count
    real(real64),intent(out)::mu,occupations(:)
    logical,intent(out)::ok
    integer::iteration,j,minimum_occupied
    real(real64)::lower,upper,mid,total,kbt
    ok=.false.;count=0;thermal_count=0;mu=0d0;occupations=0d0
    if(size(energies)<1.or.size(occupations)/=size(energies).or.any(energies(2:)<energies(:size(energies)-1)))return
    if(target>wspin*real(size(energies),real64))return
    minimum_occupied=ceiling(target/wspin-64d0*epsilon(1d0))
    if(temperature==0d0)then
      mu=energies(min(max(1,minimum_occupied),size(energies)))
      occupations(:minimum_occupied)=1d0
    else
      kbt=boltzmann_hartree_per_kelvin*temperature
      lower=energies(1)-max(1d0,64d0*kbt)
      upper=energies(size(energies))+max(1d0,64d0*kbt)
      do iteration=1,256
        mid=0.5d0*(lower+upper)
        total=wspin*sum([(fermi_value((energies(j)-mid)/kbt),j=1,size(energies))])
        if(total<target)then;lower=mid;else;upper=mid;endif
      enddo
      mu=0.5d0*(lower+upper)
      occupations=[(fermi_value((energies(j)-mu)/kbt),j=1,size(energies))]
    endif
    do j=1,size(occupations)
      if(occupations(j)>tail_tolerance)thermal_count=j
    enddo
    if(occupations(size(occupations))>tail_tolerance)return
    count=max(minimum_occupied+guard,thermal_count)
    if(count<1.or.count>size(energies))return
    do while(count<size(energies))
      if(abs(energies(count+1)-energies(count))>degeneracy_tolerance)exit
      count=count+1
    enddo
    ok=count<=size(energies).and.all(ieee_is_finite(occupations)).and.ieee_is_finite(mu)
  end subroutine thermal_inventory

  pure real(real64) function fermi_value(x)result(value)
    real(real64),intent(in)::x
    if(x>=50d0)then;value=exp(-x)
    elseif(x<=-50d0)then;value=1d0
    else;value=1d0/(1d0+exp(x))
    endif
  end function fermi_value

  subroutine build_common_column_ids(comm,local_ids,count,common_ids,ok,message)
    integer,intent(in)::comm,count
    integer(int64),intent(in)::local_ids(:)
    integer(int64),allocatable,intent(out)::common_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nproc,ierr,stat,j,p
    integer(int64),allocatable::gathered(:),catalog(:,:)
    integer(int64)::hash
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='Schwarz column size query failed';return;endif
    allocate(gathered(size(local_ids)*nproc),catalog(size(local_ids),nproc),common_ids(count),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz column workspace allocation failed',ok,message)
    if(.not.ok)return
    call MPI_Allgather(local_ids,size(local_ids),MPI_INTEGER8,gathered,size(local_ids),MPI_INTEGER8,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS,'Schwarz column source exchange failed',ok,message)
    if(.not.ok)return
    catalog=reshape(gathered,shape(catalog))
    do j=1,count
      hash=mix_hash(int(j,int64),int(z'3C6EF372FE94F82B',int64))
      do p=1,nproc
        hash=mix_hash(hash,catalog(j,p))
      enddo
      if(hash==0_int64)hash=int(j,int64)
      common_ids(j)=hash
    enddo
    call collective_gate(comm,unique_ids(common_ids),'Schwarz common column IDs collide',ok,message)
  end subroutine build_common_column_ids

  integer(int64) function state_fingerprint(state)result(hash)
    type(s_dg_hybrid_schwarz_state),intent(in)::state
    integer::j
    hash=int(z'BB67AE8584CAA73B',int64)
    hash=mix_hash(hash,int(state%basis_generation,int64));hash=mix_hash(hash,int(state%trial_count,int64))
    hash=mix_hash(hash,state%mapping_fingerprint);hash=mix_hash(hash,state%candidate_fingerprint)
    hash=mix_hash(hash,transfer(state%electron_target,0_int64));hash=mix_hash(hash,transfer(state%temperature,0_int64))
    do j=1,state%trial_count;hash=mix_hash(hash,state%column_ids(j));enddo
    if(hash==0_int64)hash=1_int64
  end function state_fingerprint

  pure integer(int64) function mix_hash(hash,value)result(mixed)
    integer(int64),intent(in)::hash,value
    mixed=ieor(ishftc(hash,17),value)
    mixed=ieor(mixed,ishftc(value,41))
  end function mix_hash

  pure logical function unique_ids(ids)result(unique)
    integer(int64),intent(in)::ids(:)
    integer::i
    unique=.true.
    do i=1,size(ids)
      if(count(ids==ids(i))/=1)then;unique=.false.;return;endif
    enddo
  end function unique_ids

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
    ok=.false.;message='Schwarz collective rank query failed'
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    failed=huge(0);if(.not.local_ok)failed=rank
    call MPI_Allreduce(failed,first_failed,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz collective status reduction failed';return;endif
    if(first_failed==huge(0))then;ok=.true.;message='';return;endif
    shared='';if(rank==first_failed)shared=detail
    call MPI_Bcast(shared,len(shared),MPI_CHARACTER,first_failed,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz collective diagnostic broadcast failed';return;endif
    message=trim(shared)
  end subroutine collective_gate
end module dg_hybrid_schwarz_state
