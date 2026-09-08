module dc_fragment_occupation
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use occupation_kernel,only:solve_weighted_state_occupations
  implicit none
  private
  public::determine_dc_fragment_occupations,assess_dc_fragment_occupation_capacity
  public::run_dc_fragment_occupation_epoch
  abstract interface
    subroutine refresh_fragment_spectrum(epoch,energies,core_weights,can_extend,ok,message)
      import real64
      integer,intent(in)::epoch
      real(real64),allocatable,intent(out)::energies(:),core_weights(:)
      logical,intent(out)::can_extend,ok
      character(*),intent(out)::message
    end subroutine
    subroutine extend_fragment_spectrum(epoch,old_count,new_count,ok,message)
      integer,intent(in)::epoch,old_count
      integer,intent(out)::new_count
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine
  end interface
contains
  subroutine run_dc_fragment_occupation_epoch(comm,fragment_count,fragment_id,representative,epoch,basis_count,&
      temperature,wspin,expected_electrons,tolerance,refresh,extend,occupations,chemical_potential,electron_count,&
      passes,extensions,ok,message)
    ! One fragment per rank, with arbitrary ranks per fragment. Callbacks use
    ! their fragment communicator. Every pass retains the same outer epoch;
    ! refresh must use its persistent bounded-update budget, never reset it.
    integer,intent(in)::comm,fragment_count,fragment_id,epoch,basis_count
    logical,intent(in)::representative
    real(real64),intent(in)::temperature,wspin,expected_electrons,tolerance
    procedure(refresh_fragment_spectrum)::refresh
    procedure(extend_fragment_spectrum)::extend
    real(real64),allocatable,intent(out)::occupations(:)
    real(real64),intent(out)::chemical_potential,electron_count
    integer,intent(out)::passes,extensions
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::values(:),weights(:),energy_table(:,:),weight_table(:,:),all_occupations(:,:)
    real(real64)::controls(4),lo(4),hi(4),mu,ne,scale
    integer::ints(2),imin(2),imax(2),ns,nmax,expected_count,new_count,stat,ierr,nproc
    integer,allocatable::representatives(:),metadata(:,:),metadata_min(:,:),metadata_max(:,:)
    logical,allocatable::mask(:),can_grow(:),needs_extension(:)
    logical::valid,extendable,sufficient
    character(512)::diagnostic
    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0;passes=0;extensions=0
    ints=[fragment_count,epoch]
    call MPI_Allreduce(ints,imin,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(ints,imax,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)return
    valid=all(imin==imax).and.fragment_count>0.and.epoch>0.and.fragment_id>=1.and.&
      fragment_id<=fragment_count.and.basis_count>0.and.fragment_count<=nproc.and.fragment_count<=huge(0)/3
    controls=[temperature,wspin,expected_electrons,tolerance]
    call epoch_status(comm,valid.and.all(ieee_is_finite(controls)),&
      'invalid occupation epoch controls',ok,message)
    if(.not.ok)return
    call MPI_Allreduce(controls,lo,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;return;endif
    call MPI_Allreduce(controls,hi,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(lo==hi).and.temperature>=0d0.and.wspin>0d0.and.&
      expected_electrons>=0d0.and.tolerance>0d0
    call epoch_status(comm,valid,'rank-disagreeing occupation epoch controls',ok,message)
    if(.not.ok)return
    allocate(representatives(fragment_count),metadata(3,fragment_count),metadata_min(3,fragment_count),&
      metadata_max(3,fragment_count),mask(fragment_count),can_grow(fragment_count),&
      needs_extension(fragment_count),stat=stat)
    call epoch_status(comm,stat==0,'cannot allocate occupation epoch directory',ok,message)
    if(.not.ok)return
    representatives=0;mask=.false.
    if(representative)then;representatives(fragment_id)=1;mask(fragment_id)=.true.;endif
    call MPI_Allreduce(MPI_IN_PLACE,representatives,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    call epoch_status(comm,ierr==MPI_SUCCESS.and.all(representatives==1),&
      'occupation epoch requires exactly one representative per fragment',ok,message)
    if(.not.ok)return
    expected_count=0
    do
      call refresh(epoch,values,weights,extendable,valid,diagnostic)
      call epoch_status(comm,valid,diagnostic,ok,message)
      if(.not.ok)return
      call epoch_status(comm,allocated(values).and.allocated(weights),&
        'fragment refresh did not publish a spectrum',ok,message)
      if(.not.ok)return
      ns=size(values)
      valid=ns>0.and.ns<=basis_count.and.size(weights)==ns
      if(expected_count>0)valid=valid.and.ns==expected_count
      call epoch_status(comm,valid,'fragment refresh changed the expected state inventory',ok,message)
      if(.not.ok)return
      valid=all(ieee_is_finite(values)).and.all(ieee_is_finite(weights))
      call epoch_status(comm,valid,'nonfinite refreshed fragment spectrum',ok,message)
      if(.not.ok)return
      call epoch_status(comm,all(weights>=0d0),'negative refreshed core weight',ok,message)
      if(.not.ok)return
      metadata=huge(0);metadata(:,fragment_id)=[ns,basis_count,merge(1,0,extendable)]
      call MPI_Allreduce(metadata,metadata_min,3*fragment_count,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;return;endif
      metadata=0;metadata(:,fragment_id)=[ns,basis_count,merge(1,0,extendable)]
      call MPI_Allreduce(metadata,metadata_max,3*fragment_count,MPI_INTEGER,MPI_MAX,comm,ierr)
      call epoch_status(comm,ierr==MPI_SUCCESS.and.all(metadata_min==metadata_max),&
        'fragment ranks disagree on refreshed inventory',ok,message)
      if(.not.ok)return
      nmax=maxval(metadata_max(1,:));can_grow=metadata_max(3,:)==1
      can_grow=can_grow.and.metadata_max(1,:)<metadata_max(2,:)
      call epoch_status(comm,int(nmax,int64)*fragment_count<=int(huge(0),int64),&
        'occupation epoch table exceeds MPI extent',ok,message)
      if(.not.ok)return
      if(allocated(energy_table))deallocate(energy_table,weight_table)
      allocate(energy_table(nmax,fragment_count),weight_table(nmax,fragment_count),stat=stat)
      call epoch_status(comm,stat==0,'cannot allocate occupation epoch tables',ok,message)
      if(.not.ok)return
      energy_table=0d0;weight_table=0d0
      if(representative)then
        energy_table(:,fragment_id)=maxval(values)
        energy_table(:ns,fragment_id)=values;weight_table(:ns,fragment_id)=weights
      endif
      call MPI_Allreduce(MPI_IN_PLACE,energy_table,nmax*fragment_count,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;return;endif
      call MPI_Allreduce(MPI_IN_PLACE,weight_table,nmax*fragment_count,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      scale=max(1d0,maxval(abs(values)),maxval(abs(energy_table(:ns,fragment_id))))
      valid=ierr==MPI_SUCCESS.and.maxval(abs(values/scale-energy_table(:ns,fragment_id)/scale))<=4096d0*epsilon(1d0)
      scale=max(1d0,maxval(weights),maxval(weight_table(:ns,fragment_id)))
      valid=valid.and.maxval(abs(weights/scale-weight_table(:ns,fragment_id)/scale))<=4096d0*epsilon(1d0)
      call epoch_status(comm,valid,'fragment ranks disagree on spectrum or core weights',ok,message)
      if(.not.ok)return
      call assess_dc_fragment_occupation_capacity(comm,weight_table,mask,can_grow,wspin,expected_electrons,&
        tolerance,sufficient,needs_extension,ok,message)
      if(.not.ok)return
      if(passes==huge(passes))then;ok=.false.;message='occupation epoch pass count overflow';return;endif
      passes=passes+1
      if(sufficient)then
        call determine_dc_fragment_occupations(comm,energy_table,weight_table,mask,temperature,wspin,&
          expected_electrons,tolerance,mu,all_occupations,ne,ok,message,needs_extension,allow_unordered=.true.)
        if(.not.ok)return
        if(.not.any(needs_extension))then
          allocate(occupations(ns),stat=stat)
          call epoch_status(comm,stat==0,'cannot publish fragment occupations',ok,message)
          if(.not.ok)then
            if(allocated(occupations))deallocate(occupations)
            return
          endif
          occupations=all_occupations(:ns,fragment_id);chemical_potential=mu;electron_count=ne
          return
        endif
      endif
      call epoch_status(comm,.not.any(needs_extension.and..not.can_grow),&
        'insufficient-spectrum tail: requested fragment is exhausted',ok,message)
      if(.not.ok)return
      expected_count=ns;valid=.true.;diagnostic=''
      if(needs_extension(fragment_id))then
        call extend(epoch,ns,new_count,valid,diagnostic)
        if(valid)then
          valid=new_count>ns.and.new_count<=basis_count
          if(.not.valid)diagnostic='fragment extension did not grow within its basis rank'
        endif
        if(valid)then;expected_count=new_count;extensions=extensions+1;endif
      endif
      call epoch_status(comm,valid,diagnostic,ok,message)
      if(.not.ok)return
    enddo
  end subroutine run_dc_fragment_occupation_epoch

  subroutine epoch_status(comm,local_ok,local_message,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    character(*),intent(in)::local_message
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,ierr,failed,first_failed
    character(512)::diagnostic
    call MPI_Comm_rank(comm,rank,ierr)
    failed=huge(0);if(.not.local_ok)failed=rank
    call MPI_Allreduce(failed,first_failed,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.first_failed==huge(0);message=''
    if(ok)return
    if(ierr/=MPI_SUCCESS)then;message='occupation epoch status reduction failed';return;endif
    diagnostic=''
    if(rank==first_failed)diagnostic=local_message
    call MPI_Bcast(diagnostic,len(diagnostic),MPI_CHARACTER,first_failed,comm,ierr)
    message=trim(diagnostic)
  end subroutine epoch_status

  subroutine assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
      wspin,expected_electrons,tolerance,capacity_sufficient,needs_extension,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::core_norms(:,:),wspin,expected_electrons,tolerance
    logical,intent(in)::representative_mask(:),can_extend(:)
    logical,intent(out)::capacity_sufficient,needs_extension(:),ok
    character(*),intent(out)::message
    integer::shape_values(5),lo(5),hi(5),ierr,bad,global_bad,f,stat
    integer,allocatable::representatives(:),extension(:)
    real(real64)::controls(3),minimum_controls(3),maximum_controls(3),capacity,local_capacity
    ok=.false.;message='invalid fragment capacity preflight';capacity_sufficient=.false.;needs_extension=.false.
    shape_values=[size(core_norms,1),size(core_norms,2),size(representative_mask),size(can_extend),size(needs_extension)]
    call MPI_Allreduce(shape_values,lo,5,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(shape_values,hi,5,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(lo/=hi))return
    if(lo(1)<1.or.lo(2)<1.or.any(lo(3:)/=lo(2)))return
    bad=0
    if(.not.all(ieee_is_finite([wspin,expected_electrons,tolerance])))bad=1
    do f=1,lo(2)
      if(.not.representative_mask(f))cycle
      if(.not.all(ieee_is_finite(core_norms(:,f))))bad=1
    enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    controls=[wspin,expected_electrons,tolerance]
    call MPI_Allreduce(controls,minimum_controls,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(controls,maximum_controls,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_controls/=maximum_controls))return
    if(wspin<=0d0.or.expected_electrons<0d0.or.tolerance<=0d0)return
    allocate(representatives(lo(2)),extension(lo(2)),stat=stat)
    bad=merge(0,1,stat==0)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    representatives=merge(1,0,representative_mask)
    call MPI_Allreduce(MPI_IN_PLACE,representatives,lo(2),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(representatives/=1))then
      message='each capacity fragment requires exactly one representative';return
    endif
    extension=0;local_capacity=0d0;bad=0
    do f=1,lo(2)
      if(.not.representative_mask(f))cycle
      extension(f)=merge(1,0,can_extend(f))
      if(any(core_norms(:,f)<0d0))bad=1
      ! Preflight must also reject finite values whose sum/product could overflow.
      if(maxval(abs(core_norms(:,f)))>huge(1d0)/real(size(core_norms),real64)/max(1d0,wspin))then
        bad=1
      else
        local_capacity=local_capacity+wspin*sum(core_norms(:,f))
      endif
    enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(local_capacity,capacity,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,extension,lo(2),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(.not.ieee_is_finite(capacity))return
    capacity_sufficient=capacity>=expected_electrons-tolerance
    if(.not.capacity_sufficient)then
      needs_extension=extension==1
      if(.not.any(needs_extension))then
        message='insufficient-spectrum capacity: every fragment is exhausted';return
      endif
    endif
    ok=.true.;message=''
  end subroutine assess_dc_fragment_occupation_capacity

  subroutine determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
      temperature,wspin,expected_electrons,tolerance,chemical_potential,occupations,&
      electron_count,ok,message,needs_extension,terminal_shell_complete,allow_unordered)
    ! Columns identify ascending fragment spectra.  A zero core norm may pad a
    ! shorter spectrum at its final energy; exactly one rank contributes each
    ! column on the total communicator.
    integer,intent(in)::comm
    real(real64),intent(in)::energies(:,:),core_norms(:,:),temperature,wspin,&
      expected_electrons,tolerance
    logical,intent(in)::representative_mask(:)
    real(real64),intent(out)::chemical_potential,electron_count
    real(real64),allocatable,intent(out)::occupations(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,optional,intent(out)::needs_extension(:)
    ! True only with an external whole-shell completeness receipt.  Without it,
    ! an occupied zero-T boundary is conservatively treated as incomplete.
    logical,optional,intent(in)::terminal_shell_complete(:)
    ! Bounded updates preserve X/P column order, not spectral order. Opt in to
    ! sorting energy/weight pairs internally and returning occupations in X order.
    logical,optional,intent(in)::allow_unordered
    integer::rank,ierr,fragment,state,boundary_state,element_count,allocation_status,local_bad,global_bad
    integer::local_shape(5),minimum_shape(5),maximum_shape(5)
    integer::optional_shape(6),minimum_optional(6),maximum_optional(6),position,original_index
    integer,allocatable::local_representatives(:),global_representatives(:)
    integer,allocatable::spectral_order(:,:)
    integer(int64)::element_count_64
    real(real64)::controls(4),minimum_controls(4),maximum_controls(4),guard_state_tail,&
      boundary_tail_charge,boundary_maximum_occupation,energy_scale,fragment_tail,fragment_max,energy_value,weight_value
    real(real64),allocatable::local_energies(:,:),global_energies(:,:),&
      local_core_norms(:,:),global_core_norms(:,:),flat_occupations(:)
    logical::kernel_ok,reorder_spectrum
    character(512)::kernel_message

    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0
    reorder_spectrum=.false.
    if(present(allow_unordered))reorder_spectrum=allow_unordered
    if(present(needs_extension))needs_extension=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation communicator query failed';return;endif

    local_shape=[size(energies,1),size(energies,2),size(core_norms,1),&
      size(core_norms,2),size(representative_mask)]
    call MPI_Allreduce(local_shape,minimum_shape,5,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation shape agreement failed';return;endif
    call MPI_Allreduce(local_shape,maximum_shape,5,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation shape agreement failed';return;endif
    if(any(minimum_shape/=maximum_shape))then
      message='rank-disagreeing fragment occupation shapes';return
    endif
    if(local_shape(1)<1.or.local_shape(2)<1.or.local_shape(3)/=local_shape(1).or.&
      local_shape(4)/=local_shape(2).or.local_shape(5)/=local_shape(2))then
      message='invalid fragment occupation shapes';return
    endif
    optional_shape=[merge(1,0,present(needs_extension)),0,merge(1,0,present(terminal_shell_complete)),0,&
      merge(1,0,present(allow_unordered)),merge(1,0,reorder_spectrum)]
    if(present(needs_extension))optional_shape(2)=size(needs_extension)
    if(present(terminal_shell_complete))optional_shape(4)=size(terminal_shell_complete)
    call MPI_Allreduce(optional_shape,minimum_optional,6,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(optional_shape,maximum_optional,6,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_optional/=maximum_optional))then
      message='rank-disagreeing fragment tail options';return
    endif
    if((optional_shape(1)==1.and.optional_shape(2)/=local_shape(2)).or.&
      (optional_shape(3)==1.and.optional_shape(4)/=local_shape(2)))then
      message='invalid fragment tail option shapes';return
    endif
    element_count_64=size(energies,kind=int64)
    if(element_count_64>int(huge(element_count),int64))then
      message='fragment occupation collective count exceeds MPI range';return
    endif
    element_count=int(element_count_64)

    local_bad=0
    if(.not.ieee_is_finite(temperature).or..not.ieee_is_finite(wspin).or.&
      .not.ieee_is_finite(expected_electrons).or..not.ieee_is_finite(tolerance))then
      local_bad=1
    else if(temperature<0d0.or.wspin<=0d0.or.expected_electrons<0d0.or.tolerance<=0d0)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation control validation failed';return;endif
    if(global_bad/=0)then;message='invalid fragment occupation controls';return;endif
    controls=[temperature,wspin,expected_electrons,tolerance]
    call MPI_Allreduce(controls,minimum_controls,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation control agreement failed';return;endif
    call MPI_Allreduce(controls,maximum_controls,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation control agreement failed';return;endif
    if(any(minimum_controls/=maximum_controls))then
      message='rank-disagreeing fragment occupation controls';return
    endif

    allocate(local_representatives(local_shape(2)),global_representatives(local_shape(2)),&
      local_energies(local_shape(1),local_shape(2)),global_energies(local_shape(1),local_shape(2)),&
      local_core_norms(local_shape(1),local_shape(2)),global_core_norms(local_shape(1),local_shape(2)),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation allocation agreement failed';return;endif
    if(global_bad/=0)then;message='cannot allocate fragment occupation collective workspace';return;endif

    local_representatives=merge(1,0,representative_mask)
    call MPI_Allreduce(local_representatives,global_representatives,local_shape(2),MPI_INTEGER,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment representative reduction failed';return;endif
    if(any(global_representatives/=1))then
      message='each fragment requires exactly one representative';return
    endif

    local_bad=0
    do fragment=1,local_shape(2)
      if(.not.representative_mask(fragment))cycle
      if(.not.all(ieee_is_finite(energies(:,fragment))).or.&
        .not.all(ieee_is_finite(core_norms(:,fragment))))then
        local_bad=1
      else if(any(core_norms(:,fragment)<0d0))then
        local_bad=1
      else if(local_shape(1)>1.and..not.reorder_spectrum)then
        if(any(energies(2:,fragment)<energies(:local_shape(1)-1,fragment)))local_bad=1
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spectrum validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid representative fragment spectrum or core norm';return;endif

    local_energies=0d0;local_core_norms=0d0
    do fragment=1,local_shape(2)
      if(.not.representative_mask(fragment))cycle
      local_energies(:,fragment)=energies(:,fragment)
      local_core_norms(:,fragment)=core_norms(:,fragment)
    enddo
    call MPI_Allreduce(local_energies,global_energies,element_count,MPI_DOUBLE_PRECISION,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spectrum gathering failed';return;endif
    call MPI_Allreduce(local_core_norms,global_core_norms,element_count,MPI_DOUBLE_PRECISION,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment core norm gathering failed';return;endif
    if(.not.all(ieee_is_finite(global_energies)).or.&
      .not.all(ieee_is_finite(global_core_norms)))then
      message='non-finite gathered fragment occupation data';return
    endif

    if(reorder_spectrum)then
      allocate(spectral_order(local_shape(1),local_shape(2)),stat=allocation_status)
      local_bad=merge(0,1,allocation_status==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate spectral permutation';return;endif
      do fragment=1,local_shape(2)
        spectral_order(:,fragment)=[(state,state=1,local_shape(1))]
        ! Stable ordering keeps exactly degenerate coefficient labels together.
        ! No arithmetic or WF rotation changes energies, weights or coefficients.
        do state=2,local_shape(1)
          energy_value=global_energies(state,fragment);weight_value=global_core_norms(state,fragment)
          original_index=spectral_order(state,fragment);position=state
          do while(position>1)
            if(global_energies(position-1,fragment)<=energy_value)exit
            global_energies(position,fragment)=global_energies(position-1,fragment)
            global_core_norms(position,fragment)=global_core_norms(position-1,fragment)
            spectral_order(position,fragment)=spectral_order(position-1,fragment)
            position=position-1
          enddo
          global_energies(position,fragment)=energy_value;global_core_norms(position,fragment)=weight_value
          spectral_order(position,fragment)=original_index
        enddo
      enddo
    endif

    kernel_ok=.false.;kernel_message='';chemical_potential=0d0;electron_count=0d0
    if(rank==0)then
      call solve_weighted_state_occupations(reshape(global_energies,[element_count]),&
        reshape(global_core_norms,[element_count]),expected_electrons,temperature,wspin,&
        flat_occupations,chemical_potential,electron_count,kernel_ok,kernel_message)
    endif
    call MPI_Bcast(kernel_ok,1,MPI_LOGICAL,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation status broadcast failed';return;endif
    call MPI_Bcast(kernel_message,len(kernel_message),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation diagnostic broadcast failed';return;endif
    if(.not.kernel_ok)then;message=trim(kernel_message);return;endif
    call MPI_Bcast(chemical_potential,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment chemical potential broadcast failed';return;endif
    call MPI_Bcast(electron_count,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment electron count broadcast failed';return;endif

    allocate(occupations(local_shape(1),local_shape(2)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation result allocation agreement failed';return;endif
    if(global_bad/=0)then;message='cannot allocate fragment occupation result';return;endif
    if(rank==0)occupations=reshape(flat_occupations,shape(occupations))
    ! The total-communicator broadcast also reaches every rank of each fragment.
    call MPI_Bcast(occupations,element_count,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupations broadcast failed';return;endif

    local_bad=0
    if(.not.ieee_is_finite(chemical_potential).or..not.ieee_is_finite(electron_count).or.&
      .not.all(ieee_is_finite(occupations)))then
      local_bad=1
    else if(any(occupations<0d0).or.any(occupations>wspin).or.&
      abs(electron_count-expected_electrons)>tolerance.or.&
      abs(sum(occupations*global_core_norms)-electron_count)>tolerance)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation result validation failed';return;endif
    if(global_bad/=0)then
      message='fragment occupations violate bounds or electron tolerance';return
    endif
    if(temperature>0d0.or.present(needs_extension))then
      boundary_tail_charge=0d0;boundary_maximum_occupation=0d0
      energy_scale=max(1d0,maxval(abs(global_energies)))
      do fragment=1,local_shape(2)
        fragment_tail=0d0;fragment_max=0d0
        boundary_state=0
        do state=local_shape(1),1,-1
          if(global_core_norms(state,fragment)>0d0)then
            boundary_state=state;exit
          endif
        enddo
        if(boundary_state>0)then
          do state=1,local_shape(1)
            if(global_core_norms(state,fragment)<=0d0)cycle
            if(abs(global_energies(state,fragment)/energy_scale-&
                global_energies(boundary_state,fragment)/energy_scale)>4096d0*epsilon(1d0))cycle
            fragment_tail=fragment_tail+&
              occupations(state,fragment)*global_core_norms(state,fragment)
            fragment_max=max(fragment_max,occupations(state,fragment))
          enddo
        endif
        boundary_tail_charge=boundary_tail_charge+fragment_tail
        boundary_maximum_occupation=max(boundary_maximum_occupation,fragment_max)
        if(present(needs_extension))then
          if(temperature>0d0)then
            needs_extension(fragment)=fragment_tail>tolerance/real(local_shape(2),real64)
          else if(boundary_state>0)then
            needs_extension(fragment)=fragment_max>0d0.or.&
              global_energies(boundary_state,fragment)/energy_scale<=chemical_potential/energy_scale+4096d0*epsilon(1d0)
            if(present(terminal_shell_complete))then
              ! Only the representative owns this fragment's completeness receipt.
              local_bad=0
              if(representative_mask(fragment))local_bad=merge(1,0,terminal_shell_complete(fragment))
              call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_SUM,comm,ierr)
              if(ierr/=MPI_SUCCESS)then;message='terminal completeness reduction failed';return;endif
              if(global_bad==1)needs_extension(fragment)=.false.
            endif
          endif
        endif
      enddo
      guard_state_tail=max(boundary_tail_charge,boundary_maximum_occupation)
      if(guard_state_tail>tolerance.and..not.present(needs_extension))then
        write(message,'(a,2(es24.16,a))')'finite-temperature guard-state tail exceeds tolerance: tail=',&
          guard_state_tail,' tolerance=',tolerance,''
        return
      endif
    endif
    if(reorder_spectrum)then
      local_energies=occupations
      do fragment=1,local_shape(2)
        occupations(spectral_order(:,fragment),fragment)=local_energies(:,fragment)
      enddo
    endif
    ok=.true.;message=''
  end subroutine determine_dc_fragment_occupations
end module dc_fragment_occupation
