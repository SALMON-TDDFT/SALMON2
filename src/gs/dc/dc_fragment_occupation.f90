module dc_fragment_occupation
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use occupation_kernel,only:solve_weighted_state_occupations
  implicit none
  private
  public::determine_dc_fragment_occupations,assess_dc_fragment_occupation_capacity
contains
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
      electron_count,ok,message,needs_extension,terminal_shell_complete)
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
    integer::rank,ierr,fragment,state,boundary_state,element_count,allocation_status,local_bad,global_bad
    integer::local_shape(5),minimum_shape(5),maximum_shape(5)
    integer::optional_shape(4),minimum_optional(4),maximum_optional(4)
    integer,allocatable::local_representatives(:),global_representatives(:)
    integer(int64)::element_count_64
    real(real64)::controls(4),minimum_controls(4),maximum_controls(4),guard_state_tail,&
      boundary_tail_charge,boundary_maximum_occupation,energy_scale,fragment_tail,fragment_max
    real(real64),allocatable::local_energies(:,:),global_energies(:,:),&
      local_core_norms(:,:),global_core_norms(:,:),flat_occupations(:)
    logical::kernel_ok
    character(512)::kernel_message

    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0
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
    optional_shape=[merge(1,0,present(needs_extension)),0,merge(1,0,present(terminal_shell_complete)),0]
    if(present(needs_extension))optional_shape(2)=size(needs_extension)
    if(present(terminal_shell_complete))optional_shape(4)=size(terminal_shell_complete)
    call MPI_Allreduce(optional_shape,minimum_optional,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(optional_shape,maximum_optional,4,MPI_INTEGER,MPI_MAX,comm,ierr)
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
      else if(local_shape(1)>1)then
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
    ok=.true.;message=''
  end subroutine determine_dc_fragment_occupations
end module dc_fragment_occupation
