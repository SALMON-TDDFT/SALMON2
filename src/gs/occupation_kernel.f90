module occupation_kernel
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::solve_spectrum_occupations,solve_weighted_state_occupations
contains
  subroutine solve_spectrum_occupations(eigenvalues,k_weights,electron_target,electronic_temperature,&
      spin_orbit,occupations,chemical_potential,electron_count,ok,message)
    real(real64),intent(in)::eigenvalues(:,:,:),k_weights(:),electron_target,electronic_temperature
    logical,intent(in)::spin_orbit
    real(real64),allocatable,intent(out)::occupations(:,:,:)
    real(real64),intent(out)::chemical_potential,electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nstate,nk,nspin,state,kpoint,spin,flat_index,flat_count,allocation_status
    real(real64)::spin_weight,spin_orbit_weight
    real(real64),allocatable::flat_eigenvalues(:),flat_weights(:),flat_occupations(:)
    logical::weighted_ok
    character(256)::weighted_message

    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0
    nstate=size(eigenvalues,1);nk=size(eigenvalues,2);nspin=size(eigenvalues,3)
    if(nstate<1.or.nk<1.or.(nspin/=1.and.nspin/=2).or.size(k_weights)/=nk)then
      message='invalid occupation spectrum shape or spin convention';return
    endif
    if(.not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(k_weights)).or.&
      .not.ieee_is_finite(electron_target).or..not.ieee_is_finite(electronic_temperature))then
      message='invalid occupation spectrum values';return
    endif
    if(any(k_weights<0d0).or.sum(k_weights)<=0d0.or.electron_target<0d0.or.electronic_temperature<0d0)then
      message='invalid occupation spectrum values';return
    endif
    spin_weight=merge(2d0,1d0,nspin==1)
    spin_orbit_weight=merge(0.5d0,1d0,spin_orbit)
    flat_count=nstate*nk*nspin
    allocate(flat_eigenvalues(flat_count),flat_weights(flat_count),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate occupation spectrum workspace';return;endif
    flat_index=0
    do spin=1,nspin;do kpoint=1,nk;do state=1,nstate
      flat_index=flat_index+1
      flat_eigenvalues(flat_index)=eigenvalues(state,kpoint,spin)
      flat_weights(flat_index)=k_weights(kpoint)*spin_orbit_weight
    enddo;enddo;enddo
    call solve_weighted_state_occupations(flat_eigenvalues,flat_weights,electron_target,&
      electronic_temperature,spin_weight,flat_occupations,chemical_potential,electron_count,&
      weighted_ok,weighted_message)
    if(.not.weighted_ok)then;message=trim(weighted_message);return;endif
    allocate(occupations(nstate,nk,nspin),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate occupation spectrum workspace';return;endif
    occupations=reshape(flat_occupations,[nstate,nk,nspin])
    ok=.true.;message=''
  end subroutine solve_spectrum_occupations

  subroutine solve_weighted_state_occupations(eigenvalues,state_weights,electron_target,&
      electronic_temperature,maximum_occupation,occupations,chemical_potential,electron_count,&
      ok,message)
    real(real64),intent(in)::eigenvalues(:),state_weights(:),electron_target,&
      electronic_temperature,maximum_occupation
    real(real64),allocatable,intent(out)::occupations(:)
    real(real64),intent(out)::chemical_potential,electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nstate,expansion,iteration,allocation_status
    real(real64)::capacity,lower_mu,upper_mu,trial_mu,lower_count,upper_count,trial_count,&
      previous_count,previous_mu,relative_change,mu_change
    real(real64),allocatable::trial_occupations(:),exact_occupations(:)
    real(real64)::exact_mu,exact_count
    logical::have_exact,fractional_required,fractional_ok

    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0
    nstate=size(eigenvalues)
    if(nstate<1.or.size(state_weights)/=nstate)then
      message='invalid weighted occupation spectrum shape';return
    endif
    if(.not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(state_weights)).or.&
      .not.ieee_is_finite(electron_target).or..not.ieee_is_finite(electronic_temperature).or.&
      .not.ieee_is_finite(maximum_occupation))then
      message='invalid weighted occupation spectrum values';return
    endif
    if(any(state_weights<0d0).or.sum(state_weights)<=0d0.or.electron_target<0d0.or.&
      electronic_temperature<0d0.or.maximum_occupation<=0d0)then
      message='invalid weighted occupation spectrum values';return
    endif
    capacity=maximum_occupation*sum(state_weights)
    if(electron_target>capacity+1d-10*max(1d0,capacity))then
      write(message,'(a,2(es24.16,a))')'occupation spectrum capacity is insufficient: requested=',&
        electron_target,' capacity=',capacity,''
      return
    endif
    allocate(occupations(nstate),trial_occupations(nstate),exact_occupations(nstate),&
      stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate occupation spectrum workspace';return;endif
    have_exact=.false.;exact_mu=0d0;exact_count=0d0
    if(electronic_temperature==0d0)then
      call fill_zero_temperature_ensemble(occupations,chemical_potential,electron_count,&
        fractional_required,fractional_ok)
      if(.not.fractional_ok)then
        message='cannot construct zero-temperature occupation ensemble';return
      endif
      if(fractional_required)then;ok=.true.;return;endif
    endif
    do expansion=0,51
      lower_mu=minval(eigenvalues)-0.2d0*real(expansion,real64)
      upper_mu=maxval(eigenvalues)+0.2d0*real(expansion,real64)
      previous_count=0d0;previous_mu=0d0
      call evaluate_count(lower_mu,occupations,lower_count)
      call remember_exact(lower_mu,occupations,lower_count)
      call evaluate_count(upper_mu,trial_occupations,upper_count)
      call remember_exact(upper_mu,trial_occupations,upper_count)
      do iteration=1,999
        trial_mu=lower_mu+0.5d0*(upper_mu-lower_mu)
        call evaluate_count(trial_mu,trial_occupations,trial_count)
        call remember_exact(trial_mu,trial_occupations,trial_count)
        if(trial_count==0d0)then
          relative_change=merge(0d0,huge(1d0),previous_count==0d0)
        else
          relative_change=abs((trial_count-previous_count)/trial_count)
        endif
        mu_change=trial_mu-previous_mu
        if(abs(trial_count-electron_target)<1d-9.and.relative_change<1d-10.and.mu_change<1d-9)then
          occupations=trial_occupations;chemical_potential=trial_mu;electron_count=trial_count
          ok=.true.;return
        endif
        if((lower_count-electron_target)*(trial_count-electron_target)>0d0)then
          lower_mu=trial_mu;lower_count=trial_count
        else
          upper_mu=trial_mu;upper_count=trial_count
        endif
        previous_count=trial_count;previous_mu=trial_mu
      enddo
    enddo
    if(electronic_temperature==0d0.and.have_exact)then
      occupations=exact_occupations;chemical_potential=exact_mu;electron_count=exact_count
      ok=.true.;return
    endif
    if(electronic_temperature==0d0)then
      call fill_zero_temperature_ensemble(occupations,chemical_potential,electron_count,&
        fractional_required,fractional_ok)
      if(fractional_ok)then;ok=.true.;return;endif
    endif
    message='constant-electron occupation solve did not converge'
  contains
    subroutine fill_zero_temperature_ensemble(values,mu,total,requires_fractional,success)
      real(real64),intent(out)::values(:),mu,total
      logical,intent(out)::requires_fractional,success
      integer,allocatable::order(:)
      integer::i,first,last,allocation_status
      real(real64)::completed_count,shell_weight,shell_capacity,remaining,shell_occupation,energy_scale

      values=0d0;mu=0d0;total=0d0;requires_fractional=.false.;success=.false.
      allocate(order(nstate),stat=allocation_status)
      if(allocation_status/=0)return
      order=[(i,i=1,nstate)]
      call sort_state_indices(order)
      energy_scale=max(1d0,maxval(abs(eigenvalues)))
      completed_count=0d0;first=1
      do while(first<=nstate)
        last=first
        do while(last<nstate)
          if(abs(eigenvalues(order(last+1))/energy_scale-&
              eigenvalues(order(first))/energy_scale)>4096d0*epsilon(1d0))exit
          last=last+1
        enddo
        shell_weight=sum(state_weights(order(first:last)))
        if(shell_weight<=0d0)then
          values(order(first:last))=maximum_occupation
          first=last+1;cycle
        endif
        shell_capacity=maximum_occupation*shell_weight
        if(electron_target<=completed_count+shell_capacity)then
          remaining=max(0d0,min(shell_capacity,electron_target-completed_count))
          if(remaining<1d-9)then
            shell_occupation=0d0
          elseif(abs(remaining-shell_capacity)<1d-9)then
            shell_occupation=maximum_occupation
          else
            shell_occupation=remaining/shell_weight
            requires_fractional=.true.
          endif
          values(order(first:last))=shell_occupation
          mu=eigenvalues(order(first))
          total=sum(values*state_weights)
          success=abs(total-electron_target)<1d-9
          return
        endif
        values(order(first:last))=maximum_occupation
        completed_count=completed_count+shell_capacity
        first=last+1
      enddo
      if(abs(completed_count-electron_target)<1d-9)then
        mu=maxval(eigenvalues);total=sum(values*state_weights);success=.true.
      endif
    end subroutine fill_zero_temperature_ensemble

    subroutine sort_state_indices(order)
      integer,intent(inout)::order(:)
      integer::start,finish,root,child,temporary
      do start=size(order)/2,1,-1
        root=start
        do while(2*root<=size(order))
          child=2*root
          if(child<size(order))then
            if(state_follows(order(child+1),order(child)))child=child+1
          endif
          if(.not.state_follows(order(child),order(root)))exit
          temporary=order(root);order(root)=order(child);order(child)=temporary
          root=child
        enddo
      enddo
      do finish=size(order),2,-1
        temporary=order(1);order(1)=order(finish);order(finish)=temporary
        root=1
        do while(2*root<=finish-1)
          child=2*root
          if(child<finish-1)then
            if(state_follows(order(child+1),order(child)))child=child+1
          endif
          if(.not.state_follows(order(child),order(root)))exit
          temporary=order(root);order(root)=order(child);order(child)=temporary
          root=child
        enddo
      enddo
    end subroutine sort_state_indices

    logical function state_follows(left,right)
      integer,intent(in)::left,right
      state_follows=eigenvalues(left)>eigenvalues(right).or.&
        (eigenvalues(left)==eigenvalues(right).and.left>right)
    end function state_follows

    subroutine remember_exact(mu,values,total)
      real(real64),intent(in)::mu,values(:),total
      if(abs(total-electron_target)>=1d-9)return
      have_exact=.true.;exact_mu=mu;exact_count=total;exact_occupations=values
    end subroutine remember_exact
    subroutine evaluate_count(mu,values,total)
      real(real64),intent(in)::mu
      real(real64),intent(out)::values(:),total
      integer::state
      real(real64)::argument
      total=0d0
      do state=1,nstate
        if(electronic_temperature==0d0)then
          values(state)=merge(0d0,maximum_occupation,eigenvalues(state)-mu>0d0)
        else
          argument=(eigenvalues(state)-mu)/electronic_temperature
          if(argument>=40d0)then
            values(state)=0d0
          else
            values(state)=maximum_occupation/(1d0+exp(argument))
          endif
        endif
        total=total+values(state)*state_weights(state)
      enddo
    end subroutine evaluate_count
  end subroutine solve_weighted_state_occupations
end module occupation_kernel
