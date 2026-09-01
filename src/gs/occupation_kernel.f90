module occupation_kernel
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::solve_spectrum_occupations
contains
  subroutine solve_spectrum_occupations(eigenvalues,k_weights,electron_target,electronic_temperature,&
      spin_orbit,occupations,chemical_potential,electron_count,ok,message)
    real(real64),intent(in)::eigenvalues(:,:,:),k_weights(:),electron_target,electronic_temperature
    logical,intent(in)::spin_orbit
    real(real64),allocatable,intent(out)::occupations(:,:,:)
    real(real64),intent(out)::chemical_potential,electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nstate,nk,nspin,expansion,iteration,allocation_status
    real(real64)::spin_weight,spin_orbit_weight,capacity,lower_mu,upper_mu,trial_mu,&
      lower_count,upper_count,trial_count,previous_count,previous_mu,relative_change,mu_change
    real(real64),allocatable::trial_occupations(:,:,:),exact_occupations(:,:,:)
    real(real64)::exact_mu,exact_count
    logical::have_exact

    ok=.false.;message='';chemical_potential=0d0;electron_count=0d0
    nstate=size(eigenvalues,1);nk=size(eigenvalues,2);nspin=size(eigenvalues,3)
    if(nstate<1.or.nk<1.or.(nspin/=1.and.nspin/=2).or.size(k_weights)/=nk)then
      message='invalid occupation spectrum shape or spin convention';return
    endif
    if(.not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(k_weights)).or.&
      .not.ieee_is_finite(electron_target).or..not.ieee_is_finite(electronic_temperature).or.&
      any(k_weights<0d0).or.sum(k_weights)<=0d0.or.electron_target<0d0.or.electronic_temperature<0d0)then
      message='invalid occupation spectrum values';return
    endif
    spin_weight=merge(2d0,1d0,nspin==1)
    spin_orbit_weight=merge(0.5d0,1d0,spin_orbit)
    capacity=spin_weight*sum(k_weights)*real(nstate*nspin,real64)*spin_orbit_weight
    if(electron_target>capacity+1d-10*max(1d0,capacity))then
      write(message,'(a,2(es24.16,a))')'occupation spectrum capacity is insufficient: requested=',&
        electron_target,' capacity=',capacity,''
      return
    endif
    allocate(occupations(nstate,nk,nspin),trial_occupations(nstate,nk,nspin),&
      exact_occupations(nstate,nk,nspin),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate occupation spectrum workspace';return;endif
    have_exact=.false.;exact_mu=0d0;exact_count=0d0
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
    message='constant-electron occupation solve did not converge'
  contains
    subroutine remember_exact(mu,values,total)
      real(real64),intent(in)::mu,values(:,:,:),total
      if(abs(total-electron_target)>=1d-9)return
      have_exact=.true.;exact_mu=mu;exact_count=total;exact_occupations=values
    end subroutine remember_exact
    subroutine evaluate_count(mu,values,total)
      real(real64),intent(in)::mu
      real(real64),intent(out)::values(:,:,:),total
      integer::state,kpoint,spin
      real(real64)::argument
      total=0d0
      do spin=1,nspin;do kpoint=1,nk;do state=1,nstate
        if(electronic_temperature==0d0)then
          values(state,kpoint,spin)=merge(0d0,spin_weight,eigenvalues(state,kpoint,spin)-mu>0d0)
        else
          argument=(eigenvalues(state,kpoint,spin)-mu)/electronic_temperature
          if(argument>=40d0)then
            values(state,kpoint,spin)=0d0
          else
            values(state,kpoint,spin)=spin_weight/(1d0+exp(argument))
          endif
        endif
        total=total+values(state,kpoint,spin)*k_weights(kpoint)
      enddo;enddo;enddo
      total=total*spin_orbit_weight
    end subroutine evaluate_count
  end subroutine solve_spectrum_occupations
end module occupation_kernel
