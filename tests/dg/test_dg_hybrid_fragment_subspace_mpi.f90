program test_dg_hybrid_fragment_subspace_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_fragment_subspace
  use dc_fragment_occupation,only:determine_dc_fragment_occupations
  implicit none
  integer,parameter::n=24,m=2
  integer::comm,rank,nproc,ierr,nlocal,i,j,k,row,iterations,precondition_calls,maximum_trial,mode
  integer(int64),allocatable::ids(:)
  integer(int64)::workspace,fingerprint
  type(s_dg_hybrid_fragment_subspace_state)::state,entry,warm,cold
  complex(real64)::h(n,n),s(n,n),g(n,n),x(n,m),res(n,m),gram(m,m)
  real(real64)::values(m),relative_residual,measured,warm_residual,cold_residual,entry_residual,norm
  logical::ok,advanced,converged
  character(64)::reason
  character(256)::message
  complex(real64)::last_h_vectors(n,m)
  real(real64)::first_shifts(m)
  integer::shifted_calls
  logical::shifts_changed
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  ! The eight-rank run includes four ranks with no local rows.
  nlocal=count([(mod(row-1,min(nproc,4))==rank,row=1,n)])
  allocate(ids(nlocal));k=0
  do row=n,1,-1
    if(mod(row-1,min(nproc,4))/=rank)cycle
    k=k+1;ids(k)=row
  enddo
  h=0d0;s=0d0
  do i=1,n
    s(i,i)=1d0+0.02d0*i;h(i,i)=0.15d0*i*s(i,i)
    if(i==n)cycle
    h(i,i+1)=cmplx(-0.04d0,0.02d0,real64);h(i+1,i)=conjg(h(i,i+1))
  enddo
  state%fragment_id=7;state%basis_generation=2;state%basis_fingerprint=101_int64
  state%metric_fingerprint=203_int64;state%state_count=m
  allocate(state%vectors(nlocal,m),state%directions(nlocal,m));state%directions=0d0
  do i=1,nlocal
    row=int(ids(i))
    state%vectors(i,1)=cmplx(1d0/(row+1),0.002d0*row,real64)
    state%vectors(i,2)=cmplx(cos(real(row,real64)),sin(0.2d0*row),real64)
  enddo
  entry=state;mode=0;precondition_calls=0;maximum_trial=0
  call advance(3,1d-14,2d0)
  call require(ok.and.advanced.and..not.converged,'bounded nonconverged update: '//trim(message))
  call require(iterations==3.and.reason=='step_cap','three-step cap must be successful')
  call require(precondition_calls==3.and.maximum_trial<=3*m,'preconditioner or trial dimension contract')
  call require(state%state_count==m.and.size(state%vectors,2)==m,'occupied plus guard count changed')
  call certify()
  call require(workspace>0_int64.and.fingerprint/=0_int64,'missing bounded-update receipt')
  warm=state;cold=state;cold%directions=0d0
  h(1,1)=h(1,1)+0.002d0;h(3,3)=h(3,3)-0.001d0
  state=warm;call advance(3,1d-14,2d0);warm=state;warm_residual=relative_residual
  call require(ok.and.state%history_rank_used>0,'next epoch did not use saved P')
  call certify()
  state=cold;call advance(3,1d-14,2d0);cold_residual=relative_residual
  call require(ok.and.warm_residual<=cold_residual+1d-12,'warm epoch worse than cold on perturbed fixture')
  state=warm
  do i=1,32
    call advance(3,1d-9,2d0)
    call require(ok,'successive bounded steps failed: '//trim(message))
    if(converged)exit
  enddo
  call require(converged,'small residual directions were discarded before the intermediate target')
  call certify()
  state=warm;call advance(3,1d0,2d0)
  call require(ok.and.converged.and.iterations==0,'intermediate target did not stop early')
  call certify()

  ! A zero preconditioner is an explicit fixture: every candidate is rejected.
  state=warm;state%directions=0d0;entry=state;mode=1
  call advance(3,1d-14,1d0)
  if(rank==0.and.(.not.ok.or.advanced.or.reason/='safe_entry_retained'))&
    write(*,'(*(g0))')'safe entry status=',ok,' advanced=',advanced,' reason=',trim(reason),&
      ' residual=',relative_residual,' message=',trim(message)
  call require(ok.and..not.advanced.and.reason=='safe_entry_retained','safe entry must remain usable')
  call require(all(state%vectors==entry%vectors),'rejected update changed safe entry coefficients')
  call certify()
  mode=2;state=warm;entry=state;call advance(3,1d-14,2d0)
  call require(.not.ok.and.all(state%vectors==entry%vectors),'callback failure must be collective/transactional')
  mode=3;call advance(3,1d-14,2d0)
  call require(.not.ok,'nonfinite callback output was accepted')
  do mode=4,7
    state=warm;entry=state;call advance(3,1d-14,2d0)
    call require(.not.ok.and.all(state%vectors==entry%vectors),'H/S callback failure changed accepted state')
  enddo
  mode=0;state=warm;state%vectors(:,2)=state%vectors(:,1)
  call advance(3,1d-14,2d0);call require(.not.ok,'rank deficient entry was accepted')
  state=warm;s=-s;call advance(3,1d-14,2d0)
  call require(.not.ok,'indefinite sampled metric was accepted');s=-s
  state=warm;state%basis_generation=1;call advance(3,1d-14,2d0)
  call require(.not.ok.and..not.allocated(state%vectors),'changed generation did not invalidate cache')
  state=warm;state%metric_fingerprint=204_int64;call advance(3,1d-14,2d0)
  call require(.not.ok.and..not.allocated(state%vectors),'changed metric did not invalidate cache')
  state=warm;call advance(0,1d-14,2d0);call require(.not.ok,'zero step count accepted')
  state=warm;call advance(3,ieee_value(0d0,ieee_quiet_nan),2d0)
  call require(.not.ok,'nonfinite tolerance accepted')
  state=warm;call advance(3,1d-14,0.5d0);call require(.not.ok,'invalid residual growth factor accepted')
  state=cold;mode=0;call advance(3,1d-14,2d0)
  cold_residual=relative_residual;call require(ok,'preconditioner scale reference failed')
  state=cold;mode=9;call advance(3,1d-14,2d0)
  call require(ok.and.abs(relative_residual-cold_residual)<1d-10,&
    'uniform preconditioner scaling changed the accepted trial space')
  mode=0
  state=warm;shifted_calls=0;shifts_changed=.false.
  call advance_shifted(1)
  call require(ok.and.shifted_calls==3.and.shifts_changed,'shifted callback did not receive changing current values')
  call certify()
  state=warm;entry=state;mode=10;call advance_shifted(1)
  call require(.not.ok.and.all(state%vectors==entry%vectors),'failed shifted callback published an update')
  mode=0;state=warm;call advance_shifted(2)
  call require(.not.ok,'both preconditioner callbacks were silently accepted')
  call advance_shifted(3);call require(.not.ok,'missing preconditioner callback was accepted')
  if(nproc>1)then
    state=warm
    if(rank==0)then
      call advance_shifted(1)
    else
      call advance(3,1d-14,2d0)
    endif
    call require(.not.ok,'rank-disagreeing callback selection accepted')
  endif
  call test_measurement_only()
  call test_seed_initialization()
  call test_growth_rollback()
  call test_extension()
  call test_metric_publication()
  if(nproc>1)then
    state=warm;call advance(merge(2,3,rank==0),1d-14,2d0)
    call require(.not.ok,'rank-disagreeing controls accepted')
  endif
  state=warm
  if(nlocal>0)ids(1)=1_int64
  call advance(3,1d-14,2d0);call require(.not.ok,'duplicate/missing rows accepted')
  if(rank==0)then
    write(*,'(a,2es14.6)')'warm/cold residuals: ',warm_residual,cold_residual
    write(*,'(a,i0,a)')'PASS hybrid fragment subspace on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine test_seed_initialization()
    type(s_dg_hybrid_fragment_subspace_state)::initial,snapshot
    complex(real64)::seeds(nlocal,6)
    complex(real64)::full_seeds(nlocal,n)
    real(real64)::energies(6),occupations(6)
    integer,allocatable::selected(:)
    integer::a,b
    energies=[-1d0,-2d0,0d0,0.5d0,0.5d0,2d0]
    occupations=[2d0,2d0,0d0,0d0,0d0,0d0]
    seeds=0d0
    do a=1,nlocal
      b=int(ids(a))
      if(b<=6)seeds(a,b)=1d0/sqrt(real(s(b,b),real64))
    enddo
    call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message)
    call require(ok,'DC seed initialization failed: '//trim(message))
    call require(initial%state_count==3.and.all(selected==[2,1,3]),'occupied-plus-guard selection wrong')
    call require(maxval(abs(initial%vectors-seeds(:,selected)))<1d-12.and.&
      all(initial%directions==0d0),'initializer lost seed coefficients or created CG history')
    call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      seeds,energies,occupations,2,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message)
    call require(ok.and.initial%state_count==5,'guard boundary split a degenerate shell')
    call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message,&
      energy_cutoff=0.5d0)
    call require(ok.and.initial%state_count==5,'energy window did not retain the whole degenerate shell')
    snapshot=initial
    full_seeds=0d0
    do a=1,nlocal
      b=int(ids(a));full_seeds(a,b)=1d0/sqrt(real(s(b,b),real64))
    enddo
    call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      full_seeds,[(real(a,real64),a=1,n)],[(2d0,a=1,n)],1,1d-8,1d-10,1d-10,&
      apply_s,initial,selected,ok,message)
    call require(.not.ok.and..not.allocated(selected),'initializer silently started with the full fragment basis')
    call require(all(initial%vectors==snapshot%vectors),'full-basis rejection changed accepted state')
    call initialize_dg_hybrid_fragment_subspace(comm,n-1,ids,7,2,101_int64,203_int64,&
      full_seeds,[(real(a,real64),a=1,n)],[(2d0,a=1,n)],1,1d-8,1d-10,1d-10,&
      apply_s,initial,selected,ok,message)
    call require(.not.ok.and..not.allocated(selected),'initializer accepted more states than basis rank')
    call require(all(initial%vectors==snapshot%vectors),'over-rank rejection changed accepted state')
    do a=1,3
      if(rank==0)then
        select case(a)
        case(1);energies(1)=ieee_value(0d0,ieee_quiet_nan)
        case(2);occupations(1)=-1d0
        case(3);energies(1)=-1.1d0
        end select
      endif
      if(a==3.and.nproc==1)cycle
      call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message)
      call require(.not.ok.and..not.allocated(selected),'invalid or rank-disagreeing seed spectrum accepted')
      call require(all(initial%vectors==snapshot%vectors),'invalid seed spectrum changed accepted state')
      energies(1)=-1d0;occupations(1)=2d0
    enddo
    energies(1)=-1d0
    if(nproc>1)then
      if(rank==0)then
        call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
          seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message,energy_cutoff=0d0)
      else
        call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
          seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message)
      endif
      call require(.not.ok.and..not.allocated(selected),'rank-disagreeing optional window accepted')
    endif
    seeds(:,2)=seeds(:,1)
    call initialize_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      seeds,energies,occupations,1,1d-8,1d-10,1d-10,apply_s,initial,selected,ok,message)
    call require(.not.ok.and..not.allocated(selected),'rank-deficient seeds were accepted')
    call require(all(initial%vectors==snapshot%vectors).and.initial%state_count==snapshot%state_count,&
      'failed initialization destroyed an accepted state')
  end subroutine

  subroutine test_measurement_only()
    type(s_dg_hybrid_fragment_subspace_state)::saved
    real(real64)::old_values(m)
    real(real64)::mu,nelectron
    real(real64),allocatable::occ(:,:)
    complex(real64)::saved_h(n,n),sx(n,m)
    logical::tail(1)
    integer::failure_mode
    state=warm;saved=state;precondition_calls=0;maximum_trial=0;mode=0
    call measure()
    call require(ok.and.fingerprint/=0_int64,'measurement-only spectrum failed: '//trim(message))
    call require(all(state%vectors==saved%vectors).and.all(state%directions==saved%directions),&
      'measurement-only call changed X or P')
    call require(precondition_calls==0.and.maximum_trial==m,'measurement-only call advanced the trial space')
    call certify();old_values=values
    h=h+0.125d0*s
    call measure()
    call require(ok.and.maxval(abs(values-old_values-0.125d0))<1d-12,'measurement reused stale energies')
    call certify();h=h-0.125d0*s
    ! A potential epoch can reverse the Rayleigh order without changing X.
    ! Measure, determine common occupations, and retain coefficient association.
    saved_h=h;call gather(state%vectors);sx=matmul(s,g(:,1:m))
    h=10d0*s-8d0*matmul(sx(:,1:1),conjg(transpose(sx(:,1:1))))-&
      11d0*matmul(sx(:,2:2),conjg(transpose(sx(:,2:2))))
    call measure()
    call require(ok.and.maxval(abs(values-[2d0,-1d0]))<1d-10,'measurement did not retain crossing state order')
    call determine_dc_fragment_occupations(comm,reshape(values,[m,1]),reshape([0.25d0,0.75d0],[m,1]),&
      [rank==0],0d0,2d0,1.5d0,1d-8,mu,occ,nelectron,ok,message,tail,allow_unordered=.true.)
    call require(ok.and.maxval(abs(occ(:,1)-[0d0,2d0]))<1d-10.and..not.any(tail),&
      'measurement-to-occupation handoff lost coefficient association')
    call require(abs(nelectron-1.5d0)<1d-8.and.all(state%vectors==saved%vectors).and.&
      all(state%directions==saved%directions),'occupation handoff changed electron count or X/P')
    h=saved_h
    do failure_mode=4,7
      mode=failure_mode;call measure()
      call require(.not.ok.and.fingerprint==0_int64.and.all(values==0d0),&
        'failed measurement published a spectrum')
      call require(all(state%vectors==saved%vectors).and.all(state%directions==saved%directions),&
        'failed measurement changed the accepted cache')
    enddo
    mode=0;state%vectors(:,2)=state%vectors(:,1)
    call measure();call require(.not.ok,'measurement accepted a rank-deficient state')
    state=saved
  end subroutine
  subroutine measure()
    call measure_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      apply_h,apply_s,1d-10,state,values,relative_residual,fingerprint,ok,message)
  end subroutine
  subroutine advance_shifted(selection)
    integer,intent(in)::selection
    if(selection==2)then
      call advance_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,precondition,3,1d-14,1d-10,2d0,state,values,iterations,relative_residual,&
        converged,advanced,reason,workspace,fingerprint,ok,message,apply_shifted_preconditioner=shifted_precondition)
    else if(selection==3)then
      call advance_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,maximum_steps=3,intermediate_tolerance=1d-14,orthogonality_tolerance=1d-10,&
        allowed_residual_growth=2d0,state=state,eigenvalues=values,iterations=iterations,&
        relative_residual=relative_residual,eigensolver_converged=converged,advanced=advanced,&
        stop_reason=reason,workspace_peak_bytes=workspace,fingerprint=fingerprint,ok=ok,message=message)
    else
      call advance_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,maximum_steps=3,intermediate_tolerance=1d-14,orthogonality_tolerance=1d-10,&
        allowed_residual_growth=2d0,state=state,eigenvalues=values,iterations=iterations,&
        relative_residual=relative_residual,eigensolver_converged=converged,advanced=advanced,&
        stop_reason=reason,workspace_peak_bytes=workspace,fingerprint=fingerprint,ok=ok,message=message,&
        apply_shifted_preconditioner=shifted_precondition)
    endif
  end subroutine
  subroutine shifted_precondition(input,shifts,output,valid)
    complex(real64),intent(in)::input(:,:)
    real(real64),intent(in)::shifts(:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    complex(real64)::hx(n,m),sx(n,m),raw(n,m)
    real(real64)::expected_shifts(m)
    integer::b
    hx=matmul(h,last_h_vectors);sx=matmul(s,last_h_vectors)
    do b=1,m
      expected_shifts(b)=real(sum(conjg(last_h_vectors(:,b))*hx(:,b)),real64)
      raw(:,b)=hx(:,b)-expected_shifts(b)*sx(:,b)
    enddo
    call require(maxval(abs(expected_shifts-shifts))<1d-12,'shifted callback received stale Rayleigh values')
    ! Validate without a rank-local collective inside unequal row loops.
    valid=.true.
    if(nlocal>0)valid=maxval(abs(input-raw(int(ids),:)))<1d-12
    call require(valid,'preconditioner did not receive raw H X - epsilon S X')
    shifted_calls=shifted_calls+1
    if(shifted_calls==1)first_shifts=shifts
    if(shifted_calls>1)shifts_changed=shifts_changed.or.maxval(abs(shifts-first_shifts))>1d-8
    call precondition(input,output,valid)
    if(mode==10.and.rank==0)valid=.false.
  end subroutine
  subroutine test_metric_publication()
    type(s_dg_hybrid_fragment_candidate_catalog)::catalog
    type(s_dg_hybrid_fragment_extension_receipt)::receipt
    type(s_dg_hybrid_fragment_subspace_state)::candidate,safe
    complex(real64)::saved_h(n,n),saved_s(n,n),u(n,n),column(n),metric(n,n)
    complex(real64),allocatable::all_rows(:,:)
    real(real64)::e(n),angle,error
    integer::a,b,exponent
    logical::used_before(3)
    saved_h=h;saved_s=s
    do a=1,n
      do b=1,n
        angle=2d0*acos(-1d0)*real((a-1)*(b-1),real64)/n
        u(a,b)=cmplx(cos(angle),sin(angle),real64)/sqrt(real(n,real64))
      enddo
    enddo
    catalog%fragment_id=7;catalog%basis_generation=2;catalog%basis_fingerprint=101_int64
    catalog%metric_fingerprint=203_int64
    allocate(catalog%coefficients(nlocal,3),catalog%energies(3),catalog%ids(3),&
      catalog%source_kind(3),catalog%used(3))
    catalog%energies=0.6d0;catalog%ids=[1_int64,2_int64,3_int64];catalog%source_kind=fragment_seed
    do exponent=3,8
      do a=1,n;e(a)=10d0**(-real(exponent,real64)*real(n-a,real64)/(n-1));enddo
      h=0d0;s=0d0
      do a=1,n
        column=u(:,a)
        do b=1,n;s(:,b)=s(:,b)+e(a)*column*conjg(column(b));enddo
        h(a,a)=1d0
      enddo
      candidate=warm;candidate%directions=0d0;catalog%used=.false.;used_before=catalog%used
      do a=1,nlocal
        b=int(ids(a));candidate%vectors(a,:)=u(b,n-m+1:)/sqrt(e(n-m+1:))
        catalog%coefficients(a,:)=u(b,:3)
      enddo
      safe=candidate
      call extend_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,1d-10,1d-13,catalog,candidate,receipt,ok,message)
      if(ok)then
        call collect_rows(candidate%vectors,all_rows)
        metric(:5,:5)=matmul(conjg(transpose(all_rows)),matmul(s,all_rows))-identity(5)
        error=maxval(abs(metric(:5,:5)))
        if(rank==0.and.error>1d-13)write(*,'(a,i0,es16.7)')'unsafe metric publication exponent/error=',exponent,error
        call require(error<=1d-13,'published extension violates requested metric orthogonality')
      else
        call require(candidate%state_count==safe%state_count.and.all(candidate%vectors==safe%vectors).and.&
          all(catalog%used.eqv.used_before),'failed metric certificate changed cache or catalog')
      endif
    enddo
    h=saved_h;s=saved_s
  end subroutine
  subroutine test_growth_rollback()
    complex(real64)::saved_h(n,n),saved_s(n,n)
    integer::a,b
    saved_h=h;saved_s=s;h=0d0;s=0d0
    do a=1,n;h(a,a)=100d0;s(a,a)=1d0;enddo
    h(1,1)=1d0;h(2,2)=2d0;h(3,3)=0.5d0
    h(1,3)=0.01d0;h(3,1)=0.01d0;h(3,4)=10d0;h(4,3)=10d0
    state=warm;state%vectors=0d0;state%directions=0d0
    do a=1,nlocal
      b=int(ids(a));if(b<=m)state%vectors(a,b)=1d0
    enddo
    entry=state;call advance(3,1d-14,1d0)
    call require(ok.and..not.advanced.and.reason=='safe_entry_retained'.and.iterations==3,&
      'residual-growing trial did not roll back successfully')
    call require(all(state%vectors==entry%vectors).and.abs(relative_residual-0.01d0)<1d-6,&
      'growth rollback did not preserve entry coefficients/residual')
    call certify();h=saved_h;s=saved_s
    ! Exact invariant subspace, but its two columns are not yet Ritz eigenvectors.
    h=0d0;s=0d0
    do a=1,n;h(a,a)=real(a,real64);s(a,a)=1d0;enddo
    h(1,2)=0.2d0;h(2,1)=0.2d0
    state=entry;call advance(3,1d-12,2d0)
    call require(ok.and.converged.and.advanced,'invariant X still needs an in-space Ritz rotation')
    call certify();h=saved_h;s=saved_s
  end subroutine
  subroutine test_extension()
    type(s_dg_hybrid_fragment_candidate_catalog)::catalog,original_catalog
    type(s_dg_hybrid_fragment_extension_receipt)::receipt
    type(s_dg_hybrid_fragment_subspace_state)::base,extended,reference,rotated
    complex(real64)::rotation(n,n),saved_h(n,n),saved_s(n,n),projector(n,n),reference_projector(n,n)
    complex(real64)::bounded_projector(n,n)
    complex(real64),allocatable::all_vectors(:,:),physical(:,:)
    real(real64)::angle,spectrum(7),reference_spectrum(7),gauge_residual,reference_residual
    real(real64)::measured_spectrum(7),measured_reference(7),measured_residual,measured_reference_residual
    real(real64)::energy_table(7,2),core_weights(7,2),mu,nelectron,density(n),reference_density(n)
    real(real64),allocatable::occupation(:,:)
    logical::representative(2),tail(2)
    integer::a,b,c,variant,old_count
    saved_h=h;saved_s=s
    catalog%fragment_id=7;catalog%basis_generation=2
    catalog%basis_fingerprint=101_int64;catalog%metric_fingerprint=203_int64
    ! Saved seed eigenvalues are deliberately not in column order.  Both degenerate states must enter.
    allocate(catalog%coefficients(nlocal,5),catalog%energies(5),catalog%ids(5),&
      catalog%source_kind(5),catalog%used(5))
    catalog%energies=[0.6d0,0.6d0,1.2d0,1.2d0,0d0]
    catalog%ids=[12_int64,11_int64,21_int64,20_int64,30_int64]
    catalog%source_kind=[fragment_seed,fragment_seed,fragment_pw,fragment_pw,fragment_projector]
    catalog%used=.false.;catalog%coefficients=0d0
    base=warm;base%vectors=0d0;base%directions=0d0
    do a=1,nlocal
      b=int(ids(a))
      if(b<=m)base%vectors(a,b)=1d0/sqrt(real(s(b,b),real64))
      if(b<=m)base%directions(a,b)=0.03d0
      if(b>=3.and.b<=7)catalog%coefficients(a,b-2)=1d0/sqrt(real(s(b,b),real64))
    enddo
    original_catalog=catalog
    do variant=0,3
      rotation=0d0
      do a=1,n;rotation(a,a)=1d0;enddo
      if(variant==1)then
        do a=1,n;rotation(a,a)=cmplx(cos(0.2d0*a),sin(0.2d0*a),real64);enddo
      else if(variant==2)then
        rotation=0d0
        do a=1,n;rotation(a,n-a+1)=1d0;enddo
      else if(variant==3)then
        ! Dense unitary within the whole construction-WF block (DFT), not just a phase rotation.
        do a=1,n
          do b=1,n
            angle=2d0*acos(-1d0)*real((a-1)*(b-1),real64)/n
            rotation(a,b)=cmplx(cos(angle),sin(angle),real64)/sqrt(real(n,real64))
          enddo
        enddo
      endif
      h=matmul(conjg(transpose(rotation)),matmul(saved_h,rotation))
      s=matmul(conjg(transpose(rotation)),matmul(saved_s,rotation))
      catalog=original_catalog;extended=base
      call rotate_rows(base%vectors,rotation,extended%vectors)
      call rotate_rows(base%directions,rotation,extended%directions)
      call rotate_rows(original_catalog%coefficients,rotation,catalog%coefficients)
      rotated=extended
      call extend_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,1d-10,1d-10,catalog,extended,receipt,ok,message)
      call require(ok,'seed extension failed: '//trim(message))
      call require(extended%state_count==4.and.receipt%source_kind==fragment_seed.and.&
        receipt%old_state_count==2.and.receipt%new_state_count==4,'next complete seed eigenspace not appended')
      call require(abs(receipt%shell_lower-0.6d0)<1d-14.and.abs(receipt%shell_upper-0.6d0)<1d-14,&
        'extension changed selected shell energy')
      call require(all(extended%vectors(:,:m)==rotated%vectors).and.&
        all(extended%directions(:,:m)==rotated%directions).and.all(extended%directions(:,m+1:)==0d0),&
        'extension did not embed old X/P bitwise with zero new history')
      call collect_rows(extended%vectors,all_vectors)
      physical=matmul(rotation,all_vectors)
      projector=matmul(physical,conjg(transpose(physical)))
      if(variant==0)reference_projector=projector
      call require(maxval(abs(projector-reference_projector))<1d-10,'seed extension depends on WF gauge')
      call extend_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,1d-10,1d-10,catalog,extended,receipt,ok,message)
      call require(ok.and.extended%state_count==7.and.receipt%source_kind==fragment_pw,&
        'PW kinetic shell plus unused projector support not appended')
      call collect_rows(extended%vectors,all_vectors)
      physical=matmul(rotation,all_vectors)
      projector=matmul(physical,conjg(transpose(physical)))
      call require(maxval(abs(matmul(conjg(transpose(all_vectors)),matmul(s,all_vectors))-&
        identity(7)))<1d-10,'extended cache is not S orthonormal')
      if(variant==0)then
        reference=extended
      else
        call collect_rows(reference%vectors,all_vectors)
        call require(maxval(abs(projector-matmul(all_vectors,conjg(transpose(all_vectors)))))<1d-10,&
          'PW/projector extension depends on WF gauge')
      endif
      old_count=extended%state_count
      call extend_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,1d-10,1d-10,catalog,extended,receipt,ok,message)
      call require(.not.ok.and.index(message,'insufficient-spectrum')>0.and.extended%state_count==old_count,&
        'exhausted candidate pool did not fail collectively')
      rotated=extended
      call measure_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,1d-10,extended,measured_spectrum,measured_residual,fingerprint,ok,message)
      call require(ok.and.all(extended%vectors==rotated%vectors).and.&
        all(extended%directions==rotated%directions),'post-extension measurement changed X/P')
      call collect_rows(extended%vectors,all_vectors)
      do a=1,7
        call require(abs(measured_spectrum(a)-real(sum(conjg(all_vectors(:,a))*&
          matmul(h,all_vectors(:,a))),real64))<1d-12,'post-extension Rayleigh value mismatch')
      enddo
      if(variant==0)then
        measured_reference=measured_spectrum;measured_reference_residual=measured_residual
      else
        call require(maxval(abs(measured_spectrum-measured_reference))<1d-10.and.&
          abs(measured_residual-measured_reference_residual)<1d-10,'measurement depends on WF gauge')
      endif
      ! Explicit identity fixture is covariant under a unitary construction-WF gauge.
      call advance_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
        apply_h,apply_s,identity_precondition,3,1d-12,1d-10,2d0,extended,spectrum,iterations,&
        gauge_residual,converged,advanced,reason,workspace,fingerprint,ok,message)
      call require(ok,'gauge-rotated bounded update failed: '//trim(message))
      call collect_rows(extended%vectors,all_vectors);physical=matmul(rotation,all_vectors)
      projector=matmul(physical,conjg(transpose(physical)))
      if(variant==0)then
        bounded_projector=projector;reference_spectrum=spectrum;reference_residual=gauge_residual
      else
        call require(maxval(abs(projector-bounded_projector))<1d-8.and.&
          maxval(abs(spectrum-reference_spectrum))<1d-10.and.abs(gauge_residual-reference_residual)<1d-10,&
          'bounded projector/spectrum/residual depends on construction-WF gauge')
      endif
      energy_table(:,1)=spectrum;energy_table(:,2)=[-0.3d0,0.4d0,0.8d0,1.2d0,1.2d0,1.2d0,1.2d0]
      core_weights=0d0;core_weights(:4,2)=0.4d0
      do a=1,7
        core_weights(a,1)=0.6d0*sum([(real(saved_s(b,b),real64)*abs(physical(b,a))**2,b=1,n)])
      enddo
      representative=[rank==0,rank==mod(1,nproc)]
      call determine_dc_fragment_occupations(comm,energy_table,core_weights,representative,0d0,2d0,2d0,1d-8,&
        mu,occupation,nelectron,ok,message,tail)
      call require(ok.and..not.any(tail),'common-mu gauge fixture failed: '//trim(message))
      density=0d0
      do a=1,7
        do b=1,n
          density(b)=density(b)+0.6d0*occupation(a,1)*real(saved_s(b,b),real64)*abs(physical(b,a))**2
        enddo
      enddo
      density(:4)=density(:4)+0.4d0*occupation(:4,2)
      call require(abs(sum(density)-2d0)<1d-8,'common-mu gauge density lost electrons')
      if(variant==0)reference_density=density
      call require(maxval(abs(density-reference_density))<1d-8,'common-mu density depends on WF gauge')
    enddo
    h=saved_h;s=saved_s
  end subroutine
  function identity(count) result(matrix)
    integer,intent(in)::count
    complex(real64)::matrix(count,count)
    integer::a
    matrix=0d0
    do a=1,count;matrix(a,a)=1d0;enddo
  end function
  subroutine collect_rows(local,all_rows)
    complex(real64),intent(in)::local(:,:)
    complex(real64),allocatable,intent(out)::all_rows(:,:)
    integer::a
    allocate(all_rows(n,size(local,2)));all_rows=0d0
    do a=1,nlocal;all_rows(int(ids(a)),:)=local(a,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,all_rows,size(all_rows),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine
  subroutine rotate_rows(local,rotation,rotated)
    complex(real64),intent(in)::local(:,:),rotation(:,:)
    complex(real64),intent(out)::rotated(:,:)
    complex(real64),allocatable::all_rows(:,:),result(:,:)
    integer::a
    call collect_rows(local,all_rows);result=matmul(conjg(transpose(rotation)),all_rows)
    do a=1,nlocal;rotated(a,:)=result(int(ids(a)),:);enddo
  end subroutine
  subroutine advance(cap,tolerance,growth)
    integer,intent(in)::cap
    real(real64),intent(in)::tolerance,growth
    call advance_dg_hybrid_fragment_subspace(comm,n,ids,7,2,101_int64,203_int64,&
      apply_h,apply_s,precondition,cap,tolerance,1d-10,growth,state,values,iterations,&
      relative_residual,converged,advanced,reason,workspace,fingerprint,ok,message)
  end subroutine
  subroutine certify()
    call gather(state%vectors);x=g(:,1:m)
    gram=matmul(conjg(transpose(x)),matmul(s,x))
    do j=1,m;gram(j,j)=gram(j,j)-1d0;enddo
    call require(maxval(abs(gram))<1d-9,'returned state is not S orthonormal')
    res=matmul(h,x);measured=0d0
    do j=1,m
      norm=sqrt(sum(abs(res(:,j))**2));res(:,j)=res(:,j)-values(j)*matmul(s,x(:,j))
      measured=max(measured,sqrt(sum(abs(res(:,j))**2))/max(1d0,abs(values(j)),norm))
    enddo
    call require(abs(measured-relative_residual)<1d-12,'residual is not distributed normalized two-norm')
  end subroutine
  subroutine gather(input)
    complex(real64),intent(in)::input(:,:)
    integer::a
    g=0d0
    do a=1,nlocal;g(int(ids(a)),1:size(input,2))=input(a,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,g,size(g),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine
  subroutine apply_h(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    integer::a
    maximum_trial=max(maximum_trial,size(input,2));call gather(input)
    if(size(input,2)==m)last_h_vectors=g(:,1:m)
    do a=1,nlocal;output(a,:)=matmul(h(int(ids(a)),:),g(:,1:size(input,2)));enddo
    valid=.true.
    if(mode==4.and.rank==0)valid=.false.
    if(mode==6.and.rank==0.and.nlocal>0)output(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  end subroutine
  subroutine apply_s(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    integer::a
    call gather(input)
    do a=1,nlocal;output(a,:)=matmul(s(int(ids(a)),:),g(:,1:size(input,2)));enddo
    valid=.true.
    if(mode==5.and.rank==0)valid=.false.
    if(mode==7.and.rank==0.and.nlocal>0)output(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  end subroutine
  subroutine precondition(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    integer::a
    precondition_calls=precondition_calls+1
    do a=1,nlocal;output(a,:)=input(a,:)/(1d0+real(h(int(ids(a)),int(ids(a))),real64));enddo
    if(mode==1)output=0d0
    if(mode==9)output=1d-12*output
    valid=.not.(mode==2.and.rank==0)
    if(mode==3.and.rank==0.and.nlocal>0)output(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  end subroutine
  subroutine identity_precondition(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    output=input;valid=.true.
  end subroutine
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::bad
    call MPI_Allreduce(merge(0,1,condition),bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(bad/=0)then
      if(rank==0)write(*,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine
end program
