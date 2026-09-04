module dg_hybrid_fragment_subspace
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,&
    ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  implicit none
  private
  type,public::s_dg_hybrid_fragment_epoch_budget
    private
    integer::epoch=0,limit=0,remaining=0,fragment_id=0,generation=0,global_count=0
    integer::comm_rank=-1,comm_size=0
    integer(int64)::basis_fp=0_int64,metric_fp=0_int64
    logical::failed=.false.
  end type
  type,public::s_dg_hybrid_fragment_subspace_state
    integer::fragment_id=0,basis_generation=0,state_count=0
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64
    complex(real64),allocatable::vectors(:,:),directions(:,:)
    integer::history_rank_used=0,maximum_trial_dimension=0
  end type
  integer,parameter,public::fragment_seed=1,fragment_pw=2,fragment_projector=3
  type,public::s_dg_hybrid_fragment_candidate_catalog
    integer::fragment_id=0,basis_generation=0
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64
    ! Distributed coefficient rows in the immutable Task 5 catalog.  Seed columns
    ! are saved DC eigenspaces mapped through seed-to-WF, never named WF columns.
    complex(real64),allocatable::coefficients(:,:)
    real(real64),allocatable::energies(:) ! DC eigenvalue or PW kinetic energy; ignored for projectors.
    integer(int64),allocatable::ids(:)
    integer,allocatable::source_kind(:)
    logical,allocatable::used(:)
  end type
  type,public::s_dg_hybrid_fragment_extension_receipt
    integer::source_kind=0,old_state_count=0,new_state_count=0,candidate_count=0
    real(real64)::shell_lower=0d0,shell_upper=0d0
  end type
  abstract interface
    subroutine fragment_apply(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine
    subroutine fragment_shifted_apply(input,rayleigh_values,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      real(real64),intent(in)::rayleigh_values(:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine
  end interface
  public::advance_dg_hybrid_fragment_subspace
  public::measure_dg_hybrid_fragment_subspace
  public::extend_dg_hybrid_fragment_subspace
  public::initialize_dg_hybrid_fragment_subspace
  public::initialize_dg_hybrid_fragment_density_checked
  public::advance_dg_hybrid_fragment_epoch
contains
  ! Single-owner adapter. The legacy initializer runs on a temporary state on
  ! MPI_COMM_SELF; only a collectively density-admissible state is published.
  ! apply_s must be fragment-local (no total-communicator collectives). Guard,
  ! cutoff presence/value and tolerances are shared controls; seed counts,
  ! spectra and occupations can differ between fragments.
  subroutine initialize_dg_hybrid_fragment_density_checked(comm,fragment_id,generation,basis_fp,metric_fp,&
      seed_coefficients,seed_energies,seed_occupations,guard_count,occupation_tolerance,energy_tolerance,&
      orthogonality_tolerance,core_basis,weights,reference_density,density_tolerance,electron_tolerance,&
      apply_s,state,selected_seeds,defects,ok,message,energy_cutoff)
    integer,intent(in)::comm,fragment_id,generation,guard_count
    integer(int64),intent(in)::basis_fp,metric_fp
    complex(real64),intent(in)::seed_coefficients(:,:),core_basis(:,:)
    real(real64),intent(in)::seed_energies(:),seed_occupations(:),occupation_tolerance,energy_tolerance,&
      orthogonality_tolerance,weights(:),reference_density(:),density_tolerance,electron_tolerance
    real(real64),optional,intent(in)::energy_cutoff
    procedure(fragment_apply)::apply_s
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    integer,allocatable,intent(out)::selected_seeds(:)
    real(real64),intent(out)::defects(2)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical::halting(3)
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      type(s_dg_hybrid_fragment_subspace_state)::working
      integer,allocatable::fragments(:),chosen(:)
      integer(int64),allocatable::local_ids(:)
      complex(real64),allocatable::physical(:,:)
      real(real64),allocatable::density(:)
      integer::np,ierr,n,p,j,stat
      real(real64)::controls(8),minimum(8),maximum(8),electrons
      logical::valid,initialized
      character(256)::why
      defects=huge(1d0);ok=.false.;message='invalid density-checked initialization'
      call MPI_Comm_size(comm,np,ierr)
      n=size(seed_coefficients,1);p=size(core_basis,1)
      controls=[occupation_tolerance,energy_tolerance,orthogonality_tolerance,&
        density_tolerance,electron_tolerance,real(guard_count,real64),&
        real(merge(1,0,present(energy_cutoff)),real64),0d0]
      if(present(energy_cutoff))controls(8)=energy_cutoff
      valid=ierr==MPI_SUCCESS.and.fragment_id>=1.and.fragment_id<=np.and.n>0.and.p>0.and.&
        size(core_basis,2)==n.and.size(weights)==p.and.size(reference_density)==p.and.&
        finite(core_basis).and.all(ieee_is_finite(weights)).and.all(weights>0d0).and.&
        all(ieee_is_finite(reference_density)).and.all(reference_density>=0d0).and.&
        all(ieee_is_finite(controls)).and.density_tolerance>0d0.and.electron_tolerance>0d0
      call check(valid,'invalid density-checked initialization input');if(.not.ok)return
      call MPI_Allreduce(controls,minimum,8,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      call check(ierr==MPI_SUCCESS,'initial density controls exchange failed');if(.not.ok)return
      call MPI_Allreduce(controls,maximum,8,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call check(ierr==MPI_SUCCESS.and.all(minimum==maximum),'initial density controls differ');if(.not.ok)return
      allocate(fragments(np),local_ids(n),stat=stat)
      call check(stat==0,'initial density directory allocation failed');if(.not.ok)return
      call MPI_Allgather(fragment_id,1,MPI_INTEGER,fragments,1,MPI_INTEGER,comm,ierr)
      valid=ierr==MPI_SUCCESS
      do j=1,np;valid=valid.and.count(fragments==j)==1;enddo
      call check(valid,'density-checked initialization requires one rank per fragment');if(.not.ok)return
      local_ids=[(int(j,int64),j=1,n)]
      call initialize_dg_hybrid_fragment_subspace(MPI_COMM_SELF,n,local_ids,fragment_id,generation,&
        basis_fp,metric_fp,seed_coefficients,seed_energies,seed_occupations,guard_count,&
        occupation_tolerance,energy_tolerance,orthogonality_tolerance,apply_s,working,chosen,initialized,why,&
        energy_cutoff)
      call check(initialized,'fragment initialization failed: '//trim(why));if(.not.ok)return
      allocate(physical(p,size(chosen)),density(p),stat=stat)
      call check(stat==0,'post-initializer density allocation failed');if(.not.ok)return
      physical=matmul(core_basis,working%vectors);density=0d0
      do j=1,size(chosen)
        density=density+seed_occupations(chosen(j))*abs(physical(:,j))**2
      enddo
      electrons=sum(weights*reference_density)
      defects(1)=sum(weights*abs(density-reference_density))/max(1d0,electrons)
      defects(2)=abs(sum(weights*(density-reference_density)))
      valid=finite(physical).and.all(ieee_is_finite(density)).and.ieee_is_finite(electrons).and.&
        all(ieee_is_finite(defects))
      call check(valid,'nonfinite post-initializer density');if(.not.ok)return
      write(why,'(a,i0,a,2es14.6)')'post-initializer density mismatch: fragment=',fragment_id,&
        ' density/electron=',defects
      call check(defects(1)<=density_tolerance.and.defects(2)<=electron_tolerance,why)
      if(.not.ok)return
      call move_alloc(working%vectors,state%vectors);call move_alloc(working%directions,state%directions)
      state%fragment_id=working%fragment_id;state%basis_generation=working%basis_generation
      state%basis_fingerprint=working%basis_fingerprint;state%metric_fingerprint=working%metric_fingerprint
      state%state_count=working%state_count;state%history_rank_used=working%history_rank_used
      state%maximum_trial_dimension=working%maximum_trial_dimension
      call move_alloc(chosen,selected_seeds)
    end subroutine
    subroutine check(valid,description)
      logical,intent(in)::valid
      character(*),intent(in)::description
      integer::rank,ierr,bad,first_bad
      character(256)::why
      call MPI_Comm_rank(comm,rank,ierr)
      bad=huge(0);if(.not.valid)bad=rank
      call MPI_Allreduce(bad,first_bad,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      ok=ierr==MPI_SUCCESS;message='initial density status exchange failed'
      if(.not.ok)return
      if(first_bad==huge(0))then;message='';return;endif
      why=description
      call MPI_Bcast(why,len(why),MPI_CHARACTER,first_bad,comm,ierr)
      ok=.false.;message=why
    end subroutine
  end subroutine initialize_dg_hybrid_fragment_density_checked

  subroutine advance_dg_hybrid_fragment_epoch(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,&
      epoch,maximum_steps,apply_h,apply_s,apply_shifted_preconditioner,intermediate_tolerance,&
      orthogonality_tolerance,allowed_residual_growth,budget,state,eigenvalues,iterations,remaining_steps,&
      relative_residual,eigensolver_converged,advanced,stop_reason,workspace_peak_bytes,fingerprint,ok,message)
    ! One persistent budget per fragment. Repeated capacity/tail passes in the
    ! same density epoch share it, including after state-count extension.
    integer,intent(in)::comm,global_count,fragment_id,generation,epoch,maximum_steps
    integer(int64),intent(in)::row_ids(:),basis_fp,metric_fp
    procedure(fragment_apply)::apply_h,apply_s
    procedure(fragment_shifted_apply)::apply_shifted_preconditioner
    real(real64),intent(in)::intermediate_tolerance,orthogonality_tolerance,allowed_residual_growth
    type(s_dg_hybrid_fragment_epoch_budget),intent(inout)::budget
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    real(real64),intent(out)::eigenvalues(:),relative_residual
    integer,intent(out)::iterations,remaining_steps
    logical,intent(out)::eigensolver_converged,advanced,ok
    character(*),intent(out)::stop_reason,message
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    integer::rank,nproc,ierr
    logical::valid,halting(3)
    ok=.false.;message='invalid fragment epoch budget';stop_reason='invalid_epoch'
    eigenvalues=0d0;iterations=0;remaining_steps=0;relative_residual=huge(1d0)
    eigensolver_converged=.false.;advanced=.false.;workspace_peak_bytes=0_int64;fingerprint=0_int64
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      valid=agree_int(comm,epoch)
      valid=agree_int(comm,maximum_steps).and.valid
      valid=agree_int(comm,budget%epoch).and.valid
      valid=agree_int(comm,budget%limit).and.valid
      valid=agree_int(comm,budget%remaining).and.valid
      valid=agree_int(comm,merge(1,0,budget%failed)).and.valid
      valid=agree_real(comm,intermediate_tolerance).and.valid
      valid=agree_real(comm,orthogonality_tolerance).and.valid
      valid=agree_real(comm,allowed_residual_growth).and.valid
      if(.not.consensus(comm,valid.and.epoch>0.and.maximum_steps>=1.and.maximum_steps<=256))return
      valid=ieee_is_finite(intermediate_tolerance).and.ieee_is_finite(orthogonality_tolerance).and.&
        ieee_is_finite(allowed_residual_growth)
      if(.not.consensus(comm,valid))return
      if(intermediate_tolerance<=0d0.or.orthogonality_tolerance<64d0*epsilon(1d0).or.&
        orthogonality_tolerance>1d-2.or.allowed_residual_growth<1d0)return
      if(budget%failed)then;message='failed fragment epoch cannot be retried';return;endif
      if(epoch<budget%epoch)return
      if(epoch==budget%epoch.and.maximum_steps/=budget%limit)return
      call validate_cache(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,state,valid,message)
      if(.not.valid)return
      call MPI_Comm_rank(comm,rank,ierr);valid=ierr==MPI_SUCCESS
      call MPI_Comm_size(comm,nproc,ierr);valid=valid.and.ierr==MPI_SUCCESS
      if(budget%epoch>0)then
        valid=valid.and.budget%fragment_id==fragment_id.and.budget%generation==generation.and.&
          budget%basis_fp==basis_fp.and.budget%metric_fp==metric_fp.and.budget%global_count==global_count.and.&
          budget%comm_rank==rank.and.budget%comm_size==nproc
      endif
      if(.not.consensus(comm,valid))then;message='fragment epoch budget provenance mismatch';return;endif
      if(epoch>budget%epoch)then
        budget%epoch=epoch;budget%limit=maximum_steps;budget%remaining=maximum_steps
        budget%fragment_id=fragment_id;budget%generation=generation
        budget%basis_fp=basis_fp;budget%metric_fp=metric_fp;budget%global_count=global_count
        budget%comm_rank=rank;budget%comm_size=nproc
      endif
      if(budget%remaining==0)then
        ! Conservative local bound, including metric certification and row
        ! ownership workspace; callback-owned H/S storage is not included.
        workspace_peak_bytes=16_int64*(4_int64*size(row_ids)*state%state_count+&
          2_int64*state%state_count*state%state_count)+8_int64*state%state_count+4_int64*global_count
        call measure_dg_hybrid_fragment_subspace(comm,global_count,row_ids,fragment_id,generation,&
          basis_fp,metric_fp,apply_h,apply_s,orthogonality_tolerance,state,eigenvalues,relative_residual,&
          fingerprint,ok,message)
        stop_reason='budget_exhausted'
        eigensolver_converged=ok.and.relative_residual<=intermediate_tolerance
      else
        call advance_dg_hybrid_fragment_subspace(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,&
          apply_h,apply_s,maximum_steps=budget%remaining,intermediate_tolerance=intermediate_tolerance,&
          orthogonality_tolerance=orthogonality_tolerance,allowed_residual_growth=allowed_residual_growth,&
          state=state,eigenvalues=eigenvalues,iterations=iterations,relative_residual=relative_residual,&
          eigensolver_converged=eigensolver_converged,advanced=advanced,stop_reason=stop_reason,&
          workspace_peak_bytes=workspace_peak_bytes,fingerprint=fingerprint,ok=ok,message=message,&
          apply_shifted_preconditioner=apply_shifted_preconditioner)
        ! Count attempted updates even when rollback retains the entry state.
        budget%remaining=budget%remaining-iterations
      endif
      remaining_steps=budget%remaining
      if(.not.ok)budget%failed=.true. ! No implicit recovery after a callback/metric failure.
    end subroutine
  end subroutine advance_dg_hybrid_fragment_epoch

  subroutine initialize_dg_hybrid_fragment_subspace(comm,global_count,row_ids,fragment_id,generation,&
      basis_fp,metric_fp,seed_coefficients,seed_energies,seed_occupations,guard_count,&
      occupation_tolerance,energy_tolerance,orthogonality_tolerance,apply_s,state,selected_seeds,ok,message,&
      energy_cutoff)
    ! Seed columns are physical DC orbitals in uncompressed WF+PW coordinates;
    ! selected-space projection may give nonzero PW components. Selection never
    ! means taking the first named WFs.
    ! The returned inventory is only an initial guess, not a tail/symmetry certificate.
    integer,intent(in)::comm,global_count,fragment_id,generation,guard_count
    integer(int64),intent(in)::row_ids(:),basis_fp,metric_fp
    complex(real64),intent(in)::seed_coefficients(:,:)
    real(real64),intent(in)::seed_energies(:),seed_occupations(:),occupation_tolerance,&
      energy_tolerance,orthogonality_tolerance
    real(real64),optional,intent(in)::energy_cutoff
    procedure(fragment_apply)::apply_s
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    integer,allocatable,intent(out)::selected_seeds(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_fragment_subspace_state)::working
    integer,allocatable::order(:),chosen(:)
    integer::nseed,nselected,occupied,j,k,tmp,stat
    real(real64)::edge,scale
    logical::valid,halting(3)
    ok=.false.;message='invalid fragment seed initialization contract'
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      nseed=size(seed_energies)
      valid=agree_int(comm,nseed)
      valid=agree_int(comm,global_count).and.valid
      valid=agree_int(comm,guard_count).and.valid
      valid=agree_real(comm,occupation_tolerance).and.valid
      valid=agree_real(comm,energy_tolerance).and.valid
      valid=agree_real(comm,orthogonality_tolerance).and.valid
      valid=agree_int(comm,merge(1,0,present(energy_cutoff))).and.valid
      if(.not.consensus(comm,valid))return
      if(present(energy_cutoff))then
        valid=agree_real(comm,energy_cutoff)
        if(.not.consensus(comm,valid.and.ieee_is_finite(energy_cutoff)))return
      endif
      valid=nseed>0.and.guard_count>=0.and.all(shape(seed_coefficients)==[size(row_ids),nseed]).and.&
        size(seed_occupations)==nseed
      valid=valid.and.ieee_is_finite(occupation_tolerance).and.ieee_is_finite(energy_tolerance).and.&
        ieee_is_finite(orthogonality_tolerance)
      if(.not.consensus(comm,valid))return
      if(occupation_tolerance<0d0.or.energy_tolerance<=0d0.or.energy_tolerance>1d-2.or.&
        orthogonality_tolerance<64d0*epsilon(1d0).or.orthogonality_tolerance>1d-2)return
      valid=finite(seed_coefficients).and.all(ieee_is_finite(seed_energies)).and.&
        all(ieee_is_finite(seed_occupations)).and.all(seed_occupations>=0d0)
      if(.not.consensus(comm,valid))return
      do j=1,nseed
        valid=agree_real(comm,seed_energies(j))
        valid=agree_real(comm,seed_occupations(j)).and.valid
        if(.not.valid)return
      enddo
      allocate(order(nseed),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate seed selection';return;endif
      order=[(j,j=1,nseed)]
      do j=2,nseed
        tmp=order(j);k=j
        do while(k>1)
          if(seed_energies(order(k-1))<=seed_energies(tmp))exit
          order(k)=order(k-1);k=k-1
        enddo
        order(k)=tmp
      enddo
      occupied=0
      do j=1,nseed
        if(seed_occupations(order(j))>occupation_tolerance)occupied=j
      enddo
      ! Saturate without integer overflow; later capacity/tail checks request
      ! more invariant shells if the saved DC inventory is insufficient.
      nselected=occupied+min(guard_count,nseed-occupied)
      if(present(energy_cutoff))then
        do j=1,nseed
          if(seed_energies(order(j))<=energy_cutoff)nselected=max(nselected,j)
        enddo
      endif
      if(nselected==0)then;message='initial DC occupied/guard inventory is empty';return;endif
      edge=seed_energies(order(nselected))
      do j=nselected+1,nseed
        scale=max(1d0,abs(edge),abs(seed_energies(order(j))))
        if(abs(seed_energies(order(j))/scale-edge/scale)>energy_tolerance)exit
        nselected=j
      enddo
      if(nselected>=global_count)then
        message='initial seed selection exhausts fragment basis; enlarge basis before initialization';return
      endif
      valid=int(nselected,int64)**2<=int(huge(0),int64)/100_int64.and.&
        int(size(row_ids),int64)*int(nselected,int64)<=int(huge(0),int64)/40_int64
      if(.not.consensus(comm,valid))then;message='initial fragment state exceeds safe workspace bounds';return;endif
      working%fragment_id=fragment_id;working%basis_generation=generation
      working%basis_fingerprint=basis_fp;working%metric_fingerprint=metric_fp
      working%state_count=nselected
      allocate(working%vectors(size(row_ids),nselected),working%directions(size(row_ids),nselected),&
        chosen(nselected),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate initial fragment state';return;endif
      chosen=order(:nselected);working%vectors=seed_coefficients(:,chosen);working%directions=0d0
      call validate_cache(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,working,valid,message)
      if(.not.valid)return
      call normalize_entry(comm,apply_s,working%vectors,orthogonality_tolerance,valid)
      if(.not.valid)then;message='physical DC seed subspace has invalid metric rank';return;endif
      call move_alloc(working%vectors,state%vectors);call move_alloc(working%directions,state%directions)
      state%fragment_id=fragment_id;state%basis_generation=generation
      state%basis_fingerprint=basis_fp;state%metric_fingerprint=metric_fp
      state%state_count=nselected;state%history_rank_used=0;state%maximum_trial_dimension=0
      call move_alloc(chosen,selected_seeds)
      ok=.true.;message=''
    end subroutine
  end subroutine initialize_dg_hybrid_fragment_subspace

  subroutine measure_dg_hybrid_fragment_subspace(comm,global_count,row_ids,fragment_id,generation,&
      basis_fp,metric_fp,apply_h,apply_s,orthogonality_tolerance,state,eigenvalues,relative_residual,&
      fingerprint,ok,message)
    ! Used after shell extension when this density epoch has no CG budget left.
    ! No normalization, Ritz rotation, preconditioning, or history update occurs.
    ! Values retain coefficient-column order; occupation packing must apply the
    ! same spectral permutation to energies, core norms and returned occupations.
    integer,intent(in)::comm,global_count,fragment_id,generation
    integer(int64),intent(in)::row_ids(:),basis_fp,metric_fp
    procedure(fragment_apply)::apply_h,apply_s
    real(real64),intent(in)::orthogonality_tolerance
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    real(real64),intent(out)::eigenvalues(:),relative_residual
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::hx(:,:),sx(:,:),residual(:,:)
    real(real64),allocatable::values(:)
    real(real64)::measured
    logical::halting(3),valid
    integer::nr,ns,stat,j
    ok=.false.;message='';eigenvalues=0d0;relative_residual=huge(1d0);fingerprint=0_int64
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      call validate_cache(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,state,valid,message)
      if(.not.valid)return
      nr=size(row_ids);ns=state%state_count
      valid=agree_real(comm,orthogonality_tolerance)
      valid=consensus(comm,valid.and.size(eigenvalues)==ns.and.ieee_is_finite(orthogonality_tolerance))
      if(.not.valid)then;message='invalid or rank-disagreeing measurement controls';return;endif
      if(orthogonality_tolerance<64d0*epsilon(1d0).or.orthogonality_tolerance>1d-2)then
        message='invalid measurement orthogonality tolerance';return
      endif
      call certify_metric(comm,apply_s,state%vectors,orthogonality_tolerance,valid)
      if(.not.valid)then;message='measurement requires a safe metric-orthonormal state';return;endif
      allocate(hx(nr,ns),sx(nr,ns),residual(nr,ns),values(ns),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate measurement workspace';return;endif
      call evaluate(comm,apply_h,apply_s,state%vectors,hx,sx,values,residual,measured,valid)
      if(.not.valid)then;message='fragment spectrum measurement failed';return;endif
      ! Intermediate diagnostic only, not an operator-epoch or final eigenpair certificate.
      fingerprint=ieor(basis_fp,ishftc(metric_fp,17))
      do j=1,ns;fingerprint=ieor(ishftc(fingerprint,7),transfer(values(j),0_int64));enddo
      if(fingerprint==0_int64)fingerprint=1_int64
      eigenvalues=values;relative_residual=measured;ok=.true.;message=''
    end subroutine
  end subroutine measure_dg_hybrid_fragment_subspace

  subroutine extend_dg_hybrid_fragment_subspace(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,&
      apply_h,apply_s,energy_tolerance,orthogonality_tolerance,catalog,state,receipt,ok,message)
    integer,intent(in)::comm,global_count,fragment_id,generation
    integer(int64),intent(in)::row_ids(:),basis_fp,metric_fp
    procedure(fragment_apply)::apply_h,apply_s
    real(real64),intent(in)::energy_tolerance,orthogonality_tolerance
    type(s_dg_hybrid_fragment_candidate_catalog),intent(inout)::catalog
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    type(s_dg_hybrid_fragment_extension_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical::halting(3),valid
    integer::np,nold,nr,stat,j,k,source,nselected,dimension,added,position,temp
    integer,allocatable::selected(:)
    logical,allocatable::used(:),mask(:)
    complex(real64),allocatable::basis(:,:),pool(:,:),newx(:,:),newp(:,:),sx(:,:),metric(:,:),rotated(:,:)
    real(real64),allocatable::values(:)
    real(real64)::shell,scale
    ok=.false.;message='invalid fragment extension catalog'
    receipt=s_dg_hybrid_fragment_extension_receipt()
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      call validate_cache(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,state,valid,message)
      if(.not.valid)return
      message='invalid fragment extension catalog'
      valid=agree_real(comm,energy_tolerance)
      valid=agree_real(comm,orthogonality_tolerance).and.valid
      valid=valid.and.ieee_is_finite(energy_tolerance).and.ieee_is_finite(orthogonality_tolerance)
      if(.not.consensus(comm,valid))return
      if(energy_tolerance<=0d0.or.energy_tolerance>1d-2.or.orthogonality_tolerance<64d0*epsilon(1d0).or.&
        orthogonality_tolerance>1d-2)return
      valid=catalog%fragment_id==fragment_id.and.catalog%basis_generation==generation.and.&
        catalog%basis_fingerprint==basis_fp.and.catalog%metric_fingerprint==metric_fp
      valid=valid.and.allocated(catalog%coefficients).and.allocated(catalog%energies).and.&
        allocated(catalog%ids).and.allocated(catalog%source_kind).and.allocated(catalog%used)
      if(.not.consensus(comm,valid))return
      np=size(catalog%ids);nr=size(row_ids);nold=state%state_count
      valid=agree_int(comm,np)
      valid=valid.and.size(catalog%coefficients,1)==nr.and.size(catalog%coefficients,2)==np.and.&
        size(catalog%energies)==np.and.size(catalog%source_kind)==np.and.size(catalog%used)==np
      if(.not.consensus(comm,valid))return
      if(int(np,int64)**2>int(huge(0),int64)/100_int64)return
      valid=finite(catalog%coefficients).and.all(ieee_is_finite(catalog%energies)).and.&
        all(catalog%ids>0_int64).and.all(catalog%source_kind>=fragment_seed).and.&
        all(catalog%source_kind<=fragment_projector)
      if(.not.consensus(comm,valid))return
      do j=1,np
        valid=agree_bits(comm,catalog%ids(j))
        valid=agree_real(comm,catalog%energies(j)).and.valid
        valid=agree_int(comm,catalog%source_kind(j)).and.valid
        valid=agree_int(comm,merge(1,0,catalog%used(j))).and.valid
        if(count(catalog%ids==catalog%ids(j))/=1)valid=.false.
        if(.not.consensus(comm,valid))return
      enddo
      allocate(selected(np),used(np),mask(np),sx(nr,nold),metric(nold,nold),&
        basis(nr,min(global_count,nold+np)),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate extension workspace';return;endif
      call checked_apply(comm,apply_s,state%vectors,sx,valid);if(.not.valid)return
      call gram_matrix(comm,state%vectors,sx,metric,valid);if(.not.valid)return
      do j=1,nold;metric(j,j)=metric(j,j)-1d0;enddo
      if(maxval(abs(metric))>orthogonality_tolerance)then
        message='extension requires a safe metric-orthonormal entry';return
      endif
      used=catalog%used
      do
        mask=.not.used.and.catalog%source_kind==fragment_seed;source=fragment_seed
        if(.not.any(mask))then
          mask=.not.used.and.catalog%source_kind==fragment_pw;source=fragment_pw
        endif
        if(any(mask))then
          shell=minval(catalog%energies,mask=mask)
          scale=max(1d0,abs(shell),maxval(abs(catalog%energies),mask=mask))
          mask=mask.and.abs(catalog%energies/scale-shell/scale)<=energy_tolerance
          receipt%shell_lower=minval(catalog%energies,mask=mask)
          receipt%shell_upper=maxval(catalog%energies,mask=mask)
          if(source==fragment_pw)mask=mask.or.(.not.used.and.catalog%source_kind==fragment_projector)
        else
          source=fragment_projector;mask=.not.used.and.catalog%source_kind==fragment_projector
          receipt%shell_lower=0d0;receipt%shell_upper=0d0
        endif
        nselected=count(mask)
        if(nselected==0.or.nold==global_count)then
          message='insufficient-spectrum: fragment candidate metric rank exhausted';return
        endif
        k=0
        do j=1,np
          if(.not.mask(j))cycle
          k=k+1;selected(k)=j
        enddo
        ! Stable IDs order the complete pool; no ID removes a degenerate direction.
        do j=2,nselected
          temp=selected(j);position=j
          do while(position>1)
            if(catalog%ids(selected(position-1))<catalog%ids(temp))exit
            selected(position)=selected(position-1);position=position-1
          enddo
          selected(position)=temp
        enddo
        if(allocated(pool))deallocate(pool)
        allocate(pool(nr,nselected),stat=stat)
        if(.not.consensus(comm,stat==0))then;message='cannot allocate next invariant shell';return;endif
        pool=catalog%coefficients(:,selected(:nselected));basis(:,:nold)=state%vectors;dimension=nold
        call append_metric_pool(comm,apply_s,basis,dimension,pool,orthogonality_tolerance,added,valid)
        if(.not.valid)then;message='invalid candidate shell metric';return;endif
        used=used.or.mask
        if(added==0)cycle ! A previously represented whole pool is consumed, never partially selected.
        exit
      enddo
      allocate(rotated(nr,added),values(added),newx(nr,dimension),newp(nr,dimension),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate expanded fragment cache';return;endif
      ! Only the new invariant pool is diagonalized; the accepted X is not rotated.
      call ritz(comm,apply_h,apply_s,basis(:,nold+1:dimension),added,rotated,values,valid)
      if(.not.valid)then;message='candidate pool Rayleigh-Ritz failed';return;endif
      newx(:,:nold)=state%vectors;newx(:,nold+1:)=rotated
      newp=0d0;newp(:,:nold)=state%directions
      call certify_metric(comm,apply_s,newx,orthogonality_tolerance,valid)
      if(.not.valid)then;message='expanded fragment metric certificate failed';return;endif
      ! Transactional publication after every rank has certified the new whole pool.
      call move_alloc(newx,state%vectors);call move_alloc(newp,state%directions)
      state%state_count=dimension;catalog%used=used
      receipt%source_kind=source;receipt%old_state_count=nold
      receipt%new_state_count=dimension;receipt%candidate_count=nselected
      ok=.true.;message=''
    end subroutine
  end subroutine extend_dg_hybrid_fragment_subspace

  subroutine advance_dg_hybrid_fragment_subspace(comm,global_count,row_ids,&
      fragment_id,basis_generation,basis_fingerprint,metric_fingerprint,&
      apply_h,apply_s,apply_preconditioner,maximum_steps,intermediate_tolerance,&
      orthogonality_tolerance,allowed_residual_growth,state,eigenvalues,&
      iterations,relative_residual,eigensolver_converged,advanced,&
      stop_reason,workspace_peak_bytes,fingerprint,ok,message,apply_shifted_preconditioner)
    ! Exactly one explicit preconditioner is mandatory. To preserve construction-WF gauge
    ! covariance, transform a physical preconditioner with the basis; a fresh
    ! coordinate-diagonal H approximation is not generally covariant. An
    ! identity callback is used only explicitly by reference fixtures.
    ! This routine never decides outer density-SCF convergence.
    integer,intent(in)::comm,global_count,fragment_id,basis_generation,maximum_steps
    integer(int64),intent(in)::row_ids(:),basis_fingerprint,metric_fingerprint
    procedure(fragment_apply)::apply_h,apply_s
    procedure(fragment_apply),optional::apply_preconditioner
    procedure(fragment_shifted_apply),optional::apply_shifted_preconditioner
    real(real64),intent(in)::intermediate_tolerance,orthogonality_tolerance,allowed_residual_growth
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    real(real64),intent(out)::eigenvalues(:),relative_residual
    integer,intent(out)::iterations
    logical,intent(out)::eigensolver_converged,advanced,ok
    character(*),intent(out)::stop_reason,message
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    complex(real64),allocatable::x(:,:),p(:,:),hx(:,:),sx(:,:),r(:,:),z(:,:),trial(:,:),&
      y(:,:),hy(:,:),sy(:,:),newp(:,:),bestx(:,:),bestp(:,:),overlap(:,:)
    real(real64),allocatable::values(:),bestvalues(:),candidate_values(:)
    logical::halting(3),valid
    real(real64)::entry_residual,best_residual,candidate_residual,history_norm
    integer::nr,ns,dim,added,step,stat,history_used,trial_peak,j,preconditioner_kind

    ok=.false.;message='';stop_reason='invalid_contract';iterations=0
    relative_residual=huge(1d0);eigensolver_converged=.false.;advanced=.false.
    eigenvalues=0d0;workspace_peak_bytes=0_int64;fingerprint=0_int64
    call suspend_traps(halting)
    call execute()
    call restore_traps(halting)
  contains
    subroutine execute()
      preconditioner_kind=merge(1,0,present(apply_preconditioner))+merge(2,0,present(apply_shifted_preconditioner))
      valid=agree_int(comm,preconditioner_kind)
      if(.not.valid.or.(preconditioner_kind/=1.and.preconditioner_kind/=2))then
        message='exactly one rank-consistent preconditioner callback is required';return
      endif
      call validate_cache(comm,global_count,row_ids,fragment_id,basis_generation,&
        basis_fingerprint,metric_fingerprint,state,valid,message)
      if(.not.valid)return
      nr=size(row_ids);ns=state%state_count
      valid=agree_int(comm,maximum_steps)
      valid=agree_real(comm,intermediate_tolerance).and.valid
      valid=agree_real(comm,orthogonality_tolerance).and.valid
      valid=agree_real(comm,allowed_residual_growth).and.valid
      valid=consensus(comm,valid.and.maximum_steps>0.and.maximum_steps<=256.and.&
        size(eigenvalues)==ns.and.ieee_is_finite(intermediate_tolerance).and.&
        ieee_is_finite(orthogonality_tolerance).and.ieee_is_finite(allowed_residual_growth))
      if(.not.valid)then;message='invalid or rank-disagreeing bounded update controls';return;endif
      if(intermediate_tolerance<=0d0.or.orthogonality_tolerance<64d0*epsilon(1d0).or.&
        orthogonality_tolerance>1d-2.or.allowed_residual_growth<1d0)then
        message='invalid bounded update tolerance or growth bound';return
      endif
      allocate(x(nr,ns),p(nr,ns),hx(nr,ns),sx(nr,ns),r(nr,ns),z(nr,ns),trial(nr,min(global_count,3*ns)),&
        y(nr,ns),hy(nr,ns),sy(nr,ns),newp(nr,ns),bestx(nr,ns),bestp(nr,ns),overlap(ns,ns),&
        values(ns),bestvalues(ns),candidate_values(ns),stat=stat)
      if(.not.consensus(comm,stat==0))then;message='cannot allocate bounded update workspace';return;endif
      ! Conservative upper bound, including nested metric/Ritz work and collective temporaries.
      workspace_peak_bytes=16_int64*(40_int64*nr*ns+100_int64*ns*ns)+&
        8_int64*32_int64*ns+4_int64*global_count
      x=state%vectors;p=state%directions
      call normalize_entry(comm,apply_s,x,orthogonality_tolerance,valid)
      if(.not.valid)then;message='no safe entry: nonfinite, indefinite or rank deficient metric';return;endif
      call evaluate(comm,apply_h,apply_s,x,hx,sx,values,r,entry_residual,valid)
      if(.not.valid)then;message='entry operator or residual failed';return;endif
      bestx=x;bestp=p;bestvalues=values;best_residual=entry_residual
      history_used=0;trial_peak=ns;stop_reason='step_cap'
      if(best_residual<=intermediate_tolerance)then
        stop_reason='intermediate_target'
      else
        do step=1,maximum_steps
          iterations=step
          if(present(apply_shifted_preconditioner))then
            z=0d0
            call apply_shifted_preconditioner(r,values,z,valid)
            valid=consensus(comm,valid.and.finite(z))
          else if(present(apply_preconditioner))then
            call checked_apply(comm,apply_preconditioner,r,z,valid)
          endif
          if(.not.valid)then;message='fragment preconditioner callback failed';return;endif
          ! P is retained across potential epochs.  It is projected afresh in the current S metric.
          call matrix_norm(comm,p,history_norm,valid)
          if(.not.valid)then;message='nonfinite direction history';return;endif
          if(history_norm>orthogonality_tolerance)history_used=ns
          trial(:,1:ns)=x;dim=ns
          call append_metric_pool(comm,apply_s,trial,dim,z,orthogonality_tolerance,added,valid)
          if(.not.valid)then;message='residual trial metric failed';return;endif
          ! A cancelled conjugate direction is roundoff, not a new normalized
          ! search direction.  R has its own scale-independent residual contract.
          if(history_norm>orthogonality_tolerance)then
            call append_metric_pool(comm,apply_s,trial,dim,p,orthogonality_tolerance,added,valid)
            if(.not.valid)then;message='history trial metric failed';return;endif
          endif
          trial_peak=max(trial_peak,dim)
          ! Even an invariant X may need an internal Ritz rotation.  This is
          ! still an nstate-sized solve, never a full fragment cold start.
          call ritz(comm,apply_h,apply_s,trial(:,:dim),ns,y,candidate_values,valid)
          if(.not.valid)then;message='fragment trial Rayleigh-Ritz failed';return;endif
          call certify_metric(comm,apply_s,y,orthogonality_tolerance,valid)
          if(.not.valid)then;message='fragment trial metric certificate failed';return;endif
          call evaluate(comm,apply_h,apply_s,y,hy,sy,candidate_values,z,candidate_residual,valid)
          if(.not.valid)then;message='nonfinite trial operator or residual';return;endif
          if(candidate_residual/allowed_residual_growth>entry_residual)cycle
          call gram_matrix(comm,x,sy,overlap,valid)
          if(.not.valid)then;message='conjugate direction projection failed';return;endif
          newp=y-matmul(x,overlap)
          x=y;p=newp;hx=hy;sx=sy;values=candidate_values;r=z
          if(candidate_residual<best_residual-64d0*epsilon(1d0)*max(1d0,best_residual).or.&
            (candidate_residual<=intermediate_tolerance.and.candidate_residual<best_residual))then
            bestx=x;bestp=p;bestvalues=values;best_residual=candidate_residual;advanced=.true.
          endif
          if(best_residual<=intermediate_tolerance)then;stop_reason='intermediate_target';exit;endif
        enddo
        if(.not.advanced)stop_reason='safe_entry_retained'
      endif
      state%vectors=bestx;state%directions=bestp
      state%history_rank_used=history_used;state%maximum_trial_dimension=trial_peak
      eigenvalues=bestvalues;relative_residual=best_residual
      eigensolver_converged=relative_residual<=intermediate_tolerance
      ! This receipt is an intermediate density-update diagnostic, not a final LCFO certificate.
      fingerprint=ieor(basis_fingerprint,ishftc(metric_fingerprint,17))
      do j=1,ns;fingerprint=ieor(ishftc(fingerprint,7),transfer(eigenvalues(j),0_int64));enddo
      if(fingerprint==0_int64)fingerprint=1_int64
      ok=.true.;message=''
    end subroutine
  end subroutine advance_dg_hybrid_fragment_subspace

  subroutine validate_cache(comm,global_count,row_ids,fragment_id,generation,basis_fp,metric_fp,state,ok,message)
    integer,intent(in)::comm,global_count,fragment_id,generation
    integer(int64),intent(in)::row_ids(:),basis_fp,metric_fp
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::owners(:)
    integer::i,stat,ierr,nstate
    logical::valid,stale
    ok=.false.;message='invalid fragment subspace cache'
    valid=agree_int(comm,global_count)
    valid=agree_int(comm,fragment_id).and.valid
    valid=agree_int(comm,generation).and.valid
    valid=agree_bits(comm,basis_fp).and.valid
    valid=agree_bits(comm,metric_fp).and.valid
    nstate=state%state_count;valid=agree_int(comm,nstate).and.valid
    if(.not.consensus(comm,valid.and.global_count>0.and.fragment_id>0.and.generation>0.and.&
      basis_fp/=0_int64.and.metric_fp/=0_int64.and.nstate>0.and.nstate<=global_count))return
    ! Guard all MPI matrix counts and the conservative byte receipt before allocating.
    if(int(nstate,int64)**2>int(huge(0),int64)/100_int64.or.&
      int(size(row_ids),int64)*int(nstate,int64)>int(huge(0),int64)/40_int64)valid=.false.
    if(.not.consensus(comm,valid))return
    stale=state%fragment_id/=fragment_id.or.state%basis_generation/=generation.or.&
      state%basis_fingerprint/=basis_fp.or.state%metric_fingerprint/=metric_fp
    if(.not.consensus(comm,.not.stale))then
      if(allocated(state%vectors))deallocate(state%vectors)
      if(allocated(state%directions))deallocate(state%directions)
      state%state_count=0;state%history_rank_used=0;state%maximum_trial_dimension=0
      message='fragment cache invalidated: basis generation or metric identity changed';return
    endif
    valid=allocated(state%vectors).and.allocated(state%directions)
    if(.not.consensus(comm,valid))return
    valid=size(state%vectors,1)==size(row_ids).and.size(state%directions,1)==size(row_ids).and.&
      size(state%vectors,2)==nstate.and.size(state%directions,2)==nstate
    valid=valid.and.finite(state%vectors).and.finite(state%directions)
    valid=valid.and.all(row_ids>=1_int64).and.all(row_ids<=int(global_count,int64))
    if(.not.consensus(comm,valid))return
    allocate(owners(global_count),stat=stat)
    if(.not.consensus(comm,stat==0))then;message='cannot allocate fragment row ownership';return;endif
    owners=0
    do i=1,size(row_ids);owners(int(row_ids(i)))=owners(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(owners/=1))then;message='fragment rows are not owned exactly once';return;endif
    ok=.true.;message=''
  end subroutine

  subroutine certify_metric(comm,apply_s,x,tolerance,ok)
    integer,intent(in)::comm
    procedure(fragment_apply)::apply_s
    complex(real64),intent(in)::x(:,:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    complex(real64),allocatable::sx(:,:),metric(:,:)
    integer::j,stat,ns
    ns=size(x,2)
    allocate(sx(size(x,1),ns),metric(ns,ns),stat=stat)
    ok=consensus(comm,stat==0);if(.not.ok)return
    call checked_apply(comm,apply_s,x,sx,ok);if(.not.ok)return
    call gram_matrix(comm,x,sx,metric,ok);if(.not.ok)return
    do j=1,ns;metric(j,j)=metric(j,j)-1d0;enddo
    ok=consensus(comm,maxval(abs(metric))<=tolerance)
  end subroutine

  subroutine normalize_entry(comm,apply_s,x,tolerance,ok)
    integer,intent(in)::comm
    procedure(fragment_apply)::apply_s
    complex(real64),intent(inout)::x(:,:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    complex(real64),allocatable::sx(:,:),gram(:,:),transform(:,:)
    real(real64),allocatable::e(:)
    integer::ns,j,stat
    real(real64)::error
    ns=size(x,2)
    allocate(sx(size(x,1),ns),gram(ns,ns),transform(ns,ns),e(ns),stat=stat)
    ok=consensus(comm,stat==0);if(.not.ok)return
    call checked_apply(comm,apply_s,x,sx,ok);if(.not.ok)return
    call gram_matrix(comm,x,sx,gram,ok);if(.not.ok)return
    transform=gram
    do j=1,ns;transform(j,j)=transform(j,j)-1d0;enddo
    error=maxval(abs(transform))
    if(error<=tolerance)return ! Preserve accepted coefficients bitwise on a safe rollback.
    transform=gram;call hermitian_eigen(comm,transform,e,ok);if(.not.ok)return
    ok=minval(e)>tolerance*max(1d0,maxval(abs(e)));if(.not.ok)return
    gram=transform
    do j=1,ns;gram(:,j)=gram(:,j)/sqrt(e(j));enddo
    x=matmul(x,matmul(gram,conjg(transpose(transform))))
    call checked_apply(comm,apply_s,x,sx,ok);if(.not.ok)return
    call gram_matrix(comm,x,sx,gram,ok);if(.not.ok)return
    do j=1,ns;gram(j,j)=gram(j,j)-1d0;enddo
    ok=consensus(comm,finite(x).and.maxval(abs(gram))<=tolerance)
  end subroutine

  subroutine append_metric_pool(comm,apply_s,basis,dimension,pool,tolerance,added,ok)
    integer,intent(in)::comm
    procedure(fragment_apply)::apply_s
    complex(real64),intent(inout)::basis(:,:)
    integer,intent(inout)::dimension
    complex(real64),intent(in)::pool(:,:)
    real(real64),intent(in)::tolerance
    integer,intent(out)::added
    logical,intent(out)::ok
    complex(real64),allocatable::q(:,:),sq(:,:),overlap(:,:),metric(:,:)
    real(real64),allocatable::e(:)
    integer::np,j,pass,stat
    real(real64)::cutoff,scale,pool_norm
    np=size(pool,2);added=0;ok=.true.
    if(np==0.or.dimension==size(basis,2))return
    allocate(q(size(pool,1),np),sq(size(pool,1),np),overlap(dimension,np),metric(np,np),e(np),stat=stat)
    ok=consensus(comm,stat==0);if(.not.ok)return
    ! A uniform preconditioner scaling must not change the trial space. Scale the
    ! entire pool (not named columns) before projecting/rank revelation.
    call matrix_norm(comm,pool,pool_norm,ok);if(.not.ok)return
    if(pool_norm==0d0)return
    q=pool/pool_norm
    do pass=1,2
      call checked_apply(comm,apply_s,q,sq,ok);if(.not.ok)return
      call gram_matrix(comm,basis(:,:dimension),sq,overlap,ok);if(.not.ok)return
      q=q-matmul(basis(:,:dimension),overlap)
    enddo
    call checked_apply(comm,apply_s,q,sq,ok);if(.not.ok)return
    call gram_matrix(comm,q,sq,metric,ok);if(.not.ok)return
    call hermitian_eigen(comm,metric,e,ok);if(.not.ok)return
    scale=maxval(abs(e))
    cutoff=max(tolerance**2,64d0*epsilon(1d0)**2*np,64d0*epsilon(1d0)*scale*np)
    ok=minval(e)>=-cutoff;if(.not.ok)return
    added=count(e>cutoff)
    ! Dependence is expected in residual/history pools. Never truncate a non-null pool by column name.
    ok=added<=size(basis,2)-dimension;if(.not.ok)return
    do j=1,np
      if(e(j)<=cutoff)cycle
      dimension=dimension+1;basis(:,dimension)=matmul(q,metric(:,j))/sqrt(e(j))
    enddo
    ok=consensus(comm,finite(basis(:,:dimension)))
  end subroutine

  subroutine ritz(comm,apply_h,apply_s,basis,nstate,x,eigenvalues,ok)
    integer,intent(in)::comm,nstate
    procedure(fragment_apply)::apply_h,apply_s
    complex(real64),intent(in)::basis(:,:)
    complex(real64),intent(out)::x(:,:)
    real(real64),intent(out)::eigenvalues(:)
    logical,intent(out)::ok
    complex(real64),allocatable::h(:,:),s(:,:),hb(:,:),sb(:,:),work(:)
    real(real64),allocatable::e(:),rwork(:)
    integer::dim,stat,info
    external::zhegv
    dim=size(basis,2)
    allocate(h(dim,dim),s(dim,dim),hb(size(basis,1),dim),sb(size(basis,1),dim),&
      work(max(1,2*dim)),e(dim),rwork(max(1,3*dim-2)),stat=stat)
    ok=consensus(comm,stat==0);if(.not.ok)return
    call checked_apply(comm,apply_h,basis,hb,ok);if(.not.ok)return
    call checked_apply(comm,apply_s,basis,sb,ok);if(.not.ok)return
    call gram_matrix(comm,basis,hb,h,ok);if(.not.ok)return
    call gram_matrix(comm,basis,sb,s,ok);if(.not.ok)return
    ok=consensus(comm,hermitian(h).and.hermitian(s));if(.not.ok)return
    call zhegv(1,'V','U',dim,h,dim,s,dim,e,work,size(work),rwork,info)
    ok=consensus(comm,info==0.and.finite(h).and.all(ieee_is_finite(e)));if(.not.ok)return
    x=matmul(basis,h(:,1:nstate));eigenvalues=e(1:nstate)
    ok=consensus(comm,finite(x))
  end subroutine

  subroutine evaluate(comm,apply_h,apply_s,x,hx,sx,values,residual,receipt,ok)
    integer,intent(in)::comm
    procedure(fragment_apply)::apply_h,apply_s
    complex(real64),intent(in)::x(:,:)
    complex(real64),intent(out)::hx(:,:),sx(:,:),residual(:,:)
    real(real64),intent(out)::values(:),receipt
    logical,intent(out)::ok
    complex(real64)::expectation
    real(real64)::hnorm,rnorm
    integer::j,ierr
    call checked_apply(comm,apply_h,x,hx,ok);if(.not.ok)return
    call checked_apply(comm,apply_s,x,sx,ok);if(.not.ok)return
    receipt=0d0
    do j=1,size(x,2)
      expectation=sum(conjg(x(:,j))*hx(:,j))
      call MPI_Allreduce(MPI_IN_PLACE,expectation,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      ok=consensus(comm,ierr==MPI_SUCCESS.and.ieee_is_finite(real(expectation,real64)).and.&
        ieee_is_finite(aimag(expectation)));if(.not.ok)return
      values(j)=real(expectation,real64)
      residual(:,j)=hx(:,j)-values(j)*sx(:,j)
      call matrix_norm(comm,hx(:,j:j),hnorm,ok);if(.not.ok)return
      call matrix_norm(comm,residual(:,j:j),rnorm,ok);if(.not.ok)return
      receipt=max(receipt,rnorm/max(1d0,abs(values(j)),hnorm))
    enddo
    ok=consensus(comm,ieee_is_finite(receipt))
  end subroutine

  subroutine matrix_norm(comm,x,norm,ok)
    integer,intent(in)::comm
    complex(real64),intent(in)::x(:,:)
    real(real64),intent(out)::norm
    logical,intent(out)::ok
    real(real64)::scale,sum_square
    integer::ierr
    ok=consensus(comm,finite(x));if(.not.ok)return
    scale=0d0;if(size(x)>0)scale=maxval(abs(x))
    call MPI_Allreduce(MPI_IN_PLACE,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS;if(.not.ok)return
    norm=0d0;if(scale==0d0)return
    sum_square=sum(abs(x/scale)**2)
    call MPI_Allreduce(MPI_IN_PLACE,sum_square,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    norm=scale*sqrt(sum_square);ok=consensus(comm,ierr==MPI_SUCCESS.and.ieee_is_finite(norm))
  end subroutine

  subroutine gram_matrix(comm,x,y,gram,ok)
    integer,intent(in)::comm
    complex(real64),intent(in)::x(:,:),y(:,:)
    complex(real64),intent(out)::gram(:,:)
    logical,intent(out)::ok
    integer::ierr
    gram=matmul(conjg(transpose(x)),y)
    call MPI_Allreduce(MPI_IN_PLACE,gram,size(gram),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    ok=consensus(comm,ierr==MPI_SUCCESS.and.finite(gram))
  end subroutine

  subroutine checked_apply(comm,apply,input,output,ok)
    integer,intent(in)::comm
    procedure(fragment_apply)::apply
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::ok
    logical::valid
    ok=consensus(comm,finite(input));if(.not.ok)return
    output=0d0;call apply(input,output,valid)
    ok=consensus(comm,valid.and.finite(output))
  end subroutine

  subroutine hermitian_eigen(comm,matrix,values,ok)
    integer,intent(in)::comm
    complex(real64),intent(inout)::matrix(:,:)
    real(real64),intent(out)::values(:)
    logical,intent(out)::ok
    complex(real64),allocatable::work(:)
    real(real64),allocatable::rwork(:)
    integer::n,info,stat
    external::zheev
    n=size(matrix,1)
    ok=consensus(comm,hermitian(matrix));if(.not.ok)return
    allocate(work(max(1,2*n)),rwork(max(1,3*n-2)),stat=stat)
    ok=consensus(comm,stat==0);if(.not.ok)return
    call zheev('V','U',n,matrix,n,values,work,size(work),rwork,info)
    ok=consensus(comm,info==0.and.finite(matrix).and.all(ieee_is_finite(values)))
  end subroutine

  logical function hermitian(matrix)
    complex(real64),intent(in)::matrix(:,:)
    hermitian=finite(matrix)
    if(hermitian)hermitian=maxval(abs(matrix-conjg(transpose(matrix))))<=&
      4096d0*epsilon(1d0)*max(1d0,maxval(abs(matrix)))*size(matrix,1)
  end function
  logical function finite(matrix)
    complex(real64),intent(in)::matrix(:,:)
    finite=all(ieee_is_finite(real(matrix,real64))).and.all(ieee_is_finite(aimag(matrix)))
  end function
  logical function consensus(comm,valid)
    integer,intent(in)::comm
    logical,intent(in)::valid
    integer::bad,ierr
    call MPI_Allreduce(merge(0,1,valid),bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    consensus=ierr==MPI_SUCCESS.and.bad==0
  end function
  logical function agree_int(comm,value)
    integer,intent(in)::comm,value
    agree_int=agree_bits(comm,int(value,int64))
  end function
  logical function agree_real(comm,value)
    integer,intent(in)::comm
    real(real64),intent(in)::value
    agree_real=agree_bits(comm,transfer(value,0_int64))
  end function
  logical function agree_bits(comm,value)
    integer,intent(in)::comm
    integer(int64),intent(in)::value
    integer(int64)::lo,hi
    integer::ierr
    call MPI_Allreduce(value,lo,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    agree_bits=ierr==MPI_SUCCESS
    call MPI_Allreduce(value,hi,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    agree_bits=agree_bits.and.ierr==MPI_SUCCESS.and.lo==hi
  end function
  subroutine suspend_traps(halting)
    logical,intent(out)::halting(3)
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
  end subroutine
  subroutine restore_traps(halting)
    logical,intent(in)::halting(3)
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1))
    call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  end subroutine
end module dg_hybrid_fragment_subspace
