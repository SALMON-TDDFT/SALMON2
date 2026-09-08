module dg_hybrid_fragment_thermal
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,&
    ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_admission,only:s_dg_hybrid_admission_report
  use dg_hybrid_projected_fragment_pipeline,only:s_dg_hybrid_projection_factorization_receipt,&
    validate_dg_hybrid_projected_basis,dg_hybrid_core_quadrature_binding
  use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,s_dg_hybrid_fragment_candidate_catalog,&
    s_dg_hybrid_fragment_epoch_budget,s_dg_hybrid_fragment_extension_receipt,&
    advance_dg_hybrid_fragment_epoch,extend_dg_hybrid_fragment_subspace
  use dg_hybrid_fragment_solver,only:measure_dg_hybrid_fragment_core_norms,reconstruct_dg_hybrid_fragment_density
  use dc_fragment_occupation,only:run_dc_fragment_occupation_epoch
  implicit none
  private
  type,public::s_dg_hybrid_thermal_state
    ! valid means a consistent unmixed thermal density, not SCF convergence.
    logical::valid=.false.
    integer::epoch=0
    type(s_dg_hybrid_fragment_subspace_state)::state
    type(s_dg_hybrid_fragment_candidate_catalog)::candidates
    real(real64),allocatable::energies(:),occupations(:),density(:)
    real(real64)::chemical_potential=0d0,electron_count=0d0,electron_defect=huge(1d0)
  end type
  abstract interface
    subroutine apply_interface(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine
    subroutine precondition_interface(input,energies,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      real(real64),intent(in)::energies(:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine
  end interface
  public::advance_dg_hybrid_thermal_state
contains
  subroutine advance_dg_hybrid_thermal_state(comm,basis,receipt,admission,core_mask,point_weights,&
      epoch,maximum_steps,temperature,wspin,target,electron_tolerance,energy_tolerance,orthogonality_tolerance,&
      intermediate_tolerance,allowed_residual_growth,apply_h,apply_s,precondition,budget,accepted,&
      iterations,passes,extensions,ok,message)
    ! One MPI rank per fragment. H/S/precondition callbacks are local only.
    ! Staging protects every published field. The persistent budget is NOT
    ! staged: attempted updates remain consumed even when another rank fails.
    integer,intent(in)::comm,epoch,maximum_steps
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    type(s_dg_hybrid_projection_factorization_receipt),intent(in)::receipt
    type(s_dg_hybrid_admission_report),intent(in)::admission
    logical,intent(in)::core_mask(:)
    real(real64),intent(in)::point_weights(:),temperature,wspin,target,electron_tolerance,energy_tolerance,&
      orthogonality_tolerance,intermediate_tolerance,allowed_residual_growth
    procedure(apply_interface)::apply_h,apply_s
    procedure(precondition_interface)::precondition
    type(s_dg_hybrid_fragment_epoch_budget),intent(inout)::budget
    type(s_dg_hybrid_thermal_state),allocatable,intent(inout)::accepted
    integer,intent(out)::iterations,passes,extensions
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_thermal_state),allocatable::work
    integer::n,np,npoint,ierr,stat,j,f,controls(2),lo(2),hi(2)
    integer,allocatable::owners(:)
    integer(int64),allocatable::rows(:)
    complex(real64),allocatable::metric(:,:)
    real(real64),allocatable::last_weights(:)
    real(real64)::parameters(8),pmin(8),pmax(8),counts(2),totals(2),local_ne
    logical::valid,halting(3)
    character(512)::diagnostic
    ok=.false.;message='';iterations=0;passes=0;extensions=0
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call execute()
    call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1))
    call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  contains
    subroutine execute()
      call MPI_Comm_size(comm,np,ierr)
      if(ierr/=MPI_SUCCESS)then;message='thermal communicator size failed';return;endif
      call status(allocated(accepted),'missing thermal entry state');if(.not.ok)return
      parameters=[temperature,wspin,target,electron_tolerance,energy_tolerance,orthogonality_tolerance,&
        intermediate_tolerance,allowed_residual_growth]
      call status(all(ieee_is_finite(parameters)),'nonfinite thermal controls');if(.not.ok)return
      call MPI_Allreduce(parameters,pmin,8,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='thermal control minimum failed';return;endif
      call MPI_Allreduce(parameters,pmax,8,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call status(ierr==MPI_SUCCESS.and.all(pmin==pmax),'rank-disagreeing thermal controls');if(.not.ok)return
      controls=[epoch,maximum_steps]
      call MPI_Allreduce(controls,lo,2,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='thermal epoch minimum failed';return;endif
      call MPI_Allreduce(controls,hi,2,MPI_INTEGER,MPI_MAX,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.all(lo==hi).and.epoch>0.and.maximum_steps>=1.and.maximum_steps<=256.and.&
        temperature>=0d0.and.wspin>0d0.and.target>=0d0.and.electron_tolerance>0d0.and.energy_tolerance>=0d0.and.&
        orthogonality_tolerance>=64d0*epsilon(1d0).and.orthogonality_tolerance<=1d-2.and.&
        intermediate_tolerance>0d0.and.allowed_residual_growth>=1d0
      call status(valid,'invalid thermal epoch controls');if(.not.ok)return
      call validate_dg_hybrid_projected_basis(comm,basis,receipt,receipt%wannier_fingerprint,ok,message)
      if(.not.ok)return
      n=size(basis%global_ids);npoint=size(basis%buffer_point_ids);f=basis%fragment_id
      valid=f>=1.and.f<=np.and.size(core_mask)==npoint.and.size(point_weights)==npoint.and.&
        admission%trial_prepared.and.admission%support_measured.and.&
        admission%basis_fingerprint==receipt%payload_fingerprint.and.admission%metric_fingerprint/=0_int64.and.&
        accepted%state%basis_fingerprint==admission%basis_fingerprint.and.&
        accepted%state%metric_fingerprint==admission%metric_fingerprint.and.&
        accepted%state%fragment_id==f.and.accepted%state%basis_generation==basis%generation.and.&
        allocated(accepted%candidates%used)
      call status(valid,'thermal admission/state context mismatch');if(.not.ok)return
      valid=all(ieee_is_finite(point_weights)).and.all(point_weights>0d0).and.any(core_mask)
      call status(valid,'invalid thermal quadrature');if(.not.ok)return
      valid=receipt%core_quadrature_fingerprint==dg_hybrid_core_quadrature_binding(&
        pack(basis%buffer_point_ids,core_mask),pack(point_weights,core_mask))
      call status(valid,'thermal core quadrature differs from admitted metric');if(.not.ok)return
      allocate(owners(np),rows(n),metric(n,n),stat=stat)
      call status(stat==0,'thermal staging allocation failed');if(.not.ok)return
      owners=0;owners(f)=1
      call MPI_Allreduce(MPI_IN_PLACE,owners,np,MPI_INTEGER,MPI_SUM,comm,ierr)
      call status(ierr==MPI_SUCCESS.and.all(owners==1),'thermal path requires one rank per fragment')
      if(.not.ok)return
      rows=[(int(j,int64),j=1,n)]
      metric=matmul(conjg(transpose(basis%buffer_values)),&
        basis%buffer_values*spread(merge(point_weights,0d0,core_mask),2,n))
      call status(all(ieee_is_finite(real(metric))).and.all(ieee_is_finite(aimag(metric))),&
        'nonfinite thermal core metric');if(.not.ok)return
      allocate(work,source=accepted,stat=stat)
      call status(stat==0,'thermal state copy allocation failed');if(.not.ok)return
      call run_dc_fragment_occupation_epoch(comm,np,f,.true.,epoch,n,temperature,wspin,target,electron_tolerance,&
        refresh,extend,work%occupations,work%chemical_potential,work%electron_count,passes,extensions,ok,message)
      if(.not.ok)return
      if(allocated(work%density))deallocate(work%density)
      allocate(work%density(npoint),stat=stat)
      call status(stat==0,'thermal density allocation failed');if(.not.ok)return
      call reconstruct_dg_hybrid_fragment_density(MPI_COMM_SELF,basis,work%state%vectors,work%occupations,wspin,&
        core_mask,point_weights,work%density,local_ne,valid,diagnostic)
      call status(valid,diagnostic);if(.not.ok)return
      counts=[sum(point_weights*work%density),sum(last_weights*work%occupations)]
      call status(all(ieee_is_finite(counts)).and.abs(counts(1)-local_ne)<=electron_tolerance,&
        'thermal local density/count inconsistency');if(.not.ok)return
      call MPI_Allreduce(counts,totals,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.all(ieee_is_finite(totals))
      valid=valid.and.all(abs(totals-work%electron_count)<=electron_tolerance)
      call status(valid,'thermal density/occupation electron inconsistency');if(.not.ok)return
      ! Target drift is diagnostic here. Outer divided SCF controls convergence.
      work%electron_defect=abs(totals(1)-target);work%electron_count=totals(1)
      work%epoch=epoch;work%valid=.true.
      call move_alloc(work,accepted)
    end subroutine
    subroutine refresh(current_epoch,energies,weights,can_grow,success,detail)
      integer,intent(in)::current_epoch
      real(real64),allocatable,intent(out)::energies(:),weights(:)
      logical,intent(out)::can_grow,success
      character(*),intent(out)::detail
      integer::steps,remaining,allocation_status
      integer(int64)::workspace,fp
      real(real64)::residual
      logical::converged,advanced
      character(256)::reason
      can_grow=.false.;success=.false.
      allocate(energies(work%state%state_count),stat=allocation_status)
      detail='thermal spectrum allocation failed';if(allocation_status/=0)return
      call advance_dg_hybrid_fragment_epoch(MPI_COMM_SELF,n,rows,f,basis%generation,&
        admission%basis_fingerprint,admission%metric_fingerprint,current_epoch,maximum_steps,&
        apply_h,checked_s,precondition,intermediate_tolerance,orthogonality_tolerance,allowed_residual_growth,&
        budget,work%state,energies,steps,remaining,residual,converged,advanced,reason,workspace,fp,success,detail)
      iterations=iterations+steps
      if(.not.success)return
      call measure_dg_hybrid_fragment_core_norms(MPI_COMM_SELF,basis,work%state%vectors,&
        core_mask,point_weights,weights,success,detail)
      if(.not.success)return
      if(allocated(work%energies))deallocate(work%energies)
      if(allocated(last_weights))deallocate(last_weights)
      allocate(work%energies(size(energies)),last_weights(size(weights)),stat=allocation_status)
      success=allocation_status==0;detail='thermal spectrum staging failed';if(.not.success)return
      work%energies=energies;last_weights=weights
      can_grow=.not.all(work%candidates%used);detail=''
    end subroutine
    subroutine extend(current_epoch,old_count,new_count,success,detail)
      integer,intent(in)::current_epoch,old_count
      integer,intent(out)::new_count
      logical,intent(out)::success
      character(*),intent(out)::detail
      type(s_dg_hybrid_fragment_extension_receipt)::extension
      success=current_epoch==epoch.and.old_count==work%state%state_count;new_count=old_count
      detail='thermal extension context mismatch';if(.not.success)return
      call extend_dg_hybrid_fragment_subspace(MPI_COMM_SELF,n,rows,f,basis%generation,&
        admission%basis_fingerprint,admission%metric_fingerprint,apply_h,checked_s,energy_tolerance,&
        orthogonality_tolerance,work%candidates,work%state,extension,success,detail)
      new_count=work%state%state_count
    end subroutine
    subroutine checked_s(input,output,success)
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::success
      complex(real64)::reference(size(output,1),size(output,2))
      real(real64)::scale
      call apply_s(input,output,success);if(.not.success)return
      reference=matmul(metric,input)
      success=all(ieee_is_finite(real(output))).and.all(ieee_is_finite(aimag(output))).and.&
        all(ieee_is_finite(real(reference))).and.all(ieee_is_finite(aimag(reference)))
      if(.not.success)return
      scale=max(1d0,maxval(abs(reference)))
      success=maxval(abs(output/scale-reference/scale))<=orthogonality_tolerance
    end subroutine
    subroutine status(local_ok,detail)
      logical,intent(in)::local_ok
      character(*),intent(in)::detail
      integer::rank,failed,first_failed,code
      character(512)::shared
      ok=.false.;message='thermal status rank lookup failed'
      call MPI_Comm_rank(comm,rank,code);if(code/=MPI_SUCCESS)return
      failed=huge(0);if(.not.local_ok)failed=rank
      call MPI_Allreduce(failed,first_failed,1,MPI_INTEGER,MPI_MIN,comm,code)
      message='thermal status reduction failed';if(code/=MPI_SUCCESS)return
      ok=first_failed==huge(0);message='';if(ok)return
      shared='';if(rank==first_failed)shared=detail
      call MPI_Bcast(shared,len(shared),MPI_CHARACTER,first_failed,comm,code)
      message='thermal diagnostic broadcast failed';if(code/=MPI_SUCCESS)return
      message=trim(shared)
    end subroutine
  end subroutine advance_dg_hybrid_thermal_state
end module dg_hybrid_fragment_thermal
