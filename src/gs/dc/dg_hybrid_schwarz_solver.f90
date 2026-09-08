module dg_hybrid_schwarz_solver
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state,extend_dg_hybrid_schwarz_state
  implicit none
  private
  abstract interface
    subroutine apply_interface(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine apply_interface
  end interface
  interface
    subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
      import real64
      character(1),intent(in)::jobz,uplo
      integer,intent(in)::n,lda,lwork
      complex(real64),intent(inout)::a(lda,*),work(*)
      real(real64),intent(out)::w(*),rwork(*)
      integer,intent(out)::info
    end subroutine zheev
  end interface
  public::advance_dg_hybrid_schwarz_epoch,assign_dg_hybrid_schwarz_occupations
contains
  subroutine assign_dg_hybrid_schwarz_occupations(comm,basis_generation,temperature,wspin,target,&
      tail_tolerance,degeneracy_tolerance,candidate_ids,candidate_energies,candidate_vectors,&
      apply_h,apply_s,state,occupations,energies,extended,ok,message)
    integer,intent(in)::comm,basis_generation
    real(real64),intent(in)::temperature,wspin,target,tail_tolerance,degeneracy_tolerance
    integer(kind=8),intent(in)::candidate_ids(:)
    real(real64),intent(in)::candidate_energies(:)
    complex(real64),intent(in)::candidate_vectors(:,:)
    procedure(apply_interface)::apply_h,apply_s
    type(s_dg_hybrid_schwarz_state),intent(inout)::state
    real(real64),allocatable,intent(out)::occupations(:),energies(:)
    logical,intent(out)::extended,ok
    character(*),intent(out)::message
    type(s_dg_hybrid_schwarz_state)::work_state
    complex(real64),allocatable::hcoeff(:,:),projected(:,:),lapack_work(:)
    real(real64),allocatable::rwork(:),minimum_candidates(:),maximum_candidates(:)
    real(real64)::orthogonality,mu,electron_count
    integer::n,lwork,stat,ierr,info,requested
    logical::success,tail_resolved

    ok=.false.;message='';extended=.false.
    success=state%valid.and.state%basis_generation==basis_generation.and.temperature>=0d0.and.&
      wspin>0d0.and.target>0d0.and.tail_tolerance>0d0.and.tail_tolerance<1d0.and.&
      degeneracy_tolerance>=0d0.and.size(candidate_ids)==state%candidate_count.and.&
      size(candidate_energies)==state%candidate_count.and.size(candidate_vectors,1)==state%local_basis_count.and.&
      size(candidate_vectors,2)==state%candidate_count.and.all(ieee_is_finite(candidate_energies)).and.&
      all(candidate_energies(2:)>=candidate_energies(:size(candidate_energies)-1))
    call collective_gate(comm,success,'invalid global Schwarz occupation context',ok,message);if(.not.ok)return
    allocate(minimum_candidates(size(candidate_energies)),maximum_candidates(size(candidate_energies)),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz occupation candidate allocation failed',ok,message);if(.not.ok)return
    call MPI_Allreduce(candidate_energies,minimum_candidates,size(candidate_energies),MPI_DOUBLE_PRECISION,&
      MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='Schwarz occupation energy minimum failed';return;endif
    call MPI_Allreduce(candidate_energies,maximum_candidates,size(candidate_energies),MPI_DOUBLE_PRECISION,&
      MPI_MAX,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.all(minimum_candidates==maximum_candidates),&
      'Schwarz occupation candidate energies differ between ranks',ok,message);if(.not.ok)return
    work_state=state;n=work_state%trial_count;lwork=max(1,2*n*n)
    call orthonormalize_columns(comm,apply_s,work_state%coefficients,max(64d0*epsilon(1d0),1d-12),&
      orthogonality,success,message)
    if(.not.success)then;ok=.false.;return;endif
    allocate(hcoeff(work_state%local_basis_count,n),projected(n,n),energies(n),occupations(n),&
      lapack_work(lwork),rwork(max(1,3*n-2)),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz occupation workspace allocation failed',ok,message);if(.not.ok)return
    call apply_h(work_state%coefficients,hcoeff,success)
    call collective_gate(comm,success.and.finite_matrix(hcoeff),'Schwarz occupation Hamiltonian failed',ok,message)
    if(.not.ok)return
    projected=matmul(conjg(transpose(work_state%coefficients)),hcoeff)
    call MPI_Allreduce(MPI_IN_PLACE,projected,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.finite_matrix(projected),&
      'Schwarz occupation Rayleigh reduction failed',ok,message);if(.not.ok)return
    projected=0.5d0*(projected+conjg(transpose(projected)))
    call zheev('V','U',n,projected,n,energies,lapack_work,lwork,rwork,info)
    call collective_gate(comm,info==0.and.all(ieee_is_finite(energies)),&
      'Schwarz occupation eigensystem failed',ok,message);if(.not.ok)return
    work_state%coefficients=matmul(work_state%coefficients,projected)
    call solve_fermi_occupations(energies,temperature,wspin,target,tail_tolerance,&
      occupations,mu,electron_count,tail_resolved,success)
    call collective_gate(comm,success,'Schwarz finite-temperature occupation solve failed',ok,message)
    if(.not.ok)return
    if(.not.tail_resolved)then
      call collective_gate(comm,n<work_state%candidate_count,&
        'Schwarz candidate capacity exhausted before resolving the 300 K occupation tail',ok,message)
      if(.not.ok)return
      requested=n+1
      do while(requested<work_state%candidate_count)
        if(abs(candidate_energies(requested+1)-candidate_energies(requested))>degeneracy_tolerance)exit
        requested=requested+1
      enddo
      call extend_dg_hybrid_schwarz_state(comm,basis_generation,requested,work_state%mapping_fingerprint,&
        work_state%candidate_fingerprint,candidate_ids,candidate_vectors,work_state,ok,message)
      if(.not.ok)return
      state=work_state;extended=.true.;deallocate(occupations,energies)
      allocate(occupations(0),energies(0));ok=.true.;message='';return
    endif
    if(allocated(work_state%energies))deallocate(work_state%energies)
    if(allocated(work_state%occupations))deallocate(work_state%occupations)
    allocate(work_state%energies(n),work_state%occupations(n),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz occupation publication allocation failed',ok,message)
    if(.not.ok)return
    work_state%energies=energies;work_state%occupations=occupations
    work_state%chemical_potential=mu;work_state%electron_count=electron_count
    work_state%electron_defect=abs(electron_count-target)
    state=work_state;ok=.true.;message=''
  end subroutine assign_dg_hybrid_schwarz_occupations

  subroutine solve_fermi_occupations(energies,temperature,wspin,target,tail_tolerance,occupations,&
      mu,electron_count,tail_resolved,ok)
    real(real64),intent(in)::energies(:),temperature,wspin,target,tail_tolerance
    real(real64),intent(out)::occupations(:),mu,electron_count
    logical,intent(out)::tail_resolved,ok
    real(real64),parameter::boltzmann_hartree_per_kelvin=3.166811563d-6
    real(real64)::lower,upper,mid,kbt,total
    integer::iteration,j,occupied
    ok=.false.;tail_resolved=.false.;mu=0d0;electron_count=0d0;occupations=0d0
    if(size(energies)<1.or.size(occupations)/=size(energies).or.target>wspin*real(size(energies),real64))return
    if(temperature==0d0)then
      occupied=ceiling(target/wspin-64d0*epsilon(1d0));occupations(:occupied)=1d0
      mu=energies(min(size(energies),max(1,occupied)))
    else
      kbt=boltzmann_hartree_per_kelvin*temperature
      lower=energies(1)-max(1d0,64d0*kbt);upper=energies(size(energies))+max(1d0,64d0*kbt)
      do iteration=1,256
        mid=0.5d0*(lower+upper)
        total=wspin*sum([(fermi_value((energies(j)-mid)/kbt),j=1,size(energies))])
        if(total<target)then;lower=mid;else;upper=mid;endif
      enddo
      mu=0.5d0*(lower+upper)
      occupations=[(fermi_value((energies(j)-mu)/kbt),j=1,size(energies))]
    endif
    electron_count=wspin*sum(occupations);tail_resolved=occupations(size(occupations))<=tail_tolerance
    ok=all(ieee_is_finite(occupations)).and.ieee_is_finite(mu).and.ieee_is_finite(electron_count).and.&
      abs(electron_count-target)<=max(1d-12,256d0*epsilon(1d0)*max(1d0,target))
  end subroutine solve_fermi_occupations

  pure real(real64) function fermi_value(x)result(value)
    real(real64),intent(in)::x
    if(x>=50d0)then;value=exp(-x)
    elseif(x<=-50d0)then;value=1d0
    else;value=1d0/(1d0+exp(x))
    endif
  end function fermi_value

  subroutine advance_dg_hybrid_schwarz_epoch(comm,basis_generation,maximum_steps,residual_tolerance,&
      orthogonality_tolerance,allowed_residual_growth,apply_h,apply_s,precondition,state,iterations,&
      maximum_residual,orthogonality_defect,converged,rolled_back,ok,message,local_publish_ok)
    integer,intent(in)::comm,basis_generation,maximum_steps
    real(real64),intent(in)::residual_tolerance,orthogonality_tolerance,allowed_residual_growth
    procedure(apply_interface)::apply_h,apply_s,precondition
    type(s_dg_hybrid_schwarz_state),intent(inout)::state
    integer,intent(out)::iterations
    real(real64),intent(out)::maximum_residual,orthogonality_defect
    logical,intent(out)::converged,rolled_back,ok
    character(*),intent(out)::message
    logical,intent(in),optional::local_publish_ok
    type(s_dg_hybrid_schwarz_state)::work_state
    complex(real64),allocatable::residual(:,:),preconditioned(:,:),direction(:,:),previous_direction(:,:),trial(:,:)
    real(real64)::current_residual,trial_residual,previous_residual,alpha,beta
    integer::step,attempt,stat
    logical::success,accepted,publish

    ok=.false.;message='';iterations=0;maximum_residual=huge(1d0)
    orthogonality_defect=huge(1d0);converged=.false.;rolled_back=.false.
    success=state%valid.and.state%basis_generation==basis_generation.and.maximum_steps>=1.and.maximum_steps<=256.and.&
      residual_tolerance>0d0.and.orthogonality_tolerance>=64d0*epsilon(1d0).and.&
      allowed_residual_growth>=1d0.and.ieee_is_finite(residual_tolerance).and.&
      ieee_is_finite(orthogonality_tolerance).and.ieee_is_finite(allowed_residual_growth).and.&
      allocated(state%coefficients).and.size(state%coefficients,1)==state%local_basis_count.and.&
      size(state%coefficients,2)==state%trial_count
    call collective_gate(comm,success,'invalid bounded Schwarz solver context',ok,message)
    if(.not.ok)return
    work_state=state
    call orthonormalize_columns(comm,apply_s,work_state%coefficients,orthogonality_tolerance,&
      orthogonality_defect,success,message)
    if(.not.success)then;ok=.false.;rolled_back=.true.;return;endif
    allocate(residual(state%local_basis_count,state%trial_count),&
      preconditioned(state%local_basis_count,state%trial_count),&
      direction(state%local_basis_count,state%trial_count),&
      previous_direction(state%local_basis_count,state%trial_count),&
      trial(state%local_basis_count,state%trial_count),stat=stat)
    call collective_gate(comm,stat==0,'bounded Schwarz solver allocation failed',ok,message)
    if(.not.ok)then;rolled_back=.true.;return;endif
    previous_direction=(0d0,0d0);previous_residual=0d0
    do step=1,maximum_steps
      call compute_residual(comm,apply_h,apply_s,work_state%coefficients,residual,current_residual,success,message)
      if(.not.success)then;ok=.false.;rolled_back=.true.;return;endif
      maximum_residual=current_residual
      if(current_residual<=residual_tolerance)then;converged=.true.;exit;endif
      call precondition(residual,preconditioned,success)
      call collective_gate(comm,success.and.finite_matrix(preconditioned),&
        'bounded Schwarz preconditioner failed',ok,message)
      if(.not.ok)then;rolled_back=.true.;return;endif
      if(step==1.or.previous_residual<=tiny(1d0))then
        direction=-preconditioned
      else
        beta=(current_residual/previous_residual)**2
        direction=-preconditioned+beta*previous_direction
      endif
      call project_direction(comm,apply_s,work_state%coefficients,direction,success,message)
      if(.not.success)then;ok=.false.;rolled_back=.true.;return;endif
      alpha=0.5d0;accepted=.false.
      do attempt=1,8
        trial=work_state%coefficients+alpha*direction
        call orthonormalize_columns(comm,apply_s,trial,orthogonality_tolerance,&
          orthogonality_defect,success,message)
        if(.not.success)then;alpha=0.5d0*alpha;cycle;endif
        call compute_residual(comm,apply_h,apply_s,trial,residual,trial_residual,success,message)
        if(.not.success)then;ok=.false.;rolled_back=.true.;return;endif
        if(trial_residual<=allowed_residual_growth*current_residual)then;accepted=.true.;exit;endif
        alpha=0.5d0*alpha
      enddo
      call collective_gate(comm,accepted,'bounded Schwarz line search rejected every trial',ok,message)
      if(.not.ok)then;rolled_back=.true.;return;endif
      publish=.true.;if(present(local_publish_ok))publish=local_publish_ok
      call collective_gate(comm,publish,'bounded Schwarz coefficient publication rolled back',ok,message)
      if(.not.ok)then;rolled_back=.true.;return;endif
      work_state%coefficients=trial;work_state%coefficient_epoch=work_state%coefficient_epoch+1
      previous_direction=direction;previous_residual=current_residual
      iterations=iterations+1;maximum_residual=trial_residual
      if(trial_residual<=residual_tolerance)then;converged=.true.;exit;endif
    enddo
    call orthonormality_defect(comm,apply_s,work_state%coefficients,orthogonality_defect,success,message)
    if(.not.success.or.orthogonality_defect>orthogonality_tolerance)then
      ok=.false.;rolled_back=.true.
      if(success)message='bounded Schwarz final orthogonality exceeds tolerance'
      return
    endif
    state=work_state;ok=.true.;message=''
  end subroutine advance_dg_hybrid_schwarz_epoch

  subroutine compute_residual(comm,apply_h,apply_s,coefficients,residual,norm,ok,message)
    integer,intent(in)::comm
    procedure(apply_interface)::apply_h,apply_s
    complex(real64),intent(in)::coefficients(:,:)
    complex(real64),intent(out)::residual(:,:)
    real(real64),intent(out)::norm
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::hcoeff(:,:),scoeff(:,:),rayleigh(:,:)
    real(real64)::local_norm,global_norm
    integer::n,stat,ierr
    logical::success
    n=size(coefficients,2);allocate(hcoeff(size(coefficients,1),n),scoeff(size(coefficients,1),n),&
      rayleigh(n,n),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz residual allocation failed',ok,message);if(.not.ok)return
    call apply_h(coefficients,hcoeff,success)
    call collective_gate(comm,success.and.finite_matrix(hcoeff),'Schwarz Hamiltonian application failed',ok,message)
    if(.not.ok)return
    call apply_s(coefficients,scoeff,success)
    call collective_gate(comm,success.and.finite_matrix(scoeff),'Schwarz metric application failed',ok,message)
    if(.not.ok)return
    rayleigh=matmul(conjg(transpose(coefficients)),hcoeff)
    call MPI_Allreduce(MPI_IN_PLACE,rayleigh,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.finite_matrix(rayleigh),&
      'Schwarz Rayleigh reduction failed',ok,message);if(.not.ok)return
    rayleigh=0.5d0*(rayleigh+conjg(transpose(rayleigh)))
    residual=hcoeff-matmul(scoeff,rayleigh)
    local_norm=sum(abs(residual)**2)
    call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    norm=sqrt(max(0d0,global_norm))
    call collective_gate(comm,ierr==MPI_SUCCESS.and.ieee_is_finite(norm).and.finite_matrix(residual),&
      'Schwarz residual norm failed',ok,message)
  end subroutine compute_residual

  subroutine project_direction(comm,apply_s,coefficients,direction,ok,message)
    integer,intent(in)::comm
    procedure(apply_interface)::apply_s
    complex(real64),intent(in)::coefficients(:,:)
    complex(real64),intent(inout)::direction(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::sdirection(:,:),overlap(:,:)
    integer::n,stat,ierr
    logical::success
    n=size(coefficients,2);allocate(sdirection(size(direction,1),n),overlap(n,n),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz direction allocation failed',ok,message);if(.not.ok)return
    call apply_s(direction,sdirection,success)
    call collective_gate(comm,success.and.finite_matrix(sdirection),'Schwarz direction metric failed',ok,message)
    if(.not.ok)return
    overlap=matmul(conjg(transpose(coefficients)),sdirection)
    call MPI_Allreduce(MPI_IN_PLACE,overlap,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.finite_matrix(overlap),&
      'Schwarz direction projection reduction failed',ok,message);if(.not.ok)return
    direction=direction-matmul(coefficients,overlap)
    ok=.true.;message=''
  end subroutine project_direction

  subroutine orthonormalize_columns(comm,apply_s,coefficients,tolerance,defect,ok,message)
    integer,intent(in)::comm
    procedure(apply_interface)::apply_s
    complex(real64),intent(inout)::coefficients(:,:)
    real(real64),intent(in)::tolerance
    real(real64),intent(out)::defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::scoeff(:,:),gram(:,:),lapack_work(:),transform(:,:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    integer::n,stat,ierr,info,lwork,j
    logical::success
    n=size(coefficients,2);lwork=max(1,2*n*n)
    allocate(scoeff(size(coefficients,1),n),gram(n,n),transform(n,n),eigenvalues(n),&
      lapack_work(lwork),rwork(max(1,3*n-2)),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz orthogonalization allocation failed',ok,message);if(.not.ok)return
    call apply_s(coefficients,scoeff,success)
    call collective_gate(comm,success.and.finite_matrix(scoeff),'Schwarz orthogonalization metric failed',ok,message)
    if(.not.ok)return
    gram=matmul(conjg(transpose(coefficients)),scoeff)
    call MPI_Allreduce(MPI_IN_PLACE,gram,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call collective_gate(comm,ierr==MPI_SUCCESS.and.finite_matrix(gram),&
      'Schwarz Gram reduction failed',ok,message);if(.not.ok)return
    gram=0.5d0*(gram+conjg(transpose(gram)))
    call zheev('V','U',n,gram,n,eigenvalues,lapack_work,lwork,rwork,info)
    call collective_gate(comm,info==0.and.all(ieee_is_finite(eigenvalues)).and.&
      minval(eigenvalues)>tolerance,'Schwarz metric rank was lost',ok,message);if(.not.ok)return
    transform=gram
    do j=1,n;transform(:,j)=transform(:,j)/sqrt(eigenvalues(j));enddo
    transform=matmul(transform,conjg(transpose(gram)))
    coefficients=matmul(coefficients,transform)
    call orthonormality_defect(comm,apply_s,coefficients,defect,ok,message)
  end subroutine orthonormalize_columns

  subroutine orthonormality_defect(comm,apply_s,coefficients,defect,ok,message)
    integer,intent(in)::comm
    procedure(apply_interface)::apply_s
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(out)::defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::scoeff(:,:),gram(:,:)
    integer::n,j,stat,ierr
    logical::success
    n=size(coefficients,2);allocate(scoeff(size(coefficients,1),n),gram(n,n),stat=stat)
    call collective_gate(comm,stat==0,'Schwarz defect allocation failed',ok,message);if(.not.ok)return
    call apply_s(coefficients,scoeff,success)
    call collective_gate(comm,success.and.finite_matrix(scoeff),'Schwarz defect metric failed',ok,message)
    if(.not.ok)return
    gram=matmul(conjg(transpose(coefficients)),scoeff)
    call MPI_Allreduce(MPI_IN_PLACE,gram,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr==MPI_SUCCESS)then
      do j=1,n;gram(j,j)=gram(j,j)-1d0;enddo
      defect=maxval(abs(gram))
    else
      defect=huge(1d0)
    endif
    call collective_gate(comm,ierr==MPI_SUCCESS.and.ieee_is_finite(defect),&
      'Schwarz orthogonality reduction failed',ok,message)
  end subroutine orthonormality_defect

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
    ok=.false.;message='Schwarz solver status rank query failed'
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    failed=huge(0);if(.not.local_ok)failed=rank
    call MPI_Allreduce(failed,first_failed,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz solver status reduction failed';return;endif
    if(first_failed==huge(0))then;ok=.true.;message='';return;endif
    shared='';if(rank==first_failed)shared=detail
    call MPI_Bcast(shared,len(shared),MPI_CHARACTER,first_failed,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Schwarz solver diagnostic broadcast failed';return;endif
    message=trim(shared)
  end subroutine collective_gate
end module dg_hybrid_schwarz_solver
