module dg_hybrid_schwarz_solver
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state
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
  public::advance_dg_hybrid_schwarz_epoch
contains
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
