program test_dg_hybrid_schwarz_solver_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state,initialize_dg_hybrid_schwarz_state
  use dg_hybrid_schwarz_solver,only:advance_dg_hybrid_schwarz_epoch
  implicit none
  integer::ierr,rank,nproc,fragment,nb,ncandidate,j,steps
  integer(int64),allocatable::candidate_ids(:)
  real(real64),allocatable::candidate_energies(:),diagonal(:)
  complex(real64),allocatable::candidate_vectors(:,:),before(:,:),metric_action(:,:),gram(:,:)
  type(s_dg_hybrid_schwarz_state)::state
  logical::ok,converged,rolled_back,inject_failure,fail_operator
  real(real64)::residual,orthogonality
  character(512)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(any(nproc==[2,4,8]),'test requires 2, 4, or 8 ranks')
  fragment=rank+1;nb=4+rank;ncandidate=4
  allocate(candidate_ids(ncandidate),candidate_energies(ncandidate),candidate_vectors(nb,ncandidate),diagonal(nb))
  candidate_ids=[(int(100*fragment+j,int64),j=1,ncandidate)]
  candidate_energies=[0d0,0.05d0,0.20d0,0.40d0]
  candidate_vectors=(0d0,0d0)
  do j=1,ncandidate
    candidate_vectors(modulo(j-1,nb)+1,j)=cmplx(1d0,0.02d0*fragment*j,real64)
  enddo
  diagonal=[(0.4d0+0.13d0*rank+0.17d0*j,j=1,nb)]
  call initialize_dg_hybrid_schwarz_state(MPI_COMM_WORLD,fragment,nproc,11,2d0,300d0,2d0,1,&
    1d-12,1d-10,77123_int64,88231_int64,candidate_ids,candidate_energies,candidate_vectors,&
    state,ok,message)
  call require(ok,'solver fixture state rejected: '//trim(message))
  inject_failure=.false.;fail_operator=.false.
  call advance_dg_hybrid_schwarz_epoch(MPI_COMM_WORLD,11,3,1d-14,1d-10,10d0,&
    apply_h,apply_s,precondition,state,steps,residual,orthogonality,converged,rolled_back,ok,message)
  call require(ok,'bounded Schwarz update failed: '//trim(message))
  call require(steps>=1.and.steps<=3,'bounded solver violated the one-to-three step cap')
  call require(.not.rolled_back,'valid bounded update reported rollback')
  call require(state%coefficient_epoch==steps,'accepted coefficient epoch does not match steps')
  call require(orthogonality<=1d-10,'accepted columns are not globally S-orthonormal')
  call require(ieee_finite(residual),'solver residual is not finite')

  allocate(before,source=state%coefficients)
  inject_failure=rank==0
  call advance_dg_hybrid_schwarz_epoch(MPI_COMM_WORLD,11,1,1d-14,1d-10,10d0,&
    apply_h,apply_s,precondition,state,steps,residual,orthogonality,converged,rolled_back,ok,message,&
    local_publish_ok=.not.inject_failure)
  call require(.not.ok.and.rolled_back,'rank-local publish failure was not rolled back collectively')
  call require(all(state%coefficients==before),'failed step changed accepted coefficients')
  call require(state%coefficient_epoch>=1,'failed step erased the accepted epoch')

  fail_operator=rank==0
  call advance_dg_hybrid_schwarz_epoch(MPI_COMM_WORLD,11,1,1d-14,1d-10,10d0,&
    apply_h,apply_s,precondition,state,steps,residual,orthogonality,converged,rolled_back,ok,message)
  call require(.not.ok.and.rolled_back,'rank-local Hamiltonian failure did not return failure and rollback')
  call require(all(state%coefficients==before),'Hamiltonian failure changed accepted coefficients')

  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid Schwarz solver on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine apply_h(input,output,success)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::success
    output=input*spread(diagonal,2,size(input,2));success=.not.fail_operator
  end subroutine apply_h
  subroutine apply_s(input,output,success)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::success
    output=input;success=.true.
  end subroutine apply_s
  subroutine precondition(input,output,success)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::success
    output=input/spread(diagonal+1d0,2,size(input,2));success=.true.
  end subroutine precondition
  logical function ieee_finite(value)
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    real(real64),intent(in)::value
    ieee_finite=ieee_is_finite(value)
  end function ieee_finite
  subroutine require(condition,detail)
    logical,intent(in)::condition
    character(*),intent(in)::detail
    integer::local_failure,global_failure,code
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,code)
    if(global_failure/=0)then
      if(.not.condition)write(0,'(a,i0,2a)')'rank ',rank,': ',trim(detail)
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif
  end subroutine require
end program test_dg_hybrid_schwarz_solver_mpi
