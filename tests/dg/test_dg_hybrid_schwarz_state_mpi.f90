program test_dg_hybrid_schwarz_state_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state,&
    initialize_dg_hybrid_schwarz_state,extend_dg_hybrid_schwarz_state,&
    validate_dg_hybrid_schwarz_dynamic_receipt
  implicit none
  integer::ierr,rank,nproc,fragment,nb,ncandidate,i,j,lo,hi,original_thermal_tail_count
  integer(int64)::mapping_fingerprint,candidate_fingerprint,accepted_fingerprint,dynamic_receipt
  integer(int64),allocatable::candidate_ids(:),minimum_ids(:),maximum_ids(:)
  real(real64),allocatable::candidate_energies(:)
  real(real64)::original_occupation
  complex(real64),allocatable::candidate_vectors(:,:)
  type(s_dg_hybrid_schwarz_state)::state,failed_state
  logical::ok
  character(512)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(any(nproc==[2,4,8]),'test requires 2, 4, or 8 ranks')
  fragment=rank+1;nb=5+rank;ncandidate=8
  mapping_fingerprint=73491_int64;candidate_fingerprint=92821_int64
  allocate(candidate_ids(ncandidate),candidate_energies(ncandidate),candidate_vectors(nb,ncandidate))
  candidate_ids=[(int(100*fragment+j,int64),j=1,ncandidate)]
  candidate_energies=[0d0,0.02d0,0.08d0,0.08d0,0.20d0,0.32d0,0.46d0,0.62d0]
  candidate_vectors=(0d0,0d0)
  do j=1,ncandidate
    i=modulo(j-1,nb)+1
    candidate_vectors(i,j)=cmplx(1d0,0.01d0*real(fragment*j,real64),real64)
  enddo

  call initialize_dg_hybrid_schwarz_state(MPI_COMM_WORLD,fragment,nproc,7,4d0,300d0,2d0,1,&
    1d-12,1d-10,mapping_fingerprint,candidate_fingerprint,candidate_ids,candidate_energies,&
    candidate_vectors,state,ok,message)
  call require(ok,'valid common inventory rejected: '//trim(message))
  call require(state%valid,'valid state was not published')
  call require(state%trial_count==4,'electron/tail/guard/degeneracy count is not four')
  call require(state%local_basis_count==nb,'unequal local basis size was not retained')
  call require(all(shape(state%coefficients)==[nb,4]),'local coefficient row block has wrong shape')
  call require(allocated(state%energies).and.allocated(state%occupations),&
    'initial state must retain common occupations for rollback diagnostics')
  call require(size(state%energies)==4.and.size(state%occupations)==4,&
    'initial diagnostic spectrum has the wrong trial size')
  call require(abs(state%electron_count-4d0)<1d-10.and.state%electron_defect<1d-10,&
    'initial common-occupation electron count is not diagnostic-ready')
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(ok.and.dynamic_receipt/=0_int64,'valid dynamic Schwarz receipt rejected: '//trim(message))
  original_occupation=state%occupations(1)
  if(rank==0)state%occupations(1)=state%occupations(1)+1d-3
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(.not.ok,'rank-local dynamic occupation perturbation was accepted')
  if(rank==0)state%occupations(1)=original_occupation
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(ok,'restored dynamic Schwarz receipt rejected: '//trim(message))
  original_thermal_tail_count=state%thermal_tail_count
  if(rank==0)state%thermal_tail_count=state%thermal_tail_count+1
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(.not.ok,'rank-local common integer perturbation was accepted')
  if(rank==0)state%thermal_tail_count=original_thermal_tail_count
  state%fragment_count=nproc+1
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(.not.ok,'collectively wrong fragment count was accepted')
  state%fragment_count=nproc
  state%fragment_id=modulo(rank+1,nproc)+1
  call validate_dg_hybrid_schwarz_dynamic_receipt(MPI_COMM_WORLD,state,dynamic_receipt,ok,message)
  call require(.not.ok,'fragment ID/rank mismatch was accepted')
  state%fragment_id=rank+1
  call require(size(state%column_ids)==4.and.all(state%column_ids/=0_int64),&
    'common column IDs are missing')
  allocate(minimum_ids(4),maximum_ids(4))
  call MPI_Allreduce(state%column_ids,minimum_ids,4,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(state%column_ids,maximum_ids,4,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierr)
  call require(all(minimum_ids==maximum_ids),'column IDs differ between fragments')
  call require(all([(count(state%column_ids==state%column_ids(j))==1,j=1,4)]),&
    'column IDs are not unique')
  call MPI_Allreduce(state%trial_count,lo,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(state%trial_count,hi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  call require(lo==hi,'trial count differs between ranks')
  accepted_fingerprint=state%fingerprint
  call require(accepted_fingerprint/=0_int64,'state fingerprint is zero')

  call extend_dg_hybrid_schwarz_state(MPI_COMM_WORLD,7,5,mapping_fingerprint,candidate_fingerprint,&
    candidate_ids,candidate_vectors,state,ok,message)
  call require(ok,'valid common extension rejected: '//trim(message))
  call require(state%trial_count==5.and.state%fingerprint/=accepted_fingerprint,&
    'extension did not publish a new five-column state')
  call require(maxval(abs(state%coefficients(:,5)-candidate_vectors(:,5)))<1d-14,&
    'extension did not use the ordered local candidate')
  accepted_fingerprint=state%fingerprint

  call extend_dg_hybrid_schwarz_state(MPI_COMM_WORLD,7,6,mapping_fingerprint,candidate_fingerprint,&
    candidate_ids,candidate_vectors,state,ok,message,local_publish_ok=rank/=0)
  call require(.not.ok,'rank-local extension failure was accepted')
  call require(state%trial_count==5.and.state%fingerprint==accepted_fingerprint,&
    'failed extension changed the accepted state')

  call extend_dg_hybrid_schwarz_state(MPI_COMM_WORLD,8,6,mapping_fingerprint,candidate_fingerprint,&
    candidate_ids,candidate_vectors,state,ok,message)
  call require(.not.ok,'stale basis generation was accepted')
  call require(state%trial_count==5.and.state%fingerprint==accepted_fingerprint,&
    'stale extension changed the accepted state')

  call extend_dg_hybrid_schwarz_state(MPI_COMM_WORLD,7,99,mapping_fingerprint,candidate_fingerprint,&
    candidate_ids,candidate_vectors,state,ok,message)
  call require(.not.ok,'candidate capacity exhaustion was accepted')
  call require(state%trial_count==5.and.state%fingerprint==accepted_fingerprint,&
    'capacity failure changed the accepted state')

  candidate_ids(2)=candidate_ids(1)
  call initialize_dg_hybrid_schwarz_state(MPI_COMM_WORLD,fragment,nproc,7,4d0,300d0,2d0,1,&
    1d-12,1d-10,mapping_fingerprint,candidate_fingerprint,candidate_ids,candidate_energies,&
    candidate_vectors,state,ok,message)
  call require(.not.ok,'duplicate local candidate IDs were accepted')
  call require(state%trial_count==5.and.state%fingerprint==accepted_fingerprint,&
    'failed reinitialization changed the accepted state')
  candidate_ids(2)=int(100*fragment+2,int64)

  call initialize_dg_hybrid_schwarz_state(MPI_COMM_WORLD,modulo(fragment,nproc)+1,nproc,7,&
    4d0,300d0,2d0,1,1d-12,1d-10,mapping_fingerprint,candidate_fingerprint,candidate_ids,&
    candidate_energies,candidate_vectors,failed_state,ok,message)
  call require(.not.ok,'changed rank-fragment mapping was accepted')
  call require(.not.failed_state%valid,'failed mapping published a state')

  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid Schwarz state on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
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
end program test_dg_hybrid_schwarz_state_mpi
