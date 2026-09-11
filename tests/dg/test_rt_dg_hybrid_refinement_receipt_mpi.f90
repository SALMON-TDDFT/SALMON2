program test_rt_dg_hybrid_refinement_receipt_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use mpi
  use rt_dg_hybrid_refinement_receipt,only:s_rt_dg_hybrid_refinement_receipt,&
    write_rt_dg_hybrid_refinement_receipt,read_rt_dg_hybrid_refinement_receipt
  implicit none
  integer::ierr,rank,nproc
  character(32)::mode
  character(512)::prefix
  character(256)::message
  logical::ok,present
  type(s_rt_dg_hybrid_refinement_receipt)::written,loaded
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call get_command_argument(1,mode);call get_command_argument(2,prefix);call fill_receipt(written)
  if(trim(mode)=='read-corrupt')then
    call read_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix),9001_int64,loaded,present,ok,message)
    call require(.not.ok.and.present,'corrupt companion receipt was accepted')
    if(rank==0)write(*,'(a,i0)')'PASS corrupt refinement receipt rejection ranks=',nproc
    call MPI_Finalize(ierr);stop
  endif
  call read_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix)//'-legacy',9001_int64,&
    loaded,present,ok,message)
  call require(ok.and..not.present,'missing legacy companion was not accepted')
  call write_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix),written,ok,message)
  call require(ok,'valid refinement receipt write failed: '//trim(message))
  call read_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix),9001_int64,loaded,present,ok,message)
  call require(ok.and.present,'valid refinement receipt read failed: '//trim(message))
  call require(loaded%version==1.and.loaded%v5_publication_fingerprint==9001_int64,'wrong receipt identity')
  call require(loaded%total_solve_count==3.and.loaded%additional_refinement_count==2,'wrong solve counts')
  call require(loaded%density_change==1.25e-6_real64.and.loaded%energy_change==2.5e-7_real64,'wrong metrics')
  call require(loaded%converged.and..not.loaded%exhausted.and.trim(loaded%exit_reason)=='refined-lcfo-converged',&
    'wrong terminal status')
  call require(any(loaded%digest/=0_int64),'missing full refinement SHA-256')
  call read_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix),9002_int64,loaded,present,ok,message)
  call require(.not.ok.and.present,'wrong v5 binding was accepted')
  if(nproc>1.and.rank==nproc-1)written%density_change=1.5e-6_real64
  call write_rt_dg_hybrid_refinement_receipt(MPI_COMM_WORLD,trim(prefix)//'-disagree',written,ok,message)
  if(nproc>1)then;call require(.not.ok,'rank-disagreeing metrics were accepted')
  else;call require(ok,'single-rank refinement receipt was rejected');endif
  if(rank==0)write(*,'(a,i0)')'PASS refinement receipt round trip ranks=',nproc
  call MPI_Finalize(ierr)
contains
  subroutine fill_receipt(receipt)
    type(s_rt_dg_hybrid_refinement_receipt),intent(out)::receipt
    receipt=s_rt_dg_hybrid_refinement_receipt();receipt%version=1
    receipt%v5_publication_fingerprint=9001_int64;receipt%fragment_id=rank+1
    receipt%total_solve_count=3;receipt%additional_refinement_count=2
    receipt%density_change=1.25e-6_real64;receipt%energy_change=2.5e-7_real64
    receipt%converged=.true.;receipt%exhausted=.false.;receipt%exit_reason='refined-lcfo-converged'
  end subroutine
  subroutine require(condition,why)
    logical,intent(in)::condition;character(*),intent(in)::why
    if(.not.condition)then;write(0,'(a,i0,2a)')'rank ',rank,': ',trim(why);call MPI_Abort(MPI_COMM_WORLD,1,ierr);endif
  end subroutine
end program
