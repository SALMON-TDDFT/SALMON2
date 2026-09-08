#include "config.h"
program test_dg_hybrid_continuation_acceptance_mpi
  use mpi, only: MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_Finalize, MPI_Init, &
    MPI_INTEGER, MPI_MAX, MPI_SUCCESS
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_continuation_acceptance,only:s_dg_hybrid_acceptance_result,evaluate_dg_hybrid_acceptance,&
    validate_dg_hybrid_acceptance_receipt
  implicit none
  integer::icomm,id_rank,nproc,ierr
  complex(real64),allocatable::d(:,:,:),dinv(:,:,:)
  complex(real64)::hvol(2,2),hface(2,2),metric(2,2),qret(2,2),gamma(2,2)
  complex(real64)::hgrid(2),sgrid_epsilon(2)
  complex(real64)::p(2,2),pinv(2,2)
  real(real64)::occupations(2),face_lambda(2),face_terms(3,2),weights(2)
  integer::clusters(2)
  type(s_dg_hybrid_acceptance_result)::result
  logical::ok
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  allocate(d(2,2,2),dinv(2,2,2))
  call passing_payload()
  call evaluate(ok,message)
  call require(ok.and..not.result%identity_only,trim(message))
  call require(.not.result%excitation_cutoff_convergence_proven,&
    'symmetry-complete retained space was mislabeled as excitation-cutoff convergence')
  call valid_receipt(result)
  call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2,1001_int64,96_int64,ok,message)
  call require(ok,'candidate-bound acceptance receipt was rejected: '//trim(message))
  result%checked_occupied_count=1
  call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2,1001_int64,96_int64,ok,message)
  call require(.not.ok,'wrong occupied count in acceptance receipt was accepted')
  call valid_receipt(result);result%operator_fingerprint=1000_int64
  call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2,1001_int64,96_int64,ok,message)
  call require(.not.ok,'stale operator fingerprint in acceptance receipt was accepted')
  call valid_receipt(result);result%checked_face_count=1
  call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2,1001_int64,96_int64,ok,message)
  call require(.not.ok,'omitted face in acceptance receipt was accepted')
  call valid_receipt(result);result%face_topology_fingerprint=0_int64
  call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2,1001_int64,0_int64,ok,message)
  call require(.not.ok,'zero face-topology fingerprint was accepted')
  if(nproc>1)then
    call valid_receipt(result);result%checked_face_count=2+id_rank
    call validate_dg_hybrid_acceptance_receipt(icomm,result,2,2+id_rank,1001_int64,96_int64,ok,message)
    call require(.not.ok,'rank-disagreeing face count was accepted')
  endif
  call nonorthogonal_payload();call evaluate(ok,message)
  call require(ok,'nonorthogonal bilinear/congruence covariance convention failed: '//trim(message))

  d(:,:,2)=reshape([cmplx(0d0,0d0,real64),cmplx(0d0,0d0,real64),&
    cmplx(1d0,0d0,real64),cmplx(1d0,0d0,real64)],[2,2])
  call evaluate(ok,message);call require(.not.ok,'basis missing a symmetry partner was accepted')
  call passing_payload();qret=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),&
    cmplx(0d0,0d0,real64),cmplx(0d0,0d0,real64)],[2,2])
  call evaluate(ok,message);call require(.not.ok,'cutoff splitting a symmetry multiplet was accepted')
  call passing_payload();face_lambda=[0.25d0,0.5d0]
  call evaluate(ok,message);call require(.not.ok,'face-local continuation lambda was accepted')
  call passing_payload();occupations=[1d0,0.5d0]
  gamma=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),&
    cmplx(0d0,0d0,real64),cmplx(0.5d0,0d0,real64)],[2,2])
  call evaluate(ok,message);call require(.not.ok,'symmetry-incompatible degenerate occupations were accepted')
  call passing_payload();hgrid(2)=1d-2
  call evaluate(ok,message);call require(.not.ok,'large reconstructed-grid DG residual was accepted')
  if(nproc>1)then
    call passing_payload();face_lambda=0.5d0+0.1d0*id_rank
    call evaluate(ok,message);call require(.not.ok,'rank-disagreeing continuation lambda was accepted')
    call passing_payload()
    call evaluate_dg_hybrid_acceptance(icomm,.true.,91_int64+id_rank,2,1,.false.,d,dinv,hvol,hface,metric,qret,gamma,&
      occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,ok,message)
    call require(.not.ok,'rank-disagreeing symmetry provenance was accepted')
    call passing_payload();if(id_rank==0)metric(2,2)=-1d0
    call evaluate(ok,message);call require(.not.ok,'rank-local indefinite DG metric was accepted')
  endif

  deallocate(d,dinv);allocate(d(2,2,1),dinv(2,2,1));call passing_payload()
  call evaluate_dg_hybrid_acceptance(icomm,.true.,91_int64,2,1,.false.,d,dinv,hvol,hface,metric,qret,gamma,&
    occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,ok,message)
  call require(.not.ok,'omitted nonidentity symmetry operation was accepted')
  call evaluate_dg_hybrid_acceptance(icomm,.true.,91_int64,1,0,.true.,d,dinv,hvol,hface,metric,qret,gamma,&
    occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,ok,message)
  call require(ok.and.result%identity_only,'identity-only acceptance path failed: '//trim(message))
  call evaluate_dg_hybrid_acceptance(icomm,.false.,91_int64,1,0,.true.,d,dinv,hvol,hface,metric,qret,gamma,&
    occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,ok,message)
  call require(.not.ok,'unfinished authoritative symmetry analysis was accepted')
  call evaluate_dg_hybrid_acceptance(icomm,.true.,0_int64,1,0,.true.,d,dinv,hvol,hface,metric,qret,gamma,&
    occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,ok,message)
  call require(.not.ok,'missing authoritative symmetry provenance was accepted')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS hybrid continuation acceptance on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine passing_payload()
    d=(0d0,0d0);dinv=(0d0,0d0)
    d(1,1,1)=1d0;d(2,2,1)=1d0;dinv(:,:,1)=d(:,:,1)
    if(size(d,3)>1)then;d(1,2,2)=1d0;d(2,1,2)=1d0;dinv(:,:,2)=d(:,:,2);endif
    hvol=reshape([cmplx(2d0,0d0,real64),cmplx(0.2d0,0d0,real64),&
      cmplx(0.2d0,0d0,real64),cmplx(2d0,0d0,real64)],[2,2])
    hface=reshape([cmplx(0.3d0,0d0,real64),cmplx(-0.1d0,0d0,real64),&
      cmplx(-0.1d0,0d0,real64),cmplx(0.3d0,0d0,real64)],[2,2])
    metric=0d0;metric(1,1)=1d0;metric(2,2)=1d0;qret=metric;gamma=metric
    occupations=[1d0,1d0];clusters=[1,1];face_lambda=0.5d0
    hgrid=[cmplx(1d0,0d0,real64),cmplx(2d0,0d0,real64)];sgrid_epsilon=hgrid
    face_terms=0d0;weights=1d0
  end subroutine passing_payload
  subroutine nonorthogonal_payload()
    call passing_payload();p=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),&
      cmplx(0.3d0,0d0,real64),cmplx(1d0,0d0,real64)],[2,2])
    pinv=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),&
      cmplx(-0.3d0,0d0,real64),cmplx(1d0,0d0,real64)],[2,2])
    d(:,:,2)=matmul(p,matmul(d(:,:,2),pinv));dinv(:,:,2)=d(:,:,2)
    metric=matmul(conjg(transpose(pinv)),pinv);gamma=matmul(p,conjg(transpose(p)))
    hvol=2d0*metric;hface=0.3d0*metric
  end subroutine nonorthogonal_payload
  subroutine valid_receipt(receipt)
    type(s_dg_hybrid_acceptance_result),intent(out)::receipt
    receipt=s_dg_hybrid_acceptance_result();receipt%valid=.true.;receipt%symmetry_complete=.true.
    receipt%grid_complete=.true.;receipt%face_complete=.true.;receipt%expected_occupied_count=2
    receipt%checked_occupied_count=2;receipt%checked_face_count=2;receipt%analysis_fingerprint=91_int64
    receipt%basis_fingerprint=92_int64;receipt%operator_fingerprint=1001_int64
    receipt%state_fingerprint=94_int64;receipt%action_fingerprint=95_int64
    receipt%action_operator_fingerprint=1001_int64;receipt%face_topology_fingerprint=96_int64
  end subroutine valid_receipt
  subroutine evaluate(passed,detail)
    logical,intent(out)::passed;character(*),intent(out)::detail
    call evaluate_dg_hybrid_acceptance(icomm,.true.,91_int64,2,1,.false.,d,dinv,hvol,hface,metric,qret,gamma,&
      occupations,clusters,face_lambda,hgrid,sgrid_epsilon,face_terms,weights,weights,1d-10,result,passed,detail)
  end subroutine evaluate
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(id_rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_continuation_acceptance_mpi
