#include "config.h"
program test_dg_hybrid_continuation_state_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_continuation_state,only:s_dg_hybrid_catalog_receipt,s_dg_hybrid_scope_receipt,&
    s_dg_hybrid_continuation_state,initialize_dg_hybrid_continuation,close_dg_hybrid_selection,&
    build_dg_hybrid_scope_receipt,validate_dg_hybrid_frozen_catalog
  implicit none
  integer,parameter::ngrid=6
  integer::comm,rank,nproc,ierr,i
  integer(int64)::fingerprint,reference_fingerprint
  integer(int64),allocatable::core_ids(:)
  integer,allocatable::effective(:),added_parent(:),added_operation(:)
  integer::universe(4),requested(1),group_action(4,2),identity_action(4,1)
  real(real64)::dc_density(ngrid),trial_density(ngrid)
  type(s_dg_hybrid_catalog_receipt)::catalog
  type(s_dg_hybrid_scope_receipt)::scope
  type(s_dg_hybrid_continuation_state)::state
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  dc_density=[(0.25d0+0.01d0*i,i=1,ngrid)]
  trial_density=dc_density+0.5d0
  allocate(core_ids(count([(mod(i-1,nproc)==rank,i=1,ngrid)])))
  core_ids=pack([(int(i,int64),i=1,ngrid)],[(mod(i-1,nproc)==rank,i=1,ngrid)])

  call valid_identity_catalog(catalog)
  call build_dg_hybrid_scope_receipt(comm,1,.true.,1,.false.,.false.,.false.,.false.,.false.,[1],&
    scope,ok,message)
  call require(ok,trim(message))
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(ok,trim(message))
  call require(all(state%seed_density==dc_density),'continuation seed is not the exact DC density')
  call require(all(state%mixed_density==dc_density),'mixed density did not start from the DC density')
  call require(any(state%seed_density/=trial_density),'trial coefficient density was used as the seed')
  call require(state%lambda==0d0.and..not.state%lambda_accepted,'lambda zero was accepted before fixed-point gates')
  call require(.not.state%accepted%valid,'accepted snapshot exists before lambda-zero convergence')
  call require(state%seed_fingerprint/=0_int64.and.fingerprint/=0_int64,'missing continuation fingerprints')
  reference_fingerprint=fingerprint
  call require(state%catalog%identity_only.and.state%catalog%nonidentity_count==0,&
    'identity-only symmetry receipt was not preserved')
  call validate_dg_hybrid_frozen_catalog(comm,state,catalog,ok,message)
  call require(ok,trim(message))
  catalog%selection_fingerprint=catalog%selection_fingerprint+1_int64
  call validate_dg_hybrid_frozen_catalog(comm,state,catalog,ok,message)
  call require(.not.ok,'catalog mutation after initialization was accepted')
  call valid_identity_catalog(catalog)

  universe=[10,20,30,40];requested=[10]
  group_action(:,1)=[1,2,3,4];group_action(:,2)=[2,1,4,3]
  call close_dg_hybrid_selection(comm,requested,universe,group_action,effective,added_parent,&
    added_operation,fingerprint,ok,message)
  call require(ok,trim(message))
  call require(all(effective==[10,20]),'selection did not expand to its complete symmetry orbit')
  call require(size(added_parent)==1.and.added_parent(1)==10.and.added_operation(1)==2,&
    'selection closure provenance is incorrect')
  identity_action(:,1)=[1,2,3,4];requested=[30]
  call close_dg_hybrid_selection(comm,requested,universe,identity_action,effective,added_parent,&
    added_operation,fingerprint,ok,message)
  call require(ok.and.all(effective==[30]).and.size(added_parent)==0,&
    'identity-only selection changed the requested catalog')
  group_action(1,2)=3
  call close_dg_hybrid_selection(comm,requested,universe,group_action,effective,added_parent,&
    added_operation,fingerprint,ok,message)
  call require(.not.ok,'malformed symmetry action was accepted')
  group_action(:,2)=[2,1,4,3]
  call close_dg_hybrid_selection(comm,requested,universe,group_action(:,2:2),effective,added_parent,&
    added_operation,fingerprint,ok,message)
  call require(.not.ok,'symmetry action without identity was accepted')

  catalog%analysis_complete=.false.
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(.not.ok,'unfinished symmetry analysis was accepted')
  call valid_identity_catalog(catalog);catalog%analysis_fingerprint=0_int64
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(.not.ok,'missing symmetry provenance was accepted')
  call valid_identity_catalog(catalog);catalog%operation_count=0
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(.not.ok,'empty symmetry operation list was accepted')
  call valid_identity_catalog(catalog);catalog%catalog_fingerprint=0_int64
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(.not.ok,'zero catalog fingerprint was accepted')

  call valid_identity_catalog(catalog)
  if(rank==0.and.size(core_ids)>0)core_ids=[core_ids,core_ids(1)]
  call initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
    state,fingerprint,ok,message)
  call require(.not.ok,'duplicate distributed core ownership was accepted')

  call build_dg_hybrid_scope_receipt(comm,1,.true.,2,.false.,.false.,.false.,.false.,.false.,[1],&
    scope,ok,message)
  call require(.not.ok,'spin-polarized continuation scope was accepted')
  call build_dg_hybrid_scope_receipt(comm,1,.false.,1,.false.,.false.,.false.,.false.,.false.,[1],&
    scope,ok,message)
  call require(.not.ok,'nonperiodic continuation scope was accepted')
  call build_dg_hybrid_scope_receipt(comm,1,.true.,1,.false.,.false.,.false.,.false.,.false.,[99],&
    scope,ok,message)
  call require(.not.ok,'unsupported XC continuation scope was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'CONTINUATION_STATE ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid continuation state on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine valid_identity_catalog(receipt)
    type(s_dg_hybrid_catalog_receipt),intent(out)::receipt
    receipt%frozen=.true.
    receipt%analysis_complete=.true.
    receipt%identity_only=.true.
    receipt%operation_count=1
    receipt%nonidentity_count=0
    receipt%analysis_fingerprint=101_int64
    receipt%catalog_fingerprint=202_int64
    receipt%selection_fingerprint=303_int64
  end subroutine valid_identity_catalog

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_continuation_state_mpi
