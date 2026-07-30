#include "config.h"
program test_dg_overlapping_wannier_scf_mpi
  use mpi
  use iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_scf
  implicit none
  integer::comm,rank,nproc,ierr,nlocal,i,rebuilds,corrupt_outer
  integer(int64),allocatable::row_ids(:),point_ids(:)
  integer(int64)::basis_fingerprint,rounded_fingerprint,changed_fingerprint
  complex(8),allocatable::srows(:,:),values(:,:)
  integer,allocatable::generations(:,:)
  real(8),allocatable::weights(:)
  real(8),allocatable::mix_current(:),mix_raw(:),mix_history(:,:),mix_result(:),mix_new_history(:,:)
  type(s_dg_overlapping_wannier_scf_state)::state,before
  type(s_dg_overlapping_wannier_scf_result)::result
  logical::ok,force_failure,force_rollback_failure,seed_used
  real(8)::external_field,external_snapshot,saved_density
  complex(8)::saved_s
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(i-1,nproc)==rank,i=1,2)])
  allocate(row_ids(nlocal),srows(nlocal,2),point_ids(nlocal),weights(nlocal),values(2,nlocal),generations(2,nlocal))
  nlocal=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;row_ids(nlocal)=i;point_ids(nlocal)=i;weights(nlocal)=1d0
    srows(nlocal,:)=(0d0,0d0);srows(nlocal,i)=1d0
    values(:,nlocal)=(0d0,0d0);values(i,nlocal)=1d0;generations(:,nlocal)=7
  enddo
  allocate(mix_current(nlocal),mix_raw(nlocal),mix_history(nlocal,2),mix_result(nlocal),&
    mix_new_history(nlocal,2))
  mix_current=0.4d0;mix_raw=1d0;mix_history=0d0
  call mix_dg_overlapping_wannier_density_history(comm,0.2d0,mix_current,mix_raw,mix_history,1,&
    mix_result,mix_new_history,i,ok,message)
  call require(ok.and.all(mix_result>0.4d0).and.all(mix_result<0.52d0),&
    'seeded Anderson mixing uses a positivity-preserving convex update')
  call compute_dg_overlapping_wannier_scf_fingerprint(comm,row_ids,srows,point_ids,weights,values,&
    generations,basis_fingerprint,77_int64)
  if(rank==0)then
    saved_s=srows(1,1)
    srows(1,1)=saved_s+cmplx(1d-13,-1d-13,8)
  endif
  call compute_dg_overlapping_wannier_scf_fingerprint(comm,row_ids,srows,point_ids,weights,values,&
    generations,rounded_fingerprint,77_int64)
  call require(rounded_fingerprint==basis_fingerprint,&
    'basis fingerprint ignores sub-tolerance floating-point roundoff')
  if(rank==0)srows(1,1)=saved_s+cmplx(1d-6,0d0,8)
  call compute_dg_overlapping_wannier_scf_fingerprint(comm,row_ids,srows,point_ids,weights,values,&
    generations,changed_fingerprint,77_int64)
  call require(changed_fingerprint/=basis_fingerprint,&
    'basis fingerprint detects physically meaningful basis changes')
  if(rank==0)srows(1,1)=saved_s
  allocate(state%density(nlocal),state%potential(nlocal),state%coefficients(2,1),state%eigenvalues(1),&
    state%density_history(nlocal,3))
  do i=1,nlocal;state%density(i)=merge(0.25d0,0.75d0,point_ids(i)==1_int64);enddo
  state%potential=state%density
  state%coefficients=reshape([cmplx(0d0,0d0,8),cmplx(1d0,0d0,8)],[2,1]);state%eigenvalues=2d0
  state%density_history=0d0;state%density_history(:,1)=state%density;state%history_count=1
  state%basis_generation=7;state%geometry_generation=3
  state%basis_fingerprint=basis_fingerprint;state%operator_fingerprint=900_int64
  rebuilds=0;force_failure=.false.;force_rollback_failure=.false.;seed_used=.false.;external_field=5d0
  if(rank==0)then
    saved_density=state%density(1)
    state%density(1)=ieee_value(0d0,ieee_quiet_nan)
  endif
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,40,100,1d-9,1d-10,build_toy,mix_toy,transaction_toy,&
    commit_toy,state,result,ok,message)
  call require(.not.ok.and.trim(message)=='nonfinite overlapping-Wannier SCF density code=116',&
    'nonfinite SCF state identifies the contaminated field collectively')
  if(rank==0)state%density(1)=saved_density
  before=state
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,1,100,1d-30,1d-10,&
    build_toy,mix_toy,transaction_toy,commit_toy,state,result,ok,message)
  call require(.not.ok.and.index(message,'did not converge')>0,&
    'outer-iteration exhaustion is an explicit failure: '//trim(message))
  call require(all(state%density==before%density),'outer-iteration exhaustion rolls back state')
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,40,100,1d-9,1d-10,build_toy,mix_toy,transaction_toy,&
    commit_toy,state,result,ok,message)
  call require(ok.and.result%converged,trim(message))
  call require(result%unmixed_density_residual<1d-9,'unmixed fixed point')
  call require(result%hamiltonian_rebuilds==result%iterations+1,'H rebuild only at outer boundaries')
  call require(state%history_count>1,'seeded mixing history retained')
  call require(seed_used,'seeded DC mixer history consumed')
  call require(abs(result%integrated_charge-1d0)<1d-10,'core-only density charge')
  before=state;force_failure=.true.;rebuilds=0
  external_field=17d0
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,4,100,1d-9,1d-10,build_toy,mix_toy,transaction_toy,&
    commit_toy,state,result,ok,message)
  call require(.not.ok,'callback failure is explicit')
  call require(all(state%density==before%density).and.&
    all(state%coefficients==before%coefficients).and.&
    state%operator_fingerprint==before%operator_fingerprint.and.&
    state%history_count==before%history_count,'atomic rollback')
  call require(external_field==17d0,'external DC field rollback')
  force_rollback_failure=.true.
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,4,100,1d-9,1d-10,&
    build_toy,mix_toy,transaction_toy,commit_toy,state,result,ok,message)
  call require(.not.ok.and.index(message,'forced toy rebuild failure')>0.and.&
    index(message,'forced toy rollback failure')>0,&
    'rollback failure preserves both original and rollback diagnostics')
  force_rollback_failure=.false.
  force_failure=.false.;corrupt_outer=merge(5,4,rank==max(0,nproc-1))
  call run_dg_overlapping_wannier_scf(comm,row_ids,srows,point_ids,weights,values,generations,7,&
    3,basis_fingerprint,[1d0],1d-12,1d-9,77_int64,2_int64,0.5d0,corrupt_outer,100,1d-9,1d-10,&
    build_toy,mix_toy,transaction_toy,&
    commit_toy,state,result,ok,message)
  if(nproc>1)call require(.not.ok,'rank-inconsistent SCF contract rejected collectively')
  if(rank==0)write(*,'(a,i0,a,i0)')'PASS overlapping-Wannier SCF on ',nproc,&
    ' ranks fingerprint=',basis_fingerprint
  call MPI_Finalize(ierr)
contains
  subroutine build_toy(comm_in,density,potential,hrows,new_potential,fingerprint,callback_ok,callback_message)
    integer,intent(in)::comm_in
    real(8),intent(in)::density(:),potential(:)
    complex(8),intent(out)::hrows(:,:)
    real(8),intent(out)::new_potential(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    integer::j
    rebuilds=rebuilds+1;hrows=(0d0,0d0)
    external_field=sum(density)+100d0
    do j=1,size(row_ids);hrows(j,int(row_ids(j)))=merge(0d0,2d0,row_ids(j)==1_int64);enddo
    new_potential=density;fingerprint=900_int64+rebuilds
    callback_ok=.not.force_failure;callback_message=''
    if(force_failure)callback_message='forced toy rebuild failure'
  end subroutine
  subroutine commit_toy()
  end subroutine
  subroutine mix_toy(comm_in,mixing,current_density,raw_density,history,history_count,mixed_density,&
      new_history,new_history_count,mix_ok,mix_message)
    integer,intent(in)::comm_in,history_count
    real(8),intent(in)::mixing,current_density(:),raw_density(:),history(:,:)
    real(8),intent(out)::mixed_density(:),new_history(:,:)
    integer,intent(out)::new_history_count
    logical,intent(out)::mix_ok
    character(*),intent(out)::mix_message
    real(8)::rate
    seed_used=seed_used.or.history_count>0
    rate=min(1d0,mixing*(1d0+0.1d0*history_count))
    mixed_density=(1d0-rate)*current_density+rate*raw_density
    new_history=history;new_history_count=min(size(history,2),history_count+1)
    if(history_count>=size(history,2).and.size(history,2)>1)new_history(:,1:-1+size(history,2))=history(:,2:)
    new_history(:,new_history_count)=raw_density
    mix_ok=.true.;mix_message=''
  end subroutine
  subroutine transaction_toy(action,transaction_ok,transaction_message)
    integer,intent(in)::action
    logical,intent(out)::transaction_ok
    character(*),intent(out)::transaction_message
    if(action==0)then
      external_snapshot=external_field
    else if(action==-1)then
      external_field=external_snapshot
    endif
    transaction_ok=.not.(action==-1.and.force_rollback_failure)
    transaction_message=''
    if(.not.transaction_ok)transaction_message='forced toy rollback failure'
  end subroutine
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    if(.not.condition)then
      write(0,'(a)')trim(text);error stop 1
    endif
  end subroutine
end program
