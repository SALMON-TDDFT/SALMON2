#include "config.h"
program test_rt_dg_hybrid_initialization_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_checkpoint,only:s_rt_dg_hybrid_ground_state_payload,&
    write_rt_dg_hybrid_ground_state_checkpoint
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,initialize_rt_dg_hybrid_from_checkpoint,&
    fingerprint_rt_dg_hybrid_scope
  use rt_dg_hybrid_density_update,only:update_rt_dg_hybrid_density
  implicit none
  integer::comm,rank,nproc,ierr,smoke_grid_count
  integer(int64)::fingerprint
  logical::ok
  character(256)::message,path,mode,grid_count_argument
  type(s_rt_dg_hybrid_ground_state_payload)::payload
  type(s_rt_dg_hybrid_state)::state
  integer(int64)::structure_before,value_before
  complex(real64),allocatable::hamiltonian_before(:)
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_command_argument(1,path);call get_command_argument(2,mode);if(len_trim(mode)==0)mode='roundtrip'
  smoke_grid_count=2
  call get_command_argument(3,grid_count_argument)
  if(len_trim(grid_count_argument)>0)read(grid_count_argument,*)smoke_grid_count
  if(trim(mode)/='read_only')then
    call build_payload(payload)
    call write_rt_dg_hybrid_ground_state_checkpoint(comm,trim(path),payload,fingerprint,ok,message)
    call require(ok,'fixture write failed: '//trim(message))
    if(trim(mode)=='write_only'.or.trim(mode)=='write_production')then
      if(rank==0)write(*,'(a,i0)')'HYBRID_RT_INIT_FINGERPRINT=',fingerprint
      call MPI_Finalize(ierr);stop
    endif
  endif
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(path),'tddft_response',.true.,1,.false.,.false.,.false.,&
    .false.,.false.,[1,0,0],state,ok,message)
  call require(ok,trim(message))
  if(trim(mode)/='read_only')call require(state%payload_fingerprint==fingerprint,&
    'RT did not retain exact serialized payload identity')
  call require(state%initial_invariants_valid,'RT startup invariant receipt is absent')
  if(trim(mode)/='read_only')call require(state%operator_structure_fingerprint==payload%operator_structure_fingerprint,&
    'RT changed the operator-union identity during redistribution')
  structure_before=state%operator_structure_fingerprint;value_before=state%operator_value_fingerprint
  allocate(hamiltonian_before,source=state%operators%hamiltonian_values)
  call update_rt_dg_hybrid_density(comm,state,state%density,project_density_local,ok,message)
  call require(ok,trim(message))
  call require(all(state%operators%hamiltonian_values==hamiltonian_before),&
    'stored initial density did not reconstruct exact H_DG(0)')
  state%density=state%density+0.125d0
  call update_rt_dg_hybrid_density(comm,state,state%density,project_density_local,ok,message)
  call require(ok,trim(message))
  call require(any_rank(state%operators%hamiltonian_values/=hamiltonian_before),&
    'density perturbation did not update the RT Hamiltonian')
  call require(state%operator_structure_fingerprint==structure_before.and.&
    state%operators%fingerprint==structure_before.and.state%operator_value_fingerprint/=value_before,&
    'density update rebuilt the operator graph or retained a stale value fingerprint')
  fingerprint=state%payload_fingerprint
  call initialize_rt_dg_hybrid_from_checkpoint(comm,trim(path),'tddft_response',.true.,2,.false.,.false.,.false.,&
    .false.,.false.,[1],state,ok,message)
  call require(.not.ok,'spinful RT scope was accepted')
  if(rank==0)then
    write(*,'(a,i0)')'HYBRID_RT_INIT_FINGERPRINT=',fingerprint
    write(*,'(a,i0,a)')'PASS hybrid RT initialization on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine build_payload(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::p
    integer::i,row,nowned
    nowned=count([(mod(row-1,nproc)==rank,row=1,2)])
    p%valid=.true.;p%final_refresh_complete=.true.;p%analysis_complete=.true.;p%identity_only=.true.
    p%global_count=2;p%global_grid_count=smoke_grid_count;p%noccupied=1;p%operation_count=1;p%nonidentity_operation_count=0
    p%position_convention_fingerprint=115_int64
    p%catalog_fingerprint=101;p%state_fingerprint=102;p%metric_fingerprint=103
    p%operator_structure_fingerprint=104;p%operator_value_fingerprint=105
    p%kinetic_fingerprint=1;p%nonlocal_fingerprint=1;p%local_fingerprint=1;p%sipg_fingerprint=1
    p%basis_fingerprint=106;p%face_fingerprint=107;p%dc_seed_fingerprint=108
    p%continuation_fingerprint=109;p%analysis_fingerprint=111
    p%selection_fingerprint=112;p%pseudopotential_fingerprint=113;p%energy_fingerprint=114
    allocate(p%position_rows(3,nowned,2),p%symmetry_representation(2,2,1));p%position_rows=(0d0,0d0)
    p%symmetry_representation=(0d0,0d0);p%symmetry_representation(1,1,1)=(1d0,0d0)
    p%symmetry_representation(2,2,1)=(1d0,0d0)
    allocate(p%row_ids(nowned),p%metric_rows(nowned,2),p%kinetic_rows(nowned,2),p%nonlocal_rows(nowned,2),&
      p%local_rows(nowned,2),p%sipg_rows(nowned,2),p%hamiltonian_rows(nowned,2),p%coefficients(nowned,1))
    i=0
    do row=1,2
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;p%row_ids(i)=row;p%metric_rows(i,:)=(0d0,0d0);p%metric_rows(i,row)=(1d0,0d0)
      p%kinetic_rows(i,:)=(0d0,0d0);p%kinetic_rows(i,row)=cmplx(real(row,8)-0.5d0,0d0,8)
      p%nonlocal_rows(i,:)=(0d0,0d0);p%local_rows(i,:)=(0d0,0d0);p%local_rows(i,row)=(0.5d0,0d0)
      p%sipg_rows(i,:)=(0d0,0d0)
      p%hamiltonian_rows(i,:)=p%kinetic_rows(i,:)+p%local_rows(i,:)
      p%coefficients(i,1)=merge((1d0,0d0),(0d0,0d0),row==1)
    enddo
    if(trim(mode)=='write_production')then
      p%kinetic_rows=(0d0,0d0);p%local_rows=(0d0,0d0);p%hamiltonian_rows=(0d0,0d0)
    endif
    call complete_payload(p)
    call component_fingerprints(p)
  end subroutine build_payload
  subroutine complete_payload(p)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer::nrow,i,point,npoint
    nrow=size(p%row_ids)
    allocate(p%metric_row_offsets(nrow+1),p%metric_column_ids(2*nrow),p%operator_row_offsets(nrow+1),&
      p%operator_column_ids(2*nrow));p%metric_row_offsets=[(2*i-1,i=1,nrow+1)]
    p%operator_row_offsets=p%metric_row_offsets
    do i=1,nrow;p%metric_column_ids(2*i-1:2*i)=[1,2];p%operator_column_ids(2*i-1:2*i)=[1,2];enddo
    npoint=count([(mod(point-1,nproc)==rank,point=1,smoke_grid_count)])
    allocate(p%grid_ids(npoint),p%grid_weights(npoint),p%partition_ids(npoint),p%basis_values(2,npoint),p%density(npoint))
    i=0
    do point=1,smoke_grid_count
      if(mod(point-1,nproc)/=rank)cycle
      i=i+1;p%grid_ids(i)=point
    enddo
    p%grid_weights=1d0;p%partition_ids=1;p%basis_values=(0d0,0d0);p%density=1d0/real(smoke_grid_count,8)
    allocate(p%occupations(1),p%eigenvalues(1));p%occupations=1d0;p%eigenvalues=merge(0d0,1d0,trim(mode)=='write_production')
    allocate(p%requested_ids(2),p%effective_ids(2),p%added_ids(0),p%closure_parent(0),p%closure_reason(0),p%closure_action(0))
    p%requested_ids=[1,2];p%effective_ids=[1,2]
    allocate(p%scope_selectors(8),p%xc_types(3));p%scope_selectors=[1,1,1,0,0,0,0,0];p%xc_types=[1,0,0]
    p%scope_fingerprint=fingerprint_rt_dg_hybrid_scope(p%scope_selectors,p%xc_types)
    allocate(p%continuation_receipt(1),p%pseudopotential_receipt(1),p%energy_receipt(1));p%continuation_receipt=0d0
    p%pseudopotential_receipt=0d0;p%energy_receipt=0d0
    allocate(p%face_ids(0),p%face_metadata(8,0),p%face_normals(3,0),p%face_offsets(1),p%face_value_offsets(1),&
      p%face_basis_ids(0),p%face_point_ids(0),p%face_weights(0),p%face_values(1,0));p%face_offsets=1;p%face_value_offsets=1
    allocate(p%nonlocal_ids(nrow),p%nonlocal_owner(nrow),p%nonlocal_values(1,nrow));p%nonlocal_ids=p%row_ids
    p%nonlocal_owner=rank;p%nonlocal_values=(0d0,0d0)
    allocate(p%interface_observables(0,0))
  end subroutine complete_payload
  subroutine component_fingerprints(p)
    use rt_dg_hybrid_checkpoint,only:fingerprint_rt_dg_hybrid_component
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%kinetic_rows,p%kinetic_fingerprint,ok)
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%nonlocal_rows,p%nonlocal_fingerprint,ok)
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%local_rows,p%local_fingerprint,ok)
    call fingerprint_rt_dg_hybrid_component(comm,p%row_ids,p%sipg_rows,p%sipg_fingerprint,ok)
  end subroutine component_fingerprints
  subroutine project_density_local(row_ids,grid_ids,density,local_rows,callback_ok,callback_message)
    integer(int64),intent(in)::row_ids(:),grid_ids(:)
    real(real64),intent(in)::density(:)
    complex(real64),intent(out)::local_rows(:,:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    integer::q;real(real64)::local_sum,global_sum;integer::local_count,global_count
    local_sum=sum(density);local_count=size(grid_ids)
    call MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    local_rows=(0d0,0d0)
    do q=1,size(row_ids)
      local_rows(q,int(row_ids(q)))=cmplx(global_sum/real(global_count,real64),0d0,real64)
      local_rows(q,3-int(row_ids(q)))=cmplx(global_sum/real(global_count,real64)-0.5d0,0d0,real64)
    enddo
    callback_ok=.true.;callback_message=''
  end subroutine project_density_local
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
  logical function any_rank(values)
    logical,intent(in)::values(:);integer::local_value,global_value
    local_value=merge(1,0,any(values));call MPI_Allreduce(local_value,global_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    any_rank=global_value==1
  end function any_rank
end program test_rt_dg_hybrid_initialization_mpi
