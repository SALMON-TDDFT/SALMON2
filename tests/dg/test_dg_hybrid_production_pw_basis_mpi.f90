#include "config.h"
program test_dg_hybrid_production_pw_basis_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog,s_dg_hybrid_production_selection
  use dg_hybrid_windowed_pw_basis,only:materialize_dg_hybrid_windowed_pw_columns
  use dg_hybrid_production_pw_basis,only:build_dg_hybrid_production_pw_basis,&
    analyze_dg_hybrid_production_selection,analyze_dg_hybrid_lcfo_selection,&
    freeze_dg_hybrid_production_selection
  use dg_hybrid_continuation_state,only:close_dg_hybrid_selection
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer::comm,rank,nproc,ierr,nowned,i,p
  integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:),mixed_row_action(:,:),&
    coset_row_action(:,:),root_row_action(:,:),root_core_fragment_ids(:)
  integer(int64),allocatable::box_ids(:),core_ids(:),root_core_ids(:)
  real(real64),allocatable::box_windows(:,:),coordinates(:,:),windows(:,:),g_vectors(:,:),root_coordinates(:,:)
  real(real64)::reciprocal_lattice(3,3),reciprocal_rotation(3,3,2)
  real(real64)::mixed_reciprocal_rotation(3,3,4)
  complex(real64),allocatable::tile(:,:)
  type(s_dg_hybrid_basis_catalog)::catalog
  type(s_dg_hybrid_production_selection)::selection
  integer,allocatable::effective_ids(:),added_parent(:),added_operation(:)
  integer(int64)::fingerprint,workspace,closure_fingerprint
  logical::ok,values_ok
  character(256)::message
#ifdef USE_MPI
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif
  nowned=count([(mod(i-1,nproc)==rank,i=1,2)])
  allocate(fragment_ids(nowned),box_ids(4),box_windows(nowned,4),row_action(4,2))
  p=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;fragment_ids(p)=i
  enddo
  box_ids=[1_int64,2_int64,3_int64,4_int64]
  do p=1,nowned
    if(fragment_ids(p)==1)box_windows(p,:)=[4d0,3d0,2d0,1d0]
    if(fragment_ids(p)==2)box_windows(p,:)=[1d0,2d0,3d0,4d0]
  enddo
  row_action(:,1)=[1,2,3,4];row_action(:,2)=[4,3,2,1]
  if(rank==0)then
    allocate(core_ids(2),source=[1_int64,2_int64]);allocate(core_fragment_ids(2),source=[1,1])
  elseif(rank==1)then
    allocate(core_ids(2),source=[3_int64,4_int64]);allocate(core_fragment_ids(2),source=[2,2])
  else
    allocate(core_ids(0),core_fragment_ids(0))
  endif
  if(nproc==1)then
    deallocate(core_ids,core_fragment_ids);allocate(core_ids(4),source=box_ids)
    allocate(core_fragment_ids(4),source=[1,1,2,2])
  endif
  allocate(coordinates(3,size(core_ids)));coordinates=0d0
  do p=1,size(core_ids);coordinates(1,p)=real(core_ids(p)-1_int64,real64);enddo
  reciprocal_lattice=0d0;reciprocal_rotation=0d0
  reciprocal_lattice(1,1)=1d0
  reciprocal_lattice(2,2)=1d0+5d-14
  reciprocal_lattice(3,3)=1d0+2d-14
  do i=1,3
    reciprocal_rotation(i,i,1)=1d0;reciprocal_rotation(i,i,2)=-1d0
  enddo
  call build_dg_hybrid_production_pw_basis(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
    windows,g_vectors,catalog,workspace,fingerprint,ok,message)
  call require(ok,'production PW basis rejected: '//trim(message))
  call require(catalog%valid.and.size(catalog%packets)==2,'production packet catalog mismatch')
  call require(size(g_vectors,2)==1,'zero-cutoff production catalog must contain only G=0')
  allocate(tile(1,size(core_ids)))
  call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,coordinates,windows,1,1,tile,ok,message)
  call require(ok,'production PW packet materialization failed: '//trim(message))
  values_ok=.true.
  do p=1,size(core_ids)
    values_ok=values_ok.and.abs(tile(1,p)-cmplx(windows(1,p),0d0,real64))<1d-12
  enddo
  call require(values_ok,'materialized production PW values mismatch')
  call require(fingerprint/=0_int64.and.workspace>0_int64,'production PW receipts are missing')
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
    windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok,'authoritative production symmetry analysis rejected: '//trim(message))
  call require(selection%analysis_complete.and..not.selection%identity_only,&
    'nontrivial production symmetry receipt is incomplete')
  call require(selection%operation_count==2.and.selection%window_operation_count==2.and.&
    selection%analysis_fingerprint/=0_int64,&
    'production symmetry receipt metadata mismatch')
  call require(all(shape(selection%packet_action)==[2,2]),'production packet action shape mismatch')
  call require(all(selection%packet_action(:,1)==[1,2]).and.all(selection%packet_action(:,2)==[2,1]),&
    'production packet action mismatch')
  call close_dg_hybrid_selection(comm,[1],selection%packet_ids,selection%packet_action,effective_ids,&
    added_parent,added_operation,closure_fingerprint,ok,message)
  call require(ok.and.all(effective_ids==[1,2]),'production packet closure failed: '//trim(message))
  call freeze_dg_hybrid_production_selection(comm,selection,effective_ids,catalog,fingerprint,ok,message)
  call require(ok.and.catalog%valid.and.size(catalog%packets)==2,&
    'closed production packet selection did not freeze: '//trim(message))
  call freeze_dg_hybrid_production_selection(comm,selection,[1],catalog,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'closed')>0,'non-closed production packet selection was accepted')
  call freeze_dg_hybrid_production_selection(comm,selection,[3],catalog,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'closed')>0,'out-of-universe production packet ID was accepted')
  selection%packet_action(1,2)=1
  call freeze_dg_hybrid_production_selection(comm,selection,[1,2],catalog,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'receipt')>0,'mutated production action receipt was accepted')
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
    windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok,'production analysis refresh failed: '//trim(message))
  selection%packets(1)%owner_rank=-99
  call freeze_dg_hybrid_production_selection(comm,selection,[1,2],catalog,fingerprint,ok,message)
  call require(ok.and.catalog%packets(1)%owner_rank==0,&
    'effective production ownership was copied instead of recomputed')
  if(nproc>1)then
    if(rank==1)selection%analysis_fingerprint=selection%analysis_fingerprint+1_int64
    call freeze_dg_hybrid_production_selection(comm,selection,[1,2],catalog,fingerprint,ok,message)
    call require(.not.ok.and.index(message,'fingerprint')>0,&
      'rank-local production receipt corruption was not rejected collectively')
    call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
      core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
      windows,g_vectors,selection,workspace,fingerprint,ok,message)
    call require(ok,'production analysis refresh after receipt corruption failed: '//trim(message))
  endif
  if(nproc>1)then
    if(rank==0)then
      allocate(root_core_ids(4),source=box_ids)
      allocate(root_core_fragment_ids(4),source=[1,1,2,2])
    else
      allocate(root_core_ids(0),root_core_fragment_ids(0))
    endif
    allocate(root_coordinates(3,size(root_core_ids)),root_row_action(size(root_core_ids),2))
    root_coordinates=0d0
    do p=1,size(root_core_ids)
      root_coordinates(1,p)=real(root_core_ids(p)-1_int64,real64)
      root_row_action(p,:)=row_action(int(root_core_ids(p)),:)
    enddo
    call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,root_core_ids,&
      root_core_fragment_ids,root_coordinates,root_row_action,reciprocal_lattice,reciprocal_rotation,&
      0d0,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
    call require(ok,'ambiguous root-owned production row layout failed: '//trim(message))
    deallocate(root_core_ids,root_core_fragment_ids,root_coordinates,root_row_action)
  endif
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action(:,1:1),reciprocal_lattice,reciprocal_rotation(:,:,1:1),&
    0d0,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok.and.selection%analysis_complete.and.selection%identity_only,&
    'explicit identity-only production analysis failed: '//trim(message))
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action(:,2:2),reciprocal_lattice,reciprocal_rotation(:,:,2:2),&
    0d0,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'group')>0,&
    'production analysis accepted an operation list without explicit identity')
  if(nproc>1)then
    if(mod(rank,2)==0)then
      effective_ids=[1,2]
    else
      effective_ids=[2,1]
    endif
    call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
      core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
      windows,g_vectors,selection,workspace,fingerprint,ok,message)
    call require(ok,'production analysis setup for distributed freeze failed: '//trim(message))
    call freeze_dg_hybrid_production_selection(comm,selection,effective_ids,catalog,fingerprint,ok,message)
    call require(.not.ok.and.index(message,'rank')>0,&
      'production freeze accepted rank-dependent effective-ID ordering')
  endif
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0.5d0,2,1d-12,&
    windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok.and.selection%shell_added==2.and.&
    all(selection%requested_packet_ids==selection%packet_ids),&
    'strict production selection omitted a cutoff-completion shell: '//trim(message))
  call freeze_dg_hybrid_production_selection(comm,selection,selection%packet_ids,catalog,fingerprint,ok,message)
  call require(ok,'strict cutoff-complete production selection did not freeze: '//trim(message))
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0.6d0,2,1d-12,&
    windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok.and.selection%effective_cutoff<selection%requested_cutoff,&
    'between-shell effective cutoff receipt was not exposed: '//trim(message))
  call freeze_dg_hybrid_production_selection(comm,selection,selection%packet_ids,catalog,fingerprint,ok,message)
  call require(ok,'between-shell production selection did not freeze: '//trim(message))
  allocate(mixed_row_action(4,4));mixed_row_action(:,1:2)=row_action
  mixed_row_action(:,3)=[1,3,2,4]
  mixed_row_action(:,4)=[4,2,3,1]
  mixed_reciprocal_rotation(:,:,1:2)=reciprocal_rotation
  mixed_reciprocal_rotation(:,:,3)=0d0
  mixed_reciprocal_rotation(1,2,3)=1d0;mixed_reciprocal_rotation(2,1,3)=1d0
  mixed_reciprocal_rotation(3,3,3)=1d0
  mixed_reciprocal_rotation(:,:,4)=-mixed_reciprocal_rotation(:,:,3)
  call analyze_dg_hybrid_production_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,mixed_row_action,reciprocal_lattice,mixed_reciprocal_rotation,0d0,2,1d-12,&
    windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'whole fragments')>0,&
    'authoritative production analysis silently downgraded a known physical group')
  call analyze_dg_hybrid_lcfo_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,mixed_row_action,reciprocal_lattice,mixed_reciprocal_rotation,&
    701_int64,0.5d0,7,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok,'LCFO-deferred production selection rejected split fragment action: '//trim(message))
  call require(selection%lcfo_symmetry_deferred,'LCFO deferral provenance is missing')
  call require(selection%wannier_symmetry_fingerprint==701_int64,'Wannier provenance was not retained')
  call require(selection%operation_count==4.and.selection%window_operation_count==1.and.&
    selection%identity_only,'fragment-local bookkeeping and physical reciprocal operations were conflated')
  call require(size(selection%row_action,2)==1.and.size(selection%reciprocal_action,2)==4.and.&
    size(selection%reciprocal_rotation,3)==4,&
    'LCFO-deferred selection discarded the authoritative reciprocal operation catalog')
  call require(selection%pw_mode_count==7.and.selection%shell_added==2.and.selection%orbit_added==2.and.&
    selection%requested_cutoff==0.5d0.and.selection%effective_cutoff>selection%requested_cutoff,&
    'LCFO-deferred PW cutoff-completion receipt is incorrect')
  call require(all(selection%requested_packet_ids==selection%packet_ids),&
    'LCFO-deferred preparation must retain the complete PW packet catalog')
  allocate(coset_row_action(4,4))
  coset_row_action(:,1)=[1,2,3,4]
  coset_row_action(:,2)=[4,3,2,1]
  coset_row_action(:,3)=[2,1,3,4]
  coset_row_action(:,4)=[1,2,4,3]
  call analyze_dg_hybrid_lcfo_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,coset_row_action,reciprocal_lattice,mixed_reciprocal_rotation,&
    701_int64,0.5d0,7,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(ok.and.selection%operation_count==4.and.selection%window_operation_count==1,&
    'LCFO-deferred selection rejected nonclosed affine coset representatives: '//trim(message))
  call analyze_dg_hybrid_lcfo_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,mixed_row_action,reciprocal_lattice,mixed_reciprocal_rotation,&
    701_int64,0.5d0,6,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'capacity')>0,&
    'wannier_pw_max clipped a completed shell/orbit instead of reporting capacity failure')
  call analyze_dg_hybrid_lcfo_selection(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,mixed_row_action,reciprocal_lattice,mixed_reciprocal_rotation,&
    0_int64,0.5d0,7,2,1d-12,windows,g_vectors,selection,workspace,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'Wannier symmetry provenance')>0,&
    'LCFO-deferred production selection accepted missing Wannier symmetry provenance')
  if(rank==0)write(*,'(a,i0,a,i0)')'PRODUCTION_PW ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid production PW basis on ',nproc,' ranks'
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::bad,global_bad
    bad=merge(0,1,condition)
#ifdef USE_MPI
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    global_bad=bad
#endif
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(text)
#ifdef USE_MPI
      call MPI_Abort(comm,1,ierr)
#else
      error stop 1
#endif
    endif
  end subroutine require
end program test_dg_hybrid_production_pw_basis_mpi
