#include "config.h"
module dg_hybrid_production_pw_basis
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog,s_dg_hybrid_production_selection
  use dg_hybrid_reciprocal_catalog,only:build_dg_hybrid_reciprocal_catalog
  use dg_hybrid_window_distribution,only:prepare_dg_hybrid_window_distribution
  use dg_hybrid_windowed_pw_basis,only:build_dg_hybrid_windowed_pw_basis
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_dg_hybrid_production_pw_basis,analyze_dg_hybrid_production_selection,&
    analyze_dg_hybrid_lcfo_selection,freeze_dg_hybrid_production_selection
contains
  subroutine analyze_dg_hybrid_production_selection(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,coordinates,row_action,&
      reciprocal_lattice,reciprocal_rotation,cutoff,tile_width,tolerance,windows,g_vectors,&
      selection,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:),core_fragment_ids(:)
    integer(int64),intent(in)::box_ids(:),core_ids(:)
    real(real64),intent(in)::box_windows(:,:),coordinates(:,:),reciprocal_lattice(3,3),&
      reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,intent(in)::row_action(:,:),tile_width
    real(real64),allocatable,intent(out)::windows(:,:),g_vectors(:,:)
    type(s_dg_hybrid_production_selection),intent(out)::selection
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_basis_catalog)::universe_catalog
    integer,allocatable::normalized_row_action(:,:),fragment_action(:,:),g_integer(:,:),g_action(:,:),&
      g_star(:),g_conjugate(:)
    real(real64),allocatable::raw_windows(:,:)
    integer(int64)::window_workspace,window_fingerprint,reciprocal_fingerprint,basis_workspace,basis_fingerprint
    real(real64)::effective_cutoff
    integer::i,j,op,packet,target_fragment,target_packet,shell_added,orbit_added
    logical::stage_ok
    character(256)::stage_message

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    selection%analysis_complete=.false.
    call normalize_dg_hybrid_production_row_action(comm,global_point_count,core_ids,row_action,&
      normalized_row_action,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call validate_dg_hybrid_production_rotation_extent(comm,size(normalized_row_action,2),&
      reciprocal_rotation,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call prepare_dg_hybrid_window_distribution(comm,global_point_count,fragment_count,fragment_ids,&
      box_ids,box_windows,core_ids,core_fragment_ids,normalized_row_action,raw_windows,fragment_action,&
      window_workspace,window_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,tolerance,&
      g_integer,g_vectors,g_action,g_star,g_conjugate,reciprocal_fingerprint,effective_cutoff,&
      shell_added,orbit_added,stage_ok,stage_message)
    if(.not.stage_ok)then
      if(index(stage_message,'identity')>0.or.index(stage_message,'closed')>0.or.&
          index(stage_message,'duplicate')>0)then
        message='authoritative production operations do not form an explicit group: '//trim(stage_message)
      else
        message=trim(stage_message)
      endif
      return
    endif
    if(.not.valid_authoritative_group(normalized_row_action,fragment_action,g_action,&
      reciprocal_rotation,tolerance))then
      message='authoritative production operations do not form an explicit group';return
    endif
    call build_dg_hybrid_windowed_pw_basis(comm,global_point_count,core_ids,coordinates,raw_windows,&
      fragment_action,normalized_row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,&
      tile_width,tolerance,windows,universe_catalog,basis_workspace,basis_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    selection%window_operation_count=size(normalized_row_action,2)
    selection%operation_count=size(g_action,2)
    selection%pw_mode_count=size(g_vectors,2)
    selection%requested_cutoff=cutoff
    selection%effective_cutoff=effective_cutoff
    selection%shell_added=shell_added
    selection%orbit_added=orbit_added
    selection%identity_only=selection%window_operation_count==1.and.selection%operation_count==1
    if(selection%identity_only)then
      selection%identity_only=all(normalized_row_action(:,1)==[(i,i=1,global_point_count)]).and.&
        all(fragment_action(:,1)==[(i,i=1,fragment_count)]).and.&
        all(g_action(:,1)==[(i,i=1,size(g_action,1))])
    endif
    allocate(selection%fragment_action,source=fragment_action)
    allocate(selection%row_action,source=normalized_row_action)
    allocate(selection%reciprocal_action,source=g_action)
    allocate(selection%reciprocal_rotation,source=reciprocal_rotation)
    allocate(selection%packet_ids(size(universe_catalog%packets)),&
      selection%packet_action(size(universe_catalog%packets),selection%operation_count))
    selection%packet_ids=[(i,i=1,size(selection%packet_ids))]
    do op=1,selection%operation_count
      do packet=1,size(universe_catalog%packets)
        target_fragment=fragment_action(universe_catalog%packets(packet)%fragment_id,op)
        target_packet=0
        do j=1,size(universe_catalog%packets)
          if(universe_catalog%packets(j)%fragment_id==target_fragment.and.&
            universe_catalog%packets(j)%star_id==universe_catalog%packets(packet)%star_id)then
            target_packet=j;exit
          endif
        enddo
        if(target_packet==0)then;message='production packet action is incomplete';return;endif
        selection%packet_action(packet,op)=target_packet
      enddo
    enddo
    allocate(selection%requested_packet_ids,source=selection%packet_ids)
    call move_alloc(universe_catalog%packets,selection%packets)
    selection%window_fingerprint=window_fingerprint
    selection%packet_fingerprint=ieor(reciprocal_fingerprint,basis_fingerprint)
    fingerprint=production_selection_fingerprint(selection)
    selection%analysis_fingerprint=fingerprint
    selection%analysis_complete=.true.
    workspace_peak_bytes=max(window_workspace,basis_workspace)
    ok=.true.
  end subroutine analyze_dg_hybrid_production_selection

  subroutine analyze_dg_hybrid_lcfo_selection(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,coordinates,row_action,&
      reciprocal_lattice,reciprocal_rotation,wannier_symmetry_fingerprint,cutoff,maximum_pw_count,tile_width,&
      tolerance,windows,g_vectors,selection,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:),core_fragment_ids(:),maximum_pw_count
    integer(int64),intent(in)::box_ids(:),core_ids(:),wannier_symmetry_fingerprint
    real(real64),intent(in)::box_windows(:,:),coordinates(:,:),reciprocal_lattice(3,3),&
      reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,intent(in)::row_action(:,:),tile_width
    real(real64),allocatable,intent(out)::windows(:,:),g_vectors(:,:)
    type(s_dg_hybrid_production_selection),intent(out)::selection
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_basis_catalog)::universe_catalog
    integer,allocatable::normalized_physical_row_action(:,:),identity_row_action(:,:),fragment_action(:,:),&
      g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable::raw_windows(:,:)
    real(real64)::effective_cutoff
    integer::i,op,packet,ierr,local_bad,global_bad,minimum_integer,maximum_integer,shell_added,orbit_added
    integer(int64)::minimum_fingerprint,maximum_fingerprint,window_workspace,window_fingerprint,&
      reciprocal_fingerprint,basis_workspace,basis_fingerprint,physical_row_fingerprint
    logical::stage_ok
    character(256)::stage_message

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    selection%analysis_complete=.false.
    local_bad=merge(0,1,wannier_symmetry_fingerprint/=0_int64)
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier symmetry provenance validation failed';return;endif
    call MPI_Allreduce(wannier_symmetry_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier symmetry provenance rank agreement failed';return;endif
    call MPI_Allreduce(wannier_symmetry_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='Wannier symmetry provenance rank agreement failed';return
    endif
    call MPI_Allreduce(maximum_pw_count,minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='PW capacity rank agreement failed';return;endif
    call MPI_Allreduce(maximum_pw_count,maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='PW capacity rank agreement failed';return
    endif
#else
    global_bad=local_bad
    minimum_integer=maximum_pw_count;maximum_integer=maximum_pw_count
#endif
    if(global_bad/=0)then;message='Wannier symmetry provenance is missing';return;endif
    if(minimum_integer<0)then;message='PW capacity must be nonnegative';return;endif
    call normalize_dg_hybrid_production_row_action(comm,global_point_count,core_ids,row_action,&
      normalized_physical_row_action,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call validate_dg_hybrid_production_rotation_extent(comm,size(normalized_physical_row_action,2),&
      reciprocal_rotation,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    if(.not.valid_row_action_catalog(normalized_physical_row_action))then
      message='authoritative physical row operations are not valid permutations';return
    endif
    physical_row_fingerprint=int(z'1F83D9ABFB41BD6B',int64)
    do op=1,size(normalized_physical_row_action,2);do i=1,global_point_count
      physical_row_fingerprint=ieor(ishftc(physical_row_fingerprint,7),&
        int(normalized_physical_row_action(i,op),int64))
    enddo;enddo
    allocate(identity_row_action(global_point_count,1))
    identity_row_action(:,1)=[(i,i=1,global_point_count)]
    call prepare_dg_hybrid_window_distribution(comm,global_point_count,fragment_count,fragment_ids,&
      box_ids,box_windows,core_ids,core_fragment_ids,identity_row_action,raw_windows,fragment_action,&
      window_workspace,window_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,tolerance,&
      g_integer,g_vectors,g_action,g_star,g_conjugate,reciprocal_fingerprint,effective_cutoff,&
      shell_added,orbit_added,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    if(maximum_pw_count>0.and.size(g_vectors,2)>maximum_pw_count)then
      write(message,'(a,i0,a,i0)')'completed reciprocal catalog exceeds PW capacity: required=',&
        size(g_vectors,2),' capacity=',maximum_pw_count
      return
    endif
    call build_dg_hybrid_windowed_pw_basis(comm,global_point_count,core_ids,coordinates,raw_windows,&
      fragment_action,identity_row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,&
      tile_width,tolerance,windows,universe_catalog,basis_workspace,basis_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    selection%window_operation_count=1
    selection%operation_count=size(g_action,2)
    selection%pw_mode_count=size(g_vectors,2)
    selection%requested_cutoff=cutoff
    selection%effective_cutoff=effective_cutoff
    selection%shell_added=shell_added
    selection%orbit_added=orbit_added
    selection%identity_only=.true.
    allocate(selection%fragment_action,source=fragment_action)
    allocate(selection%row_action,source=identity_row_action)
    allocate(selection%reciprocal_action,source=g_action)
    allocate(selection%reciprocal_rotation,source=reciprocal_rotation)
    allocate(selection%packet_ids(size(universe_catalog%packets)),&
      selection%packet_action(size(universe_catalog%packets),selection%operation_count))
    selection%packet_ids=[(i,i=1,size(selection%packet_ids))]
    do op=1,selection%operation_count;do packet=1,size(selection%packet_ids)
      selection%packet_action(packet,op)=packet
    enddo;enddo
    allocate(selection%requested_packet_ids,source=selection%packet_ids)
    call move_alloc(universe_catalog%packets,selection%packets)
    selection%lcfo_symmetry_deferred=.true.
    selection%wannier_symmetry_fingerprint=wannier_symmetry_fingerprint
    selection%window_fingerprint=window_fingerprint
    selection%packet_fingerprint=ieor(ieor(reciprocal_fingerprint,basis_fingerprint),physical_row_fingerprint)
    fingerprint=production_selection_fingerprint(selection)
    selection%analysis_fingerprint=fingerprint
    selection%analysis_complete=.true.
    workspace_peak_bytes=max(window_workspace,basis_workspace)
    ok=.true.
  end subroutine analyze_dg_hybrid_lcfo_selection

  subroutine freeze_dg_hybrid_production_selection(comm,selection,effective_ids,catalog,fingerprint,ok,message)
    integer,intent(in)::comm,effective_ids(:)
    type(s_dg_hybrid_production_selection),intent(in)::selection
    type(s_dg_hybrid_basis_catalog),intent(out)::catalog
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,op,position,local_bad,global_bad,ierr,minimum_integer,maximum_integer,nproc
    integer(int64)::receipt_fingerprint,minimum_fingerprint,maximum_fingerprint

    ok=.false.;message='';fingerprint=0_int64;catalog%valid=.false.;local_bad=0
#ifdef USE_MPI
    call MPI_Allreduce(size(effective_ids),minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production effective-ID rank agreement failed';return;endif
    call MPI_Allreduce(size(effective_ids),maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='production effective-ID rank agreement failed';return
    endif
    do i=1,size(effective_ids)
      call MPI_Allreduce(effective_ids(i),minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='production effective-ID rank agreement failed';return;endif
      call MPI_Allreduce(effective_ids(i),maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='production effective-ID rank agreement failed';return
      endif
    enddo
#endif
    if(.not.selection%analysis_complete.or.selection%analysis_fingerprint==0_int64.or.&
      (selection%lcfo_symmetry_deferred.and.selection%wannier_symmetry_fingerprint==0_int64).or.&
      size(effective_ids)<1.or..not.allocated(selection%packet_ids).or.&
      .not.allocated(selection%requested_packet_ids).or.&
      .not.allocated(selection%packet_action).or..not.allocated(selection%packets).or.&
      .not.allocated(selection%fragment_action).or..not.allocated(selection%row_action).or.&
      .not.allocated(selection%reciprocal_action).or..not.allocated(selection%reciprocal_rotation))local_bad=1
    if(local_bad==0)then
      if(selection%operation_count<1.or.selection%window_operation_count<1.or.&
        selection%pw_mode_count<1.or.selection%shell_added<0.or.selection%orbit_added<0.or.&
        .not.ieee_is_finite(selection%requested_cutoff).or..not.ieee_is_finite(selection%effective_cutoff).or.&
        selection%requested_cutoff<0d0.or.selection%effective_cutoff<0d0.or.&
        size(selection%packets)/=size(selection%packet_ids).or.&
        any(shape(selection%packet_action)/=[size(selection%packet_ids),selection%operation_count]).or.&
        size(selection%fragment_action,2)/=selection%window_operation_count.or.&
        size(selection%row_action,2)/=selection%window_operation_count.or.&
        any(shape(selection%reciprocal_action)/=[selection%pw_mode_count,selection%operation_count]).or.&
        any(shape(selection%reciprocal_rotation)/=[3,3,selection%operation_count]))local_bad=1
    endif
    if(local_bad==0)then
      do i=1,size(selection%packets)
        if(.not.allocated(selection%packets(i)%g_indices))local_bad=1
        if(allocated(selection%packets(i)%g_indices))then
          if(any(selection%packets(i)%g_indices<1).or.&
            any(selection%packets(i)%g_indices>selection%pw_mode_count))local_bad=1
        endif
      enddo
    endif
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production receipt validation reduction failed';return;endif
#else
    global_bad=local_bad
#endif
    if(global_bad/=0)then;message='invalid production analysis receipt';return;endif
    receipt_fingerprint=production_selection_fingerprint(selection)
#ifdef USE_MPI
    call MPI_Allreduce(receipt_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production receipt rank agreement failed';return;endif
    call MPI_Allreduce(receipt_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='production receipt rank agreement failed';return
    endif
#endif
    local_bad=merge(0,1,receipt_fingerprint==selection%analysis_fingerprint)
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production receipt comparison reduction failed';return;endif
#else
    global_bad=local_bad
#endif
    if(global_bad/=0)then;message='production analysis receipt fingerprint mismatch';return;endif
    local_bad=0
    do i=1,size(effective_ids)
      if(count(selection%packet_ids==effective_ids(i))/=1.or.count(effective_ids==effective_ids(i))/=1)local_bad=1
    enddo
    if(local_bad==0)then
      do i=1,size(effective_ids)
        position=findloc(selection%packet_ids,effective_ids(i),dim=1)
        do op=1,selection%operation_count
          if(count(effective_ids==selection%packet_action(position,op))/=1)local_bad=1
        enddo
      enddo
    endif
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production selection validation reduction failed';return;endif
#else
    global_bad=local_bad
#endif
    if(global_bad/=0)then;message='effective production packet selection is not closed';return;endif
    allocate(catalog%packets(size(effective_ids)))
#ifdef USE_MPI
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production effective ownership communicator failed';return;endif
#else
    nproc=1
#endif
    fingerprint=selection%analysis_fingerprint
    do i=1,size(effective_ids)
      position=findloc(selection%packet_ids,effective_ids(i),dim=1)
      catalog%packets(i)=selection%packets(position)
      catalog%packets(i)%owner_rank=mod(catalog%packets(i)%fragment_id-1,nproc)
      fingerprint=ieor(ishftc(fingerprint,11),int(effective_ids(i),int64))
      fingerprint=ieor(ishftc(fingerprint,11),int(catalog%packets(i)%owner_rank,int64))
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    catalog%window_fingerprint=selection%window_fingerprint
    catalog%packet_fingerprint=fingerprint
    catalog%catalog_fingerprint=ieor(selection%analysis_fingerprint,ishftc(fingerprint,17))
    if(catalog%catalog_fingerprint==0_int64)catalog%catalog_fingerprint=1_int64
    catalog%operation_count=selection%operation_count
    catalog%window_operation_count=selection%window_operation_count
    catalog%pw_mode_count=selection%pw_mode_count
    catalog%shell_added=selection%shell_added
    catalog%orbit_added=selection%orbit_added
    catalog%requested_cutoff=selection%requested_cutoff
    catalog%effective_cutoff=selection%effective_cutoff
    catalog%wannier_fingerprint=selection%wannier_symmetry_fingerprint
    catalog%valid=.true.;ok=.true.
  end subroutine freeze_dg_hybrid_production_selection

  integer(int64) function production_selection_fingerprint(selection) result(fingerprint)
    type(s_dg_hybrid_production_selection),intent(in)::selection
    integer::i,j,op
    integer(int64)::bits

    fingerprint=ieor(int(z'510E527FADE682D1',int64),selection%window_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,7),selection%packet_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,7),int(selection%operation_count,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(selection%window_operation_count,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(selection%pw_mode_count,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(selection%shell_added,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(selection%orbit_added,int64))
    bits=transfer(selection%requested_cutoff,bits)
    fingerprint=ieor(ishftc(fingerprint,7),bits)
    bits=transfer(selection%effective_cutoff,bits)
    fingerprint=ieor(ishftc(fingerprint,7),bits)
    fingerprint=ieor(ishftc(fingerprint,7),merge(1_int64,0_int64,selection%identity_only))
    fingerprint=ieor(ishftc(fingerprint,7),merge(1_int64,0_int64,selection%lcfo_symmetry_deferred))
    fingerprint=ieor(ishftc(fingerprint,7),selection%wannier_symmetry_fingerprint)
    do op=1,selection%window_operation_count
      do i=1,size(selection%fragment_action,1)
        fingerprint=ieor(ishftc(fingerprint,7),int(selection%fragment_action(i,op),int64))
      enddo
      do i=1,size(selection%row_action,1)
        fingerprint=ieor(ishftc(fingerprint,7),int(selection%row_action(i,op),int64))
      enddo
    enddo
    do op=1,selection%operation_count
      do i=1,size(selection%reciprocal_action,1)
        fingerprint=ieor(ishftc(fingerprint,7),int(selection%reciprocal_action(i,op),int64))
      enddo
      do j=1,3;do i=1,3
        bits=transfer(selection%reciprocal_rotation(i,j,op),bits)
        fingerprint=ieor(ishftc(fingerprint,7),bits)
      enddo;enddo
    enddo
    do i=1,size(selection%packet_ids)
      fingerprint=ieor(ishftc(fingerprint,11),int(selection%packet_ids(i),int64))
      do op=1,selection%operation_count
        fingerprint=ieor(ishftc(fingerprint,11),int(selection%packet_action(i,op),int64))
      enddo
      fingerprint=ieor(ishftc(fingerprint,11),int(selection%packets(i)%fragment_id,int64))
      fingerprint=ieor(ishftc(fingerprint,11),int(selection%packets(i)%star_id,int64))
      do j=1,size(selection%packets(i)%g_indices)
        fingerprint=ieor(ishftc(fingerprint,11),int(selection%packets(i)%g_indices(j),int64))
      enddo
    enddo
    do i=1,size(selection%requested_packet_ids)
      fingerprint=ieor(ishftc(fingerprint,11),int(selection%requested_packet_ids(i),int64))
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
  end function production_selection_fingerprint

  subroutine normalize_dg_hybrid_production_row_action(comm,global_point_count,core_ids,row_action,&
      normalized_row_action,ok,message)
    integer,intent(in)::comm,global_point_count,row_action(:,:)
    integer(int64),intent(in)::core_ids(:)
    integer,allocatable,intent(out)::normalized_row_action(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,op,ierr,local_full,all_full,local_distributed,all_distributed,local_bad,global_bad,&
      minimum_integer,maximum_integer

    ok=.false.;message=''
#ifdef USE_MPI
    call MPI_Allreduce(size(row_action,2),minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production operation-count agreement failed';return;endif
    call MPI_Allreduce(size(row_action,2),maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer.or.minimum_integer<1)then
      message='empty or inconsistent authoritative production operation list';return
    endif
    local_full=merge(1,0,size(row_action,1)==global_point_count)
    local_distributed=merge(1,0,size(row_action,1)==size(core_ids))
    call MPI_Allreduce(local_full,all_full,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production row-action layout agreement failed';return;endif
    call MPI_Allreduce(local_distributed,all_distributed,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.(all_full==0.and.all_distributed==0))then
      message='invalid or rank-inconsistent production row-action extent';return
    endif
    allocate(normalized_row_action(global_point_count,minimum_integer));normalized_row_action=0
    if(all_full==1)then
      normalized_row_action=row_action
    else
      local_bad=merge(0,1,all(core_ids>=1_int64.and.core_ids<=int(global_point_count,int64)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        message='invalid distributed production row IDs';return
      endif
      do op=1,minimum_integer;do i=1,size(core_ids)
        normalized_row_action(int(core_ids(i)),op)=row_action(i,op)
      enddo;enddo
      call MPI_Allreduce(MPI_IN_PLACE,normalized_row_action,size(normalized_row_action),MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='distributed production row-action reduction failed';return;endif
    endif
#else
    if(size(row_action,2)<1.or.size(row_action,1)/=global_point_count)then
      message='non-MPI production row action must be complete';return
    endif
    allocate(normalized_row_action,source=row_action)
#endif
    ok=.true.
  end subroutine normalize_dg_hybrid_production_row_action

  subroutine validate_dg_hybrid_production_rotation_extent(comm,operation_count,rotation,ok,message)
    integer,intent(in)::comm,operation_count
    real(real64),intent(in)::rotation(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::local_bad,global_bad,ierr

    ok=.false.;message=''
    local_bad=merge(0,1,all(shape(rotation)==[3,3,operation_count]))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='rank-inconsistent production reciprocal-operation extent';return
    endif
#else
    if(local_bad/=0)then;message='invalid production reciprocal-operation extent';return;endif
#endif
    ok=.true.
  end subroutine validate_dg_hybrid_production_rotation_extent

  logical function valid_row_action_catalog(row_action)
    integer,intent(in)::row_action(:,:)
    integer::operation,i

    valid_row_action_catalog=.false.
    if(size(row_action,1)<1.or.size(row_action,2)<1.or.any(row_action<1).or.&
      any(row_action>size(row_action,1)))return
    do operation=1,size(row_action,2)
      do i=1,size(row_action,1)
        if(count(row_action(:,operation)==i)/=1)return
      enddo
    enddo
    valid_row_action_catalog=.true.
  end function valid_row_action_catalog

  logical function valid_authoritative_group(row_action,fragment_action,g_action,rotation,tolerance)
    integer,intent(in)::row_action(:,:),fragment_action(:,:),g_action(:,:)
    real(real64),intent(in)::rotation(:,:,:),tolerance
    integer::a,b,c,i,identity_count,match_count
    logical::matches

    valid_authoritative_group=.false.
    if(size(row_action,2)<1.or.size(fragment_action,2)/=size(row_action,2).or.&
      size(g_action,2)/=size(row_action,2).or.size(rotation,3)/=size(row_action,2))return
    identity_count=0
    do a=1,size(row_action,2)
      if(all(row_action(:,a)==[(i,i=1,size(row_action,1))]).and.&
        all(fragment_action(:,a)==[(i,i=1,size(fragment_action,1))]).and.&
        all(g_action(:,a)==[(i,i=1,size(g_action,1))]).and.&
        maxval(abs(rotation(:,:,a)-identity_matrix()))<=100d0*tolerance)identity_count=identity_count+1
      do b=a+1,size(row_action,2)
        if(all(row_action(:,a)==row_action(:,b)).and.&
          maxval(abs(rotation(:,:,a)-rotation(:,:,b)))<=100d0*tolerance)return
      enddo
    enddo
    if(identity_count/=1)return
    do a=1,size(row_action,2);do b=1,size(row_action,2)
      match_count=0
      do c=1,size(row_action,2)
        matches=all(row_action(row_action(:,a),b)==row_action(:,c)).and.&
          all(fragment_action(fragment_action(:,a),b)==fragment_action(:,c)).and.&
          all(g_action(g_action(:,a),b)==g_action(:,c)).and.&
          maxval(abs(matmul(rotation(:,:,b),rotation(:,:,a))-rotation(:,:,c)))<=100d0*tolerance
        if(matches)match_count=match_count+1
      enddo
      if(match_count/=1)return
    enddo;enddo
    valid_authoritative_group=.true.
  contains
    pure function identity_matrix() result(identity)
      real(real64)::identity(3,3)
      integer::axis
      identity=0d0
      do axis=1,3;identity(axis,axis)=1d0;enddo
    end function identity_matrix
  end function valid_authoritative_group

  subroutine build_dg_hybrid_production_pw_basis(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,coordinates,row_action,&
      reciprocal_lattice,reciprocal_rotation,cutoff,tile_width,tolerance,windows,g_vectors,&
      catalog,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:),core_fragment_ids(:)
    integer(int64),intent(in)::box_ids(:),core_ids(:)
    real(real64),intent(in)::box_windows(:,:),coordinates(:,:),reciprocal_lattice(3,3),&
      reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,intent(in)::row_action(:,:),tile_width
    real(real64),allocatable,intent(out)::windows(:,:),g_vectors(:,:)
    type(s_dg_hybrid_basis_catalog),intent(out)::catalog
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::fragment_action(:,:),g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable::raw_windows(:,:)
    integer(int64)::window_workspace,window_fingerprint,reciprocal_fingerprint,basis_workspace,basis_fingerprint
    real(real64)::effective_cutoff
    integer::shell_added,orbit_added
    logical::stage_ok
    character(256)::stage_message
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    catalog%valid=.false.
    call prepare_dg_hybrid_window_distribution(comm,global_point_count,fragment_count,fragment_ids,&
      box_ids,box_windows,core_ids,core_fragment_ids,row_action,raw_windows,fragment_action,&
      window_workspace,window_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,tolerance,&
      g_integer,g_vectors,g_action,g_star,g_conjugate,reciprocal_fingerprint,effective_cutoff,&
      shell_added,orbit_added,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_windowed_pw_basis(comm,global_point_count,core_ids,coordinates,raw_windows,&
      fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,&
      tile_width,tolerance,windows,catalog,basis_workspace,basis_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);catalog%valid=.false.;return;endif
    catalog%window_fingerprint=window_fingerprint
    catalog%packet_fingerprint=ieor(reciprocal_fingerprint,basis_fingerprint)
    catalog%requested_cutoff=cutoff
    catalog%effective_cutoff=effective_cutoff
    catalog%shell_added=shell_added
    catalog%orbit_added=orbit_added
    catalog%catalog_fingerprint=ieor(ieor(ishftc(window_fingerprint,11),&
      ishftc(reciprocal_fingerprint,23)),basis_fingerprint)
    if(catalog%catalog_fingerprint==0_int64)catalog%catalog_fingerprint=1_int64
    fingerprint=catalog%catalog_fingerprint
    workspace_peak_bytes=max(window_workspace,basis_workspace)
    deallocate(raw_windows,fragment_action,g_integer,g_action,g_star,g_conjugate)
    ok=.true.
  end subroutine build_dg_hybrid_production_pw_basis
end module dg_hybrid_production_pw_basis
