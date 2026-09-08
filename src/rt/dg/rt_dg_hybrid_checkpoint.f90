#include "config.h"
module rt_dg_hybrid_checkpoint
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter::checkpoint_version=2,legacy_checkpoint_version=1
  integer,parameter,public::rt_dg_hybrid_checkpoint_version=2
  character(16),parameter::checkpoint_magic='SALMON_DG_HYB01 '
  character(16),parameter::occupied_magic='SALMON_DG_OCC02 '
  integer,parameter::occupied_version=2
  integer,parameter,public::rt_dg_hybrid_occupied_checkpoint_version=2
  integer,parameter::ground_state_version=3
  integer,parameter,public::rt_dg_hybrid_ground_state_checkpoint_version=3
  integer,parameter,public::rt_dg_hybrid_energy_window_explicit=1,&
    rt_dg_hybrid_energy_window_legacy_dynamic=2
  ! Version-3 vector operator item one is the homogeneous velocity-gauge
  ! coupling.  Cell-wrapped position remains a separately named payload.
  integer,parameter,public::rt_dg_hybrid_vector_canonical_momentum=1
  character(16),parameter::ground_state_magic='SALMON_DG_GS001 '
  integer,parameter::ground_state_logical_count=15,ground_state_integer_count=25,&
    ground_state_fingerprint_count=56,ground_state_real_count=39

  type,public::s_rt_dg_hybrid_construction_catalog
    logical::valid=.false.
    integer::global_count=0
    integer(int64),allocatable::ids(:)
    integer,allocatable::generations(:),ordering(:),ownership(:)
    integer(int64)::ids_fingerprint=0_int64,generation_fingerprint=0_int64,&
      ordering_fingerprint=0_int64,ownership_fingerprint=0_int64,&
      provenance_fingerprint=0_int64,catalog_fingerprint=0_int64
  end type s_rt_dg_hybrid_construction_catalog

  type,public::s_rt_dg_hybrid_certified_basis_payload
    logical::valid=.false.,localization_converged=.false.,localization_symmetry_constrained=.false.
    integer::construction_count=0,certified_count=0,occupied_count=0,localization_iterations=0
    integer(int64),allocatable::construction_row_ids(:),transformation_row_ids(:)
    complex(real64),allocatable::c_cert(:,:),u_rt(:,:),b_rt(:,:),initial_occupied_amplitudes(:,:)
    real(real64),allocatable::certified_eigenvalues(:),occupations(:),centers(:,:),&
      spreads_before(:),spreads_after(:)
    real(real64)::spread_before_total=0d0,spread_after_total=0d0,spread_improvement=0d0
    real(real64)::transform_unitarity_defect=huge(0d0),certified_metric_defect=huge(0d0),&
      rt_metric_defect=huge(0d0),embedding_defect=huge(0d0),projector_invariance_defect=huge(0d0)
    real(real64)::target_symmetry_defect_before=huge(0d0),target_symmetry_defect_after=huge(0d0),&
      energy_symmetry_defect_before=huge(0d0),energy_symmetry_defect_after=huge(0d0),&
      symmetry_defect_invariance=huge(0d0),scalar_covariance_defect=huge(0d0),&
      vector_covariance_defect=huge(0d0),tensor_covariance_defect=huge(0d0)
    integer(int64)::c_cert_fingerprint=0_int64,u_rt_fingerprint=0_int64,b_rt_fingerprint=0_int64,&
      initial_state_fingerprint=0_int64,transformation_fingerprint=0_int64,&
      operator_fingerprint=0_int64,fingerprint=0_int64
  end type s_rt_dg_hybrid_certified_basis_payload

  type,public::s_rt_dg_hybrid_electron_count_receipt
    logical::valid=.false.
    real(real64)::expected_count=0d0,actual_count=0d0,tolerance=0d0,defect=huge(0d0),&
      omitted_tail=huge(0d0),chemical_potential=0d0
    integer(int64)::fingerprint=0_int64
  end type s_rt_dg_hybrid_electron_count_receipt

  type,public::s_rt_dg_hybrid_rt_space_payload
    logical::valid=.false.
    integer::rank=0,operation_count=0,scalar_count=0,vector_count=0,tensor_count=0
    integer(int64),allocatable::row_ids(:)
    integer,allocatable::row_owner_keys(:),grid_owner_keys(:)
    complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      local_rows(:,:),sipg_rows(:,:),hamiltonian_rows(:,:),representation(:,:,:),&
      scalar_operator_rows(:,:,:),vector_operator_rows(:,:,:,:),tensor_operator_rows(:,:,:,:,:),&
      basis_values(:,:)
    real(real64),allocatable::cartesian_rotations(:,:,:),density(:)
    integer(int64)::metric_fingerprint=0_int64,kinetic_fingerprint=0_int64,&
      nonlocal_fingerprint=0_int64,local_fingerprint=0_int64,sipg_fingerprint=0_int64,&
      hamiltonian_fingerprint=0_int64,basis_fingerprint=0_int64,density_fingerprint=0_int64,&
      ownership_fingerprint=0_int64,scalar_fingerprint=0_int64,vector_fingerprint=0_int64,&
      tensor_fingerprint=0_int64,representation_fingerprint=0_int64,fingerprint=0_int64
  end type s_rt_dg_hybrid_rt_space_payload

  type,public::s_rt_dg_hybrid_energy_window_receipt
    logical::valid=.false.,compatibility_dynamic_rank=.false.,proof_state_present=.false.
    integer::mode=0,construction_rank=0,solved_rank=0,occupied_rank=0,requested_rank=0,&
      certified_rank=0,extension_states=0,boundary_cluster_rank=0,proof_status=0
    real(real64)::window_size=0d0,e_homo=0d0,requested_cutoff=0d0,certified_cutoff=0d0,&
      extension_energy=0d0,proof_energy=0d0
    integer(int64)::fingerprint=0_int64
  end type s_rt_dg_hybrid_energy_window_receipt

  type,public::s_rt_dg_hybrid_symmetry_receipt
    logical::valid=.false.
    integer::worst_operation=0
    real(real64)::occupied_subspace_defect=huge(0d0),occupied_projector_defect=huge(0d0),&
      target_subspace_defect=huge(0d0),target_energy_defect=huge(0d0),density_defect=huge(0d0),&
      scalar_covariance_defect=huge(0d0),vector_covariance_defect=huge(0d0),&
      tensor_covariance_defect=huge(0d0),final_basis_defect=huge(0d0),&
      worst_operation_defect=huge(0d0),maximum_physical_defect=huge(0d0)
    integer(int64)::fingerprint=0_int64
  end type s_rt_dg_hybrid_symmetry_receipt

  type,public::s_rt_dg_hybrid_handoff_receipts
    logical::valid=.false.
    integer(int64)::position_fingerprint=0_int64,nonlocal_fingerprint=0_int64,&
      face_fingerprint=0_int64,pseudopotential_fingerprint=0_int64,&
      transformation_fingerprint=0_int64,fingerprint=0_int64
  end type s_rt_dg_hybrid_handoff_receipts

  type,public::s_rt_dg_hybrid_ground_state_payload
    logical::valid=.false.,final_refresh_complete=.false.,analysis_complete=.false.,identity_only=.false.
    integer::global_count=0,global_grid_count=0,noccupied=0,operation_count=0,nonidentity_operation_count=0
    integer(int64)::catalog_fingerprint=0_int64,state_fingerprint=0_int64,metric_fingerprint=0_int64,&
      operator_structure_fingerprint=0_int64,operator_value_fingerprint=0_int64,&
      kinetic_fingerprint=0_int64,nonlocal_fingerprint=0_int64,local_fingerprint=0_int64,&
      sipg_fingerprint=0_int64,basis_fingerprint=0_int64,face_fingerprint=0_int64,&
      dc_seed_fingerprint=0_int64,continuation_fingerprint=0_int64,scope_fingerprint=0_int64,&
      selection_fingerprint=0_int64,&
      analysis_fingerprint=0_int64,pseudopotential_fingerprint=0_int64,energy_fingerprint=0_int64,&
      position_convention_fingerprint=0_int64,payload_fingerprint=0_int64
    integer(int64),allocatable::row_ids(:),grid_ids(:),face_ids(:),face_point_ids(:),nonlocal_ids(:)
    integer,allocatable::metric_row_offsets(:),metric_column_ids(:),operator_row_offsets(:),operator_column_ids(:),&
      partition_ids(:),face_metadata(:,:),face_offsets(:),face_weight_offsets(:),face_basis_offsets(:),&
      face_value_offsets(:),face_observable_offsets(:),face_basis_ids(:),nonlocal_owner(:),&
      requested_ids(:),effective_ids(:),added_ids(:),closure_parent(:),closure_reason(:),closure_action(:),&
      scope_selectors(:),xc_types(:)
    real(real64),allocatable::grid_weights(:),face_normals(:,:),face_weights(:),density(:),occupations(:),&
      eigenvalues(:),continuation_receipt(:),pseudopotential_receipt(:),energy_receipt(:)
    complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),local_rows(:,:),&
      sipg_rows(:,:),hamiltonian_rows(:,:),basis_values(:,:),face_values(:,:),nonlocal_values(:,:),&
      coefficients(:,:),interface_observables(:,:)
    complex(real64),allocatable::position_rows(:,:,:),symmetry_representation(:,:,:)
    type(s_rt_dg_hybrid_construction_catalog)::construction_catalog
    type(s_rt_dg_hybrid_certified_basis_payload)::certified_basis
    type(s_rt_dg_hybrid_electron_count_receipt)::electron_count
    type(s_rt_dg_hybrid_rt_space_payload)::rt_space
    type(s_rt_dg_hybrid_energy_window_receipt)::energy_window
    type(s_rt_dg_hybrid_symmetry_receipt)::symmetry_receipt
    type(s_rt_dg_hybrid_handoff_receipts)::handoff_receipts
  end type s_rt_dg_hybrid_ground_state_payload
  public::write_rt_dg_hybrid_checkpoint,read_rt_dg_hybrid_checkpoint,&
    write_rt_dg_hybrid_occupied_checkpoint,read_rt_dg_hybrid_occupied_checkpoint,&
    write_rt_dg_hybrid_ground_state_checkpoint,read_rt_dg_hybrid_ground_state_checkpoint,&
    read_rt_dg_hybrid_ground_state_checkpoint_coalesced,&
    fingerprint_rt_dg_hybrid_ground_state_payload,authenticate_rt_dg_hybrid_ground_state_payload,&
    fingerprint_rt_dg_hybrid_component
  interface
    function c_rename(old_path,new_path) bind(C,name='rename') result(status)
      import::c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
      integer(c_int)::status
    end function c_rename
  end interface
contains
  subroutine fingerprint_rt_dg_hybrid_component(comm,row_ids,values,fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::values(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
#ifdef USE_MPI
    integer::i,j,ierr;integer(int64)::local_hash,bits
    local_hash=0_int64
    if(size(values,1)/=size(row_ids).or..not.finite_matrix(values))then;fingerprint=0_int64;ok=.false.;return;endif
    do i=1,size(row_ids);do j=1,size(values,2)
      bits=transfer(real(values(i,j)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,row_ids(i)),mod(7*j,63)))
      bits=transfer(aimag(values(i,j)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,ishftc(row_ids(i),17)),mod(11*j,63)))
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(fingerprint,ishftc(int(size(values,2),int64),31))
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=ierr==MPI_SUCCESS
#else
    fingerprint=0_int64;ok=.false.
#endif
  end subroutine fingerprint_rt_dg_hybrid_component

  subroutine fingerprint_rt_dg_hybrid_ground_state_payload(comm,payload,fingerprint,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    integer(int64)::common_hash,minimum_common_hash,maximum_common_hash,local_hash,global_hash
    ok=.false.;message='';fingerprint=0_int64
    call validate_ground_state_payload(payload,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid named Hybrid ground-state payload';return
    endif
    call hash_ground_state_common(payload,common_hash)
    call MPI_Allreduce(common_hash,minimum_common_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(common_hash,maximum_common_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_common_hash/=maximum_common_hash)then
      message='rank-disagreeing named Hybrid ground-state metadata';return
    endif
    call validate_ground_state_global_ownership(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid named Hybrid ground-state ownership';return
    endif
    call validate_ground_state_distributed_relations(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid named Hybrid ground-state distributed relation';return
    endif
    call hash_ground_state_payload(payload,local_hash)
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='named Hybrid ground-state digest reduction failed';return;endif
    fingerprint=ground_state_mix_hash(common_hash,global_hash);if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='named Hybrid ground-state fingerprint requires MPI';fingerprint=0_int64
#endif
  end subroutine fingerprint_rt_dg_hybrid_ground_state_payload

  subroutine authenticate_rt_dg_hybrid_ground_state_payload(comm,payload,expected_fingerprint,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::expected_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::observed_fingerprint
    call fingerprint_rt_dg_hybrid_ground_state_payload(comm,payload,observed_fingerprint,ok,message)
    if(.not.ok)return
    if(expected_fingerprint==0_int64.or.observed_fingerprint/=expected_fingerprint.or.&
      payload%payload_fingerprint/=expected_fingerprint)then
      ok=.false.;message='named Hybrid ground-state authentication failed'
    endif
  end subroutine authenticate_rt_dg_hybrid_ground_state_payload

#ifdef USE_MPI
  subroutine validate_ground_state_global_ownership(comm,payload,bad,ierr)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer,intent(out)::bad,ierr
    integer::i,nproc
    integer,allocatable::rows(:),grid(:),rt_rows(:),transformation_rows(:)
    bad=0;ierr=MPI_SUCCESS
    allocate(rows(payload%global_count),grid(payload%global_grid_count),rt_rows(payload%rt_space%rank),&
      transformation_rows(payload%certified_basis%certified_count))
    rows=0;grid=0;rt_rows=0;transformation_rows=0
    do i=1,size(payload%row_ids);rows(int(payload%row_ids(i)))=rows(int(payload%row_ids(i)))+1;enddo
    do i=1,size(payload%grid_ids);grid(int(payload%grid_ids(i)))=grid(int(payload%grid_ids(i)))+1;enddo
    do i=1,size(payload%rt_space%row_ids)
      rt_rows(int(payload%rt_space%row_ids(i)))=rt_rows(int(payload%rt_space%row_ids(i)))+1
    enddo
    do i=1,size(payload%certified_basis%transformation_row_ids)
      transformation_rows(int(payload%certified_basis%transformation_row_ids(i)))=&
        transformation_rows(int(payload%certified_basis%transformation_row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,rows,size(rows),MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,grid,size(grid),MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,rt_rows,size(rt_rows),MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,transformation_rows,size(transformation_rows),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(any(rows/=1).or.any(grid/=1).or.any(rt_rows/=1).or.any(transformation_rows/=1))bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call check_unique(payload%face_ids)
    if(ierr/=MPI_SUCCESS)return
    call check_unique(payload%nonlocal_ids)
  contains
    subroutine check_unique(local_ids)
      integer(int64),intent(in)::local_ids(:)
      integer::destination,j,local_duplicate,local_invalid,global_duplicate,total_received
      integer,allocatable::send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:),next(:)
      integer(int64),allocatable::send_ids(:),receive_ids(:)
      allocate(send_counts(nproc),receive_counts(nproc),send_displacements(nproc),&
        receive_displacements(nproc),next(nproc),send_ids(size(local_ids)))
      send_counts=0;local_invalid=0
      do j=1,size(local_ids)
        if(local_ids(j)<1_int64)then
          local_invalid=1
        else
          destination=int(mod(local_ids(j)-1_int64,int(nproc,int64)))+1
          send_counts(destination)=send_counts(destination)+1
        endif
      enddo
      call MPI_Alltoall(send_counts,1,MPI_INTEGER,receive_counts,1,MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      send_displacements(1)=0;receive_displacements(1)=0
      do j=2,nproc
        send_displacements(j)=send_displacements(j-1)+send_counts(j-1)
        receive_displacements(j)=receive_displacements(j-1)+receive_counts(j-1)
      enddo
      next=send_displacements+1
      do j=1,size(local_ids)
        if(local_ids(j)<1_int64)cycle
        destination=int(mod(local_ids(j)-1_int64,int(nproc,int64)))+1
        send_ids(next(destination))=local_ids(j);next(destination)=next(destination)+1
      enddo
      total_received=sum(receive_counts);allocate(receive_ids(total_received))
      call MPI_Alltoallv(send_ids,send_counts,send_displacements,MPI_INTEGER8,receive_ids,receive_counts,&
        receive_displacements,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call sort_integer64(receive_ids)
      local_duplicate=local_invalid
      do j=2,size(receive_ids)
        if(receive_ids(j)==receive_ids(j-1))local_duplicate=1
      enddo
      call MPI_Allreduce(local_duplicate,global_duplicate,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr==MPI_SUCCESS.and.global_duplicate/=0)bad=1
    end subroutine check_unique

    subroutine sort_integer64(values)
      integer(int64),intent(inout)::values(:)
      integer(int64)::left,middle,right,width,source,destination,other,nvalues
      integer(int64),allocatable::scratch(:)
      nvalues=size(values,kind=int64);if(nvalues<2_int64)return
      allocate(scratch(size(values)));width=1_int64
      do
        scratch=values;left=1_int64
        do while(left<=nvalues)
          middle=min(left+width-1_int64,nvalues);right=min(left+2_int64*width-1_int64,nvalues)
          if(middle<right)then
            source=left;other=middle+1_int64
            do destination=left,right
              if(source>middle)then
                scratch(destination)=values(other);other=other+1_int64
              else if(other>right)then
                scratch(destination)=values(source);source=source+1_int64
              else if(values(source)<=values(other))then
                scratch(destination)=values(source);source=source+1_int64
              else
                scratch(destination)=values(other);other=other+1_int64
              endif
            enddo
          endif
          left=left+2_int64*width
        enddo
        values=scratch
        if(width>=nvalues-width)exit
        width=2_int64*width
      enddo
    end subroutine sort_integer64
  end subroutine validate_ground_state_global_ownership

  subroutine validate_ground_state_distributed_relations(comm,payload,bad,ierr)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer,intent(out)::bad,ierr
    integer::i,row,r,local_bad,global_bad
    real(real64)::safe_limit
    complex(real64),allocatable::local_u(:,:),full_u(:,:),expected_b(:,:),expected_a(:,:)
    r=payload%certified_basis%certified_count;local_bad=0;bad=0;ierr=MPI_SUCCESS
    allocate(local_u(r,r),full_u(r,r));local_u=(0d0,0d0)
    do i=1,size(payload%certified_basis%transformation_row_ids)
      row=int(payload%certified_basis%transformation_row_ids(i))
      local_u(row,:)=payload%certified_basis%u_rt(i,:)
    enddo
    call MPI_Allreduce(local_u,full_u,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    safe_limit=sqrt(huge(0d0))/max(4d0,4d0*real(r,real64))
    if(maximum_complex_component(full_u)>safe_limit.or.&
      maximum_complex_component(payload%certified_basis%c_cert)>safe_limit)then
      local_bad=1
    else
      allocate(expected_b(size(payload%certified_basis%c_cert,1),r),&
        expected_a(r,payload%noccupied))
      expected_b=matmul(payload%certified_basis%c_cert,full_u)
      expected_a=conjg(transpose(full_u(1:payload%noccupied,:)))
      if(maximum_complex_component(payload%certified_basis%b_rt-expected_b)>1d-11.or.&
        maximum_complex_component(payload%certified_basis%initial_occupied_amplitudes-expected_a)>1d-11)&
        local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)bad=global_bad
  end subroutine validate_ground_state_distributed_relations
#endif

  subroutine write_rt_dg_hybrid_checkpoint(comm,path,catalog_fingerprint,metric,operators,coefficients_owned,&
      state_fingerprint,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::catalog_fingerprint,state_fingerprint
    type(s_dg_hybrid_sparse_metric),intent(in)::metric
    type(s_dg_hybrid_sparse_operators),intent(in)::operators
    complex(real64),intent(in)::coefficients_owned(:)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,n,nowned,row,i,k,root,position,unit,io_status,close_status,local_bad,global_bad,allocation_status
    integer::metric_degree,operator_degree,max_metric_degree,max_operator_degree
    integer(int64)::metadata_hash,minimum_metadata_hash,maximum_metadata_hash
    integer,allocatable::ownership(:),owner(:),owner_position(:),metric_degrees(:),operator_degrees(:),&
      metric_columns(:),operator_columns(:)
    complex(real64),allocatable::metric_values(:),operator_hamiltonian(:),operator_position(:,:)
    complex(real64)::coefficient
    logical::file_opened
    character(16)::path_probe
    character(:),allocatable::temporary_path
    ok=.false.;message='';payload_fingerprint=0_int64;n=metric%global_count;nowned=size(metric%owned_row_ids);local_bad=0
    io_status=-1;file_opened=.false.
    temporary_path=trim(path)//'.tmp.'//trim(int64_string(catalog_fingerprint))//'.'//trim(int64_string(state_fingerprint))
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,path_probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent hybrid checkpoint path';return;endif
    if(.not.metric%valid.or..not.operators%valid.or.n<1.or.operators%global_count/=n)local_bad=1
    if(.not.ieee_is_finite(metric%condition_estimate).or.metric%condition_estimate<1d0.or.&
      .not.ieee_is_finite(metric%maximum_value).or.metric%maximum_value<=0d0)local_bad=1
    if(catalog_fingerprint==0_int64.or.state_fingerprint==0_int64.or.metric%fingerprint==0_int64.or.&
      operators%fingerprint==0_int64.or.operators%selection_fingerprint==0_int64.or.&
      operators%window_fingerprint==0_int64.or.operators%packet_fingerprint==0_int64.or.&
      operators%complement_fingerprint==0_int64.or.operators%metric_fingerprint==0_int64.or.&
      operators%position_convention_fingerprint==0_int64)local_bad=1
    if(size(coefficients_owned)/=nowned.or.size(operators%owned_row_ids)/=nowned)local_bad=1
    if(nowned==huge(0))local_bad=1
    if(size(metric%active_rows)/=n.or.size(metric%packet_ids)/=n)local_bad=1
    if(local_bad==0)then
      if(size(metric%row_offsets)/=nowned+1.or.size(operators%row_offsets)/=nowned+1)local_bad=1
    endif
    if(local_bad==0)then
      if(any(metric%owned_row_ids<1_int64).or.any(metric%owned_row_ids>int(n,int64)))local_bad=1
      if(any(metric%owned_row_ids/=operators%owned_row_ids).or..not.finite_vector(coefficients_owned))local_bad=1
      if(metric%row_offsets(1)/=1.or.operators%row_offsets(1)/=1)local_bad=1
      if(any(metric%row_offsets(2:)<metric%row_offsets(:nowned)).or.&
        any(operators%row_offsets(2:)<operators%row_offsets(:nowned)))local_bad=1
      if(metric%row_offsets(nowned+1)-1/=size(metric%column_ids).or.&
        operators%row_offsets(nowned+1)-1/=size(operators%column_ids))local_bad=1
      if(size(metric%values)/=size(metric%column_ids).or.&
        size(operators%hamiltonian_values)/=size(operators%column_ids).or.&
        size(operators%position_values,1)/=3.or.size(operators%position_values,2)/=size(operators%column_ids))local_bad=1
      if(any(metric%column_ids<1).or.any(metric%column_ids>n).or.any(operators%column_ids<1).or.&
        any(operators%column_ids>n))local_bad=1
      if(.not.finite_vector(metric%values).or..not.finite_vector(operators%hamiltonian_values).or.&
        .not.finite_matrix(operators%position_values))local_bad=1
      if(metric%numerical_rank/=count(metric%active_rows).or.&
        .not.valid_packet_activity(metric%packet_ids,metric%active_rows))local_bad=1
      do i=1,nowned
        if(.not.strictly_increasing(metric%column_ids(metric%row_offsets(i):metric%row_offsets(i+1)-1)).or.&
          .not.strictly_increasing(operators%column_ids(operators%row_offsets(i):operators%row_offsets(i+1)-1)))local_bad=1
      enddo
    endif
    if(local_bad==0)then
      if(any(metric%packet_ids<1))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid checkpoint write contract';return;endif
    metadata_hash=catalog_fingerprint
    metadata_hash=mix_hash(metadata_hash,state_fingerprint);metadata_hash=mix_hash(metadata_hash,metric%fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%fingerprint);metadata_hash=mix_hash(metadata_hash,int(n,int64))
    metadata_hash=mix_hash(metadata_hash,int(metric%numerical_rank,int64))
    metadata_hash=mix_hash(metadata_hash,transfer(metric%condition_estimate,metadata_hash))
    metadata_hash=mix_hash(metadata_hash,transfer(metric%maximum_value,metadata_hash))
    metadata_hash=mix_hash(metadata_hash,operators%selection_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%window_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%packet_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%complement_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%metric_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%position_convention_fingerprint)
    do row=1,n
      metadata_hash=mix_hash(metadata_hash,merge(1_int64,0_int64,metric%active_rows(row)))
      metadata_hash=mix_hash(metadata_hash,int(metric%packet_ids(row),int64))
    enddo
    call MPI_Allreduce(metadata_hash,minimum_metadata_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(metadata_hash,maximum_metadata_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_metadata_hash/=maximum_metadata_hash)then
      message='rank-disagreeing hybrid checkpoint metadata';return
    endif
    allocate(ownership(n),owner(n),owner_position(n),metric_degrees(n),operator_degrees(n),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid checkpoint metadata';return;endif
    ownership=0;owner=-1;owner_position=0;metric_degrees=0;operator_degrees=0
    do i=1,nowned
      row=int(metric%owned_row_ids(i));ownership(row)=ownership(row)+1;owner(row)=rank;owner_position(row)=i
      metric_degrees(row)=metric%row_offsets(i+1)-metric%row_offsets(i)
      operator_degrees(row)=operators%row_offsets(i+1)-operators%row_offsets(i)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owner,n,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,n,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,metric_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,operator_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    if(any(ownership/=1).or.any(metric_degrees<0).or.any(operator_degrees<0))local_bad=1
    max_metric_degree=maxval(metric_degrees);max_operator_degree=maxval(operator_degrees)
    if(max_operator_degree>huge(0)/3)local_bad=1
    allocate(metric_columns(max(1,max_metric_degree)),metric_values(max(1,max_metric_degree)),&
      operator_columns(max(1,max_operator_degree)),&
      operator_hamiltonian(max(1,max_operator_degree)),operator_position(3,max(1,max_operator_degree)),stat=allocation_status)
    if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='invalid or unallocatable checkpoint row workspace';return;endif
    io_status=0
    if(rank==0)then
      open(newunit=unit,file=temporary_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status)
      file_opened=io_status==0
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='cannot open hybrid checkpoint for writing';return;endif
    payload_fingerprint=catalog_fingerprint;call hash_int(state_fingerprint);call hash_int(metric%fingerprint)
    call hash_int(operators%fingerprint);call hash_int(int(n,int64));call hash_int(int(metric%numerical_rank,int64))
    call hash_int(transfer(metric%condition_estimate,payload_fingerprint));call hash_int(transfer(metric%maximum_value,payload_fingerprint))
    call hash_int(operators%selection_fingerprint);call hash_int(operators%window_fingerprint)
    call hash_int(operators%packet_fingerprint);call hash_int(operators%complement_fingerprint)
    call hash_int(operators%metric_fingerprint);call hash_int(operators%position_convention_fingerprint)
    do row=1,n
      call hash_int(merge(1_int64,0_int64,metric%active_rows(row)))
      call hash_int(int(metric%packet_ids(row),int64));call hash_int(int(metric_degrees(row),int64))
      call hash_int(int(operator_degrees(row),int64))
    enddo
    if(rank==0)then
      write(unit,iostat=io_status)checkpoint_magic,checkpoint_version,n,metric%numerical_rank,catalog_fingerprint,&
        state_fingerprint,metric%fingerprint,metric%condition_estimate,metric%maximum_value,&
        operators%selection_fingerprint,operators%window_fingerprint,operators%packet_fingerprint,&
        operators%complement_fingerprint,operators%metric_fingerprint,operators%position_convention_fingerprint,&
        operators%fingerprint,metric%active_rows,metric%packet_ids,metric_degrees,operator_degrees
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do row=1,n
      root=owner(row);position=owner_position(row);metric_degree=metric_degrees(row);operator_degree=operator_degrees(row)
      if(rank==root)then
        metric_columns(1:metric_degree)=metric%column_ids(metric%row_offsets(position):metric%row_offsets(position+1)-1)
        metric_values(1:metric_degree)=metric%values(metric%row_offsets(position):metric%row_offsets(position+1)-1)
        operator_columns(1:operator_degree)=operators%column_ids(operators%row_offsets(position):operators%row_offsets(position+1)-1)
        operator_hamiltonian(1:operator_degree)=&
          operators%hamiltonian_values(operators%row_offsets(position):operators%row_offsets(position+1)-1)
        operator_position(:,1:operator_degree)=&
          operators%position_values(:,operators%row_offsets(position):operators%row_offsets(position+1)-1)
        coefficient=coefficients_owned(position)
      endif
      call MPI_Bcast(metric_columns,metric_degree,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(metric_values,metric_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_columns,operator_degree,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_hamiltonian,operator_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_position,3*operator_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(coefficient,1,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call hash_row(row,metric_degree,operator_degree,metric_columns,metric_values,operator_columns,&
        operator_hamiltonian,operator_position,coefficient)
      if(rank==0)write(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
        operator_columns(1:operator_degree),operator_hamiltonian(1:operator_degree),&
        operator_position(:,1:operator_degree),coefficient
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)then
      write(unit,iostat=io_status)payload_fingerprint
      close_status=0;close(unit,iostat=close_status);file_opened=.false.
      if(io_status==0)io_status=close_status
    endif
    call sync_io(io_status,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='hybrid checkpoint final write failed';return;endif
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='hybrid checkpoint atomic publication failed';return;endif
    call cleanup();ok=.true.;return
900 message='hybrid checkpoint MPI stream failed';if(rank==0.and.file_opened)close(unit);call cleanup();return
910 message='hybrid checkpoint file write failed';if(rank==0.and.file_opened)close(unit);call cleanup();return
#else
    ok=.false.;message='hybrid checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine hash_int(value)
      integer(int64),intent(in)::value;payload_fingerprint=ieor(ishftc(payload_fingerprint,9),value)
    end subroutine hash_int
    subroutine hash_complex(value)
      complex(real64),intent(in)::value;call hash_int(transfer(real(value),payload_fingerprint));&
        call hash_int(transfer(aimag(value),payload_fingerprint))
    end subroutine hash_complex
    subroutine hash_row(global_row,md,od,mc,mv,oc,oh,op,c)
      integer,intent(in)::global_row,md,od,mc(:),oc(:);complex(real64),intent(in)::mv(:),oh(:),op(:,:),c
      integer::a,b
      call hash_int(int(global_row,int64));call hash_int(int(md,int64));call hash_int(int(od,int64))
      do a=1,md;call hash_int(int(mc(a),int64));call hash_complex(mv(a));enddo
      do a=1,od
        call hash_int(int(oc(a),int64));call hash_complex(oh(a))
        do b=1,3;call hash_complex(op(b,a));enddo
      enddo
      call hash_complex(c)
    end subroutine hash_row
    subroutine cleanup()
      if(allocated(ownership))deallocate(ownership);if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position);if(allocated(metric_degrees))deallocate(metric_degrees)
      if(allocated(operator_degrees))deallocate(operator_degrees);if(allocated(metric_columns))deallocate(metric_columns)
      if(allocated(metric_values))deallocate(metric_values);if(allocated(operator_columns))deallocate(operator_columns)
      if(allocated(operator_hamiltonian))deallocate(operator_hamiltonian)
      if(allocated(operator_position))deallocate(operator_position)
    end subroutine cleanup
#endif
  end subroutine write_rt_dg_hybrid_checkpoint

  subroutine read_rt_dg_hybrid_checkpoint(comm,path,expected_catalog,expected_state,expected_selection,expected_window,&
      expected_packet,expected_complement,expected_metric,expected_position,expected_operator,metric,operators,&
      coefficients_owned,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
      expected_complement,expected_metric,expected_position,expected_operator
    type(s_dg_hybrid_sparse_metric),intent(out)::metric
    type(s_dg_hybrid_sparse_operators),intent(out)::operators
    complex(real64),allocatable,intent(out)::coefficients_owned(:)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,version,n,numerical_rank,row,i,k,nowned,metric_nnz,operator_nnz,&
      metric_degree,operator_degree,max_metric_degree,max_operator_degree,allocation_status,position
    integer::integer_header(3),allocation_bad
    integer(int64)::fingerprint_header(3),operator_header(7)
    integer(int64)::catalog,state_fp,metric_fp,selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,&
      position_fp,operator_fp,stored_fingerprint
    integer(int64)::file_size
    integer(int64)::expected_values(9),minimum_expected(9),maximum_expected(9)
    real(real64)::condition,maximum_value
    logical,allocatable::active_rows(:)
    integer,allocatable::packet_ids(:),metric_degrees(:),operator_degrees(:),metric_columns(:),operator_columns(:)
    complex(real64),allocatable::metric_values(:),operator_metric(:),operator_hamiltonian(:),operator_position(:,:)
    complex(real64)::coefficient
    logical::file_opened
    character(16)::magic,path_probe
    ok=.false.;message='';payload_fingerprint=0_int64;file_opened=.false.
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    expected_values=[expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
      expected_complement,expected_metric,expected_position,expected_operator]
    call MPI_Allreduce(expected_values,minimum_expected,9,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(expected_values,maximum_expected,9,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(minimum_expected/=maximum_expected).or.any(expected_values==0_int64))then
      message='rank-disagreeing expected checkpoint provenance';return
    endif
    call validate_path(path,path_probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent hybrid restart path';return;endif
    io_status=0;file_size=0_int64
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status)
      file_opened=io_status==0
      if(io_status==0)inquire(unit=unit,size=file_size,iostat=io_status)
      if(io_status==0)read(unit,iostat=io_status)magic,version,n,numerical_rank,catalog,state_fp,metric_fp,condition,maximum_value,&
        selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,position_fp,operator_fp
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then
      if(rank==0.and.file_opened)close(unit);message='cannot read hybrid checkpoint header';return
    endif
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(magic,len(magic),MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)then
      integer_header=[version,n,numerical_rank]
      fingerprint_header=[catalog,state_fp,metric_fp]
      operator_header=[selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,position_fp,operator_fp]
    endif
    call MPI_Bcast(integer_header,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(fingerprint_header,3,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(condition,1,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(maximum_value,1,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(operator_header,7,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    version=integer_header(1);n=integer_header(2);numerical_rank=integer_header(3)
    catalog=fingerprint_header(1);state_fp=fingerprint_header(2);metric_fp=fingerprint_header(3)
    selection_fp=operator_header(1);window_fp=operator_header(2);packet_fp=operator_header(3)
    complement_fp=operator_header(4);operator_metric_fp=operator_header(5);position_fp=operator_header(6);operator_fp=operator_header(7)
    if(magic/=checkpoint_magic.or.(version/=checkpoint_version.and.version/=legacy_checkpoint_version).or.&
      n<1.or.numerical_rank<1.or.numerical_rank>n.or.&
      n>huge(0)/3.or.int(n,int64)>file_size/4_int64)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif;message='incompatible hybrid checkpoint version';return
    endif
    if(catalog/=expected_catalog.or.state_fp/=expected_state.or.selection_fp/=expected_selection.or.&
      window_fp/=expected_window.or.packet_fp/=expected_packet.or.complement_fp/=expected_complement.or.&
      metric_fp/=expected_metric.or.position_fp/=expected_position.or.operator_fp/=expected_operator.or.&
      operator_metric_fp/=metric_fp)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif;message='stale hybrid checkpoint provenance';return
    endif
    allocate(active_rows(n),packet_ids(n),metric_degrees(n),operator_degrees(n),stat=allocation_status)
    if(rank==0.and.allocation_status==0)read(unit,iostat=io_status)active_rows,packet_ids,metric_degrees,operator_degrees
    allocation_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(allocation_bad,k,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.k/=0)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif
      call cleanup_buffers();message='cannot allocate checkpoint metadata';return
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    call MPI_Bcast(active_rows,n,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(packet_ids,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(metric_degrees,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(operator_degrees,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    if(any(packet_ids<1).or.any(packet_ids>n).or.any(metric_degrees<0).or.any(operator_degrees<0).or.n>huge(0)/3.or.&
      numerical_rank/=count(active_rows).or..not.valid_packet_activity(packet_ids,active_rows).or.&
      .not.ieee_is_finite(condition).or.condition<1d0.or..not.ieee_is_finite(maximum_value).or.maximum_value<=0d0)goto 920
    if(rank>=n)then;nowned=0;else;nowned=(n-1-rank)/nproc+1;endif
    metric_nnz=0;operator_nnz=0
    do row=1,n
      if(metric_degrees(row)>n.or.operator_degrees(row)>n)goto 920
      if(mod(row-1,nproc)==rank)then
        if(metric_degrees(row)>huge(0)-metric_nnz.or.operator_degrees(row)>huge(0)-operator_nnz)goto 920
        metric_nnz=metric_nnz+metric_degrees(row);operator_nnz=operator_nnz+operator_degrees(row)
      endif
    enddo
    max_metric_degree=maxval(metric_degrees);max_operator_degree=maxval(operator_degrees)
    allocate(metric_columns(max(1,max_metric_degree)),metric_values(max(1,max_metric_degree)),&
      operator_columns(max(1,max_operator_degree)),operator_metric(max(1,max_operator_degree)),&
      operator_hamiltonian(max(1,max_operator_degree)),operator_position(3,max(1,max_operator_degree)),&
      metric%owned_row_ids(nowned),metric%row_offsets(nowned+1),metric%column_ids(metric_nnz),metric%values(metric_nnz),&
      metric%active_rows(n),metric%packet_ids(n),operators%owned_row_ids(nowned),operators%row_offsets(nowned+1),&
      operators%column_ids(operator_nnz),operators%metric_values(operator_nnz),operators%hamiltonian_values(operator_nnz),&
      operators%position_values(3,operator_nnz),coefficients_owned(nowned),stat=allocation_status)
    allocation_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(allocation_bad,k,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.k/=0)goto 920
    metric%row_offsets(1)=1;operators%row_offsets(1)=1;position=0;metric_nnz=0;operator_nnz=0
    payload_fingerprint=catalog;call hash_int_read(state_fp);call hash_int_read(metric_fp);call hash_int_read(operator_fp)
    call hash_int_read(int(n,int64));call hash_int_read(int(numerical_rank,int64))
    call hash_int_read(transfer(condition,payload_fingerprint));call hash_int_read(transfer(maximum_value,payload_fingerprint))
    call hash_int_read(selection_fp);call hash_int_read(window_fp);call hash_int_read(packet_fp);call hash_int_read(complement_fp)
    call hash_int_read(operator_metric_fp);call hash_int_read(position_fp)
    do row=1,n
      call hash_int_read(merge(1_int64,0_int64,active_rows(row)))
      call hash_int_read(int(packet_ids(row),int64));call hash_int_read(int(metric_degrees(row),int64))
      call hash_int_read(int(operator_degrees(row),int64))
    enddo
    do row=1,n
      metric_degree=metric_degrees(row);operator_degree=operator_degrees(row)
      if(rank==0)then
        if(version==legacy_checkpoint_version)then
          read(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
            operator_columns(1:operator_degree),operator_metric(1:operator_degree),operator_hamiltonian(1:operator_degree),&
            operator_position(:,1:operator_degree),coefficient
        else
          read(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
            operator_columns(1:operator_degree),operator_hamiltonian(1:operator_degree),&
            operator_position(:,1:operator_degree),coefficient
        endif
      endif
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
      call MPI_Bcast(metric_columns,metric_degree,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(metric_values,metric_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(operator_columns,operator_degree,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      if(version==legacy_checkpoint_version)then
        call MPI_Bcast(operator_metric,operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      endif
      call MPI_Bcast(operator_hamiltonian,operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(operator_position,3*operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(coefficient,1,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      if(any(metric_columns(1:metric_degree)<1).or.any(metric_columns(1:metric_degree)>n).or.&
        any(operator_columns(1:operator_degree)<1).or.any(operator_columns(1:operator_degree)>n))goto 920
      if(.not.strictly_increasing(metric_columns(1:metric_degree)).or.&
        .not.strictly_increasing(operator_columns(1:operator_degree)))goto 920
      if(.not.finite_vector(metric_values(1:metric_degree)).or.&
        .not.finite_vector(operator_hamiltonian(1:operator_degree)).or.&
        .not.finite_matrix(operator_position(:,1:operator_degree)).or..not.finite_vector([coefficient]))goto 920
      if(version==legacy_checkpoint_version)then
        if(.not.finite_vector(operator_metric(1:operator_degree)))goto 920
      else
        operator_metric(1:operator_degree)=(0d0,0d0)
        do i=1,operator_degree
          k=findloc(metric_columns(1:metric_degree),operator_columns(i),dim=1)
          if(k>0)operator_metric(i)=metric_values(k)
        enddo
      endif
      call hash_row_read(row,metric_degree,operator_degree,metric_columns,metric_values,operator_columns,operator_metric,&
        operator_hamiltonian,operator_position,coefficient)
      if(mod(row-1,nproc)==rank)then
        position=position+1;metric%owned_row_ids(position)=row;operators%owned_row_ids(position)=row
        metric%column_ids(metric_nnz+1:metric_nnz+metric_degree)=metric_columns(1:metric_degree)
        metric%values(metric_nnz+1:metric_nnz+metric_degree)=metric_values(1:metric_degree);metric_nnz=metric_nnz+metric_degree
        metric%row_offsets(position+1)=metric_nnz+1
        operators%column_ids(operator_nnz+1:operator_nnz+operator_degree)=operator_columns(1:operator_degree)
        operators%metric_values(operator_nnz+1:operator_nnz+operator_degree)=operator_metric(1:operator_degree)
        operators%hamiltonian_values(operator_nnz+1:operator_nnz+operator_degree)=operator_hamiltonian(1:operator_degree)
        operators%position_values(:,operator_nnz+1:operator_nnz+operator_degree)=operator_position(:,1:operator_degree)
        operator_nnz=operator_nnz+operator_degree;operators%row_offsets(position+1)=operator_nnz+1;coefficients_owned(position)=coefficient
      endif
    enddo
    if(rank==0)then;read(unit,iostat=io_status)stored_fingerprint;close(unit);file_opened=.false.;endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    call MPI_Bcast(stored_fingerprint,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(stored_fingerprint/=payload_fingerprint)then;call cleanup_read();message='corrupt hybrid checkpoint payload';return;endif
    metric%valid=.true.;metric%global_count=n;metric%numerical_rank=numerical_rank;metric%condition_estimate=condition
    metric%maximum_value=maximum_value;metric%fingerprint=metric_fp;metric%active_rows=active_rows;metric%packet_ids=packet_ids
    metric%max_row_nnz=max_metric_degree
    operators%valid=.true.;operators%global_count=n;operators%selection_fingerprint=selection_fp
    operators%window_fingerprint=window_fp;operators%packet_fingerprint=packet_fp
    operators%complement_fingerprint=complement_fp
    operators%metric_fingerprint=operator_metric_fp;operators%position_convention_fingerprint=position_fp
    operators%fingerprint=operator_fp;ok=.true.;call cleanup_buffers();return
920 if(rank==0.and.file_opened)close(unit);call cleanup_read();message='hybrid checkpoint read failed';return
#else
    ok=.false.;message='hybrid checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine hash_int_read(value)
      integer(int64),intent(in)::value;payload_fingerprint=ieor(ishftc(payload_fingerprint,9),value)
    end subroutine hash_int_read
    subroutine hash_complex_read(value)
      complex(real64),intent(in)::value;call hash_int_read(transfer(real(value),payload_fingerprint));&
        call hash_int_read(transfer(aimag(value),payload_fingerprint))
    end subroutine hash_complex_read
    subroutine hash_row_read(global_row,md,od,mc,mv,oc,om,oh,op,c)
      integer,intent(in)::global_row,md,od,mc(:),oc(:);complex(real64),intent(in)::mv(:),om(:),oh(:),op(:,:),c
      integer::a,b
      call hash_int_read(int(global_row,int64));call hash_int_read(int(md,int64));call hash_int_read(int(od,int64))
      do a=1,md;call hash_int_read(int(mc(a),int64));call hash_complex_read(mv(a));enddo
      do a=1,od
        call hash_int_read(int(oc(a),int64))
        if(version==legacy_checkpoint_version)call hash_complex_read(om(a))
        call hash_complex_read(oh(a))
        do b=1,3;call hash_complex_read(op(b,a));enddo
      enddo
      call hash_complex_read(c)
    end subroutine hash_row_read
    subroutine cleanup_buffers()
      if(allocated(active_rows))deallocate(active_rows);if(allocated(packet_ids))deallocate(packet_ids)
      if(allocated(metric_degrees))deallocate(metric_degrees);if(allocated(operator_degrees))deallocate(operator_degrees)
      if(allocated(metric_columns))deallocate(metric_columns);if(allocated(metric_values))deallocate(metric_values)
      if(allocated(operator_columns))deallocate(operator_columns);if(allocated(operator_metric))deallocate(operator_metric)
      if(allocated(operator_hamiltonian))deallocate(operator_hamiltonian);if(allocated(operator_position))deallocate(operator_position)
    end subroutine cleanup_buffers
    subroutine cleanup_read()
      call cleanup_buffers()
      if(allocated(coefficients_owned))deallocate(coefficients_owned)
      if(allocated(metric%owned_row_ids))deallocate(metric%owned_row_ids)
      if(allocated(metric%row_offsets))deallocate(metric%row_offsets)
      if(allocated(metric%column_ids))deallocate(metric%column_ids)
      if(allocated(metric%values))deallocate(metric%values)
      if(allocated(metric%active_rows))deallocate(metric%active_rows)
      if(allocated(metric%packet_ids))deallocate(metric%packet_ids)
      if(allocated(operators%owned_row_ids))deallocate(operators%owned_row_ids)
      if(allocated(operators%row_offsets))deallocate(operators%row_offsets)
      if(allocated(operators%column_ids))deallocate(operators%column_ids)
      if(allocated(operators%metric_values))deallocate(operators%metric_values)
      if(allocated(operators%hamiltonian_values))deallocate(operators%hamiltonian_values)
      if(allocated(operators%position_values))deallocate(operators%position_values)
      metric%valid=.false.;operators%valid=.false.
    end subroutine cleanup_read
#endif
  end subroutine read_rt_dg_hybrid_checkpoint

  subroutine write_rt_dg_hybrid_ground_state_checkpoint(comm,path,payload,payload_fingerprint,ok,message,interrupt_after_write)
    integer,intent(in)::comm
    character(*),intent(in)::path
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,optional,intent(in)::interrupt_after_write
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,owner,local_bad,global_bad,&
      header_i(ground_state_integer_count),attempt
    integer(int64)::header_fp(ground_state_fingerprint_count),local_hash,global_hash,common_hash,&
      minimum_common_hash,maximum_common_hash,nonce,&
      computed_component_fingerprints(4)
    real(real64)::header_r(ground_state_real_count)
    logical::header_l(ground_state_logical_count),opened,verified_ok,fingerprint_ok,created
    character(256)::verified_message
    character(:),allocatable::temporary_path
    type(s_rt_dg_hybrid_ground_state_payload)::verified
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_ground_state_payload(payload,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid complete DG ground-state payload';return;endif
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%kinetic_rows,&
      computed_component_fingerprints(1),fingerprint_ok);if(.not.fingerprint_ok)return
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%nonlocal_rows,&
      computed_component_fingerprints(2),fingerprint_ok);if(.not.fingerprint_ok)return
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%local_rows,&
      computed_component_fingerprints(3),fingerprint_ok);if(.not.fingerprint_ok)return
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%sipg_rows,&
      computed_component_fingerprints(4),fingerprint_ok);if(.not.fingerprint_ok)return
    if(any(computed_component_fingerprints/=[payload%kinetic_fingerprint,payload%nonlocal_fingerprint,&
      payload%local_fingerprint,payload%sipg_fingerprint]))then
      message='complete DG ground-state component fingerprint mismatch';return
    endif
    call hash_ground_state_common(payload,common_hash)
    call MPI_Allreduce(common_hash,minimum_common_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(common_hash,maximum_common_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_common_hash/=maximum_common_hash)then
      message='rank-disagreeing complete DG ground-state metadata';return
    endif
    call validate_ground_state_global_ownership(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid complete DG ground-state ownership';return
    endif
    call validate_ground_state_distributed_relations(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid complete DG ground-state distributed relation';return
    endif
    call hash_ground_state_payload(payload,local_hash)
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='complete DG ground-state fingerprint failed';return;endif
    global_hash=ground_state_mix_hash(common_hash,global_hash);if(global_hash==0_int64)global_hash=1_int64
    payload_fingerprint=global_hash
    nonce=0_int64;if(rank==0)call system_clock(count=nonce)
    call MPI_Bcast(nonce,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='complete DG ground-state temporary-name broadcast failed';return;endif
    created=.false.
    do attempt=0,31
      temporary_path=trim(path)//'.tmp.'//trim(int64_string(payload%catalog_fingerprint))//'.'//&
        trim(int64_string(payload%state_fingerprint))//'.'//trim(int64_string(nonce+int(attempt,int64)))
      io_status=0
      if(rank==0)then
        open(newunit=unit,file=temporary_path,status='new',access='stream',form='unformatted',action='write',iostat=io_status)
        opened=io_status==0
      endif
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(io_status==0)then;created=.true.;exit;endif
    enddo
    if(.not.created)then;message='cannot create unique complete DG ground-state temporary file';return;endif
    if(rank==0)write(unit,iostat=io_status)ground_state_magic,ground_state_version,nproc
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do owner=0,nproc-1
      if(rank==owner)then
        call pack_ground_state_header(payload,header_l,header_i,header_fp,header_r)
        header_fp(20)=payload_fingerprint
      endif
      call MPI_Bcast(header_l,ground_state_logical_count,MPI_LOGICAL,owner,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(header_i,ground_state_integer_count,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(header_fp,ground_state_fingerprint_count,MPI_INTEGER8,owner,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(header_r,ground_state_real_count,MPI_DOUBLE_PRECISION,owner,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      if(rank==0)write(unit,iostat=io_status)header_l,header_i,header_fp,header_r
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
      call write_ground_state_arrays(comm,owner,unit,rank,payload,io_status,ierr)
      if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(rank==0)then;close(unit,iostat=io_status);opened=.false.;endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    if(present(interrupt_after_write))then
      if(interrupt_after_write)then;message='injected interruption after complete checkpoint write';return;endif
    endif
    call read_rt_dg_hybrid_ground_state_checkpoint(comm,temporary_path,verified,global_hash,verified_ok,verified_message)
    if(.not.verified_ok.or.global_hash/=payload_fingerprint)then
      message='complete DG ground-state checkpoint verification failed';return
    endif
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;message='complete DG ground-state checkpoint publication failed';return;endif
    ok=.true.;return
900 message='complete DG ground-state checkpoint MPI stream failed';if(rank==0.and.opened)close(unit);return
910 message='complete DG ground-state checkpoint write failed';if(rank==0.and.opened)close(unit);return
#else
    ok=.false.;message='complete DG ground-state checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  end subroutine write_rt_dg_hybrid_ground_state_checkpoint

  subroutine read_rt_dg_hybrid_ground_state_checkpoint(comm,path,payload,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::payload
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,owner,file_nproc,version,local_bad,global_bad,&
      header_i(ground_state_integer_count),i
    integer,allocatable::grid_ownership(:)
    integer(int64)::header_fp(ground_state_fingerprint_count),local_hash,global_hash,&
      computed_component_fingerprints(4),common_hash,&
      minimum_common_hash,maximum_common_hash
    real(real64)::header_r(ground_state_real_count)
    logical::header_l(ground_state_logical_count),opened,fingerprint_ok
    character(16)::magic
    payload=s_rt_dg_hybrid_ground_state_payload()
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status)
      opened=io_status==0
      if(io_status==0)read(unit,iostat=io_status)magic,version,file_nproc
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    call MPI_Bcast(magic,16,MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(version,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(file_nproc,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    if(magic/=ground_state_magic.or.version/=ground_state_version.or.file_nproc/=nproc)then
      message='incompatible complete DG ground-state checkpoint';goto 920
    endif
    do owner=0,nproc-1
      if(rank==0)read(unit,iostat=io_status)header_l,header_i,header_fp,header_r
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
      call MPI_Bcast(header_l,ground_state_logical_count,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(header_i,ground_state_integer_count,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(header_fp,ground_state_fingerprint_count,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(header_r,ground_state_real_count,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      if(rank==owner)then
        call set_ground_state_header(payload,header_l,header_i,header_fp,header_r)
      endif
      call read_ground_state_arrays(comm,owner,unit,rank,payload,io_status,ierr)
      if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    enddo
    if(rank==0)then;close(unit);opened=.false.;endif
    call validate_ground_state_payload(payload,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='corrupt complete DG ground-state payload';goto 920;endif
    call hash_ground_state_common(payload,common_hash)
    call MPI_Allreduce(common_hash,minimum_common_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(common_hash,maximum_common_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_common_hash/=maximum_common_hash)then
      message='rank-disagreeing complete DG ground-state metadata';goto 920
    endif
    allocate(grid_ownership(payload%global_grid_count));grid_ownership=0
    do i=1,size(payload%grid_ids)
      grid_ownership(int(payload%grid_ids(i)))=grid_ownership(int(payload%grid_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,grid_ownership,payload%global_grid_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(grid_ownership/=1))then
      message='corrupt complete DG ground-state grid catalog';goto 920
    endif
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%kinetic_rows,&
      computed_component_fingerprints(1),fingerprint_ok);if(.not.fingerprint_ok)goto 920
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%nonlocal_rows,&
      computed_component_fingerprints(2),fingerprint_ok);if(.not.fingerprint_ok)goto 920
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%local_rows,&
      computed_component_fingerprints(3),fingerprint_ok);if(.not.fingerprint_ok)goto 920
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%sipg_rows,&
      computed_component_fingerprints(4),fingerprint_ok);if(.not.fingerprint_ok)goto 920
    if(any(computed_component_fingerprints/=[payload%kinetic_fingerprint,payload%nonlocal_fingerprint,&
      payload%local_fingerprint,payload%sipg_fingerprint]))then
      message='corrupt complete DG ground-state component fingerprint';goto 920
    endif
    call validate_ground_state_global_ownership(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='corrupt complete DG ground-state ownership';goto 920
    endif
    call validate_ground_state_distributed_relations(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='corrupt complete DG ground-state distributed relation';goto 920
    endif
    call hash_ground_state_payload(payload,local_hash)
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 920
    global_hash=ground_state_mix_hash(common_hash,global_hash);if(global_hash==0_int64)global_hash=1_int64
    local_bad=merge(0,1,global_hash==payload%payload_fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='corrupt complete DG ground-state fingerprint';goto 920;endif
    payload_fingerprint=global_hash;ok=.true.;return
920 if(rank==0.and.opened)close(unit);payload=s_rt_dg_hybrid_ground_state_payload()
    payload_fingerprint=0_int64
    if(len_trim(message)==0)message='complete DG ground-state checkpoint read failed'
    return
#else
    ok=.false.;message='complete DG ground-state checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  end subroutine read_rt_dg_hybrid_ground_state_checkpoint

  subroutine read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,path,payload,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    type(s_rt_dg_hybrid_ground_state_payload),intent(out)::payload
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_rt_dg_hybrid_ground_state_payload)::shard
    integer::rank,nproc,ierr,unit,io_status,owner,receiver,file_nproc,version,bad,global_bad,&
      header_i(ground_state_integer_count)
    integer(int64)::header_fp(ground_state_fingerprint_count),local_hash,global_hash,&
      shard_common_hash,reference_common_hash,assembled_common_hash,min_common,max_common,&
      shard_payload_fingerprint,reference_payload_fingerprint
    real(real64)::header_r(ground_state_real_count)
    logical::header_l(ground_state_logical_count),opened,have_common
    character(16)::magic
    payload=s_rt_dg_hybrid_ground_state_payload()
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=0
    local_hash=0_int64;have_common=.false.;reference_common_hash=0_int64;reference_payload_fingerprint=0_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status)
      opened=io_status==0;if(io_status==0)read(unit,iostat=io_status)magic,version,file_nproc
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 930
    call MPI_Bcast(magic,16,MPI_CHARACTER,0,comm,ierr);call MPI_Bcast(version,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(file_nproc,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.magic/=ground_state_magic.or.version/=ground_state_version.or.file_nproc<1)goto 930
    do owner=0,file_nproc-1
      receiver=mod(owner,nproc)
      if(rank==0)read(unit,iostat=io_status)header_l,header_i,header_fp,header_r
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 930
      call MPI_Bcast(header_l,ground_state_logical_count,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      call MPI_Bcast(header_i,ground_state_integer_count,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      call MPI_Bcast(header_fp,ground_state_fingerprint_count,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      call MPI_Bcast(header_r,ground_state_real_count,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      if(rank==receiver)then
        shard=s_rt_dg_hybrid_ground_state_payload()
        call set_ground_state_header(shard,header_l,header_i,header_fp,header_r)
      endif
      call read_ground_state_arrays(comm,receiver,unit,rank,shard,io_status,ierr)
      if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 930
      bad=0
      if(rank==receiver)then
        call validate_ground_state_payload(shard,bad)
        if(bad==0)then
          call hash_ground_state_common(shard,shard_common_hash)
          shard_payload_fingerprint=shard%payload_fingerprint
        endif
      endif
      call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
      call MPI_Bcast(shard_common_hash,1,MPI_INTEGER8,receiver,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      call MPI_Bcast(shard_payload_fingerprint,1,MPI_INTEGER8,receiver,comm,ierr);if(ierr/=MPI_SUCCESS)goto 930
      bad=0
      if(owner==0)then
        reference_common_hash=shard_common_hash
        reference_payload_fingerprint=shard_payload_fingerprint
      else if(shard_common_hash/=reference_common_hash.or.&
          shard_payload_fingerprint/=reference_payload_fingerprint)then
        bad=1
      endif
      call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
      if(rank==receiver)then
        if(.not.have_common)then
          call copy_ground_state_common(shard,payload);have_common=.true.
        endif
        call append_initialization_shard(payload,shard)
        shard=s_rt_dg_hybrid_ground_state_payload()
      endif
    enddo
    call broadcast_coalesced_ground_state_common(comm,rank,payload,have_common,ierr)
    if(ierr/=MPI_SUCCESS)goto 930
    if(file_nproc/=nproc)then
      call redistribute_coalesced_ground_state(comm,rank,nproc,payload,ierr)
      if(ierr/=MPI_SUCCESS)goto 930
    endif
    if(rank==0)then;close(unit);opened=.false.;endif
    call validate_ground_state_payload(payload,bad)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    call hash_ground_state_common(payload,assembled_common_hash)
    call MPI_Allreduce(assembled_common_hash,min_common,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(assembled_common_hash,max_common,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.min_common/=max_common.or.&
      assembled_common_hash/=reference_common_hash)goto 930
    call validate_ground_state_global_ownership(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    call validate_ground_state_distributed_relations(comm,payload,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    call hash_ground_state_payload(payload,local_hash)
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 930
    global_hash=ground_state_mix_hash(assembled_common_hash,global_hash);if(global_hash==0_int64)global_hash=1_int64
    bad=merge(0,1,global_hash==reference_payload_fingerprint.and.&
      payload%payload_fingerprint==reference_payload_fingerprint)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    payload_fingerprint=global_hash;ok=.true.;return
930 if(rank==0.and.opened)close(unit);payload=s_rt_dg_hybrid_ground_state_payload()
    payload_fingerprint=0_int64
    if(len_trim(message)==0)message='complete DG ground-state checkpoint coalescing failed'
    return
#else
    ok=.false.;message='complete DG ground-state checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  end subroutine read_rt_dg_hybrid_ground_state_checkpoint_coalesced

  subroutine pack_ground_state_header(p,l,h,f,r)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::p
    logical,intent(out)::l(ground_state_logical_count)
    integer,intent(out)::h(ground_state_integer_count)
    integer(int64),intent(out)::f(ground_state_fingerprint_count)
    real(real64),intent(out)::r(ground_state_real_count)
    l=[p%valid,p%final_refresh_complete,p%analysis_complete,p%identity_only,&
      p%construction_catalog%valid,p%certified_basis%valid,p%certified_basis%localization_converged,&
      p%certified_basis%localization_symmetry_constrained,p%electron_count%valid,p%rt_space%valid,&
      p%energy_window%valid,p%energy_window%compatibility_dynamic_rank,p%energy_window%proof_state_present,&
      p%symmetry_receipt%valid,p%handoff_receipts%valid]
    h=[p%global_count,p%global_grid_count,p%noccupied,p%operation_count,p%nonidentity_operation_count,&
      p%construction_catalog%global_count,p%certified_basis%construction_count,&
      p%certified_basis%certified_count,p%certified_basis%occupied_count,&
      p%certified_basis%localization_iterations,p%rt_space%rank,p%rt_space%operation_count,&
      p%rt_space%scalar_count,p%rt_space%vector_count,p%rt_space%tensor_count,p%energy_window%mode,&
      p%energy_window%construction_rank,p%energy_window%solved_rank,p%energy_window%occupied_rank,&
      p%energy_window%requested_rank,p%energy_window%certified_rank,p%energy_window%extension_states,&
      p%energy_window%boundary_cluster_rank,p%energy_window%proof_status,p%symmetry_receipt%worst_operation]
    f=[p%catalog_fingerprint,p%state_fingerprint,p%metric_fingerprint,p%operator_structure_fingerprint,&
      p%operator_value_fingerprint,p%kinetic_fingerprint,p%nonlocal_fingerprint,p%local_fingerprint,&
      p%sipg_fingerprint,p%basis_fingerprint,p%face_fingerprint,p%dc_seed_fingerprint,&
      p%continuation_fingerprint,p%scope_fingerprint,p%analysis_fingerprint,p%selection_fingerprint,&
      p%pseudopotential_fingerprint,p%energy_fingerprint,p%position_convention_fingerprint,p%payload_fingerprint,&
      p%construction_catalog%ids_fingerprint,p%construction_catalog%generation_fingerprint,&
      p%construction_catalog%ordering_fingerprint,p%construction_catalog%ownership_fingerprint,&
      p%construction_catalog%provenance_fingerprint,p%construction_catalog%catalog_fingerprint,&
      p%certified_basis%c_cert_fingerprint,p%certified_basis%u_rt_fingerprint,&
      p%certified_basis%b_rt_fingerprint,p%certified_basis%initial_state_fingerprint,&
      p%certified_basis%transformation_fingerprint,p%certified_basis%operator_fingerprint,&
      p%certified_basis%fingerprint,p%electron_count%fingerprint,p%rt_space%metric_fingerprint,&
      p%rt_space%kinetic_fingerprint,p%rt_space%nonlocal_fingerprint,p%rt_space%local_fingerprint,&
      p%rt_space%sipg_fingerprint,p%rt_space%hamiltonian_fingerprint,p%rt_space%basis_fingerprint,&
      p%rt_space%density_fingerprint,p%rt_space%ownership_fingerprint,p%rt_space%scalar_fingerprint,&
      p%rt_space%vector_fingerprint,p%rt_space%tensor_fingerprint,p%rt_space%representation_fingerprint,&
      p%rt_space%fingerprint,p%energy_window%fingerprint,p%symmetry_receipt%fingerprint,&
      p%handoff_receipts%position_fingerprint,p%handoff_receipts%nonlocal_fingerprint,&
      p%handoff_receipts%face_fingerprint,p%handoff_receipts%pseudopotential_fingerprint,&
      p%handoff_receipts%transformation_fingerprint,p%handoff_receipts%fingerprint]
    r=[p%certified_basis%spread_before_total,p%certified_basis%spread_after_total,&
      p%certified_basis%spread_improvement,p%certified_basis%transform_unitarity_defect,&
      p%certified_basis%certified_metric_defect,p%certified_basis%rt_metric_defect,&
      p%certified_basis%embedding_defect,p%certified_basis%projector_invariance_defect,&
      p%certified_basis%target_symmetry_defect_before,p%certified_basis%target_symmetry_defect_after,&
      p%certified_basis%energy_symmetry_defect_before,p%certified_basis%energy_symmetry_defect_after,&
      p%certified_basis%symmetry_defect_invariance,p%certified_basis%scalar_covariance_defect,&
      p%certified_basis%vector_covariance_defect,p%certified_basis%tensor_covariance_defect,&
      p%electron_count%expected_count,p%electron_count%actual_count,p%electron_count%tolerance,&
      p%electron_count%defect,p%electron_count%omitted_tail,p%electron_count%chemical_potential,&
      p%energy_window%window_size,p%energy_window%e_homo,p%energy_window%requested_cutoff,&
      p%energy_window%certified_cutoff,p%energy_window%extension_energy,p%energy_window%proof_energy,&
      p%symmetry_receipt%occupied_subspace_defect,p%symmetry_receipt%occupied_projector_defect,&
      p%symmetry_receipt%target_subspace_defect,p%symmetry_receipt%target_energy_defect,&
      p%symmetry_receipt%density_defect,p%symmetry_receipt%scalar_covariance_defect,&
      p%symmetry_receipt%vector_covariance_defect,p%symmetry_receipt%tensor_covariance_defect,&
      p%symmetry_receipt%final_basis_defect,p%symmetry_receipt%worst_operation_defect,&
      p%symmetry_receipt%maximum_physical_defect]
  end subroutine pack_ground_state_header

  subroutine set_ground_state_header(p,l,h,f,r)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    logical,intent(in)::l(ground_state_logical_count)
    integer,intent(in)::h(ground_state_integer_count)
    integer(int64),intent(in)::f(ground_state_fingerprint_count)
    real(real64),intent(in)::r(ground_state_real_count)
    p%valid=l(1);p%final_refresh_complete=l(2);p%analysis_complete=l(3);p%identity_only=l(4)
    p%construction_catalog%valid=l(5);p%certified_basis%valid=l(6)
    p%certified_basis%localization_converged=l(7);p%certified_basis%localization_symmetry_constrained=l(8)
    p%electron_count%valid=l(9);p%rt_space%valid=l(10);p%energy_window%valid=l(11)
    p%energy_window%compatibility_dynamic_rank=l(12);p%energy_window%proof_state_present=l(13)
    p%symmetry_receipt%valid=l(14);p%handoff_receipts%valid=l(15)
    p%global_count=h(1);p%global_grid_count=h(2);p%noccupied=h(3);p%operation_count=h(4)
    p%nonidentity_operation_count=h(5);p%construction_catalog%global_count=h(6)
    p%certified_basis%construction_count=h(7);p%certified_basis%certified_count=h(8)
    p%certified_basis%occupied_count=h(9);p%certified_basis%localization_iterations=h(10)
    p%rt_space%rank=h(11);p%rt_space%operation_count=h(12);p%rt_space%scalar_count=h(13)
    p%rt_space%vector_count=h(14);p%rt_space%tensor_count=h(15);p%energy_window%mode=h(16)
    p%energy_window%construction_rank=h(17);p%energy_window%solved_rank=h(18)
    p%energy_window%occupied_rank=h(19);p%energy_window%requested_rank=h(20)
    p%energy_window%certified_rank=h(21);p%energy_window%extension_states=h(22)
    p%energy_window%boundary_cluster_rank=h(23);p%energy_window%proof_status=h(24)
    p%symmetry_receipt%worst_operation=h(25)
    p%catalog_fingerprint=f(1);p%state_fingerprint=f(2);p%metric_fingerprint=f(3)
    p%operator_structure_fingerprint=f(4);p%operator_value_fingerprint=f(5);p%kinetic_fingerprint=f(6)
    p%nonlocal_fingerprint=f(7);p%local_fingerprint=f(8);p%sipg_fingerprint=f(9);p%basis_fingerprint=f(10)
    p%face_fingerprint=f(11);p%dc_seed_fingerprint=f(12);p%continuation_fingerprint=f(13);p%scope_fingerprint=f(14)
    p%analysis_fingerprint=f(15);p%selection_fingerprint=f(16);p%pseudopotential_fingerprint=f(17)
    p%energy_fingerprint=f(18);p%position_convention_fingerprint=f(19);p%payload_fingerprint=f(20)
    p%construction_catalog%ids_fingerprint=f(21);p%construction_catalog%generation_fingerprint=f(22)
    p%construction_catalog%ordering_fingerprint=f(23);p%construction_catalog%ownership_fingerprint=f(24)
    p%construction_catalog%provenance_fingerprint=f(25);p%construction_catalog%catalog_fingerprint=f(26)
    p%certified_basis%c_cert_fingerprint=f(27);p%certified_basis%u_rt_fingerprint=f(28)
    p%certified_basis%b_rt_fingerprint=f(29);p%certified_basis%initial_state_fingerprint=f(30)
    p%certified_basis%transformation_fingerprint=f(31);p%certified_basis%operator_fingerprint=f(32)
    p%certified_basis%fingerprint=f(33);p%electron_count%fingerprint=f(34)
    p%rt_space%metric_fingerprint=f(35);p%rt_space%kinetic_fingerprint=f(36)
    p%rt_space%nonlocal_fingerprint=f(37);p%rt_space%local_fingerprint=f(38);p%rt_space%sipg_fingerprint=f(39)
    p%rt_space%hamiltonian_fingerprint=f(40);p%rt_space%basis_fingerprint=f(41)
    p%rt_space%density_fingerprint=f(42);p%rt_space%ownership_fingerprint=f(43)
    p%rt_space%scalar_fingerprint=f(44);p%rt_space%vector_fingerprint=f(45)
    p%rt_space%tensor_fingerprint=f(46);p%rt_space%representation_fingerprint=f(47);p%rt_space%fingerprint=f(48)
    p%energy_window%fingerprint=f(49);p%symmetry_receipt%fingerprint=f(50)
    p%handoff_receipts%position_fingerprint=f(51);p%handoff_receipts%nonlocal_fingerprint=f(52)
    p%handoff_receipts%face_fingerprint=f(53);p%handoff_receipts%pseudopotential_fingerprint=f(54)
    p%handoff_receipts%transformation_fingerprint=f(55);p%handoff_receipts%fingerprint=f(56)
    p%certified_basis%spread_before_total=r(1);p%certified_basis%spread_after_total=r(2)
    p%certified_basis%spread_improvement=r(3);p%certified_basis%transform_unitarity_defect=r(4)
    p%certified_basis%certified_metric_defect=r(5);p%certified_basis%rt_metric_defect=r(6)
    p%certified_basis%embedding_defect=r(7);p%certified_basis%projector_invariance_defect=r(8)
    p%certified_basis%target_symmetry_defect_before=r(9);p%certified_basis%target_symmetry_defect_after=r(10)
    p%certified_basis%energy_symmetry_defect_before=r(11);p%certified_basis%energy_symmetry_defect_after=r(12)
    p%certified_basis%symmetry_defect_invariance=r(13);p%certified_basis%scalar_covariance_defect=r(14)
    p%certified_basis%vector_covariance_defect=r(15);p%certified_basis%tensor_covariance_defect=r(16)
    p%electron_count%expected_count=r(17);p%electron_count%actual_count=r(18);p%electron_count%tolerance=r(19)
    p%electron_count%defect=r(20);p%electron_count%omitted_tail=r(21);p%electron_count%chemical_potential=r(22)
    p%energy_window%window_size=r(23);p%energy_window%e_homo=r(24);p%energy_window%requested_cutoff=r(25)
    p%energy_window%certified_cutoff=r(26);p%energy_window%extension_energy=r(27);p%energy_window%proof_energy=r(28)
    p%symmetry_receipt%occupied_subspace_defect=r(29);p%symmetry_receipt%occupied_projector_defect=r(30)
    p%symmetry_receipt%target_subspace_defect=r(31);p%symmetry_receipt%target_energy_defect=r(32)
    p%symmetry_receipt%density_defect=r(33);p%symmetry_receipt%scalar_covariance_defect=r(34)
    p%symmetry_receipt%vector_covariance_defect=r(35);p%symmetry_receipt%tensor_covariance_defect=r(36)
    p%symmetry_receipt%final_basis_defect=r(37);p%symmetry_receipt%worst_operation_defect=r(38)
    p%symmetry_receipt%maximum_physical_defect=r(39)
  end subroutine set_ground_state_header

#ifdef USE_MPI
  subroutine broadcast_coalesced_ground_state_common(comm,rank,p,have_common,ierr)
    integer,intent(in)::comm,rank
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    logical,intent(inout)::have_common
    integer,intent(out)::ierr
    logical::l(ground_state_logical_count)
    integer::h(ground_state_integer_count),layout(5)
    integer(int64)::f(ground_state_fingerprint_count)
    real(real64)::r(ground_state_real_count)
    if(rank==0)then
      call pack_ground_state_header(p,l,h,f,r)
      layout=[size(p%face_metadata,1),size(p%face_normals,1),size(p%face_values,1),&
        size(p%interface_observables,1),size(p%nonlocal_values,1)]
    endif
    call MPI_Bcast(l,ground_state_logical_count,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(h,ground_state_integer_count,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(f,ground_state_fingerprint_count,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(r,ground_state_real_count,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(layout,size(layout),MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call set_ground_state_header(p,l,h,f,r)
    call bcast_i1(p%scope_selectors);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%xc_types);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%requested_ids);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%effective_ids);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%added_ids);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%closure_parent);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%closure_reason);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%closure_action);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%occupations);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%eigenvalues);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%continuation_receipt);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%pseudopotential_receipt);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%energy_receipt);if(ierr/=MPI_SUCCESS)return
    call bcast_z3(p%symmetry_representation);if(ierr/=MPI_SUCCESS)return
    call bcast_i64_1(p%construction_catalog%ids);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%construction_catalog%generations);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%construction_catalog%ordering);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%construction_catalog%ownership);if(ierr/=MPI_SUCCESS)return
    call bcast_z2(p%certified_basis%initial_occupied_amplitudes);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%certified_basis%certified_eigenvalues);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%certified_basis%occupations);if(ierr/=MPI_SUCCESS)return
    call bcast_r2(p%certified_basis%centers);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%certified_basis%spreads_before);if(ierr/=MPI_SUCCESS)return
    call bcast_r1(p%certified_basis%spreads_after);if(ierr/=MPI_SUCCESS)return
    call bcast_i1(p%rt_space%row_owner_keys);if(ierr/=MPI_SUCCESS)return
    call bcast_z3(p%rt_space%representation);if(ierr/=MPI_SUCCESS)return
    call bcast_r3(p%rt_space%cartesian_rotations);if(ierr/=MPI_SUCCESS)return
    if(.not.have_common)call initialize_empty_distributed(p,layout)
    have_common=.true.
  contains
    subroutine bcast_i64_1(a)
      integer(int64),allocatable,intent(inout)::a(:)
      integer::count
      if(rank==0)count=size(a)
      call MPI_Bcast(count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(count));endif
      call MPI_Bcast(a,count,MPI_INTEGER8,0,comm,ierr)
    end subroutine bcast_i64_1
    subroutine bcast_i1(a)
      integer,allocatable,intent(inout)::a(:)
      integer::count
      if(rank==0)count=size(a)
      call MPI_Bcast(count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(count));endif
      call MPI_Bcast(a,count,MPI_INTEGER,0,comm,ierr)
    end subroutine bcast_i1
    subroutine bcast_r1(a)
      real(real64),allocatable,intent(inout)::a(:)
      integer::count
      if(rank==0)count=size(a)
      call MPI_Bcast(count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(count));endif
      call MPI_Bcast(a,count,MPI_DOUBLE_PRECISION,0,comm,ierr)
    end subroutine bcast_r1
    subroutine bcast_r2(a)
      real(real64),allocatable,intent(inout)::a(:,:)
      integer::dims(2)
      if(rank==0)dims=shape(a)
      call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(dims(1),dims(2)));endif
      call MPI_Bcast(a,product(dims),MPI_DOUBLE_PRECISION,0,comm,ierr)
    end subroutine bcast_r2
    subroutine bcast_r3(a)
      real(real64),allocatable,intent(inout)::a(:,:,:)
      integer::dims(3)
      if(rank==0)dims=shape(a)
      call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(dims(1),dims(2),dims(3)));endif
      call MPI_Bcast(a,product(dims),MPI_DOUBLE_PRECISION,0,comm,ierr)
    end subroutine bcast_r3
    subroutine bcast_z2(a)
      complex(real64),allocatable,intent(inout)::a(:,:)
      integer::dims(2)
      if(rank==0)dims=shape(a)
      call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(dims(1),dims(2)));endif
      call MPI_Bcast(a,product(dims),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    end subroutine bcast_z2
    subroutine bcast_z3(a)
      complex(real64),allocatable,intent(inout)::a(:,:,:)
      integer::dims(3)
      if(rank==0)dims=shape(a)
      call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(rank/=0)then;if(allocated(a))deallocate(a);allocate(a(dims(1),dims(2),dims(3)));endif
      call MPI_Bcast(a,product(dims),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    end subroutine bcast_z3
    subroutine initialize_empty_distributed(q,distributed_layout)
      type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::q
      integer,intent(in)::distributed_layout(5)
      integer::n,r_rt
      n=q%global_count;r_rt=q%certified_basis%certified_count
      allocate(q%row_ids(0),q%grid_ids(0),q%face_ids(0),q%face_point_ids(0),q%nonlocal_ids(0),&
        q%partition_ids(0),q%metric_row_offsets(1),q%metric_column_ids(0),q%operator_row_offsets(1),&
        q%operator_column_ids(0),q%face_offsets(1),q%face_weight_offsets(1),q%face_basis_offsets(1),&
        q%face_value_offsets(1),q%face_observable_offsets(1),q%face_basis_ids(0),&
        q%nonlocal_owner(0),q%grid_weights(0),q%face_weights(0),q%density(0))
      q%metric_row_offsets=1;q%operator_row_offsets=1;q%face_offsets=1;q%face_weight_offsets=1
      q%face_basis_offsets=1;q%face_value_offsets=1;q%face_observable_offsets=1
      allocate(q%face_metadata(distributed_layout(1),0),q%face_normals(distributed_layout(2),0),&
        q%face_values(distributed_layout(3),0),q%interface_observables(distributed_layout(4),0),&
        q%nonlocal_values(distributed_layout(5),0),q%metric_rows(0,n),q%kinetic_rows(0,n),&
        q%nonlocal_rows(0,n),q%local_rows(0,n),q%sipg_rows(0,n),q%hamiltonian_rows(0,n),&
        q%basis_values(n,0),q%coefficients(0,q%noccupied),q%position_rows(3,0,n))
      allocate(q%certified_basis%construction_row_ids(0),q%certified_basis%transformation_row_ids(0),&
        q%certified_basis%c_cert(0,r_rt),q%certified_basis%u_rt(0,r_rt),q%certified_basis%b_rt(0,r_rt))
      allocate(q%rt_space%row_ids(0),q%rt_space%grid_owner_keys(0),q%rt_space%metric_rows(0,r_rt),&
        q%rt_space%kinetic_rows(0,r_rt),q%rt_space%nonlocal_rows(0,r_rt),q%rt_space%local_rows(0,r_rt),&
        q%rt_space%sipg_rows(0,r_rt),q%rt_space%hamiltonian_rows(0,r_rt),&
        q%rt_space%scalar_operator_rows(0,r_rt,q%rt_space%scalar_count),&
        q%rt_space%vector_operator_rows(0,r_rt,3,q%rt_space%vector_count),&
        q%rt_space%tensor_operator_rows(0,r_rt,3,3,q%rt_space%tensor_count),&
        q%rt_space%basis_values(r_rt,0),q%rt_space%density(0))
    end subroutine initialize_empty_distributed
  end subroutine broadcast_coalesced_ground_state_common

  subroutine redistribute_coalesced_ground_state(comm,rank,nproc,p,ierr)
    integer,intent(in)::comm,rank,nproc
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer,intent(out)::ierr
    type(s_rt_dg_hybrid_ground_state_payload)::target
    integer::n,g,r_rt,nocc,id,i,source_index,owner,destination,local_index,first,last,&
      nlocal,metric_nnz,operator_nnz,layout(5),record_owner_buffer,transfer_status(MPI_STATUS_SIZE)
    integer,allocatable::local_counts(:),global_counts(:),local_owners(:),global_owners(:),&
      local_metric_degrees(:),global_metric_degrees(:),local_operator_degrees(:),global_operator_degrees(:),&
      metric_columns_buffer(:),operator_columns_buffer(:),local_positions(:)
    complex(real64),allocatable::row_components(:,:),coefficient_buffer(:),position_buffer(:,:),&
      c_buffer(:),b_buffer(:),grid_basis_buffer(:),rt_grid_basis_buffer(:),u_buffer(:),&
      rt_components(:,:),rt_scalar_buffer(:,:),rt_vector_buffer(:,:,:),rt_tensor_buffer(:,:,:,:)
    real(real64)::grid_real_buffer(3)
    integer::grid_integer_buffer(2)
    n=p%global_count;g=p%global_grid_count;r_rt=p%certified_basis%certified_count;nocc=p%noccupied
    ierr=MPI_SUCCESS
    call copy_ground_state_common(p,target)

    allocate(local_counts(n),global_counts(n),local_owners(n),global_owners(n),&
      local_metric_degrees(n),global_metric_degrees(n),local_operator_degrees(n),global_operator_degrees(n))
    allocate(local_positions(n));local_positions=0
    local_counts=0;local_owners=0;local_metric_degrees=0;local_operator_degrees=0
    do i=1,size(p%row_ids)
      id=int(p%row_ids(i));local_counts(id)=local_counts(id)+1;local_owners(id)=rank+1
      local_positions(id)=i
      local_metric_degrees(id)=p%metric_row_offsets(i+1)-p%metric_row_offsets(i)
      local_operator_degrees(id)=p%operator_row_offsets(i+1)-p%operator_row_offsets(i)
    enddo
    call MPI_Allreduce(local_counts,global_counts,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_owners,global_owners,n,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_metric_degrees,global_metric_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_operator_degrees,global_operator_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(any(global_counts/=1))then;ierr=1;return;endif
    nlocal=count([(mod(id-1,nproc)==rank,id=1,n)]);metric_nnz=0;operator_nnz=0
    do id=1,n
      if(mod(id-1,nproc)/=rank)cycle
      metric_nnz=metric_nnz+global_metric_degrees(id)
      operator_nnz=operator_nnz+global_operator_degrees(id)
    enddo
    deallocate(target%row_ids,target%metric_rows,target%kinetic_rows,target%nonlocal_rows,target%local_rows,&
      target%sipg_rows,target%hamiltonian_rows,target%coefficients,target%position_rows,&
      target%metric_row_offsets,target%metric_column_ids,target%operator_row_offsets,target%operator_column_ids,&
      target%certified_basis%construction_row_ids,target%certified_basis%c_cert,target%certified_basis%b_rt)
    allocate(target%row_ids(nlocal),target%metric_rows(nlocal,n),target%kinetic_rows(nlocal,n),&
      target%nonlocal_rows(nlocal,n),target%local_rows(nlocal,n),target%sipg_rows(nlocal,n),&
      target%hamiltonian_rows(nlocal,n),target%coefficients(nlocal,nocc),target%position_rows(3,nlocal,n),&
      target%metric_row_offsets(nlocal+1),target%metric_column_ids(metric_nnz),&
      target%operator_row_offsets(nlocal+1),target%operator_column_ids(operator_nnz),&
      target%certified_basis%construction_row_ids(nlocal),target%certified_basis%c_cert(nlocal,r_rt),&
      target%certified_basis%b_rt(nlocal,r_rt))
    target%metric_row_offsets(1)=1;target%operator_row_offsets(1)=1;local_index=0
    do id=1,n
      if(mod(id-1,nproc)/=rank)cycle
      local_index=local_index+1;target%row_ids(local_index)=int(id,int64)
      target%certified_basis%construction_row_ids(local_index)=int(id,int64)
      target%metric_row_offsets(local_index+1)=target%metric_row_offsets(local_index)+global_metric_degrees(id)
      target%operator_row_offsets(local_index+1)=target%operator_row_offsets(local_index)+global_operator_degrees(id)
    enddo
    allocate(metric_columns_buffer(max(1,maxval(global_metric_degrees))),&
      operator_columns_buffer(max(1,maxval(global_operator_degrees))),row_components(6,n),&
      coefficient_buffer(nocc),position_buffer(3,n),c_buffer(r_rt),b_buffer(r_rt))
    do id=1,n
      owner=global_owners(id)-1;source_index=0
      destination=mod(id-1,nproc)
      if(rank==owner)source_index=local_positions(id)
      if(rank==owner)then
        first=p%metric_row_offsets(source_index);last=p%metric_row_offsets(source_index+1)-1
        if(last>=first)metric_columns_buffer(:last-first+1)=p%metric_column_ids(first:last)
        first=p%operator_row_offsets(source_index);last=p%operator_row_offsets(source_index+1)-1
        if(last>=first)operator_columns_buffer(:last-first+1)=p%operator_column_ids(first:last)
        row_components(1,:)=p%metric_rows(source_index,:);row_components(2,:)=p%kinetic_rows(source_index,:)
        row_components(3,:)=p%nonlocal_rows(source_index,:);row_components(4,:)=p%local_rows(source_index,:)
        row_components(5,:)=p%sipg_rows(source_index,:);row_components(6,:)=p%hamiltonian_rows(source_index,:)
        coefficient_buffer=p%coefficients(source_index,:);position_buffer=p%position_rows(:,source_index,:)
        c_buffer=p%certified_basis%c_cert(source_index,:);b_buffer=p%certified_basis%b_rt(source_index,:)
      endif
      if(owner/=destination)then
        if(rank==owner)then
          call MPI_Send(metric_columns_buffer,global_metric_degrees(id),MPI_INTEGER,destination,29101,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(operator_columns_buffer,global_operator_degrees(id),MPI_INTEGER,destination,29102,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(row_components,6*n,MPI_DOUBLE_COMPLEX,destination,29103,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(coefficient_buffer,nocc,MPI_DOUBLE_COMPLEX,destination,29104,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(position_buffer,3*n,MPI_DOUBLE_COMPLEX,destination,29105,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(c_buffer,r_rt,MPI_DOUBLE_COMPLEX,destination,29106,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(b_buffer,r_rt,MPI_DOUBLE_COMPLEX,destination,29107,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
        else if(rank==destination)then
          call MPI_Recv(metric_columns_buffer,global_metric_degrees(id),MPI_INTEGER,owner,29101,comm,&
            transfer_status,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(operator_columns_buffer,global_operator_degrees(id),MPI_INTEGER,owner,29102,comm,&
            transfer_status,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(row_components,6*n,MPI_DOUBLE_COMPLEX,owner,29103,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(coefficient_buffer,nocc,MPI_DOUBLE_COMPLEX,owner,29104,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(position_buffer,3*n,MPI_DOUBLE_COMPLEX,owner,29105,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(c_buffer,r_rt,MPI_DOUBLE_COMPLEX,owner,29106,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(b_buffer,r_rt,MPI_DOUBLE_COMPLEX,owner,29107,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
        endif
      endif
      if(rank==destination)then
        local_index=(id-1)/nproc+1
        target%metric_rows(local_index,:)=row_components(1,:)
        target%kinetic_rows(local_index,:)=row_components(2,:)
        target%nonlocal_rows(local_index,:)=row_components(3,:)
        target%local_rows(local_index,:)=row_components(4,:)
        target%sipg_rows(local_index,:)=row_components(5,:)
        target%hamiltonian_rows(local_index,:)=row_components(6,:)
        target%coefficients(local_index,:)=coefficient_buffer;target%position_rows(:,local_index,:)=position_buffer
        target%certified_basis%c_cert(local_index,:)=c_buffer;target%certified_basis%b_rt(local_index,:)=b_buffer
        first=target%metric_row_offsets(local_index);last=target%metric_row_offsets(local_index+1)-1
        if(last>=first)target%metric_column_ids(first:last)=metric_columns_buffer(:last-first+1)
        first=target%operator_row_offsets(local_index);last=target%operator_row_offsets(local_index+1)-1
        if(last>=first)target%operator_column_ids(first:last)=operator_columns_buffer(:last-first+1)
      endif
    enddo
    deallocate(local_counts,global_counts,local_owners,global_owners,local_metric_degrees,global_metric_degrees,&
      local_operator_degrees,global_operator_degrees,metric_columns_buffer,operator_columns_buffer,&
      row_components,coefficient_buffer,position_buffer,c_buffer,b_buffer,local_positions)

    allocate(local_counts(g),global_counts(g),local_owners(g),global_owners(g),local_positions(g))
    local_counts=0;local_owners=0;local_positions=0
    do i=1,size(p%grid_ids)
      id=int(p%grid_ids(i));local_counts(id)=local_counts(id)+1;local_owners(id)=rank+1
      local_positions(id)=i
    enddo
    call MPI_Allreduce(local_counts,global_counts,g,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_owners,global_owners,g,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(global_counts/=1))then;ierr=1;return;endif
    nlocal=count([(mod(id-1,nproc)==rank,id=1,g)])
    deallocate(target%grid_ids,target%partition_ids,target%grid_weights,target%density,target%basis_values,&
      target%rt_space%grid_owner_keys,target%rt_space%basis_values,target%rt_space%density)
    allocate(target%grid_ids(nlocal),target%partition_ids(nlocal),target%grid_weights(nlocal),target%density(nlocal),&
      target%basis_values(n,nlocal),target%rt_space%grid_owner_keys(nlocal),&
      target%rt_space%basis_values(r_rt,nlocal),target%rt_space%density(nlocal),&
      grid_basis_buffer(n),rt_grid_basis_buffer(r_rt))
    do id=1,g
      owner=global_owners(id)-1;source_index=0
      destination=mod(id-1,nproc)
      if(rank==owner)source_index=local_positions(id)
      if(rank==owner)then
        grid_integer_buffer=[p%partition_ids(source_index),p%rt_space%grid_owner_keys(source_index)]
        grid_real_buffer=[p%grid_weights(source_index),p%density(source_index),p%rt_space%density(source_index)]
        grid_basis_buffer=p%basis_values(:,source_index);rt_grid_basis_buffer=p%rt_space%basis_values(:,source_index)
      endif
      if(owner/=destination)then
        if(rank==owner)then
          call MPI_Send(grid_integer_buffer,2,MPI_INTEGER,destination,29111,comm,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Send(grid_real_buffer,3,MPI_DOUBLE_PRECISION,destination,29112,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(grid_basis_buffer,n,MPI_DOUBLE_COMPLEX,destination,29113,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(rt_grid_basis_buffer,r_rt,MPI_DOUBLE_COMPLEX,destination,29114,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
        else if(rank==destination)then
          call MPI_Recv(grid_integer_buffer,2,MPI_INTEGER,owner,29111,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(grid_real_buffer,3,MPI_DOUBLE_PRECISION,owner,29112,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(grid_basis_buffer,n,MPI_DOUBLE_COMPLEX,owner,29113,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(rt_grid_basis_buffer,r_rt,MPI_DOUBLE_COMPLEX,owner,29114,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
        endif
      endif
      if(rank==destination)then
        local_index=(id-1)/nproc+1;target%grid_ids(local_index)=int(id,int64)
        target%partition_ids(local_index)=grid_integer_buffer(1)
        target%rt_space%grid_owner_keys(local_index)=grid_integer_buffer(2)
        target%grid_weights(local_index)=grid_real_buffer(1);target%density(local_index)=grid_real_buffer(2)
        target%rt_space%density(local_index)=grid_real_buffer(3)
        target%basis_values(:,local_index)=grid_basis_buffer
        target%rt_space%basis_values(:,local_index)=rt_grid_basis_buffer
      endif
    enddo
    deallocate(local_counts,global_counts,local_owners,global_owners,grid_basis_buffer,rt_grid_basis_buffer,&
      local_positions)

    allocate(local_counts(r_rt),global_counts(r_rt),local_owners(r_rt),global_owners(r_rt),u_buffer(r_rt),&
      local_positions(r_rt))
    local_counts=0;local_owners=0;local_positions=0
    do i=1,size(p%certified_basis%transformation_row_ids)
      id=int(p%certified_basis%transformation_row_ids(i));local_counts(id)=local_counts(id)+1
      local_owners(id)=rank+1;local_positions(id)=i
    enddo
    call MPI_Allreduce(local_counts,global_counts,r_rt,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_owners,global_owners,r_rt,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(global_counts/=1))then;ierr=1;return;endif
    nlocal=count([(mod(id-1,nproc)==rank,id=1,r_rt)])
    deallocate(target%certified_basis%transformation_row_ids,target%certified_basis%u_rt)
    allocate(target%certified_basis%transformation_row_ids(nlocal),target%certified_basis%u_rt(nlocal,r_rt))
    do id=1,r_rt
      owner=global_owners(id)-1;source_index=0
      destination=mod(id-1,nproc)
      if(rank==owner)then
        source_index=local_positions(id)
        u_buffer=p%certified_basis%u_rt(source_index,:)
      endif
      if(owner/=destination)then
        if(rank==owner)call MPI_Send(u_buffer,r_rt,MPI_DOUBLE_COMPLEX,destination,29121,comm,ierr)
        if(rank==destination)call MPI_Recv(u_buffer,r_rt,MPI_DOUBLE_COMPLEX,owner,29121,comm,transfer_status,ierr)
        if((rank==owner.or.rank==destination).and.ierr/=MPI_SUCCESS)return
      endif
      if(rank==destination)then
        local_index=(id-1)/nproc+1;target%certified_basis%transformation_row_ids(local_index)=int(id,int64)
        target%certified_basis%u_rt(local_index,:)=u_buffer
      endif
    enddo

    local_counts=0;local_owners=0;local_positions=0
    do i=1,size(p%rt_space%row_ids)
      id=int(p%rt_space%row_ids(i));local_counts(id)=local_counts(id)+1;local_owners(id)=rank+1
      local_positions(id)=i
    enddo
    call MPI_Allreduce(local_counts,global_counts,r_rt,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_owners,global_owners,r_rt,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(global_counts/=1))then;ierr=1;return;endif
    deallocate(target%rt_space%row_ids,target%rt_space%metric_rows,target%rt_space%kinetic_rows,&
      target%rt_space%nonlocal_rows,target%rt_space%local_rows,target%rt_space%sipg_rows,&
      target%rt_space%hamiltonian_rows,target%rt_space%scalar_operator_rows,&
      target%rt_space%vector_operator_rows,target%rt_space%tensor_operator_rows)
    allocate(target%rt_space%row_ids(nlocal),target%rt_space%metric_rows(nlocal,r_rt),&
      target%rt_space%kinetic_rows(nlocal,r_rt),target%rt_space%nonlocal_rows(nlocal,r_rt),&
      target%rt_space%local_rows(nlocal,r_rt),target%rt_space%sipg_rows(nlocal,r_rt),&
      target%rt_space%hamiltonian_rows(nlocal,r_rt),&
      target%rt_space%scalar_operator_rows(nlocal,r_rt,target%rt_space%scalar_count),&
      target%rt_space%vector_operator_rows(nlocal,r_rt,3,target%rt_space%vector_count),&
      target%rt_space%tensor_operator_rows(nlocal,r_rt,3,3,target%rt_space%tensor_count),&
      rt_components(6,r_rt),rt_scalar_buffer(r_rt,target%rt_space%scalar_count),&
      rt_vector_buffer(r_rt,3,target%rt_space%vector_count),&
      rt_tensor_buffer(r_rt,3,3,target%rt_space%tensor_count))
    do id=1,r_rt
      owner=global_owners(id)-1;source_index=0
      destination=mod(id-1,nproc)
      if(rank==owner)then
        source_index=local_positions(id)
        rt_components(1,:)=p%rt_space%metric_rows(source_index,:)
        rt_components(2,:)=p%rt_space%kinetic_rows(source_index,:)
        rt_components(3,:)=p%rt_space%nonlocal_rows(source_index,:)
        rt_components(4,:)=p%rt_space%local_rows(source_index,:)
        rt_components(5,:)=p%rt_space%sipg_rows(source_index,:)
        rt_components(6,:)=p%rt_space%hamiltonian_rows(source_index,:)
        rt_scalar_buffer=p%rt_space%scalar_operator_rows(source_index,:,:)
        rt_vector_buffer=p%rt_space%vector_operator_rows(source_index,:,:,:)
        rt_tensor_buffer=p%rt_space%tensor_operator_rows(source_index,:,:,:,:)
      endif
      if(owner/=destination)then
        if(rank==owner)then
          call MPI_Send(rt_components,6*r_rt,MPI_DOUBLE_COMPLEX,destination,29131,comm,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Send(rt_scalar_buffer,r_rt*target%rt_space%scalar_count,MPI_DOUBLE_COMPLEX,destination,29132,&
            comm,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Send(rt_vector_buffer,3*r_rt*target%rt_space%vector_count,MPI_DOUBLE_COMPLEX,destination,29133,&
            comm,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Send(rt_tensor_buffer,9*r_rt*target%rt_space%tensor_count,MPI_DOUBLE_COMPLEX,destination,29134,&
            comm,ierr);if(ierr/=MPI_SUCCESS)return
        else if(rank==destination)then
          call MPI_Recv(rt_components,6*r_rt,MPI_DOUBLE_COMPLEX,owner,29131,comm,transfer_status,ierr)
          if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(rt_scalar_buffer,r_rt*target%rt_space%scalar_count,MPI_DOUBLE_COMPLEX,owner,29132,comm,&
            transfer_status,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(rt_vector_buffer,3*r_rt*target%rt_space%vector_count,MPI_DOUBLE_COMPLEX,owner,29133,comm,&
            transfer_status,ierr);if(ierr/=MPI_SUCCESS)return
          call MPI_Recv(rt_tensor_buffer,9*r_rt*target%rt_space%tensor_count,MPI_DOUBLE_COMPLEX,owner,29134,comm,&
            transfer_status,ierr);if(ierr/=MPI_SUCCESS)return
        endif
      endif
      if(rank==destination)then
        local_index=(id-1)/nproc+1;target%rt_space%row_ids(local_index)=int(id,int64)
        target%rt_space%metric_rows(local_index,:)=rt_components(1,:)
        target%rt_space%kinetic_rows(local_index,:)=rt_components(2,:)
        target%rt_space%nonlocal_rows(local_index,:)=rt_components(3,:)
        target%rt_space%local_rows(local_index,:)=rt_components(4,:)
        target%rt_space%sipg_rows(local_index,:)=rt_components(5,:)
        target%rt_space%hamiltonian_rows(local_index,:)=rt_components(6,:)
        target%rt_space%scalar_operator_rows(local_index,:,:)=rt_scalar_buffer
        target%rt_space%vector_operator_rows(local_index,:,:,:)=rt_vector_buffer
        target%rt_space%tensor_operator_rows(local_index,:,:,:,:)=rt_tensor_buffer
      endif
    enddo
    deallocate(local_counts,global_counts,local_owners,global_owners,u_buffer,rt_components,local_positions,&
      rt_scalar_buffer,rt_vector_buffer,rt_tensor_buffer)

    if(rank==0)layout=[size(p%face_metadata,1),size(p%face_normals,1),size(p%face_values,1),&
      size(p%interface_observables,1),size(p%nonlocal_values,1)]
    call MPI_Bcast(layout,size(layout),MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    record_owner_buffer=0
    if(size(p%face_ids)>0)then
      if(any([size(p%face_metadata,1),size(p%face_normals,1),size(p%face_values,1),&
        size(p%interface_observables,1)]/=layout(:4)))record_owner_buffer=1
    endif
    if(size(p%nonlocal_ids)>0)then
      if(size(p%nonlocal_values,1)/=layout(5))record_owner_buffer=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,record_owner_buffer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(record_owner_buffer/=0)then;ierr=1;return;endif

    call redistribute_face_payload(comm,rank,nproc,p,target,layout,ierr)
    if(ierr/=MPI_SUCCESS)return
    call redistribute_nonlocal_payload(comm,rank,nproc,p,target,layout(5),ierr)
    if(ierr/=MPI_SUCCESS)return
    call replace_distributed_payload(p,target)
  contains
    subroutine redistribute_face_payload(comm,rank,nproc,source,target,face_layout,ierr)
      integer,intent(in)::comm,rank,nproc,face_layout(5)
      type(s_rt_dg_hybrid_ground_state_payload),intent(in)::source
      type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::target
      integer,intent(out)::ierr
      integer::i,k,category,destination,packed,source_index,original_index,nlocal,nreceived,&
        first,last,entry_count,position,local_bad,global_bad,transfer_total
      integer,allocatable::send_record_counts(:),receive_record_counts(:),send_record_displacements(:),&
        receive_record_displacements(:),next_record(:),send_positions(:),record_order(:),&
        scaled_send_counts(:),scaled_receive_counts(:),scaled_send_displacements(:),&
        scaled_receive_displacements(:),send_data_counts(:,:),receive_data_counts(:,:),&
        send_data_displacements(:,:),receive_data_displacements(:,:),next_data(:,:),&
        send_degrees(:,:),receive_degrees(:,:),receive_data_offsets(:,:),&
        send_metadata(:,:),receive_metadata(:,:),send_basis_ids(:),receive_basis_ids(:)
      integer(int64),allocatable::send_face_ids(:),receive_face_ids(:),send_point_ids(:),receive_point_ids(:)
      integer(int64),allocatable::send_data_counts64(:,:)
      integer(int64)::degree_totals(5),count_limit
      real(real64),allocatable::send_normals(:,:),receive_normals(:,:),send_weights(:),receive_weights(:)
      complex(real64),allocatable::send_values(:),receive_values(:),send_observables(:),receive_observables(:)
      ierr=MPI_SUCCESS;nlocal=size(source%face_ids);local_bad=0
      if(nlocal>0)then
        if(any(source%face_ids<=0_int64))local_bad=1
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then;ierr=1;return;endif
      allocate(send_record_counts(nproc),receive_record_counts(nproc),send_record_displacements(nproc),&
        receive_record_displacements(nproc),next_record(nproc))
      send_record_counts=0
      do i=1,nlocal
        destination=int(mod(source%face_ids(i)-1_int64,int(nproc,int64)))+1
        send_record_counts(destination)=send_record_counts(destination)+1
      enddo
      call MPI_Alltoall(send_record_counts,1,MPI_INTEGER,receive_record_counts,1,MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      allocate(scaled_send_counts(nproc),scaled_receive_counts(nproc),scaled_send_displacements(nproc),&
        scaled_receive_displacements(nproc))
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,1,scaled_send_counts,&
        scaled_receive_counts,send_record_displacements,receive_record_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      nreceived=sum(scaled_receive_counts);next_record=send_record_displacements+1
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,5,scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,face_layout(1),scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,face_layout(2),scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      allocate(send_positions(nlocal),send_face_ids(nlocal),send_degrees(5,nlocal),&
        send_metadata(face_layout(1),nlocal),send_normals(face_layout(2),nlocal),&
        receive_face_ids(nreceived),receive_degrees(5,nreceived),receive_metadata(face_layout(1),nreceived),&
        receive_normals(face_layout(2),nreceived),record_order(nreceived))
      do source_index=1,nlocal
        destination=int(mod(source%face_ids(source_index)-1_int64,int(nproc,int64)))+1
        packed=next_record(destination);next_record(destination)=packed+1
        send_positions(packed)=source_index;send_face_ids(packed)=source%face_ids(source_index)
        send_degrees(:,packed)=[source%face_offsets(source_index+1)-source%face_offsets(source_index),&
          source%face_weight_offsets(source_index+1)-source%face_weight_offsets(source_index),&
          source%face_basis_offsets(source_index+1)-source%face_basis_offsets(source_index),&
          source%face_value_offsets(source_index+1)-source%face_value_offsets(source_index),&
          source%face_observable_offsets(source_index+1)-source%face_observable_offsets(source_index)]
        send_metadata(:,packed)=source%face_metadata(:,source_index)
        send_normals(:,packed)=source%face_normals(:,source_index)
      enddo
      call MPI_Alltoallv(send_face_ids,send_record_counts,send_record_displacements,MPI_INTEGER8,&
        receive_face_ids,receive_record_counts,receive_record_displacements,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,5,scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_degrees,scaled_send_counts,scaled_send_displacements,MPI_INTEGER,&
        receive_degrees,scaled_receive_counts,scaled_receive_displacements,MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,face_layout(1),scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_metadata,scaled_send_counts,scaled_send_displacements,MPI_INTEGER,&
        receive_metadata,scaled_receive_counts,scaled_receive_displacements,MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_record_counts,receive_record_counts,face_layout(2),scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_normals,scaled_send_counts,scaled_send_displacements,MPI_DOUBLE_PRECISION,&
        receive_normals,scaled_receive_counts,scaled_receive_displacements,MPI_DOUBLE_PRECISION,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      do i=1,nreceived;record_order(i)=i;enddo
      call sort_id_positions(receive_face_ids,record_order)
      local_bad=0;degree_totals=0_int64;count_limit=int(huge(0),int64)
      do i=1,nreceived
        if(any(receive_degrees(:,i)<0))local_bad=1
        degree_totals=degree_totals+int(receive_degrees(:,i),int64)
        if(receive_face_ids(i)<=0_int64)then
          local_bad=1
        else if(mod(receive_face_ids(i)-1_int64,int(nproc,int64))/=int(rank,int64))then
          local_bad=1
        endif
        if(i>1)then
          if(receive_face_ids(i)==receive_face_ids(i-1))local_bad=1
        endif
      enddo
      if(any(degree_totals>count_limit))local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then;ierr=1;return;endif

      allocate(send_data_counts(nproc,5),receive_data_counts(nproc,5),send_data_displacements(nproc,5),&
        receive_data_displacements(nproc,5),next_data(nproc,5),send_data_counts64(nproc,5))
      send_data_counts64=0_int64
      do packed=1,nlocal
        destination=int(mod(send_face_ids(packed)-1_int64,int(nproc,int64)))+1
        send_data_counts64(destination,1)=send_data_counts64(destination,1)+int(send_degrees(1,packed),int64)
        send_data_counts64(destination,2)=send_data_counts64(destination,2)+int(send_degrees(2,packed),int64)
        send_data_counts64(destination,3)=send_data_counts64(destination,3)+int(send_degrees(3,packed),int64)
        send_data_counts64(destination,4)=send_data_counts64(destination,4)+&
          int(face_layout(3),int64)*int(send_degrees(4,packed),int64)
        send_data_counts64(destination,5)=send_data_counts64(destination,5)+&
          int(face_layout(4),int64)*int(send_degrees(5,packed),int64)
      enddo
      local_bad=0;if(any(send_data_counts64<0_int64).or.any(send_data_counts64>count_limit))local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then;ierr=1;return;endif
      send_data_counts=int(send_data_counts64)
      do category=1,5
        call MPI_Alltoall(send_data_counts(:,category),1,MPI_INTEGER,receive_data_counts(:,category),1,&
          MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
        call scale_counts_checked(comm,send_data_counts(:,category),receive_data_counts(:,category),1,&
          scaled_send_counts,scaled_receive_counts,send_data_displacements(:,category),&
          receive_data_displacements(:,category),ierr)
        if(ierr/=MPI_SUCCESS)return
      enddo
      next_data=send_data_displacements+1
      allocate(send_point_ids(sum(send_data_counts(:,1))),receive_point_ids(sum(receive_data_counts(:,1))),&
        send_weights(sum(send_data_counts(:,2))),receive_weights(sum(receive_data_counts(:,2))),&
        send_basis_ids(sum(send_data_counts(:,3))),receive_basis_ids(sum(receive_data_counts(:,3))),&
        send_values(sum(send_data_counts(:,4))),receive_values(sum(receive_data_counts(:,4))),&
        send_observables(sum(send_data_counts(:,5))),receive_observables(sum(receive_data_counts(:,5))))
      do packed=1,nlocal
        source_index=send_positions(packed)
        destination=int(mod(send_face_ids(packed)-1_int64,int(nproc,int64)))+1
        first=source%face_offsets(source_index);last=source%face_offsets(source_index+1)-1
        entry_count=max(0,last-first+1);position=next_data(destination,1)
        if(entry_count>0)send_point_ids(position:position+entry_count-1)=source%face_point_ids(first:last)
        next_data(destination,1)=position+entry_count
        first=source%face_weight_offsets(source_index);last=source%face_weight_offsets(source_index+1)-1
        entry_count=max(0,last-first+1);position=next_data(destination,2)
        if(entry_count>0)send_weights(position:position+entry_count-1)=source%face_weights(first:last)
        next_data(destination,2)=position+entry_count
        first=source%face_basis_offsets(source_index);last=source%face_basis_offsets(source_index+1)-1
        entry_count=max(0,last-first+1);position=next_data(destination,3)
        if(entry_count>0)send_basis_ids(position:position+entry_count-1)=source%face_basis_ids(first:last)
        next_data(destination,3)=position+entry_count
        first=source%face_value_offsets(source_index);last=source%face_value_offsets(source_index+1)-1
        entry_count=face_layout(3)*max(0,last-first+1);position=next_data(destination,4)
        if(entry_count>0)send_values(position:position+entry_count-1)=&
          reshape(source%face_values(:,first:last),[entry_count])
        next_data(destination,4)=position+entry_count
        first=source%face_observable_offsets(source_index);last=source%face_observable_offsets(source_index+1)-1
        entry_count=face_layout(4)*max(0,last-first+1);position=next_data(destination,5)
        if(entry_count>0)send_observables(position:position+entry_count-1)=&
          reshape(source%interface_observables(:,first:last),[entry_count])
        next_data(destination,5)=position+entry_count
      enddo
      call MPI_Alltoallv(send_point_ids,send_data_counts(:,1),send_data_displacements(:,1),MPI_INTEGER8,&
        receive_point_ids,receive_data_counts(:,1),receive_data_displacements(:,1),MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_weights,send_data_counts(:,2),send_data_displacements(:,2),MPI_DOUBLE_PRECISION,&
        receive_weights,receive_data_counts(:,2),receive_data_displacements(:,2),MPI_DOUBLE_PRECISION,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_basis_ids,send_data_counts(:,3),send_data_displacements(:,3),MPI_INTEGER,&
        receive_basis_ids,receive_data_counts(:,3),receive_data_displacements(:,3),MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_values,send_data_counts(:,4),send_data_displacements(:,4),MPI_DOUBLE_COMPLEX,&
        receive_values,receive_data_counts(:,4),receive_data_displacements(:,4),MPI_DOUBLE_COMPLEX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_observables,send_data_counts(:,5),send_data_displacements(:,5),MPI_DOUBLE_COMPLEX,&
        receive_observables,receive_data_counts(:,5),receive_data_displacements(:,5),MPI_DOUBLE_COMPLEX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return

      allocate(receive_data_offsets(5,nreceived+1));receive_data_offsets(:,1)=1
      do i=1,nreceived
        receive_data_offsets(1,i+1)=receive_data_offsets(1,i)+receive_degrees(1,i)
        receive_data_offsets(2,i+1)=receive_data_offsets(2,i)+receive_degrees(2,i)
        receive_data_offsets(3,i+1)=receive_data_offsets(3,i)+receive_degrees(3,i)
        receive_data_offsets(4,i+1)=receive_data_offsets(4,i)+face_layout(3)*receive_degrees(4,i)
        receive_data_offsets(5,i+1)=receive_data_offsets(5,i)+face_layout(4)*receive_degrees(5,i)
      enddo
      deallocate(target%face_ids,target%face_point_ids,target%face_metadata,target%face_offsets,&
        target%face_weight_offsets,target%face_basis_offsets,target%face_value_offsets,&
        target%face_observable_offsets,target%face_basis_ids,target%face_normals,target%face_weights,&
        target%face_values,target%interface_observables)
      allocate(target%face_ids(nreceived),target%face_metadata(face_layout(1),nreceived),&
        target%face_normals(face_layout(2),nreceived),target%face_offsets(nreceived+1),&
        target%face_weight_offsets(nreceived+1),target%face_basis_offsets(nreceived+1),&
        target%face_value_offsets(nreceived+1),target%face_observable_offsets(nreceived+1))
      target%face_offsets(1)=1;target%face_weight_offsets(1)=1;target%face_basis_offsets(1)=1
      target%face_value_offsets(1)=1;target%face_observable_offsets(1)=1
      do k=1,nreceived
        original_index=record_order(k);target%face_ids(k)=receive_face_ids(k)
        target%face_offsets(k+1)=target%face_offsets(k)+receive_degrees(1,original_index)
        target%face_weight_offsets(k+1)=target%face_weight_offsets(k)+receive_degrees(2,original_index)
        target%face_basis_offsets(k+1)=target%face_basis_offsets(k)+receive_degrees(3,original_index)
        target%face_value_offsets(k+1)=target%face_value_offsets(k)+receive_degrees(4,original_index)
        target%face_observable_offsets(k+1)=&
          target%face_observable_offsets(k)+receive_degrees(5,original_index)
      enddo
      allocate(target%face_point_ids(target%face_offsets(nreceived+1)-1),&
        target%face_weights(target%face_weight_offsets(nreceived+1)-1),&
        target%face_basis_ids(target%face_basis_offsets(nreceived+1)-1),&
        target%face_values(face_layout(3),target%face_value_offsets(nreceived+1)-1),&
        target%interface_observables(face_layout(4),target%face_observable_offsets(nreceived+1)-1))
      do k=1,nreceived
        original_index=record_order(k);target%face_metadata(:,k)=receive_metadata(:,original_index)
        target%face_normals(:,k)=receive_normals(:,original_index)
        transfer_total=receive_degrees(1,original_index);first=target%face_offsets(k)
        position=receive_data_offsets(1,original_index)
        if(transfer_total>0)target%face_point_ids(first:first+transfer_total-1)=&
          receive_point_ids(position:position+transfer_total-1)
        transfer_total=receive_degrees(2,original_index);first=target%face_weight_offsets(k)
        position=receive_data_offsets(2,original_index)
        if(transfer_total>0)target%face_weights(first:first+transfer_total-1)=&
          receive_weights(position:position+transfer_total-1)
        transfer_total=receive_degrees(3,original_index);first=target%face_basis_offsets(k)
        position=receive_data_offsets(3,original_index)
        if(transfer_total>0)target%face_basis_ids(first:first+transfer_total-1)=&
          receive_basis_ids(position:position+transfer_total-1)
        transfer_total=receive_degrees(4,original_index);first=target%face_value_offsets(k)
        position=receive_data_offsets(4,original_index);entry_count=face_layout(3)*transfer_total
        if(entry_count>0)target%face_values(:,first:first+transfer_total-1)=&
          reshape(receive_values(position:position+entry_count-1),[face_layout(3),transfer_total])
        transfer_total=receive_degrees(5,original_index);first=target%face_observable_offsets(k)
        position=receive_data_offsets(5,original_index);entry_count=face_layout(4)*transfer_total
        if(entry_count>0)target%interface_observables(:,first:first+transfer_total-1)=&
          reshape(receive_observables(position:position+entry_count-1),[face_layout(4),transfer_total])
      enddo
    end subroutine redistribute_face_payload

    subroutine redistribute_nonlocal_payload(comm,rank,nproc,source,target,value_layout,ierr)
      integer,intent(in)::comm,rank,nproc,value_layout
      type(s_rt_dg_hybrid_ground_state_payload),intent(in)::source
      type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::target
      integer,intent(out)::ierr
      integer::i,k,destination,packed,source_index,nlocal,nreceived,local_bad,global_bad
      integer,allocatable::send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:),&
        next_record(:),send_positions(:),record_order(:),send_owner_keys(:),receive_owner_keys(:),&
        scaled_send_counts(:),scaled_receive_counts(:),scaled_send_displacements(:),scaled_receive_displacements(:)
      integer(int64),allocatable::send_ids(:),receive_ids(:)
      complex(real64),allocatable::send_values(:,:),receive_values(:,:)
      ierr=MPI_SUCCESS;nlocal=size(source%nonlocal_ids);local_bad=0
      if(nlocal>0)then
        if(any(source%nonlocal_ids<=0_int64))local_bad=1
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then;ierr=1;return;endif
      allocate(send_counts(nproc),receive_counts(nproc),send_displacements(nproc),receive_displacements(nproc),&
        next_record(nproc));send_counts=0
      do i=1,nlocal
        destination=int(mod(source%nonlocal_ids(i)-1_int64,int(nproc,int64)))+1
        send_counts(destination)=send_counts(destination)+1
      enddo
      call MPI_Alltoall(send_counts,1,MPI_INTEGER,receive_counts,1,MPI_INTEGER,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      allocate(scaled_send_counts(nproc),scaled_receive_counts(nproc),scaled_send_displacements(nproc),&
        scaled_receive_displacements(nproc))
      call scale_counts_checked(comm,send_counts,receive_counts,1,scaled_send_counts,scaled_receive_counts,&
        send_displacements,receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      nreceived=sum(scaled_receive_counts);next_record=send_displacements+1
      call scale_counts_checked(comm,send_counts,receive_counts,value_layout,scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      allocate(send_positions(nlocal),send_ids(nlocal),send_owner_keys(nlocal),send_values(value_layout,nlocal),&
        receive_ids(nreceived),receive_owner_keys(nreceived),receive_values(value_layout,nreceived),&
        record_order(nreceived))
      do source_index=1,nlocal
        destination=int(mod(source%nonlocal_ids(source_index)-1_int64,int(nproc,int64)))+1
        packed=next_record(destination);next_record(destination)=packed+1
        send_positions(packed)=source_index;send_ids(packed)=source%nonlocal_ids(source_index)
        send_owner_keys(packed)=source%nonlocal_owner(source_index)
        send_values(:,packed)=source%nonlocal_values(:,source_index)
      enddo
      call MPI_Alltoallv(send_ids,send_counts,send_displacements,MPI_INTEGER8,receive_ids,receive_counts,&
        receive_displacements,MPI_INTEGER8,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_owner_keys,send_counts,send_displacements,MPI_INTEGER,receive_owner_keys,&
        receive_counts,receive_displacements,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call scale_counts_checked(comm,send_counts,receive_counts,value_layout,scaled_send_counts,&
        scaled_receive_counts,scaled_send_displacements,scaled_receive_displacements,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Alltoallv(send_values,scaled_send_counts,scaled_send_displacements,MPI_DOUBLE_COMPLEX,&
        receive_values,scaled_receive_counts,scaled_receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      do i=1,nreceived;record_order(i)=i;enddo
      call sort_id_positions(receive_ids,record_order)
      local_bad=0
      do i=1,nreceived
        if(receive_ids(i)<=0_int64)then
          local_bad=1
        else if(mod(receive_ids(i)-1_int64,int(nproc,int64))/=int(rank,int64))then
          local_bad=1
        endif
        if(i>1)then
          if(receive_ids(i)==receive_ids(i-1))local_bad=1
        endif
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_bad/=0)then;ierr=1;return;endif
      deallocate(target%nonlocal_ids,target%nonlocal_owner,target%nonlocal_values)
      allocate(target%nonlocal_ids(nreceived),target%nonlocal_owner(nreceived),&
        target%nonlocal_values(value_layout,nreceived))
      do k=1,nreceived
        source_index=record_order(k);target%nonlocal_ids(k)=receive_ids(k)
        target%nonlocal_owner(k)=receive_owner_keys(source_index)
        target%nonlocal_values(:,k)=receive_values(:,source_index)
      enddo
    end subroutine redistribute_nonlocal_payload

    subroutine scale_counts_checked(comm,base_send_counts,base_receive_counts,scale,send_counts,receive_counts,&
        send_displacements,receive_displacements,ierr)
      integer,intent(in)::comm,base_send_counts(:),base_receive_counts(:),scale
      integer,intent(out)::send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:),ierr
      integer::i,local_bad,global_bad
      integer(int64)::send_prefix,receive_prefix,scaled,limit
      send_counts=0;receive_counts=0;send_displacements=0;receive_displacements=0
      local_bad=0;send_prefix=0_int64;receive_prefix=0_int64;limit=int(huge(0),int64)
      if(scale<0.or.size(base_send_counts)/=size(send_counts).or.&
        size(base_receive_counts)/=size(receive_counts).or.&
        size(send_displacements)/=size(send_counts).or.&
        size(receive_displacements)/=size(receive_counts))local_bad=1
      if(local_bad==0)then
        do i=1,size(send_counts)
          if(base_send_counts(i)<0.or.base_receive_counts(i)<0)then
            local_bad=1;cycle
          endif
          scaled=int(base_send_counts(i),int64)*int(scale,int64)
          if(scaled>limit.or.send_prefix>limit-scaled)then
            local_bad=1
          else
            send_displacements(i)=int(send_prefix);send_counts(i)=int(scaled);send_prefix=send_prefix+scaled
          endif
          scaled=int(base_receive_counts(i),int64)*int(scale,int64)
          if(scaled>limit.or.receive_prefix>limit-scaled)then
            local_bad=1
          else
            receive_displacements(i)=int(receive_prefix);receive_counts(i)=int(scaled)
            receive_prefix=receive_prefix+scaled
          endif
        enddo
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr==MPI_SUCCESS.and.global_bad/=0)ierr=1
    end subroutine scale_counts_checked

    subroutine sort_id_positions(ids,positions)
      integer(int64),intent(inout)::ids(:)
      integer,intent(inout)::positions(:)
      integer(int64),allocatable::scratch_ids(:)
      integer,allocatable::scratch_positions(:)
      integer(int64)::left,middle,right,width,source_index,other_index,destination_index,nvalues
      nvalues=size(ids,kind=int64);if(nvalues<2_int64)return
      allocate(scratch_ids(size(ids)),scratch_positions(size(positions)));width=1_int64
      do
        scratch_ids=ids;scratch_positions=positions;left=1_int64
        do while(left<=nvalues)
          middle=min(left+width-1_int64,nvalues);right=min(left+2_int64*width-1_int64,nvalues)
          if(middle<right)then
            source_index=left;other_index=middle+1_int64
            do destination_index=left,right
              if(source_index>middle)then
                scratch_ids(destination_index)=ids(other_index)
                scratch_positions(destination_index)=positions(other_index);other_index=other_index+1_int64
              else if(other_index>right)then
                scratch_ids(destination_index)=ids(source_index)
                scratch_positions(destination_index)=positions(source_index);source_index=source_index+1_int64
              else if(ids(source_index)<=ids(other_index))then
                scratch_ids(destination_index)=ids(source_index)
                scratch_positions(destination_index)=positions(source_index);source_index=source_index+1_int64
              else
                scratch_ids(destination_index)=ids(other_index)
                scratch_positions(destination_index)=positions(other_index);other_index=other_index+1_int64
              endif
            enddo
          endif
          left=left+2_int64*width
        enddo
        ids=scratch_ids;positions=scratch_positions
        if(width>=nvalues-width)exit
        width=2_int64*width
      enddo
    end subroutine sort_id_positions

    subroutine replace_distributed_payload(destination_payload,source_payload)
      type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::destination_payload,source_payload
      call move_alloc(source_payload%row_ids,destination_payload%row_ids)
      call move_alloc(source_payload%grid_ids,destination_payload%grid_ids)
      call move_alloc(source_payload%face_ids,destination_payload%face_ids)
      call move_alloc(source_payload%face_point_ids,destination_payload%face_point_ids)
      call move_alloc(source_payload%nonlocal_ids,destination_payload%nonlocal_ids)
      call move_alloc(source_payload%metric_row_offsets,destination_payload%metric_row_offsets)
      call move_alloc(source_payload%metric_column_ids,destination_payload%metric_column_ids)
      call move_alloc(source_payload%operator_row_offsets,destination_payload%operator_row_offsets)
      call move_alloc(source_payload%operator_column_ids,destination_payload%operator_column_ids)
      call move_alloc(source_payload%partition_ids,destination_payload%partition_ids)
      call move_alloc(source_payload%face_metadata,destination_payload%face_metadata)
      call move_alloc(source_payload%face_offsets,destination_payload%face_offsets)
      call move_alloc(source_payload%face_weight_offsets,destination_payload%face_weight_offsets)
      call move_alloc(source_payload%face_basis_offsets,destination_payload%face_basis_offsets)
      call move_alloc(source_payload%face_value_offsets,destination_payload%face_value_offsets)
      call move_alloc(source_payload%face_observable_offsets,destination_payload%face_observable_offsets)
      call move_alloc(source_payload%face_basis_ids,destination_payload%face_basis_ids)
      call move_alloc(source_payload%nonlocal_owner,destination_payload%nonlocal_owner)
      call move_alloc(source_payload%grid_weights,destination_payload%grid_weights)
      call move_alloc(source_payload%face_normals,destination_payload%face_normals)
      call move_alloc(source_payload%face_weights,destination_payload%face_weights)
      call move_alloc(source_payload%density,destination_payload%density)
      call move_alloc(source_payload%metric_rows,destination_payload%metric_rows)
      call move_alloc(source_payload%kinetic_rows,destination_payload%kinetic_rows)
      call move_alloc(source_payload%nonlocal_rows,destination_payload%nonlocal_rows)
      call move_alloc(source_payload%local_rows,destination_payload%local_rows)
      call move_alloc(source_payload%sipg_rows,destination_payload%sipg_rows)
      call move_alloc(source_payload%hamiltonian_rows,destination_payload%hamiltonian_rows)
      call move_alloc(source_payload%basis_values,destination_payload%basis_values)
      call move_alloc(source_payload%face_values,destination_payload%face_values)
      call move_alloc(source_payload%nonlocal_values,destination_payload%nonlocal_values)
      call move_alloc(source_payload%coefficients,destination_payload%coefficients)
      call move_alloc(source_payload%interface_observables,destination_payload%interface_observables)
      call move_alloc(source_payload%position_rows,destination_payload%position_rows)
      call move_alloc(source_payload%certified_basis%construction_row_ids,&
        destination_payload%certified_basis%construction_row_ids)
      call move_alloc(source_payload%certified_basis%transformation_row_ids,&
        destination_payload%certified_basis%transformation_row_ids)
      call move_alloc(source_payload%certified_basis%c_cert,destination_payload%certified_basis%c_cert)
      call move_alloc(source_payload%certified_basis%u_rt,destination_payload%certified_basis%u_rt)
      call move_alloc(source_payload%certified_basis%b_rt,destination_payload%certified_basis%b_rt)
      call move_alloc(source_payload%rt_space%row_ids,destination_payload%rt_space%row_ids)
      call move_alloc(source_payload%rt_space%grid_owner_keys,destination_payload%rt_space%grid_owner_keys)
      call move_alloc(source_payload%rt_space%metric_rows,destination_payload%rt_space%metric_rows)
      call move_alloc(source_payload%rt_space%kinetic_rows,destination_payload%rt_space%kinetic_rows)
      call move_alloc(source_payload%rt_space%nonlocal_rows,destination_payload%rt_space%nonlocal_rows)
      call move_alloc(source_payload%rt_space%local_rows,destination_payload%rt_space%local_rows)
      call move_alloc(source_payload%rt_space%sipg_rows,destination_payload%rt_space%sipg_rows)
      call move_alloc(source_payload%rt_space%hamiltonian_rows,destination_payload%rt_space%hamiltonian_rows)
      call move_alloc(source_payload%rt_space%scalar_operator_rows,destination_payload%rt_space%scalar_operator_rows)
      call move_alloc(source_payload%rt_space%vector_operator_rows,destination_payload%rt_space%vector_operator_rows)
      call move_alloc(source_payload%rt_space%tensor_operator_rows,destination_payload%rt_space%tensor_operator_rows)
      call move_alloc(source_payload%rt_space%basis_values,destination_payload%rt_space%basis_values)
      call move_alloc(source_payload%rt_space%density,destination_payload%rt_space%density)
      source_payload=s_rt_dg_hybrid_ground_state_payload()
    end subroutine replace_distributed_payload

  end subroutine redistribute_coalesced_ground_state

  subroutine copy_ground_state_common(source,target)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::source
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::target
    logical::l(ground_state_logical_count);integer::h(ground_state_integer_count)
    integer(int64)::f(ground_state_fingerprint_count);real(real64)::r(ground_state_real_count)
    call pack_ground_state_header(source,l,h,f,r);call set_ground_state_header(target,l,h,f,r)
    allocate(target%construction_catalog%ids,source=source%construction_catalog%ids)
    allocate(target%construction_catalog%generations,source=source%construction_catalog%generations)
    allocate(target%construction_catalog%ordering,source=source%construction_catalog%ordering)
    allocate(target%construction_catalog%ownership,source=source%construction_catalog%ownership)
    allocate(target%certified_basis%construction_row_ids(0),&
      target%certified_basis%transformation_row_ids(0),&
      target%certified_basis%c_cert(0,source%certified_basis%certified_count),&
      target%certified_basis%u_rt(0,source%certified_basis%certified_count),&
      target%certified_basis%b_rt(0,source%certified_basis%certified_count))
    allocate(target%certified_basis%initial_occupied_amplitudes,&
      source=source%certified_basis%initial_occupied_amplitudes)
    allocate(target%certified_basis%certified_eigenvalues,source=source%certified_basis%certified_eigenvalues)
    allocate(target%certified_basis%occupations,source=source%certified_basis%occupations)
    allocate(target%certified_basis%centers,source=source%certified_basis%centers)
    allocate(target%certified_basis%spreads_before,source=source%certified_basis%spreads_before)
    allocate(target%certified_basis%spreads_after,source=source%certified_basis%spreads_after)
    allocate(target%rt_space%row_ids(0),target%rt_space%grid_owner_keys(0),&
      target%rt_space%metric_rows(0,source%rt_space%rank),&
      target%rt_space%kinetic_rows(0,source%rt_space%rank),&
      target%rt_space%nonlocal_rows(0,source%rt_space%rank),&
      target%rt_space%local_rows(0,source%rt_space%rank),&
      target%rt_space%sipg_rows(0,source%rt_space%rank),&
      target%rt_space%hamiltonian_rows(0,source%rt_space%rank),&
      target%rt_space%scalar_operator_rows(0,source%rt_space%rank,source%rt_space%scalar_count),&
      target%rt_space%vector_operator_rows(0,source%rt_space%rank,3,source%rt_space%vector_count),&
      target%rt_space%tensor_operator_rows(0,source%rt_space%rank,3,3,source%rt_space%tensor_count),&
      target%rt_space%basis_values(source%rt_space%rank,0),target%rt_space%density(0))
    allocate(target%rt_space%row_owner_keys,source=source%rt_space%row_owner_keys)
    allocate(target%rt_space%representation,source=source%rt_space%representation)
    allocate(target%rt_space%cartesian_rotations,source=source%rt_space%cartesian_rotations)
    allocate(target%scope_selectors,source=source%scope_selectors);allocate(target%xc_types,source=source%xc_types)
    allocate(target%requested_ids,source=source%requested_ids);allocate(target%effective_ids,source=source%effective_ids)
    allocate(target%added_ids,source=source%added_ids);allocate(target%closure_parent,source=source%closure_parent)
    allocate(target%closure_reason,source=source%closure_reason);allocate(target%closure_action,source=source%closure_action)
    allocate(target%occupations,source=source%occupations);allocate(target%eigenvalues,source=source%eigenvalues)
    allocate(target%continuation_receipt,source=source%continuation_receipt)
    allocate(target%pseudopotential_receipt,source=source%pseudopotential_receipt)
    allocate(target%energy_receipt,source=source%energy_receipt)
    allocate(target%symmetry_representation,source=source%symmetry_representation)
    allocate(target%row_ids(0),target%grid_ids(0),target%grid_weights(0),target%density(0))
    allocate(target%partition_ids(0))
    allocate(target%face_ids(0),target%face_point_ids(0),target%nonlocal_ids(0),target%nonlocal_owner(0))
    allocate(target%face_metadata(size(source%face_metadata,1),0),target%face_normals(size(source%face_normals,1),0),&
      target%face_weights(0),target%face_values(size(source%face_values,1),0),&
      target%interface_observables(size(source%interface_observables,1),0),&
      target%nonlocal_values(size(source%nonlocal_values,1),0),target%face_offsets(1),&
      target%face_weight_offsets(1),target%face_basis_offsets(1),target%face_value_offsets(1),&
      target%face_observable_offsets(1),target%face_basis_ids(0))
    target%face_offsets=1;target%face_weight_offsets=1;target%face_basis_offsets=1
    target%face_value_offsets=1;target%face_observable_offsets=1
    allocate(target%metric_rows(0,source%global_count),target%kinetic_rows(0,source%global_count),&
      target%nonlocal_rows(0,source%global_count),target%local_rows(0,source%global_count),&
      target%sipg_rows(0,source%global_count),target%hamiltonian_rows(0,source%global_count),&
      target%coefficients(0,source%noccupied),target%position_rows(3,0,source%global_count),&
      target%basis_values(source%global_count,0),target%metric_row_offsets(1),target%metric_column_ids(0),&
      target%operator_row_offsets(1),target%operator_column_ids(0))
    target%metric_row_offsets=1;target%operator_row_offsets=1
  end subroutine copy_ground_state_common

  subroutine append_initialization_shard(target,source)
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::target
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::source
    call append_i64(target%row_ids,source%row_ids);call append_i64(target%grid_ids,source%grid_ids)
    call append_i1(target%partition_ids,source%partition_ids)
    call append_i64(target%face_ids,source%face_ids)
    call append_i64(target%nonlocal_ids,source%nonlocal_ids);call append_i1(target%nonlocal_owner,source%nonlocal_owner)
    call append_i2_columns(target%face_metadata,source%face_metadata)
    call append_r2_columns(target%face_normals,source%face_normals)
    call append_z2_columns(target%nonlocal_values,source%nonlocal_values)
    call append_csr_i64(target%face_offsets,target%face_point_ids,source%face_offsets,source%face_point_ids)
    call append_csr_r(target%face_weight_offsets,target%face_weights,&
      source%face_weight_offsets,source%face_weights)
    call append_csr(target%face_basis_offsets,target%face_basis_ids,&
      source%face_basis_offsets,source%face_basis_ids)
    call append_csr_z2(target%face_value_offsets,target%face_values,&
      source%face_value_offsets,source%face_values)
    call append_csr_z2(target%face_observable_offsets,target%interface_observables,&
      source%face_observable_offsets,source%interface_observables)
    call append_r1(target%grid_weights,source%grid_weights);call append_r1(target%density,source%density)
    call append_z2_rows(target%metric_rows,source%metric_rows);call append_z2_rows(target%kinetic_rows,source%kinetic_rows)
    call append_z2_rows(target%nonlocal_rows,source%nonlocal_rows);call append_z2_rows(target%local_rows,source%local_rows)
    call append_z2_rows(target%sipg_rows,source%sipg_rows);call append_z2_rows(target%hamiltonian_rows,source%hamiltonian_rows)
    call append_z2_rows(target%coefficients,source%coefficients);call append_z2_columns(target%basis_values,source%basis_values)
    call append_z3_middle(target%position_rows,source%position_rows)
    call append_csr(target%metric_row_offsets,target%metric_column_ids,source%metric_row_offsets,source%metric_column_ids)
    call append_csr(target%operator_row_offsets,target%operator_column_ids,source%operator_row_offsets,source%operator_column_ids)
    call append_i64(target%certified_basis%construction_row_ids,source%certified_basis%construction_row_ids)
    call append_z2_rows(target%certified_basis%c_cert,source%certified_basis%c_cert)
    call append_i64(target%certified_basis%transformation_row_ids,source%certified_basis%transformation_row_ids)
    call append_z2_rows(target%certified_basis%u_rt,source%certified_basis%u_rt)
    call append_z2_rows(target%certified_basis%b_rt,source%certified_basis%b_rt)
    call append_i64(target%rt_space%row_ids,source%rt_space%row_ids)
    call append_i1(target%rt_space%grid_owner_keys,source%rt_space%grid_owner_keys)
    call append_z2_rows(target%rt_space%metric_rows,source%rt_space%metric_rows)
    call append_z2_rows(target%rt_space%kinetic_rows,source%rt_space%kinetic_rows)
    call append_z2_rows(target%rt_space%nonlocal_rows,source%rt_space%nonlocal_rows)
    call append_z2_rows(target%rt_space%local_rows,source%rt_space%local_rows)
    call append_z2_rows(target%rt_space%sipg_rows,source%rt_space%sipg_rows)
    call append_z2_rows(target%rt_space%hamiltonian_rows,source%rt_space%hamiltonian_rows)
    call append_z3_rows(target%rt_space%scalar_operator_rows,source%rt_space%scalar_operator_rows)
    call append_z4_rows(target%rt_space%vector_operator_rows,source%rt_space%vector_operator_rows)
    call append_z5_rows(target%rt_space%tensor_operator_rows,source%rt_space%tensor_operator_rows)
    call append_z2_columns(target%rt_space%basis_values,source%rt_space%basis_values)
    call append_r1(target%rt_space%density,source%rt_space%density)
  contains
    subroutine append_i64(a,b)
      integer(int64),allocatable,intent(inout)::a(:);integer(int64),intent(in)::b(:);integer(int64),allocatable::t(:)
      allocate(t(size(a)+size(b)));t(:size(a))=a;t(size(a)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_i1(a,b)
      integer,allocatable,intent(inout)::a(:);integer,intent(in)::b(:);integer,allocatable::t(:)
      allocate(t(size(a)+size(b)));t(:size(a))=a;t(size(a)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_i2_columns(a,b)
      integer,allocatable,intent(inout)::a(:,:);integer,intent(in)::b(:,:);integer,allocatable::t(:,:)
      allocate(t(size(a,1),size(a,2)+size(b,2)));t(:,:size(a,2))=a;t(:,size(a,2)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_r2_columns(a,b)
      real(real64),allocatable,intent(inout)::a(:,:);real(real64),intent(in)::b(:,:);real(real64),allocatable::t(:,:)
      allocate(t(size(a,1),size(a,2)+size(b,2)));t(:,:size(a,2))=a;t(:,size(a,2)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_csr_i64(offsets,columns,new_offsets,new_columns)
      integer,allocatable,intent(inout)::offsets(:)
      integer(int64),allocatable,intent(inout)::columns(:)
      integer,intent(in)::new_offsets(:)
      integer(int64),intent(in)::new_columns(:)
      integer,allocatable::next_offsets(:)
      integer(int64),allocatable::next_columns(:)
      integer::old_rows,old_columns
      old_rows=size(offsets)-1;old_columns=size(columns)
      allocate(next_offsets(old_rows+size(new_offsets)),next_columns(old_columns+size(new_columns)))
      next_offsets(:old_rows+1)=offsets
      next_offsets(old_rows+2:)=old_columns+new_offsets(2:)
      next_columns(:old_columns)=columns;next_columns(old_columns+1:)=new_columns
      call move_alloc(next_offsets,offsets);call move_alloc(next_columns,columns)
    end subroutine
    subroutine append_r1(a,b)
      real(real64),allocatable,intent(inout)::a(:);real(real64),intent(in)::b(:);real(real64),allocatable::t(:)
      allocate(t(size(a)+size(b)));t(:size(a))=a;t(size(a)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_csr_r(offsets,values,new_offsets,new_values)
      integer,allocatable,intent(inout)::offsets(:)
      real(real64),allocatable,intent(inout)::values(:)
      integer,intent(in)::new_offsets(:)
      real(real64),intent(in)::new_values(:)
      integer,allocatable::next_offsets(:)
      real(real64),allocatable::next_values(:)
      integer::old_rows,old_nnz
      old_rows=size(offsets)-1;old_nnz=size(values)
      allocate(next_offsets(old_rows+size(new_offsets)),next_values(old_nnz+size(new_values)))
      next_offsets(:old_rows+1)=offsets;next_offsets(old_rows+2:)=old_nnz+new_offsets(2:)
      next_values(:old_nnz)=values;next_values(old_nnz+1:)=new_values
      call move_alloc(next_offsets,offsets);call move_alloc(next_values,values)
    end subroutine append_csr_r
    subroutine append_csr_z2(offsets,values,new_offsets,new_values)
      integer,allocatable,intent(inout)::offsets(:)
      complex(real64),allocatable,intent(inout)::values(:,:)
      integer,intent(in)::new_offsets(:)
      complex(real64),intent(in)::new_values(:,:)
      integer,allocatable::next_offsets(:)
      complex(real64),allocatable::next_values(:,:)
      integer::old_rows,old_nnz
      old_rows=size(offsets)-1;old_nnz=size(values,2)
      allocate(next_offsets(old_rows+size(new_offsets)),next_values(size(values,1),old_nnz+size(new_values,2)))
      next_offsets(:old_rows+1)=offsets;next_offsets(old_rows+2:)=old_nnz+new_offsets(2:)
      next_values(:,:old_nnz)=values;next_values(:,old_nnz+1:)=new_values
      call move_alloc(next_offsets,offsets);call move_alloc(next_values,values)
    end subroutine append_csr_z2
    subroutine append_z2_rows(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:);complex(real64),intent(in)::b(:,:);complex(real64),allocatable::t(:,:)
      allocate(t(size(a,1)+size(b,1),size(a,2)));t(:size(a,1),:)=a;t(size(a,1)+1:,:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_z2_columns(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:);complex(real64),intent(in)::b(:,:);complex(real64),allocatable::t(:,:)
      allocate(t(size(a,1),size(a,2)+size(b,2)));t(:,:size(a,2))=a;t(:,size(a,2)+1:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_z3_middle(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:,:);complex(real64),intent(in)::b(:,:,:);complex(real64),allocatable::t(:,:,:)
      allocate(t(size(a,1),size(a,2)+size(b,2),size(a,3)));t(:,:size(a,2),:)=a;t(:,size(a,2)+1:,:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_z3_rows(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:,:);complex(real64),intent(in)::b(:,:,:)
      complex(real64),allocatable::t(:,:,:)
      allocate(t(size(a,1)+size(b,1),size(a,2),size(a,3)))
      t(:size(a,1),:,:)=a;t(size(a,1)+1:,:,:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_z4_rows(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:,:,:);complex(real64),intent(in)::b(:,:,:,:)
      complex(real64),allocatable::t(:,:,:,:)
      allocate(t(size(a,1)+size(b,1),size(a,2),size(a,3),size(a,4)))
      t(:size(a,1),:,:,:)=a;t(size(a,1)+1:,:,:,:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_z5_rows(a,b)
      complex(real64),allocatable,intent(inout)::a(:,:,:,:,:);complex(real64),intent(in)::b(:,:,:,:,:)
      complex(real64),allocatable::t(:,:,:,:,:)
      allocate(t(size(a,1)+size(b,1),size(a,2),size(a,3),size(a,4),size(a,5)))
      t(:size(a,1),:,:,:,:)=a;t(size(a,1)+1:,:,:,:,:)=b;call move_alloc(t,a)
    end subroutine
    subroutine append_csr(offsets,columns,new_offsets,new_columns)
      integer,allocatable,intent(inout)::offsets(:),columns(:);integer,intent(in)::new_offsets(:),new_columns(:)
      integer,allocatable::to(:),tc(:);integer::old_rows,old_nnz
      old_rows=size(offsets)-1;old_nnz=size(columns);allocate(to(old_rows+size(new_offsets)),tc(old_nnz+size(new_columns)))
      to(:old_rows+1)=offsets;to(old_rows+2:)=new_offsets(2:)+old_nnz
      tc(:old_nnz)=columns;tc(old_nnz+1:)=new_columns;call move_alloc(to,offsets);call move_alloc(tc,columns)
    end subroutine
  end subroutine append_initialization_shard
#endif

#ifdef USE_MPI
  subroutine write_ground_state_arrays(comm,owner,unit,rank,payload,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer,intent(inout)::io_status
    integer,intent(out)::ierr
    call stream_write_i64_1(comm,owner,unit,rank,payload%row_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%grid_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%face_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%face_point_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%nonlocal_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%partition_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%metric_row_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%metric_column_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%operator_row_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%operator_column_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i2(comm,owner,unit,rank,payload%face_metadata,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_weight_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_basis_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_value_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_observable_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%face_basis_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%nonlocal_owner,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%requested_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%effective_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%added_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%closure_parent,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%closure_reason,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%closure_action,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%scope_selectors,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%xc_types,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%grid_weights,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r2(comm,owner,unit,rank,payload%face_normals,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%face_weights,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%density,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%occupations,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%eigenvalues,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%continuation_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%pseudopotential_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%energy_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%metric_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%kinetic_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%nonlocal_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%local_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%sipg_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%hamiltonian_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%basis_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%face_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%nonlocal_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%coefficients,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%interface_observables,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z3(comm,owner,unit,rank,payload%position_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z3(comm,owner,unit,rank,payload%symmetry_representation,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%construction_catalog%ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%construction_catalog%generations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%construction_catalog%ordering,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%construction_catalog%ownership,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%certified_basis%construction_row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%certified_basis%c_cert,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%certified_basis%transformation_row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%certified_basis%u_rt,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%certified_basis%b_rt,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%certified_basis%initial_occupied_amplitudes,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%certified_basis%certified_eigenvalues,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%certified_basis%occupations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r2(comm,owner,unit,rank,payload%certified_basis%centers,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%certified_basis%spreads_before,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%certified_basis%spreads_after,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i64_1(comm,owner,unit,rank,payload%rt_space%row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%rt_space%row_owner_keys,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_i1(comm,owner,unit,rank,payload%rt_space%grid_owner_keys,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%metric_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%kinetic_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%nonlocal_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%local_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%sipg_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%hamiltonian_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z3(comm,owner,unit,rank,payload%rt_space%representation,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r3(comm,owner,unit,rank,payload%rt_space%cartesian_rotations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z3(comm,owner,unit,rank,payload%rt_space%scalar_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z4(comm,owner,unit,rank,payload%rt_space%vector_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z5(comm,owner,unit,rank,payload%rt_space%tensor_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_z2(comm,owner,unit,rank,payload%rt_space%basis_values,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_write_r1(comm,owner,unit,rank,payload%rt_space%density,io_status,ierr)
  end subroutine write_ground_state_arrays

  subroutine read_ground_state_arrays(comm,owner,unit,rank,payload,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::payload
    integer,intent(inout)::io_status
    integer,intent(out)::ierr
    call stream_read_i64_1(comm,owner,unit,rank,payload%row_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%grid_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%face_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%face_point_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%nonlocal_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%partition_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%metric_row_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%metric_column_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%operator_row_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%operator_column_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i2(comm,owner,unit,rank,payload%face_metadata,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_weight_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_basis_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_value_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_observable_offsets,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%face_basis_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%nonlocal_owner,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%requested_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%effective_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%added_ids,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%closure_parent,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%closure_reason,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%closure_action,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%scope_selectors,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%xc_types,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%grid_weights,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r2(comm,owner,unit,rank,payload%face_normals,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%face_weights,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%density,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%occupations,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%eigenvalues,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%continuation_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%pseudopotential_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%energy_receipt,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%metric_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%kinetic_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%nonlocal_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%local_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%sipg_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%hamiltonian_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%basis_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%face_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%nonlocal_values,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%coefficients,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%interface_observables,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z3(comm,owner,unit,rank,payload%position_rows,io_status,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z3(comm,owner,unit,rank,payload%symmetry_representation,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%construction_catalog%ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%construction_catalog%generations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%construction_catalog%ordering,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%construction_catalog%ownership,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%certified_basis%construction_row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%certified_basis%c_cert,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%certified_basis%transformation_row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%certified_basis%u_rt,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%certified_basis%b_rt,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%certified_basis%initial_occupied_amplitudes,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%certified_basis%certified_eigenvalues,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%certified_basis%occupations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r2(comm,owner,unit,rank,payload%certified_basis%centers,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%certified_basis%spreads_before,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%certified_basis%spreads_after,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i64_1(comm,owner,unit,rank,payload%rt_space%row_ids,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%rt_space%row_owner_keys,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_i1(comm,owner,unit,rank,payload%rt_space%grid_owner_keys,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%metric_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%kinetic_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%nonlocal_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%local_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%sipg_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%hamiltonian_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z3(comm,owner,unit,rank,payload%rt_space%representation,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r3(comm,owner,unit,rank,payload%rt_space%cartesian_rotations,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z3(comm,owner,unit,rank,payload%rt_space%scalar_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z4(comm,owner,unit,rank,payload%rt_space%vector_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z5(comm,owner,unit,rank,payload%rt_space%tensor_operator_rows,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_z2(comm,owner,unit,rank,payload%rt_space%basis_values,io_status,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call stream_read_r1(comm,owner,unit,rank,payload%rt_space%density,io_status,ierr)
  end subroutine read_ground_state_arrays
#endif

  subroutine write_rt_dg_hybrid_occupied_checkpoint(comm,path,global_count,row_ids,coefficients,occupations,eigenvalues,&
      catalog_fingerprint,basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,scf_receipts,&
      maximum_scf_residual,payload_fingerprint,ok,message)
    integer,intent(in)::comm,global_count
    character(*),intent(in)::path
    integer(int64),intent(in)::row_ids(:),catalog_fingerprint,basis_fingerprint,provenance_fingerprints(6),&
      operator_fingerprint,state_fingerprint
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),eigenvalues(:),scf_receipts(5),maximum_scf_residual
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,noccupied,nlocal,i,j,row,root,position,unit,io_status,close_status,bad,global_bad,status
    integer,allocatable::counts(:),owners(:),positions(:)
    complex(real64),allocatable::row_values(:)
    integer(int64)::minimum_i,maximum_i,bits
    real(real64)::minimum_r,maximum_r
    character(16)::probe
    character(:),allocatable::temporary_path
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    nlocal=size(row_ids);noccupied=size(coefficients,2);bad=0
    if(global_count<1.or.noccupied<1.or.size(coefficients,1)/=nlocal.or.size(occupations)/=noccupied.or.&
      size(eigenvalues)/=noccupied)bad=1
    call agree_int(global_count);call agree_int(noccupied)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint dimensions';return;endif
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or..not.finite_matrix(coefficients).or.&
      any(.not.ieee_is_finite(occupations)).or.any(.not.ieee_is_finite(eigenvalues)).or.&
      any(.not.ieee_is_finite(scf_receipts)).or.any(occupations<0d0))bad=1
    if(catalog_fingerprint==0_int64.or.basis_fingerprint==0_int64.or.any(provenance_fingerprints==0_int64).or.&
      operator_fingerprint==0_int64.or.state_fingerprint==0_int64)bad=1
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0.or.any(scf_receipts<0d0).or.&
      any(scf_receipts>maximum_scf_residual))bad=1
    call agree_i64(catalog_fingerprint);call agree_i64(basis_fingerprint);call agree_i64(operator_fingerprint);call agree_i64(state_fingerprint)
    do i=1,6;call agree_i64(provenance_fingerprints(i));enddo
    call agree_real(maximum_scf_residual)
    do i=1,noccupied;call agree_real(occupations(i));call agree_real(eigenvalues(i));enddo
    do i=1,5;call agree_real(scf_receipts(i));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint write contract';return;endif
    allocate(counts(global_count),owners(global_count),positions(global_count),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='cannot allocate occupied checkpoint workspace';return;endif
    counts=0;owners=-1;positions=0
    do i=1,nlocal;row=int(row_ids(i));counts(row)=counts(row)+1;owners(row)=rank;positions(row)=i;enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,positions,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(counts==1));call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='occupied checkpoint rows are not exactly once';return;endif
    temporary_path=trim(path)//'.occupied.tmp.'//trim(int64_string(state_fingerprint))
    if(rank==0)then
      open(newunit=unit,file=temporary_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status);opened=io_status==0
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    payload_fingerprint=catalog_fingerprint
    call occupied_hash(basis_fingerprint);do i=1,6;call occupied_hash(provenance_fingerprints(i));enddo
    call occupied_hash(operator_fingerprint);call occupied_hash(state_fingerprint)
    call occupied_hash(int(global_count,int64));call occupied_hash(int(noccupied,int64))
    do i=1,noccupied;call occupied_hash(transfer(occupations(i),bits));call occupied_hash(transfer(eigenvalues(i),bits));enddo
    do i=1,5;call occupied_hash(transfer(scf_receipts(i),bits));enddo
    if(rank==0)write(unit,iostat=io_status)occupied_magic,occupied_version,global_count,noccupied,catalog_fingerprint,&
      basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,occupations,eigenvalues,scf_receipts
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do row=1,global_count
      root=owners(row);position=positions(row);if(rank==root)row_values=coefficients(position,:)
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call occupied_hash(int(row,int64));do j=1,noccupied;call occupied_hash(transfer(real(row_values(j)),bits));call occupied_hash(transfer(aimag(row_values(j)),bits));enddo
      if(rank==0)write(unit,iostat=io_status)row_values
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)then
      write(unit,iostat=io_status)payload_fingerprint;close_status=0;close(unit,iostat=close_status);opened=.false.;if(io_status==0)io_status=close_status
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call clean;ok=.true.;return
900 message='occupied checkpoint MPI stream failed';if(rank==0.and.opened)close(unit);call clean;return
910 message='occupied checkpoint publication failed';if(rank==0.and.opened)close(unit);call clean;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_int(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_i64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,minimum_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minimum_i/=maximum_i)bad=1
    end subroutine
    subroutine agree_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,minimum_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.transfer(minimum_r,bits)/=transfer(maximum_r,bits))bad=1
    end subroutine
    subroutine occupied_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine clean
      if(allocated(counts))deallocate(counts);if(allocated(owners))deallocate(owners)
      if(allocated(positions))deallocate(positions);if(allocated(row_values))deallocate(row_values)
    end subroutine
#endif
  end subroutine write_rt_dg_hybrid_occupied_checkpoint

  subroutine read_rt_dg_hybrid_occupied_checkpoint(comm,path,expected_catalog,expected_basis,expected_operator,&
      expected_state,expected_provenance,expected_occupations,maximum_scf_residual,global_count,row_ids,coefficients,&
      occupations,eigenvalues,scf_receipts,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::expected_catalog,expected_basis,expected_operator,expected_state,expected_provenance(6)
    real(real64),intent(in)::expected_occupations(:),maximum_scf_residual
    integer,intent(out)::global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),allocatable,intent(out)::occupations(:),eigenvalues(:)
    real(real64),intent(out)::scf_receipts(5)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,bad,global_bad,noccupied,version,nlocal,row,position,j,status
    integer(int64)::catalog,basis,provenance(6),operator_receipt,state,stored_fingerprint,bits,lo_i,hi_i,file_size
    real(real64)::lo_r,hi_r
    complex(real64),allocatable::row_values(:)
    character(16)::magic,probe
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    bad=0
    call agree_expected(expected_catalog);call agree_expected(expected_basis);call agree_expected(expected_operator);call agree_expected(expected_state)
    do j=1,6;call agree_expected(expected_provenance(j));enddo
    call agree_expected_count(size(expected_occupations))
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    call agree_expected_real(maximum_scf_residual)
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0)bad=1
    do j=1,size(expected_occupations);call agree_expected_real(expected_occupations(j));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status);opened=io_status==0
      if(io_status==0)inquire(unit=unit,size=file_size,iostat=io_status)
      if(io_status==0)read(unit,iostat=io_status)magic,version,global_count,noccupied,catalog,basis,provenance,operator_receipt,state
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(magic,len(magic),MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(version,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(global_count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(noccupied,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(catalog,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(basis,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(operator_receipt,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(state,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(provenance,6,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=0
    if(magic/=occupied_magic.or.version/=occupied_version.or.global_count<1.or.noccupied<1)bad=1
    if(catalog/=expected_catalog.or.basis/=expected_basis.or.any(provenance/=expected_provenance).or.&
      operator_receipt/=expected_operator.or.state/=expected_state)bad=1
    if(global_count>huge(0)-nproc.or.noccupied>huge(0)/max(1,global_count))bad=1
    if(file_size<1_int64.or.int(noccupied,int64)>file_size/16_int64.or.&
      int(global_count,int64)>file_size/max(16_int64,16_int64*int(noccupied,int64)))bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='stale or incompatible occupied checkpoint';goto 920;endif
    nlocal=(global_count+nproc-1-rank)/nproc
    allocate(row_ids(nlocal),coefficients(nlocal,noccupied),occupations(noccupied),eigenvalues(noccupied),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate occupied checkpoint output';goto 920;endif
    if(rank==0)read(unit,iostat=io_status)occupations,eigenvalues,scf_receipts
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(occupations,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(eigenvalues,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(scf_receipts,5,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(ieee_is_finite(occupations)).and.all(occupations>=0d0).and.&
      all(ieee_is_finite(eigenvalues)).and.all(ieee_is_finite(scf_receipts)).and.all(scf_receipts>=0d0).and.&
      all(scf_receipts<=maximum_scf_residual).and.size(expected_occupations)==noccupied)
    if(bad==0)then;if(any(occupations/=expected_occupations))bad=1;endif
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    payload_fingerprint=catalog;call read_hash(basis);do j=1,6;call read_hash(provenance(j));enddo
    call read_hash(operator_receipt);call read_hash(state)
    call read_hash(int(global_count,int64));call read_hash(int(noccupied,int64))
    do j=1,noccupied;call read_hash(transfer(occupations(j),bits));call read_hash(transfer(eigenvalues(j),bits));enddo
    do j=1,5;call read_hash(transfer(scf_receipts(j),bits));enddo
    position=0
    do row=1,global_count
      if(rank==0)read(unit,iostat=io_status)row_values
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      if(.not.finite_vector(row_values))goto 930
      call read_hash(int(row,int64));do j=1,noccupied;call read_hash(transfer(real(row_values(j)),bits));call read_hash(transfer(aimag(row_values(j)),bits));enddo
      if(mod(row-1,nproc)==rank)then;position=position+1;row_ids(position)=row;coefficients(position,:)=row_values;endif
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)read(unit,iostat=io_status)stored_fingerprint
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(stored_fingerprint,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    if(rank==0)then;close(unit,iostat=io_status);opened=.false.;endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0.or.stored_fingerprint/=payload_fingerprint)goto 930
    if(allocated(row_values))deallocate(row_values);ok=.true.;return
900 message='occupied checkpoint MPI read failed';goto 920
910 message='occupied checkpoint file read failed';goto 920
930 message='corrupt occupied checkpoint payload'
920 if(rank==0.and.opened)close(unit);call cleanup_output;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0
#endif
  contains
#ifdef USE_MPI
    subroutine agree_expected(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,lo_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo_i/=hi_i.or.value==0_int64)bad=1
    end subroutine
    subroutine agree_expected_count(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_expected_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,lo_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.transfer(lo_r,bits)/=transfer(hi_r,bits))bad=1
    end subroutine
    subroutine read_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine cleanup_output
      if(allocated(row_ids))deallocate(row_ids);if(allocated(coefficients))deallocate(coefficients)
      if(allocated(occupations))deallocate(occupations);if(allocated(eigenvalues))deallocate(eigenvalues)
      if(allocated(row_values))deallocate(row_values);global_count=0;payload_fingerprint=0_int64
    end subroutine
#endif
  end subroutine read_rt_dg_hybrid_occupied_checkpoint

#ifdef USE_MPI
  subroutine valid_stream_extent(unit,rank,count,element_bytes,comm,valid,ierr)
    integer,intent(in)::unit,rank,element_bytes,comm
    integer(int64),intent(in)::count
    logical,intent(out)::valid
    integer,intent(out)::ierr
    integer(int64)::file_size,file_position,available_bytes
    integer::inquire_status
    file_size=0_int64;file_position=1_int64;inquire_status=0
    if(rank==0)inquire(unit=unit,size=file_size,pos=file_position,iostat=inquire_status)
    call MPI_Bcast(inquire_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(file_position,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    available_bytes=max(0_int64,file_size-file_position+1_int64)
    valid=inquire_status==0.and.count>=0_int64.and.&
      count<=available_bytes/max(1_int64,int(element_bytes,int64))
  end subroutine valid_stream_extent

  subroutine valid_stream_shape(unit,rank,dims,element_bytes,comm,count,valid,ierr)
    integer,intent(in)::unit,rank,dims(:),element_bytes,comm
    integer,intent(out)::count
    logical,intent(out)::valid
    integer,intent(out)::ierr
    integer::i
    integer(int64)::wide_count
    valid=.false.;count=0;ierr=MPI_SUCCESS
    if(any(dims<0))return
    wide_count=1_int64
    do i=1,size(dims)
      if(dims(i)==0)then
        wide_count=0_int64;exit
      endif
      if(wide_count>huge(0_int64)/int(dims(i),int64))return
      wide_count=wide_count*int(dims(i),int64)
    enddo
    if(wide_count>int(huge(0),int64))return
    call valid_stream_extent(unit,rank,wide_count,element_bytes,comm,valid,ierr)
    if(valid)count=int(wide_count)
  end subroutine valid_stream_shape

  subroutine stream_write_i64_1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer(int64),allocatable,intent(in)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    integer(int64),allocatable::buffer(:);logical::extent_ok
    n=0;if(rank==owner.and.allocated(a))n=size(a);call MPI_Bcast(n,1,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.n>0)then
      if(rank==0)then;allocate(buffer(n));call MPI_Recv(buffer,n,MPI_INTEGER8,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,n,MPI_INTEGER8,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)n
      if(io_status==0.and.n>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_i64_1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer(int64),allocatable,intent(inout)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    integer(int64),allocatable::buffer(:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)n;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(n,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_extent(unit,rank,int(n,int64),8,comm,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(n));if(n>0)read(unit,iostat=io_status)a
      else;allocate(buffer(n));if(n>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(n));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.n>0)then
      if(rank==0)call MPI_Send(buffer,n,MPI_INTEGER8,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,n,MPI_INTEGER8,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_i1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer,allocatable,intent(in)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    integer,allocatable::buffer(:);logical::extent_ok
    n=0;if(rank==owner.and.allocated(a))n=size(a);call MPI_Bcast(n,1,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.n>0)then
      if(rank==0)then;allocate(buffer(n));call MPI_Recv(buffer,n,MPI_INTEGER,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,n,MPI_INTEGER,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)n
      if(io_status==0.and.n>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_i1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer,allocatable,intent(inout)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    integer,allocatable::buffer(:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)n;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(n,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_extent(unit,rank,int(n,int64),4,comm,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(n));if(n>0)read(unit,iostat=io_status)a
      else;allocate(buffer(n));if(n>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(n));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.n>0)then
      if(rank==0)call MPI_Send(buffer,n,MPI_INTEGER,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,n,MPI_INTEGER,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_i2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer,allocatable,intent(in)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),status(MPI_STATUS_SIZE)
    integer,allocatable::buffer(:,:);logical::extent_ok
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a);call MPI_Bcast(dims,2,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2)));call MPI_Recv(buffer,product(dims),MPI_INTEGER,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_INTEGER,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_i2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;integer,allocatable,intent(inout)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),count,status(MPI_STATUS_SIZE)
    integer,allocatable::buffer(:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,4,comm,count,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_INTEGER,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_INTEGER,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_r1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(in)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:);logical::extent_ok
    n=0;if(rank==owner.and.allocated(a))n=size(a);call MPI_Bcast(n,1,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.n>0)then
      if(rank==0)then;allocate(buffer(n));call MPI_Recv(buffer,n,MPI_DOUBLE_PRECISION,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,n,MPI_DOUBLE_PRECISION,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)n
      if(io_status==0.and.n>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_r1(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(inout)::a(:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::n,status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)n;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(n,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_extent(unit,rank,int(n,int64),8,comm,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(n));if(n>0)read(unit,iostat=io_status)a
      else;allocate(buffer(n));if(n>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(n));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.n>0)then
      if(rank==0)call MPI_Send(buffer,n,MPI_DOUBLE_PRECISION,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,n,MPI_DOUBLE_PRECISION,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_r2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(in)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:,:);logical::extent_ok
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a);call MPI_Bcast(dims,2,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2)))
        call MPI_Recv(buffer,product(dims),MPI_DOUBLE_PRECISION,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_PRECISION,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_r2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(inout)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),count,status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,8,comm,count,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_PRECISION,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_PRECISION,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_z2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(in)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:);logical::extent_ok
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a);call MPI_Bcast(dims,2,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2)));call MPI_Recv(buffer,product(dims),MPI_DOUBLE_COMPLEX,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_COMPLEX,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine
  subroutine stream_read_z2(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(inout)::a(:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(2),count,status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,16,comm,count,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_COMPLEX,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_COMPLEX,0,29031,comm,status,ierr)
    endif
  end subroutine
  subroutine stream_write_z3(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(in)::a(:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(3),status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:)
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a);call MPI_Bcast(dims,3,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2),dims(3)))
        call MPI_Recv(buffer,product(dims),MPI_DOUBLE_COMPLEX,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_COMPLEX,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine stream_write_z3
  subroutine stream_read_z3(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(inout)::a(:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(3),count,status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:)
    logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims;call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,16,comm,count,extent_ok,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2),dims(3)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2),dims(3)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2),dims(3)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_COMPLEX,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_COMPLEX,0,29031,comm,status,ierr)
    endif
  end subroutine stream_read_z3
  subroutine stream_write_r3(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(in)::a(:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(3),status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:,:,:)
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a)
    call MPI_Bcast(dims,3,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2),dims(3)))
        call MPI_Recv(buffer,product(dims),MPI_DOUBLE_PRECISION,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_PRECISION,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine stream_write_r3
  subroutine stream_read_r3(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;real(real64),allocatable,intent(inout)::a(:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(3),count,status(MPI_STATUS_SIZE)
    real(real64),allocatable::buffer(:,:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,8,comm,count,extent_ok,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2),dims(3)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2),dims(3)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2),dims(3)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_PRECISION,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_PRECISION,0,29031,comm,status,ierr)
    endif
  end subroutine stream_read_r3
  subroutine stream_write_z4(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(in)::a(:,:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(4),status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:,:)
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a)
    call MPI_Bcast(dims,4,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2),dims(3),dims(4)))
        call MPI_Recv(buffer,product(dims),MPI_DOUBLE_COMPLEX,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_COMPLEX,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine stream_write_z4
  subroutine stream_read_z4(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(inout)::a(:,:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(4),count,status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,4,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,16,comm,count,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2),dims(3),dims(4)));if(count>0)read(unit,iostat=io_status)a
      else;allocate(buffer(dims(1),dims(2),dims(3),dims(4)));if(count>0)read(unit,iostat=io_status)buffer;endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2),dims(3),dims(4)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_COMPLEX,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_COMPLEX,0,29031,comm,status,ierr)
    endif
  end subroutine stream_read_z4
  subroutine stream_write_z5(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(in)::a(:,:,:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(5),status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:,:,:)
    dims=0;if(rank==owner.and.allocated(a))dims=shape(a)
    call MPI_Bcast(dims,5,MPI_INTEGER,owner,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(owner/=0.and.product(dims)>0)then
      if(rank==0)then
        allocate(buffer(dims(1),dims(2),dims(3),dims(4),dims(5)))
        call MPI_Recv(buffer,product(dims),MPI_DOUBLE_COMPLEX,owner,29032,comm,status,ierr)
      else if(rank==owner)then;call MPI_Send(a,product(dims),MPI_DOUBLE_COMPLEX,0,29032,comm,ierr);endif
      if((rank==0.or.rank==owner).and.ierr/=MPI_SUCCESS)return
    endif
    if(rank==0)then
      write(unit,iostat=io_status)dims
      if(io_status==0.and.product(dims)>0)then
        if(owner==0)then;write(unit,iostat=io_status)a;else;write(unit,iostat=io_status)buffer;endif
      endif
    endif
    call sync_io(io_status,comm,ierr)
  end subroutine stream_write_z5
  subroutine stream_read_z5(comm,owner,unit,rank,a,io_status,ierr)
    integer,intent(in)::comm,owner,unit,rank;complex(real64),allocatable,intent(inout)::a(:,:,:,:,:)
    integer,intent(inout)::io_status;integer,intent(out)::ierr;integer::dims(5),count,status(MPI_STATUS_SIZE)
    complex(real64),allocatable::buffer(:,:,:,:,:);logical::extent_ok
    if(rank==0)read(unit,iostat=io_status)dims
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    call MPI_Bcast(dims,5,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call valid_stream_shape(unit,rank,dims,16,comm,count,extent_ok,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.extent_ok)then;io_status=1;call sync_io(io_status,comm,ierr);return;endif
    if(rank==0)then
      if(owner==0)then;allocate(a(dims(1),dims(2),dims(3),dims(4),dims(5)));if(count>0)read(unit,iostat=io_status)a
      else
        allocate(buffer(dims(1),dims(2),dims(3),dims(4),dims(5)))
        if(count>0)read(unit,iostat=io_status)buffer
      endif
    else if(rank==owner)then;allocate(a(dims(1),dims(2),dims(3),dims(4),dims(5)));endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)return
    if(owner/=0.and.count>0)then
      if(rank==0)call MPI_Send(buffer,count,MPI_DOUBLE_COMPLEX,owner,29031,comm,ierr)
      if(rank==owner)call MPI_Recv(a,count,MPI_DOUBLE_COMPLEX,0,29031,comm,status,ierr)
    endif
  end subroutine stream_read_z5
#endif

  subroutine validate_ground_state_payload(payload,bad)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer,intent(out)::bad
    integer::nrow,npoint,nrtrow,r,i
    real(real64)::hamiltonian_scale,electron_difference,electron_scale
    bad=0;nrow=0;npoint=0
    if(allocated(payload%row_ids))nrow=size(payload%row_ids)
    if(allocated(payload%grid_ids))npoint=size(payload%grid_ids)
    if(.not.payload%valid.or..not.payload%final_refresh_complete.or..not.payload%analysis_complete.or.&
      payload%global_count<1.or.payload%global_grid_count<1.or.payload%noccupied<1.or.&
      payload%noccupied>payload%global_count.or.&
      payload%operation_count<1.or.payload%nonidentity_operation_count<0.or.&
      payload%nonidentity_operation_count>=payload%operation_count)bad=1
    if(any([payload%catalog_fingerprint,payload%state_fingerprint,payload%metric_fingerprint,&
      payload%operator_structure_fingerprint,payload%operator_value_fingerprint,payload%kinetic_fingerprint,&
      payload%nonlocal_fingerprint,payload%local_fingerprint,payload%sipg_fingerprint,payload%basis_fingerprint,&
      payload%face_fingerprint,payload%dc_seed_fingerprint,payload%continuation_fingerprint,payload%scope_fingerprint,&
      payload%analysis_fingerprint,payload%selection_fingerprint,payload%pseudopotential_fingerprint,&
      payload%energy_fingerprint,payload%position_convention_fingerprint]==0_int64))bad=1
    if(bad/=0)return
    if(.not.allocated(payload%row_ids).or..not.allocated(payload%metric_rows).or.&
      .not.allocated(payload%kinetic_rows).or..not.allocated(payload%nonlocal_rows).or.&
      .not.allocated(payload%local_rows).or..not.allocated(payload%sipg_rows).or.&
      .not.allocated(payload%hamiltonian_rows).or..not.allocated(payload%coefficients).or.&
      .not.allocated(payload%occupations).or..not.allocated(payload%eigenvalues).or.&
      .not.allocated(payload%position_rows).or..not.allocated(payload%symmetry_representation))then;bad=1;return;endif
    if(.not.allocated(payload%metric_row_offsets).or..not.allocated(payload%metric_column_ids).or.&
      .not.allocated(payload%operator_row_offsets).or..not.allocated(payload%operator_column_ids))then;bad=1;return;endif
    if(size(payload%metric_row_offsets)/=nrow+1.or.size(payload%operator_row_offsets)/=nrow+1)then;bad=1;return;endif
    if(payload%metric_row_offsets(1)/=1.or.payload%operator_row_offsets(1)/=1.or.&
      payload%metric_row_offsets(nrow+1)-1/=size(payload%metric_column_ids).or.&
      payload%operator_row_offsets(nrow+1)-1/=size(payload%operator_column_ids).or.&
      any(payload%metric_row_offsets<1).or.any(payload%metric_row_offsets>size(payload%metric_column_ids)+1).or.&
      any(payload%operator_row_offsets<1).or.any(payload%operator_row_offsets>size(payload%operator_column_ids)+1).or.&
      any(payload%metric_row_offsets(2:)<payload%metric_row_offsets(:nrow)).or.&
      any(payload%operator_row_offsets(2:)<payload%operator_row_offsets(:nrow)))then;bad=1;return;endif
    if(bad==0.and.nrow>0)then
      if(any(payload%row_ids<1_int64).or.any(payload%row_ids>int(payload%global_count,int64)).or.&
        any(payload%metric_column_ids<1).or.any(payload%metric_column_ids>payload%global_count).or.&
        any(payload%operator_column_ids<1).or.any(payload%operator_column_ids>payload%global_count))bad=1
    endif
    if(any(shape(payload%metric_rows)/=[nrow,payload%global_count]).or.&
      any(shape(payload%kinetic_rows)/=shape(payload%metric_rows)).or.&
      any(shape(payload%nonlocal_rows)/=shape(payload%metric_rows)).or.&
      any(shape(payload%local_rows)/=shape(payload%metric_rows)).or.&
      any(shape(payload%sipg_rows)/=shape(payload%metric_rows)).or.&
      any(shape(payload%hamiltonian_rows)/=shape(payload%metric_rows)).or.&
      any(shape(payload%position_rows)/=[3,nrow,payload%global_count]).or.&
      any(shape(payload%symmetry_representation)/=[payload%global_count,payload%global_count,payload%operation_count]).or.&
      any(shape(payload%coefficients)/=[nrow,payload%noccupied]).or.&
      size(payload%occupations)/=payload%noccupied.or.size(payload%eigenvalues)/=payload%noccupied)then;bad=1;return;endif
    if(bad==0)then
      if(.not.finite_matrix(payload%metric_rows).or..not.finite_matrix(payload%kinetic_rows).or.&
        .not.finite_matrix(payload%nonlocal_rows).or..not.finite_matrix(payload%local_rows).or.&
        .not.finite_matrix(payload%sipg_rows).or..not.finite_matrix(payload%hamiltonian_rows).or.&
        .not.finite_matrix(payload%coefficients))bad=1
      if(any(.not.ieee_is_finite(real(payload%position_rows))).or.any(.not.ieee_is_finite(aimag(payload%position_rows))).or.&
        any(.not.ieee_is_finite(real(payload%symmetry_representation))).or.&
        any(.not.ieee_is_finite(aimag(payload%symmetry_representation))))bad=1
    endif
    if(bad/=0)return
    if(maximum_complex_component(payload%kinetic_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%nonlocal_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%local_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%sipg_rows)>huge(0d0)/8d0)then;bad=1;return;endif
    hamiltonian_scale=max(1d0,maximum_complex_component(payload%hamiltonian_rows),&
      maximum_complex_component(payload%kinetic_rows),maximum_complex_component(payload%nonlocal_rows),&
      maximum_complex_component(payload%local_rows),maximum_complex_component(payload%sipg_rows))
    if(maximum_complex_component(payload%hamiltonian_rows-&
      (payload%kinetic_rows+payload%nonlocal_rows+payload%local_rows+payload%sipg_rows))>&
      64d0*epsilon(1d0)*hamiltonian_scale)then;bad=1;return;endif
    if(.not.allocated(payload%grid_ids).or..not.allocated(payload%grid_weights).or.&
      .not.allocated(payload%partition_ids).or..not.allocated(payload%basis_values).or.&
      .not.allocated(payload%density))then;bad=1;return;endif
    if(size(payload%grid_weights)/=npoint.or.size(payload%partition_ids)/=npoint.or.&
      any(shape(payload%basis_values)/=[payload%global_count,npoint]).or.size(payload%density)/=npoint)then;bad=1;return;endif
    if(npoint>0)then
      if(any(payload%grid_ids<1_int64).or.any(payload%grid_ids>int(payload%global_grid_count,int64)))bad=1
    endif
    if(.not.allocated(payload%requested_ids).or..not.allocated(payload%effective_ids).or.&
      .not.allocated(payload%added_ids).or..not.allocated(payload%closure_parent).or.&
      .not.allocated(payload%closure_reason).or..not.allocated(payload%closure_action).or.&
      .not.allocated(payload%scope_selectors).or..not.allocated(payload%xc_types).or.&
      .not.allocated(payload%continuation_receipt).or..not.allocated(payload%pseudopotential_receipt).or.&
      .not.allocated(payload%energy_receipt))then;bad=1;return;endif
    if(.not.allocated(payload%face_ids).or..not.allocated(payload%face_point_ids).or.&
      .not.allocated(payload%face_metadata).or..not.allocated(payload%face_offsets).or.&
      .not.allocated(payload%face_weight_offsets).or..not.allocated(payload%face_basis_offsets).or.&
      .not.allocated(payload%face_value_offsets).or..not.allocated(payload%face_observable_offsets).or.&
      .not.allocated(payload%face_basis_ids).or.&
      .not.allocated(payload%face_normals).or..not.allocated(payload%face_weights).or.&
      .not.allocated(payload%face_values).or..not.allocated(payload%interface_observables).or.&
      .not.allocated(payload%nonlocal_ids).or..not.allocated(payload%nonlocal_owner).or.&
      .not.allocated(payload%nonlocal_values))then;bad=1;return;endif
    if(bad==0)then
      if(size(payload%face_metadata,2)/=size(payload%face_ids).or.&
        size(payload%face_normals,2)/=size(payload%face_ids).or.size(payload%face_offsets)/=size(payload%face_ids)+1.or.&
        size(payload%face_weight_offsets)/=size(payload%face_ids)+1.or.&
        size(payload%face_basis_offsets)/=size(payload%face_ids)+1.or.&
        size(payload%face_value_offsets)/=size(payload%face_ids)+1.or.&
        size(payload%face_observable_offsets)/=size(payload%face_ids)+1.or.&
        size(payload%interface_observables,1)/=3.or.&
        size(payload%nonlocal_owner)/=size(payload%nonlocal_ids).or.&
        size(payload%nonlocal_values,2)/=size(payload%nonlocal_ids).or.&
        size(payload%added_ids)/=size(payload%closure_parent).or.size(payload%added_ids)/=size(payload%closure_reason).or.&
        size(payload%added_ids)/=size(payload%closure_action).or.size(payload%requested_ids)<1.or.&
        size(payload%effective_ids)<size(payload%requested_ids).or.&
        size(payload%added_ids)/=size(payload%effective_ids)-size(payload%requested_ids))bad=1
      if(bad/=0)return
      if(any(payload%requested_ids<=0).or.any(payload%effective_ids<=0).or.any(payload%added_ids<=0).or.&
          any(payload%closure_parent<=0).or.any(payload%closure_reason<=0).or.any(payload%closure_action<=0))then
        bad=1;return
      endif
      do i=1,size(payload%requested_ids)
        if(count(payload%requested_ids==payload%requested_ids(i))/=1.or.&
            count(payload%effective_ids==payload%requested_ids(i))/=1)bad=1
      enddo
      do i=1,size(payload%effective_ids)
        if(count(payload%effective_ids==payload%effective_ids(i))/=1.or.&
            count(payload%requested_ids==payload%effective_ids(i))+&
              count(payload%added_ids==payload%effective_ids(i))/=1)bad=1
      enddo
      do i=1,size(payload%added_ids)
        if(count(payload%added_ids==payload%added_ids(i))/=1.or.&
            count(payload%requested_ids==payload%added_ids(i))/=0.or.&
            count(payload%effective_ids==payload%added_ids(i))/=1.or.&
            count(payload%effective_ids==payload%closure_parent(i))/=1)bad=1
      enddo
      if(bad/=0)return
      if(size(payload%face_ids)>0)then
        if(any(payload%face_ids<=0_int64))bad=1
      endif
      if(size(payload%nonlocal_ids)>0)then
        if(any(payload%nonlocal_ids<=0_int64))bad=1
      endif
      if(payload%face_offsets(1)/=1.or.payload%face_weight_offsets(1)/=1.or.&
        payload%face_basis_offsets(1)/=1.or.payload%face_value_offsets(1)/=1.or.&
        payload%face_observable_offsets(1)/=1.or.&
        payload%face_offsets(size(payload%face_offsets))-1/=size(payload%face_point_ids).or.&
        payload%face_weight_offsets(size(payload%face_weight_offsets))-1/=size(payload%face_weights).or.&
        payload%face_basis_offsets(size(payload%face_basis_offsets))-1/=size(payload%face_basis_ids).or.&
        payload%face_value_offsets(size(payload%face_value_offsets))-1/=size(payload%face_values,2).or.&
        payload%face_observable_offsets(size(payload%face_observable_offsets))-1/=&
          size(payload%interface_observables,2).or.&
        any(payload%face_offsets<1).or.any(payload%face_offsets>size(payload%face_point_ids)+1).or.&
        any(payload%face_weight_offsets<1).or.any(payload%face_weight_offsets>size(payload%face_weights)+1).or.&
        any(payload%face_basis_offsets<1).or.any(payload%face_basis_offsets>size(payload%face_basis_ids)+1).or.&
        any(payload%face_value_offsets<1).or.any(payload%face_value_offsets>size(payload%face_values,2)+1).or.&
        any(payload%face_observable_offsets<1).or.&
          any(payload%face_observable_offsets>size(payload%interface_observables,2)+1).or.&
        any(payload%face_offsets(2:)<payload%face_offsets(:size(payload%face_offsets)-1)).or.&
        any(payload%face_weight_offsets(2:)<&
          payload%face_weight_offsets(:size(payload%face_weight_offsets)-1)).or.&
        any(payload%face_basis_offsets(2:)<&
          payload%face_basis_offsets(:size(payload%face_basis_offsets)-1)).or.&
        any(payload%face_value_offsets(2:)<&
          payload%face_value_offsets(:size(payload%face_value_offsets)-1)).or.&
        any(payload%face_observable_offsets(2:)<&
          payload%face_observable_offsets(:size(payload%face_observable_offsets)-1)))bad=1
      if(size(payload%face_point_ids)>0)then
        if(any(payload%face_point_ids<=0_int64))bad=1
      endif
      if(size(payload%face_basis_ids)>0)then
        if(any(payload%face_basis_ids<1).or.any(payload%face_basis_ids>payload%global_count))bad=1
      endif
      if(.not.finite_matrix(payload%basis_values).or..not.finite_matrix(payload%face_values).or.&
        .not.finite_matrix(payload%nonlocal_values).or..not.finite_matrix(payload%interface_observables).or.&
        any(.not.ieee_is_finite(payload%grid_weights)).or.any(.not.ieee_is_finite(payload%face_normals)).or.&
        any(.not.ieee_is_finite(payload%face_weights)).or.any(.not.ieee_is_finite(payload%density)).or.&
        any(.not.ieee_is_finite(payload%occupations)).or.any(.not.ieee_is_finite(payload%eigenvalues)).or.&
        any(.not.ieee_is_finite(payload%continuation_receipt)).or.&
        any(.not.ieee_is_finite(payload%pseudopotential_receipt)).or.any(.not.ieee_is_finite(payload%energy_receipt)))bad=1
    endif
    if(.not.payload%construction_catalog%valid.or.&
        payload%construction_catalog%global_count/=payload%global_count.or.&
        payload%construction_catalog%catalog_fingerprint/=payload%catalog_fingerprint)bad=1
    if(.not.allocated(payload%construction_catalog%ids).or.&
      .not.allocated(payload%construction_catalog%generations).or.&
      .not.allocated(payload%construction_catalog%ordering).or.&
      .not.allocated(payload%construction_catalog%ownership))then;bad=1;return;endif
    if(size(payload%construction_catalog%ids)/=payload%global_count.or.&
      size(payload%construction_catalog%generations)/=payload%global_count.or.&
      size(payload%construction_catalog%ordering)/=payload%global_count.or.&
      size(payload%construction_catalog%ownership)/=payload%global_count)then;bad=1;return;endif
    if(any(payload%construction_catalog%ids<=0_int64).or.any(payload%construction_catalog%generations<0).or.&
      any(payload%construction_catalog%ownership<=0))bad=1
    do i=1,payload%global_count
      if(count(payload%construction_catalog%ids==payload%construction_catalog%ids(i))/=1.or.&
        count(payload%construction_catalog%ordering==i)/=1)bad=1
    enddo
    if(any([payload%construction_catalog%ids_fingerprint,&
      payload%construction_catalog%generation_fingerprint,payload%construction_catalog%ordering_fingerprint,&
      payload%construction_catalog%ownership_fingerprint,payload%construction_catalog%provenance_fingerprint,&
      payload%construction_catalog%catalog_fingerprint]==0_int64))bad=1

    r=payload%certified_basis%certified_count
    if(.not.payload%certified_basis%valid.or..not.payload%certified_basis%localization_converged.or.&
      payload%certified_basis%localization_symmetry_constrained.or.&
      payload%certified_basis%construction_count/=payload%global_count.or.r<1.or.r>payload%global_count.or.&
      r<payload%noccupied.or.payload%certified_basis%occupied_count/=payload%noccupied.or.&
      payload%certified_basis%localization_iterations<0)bad=1
    if(bad/=0)return
    if(.not.allocated(payload%certified_basis%construction_row_ids).or.&
      .not.allocated(payload%certified_basis%transformation_row_ids).or.&
      .not.allocated(payload%certified_basis%c_cert).or..not.allocated(payload%certified_basis%u_rt).or.&
      .not.allocated(payload%certified_basis%b_rt).or.&
      .not.allocated(payload%certified_basis%initial_occupied_amplitudes).or.&
      .not.allocated(payload%certified_basis%certified_eigenvalues).or.&
      .not.allocated(payload%certified_basis%occupations).or..not.allocated(payload%certified_basis%centers).or.&
      .not.allocated(payload%certified_basis%spreads_before).or.&
      .not.allocated(payload%certified_basis%spreads_after))then;bad=1;return;endif
    if(size(payload%certified_basis%construction_row_ids)/=nrow.or.&
      any(shape(payload%certified_basis%c_cert)/=[nrow,r]).or.&
      any(shape(payload%certified_basis%b_rt)/=[nrow,r]).or.&
      size(payload%certified_basis%transformation_row_ids)/=size(payload%certified_basis%u_rt,1).or.&
      size(payload%certified_basis%u_rt,2)/=r.or.&
      any(shape(payload%certified_basis%initial_occupied_amplitudes)/=[r,payload%noccupied]).or.&
      size(payload%certified_basis%certified_eigenvalues)/=r.or.&
      size(payload%certified_basis%occupations)/=payload%noccupied.or.&
      any(shape(payload%certified_basis%centers)/=[3,r]).or.&
      size(payload%certified_basis%spreads_before)/=r.or.size(payload%certified_basis%spreads_after)/=r)then;bad=1;return;endif
    if(any(payload%certified_basis%construction_row_ids/=payload%row_ids))then;bad=1;return;endif
    if(size(payload%certified_basis%transformation_row_ids)>0)then
      if(any(payload%certified_basis%transformation_row_ids<1_int64).or.&
        any(payload%certified_basis%transformation_row_ids>int(r,int64)))then;bad=1;return;endif
    endif
    if(bad==0)then
      if(.not.finite_matrix(payload%certified_basis%c_cert).or..not.finite_matrix(payload%certified_basis%u_rt).or.&
        .not.finite_matrix(payload%certified_basis%b_rt).or.&
        .not.finite_matrix(payload%certified_basis%initial_occupied_amplitudes).or.&
        any(.not.ieee_is_finite(payload%certified_basis%certified_eigenvalues)).or.&
        any(.not.ieee_is_finite(payload%certified_basis%occupations)).or.&
        any(.not.ieee_is_finite(payload%certified_basis%centers)).or.&
        any(.not.ieee_is_finite(payload%certified_basis%spreads_before)).or.&
        any(.not.ieee_is_finite(payload%certified_basis%spreads_after)))bad=1
    endif
    if(bad/=0)return
    if(r>1)then
      if(any(payload%certified_basis%certified_eigenvalues(2:)<&
        payload%certified_basis%certified_eigenvalues(:r-1)))then;bad=1;return;endif
    endif
    if(any([payload%certified_basis%c_cert_fingerprint,payload%certified_basis%u_rt_fingerprint,&
      payload%certified_basis%b_rt_fingerprint,payload%certified_basis%initial_state_fingerprint,&
      payload%certified_basis%transformation_fingerprint,payload%certified_basis%operator_fingerprint,&
      payload%certified_basis%fingerprint]==0_int64))bad=1
    if(any(.not.ieee_is_finite([payload%certified_basis%spread_before_total,&
      payload%certified_basis%spread_after_total,payload%certified_basis%spread_improvement,&
      payload%certified_basis%transform_unitarity_defect,payload%certified_basis%certified_metric_defect,&
      payload%certified_basis%rt_metric_defect,payload%certified_basis%embedding_defect,&
      payload%certified_basis%projector_invariance_defect,payload%certified_basis%target_symmetry_defect_before,&
      payload%certified_basis%target_symmetry_defect_after,payload%certified_basis%energy_symmetry_defect_before,&
      payload%certified_basis%energy_symmetry_defect_after,payload%certified_basis%symmetry_defect_invariance,&
      payload%certified_basis%scalar_covariance_defect,payload%certified_basis%vector_covariance_defect,&
      payload%certified_basis%tensor_covariance_defect])))then;bad=1;return;endif
    if(min(payload%certified_basis%spread_before_total,payload%certified_basis%spread_after_total,&
      payload%certified_basis%transform_unitarity_defect,&
      payload%certified_basis%certified_metric_defect,payload%certified_basis%rt_metric_defect,&
      payload%certified_basis%embedding_defect,payload%certified_basis%projector_invariance_defect,&
      payload%certified_basis%target_symmetry_defect_before,payload%certified_basis%target_symmetry_defect_after,&
      payload%certified_basis%energy_symmetry_defect_before,payload%certified_basis%energy_symmetry_defect_after,&
      payload%certified_basis%symmetry_defect_invariance,payload%certified_basis%scalar_covariance_defect,&
      payload%certified_basis%vector_covariance_defect,payload%certified_basis%tensor_covariance_defect)<0d0)bad=1

    if(.not.payload%electron_count%valid.or.payload%electron_count%fingerprint==0_int64)then;bad=1;return;endif
    if(any(.not.ieee_is_finite([payload%electron_count%expected_count,payload%electron_count%actual_count,&
      payload%electron_count%tolerance,payload%electron_count%defect,payload%electron_count%omitted_tail,&
      payload%electron_count%chemical_potential])))then;bad=1;return;endif
    electron_difference=abs(payload%electron_count%expected_count-payload%electron_count%actual_count)
    electron_scale=max(1d0,abs(payload%electron_count%expected_count),abs(payload%electron_count%actual_count))
    if(payload%electron_count%expected_count<0d0.or.&
      payload%electron_count%actual_count<0d0.or.payload%electron_count%tolerance<=0d0.or.&
      payload%electron_count%defect<0d0.or.payload%electron_count%omitted_tail<0d0.or.&
      payload%electron_count%defect>payload%electron_count%tolerance.or.&
      payload%electron_count%omitted_tail>payload%electron_count%tolerance.or.&
      abs(payload%electron_count%defect-electron_difference)>64d0*epsilon(1d0)*electron_scale.or.&
      abs(sum(payload%certified_basis%occupations)-payload%electron_count%actual_count)>&
        payload%electron_count%tolerance.or.&
      electron_difference>payload%electron_count%tolerance)then;bad=1;return;endif

    nrtrow=0;if(allocated(payload%rt_space%row_ids))nrtrow=size(payload%rt_space%row_ids)
    if(.not.payload%rt_space%valid.or.payload%rt_space%rank/=r.or.&
      payload%rt_space%operation_count/=payload%operation_count.or.payload%rt_space%scalar_count<0.or.&
      payload%rt_space%vector_count<rt_dg_hybrid_vector_canonical_momentum.or.&
      payload%rt_space%tensor_count<0)bad=1
    if(bad/=0)return
    if(.not.allocated(payload%rt_space%row_ids).or..not.allocated(payload%rt_space%row_owner_keys).or.&
      .not.allocated(payload%rt_space%grid_owner_keys).or..not.allocated(payload%rt_space%metric_rows).or.&
      .not.allocated(payload%rt_space%kinetic_rows).or..not.allocated(payload%rt_space%nonlocal_rows).or.&
      .not.allocated(payload%rt_space%local_rows).or..not.allocated(payload%rt_space%sipg_rows).or.&
      .not.allocated(payload%rt_space%hamiltonian_rows).or..not.allocated(payload%rt_space%representation).or.&
      .not.allocated(payload%rt_space%cartesian_rotations).or.&
      .not.allocated(payload%rt_space%scalar_operator_rows).or.&
      .not.allocated(payload%rt_space%vector_operator_rows).or.&
      .not.allocated(payload%rt_space%tensor_operator_rows).or..not.allocated(payload%rt_space%basis_values).or.&
      .not.allocated(payload%rt_space%density))then;bad=1;return;endif
    if(size(payload%rt_space%row_owner_keys)/=r.or.size(payload%rt_space%grid_owner_keys)/=npoint.or.&
      any(shape(payload%rt_space%metric_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%kinetic_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%nonlocal_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%local_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%sipg_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%hamiltonian_rows)/=[nrtrow,r]).or.&
      any(shape(payload%rt_space%representation)/=[r,r,payload%rt_space%operation_count]).or.&
      any(shape(payload%rt_space%cartesian_rotations)/=[3,3,payload%rt_space%operation_count]).or.&
      any(shape(payload%rt_space%scalar_operator_rows)/=[nrtrow,r,payload%rt_space%scalar_count]).or.&
      any(shape(payload%rt_space%vector_operator_rows)/=[nrtrow,r,3,payload%rt_space%vector_count]).or.&
      any(shape(payload%rt_space%tensor_operator_rows)/=[nrtrow,r,3,3,payload%rt_space%tensor_count]).or.&
      any(shape(payload%rt_space%basis_values)/=[r,npoint]).or.size(payload%rt_space%density)/=npoint)then;bad=1;return;endif
    if(nrtrow>0)then
      if(any(payload%rt_space%row_ids<1_int64).or.any(payload%rt_space%row_ids>int(r,int64)))bad=1
    endif
    if(any(payload%rt_space%row_owner_keys<=0).or.any(payload%rt_space%grid_owner_keys<=0))bad=1
    if(bad==0)then
      if(.not.finite_matrix(payload%rt_space%metric_rows).or.&
        .not.finite_matrix(payload%rt_space%kinetic_rows).or.&
        .not.finite_matrix(payload%rt_space%nonlocal_rows).or.&
        .not.finite_matrix(payload%rt_space%local_rows).or.&
        .not.finite_matrix(payload%rt_space%sipg_rows).or.&
        .not.finite_matrix(payload%rt_space%hamiltonian_rows).or.&
        any(.not.ieee_is_finite(real(payload%rt_space%representation))).or.&
        any(.not.ieee_is_finite(aimag(payload%rt_space%representation))).or.&
        any(.not.ieee_is_finite(payload%rt_space%cartesian_rotations)).or.&
        any(.not.ieee_is_finite(real(payload%rt_space%scalar_operator_rows))).or.&
        any(.not.ieee_is_finite(aimag(payload%rt_space%scalar_operator_rows))).or.&
        any(.not.ieee_is_finite(real(payload%rt_space%vector_operator_rows))).or.&
        any(.not.ieee_is_finite(aimag(payload%rt_space%vector_operator_rows))).or.&
        any(.not.ieee_is_finite(real(payload%rt_space%tensor_operator_rows))).or.&
        any(.not.ieee_is_finite(aimag(payload%rt_space%tensor_operator_rows))).or.&
        .not.finite_matrix(payload%rt_space%basis_values).or.&
        any(.not.ieee_is_finite(payload%rt_space%density)))bad=1
    endif
    if(bad/=0)return
    if(maximum_complex_component(payload%rt_space%kinetic_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%rt_space%nonlocal_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%rt_space%local_rows)>huge(0d0)/8d0.or.&
      maximum_complex_component(payload%rt_space%sipg_rows)>huge(0d0)/8d0)then;bad=1;return;endif
    hamiltonian_scale=max(1d0,maximum_complex_component(payload%rt_space%hamiltonian_rows),&
      maximum_complex_component(payload%rt_space%kinetic_rows),&
      maximum_complex_component(payload%rt_space%nonlocal_rows),&
      maximum_complex_component(payload%rt_space%local_rows),&
      maximum_complex_component(payload%rt_space%sipg_rows))
    if(maximum_complex_component(payload%rt_space%hamiltonian_rows-&
      (payload%rt_space%kinetic_rows+payload%rt_space%nonlocal_rows+payload%rt_space%local_rows+&
      payload%rt_space%sipg_rows))>64d0*epsilon(1d0)*hamiltonian_scale)then;bad=1;return;endif
    if(any([payload%rt_space%metric_fingerprint,payload%rt_space%kinetic_fingerprint,&
      payload%rt_space%nonlocal_fingerprint,payload%rt_space%local_fingerprint,&
      payload%rt_space%sipg_fingerprint,payload%rt_space%hamiltonian_fingerprint,&
      payload%rt_space%basis_fingerprint,payload%rt_space%density_fingerprint,&
      payload%rt_space%ownership_fingerprint,payload%rt_space%scalar_fingerprint,&
      payload%rt_space%vector_fingerprint,payload%rt_space%tensor_fingerprint,&
      payload%rt_space%representation_fingerprint,payload%rt_space%fingerprint]==0_int64))bad=1

    if(.not.payload%energy_window%valid.or.payload%energy_window%fingerprint==0_int64.or.&
      payload%energy_window%construction_rank/=payload%global_count.or.&
      payload%energy_window%solved_rank/=payload%energy_window%construction_rank.or.&
      payload%energy_window%occupied_rank/=payload%noccupied.or.&
      payload%energy_window%requested_rank<payload%noccupied.or.&
      payload%energy_window%requested_rank>r.or.payload%energy_window%certified_rank/=r.or.&
      payload%energy_window%boundary_cluster_rank/=r.or.&
      payload%energy_window%extension_states/=r-payload%energy_window%requested_rank.or.&
      (payload%energy_window%proof_state_present.neqv.r<payload%energy_window%solved_rank))then
      bad=1;return
    endif
    if(any(.not.ieee_is_finite([payload%energy_window%window_size,payload%energy_window%e_homo,&
      payload%energy_window%requested_cutoff,payload%energy_window%certified_cutoff,&
      payload%energy_window%extension_energy,payload%energy_window%proof_energy])))then;bad=1;return;endif
    if(maxval(abs([payload%energy_window%window_size,payload%energy_window%e_homo,&
      payload%energy_window%requested_cutoff,payload%energy_window%certified_cutoff,&
      payload%energy_window%extension_energy,payload%energy_window%proof_energy]))>huge(0d0)/4d0)then
      bad=1;return
    endif
    if(payload%energy_window%e_homo/=payload%certified_basis%certified_eigenvalues(payload%noccupied).or.&
      payload%energy_window%certified_cutoff/=payload%certified_basis%certified_eigenvalues(r).or.&
      payload%energy_window%extension_energy/=&
        max(0d0,payload%energy_window%certified_cutoff-payload%energy_window%requested_cutoff))then
      bad=1;return
    endif
    select case(payload%energy_window%mode)
    case(rt_dg_hybrid_energy_window_explicit)
      if(payload%energy_window%compatibility_dynamic_rank.or.payload%energy_window%window_size<0d0.or.&
        .not.payload%energy_window%proof_state_present.or.&
        payload%energy_window%requested_cutoff/=&
          payload%energy_window%e_homo+payload%energy_window%window_size)then
        bad=1;return
      endif
    case(rt_dg_hybrid_energy_window_legacy_dynamic)
      if(.not.payload%energy_window%compatibility_dynamic_rank.or.payload%energy_window%window_size/=-1d0.or.&
        payload%energy_window%requested_cutoff/=&
          payload%certified_basis%certified_eigenvalues(payload%energy_window%requested_rank))then
        bad=1;return
      endif
    case default
      bad=1;return
    end select
    if(payload%energy_window%proof_state_present)then
      if(payload%energy_window%proof_status<=0.or.&
        payload%energy_window%proof_energy<=max(payload%energy_window%certified_cutoff,&
          payload%energy_window%requested_cutoff))then;bad=1;return;endif
    else if(payload%energy_window%proof_status/=0.or.payload%energy_window%proof_energy/=0d0)then
      bad=1
    endif
    if(.not.payload%symmetry_receipt%valid.or.payload%symmetry_receipt%worst_operation<1.or.&
      payload%symmetry_receipt%worst_operation>payload%operation_count.or.&
      payload%symmetry_receipt%fingerprint==0_int64)then;bad=1;return;endif
    if(any(.not.ieee_is_finite([payload%symmetry_receipt%occupied_subspace_defect,&
      payload%symmetry_receipt%occupied_projector_defect,payload%symmetry_receipt%target_subspace_defect,&
      payload%symmetry_receipt%target_energy_defect,payload%symmetry_receipt%density_defect,&
      payload%symmetry_receipt%scalar_covariance_defect,payload%symmetry_receipt%vector_covariance_defect,&
      payload%symmetry_receipt%tensor_covariance_defect,payload%symmetry_receipt%final_basis_defect,&
      payload%symmetry_receipt%worst_operation_defect,payload%symmetry_receipt%maximum_physical_defect])))then
      bad=1;return
    endif
    if(min(payload%symmetry_receipt%occupied_subspace_defect,payload%symmetry_receipt%occupied_projector_defect,&
      payload%symmetry_receipt%target_subspace_defect,payload%symmetry_receipt%target_energy_defect,&
      payload%symmetry_receipt%density_defect,payload%symmetry_receipt%scalar_covariance_defect,&
      payload%symmetry_receipt%vector_covariance_defect,payload%symmetry_receipt%tensor_covariance_defect,&
      payload%symmetry_receipt%final_basis_defect,payload%symmetry_receipt%worst_operation_defect,&
      payload%symmetry_receipt%maximum_physical_defect)<0d0)bad=1
    if(.not.payload%handoff_receipts%valid.or.payload%handoff_receipts%fingerprint==0_int64.or.&
      payload%handoff_receipts%position_fingerprint/=payload%position_convention_fingerprint.or.&
      payload%handoff_receipts%nonlocal_fingerprint/=payload%nonlocal_fingerprint.or.&
      payload%handoff_receipts%face_fingerprint/=payload%face_fingerprint.or.&
      payload%handoff_receipts%pseudopotential_fingerprint/=payload%pseudopotential_fingerprint.or.&
      payload%handoff_receipts%transformation_fingerprint/=&
        payload%certified_basis%transformation_fingerprint)bad=1
  end subroutine validate_ground_state_payload

  subroutine hash_ground_state_payload(payload,hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(out)::hash
    integer::i,first,last
    integer(int64)::record
    hash=0_int64
    do i=1,size(payload%row_ids)
      record=1001_int64;call add_i64(record,payload%row_ids(i))
      first=payload%metric_row_offsets(i);last=payload%metric_row_offsets(i+1)-1
      call add_i1(record,payload%metric_column_ids(first:last))
      first=payload%operator_row_offsets(i);last=payload%operator_row_offsets(i+1)-1
      call add_i1(record,payload%operator_column_ids(first:last))
      call add_z1(record,payload%metric_rows(i,:));call add_z1(record,payload%kinetic_rows(i,:))
      call add_z1(record,payload%nonlocal_rows(i,:));call add_z1(record,payload%local_rows(i,:))
      call add_z1(record,payload%sipg_rows(i,:));call add_z1(record,payload%hamiltonian_rows(i,:))
      call add_z1(record,payload%coefficients(i,:));call add_z2(record,payload%position_rows(:,i,:))
      call add_i64(record,payload%certified_basis%construction_row_ids(i))
      call add_z1(record,payload%certified_basis%c_cert(i,:))
      call add_z1(record,payload%certified_basis%b_rt(i,:))
      hash=ieor(hash,record)
    enddo
    do i=1,size(payload%grid_ids)
      record=1002_int64;call add_i64(record,payload%grid_ids(i))
      call add_i64(record,int(payload%partition_ids(i),int64));call add_r(record,payload%grid_weights(i))
      call add_r(record,payload%density(i));call add_z1(record,payload%basis_values(:,i))
      call add_i64(record,int(payload%rt_space%grid_owner_keys(i),int64))
      call add_r(record,payload%rt_space%density(i));call add_z1(record,payload%rt_space%basis_values(:,i))
      hash=ieor(hash,record)
    enddo
    do i=1,size(payload%face_ids)
      record=1003_int64;call add_i64(record,payload%face_ids(i))
      call add_i1(record,payload%face_metadata(:,i))
      first=payload%face_offsets(i);last=payload%face_offsets(i+1)-1
      call add_i64_1(record,payload%face_point_ids(first:last))
      first=payload%face_weight_offsets(i);last=payload%face_weight_offsets(i+1)-1
      call add_r1(record,payload%face_weights(first:last))
      first=payload%face_basis_offsets(i);last=payload%face_basis_offsets(i+1)-1
      call add_i1(record,payload%face_basis_ids(first:last))
      first=payload%face_value_offsets(i);last=payload%face_value_offsets(i+1)-1
      call add_z2(record,payload%face_values(:,first:last))
      first=payload%face_observable_offsets(i);last=payload%face_observable_offsets(i+1)-1
      call add_z2(record,payload%interface_observables(:,first:last))
      call add_r1(record,payload%face_normals(:,i))
      hash=ieor(hash,record)
    enddo
    do i=1,size(payload%nonlocal_ids)
      record=1004_int64;call add_i64(record,payload%nonlocal_ids(i))
      call add_i64(record,int(payload%nonlocal_owner(i),int64));call add_z1(record,payload%nonlocal_values(:,i))
      hash=ieor(hash,record)
    enddo
    do i=1,size(payload%rt_space%row_ids)
      record=1005_int64;call add_i64(record,payload%rt_space%row_ids(i))
      call add_z1(record,payload%rt_space%metric_rows(i,:));call add_z1(record,payload%rt_space%kinetic_rows(i,:))
      call add_z1(record,payload%rt_space%nonlocal_rows(i,:));call add_z1(record,payload%rt_space%local_rows(i,:))
      call add_z1(record,payload%rt_space%sipg_rows(i,:));call add_z1(record,payload%rt_space%hamiltonian_rows(i,:))
      call add_z2(record,payload%rt_space%scalar_operator_rows(i,:,:))
      call add_z3(record,payload%rt_space%vector_operator_rows(i,:,:,:))
      call add_z4(record,payload%rt_space%tensor_operator_rows(i,:,:,:,:))
      hash=ieor(hash,record)
    enddo
    do i=1,size(payload%certified_basis%transformation_row_ids)
      record=1006_int64;call add_i64(record,payload%certified_basis%transformation_row_ids(i))
      call add_z1(record,payload%certified_basis%u_rt(i,:));hash=ieor(hash,record)
    enddo
  contains
    subroutine add_i64(seed,value)
      integer(int64),intent(inout)::seed;integer(int64),intent(in)::value
      seed=ground_state_mix_hash(seed,value)
    end subroutine
    subroutine add_i1(seed,values)
      integer(int64),intent(inout)::seed;integer,intent(in)::values(:);integer::j
      seed=ground_state_mix_hash(seed,int(size(values),int64));do j=1,size(values);call add_i64(seed,int(values(j),int64));enddo
    end subroutine
    subroutine add_i64_1(seed,values)
      integer(int64),intent(inout)::seed;integer(int64),intent(in)::values(:);integer::j
      seed=ground_state_mix_hash(seed,int(size(values),int64));do j=1,size(values);call add_i64(seed,values(j));enddo
    end subroutine
    subroutine add_r(seed,value)
      integer(int64),intent(inout)::seed;real(real64),intent(in)::value;integer(int64)::bits
      bits=transfer(value,bits);seed=ground_state_mix_hash(seed,bits)
    end subroutine
    subroutine add_r1(seed,values)
      integer(int64),intent(inout)::seed;real(real64),intent(in)::values(:);integer::j
      seed=ground_state_mix_hash(seed,int(size(values),int64));do j=1,size(values);call add_r(seed,values(j));enddo
    end subroutine
    subroutine add_z(seed,value)
      integer(int64),intent(inout)::seed;complex(real64),intent(in)::value
      call add_r(seed,real(value,real64));call add_r(seed,aimag(value))
    end subroutine
    subroutine add_z1(seed,values)
      integer(int64),intent(inout)::seed;complex(real64),intent(in)::values(:);integer::j
      seed=ground_state_mix_hash(seed,int(size(values),int64));do j=1,size(values);call add_z(seed,values(j));enddo
    end subroutine
    subroutine add_z2(seed,values)
      integer(int64),intent(inout)::seed;complex(real64),intent(in)::values(:,:);integer::j,k
      seed=ground_state_mix_hash(seed,int(size(values,1),int64));seed=ground_state_mix_hash(seed,int(size(values,2),int64))
      do k=1,size(values,2);do j=1,size(values,1);call add_z(seed,values(j,k));enddo;enddo
    end subroutine
    subroutine add_z3(seed,values)
      integer(int64),intent(inout)::seed;complex(real64),intent(in)::values(:,:,:);integer::j,k,l
      seed=ground_state_mix_hash(seed,int(size(values,1),int64));seed=ground_state_mix_hash(seed,int(size(values,2),int64))
      seed=ground_state_mix_hash(seed,int(size(values,3),int64))
      do l=1,size(values,3);do k=1,size(values,2);do j=1,size(values,1)
        call add_z(seed,values(j,k,l))
      enddo;enddo;enddo
    end subroutine
    subroutine add_z4(seed,values)
      integer(int64),intent(inout)::seed;complex(real64),intent(in)::values(:,:,:,:);integer::j,k,l,m
      seed=ground_state_mix_hash(seed,int(size(values,1),int64));seed=ground_state_mix_hash(seed,int(size(values,2),int64))
      seed=ground_state_mix_hash(seed,int(size(values,3),int64));seed=ground_state_mix_hash(seed,int(size(values,4),int64))
      do m=1,size(values,4);do l=1,size(values,3);do k=1,size(values,2);do j=1,size(values,1)
        call add_z(seed,values(j,k,l,m))
      enddo;enddo;enddo;enddo
    end subroutine
  end subroutine hash_ground_state_payload

  subroutine hash_z_cube(hash,a)
    integer(int64),intent(inout)::hash;complex(real64),intent(in)::a(:,:,:)
    integer::i,j,k;integer(int64)::bits
    hash=ground_state_mix_hash(hash,int(size(a,1),int64))
    hash=ground_state_mix_hash(hash,int(size(a,2),int64))
    hash=ground_state_mix_hash(hash,int(size(a,3),int64))
    do k=1,size(a,3);do j=1,size(a,2);do i=1,size(a,1)
      bits=transfer(real(a(i,j,k)),bits);hash=ground_state_mix_hash(hash,bits)
      bits=transfer(aimag(a(i,j,k)),bits);hash=ground_state_mix_hash(hash,bits)
    enddo;enddo;enddo
  end subroutine hash_z_cube

  subroutine hash_r_cube(hash,a)
    integer(int64),intent(inout)::hash;real(real64),intent(in)::a(:,:,:)
    integer::i,j,k;integer(int64)::bits
    hash=ground_state_mix_hash(hash,int(size(a,1),int64))
    hash=ground_state_mix_hash(hash,int(size(a,2),int64))
    hash=ground_state_mix_hash(hash,int(size(a,3),int64))
    do k=1,size(a,3);do j=1,size(a,2);do i=1,size(a,1)
      bits=transfer(a(i,j,k),bits);hash=ground_state_mix_hash(hash,bits)
    enddo;enddo;enddo
  end subroutine hash_r_cube

  subroutine hash_ground_state_common(payload,hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(out)::hash
    logical::l(ground_state_logical_count);integer::h(ground_state_integer_count),i
    integer(int64)::f(ground_state_fingerprint_count),bits
    real(real64)::r(ground_state_real_count)
    call pack_ground_state_header(payload,l,h,f,r);f(20)=0_int64
    hash=3001_int64
    do i=1,size(l);hash=ground_state_mix_hash(hash,merge(1_int64,0_int64,l(i)));enddo
    do i=1,size(h);hash=ground_state_mix_hash(hash,int(h(i),int64));enddo
    do i=1,size(f);hash=ground_state_mix_hash(hash,f(i));enddo
    do i=1,size(r);bits=transfer(r(i),bits);hash=ground_state_mix_hash(hash,bits);enddo
    call hash_i_array(hash,payload%requested_ids);call hash_i_array(hash,payload%effective_ids)
    call hash_i_array(hash,payload%added_ids);call hash_i_array(hash,payload%closure_parent)
    call hash_i_array(hash,payload%closure_reason);call hash_i_array(hash,payload%closure_action)
    call hash_i_array(hash,payload%scope_selectors);call hash_i_array(hash,payload%xc_types)
    call hash_r_array(hash,payload%occupations);call hash_r_array(hash,payload%eigenvalues)
    call hash_r_array(hash,payload%continuation_receipt);call hash_r_array(hash,payload%pseudopotential_receipt)
    call hash_r_array(hash,payload%energy_receipt)
    call hash_i64_array(hash,payload%construction_catalog%ids)
    call hash_i_array(hash,payload%construction_catalog%generations)
    call hash_i_array(hash,payload%construction_catalog%ordering)
    call hash_i_array(hash,payload%construction_catalog%ownership)
    call hash_r_array(hash,payload%certified_basis%certified_eigenvalues)
    call hash_r_array(hash,payload%certified_basis%occupations)
    call hash_z_matrix(hash,payload%certified_basis%initial_occupied_amplitudes)
    call hash_r_matrix(hash,payload%certified_basis%centers)
    call hash_r_array(hash,payload%certified_basis%spreads_before)
    call hash_r_array(hash,payload%certified_basis%spreads_after)
    call hash_i_array(hash,payload%rt_space%row_owner_keys)
    call hash_z_cube(hash,payload%rt_space%representation)
    call hash_r_cube(hash,payload%rt_space%cartesian_rotations)
  end subroutine hash_ground_state_common

  subroutine hash_i64_array(hash,a)
    integer(int64),intent(inout)::hash;integer(int64),allocatable,intent(in)::a(:);integer::i
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a),int64))
    do i=1,size(a);hash=ground_state_mix_hash(hash,a(i));enddo
  end subroutine
  subroutine hash_i_array(hash,a)
    integer(int64),intent(inout)::hash;integer,allocatable,intent(in)::a(:);integer::i
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a),int64))
    do i=1,size(a);hash=ground_state_mix_hash(hash,int(a(i),int64));enddo
  end subroutine
  subroutine hash_i_matrix(hash,a)
    integer(int64),intent(inout)::hash;integer,allocatable,intent(in)::a(:,:);integer::i,j
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a,1),int64));hash=ground_state_mix_hash(hash,int(size(a,2),int64))
    do j=1,size(a,2);do i=1,size(a,1)
      hash=ground_state_mix_hash(hash,int(a(i,j),int64))
    enddo;enddo
  end subroutine
  subroutine hash_r_array(hash,a)
    integer(int64),intent(inout)::hash;real(real64),allocatable,intent(in)::a(:);integer::i;integer(int64)::bits
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a),int64))
    do i=1,size(a);bits=transfer(a(i),bits);hash=ground_state_mix_hash(hash,bits);enddo
  end subroutine
  subroutine hash_r_matrix(hash,a)
    integer(int64),intent(inout)::hash;real(real64),allocatable,intent(in)::a(:,:);integer::i,j;integer(int64)::bits
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a,1),int64));hash=ground_state_mix_hash(hash,int(size(a,2),int64))
    do j=1,size(a,2);do i=1,size(a,1)
      bits=transfer(a(i,j),bits);hash=ground_state_mix_hash(hash,bits)
    enddo;enddo
  end subroutine
  subroutine hash_z_matrix(hash,a)
    integer(int64),intent(inout)::hash;complex(real64),allocatable,intent(in)::a(:,:);integer::i,j;integer(int64)::bits
    if(.not.allocated(a))then;hash=ground_state_mix_hash(hash,-1_int64);return;endif
    hash=ground_state_mix_hash(hash,int(size(a,1),int64));hash=ground_state_mix_hash(hash,int(size(a,2),int64))
    do j=1,size(a,2);do i=1,size(a,1)
      bits=transfer(real(a(i,j)),bits);hash=ground_state_mix_hash(hash,bits)
      bits=transfer(aimag(a(i,j)),bits);hash=ground_state_mix_hash(hash,bits)
    enddo;enddo
  end subroutine

  subroutine validate_path(path,probe,comm,ierr)
    character(*),intent(in)::path;character(16),intent(out)::probe;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    integer::i,status
    integer(int64)::local_hash,minimum_hash,maximum_hash
    local_hash=int(len_trim(path),int64)
    do i=1,len_trim(path);local_hash=ieor(ishftc(local_hash,7),int(iachar(path(i:i)),int64));enddo
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.minimum_hash/=maximum_hash)ierr=MPI_ERR_OTHER
    probe=''
#else
    probe='';ierr=1
#endif
  end subroutine validate_path
  subroutine sync_io(io_status,comm,ierr)
    integer,intent(inout)::io_status;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
#else
    ierr=1
#endif
  end subroutine sync_io
  subroutine atomic_rename(old_path,new_path,status)
    character(*),intent(in)::old_path,new_path;integer,intent(out)::status
    character(c_char),allocatable::old_c(:),new_c(:)
    integer::i
    allocate(old_c(len_trim(old_path)+1),new_c(len_trim(new_path)+1))
    do i=1,len_trim(old_path);old_c(i)=old_path(i:i);enddo;old_c(size(old_c))=c_null_char
    do i=1,len_trim(new_path);new_c(i)=new_path(i:i);enddo;new_c(size(new_c))=c_null_char
    status=int(c_rename(old_c,new_c))
  end subroutine atomic_rename
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
  real(real64) function maximum_complex_component(values)
    complex(real64),intent(in)::values(:,:)
    if(size(values)==0)then
      maximum_complex_component=0d0
    else
      maximum_complex_component=max(maxval(abs(real(values))),maxval(abs(aimag(values))))
    endif
  end function maximum_complex_component
  logical function strictly_increasing(values)
    integer,intent(in)::values(:);integer::i
    strictly_increasing=.true.
    do i=2,size(values);if(values(i)<=values(i-1))then;strictly_increasing=.false.;return;endif;enddo
  end function strictly_increasing
  logical function valid_packet_activity(packet_ids,active_rows)
    integer,intent(in)::packet_ids(:);logical,intent(in)::active_rows(:)
    integer,allocatable::states(:);integer::i,status
    valid_packet_activity=.false.
    if(size(packet_ids)/=size(active_rows).or.any(packet_ids<1).or.any(packet_ids>size(packet_ids)))return
    allocate(states(size(packet_ids)),stat=status);if(status/=0)return;states=-1
    do i=1,size(packet_ids)
      if(states(packet_ids(i))<0)then
        states(packet_ids(i))=merge(1,0,active_rows(i))
      else if(states(packet_ids(i))/=merge(1,0,active_rows(i)))then
        deallocate(states);return
      endif
    enddo
    deallocate(states);valid_packet_activity=.true.
  end function valid_packet_activity
  pure integer(int64) function mix_hash(seed,value)
    integer(int64),intent(in)::seed,value
    mix_hash=ieor(ishftc(seed,9),value)
  end function mix_hash
  pure integer(int64) function ground_state_mix_hash(seed,value)
    integer(int64),intent(in)::seed,value
    ground_state_mix_hash=multiply_ground_state_prime(ieor(seed,value))
  end function ground_state_mix_hash
  pure integer(int64) function multiply_ground_state_prime(value)
    integer(int64),intent(in)::value
    multiply_ground_state_prime=shiftl(value,40)
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,shiftl(value,8))
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,shiftl(value,7))
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,shiftl(value,5))
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,shiftl(value,4))
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,shiftl(value,1))
    multiply_ground_state_prime=add_modulo_64(multiply_ground_state_prime,value)
  end function multiply_ground_state_prime
  pure integer(int64) function add_modulo_64(left,right)
    integer(int64),intent(in)::left,right
    integer(int64)::x,y,carry
    x=left;y=right
    do while(y/=0_int64)
      carry=iand(x,y)
      x=ieor(x,y)
      y=shiftl(carry,1)
    enddo
    add_modulo_64=x
  end function add_modulo_64
  function int64_string(value) result(text)
    integer(int64),intent(in)::value;character(32)::text
    write(text,'(i0)')value
  end function int64_string
end module rt_dg_hybrid_checkpoint
