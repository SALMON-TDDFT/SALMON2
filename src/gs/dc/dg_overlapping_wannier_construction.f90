#include "config.h"
module dg_overlapping_wannier_construction
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_overlapping_wannier_construction
    integer::candidate_rank=0,target_rank=0,retained_rank=0,generation=0
    integer(int64)::transform_fingerprint=0_int64
    integer,allocatable::center_owner_rank(:),center_owner_fragment(:)
    integer(int64),allocatable::physical_grid_ids(:),center_box_point_ids(:)
    complex(real64),allocatable::value(:,:),gradient(:,:,:),transform(:,:)
    complex(real64),allocatable::symmetry_representation(:,:,:)
    real(real64)::occupied_inclusion_residual=huge(1d0)
    real(real64)::symmetry_closure_residual=huge(1d0)
    real(real64)::boundary_value_max=huge(1d0),boundary_gradient_max=huge(1d0)
    real(real64)::metric_minimum_eigenvalue=0d0,metric_condition_number=huge(1d0)
  end type
  public::construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction
  public::verify_dg_overlapping_wannier_periodic_closure
contains
  subroutine verify_dg_overlapping_wannier_periodic_closure(comm,box_ids,symmetry_target_box_ids,&
      values,gradients,symmetry_representation,gradient_transform,expected_box_count,tolerance,&
      residual,fingerprint,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::box_ids(:),symmetry_target_box_ids(:,:),expected_box_count
    complex(real64),intent(in)::values(:,:),gradients(:,:,:),symmetry_representation(:,:,:)
    real(real64),intent(in)::gradient_transform(:,:,:),tolerance
    real(real64),intent(out)::residual
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::all_values(:,:),all_gradients(:,:,:),local_values(:,:),local_gradients(:,:,:)
    complex(real64),allocatable::mapped_value(:),mapped_gradient(:,:)
    integer,allocatable::owners(:)
    integer::nwann,nsym,nbox,p,j,isym,target,ierr,bad,global_bad,rank
    integer(int64)::local_hash,bits,payload_count64
    ok=.false.;message='';residual=huge(1d0);fingerprint=0_int64
    call MPI_Comm_rank(comm,rank,ierr)
    nwann=size(values,1);nsym=size(symmetry_target_box_ids,2)
    bad=0
    if(expected_box_count<1_int64.or.expected_box_count>10000000_int64.or.nwann<1.or.nsym<1)bad=1
    if(size(values,2)/=size(box_ids).or.any(shape(gradients)/=[3,nwann,size(box_ids)]).or.&
       size(symmetry_target_box_ids,1)/=size(box_ids).or.&
       any(shape(symmetry_representation)/=[nwann,nwann,nsym]).or.&
       any(shape(gradient_transform)/=[3,3,nsym]).or.tolerance<=0d0)bad=1
    if(any(box_ids<1_int64).or.any(box_ids>expected_box_count))bad=1
    if(nwann>0.and.expected_box_count<=huge(1_int64)/int(nwann,int64))then
      payload_count64=int(nwann,int64)*expected_box_count
      if(payload_count64>int(huge(nbox),int64)/3_int64.or.&
         payload_count64>268435456_int64)bad=1
    else
      payload_count64=0_int64;bad=1
    endif
    if(.not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
       .not.all(ieee_is_finite(real(gradients))).or..not.all(ieee_is_finite(aimag(gradients))).or.&
       .not.all(ieee_is_finite(gradient_transform)).or..not.ieee_is_finite(tolerance))bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid authoritative periodic closure payload';return;endif
    nbox=int(expected_box_count)
    allocate(local_values(nwann,nbox),all_values(nwann,nbox),local_gradients(3,nwann,nbox),&
      all_gradients(3,nwann,nbox),owners(nbox))
    local_values=(0d0,0d0);local_gradients=(0d0,0d0);owners=0
    do p=1,size(box_ids)
      local_values(:,int(box_ids(p)))=values(:,p)
      local_gradients(:,:,int(box_ids(p)))=gradients(:,:,p)
      owners(int(box_ids(p)))=owners(int(box_ids(p)))+1
    enddo
    call MPI_Allreduce(local_values,all_values,nwann*nbox,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_gradients,all_gradients,3*nwann*nbox,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,owners,nbox,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(owners/=1))then;message='periodic closure box ownership is incomplete';return;endif
    allocate(mapped_value(nwann),mapped_gradient(3,nwann));residual=0d0;local_hash=0_int64
    do isym=1,nsym;do p=1,size(box_ids)
      target=int(symmetry_target_box_ids(p,isym))
      if(target<1.or.target>nbox)then;bad=1;cycle;endif
      mapped_value=matmul(transpose(symmetry_representation(:,:,isym)),all_values(:,target))
      mapped_gradient=matmul(gradient_transform(:,:,isym),&
        matmul(all_gradients(:,:,target),symmetry_representation(:,:,isym)))
      residual=max(residual,maxval(abs(values(:,p)-mapped_value)),&
        maxval(abs(gradients(:,:,p)-mapped_gradient)))
      do j=1,nwann
        bits=transfer(real(mapped_value(j),real64),bits)
        local_hash=ieor(local_hash,ishftc(bits,mod(j+7*isym+int(modulo(box_ids(p),63_int64)),63)))
      enddo
      local_hash=ieor(local_hash,ishftc(symmetry_target_box_ids(p,isym),&
        mod(11*isym+int(modulo(box_ids(p),53_int64)),63)))
    enddo;enddo
    if(rank==0)then
      do isym=1,nsym;do j=1,3;do p=1,3
        bits=transfer(gradient_transform(j,p,isym),bits)
        local_hash=ieor(local_hash,ishftc(bits,mod(13*j+17*p+19*isym,63)))
      enddo;enddo;enddo
    endif
    call MPI_Allreduce(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(fingerprint,int(z'510E527FADE682D1',int64))
    if(fingerprint==0_int64)fingerprint=1_int64
    if(bad/=0.or.residual>tolerance)then
      message='authoritative periodic value/gradient closure failed';return
    endif
    ok=.true.
#else
    residual=huge(1d0);fingerprint=0_int64;ok=.false.
    message='authoritative periodic closure requires MPI'
#endif
  end subroutine

  subroutine construct_dg_overlapping_wannier_basis(comm,ncandidate,ntarget,noccupied,physical_ids,&
      core_fragment,weights,localization_coordinate,boundary_mask,candidate_value,candidate_gradient,&
      occupied_coefficients,expected_core_count,generation,boundary_value_tolerance,&
      boundary_gradient_tolerance,rank_tolerance,result,ok,message,core_mask,&
      box_point_ids,symmetry_target_box_ids,expected_box_count,symmetry_tolerance,&
      periodic_localization_phase)
    integer,intent(in)::comm,ncandidate,ntarget,noccupied,generation
    integer(int64),intent(in)::physical_ids(:),expected_core_count
    integer,intent(in)::core_fragment(:)
    real(real64),intent(in)::weights(:),localization_coordinate(:)
    logical,intent(in)::boundary_mask(:)
    complex(real64),intent(in)::candidate_value(:,:),candidate_gradient(:,:,:)
    complex(real64),intent(in)::occupied_coefficients(:,:)
    real(real64),intent(in)::boundary_value_tolerance,boundary_gradient_tolerance,rank_tolerance
    type(s_dg_overlapping_wannier_construction),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::core_mask(:)
    integer(int64),intent(in),optional::box_point_ids(:),symmetry_target_box_ids(:,:)
    integer(int64),intent(in),optional::expected_box_count
    real(real64),intent(in),optional::symmetry_tolerance
    complex(real64),intent(in),optional::periodic_localization_phase(:,:)
#ifdef USE_MPI
    complex(real64),allocatable::s_local(:,:),s(:,:),l_local(:,:),localizer(:,:),x(:,:),&
      occ_gram(:,:),occ_vectors(:,:),a_occ(:,:),overlap_occ_x(:,:),residual(:,:),&
      residual_gram(:,:),residual_vectors(:,:),q_comp(:,:),comp_localizer(:,:),&
      comp_vectors(:,:),transform(:,:),final_metric(:,:),metric_vectors(:,:),&
      occupied_reference(:,:),symmetry_product(:,:),orthogonal_transform(:,:),mapped_transform(:,:),&
      candidate_symmetry(:,:,:),all_candidate_flat(:),all_candidate(:,:),spatial_overlap(:,:),&
      link_local(:,:,:),link_global(:,:,:),symmetrized_localizer(:,:),all_phase_flat(:),&
      canonical_phase(:,:),all_wannier(:,:)
    complex(real64),allocatable::symmetry_mapped_wannier(:)
    real(real64),allocatable::spectrum(:),occ_spectrum(:),residual_spectrum(:),comp_spectrum(:),&
      metric_spectrum(:)
    logical,allocatable::integration_core(:),all_core_unsorted(:),canonical_core(:)
    real(real64),allocatable::all_box_weights_unsorted(:),all_box_weights(:)
    integer(int64),allocatable::all_box_ids(:),all_target_ids(:,:),canonical_target_ids(:,:),&
      sorted_target_ids(:)
    integer,allocatable::box_counts(:),box_displs(:),value_counts(:),value_displs(:),all_fragments(:),&
      canonical_fragments(:),canonical_ranks(:),&
      phase_counts(:),phase_displs(:)
    integer::nlocal,i,j,p,ierr,nproc,rank,local_bad,global_bad,&
      ncomp,nneed,start_index,matrix_count,scalar_min(4),scalar_max(4),scalar_local(4),&
      symmetry_present,symmetry_present_min,symmetry_present_max,nsym,isym,jsym,ksym,&
      nsym_min,nsym_max,core_mask_present,core_mask_present_min,core_mask_present_max,&
      phase_present,phase_present_min,phase_present_max,fragment_min,fragment_max
    integer::total_box,source,target,axis,source_axis,conjugate_flag,center_index
    real(real64)::largest,local_boundary_value,local_boundary_gradient
    real(real64)::tolerance_local(3),tolerance_min(3),tolerance_max(3)
    real(real64)::active_symmetry_tolerance,symmetry_tolerance_min,symmetry_tolerance_max,&
      symmetry_defect,product_residual
    complex(real64)::phase_ratio,trial_ratio
    real(real64),parameter::periodic_real_weight(3)=[sqrt(2d0),sqrt(3d0),sqrt(5d0)]
    real(real64),parameter::periodic_imag_weight(3)=[sqrt(7d0),sqrt(11d0),sqrt(13d0)]
    integer(int64)::local_fingerprint,global_fingerprint,quantized_density,matrix_count64,&
      expected_min,expected_max,total_box64,value_total64,running64
    integer(int64)::expected_box_min,expected_box_max
    integer(int64)::local_core_count,global_core_count
    logical::matched_product,phase_covariant,source_axis_used(3)

    ok=.false.;message='';nlocal=size(physical_ids)
    call MPI_Comm_size(comm,nproc,ierr)
    symmetry_present=merge(1,0,present(symmetry_target_box_ids))
    call MPI_Allreduce(symmetry_present,symmetry_present_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(symmetry_present,symmetry_present_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(symmetry_present_min/=symmetry_present_max)then
      message='inconsistent periodic-box symmetry contract across ranks';return
    endif
    core_mask_present=merge(1,0,present(core_mask))
    call MPI_Allreduce(core_mask_present,core_mask_present_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(core_mask_present,core_mask_present_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(core_mask_present_min/=core_mask_present_max)then
      message='inconsistent periodic-box core-mask contract across ranks';return
    endif
    phase_present=merge(1,0,present(periodic_localization_phase))
    call MPI_Allreduce(phase_present,phase_present_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(phase_present,phase_present_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(phase_present_min/=phase_present_max)then
      message='inconsistent periodic localization phase contract across ranks';return
    endif
    nsym=0
    if(present(symmetry_target_box_ids))nsym=size(symmetry_target_box_ids,2)
    call MPI_Allreduce(nsym,nsym_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nsym,nsym_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(nsym_min/=nsym_max)then
      message='inconsistent periodic-box symmetry count across ranks';return
    endif
    active_symmetry_tolerance=rank_tolerance
    if(present(symmetry_tolerance))active_symmetry_tolerance=symmetry_tolerance
    call MPI_Allreduce(active_symmetry_tolerance,symmetry_tolerance_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(active_symmetry_tolerance,symmetry_tolerance_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scalar_local=[ncandidate,ntarget,noccupied,generation]
    call MPI_Allreduce(scalar_local,scalar_min,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(scalar_local,scalar_max,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    tolerance_local=[boundary_value_tolerance,boundary_gradient_tolerance,rank_tolerance]
    call MPI_Allreduce(tolerance_local,tolerance_min,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(tolerance_local,tolerance_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(any(scalar_min/=scalar_max).or.any(tolerance_min/=tolerance_max).or.&
        expected_min/=expected_max.or.symmetry_tolerance_min/=symmetry_tolerance_max)then
      message='inconsistent overlapping-Wannier construction contract across ranks';return
    endif
    local_bad=0
    if(ncandidate<=0.or.noccupied<=0.or.ntarget<noccupied.or.ntarget>ncandidate) local_bad=1
    if(generation<=0.or.expected_core_count<=0_int64.or.rank_tolerance<=0d0) local_bad=1
    if(boundary_value_tolerance<=0d0.or.boundary_gradient_tolerance<=0d0)local_bad=1
    if(size(core_fragment)/=nlocal.or.size(weights)/=nlocal.or.&
        size(localization_coordinate)/=nlocal.or.size(boundary_mask)/=nlocal)local_bad=1
    if(size(candidate_value,1)/=ncandidate.or.size(candidate_value,2)/=nlocal)local_bad=1
    if(size(candidate_gradient,1)/=3.or.size(candidate_gradient,2)/=ncandidate.or.&
        size(candidate_gradient,3)/=nlocal)local_bad=1
    if(size(occupied_coefficients,1)/=ncandidate.or.size(occupied_coefficients,2)/=noccupied)local_bad=1
    if(present(core_mask))then
      if(size(core_mask)/=nlocal)local_bad=1
    endif
    if(present(symmetry_target_box_ids))then
      if(.not.present(box_point_ids).or..not.present(expected_box_count).or..not.present(core_mask))then
        local_bad=1
      else
        if(size(box_point_ids)/=nlocal.or.expected_box_count<=0_int64)local_bad=1
      endif
      if(size(symmetry_target_box_ids,1)/=nlocal.or.size(symmetry_target_box_ids,2)<=0)local_bad=1
      if(.not.present(periodic_localization_phase))then
        local_bad=1
      else
        if(size(periodic_localization_phase,1)/=3.or.size(periodic_localization_phase,2)/=nlocal)local_bad=1
        if(.not.all(ieee_is_finite(real(periodic_localization_phase))).or.&
            .not.all(ieee_is_finite(aimag(periodic_localization_phase))))local_bad=1
        if(any(abs(abs(periodic_localization_phase)-1d0)>active_symmetry_tolerance))local_bad=1
      endif
    endif
    if(active_symmetry_tolerance<=0d0.or..not.ieee_is_finite(active_symmetry_tolerance))local_bad=1
    if(any(physical_ids<=0_int64).or.any(core_fragment<=0).or.any(weights<=0d0))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(localization_coordinate)))local_bad=1
    if(.not.all(ieee_is_finite(real(candidate_value))).or.&
        .not.all(ieee_is_finite(aimag(candidate_value))))local_bad=1
    if(.not.all(ieee_is_finite(real(candidate_gradient))).or.&
        .not.all(ieee_is_finite(aimag(candidate_gradient))))local_bad=1
    if(.not.all(ieee_is_finite(real(occupied_coefficients))).or.&
        .not.all(ieee_is_finite(aimag(occupied_coefficients))))local_bad=1
    if(ncandidate>0)then
      if(int(ncandidate,int64)>huge(1_int64)/int(ncandidate,int64))then
        local_bad=1;matrix_count64=0_int64
      else
        matrix_count64=int(ncandidate,int64)*int(ncandidate,int64)
        if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
        if(nsym>0.and.matrix_count64>int(huge(matrix_count),int64)/3_int64)local_bad=1
      endif
    else
      matrix_count64=0_int64
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid overlapping-Wannier construction metadata';return;endif
    if(nsym>0)then
      call MPI_Allreduce(expected_box_count,expected_box_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      call MPI_Allreduce(expected_box_count,expected_box_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(expected_box_min/=expected_box_max)then
        message='inconsistent periodic-box point count across ranks';return
      endif
    endif
    allocate(integration_core(nlocal));integration_core=.true.
    if(present(core_mask))integration_core=core_mask
    local_core_count=int(count(integration_core),int64)
    call MPI_Allreduce(local_core_count,global_core_count,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if(global_core_count/=expected_core_count)then
      message='periodic-box core point count does not match construction contract';return
    endif
    fragment_min=huge(fragment_min);fragment_max=0
    if(nlocal>0)then
      fragment_min=minval(core_fragment);fragment_max=maxval(core_fragment)
    endif
    if(fragment_min/=fragment_max)then
      message='construction call must contain exactly one fragment box';return
    endif
    allocate(occupied_reference(ncandidate,noccupied));occupied_reference=occupied_coefficients
    call MPI_Bcast(occupied_reference,ncandidate*noccupied,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    local_bad=merge(0,1,maxval(abs(occupied_reference-occupied_coefficients))<=rank_tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='inconsistent occupied candidate subspace across ranks';return;endif
    if(nsym>0)then
      allocate(box_counts(nproc),box_displs(nproc),value_counts(nproc),value_displs(nproc),&
        phase_counts(nproc),phase_displs(nproc))
      call MPI_Allgather(nlocal,1,MPI_INTEGER,box_counts,1,MPI_INTEGER,comm,ierr)
      local_bad=0;total_box64=0_int64
      do i=1,nproc
        if(box_counts(i)<0.or.total_box64>int(huge(total_box),int64)-int(box_counts(i),int64))local_bad=1
        if(local_bad==0)total_box64=total_box64+int(box_counts(i),int64)
      enddo
      if(total_box64>0_int64.and.int(ncandidate,int64)>int(huge(total_box),int64)/total_box64)local_bad=1
      if(total_box64>0_int64.and.int(nsym,int64)>int(huge(total_box),int64)/total_box64)local_bad=1
      if(total_box64>int(huge(total_box),int64)/3_int64)local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0)then;message='periodic-box collective extent overflow';return;endif
      total_box=int(total_box64);value_total64=int(ncandidate,int64)*total_box64
      box_displs(1)=0;value_displs(1)=0;phase_displs(1)=0;running64=0_int64
      do i=1,nproc
        if(int(ncandidate,int64)>0_int64.and.int(box_counts(i),int64)>&
            int(huge(total_box),int64)/int(ncandidate,int64))then
          message='periodic-box collective count overflow';return
        endif
        value_counts(i)=ncandidate*box_counts(i)
        phase_counts(i)=3*box_counts(i)
        if(i>1)then
          box_displs(i)=int(running64)
          value_displs(i)=int(int(ncandidate,int64)*running64)
          phase_displs(i)=int(3_int64*running64)
        endif
        running64=running64+int(box_counts(i),int64)
      enddo
      allocate(all_box_ids(total_box),all_box_weights_unsorted(total_box),all_fragments(total_box),&
        all_core_unsorted(total_box))
      allocate(all_candidate_flat(int(value_total64)),all_target_ids(total_box,nsym),&
        all_phase_flat(3*total_box))
      call MPI_Allgatherv(box_point_ids,nlocal,MPI_INTEGER8,all_box_ids,box_counts,box_displs,MPI_INTEGER8,comm,ierr)
      call MPI_Allgatherv(weights,nlocal,MPI_DOUBLE_PRECISION,all_box_weights_unsorted,box_counts,&
        box_displs,MPI_DOUBLE_PRECISION,comm,ierr)
      call MPI_Allgatherv(core_fragment,nlocal,MPI_INTEGER,all_fragments,box_counts,box_displs,MPI_INTEGER,comm,ierr)
      call MPI_Allgatherv(integration_core,nlocal,MPI_LOGICAL,all_core_unsorted,box_counts,box_displs,&
        MPI_LOGICAL,comm,ierr)
      call MPI_Allgatherv(candidate_value,ncandidate*nlocal,MPI_DOUBLE_COMPLEX,all_candidate_flat,&
        value_counts,value_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
      call MPI_Allgatherv(periodic_localization_phase,3*nlocal,MPI_DOUBLE_COMPLEX,all_phase_flat,&
        phase_counts,phase_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
      do isym=1,nsym
        call MPI_Allgatherv(symmetry_target_box_ids(:,isym),nlocal,MPI_INTEGER8,all_target_ids(:,isym),&
          box_counts,box_displs,MPI_INTEGER8,comm,ierr)
      enddo
      if(int(total_box,int64)/=expected_box_count)then
        message='periodic-box construction call must contain complete fragment boxes';return
      endif
      allocate(all_candidate(ncandidate,total_box),all_box_weights(total_box),canonical_core(total_box),&
        canonical_phase(3,total_box),canonical_fragments(total_box),canonical_ranks(total_box))
      allocate(canonical_target_ids(total_box,nsym))
      do source=1,total_box
        if(all_box_ids(source)<1_int64.or.all_box_ids(source)>int(total_box,int64))then
          message='invalid periodic-box point id';return
        endif
        target=int(all_box_ids(source))
        all_candidate(:,target)=all_candidate_flat((source-1)*ncandidate+1:source*ncandidate)
        all_box_weights(target)=all_box_weights_unsorted(source)
        canonical_core(target)=all_core_unsorted(source)
        canonical_fragments(target)=all_fragments(source)
        canonical_ranks(target)=count(box_displs<=source-1)-1
        canonical_phase(:,target)=all_phase_flat((source-1)*3+1:source*3)
        canonical_target_ids(target,:)=all_target_ids(source,:)
      enddo
      allocate(sorted_target_ids,source=all_box_ids);call sort_ids(sorted_target_ids)
      if(any(sorted_target_ids/=[(int(source,int64),source=1,total_box)]))then
        message='duplicate or missing periodic-box point';return
      endif
      deallocate(sorted_target_ids)
      do isym=1,nsym
        allocate(sorted_target_ids,source=canonical_target_ids(:,isym));call sort_ids(sorted_target_ids)
        if(any(sorted_target_ids/=[(int(source,int64),source=1,total_box)]))then
          message='periodic-box symmetry point map is not a permutation';return
        endif
        deallocate(sorted_target_ids)
        do source=1,total_box
          target=int(canonical_target_ids(source,isym))
          if(canonical_core(target).neqv.canonical_core(source))then
            message='periodic-box symmetry does not preserve the core-buffer partition';return
          endif
        enddo
        source_axis_used=.false.
        do axis=1,3
          phase_covariant=.false.
          do source_axis=1,3
            if(source_axis_used(source_axis))cycle
            do conjugate_flag=0,1
              target=int(canonical_target_ids(1,isym))
              if(conjugate_flag==0)then
                phase_ratio=canonical_phase(axis,target)/canonical_phase(source_axis,1)
              else
                phase_ratio=canonical_phase(axis,target)/conjg(canonical_phase(source_axis,1))
              endif
              product_residual=0d0
              do source=1,total_box
                target=int(canonical_target_ids(source,isym))
                if(conjugate_flag==0)then
                  trial_ratio=canonical_phase(axis,target)/canonical_phase(source_axis,source)
                else
                  trial_ratio=canonical_phase(axis,target)/conjg(canonical_phase(source_axis,source))
                endif
                product_residual=max(product_residual,abs(trial_ratio-phase_ratio))
              enddo
              if(product_residual<=active_symmetry_tolerance)then
                phase_covariant=.true.;source_axis_used(source_axis)=.true.;exit
              endif
            enddo
            if(phase_covariant)exit
          enddo
          if(.not.phase_covariant)then
            message='periodic localization phases are not covariant under box symmetry';return
          endif
        enddo
      enddo
    endif

    allocate(s_local(ncandidate,ncandidate),l_local(ncandidate,ncandidate))
    s_local=(0d0,0d0);l_local=(0d0,0d0)
    if(nsym>0)then
      allocate(link_local(3,ncandidate,ncandidate),link_global(3,ncandidate,ncandidate))
      link_local=(0d0,0d0)
    endif
    do p=1,nlocal
      do j=1,ncandidate;do i=1,ncandidate
        s_local(i,j)=s_local(i,j)+weights(p)*conjg(candidate_value(i,p))*candidate_value(j,p)
        if(nsym>0)then
          do axis=1,3
            link_local(axis,i,j)=link_local(axis,i,j)+weights(p)*conjg(candidate_value(i,p))*&
              periodic_localization_phase(axis,p)*candidate_value(j,p)
          enddo
        else
          l_local(i,j)=l_local(i,j)+weights(p)*localization_coordinate(p)*&
            conjg(candidate_value(i,p))*candidate_value(j,p)
        endif
      enddo;enddo
    enddo
    allocate(s(ncandidate,ncandidate),localizer(ncandidate,ncandidate))
    matrix_count=int(matrix_count64)
    call MPI_Allreduce(s_local,s,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(nsym>0)then
      call MPI_Allreduce(link_local,link_global,3*matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      localizer=(0d0,0d0)
      do axis=1,3
        localizer=localizer+periodic_real_weight(axis)*&
          0.5d0*(link_global(axis,:,:)+conjg(transpose(link_global(axis,:,:))))+&
          periodic_imag_weight(axis)*cmplx(0d0,-0.5d0,real64)*&
          (link_global(axis,:,:)-conjg(transpose(link_global(axis,:,:))))
      enddo
    else
      call MPI_Allreduce(l_local,localizer,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif
    s=0.5d0*(s+conjg(transpose(s)));localizer=0.5d0*(localizer+conjg(transpose(localizer)))
    if(.not.finite_complex_matrix(s).or..not.finite_complex_matrix(localizer))then
      ok=.false.;message='nonfinite periodic-box overlap or localization matrix';return
    endif
    if(.not.positive_definite_above(s,rank_tolerance))then
      ok=.false.;message='candidate rank loss in overlapping-Wannier construction';return
    endif
    call hermitian_eigensystem(s,spectrum,x,ok,message)
    if(.not.ok)return
    largest=maxval(spectrum)
    if(largest<=0d0.or.count(spectrum>rank_tolerance*largest)<ncandidate)then
      ok=.false.;message='candidate rank loss in overlapping-Wannier construction';return
    endif
    do j=1,ncandidate
      x(:,j)=x(:,j)/sqrt(spectrum(j))
    enddo
    if(nsym>0)then
      allocate(candidate_symmetry(ncandidate,ncandidate,nsym),spatial_overlap(ncandidate,ncandidate))
      do isym=1,nsym
        spatial_overlap=(0d0,0d0)
        do source=1,total_box
          target=int(canonical_target_ids(source,isym))
          if(abs(all_box_weights(target)-all_box_weights(source))>&
              active_symmetry_tolerance*max(1d0,all_box_weights(source)))then
            ok=.false.;message='periodic-box symmetry does not preserve quadrature weights';return
          endif
          do j=1,ncandidate;do i=1,ncandidate
            spatial_overlap(i,j)=spatial_overlap(i,j)+all_box_weights(target)*&
              conjg(all_candidate(i,target))*all_candidate(j,source)
          enddo;enddo
        enddo
        candidate_symmetry(:,:,isym)=matmul(conjg(transpose(x)),matmul(spatial_overlap,x))
      enddo
      if(.not.all(ieee_is_finite(real(candidate_symmetry))).or.&
          .not.all(ieee_is_finite(aimag(candidate_symmetry))))then
        ok=.false.;message='nonfinite periodic-box candidate symmetry representation';return
      endif
      allocate(symmetry_product(ncandidate,ncandidate))
      symmetry_defect=0d0
      do isym=1,nsym
        symmetry_defect=max(symmetry_defect,maxval(abs(matmul(conjg(transpose(candidate_symmetry(:,:,isym))),&
          candidate_symmetry(:,:,isym))-identity_complex(ncandidate))))
      enddo
      if(symmetry_defect>active_symmetry_tolerance)then
        ok=.false.;message='periodic-box candidate symmetry representation is not unitary';return
      endif
      do isym=1,nsym;do jsym=1,nsym
        symmetry_product=matmul(candidate_symmetry(:,:,isym),candidate_symmetry(:,:,jsym))
        matched_product=.false.
        do ksym=1,nsym
          product_residual=maxval(abs(symmetry_product-candidate_symmetry(:,:,ksym)))
          do source=1,total_box
            target=int(canonical_target_ids(int(canonical_target_ids(source,jsym)),isym))
            if(target/=int(canonical_target_ids(source,ksym)))product_residual=huge(1d0)
          enddo
          if(product_residual<=active_symmetry_tolerance)matched_product=.true.
        enddo
        if(.not.matched_product)then
          ok=.false.;message='periodic-box spatial and candidate symmetry representations are not homomorphic';return
        endif
      enddo;enddo
      orthogonal_transform=matmul(conjg(transpose(x)),matmul(localizer,x))
      allocate(symmetrized_localizer(ncandidate,ncandidate));symmetrized_localizer=(0d0,0d0)
      do isym=1,nsym
        symmetrized_localizer=symmetrized_localizer+matmul(conjg(transpose(candidate_symmetry(:,:,isym))),&
          matmul(orthogonal_transform,candidate_symmetry(:,:,isym)))
      enddo
      symmetrized_localizer=symmetrized_localizer/real(nsym,real64)
      localizer=matmul(s,matmul(x,matmul(symmetrized_localizer,&
        matmul(conjg(transpose(x)),s))))
      localizer=0.5d0*(localizer+conjg(transpose(localizer)))
      if(.not.finite_complex_matrix(localizer))then
        ok=.false.;message='nonfinite symmetry-averaged periodic localization matrix';return
      endif
    endif

    occ_gram=matmul(conjg(transpose(occupied_coefficients)),matmul(s,occupied_coefficients))
    occ_gram=0.5d0*(occ_gram+conjg(transpose(occ_gram)))
    if(.not.finite_complex_matrix(occ_gram))then
      ok=.false.;message='nonfinite occupied periodic-box Gram matrix';return
    endif
    if(.not.positive_definite_above(occ_gram,rank_tolerance))then
      ok=.false.;message='occupied candidate rank loss';return
    endif
    call hermitian_eigensystem(occ_gram,occ_spectrum,occ_vectors,ok,message)
    if(.not.ok)return
    if(minval(occ_spectrum)<=rank_tolerance*maxval(occ_spectrum))then
      ok=.false.;message='occupied candidate rank loss';return
    endif
    a_occ=matmul(occupied_coefficients,occ_vectors)
    do j=1,noccupied
      a_occ(:,j)=a_occ(:,j)/sqrt(occ_spectrum(j))
    enddo
    occ_gram=matmul(conjg(transpose(a_occ)),matmul(localizer,a_occ))
    occ_gram=0.5d0*(occ_gram+conjg(transpose(occ_gram)))
    if(.not.finite_complex_matrix(occ_gram))then
      ok=.false.;message='nonfinite occupied periodic localization matrix';return
    endif
    call hermitian_eigensystem(occ_gram,occ_spectrum,occ_vectors,ok,message)
    if(.not.ok)return
    a_occ=matmul(a_occ,occ_vectors)

    overlap_occ_x=matmul(conjg(transpose(a_occ)),matmul(s,x))
    residual=x-matmul(a_occ,overlap_occ_x)
    residual_gram=matmul(conjg(transpose(residual)),matmul(s,residual))
    residual_gram=0.5d0*(residual_gram+conjg(transpose(residual_gram)))
    if(.not.finite_complex_matrix(residual_gram))then
      ok=.false.;message='nonfinite periodic localization complement Gram matrix';return
    endif
    call hermitian_eigensystem(residual_gram,residual_spectrum,residual_vectors,ok,message)
    if(.not.ok)return
    ncomp=count(residual_spectrum>rank_tolerance*max(1d0,maxval(residual_spectrum)))
    nneed=ntarget-noccupied
    if(ncomp<nneed)then;ok=.false.;message='target rank loss in localization complement';return;endif
    if(nneed>0)then
      start_index=ncandidate-ncomp+1
      q_comp=matmul(residual,residual_vectors(:,start_index:ncandidate))
      do j=1,ncomp
        q_comp(:,j)=q_comp(:,j)/sqrt(residual_spectrum(start_index+j-1))
      enddo
      comp_localizer=matmul(conjg(transpose(q_comp)),matmul(localizer,q_comp))
      comp_localizer=0.5d0*(comp_localizer+conjg(transpose(comp_localizer)))
      if(.not.finite_complex_matrix(comp_localizer))then
        ok=.false.;message='nonfinite periodic localization complement matrix';return
      endif
      call hermitian_eigensystem(comp_localizer,comp_spectrum,comp_vectors,ok,message)
      if(.not.ok)return
      allocate(transform(ncandidate,ntarget))
      transform(:,1:noccupied)=a_occ
      transform(:,noccupied+1:ntarget)=matmul(q_comp,comp_vectors(:,1:nneed))
    else
      allocate(transform,source=a_occ)
    endif
    result%symmetry_closure_residual=0d0
    if(nsym>0)then
      orthogonal_transform=matmul(conjg(transpose(x)),matmul(s,transform))
      allocate(result%symmetry_representation(ntarget,ntarget,nsym))
      do isym=1,nsym
        mapped_transform=matmul(candidate_symmetry(:,:,isym),orthogonal_transform)
        result%symmetry_representation(:,:,isym)=matmul(conjg(transpose(orthogonal_transform)),&
          mapped_transform)
        result%symmetry_closure_residual=max(result%symmetry_closure_residual,&
          maxval(abs(mapped_transform-matmul(orthogonal_transform,&
          result%symmetry_representation(:,:,isym)))))
      enddo
      if(result%symmetry_closure_residual>active_symmetry_tolerance)then
        ok=.false.;message='retained periodic-box Wannier space is not symmetry closed';return
      endif
    endif

    result%candidate_rank=ncandidate;result%target_rank=ntarget
    result%retained_rank=ntarget;result%generation=generation
    allocate(result%physical_grid_ids,source=physical_ids)
    allocate(result%transform,source=transform)
    allocate(result%value(ntarget,nlocal),result%gradient(3,ntarget,nlocal))
    result%value=matmul(transpose(transform),candidate_value)
    do p=1,nlocal;do i=1,3
      result%gradient(i,:,p)=matmul(transpose(transform),candidate_gradient(i,:,p))
    enddo;enddo
    if(.not.finite_complex_matrix(result%value).or.&
        .not.all(ieee_is_finite(real(result%gradient))).or.&
        .not.all(ieee_is_finite(aimag(result%gradient))))then
      ok=.false.;message='nonfinite periodic-box Wannier value or gradient tails';return
    endif
    if(nsym>0)then
      all_wannier=matmul(transpose(transform),all_candidate)
      if(.not.finite_complex_matrix(all_wannier))then
        ok=.false.;message='nonfinite periodic-box Wannier tails';return
      endif
      allocate(result%center_box_point_ids(ntarget))
      do j=1,ntarget
        largest=maxval(abs(all_wannier(j,:))**2);center_index=0
        if(largest<=0d0.or..not.ieee_is_finite(largest))then
          ok=.false.;message='periodic-box Wannier has no finite nonzero center density';return
        endif
        do source=1,total_box
          if(largest-abs(all_wannier(j,source))**2<=&
              active_symmetry_tolerance*largest)then
            center_index=source;exit
          endif
        enddo
        if(.not.canonical_core(center_index))then
          ok=.false.;message='periodic-box Wannier center is not owned by the fragment core';return
        endif
        result%center_box_point_ids(j)=int(center_index,int64)
      enddo
      allocate(symmetry_mapped_wannier(total_box))
      do isym=1,nsym;do j=1,ntarget
        symmetry_mapped_wannier=matmul(result%symmetry_representation(:,j,isym),all_wannier)
        largest=maxval(abs(symmetry_mapped_wannier)**2)
        if(largest<=0d0.or..not.ieee_is_finite(largest))then
          ok=.false.;message='symmetry-mapped periodic-box Wannier has no finite center density';return
        endif
        target=int(canonical_target_ids(int(result%center_box_point_ids(j)),isym))
        if(largest-abs(symmetry_mapped_wannier(target))**2 .gt. &
            active_symmetry_tolerance*largest)then
          ok=.false.;message='periodic-box Wannier centers do not form the required symmetry orbit';return
        endif
      enddo;enddo
    endif

    final_metric=matmul(conjg(transpose(transform)),matmul(s,transform))
    final_metric=0.5d0*(final_metric+conjg(transpose(final_metric)))
    if(.not.finite_complex_matrix(final_metric))then
      ok=.false.;message='nonfinite retained periodic-box Wannier metric';return
    endif
    call hermitian_eigensystem(final_metric,metric_spectrum,metric_vectors,ok,message)
    if(.not.ok)return
    result%metric_minimum_eigenvalue=minval(metric_spectrum)
    if(result%metric_minimum_eigenvalue<=rank_tolerance*maxval(metric_spectrum))then
      ok=.false.;message='retained periodic-box Wannier metric rank loss';return
    endif
    result%metric_condition_number=maxval(metric_spectrum)/result%metric_minimum_eigenvalue

    occ_gram=matmul(conjg(transpose(occupied_coefficients)),matmul(s,occupied_coefficients))
    overlap_occ_x=matmul(conjg(transpose(transform)),matmul(s,occupied_coefficients))
    result%occupied_inclusion_residual=maxval(abs(occ_gram-&
      matmul(conjg(transpose(overlap_occ_x)),overlap_occ_x)))/max(1d0,maxval(abs(occ_gram)))
    if(.not.ieee_is_finite(result%occupied_inclusion_residual))then
      ok=.false.;message='nonfinite occupied inclusion residual';return
    endif
    if(result%occupied_inclusion_residual>rank_tolerance)then
      ok=.false.;message='occupied inclusion tolerance exceeded';return
    endif
    local_boundary_value=0d0;local_boundary_gradient=0d0
    do p=1,nlocal
      if(.not.boundary_mask(p))cycle
      local_boundary_value=max(local_boundary_value,maxval(abs(result%value(:,p))))
      local_boundary_gradient=max(local_boundary_gradient,maxval(abs(result%gradient(:,:,p))))
    enddo
    call MPI_Allreduce(local_boundary_value,result%boundary_value_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_boundary_gradient,result%boundary_gradient_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(result%boundary_value_max>boundary_value_tolerance.or.&
        result%boundary_gradient_max>boundary_gradient_tolerance)then
      ok=.false.;message='buffer-boundary value or gradient tolerance exceeded';return
    endif

    call MPI_Comm_rank(comm,rank,ierr)
    allocate(result%center_owner_rank(ntarget),result%center_owner_fragment(ntarget))
    if(nsym>0)then
      do j=1,ntarget
        result%center_owner_rank(j)=canonical_ranks(int(result%center_box_point_ids(j)))
        result%center_owner_fragment(j)=canonical_fragments(int(result%center_box_point_ids(j)))
      enddo
    else
      result%center_owner_rank=rank
      result%center_owner_fragment=fragment_min
    endif

    local_fingerprint=0_int64
    if(rank==0)local_fingerprint=ieor(int(generation,int64),ishftc(int(ntarget,int64),11))
    do p=1,nlocal
      quantized_density=nint(sum(abs(result%value(:,p))**2)*1d10,int64)
      if(present(box_point_ids))then
        local_fingerprint=ieor(local_fingerprint,ieor(ishftc(box_point_ids(p),17),quantized_density))
      else
        local_fingerprint=ieor(local_fingerprint,ieor(ishftc(physical_ids(p),17),quantized_density))
      endif
    enddo
    call MPI_Allreduce(local_fingerprint,global_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    result%transform_fingerprint=global_fingerprint
    if(result%transform_fingerprint==0_int64)result%transform_fingerprint=1_int64
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier construction requires MPI'
#endif
  end subroutine

  function identity_complex(n) result(identity)
    integer,intent(in)::n
    complex(real64)::identity(n,n)
    integer::i
    identity=(0d0,0d0)
    do i=1,n;identity(i,i)=1d0;enddo
  end function identity_complex

  logical function finite_complex_matrix(matrix)
    complex(real64),intent(in)::matrix(:,:)
    finite_complex_matrix=all(ieee_is_finite(real(matrix))).and.all(ieee_is_finite(aimag(matrix)))
  end function finite_complex_matrix

  subroutine sort_ids(ids)
    integer(int64),intent(inout)::ids(:)
    integer::i,j
    integer(int64)::key
    do i=2,size(ids)
      key=ids(i);j=i-1
      do while(j>=1)
        if(ids(j)<=key)exit
        ids(j+1)=ids(j);j=j-1
      enddo
      ids(j+1)=key
    enddo
  end subroutine sort_ids

  logical function positive_definite_above(matrix,relative_tolerance)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),intent(in)::relative_tolerance
    complex(real64),allocatable::factor(:,:)
    real(real64)::pivot,scale
    integer::i,j,n
    n=size(matrix,1);positive_definite_above=.false.
    if(n<=0.or.size(matrix,2)/=n)return
    scale=max(1d0,maxval(abs([(real(matrix(i,i)),i=1,n)])))
    allocate(factor(n,n));factor=(0d0,0d0)
    do j=1,n
      pivot=real(matrix(j,j))-sum(abs(factor(j,1:j-1))**2)
      if(pivot<=relative_tolerance*scale)return
      factor(j,j)=sqrt(pivot)
      do i=j+1,n
        factor(i,j)=(matrix(i,j)-sum(factor(i,1:j-1)*conjg(factor(j,1:j-1))))/factor(j,j)
      enddo
    enddo
    positive_definite_above=.true.
  end function positive_definite_above

  subroutine hermitian_eigensystem(matrix,eigenvalues,eigenvectors,ok,message)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable,intent(out)::eigenvectors(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::work(:)
    real(real64),allocatable::rwork(:)
    integer::n,lwork,info
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character(1),intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface
    n=size(matrix,1);ok=.false.;message=''
    if(n<=0.or.size(matrix,2)/=n)then;message='invalid Hermitian eigensystem shape';return;endif
    allocate(eigenvectors,source=matrix);allocate(eigenvalues(n),rwork(max(1,3*n-2)),work(1))
    lwork=-1;call zheev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0)then;message='Hermitian workspace query failed';return;endif
    lwork=max(1,int(real(work(1))));deallocate(work);allocate(work(lwork))
    call zheev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0)then;message='Hermitian eigensystem failed';return;endif
    ok=.true.
  end subroutine

  subroutine release_dg_overlapping_wannier_construction(result)
    type(s_dg_overlapping_wannier_construction),intent(inout)::result
    if(allocated(result%center_owner_rank))deallocate(result%center_owner_rank)
    if(allocated(result%center_owner_fragment))deallocate(result%center_owner_fragment)
    if(allocated(result%physical_grid_ids))deallocate(result%physical_grid_ids)
    if(allocated(result%center_box_point_ids))deallocate(result%center_box_point_ids)
    if(allocated(result%value))deallocate(result%value)
    if(allocated(result%gradient))deallocate(result%gradient)
    if(allocated(result%transform))deallocate(result%transform)
    if(allocated(result%symmetry_representation))deallocate(result%symmetry_representation)
  end subroutine
end module dg_overlapping_wannier_construction
