#include "config.h"
module dg_overlapping_wannier_construction
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
#ifdef USE_EIGENEXA
  use structures,only:s_parallel_info
  use dg_overlapping_wannier_metric,only:assemble_dg_eigenexa_cyclic_metric_block
  use eigen_eigenexa,only:eigen_pdsyevd_ex_distributed_blocks
  use eigen_libs_mod,only:eigen_owner_node,eigen_translate_g2l,eigen_translate_l2g,&
    eigen_loop_start,eigen_loop_end,eigen_get_matdims
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
    real(real64)::projection_inclusion_residual=huge(1d0)
    real(real64)::symmetry_closure_residual=huge(1d0)
    real(real64)::boundary_value_max=huge(1d0),boundary_gradient_max=huge(1d0)
    real(real64)::metric_minimum_eigenvalue=0d0,metric_condition_number=huge(1d0)
  end type
  type,public::s_dg_translation_orbit_accumulator
    logical::initialized=.false.
    logical,allocatable::visited(:)
    integer(int64)::catalog_fingerprint=0_int64
    integer(int64)::table_fingerprint=0_int64
  end type
  type,public::s_dg_prepared_translation_action
    integer::global_row_count=0,ntranslation=0,ngenerator=0,identity_operation=0
    integer::construction_collective_count=0
    integer(int64)::catalog_fingerprint=0_int64,workspace_peak_bytes=0_int64
    integer(int64),allocatable::row_ids(:)
    integer,allocatable::generator_orders(:),element_words(:,:),product_table(:,:),&
      generator_maps(:,:),element_maps(:,:)
  end type
  type,public::s_dg_prepared_spectral_basins
    logical::initialized=.false.
    integer::global_row_count=0,nlocal=0,nstate=0,basin_count=0
    integer(int64)::state_frame_fingerprint=0_int64,basin_fingerprint=0_int64
    integer(int64)::workspace_peak_bytes=0_int64
    real(real64)::state_frame_defect=huge(1d0),tolerance=0d0
    integer,allocatable::basin_offsets(:),point_indices(:)
  end type
  public::construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction
  public::verify_dg_overlapping_wannier_periodic_closure
  public::assemble_dg_distributed_candidate_symmetry
  public::assemble_dg_distributed_basis_symmetry_overlap
  public::assemble_dg_distributed_basis_symmetry_overlap_rows
  public::gather_dg_single_symmetry_representation
  public::validate_dg_row_owned_group_representation
  public::validate_dg_streamed_affine_representation
  public::build_dg_pointwise_affine_owner_map
  public::select_dg_fixed_rank_symmetry_closed_subspace
  public::build_dg_distributed_symmetry_closed_basis
  public::build_dg_group_averaged_occupied_candidates_reference
  public::orthonormalize_dg_distributed_seed_space
  public::align_dg_fragment_wannier_gauge
  public::replicate_dg_fragment_wannier_representative
  public::verify_dg_fragment_center_orbit
  public::verify_dg_fragment_subspace_density_covariance
  public::verify_dg_fragment_wannier_streaming_closure
  public::build_dg_core_owned_occupied_subspace
  public::verify_dg_uniform_fragment_target_rank
  public::assign_dg_overlapping_wannier_occupations
  public::find_dg_group_identity
  public::select_dg_group_generators
  public::build_dg_smooth_partition_of_unity
  public::compose_dg_buffered_orbital_tile_to_physical_grid
  public::accumulate_dg_lcfo_buffer_contributions_to_core
  public::measure_dg_rank_fixed_symmetry_residuals
#ifdef USE_EIGENEXA
  public::measure_dg_rank_fixed_symmetry_residuals_eigenexa
  public::build_dg_group_averaged_occupied_candidates_eigenexa
  public::build_dg_cocycle_averaged_occupied_candidates_eigenexa
  public::split_dg_translation_character_sector_eigenexa
#endif
  public::exchange_dg_point_permuted_orbital_rows
  public::accept_dg_boundary_calibrated_symmetry
  public::solve_dg_affine_common_fixed_point
  public::compute_dg_periodic_wannier_centers
  public::verify_dg_wannier_center_affine_orbits
  public::diagnose_dg_point_center_gauge
  public::build_dg_finite_abelian_character_table
  public::inverse_dg_translation_character_orbits
  public::accumulate_dg_translation_character_orbit_sector
  public::accumulate_dg_translation_character_orbit_sector_values
  public::apply_dg_row_owned_orbital_transform_streamed
  public::validate_dg_factored_point_cogroup_gauge
  public::build_dg_translation_character_intertwining_phase
  public::prepare_dg_translation_character_action
  public::build_dg_translation_character_intertwining_phase_prepared
  public::release_dg_prepared_translation_action
  public::materialize_dg_row_owned_sector_on_spatial_grid
  public::validate_dg_translation_sector_cluster
  public::build_dg_balanced_orbital_ownership
  public::transpose_dg_spatial_cores_to_orbital_owners
  public::redistribute_dg_owned_orbitals_to_center_fragments
  public::redistribute_dg_buffer_orbitals_to_center_fragments
  public::assign_dg_periodic_centers_to_fragments
  public::build_dg_equal_count_spectral_windows
  public::build_dg_spectral_density_descriptors
  public::build_dg_occupied_empty_moment_descriptors
  public::build_dg_periodic_spectral_basins
  public::project_dg_single_spectral_basin_operator
  public::prepare_dg_spectral_basin_operators,project_dg_prepared_spectral_basin_operator,&
    release_dg_prepared_spectral_basins
  public::diagonalize_dg_spectral_basin_operator
  public::select_dg_spectral_basin_channel_ranks
contains

  subroutine select_dg_spectral_basin_channel_ranks(comm,spectra,block_ends,basin_generator_maps,&
      retained_rank,tolerance,selected_ranks,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,basin_generator_maps(:,:),retained_rank
    real(real64),intent(in)::spectra(:,:),tolerance
    logical,intent(in)::block_ends(:,:)
    integer,intent(out)::selected_ranks(:)
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,nbasin,ngenerator,i,j,g,b,target,head,tail,norbit,o,r,total,new_total,status
    integer::ierr,local_bad,global_bad,minint,maxint
    integer,allocatable::orbit_id(:),queue(:),orbit_size(:),representative(:),choice(:,:),prior(:,:)
    real(real64),allocatable::score(:),next_score(:)
    integer(int64)::elements,bytes,minhash,maxhash,hash_value,quantized
    real(real64)::minreal,maxreal,scale,candidate,quantum

    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    selected_ranks=0;nstate=size(spectra,1);nbasin=size(spectra,2);ngenerator=size(basin_generator_maps,2)
    local_bad=0
    if(nstate<1.or.nbasin<1.or.retained_rank/=nstate.or.size(selected_ranks)/=nbasin.or.&
        any(shape(block_ends)/=shape(spectra)).or.size(basin_generator_maps,1)/=nbasin.or.&
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.tolerance>huge(1d0)/100d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(spectra)).or.any(basin_generator_maps<1).or.&
        any(basin_generator_maps>nbasin))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid spectral basin channel catalog contract';return;endif
    do i=1,4
      select case(i)
      case(1);local_bad=nstate
      case(2);local_bad=nbasin
      case(3);local_bad=ngenerator
      case default;local_bad=retained_rank
      end select
      call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin channel catalog dimensions disagree';return;endif
    enddo
    call MPI_Allreduce(tolerance,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
      message='spectral basin channel catalog tolerance disagrees';return
    endif
    do b=1,nbasin
      do i=1,nstate
        call MPI_Allreduce(spectra(i,b),minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
        if(ierr/=MPI_SUCCESS)return
        call MPI_Allreduce(spectra(i,b),maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS)return
        if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
          message='spectral basin spectra disagree across ranks';return
        endif
        local_bad=merge(1,0,block_ends(i,b))
        call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
        call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
        if(minint/=maxint)then;message='spectral basin block boundaries disagree across ranks';return;endif
      enddo
    enddo
    do g=1,ngenerator;do b=1,nbasin
      call MPI_Allreduce(basin_generator_maps(b,g),minint,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(basin_generator_maps(b,g),maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin orbit maps disagree across ranks';return;endif
    enddo;enddo
    local_bad=0
    do b=1,nbasin
      scale=max(1d0,maxval(abs(spectra(:,b))))
      if(.not.block_ends(nstate,b).or.any(spectra(2:nstate,b)>spectra(1:nstate-1,b)+tolerance*scale).or.&
          minval(spectra(:,b))< -10d0*tolerance*scale)local_bad=1
    enddo
    do g=1,ngenerator
      do b=1,nbasin
        if(count(basin_generator_maps(:,g)==b)/=1)local_bad=1
      enddo
    enddo
    if(nstate>=huge(0))then
      local_bad=1
    elseif(int(nstate,int64)>huge(0_int64)/int(nbasin,int64))then
      local_bad=1
    elseif((int(nstate,int64)+1_int64)*int(nbasin,int64)>int(huge(0),int64))then
      local_bad=1
    endif
    scale=huge(1d0)/16d0/real(nbasin,real64)/real(nstate,real64)
    if(maxval(abs(spectra))>scale)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid spectral basin spectra or orbit action';return;endif
    elements=int(nstate,int64)*int(nbasin,int64)
    allocate(orbit_id(nbasin),queue(nbasin),orbit_size(nbasin),representative(nbasin),&
      choice(nbasin,0:nstate),prior(nbasin,0:nstate),score(0:nstate),next_score(0:nstate),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(orbit_id))deallocate(orbit_id)
      if(allocated(queue))deallocate(queue)
      if(allocated(orbit_size))deallocate(orbit_size)
      if(allocated(representative))deallocate(representative)
      if(allocated(choice))deallocate(choice)
      if(allocated(prior))deallocate(prior)
      if(allocated(score))deallocate(score)
      if(allocated(next_score))deallocate(next_score)
      message='spectral basin channel catalog workspace allocation failed';return
    endif
    orbit_id=0;orbit_size=0;representative=0;norbit=0;local_bad=0
    do b=1,nbasin
      if(orbit_id(b)/=0)cycle
      norbit=norbit+1;representative(norbit)=b;head=1;tail=1;queue(1)=b;orbit_id(b)=norbit
      do while(head<=tail)
        target=queue(head);head=head+1
        do g=1,ngenerator
          i=basin_generator_maps(target,g)
          if(orbit_id(i)==0)then
            tail=tail+1;queue(tail)=i;orbit_id(i)=norbit
          elseif(orbit_id(i)/=norbit)then
            local_bad=1
          endif
        enddo
      enddo
      orbit_size(norbit)=tail
    enddo
    do o=1,norbit
      b=representative(o);scale=max(1d0,maxval(abs(spectra(:,b))))
      do j=1,nbasin
        if(orbit_id(j)/=o)cycle
        if(any(block_ends(:,j).neqv.block_ends(:,b)).or.&
            maxval(abs(spectra(:,j)-spectra(:,b)))>10d0*tolerance*scale)local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='symmetry-related spectral basins have inconsistent eigenspaces';return
    endif
    score=-huge(1d0);score(0)=0d0;choice=-1;prior=-1
    do o=1,norbit
      next_score=-huge(1d0);b=representative(o);scale=max(1d0,maxval(abs(spectra(:,b))))
      do total=0,nstate
        if(score(total)<=-0.5d0*huge(1d0))cycle
        do r=0,nstate
          if(r>0)then
            if(.not.block_ends(r,b).or.spectra(r,b)<=tolerance*scale)cycle
          endif
          new_total=total+orbit_size(o)*r;if(new_total>nstate)cycle
          candidate=score(total)+real(orbit_size(o),real64)*sum(spectra(1:r,b))
          if(candidate>next_score(new_total)+tolerance*max(1d0,abs(candidate)))then
            next_score(new_total)=candidate;choice(o,new_total)=r;prior(o,new_total)=total
          endif
        enddo
      enddo
      score=next_score
    enddo
    local_bad=merge(0,1,choice(norbit,nstate)>=0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='spectral basin channel catalog cannot span the retained rank without splitting a block';return
    endif
    total=nstate
    do o=norbit,1,-1
      r=choice(o,total)
      do b=1,nbasin;if(orbit_id(b)==o)selected_ranks(b)=r;enddo
      total=prior(o,total)
    enddo
    if(total/=0.or.sum(selected_ranks)/=nstate)then
      message='spectral basin channel catalog rank reconstruction failed';return
    endif
    quantum=100d0*tolerance;hash_value=int(nstate,int64)
    do b=1,nbasin
      hash_value=ieor(ishftc(hash_value,7),int(selected_ranks(b),int64))
      do i=1,nstate
        if(abs(spectra(i,b))>0.25d0*real(huge(0_int64),real64)*quantum)then
          message='spectral basin channel catalog fingerprint range is unsafe';return
        endif
        quantized=nint(spectra(i,b)/quantum,int64);hash_value=ieor(ishftc(hash_value,11),quantized)
        hash_value=ieor(ishftc(hash_value,5),merge(1_int64,0_int64,block_ends(i,b)))
      enddo
    enddo
    do g=1,ngenerator;do b=1,nbasin
      hash_value=ieor(ishftc(hash_value,3),int(basin_generator_maps(b,g),int64))
    enddo;enddo
    if(hash_value==0_int64)hash_value=1_int64
    bytes=4_int64*(4_int64*int(nbasin,int64)+2_int64*int(nbasin,int64)*int(nstate+1,int64))+&
      16_int64*int(nstate+1,int64)
    call MPI_Allreduce(bytes,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(bytes,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minhash/=maxhash)then;message='spectral basin channel catalog workspace disagrees';return;endif
    fingerprint=hash_value;workspace_peak_bytes=maxhash;ok=.true.
    deallocate(orbit_id,queue,orbit_size,representative,choice,prior,score,next_score)
#else
    ok=.false.;message='spectral basin channel catalog requires MPI';fingerprint=0_int64
    workspace_peak_bytes=0_int64;selected_ranks=0
#endif
  end subroutine select_dg_spectral_basin_channel_ranks

  subroutine diagonalize_dg_spectral_basin_operator(comm,basin_operator,operator_fingerprint,tolerance,&
      spectrum,block_offsets,eigensystem_residual,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(inout)::basin_operator(:,:)
    integer(int64),intent(in)::operator_fingerprint
    real(real64),intent(in)::tolerance
    real(real64),allocatable,intent(out)::spectrum(:)
    integer,allocatable,intent(out)::block_offsets(:)
    real(real64),intent(out)::eigensystem_residual
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::original_operator(:,:),eigenvectors(:,:),work(:),residual_matrix(:,:),column(:)
    real(real64),allocatable::rwork(:)
    integer,allocatable::offset_workspace(:)
    complex(real64)::work_query(1)
    integer::n,rank,i,j,lwork,lapack_info,ierr,status,local_bad,global_bad,minint,maxint,nblock
    integer(int64)::minhash,maxhash,elements,bytes,root_bytes,hash_value,quantized
    real(real64)::minreal,maxreal,scale,quantum
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character(1),intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine zheev
    end interface

    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    eigensystem_residual=huge(1d0);n=size(basin_operator,1)
    local_bad=0
    if(n<1.or.size(basin_operator,2)/=n.or.operator_fingerprint==0_int64.or.&
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.tolerance>huge(1d0)/100d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(real(basin_operator))).or.&
        .not.all(ieee_is_finite(aimag(basin_operator))))then
      local_bad=1
    elseif(maxval(abs(basin_operator-conjg(transpose(basin_operator))))>10d0*tolerance)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid spectral basin eigensystem contract';return;endif
    call MPI_Allreduce(n,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(n,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minint/=maxint)then;message='spectral basin eigensystem dimension disagrees';return;endif
    call MPI_Allreduce(operator_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(operator_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minhash/=maxhash)then;message='spectral basin operator provenance disagrees';return;endif
    call MPI_Allreduce(tolerance,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
      message='spectral basin eigensystem tolerance disagrees';return
    endif
    local_bad=0
    if(int(n,int64)>huge(0_int64)/int(n,int64))then
      local_bad=1
    else
      elements=int(n,int64)*int(n,int64)
      if(elements>int(huge(0),int64).or.int(n,int64)>int(huge(0),int64)/5_int64)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin eigensystem extent overflow';return;endif
    allocate(spectrum(n),offset_workspace(n+1),column(n),stat=status)
    local_bad=merge(0,1,status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(spectrum))deallocate(spectrum)
      if(allocated(offset_workspace))deallocate(offset_workspace)
      if(allocated(column))deallocate(column)
      message='spectral basin eigensystem output allocation failed';return
    endif
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='spectral basin eigensystem communicator query failed';return;endif
    lapack_info=0;lwork=1
    if(rank==0)then
      allocate(original_operator(n,n),eigenvectors(n,n),rwork(max(1,3*n-2)),stat=status)
      if(status==0)then
        original_operator=basin_operator;eigenvectors=basin_operator
        call zheev('V','U',n,eigenvectors,n,spectrum,work_query,-1,rwork,lapack_info)
        if(lapack_info==0.and.ieee_is_finite(real(work_query(1))).and.real(work_query(1))>=1d0.and.&
            real(work_query(1))<=real(huge(0),real64))then
          lwork=ceiling(real(work_query(1)))
          allocate(work(lwork),residual_matrix(n,n),stat=status)
        else
          status=1
        endif
      endif
    else
      status=0
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      if(allocated(spectrum))deallocate(spectrum)
      if(allocated(offset_workspace))deallocate(offset_workspace)
      if(allocated(column))deallocate(column)
      if(allocated(original_operator))deallocate(original_operator)
      if(allocated(eigenvectors))deallocate(eigenvectors)
      if(allocated(rwork))deallocate(rwork)
      if(allocated(work))deallocate(work)
      if(allocated(residual_matrix))deallocate(residual_matrix)
      message='spectral basin eigensystem workspace allocation failed';return
    endif
    if(rank==0)then
      call zheev('V','U',n,eigenvectors,n,spectrum,work,lwork,rwork,lapack_info)
      if(lapack_info==0)then
        residual_matrix=matmul(original_operator,eigenvectors)
        do j=1,n;residual_matrix(:,j)=residual_matrix(:,j)-spectrum(j)*eigenvectors(:,j);enddo
        eigensystem_residual=maxval(abs(residual_matrix))/max(1d0,maxval(abs(original_operator)))
        do j=1,n/2
          column=eigenvectors(:,j);eigenvectors(:,j)=eigenvectors(:,n-j+1);eigenvectors(:,n-j+1)=column
          scale=spectrum(j);spectrum(j)=spectrum(n-j+1);spectrum(n-j+1)=scale
        enddo
        basin_operator=eigenvectors
      endif
    endif
    call MPI_Bcast(lapack_info,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(lapack_info/=0)then;message='spectral basin Hermitian eigensystem failed';return;endif
    call MPI_Bcast(spectrum,n,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(basin_operator,n*n,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(eigensystem_residual,1,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,all(ieee_is_finite(spectrum)).and.&
      all(ieee_is_finite(real(basin_operator))).and.all(ieee_is_finite(aimag(basin_operator))).and.&
      ieee_is_finite(eigensystem_residual).and.eigensystem_residual<=10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin eigensystem residual failed';return;endif
    scale=max(1d0,maxval(abs(spectrum)));nblock=1;offset_workspace(1)=1
    do i=2,n
      if(abs(spectrum(i)-spectrum(i-1))>tolerance*scale)then
        nblock=nblock+1;offset_workspace(nblock)=i
      endif
    enddo
    nblock=nblock+1;offset_workspace(nblock)=n+1
    allocate(block_offsets(nblock),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(block_offsets))deallocate(block_offsets)
      message='spectral basin block allocation failed';return
    endif
    block_offsets=offset_workspace(1:nblock)
    quantum=100d0*tolerance;local_bad=0
    if(maxval(abs(spectrum))>0.25d0*real(huge(0_int64),real64)*quantum)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin spectrum fingerprint range is unsafe';return;endif
    hash_value=operator_fingerprint
    do i=1,n
      quantized=nint(spectrum(i)/quantum,int64);hash_value=ieor(ishftc(hash_value,7),quantized)
    enddo
    do i=1,nblock;hash_value=ieor(ishftc(hash_value,11),int(block_offsets(i),int64));enddo
    if(hash_value==0_int64)hash_value=1_int64
    root_bytes=48_int64*elements+16_int64*int(max(1,lwork),int64)+&
      8_int64*int(max(1,3*n-2),int64)+16_int64*int(n,int64)+&
      8_int64*int(n,int64)+8_int64*int(n+1,int64)
    bytes=16_int64*int(n,int64)+8_int64*int(n,int64)+8_int64*int(n+1,int64)
    call MPI_Allreduce(root_bytes,workspace_peak_bytes,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    workspace_peak_bytes=max(workspace_peak_bytes,bytes)
    fingerprint=hash_value;ok=.true.
    if(allocated(original_operator))deallocate(original_operator)
    if(allocated(eigenvectors))deallocate(eigenvectors)
    if(allocated(work))deallocate(work)
    if(allocated(rwork))deallocate(rwork)
    if(allocated(residual_matrix))deallocate(residual_matrix)
    deallocate(offset_workspace,column)
#else
    ok=.false.;message='spectral basin eigensystem requires MPI';fingerprint=0_int64
    workspace_peak_bytes=0_int64;eigensystem_residual=huge(1d0)
#endif
  end subroutine diagonalize_dg_spectral_basin_operator

  subroutine release_dg_prepared_spectral_basins(prepared)
    type(s_dg_prepared_spectral_basins),intent(inout)::prepared
    if(allocated(prepared%basin_offsets))deallocate(prepared%basin_offsets)
    if(allocated(prepared%point_indices))deallocate(prepared%point_indices)
    prepared%initialized=.false.;prepared%global_row_count=0;prepared%nlocal=0;prepared%nstate=0
    prepared%basin_count=0;prepared%state_frame_fingerprint=0_int64;prepared%basin_fingerprint=0_int64
    prepared%workspace_peak_bytes=0_int64;prepared%state_frame_defect=huge(1d0);prepared%tolerance=0d0
  end subroutine release_dg_prepared_spectral_basins

  subroutine prepare_dg_spectral_basin_operators(comm,row_ids,global_row_count,state_values,point_weights,&
      basin_labels,basin_count,state_frame_fingerprint,state_frame_defect,basin_fingerprint,tolerance,&
      prepared,ok,message)
    integer,intent(in)::comm,global_row_count,basin_labels(:),basin_count
    integer(int64),intent(in)::row_ids(:),state_frame_fingerprint,basin_fingerprint
    complex(real64),intent(in)::state_values(:,:)
    real(real64),intent(in)::point_weights(:),state_frame_defect,tolerance
    type(s_dg_prepared_spectral_basins),intent(inout)::prepared
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nstate,i,b,p,ierr,local_bad,global_bad,minint,maxint,status
    integer,allocatable::ownership(:),counts(:),cursor(:)
    integer(int64)::minhash,maxhash,bytes
    real(real64)::minreal,maxreal,local_max,global_max,local_weight_max,global_weight_max,safe_coefficient
    nlocal=size(row_ids);nstate=size(state_values,1);ok=.false.;message=''
    call release_dg_prepared_spectral_basins(prepared)
    local_bad=0
    if(global_row_count<1.or.nstate<1.or.basin_count<1.or.size(state_values,2)/=nlocal.or.&
        size(point_weights)/=nlocal.or.size(basin_labels)/=nlocal)then
      local_bad=1
    elseif(state_frame_fingerprint==0_int64.or.basin_fingerprint==0_int64)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.tolerance>huge(1d0)/100d0.or.&
        .not.ieee_is_finite(state_frame_defect).or.state_frame_defect<0d0.or.state_frame_defect>tolerance)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(real(state_values))).or..not.all(ieee_is_finite(aimag(state_values))).or.&
        .not.all(ieee_is_finite(point_weights)))then
      local_bad=1
    elseif(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)).or.any(point_weights<0d0).or.&
        any(basin_labels<1).or.any(basin_labels>basin_count))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid prepared spectral basin contract';return;endif
    do i=1,3
      if(i==1)then;local_bad=global_row_count;elseif(i==2)then;local_bad=nstate;else;local_bad=basin_count;endif
      call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='prepared spectral basin metadata disagrees';return;endif
    enddo
    do i=1,2
      if(i==1)then;bytes=state_frame_fingerprint;else;bytes=basin_fingerprint;endif
      call MPI_Allreduce(bytes,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bytes,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minhash/=maxhash)then;message='prepared spectral basin provenance disagrees';return;endif
    enddo
    do i=1,2
      if(i==1)then;local_max=tolerance;else;local_max=state_frame_defect;endif
      call MPI_Allreduce(local_max,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_max,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
        message='prepared spectral basin receipts disagree';return
      endif
    enddo
    local_bad=0
    if(int(nstate,int64)>huge(0_int64)/int(nstate,int64).or.&
        int(nstate,int64)*int(nstate,int64)>int(huge(0),int64))local_bad=1
    if(int(global_row_count,int64)>huge(0_int64)/4_int64.or.&
        int(nlocal,int64)>huge(0_int64)/4_int64.or.&
        int(basin_count,int64)+1_int64>huge(0_int64)/4_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared spectral basin extent overflow';return;endif
    allocate(ownership(global_row_count),counts(basin_count),cursor(basin_count),&
      prepared%basin_offsets(basin_count+1),prepared%point_indices(nlocal),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(ownership))deallocate(ownership)
      if(allocated(counts))deallocate(counts)
      if(allocated(cursor))deallocate(cursor)
      call release_dg_prepared_spectral_basins(prepared)
      message='prepared spectral basin allocation failed';return
    endif
    ownership=0
    do p=1,nlocal;ownership(int(row_ids(p)))=ownership(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      deallocate(ownership,counts,cursor);call release_dg_prepared_spectral_basins(prepared)
      message='prepared spectral basin rows are not owned exactly once';return
    endif
    local_max=0d0;local_weight_max=0d0
    if(nlocal>0)then;local_max=maxval(abs(state_values));local_weight_max=maxval(point_weights);endif
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_weight_max,global_weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_weight_max>0d0)then
      safe_coefficient=sqrt((huge(1d0)/16d0/real(global_row_count,real64))/global_weight_max)
      local_bad=merge(0,1,global_max<=safe_coefficient)
    else
      local_bad=0
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(ownership,counts,cursor);call release_dg_prepared_spectral_basins(prepared)
      message='prepared spectral basin input magnitude is unsafe';return
    endif
    counts=0
    do p=1,nlocal;counts(basin_labels(p))=counts(basin_labels(p))+1;enddo
    prepared%basin_offsets(1)=1
    do b=1,basin_count;prepared%basin_offsets(b+1)=prepared%basin_offsets(b)+counts(b);enddo
    cursor=prepared%basin_offsets(1:basin_count)
    do p=1,nlocal
      b=basin_labels(p);prepared%point_indices(cursor(b))=p;cursor(b)=cursor(b)+1
    enddo
    prepared%initialized=.true.;prepared%global_row_count=global_row_count;prepared%nlocal=nlocal
    prepared%nstate=nstate;prepared%basin_count=basin_count
    prepared%state_frame_fingerprint=state_frame_fingerprint;prepared%basin_fingerprint=basin_fingerprint
    prepared%state_frame_defect=state_frame_defect;prepared%tolerance=tolerance
    bytes=int(global_row_count,int64)+2_int64*int(basin_count,int64)+1_int64+int(nlocal,int64)
    bytes=4_int64*bytes
    prepared%workspace_peak_bytes=bytes
    deallocate(ownership,counts,cursor);ok=.true.
#else
    call release_dg_prepared_spectral_basins(prepared);ok=.false.
    message='prepared spectral basins require MPI'
#endif
  end subroutine prepare_dg_spectral_basin_operators

  subroutine project_dg_prepared_spectral_basin_operator(comm,prepared,state_values,point_weights,&
      basin_index,basin_operator,hermiticity_defect,operator_trace,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,basin_index
    type(s_dg_prepared_spectral_basins),intent(in)::prepared
    complex(real64),intent(in)::state_values(:,:)
    real(real64),intent(in)::point_weights(:)
    complex(real64),allocatable,intent(inout)::basin_operator(:,:)
    real(real64),intent(out)::hermiticity_defect,operator_trace
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,q,p,ierr,local_bad,global_bad,minint,maxint,status
    integer(int64)::hash_value,quantized,operator_bytes
    real(real64)::quantum
    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    hermiticity_defect=huge(1d0);operator_trace=0d0
    local_bad=0
    if(.not.prepared%initialized.or.prepared%nstate<1.or.prepared%basin_count<1.or.&
        prepared%nlocal/=size(point_weights).or.&
        .not.all(shape(state_values)==[prepared%nstate,prepared%nlocal]))local_bad=1
    if(.not.allocated(prepared%basin_offsets).or..not.allocated(prepared%point_indices))then
      local_bad=1
    else
      if(size(prepared%basin_offsets)/=prepared%basin_count+1.or.&
          size(prepared%point_indices)/=prepared%nlocal)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid prepared spectral basin state';return;endif
    call MPI_Allreduce(basin_index,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(basin_index,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint.or.basin_index<1.or.basin_index>prepared%basin_count)then
      message='prepared spectral basin index disagrees';return
    endif
    local_bad=0
    if(allocated(basin_operator))then
      if(any(shape(basin_operator)/=[prepared%nstate,prepared%nstate]))local_bad=1
    else
      allocate(basin_operator(prepared%nstate,prepared%nstate),stat=status);if(status/=0)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(basin_operator))deallocate(basin_operator)
      message='prepared spectral basin operator allocation failed';return
    endif
    basin_operator=(0d0,0d0)
    do q=prepared%basin_offsets(basin_index),prepared%basin_offsets(basin_index+1)-1
      p=prepared%point_indices(q)
      do j=1,prepared%nstate;do i=1,prepared%nstate
        basin_operator(i,j)=basin_operator(i,j)+point_weights(p)*conjg(state_values(i,p))*state_values(j,p)
      enddo;enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,basin_operator,prepared%nstate*prepared%nstate,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,all(ieee_is_finite(real(basin_operator))).and.all(ieee_is_finite(aimag(basin_operator))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared spectral basin operator is nonfinite';return;endif
    hermiticity_defect=maxval(abs(basin_operator-conjg(transpose(basin_operator))))
    operator_trace=0d0
    do i=1,prepared%nstate;operator_trace=operator_trace+real(basin_operator(i,i),real64);enddo
    if(hermiticity_defect>10d0*prepared%tolerance.or.operator_trace< -10d0*prepared%tolerance)then
      message='prepared spectral basin operator failed Hermiticity';return
    endif
    basin_operator=0.5d0*(basin_operator+conjg(transpose(basin_operator)))
    quantum=100d0*prepared%tolerance;local_bad=0
    if(maxval(abs(real(basin_operator)))>0.25d0*real(huge(0_int64),real64)*quantum.or.&
        maxval(abs(aimag(basin_operator)))>0.25d0*real(huge(0_int64),real64)*quantum)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared spectral basin fingerprint range is unsafe';return;endif
    hash_value=ieor(prepared%state_frame_fingerprint,ishftc(prepared%basin_fingerprint,13))
    hash_value=ieor(hash_value,int(basin_index,int64))
    do j=1,prepared%nstate;do i=1,prepared%nstate
      quantized=nint(real(basin_operator(i,j))/quantum,int64);hash_value=ieor(ishftc(hash_value,7),quantized)
      quantized=nint(aimag(basin_operator(i,j))/quantum,int64);hash_value=ieor(ishftc(hash_value,11),quantized)
    enddo;enddo
    if(hash_value==0_int64)hash_value=1_int64
    operator_bytes=16_int64*int(prepared%nstate,int64)*int(prepared%nstate,int64)
    if(prepared%workspace_peak_bytes>huge(0_int64)-operator_bytes)then
      message='prepared spectral basin workspace receipt overflow';return
    endif
    fingerprint=hash_value;workspace_peak_bytes=prepared%workspace_peak_bytes+operator_bytes;ok=.true.
#else
    ok=.false.;message='prepared spectral basin operator requires MPI';fingerprint=0_int64
    workspace_peak_bytes=0_int64;hermiticity_defect=huge(1d0);operator_trace=0d0
#endif
  end subroutine project_dg_prepared_spectral_basin_operator

  subroutine project_dg_single_spectral_basin_operator(comm,row_ids,global_row_count,state_values,&
      point_weights,basin_labels,basin_count,basin_index,state_frame_fingerprint,state_frame_defect,&
      basin_fingerprint,tolerance,basin_operator,hermiticity_defect,operator_trace,fingerprint,&
      workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count,basin_labels(:),basin_count,basin_index
    integer(int64),intent(in)::row_ids(:),state_frame_fingerprint,basin_fingerprint
    complex(real64),intent(in)::state_values(:,:)
    real(real64),intent(in)::point_weights(:),state_frame_defect,tolerance
    complex(real64),allocatable,intent(out)::basin_operator(:,:)
    real(real64),intent(out)::hermiticity_defect,operator_trace
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,nlocal,i,j,p,ierr,local_bad,global_bad,minint,maxint,status
    integer,allocatable::ownership(:)
    integer(int64)::minhash,maxhash,elements,bytes,hash_value,quantized
    complex(real64),allocatable::local_operator(:,:)
    real(real64)::minreal,maxreal,local_max,global_max,local_weight_max,global_weight_max,&
      safe_coefficient,quantum
    logical::receipt_valid
    nstate=size(state_values,1);nlocal=size(state_values,2);ok=.false.;message=''
    fingerprint=0_int64;workspace_peak_bytes=0_int64;hermiticity_defect=huge(1d0);operator_trace=0d0
    local_bad=0
    if(global_row_count<1.or.nstate<1.or.size(row_ids)/=nlocal.or.size(point_weights)/=nlocal.or.&
        size(basin_labels)/=nlocal.or.basin_count<1.or.basin_index<1.or.basin_index>basin_count)then
      local_bad=1
    elseif(state_frame_fingerprint==0_int64.or.basin_fingerprint==0_int64)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.tolerance>huge(1d0)/100d0.or.&
        .not.ieee_is_finite(state_frame_defect).or.state_frame_defect<0d0.or.state_frame_defect>tolerance)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(real(state_values))).or..not.all(ieee_is_finite(aimag(state_values))).or.&
        .not.all(ieee_is_finite(point_weights)))then
      local_bad=1
    elseif(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)).or.any(point_weights<0d0).or.&
        any(basin_labels<1).or.any(basin_labels>basin_count))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid single spectral basin operator contract';return;endif
    do i=1,4
      select case(i)
      case(1);local_bad=global_row_count
      case(2);local_bad=nstate
      case(3);local_bad=basin_count
      case default;local_bad=basin_index
      end select
      call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin operator metadata disagrees across ranks';return;endif
    enddo
    do i=1,2
      if(i==1)then;hash_value=state_frame_fingerprint;else;hash_value=basin_fingerprint;endif
      call MPI_Allreduce(hash_value,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(hash_value,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minhash/=maxhash)then;message='spectral basin operator provenance disagrees across ranks';return;endif
    enddo
    do i=1,2
      if(i==1)then;local_max=tolerance;else;local_max=state_frame_defect;endif
      call MPI_Allreduce(local_max,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_max,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
        message='spectral basin operator tolerance receipt disagrees across ranks';return
      endif
    enddo
    receipt_valid=int(nstate,int64)<=huge(0_int64)/int(nstate,int64)
    if(receipt_valid)then
      elements=int(nstate,int64)*int(nstate,int64)
      receipt_valid=elements<=int(huge(0),int64).and.elements<=huge(0_int64)/32_int64
    else
      elements=0_int64
    endif
    if(receipt_valid)then
      bytes=32_int64*elements
      receipt_valid=int(global_row_count,int64)<=huge(0_int64)/4_int64.and.&
        bytes<=huge(0_int64)-4_int64*int(global_row_count,int64)
      if(receipt_valid)bytes=bytes+4_int64*int(global_row_count,int64)
    else
      bytes=0_int64
    endif
    call MPI_Allreduce(merge(0,1,receipt_valid),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin operator workspace overflow';return;endif
    allocate(ownership(global_row_count),local_operator(nstate,nstate),basin_operator(nstate,nstate),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(ownership))deallocate(ownership)
      if(allocated(local_operator))deallocate(local_operator)
      if(allocated(basin_operator))deallocate(basin_operator)
      message='spectral basin operator allocation failed';return
    endif
    ownership=0
    do p=1,nlocal;ownership(int(row_ids(p)))=ownership(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      deallocate(ownership,local_operator,basin_operator)
      message='spectral basin operator rows are not owned exactly once';return
    endif
    local_max=0d0;local_weight_max=0d0
    if(nlocal>0)then;local_max=maxval(abs(state_values));local_weight_max=maxval(point_weights);endif
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_weight_max,global_weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_weight_max>0d0)then
      safe_coefficient=sqrt((huge(1d0)/16d0/real(global_row_count,real64))/global_weight_max)
      local_bad=merge(0,1,global_max<=safe_coefficient)
    else
      local_bad=0
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(ownership,local_operator,basin_operator)
      message='spectral basin operator input magnitude is unsafe';return
    endif
    local_operator=(0d0,0d0)
    do p=1,nlocal
      if(basin_labels(p)/=basin_index)cycle
      do j=1,nstate;do i=1,nstate
        local_operator(i,j)=local_operator(i,j)+point_weights(p)*conjg(state_values(i,p))*state_values(j,p)
      enddo;enddo
    enddo
    call MPI_Allreduce(local_operator,basin_operator,nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;deallocate(ownership,local_operator,basin_operator);return;endif
    local_bad=merge(0,1,all(ieee_is_finite(real(basin_operator))).and.&
      all(ieee_is_finite(aimag(basin_operator))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(ownership,local_operator,basin_operator);message='spectral basin operator is nonfinite';return
    endif
    hermiticity_defect=maxval(abs(basin_operator-conjg(transpose(basin_operator))))
    operator_trace=0d0
    do i=1,nstate;operator_trace=operator_trace+real(basin_operator(i,i),real64);enddo
    if(hermiticity_defect>10d0*tolerance.or.operator_trace< -10d0*tolerance)then
      deallocate(ownership,local_operator,basin_operator);message='spectral basin operator failed Hermiticity';return
    endif
    basin_operator=0.5d0*(basin_operator+conjg(transpose(basin_operator)))
    quantum=100d0*tolerance;local_bad=0
    if(.not.ieee_is_finite(quantum).or.quantum<=0d0)then
      local_bad=1
    elseif(maxval(abs(real(basin_operator)))>0.25d0*real(huge(0_int64),real64)*quantum.or.&
        maxval(abs(aimag(basin_operator)))>0.25d0*real(huge(0_int64),real64)*quantum)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(ownership,local_operator,basin_operator);message='spectral basin operator fingerprint range is unsafe';return
    endif
    hash_value=ieor(state_frame_fingerprint,ishftc(basin_fingerprint,13))
    hash_value=ieor(hash_value,int(basin_index,int64))
    do j=1,nstate;do i=1,nstate
      quantized=nint(real(basin_operator(i,j))/quantum,int64);hash_value=ieor(ishftc(hash_value,7),quantized)
      quantized=nint(aimag(basin_operator(i,j))/quantum,int64);hash_value=ieor(ishftc(hash_value,11),quantized)
    enddo;enddo
    if(hash_value==0_int64)hash_value=1_int64
    fingerprint=hash_value;workspace_peak_bytes=bytes
    deallocate(ownership,local_operator);ok=.true.
#else
    ok=.false.;message='single spectral basin operator requires MPI';fingerprint=0_int64
    workspace_peak_bytes=0_int64;hermiticity_defect=huge(1d0);operator_trace=0d0
#endif
  end subroutine project_dg_single_spectral_basin_operator

  subroutine build_dg_periodic_spectral_basins(comm,row_ids,grid_shape,occupied_density,&
      empty_moment_density,shared_density,generator_maps,tolerance,basin_labels,basin_count,&
      basin_orbit_map,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,grid_shape(3),generator_maps(:,:)
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::occupied_density(:),empty_moment_density(:,:),shared_density(:,:),tolerance
    integer,allocatable,intent(out)::basin_labels(:),basin_orbit_map(:,:)
    integer,intent(out)::basin_count
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nfeature,nempty,nshared,ngenerator,npoint,i,j,g,p,x,y,z,neighbor,choice,root,&
      ierr,local_bad,global_bad,minint,maxint,status,target,target_basin
    integer,allocatable::ownership(:),parent(:),roots(:),global_labels(:),root_ids(:),target_counts(:)
    integer(int64)::npoint8,bits,minbits,maxbits,bytes,hash_value
    real(real64),allocatable::local_score(:),global_score(:)
    real(real64)::local_max,global_max,candidate,comparison_scale
    logical::receipt_valid
    nlocal=size(row_ids);nempty=size(empty_moment_density,2);nshared=size(shared_density,2);nfeature=nempty+nshared
    ngenerator=size(generator_maps,2);ok=.false.;message='';basin_count=0
    fingerprint=0_int64;workspace_peak_bytes=0_int64;local_bad=0
    if(any(grid_shape<1).or.size(occupied_density)/=nlocal.or.&
        size(empty_moment_density,1)/=nlocal.or.size(shared_density,1)/=nlocal)then
      local_bad=1
    elseif(size(generator_maps,1)<1.or.nfeature<1)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.tolerance>1d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(occupied_density)).or.&
        .not.all(ieee_is_finite(empty_moment_density)).or..not.all(ieee_is_finite(shared_density)))then
      local_bad=1
    elseif(any(occupied_density<0d0).or.any(empty_moment_density<0d0).or.any(shared_density<0d0))then
      local_bad=1
    endif
    npoint8=1_int64
    do i=1,3
      if(int(grid_shape(i),int64)>huge(0_int64)/npoint8)then;local_bad=1;exit;endif
      npoint8=npoint8*int(grid_shape(i),int64)
    enddo
    if(npoint8>int(huge(0),int64))local_bad=1
    if(local_bad==0)then
      npoint=int(npoint8)
      if(size(generator_maps,1)/=npoint.or.any(row_ids<1_int64).or.any(row_ids>npoint8))local_bad=1
    else
      npoint=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid periodic spectral basin contract';return;endif
    do i=1,3
      call MPI_Allreduce(grid_shape(i),minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(grid_shape(i),maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin grid shape disagrees across ranks';return;endif
    enddo
    call MPI_Allreduce(ngenerator,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(ngenerator,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minint/=maxint)then;message='spectral basin generator count disagrees across ranks';return;endif
    do i=1,2
      if(i==1)then;local_bad=nempty;else;local_bad=nshared;endif
      call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin feature count disagrees across ranks';return;endif
    enddo
    bits=transfer(tolerance,0_int64)
    call MPI_Allreduce(bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minbits/=maxbits)then;message='spectral basin tolerance disagrees across ranks';return;endif
    do g=1,ngenerator;do i=1,npoint
      call MPI_Allreduce(generator_maps(i,g),minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(generator_maps(i,g),maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral basin generator maps disagree across ranks';return;endif
    enddo;enddo
    local_bad=merge(0,1,all(generator_maps>=1).and.all(generator_maps<=npoint))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin generator map is out of range';return;endif
    receipt_valid=npoint8<=huge(0_int64)/56_int64
    if(receipt_valid)then;bytes=56_int64*npoint8;else;bytes=0_int64;endif
    if(receipt_valid.and.ngenerator>0)then
      if(npoint8>huge(0_int64)/(4_int64*int(ngenerator,int64)))then
        receipt_valid=.false.
      elseif(bytes>huge(0_int64)-4_int64*npoint8*int(ngenerator,int64))then
        receipt_valid=.false.
      else
        bytes=bytes+4_int64*npoint8*int(ngenerator,int64)
      endif
    endif
    call MPI_Allreduce(merge(0,1,receipt_valid),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral basin workspace extent overflow';return;endif
    allocate(ownership(npoint),parent(npoint),roots(npoint),global_labels(npoint),root_ids(npoint),&
      local_score(npoint),global_score(npoint),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(ownership))deallocate(ownership)
      if(allocated(parent))deallocate(parent)
      if(allocated(roots))deallocate(roots)
      if(allocated(global_labels))deallocate(global_labels)
      if(allocated(root_ids))deallocate(root_ids)
      if(allocated(local_score))deallocate(local_score)
      if(allocated(global_score))deallocate(global_score)
      message='spectral basin workspace allocation failed';return
    endif
    ownership=0
    do p=1,nlocal;ownership(int(row_ids(p)))=ownership(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,npoint,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      deallocate(ownership,parent,roots,global_labels,root_ids,local_score,global_score)
      message='spectral basin rows are not owned exactly once';return
    endif
    local_score=0d0
    local_max=0d0;if(nlocal>0)local_max=maxval(occupied_density)
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(global_max>0d0)then
      do p=1,nlocal;local_score(int(row_ids(p)))=occupied_density(p)/global_max;enddo
    endif
    do j=1,size(empty_moment_density,2)
      local_max=0d0;if(nlocal>0)local_max=maxval(empty_moment_density(:,j))
      call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(global_max>0d0)then
        do p=1,nlocal
          local_score(int(row_ids(p)))=max(local_score(int(row_ids(p))),empty_moment_density(p,j)/global_max)
        enddo
      endif
    enddo
    do j=1,size(shared_density,2)
      local_max=0d0;if(nlocal>0)local_max=maxval(shared_density(:,j))
      call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(global_max>0d0)then
        do p=1,nlocal
          local_score(int(row_ids(p)))=max(local_score(int(row_ids(p))),shared_density(p,j)/global_max)
        enddo
      endif
    enddo
    call MPI_Allreduce(local_score,global_score,npoint,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(global_score)).or.maxval(global_score)<=0d0)then
      deallocate(ownership,parent,roots,global_labels,root_ids,local_score,global_score)
      message='spectral basin score is invalid';return
    endif
    do i=1,npoint
      x=modulo(i-1,grid_shape(1))+1
      y=modulo((i-1)/grid_shape(1),grid_shape(2))+1
      z=(i-1)/(grid_shape(1)*grid_shape(2))+1
      choice=i
      do j=1,6
        select case(j)
        case(1);neighbor=modulo(x,grid_shape(1))+1+grid_shape(1)*(y-1+grid_shape(2)*(z-1))
        case(2);neighbor=modulo(x-2,grid_shape(1))+1+grid_shape(1)*(y-1+grid_shape(2)*(z-1))
        case(3);neighbor=x+grid_shape(1)*(modulo(y,grid_shape(2))+grid_shape(2)*(z-1))
        case(4);neighbor=x+grid_shape(1)*(modulo(y-2,grid_shape(2))+grid_shape(2)*(z-1))
        case(5);neighbor=x+grid_shape(1)*(y-1+grid_shape(2)*modulo(z,grid_shape(3)))
        case default;neighbor=x+grid_shape(1)*(y-1+grid_shape(2)*modulo(z-2,grid_shape(3)))
        end select
        comparison_scale=max(1d0,abs(global_score(choice)),abs(global_score(neighbor)))
        if(global_score(neighbor)>global_score(choice)+tolerance*comparison_scale.or.&
            (abs(global_score(neighbor)-global_score(choice))<=tolerance*comparison_scale.and.neighbor<choice))choice=neighbor
      enddo
      parent(i)=choice
    enddo
    do i=1,npoint
      root=i
      do j=1,npoint
        if(parent(root)==root)exit
        root=parent(root)
      enddo
      if(parent(root)/=root)then
        deallocate(ownership,parent,roots,global_labels,root_ids,local_score,global_score)
        message='spectral basin watershed did not terminate';return
      endif
      roots(i)=root
    enddo
    basin_count=0;global_labels=0
    do i=1,npoint
      if(roots(i)/=i)cycle
      basin_count=basin_count+1;root_ids(basin_count)=i
    enddo
    do i=1,npoint
      do j=1,basin_count
        if(roots(i)==root_ids(j))then;global_labels(i)=j;exit;endif
      enddo
    enddo
    allocate(basin_labels(nlocal),basin_orbit_map(basin_count,ngenerator),target_counts(basin_count),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(basin_labels))deallocate(basin_labels)
      if(allocated(basin_orbit_map))deallocate(basin_orbit_map)
      if(allocated(target_counts))deallocate(target_counts)
      deallocate(ownership,parent,roots,global_labels,root_ids,local_score,global_score)
      message='spectral basin output allocation failed';return
    endif
    do p=1,nlocal;basin_labels(p)=global_labels(int(row_ids(p)));enddo
    do g=1,ngenerator
      basin_orbit_map(:,g)=0
      do i=1,npoint
        target_basin=global_labels(generator_maps(i,g));j=global_labels(i)
        if(basin_orbit_map(j,g)==0)basin_orbit_map(j,g)=target_basin
        if(basin_orbit_map(j,g)/=target_basin)local_bad=1
      enddo
      target_counts=0
      do j=1,basin_count
        if(basin_orbit_map(j,g)>=1.and.basin_orbit_map(j,g)<=basin_count)&
          target_counts(basin_orbit_map(j,g))=target_counts(basin_orbit_map(j,g))+1
      enddo
      if(any(target_counts/=1))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(basin_labels,basin_orbit_map,target_counts,ownership,parent,roots,global_labels,root_ids,&
        local_score,global_score);message='spectral basins are not closed under the generators';return
    endif
    hash_value=ieor(int(npoint,int64),ishftc(int(basin_count,int64),19))
    do i=1,npoint;hash_value=ieor(ishftc(hash_value,7),int(global_labels(i),int64));enddo
    do g=1,ngenerator;do j=1,basin_count
      hash_value=ieor(ishftc(hash_value,9),int(basin_orbit_map(j,g),int64))
    enddo;enddo
    if(hash_value==0_int64)hash_value=1_int64
    fingerprint=hash_value;workspace_peak_bytes=bytes
    deallocate(target_counts,ownership,parent,roots,global_labels,root_ids,local_score,global_score);ok=.true.
#else
    ok=.false.;message='periodic spectral basins require MPI';basin_count=0
    fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_periodic_spectral_basins

  subroutine build_dg_occupied_empty_moment_descriptors(comm,row_ids,global_row_count,state_values,&
      eigenvalues,occupations,maximum_moment,tolerance,occupied_density,empty_moment_density,&
      shared_density,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count,maximum_moment
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::state_values(:,:)
    real(real64),intent(in)::eigenvalues(:),occupations(:),tolerance
    real(real64),allocatable,intent(out)::occupied_density(:),empty_moment_density(:,:),shared_density(:,:)
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,nlocal,i,k,p,first_empty,ierr,local_bad,global_bad,minint,maxint,status
    integer(int64)::bits,minbits,maxbits,base_fingerprint,base_workspace,elements
    real(real64)::energy_origin,energy_scale,x
    real(real64),allocatable::base_weights(:,:),base_empty(:,:),total_empty(:)
    nstate=size(state_values,1);nlocal=size(state_values,2)
    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=0
    if(nstate<1.or.size(eigenvalues)/=nstate.or.size(occupations)/=nstate.or.&
        maximum_moment<0.or.maximum_moment>8)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(occupations)))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied/empty moment contract';return;endif
    call MPI_Allreduce(maximum_moment,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(maximum_moment,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minint/=maxint)then;message='empty moment order disagrees across ranks';return;endif
    do i=1,nstate
      bits=transfer(eigenvalues(i),0_int64)
      call MPI_Allreduce(bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minbits/=maxbits)then;message='moment eigenvalues disagree across ranks';return;endif
    enddo
    first_empty=0;local_bad=0
    do i=2,nstate
      if(eigenvalues(i)<eigenvalues(i-1)-tolerance*max(1d0,abs(eigenvalues(i)),abs(eigenvalues(i-1))))local_bad=1
    enddo
    do i=1,nstate
      if(occupations(i)<=tolerance)then;first_empty=i;exit;endif
    enddo
    if(first_empty<=1)then
      local_bad=1
    elseif(any(occupations(first_empty:nstate)>tolerance))then
      local_bad=1
    elseif(maxval(abs(occupations(1:first_empty-1)-occupations(1)))>10d0*tolerance.or.&
        occupations(1)<=tolerance)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='occupied/empty moment descriptors require a sorted gapped integer-occupation spectrum';return
    endif
    if(int(nlocal,int64)>huge(0_int64)/int(maximum_moment+1,int64))then
      local_bad=1;elements=0_int64
    else
      local_bad=0;elements=int(nlocal,int64)*int(maximum_moment+1,int64)
      if(elements>huge(0_int64)/8_int64)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='empty moment extent overflow';return;endif
    allocate(base_weights(nstate,1),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(base_weights))deallocate(base_weights)
      message='empty moment weight allocation failed';return
    endif
    base_weights=0d0;base_weights(first_empty:nstate,1)=1d0
    call build_dg_spectral_density_descriptors(comm,row_ids,global_row_count,state_values,occupations,&
      base_weights,tolerance,occupied_density,base_empty,total_empty,shared_density,base_fingerprint,&
      base_workspace,ok,message)
    deallocate(base_weights)
    if(.not.ok)return
    allocate(empty_moment_density(nlocal,maximum_moment+1),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(empty_moment_density))deallocate(empty_moment_density)
      deallocate(occupied_density,base_empty,total_empty,shared_density)
      message='empty moment density allocation failed';ok=.false.;return
    endif
    empty_moment_density=0d0;empty_moment_density(:,1)=total_empty
    energy_origin=eigenvalues(first_empty);energy_scale=eigenvalues(nstate)-energy_origin
    if(energy_scale>tolerance*max(1d0,abs(energy_origin),abs(eigenvalues(nstate))))then
      do i=first_empty,nstate
        x=max(0d0,min(1d0,(eigenvalues(i)-energy_origin)/energy_scale))
        do k=1,maximum_moment
          do p=1,nlocal
            empty_moment_density(p,k+1)=empty_moment_density(p,k+1)+x**k*abs(state_values(i,p))**2
          enddo
        enddo
      enddo
    endif
    local_bad=merge(0,1,all(ieee_is_finite(empty_moment_density)).and.all(empty_moment_density>=0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(occupied_density,base_empty,total_empty,shared_density,empty_moment_density)
      message='empty moment density is nonfinite';ok=.false.;return
    endif
    fingerprint=base_fingerprint
    do i=1,nstate;fingerprint=ieor(ishftc(fingerprint,9),transfer(eigenvalues(i),0_int64));enddo
    fingerprint=ieor(fingerprint,int(maximum_moment,int64));if(fingerprint==0_int64)fingerprint=1_int64
    if(base_workspace>huge(0_int64)-8_int64*elements)then
      deallocate(occupied_density,base_empty,total_empty,shared_density,empty_moment_density)
      message='empty moment workspace receipt overflow';ok=.false.;return
    endif
    workspace_peak_bytes=base_workspace+8_int64*elements
    deallocate(base_empty,total_empty);ok=.true.
#else
    ok=.false.;message='occupied/empty moment descriptors require MPI';fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_occupied_empty_moment_descriptors

  subroutine build_dg_spectral_density_descriptors(comm,row_ids,global_row_count,state_values,&
      occupations,window_weights,tolerance,occupied_density,unoccupied_density,total_unoccupied_density,shared_density,&
      fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::state_values(:,:)
    real(real64),intent(in)::occupations(:),window_weights(:,:),tolerance
    real(real64),allocatable,intent(out)::occupied_density(:),unoccupied_density(:,:),&
      total_unoccupied_density(:),shared_density(:,:)
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,nlocal,nwindow,i,j,p,ierr,local_bad,global_bad,minint,maxint,status
    integer,allocatable::ownership_count(:)
    integer(int64)::bits,minbits,maxbits,elements,local_hash,global_hash,row_hash,quantized
    real(real64)::minreal,maxreal,maxcoefficient,global_maxcoefficient,safe_coefficient,denominator
    logical::receipt_valid
    nstate=size(state_values,1);nlocal=size(state_values,2);nwindow=size(window_weights,2)
    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=0
    if(global_row_count<1.or.nstate<1.or.nwindow<1.or.size(row_ids)/=nlocal)then
      local_bad=1
    elseif(size(occupations)/=nstate.or.size(window_weights,1)/=nstate)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(real(state_values))).or..not.all(ieee_is_finite(aimag(state_values))))then
      local_bad=1
    elseif(.not.all(ieee_is_finite(occupations)).or..not.all(ieee_is_finite(window_weights)))then
      local_bad=1
    elseif(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)).or.&
        any(occupations<0d0).or.any(window_weights<0d0))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid spectral-density contract';return;endif
    do i=1,3
      select case(i)
      case(1);local_bad=global_row_count
      case(2);local_bad=nstate
      case default;local_bad=nwindow
      end select
      call MPI_Allreduce(local_bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minint/=maxint)then;message='spectral-density metadata disagrees across ranks';return;endif
    enddo
    call MPI_Allreduce(tolerance,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
      message='spectral-density tolerance disagrees across ranks';return
    endif
    do i=1,nstate
      bits=transfer(occupations(i),0_int64)
      call MPI_Allreduce(bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minbits/=maxbits)then;message='spectral occupations disagree across ranks';return;endif
      do j=1,nwindow
        bits=transfer(window_weights(i,j),0_int64)
        call MPI_Allreduce(bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
        call MPI_Allreduce(bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
        if(minbits/=maxbits)then;message='spectral weights disagree across ranks';return;endif
      enddo
    enddo
    local_bad=merge(0,1,maxval(abs(sum(window_weights,dim=2)-merge(0d0,1d0,occupations>tolerance)))<=&
      10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral weights do not cover the unoccupied states';return;endif
    allocate(ownership_count(global_row_count),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(ownership_count))deallocate(ownership_count)
      message='spectral-density ownership allocation failed';return
    endif
    ownership_count=0
    do p=1,nlocal;ownership_count(int(row_ids(p)))=ownership_count(int(row_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      deallocate(ownership_count);message='spectral-density rows are not owned exactly once';return
    endif
    maxcoefficient=0d0
    if(nlocal>0)maxcoefficient=maxval(abs(state_values))
    call MPI_Allreduce(maxcoefficient,global_maxcoefficient,1,&
      MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    denominator=max(1d0,sum(occupations)+sum(window_weights))
    safe_coefficient=sqrt(huge(1d0)/(16d0*denominator))
    local_bad=merge(0,1,global_maxcoefficient<=safe_coefficient)
    receipt_valid=int(nlocal,int64)<=huge(0_int64)/int(2+2*nwindow,int64)
    if(receipt_valid)then
      elements=int(nlocal,int64)*int(2+2*nwindow,int64)
      receipt_valid=elements<=huge(0_int64)/8_int64
    else
      elements=0_int64
    endif
    if(.not.receipt_valid)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(ownership_count);message='spectral-density magnitude or extent is unsafe';return
    endif
    allocate(occupied_density(nlocal),unoccupied_density(nlocal,nwindow),total_unoccupied_density(nlocal),&
      shared_density(nlocal,nwindow),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(occupied_density))deallocate(occupied_density)
      if(allocated(unoccupied_density))deallocate(unoccupied_density)
      if(allocated(total_unoccupied_density))deallocate(total_unoccupied_density)
      if(allocated(shared_density))deallocate(shared_density)
      deallocate(ownership_count);message='spectral-density allocation failed';return
    endif
    occupied_density=0d0;unoccupied_density=0d0
    do p=1,nlocal
      do i=1,nstate
        occupied_density(p)=occupied_density(p)+occupations(i)*abs(state_values(i,p))**2
        do j=1,nwindow
          unoccupied_density(p,j)=unoccupied_density(p,j)+window_weights(i,j)*abs(state_values(i,p))**2
        enddo
      enddo
    enddo
    total_unoccupied_density=sum(unoccupied_density,dim=2)
    do j=1,nwindow
      do p=1,nlocal
        denominator=occupied_density(p)+unoccupied_density(p,j)
        if(denominator>tiny(1d0))then
          shared_density(p,j)=2d0*sqrt(occupied_density(p)*unoccupied_density(p,j))/denominator
        else
          shared_density(p,j)=0d0
        endif
      enddo
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(occupied_density)).and.&
      all(ieee_is_finite(unoccupied_density)).and.all(ieee_is_finite(total_unoccupied_density)).and.&
      all(ieee_is_finite(shared_density)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(occupied_density,unoccupied_density,total_unoccupied_density,shared_density,ownership_count)
      message='spectral-density accumulation is nonfinite';return
    endif
    local_hash=0_int64;local_bad=0
    do p=1,nlocal
      row_hash=ieor(row_ids(p),ishftc(row_ids(p),23))
      denominator=100d0*tolerance
      if(occupied_density(p)>0.25d0*real(huge(0_int64),real64)*denominator)local_bad=1
      if(local_bad==0)then
        quantized=nint(occupied_density(p)/denominator,int64);row_hash=ieor(ishftc(row_hash,7),quantized)
        do j=1,nwindow
          if(unoccupied_density(p,j)>0.25d0*real(huge(0_int64),real64)*denominator)then
            local_bad=1;exit
          endif
          quantized=nint(unoccupied_density(p,j)/denominator,int64)
          row_hash=ieor(ishftc(row_hash,7),quantized)
          quantized=nint(shared_density(p,j)/denominator,int64)
          row_hash=ieor(ishftc(row_hash,7),quantized)
        enddo
      endif
      local_hash=ieor(local_hash,row_hash)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(occupied_density,unoccupied_density,total_unoccupied_density,shared_density,ownership_count)
      message='spectral-density fingerprint range is unsafe';return
    endif
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      deallocate(occupied_density,unoccupied_density,total_unoccupied_density,shared_density,ownership_count);return
    endif
    if(global_hash==0_int64)global_hash=1_int64
    fingerprint=global_hash;workspace_peak_bytes=8_int64*elements+4_int64*int(global_row_count,int64)
    deallocate(ownership_count);ok=.true.
#else
    ok=.false.;message='spectral density descriptors require MPI';fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_spectral_density_descriptors

  subroutine build_dg_equal_count_spectral_windows(comm,eigenvalues,occupations,nwindow,tolerance,&
      window_weights,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,nwindow
    real(real64),intent(in)::eigenvalues(:),occupations(:),tolerance
    real(real64),allocatable,intent(out)::window_weights(:,:)
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,nunoccupied,first_unoccupied,i,j,k,cluster_last,cluster_size,window,&
      ierr,local_bad,global_bad,status,minint,maxint
    integer(int64)::local_bits,minbits,maxbits,nelement,hash_value
    real(real64)::scale,minreal,maxreal,cluster_coordinate,window_coordinate,fraction
    nstate=size(eigenvalues);ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=0
    if(nstate<1.or.size(occupations)/=nstate.or.nwindow<1)then
      local_bad=1
    elseif(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      local_bad=1
    elseif(.not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(occupations)))then
      local_bad=1
    elseif(any(occupations<0d0))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid spectral-window contract';return;endif
    call MPI_Allreduce(nstate,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nstate,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minint/=maxint)then;message='spectral state count disagrees across ranks';return;endif
    call MPI_Allreduce(nwindow,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nwindow,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minint/=maxint)then;message='spectral window count disagrees across ranks';return;endif
    call MPI_Allreduce(tolerance,minreal,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxreal,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(transfer(minreal,0_int64)/=transfer(maxreal,0_int64))then
      message='spectral tolerance disagrees across ranks';return
    endif
    do i=1,nstate
      local_bits=transfer(eigenvalues(i),0_int64)
      call MPI_Allreduce(local_bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minbits/=maxbits)then;message='spectral eigenvalues disagree across ranks';return;endif
      local_bits=transfer(occupations(i),0_int64)
      call MPI_Allreduce(local_bits,minbits,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(local_bits,maxbits,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(minbits/=maxbits)then;message='spectral occupations disagree across ranks';return;endif
    enddo
    local_bad=0
    do i=2,nstate
      scale=max(1d0,abs(eigenvalues(i-1)),abs(eigenvalues(i)))
      if(eigenvalues(i)<eigenvalues(i-1)-tolerance*scale)local_bad=1
    enddo
    first_unoccupied=0
    do i=1,nstate
      if(occupations(i)<=tolerance)then;first_unoccupied=i;exit;endif
    enddo
    if(first_unoccupied==0)local_bad=1
    if(first_unoccupied>0)then
      if(any(occupations(first_unoccupied:nstate)>tolerance))local_bad=1
      nunoccupied=nstate-first_unoccupied+1
      if(nwindow>nunoccupied)local_bad=1
    else
      nunoccupied=0
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid ordered occupied/unoccupied spectrum';return;endif
    if(int(nstate,int64)>huge(0_int64)/int(nwindow,int64))then
      local_bad=1;nelement=0_int64
    else
      local_bad=0;nelement=int(nstate,int64)*int(nwindow,int64)
      if(nelement>huge(0_int64)/8_int64)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spectral-window extent overflow';return;endif
    allocate(window_weights(nstate,nwindow),stat=status)
    call MPI_Allreduce(merge(0,1,status==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(window_weights))deallocate(window_weights)
      message='spectral-window allocation failed';return
    endif
    window_weights=0d0;i=first_unoccupied
    do while(i<=nstate)
      cluster_last=i
      do while(cluster_last<nstate)
        scale=max(1d0,abs(eigenvalues(cluster_last)),abs(eigenvalues(cluster_last+1)))
        if(abs(eigenvalues(cluster_last+1)-eigenvalues(cluster_last))>tolerance*scale)exit
        cluster_last=cluster_last+1
      enddo
      cluster_size=cluster_last-i+1
      j=i-first_unoccupied
      cluster_coordinate=(real(j,real64)+0.5d0*real(cluster_size,real64))/real(nunoccupied,real64)
      window_coordinate=cluster_coordinate*real(nwindow,real64)+0.5d0
      if(window_coordinate<=1d0)then
        window_weights(i:cluster_last,1)=1d0
      elseif(window_coordinate>=real(nwindow,real64))then
        window_weights(i:cluster_last,nwindow)=1d0
      else
        window=int(floor(window_coordinate));fraction=window_coordinate-real(window,real64)
        ! Smoothstep keeps a tolerance-degenerate cluster intact while avoiding
        ! a discontinuous transfer at an equal-count window boundary.
        fraction=fraction*fraction*(3d0-2d0*fraction)
        window_weights(i:cluster_last,window)=1d0-fraction
        window_weights(i:cluster_last,window+1)=fraction
      endif
      i=cluster_last+1
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(window_weights)).and.all(window_weights>=0d0).and.&
      maxval(abs(sum(window_weights(first_unoccupied:nstate,:),dim=2)-1d0))<=10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      deallocate(window_weights);message='spectral-window partition failed';return
    endif
    hash_value=ieor(int(nstate,int64),ishftc(int(nwindow,int64),17))
    do i=1,nstate
      hash_value=ieor(ishftc(hash_value,7),transfer(eigenvalues(i),0_int64))
      hash_value=ieor(ishftc(hash_value,11),transfer(occupations(i),0_int64))
      do k=1,nwindow
        hash_value=ieor(ishftc(hash_value,5),int(nint(window_weights(i,k)*1024d0),int64))
      enddo
    enddo
    if(hash_value==0_int64)hash_value=1_int64
    fingerprint=hash_value;workspace_peak_bytes=8_int64*nelement;ok=.true.
#else
    ok=.false.;message='spectral windows require MPI';fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_equal_count_spectral_windows

  subroutine materialize_dg_row_owned_sector_on_spatial_grid(comm,row_ids,global_state_count,&
      sector_rows,local_basis,basis_fingerprint,spatial_sector,output_fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_state_count
    integer(int64),intent(in)::row_ids(:)
    integer(int64),intent(in)::basis_fingerprint
    complex(real64),intent(in)::sector_rows(:,:),local_basis(:,:)
    complex(real64),allocatable,intent(out)::spatial_sector(:,:)
    integer(int64),intent(out)::output_fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nrow,m,npoint,i,p,rank,ierr,local_bad,global_bad,status,minint,maxint
    integer,allocatable::owner(:),position(:),ownership_count(:)
    complex(real64),allocatable::stream_row(:)
    integer(int64)::elements,bytes,term,minhash,maxhash,bits
    real(real64)::local_magnitude,global_magnitude,safe_magnitude
    logical::receipt_valid
    nrow=size(row_ids);m=size(sector_rows,2);npoint=size(local_basis,2)
    ok=.false.;message='';output_fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=merge(0,1,global_state_count>=1.and.m>=1.and.npoint>=0.and.&
      all(shape(sector_rows)==[nrow,m]).and.size(local_basis,1)==global_state_count.and.&
      basis_fingerprint/=0_int64.and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_state_count,int64)).and.&
      all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
      all(ieee_is_finite(real(local_basis))).and.all(ieee_is_finite(aimag(local_basis))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid row-owned spatial sector contract';return;endif
    call MPI_Allreduce(global_state_count,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(global_state_count,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='spatial sector state extent disagrees';return;endif
    call MPI_Allreduce(m,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='spatial sector rank disagrees';return;endif
    call MPI_Allreduce(basis_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(basis_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='spatial sector basis provenance disagrees';return;endif
    local_bad=merge(1,0,int(npoint,int64)>huge(0_int64)/int(m,int64))
    elements=0_int64;bytes=0_int64;receipt_valid=local_bad==0
    if(receipt_valid)then
      elements=int(npoint,int64)*int(m,int64)
      if(elements>huge(0_int64)-int(m,int64))then
        receipt_valid=.false.
      else
        elements=elements+int(m,int64)
      endif
      if(elements>huge(0_int64)/16_int64)receipt_valid=.false.
      if(receipt_valid)bytes=16_int64*elements
      if(int(global_state_count,int64)>huge(0_int64)/12_int64)receipt_valid=.false.
      if(receipt_valid)then
        term=12_int64*int(global_state_count,int64)
        if(bytes>huge(0_int64)-term)then
          receipt_valid=.false.
        else
          bytes=bytes+term
        endif
      endif
    endif
    local_bad=merge(0,1,receipt_valid)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spatial sector workspace overflows';return;endif
    workspace_peak_bytes=bytes
    local_magnitude=maxval(abs(local_basis))
    if(nrow>0)local_magnitude=max(local_magnitude,maxval(abs(sector_rows)))
    call MPI_Allreduce(local_magnitude,global_magnitude,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_magnitude=sqrt(huge(1d0))/(4d0*real(global_state_count,real64))
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.global_magnitude<=safe_magnitude)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='spatial sector input magnitude is unsafe';return;endif
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(spatial_sector(npoint,m),stream_row(m),owner(global_state_count),position(global_state_count),&
      ownership_count(global_state_count),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(spatial_sector))deallocate(spatial_sector)
      message='spatial sector allocation failed';return
    endif
    owner=0;position=0;ownership_count=0
    do i=1,nrow
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;ownership_count(int(row_ids(i)))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_state_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,position,global_state_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_state_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;message='spatial sector coefficient rows are not uniquely owned';return;endif
    spatial_sector=(0d0,0d0)
    do i=1,global_state_count
      stream_row=(0d0,0d0)
      if(rank==owner(i)-1)stream_row=sector_rows(position(i),:)
      call MPI_Bcast(stream_row,m,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      do p=1,npoint;spatial_sector(p,:)=spatial_sector(p,:)+local_basis(i,p)*stream_row;enddo
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(real(spatial_sector))).and.all(ieee_is_finite(aimag(spatial_sector))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='materialized spatial sector is nonfinite';return;endif
    output_fingerprint=ieor(basis_fingerprint,int(global_state_count,int64))
    do i=1,global_state_count
      stream_row=(0d0,0d0)
      if(rank==owner(i)-1)stream_row=sector_rows(position(i),:)
      call MPI_Bcast(stream_row,m,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      output_fingerprint=ieor(ishftc(output_fingerprint,9),int(i,int64))
      do p=1,m
        bits=transfer(real(stream_row(p),real64),bits);output_fingerprint=ieor(ishftc(output_fingerprint,9),bits)
        bits=transfer(aimag(stream_row(p)),bits);output_fingerprint=ieor(ishftc(output_fingerprint,9),bits)
      enddo
    enddo
    if(output_fingerprint==0_int64)output_fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='row-owned spatial sector materialization requires MPI';output_fingerprint=0_int64
    workspace_peak_bytes=0_int64
#endif
  end subroutine materialize_dg_row_owned_sector_on_spatial_grid

  subroutine build_dg_translation_character_intertwining_phase(comm,row_ids,global_row_count,&
      generator_maps,generator_orders,element_words,product_table,identity_operation,reference_character,target_character,&
      catalog_fingerprint,tolerance,&
      local_phase,fingerprint,phase_payload_fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count,generator_orders(:),element_words(:,:),product_table(:,:),identity_operation
    integer(int64),intent(in)::row_ids(:),generator_maps(:,:)
    complex(real64),intent(in)::reference_character(:),target_character(:)
    integer(int64),intent(in)::catalog_fingerprint
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::local_phase(:)
    integer(int64),intent(out)::fingerprint,phase_payload_fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nt,ng,nlocal,x,g,h,target,ierr,local_bad,global_bad,allocation_status,minint,maxint
    integer,allocatable::permutation_count(:),ownership_count(:),generator_map(:),element_map(:),right_map(:),product_map(:)
    complex(real64),allocatable::global_phase(:),ratio(:)
    logical,allocatable::assigned(:)
    integer(int64)::bits,metadata_hash,minhash,maxhash,elements
    real(real64)::mintol,maxtol
    nt=size(reference_character);ng=size(generator_maps,2);nlocal=size(row_ids);ok=.false.;message=''
    fingerprint=0_int64;phase_payload_fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=merge(0,1,global_row_count>=1.and.nt>=1.and.size(target_character)==nt.and.&
      ng>=1.and.size(generator_orders)==ng.and.all(generator_orders>=1).and.&
      all(generator_orders<=nt).and.all(mod(nt,generator_orders)==0).and.&
      size(generator_maps,1)==nlocal.and.all(shape(element_words)==[nt,ng]).and.&
      all(shape(product_table)==[nt,nt]).and.all(element_words>=0).and.&
      identity_operation>=1.and.identity_operation<=nt.and.&
      all(generator_maps>=1_int64).and.all(generator_maps<=int(global_row_count,int64)).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.catalog_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(reference_character))).and.all(ieee_is_finite(aimag(reference_character))).and.&
      all(ieee_is_finite(real(target_character))).and.all(ieee_is_finite(aimag(target_character))).and.&
      maxval(abs(abs(reference_character)-1d0))<=10d0*tolerance.and.&
      maxval(abs(abs(target_character)-1d0))<=10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid translation intertwining-phase contract';return;endif
    call MPI_Allreduce(global_row_count,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(global_row_count,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='translation phase row extent disagrees';return;endif
    call MPI_Allreduce(nt,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nt,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='translation phase order disagrees';return;endif
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='translation phase tolerance disagrees';return;endif
    do h=1,ng
      local_bad=merge(1,0,any(element_words(:,h)>=generator_orders(h)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='translation generator words exceed their orders';return;endif
    enddo
    do h=1,ng;call agree_integer(generator_orders(h));if(global_bad/=0)return;enddo
    do g=1,nt
      do h=1,ng;call agree_integer(element_words(g,h));if(global_bad/=0)return;enddo
      do h=1,nt;call agree_integer(product_table(h,g));if(global_bad/=0)return;enddo
      call agree_complex(reference_character(g));if(global_bad/=0)return
      call agree_complex(target_character(g));if(global_bad/=0)return
    enddo
    metadata_hash=catalog_fingerprint
    metadata_hash=ieor(ishftc(metadata_hash,7),int(identity_operation,int64))
    do g=1,ng
      metadata_hash=ieor(ishftc(metadata_hash,7),int(generator_orders(g),int64))
    enddo
    do g=1,nt
      bits=transfer(real(reference_character(g),real64),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(aimag(reference_character(g)),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(real(target_character(g),real64),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(aimag(target_character(g)),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      do h=1,ng;metadata_hash=ieor(ishftc(metadata_hash,7),int(element_words(g,h),int64));enddo
      do h=1,nt;metadata_hash=ieor(ishftc(metadata_hash,7),int(product_table(h,g),int64));enddo
    enddo
    allocate(generator_map(global_row_count),element_map(global_row_count),right_map(global_row_count),&
      product_map(global_row_count),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='translation map streaming allocation failed';return;endif
    do g=1,ng
      call gather_generator_map(g,generator_map);if(global_bad/=0)return
      do x=1,global_row_count;metadata_hash=ieor(ishftc(metadata_hash,7),int(generator_map(x),int64));enddo
    enddo
    call MPI_Allreduce(metadata_hash,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(metadata_hash,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='translation phase metadata disagree';return;endif
    local_bad=merge(1,0,int(global_row_count,int64)>huge(0_int64)/64_int64.or.&
        int(nlocal,int64)>huge(0_int64)/16_int64.or.int(nt,int64)>huge(0_int64)/16_int64)
    if(local_bad==0)then
      elements=64_int64*int(global_row_count,int64)
      if(elements>huge(0_int64)-16_int64*int(nlocal,int64)-16_int64*int(nt,int64))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='translation phase workspace receipt overflows';return
    endif
    workspace_peak_bytes=elements+16_int64*int(nlocal,int64)+16_int64*int(nt,int64)
    allocate(global_phase(global_row_count),ratio(nt),assigned(global_row_count),local_phase(nlocal),&
      permutation_count(global_row_count),ownership_count(global_row_count),&
      stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(local_phase))deallocate(local_phase)
      message='translation intertwining-phase allocation failed';return
    endif
    ownership_count=0;do x=1,nlocal;ownership_count(int(row_ids(x)))=ownership_count(int(row_ids(x)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;message='translation phase rows are not uniquely owned';return;endif
    local_bad=merge(1,0,any(product_table<1).or.any(product_table>nt))
    if(local_bad==0)then
      do g=1,nt
        permutation_count=0
        call build_element_map(g,element_map);call build_element_map(identity_operation,product_map)
        do x=1,global_row_count
          target=element_map(x);permutation_count(target)=permutation_count(target)+1
          if(product_map(x)/=x)local_bad=1
        enddo
        if(any(permutation_count/=1))local_bad=1
      enddo
      do g=1,nt;do h=1,nt
        call build_element_map(g,element_map);call build_element_map(h,right_map)
        call build_element_map(product_table(h,g),product_map)
        do x=1,global_row_count
          if(element_map(right_map(x))/=product_map(x))local_bad=1
        enddo
      enddo;enddo
      do g=1,ng
        call gather_generator_map(g,generator_map)
        do x=1,global_row_count
          target=x
          do h=1,generator_orders(g);target=generator_map(target);enddo
          if(target/=x)local_bad=1
        enddo
      enddo
      do g=1,nt
        call build_element_map(g,element_map)
        call build_element_map(identity_operation,product_map)
        if(g/=identity_operation.and.any(element_map==product_map))local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='translation maps do not realize the supplied finite-group action';return
    endif
    ratio=target_character*conjg(reference_character);global_phase=(0d0,0d0);assigned=.false.
    do x=1,global_row_count;permutation_count(x)=x;enddo
    do g=1,nt
      call build_element_map(g,element_map)
      permutation_count=min(permutation_count,element_map)
    enddo
    do x=1,global_row_count
      if(permutation_count(x)/=x)cycle
      do g=1,nt
        call build_element_map(g,element_map);target=element_map(x)
        if(assigned(target))then
          if(abs(global_phase(target)-ratio(g))>10d0*tolerance)local_bad=1
        else
          global_phase(target)=ratio(g);assigned(target)=.true.
        endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='translation phase orbit is inconsistent with character ratio';return
    endif
    local_bad=0
    do x=1,global_row_count;do g=1,nt
      call build_element_map(g,element_map);target=element_map(x)
      if(abs(global_phase(target)-ratio(g)*global_phase(x))>10d0*tolerance)local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='translation phase fails the complete character covariance action';return
    endif
    do x=1,nlocal;local_phase(x)=global_phase(int(row_ids(x)));enddo
    fingerprint=metadata_hash;phase_payload_fingerprint=int(z'243F6A8885A308D3',int64)
    do x=1,global_row_count
      bits=transfer(real(global_phase(x),real64),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
      phase_payload_fingerprint=ieor(ishftc(phase_payload_fingerprint,11),bits)
      bits=transfer(aimag(global_phase(x)),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
      phase_payload_fingerprint=ieor(ishftc(phase_payload_fingerprint,11),bits)
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    if(phase_payload_fingerprint==0_int64)phase_payload_fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='translation intertwining phase requires MPI';fingerprint=0_int64
    phase_payload_fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_integer(value)
      integer,intent(in)::value
      call MPI_Allreduce(value,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;global_bad=1;return;endif
      call MPI_Allreduce(value,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minint/=maxint)
      if(global_bad/=0)message='translation phase integer metadata disagree'
    end subroutine
    subroutine agree_int64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;global_bad=1;return;endif
      call MPI_Allreduce(value,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minhash/=maxhash)
      if(global_bad/=0)message='translation phase map metadata disagree'
    end subroutine
    subroutine agree_complex(value)
      complex(real64),intent(in)::value
      bits=transfer(real(value,real64),bits);call agree_int64(bits);if(global_bad/=0)return
      bits=transfer(aimag(value),bits);call agree_int64(bits)
    end subroutine
    subroutine gather_generator_map(generator,map)
      integer,intent(in)::generator
      integer,intent(out)::map(:)
      integer::index
      map=0
      do index=1,nlocal;map(int(row_ids(index)))=int(generator_maps(index,generator));enddo
      call MPI_Allreduce(MPI_IN_PLACE,map,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.any(map<1).or.any(map>global_row_count))
      if(global_bad/=0)message='translation generator map stream is invalid'
    end subroutine
    subroutine build_element_map(element,map)
      integer,intent(in)::element
      integer,intent(out)::map(:)
      integer::generator,power,index
      do index=1,global_row_count;map(index)=index;enddo
      do generator=1,ng
        call gather_generator_map(generator,generator_map);if(global_bad/=0)return
        do power=1,element_words(element,generator);map=generator_map(map);enddo
      enddo
    end subroutine
#endif
  end subroutine build_dg_translation_character_intertwining_phase

  subroutine release_dg_prepared_translation_action(action)
    type(s_dg_prepared_translation_action),intent(inout)::action
    if(allocated(action%row_ids))deallocate(action%row_ids)
    if(allocated(action%generator_orders))deallocate(action%generator_orders)
    if(allocated(action%element_words))deallocate(action%element_words)
    if(allocated(action%product_table))deallocate(action%product_table)
    if(allocated(action%generator_maps))deallocate(action%generator_maps)
    if(allocated(action%element_maps))deallocate(action%element_maps)
    action%global_row_count=0;action%ntranslation=0;action%ngenerator=0;action%identity_operation=0
    action%construction_collective_count=0;action%catalog_fingerprint=0_int64
    action%workspace_peak_bytes=0_int64
  end subroutine release_dg_prepared_translation_action

  subroutine prepare_dg_translation_character_action(comm,row_ids,global_row_count,generator_maps,&
      generator_orders,element_words,product_table,identity_operation,catalog_fingerprint,tolerance,action,ok,message)
    integer,intent(in)::comm,global_row_count,generator_orders(:),element_words(:,:),product_table(:,:),identity_operation
    integer(int64),intent(in)::row_ids(:),generator_maps(:,:)
    integer(int64),intent(in)::catalog_fingerprint
    real(real64),intent(in)::tolerance
    type(s_dg_prepared_translation_action),intent(inout)::action
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nt,ng,nlocal,status,ierr,local_bad,global_bad,g,h,k,x,p,power,target,metadata_count,index,&
      minint,maxint
    integer,allocatable::ownership_count(:),permutation_count(:),metadata(:),metadata_minimum(:),metadata_maximum(:)
    integer(int64)::integer_elements,term,minhash,maxhash
    real(real64)::mintol,maxtol
    call release_dg_prepared_translation_action(action)
    nt=size(element_words,1);ng=size(generator_maps,2);nlocal=size(row_ids)
    ok=.false.;message=''
    local_bad=merge(0,1,global_row_count>=1.and.nt>=1.and.ng>=1.and.nlocal>=1.and.&
      size(generator_maps,1)==nlocal.and.size(generator_orders)==ng.and.&
      size(element_words,2)==ng.and.all(shape(product_table)==[nt,nt]).and.&
      identity_operation>=1.and.identity_operation<=nt.and.catalog_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      all(generator_maps>=1_int64).and.all(generator_maps<=int(global_row_count,int64)).and.&
      all(generator_orders>=1).and.all(generator_orders<=nt).and.all(mod(nt,generator_orders)==0).and.&
      all(element_words>=0).and.all(product_table>=1).and.all(product_table<=nt))
    if(local_bad==0)then
      do g=1,ng;if(any(element_words(:,g)>=generator_orders(g)))local_bad=1;enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid prepared translation action contract';return;endif
    call agree_integer(global_row_count);if(global_bad/=0)return
    call agree_integer(nt);if(global_bad/=0)return
    call agree_integer(ng);if(global_bad/=0)return
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='prepared translation tolerance disagrees';return;endif
    call MPI_Allreduce(catalog_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(catalog_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='prepared translation catalog disagrees';return;endif
    local_bad=0;integer_elements=0_int64
    if(int(nt,int64)>huge(0_int64)/(int(ng,int64)+int(nt,int64)))local_bad=1
    if(local_bad==0)then
      term=int(ng,int64)+int(nt,int64)*int(ng+nt,int64)+1_int64
      if(term>int(huge(0),int64))local_bad=1
    endif
    if(local_bad==0)then
      metadata_count=int(term)
      if(int(global_row_count,int64)>huge(0_int64)/(int(nt,int64)+int(ng,int64)+2_int64))local_bad=1
    endif
    if(local_bad==0)integer_elements=int(global_row_count,int64)*(int(nt,int64)+int(ng,int64)+2_int64)
    if(local_bad==0)then
      term=int(nlocal,int64)+int(ng,int64)+int(nt,int64)*int(ng,int64)+&
        int(nt,int64)*int(nt,int64)+3_int64*int(metadata_count,int64)
      if(integer_elements>huge(0_int64)-term)local_bad=1
    endif
    if(local_bad==0.and.integer_elements+term>huge(0_int64)/4_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared translation action workspace overflows';ok=.false.;return;endif
    allocate(action%row_ids(nlocal),action%generator_orders(ng),action%element_words(nt,ng),&
      action%product_table(nt,nt),action%generator_maps(global_row_count,ng),&
      action%element_maps(global_row_count,nt),ownership_count(global_row_count),&
      permutation_count(global_row_count),metadata(metadata_count),metadata_minimum(metadata_count),&
      metadata_maximum(metadata_count),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_prepared_translation_action(action);call cleanup_temporary()
      message='prepared translation action allocation failed';ok=.false.;return
    endif
    ownership_count=0
    do x=1,nlocal;ownership_count(int(row_ids(x)))=ownership_count(int(row_ids(x)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      call release_dg_prepared_translation_action(action);call cleanup_temporary()
      message='prepared translation rows are not uniquely owned';return
    endif
    index=1;metadata(index)=identity_operation
    do g=1,ng;index=index+1;metadata(index)=generator_orders(g);enddo
    do g=1,ng;do h=1,nt;index=index+1;metadata(index)=element_words(h,g);enddo;enddo
    do g=1,nt;do h=1,nt;index=index+1;metadata(index)=product_table(h,g);enddo;enddo
    call MPI_Allreduce(metadata,metadata_minimum,metadata_count,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call release_dg_prepared_translation_action(action);call cleanup_temporary();return;endif
    call MPI_Allreduce(metadata,metadata_maximum,metadata_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(metadata_minimum/=metadata_maximum))then
      call release_dg_prepared_translation_action(action);call cleanup_temporary()
      message='prepared translation integer metadata disagree';return
    endif
    action%row_ids=row_ids;action%generator_orders=generator_orders
    action%element_words=element_words;action%product_table=product_table
    action%global_row_count=global_row_count;action%ntranslation=nt;action%ngenerator=ng
    action%identity_operation=identity_operation;action%catalog_fingerprint=catalog_fingerprint
    do g=1,ng
      action%generator_maps(:,g)=0
      do x=1,nlocal
        action%generator_maps(int(row_ids(x)),g)=int(generator_maps(x,g))
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,action%generator_maps(:,g),global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then
        call release_dg_prepared_translation_action(action);call cleanup_temporary()
        message='prepared translation generator gather failed';ok=.false.;return
      endif
    enddo
    do g=1,nt
      do x=1,global_row_count;action%element_maps(x,g)=x;enddo
      do p=1,ng
        do power=1,element_words(g,p)
          action%element_maps(:,g)=action%generator_maps(action%element_maps(:,g),p)
        enddo
      enddo
    enddo
    local_bad=0
    do g=1,ng
      permutation_count=0
      do x=1,global_row_count
        target=action%generator_maps(x,g)
        if(target<1.or.target>global_row_count)then;local_bad=1;cycle;endif
        permutation_count(target)=permutation_count(target)+1
      enddo
      if(any(permutation_count/=1))local_bad=1
      do x=1,global_row_count;target=x
        do power=1,generator_orders(g);target=action%generator_maps(target,g);enddo
        if(target/=x)local_bad=1
      enddo
    enddo
    do g=1,nt
      permutation_count=0
      do x=1,global_row_count
        target=action%element_maps(x,g);permutation_count(target)=permutation_count(target)+1
      enddo
      if(any(permutation_count/=1))local_bad=1
      if(g==identity_operation)then
        do x=1,global_row_count;if(action%element_maps(x,g)/=x)local_bad=1;enddo
      else
        do x=1,global_row_count;if(action%element_maps(x,g)==x)local_bad=1;enddo
      endif
    enddo
    do g=1,nt;do h=1,nt
      k=product_table(h,g)
      do x=1,global_row_count
        if(action%element_maps(action%element_maps(x,h),g)/=action%element_maps(x,k))local_bad=1
      enddo
    enddo;enddo
    do g=1,nt;do h=1,nt;do k=1,nt
      if(product_table(product_table(g,h),k)/=product_table(g,product_table(h,k)))local_bad=1
    enddo;enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call cleanup_temporary()
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_prepared_translation_action(action)
      message='prepared translation maps do not realize the supplied free group action';return
    endif
    action%workspace_peak_bytes=4_int64*(integer_elements+term)
    action%construction_collective_count=ng+10
    ok=.true.
#else
    call release_dg_prepared_translation_action(action)
    ok=.false.;message='prepared translation action requires MPI'
#endif
  contains
#ifdef USE_MPI
    subroutine agree_integer(value)
      integer,intent(in)::value
      call MPI_Allreduce(value,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;global_bad=1;return;endif
      call MPI_Allreduce(value,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minint/=maxint)
      if(global_bad/=0)message='prepared translation dimensions disagree'
    end subroutine
    subroutine cleanup_temporary()
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(permutation_count))deallocate(permutation_count)
      if(allocated(metadata))deallocate(metadata)
      if(allocated(metadata_minimum))deallocate(metadata_minimum)
      if(allocated(metadata_maximum))deallocate(metadata_maximum)
    end subroutine
#endif
  end subroutine prepare_dg_translation_character_action

  subroutine build_dg_translation_character_intertwining_phase_prepared(comm,action,reference_character,&
      target_character,tolerance,local_phase,fingerprint,phase_payload_fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    type(s_dg_prepared_translation_action),intent(in)::action
    complex(real64),intent(in)::reference_character(:),target_character(:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::local_phase(:)
    integer(int64),intent(out)::fingerprint,phase_payload_fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nt,ng,nlocal,status,ierr,local_bad,global_bad,minint,maxint,g,h,x,target
    integer,allocatable::orbit_representative(:)
    complex(real64),allocatable::global_phase(:),ratio(:)
    logical,allocatable::assigned(:)
    integer(int64)::bits,metadata_hash,minhash,maxhash,elements
    real(real64)::mintol,maxtol
    nt=action%ntranslation;ng=action%ngenerator
    if(allocated(action%row_ids))then;nlocal=size(action%row_ids);else;nlocal=0;endif
    ok=.false.;message='';fingerprint=0_int64;phase_payload_fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=merge(0,1,action%global_row_count>=1.and.nt>=1.and.ng>=1.and.nlocal>=1.and.&
      action%catalog_fingerprint/=0_int64.and.size(reference_character)==nt.and.size(target_character)==nt.and.&
      allocated(action%generator_orders).and.allocated(action%element_words).and.allocated(action%product_table).and.&
      allocated(action%generator_maps).and.allocated(action%element_maps).and.&
      all(shape(action%element_maps)==[action%global_row_count,nt]).and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(reference_character))).and.all(ieee_is_finite(aimag(reference_character))).and.&
      all(ieee_is_finite(real(target_character))).and.all(ieee_is_finite(aimag(target_character))).and.&
      maxval(abs(abs(reference_character)-1d0))<=10d0*tolerance.and.&
      maxval(abs(abs(target_character)-1d0))<=10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid prepared translation phase contract';return;endif
    call MPI_Allreduce(nt,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nt,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='prepared translation order disagrees';return;endif
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='prepared translation tolerance disagrees';return;endif
    metadata_hash=action%catalog_fingerprint
    metadata_hash=ieor(ishftc(metadata_hash,7),int(action%identity_operation,int64))
    do g=1,ng;metadata_hash=ieor(ishftc(metadata_hash,7),int(action%generator_orders(g),int64));enddo
    do g=1,nt
      call agree_complex(reference_character(g));if(global_bad/=0)return
      call agree_complex(target_character(g));if(global_bad/=0)return
      bits=transfer(real(reference_character(g),real64),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(aimag(reference_character(g)),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(real(target_character(g),real64),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      bits=transfer(aimag(target_character(g)),bits);metadata_hash=ieor(ishftc(metadata_hash,7),bits)
      do h=1,ng;metadata_hash=ieor(ishftc(metadata_hash,7),int(action%element_words(g,h),int64));enddo
      do h=1,nt;metadata_hash=ieor(ishftc(metadata_hash,7),int(action%product_table(h,g),int64));enddo
    enddo
    do g=1,ng;do x=1,action%global_row_count
      metadata_hash=ieor(ishftc(metadata_hash,7),int(action%generator_maps(x,g),int64))
    enddo;enddo
    call MPI_Allreduce(metadata_hash,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(metadata_hash,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='prepared translation metadata disagree';return;endif
    local_bad=merge(1,0,int(action%global_row_count,int64)>huge(0_int64)/29_int64.or.&
      int(nlocal+nt,int64)>huge(0_int64)/16_int64)
    if(local_bad==0)then
      elements=29_int64*int(action%global_row_count,int64)
      if(elements>huge(0_int64)-16_int64*int(nlocal+nt,int64))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared translation phase workspace overflows';return;endif
    workspace_peak_bytes=elements+16_int64*int(nlocal+nt,int64)
    allocate(global_phase(action%global_row_count),ratio(nt),assigned(action%global_row_count),&
      orbit_representative(action%global_row_count),local_phase(nlocal),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(local_phase))deallocate(local_phase)
      message='prepared translation phase allocation failed';return
    endif
    ratio=target_character*conjg(reference_character);global_phase=(0d0,0d0);assigned=.false.
    do x=1,action%global_row_count;orbit_representative(x)=x;enddo
    do g=1,nt;orbit_representative=min(orbit_representative,action%element_maps(:,g));enddo
    local_bad=0
    do x=1,action%global_row_count
      if(orbit_representative(x)/=x)cycle
      do g=1,nt
        target=action%element_maps(x,g)
        if(assigned(target))then
          if(abs(global_phase(target)-ratio(g))>10d0*tolerance)local_bad=1
        else
          global_phase(target)=ratio(g);assigned(target)=.true.
        endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared translation phase orbit is inconsistent';return;endif
    local_bad=0
    do x=1,action%global_row_count;do g=1,nt
      target=action%element_maps(x,g)
      if(abs(global_phase(target)-ratio(g)*global_phase(x))>10d0*tolerance)local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='prepared translation phase fails covariance';return;endif
    do x=1,nlocal;local_phase(x)=global_phase(int(action%row_ids(x)));enddo
    fingerprint=metadata_hash;phase_payload_fingerprint=int(z'243F6A8885A308D3',int64)
    do x=1,action%global_row_count
      bits=transfer(real(global_phase(x),real64),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
      phase_payload_fingerprint=ieor(ishftc(phase_payload_fingerprint,11),bits)
      bits=transfer(aimag(global_phase(x)),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
      phase_payload_fingerprint=ieor(ishftc(phase_payload_fingerprint,11),bits)
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    if(phase_payload_fingerprint==0_int64)phase_payload_fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='prepared translation phase requires MPI';fingerprint=0_int64
    phase_payload_fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_complex(value)
      complex(real64),intent(in)::value
      bits=transfer(real(value,real64),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;global_bad=1;return;endif
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;global_bad=1;message='prepared character metadata disagree';return;endif
      bits=transfer(aimag(value),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;global_bad=1;return;endif
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minhash/=maxhash)
      if(global_bad/=0)message='prepared character metadata disagree'
    end subroutine agree_complex
#endif
  end subroutine build_dg_translation_character_intertwining_phase_prepared

  subroutine validate_dg_factored_point_cogroup_gauge(comm,local_basis,weights,point_maps,translation_maps,&
      point_product,point_identity,translation_cocycle,ntranslation,subspace_defect,tolerance,identity_defect,&
      unitarity_defect,closure_defect,workspace_peak_bytes,ok,message,generator_count,checked_pair_count,&
      prepared_operation_count)
    integer,intent(in)::comm,point_product(:,:),point_identity,translation_cocycle(:,:),ntranslation
    complex(real64),intent(in)::local_basis(:,:)
    real(real64),intent(in)::weights(:),subspace_defect,tolerance
    integer(int64),intent(in)::point_maps(:,:),translation_maps(:,:)
    real(real64),intent(out)::identity_defect,unitarity_defect,closure_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(out),optional::generator_count,checked_pair_count,prepared_operation_count
#ifdef USE_MPI
    integer::npoint,nstate,nlocal,left,right,product,cocycle,rank,nproc,owner,base,remainder,first,count,&
      local_bad,global_bad,ierr,status,minvalue,maxvalue,operation,index,nchecked,a,b,c,step
    integer,allocatable::point_generators(:)
    logical,allocatable::generator_mask(:)
    logical::generator_ok
    integer(int64),allocatable::global_left(:),global_right(:),global_product(:),global_translation(:),expected_map(:,:)
    integer(int64),allocatable::row_ids(:),expected_row_ids(:)
    complex(real64),allocatable::point_rows(:,:,:),expected_rows(:,:,:),product_rows(:,:),remote_rows(:,:)
    real(real64)::local_closure,global_closure,norm_bound,propagation_factor,propagated_closure
    integer(int64)::expected_workspace,global_extent,product_elements,persistent_bytes,term_bytes,&
      cache_bytes,live_bytes,maximum_remote_bytes,global_workspace_peak
    npoint=size(point_maps,2)
    nstate=size(local_basis,1);nlocal=size(local_basis,2)
    ok=.false.;message='';identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);workspace_peak_bytes=0_int64
    if(present(generator_count))generator_count=0
    if(present(checked_pair_count))checked_pair_count=0
    if(present(prepared_operation_count))prepared_operation_count=0
    local_bad=merge(0,1,npoint>=1.and.npoint<=48.and.nstate>=1.and.nlocal>=1.and.&
      point_identity>=1.and.point_identity<=npoint.and.&
      ntranslation>=1.and.all(shape(translation_cocycle)==[npoint,npoint]).and.&
      all(shape(point_product)==[npoint,npoint]).and.all(point_product>=1).and.all(point_product<=npoint).and.&
      size(point_maps,1)==nlocal.and.all(shape(translation_maps)==[nlocal,ntranslation]).and.&
      size(weights)==nlocal.and.all(translation_cocycle>=1).and.all(translation_cocycle<=ntranslation))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='invalid factored point-cogroup/cocycle contract'
      return
    endif
    call agree_integer(npoint);if(global_bad/=0)return
    call agree_integer(ntranslation);if(global_bad/=0)return
    call agree_integer(nstate);if(global_bad/=0)return
    call agree_integer(nlocal);if(global_bad/=0)return
    call agree_integer(point_identity);if(global_bad/=0)return
    do right=1,npoint;do left=1,npoint
      call agree_integer(point_product(left,right));if(global_bad/=0)return
      call agree_integer(translation_cocycle(left,right));if(global_bad/=0)return
    enddo;enddo
    local_bad=0
    do a=1,npoint;do b=1,npoint;do c=1,npoint
      if(point_product(point_product(a,b),c)/=point_product(a,point_product(b,c)))local_bad=1
    enddo;enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='point-cogroup product is not associative';return
    endif
    call select_dg_group_generators(point_product,point_identity,point_generators,generator_ok,message)
    local_bad=merge(0,1,generator_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='point-cogroup generator selection failed collectively';return
    endif
    allocate(generator_mask(npoint),stat=status)
    local_bad=merge(0,1,status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(generator_mask))deallocate(generator_mask)
      message='point-cogroup generator mask allocation failed';return
    endif
    generator_mask=.false.
    if(size(point_generators)>0)generator_mask(point_generators)=.true.
    if(present(generator_count))generator_count=size(point_generators)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='factored proof rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='factored proof size query failed';return;endif
    local_bad=0
    if(int(nlocal,int64)>huge(0_int64)/int(nproc,int64))local_bad=1
    if(local_bad==0)then
      global_extent=int(nlocal,int64)*int(nproc,int64)
      if(global_extent>int(huge(0),int64))local_bad=1
    endif
    if(int(max(1,nstate/nproc+1),int64)>huge(0_int64)/int(nstate,int64))local_bad=1
    if(local_bad==0)then
      product_elements=int(max(1,nstate/nproc+1),int64)*int(nstate,int64)
      if(product_elements>int(huge(0),int64))local_bad=1
      if(global_extent>huge(0_int64)/40_int64)then
        local_bad=1
      else
        persistent_bytes=40_int64*global_extent
        if(product_elements>huge(0_int64)/16_int64)then
          local_bad=1
        else
          term_bytes=16_int64*product_elements
          if(persistent_bytes>huge(0_int64)-term_bytes)local_bad=1
        endif
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='factored cocycle proof extent overflows';return;endif
    allocate(global_left(int(global_extent)),global_right(int(global_extent)),global_product(int(global_extent)),&
      global_translation(int(global_extent)),expected_map(nlocal,1),&
      product_rows(max(1,nstate/nproc+1),nstate),stat=status)
    local_bad=merge(0,1,status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup_factored_workspace();message='factored cocycle proof allocation failed';return
    endif
    ! Every spatial action column must be a permutation before any value is used as an index.
    do operation=1,npoint
      call MPI_Allgather(point_maps(:,operation),nlocal,MPI_INTEGER8,global_left,nlocal,MPI_INTEGER8,comm,ierr)
      local_bad=merge(0,1,ierr==MPI_SUCCESS)
      if(local_bad==0)local_bad=merge(0,1,all(global_left>=1_int64).and.all(global_left<=global_extent))
      if(local_bad==0)then
        global_right=0_int64
        do index=1,int(global_extent);global_right(int(global_left(index)))=global_right(int(global_left(index)))+1_int64;enddo
        local_bad=merge(0,1,all(global_right==1_int64))
      endif
      call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.local_bad/=0)then
        call cleanup_factored_workspace();message='point representative map is not a global permutation';return
      endif
    enddo
    do operation=1,ntranslation
      call MPI_Allgather(translation_maps(:,operation),nlocal,MPI_INTEGER8,global_left,nlocal,MPI_INTEGER8,comm,ierr)
      local_bad=merge(0,1,ierr==MPI_SUCCESS)
      if(local_bad==0)local_bad=merge(0,1,all(global_left>=1_int64).and.all(global_left<=global_extent))
      if(local_bad==0)then
        global_right=0_int64
        do index=1,int(global_extent);global_right(int(global_left(index)))=global_right(int(global_left(index)))+1_int64;enddo
        local_bad=merge(0,1,all(global_right==1_int64))
      endif
      call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.local_bad/=0)then
        call cleanup_factored_workspace();message='translation cocycle map is not a global permutation';return
      endif
    enddo
    call validate_dg_streamed_affine_representation(comm,local_basis,weights,point_maps,&
      point_identity,subspace_defect,tolerance,identity_defect,unitarity_defect,closure_defect,&
      workspace_peak_bytes,ok,message,row_ids,point_rows)
    if(.not.ok)then;call cleanup_factored_workspace();return;endif
    if(present(prepared_operation_count))prepared_operation_count=npoint
    if(persistent_bytes>huge(0_int64)-term_bytes)then
      message='factored persistent workspace peak overflows';call cleanup_factored_workspace();return
    endif
    persistent_bytes=persistent_bytes+term_bytes
    if(int(size(point_rows),int64)>huge(0_int64)/16_int64)then
      message='prepared point representation receipt overflows';call cleanup_factored_workspace();return
    endif
    cache_bytes=16_int64*int(size(point_rows),int64)
    if(persistent_bytes>huge(0_int64)-cache_bytes)then
      message='prepared point representation peak overflows';call cleanup_factored_workspace();return
    endif
    persistent_bytes=persistent_bytes+cache_bytes
    workspace_peak_bytes=max(workspace_peak_bytes,persistent_bytes)
    base=nstate/nproc;remainder=mod(nstate,nproc);local_closure=0d0
    maximum_remote_bytes=0_int64;local_bad=0
    do owner=0,nproc-1
      count=base+merge(1,0,owner<remainder)
      if(int(count,int64)>huge(0_int64)/int(nstate,int64))then
        local_bad=1
      else
        term_bytes=int(count,int64)*int(nstate,int64)
        if(term_bytes>huge(0_int64)/16_int64)then
          local_bad=1
        else
          maximum_remote_bytes=max(maximum_remote_bytes,16_int64*term_bytes)
        endif
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='factored remote tile preallocation receipt overflows';call cleanup_factored_workspace();return
    endif
    nchecked=0
    do left=1,npoint;do right=1,npoint
      product=point_product(right,left);cocycle=translation_cocycle(right,left)
      call MPI_Allgather(point_maps(:,left),nlocal,MPI_INTEGER8,global_left,nlocal,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='factored left map gather failed';call cleanup_factored_workspace();return;endif
      call MPI_Allgather(point_maps(:,right),nlocal,MPI_INTEGER8,global_right,nlocal,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='factored right map gather failed';call cleanup_factored_workspace();return;endif
      call MPI_Allgather(point_maps(:,product),nlocal,MPI_INTEGER8,global_product,nlocal,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='factored product map gather failed';call cleanup_factored_workspace();return;endif
      call MPI_Allgather(translation_maps(:,cocycle),nlocal,MPI_INTEGER8,global_translation,nlocal,MPI_INTEGER8,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='factored translation map gather failed';call cleanup_factored_workspace();return;endif
      local_bad=0
      do owner=1,nlocal
        ! factor_dg_affine_translation_cocycle defines r_left r_right = t_cocycle r_product.
        expected_map(owner,1)=global_translation(int(global_product(rank*nlocal+owner)))
        if(global_right(int(global_left(rank*nlocal+owner)))/=expected_map(owner,1))local_bad=1
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        ok=.false.;message='point representatives do not realize the supplied translation cocycle';&
        call cleanup_factored_workspace();return
      endif
      if(.not.generator_mask(left).and..not.generator_mask(right))cycle
      nchecked=nchecked+1
      product_rows=(0d0,0d0)
      live_bytes=persistent_bytes;local_bad=0
      if(local_bad==0)then
        if(maximum_remote_bytes>huge(0_int64)-live_bytes)then
          local_bad=1
        else
          live_bytes=live_bytes+maximum_remote_bytes
        endif
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        ok=.false.;message='factored cocycle live workspace preallocation receipt overflows';&
        call cleanup_factored_workspace();return
      endif
      workspace_peak_bytes=max(workspace_peak_bytes,live_bytes)
      do owner=0,nproc-1
        count=base+merge(1,0,owner<remainder);first=owner*base+min(owner,remainder)+1
        allocate(remote_rows(count,nstate),stat=status)
        call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
          ok=.false.;message='cocycle product tile allocation failed';call cleanup_factored_workspace();return
        endif
        if(rank==owner)remote_rows=point_rows(:,:,right)
        call MPI_Bcast(remote_rows,count*nstate,MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup_factored_workspace();message='cocycle product tile broadcast failed';return;endif
        product_rows(1:size(row_ids),:)=product_rows(1:size(row_ids),:)+&
          matmul(point_rows(:,first:first+count-1,left),remote_rows)
        deallocate(remote_rows)
      enddo
      call assemble_dg_distributed_basis_symmetry_overlap_rows(comm,local_basis,weights,&
        expected_map,expected_row_ids,expected_rows,expected_workspace,ok,message)
      if(.not.ok)then;call cleanup_factored_workspace();return;endif
      if(any(row_ids/=expected_row_ids))then
        ok=.false.;message='factored expected row ownership changed';call cleanup_factored_workspace();return
      endif
      local_closure=max(local_closure,maxval(abs(product_rows(1:size(row_ids),:)-expected_rows(:,:,1))))
      live_bytes=persistent_bytes;local_bad=0
      if(local_bad==0.and.expected_workspace>huge(0_int64)-live_bytes)then
        local_bad=1
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        ok=.false.;message='factored expected-action workspace receipt overflows';call cleanup_factored_workspace();return
      endif
      workspace_peak_bytes=max(workspace_peak_bytes,live_bytes+expected_workspace)
      deallocate(expected_row_ids,expected_rows)
    enddo;enddo
    call MPI_Allreduce(local_closure,global_closure,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(workspace_peak_bytes,global_workspace_peak,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)workspace_peak_bytes=global_workspace_peak
    ! Convert the measured max-entry unitarity defect into an operator-norm
    ! bound, then propagate generator residuals over both word directions.
    if(unitarity_defect>(huge(1d0)-1d0)/real(nstate,real64))then
      propagated_closure=huge(1d0)
    else
      norm_bound=sqrt(1d0+real(nstate,real64)*unitarity_defect)
      propagation_factor=0d0
      do step=1,2*npoint
        if(propagation_factor>(huge(1d0)-1d0)/max(1d0,norm_bound))then
          propagation_factor=huge(1d0);exit
        endif
        propagation_factor=1d0+norm_bound*propagation_factor
      enddo
      if(global_closure>0d0)then
        if(propagation_factor>huge(1d0)/global_closure)then
          propagated_closure=huge(1d0)
        else
          propagated_closure=propagation_factor*global_closure
        endif
      else
        propagated_closure=0d0
      endif
      if(propagated_closure>huge(1d0)/2d0)then
        propagated_closure=huge(1d0)
      else
        propagated_closure=2d0*propagated_closure
      endif
    endif
    closure_defect=max(closure_defect,propagated_closure)
    if(present(checked_pair_count))checked_pair_count=nchecked
    ok=ierr==MPI_SUCCESS.and.closure_defect<=tolerance
    if(ok)then;message='';else;message='factored point-cogroup internal action violates cocycle closure';endif
    call cleanup_factored_workspace()
  contains
    subroutine agree_integer(value)
      integer,intent(in)::value
      call MPI_Allreduce(value,minvalue,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;global_bad=1;message='factored point-cogroup metadata reduction failed';return;endif
      call MPI_Allreduce(value,maxvalue,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minvalue/=maxvalue)
      if(global_bad/=0)message='factored point-cogroup metadata disagree across ranks'
    end subroutine
    subroutine cleanup_factored_workspace()
      if(allocated(global_left))deallocate(global_left)
      if(allocated(global_right))deallocate(global_right)
      if(allocated(global_product))deallocate(global_product)
      if(allocated(global_translation))deallocate(global_translation)
      if(allocated(expected_map))deallocate(expected_map)
      if(allocated(product_rows))deallocate(product_rows)
      if(allocated(remote_rows))deallocate(remote_rows)
      if(allocated(row_ids))deallocate(row_ids)
      if(allocated(expected_row_ids))deallocate(expected_row_ids)
      if(allocated(point_rows))deallocate(point_rows)
      if(allocated(expected_rows))deallocate(expected_rows)
      if(allocated(point_generators))deallocate(point_generators)
      if(allocated(generator_mask))deallocate(generator_mask)
    end subroutine
#else
    ok=.false.;message='factored point-cogroup validation requires MPI'
    identity_defect=huge(1d0);unitarity_defect=huge(1d0);closure_defect=huge(1d0);workspace_peak_bytes=0_int64
#endif
  end subroutine validate_dg_factored_point_cogroup_gauge

  subroutine accumulate_dg_translation_character_orbit_sector(comm,state,character_index,characters,&
      catalog_fingerprint,sector_values,sector_gradients,initialize,finalize,tolerance,&
      orbit_values,orbit_gradients,workspace_bytes,ok,message)
    integer,intent(in)::comm,character_index
    type(s_dg_translation_orbit_accumulator),intent(inout)::state
    complex(real64),intent(in)::characters(:,:),sector_values(:,:),sector_gradients(:,:,:)
    integer(int64),intent(in)::catalog_fingerprint
    logical,intent(in)::initialize,finalize
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(inout)::orbit_values(:,:,:),orbit_gradients(:,:,:,:)
    integer(int64),intent(out)::workspace_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ntranslation,nlocal,ninternal,t,c,ierr,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer,initialize_integer,finalize_integer
    integer(int64)::minimum_fingerprint,maximum_fingerprint,elements,gradient_elements,&
      table_fingerprint,raw_bits,state_signature,minimum_state_signature,maximum_state_signature
    real(real64)::normalization,safe_magnitude
    workspace_bytes=0_int64
    ntranslation=size(characters,1);nlocal=size(sector_values,1);ninternal=size(sector_values,2)
    local_bad=merge(0,1,ntranslation>=1.and.ninternal>=1.and.size(characters,2)==ntranslation.and.&
      character_index>=1.and.character_index<=ntranslation.and.&
      catalog_fingerprint/=0_int64.and.tolerance>0d0.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(characters))).and.all(ieee_is_finite(aimag(characters))).and.&
      maxval(abs(abs(characters)-1d0))<=10d0*tolerance.and.&
      all(shape(sector_gradients)==[3,nlocal,ninternal]).and.&
      all(ieee_is_finite(real(sector_values))).and.all(ieee_is_finite(aimag(sector_values))).and.&
      all(ieee_is_finite(real(sector_gradients))).and.all(ieee_is_finite(aimag(sector_gradients))))
    if(.not.initialize)then
      if(.not.allocated(orbit_values).or..not.allocated(orbit_gradients))then
        local_bad=1
      elseif(any(shape(orbit_values)/=[nlocal,ninternal,ntranslation]).or.&
          any(shape(orbit_gradients)/=[3,nlocal,ninternal,ntranslation]))then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='invalid streamed inverse-character accumulation contract';return
    endif
    safe_magnitude=sqrt(huge(1d0))/max(4d0,4d0*real(ntranslation,real64))
    local_bad=merge(0,1,maxval(abs(sector_values))<=safe_magnitude.and.&
      maxval(abs(sector_gradients))<=safe_magnitude)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='streamed inverse sector magnitude is unsafe';return
    endif
    call agree_integer(character_index);if(global_bad/=0)return
    call agree_integer(ntranslation);if(global_bad/=0)return
    call agree_integer(ninternal);if(global_bad/=0)return
    initialize_integer=merge(1,0,initialize);call agree_integer(initialize_integer);if(global_bad/=0)return
    finalize_integer=merge(1,0,finalize);call agree_integer(finalize_integer);if(global_bad/=0)return
    call MPI_Allreduce(catalog_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='streamed inverse catalog reduction failed';return;endif
    call MPI_Allreduce(catalog_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      ok=.false.;message='streamed inverse catalogs disagree across ranks';return
    endif
    table_fingerprint=catalog_fingerprint
    do c=1,ntranslation;do t=1,ntranslation
      raw_bits=transfer(real(characters(c,t),real64),raw_bits)
      table_fingerprint=ieor(ishftc(table_fingerprint,11),raw_bits)
      raw_bits=transfer(aimag(characters(c,t)),raw_bits)
      table_fingerprint=ieor(ishftc(table_fingerprint,11),raw_bits)
    enddo;enddo
    call MPI_Allreduce(table_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='streamed inverse table reduction failed';return;endif
    call MPI_Allreduce(table_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      ok=.false.;message='streamed inverse character tables disagree across ranks';return
    endif
    state_signature=merge(1_int64,0_int64,state%initialized)
    state_signature=ieor(ishftc(state_signature,7),state%catalog_fingerprint)
    state_signature=ieor(ishftc(state_signature,7),state%table_fingerprint)
    if(allocated(state%visited))then
      state_signature=ieor(ishftc(state_signature,7),int(size(state%visited),int64))
      do t=1,size(state%visited)
        state_signature=ieor(ishftc(state_signature,7),merge(int(t,int64),0_int64,state%visited(t)))
      enddo
    endif
    call MPI_Allreduce(state_signature,minimum_state_signature,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='streamed inverse state reduction failed';return;endif
    call MPI_Allreduce(state_signature,maximum_state_signature,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_state_signature/=maximum_state_signature)then
      ok=.false.;message='streamed inverse state disagrees across ranks';return
    endif
    local_bad=0
    if(state%initialized.neqv.allocated(state%visited))local_bad=1
    if(state%initialized)then
      if(allocated(state%visited))then
        if(size(state%visited)/=ntranslation)local_bad=1
      endif
      if(state%catalog_fingerprint==0_int64.or.state%table_fingerprint==0_int64)local_bad=1
    elseif(state%catalog_fingerprint/=0_int64.or.state%table_fingerprint/=0_int64)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='streamed inverse state contract is invalid';return
    endif
    local_bad=0
    if(int(nlocal,int64)>huge(0_int64)/int(ninternal,int64))then
      local_bad=1
      elements=0_int64
    else
      elements=int(nlocal,int64)*int(ninternal,int64)
    endif
    if(local_bad==0.and.elements>huge(0_int64)/int(ntranslation,int64))then
      local_bad=1
    else if(local_bad==0)then
      elements=elements*int(ntranslation,int64)
    endif
    if(local_bad==0.and.elements>huge(0_int64)/3_int64)then
      local_bad=1
    else if(local_bad==0)then
      gradient_elements=3_int64*elements
    endif
    if(local_bad==0)then
      if(elements>huge(0_int64)-gradient_elements)then
        local_bad=1
      elseif(elements+gradient_elements>(huge(0_int64)-int(ntranslation,int64))/16_int64)then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='streamed inverse output extent or receipt overflows';return
    endif
    workspace_bytes=16_int64*(elements+gradient_elements)+int(ntranslation,int64)
    if(initialize)then
      local_bad=merge(1,0,state%initialized.or.allocated(state%visited))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        ok=.false.;message='streamed inverse cannot reinitialize an active sequence';return
      endif
      if(allocated(orbit_values))deallocate(orbit_values)
      if(allocated(orbit_gradients))deallocate(orbit_gradients)
      allocate(orbit_values(nlocal,ninternal,ntranslation),&
        orbit_gradients(3,nlocal,ninternal,ntranslation),stat=allocation_status)
      call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        if(allocated(orbit_values))deallocate(orbit_values)
        if(allocated(orbit_gradients))deallocate(orbit_gradients)
        ok=.false.;message='streamed inverse-character output allocation failed';return
      endif
      orbit_values=(0d0,0d0);orbit_gradients=(0d0,0d0)
      allocate(state%visited(ntranslation),stat=allocation_status)
      call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        if(allocated(state%visited))deallocate(state%visited)
        if(allocated(orbit_values))deallocate(orbit_values)
        if(allocated(orbit_gradients))deallocate(orbit_gradients)
        ok=.false.;message='streamed inverse state allocation failed';return
      endif
      state%visited=.false.;state%initialized=.true.;state%catalog_fingerprint=catalog_fingerprint
      state%table_fingerprint=table_fingerprint
    elseif(.not.state%initialized.or.state%catalog_fingerprint/=catalog_fingerprint.or.&
        state%table_fingerprint/=table_fingerprint)then
      ok=.false.;message='streamed inverse sequence is not initialized or catalog-bound';return
    endif
    local_bad=merge(1,0,state%visited(character_index))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='streamed inverse character sector is duplicated';return
    endif
    normalization=1d0/sqrt(real(ntranslation,real64))
    do t=1,ntranslation
      orbit_values(:,:,t)=orbit_values(:,:,t)+normalization*&
        conjg(characters(character_index,t))*sector_values
      orbit_gradients(:,:,:,t)=orbit_gradients(:,:,:,t)+normalization*&
        conjg(characters(character_index,t))*sector_gradients
    enddo
    state%visited(character_index)=.true.
    if(finalize)then
      local_bad=merge(1,0,.not.all(state%visited))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        ok=.false.;message='streamed inverse character sequence is incomplete';return
      endif
      state%initialized=.false.;state%catalog_fingerprint=0_int64;state%table_fingerprint=0_int64
      deallocate(state%visited)
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='streamed inverse-character accumulation requires MPI'
#endif
  contains
#ifdef USE_MPI
    subroutine agree_integer(value)
      integer,intent(in)::value
      call MPI_Allreduce(value,minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;global_bad=1;ok=.false.;message='streamed inverse metadata reduction failed';return;endif
      call MPI_Allreduce(value,maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)
      if(global_bad/=0)then;ok=.false.;message='streamed inverse metadata disagree across ranks';endif
    end subroutine agree_integer
#endif
  end subroutine accumulate_dg_translation_character_orbit_sector

  subroutine accumulate_dg_translation_character_orbit_sector_values(comm,state,character_index,characters,&
      catalog_fingerprint,sector_values,initialize,finalize,tolerance,orbit_values,workspace_bytes,ok,message)
    integer,intent(in)::comm,character_index
    type(s_dg_translation_orbit_accumulator),intent(inout)::state
    complex(real64),intent(in)::characters(:,:),sector_values(:,:)
    integer(int64),intent(in)::catalog_fingerprint
    logical,intent(in)::initialize,finalize
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(inout)::orbit_values(:,:)
    integer(int64),intent(out)::workspace_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nt,nlocal,m,t,c,ierr,local_bad,global_bad,allocation_status,flag,minflag,maxflag,minnt,maxnt,minm,maxm
    integer(int64)::elements,raw_bits,table_hash,minhash,maxhash,state_hash,minstate,maxstate
    real(real64)::normalization,safe_magnitude
    nt=size(characters,1);nlocal=size(sector_values,1);m=size(sector_values,2)
    ok=.false.;message='';workspace_bytes=0_int64
    local_bad=merge(0,1,nt>=1.and.m>=1.and.size(characters,2)==nt.and.character_index>=1.and.&
      character_index<=nt.and.catalog_fingerprint/=0_int64.and.tolerance>0d0.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(characters))).and.all(ieee_is_finite(aimag(characters))).and.&
      maxval(abs(abs(characters)-1d0))<=10d0*tolerance.and.&
      all(ieee_is_finite(real(sector_values))).and.all(ieee_is_finite(aimag(sector_values))))
    if(.not.initialize.and.allocated(orbit_values))then
      if(size(orbit_values,2)/=nlocal.or.&
          int(size(orbit_values,1),int64)/=int(m,int64)*int(nt,int64))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid values-only streamed inverse contract';return;endif
    call MPI_Allreduce(nt,minnt,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nt,maxnt,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minnt/=maxnt)then;message='values-only streamed translation extent disagrees';return;endif
    call MPI_Allreduce(m,minm,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxm,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minm/=maxm)then;message='values-only streamed internal extent disagrees';return;endif
    safe_magnitude=sqrt(huge(1d0))/(4d0*sqrt(real(nt,real64)))
    local_bad=merge(0,1,maxval(abs(sector_values))<=safe_magnitude)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='values-only streamed sector magnitude is unsafe';return;endif
    flag=character_index;call MPI_Allreduce(flag,minflag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(flag,maxflag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minflag/=maxflag)then;message='values-only streamed index disagrees';return;endif
    flag=merge(1,0,initialize);call MPI_Allreduce(flag,minflag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(flag,maxflag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minflag/=maxflag)then;message='values-only streamed initialize disagrees';return;endif
    flag=merge(1,0,finalize);call MPI_Allreduce(flag,minflag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(flag,maxflag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minflag/=maxflag)then;message='values-only streamed finalize disagrees';return;endif
    table_hash=catalog_fingerprint
    do c=1,nt;do t=1,nt
      raw_bits=transfer(real(characters(c,t),real64),raw_bits);table_hash=ieor(ishftc(table_hash,11),raw_bits)
      raw_bits=transfer(aimag(characters(c,t)),raw_bits);table_hash=ieor(ishftc(table_hash,11),raw_bits)
    enddo;enddo
    call MPI_Allreduce(table_hash,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(table_hash,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='values-only streamed catalog disagrees';return;endif
    state_hash=merge(1_int64,0_int64,state%initialized)
    state_hash=ieor(ishftc(state_hash,7),state%catalog_fingerprint)
    state_hash=ieor(ishftc(state_hash,7),state%table_fingerprint)
    if(allocated(state%visited))then
      state_hash=ieor(ishftc(state_hash,7),int(size(state%visited),int64))
      do t=1,size(state%visited);state_hash=ieor(ishftc(state_hash,7),merge(int(t,int64),0_int64,state%visited(t)));enddo
    endif
    call MPI_Allreduce(state_hash,minstate,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(state_hash,maxstate,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minstate/=maxstate)then;message='values-only streamed state disagrees';return;endif
    local_bad=0
    if(state%initialized.neqv.allocated(state%visited))local_bad=1
    if(state%initialized.and.allocated(state%visited))then;if(size(state%visited)/=nt)local_bad=1;endif
    if(int(m,int64)*int(nt,int64)>int(huge(0),int64))local_bad=1
    if(int(nlocal,int64)>huge(0_int64)/int(m,int64))local_bad=1
    if(local_bad==0)then;elements=int(nlocal,int64)*int(m,int64);else;elements=0_int64;endif
    if(local_bad==0.and.elements>huge(0_int64)/int(nt,int64))local_bad=1
    if(local_bad==0)elements=elements*int(nt,int64)
    if(local_bad==0.and.elements>(huge(0_int64)-int(nt,int64))/16_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='values-only streamed extent/state invalid';return;endif
    workspace_bytes=16_int64*elements+int(nt,int64)
    if(initialize)then
      local_bad=merge(1,0,state%initialized.or.allocated(state%visited))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='values-only streamed reinitialize rejected';return;endif
      if(allocated(orbit_values))deallocate(orbit_values)
      allocate(orbit_values(m*nt,nlocal),state%visited(nt),stat=allocation_status)
      call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        if(allocated(orbit_values))deallocate(orbit_values);if(allocated(state%visited))deallocate(state%visited)
        message='values-only streamed allocation failed';return
      endif
      orbit_values=(0d0,0d0);state%visited=.false.;state%initialized=.true.
      state%catalog_fingerprint=catalog_fingerprint;state%table_fingerprint=table_hash
    elseif(.not.state%initialized.or.state%catalog_fingerprint/=catalog_fingerprint.or.&
        state%table_fingerprint/=table_hash)then
      message='values-only streamed sequence is not catalog-bound';return
    endif
    local_bad=merge(1,0,state%visited(character_index))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='values-only streamed duplicate sector';return;endif
    normalization=1d0/sqrt(real(nt,real64))
    do t=1,nt
      orbit_values((t-1)*m+1:t*m,:)=orbit_values((t-1)*m+1:t*m,:)+&
        normalization*conjg(characters(character_index,t))*transpose(sector_values)
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(real(orbit_values))).and.all(ieee_is_finite(aimag(orbit_values))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='values-only streamed inverse accumulation became nonfinite';return
    endif
    state%visited(character_index)=.true.
    if(finalize)then
      local_bad=merge(1,0,.not.all(state%visited));call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='values-only streamed sequence incomplete';return;endif
      state%initialized=.false.;state%catalog_fingerprint=0_int64;state%table_fingerprint=0_int64
      deallocate(state%visited)
    endif
    ok=.true.
#else
    ok=.false.;message='values-only streamed inverse requires MPI';workspace_bytes=0_int64
#endif
  end subroutine accumulate_dg_translation_character_orbit_sector_values

  subroutine apply_dg_row_owned_orbital_transform_streamed(comm,row_ids,global_row_count,transform_rows,&
      input_values,input_gradients,output_values,output_gradients,workspace_bytes,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::transform_rows(:,:),input_values(:,:),input_gradients(:,:,:)
    complex(real64),allocatable,intent(out)::output_values(:,:),output_gradients(:,:,:)
    integer(int64),intent(out)::workspace_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nout,npoint,i,j,p,rank,ierr,local_bad,global_bad,allocation_status
    integer,allocatable::owner(:),position(:),count(:)
    complex(real64),allocatable::stream_row(:)
    integer(int64)::elements
    nlocal=size(row_ids);nout=size(transform_rows,2);npoint=size(input_values,2)
    ok=.false.;message='';workspace_bytes=0_int64
    local_bad=merge(0,1,global_row_count>=1.and.nout>=1.and.npoint>=0.and.&
      all(shape(transform_rows)==[nlocal,nout]).and.size(input_values,1)==global_row_count.and.&
      all(shape(input_gradients)==[3,global_row_count,npoint]).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      all(ieee_is_finite(real(transform_rows))).and.all(ieee_is_finite(aimag(transform_rows))).and.&
      all(ieee_is_finite(real(input_values))).and.all(ieee_is_finite(aimag(input_values))).and.&
      all(ieee_is_finite(real(input_gradients))).and.all(ieee_is_finite(aimag(input_gradients))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid row-streamed orbital transform contract';return;endif
    if(int(nout,int64)>huge(0_int64)/max(1_int64,int(npoint,int64)))then
      local_bad=1;elements=0_int64
    else
      local_bad=0;elements=int(nout,int64)*int(npoint,int64)
    endif
    if(local_bad==0.and.elements>(huge(0_int64)-int(nout,int64))/64_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='row-streamed orbital transform extent overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(owner(global_row_count),position(global_row_count),count(global_row_count),stream_row(nout),&
      output_values(nout,npoint),output_gradients(3,nout,npoint),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='row-streamed orbital transform allocation failed';return;endif
    owner=0;position=0;count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;count(int(row_ids(i)))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(count/=1))then;message='row-streamed transform rows are not uniquely owned';return;endif
    output_values=(0d0,0d0);output_gradients=(0d0,0d0)
    do i=1,global_row_count
      stream_row=(0d0,0d0)
      if(rank==owner(i)-1)stream_row=transform_rows(position(i),:)
      call MPI_Bcast(stream_row,nout,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='row-streamed transform broadcast failed';return;endif
      do p=1,npoint
        output_values(:,p)=output_values(:,p)+stream_row*input_values(i,p)
        do j=1,3;output_gradients(j,:,p)=output_gradients(j,:,p)+stream_row*input_gradients(j,i,p);enddo
      enddo
    enddo
    workspace_bytes=64_int64*elements+16_int64*int(nout,int64)
    ok=.true.
#else
    ok=.false.;message='row-streamed orbital transform requires MPI';workspace_bytes=0_int64
#endif
  end subroutine apply_dg_row_owned_orbital_transform_streamed

  subroutine validate_dg_translation_sector_cluster(eigenvalues,sector_real_dimension,tolerance,ok,message)
    ! EigenExa supplies eigenvalues in ascending order; reject unsorted external callers explicitly.
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer,intent(in)::sector_real_dimension
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64)::scale
    integer::i
    ok=.false.;message=''
    if(size(eigenvalues)<1.or.sector_real_dimension<1.or.sector_real_dimension>size(eigenvalues).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(eigenvalues)))then
      message='invalid translation-sector cluster contract';return
    endif
    scale=max(1d0,maxval(abs(eigenvalues)))
    do i=2,size(eigenvalues)
      if(eigenvalues(i)<eigenvalues(i-1))then
        message='translation-sector eigenvalues are not ordered';return
      endif
    enddo
    if(maxval(abs(eigenvalues(1:sector_real_dimension)))>10d0*tolerance*scale)then
      message='translation-sector selected cluster is rank losing';return
    endif
    if(sector_real_dimension<size(eigenvalues))then
      if(eigenvalues(sector_real_dimension+1)<=100d0*tolerance*scale)then
        message='translation-sector boundary splits a spectral cluster';return
      endif
    endif
    ok=.true.
  end subroutine validate_dg_translation_sector_cluster

  subroutine build_dg_smooth_partition_of_unity(comm,physical_ids,raw_weight,raw_gradient,&
      partition_weight,partition_gradient,sum_defect,gradient_defect,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::physical_ids(:)
    real(real64),intent(in)::raw_weight(:),raw_gradient(:,:)
    real(real64),intent(out)::partition_weight(:),partition_gradient(:,:)
    real(real64),intent(out)::sum_defect,gradient_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,nlocal,nlocal_max,source,source_count,p,q,ierr,local_bad,global_bad
    integer(int64),allocatable::source_ids(:)
    integer,allocatable::local_order(:)
    real(real64),allocatable::source_weight(:),source_gradient(:,:),denominator(:),&
      denominator_gradient(:,:),check_sum(:),check_gradient(:,:)
    real(real64)::local_sum_defect,local_gradient_defect
    ok=.false.;message='';sum_defect=huge(1d0);gradient_defect=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nlocal=size(physical_ids)
    call MPI_Allreduce(nlocal,nlocal_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.nlocal>0.and.&
      size(raw_weight)==nlocal.and.all(shape(raw_gradient)==[3,nlocal]).and.&
      size(partition_weight)==nlocal.and.all(shape(partition_gradient)==[3,nlocal]).and.&
      all(physical_ids>0_int64).and.all(raw_weight>=0d0).and.&
      all(ieee_is_finite(raw_weight)).and.all(ieee_is_finite(raw_gradient)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid smooth partition-of-unity contract';return
    endif
    allocate(source_ids(nlocal_max),source_weight(nlocal_max),source_gradient(3,nlocal_max),&
      denominator(nlocal),denominator_gradient(3,nlocal),check_sum(nlocal),check_gradient(3,nlocal),&
      local_order(nlocal))
    call sort_dg_int64_index(physical_ids,local_order)
    local_bad=0
    do p=2,nlocal
      if(physical_ids(local_order(p))==physical_ids(local_order(p-1)))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='smooth partition requires unique physical-grid ids within each fragment';return
    endif
    denominator=0d0;denominator_gradient=0d0
    do source=0,nproc-1
      source_count=nlocal
      call MPI_Bcast(source_count,1,MPI_INTEGER,source,comm,ierr)
      if(rank==source)then
        source_ids(1:source_count)=physical_ids
        source_weight(1:source_count)=raw_weight
        source_gradient(:,1:source_count)=raw_gradient
      endif
      call MPI_Bcast(source_ids,source_count,MPI_INTEGER8,source,comm,ierr)
      call MPI_Bcast(source_weight,source_count,MPI_DOUBLE_PRECISION,source,comm,ierr)
      call MPI_Bcast(source_gradient,3*source_count,MPI_DOUBLE_PRECISION,source,comm,ierr)
      do q=1,source_count
        p=find_dg_int64_index(physical_ids,local_order,source_ids(q))
        if(p==0)cycle
        denominator(p)=denominator(p)+source_weight(q)
        denominator_gradient(:,p)=denominator_gradient(:,p)+source_gradient(:,q)
      enddo
    enddo
    local_bad=merge(0,1,all(denominator>tiny(1d0)).and.all(ieee_is_finite(denominator)).and.&
      all(ieee_is_finite(denominator_gradient)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='smooth partition has missing or nonfinite physical-grid coverage';return
    endif
    partition_weight=raw_weight/denominator
    do p=1,nlocal
      partition_gradient(:,p)=(raw_gradient(:,p)*denominator(p)-&
        raw_weight(p)*denominator_gradient(:,p))/denominator(p)**2
    enddo
    check_sum=0d0;check_gradient=0d0
    do source=0,nproc-1
      source_count=nlocal
      call MPI_Bcast(source_count,1,MPI_INTEGER,source,comm,ierr)
      if(rank==source)then
        source_ids(1:source_count)=physical_ids
        source_weight(1:source_count)=partition_weight
        source_gradient(:,1:source_count)=partition_gradient
      endif
      call MPI_Bcast(source_ids,source_count,MPI_INTEGER8,source,comm,ierr)
      call MPI_Bcast(source_weight,source_count,MPI_DOUBLE_PRECISION,source,comm,ierr)
      call MPI_Bcast(source_gradient,3*source_count,MPI_DOUBLE_PRECISION,source,comm,ierr)
      do q=1,source_count
        p=find_dg_int64_index(physical_ids,local_order,source_ids(q))
        if(p==0)cycle
        check_sum(p)=check_sum(p)+source_weight(q)
        check_gradient(:,p)=check_gradient(:,p)+source_gradient(:,q)
      enddo
    enddo
    local_sum_defect=maxval(abs(check_sum-1d0));local_gradient_defect=maxval(abs(check_gradient))
    call MPI_Allreduce(local_sum_defect,sum_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_gradient_defect,gradient_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.ieee_is_finite(sum_defect).and.ieee_is_finite(gradient_defect)
    if(ok)then;message='';else;message='smooth partition normalization receipt is nonfinite';endif
#else
    ok=.false.;message='smooth partition-of-unity construction requires MPI'
    sum_defect=huge(1d0);gradient_defect=huge(1d0)
#endif
  end subroutine build_dg_smooth_partition_of_unity

  subroutine compose_dg_buffered_orbital_tile_to_physical_grid(comm,physical_ids,partition_weight,&
      buffer_values,owned_ids,owned_values,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::physical_ids(:)
    real(real64),intent(in)::partition_weight(:)
    complex(real64),intent(in)::buffer_values(:,:)
    integer(int64),allocatable,intent(out)::owned_ids(:)
    complex(real64),allocatable,intent(out)::owned_values(:,:)
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,norb,norb_min,norb_max,nlocal,local_bad,global_bad
    integer::p,owner,cursor,total_send,total_recv,nowned,i,j
    integer(int64)::scaled_count,local_hash,bits,real_bytes,complex_bytes,int_bytes,&
      local_max_id,global_point_count,points_per_owner
    integer,allocatable::send_counts(:),recv_counts(:),send_displs(:),recv_displs(:),fill(:),order(:),&
      value_send_counts(:),value_recv_counts(:),value_send_displs(:),value_recv_displs(:)
    integer(int64),allocatable::send_ids(:),recv_ids(:)
    real(real64),allocatable::send_weights(:),recv_weights(:),owned_weight_sum(:)
    complex(real64),allocatable::send_values(:,:),recv_values(:,:)
    logical::counts_ok

    ok=.false.;message='';fingerprint=0_int64;workspace_peak_bytes=0_int64
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    norb=size(buffer_values,1);nlocal=size(physical_ids)
    call MPI_Allreduce(norb,norb_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(norb,norb_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.nproc>0.and.norb>0.and.norb_min==norb_max.and.nlocal>0.and.&
      size(partition_weight)==nlocal.and.size(buffer_values,2)==nlocal.and.all(physical_ids>0_int64).and.&
      all(partition_weight>=0d0).and.all(ieee_is_finite(partition_weight)).and.&
      all(ieee_is_finite(real(buffer_values))).and.all(ieee_is_finite(aimag(buffer_values))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid buffer-composed orbital tile contract';return
    endif
    allocate(order(nlocal));call sort_dg_int64_index(physical_ids,order);local_bad=0
    do p=2,nlocal
      if(physical_ids(order(p))==physical_ids(order(p-1)))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr);deallocate(order)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='buffer composition requires unique physical-grid ids within each fragment';return
    endif
    local_max_id=maxval(physical_ids)
    call MPI_Allreduce(local_max_id,global_point_count,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.global_point_count>0_int64.and.&
      modulo(global_point_count,int(nproc,int64))==0_int64)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='buffer composition physical-grid ownership is not rank balanced';return
    endif
    points_per_owner=global_point_count/int(nproc,int64)
    allocate(send_counts(nproc),recv_counts(nproc),send_displs(nproc),recv_displs(nproc),fill(nproc))
    send_counts=0
    do p=1,nlocal
      owner=int((physical_ids(p)-1_int64)/points_per_owner)+1
      if(send_counts(owner)==huge(send_counts(owner)))then;local_bad=1;exit;endif
      send_counts(owner)=send_counts(owner)+1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='buffer composition send count overflow';return;endif
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr)
    call build_checked_mpi_displacements(send_counts,send_displs,total_send,counts_ok)
    if(counts_ok)call build_checked_mpi_displacements(recv_counts,recv_displs,total_recv,counts_ok)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.counts_ok.and.total_send==nlocal)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='buffer composition displacement overflow';return;endif
    allocate(send_ids(total_send),recv_ids(total_recv),send_weights(total_send),recv_weights(total_recv),&
      send_values(norb,total_send),recv_values(norb,total_recv));fill=send_displs
    do p=1,nlocal
      owner=int((physical_ids(p)-1_int64)/points_per_owner)+1
      cursor=fill(owner)+1;fill(owner)=cursor
      send_ids(cursor)=physical_ids(p);send_weights(cursor)=partition_weight(p)
      send_values(:,cursor)=buffer_values(:,p)
    enddo
    call MPI_Alltoallv(send_ids,send_counts,send_displs,MPI_INTEGER8,&
      recv_ids,recv_counts,recv_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='buffer composition ID exchange failed';return;endif
    call MPI_Alltoallv(send_weights,send_counts,send_displs,MPI_DOUBLE_PRECISION,&
      recv_weights,recv_counts,recv_displs,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='buffer composition weight exchange failed';return;endif
    allocate(value_send_counts(nproc),value_recv_counts(nproc),value_send_displs(nproc),value_recv_displs(nproc))
    local_bad=0
    do owner=1,nproc
      scaled_count=int(norb,int64)*int(send_counts(owner),int64)
      if(scaled_count>int(huge(0),int64))local_bad=1
      value_send_counts(owner)=int(min(scaled_count,int(huge(0),int64)))
      scaled_count=int(norb,int64)*int(recv_counts(owner),int64)
      if(scaled_count>int(huge(0),int64))local_bad=1
      value_recv_counts(owner)=int(min(scaled_count,int(huge(0),int64)))
    enddo
    call build_checked_mpi_displacements(value_send_counts,value_send_displs,cursor,counts_ok)
    if(counts_ok)call build_checked_mpi_displacements(value_recv_counts,value_recv_displs,cursor,counts_ok)
    call MPI_Allreduce(local_bad+merge(0,1,counts_ok),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='buffer composition value count overflow';return;endif
    call MPI_Alltoallv(send_values,value_send_counts,value_send_displs,MPI_DOUBLE_COMPLEX,&
      recv_values,value_recv_counts,value_recv_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='buffer composition exchange failed';return;endif
    allocate(order(total_recv));call sort_dg_int64_index(recv_ids,order)
    nowned=1
    do i=2,total_recv
      if(recv_ids(order(i))/=recv_ids(order(i-1)))nowned=nowned+1
    enddo
    allocate(owned_ids(nowned),owned_values(norb,nowned),owned_weight_sum(nowned))
    owned_values=(0d0,0d0);owned_weight_sum=0d0;j=0
    do i=1,total_recv
      if(i==1)then
        j=1;owned_ids(j)=recv_ids(order(i))
      elseif(recv_ids(order(i))/=recv_ids(order(i-1)))then
        j=j+1;owned_ids(j)=recv_ids(order(i))
      endif
      owned_values(:,j)=owned_values(:,j)+recv_weights(order(i))*recv_values(:,order(i))
      owned_weight_sum(j)=owned_weight_sum(j)+recv_weights(order(i))
    enddo
    local_bad=merge(0,1,nowned>0.and.all(abs(owned_weight_sum-1d0)<=1d3*epsilon(1d0)).and.&
      all(ieee_is_finite(real(owned_values))).and.all(ieee_is_finite(aimag(owned_values))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='buffer composition has missing, excess, or nonfinite physical-grid coverage';return
    endif
    local_hash=0_int64
    do p=1,nowned
      local_hash=ieor(local_hash,ishftc(owned_ids(p),int(modulo(owned_ids(p),63_int64))))
      do i=1,norb
        bits=transfer(real(owned_values(i,p)),bits)
        local_hash=ieor(local_hash,ishftc(bits,mod(7*i+11*int(modulo(owned_ids(p),63_int64)),63)))
        bits=transfer(aimag(owned_values(i,p)),bits)
        local_hash=ieor(local_hash,ishftc(bits,mod(13*i+17*int(modulo(owned_ids(p),63_int64)),63)))
      enddo
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(fingerprint,int(z'3C6EF372FE94F82B',int64))
    if(fingerprint==0_int64)fingerprint=1_int64
    real_bytes=int(storage_size(0d0)/8,int64);complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    int_bytes=int(storage_size(0_int64)/8,int64)
    workspace_peak_bytes=int_bytes*int(size(send_ids)+size(recv_ids)+size(owned_ids),int64)+&
      real_bytes*int(size(send_weights)+size(recv_weights)+size(owned_weight_sum),int64)+&
      complex_bytes*int(size(send_values)+size(recv_values)+size(owned_values),int64)+&
      int(storage_size(0)/8,int64)*int(size(send_counts)+size(recv_counts)+size(send_displs)+&
      size(recv_displs)+size(fill)+size(order)+size(value_send_counts)+size(value_recv_counts)+&
      size(value_send_displs)+size(value_recv_displs),int64)
    call MPI_Allreduce(MPI_IN_PLACE,workspace_peak_bytes,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.workspace_peak_bytes>0_int64
    if(ok)then;message='';else;message='buffer composition workspace receipt reduction failed';endif
#else
    ok=.false.;message='buffer-composed orbital tiles require MPI'
    fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine compose_dg_buffered_orbital_tile_to_physical_grid

  subroutine sort_dg_int64_index(values,order)
    integer(int64),intent(in)::values(:)
    integer,intent(out)::order(:)
    integer::i
    order=[(i,i=1,size(values))]
    call sort_range(1,size(order))
  contains
    recursive subroutine sort_range(left,right)
      integer,intent(in)::left,right
      integer::i,j,temporary
      integer(int64)::pivot
      if(left>=right)return
      i=left;j=right;pivot=values(order((left+right)/2))
      do
        do while(values(order(i))<pivot);i=i+1;enddo
        do while(values(order(j))>pivot);j=j-1;enddo
        if(i>j)exit
        temporary=order(i);order(i)=order(j);order(j)=temporary
        i=i+1;j=j-1
        if(i>j)exit
      enddo
      if(left<j)call sort_range(left,j)
      if(i<right)call sort_range(i,right)
    end subroutine sort_range
  end subroutine sort_dg_int64_index

  integer function find_dg_int64_index(values,order,target) result(location)
    integer(int64),intent(in)::values(:),target
    integer,intent(in)::order(:)
    integer::left,right,middle
    location=0;left=1;right=size(order)
    do while(left<=right)
      middle=(left+right)/2
      if(values(order(middle))<target)then
        left=middle+1
      elseif(values(order(middle))>target)then
        right=middle-1
      else
        location=order(middle);return
      endif
    enddo
  end function find_dg_int64_index

  subroutine select_dg_group_generators(product_table,identity_operation,generators,ok,message)
    integer,intent(in)::product_table(:,:),identity_operation
    integer,allocatable,intent(out)::generators(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::reached(:)
    integer,allocatable::work_generators(:)
    integer::n,operation,generator_index,ngenerator,status
    logical::changed,has_inverse

    ok=.false.;message='';n=size(product_table,1)
    if(n<1.or.size(product_table,2)/=n.or.identity_operation<1.or.identity_operation>n.or.&
        any(product_table<1).or.any(product_table>n))then
      message='invalid group-generator product table';return
    endif
    do operation=1,n
      if(product_table(identity_operation,operation)/=operation.or.&
          product_table(operation,identity_operation)/=operation)then
        message='group-generator identity operation is invalid';return
      endif
      has_inverse=any(product_table(operation,:)==identity_operation.and.&
        product_table(:,operation)==identity_operation)
      if(.not.has_inverse)then;message='group-generator product table lacks an inverse';return;endif
    enddo
    allocate(reached(n),work_generators(max(0,n-1)),stat=status)
    if(status/=0)then;message='group-generator workspace allocation failed';return;endif
    reached=.false.;reached(identity_operation)=.true.;ngenerator=0
    do while(.not.all(reached))
      do operation=1,n
        if(.not.reached(operation))exit
      enddo
      ngenerator=ngenerator+1;work_generators(ngenerator)=operation;reached(operation)=.true.
      changed=.true.
      do while(changed)
        changed=.false.
        do operation=1,n
          if(.not.reached(operation))cycle
          do generator_index=1,ngenerator
            if(.not.reached(product_table(operation,work_generators(generator_index))))then
              reached(product_table(operation,work_generators(generator_index)))=.true.;changed=.true.
            endif
            if(.not.reached(product_table(work_generators(generator_index),operation)))then
              reached(product_table(work_generators(generator_index),operation))=.true.;changed=.true.
            endif
          enddo
        enddo
      enddo
    enddo
    allocate(generators(ngenerator),stat=status)
    if(status/=0)then;message='group-generator result allocation failed';return;endif
    generators=work_generators(1:ngenerator)
    ok=.true.
  end subroutine select_dg_group_generators

#ifdef USE_EIGENEXA
  subroutine measure_dg_rank_fixed_symmetry_residuals_eigenexa(info,comm,basis,weights,&
      symmetry_target_box_ids,boundary_mask,total_residual,boundary_residual,interior_residual,&
      ok,message,workspace_peak_bytes)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::comm
    complex(real64),intent(in)::basis(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    logical,intent(in)::boundary_mask(:)
    real(real64),intent(out)::total_residual(:),boundary_residual(:),interior_residual(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(out),optional::workspace_peak_bytes
    integer,parameter::orbital_tile_size=32
    real(real64),allocatable::local_cyclic_metric(:,:),local_cyclic_vectors(:,:),metric_spectrum(:),&
      eigenvector_tile(:,:),eigenvector_row(:),inverse_root_tile(:,:),local_norms(:),global_norms(:)
    complex(real64),allocatable::orthonormal_basis(:,:),image_tile(:,:),overlap_tile(:,:),&
      reduced_overlap_tile(:,:)
    integer(int64)::metric_peak,peak_bytes,current_bytes,real_bytes,complex_bytes
    integer::nstate,nlocal,noperation,tile_first,tile_count,i,j,operation,ierr
    logical::eigen_ok
    character(256)::detail

    ok=.false.;message='';metric_peak=0_int64;peak_bytes=0_int64;current_bytes=0_int64
    if(present(workspace_peak_bytes))workspace_peak_bytes=0_int64
    if(.not.info%flag_eigenexa_init.or.size(basis,1)<=0)then
      message='OW-sized EigenExa descriptor is not initialized';return
    endif
    call assemble_dg_eigenexa_cyclic_metric_block(comm,info%nprow,info%npcol,info%myrow,info%mycol,&
      info%nrow_local,info%ncol_local,basis,weights,local_cyclic_metric,metric_peak,ok,detail)
    if(.not.ok)then;message='distributed rank-fixed metric: '//trim(detail);return;endif
    allocate(local_cyclic_vectors(info%nrow_local,info%ncol_local),metric_spectrum(size(basis,1)))
    call eigen_pdsyevd_ex_distributed_blocks(info,size(basis,1),local_cyclic_metric,metric_spectrum,&
      local_cyclic_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then
      ok=.false.;peak_bytes=metric_peak*int(storage_size(0d0)/8,int64)
      if(present(workspace_peak_bytes))workspace_peak_bytes=peak_bytes
      message='distributed rank-fixed eigensystem: '//trim(detail);return
    endif
    if(minval(metric_spectrum)<=epsilon(1d0)*max(1d0,maxval(metric_spectrum)))then
      ok=.false.;peak_bytes=metric_peak*int(storage_size(0d0)/8,int64)
      if(present(workspace_peak_bytes))workspace_peak_bytes=peak_bytes
      message='distributed rank-fixed occupied metric is singular';return
    endif
    nstate=size(basis,1);nlocal=size(basis,2);noperation=size(symmetry_target_box_ids,2)
    ok=nlocal>0.and.noperation>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.size(boundary_mask)==nlocal.and.&
      size(total_residual)==noperation.and.size(boundary_residual)==noperation.and.&
      size(interior_residual)==noperation.and.all(weights>=0d0)
    if(.not.ok)then;message='invalid distributed rank-fixed residual contract';return;endif
    real_bytes=int(storage_size(0d0)/8,int64);complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    allocate(orthonormal_basis(nstate,nlocal),eigenvector_row(nstate),local_norms(3),global_norms(3))
    current_bytes=real_bytes*int(size(local_cyclic_metric)+size(local_cyclic_vectors)+&
      size(metric_spectrum)+size(eigenvector_row)+size(local_norms)+size(global_norms),int64)+&
      complex_bytes*int(size(orthonormal_basis),int64)
    peak_bytes=max(metric_peak*real_bytes,current_bytes)
    do tile_first=1,nstate,orbital_tile_size
      tile_count=min(orbital_tile_size,nstate-tile_first+1)
      allocate(eigenvector_tile(tile_count,nstate),inverse_root_tile(tile_count,nstate))
      call gather_cyclic_eigenvector_rows(tile_first,tile_count,eigenvector_tile,ok,detail)
      if(.not.ok)then;message=trim(detail);return;endif
      do j=1,nstate
        call gather_cyclic_eigenvector_rows(j,1,eigenvector_row,ok,detail)
        if(.not.ok)then;message=trim(detail);return;endif
        do i=1,tile_count
          inverse_root_tile(i,j)=sum(eigenvector_tile(i,:)*eigenvector_row/sqrt(metric_spectrum))
        enddo
      enddo
      orthonormal_basis(tile_first:tile_first+tile_count-1,:)=matmul(inverse_root_tile,basis)
      peak_bytes=max(peak_bytes,current_bytes+real_bytes*&
        int(size(eigenvector_tile)+size(inverse_root_tile),int64))
      deallocate(eigenvector_tile,inverse_root_tile)
    enddo
    total_residual=0d0;boundary_residual=0d0;interior_residual=0d0
    do operation=1,noperation
      local_norms=0d0
      do tile_first=1,nstate,orbital_tile_size
        tile_count=min(orbital_tile_size,nstate-tile_first+1)
        allocate(image_tile(tile_count,nlocal),overlap_tile(nstate,tile_count),&
          reduced_overlap_tile(nstate,tile_count));overlap_tile=(0d0,0d0)
        call exchange_dg_point_permuted_orbital_rows(comm,&
          orthonormal_basis(tile_first:tile_first+tile_count-1,:),&
          symmetry_target_box_ids(:,operation),image_tile,ok,detail)
        if(.not.ok)then;message=trim(detail);return;endif
        do i=1,nstate
          overlap_tile(i,:)=matmul(conjg(orthonormal_basis(i,:))*weights,transpose(image_tile))
        enddo
        call MPI_Allreduce(overlap_tile,reduced_overlap_tile,nstate*tile_count,&
          MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;ok=.false.;message='distributed affine overlap reduction failed';return;endif
        image_tile=image_tile-matmul(transpose(reduced_overlap_tile),orthonormal_basis)
        local_norms(1)=local_norms(1)+sum(spread(weights,1,tile_count)*abs(image_tile)**2)
        local_norms(2)=local_norms(2)+sum(spread(weights*merge(1d0,0d0,boundary_mask),1,tile_count)*&
          abs(image_tile)**2)
        local_norms(3)=local_norms(3)+sum(spread(weights*merge(0d0,1d0,boundary_mask),1,tile_count)*&
          abs(image_tile)**2)
        peak_bytes=max(peak_bytes,current_bytes+complex_bytes*&
          int(size(image_tile)+size(overlap_tile)+size(reduced_overlap_tile),int64))
        deallocate(image_tile,overlap_tile,reduced_overlap_tile)
      enddo
      call MPI_Allreduce(local_norms,global_norms,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='distributed affine norm reduction failed';return;endif
      total_residual(operation)=sqrt(max(0d0,global_norms(1)))
      boundary_residual(operation)=sqrt(max(0d0,global_norms(2)))
      interior_residual(operation)=sqrt(max(0d0,global_norms(3)))
    enddo
    ok=all(ieee_is_finite(total_residual)).and.all(ieee_is_finite(boundary_residual)).and.&
      all(ieee_is_finite(interior_residual))
    if(.not.ok)then;message='distributed affine residual is nonfinite';return;endif
    if(present(workspace_peak_bytes))workspace_peak_bytes=peak_bytes
    message=''
  contains
    subroutine gather_cyclic_eigenvector_rows(first_row,row_count,rows,rows_ok,rows_message)
      integer,intent(in)::first_row,row_count
      real(real64),intent(out)::rows(..)
      logical,intent(out)::rows_ok
      character(*),intent(out)::rows_message
      real(real64),allocatable::local_rows(:,:)
      integer::global_row,global_col,local_row,local_col,row_offset,collective_error
      allocate(local_rows(row_count,nstate));local_rows=0d0
      do row_offset=1,row_count
        global_row=first_row+row_offset-1
        if(mod(global_row-1,info%nprow)/=info%myrow-1)cycle
        local_row=(global_row-1)/info%nprow+1
        do global_col=1,nstate
          if(mod(global_col-1,info%npcol)/=info%mycol-1)cycle
          local_col=(global_col-1)/info%npcol+1
          local_rows(row_offset,global_col)=local_cyclic_vectors(local_row,local_col)
        enddo
      enddo
      select rank(rows)
      rank(1)
        call MPI_Allreduce(local_rows(1,:),rows,nstate,MPI_DOUBLE_PRECISION,MPI_SUM,comm,collective_error)
      rank(2)
        call MPI_Allreduce(local_rows,rows,row_count*nstate,MPI_DOUBLE_PRECISION,MPI_SUM,comm,collective_error)
      end select
      rows_ok=collective_error==MPI_SUCCESS
      if(rows_ok)then;rows_message='';else;rows_message='cyclic eigenvector row gather failed';endif
    end subroutine gather_cyclic_eigenvector_rows
  end subroutine measure_dg_rank_fixed_symmetry_residuals_eigenexa
#endif

  subroutine build_checked_mpi_displacements(counts,displacements,total_count,ok)
    integer,intent(in)::counts(:)
    integer,intent(out)::displacements(:),total_count
    logical,intent(out)::ok
    integer(int64)::running
    integer::i
    ok=size(counts)>0.and.size(displacements)==size(counts).and.all(counts>=0)
    running=0_int64;total_count=0
    if(.not.ok)return
    do i=1,size(counts)
      if(running>int(huge(total_count),int64))then;ok=.false.;return;end if
      displacements(i)=int(running)
      running=running+int(counts(i),int64)
    end do
    if(running>int(huge(total_count),int64))then;ok=.false.;return;end if
    total_count=int(running)
  end subroutine build_checked_mpi_displacements

  subroutine assign_dg_periodic_centers_to_fragments(global_grid,centers,all_core_ids,fragment_ids,&
      tolerance,center_ids,center_owners,center_fragments,ok,message)
    integer,intent(in)::global_grid(3),fragment_ids(:)
    real(real64),intent(in)::centers(:,:),tolerance
    integer(int64),intent(in)::all_core_ids(:,:)
    integer(int64),allocatable,intent(out)::center_ids(:)
    integer,allocatable,intent(out)::center_owners(:),center_fragments(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::ncenter,center,axis,owner,location(2),grid_index(3)
    real(real64)::scaled,boundary_distance
    integer(int64)::global_count

    ncenter=size(centers,2);global_count=int(global_grid(1),int64)*int(global_grid(2),int64)*&
      int(global_grid(3),int64)
    ok=all(global_grid>0).and.ncenter>0.and.size(centers,1)==3.and.&
      size(all_core_ids,2)==size(fragment_ids).and.size(all_core_ids)==int(global_count).and.&
      tolerance>0d0.and.all(ieee_is_finite(centers)).and.all(fragment_ids>0)
    if(.not.ok)then;message='invalid periodic center-to-fragment ownership contract';return;end if
    if(any(all_core_ids<1_int64).or.any(all_core_ids>global_count))then
      ok=.false.;message='center ownership core ID is outside the global grid';return
    end if
    allocate(center_ids(ncenter),center_owners(ncenter),center_fragments(ncenter))
    do center=1,ncenter
      do axis=1,3
        scaled=modulo(centers(axis,center),1d0)*real(global_grid(axis),real64)
        boundary_distance=abs(scaled+0.5d0-anint(scaled+0.5d0))
        if(boundary_distance<=tolerance*real(global_grid(axis),real64))then
          grid_index(axis)=modulo(ceiling(scaled-0.5d0),global_grid(axis))
        else
          grid_index(axis)=modulo(floor(scaled+0.5d0),global_grid(axis))
        end if
      end do
      center_ids(center)=1_int64+int(grid_index(1),int64)+int(global_grid(1),int64)*(&
        int(grid_index(2),int64)+int(global_grid(2),int64)*int(grid_index(3),int64))
      location=findloc(all_core_ids,center_ids(center))
      if(any(location<1))then
        ok=.false.;message='periodic center is not covered by a unique fragment core';return
      end if
      if(count(all_core_ids==center_ids(center))/=1)then
        ok=.false.;message='periodic center core ownership is not unique';return
      end if
      owner=location(2)-1;center_owners(center)=owner;center_fragments(center)=fragment_ids(owner+1)
    end do
    ok=.true.;message=''
  end subroutine assign_dg_periodic_centers_to_fragments


  subroutine redistribute_dg_buffer_orbitals_to_center_fragments(comm,buffer_values,core_positions,&
      core_ids,center_owners,buffer_ids,local_orbitals,local_values,ok,message)
    integer,intent(in)::comm,core_positions(:),center_owners(:)
    complex(real64),intent(in)::buffer_values(:,:)
    integer(int64),intent(in)::core_ids(:),buffer_ids(:)
    integer,allocatable,intent(out)::local_orbitals(:)
    complex(real64),allocatable,intent(out)::local_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,norbital,ncore,nbuffer,nglobal,source,destination,orbital,point,&
      source_position,position,total_send,total_receive,local_bad,global_bad,status,minval_i,maxval_i
    integer,allocatable::send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:),&
      buffer_owner(:),occurrence(:)
    integer(int64),allocatable::all_core_ids(:,:),all_buffer_ids(:,:)
    complex(real64),allocatable::send_values(:),receive_values(:)
    integer(int64)::count64

    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    norbital=size(center_owners);ncore=size(core_ids);nbuffer=size(buffer_ids)
    if(norbital<1.or.ncore<1.or.nbuffer<1.or.size(core_positions)/=ncore.or.&
        size(buffer_values,1)/=norbital.or.any(core_positions<1).or.&
        any(core_positions>size(buffer_values,2)).or.any(center_owners<0).or.&
        any(center_owners>=nproc).or.any(core_ids<1_int64).or.any(buffer_ids<1_int64).or.&
        .not.all(ieee_is_finite(real(buffer_values))).or.&
        .not.all(ieee_is_finite(aimag(buffer_values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid direct buffer redistribution contract';return;endif
    call MPI_Allreduce(norbital,minval_i,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(norbital,maxval_i,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minval_i/=maxval_i)then;message='direct redistribution orbital extent disagrees';return;endif
    call MPI_Allreduce(ncore,minval_i,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(ncore,maxval_i,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minval_i/=maxval_i)then;message='direct redistribution core extent disagrees';return;endif
    call MPI_Allreduce(nbuffer,minval_i,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(nbuffer,maxval_i,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minval_i/=maxval_i.or.int(ncore,int64)>int(huge(1),int64)/int(nproc,int64))then
      message='direct redistribution buffer extent disagrees or overflows';return
    endif
    nglobal=ncore*nproc
    allocate(all_core_ids(ncore,nproc),all_buffer_ids(nbuffer,nproc),send_counts(nproc),&
      receive_counts(nproc),send_displacements(nproc),receive_displacements(nproc),&
      buffer_owner(nbuffer),occurrence(nglobal),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='direct redistribution metadata allocation failed';return;endif
    call MPI_Allgather(core_ids,ncore,MPI_INTEGER8,all_core_ids,ncore,MPI_INTEGER8,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allgather(buffer_ids,nbuffer,MPI_INTEGER8,all_buffer_ids,nbuffer,MPI_INTEGER8,comm,ierr);if(ierr/=MPI_SUCCESS)return
    occurrence=0
    do source=1,nproc;do point=1,ncore
      if(all_core_ids(point,source)<1_int64.or.all_core_ids(point,source)>int(nglobal,int64))then
        local_bad=1
      else
        occurrence(int(all_core_ids(point,source)))=occurrence(int(all_core_ids(point,source)))+1
      endif
    enddo;enddo
    if(any(occurrence/=1))local_bad=1
    do point=1,nbuffer
      buffer_owner(point)=0
      do source=1,nproc
        if(any(all_core_ids(:,source)==buffer_ids(point)))then;buffer_owner(point)=source;exit;endif
      enddo
      if(buffer_owner(point)==0)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='direct redistribution physical IDs are incomplete';return;endif
    do destination=0,nproc-1
      count64=int(count(center_owners==destination),int64)*int(count([(any(core_ids==&
        all_buffer_ids(point,destination+1)),point=1,nbuffer)]),int64)
      if(count64>int(huge(1),int64))then;local_bad=1;send_counts(destination+1)=0
      else;send_counts(destination+1)=int(count64);endif
      count64=int(count(center_owners==rank),int64)*int(count(buffer_owner==destination+1),int64)
      if(count64>int(huge(1),int64))then;local_bad=1;receive_counts(destination+1)=0
      else;receive_counts(destination+1)=int(count64);endif
    enddo
    call build_checked_mpi_displacements(send_counts,send_displacements,total_send,ok);if(.not.ok)local_bad=1
    call build_checked_mpi_displacements(receive_counts,receive_displacements,total_receive,ok);if(.not.ok)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='direct redistribution MPI counts overflow';return;endif
    allocate(send_values(total_send),receive_values(total_receive),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='direct redistribution value allocation failed';return;endif
    position=0
    do destination=0,nproc-1;do orbital=1,norbital
      if(center_owners(orbital)/=destination)cycle
      do point=1,nbuffer
        source_position=findloc(core_ids,all_buffer_ids(point,destination+1),dim=1)
        if(source_position<1)cycle
        position=position+1;send_values(position)=buffer_values(orbital,core_positions(source_position))
      enddo
    enddo;enddo
    if(position/=total_send)then;message='direct redistribution send packing is incomplete';return;endif
    call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,&
      receive_values,receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='direct redistribution MPI Alltoallv failed';return;endif
    local_orbitals=pack([(orbital,orbital=1,norbital)],center_owners==rank)
    allocate(local_values(size(local_orbitals),nbuffer),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='direct redistribution output allocation failed';return;endif
    local_values=(0d0,0d0);position=0
    do source=1,nproc;do orbital=1,size(local_orbitals);do point=1,nbuffer
      if(buffer_owner(point)/=source)cycle
      position=position+1;local_values(orbital,point)=receive_values(position)
    enddo;enddo;enddo
    ok=position==total_receive.and.all(ieee_is_finite(real(local_values))).and.&
      all(ieee_is_finite(aimag(local_values)))
    if(.not.ok)message='direct redistribution produced incomplete or nonfinite values'
#else
    ok=.false.;message='direct buffer redistribution requires MPI'
#endif
  end subroutine redistribute_dg_buffer_orbitals_to_center_fragments

  subroutine redistribute_dg_owned_orbitals_to_center_fragments(comm,owned_orbitals,global_ids,&
      owned_values,center_owners,local_buffer_ids,local_orbitals,local_values,ok,message)
    integer,intent(in)::comm,owned_orbitals(:),center_owners(:)
    integer(int64),intent(in)::global_ids(:),local_buffer_ids(:)
    complex(real64),intent(in)::owned_values(:,:)
    integer,allocatable,intent(out)::local_orbitals(:)
    complex(real64),allocatable,intent(out)::local_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,allocatable::orbital_counts(:),orbital_displacements(:),orbital_owners(:),buffer_counts(:),&
      buffer_displacements(:),send_counts(:),receive_counts(:),send_displacements(:),&
      receive_displacements(:),position_by_id(:),source_orbitals(:)
    integer(int64),allocatable::all_buffer_ids(:)
    complex(real64),allocatable::send_values(:),receive_values(:)
    integer::rank,nproc,ierr,norbital,nglobal,destination,source,orbital,point,position,local_row,&
      local_bad,global_bad,total_send,total_receive,total_buffer,min_norbital,max_norbital
    integer(int64)::count64

    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    norbital=size(center_owners);nglobal=size(global_ids);local_bad=0
    if(norbital<1.or.nglobal<1.or.size(owned_values,1)/=size(owned_orbitals).or.&
        size(owned_values,2)/=nglobal.or.size(local_buffer_ids)<1.or.&
        any(center_owners<0).or.any(center_owners>=nproc))local_bad=1
    call MPI_Allreduce(norbital,min_norbital,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(norbital,max_norbital,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(min_norbital/=max_norbital.or.ierr/=MPI_SUCCESS)local_bad=1
    call build_dg_balanced_orbital_ownership(norbital,nproc,orbital_counts,orbital_displacements,&
      orbital_owners,ok,message)
    if(.not.ok)local_bad=1
    if(ok)then
      if(size(owned_orbitals)/=orbital_counts(rank+1))local_bad=1
      if(size(owned_orbitals)>0)then
        if(any(owned_orbitals<1).or.any(owned_orbitals>norbital))then
          local_bad=1
        else if(any(orbital_owners(owned_orbitals)/=rank))then
          local_bad=1
        end if
      end if
    end if
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='invalid orbital-to-center-fragment redistribution contract';return
    end if
    allocate(buffer_counts(nproc),buffer_displacements(nproc),send_counts(nproc),receive_counts(nproc),&
      send_displacements(nproc),receive_displacements(nproc),position_by_id(nglobal))
    call MPI_Allgather(size(local_buffer_ids),1,MPI_INTEGER,buffer_counts,1,MPI_INTEGER,comm,ierr)
    call build_checked_mpi_displacements(buffer_counts,buffer_displacements,total_buffer,ok)
    if(.not.ok)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='orbital redistribution buffer displacement overflow';return
    end if
    allocate(all_buffer_ids(total_buffer))
    call MPI_Allgatherv(local_buffer_ids,size(local_buffer_ids),MPI_INTEGER8,all_buffer_ids,&
      buffer_counts,buffer_displacements,MPI_INTEGER8,comm,ierr)
    position_by_id=0
    do point=1,nglobal
      if(global_ids(point)<1_int64.or.global_ids(point)>int(nglobal,int64))then;local_bad=1;cycle;end if
      if(position_by_id(int(global_ids(point)))/=0)then;local_bad=1;cycle;end if
      position_by_id(int(global_ids(point)))=point
    end do
    if(any(position_by_id==0))local_bad=1
    do destination=0,nproc-1
      count64=int(count(center_owners(owned_orbitals)==destination),int64)*&
        int(buffer_counts(destination+1),int64)
      if(count64>int(huge(1),int64))then
        local_bad=1;send_counts(destination+1)=0
      else
        send_counts(destination+1)=int(count64)
      end if
      source_orbitals=pack([(orbital,orbital=1,norbital)],&
        orbital_owners==destination.and.center_owners==rank)
      count64=int(size(source_orbitals),int64)*int(size(local_buffer_ids),int64)
      if(count64>int(huge(1),int64))then
        local_bad=1;receive_counts(destination+1)=0
      else
        receive_counts(destination+1)=int(count64)
      end if
    end do
    if(any(all_buffer_ids<1_int64).or.any(all_buffer_ids>int(nglobal,int64)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='orbital-to-fragment MPI count or physical ID is invalid';return
    end if
    call build_checked_mpi_displacements(send_counts,send_displacements,total_send,ok)
    if(.not.ok)local_bad=1
    call build_checked_mpi_displacements(receive_counts,receive_displacements,total_receive,ok)
    if(.not.ok)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='orbital redistribution MPI displacement overflow';return
    end if
    allocate(send_values(total_send),receive_values(total_receive))
    position=0
    do destination=0,nproc-1
      do local_row=1,size(owned_orbitals)
        orbital=owned_orbitals(local_row)
        if(center_owners(orbital)/=destination)cycle
        do point=1,buffer_counts(destination+1)
          position=position+1
          send_values(position)=owned_values(local_row,position_by_id(int(all_buffer_ids(&
            buffer_displacements(destination+1)+point))))
        end do
      end do
    end do
    call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,&
      receive_values,receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='orbital-to-center-fragment MPI Alltoallv failed';return;end if
    local_orbitals=pack([(orbital,orbital=1,norbital)],center_owners==rank)
    allocate(local_values(size(local_orbitals),size(local_buffer_ids)));local_values=(0d0,0d0)
    position=0
    do source=0,nproc-1
      source_orbitals=pack([(orbital,orbital=1,norbital)],&
        orbital_owners==source.and.center_owners==rank)
      do orbital=1,size(source_orbitals)
        local_row=findloc(local_orbitals,source_orbitals(orbital),dim=1)
        do point=1,size(local_buffer_ids)
          position=position+1;local_values(local_row,point)=receive_values(position)
        end do
      end do
    end do
    ok=all(ieee_is_finite(real(local_values))).and.all(ieee_is_finite(aimag(local_values)))
    if(ok)then;message='';else;message='orbital-to-fragment redistribution produced nonfinite values';end if
#else
    ok=.false.;message='orbital-to-center-fragment redistribution requires MPI'
#endif
  end subroutine redistribute_dg_owned_orbitals_to_center_fragments


  subroutine build_dg_balanced_orbital_ownership(norbital,nproc,counts,displacements,owners,ok,message)
    integer,intent(in)::norbital,nproc
    integer,allocatable,intent(out)::counts(:),displacements(:),owners(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,first,last
    ok=norbital>0.and.nproc>0
    if(.not.ok)then;message='invalid balanced orbital ownership contract';return;end if
    allocate(counts(nproc),displacements(nproc),owners(norbital))
    counts=norbital/nproc
    counts(1:modulo(norbital,nproc))=counts(1:modulo(norbital,nproc))+1
    displacements(1)=0
    do rank=2,nproc;displacements(rank)=displacements(rank-1)+counts(rank-1);end do
    owners=-1
    do rank=0,nproc-1
      first=displacements(rank+1)+1;last=first+counts(rank+1)-1
      if(last>=first)owners(first:last)=rank
    end do
    ok=sum(counts)==norbital.and.maxval(counts)-minval(counts)<=1.and.all(owners>=0)
    if(ok)then;message='';else;message='balanced orbital ownership construction failed';end if
  end subroutine build_dg_balanced_orbital_ownership

  subroutine transpose_dg_spatial_cores_to_orbital_owners(comm,local_values,local_ids,batch_size,&
      owned_orbitals,global_ids,owned_values,ok,message)
    integer,intent(in)::comm,batch_size
    complex(real64),intent(in)::local_values(:,:)
    integer(int64),intent(in)::local_ids(:)
    integer,allocatable,intent(out)::owned_orbitals(:)
    integer(int64),allocatable,intent(out)::global_ids(:)
    complex(real64),allocatable,intent(out)::owned_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,allocatable::orbital_counts(:),orbital_displacements(:),orbital_owners(:),core_counts(:),&
      core_displacements(:),send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:),&
      batch_orbitals(:),local_batch_orbitals(:)
    complex(real64),allocatable::send_values(:),receive_values(:)
    integer::rank,nproc,ierr,norbital,nlocal,nglobal,batch_first,batch_last,destination,source,&
      orbital,point,position,local_row,global_point,local_bad,global_bad,total_core,total_send,total_receive,&
      min_norbital,max_norbital
    integer(int64)::count64

    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    norbital=size(local_values,1);nlocal=size(local_ids);local_bad=0
    if(norbital<1.or.nlocal<1.or.size(local_values,2)/=nlocal.or.batch_size<1.or.&
        any(local_ids<1_int64))local_bad=1
    call MPI_Allreduce(norbital,min_norbital,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(norbital,max_norbital,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(min_norbital/=max_norbital.or.ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='invalid spatial-to-orbital transpose contract';return
    end if
    call build_dg_balanced_orbital_ownership(norbital,nproc,orbital_counts,orbital_displacements,&
      orbital_owners,ok,message)
    if(.not.ok)return
    allocate(core_counts(nproc),core_displacements(nproc),send_counts(nproc),receive_counts(nproc),&
      send_displacements(nproc),receive_displacements(nproc))
    call MPI_Allgather(nlocal,1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    call build_checked_mpi_displacements(core_counts,core_displacements,total_core,ok)
    if(.not.ok)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='spatial core MPI displacement overflow';return
    end if
    nglobal=total_core;allocate(global_ids(nglobal))
    call MPI_Allgatherv(local_ids,nlocal,MPI_INTEGER8,global_ids,core_counts,core_displacements,&
      MPI_INTEGER8,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.ids_cover_unique_range(global_ids,nglobal))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      ok=.false.;message='spatial core IDs are not globally unique';return
    end if
    owned_orbitals=[(orbital,orbital=orbital_displacements(rank+1)+1,&
      orbital_displacements(rank+1)+orbital_counts(rank+1))]
    allocate(owned_values(size(owned_orbitals),nglobal));owned_values=(0d0,0d0)
    do batch_first=1,norbital,batch_size
      batch_last=min(norbital,batch_first+batch_size-1)
      batch_orbitals=[(orbital,orbital=batch_first,batch_last)]
      local_batch_orbitals=pack(batch_orbitals,orbital_owners(batch_orbitals)==rank)
      do destination=0,nproc-1
        count64=int(count(orbital_owners(batch_orbitals)==destination),int64)*int(nlocal,int64)
        if(count64>int(huge(1),int64))then
          local_bad=1;send_counts(destination+1)=0
        else
          send_counts(destination+1)=int(count64)
        end if
        count64=int(size(local_batch_orbitals),int64)*int(core_counts(destination+1),int64)
        if(count64>int(huge(1),int64))then
          local_bad=1;receive_counts(destination+1)=0
        else
          receive_counts(destination+1)=int(count64)
        end if
      end do
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
        ok=.false.;message='spatial-to-orbital MPI count overflow';return
      end if
      call build_checked_mpi_displacements(send_counts,send_displacements,total_send,ok)
      if(.not.ok)local_bad=1
      call build_checked_mpi_displacements(receive_counts,receive_displacements,total_receive,ok)
      if(.not.ok)local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
        ok=.false.;message='spatial-to-orbital MPI displacement overflow';return
      end if
      allocate(send_values(total_send),receive_values(total_receive))
      position=0
      do destination=0,nproc-1;do point=1,nlocal
        do orbital=batch_first,batch_last
          if(orbital_owners(orbital)/=destination)cycle
          position=position+1;send_values(position)=local_values(orbital,point)
        end do
      end do;end do
      call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,&
        receive_values,receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='spatial-to-orbital MPI Alltoallv failed';return;end if
      position=0
      do source=0,nproc-1;do point=1,core_counts(source+1)
        global_point=core_displacements(source+1)+point
        do orbital=1,size(local_batch_orbitals)
          position=position+1
          local_row=findloc(owned_orbitals,local_batch_orbitals(orbital),dim=1)
          owned_values(local_row,global_point)=receive_values(position)
        end do
      end do;end do
      deallocate(send_values,receive_values,batch_orbitals,local_batch_orbitals)
    end do
    ok=all(ieee_is_finite(real(owned_values))).and.all(ieee_is_finite(aimag(owned_values)))
    if(ok)then;message='';else;message='spatial-to-orbital transpose produced nonfinite values';end if
#else
    ok=.false.;message='spatial-to-orbital transpose requires MPI'
#endif
  contains
    logical function ids_cover_unique_range(ids,extent) result(valid)
      integer(int64),intent(in)::ids(:)
      integer,intent(in)::extent
      logical,allocatable::seen(:)
      integer::i
      valid=size(ids)==extent
      if(.not.valid)return
      if(any(ids<1_int64).or.any(ids>int(extent,int64)))then;valid=.false.;return;end if
      allocate(seen(extent));seen=.false.
      do i=1,size(ids)
        if(seen(int(ids(i))))then;valid=.false.;return;end if
        seen(int(ids(i)))=.true.
      end do
      valid=all(seen)
    end function ids_cover_unique_range
  end subroutine transpose_dg_spatial_cores_to_orbital_owners

  subroutine verify_dg_wannier_center_affine_orbits(centers,integer_rotations,&
      fractional_translations,tolerance,ok,message,moment_magnitudes,failed_operation)
    real(real64),intent(in)::centers(:,:),fractional_translations(:,:),tolerance
    integer,intent(in)::integer_rotations(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),intent(in),optional::moment_magnitudes(:,:)
    integer,intent(out),optional::failed_operation
    real(real64),allocatable::mapped_centers(:,:)
    integer,allocatable::matched_target(:)
    logical,allocatable::seen(:)
    integer::nwann,noperation,operation,source,target
    real(real64)::difference(3),nearest_residual,moment_min,moment_max

    if(present(failed_operation))failed_operation=0
    nwann=size(centers,2);noperation=size(integer_rotations,3)
    ok=nwann>0.and.size(centers,1)==3.and.noperation>0.and.size(integer_rotations,1)==3.and.&
      size(integer_rotations,2)==3.and.all(shape(fractional_translations)==[3,noperation]).and.&
      tolerance>0d0.and.all(ieee_is_finite(centers)).and.&
      all(ieee_is_finite(fractional_translations))
    if(ok.and.present(moment_magnitudes))then
      ok=all(shape(moment_magnitudes)==[3,nwann]).and.all(moment_magnitudes>=0d0).and.&
        all(ieee_is_finite(moment_magnitudes))
    endif
    if(.not.ok)then;message='invalid Wannier center affine-orbit contract';return;end if
    allocate(mapped_centers(3,nwann),matched_target(nwann),seen(nwann))
    do operation=1,noperation
      mapped_centers=modulo(matmul(real(integer_rotations(:,:,operation),real64),centers)+&
        spread(fractional_translations(:,operation),2,nwann),1d0)
      matched_target=0
      do source=1,nwann
        seen=.false.
        if(.not.augment_center_match(source))then
          nearest_residual=huge(1d0)
          do target=1,nwann
            difference=mapped_centers(:,source)-centers(:,target)
            difference=difference-anint(difference)
            nearest_residual=min(nearest_residual,maxval(abs(difference)))
          enddo
          if(present(moment_magnitudes))then
            moment_min=minval(moment_magnitudes(:,source))
            moment_max=maxval(moment_magnitudes(:,source))
          else
            moment_min=-1d0;moment_max=-1d0
          endif
          ok=.false.
          if(present(failed_operation))failed_operation=operation
          write(message,'(a,i0,a,i0,4(a,es12.4))')'localized Wannier center orbit mismatch operation=',&
            operation,' source=',source,' nearest_residual=',nearest_residual,' tolerance=',tolerance,&
            ' moment_min=',moment_min,' moment_max=',moment_max
          return
        end if
      end do
    end do
    ok=.true.;message=''
  contains
    recursive logical function augment_center_match(source_index) result(found)
      integer,intent(in)::source_index
      integer::target,previous_source
      real(real64)::difference(3)
      found=.false.
      do target=1,nwann
        if(seen(target))cycle
        difference=mapped_centers(:,source_index)-centers(:,target)
        difference=difference-anint(difference)
        if(maxval(abs(difference))>tolerance)cycle
        seen(target)=.true.;previous_source=matched_target(target)
        if(previous_source==0)then
          matched_target(target)=source_index;found=.true.;return
        end if
        if(augment_center_match(previous_source))then
          matched_target(target)=source_index;found=.true.;return
        end if
      end do
    end function augment_center_match
  end subroutine verify_dg_wannier_center_affine_orbits

  subroutine diagnose_dg_point_center_gauge(comm,local_basis,weights,point_map,&
      integer_rotation,fractional_translation,centers,tolerance,monomial_defect,&
      center_block_leakage,unitarity_defect,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,integer_rotation(:,:)
    complex(real64),intent(in)::local_basis(:,:)
    real(real64),intent(in)::weights(:),fractional_translation(:),centers(:,:),tolerance
    integer(int64),intent(in)::point_map(:)
    real(real64),intent(out)::monomial_defect,center_block_leakage,unitarity_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,nstate,nlocal,local_bad,global_bad,status,i,j,owner,base,remainder,&
      owner_count,owner_first,target,source
    integer(int64),allocatable::row_ids(:),map_rows(:,:),global_map(:)
    complex(real64),allocatable::rows(:,:,:),remote_rows(:,:),unitarity_tile(:,:)
    real(real64),allocatable::column_max(:),global_column_max(:),local_leakage(:),global_leakage(:)
    integer,allocatable::map_count(:)
    real(real64)::mapped(3),difference(3),row_max,local_monomial,local_unitarity,expected
    integer(int64)::overlap_workspace,complex_bytes,real_bytes,integer_bytes,extra_bytes,term
    nstate=size(local_basis,1);nlocal=size(local_basis,2)
    ok=.false.;message='';monomial_defect=huge(1d0);center_block_leakage=huge(1d0)
    unitarity_defect=huge(1d0);workspace_peak_bytes=0_int64
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.nstate>0.and.nlocal>0.and.size(weights)==nlocal.and.&
      size(point_map)==nlocal.and.all(shape(integer_rotation)==[3,3]).and.&
      size(fractional_translation)==3.and.all(shape(centers)==[3,nstate]).and.&
      tolerance>0d0.and.ieee_is_finite(tolerance).and.all(ieee_is_finite(weights)).and.&
      all(weights>=0d0).and.all(ieee_is_finite(real(local_basis))).and.&
      all(ieee_is_finite(aimag(local_basis))).and.all(ieee_is_finite(fractional_translation)).and.&
      all(ieee_is_finite(centers)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid point center-gauge diagnostic contract';return;endif
    call agree_integer(nstate);if(global_bad/=0)return
    call agree_integer(nlocal);if(global_bad/=0)return
    do i=1,3;do j=1,3;call agree_integer(integer_rotation(i,j));if(global_bad/=0)return;enddo;enddo
    allocate(global_map(nlocal*nproc),map_count(nlocal*nproc),map_rows(nlocal,1),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='point center-gauge map allocation failed';return;endif
    call MPI_Allgather(point_map,nlocal,MPI_INTEGER8,global_map,nlocal,MPI_INTEGER8,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.all(global_map>=1_int64).and.&
      all(global_map<=int(nlocal*nproc,int64)))
    map_count=0
    if(local_bad==0)then
      do i=1,size(global_map);map_count(int(global_map(i)))=map_count(int(global_map(i)))+1;enddo
      if(any(map_count/=1))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='point center-gauge map is not a permutation';return;endif
    map_rows(:,1)=point_map
    call assemble_dg_distributed_basis_symmetry_overlap_rows(comm,local_basis,weights,map_rows,&
      row_ids,rows,overlap_workspace,ok,message)
    if(.not.ok)then;call cleanup();return;endif
    allocate(column_max(nstate),global_column_max(nstate),local_leakage(nstate),global_leakage(nstate),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='point center-gauge receipt allocation failed';return;endif
    column_max=0d0;local_leakage=0d0;local_monomial=0d0
    do i=1,size(row_ids)
      target=int(row_ids(i));row_max=0d0
      do source=1,nstate
        row_max=max(row_max,abs(rows(i,source,1))**2)
        column_max(source)=max(column_max(source),abs(rows(i,source,1))**2)
        mapped=modulo(matmul(real(integer_rotation,real64),centers(:,source))+fractional_translation,1d0)
        difference=mapped-centers(:,target);difference=difference-anint(difference)
        if(maxval(abs(difference))>tolerance)&
          local_leakage(source)=local_leakage(source)+abs(rows(i,source,1))**2
      enddo
      local_monomial=max(local_monomial,abs(1d0-row_max))
    enddo
    call MPI_Allreduce(column_max,global_column_max,nstate,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_leakage,global_leakage,nstate,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    monomial_defect=max(local_monomial,maxval(abs(1d0-global_column_max)))
    call MPI_Allreduce(MPI_IN_PLACE,monomial_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    center_block_leakage=maxval(global_leakage)
    base=nstate/nproc;remainder=mod(nstate,nproc);local_unitarity=0d0
    do owner=0,nproc-1
      owner_count=base+merge(1,0,owner<remainder);owner_first=owner*base+min(owner,remainder)+1
      allocate(remote_rows(owner_count,nstate),unitarity_tile(size(row_ids),owner_count),stat=status)
      call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='point center-gauge unitarity allocation failed';return;endif
      if(rank==owner)remote_rows=rows(:,:,1)
      call MPI_Bcast(remote_rows,size(remote_rows),MPI_DOUBLE_COMPLEX,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='point center-gauge row broadcast failed';return;endif
      unitarity_tile=matmul(rows(:,:,1),conjg(transpose(remote_rows)))
      do j=1,owner_count;do i=1,size(row_ids)
        expected=merge(1d0,0d0,int(row_ids(i))==owner_first+j-1)
        local_unitarity=max(local_unitarity,abs(unitarity_tile(i,j)-expected))
      enddo;enddo
      deallocate(remote_rows,unitarity_tile)
    enddo
    call MPI_Allreduce(local_unitarity,unitarity_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    complex_bytes=16_int64;real_bytes=8_int64;integer_bytes=8_int64
    extra_bytes=integer_bytes*int(size(global_map)+size(row_ids),int64)
    term=4_int64*real_bytes*int(nstate,int64)
    if(extra_bytes>huge(0_int64)-term)then;call cleanup();message='point center-gauge receipt overflows';return;endif
    workspace_peak_bytes=max(overlap_workspace,extra_bytes+term+complex_bytes*int(size(rows),int64))
    call MPI_Allreduce(MPI_IN_PLACE,workspace_peak_bytes,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.unitarity_defect<=tolerance
    if(ok)then;message='';else;message='point center-gauge representation is not unitary';endif
    call cleanup()
  contains
    subroutine agree_integer(value)
      integer,intent(in)::value
      integer::minimum,maximum
      call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;global_bad=1;message='point center-gauge metadata reduction failed';return;endif
      call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minimum/=maximum)
      if(global_bad/=0)message='point center-gauge metadata disagree across ranks'
    end subroutine
    subroutine cleanup()
      if(allocated(row_ids))deallocate(row_ids)
      if(allocated(map_rows))deallocate(map_rows)
      if(allocated(global_map))deallocate(global_map)
      if(allocated(map_count))deallocate(map_count)
      if(allocated(rows))deallocate(rows)
      if(allocated(remote_rows))deallocate(remote_rows)
      if(allocated(unitarity_tile))deallocate(unitarity_tile)
      if(allocated(column_max))deallocate(column_max)
      if(allocated(global_column_max))deallocate(global_column_max)
      if(allocated(local_leakage))deallocate(local_leakage)
      if(allocated(global_leakage))deallocate(global_leakage)
    end subroutine
#else
    ok=.false.;message='point center-gauge diagnostic requires MPI'
    monomial_defect=huge(1d0);center_block_leakage=huge(1d0);unitarity_defect=huge(1d0)
    workspace_peak_bytes=0_int64
#endif
  end subroutine diagnose_dg_point_center_gauge

  subroutine build_dg_finite_abelian_character_table(translations,product_table,identity_operation,&
      tolerance,canonical_operations,inverse_operations,generator_count,generators,element_words,&
      characters,conjugate_characters,fingerprint,ok,message)
    real(real64),intent(in)::translations(:,:),tolerance
    integer,intent(in)::product_table(:,:),identity_operation
    integer,intent(out)::canonical_operations(:),inverse_operations(:),generator_count
    integer,allocatable,intent(out)::generators(:),element_words(:,:)
    complex(real64),intent(out)::characters(:,:)
    integer,intent(out)::conjugate_characters(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::n,i,j,k,g,candidate,changed,target,head,tail,character_count,allocation_status,&
      greedy_count,desired_count
    integer,allocatable::input_to_canonical(:),generator_buffer(:),generator_trial(:),&
      generator_best(:),generator_orders(:),&
      phase_index(:),queue(:)
    integer(int64),allocatable::translation_key(:,:)
    logical,allocatable::reached(:),word_known(:)
    logical::generator_subset_found
    complex(real64),allocatable::trial_character(:)
    real(real64)::delta(3),angle,quantum
    integer(int64),parameter::generator_subset_search_budget=10000000_int64
    integer(int64)::quantized,extent64,combination_count,generator_subset_trials

    ok=.false.;message='';fingerprint=0_int64;generator_count=0
    n=size(translations,2)
    if(n<=0)then;message='finite translation group is empty';return;endif
    if(tolerance<16d0*acos(-1d0)/real(huge(0_int64),real64))then
      message='finite translation tolerance is too small for canonical int64 keys';return
    endif
    if(size(translations,1)/=3.or.any(shape(product_table)/=[n,n]).or.&
        size(canonical_operations)/=n.or.size(inverse_operations)/=n.or.&
        any(shape(characters)/=[n,n]).or.size(conjugate_characters)/=n.or.&
        identity_operation<1.or.identity_operation>n.or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(translations)).or.any(product_table<1).or.any(product_table>n))then
      message='invalid finite-abelian character-table contract';return
    endif
    characters=(0d0,0d0);inverse_operations=0;conjugate_characters=0
    quantum=tolerance/4d0
    allocate(input_to_canonical(n),generator_buffer(n),generator_trial(n),generator_best(n),&
      reached(n),translation_key(3,n),stat=allocation_status)
    if(allocation_status/=0)then;message='finite-group canonical metadata allocation failed';return;endif
    do i=1,n;do j=1,3
      delta(j)=modulo(translations(j,i),1d0)
      if(abs(delta(j)-1d0)<=tolerance.or.abs(delta(j))<=tolerance)delta(j)=0d0
      translation_key(j,i)=nint(delta(j)/quantum,int64)
    enddo;enddo
    canonical_operations=[(i,i=1,n)]
    do i=1,n-1;do j=i+1,n
      if(key_less(canonical_operations(j),canonical_operations(i)))then
        k=canonical_operations(i);canonical_operations(i)=canonical_operations(j);canonical_operations(j)=k
      endif
    enddo;enddo
    do i=1,n;input_to_canonical(canonical_operations(i))=i;enddo
    if(canonical_operations(1)/=identity_operation)then
      message='designated identity is not the zero translation';return
    endif
    do i=1,n-1;do j=i+1,n
      delta=translations(:,i)-translations(:,j);delta=delta-anint(delta)
      if(maxval(abs(delta))<=tolerance)then;message='finite translation catalog is nonfaithful';return;endif
    enddo;enddo
    do i=1,n
      if(canonical_product(1,i)/=i.or.canonical_product(i,1)/=i)then
        message='designated translation identity violates the product table';return
      endif
      do j=1,n
        if(canonical_product(i,j)/=canonical_product(j,i))then
          message='finite translation group is nonabelian';return
        endif
        delta=translations(:,canonical_operations(i))+translations(:,canonical_operations(j))-&
          translations(:,canonical_operations(canonical_product(i,j)))
        delta=delta-anint(delta)
        if(maxval(abs(delta))>tolerance)then
          message='translation product table disagrees with geometry';return
        endif
        do k=1,n
          if(canonical_product(canonical_product(i,j),k)/=&
              canonical_product(i,canonical_product(j,k)))then
            message='finite translation product table is nonassociative';return
          endif
        enddo
      enddo
      do j=1,n
        if(canonical_product(i,j)==1.and.canonical_product(j,i)==1)inverse_operations(i)=j
      enddo
      if(inverse_operations(i)==0)then;message='finite translation element has no inverse';return;endif
    enddo
    reached=.false.;reached(1)=.true.;generator_count=0
    do candidate=2,n
      if(reached(candidate))cycle
      generator_count=generator_count+1;generator_buffer(generator_count)=candidate
      changed=1
      do while(changed/=0)
        changed=0
        do i=1,n
          if(.not.reached(i))cycle
          target=canonical_product(i,candidate)
          if(.not.reached(target))then;reached(target)=.true.;changed=1;endif
        enddo
      enddo
    enddo
    greedy_count=generator_count;generator_subset_found=.false.;generator_subset_trials=0_int64
    do desired_count=1,greedy_count
      call search_generator_subsets(1,2,desired_count)
      if(generator_subset_found)then;generator_count=desired_count;generator_buffer(:generator_count)=&
        generator_best(:generator_count);exit;endif
      if(generator_subset_trials>=generator_subset_search_budget)then
        message='minimum generator subset search exceeds supported deterministic budget';return
      endif
    enddo
    reached=.false.;reached(1)=.true.;changed=1
    do while(changed/=0)
      changed=0
      do i=1,n
        if(.not.reached(i))cycle
        do g=1,generator_count
          target=canonical_product(i,generator_buffer(g))
          if(.not.reached(target))then;reached(target)=.true.;changed=1;endif
        enddo
      enddo
    enddo
    if(.not.all(reached))then;message='canonical generators do not span translation group';return;endif
    extent64=int(n,int64)*int(generator_count,int64)
    if(extent64>int(huge(0),int64))then;message='finite-group word extent overflows';return;endif
    allocate(generators(generator_count),element_words(n,generator_count),generator_orders(generator_count),&
      phase_index(generator_count),queue(n),word_known(n),&
      trial_character(n),stat=allocation_status)
    if(allocation_status/=0)then;message='finite-group word allocation failed';return;endif
    generators=generator_buffer(:generator_count);element_words=0
    word_known=.false.;word_known(1)=.true.;head=1;tail=1;queue(1)=1
    do while(head<=tail)
      i=queue(head);head=head+1
      do g=1,generator_count
        target=canonical_product(i,generators(g))
        if(.not.word_known(target))then
          element_words(target,:)=element_words(i,:)
          if(element_words(target,g)==huge(0))then;message='generator word exponent overflows';return;endif
          element_words(target,g)=element_words(target,g)+1
          word_known(target)=.true.;tail=tail+1;queue(tail)=target
        endif
      enddo
    enddo
    if(.not.all(word_known))then;message='generator words do not cover translation group';return;endif
    combination_count=1_int64
    do g=1,generator_count
      generator_orders(g)=1;target=generators(g)
      do while(target/=1.and.generator_orders(g)<n)
        target=canonical_product(target,generators(g));generator_orders(g)=generator_orders(g)+1
      enddo
      if(target/=1)then;message='generator order exceeds finite group';return;endif
      if(combination_count>huge(combination_count)/int(generator_orders(g),int64))then
        message='character phase enumeration count overflows';return
      endif
      combination_count=combination_count*int(generator_orders(g),int64)
    enddo
    phase_index=0;character_count=0
    call enumerate_phases(1)
    if(character_count/=n)then;message='finite-group character enumeration is incomplete';return;endif
    do i=1,n
      do j=1,n
        if(maxval(abs(characters(j,:)-conjg(characters(i,:))))<=100d0*tolerance)&
          conjugate_characters(i)=j
      enddo
      if(conjugate_characters(i)==0)then;message='character has no conjugate partner';return;endif
    enddo
    fingerprint=int(z'6A09E667F3BCC909',int64)
    do i=1,n;do j=1,3
      fingerprint=ieor(fingerprint,ishftc(translation_key(j,canonical_operations(i)),mod(5*i+7*j,63)))
    enddo;enddo
    do i=1,n;do j=1,n
      angle=modulo(atan2(aimag(characters(i,j)),real(characters(i,j))),2d0*acos(-1d0))
      quantized=nint(angle/quantum,int64)
      fingerprint=ieor(fingerprint,ishftc(quantized,mod(11*i+13*j,63)))
    enddo;enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.;message=''
  contains
    integer function canonical_product(left,right) result(product)
      integer,intent(in)::left,right
      product=input_to_canonical(product_table(canonical_operations(left),canonical_operations(right)))
    end function
    logical function key_less(left,right) result(less)
      integer,intent(in)::left,right
      integer::axis
      less=.false.
      do axis=1,3
        if(translation_key(axis,left)<translation_key(axis,right))then;less=.true.;return;endif
        if(translation_key(axis,left)>translation_key(axis,right))return
      enddo
      less=left<right
    end function
    recursive subroutine enumerate_phases(level)
      integer,intent(in)::level
      integer::choice,element,relation_left,relation_right,relation_product,existing
      complex(real64)::value,root
      logical::valid,duplicate
      if(character_count>=n)return
      if(level<=generator_count)then
        do choice=0,generator_orders(level)-1
          phase_index(level)=choice;call enumerate_phases(level+1)
        enddo
        return
      endif
      do element=1,n
        value=(1d0,0d0)
        do g=1,generator_count
          root=exp(cmplx(0d0,2d0*acos(-1d0)*real(phase_index(g),real64)/&
            real(generator_orders(g),real64),real64))
          value=value*root**element_words(element,g)
        enddo
        trial_character(element)=value
      enddo
      valid=.true.
      do relation_left=1,n;do relation_right=1,n
        relation_product=canonical_product(relation_left,relation_right)
        if(abs(trial_character(relation_product)-trial_character(relation_left)*&
            trial_character(relation_right))>100d0*tolerance)valid=.false.
      enddo;enddo
      if(.not.valid)return
      duplicate=.false.
      do existing=1,character_count
        if(maxval(abs(characters(existing,:)-trial_character))<=100d0*tolerance)duplicate=.true.
      enddo
      if(duplicate)return
      character_count=character_count+1;characters(character_count,:)=trial_character
    end subroutine
    recursive subroutine search_generator_subsets(level,start,needed)
      integer,intent(in)::level,start,needed
      integer::choice
      if(generator_subset_found.or.generator_subset_trials>=generator_subset_search_budget)return
      if(level>needed)then
        if(generator_subset_trials>=generator_subset_search_budget)return
        generator_subset_trials=generator_subset_trials+1_int64
        if(subset_spans(needed))then
          generator_best(:needed)=generator_trial(:needed);generator_subset_found=.true.
        endif
        return
      endif
      do choice=start,n-(needed-level)
        generator_trial(level)=choice
        call search_generator_subsets(level+1,choice+1,needed)
        if(generator_subset_found.or.generator_subset_trials>=generator_subset_search_budget)return
      enddo
    end subroutine
    logical function subset_spans(count) result(spans)
      integer,intent(in)::count
      integer::element,igen,mapped,progress
      reached=.false.;reached(1)=.true.;progress=1
      do while(progress/=0)
        progress=0
        do element=1,n
          if(.not.reached(element))cycle
          do igen=1,count
            mapped=canonical_product(element,generator_trial(igen))
            if(.not.reached(mapped))then;reached(mapped)=.true.;progress=1;endif
          enddo
        enddo
      enddo
      spans=all(reached)
    end function
  end subroutine build_dg_finite_abelian_character_table

  subroutine compute_dg_periodic_wannier_centers(comm,values,weights,periodic_phases,centers,&
      moment_magnitudes,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::values(:,:),periodic_phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::centers(:,:),moment_magnitudes(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::local_moments(:,:),global_moments(:,:)
    real(real64),allocatable::local_norm(:),global_norm(:)
    integer::nwann,npoint,iw,axis,ierr
    real(real64),parameter::two_pi=2d0*acos(-1d0)
    nwann=size(values,1);npoint=size(values,2)
    ok=nwann>0.and.npoint>0.and.size(weights)==npoint.and.all(shape(periodic_phases)==[3,npoint]).and.&
      all(shape(centers)==[3,nwann]).and.all(shape(moment_magnitudes)==[3,nwann]).and.&
      all(weights>=0d0).and.all(ieee_is_finite(weights))
    if(.not.ok)then;message='invalid periodic Wannier center contract';return;end if
    allocate(local_moments(3,nwann),global_moments(3,nwann),local_norm(nwann),global_norm(nwann))
    do iw=1,nwann
      local_norm(iw)=sum(weights*abs(values(iw,:))**2)
      do axis=1,3
        local_moments(axis,iw)=sum(weights*abs(values(iw,:))**2*periodic_phases(axis,:))
      end do
    end do
    call MPI_Allreduce(local_norm,global_norm,nwann,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moments,global_moments,3*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_norm<=epsilon(1d0)))then
      ok=.false.;message='periodic Wannier center normalization failed';return
    end if
    do iw=1,nwann;do axis=1,3
      moment_magnitudes(axis,iw)=abs(global_moments(axis,iw))/global_norm(iw)
      centers(axis,iw)=modulo(atan2(aimag(global_moments(axis,iw)),&
        real(global_moments(axis,iw)))/two_pi,1d0)
    end do;end do
    ok=all(ieee_is_finite(centers)).and.all(ieee_is_finite(moment_magnitudes))
    if(ok)then;message='';else;message='periodic Wannier center is not finite';end if
#else
    ok=.false.;message='periodic Wannier center measurement requires MPI'
#endif
  end subroutine compute_dg_periodic_wannier_centers

  subroutine solve_dg_affine_common_fixed_point(integer_rotations,fractional_translations,tolerance,&
      has_common_center,center,maximum_residual,ok,message)
    integer,intent(in)::integer_rotations(:,:,:)
    real(real64),intent(in)::fractional_translations(:,:),tolerance
    logical,intent(out)::has_common_center,ok
    real(real64),intent(out)::center(3),maximum_residual
    character(*),intent(out)::message
    complex(real64)::normal_complex(3,3)
    complex(real64),allocatable::eigenvectors(:,:)
    real(real64),allocatable::eigenvalues(:)
    real(real64)::normal(3,3),rhs(3),trial(3),updated(3),a(3,3),delta(3),residual,&
      best_residual,best_norm,trial_norm,best_lex,trial_lex
    integer::sx,sy,sz,iteration,operation,i,j,k,noperation
    logical::eigen_ok
    character(256)::detail

    noperation=size(integer_rotations,3);center=0d0;maximum_residual=huge(1d0)
    has_common_center=.false.;ok=noperation>0.and.size(integer_rotations,1)==3.and.&
      size(integer_rotations,2)==3.and.all(shape(fractional_translations)==[3,noperation]).and.&
      tolerance>0d0.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(fractional_translations))
    if(.not.ok)then;message='invalid affine common-center contract';return;end if
    best_residual=huge(1d0);best_norm=huge(1d0);best_lex=huge(1d0)
    do sx=0,3;do sy=0,3;do sz=0,3
      trial=0.25d0*[real(sx,real64),real(sy,real64),real(sz,real64)]
      do iteration=1,16
        normal=0d0;rhs=0d0
        do operation=1,noperation
          a=0d0;do i=1,3;a(i,i)=1d0;end do
          a=a-real(integer_rotations(:,:,operation),real64)
          delta=matmul(a,trial)-fractional_translations(:,operation)
          normal=normal+matmul(transpose(a),a)
          rhs=rhs+matmul(transpose(a),fractional_translations(:,operation)+anint(delta))
        end do
        normal_complex=cmplx(normal,0d0,real64)
        call hermitian_eigensystem(normal_complex,eigenvalues,eigenvectors,eigen_ok,detail)
        if(.not.eigen_ok)then;ok=.false.;message='affine common-center normal solve failed';return;end if
        updated=0d0
        do k=1,3
          if(eigenvalues(k)<=tolerance*max(1d0,maxval(eigenvalues)))cycle
          updated=updated+real(eigenvectors(:,k),real64)*&
            dot_product(real(eigenvectors(:,k),real64),rhs)/eigenvalues(k)
        end do
        updated=modulo(updated,1d0)
        if(maxval(abs(modulo(updated-trial+0.5d0,1d0)-0.5d0))<=tolerance)exit
        trial=updated
      end do
      residual=0d0
      do operation=1,noperation
        a=0d0;do i=1,3;a(i,i)=1d0;end do
        a=a-real(integer_rotations(:,:,operation),real64)
        delta=matmul(a,updated)-fractional_translations(:,operation)
        residual=max(residual,maxval(abs(delta-anint(delta))))
      end do
      trial_norm=sum(min(updated,1d0-updated)**2)
      trial_lex=updated(1)+1d-3*updated(2)+1d-6*updated(3)
      if(residual<best_residual-tolerance.or.&
          (abs(residual-best_residual)<=tolerance.and.(trial_norm<best_norm-tolerance.or.&
          (abs(trial_norm-best_norm)<=tolerance.and.trial_lex<best_lex))))then
        best_residual=residual;best_norm=trial_norm;best_lex=trial_lex;center=updated
      end if
    end do;end do;end do
    maximum_residual=best_residual;has_common_center=best_residual<=tolerance
    ok=.true.;message=''
  end subroutine solve_dg_affine_common_fixed_point

  subroutine accept_dg_boundary_calibrated_symmetry(boundary_residual,interior_residual,&
      boundary_allowance,interior_tolerance,ok,message)
    real(real64),intent(in)::boundary_residual(:),interior_residual(:),boundary_allowance,interior_tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    ok=size(boundary_residual)>0.and.size(interior_residual)==size(boundary_residual).and.&
      boundary_allowance>=0d0.and.interior_tolerance>0d0.and.ieee_is_finite(boundary_allowance).and.&
      ieee_is_finite(interior_tolerance).and.all(ieee_is_finite(boundary_residual)).and.&
      all(ieee_is_finite(interior_residual)).and.all(boundary_residual>=0d0).and.&
      all(interior_residual>=0d0)
    if(.not.ok)then;message='invalid boundary-calibrated symmetry gate';return;end if
    if(maxval(interior_residual)>interior_tolerance)then
      ok=.false.;message='LCFO interior symmetry residual exceeds strict tolerance';return
    end if
    if(maxval(boundary_residual)>boundary_allowance)then
      ok=.false.;message='LCFO boundary symmetry residual exceeds measured stitching allowance';return
    end if
    message=''
  end subroutine accept_dg_boundary_calibrated_symmetry

  subroutine measure_dg_rank_fixed_symmetry_residuals(comm,basis,weights,symmetry_target_box_ids,&
      boundary_mask,representation,total_residual,boundary_residual,interior_residual,ok,message,&
      workspace_peak_bytes)
    integer,intent(in)::comm
    complex(real64),intent(in)::basis(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    logical,intent(in)::boundary_mask(:)
    complex(real64),intent(out),optional::representation(:,:,:)
    real(real64),intent(out)::total_residual(:),boundary_residual(:),interior_residual(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(out),optional::workspace_peak_bytes
#ifdef USE_MPI
    integer,parameter::orbital_tile_size=32
    complex(real64),allocatable::local_metric(:,:),metric(:,:),metric_vectors(:,:),metric_inverse_sqrt(:,:),&
      orthonormal_basis(:,:),image_tile(:,:),local_overlap(:,:),global_overlap(:,:)
    real(real64),allocatable::metric_spectrum(:),local_norms(:),global_norms(:)
    integer::ierr,nstate,nlocal,noperation,operation,i,tile_first,tile_count
    integer(int64)::tile_bytes,peak_bytes,base_workspace_bytes,complex_bytes,real_bytes
    logical::eigen_ok,representation_shape_ok
    character(256)::detail

    nstate=size(basis,1);nlocal=size(basis,2);noperation=size(symmetry_target_box_ids,2)
    representation_shape_ok=.true.
    if(present(representation))representation_shape_ok=all(shape(representation)==[nstate,nstate,noperation])
    ok=nstate>0.and.nlocal>0.and.noperation>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.size(boundary_mask)==nlocal.and.&
      representation_shape_ok.and.&
      size(total_residual)==noperation.and.size(boundary_residual)==noperation.and.&
      size(interior_residual)==noperation.and.all(weights>=0d0)
    peak_bytes=0_int64;if(present(workspace_peak_bytes))workspace_peak_bytes=0_int64
    if(.not.ok)then;message='invalid rank-fixed symmetry residual contract';return;end if
    allocate(local_metric(nstate,nstate),metric(nstate,nstate),metric_inverse_sqrt(nstate,nstate),&
      orthonormal_basis(nstate,nlocal),local_overlap(nstate,nstate),global_overlap(nstate,nstate),&
      local_norms(3),global_norms(3))
    do i=1,nstate
      local_metric(i,:)=matmul(conjg(basis(i,:))*weights,transpose(basis))
    end do
    call MPI_Allreduce(local_metric,metric,nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call hermitian_eigensystem(metric,metric_spectrum,metric_vectors,eigen_ok,detail)
    if(.not.eigen_ok.or.minval(metric_spectrum)<=epsilon(1d0)*maxval(metric_spectrum))then
      ok=.false.;message='rank-fixed occupied metric is singular';return
    end if
    metric_inverse_sqrt=metric_vectors
    do i=1,nstate;metric_inverse_sqrt(:,i)=metric_inverse_sqrt(:,i)/sqrt(metric_spectrum(i));end do
    metric_inverse_sqrt=matmul(metric_inverse_sqrt,conjg(transpose(metric_vectors)))
    orthonormal_basis=matmul(metric_inverse_sqrt,basis)
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    real_bytes=int(storage_size(0d0)/8,int64)
    base_workspace_bytes=complex_bytes*int(size(local_metric)+size(metric)+size(metric_vectors)+&
      size(metric_inverse_sqrt)+size(orthonormal_basis)+size(local_overlap)+size(global_overlap),int64)+&
      real_bytes*int(size(metric_spectrum)+size(local_norms)+size(global_norms),int64)
    peak_bytes=base_workspace_bytes
    do operation=1,noperation
      local_overlap=(0d0,0d0)
      do tile_first=1,nstate,orbital_tile_size
        tile_count=min(orbital_tile_size,nstate-tile_first+1)
        allocate(image_tile(tile_count,nlocal))
        tile_bytes=int(storage_size((0d0,0d0))/8,int64)*int(size(image_tile),int64)
        peak_bytes=max(peak_bytes,base_workspace_bytes+tile_bytes)
        call exchange_dg_point_permuted_orbital_rows(comm,&
          orthonormal_basis(tile_first:tile_first+tile_count-1,:),&
          symmetry_target_box_ids(:,operation),image_tile,ok,message)
        if(.not.ok)return
        do i=1,nstate
          local_overlap(i,tile_first:tile_first+tile_count-1)=&
            matmul(conjg(orthonormal_basis(i,:))*weights,transpose(image_tile))
        enddo
        deallocate(image_tile)
      enddo
      call MPI_Allreduce(local_overlap,global_overlap,nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(present(representation))representation(:,:,operation)=global_overlap
      local_norms=0d0
      do tile_first=1,nstate,orbital_tile_size
        tile_count=min(orbital_tile_size,nstate-tile_first+1)
        allocate(image_tile(tile_count,nlocal))
        tile_bytes=int(storage_size((0d0,0d0))/8,int64)*int(size(image_tile),int64)
        peak_bytes=max(peak_bytes,base_workspace_bytes+tile_bytes)
        call exchange_dg_point_permuted_orbital_rows(comm,&
          orthonormal_basis(tile_first:tile_first+tile_count-1,:),&
          symmetry_target_box_ids(:,operation),image_tile,ok,message)
        if(.not.ok)return
        image_tile=image_tile-matmul(transpose(global_overlap(:,tile_first:tile_first+tile_count-1)),&
          orthonormal_basis)
        local_norms(1)=local_norms(1)+sum(spread(weights,1,tile_count)*abs(image_tile)**2)
        local_norms(2)=local_norms(2)+sum(spread(weights*merge(1d0,0d0,boundary_mask),1,tile_count)*&
          abs(image_tile)**2)
        local_norms(3)=local_norms(3)+sum(spread(weights*merge(0d0,1d0,boundary_mask),1,tile_count)*&
          abs(image_tile)**2)
        deallocate(image_tile)
      enddo
      call MPI_Allreduce(local_norms,global_norms,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      total_residual(operation)=sqrt(max(0d0,global_norms(1)))
      boundary_residual(operation)=sqrt(max(0d0,global_norms(2)))
      interior_residual(operation)=sqrt(max(0d0,global_norms(3)))
    end do
    ok=all(ieee_is_finite(total_residual)).and.all(ieee_is_finite(boundary_residual)).and.&
      all(ieee_is_finite(interior_residual))
    if(present(workspace_peak_bytes))workspace_peak_bytes=peak_bytes
    if(ok)then;message='';else;message='rank-fixed symmetry residual is not finite';end if
#else
    ok=.false.;message='rank-fixed symmetry residual measurement requires MPI'
    if(present(workspace_peak_bytes))workspace_peak_bytes=0_int64
#endif
  end subroutine measure_dg_rank_fixed_symmetry_residuals

  subroutine exchange_dg_point_permuted_orbital_rows(comm,basis,target_global_ids,image,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::basis(:,:)
    integer(int64),intent(in)::target_global_ids(:)
    complex(real64),intent(out)::image(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,nstate,nlocal,nlocal_min,nlocal_max,local_bad,global_bad
    integer::point,owner,index,cursor,total_send,total_recv
    integer(int64)::scaled_count
    logical::counts_ok
    integer,allocatable::send_counts(:),recv_counts(:),send_displs(:),recv_displs(:),fill(:)
    integer,allocatable::request_indices(:),request_destinations(:),received_requests(:)
    integer,allocatable::target_hits(:)
    integer,allocatable::value_send_counts(:),value_recv_counts(:),value_send_displs(:),value_recv_displs(:)
    complex(real64),allocatable::response_send(:,:),response_recv(:,:)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nstate=size(basis,1);nlocal=size(basis,2);ok=.false.;message=''
    call MPI_Allreduce(nlocal,nlocal_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nlocal,nlocal_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,nstate>0.and.nlocal>0.and.nlocal_min==nlocal_max.and.&
      size(target_global_ids)==nlocal.and.all(shape(image)==[nstate,nlocal]).and.&
      all(target_global_ids>=1_int64).and.&
      all(target_global_ids<=int(nproc,int64)*int(nlocal,int64)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='invalid distributed symmetry point target';return
    endif
    allocate(send_counts(nproc),recv_counts(nproc),send_displs(nproc),recv_displs(nproc),fill(nproc))
    send_counts=0
    do point=1,nlocal
      owner=int((target_global_ids(point)-1_int64)/int(nlocal,int64))
      send_counts(owner+1)=send_counts(owner+1)+1
    enddo
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr)
    call build_checked_mpi_displacements(send_counts,send_displs,total_send,counts_ok)
    if(counts_ok)call build_checked_mpi_displacements(recv_counts,recv_displs,total_recv,counts_ok)
    local_bad=merge(0,1,counts_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='distributed symmetry request count overflow';return
    endif
    allocate(request_indices(total_send),request_destinations(total_send),received_requests(total_recv))
    fill=send_displs
    do point=1,nlocal
      owner=int((target_global_ids(point)-1_int64)/int(nlocal,int64))+1
      cursor=fill(owner)+1;fill(owner)=cursor
      request_indices(cursor)=int(modulo(target_global_ids(point)-1_int64,int(nlocal,int64)))+1
      request_destinations(cursor)=point
    enddo
    call MPI_Alltoallv(request_indices,send_counts,send_displs,MPI_INTEGER,&
      received_requests,recv_counts,recv_displs,MPI_INTEGER,comm,ierr)
    allocate(target_hits(nlocal));target_hits=0
    do index=1,total_recv
      if(received_requests(index)>=1.and.received_requests(index)<=nlocal)&
        target_hits(received_requests(index))=target_hits(received_requests(index))+1
    enddo
    local_bad=merge(0,1,all(received_requests>=1).and.all(received_requests<=nlocal).and.&
      all(target_hits==1))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='distributed symmetry targets are not a complete point permutation';return
    endif
    allocate(response_send(nstate,total_recv),response_recv(nstate,total_send))
    do index=1,total_recv
      response_send(:,index)=basis(:,received_requests(index))
    enddo
    allocate(value_send_counts(nproc),value_recv_counts(nproc),value_send_displs(nproc),value_recv_displs(nproc))
    local_bad=0
    do owner=1,nproc
      scaled_count=int(nstate,int64)*int(recv_counts(owner),int64)
      if(scaled_count>int(huge(0),int64))local_bad=1
      value_send_counts(owner)=int(min(scaled_count,int(huge(0),int64)))
      scaled_count=int(nstate,int64)*int(send_counts(owner),int64)
      if(scaled_count>int(huge(0),int64))local_bad=1
      value_recv_counts(owner)=int(min(scaled_count,int(huge(0),int64)))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='distributed symmetry value count overflow';return
    endif
    call build_checked_mpi_displacements(value_send_counts,value_send_displs,cursor,counts_ok)
    if(counts_ok)call build_checked_mpi_displacements(value_recv_counts,value_recv_displs,cursor,counts_ok)
    local_bad=merge(0,1,counts_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='distributed symmetry value displacement overflow';return
    endif
    call MPI_Alltoallv(response_send,value_send_counts,value_send_displs,MPI_DOUBLE_COMPLEX,&
      response_recv,value_recv_counts,value_recv_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
    image=(0d0,0d0)
    do index=1,total_send
      image(:,request_destinations(index))=response_recv(:,index)
    enddo
    ok=.true.
#else
    ok=.false.;message='distributed symmetry point exchange requires MPI'
#endif
  end subroutine exchange_dg_point_permuted_orbital_rows

  subroutine accumulate_dg_lcfo_buffer_contributions_to_core(comm,buffer_ids,buffer_contributions,&
      core_ids,core_values,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::buffer_ids(:),core_ids(:)
    complex(real64),intent(in)::buffer_contributions(:,:)
    complex(real64),intent(out)::core_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer(int64),allocatable::source_ids(:),sorted_core_ids(:)
    integer,allocatable::sorted_core_positions(:)
    complex(real64),allocatable::source_values(:,:)
    integer::rank,nproc,ierr,source,point,position,nstate,nbox,ncore,i,j,key_position
    integer::local_shape(3),minimum_shape(3),maximum_shape(3)
    integer(int64)::key_id

    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nstate=size(buffer_contributions,1);nbox=size(buffer_ids);ncore=size(core_ids)
    local_shape=[nstate,nbox,ncore]
    call MPI_Allreduce(local_shape,minimum_shape,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_shape,maximum_shape,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.all(minimum_shape==maximum_shape).and.nstate>0.and.nbox>0.and.ncore>0.and.&
      size(buffer_contributions,2)==nbox.and.all(shape(core_values)==[nstate,ncore])
    if(.not.ok)then;message='invalid or rank-inconsistent LCFO contribution shape';return;end if
    allocate(source_ids(nbox),source_values(nstate,nbox),sorted_core_ids(ncore),&
      sorted_core_positions(ncore))
    sorted_core_ids=core_ids;sorted_core_positions=[(i,i=1,ncore)]
    do i=2,ncore
      key_id=sorted_core_ids(i);key_position=sorted_core_positions(i);j=i-1
      do while(j>=1)
        if(sorted_core_ids(j)<=key_id)exit
        sorted_core_ids(j+1)=sorted_core_ids(j);sorted_core_positions(j+1)=sorted_core_positions(j);j=j-1
      end do
      sorted_core_ids(j+1)=key_id;sorted_core_positions(j+1)=key_position
    end do
    if(any(sorted_core_ids(2:ncore)==sorted_core_ids(1:ncore-1)))then
      ok=.false.;message='LCFO core physical IDs are not uniquely owned';return
    end if
    core_values=(0d0,0d0)
    do source=0,nproc-1
      if(rank==source)then;source_ids=buffer_ids;source_values=buffer_contributions;end if
      call MPI_Bcast(source_ids,nbox,MPI_INTEGER8,source,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='LCFO contribution ID stream failed';return;end if
      call MPI_Bcast(source_values,nstate*nbox,MPI_DOUBLE_COMPLEX,source,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='LCFO contribution value stream failed';return;end if
      do point=1,nbox
        position=find_core_position(source_ids(point),sorted_core_ids,sorted_core_positions)
        if(position>0)core_values(:,position)=core_values(:,position)+source_values(:,point)
      end do
    end do
    ok=all(ieee_is_finite(real(core_values))).and.all(ieee_is_finite(aimag(core_values)))
    if(ok)then;message='';else;message='LCFO accumulated core values are not finite';end if
#else
    ok=.false.;message='LCFO distributed contribution accumulation requires MPI'
#endif
  contains
    integer function find_core_position(id,sorted_ids,sorted_positions) result(position)
      integer(int64),intent(in)::id,sorted_ids(:)
      integer,intent(in)::sorted_positions(:)
      integer::left,right,middle
      position=0;left=1;right=size(sorted_ids)
      do while(left<=right)
        middle=(left+right)/2
        if(sorted_ids(middle)==id)then;position=sorted_positions(middle);return;end if
        if(sorted_ids(middle)<id)then;left=middle+1;else;right=middle-1;end if
      end do
    end function find_core_position
  end subroutine accumulate_dg_lcfo_buffer_contributions_to_core

  subroutine find_dg_group_identity(product_table,identity_operation,ok,message)
    integer,intent(in)::product_table(:,:)
    integer,intent(out)::identity_operation
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::candidate,operation,n,match_count
    logical::is_identity
    ok=.false.;message='';identity_operation=0;n=size(product_table,1)
    if(n<1.or.size(product_table,2)/=n.or.any(product_table<1).or.any(product_table>n))then
      message='invalid group product table';return
    end if
    match_count=0
    do candidate=1,n
      is_identity=.true.
      do operation=1,n
        if(product_table(candidate,operation)/=operation.or.&
            product_table(operation,candidate)/=operation)is_identity=.false.
      end do
      if(is_identity)then
        match_count=match_count+1;identity_operation=candidate
      end if
    end do
    if(match_count/=1)then;message='group product table has no unique identity';return;end if
    ok=.true.
  end subroutine find_dg_group_identity

  subroutine orthonormalize_dg_distributed_seed_space(comm,seed_values,weights,tolerance,basis,&
      retained_rank,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::seed_values(:,:)
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),allocatable,intent(out)::basis(:,:)
    integer,intent(out)::retained_rank
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer(int64),allocatable::identity_map(:,:)
    integer::rank,nlocal,nseed,ierr,i
    integer::identity_product(1,1)
    call MPI_Comm_rank(comm,rank,ierr)
    nlocal=size(seed_values,2);nseed=size(seed_values,1)
    allocate(identity_map(nlocal,1))
    identity_map(:,1)=int(rank,int64)*int(nlocal,int64)+[(int(i,int64),i=1,nlocal)]
    identity_product=1
    call build_dg_distributed_symmetry_closed_basis(comm,seed_values,weights,identity_map,&
      identity_product,nseed,nseed,tolerance,basis,retained_rank,ok,message,&
      minimum_rank=nseed)
#else
    retained_rank=0;ok=.false.;message='distributed seed orthonormalization requires MPI'
#endif
  end subroutine orthonormalize_dg_distributed_seed_space

  subroutine build_dg_distributed_symmetry_closed_basis(comm,seed_values,weights,&
      symmetry_target_box_ids,product_table,required_seed_count,target_rank,tolerance,basis,retained_rank,&
      ok,message,minimum_rank,required_retained_rank)
    integer,intent(in)::comm,required_seed_count,target_rank
    complex(real64),intent(in)::seed_values(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    integer,intent(in)::product_table(:,:)
    complex(real64),allocatable,intent(out)::basis(:,:)
    integer,intent(out)::retained_rank
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::minimum_rank
    integer,intent(out),optional::required_retained_rank
#ifdef USE_MPI
    complex(real64),allocatable::owner_seed(:),image(:),local_overlap(:),global_overlap(:)
    integer(int64),allocatable::all_maps(:,:,:)
    logical,allocatable::seen(:)
    real(real64)::local_norm,global_norm
    integer::rank,nproc,ierr,nseed,nlocal,noperation,seed,operation,owner,point,&
      target_owner,target_point,pass,iw,rank_before,local_bad,global_bad,&
      source_owner,source_point,effective_minimum_rank
    integer(int64)::final_target
    logical::orbit_exceeds

    ok=.false.;message='';retained_rank=0;if(present(required_retained_rank))required_retained_rank=0
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nseed=size(seed_values,1);nlocal=size(seed_values,2);noperation=size(symmetry_target_box_ids,2)
    effective_minimum_rank=target_rank;if(present(minimum_rank))effective_minimum_rank=minimum_rank
    local_bad=merge(0,1,nseed>0.and.nlocal>0.and.noperation>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.required_seed_count>=0.and.&
      required_seed_count<=nseed.and.target_rank>0.and.effective_minimum_rank>0.and.&
      effective_minimum_rank<=target_rank.and.tolerance>0d0.and.&
      all(shape(product_table)==[noperation,noperation]).and.&
      all(ieee_is_finite(weights)).and.all(weights>=0d0).and.&
      all(ieee_is_finite(real(seed_values))).and.all(ieee_is_finite(aimag(seed_values))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed symmetry-closed basis contract';return;end if
    do operation=1,noperation;do point=1,nlocal
      target_owner=int((symmetry_target_box_ids(point,operation)-1_int64)/int(nlocal,int64))
      target_point=int(modulo(symmetry_target_box_ids(point,operation)-1_int64,int(nlocal,int64)))+1
      if(target_owner<0.or.target_owner>=nproc.or.target_point<1.or.target_point>nlocal)local_bad=1
    end do;end do
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='distributed symmetry-closed basis point map is invalid';return;end if
    allocate(all_maps(nlocal,noperation,nproc),seen(nlocal*nproc))
    call MPI_Allgather(symmetry_target_box_ids,nlocal*noperation,MPI_INTEGER8,all_maps,&
      nlocal*noperation,MPI_INTEGER8,comm,ierr)
    local_bad=0
    do operation=1,noperation
      seen=.false.
      do source_owner=1,nproc;do source_point=1,nlocal
        final_target=all_maps(source_point,operation,source_owner)
        if(final_target<1_int64.or.final_target>int(nlocal*nproc,int64))then
          local_bad=1
        elseif(seen(int(final_target)))then
          local_bad=1
        else
          seen(int(final_target))=.true.
        end if
      end do;end do
      if(.not.all(seen))local_bad=1
    end do
    if(any(product_table<1).or.any(product_table>noperation))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='point maps are not a closed permutation group';return;end if
    allocate(basis(target_rank,nlocal),owner_seed(nlocal),image(nlocal),&
      local_overlap(target_rank),global_overlap(target_rank));basis=(0d0,0d0)
    do seed=1,nseed
      rank_before=retained_rank;orbit_exceeds=.false.
      do operation=1,noperation
        image=(0d0,0d0)
        do owner=0,nproc-1
          if(rank==owner)owner_seed=seed_values(seed,:)
          call MPI_Bcast(owner_seed,nlocal,MPI_DOUBLE_COMPLEX,owner,comm,ierr)
          do point=1,nlocal
            target_owner=int((symmetry_target_box_ids(point,operation)-1_int64)/int(nlocal,int64))
            if(target_owner/=owner)cycle
            target_point=int(modulo(symmetry_target_box_ids(point,operation)-1_int64,&
              int(nlocal,int64)))+1
            image(point)=owner_seed(target_point)
          end do
        end do
        do pass=1,2
          local_overlap=(0d0,0d0)
          do iw=1,retained_rank
            local_overlap(iw)=sum(weights*conjg(basis(iw,:))*image)
          end do
          call MPI_Allreduce(local_overlap,global_overlap,target_rank,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,comm,ierr)
          do iw=1,retained_rank
            image=image-global_overlap(iw)*basis(iw,:)
          end do
        end do
        local_norm=sum(weights*abs(image)**2)
        call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        if(global_norm<=tolerance**2)cycle
        if(retained_rank==target_rank)then;orbit_exceeds=.true.;exit;end if
        retained_rank=retained_rank+1;basis(retained_rank,:)=image/sqrt(global_norm)
      end do
      if(orbit_exceeds)then
        basis(rank_before+1:retained_rank,:)=(0d0,0d0);retained_rank=rank_before
        if(seed<=required_seed_count)then
          message='required symmetry orbit exceeds target rank';return
        end if
      end if
      if(seed==required_seed_count.and.present(required_retained_rank))required_retained_rank=retained_rank
      if(seed>=required_seed_count.and.retained_rank>=effective_minimum_rank)exit
    end do
    if(required_seed_count==0.and.present(required_retained_rank))required_retained_rank=0
    if(required_seed_count>0.and.present(required_retained_rank))then
      if(required_retained_rank<1)then;message='required seed closure was not completed';return;end if
    end if
    if(retained_rank<effective_minimum_rank)then
      message='symmetry-closed seeds do not fill minimum rank';return
    end if
    ok=.true.
#else
    retained_rank=0;ok=.false.;message='distributed symmetry-closed basis requires MPI'
#endif
  end subroutine

  subroutine build_dg_group_averaged_occupied_candidates_reference(comm,occupied,weights,&
      symmetry_target_box_ids,product_table,identity_operation,tolerance,candidates,spectrum,&
      candidate_rank,projector_trace,closure_residual,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,product_table(:,:),identity_operation
    complex(real64),intent(in)::occupied(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    complex(real64),allocatable,intent(out)::candidates(:,:)
    real(real64),allocatable,intent(out)::spectrum(:)
    integer,intent(out)::candidate_rank
    real(real64),intent(out)::projector_trace,closure_residual
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::left_image(:,:),right_image(:,:),local_block(:,:),global_block(:,:),&
      local_occupied_metric(:,:),occupied_metric(:,:),&
      orbit_gram(:,:),orbit_vectors(:,:),label(:,:),composed_label(:,:),expected_label(:,:)
    real(real64),allocatable::all_spectrum(:),total_residual(:),boundary_residual(:),&
      interior_residual(:)
    logical,allocatable::no_boundary(:)
    integer::noccupied,nlocal,noperation,orbit_rank,left_operation,right_operation,&
      left_first,right_first,i,j,k,ierr,rank
    real(real64)::local_group_defect,global_group_defect,occupied_metric_defect,occupied_metric_scale
    integer(int64)::complex_bytes,real_bytes
    logical::eigen_ok
    character(256)::detail

    ok=.false.;message='';candidate_rank=0;projector_trace=huge(1d0)
    closure_residual=huge(1d0);workspace_peak_bytes=0_int64
    noccupied=size(occupied,1);nlocal=size(occupied,2)
    noperation=size(symmetry_target_box_ids,2)
    if(noccupied<1.or.nlocal<1.or.noperation<1.or.size(weights)/=nlocal.or.&
        size(symmetry_target_box_ids,1)/=nlocal.or.&
        any(shape(product_table)/=[noperation,noperation]).or.&
        identity_operation<1.or.identity_operation>noperation.or.tolerance<=0d0.or.&
        .not.all(ieee_is_finite(weights)).or.any(weights<0d0).or.&
        .not.all(ieee_is_finite(real(occupied))).or.&
        .not.all(ieee_is_finite(aimag(occupied))))then
      message='invalid group-averaged occupied-projector contract';return
    endif
    if(noccupied>huge(orbit_rank)/noperation)then
      message='group-averaged occupied orbit rank overflow';return
    endif
    if(any(product_table<1).or.any(product_table>noperation))then
      message='group-averaged occupied product table is invalid';return
    endif
    allocate(local_occupied_metric(noccupied,noccupied),occupied_metric(noccupied,noccupied))
    do j=1,noccupied;do i=1,noccupied
      local_occupied_metric(i,j)=sum(weights*conjg(occupied(i,:))*occupied(j,:))
    enddo;enddo
    call MPI_Allreduce(local_occupied_metric,occupied_metric,noccupied*noccupied,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='group-averaged occupied metric reduction failed';return
    endif
    occupied_metric_scale=max(1d0,maxval(abs(occupied_metric)))
    do i=1,noccupied;occupied_metric(i,i)=occupied_metric(i,i)-1d0;enddo
    occupied_metric_defect=maxval(abs(occupied_metric))/occupied_metric_scale
    if(occupied_metric_defect>tolerance)then
      message='group-averaged occupied input is not metric orthonormal';return
    endif
    deallocate(local_occupied_metric,occupied_metric)
    call MPI_Comm_rank(comm,rank,ierr)
    allocate(label(1,nlocal),composed_label(1,nlocal),expected_label(1,nlocal))
    do i=1,nlocal;label(1,i)=cmplx(real(rank*nlocal+i,real64),0d0,real64);enddo
    do left_operation=1,noperation;do right_operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,label,&
        symmetry_target_box_ids(:,right_operation),expected_label,ok,detail)
      if(.not.ok)then;message='group-averaged right point action: '//trim(detail);return;endif
      call exchange_dg_point_permuted_orbital_rows(comm,expected_label,&
        symmetry_target_box_ids(:,left_operation),composed_label,ok,detail)
      if(.not.ok)then;message='group-averaged composed point action: '//trim(detail);return;endif
      call exchange_dg_point_permuted_orbital_rows(comm,label,&
        symmetry_target_box_ids(:,product_table(right_operation,left_operation)),expected_label,ok,detail)
      if(.not.ok)then;message='group-averaged product point action: '//trim(detail);return;endif
      local_group_defect=maxval(abs(composed_label-expected_label))
      call MPI_Allreduce(local_group_defect,global_group_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_group_defect>0d0)then
        ok=.false.;message='group-averaged point actions do not realize the product table';return
      endif
    enddo;enddo
    call exchange_dg_point_permuted_orbital_rows(comm,label,&
      symmetry_target_box_ids(:,identity_operation),expected_label,ok,detail)
    if(.not.ok.or.maxval(abs(expected_label-label))>0d0)then
      ok=.false.;message='group-averaged identity point action is invalid';return
    endif
    deallocate(label,composed_label,expected_label)
    orbit_rank=noccupied*noperation
    allocate(orbit_gram(orbit_rank,orbit_rank),left_image(noccupied,nlocal),&
      right_image(noccupied,nlocal),local_block(noccupied,noccupied),&
      global_block(noccupied,noccupied));orbit_gram=(0d0,0d0)
    do left_operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
        symmetry_target_box_ids(:,left_operation),left_image,ok,detail)
      if(.not.ok)then;message='group-averaged occupied left image: '//trim(detail);return;endif
      left_first=(left_operation-1)*noccupied+1
      do right_operation=1,noperation
        call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
          symmetry_target_box_ids(:,right_operation),right_image,ok,detail)
        if(.not.ok)then;message='group-averaged occupied right image: '//trim(detail);return;endif
        do j=1,noccupied;do i=1,noccupied
          local_block(i,j)=sum(weights*conjg(left_image(i,:))*right_image(j,:))/&
            real(noperation,real64)
        enddo;enddo
        call MPI_Allreduce(local_block,global_block,noccupied*noccupied,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;ok=.false.;message='group-averaged orbit-Gram reduction failed';return;endif
        right_first=(right_operation-1)*noccupied+1
        orbit_gram(left_first:left_first+noccupied-1,right_first:right_first+noccupied-1)=global_block
      enddo
    enddo
    orbit_gram=0.5d0*(orbit_gram+conjg(transpose(orbit_gram)))
    call hermitian_eigensystem(orbit_gram,all_spectrum,orbit_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;ok=.false.;message='group-averaged occupied eigensystem: '//trim(detail);return;endif
    projector_trace=sum(all_spectrum)
    candidate_rank=count(all_spectrum>tolerance*max(1d0,maxval(abs(all_spectrum))))
    if(candidate_rank<1)then;ok=.false.;message='group-averaged occupied projector has zero rank';return;endif
    allocate(candidates(candidate_rank,nlocal),spectrum(candidate_rank));candidates=(0d0,0d0)
    do i=1,candidate_rank
      j=orbit_rank-i+1
      spectrum(i)=all_spectrum(j)
    enddo
    do left_operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
        symmetry_target_box_ids(:,left_operation),left_image,ok,detail)
      if(.not.ok)then;message='group-averaged occupied reconstruction image: '//trim(detail);return;endif
      left_first=(left_operation-1)*noccupied+1
      do i=1,candidate_rank
        j=orbit_rank-i+1
        do k=1,noccupied
          candidates(i,:)=candidates(i,:)+orbit_vectors(left_first+k-1,j)*left_image(k,:)/&
            sqrt(real(noperation,real64)*spectrum(i))
        enddo
      enddo
    enddo
    allocate(total_residual(noperation),boundary_residual(noperation),&
      interior_residual(noperation),no_boundary(nlocal));no_boundary=.false.
    call measure_dg_rank_fixed_symmetry_residuals(comm,candidates,weights,&
      symmetry_target_box_ids,no_boundary,total_residual=total_residual,&
      boundary_residual=boundary_residual,interior_residual=interior_residual,&
      ok=ok,message=detail)
    if(.not.ok)then;message='group-averaged occupied closure: '//trim(detail);return;endif
    closure_residual=maxval(total_residual)
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    real_bytes=int(storage_size(0d0)/8,int64)
    workspace_peak_bytes=complex_bytes*int(size(left_image)+size(right_image)+size(local_block)+&
      size(global_block)+size(orbit_gram)+size(orbit_vectors)+size(candidates),int64)+&
      real_bytes*int(size(all_spectrum)+size(spectrum)+size(total_residual)+&
      size(boundary_residual)+size(interior_residual),int64)
    ok=ieee_is_finite(projector_trace).and.ieee_is_finite(closure_residual).and.&
      workspace_peak_bytes>0_int64
    if(ok)then;message='';else;message='nonfinite group-averaged occupied-projector receipt';endif
#else
    candidate_rank=0;projector_trace=huge(1d0);closure_residual=huge(1d0)
    workspace_peak_bytes=0_int64;ok=.false.
    message='group-averaged occupied projector requires MPI'
#endif
  end subroutine

#if defined(USE_MPI) && defined(USE_EIGENEXA)
  subroutine build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      symmetry_target_box_ids,product_table,identity_operation,requested_count,tolerance,&
      candidates,spectrum,candidate_rank,projector_trace,closure_residual,gamma_real_defect,&
      workspace_peak_bytes,ok,message,selected_edge,rejected_edge,cluster_gap,&
      cocycle_translation_target_box_ids,translation_cocycle)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::comm,product_table(:,:),identity_operation,requested_count
    complex(real64),intent(in)::occupied(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    complex(real64),allocatable,intent(out)::candidates(:,:)
    real(real64),allocatable,intent(out)::spectrum(:)
    integer,intent(out)::candidate_rank
    real(real64),intent(out)::projector_trace,closure_residual,gamma_real_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),intent(out),optional::selected_edge,rejected_edge,cluster_gap
    integer(int64),intent(in),optional::cocycle_translation_target_box_ids(:,:)
    integer,intent(in),optional::translation_cocycle(:,:)
    complex(real64),allocatable::left_image(:,:),right_image(:,:),local_block(:,:),global_block(:,:),&
      local_occupied_metric(:,:),occupied_metric(:,:),&
      label(:,:),composed_label(:,:),expected_label(:,:),cocycle_label(:,:)
    real(real64),allocatable::cyclic_gram(:,:),cyclic_vectors(:,:),all_spectrum(:),&
      eigenvector(:),total_residual(:),boundary_residual(:),interior_residual(:)
    logical,allocatable::no_boundary(:)
    integer::noccupied,nlocal,noperation,orbit_rank,left_operation,right_operation,&
      left_first,right_first,i,j,k,global_row,global_column,local_row,local_column,ierr,rank
    real(real64)::local_imaginary,global_imaginary,scale,local_group_defect,global_group_defect,&
      occupied_metric_defect,occupied_metric_scale
    integer(int64)::complex_bytes,real_bytes
    logical::eigen_ok
    character(256)::detail

    ok=.false.;message='';candidate_rank=0;projector_trace=huge(1d0)
    closure_residual=huge(1d0);gamma_real_defect=huge(1d0);workspace_peak_bytes=0_int64
    if(present(selected_edge))selected_edge=huge(1d0)
    if(present(rejected_edge))rejected_edge=huge(1d0)
    if(present(cluster_gap))cluster_gap=huge(1d0)
    noccupied=size(occupied,1);nlocal=size(occupied,2);noperation=size(symmetry_target_box_ids,2)
    if(present(cocycle_translation_target_box_ids).neqv.present(translation_cocycle))then
      message='distributed group-average cocycle inputs must be requested together';return
    endif
    if(.not.info%flag_eigenexa_init.or.noccupied<1.or.nlocal<1.or.noperation<1.or.&
        size(weights)/=nlocal.or.size(symmetry_target_box_ids,1)/=nlocal.or.&
        any(shape(product_table)/=[noperation,noperation]).or.any(product_table<1).or.&
        any(product_table>noperation).or.identity_operation<1.or.identity_operation>noperation.or.&
        requested_count<1.or.tolerance<=0d0.or..not.all(ieee_is_finite(weights)).or.&
        any(weights<0d0).or..not.all(ieee_is_finite(real(occupied))).or.&
        .not.all(ieee_is_finite(aimag(occupied))))then
      message='invalid distributed group-averaged occupied-projector contract';return
    endif
    if(present(translation_cocycle))then
      if(size(cocycle_translation_target_box_ids,1)/=nlocal.or.&
          size(cocycle_translation_target_box_ids,2)<1.or.&
          any(shape(translation_cocycle)/=[noperation,noperation]).or.&
          any(translation_cocycle<1).or.&
          any(translation_cocycle>size(cocycle_translation_target_box_ids,2)))then
        message='invalid distributed group-average translation cocycle';return
      endif
    endif
    if(noccupied>huge(orbit_rank)/noperation)then
      message='distributed group-averaged occupied orbit rank overflow';return
    endif
    orbit_rank=noccupied*noperation
    if(requested_count>orbit_rank.or.info%nrow_local<1.or.info%ncol_local<1)then
      message='invalid distributed group-averaged requested rank or EigenExa layout';return
    endif
    allocate(local_occupied_metric(noccupied,noccupied),occupied_metric(noccupied,noccupied))
    do j=1,noccupied;do i=1,noccupied
      local_occupied_metric(i,j)=sum(weights*conjg(occupied(i,:))*occupied(j,:))
    enddo;enddo
    call MPI_Allreduce(local_occupied_metric,occupied_metric,noccupied*noccupied,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='distributed group-average occupied metric reduction failed';return
    endif
    occupied_metric_scale=max(1d0,maxval(abs(occupied_metric)))
    do i=1,noccupied;occupied_metric(i,i)=occupied_metric(i,i)-1d0;enddo
    occupied_metric_defect=maxval(abs(occupied_metric))/occupied_metric_scale
    if(occupied_metric_defect>tolerance)then
      message='distributed group-average occupied input is not metric orthonormal';return
    endif
    deallocate(local_occupied_metric,occupied_metric)
    call MPI_Comm_rank(comm,rank,ierr)
    allocate(label(1,nlocal),composed_label(1,nlocal),expected_label(1,nlocal),cocycle_label(1,nlocal))
    do i=1,nlocal;label(1,i)=cmplx(real(rank*nlocal+i,real64),0d0,real64);enddo
    do left_operation=1,noperation;do right_operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,label,&
        symmetry_target_box_ids(:,right_operation),expected_label,ok,detail)
      if(.not.ok)then;message='distributed group-average right point action: '//trim(detail);return;endif
      call exchange_dg_point_permuted_orbital_rows(comm,expected_label,&
        symmetry_target_box_ids(:,left_operation),composed_label,ok,detail)
      if(.not.ok)then;message='distributed group-average composed point action: '//trim(detail);return;endif
      if(present(translation_cocycle))then
        call exchange_dg_point_permuted_orbital_rows(comm,label,&
          cocycle_translation_target_box_ids(:,translation_cocycle(right_operation,left_operation)),&
          cocycle_label,ok,detail)
        if(.not.ok)then;message='distributed group-average cocycle translation: '//trim(detail);return;endif
        call exchange_dg_point_permuted_orbital_rows(comm,cocycle_label,&
          symmetry_target_box_ids(:,product_table(right_operation,left_operation)),expected_label,ok,detail)
        if(.not.ok)then;message='distributed group-average cocycle representative: '//trim(detail);return;endif
      else
        call exchange_dg_point_permuted_orbital_rows(comm,label,&
          symmetry_target_box_ids(:,product_table(right_operation,left_operation)),expected_label,ok,detail)
        if(.not.ok)then;message='distributed group-average product point action: '//trim(detail);return;endif
      endif
      local_group_defect=maxval(abs(composed_label-expected_label))
      call MPI_Allreduce(local_group_defect,global_group_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_group_defect>0d0)then
        ok=.false.;message='distributed group-average point actions do not realize product table';return
      endif
    enddo;enddo
    call exchange_dg_point_permuted_orbital_rows(comm,label,&
      symmetry_target_box_ids(:,identity_operation),expected_label,ok,detail)
    if(.not.ok.or.maxval(abs(expected_label-label))>0d0)then
      ok=.false.;message='distributed group-average identity point action is invalid';return
    endif
    deallocate(label,composed_label,expected_label,cocycle_label)
    allocate(left_image(noccupied,nlocal),right_image(noccupied,nlocal),&
      local_block(noccupied,noccupied),global_block(noccupied,noccupied),&
      cyclic_gram(info%nrow_local,info%ncol_local),&
      cyclic_vectors(info%nrow_local,info%ncol_local),all_spectrum(orbit_rank))
    cyclic_gram=0d0;local_imaginary=0d0;scale=0d0
    do left_operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
        symmetry_target_box_ids(:,left_operation),left_image,ok,detail)
      if(.not.ok)then;message='distributed group-average left image: '//trim(detail);return;endif
      left_first=(left_operation-1)*noccupied+1
      do right_operation=left_operation,noperation
        call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
          symmetry_target_box_ids(:,right_operation),right_image,ok,detail)
        if(.not.ok)then;message='distributed group-average right image: '//trim(detail);return;endif
        do j=1,noccupied;do i=1,noccupied
          local_block(i,j)=sum(weights*conjg(left_image(i,:))*right_image(j,:))/&
            real(noperation,real64)
        enddo;enddo
        call MPI_Allreduce(local_block,global_block,noccupied*noccupied,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;ok=.false.;message='distributed group-average block reduction failed';return;endif
        local_imaginary=max(local_imaginary,maxval(abs(aimag(global_block))))
        scale=max(scale,maxval(abs(global_block)))
        right_first=(right_operation-1)*noccupied+1
        do j=1,noccupied
          global_column=right_first+j-1
          if(eigen_owner_node(global_column,info%npcol,info%mycol)/=info%mycol)cycle
          local_column=eigen_translate_g2l(global_column,info%npcol,info%mycol)
          do i=1,noccupied
            global_row=left_first+i-1
            if(eigen_owner_node(global_row,info%nprow,info%myrow)/=info%myrow)cycle
            local_row=eigen_translate_g2l(global_row,info%nprow,info%myrow)
            if(left_operation==right_operation)then
              cyclic_gram(local_row,local_column)=0.5d0*&
                (real(global_block(i,j),real64)+real(global_block(j,i),real64))
            else
              cyclic_gram(local_row,local_column)=real(global_block(i,j),real64)
            endif
          enddo
        enddo
        if(left_operation/=right_operation)then
          do j=1,noccupied
            global_column=left_first+j-1
            if(eigen_owner_node(global_column,info%npcol,info%mycol)/=info%mycol)cycle
            local_column=eigen_translate_g2l(global_column,info%npcol,info%mycol)
            do i=1,noccupied
              global_row=right_first+i-1
              if(eigen_owner_node(global_row,info%nprow,info%myrow)/=info%myrow)cycle
              local_row=eigen_translate_g2l(global_row,info%nprow,info%myrow)
              cyclic_gram(local_row,local_column)=real(global_block(j,i),real64)
            enddo
          enddo
        endif
      enddo
    enddo
    call MPI_Allreduce(local_imaginary,global_imaginary,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    gamma_real_defect=global_imaginary/max(1d0,scale)
    if(ierr/=MPI_SUCCESS.or.gamma_real_defect>tolerance)then
      ok=.false.;message='distributed group-average is not Gamma real';return
    endif
    call eigen_pdsyevd_ex_distributed_blocks(info,orbit_rank,cyclic_gram,all_spectrum,&
      cyclic_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;ok=.false.;message='distributed group-average eigensystem: '//trim(detail);return;endif
    if(present(selected_edge))selected_edge=all_spectrum(orbit_rank-requested_count+1)
    if(requested_count<orbit_rank)then
      if(present(rejected_edge))rejected_edge=all_spectrum(orbit_rank-requested_count)
      if(present(cluster_gap))cluster_gap=&
        all_spectrum(orbit_rank-requested_count+1)-all_spectrum(orbit_rank-requested_count)
    else
      if(present(rejected_edge))rejected_edge=0d0
      if(present(cluster_gap))cluster_gap=all_spectrum(1)
    endif
    if(requested_count<orbit_rank)then
      if(abs(all_spectrum(orbit_rank-requested_count+1)-all_spectrum(orbit_rank-requested_count))<=&
          tolerance*max(1d0,maxval(abs(all_spectrum))))then
        ok=.false.;message='requested occupied rank cuts a group-averaged degenerate block';return
      endif
    endif
    projector_trace=sum(all_spectrum);candidate_rank=requested_count
    allocate(candidates(candidate_rank,nlocal),spectrum(candidate_rank),eigenvector(orbit_rank))
    candidates=(0d0,0d0)
    do i=1,candidate_rank
      j=orbit_rank-i+1;spectrum(i)=all_spectrum(j)
      if(spectrum(i)<=tolerance*max(1d0,maxval(abs(all_spectrum))))then
        ok=.false.;message='distributed group-average requested candidate has zero weight';return
      endif
      call gather_eigenvector(j,eigenvector,ok,detail)
      if(.not.ok)then;message=trim(detail);return;endif
      do left_operation=1,noperation
        call exchange_dg_point_permuted_orbital_rows(comm,occupied,&
          symmetry_target_box_ids(:,left_operation),left_image,ok,detail)
        if(.not.ok)then;message='distributed group-average reconstruction: '//trim(detail);return;endif
        left_first=(left_operation-1)*noccupied+1
        do k=1,noccupied
          candidates(i,:)=candidates(i,:)+eigenvector(left_first+k-1)*left_image(k,:)/&
            sqrt(real(noperation,real64)*spectrum(i))
        enddo
      enddo
    enddo
    allocate(total_residual(noperation),boundary_residual(noperation),interior_residual(noperation),&
      no_boundary(nlocal));no_boundary=.false.
    call measure_dg_rank_fixed_symmetry_residuals(comm,candidates,weights,symmetry_target_box_ids,&
      no_boundary,total_residual=total_residual,boundary_residual=boundary_residual,&
      interior_residual=interior_residual,ok=ok,message=detail)
    if(.not.ok)then;message='distributed group-average closure: '//trim(detail);return;endif
    closure_residual=maxval(total_residual)
    complex_bytes=int(storage_size((0d0,0d0))/8,int64);real_bytes=int(storage_size(0d0)/8,int64)
    workspace_peak_bytes=complex_bytes*int(size(left_image)+size(right_image)+size(local_block)+&
      size(global_block)+size(candidates),int64)+real_bytes*int(size(cyclic_gram)+&
      size(cyclic_vectors)+size(all_spectrum)+size(eigenvector)+size(spectrum)+&
      size(total_residual)+size(boundary_residual)+size(interior_residual),int64)
    ok=ieee_is_finite(projector_trace).and.ieee_is_finite(closure_residual).and.&
      workspace_peak_bytes>0_int64
    if(ok)then;message='';else;message='nonfinite distributed group-average receipt';endif
  contains
    subroutine gather_eigenvector(column,values,gather_ok,gather_message)
      integer,intent(in)::column
      real(real64),intent(out)::values(:)
      logical,intent(out)::gather_ok
      character(*),intent(out)::gather_message
      integer::lr,lc,gr,gc,error,row_start,row_end,column_start,column_end
      values=0d0;row_start=eigen_loop_start(1,info%nprow,info%myrow)
      row_end=eigen_loop_end(orbit_rank,info%nprow,info%myrow)
      column_start=eigen_loop_start(1,info%npcol,info%mycol)
      column_end=eigen_loop_end(orbit_rank,info%npcol,info%mycol)
      do lc=column_start,column_end
        gc=eigen_translate_l2g(lc,info%npcol,info%mycol);if(gc/=column)cycle
        do lr=row_start,row_end
          gr=eigen_translate_l2g(lr,info%nprow,info%myrow);values(gr)=cyclic_vectors(lr,lc)
        enddo
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,values,orbit_rank,MPI_DOUBLE_PRECISION,MPI_SUM,comm,error)
      gather_ok=error==MPI_SUCCESS
      if(gather_ok)then;gather_message='';else;gather_message='distributed group-average vector gather failed';endif
    end subroutine
  end subroutine

#endif

  subroutine inverse_dg_translation_character_orbits(comm,row_ids,global_row_count,characters,product_table,&
      identity_operation,catalog_fingerprint,sector_values,sector_gradients,tolerance,orbit_values,orbit_gradients,density_defect,&
      orthogonality_defect,gamma_real_defect,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count,product_table(:,:),identity_operation
    integer(int64),intent(in)::row_ids(:),catalog_fingerprint
    complex(real64),intent(in)::characters(:,:),sector_values(:,:,:),sector_gradients(:,:,:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::orbit_values(:,:,:),orbit_gradients(:,:,:,:)
    real(real64),intent(out)::density_defect,orthogonality_defect,gamma_real_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,ninternal,ntranslation,i,j,a,b,t,c,ierr,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer
    integer,allocatable::ownership_count(:)
    real(real64)::normalization,local_value,global_value,sector_density,orbit_density,scale
    complex(real64)::overlap,local_overlap,expected
    integer(int64)::local_fingerprint,metadata_fingerprint,minimum_fingerprint,maximum_fingerprint,&
      output_elements,gradient_elements,total_elements,byte_count,quantized,element_hash,raw_bits
    complex(real64)::bilinear
    real(real64)::minimum_tolerance,maximum_tolerance,safe_input_magnitude
    ok=.false.;message='';density_defect=huge(1d0);orthogonality_defect=huge(1d0)
    gamma_real_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    nlocal=size(row_ids);ninternal=size(sector_values,2);ntranslation=size(characters,1)
    local_bad=0
    if(global_row_count<1.or.ntranslation<1.or.size(characters,2)/=ntranslation.or.ninternal<1.or.&
        size(sector_values,1)/=nlocal.or.size(sector_values,3)/=ntranslation)then
      local_bad=1
    elseif(any(shape(product_table)/=[ntranslation,ntranslation]).or.&
        identity_operation<1.or.identity_operation>ntranslation)then
      local_bad=1
    elseif(any(shape(sector_gradients)/=[3,nlocal,ninternal,ntranslation]).or.&
        tolerance<16d0*acos(-1d0)/real(huge(0_int64),real64).or.tolerance>1d-2.or.&
        .not.ieee_is_finite(tolerance))then
      local_bad=1
    elseif(.not.all(ieee_is_finite(real(characters))).or.&
        .not.all(ieee_is_finite(aimag(characters))).or.&
        .not.all(ieee_is_finite(real(sector_values))).or.&
        .not.all(ieee_is_finite(aimag(sector_values))).or.&
        .not.all(ieee_is_finite(real(sector_gradients))).or.&
        .not.all(ieee_is_finite(aimag(sector_gradients))).or.any(row_ids<1_int64).or.&
        any(row_ids>int(global_row_count,int64)).or.catalog_fingerprint==0_int64)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid inverse translation-character transform contract';return
    endif
    call agree_integer(ntranslation);if(global_bad/=0)return
    call agree_integer(ninternal);if(global_bad/=0)return
    call agree_integer(global_row_count);if(global_bad/=0)return
    call agree_integer(identity_operation);if(global_bad/=0)return
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='inverse character metadata reduction failed';return;endif
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_tolerance/=maximum_tolerance)then
      message='inverse character metadata disagree across ranks';return
    endif
    ! Parseval preserves the all-orbit norm; the quarter-range margin keeps its
    ! tolerance-quantized density/gradient receipt inside signed int64.
    safe_input_magnitude=sqrt(0.25d0*real(huge(0_int64),real64)*100d0*tolerance/&
      (real(ntranslation,real64)*real(ninternal,real64)))
    local_bad=merge(1,0,maxval(abs(sector_values))>safe_input_magnitude.or.&
      maxval(abs(sector_gradients))>safe_input_magnitude)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='inverse character sector magnitude exceeds safe receipt range';return
    endif
    metadata_fingerprint=catalog_fingerprint
    do i=1,ntranslation;do j=1,ntranslation
      metadata_fingerprint=ieor(ishftc(metadata_fingerprint,7),int(product_table(i,j),int64))
      raw_bits=transfer(real(characters(i,j),real64),raw_bits)
      metadata_fingerprint=ieor(ishftc(metadata_fingerprint,7),raw_bits)
      raw_bits=transfer(aimag(characters(i,j)),raw_bits)
      metadata_fingerprint=ieor(ishftc(metadata_fingerprint,7),raw_bits)
    enddo;enddo
    call MPI_Allreduce(metadata_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='inverse character metadata reduction failed';return;endif
    call MPI_Allreduce(metadata_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='inverse character catalogs disagree across ranks';return
    endif
    allocate(ownership_count(global_row_count),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='inverse character ownership allocation failed';return;endif
    ownership_count=0
    do i=1,nlocal;ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      message='inverse character rows are not uniquely and completely owned';return
    endif
    local_bad=merge(1,0,any(product_table<1).or.any(product_table>ntranslation))
    do c=1,ntranslation
      if(abs(characters(c,identity_operation)-1d0)>10d0*tolerance.or.&
          maxval(abs(abs(characters(c,:))-1d0))>10d0*tolerance)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='inverse character catalog fails collective range or unit validation';return
    endif
    local_bad=0
    if(local_bad==0)then;do i=1,ntranslation
      if(product_table(identity_operation,i)/=i.or.product_table(i,identity_operation)/=i)then
        local_bad=1
      endif
      do j=1,ntranslation
        do c=1,ntranslation
          expected=characters(c,i)*characters(c,j)
          if(abs(characters(c,product_table(i,j))-expected)>10d0*tolerance)then
            local_bad=1
          endif
        enddo
      enddo
    enddo;endif
    do c=1,ntranslation
      do j=1,ntranslation
        overlap=sum(characters(c,:)*conjg(characters(j,:)))
        expected=merge(cmplx(real(ntranslation,real64),0d0,real64),(0d0,0d0),c==j)
        if(abs(overlap-expected)>10d0*tolerance*real(ntranslation,real64))then
          local_bad=1
        endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='inverse character catalog fails collective algebra validation';return
    endif
    allocate(orbit_values(nlocal,ninternal,ntranslation),&
      orbit_gradients(3,nlocal,ninternal,ntranslation),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='inverse character orbit allocation failed';return
    endif
    normalization=1d0/sqrt(real(ntranslation,real64));orbit_values=(0d0,0d0)
    orbit_gradients=(0d0,0d0)
    do t=1,ntranslation;do c=1,ntranslation
      orbit_values(:,:,t)=orbit_values(:,:,t)+normalization*conjg(characters(c,t))*sector_values(:,:,c)
      orbit_gradients(:,:,:,t)=orbit_gradients(:,:,:,t)+&
        normalization*conjg(characters(c,t))*sector_gradients(:,:,:,c)
    enddo;enddo
    local_value=0d0
    do i=1,nlocal
      sector_density=sum(abs(sector_values(i,:,:))**2);orbit_density=sum(abs(orbit_values(i,:,:))**2)
      local_value=max(local_value,abs(orbit_density-sector_density))
    enddo
    call MPI_Allreduce(local_value,density_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='inverse character density reduction failed';return;endif
    orthogonality_defect=0d0
    do t=1,ntranslation;do a=1,ninternal;do j=1,ntranslation;do b=1,ninternal
      local_overlap=sum(conjg(orbit_values(:,a,t))*orbit_values(:,b,j))
      call MPI_Allreduce(local_overlap,overlap,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='inverse character Gram reduction failed';return;endif
      expected=merge((1d0,0d0),(0d0,0d0),t==j.and.a==b)
      orthogonality_defect=max(orthogonality_defect,abs(overlap-expected))
    enddo;enddo;enddo;enddo
    local_value=max(maxval(abs(aimag(orbit_values))),maxval(abs(aimag(orbit_gradients))))
    call MPI_Allreduce(local_value,gamma_real_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='inverse character Gamma-real reduction failed';return;endif
    scale=max(1d0,real(ntranslation,real64)*real(ninternal,real64))
    if(density_defect>10d0*tolerance*scale.or.orthogonality_defect>10d0*tolerance*scale.or.&
        gamma_real_defect>10d0*tolerance*scale)then
      message='inverse character orbit validation failed';return
    endif
    ! Physical invariant fingerprint: density and three gradient norms, not raw coefficients.
    local_fingerprint=0_int64
    do i=1,nlocal
      element_hash=row_ids(i)
      bilinear=sum(orbit_values(i,:,:)*conjg(orbit_values(i,:,:)))
      quantized=nint(real(bilinear,real64)/(100d0*tolerance),int64)
      element_hash=ieor(ishftc(element_hash,9),quantized)
      do j=1,3
        bilinear=sum(orbit_gradients(j,i,:,:)*conjg(orbit_gradients(j,i,:,:)))
        quantized=nint(real(bilinear,real64)/(100d0*tolerance),int64)
        element_hash=ieor(ishftc(element_hash,9),quantized)
      enddo
      local_fingerprint=ieor(local_fingerprint,element_hash)
    enddo
    call MPI_Allreduce(local_fingerprint,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='inverse character fingerprint reduction failed';return;endif
    fingerprint=ieor(fingerprint,metadata_fingerprint)
    output_elements=size(orbit_values,kind=int64);gradient_elements=size(orbit_gradients,kind=int64)
    if(output_elements>huge(0_int64)-gradient_elements)then
      message='inverse character workspace receipt overflows';return
    endif
    total_elements=output_elements+gradient_elements
    if(total_elements>huge(0_int64)/16_int64)then
      message='inverse character workspace receipt overflows';return
    endif
    ! Accounted output allocation; this is not allocator/RSS telemetry.
    byte_count=16_int64*total_elements;workspace_peak_bytes=byte_count
    ok=workspace_peak_bytes>0_int64
    if(ok)then;message='';else;message='inverse character workspace receipt overflows';endif
#else
    ok=.false.;message='inverse character transform requires MPI';density_defect=huge(1d0)
    orthogonality_defect=huge(1d0);gamma_real_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    allocate(orbit_values(0,0,0),orbit_gradients(0,0,0,0))
#endif
  contains
#ifdef USE_MPI
    subroutine agree_integer(value)
      integer,intent(in)::value
      call MPI_Allreduce(value,minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then
        global_bad=1;message='inverse character metadata reduction failed';return
      endif
      call MPI_Allreduce(value,maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      global_bad=merge(1,0,ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)
      if(global_bad/=0)message='inverse character metadata disagree across ranks'
    end subroutine agree_integer
#endif
  end subroutine inverse_dg_translation_character_orbits

#if defined(USE_MPI) && defined(USE_EIGENEXA)
  subroutine split_dg_translation_character_sector_eigenexa(info,comm,row_ids,generator_rows,gamma_rows,&
      characters,generator_orders,element_words,character_conjugates,requested_character,tolerance,&
      catalog_fingerprint,sector_vectors,sector_rank,&
      identity_defect,unitarity_defect,commutator_defect,order_defect,gamma_pairing_defect,&
      fingerprint,workspace_peak_bytes,ok,message)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::comm,generator_orders(:),element_words(:,:),character_conjugates(:),requested_character
    integer(int64),intent(in)::row_ids(:)
    ! Opaque provenance label produced and validated by the canonical Task1 catalog builder.
    integer(int64),intent(in)::catalog_fingerprint
    complex(real64),intent(in)::generator_rows(:,:,:),gamma_rows(:,:),characters(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::sector_vectors(:,:)
    integer,intent(out)::sector_rank
    real(real64),intent(out)::identity_defect,unitarity_defect,commutator_defect,order_defect,&
      gamma_pairing_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_parallel_info)::real_info
    complex(real64),allocatable::discriminator_rows(:,:),adjoint_rows(:,:),product_rows(:,:),&
      left_product_rows(:,:),power_rows(:,:),identity_rows(:,:),trial(:),candidate_rows(:,:),&
      gram(:,:),gamma_projector_rows(:,:),sector_projector_rows(:,:),stream_row(:),remote_sector(:)
    complex(real64),allocatable::sector_multiplicities(:)
    real(real64),allocatable::real_matrix(:,:),real_vectors(:,:),eigenvalues(:),real_column(:)
    integer,allocatable::ownership(:),row_owner(:),row_position(:)
    integer::n,nlocal,ngenerator,ncharacter,multiplicity,i,j,g,h,k,ierr,rank,lr,lc,gr,gc,&
      real_n,accepted,local_bad,global_bad,allocation_status
    real(real64)::scale,norm_value,local_defect
    complex(real64)::overlap,phase
    logical::eigen_ok,receipt_ok
    character(256)::detail

    ok=.false.;message='';sector_rank=0;identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    commutator_defect=huge(1d0);order_defect=huge(1d0);gamma_pairing_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64
    n=size(generator_rows,2);nlocal=size(row_ids);ngenerator=size(generator_rows,3)
    ncharacter=size(characters,1)
    local_bad=merge(0,1,info%flag_eigenexa_init.and.n>=1.and.nlocal>=0.and.&
        (ngenerator>=1.or.ncharacter==1).and.&
        ncharacter>=1.and.size(generator_rows,1)==nlocal.and.all(shape(gamma_rows)==[nlocal,n]).and.&
        size(characters,2)==ngenerator.and.size(generator_orders)==ngenerator.and.&
        all(shape(element_words)==[ncharacter,ngenerator]).and.&
        size(character_conjugates)==ncharacter.and.requested_character>=1.and.&
        requested_character<=ncharacter.and.catalog_fingerprint/=0_int64.and.&
        ieee_is_finite(tolerance).and.tolerance<=1d-2.and.&
        tolerance>=16d0*acos(-1d0)/real(huge(0_int64),real64).and.&
        all(row_ids>=1_int64).and.all(row_ids<=int(n,int64)).and.all(generator_orders>=1).and.&
        all(element_words>=0).and.&
        all(ieee_is_finite(real(generator_rows))).and.&
        all(ieee_is_finite(aimag(generator_rows))).and.&
        all(ieee_is_finite(real(gamma_rows))).and.all(ieee_is_finite(aimag(gamma_rows))).and.&
        all(ieee_is_finite(real(characters))).and.&
        all(ieee_is_finite(aimag(characters))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid translation-character sector contract';return
    endif
    local_bad=merge(0,1,mod(n,ncharacter)==0.and.(ngenerator>=1.or.&
      (ncharacter==1.and.ngenerator==0)).and.all(character_conjugates>=1).and.&
      all(character_conjugates<=ncharacter))
    do g=1,ngenerator
      if(any(element_words(:,g)>=generator_orders(g)))local_bad=1
      if(maxval(abs(abs(characters(:,g))-1d0))>tolerance)local_bad=1
      if(maxval(abs(characters(:,g)**generator_orders(g)-1d0))>10d0*tolerance)local_bad=1
    enddo
    if(local_bad==0)then
      do i=1,ncharacter
        if(maxval(abs(characters(character_conjugates(i),:)-conjg(characters(i,:))))>10d0*tolerance)&
          local_bad=1
        if(character_conjugates(character_conjugates(i))/=i)local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='translation characters violate finite-order conjugate pairing';return
    endif
    if(n>huge(0)/2)then;message='translation-sector realification extent overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr)
    allocate(ownership(n),row_owner(n),row_position(n),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector ownership allocation failed';return;endif
    ownership=0;row_owner=0;row_position=0
    do i=1,nlocal
      ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1
      row_owner(int(row_ids(i)))=rank+1;row_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='translation ownership count reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,row_owner,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='translation ownership rank reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,row_position,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      message='translation generator rows do not uniquely partition the global matrix';return
    endif
    multiplicity=n/ncharacter
    allocate(identity_rows(nlocal,n),power_rows(nlocal,n),product_rows(nlocal,n),&
      left_product_rows(nlocal,n),adjoint_rows(nlocal,n),discriminator_rows(nlocal,n),stream_row(n),&
      stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector row workspace allocation failed';return;endif
    identity_rows=(0d0,0d0)
    do i=1,nlocal;identity_rows(i,int(row_ids(i)))=1d0;enddo
    identity_defect=0d0;unitarity_defect=0d0;commutator_defect=0d0;order_defect=0d0
    if(ncharacter==1.and.ngenerator==0)then
      call distributed_adjoint(gamma_rows,adjoint_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='trivial Gamma adjoint stream failed';return;endif
      call distributed_product(adjoint_rows,gamma_rows,product_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='trivial Gamma unitarity stream failed';return;endif
      local_defect=maxval(abs(product_rows-identity_rows))
      call MPI_Allreduce(local_defect,gamma_pairing_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='trivial Gamma unitarity reduction failed';return;endif
      call distributed_product(gamma_rows,conjg(gamma_rows),product_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='trivial Gamma involution stream failed';return;endif
      local_defect=maxval(abs(product_rows-identity_rows))
      call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='trivial Gamma involution reduction failed';return;endif
      gamma_pairing_defect=max(gamma_pairing_defect,norm_value)
      if(gamma_pairing_defect>10d0*tolerance)then
        message='trivial translation Gamma sewing is not unitary involutory';return
      endif
      allocate(sector_vectors(nlocal,n),stat=allocation_status)
      call allocation_consensus(allocation_status,global_bad)
      if(global_bad/=0)then;message='trivial translation-sector allocation failed';return;endif
      sector_vectors=identity_rows;sector_rank=n
      fingerprint=ieor(ieor(1469598103934665603_int64,catalog_fingerprint),int(n,int64))
      call compute_workspace_receipt(receipt_ok)
      if(.not.receipt_ok)then;message='trivial translation-sector workspace receipt overflows';return;endif
      ok=.true.;message='';return
    endif
    discriminator_rows=(0d0,0d0);scale=max(1d0,maxval(abs(generator_rows)))
    call MPI_Allreduce(MPI_IN_PLACE,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='translation generator scale reduction failed';return;endif
    do g=1,ngenerator
      call distributed_adjoint(generator_rows(:,:,g),adjoint_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation generator adjoint stream failed';return;endif
      call distributed_product(adjoint_rows,generator_rows(:,:,g),product_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation generator unitarity stream failed';return;endif
      local_defect=maxval(abs(product_rows-identity_rows))
      call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation unitarity defect reduction failed';return;endif
      unitarity_defect=max(unitarity_defect,norm_value)
      power_rows=identity_rows
      do k=1,generator_orders(g)
        call distributed_product(power_rows,generator_rows(:,:,g),product_rows,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation generator order stream failed';return;endif
        power_rows=product_rows
      enddo
      local_defect=maxval(abs(power_rows-identity_rows))
      call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation order defect reduction failed';return;endif
      order_defect=max(order_defect,norm_value)
      phase=characters(requested_character,g)
      discriminator_rows=discriminator_rows+2d0*identity_rows-conjg(phase)*generator_rows(:,:,g)-&
        phase*adjoint_rows
      do h=1,g-1
        call distributed_product(generator_rows(:,:,g),generator_rows(:,:,h),product_rows,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation left commutator stream failed';return;endif
        call distributed_product(generator_rows(:,:,h),generator_rows(:,:,g),left_product_rows,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation generator commutator stream failed';return;endif
        local_defect=maxval(abs(product_rows-left_product_rows))
        call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation commutator defect reduction failed';return;endif
        commutator_defect=max(commutator_defect,norm_value)
      enddo
    enddo
    if(max(unitarity_defect,commutator_defect,order_defect)>tolerance*scale)then
      message='translation generators fail unitary commuting finite-order gates';return
    endif
    gamma_pairing_defect=0d0
    call distributed_adjoint(gamma_rows,adjoint_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma sewing adjoint stream failed';return;endif
    call distributed_product(adjoint_rows,gamma_rows,product_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma sewing unitarity stream failed';return;endif
    local_defect=maxval(abs(product_rows-identity_rows))
    call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma sewing unitarity reduction failed';return;endif
    gamma_pairing_defect=max(gamma_pairing_defect,norm_value)
    call distributed_product(gamma_rows,conjg(gamma_rows),product_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma sewing involution stream failed';return;endif
    local_defect=maxval(abs(product_rows-identity_rows))
    call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma sewing involution reduction failed';return;endif
    gamma_pairing_defect=max(gamma_pairing_defect,norm_value)
    do g=1,ngenerator
      call distributed_product(gamma_rows,conjg(generator_rows(:,:,g)),product_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Gamma left covariance stream failed';return;endif
      call distributed_product(product_rows,adjoint_rows,left_product_rows,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Gamma right covariance stream failed';return;endif
      local_defect=maxval(abs(left_product_rows-generator_rows(:,:,g)))
      call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Gamma covariance defect reduction failed';return;endif
      gamma_pairing_defect=max(gamma_pairing_defect,norm_value)
    enddo
    if(gamma_pairing_defect>10d0*tolerance*scale)then
      message='Gamma sewing does not pair conjugate translation characters';return
    endif
    allocate(sector_multiplicities(ncharacter),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation multiplicity allocation failed';return;endif
    sector_multiplicities=(0d0,0d0)
    do h=1,ncharacter
      power_rows=identity_rows
      do g=1,ngenerator
        do k=1,element_words(h,g)
          call distributed_product(power_rows,generator_rows(:,:,g),product_rows,ierr)
          if(ierr/=MPI_SUCCESS)then;message='translation element-word stream failed';return;endif
          power_rows=product_rows
        enddo
      enddo
      overlap=(0d0,0d0)
      do i=1,nlocal;overlap=overlap+power_rows(i,int(row_ids(i)));enddo
      call MPI_Allreduce(MPI_IN_PLACE,overlap,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation element trace reduction failed';return;endif
      if(all(element_words(h,:)==0))then
        local_defect=maxval(abs(power_rows-identity_rows))
        call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation identity defect reduction failed';return;endif
        identity_defect=max(identity_defect,norm_value)
      endif
      do i=1,ncharacter
        phase=(1d0,0d0)
        do g=1,ngenerator;phase=phase*characters(i,g)**element_words(h,g);enddo
        sector_multiplicities(i)=sector_multiplicities(i)+conjg(phase)*overlap/real(ncharacter,real64)
      enddo
    enddo
    if(maxval(abs(real(sector_multiplicities,real64)-real(multiplicity,real64)))>10d0*tolerance*scale.or.&
        maxval(abs(aimag(sector_multiplicities)))>10d0*tolerance*scale)then
      message='translation characters do not have equal complete multiplicity';return
    endif
    real_n=2*n;real_info=info
    call eigen_get_matdims(real_n,real_info%nrow_local,real_info%ncol_local)
    allocate(real_matrix(real_info%nrow_local,real_info%ncol_local),&
      real_vectors(real_info%nrow_local,real_info%ncol_local),eigenvalues(real_n),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector EigenExa allocation failed';return;endif
    real_matrix=0d0
    do i=1,n
      call broadcast_row(discriminator_rows,i,stream_row,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation discriminator row stream failed';return;endif
      do j=1,n
        call put_realified_entry(i,j,real(stream_row(j),real64))
        call put_realified_entry(i,n+j,-aimag(stream_row(j)))
        call put_realified_entry(n+i,j,aimag(stream_row(j)))
        call put_realified_entry(n+i,n+j,real(stream_row(j),real64))
      enddo
    enddo
    call eigen_pdsyevd_ex_distributed_blocks(real_info,real_n,real_matrix,eigenvalues,&
      real_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;message='translation-sector EigenExa solve: '//trim(detail);return;endif
    call validate_dg_translation_sector_cluster(eigenvalues,2*multiplicity,tolerance,eigen_ok,detail)
    if(.not.eigen_ok)then;message=trim(detail);return;endif
    allocate(candidate_rows(nlocal,2*multiplicity),real_column(real_n),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector candidate allocation failed';return;endif
    candidate_rows=(0d0,0d0)
    do j=1,2*multiplicity
      real_column=0d0
      do lc=eigen_loop_start(1,real_info%npcol,real_info%mycol),&
          eigen_loop_end(real_n,real_info%npcol,real_info%mycol)
        gc=eigen_translate_l2g(lc,real_info%npcol,real_info%mycol);if(gc/=j)cycle
        do lr=eigen_loop_start(1,real_info%nprow,real_info%myrow),&
            eigen_loop_end(real_n,real_info%nprow,real_info%myrow)
          gr=eigen_translate_l2g(lr,real_info%nprow,real_info%myrow)
          real_column(gr)=real_vectors(lr,lc)
        enddo
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,real_column,real_n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation-sector eigenvector gather failed';return;endif
      do i=1,nlocal
        candidate_rows(i,j)=cmplx(real_column(int(row_ids(i))),real_column(n+int(row_ids(i))),real64)
      enddo
    enddo
    allocate(sector_vectors(nlocal,multiplicity),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector vector allocation failed';return;endif
    sector_vectors=(0d0,0d0);accepted=0
    do j=1,2*multiplicity
      trial=candidate_rows(:,j)
      do k=1,accepted
        overlap=sum(conjg(sector_vectors(:,k))*trial)
        call MPI_Allreduce(MPI_IN_PLACE,overlap,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='translation-sector orthogonalization reduction failed';return;endif
        trial=trial-overlap*sector_vectors(:,k)
      enddo
      norm_value=sum(abs(trial)**2)
      call MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation-sector norm reduction failed';return;endif
      norm_value=sqrt(norm_value);if(norm_value<=100d0*tolerance)cycle
      accepted=accepted+1;sector_vectors(:,accepted)=trial/norm_value
      if(accepted==multiplicity)exit
    enddo
    if(accepted/=multiplicity)then;message='translation-sector complex reconstruction lost rank';return;endif
    sector_rank=multiplicity
    if(multiplicity>0.and.multiplicity>huge(0)/multiplicity)then
      message='translation-sector Gram MPI count overflows';return
    endif
    allocate(gram(multiplicity,multiplicity),stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector Gram allocation failed';return;endif
    gram=matmul(conjg(transpose(sector_vectors)),sector_vectors)
    call MPI_Allreduce(MPI_IN_PLACE,gram,multiplicity*multiplicity,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='translation-sector Gram reduction failed';return;endif
    local_defect=0d0
    do j=1,multiplicity;do i=1,multiplicity
      local_defect=max(local_defect,abs(gram(i,j)-merge(1d0,0d0,i==j)))
    enddo;enddo
    unitarity_defect=max(unitarity_defect,local_defect)
    if(local_defect>10d0*tolerance)then
      message='translation-sector frame is not orthonormal';return
    endif
    allocate(sector_projector_rows(nlocal,n),gamma_projector_rows(nlocal,n),remote_sector(multiplicity),&
      stat=allocation_status)
    call allocation_consensus(allocation_status,global_bad)
    if(global_bad/=0)then;message='translation-sector projector allocation failed';return;endif
    do j=1,n
      call broadcast_sector_row(j,remote_sector,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation sector projector stream failed';return;endif
      do i=1,nlocal
        sector_projector_rows(i,j)=sum(sector_vectors(i,:)*conjg(remote_sector))
      enddo
    enddo
    call distributed_product(gamma_rows,conjg(sector_projector_rows),product_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma projector left stream failed';return;endif
    call distributed_adjoint(gamma_rows,adjoint_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma projector adjoint stream failed';return;endif
    call distributed_product(product_rows,adjoint_rows,gamma_projector_rows,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Gamma projector right stream failed';return;endif
    if(character_conjugates(requested_character)==requested_character)then
      local_defect=maxval(abs(gamma_projector_rows-sector_projector_rows))
      call MPI_Allreduce(local_defect,norm_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Gamma projector defect reduction failed';return;endif
      gamma_pairing_defect=max(gamma_pairing_defect,norm_value)
    endif
    if(gamma_pairing_defect>10d0*tolerance)then
      message='translation-sector Gamma projector is not self-conjugate';return
    endif
    fingerprint=ieor(1469598103934665603_int64,catalog_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,13),int(requested_character,int64))
    call MPI_Comm_rank(comm,rank,ierr)
    do i=1,n
      call broadcast_row(sector_projector_rows,i,stream_row,ierr)
      if(ierr/=MPI_SUCCESS)then;message='translation-sector fingerprint row stream failed';return;endif
      if(rank==0)then
        do j=1,n
          fingerprint=ieor(fingerprint,nint(real(stream_row(j),real64)/(100d0*tolerance),int64))
          fingerprint=ishftc(fingerprint,13)
          fingerprint=ieor(fingerprint,nint(aimag(stream_row(j))/(100d0*tolerance),int64))
          fingerprint=ishftc(fingerprint,13)
        enddo
      endif
    enddo
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='translation-sector fingerprint broadcast failed';return;endif
    call compute_workspace_receipt(receipt_ok)
    if(.not.receipt_ok)then;message='translation-sector workspace receipt overflows';return;endif
    ok=.true.;message=''
  contains
    subroutine broadcast_row(rows,global_row,values,error)
      complex(real64),intent(in)::rows(:,:)
      integer,intent(in)::global_row
      complex(real64),intent(out)::values(:)
      integer,intent(out)::error
      values=(0d0,0d0)
      if(rank==row_owner(global_row)-1)values=rows(row_position(global_row),:)
      call MPI_Bcast(values,size(values),MPI_DOUBLE_COMPLEX,row_owner(global_row)-1,comm,error)
    end subroutine

    subroutine distributed_adjoint(rows,adjoint,error)
      complex(real64),intent(in)::rows(:,:)
      complex(real64),intent(out)::adjoint(:,:)
      integer,intent(out)::error
      complex(real64),allocatable::local_column(:)
      integer::ii,jj,global_target
      allocate(local_column(n),stat=allocation_status)
      call allocation_consensus(allocation_status,global_bad)
      if(global_bad/=0)then;error=1;return;endif
      error=MPI_SUCCESS
      do global_target=1,n
        local_column=(0d0,0d0)
        do jj=1,nlocal;local_column(int(row_ids(jj)))=conjg(rows(jj,global_target));enddo
        call MPI_Allreduce(MPI_IN_PLACE,local_column,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,error)
        if(error/=MPI_SUCCESS)return
        if(rank==row_owner(global_target)-1)then
          ii=row_position(global_target);adjoint(ii,:)=local_column
        endif
      enddo
    end subroutine

    subroutine distributed_product(left_rows,right_rows,result_rows,error)
      complex(real64),intent(in)::left_rows(:,:),right_rows(:,:)
      complex(real64),intent(out)::result_rows(:,:)
      integer,intent(out)::error
      integer::ii,kk
      result_rows=(0d0,0d0);error=MPI_SUCCESS
      do kk=1,n
        call broadcast_row(right_rows,kk,stream_row,error);if(error/=MPI_SUCCESS)return
        do ii=1,nlocal;result_rows(ii,:)=result_rows(ii,:)+left_rows(ii,kk)*stream_row;enddo
      enddo
    end subroutine

    subroutine put_realified_entry(global_row,global_column,value)
      integer,intent(in)::global_row,global_column
      real(real64),intent(in)::value
      integer::local_row,local_column
      if(eigen_owner_node(global_row,real_info%nprow,real_info%myrow)/=real_info%myrow)return
      if(eigen_owner_node(global_column,real_info%npcol,real_info%mycol)/=real_info%mycol)return
      local_row=eigen_translate_g2l(global_row,real_info%nprow,real_info%myrow)
      local_column=eigen_translate_g2l(global_column,real_info%npcol,real_info%mycol)
      real_matrix(local_row,local_column)=value
    end subroutine

    subroutine broadcast_sector_row(global_row,values,error)
      integer,intent(in)::global_row
      complex(real64),intent(out)::values(:)
      integer,intent(out)::error
      values=(0d0,0d0)
      if(rank==row_owner(global_row)-1)values=sector_vectors(row_position(global_row),:)
      call MPI_Bcast(values,size(values),MPI_DOUBLE_COMPLEX,row_owner(global_row)-1,comm,error)
    end subroutine

    subroutine allocation_consensus(local_status,global_status)
      integer,intent(in)::local_status
      integer,intent(out)::global_status
      integer::error
      call MPI_Allreduce(local_status,global_status,1,MPI_INTEGER,MPI_MAX,comm,error)
      if(error/=MPI_SUCCESS)global_status=max(1,global_status)
    end subroutine

    subroutine compute_workspace_receipt(receipt_valid)
      logical,intent(out)::receipt_valid
      integer(int64)::complex_elements,real_elements,integer_elements,complex_bytes,real_bytes,integer_bytes
      complex_elements=0_int64;real_elements=0_int64;integer_elements=0_int64
      receipt_valid=.true.
      call add_count(complex_elements,size(generator_rows,kind=int64),receipt_valid)
      call add_count(complex_elements,size(gamma_rows,kind=int64),receipt_valid)
      call add_count(complex_elements,size(characters,kind=int64),receipt_valid)
      if(allocated(identity_rows))call add_count(complex_elements,size(identity_rows,kind=int64),receipt_valid)
      if(allocated(power_rows))call add_count(complex_elements,size(power_rows,kind=int64),receipt_valid)
      if(allocated(product_rows))call add_count(complex_elements,size(product_rows,kind=int64),receipt_valid)
      if(allocated(left_product_rows))call add_count(complex_elements,size(left_product_rows,kind=int64),receipt_valid)
      if(allocated(adjoint_rows))call add_count(complex_elements,size(adjoint_rows,kind=int64),receipt_valid)
      if(allocated(discriminator_rows))call add_count(complex_elements,size(discriminator_rows,kind=int64),receipt_valid)
      if(allocated(candidate_rows))call add_count(complex_elements,size(candidate_rows,kind=int64),receipt_valid)
      if(allocated(sector_vectors))call add_count(complex_elements,size(sector_vectors,kind=int64),receipt_valid)
      if(allocated(sector_projector_rows))call add_count(complex_elements,size(sector_projector_rows,kind=int64),receipt_valid)
      if(allocated(gamma_projector_rows))call add_count(complex_elements,size(gamma_projector_rows,kind=int64),receipt_valid)
      if(allocated(gram))call add_count(complex_elements,size(gram,kind=int64),receipt_valid)
      if(allocated(sector_multiplicities))&
        call add_count(complex_elements,size(sector_multiplicities,kind=int64),receipt_valid)
      if(allocated(stream_row))call add_count(complex_elements,size(stream_row,kind=int64),receipt_valid)
      if(allocated(remote_sector))call add_count(complex_elements,size(remote_sector,kind=int64),receipt_valid)
      ! One N-element column is allocated transiently inside distributed_adjoint.
      call add_count(complex_elements,int(n,int64),receipt_valid)
      if(allocated(real_matrix))then
        call add_count(real_elements,size(real_matrix,kind=int64),receipt_valid)
        ! Conservative allowance for two same-sized EigenExa internal distributed blocks.
        call add_count(real_elements,size(real_matrix,kind=int64),receipt_valid)
        call add_count(real_elements,size(real_matrix,kind=int64),receipt_valid)
      endif
      if(allocated(real_vectors))call add_count(real_elements,size(real_vectors,kind=int64),receipt_valid)
      if(allocated(eigenvalues))call add_count(real_elements,size(eigenvalues,kind=int64),receipt_valid)
      if(allocated(real_column))call add_count(real_elements,size(real_column,kind=int64),receipt_valid)
      call add_count(integer_elements,size(row_ids,kind=int64),receipt_valid)
      call add_count(integer_elements,size(generator_orders,kind=int64),receipt_valid)
      call add_count(integer_elements,size(element_words,kind=int64),receipt_valid)
      call add_count(integer_elements,size(character_conjugates,kind=int64),receipt_valid)
      if(allocated(ownership))call add_count(integer_elements,size(ownership,kind=int64),receipt_valid)
      if(allocated(row_owner))call add_count(integer_elements,size(row_owner,kind=int64),receipt_valid)
      if(allocated(row_position))call add_count(integer_elements,size(row_position,kind=int64),receipt_valid)
      if(.not.receipt_valid)return
      if(complex_elements>huge(0_int64)/16_int64.or.real_elements>huge(0_int64)/8_int64.or.&
          integer_elements>huge(0_int64)/8_int64)then;receipt_valid=.false.;return;endif
      complex_bytes=16_int64*complex_elements;real_bytes=8_int64*real_elements
      integer_bytes=8_int64*integer_elements
      if(complex_bytes>huge(0_int64)-real_bytes)then;receipt_valid=.false.;return;endif
      workspace_peak_bytes=complex_bytes+real_bytes
      if(workspace_peak_bytes>huge(0_int64)-integer_bytes)then;receipt_valid=.false.;return;endif
      workspace_peak_bytes=workspace_peak_bytes+integer_bytes
      receipt_valid=workspace_peak_bytes>0_int64
    end subroutine

    subroutine add_count(total,count,valid)
      integer(int64),intent(inout)::total
      integer(int64),intent(in)::count
      logical,intent(inout)::valid
      if(.not.valid)return
      if(count<0_int64.or.total>huge(0_int64)-count)then;valid=.false.;return;endif
      total=total+count
    end subroutine
  end subroutine split_dg_translation_character_sector_eigenexa

  subroutine build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      translation_target_box_ids,representative_target_box_ids,point_product,translation_cocycle,&
      identity_operation,requested_count,tolerance,candidates,spectrum,candidate_rank,projector_trace,&
      closure_residual,gamma_real_defect,workspace_peak_bytes,ok,message,selected_edge,rejected_edge,cluster_gap)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::comm,point_product(:,:),translation_cocycle(:,:),identity_operation,requested_count
    complex(real64),intent(in)::occupied(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),intent(in)::translation_target_box_ids(:,:),representative_target_box_ids(:,:)
    complex(real64),allocatable,intent(out)::candidates(:,:)
    real(real64),allocatable,intent(out)::spectrum(:)
    integer,intent(out)::candidate_rank
    real(real64),intent(out)::projector_trace,closure_residual,gamma_real_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),intent(out),optional::selected_edge,rejected_edge,cluster_gap

    call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      representative_target_box_ids,point_product,identity_operation,requested_count,tolerance,&
      candidates,spectrum,candidate_rank,projector_trace,closure_residual,gamma_real_defect,&
      workspace_peak_bytes,ok,message,selected_edge,rejected_edge,cluster_gap,&
      translation_target_box_ids,translation_cocycle)
  end subroutine
#endif

  subroutine select_dg_fixed_rank_symmetry_closed_subspace(metric,occupied,localizer,&
      representation,product_table,target_rank,tolerance,transform,occupied_inclusion,&
      subspace_leakage,ok,message)
    complex(real64),intent(in)::metric(:,:),occupied(:,:),localizer(:,:),representation(:,:,:)
    integer,intent(in)::product_table(:,:),target_rank
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::transform(:,:)
    real(real64),intent(out)::occupied_inclusion,subspace_leakage
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::metric_vectors(:,:),metric_sqrt(:,:),metric_inverse_sqrt(:,:),&
      orthogonal_representation(:,:,:),orthogonal_occupied(:,:),occupied_projector(:,:),&
      identity(:,:),difference(:,:),orthogonal_localizer(:,:),averaged_localizer(:,:),&
      complement_projector(:,:),complement_vectors(:,:),complement_basis(:,:),&
      complement_localizer(:,:),complement_eigenvectors(:,:),selected_orthogonal(:,:),&
      selected_projector(:,:),work(:,:)
    real(real64),allocatable::metric_spectrum(:),complement_projector_spectrum(:),&
      complement_spectrum(:),occupied_spectrum(:)
    complex(real64),allocatable::occupied_vectors(:,:),occupied_inverse_sqrt(:,:)
    real(real64)::scale,defect,boundary_scale
    integer::n,noccupied,noperation,ncomplement,nselect,i,j,operation,left,right,product,column
    logical::eigen_ok
    character(256)::detail

    ok=.false.;message='';occupied_inclusion=huge(1d0);subspace_leakage=huge(1d0)
    n=size(metric,1);noccupied=size(occupied,2);noperation=size(representation,3)
    if(n<1.or.size(metric,2)/=n.or.size(localizer,1)/=n.or.size(localizer,2)/=n.or.&
        size(occupied,1)/=n.or.noccupied<1.or.target_rank<noccupied.or.target_rank>n.or.&
        size(representation,1)/=n.or.size(representation,2)/=n.or.noperation<1.or.&
        any(shape(product_table)/=[noperation,noperation]).or.tolerance<=0d0.or.&
        .not.ieee_is_finite(tolerance))then
      message='invalid fixed-rank symmetry-closed subspace contract';return
    end if
    if(.not.all(ieee_is_finite(real(metric))).or..not.all(ieee_is_finite(aimag(metric))).or.&
        .not.all(ieee_is_finite(real(occupied))).or..not.all(ieee_is_finite(aimag(occupied))).or.&
        .not.all(ieee_is_finite(real(localizer))).or..not.all(ieee_is_finite(aimag(localizer))).or.&
        .not.all(ieee_is_finite(real(representation))).or.&
        .not.all(ieee_is_finite(aimag(representation))))then
      message='nonfinite fixed-rank symmetry-closed subspace input';return
    end if
    scale=max(1d0,maxval(abs(metric)))
    if(maxval(abs(metric-conjg(transpose(metric))))>tolerance*scale.or.&
        maxval(abs(localizer-conjg(transpose(localizer))))>&
        tolerance*max(1d0,maxval(abs(localizer))))then
      message='fixed-rank metric/localizer is not Hermitian';return
    end if
    call hermitian_eigensystem(metric,metric_spectrum,metric_vectors,eigen_ok,detail)
    if(.not.eigen_ok.or.minval(metric_spectrum)<=tolerance*maxval(metric_spectrum))then
      message='fixed-rank candidate metric is not positive definite';return
    end if
    allocate(metric_sqrt(n,n),metric_inverse_sqrt(n,n),identity(n,n),difference(n,n),&
      orthogonal_representation(n,n,noperation),orthogonal_occupied(n,noccupied),&
      occupied_projector(n,n),orthogonal_localizer(n,n),averaged_localizer(n,n),&
      complement_projector(n,n),selected_projector(n,n),work(n,n))
    metric_sqrt=metric_vectors;metric_inverse_sqrt=metric_vectors
    do i=1,n
      metric_sqrt(:,i)=sqrt(metric_spectrum(i))*metric_sqrt(:,i)
      metric_inverse_sqrt(:,i)=metric_inverse_sqrt(:,i)/sqrt(metric_spectrum(i))
    end do
    metric_sqrt=matmul(metric_sqrt,conjg(transpose(metric_vectors)))
    metric_inverse_sqrt=matmul(metric_inverse_sqrt,conjg(transpose(metric_vectors)))
    identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;end do
    orthogonal_occupied=matmul(metric_sqrt,occupied)
    difference(1:noccupied,1:noccupied)=&
      matmul(conjg(transpose(orthogonal_occupied)),orthogonal_occupied)
    call hermitian_eigensystem(difference(1:noccupied,1:noccupied),occupied_spectrum,&
      occupied_vectors,eigen_ok,detail)
    if(.not.eigen_ok.or.minval(occupied_spectrum)<=tolerance*maxval(occupied_spectrum))then
      message='fixed-rank occupied coefficients are metric rank deficient';return
    end if
    occupied_inverse_sqrt=occupied_vectors
    do i=1,noccupied
      occupied_inverse_sqrt(:,i)=occupied_inverse_sqrt(:,i)/sqrt(occupied_spectrum(i))
    end do
    occupied_inverse_sqrt=matmul(occupied_inverse_sqrt,conjg(transpose(occupied_vectors)))
    orthogonal_occupied=matmul(orthogonal_occupied,occupied_inverse_sqrt)
    occupied_projector=matmul(orthogonal_occupied,conjg(transpose(orthogonal_occupied)))
    do operation=1,noperation
      orthogonal_representation(:,:,operation)=matmul(metric_sqrt,&
        matmul(representation(:,:,operation),metric_inverse_sqrt))
      defect=maxval(abs(matmul(conjg(transpose(orthogonal_representation(:,:,operation))),&
        orthogonal_representation(:,:,operation))-identity))
      if(defect>tolerance)then;message='fixed-rank symmetry representation is not metric unitary';return;end if
      defect=maxval(abs(matmul(orthogonal_representation(:,:,operation),occupied_projector)-&
        matmul(occupied_projector,orthogonal_representation(:,:,operation))))
      if(defect>tolerance)then;message='occupied subspace is not closed under full-system symmetry';return;end if
    end do
    do left=1,noperation;do right=1,noperation
      product=product_table(left,right)
      if(product<1.or.product>noperation)then;message='fixed-rank symmetry product is invalid';return;end if
      defect=maxval(abs(matmul(orthogonal_representation(:,:,left),&
        orthogonal_representation(:,:,right))-orthogonal_representation(:,:,product)))
      if(defect>tolerance)then;message='fixed-rank symmetry representation is not group closed';return;end if
    end do;end do
    orthogonal_localizer=matmul(metric_inverse_sqrt,matmul(localizer,metric_inverse_sqrt))
    averaged_localizer=(0d0,0d0)
    do operation=1,noperation
      averaged_localizer=averaged_localizer+matmul(orthogonal_representation(:,:,operation),&
        matmul(orthogonal_localizer,conjg(transpose(orthogonal_representation(:,:,operation)))))
    end do
    averaged_localizer=averaged_localizer/real(noperation,real64)
    averaged_localizer=0.5d0*(averaged_localizer+conjg(transpose(averaged_localizer)))
    complement_projector=identity-occupied_projector
    call hermitian_eigensystem(complement_projector,complement_projector_spectrum,&
      complement_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;message='fixed-rank occupied complement diagonalization failed';return;end if
    ncomplement=n-noccupied
    if(count(complement_projector_spectrum>0.5d0)/=ncomplement)then
      message='fixed-rank occupied complement has inconsistent dimension';return
    end if
    allocate(complement_basis(n,ncomplement))
    column=0
    do i=1,n
      if(complement_projector_spectrum(i)<=0.5d0)cycle
      column=column+1;complement_basis(:,column)=complement_vectors(:,i)
    end do
    complement_localizer=matmul(conjg(transpose(complement_basis)),&
      matmul(averaged_localizer,complement_basis))
    call hermitian_eigensystem(complement_localizer,complement_spectrum,&
      complement_eigenvectors,eigen_ok,detail)
    if(.not.eigen_ok)then;message='fixed-rank complement localizer diagonalization failed';return;end if
    nselect=target_rank-noccupied
    if(nselect>0.and.nselect<ncomplement)then
      boundary_scale=max(1d0,maxval(abs(complement_spectrum)))
      if(abs(complement_spectrum(nselect+1)-complement_spectrum(nselect))<=tolerance*boundary_scale)then
        message='target rank cuts a symmetry-degenerate block';return
      end if
    end if
    allocate(selected_orthogonal(n,target_rank));selected_orthogonal(:,1:noccupied)=orthogonal_occupied
    if(nselect>0)selected_orthogonal(:,noccupied+1:target_rank)=matmul(complement_basis,&
      complement_eigenvectors(:,1:nselect))
    allocate(transform(n,target_rank));transform=matmul(metric_inverse_sqrt,selected_orthogonal)
    selected_projector=matmul(selected_orthogonal,conjg(transpose(selected_orthogonal)))
    occupied_inclusion=maxval(abs(matmul(identity-selected_projector,orthogonal_occupied)))
    subspace_leakage=0d0
    do operation=1,noperation
      work=matmul(identity-selected_projector,&
        matmul(orthogonal_representation(:,:,operation),selected_projector))
      subspace_leakage=max(subspace_leakage,maxval(abs(work)))
    end do
    if(occupied_inclusion>tolerance)then;message='selected subspace lost occupied inclusion';return;end if
    if(subspace_leakage>tolerance)then;message='selected subspace is not symmetry closed';return;end if
    ok=.true.
  end subroutine

  subroutine build_dg_pointwise_affine_owner_map(global_grid,local_physical_ids,all_physical_ids,&
      integer_rotation,fractional_translation,tolerance,target_physical_ids,target_owner,&
      target_local_index,lattice_wrap,ok,message)
    integer,intent(in)::global_grid(3),integer_rotation(3,3)
    integer(int64),intent(in)::local_physical_ids(:),all_physical_ids(:,:)
    real(real64),intent(in)::fractional_translation(3),tolerance
    integer(int64),allocatable,intent(out)::target_physical_ids(:)
    integer,allocatable,intent(out)::target_owner(:),target_local_index(:),lattice_wrap(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::global_count,source_id,target_id,rotation_determinant
    integer::nlocal,nowner,point,axis,input_axis,source_grid(3),mapped_grid(3),location(2)
    real(real64)::translation_grid(3),mapped_coordinate,nearest

    ok=.false.;message='';nlocal=size(local_physical_ids);nowner=size(all_physical_ids,2)
    if(any(global_grid<1).or.nlocal<1.or.nowner<1.or.size(all_physical_ids,1)<1.or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(fractional_translation)))then
      message='invalid pointwise affine owner-map contract';return
    end if
    rotation_determinant=int(integer_rotation(1,1),int64)*(&
      int(integer_rotation(2,2),int64)*int(integer_rotation(3,3),int64)-&
      int(integer_rotation(2,3),int64)*int(integer_rotation(3,2),int64))-&
      int(integer_rotation(1,2),int64)*(&
      int(integer_rotation(2,1),int64)*int(integer_rotation(3,3),int64)-&
      int(integer_rotation(2,3),int64)*int(integer_rotation(3,1),int64))+&
      int(integer_rotation(1,3),int64)*(&
      int(integer_rotation(2,1),int64)*int(integer_rotation(3,2),int64)-&
      int(integer_rotation(2,2),int64)*int(integer_rotation(3,1),int64))
    if(abs(rotation_determinant)/=1_int64)then
      message='affine rotation must be unimodular';return
    end if
    global_count=int(global_grid(1),int64)*int(global_grid(2),int64)*int(global_grid(3),int64)
    if(global_count<1_int64.or.any(local_physical_ids<1_int64).or.&
        any(local_physical_ids>global_count).or.any(all_physical_ids<1_int64).or.&
        any(all_physical_ids>global_count))then
      message='pointwise affine owner-map physical ID is outside the global grid';return
    end if
    if(size(all_physical_ids)/=int(global_count))then
      message='pointwise affine owner table does not cover the global grid';return
    end if
    do axis=1,3
      translation_grid(axis)=fractional_translation(axis)*real(global_grid(axis),real64)
      if(abs(translation_grid(axis)-anint(translation_grid(axis)))>tolerance)then
        message='affine translation is incommensurate with the global grid';return
      end if
    end do
    allocate(target_physical_ids(nlocal),target_owner(nlocal),target_local_index(nlocal),&
      lattice_wrap(3,nlocal))
    do point=1,nlocal
      source_id=local_physical_ids(point)-1_int64
      source_grid(1)=int(modulo(source_id,int(global_grid(1),int64)))
      source_grid(2)=int(modulo(source_id/int(global_grid(1),int64),int(global_grid(2),int64)))
      source_grid(3)=int(source_id/(int(global_grid(1),int64)*int(global_grid(2),int64)))
      do axis=1,3
        mapped_coordinate=translation_grid(axis)
        do input_axis=1,3
          mapped_coordinate=mapped_coordinate+real(integer_rotation(axis,input_axis),real64)*&
            real(source_grid(input_axis),real64)*real(global_grid(axis),real64)/&
            real(global_grid(input_axis),real64)
        end do
        nearest=anint(mapped_coordinate)
        if(abs(mapped_coordinate-nearest)>tolerance)then
          message='affine rotation is incommensurate with the global grid';return
        end if
        lattice_wrap(axis,point)=floor(nearest/real(global_grid(axis),real64))
        mapped_grid(axis)=modulo(int(nearest),global_grid(axis))
      end do
      target_id=int(mapped_grid(1),int64)+int(global_grid(1),int64)*(&
        int(mapped_grid(2),int64)+int(global_grid(2),int64)*int(mapped_grid(3),int64))+1_int64
      if(count(all_physical_ids==target_id)/=1)then
        message='mapped global grid point does not have exactly one owner';return
      end if
      location=findloc(all_physical_ids,target_id)
      target_physical_ids(point)=target_id
      target_local_index(point)=location(1);target_owner(point)=location(2)-1
    end do
    ok=.true.
  end subroutine
  subroutine verify_dg_fragment_subspace_density_covariance(values,target_ids,tolerance,ok,message)
    complex(real64),intent(in)::values(:,:)
    integer(int64),intent(in)::target_ids(:,:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::density(:)
    real(real64)::scale
    integer::operation,point,target
    ok=.false.;message=''
    if(size(values,1)<1.or.size(values,2)<1.or.size(target_ids,1)/=size(values,2).or. &
        size(target_ids,2)<1.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance).or. &
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))))then
      message='invalid Wannier subspace-density covariance contract';return
    end if
    if(any(target_ids<1_int64).or.any(target_ids>int(size(values,2),int64)))then
      message='Wannier subspace-density symmetry target is out of range';return
    end if
    allocate(density(size(values,2)));density=sum(abs(values)**2,dim=1)
    scale=max(1d0,maxval(density))
    do operation=1,size(target_ids,2);do point=1,size(values,2)
      target=int(target_ids(point,operation))
      if(abs(density(target)-density(point))>tolerance*scale)then
        message='periodic-box Wannier subspace density is not symmetry covariant';return
      end if
    end do;end do
    ok=.true.
  end subroutine verify_dg_fragment_subspace_density_covariance

  subroutine assign_dg_overlapping_wannier_occupations(electron_count,occupations,ok,message)
    real(real64),intent(in)::electron_count
    real(real64),intent(out)::occupations(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::fully_occupied
    real(real64)::remainder

    occupations=0d0
    ok=size(occupations)>0.and.electron_count>=0d0.and.ieee_is_finite(electron_count).and.&
      electron_count<=2d0*real(size(occupations),real64)+10d0*epsilon(1d0)
    if(.not.ok)then
      message='overlapping-Wannier: insufficient global bands for electron count'
      return
    endif
    fully_occupied=min(size(occupations),int(electron_count/2d0))
    if(fully_occupied>0)occupations(1:fully_occupied)=2d0
    remainder=electron_count-2d0*real(fully_occupied,real64)
    if(remainder>10d0*epsilon(1d0).and.fully_occupied<size(occupations))&
      occupations(fully_occupied+1)=remainder
    ok=all(occupations>=0d0).and.all(occupations<=2d0).and.&
      abs(sum(occupations)-electron_count)<=10d0*epsilon(1d0)*max(1d0,electron_count)
    if(ok)then
      message=''
    else
      occupations=0d0
      message='overlapping-Wannier: invalid global occupations'
    endif
  end subroutine

  subroutine verify_dg_uniform_fragment_target_rank(comm,local_target_rank,ok,message)
    integer,intent(in)::comm,local_target_rank
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::minimum_rank,maximum_rank,ierr
    call MPI_Allreduce(local_target_rank,minimum_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      ok=.false.;message='cannot reduce minimum fragment target rank';return
    endif
    call MPI_Allreduce(local_target_rank,maximum_rank,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      ok=.false.;message='cannot reduce maximum fragment target rank';return
    endif
    ok=local_target_rank>0.and.minimum_rank==maximum_rank
#else
    ok=local_target_rank>0
#endif
    if(ok)then
      message=''
    else
      message='fragment target rank differs across DC fragments'
    endif
  end subroutine

  subroutine verify_dg_fragment_center_orbit(local_center_box_ids,global_center_box_ids,&
      symmetry_target_box_ids,ok,message)
    integer(int64),intent(in)::local_center_box_ids(:),global_center_box_ids(:),&
      symmetry_target_box_ids(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::center,operation
    ok=size(local_center_box_ids)>0.and.size(symmetry_target_box_ids,1)>0.and.&
      size(symmetry_target_box_ids,2)>0.and.&
      size(global_center_box_ids)==size(local_center_box_ids)*size(symmetry_target_box_ids,2)
    if(ok)ok=all(local_center_box_ids>=1_int64).and.&
      all(local_center_box_ids<=int(size(symmetry_target_box_ids,1),int64))
    if(ok)then
      do operation=1,size(symmetry_target_box_ids,2)
        do center=1,size(local_center_box_ids)
          if(count(global_center_box_ids(center::size(local_center_box_ids))==symmetry_target_box_ids(&
              int(local_center_box_ids(center)),operation))/=1)ok=.false.
        enddo
      enddo
    endif
    if(ok)then
      message=''
    else
      message='translated Wannier centers do not form a complete bond-center orbit'
    endif
  end subroutine

  subroutine build_dg_core_owned_occupied_subspace(candidate,core_mask,weights,occupations,&
      owned_electron_count,coefficients,core_electron_count,ok,message)
    complex(real64),intent(in)::candidate(:,:)
    logical,intent(in)::core_mask(:)
    real(real64),intent(in)::weights(:),occupations(:),owned_electron_count
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::core_electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::occupied_index(:)
    complex(real64),allocatable::core_gram(:,:),vectors(:,:)
    real(real64),allocatable::spectrum(:),core_norm(:)
    integer::ncandidate,nbox,nband,nowned,i,j,p

    ok=.false.;message='';core_electron_count=0d0
    ncandidate=size(candidate,1);nbox=size(candidate,2)
    if(ncandidate<1.or.nbox<1.or.size(core_mask)/=nbox.or.size(weights)/=nbox.or.&
        size(occupations)/=ncandidate.or.any(weights<=0d0).or.&
        .not.ieee_is_finite(owned_electron_count).or.owned_electron_count<=0d0)then
      message='invalid core-owned occupied-subspace contract';return
    endif
    occupied_index=pack([(i,i=1,ncandidate)],occupations>1d-12)
    nband=size(occupied_index)
    if(nband<1)then;message='DC core-owned occupied subspace has no occupied bands';return;endif
    allocate(core_norm(ncandidate));core_norm=0d0
    do i=1,ncandidate
      do p=1,nbox
        if(core_mask(p))core_norm(i)=core_norm(i)+weights(p)*abs(candidate(i,p))**2
      enddo
      core_electron_count=core_electron_count+occupations(i)*core_norm(i)
    enddo
    nowned=nint(0.5d0*owned_electron_count)
    if(nowned<1.or.nowned>nband.or.abs(owned_electron_count-2d0*real(nowned,real64))>1d-10)then
      message='core-owned ionic valence does not define an integral occupied rank';return
    endif
    allocate(core_gram(nband,nband));core_gram=(0d0,0d0)
    do j=1,nband;do i=1,nband
      do p=1,nbox
        if(core_mask(p))core_gram(i,j)=core_gram(i,j)+weights(p)*&
          conjg(candidate(occupied_index(i),p))*candidate(occupied_index(j),p)
      enddo
    enddo;enddo
    core_gram=0.5d0*(core_gram+conjg(transpose(core_gram)))
    call hermitian_eigensystem(core_gram,spectrum,vectors,ok,message)
    if(.not.ok)return
    if(spectrum(nband-nowned+1)<=epsilon(1d0)*max(1d0,spectrum(nband)))then
      ok=.false.;message='DC core-owned occupied subspace is rank deficient';return
    endif
    allocate(coefficients(ncandidate,nowned));coefficients=(0d0,0d0)
    do j=1,nowned
      do i=1,nband
        coefficients(occupied_index(i),j)=vectors(i,nband-nowned+j)
      enddo
    enddo
    ok=.true.
  end subroutine

  subroutine verify_dg_fragment_wannier_streaming_closure(comm,fragment_id,local_target_count,&
      box_ids,symmetry_target_box_ids,values,gradients,tolerance,residual,fingerprint,ok,message)
    integer,intent(in)::comm,fragment_id,local_target_count
    integer(int64),intent(in)::box_ids(:),symmetry_target_box_ids(:,:)
    complex(real64),intent(in)::values(:,:),gradients(:,:,:)
    real(real64),intent(in)::tolerance
    real(real64),intent(out)::residual
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,allocatable::rank_fragment(:),target_rank_all(:,:)
    complex(real64),allocatable::send_buffer(:,:),receive_buffer(:,:),receive_gradients(:,:,:)
    integer::rank,nproc,ierr,nbox,nsym,ntarget,operation,owner,target_rank,preimage_rank,&
      target_owner,target_point,p,iw,axis,local_bad,global_bad,tag
    integer(int64)::local_hash,bits
    real(real64)::local_residual,source_density,target_density,scale

    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nbox=size(box_ids);nsym=size(symmetry_target_box_ids,2);ntarget=size(values,1)
    local_bad=merge(0,1,nproc>0.and.local_target_count>0.and.ntarget==local_target_count*nproc.and.&
      size(values,2)==nbox.and.all(shape(gradients)==[3,ntarget,nbox]).and.&
      size(symmetry_target_box_ids,1)==nbox.and.nsym==nproc.and.tolerance>0d0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      ok=.false.;message='invalid streaming fragment Wannier closure contract'
      residual=huge(1d0);fingerprint=0_int64;return
    endif
    allocate(rank_fragment(nproc),target_rank_all(nsym,nproc),&
      send_buffer(local_target_count,nbox),receive_buffer(local_target_count,nbox),&
      receive_gradients(3,local_target_count,nbox))
    call MPI_Allgather(fragment_id,1,MPI_INTEGER,rank_fragment,1,MPI_INTEGER,comm,ierr)
    do operation=1,nsym
      target_rank=int((symmetry_target_box_ids(1,operation)-1_int64)/int(nbox,int64))+1
      target_rank=findloc(rank_fragment,target_rank,dim=1)-1
      call MPI_Allgather(target_rank,1,MPI_INTEGER,target_rank_all(operation,:),1,MPI_INTEGER,comm,ierr)
    enddo
    local_bad=merge(0,1,all(target_rank_all>=0).and.all(target_rank_all<nproc))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      ok=.false.;message='fragment symmetry is not a rank permutation'
      residual=huge(1d0);fingerprint=0_int64;return
    endif
    local_residual=0d0;local_hash=0_int64
    do operation=1,nsym
      target_rank=target_rank_all(operation,rank+1)
      preimage_rank=findloc(target_rank_all(operation,:)==rank,.true.,dim=1)-1
      if(preimage_rank<0)then
        ok=.false.;message='fragment symmetry rank permutation has no inverse'
        residual=huge(1d0);fingerprint=0_int64;return
      endif
      do owner=0,nproc-1
        target_owner=target_rank_all(operation,owner+1)
        send_buffer=values(target_owner*local_target_count+1:(target_owner+1)*local_target_count,:)
        tag=operation*nproc+owner
        call MPI_Sendrecv(send_buffer,local_target_count*nbox,MPI_DOUBLE_COMPLEX,preimage_rank,tag,&
          receive_buffer,local_target_count*nbox,MPI_DOUBLE_COMPLEX,target_rank,tag,comm,MPI_STATUS_IGNORE,ierr)
        do p=1,nbox
          target_point=int(modulo(symmetry_target_box_ids(p,operation)-1_int64,int(nbox,int64)))+1
          source_density=sum(abs(values(owner*local_target_count+1:&
            (owner+1)*local_target_count,p))**2)
          target_density=sum(abs(receive_buffer(:,target_point))**2)
          scale=max(1d0,source_density,target_density)
          local_residual=max(local_residual,abs(source_density-target_density)/scale)
        enddo
        do axis=1,3
          send_buffer=gradients(axis,target_owner*local_target_count+1:&
            (target_owner+1)*local_target_count,:)
          tag=nsym*nproc+axis*nsym*nproc+operation*nproc+owner
          call MPI_Sendrecv(send_buffer,local_target_count*nbox,MPI_DOUBLE_COMPLEX,preimage_rank,tag,&
            receive_buffer,local_target_count*nbox,MPI_DOUBLE_COMPLEX,target_rank,tag,comm,MPI_STATUS_IGNORE,ierr)
          receive_gradients(axis,:,:)=receive_buffer
        enddo
        do p=1,nbox
          target_point=int(modulo(symmetry_target_box_ids(p,operation)-1_int64,int(nbox,int64)))+1
          source_density=sum(abs(gradients(:,owner*local_target_count+1:&
            (owner+1)*local_target_count,p))**2)
          target_density=sum(abs(receive_gradients(:,:,target_point))**2)
          scale=max(1d0,source_density,target_density)
          local_residual=max(local_residual,abs(source_density-target_density)/scale)
        enddo
        bits=transfer(real(receive_buffer(1,1),real64),bits)
        local_hash=ieor(local_hash,ishftc(bits,mod(11*operation+7*owner+rank,63)))
      enddo
    enddo
    call MPI_Allreduce(local_residual,residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(fingerprint,int(z'6A09E667F3BCC909',int64))
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=ieee_is_finite(residual).and.residual<=tolerance
    if(ok)then;message='';else;message='streaming fragment Wannier subspace-density closure failed';endif
#else
    residual=huge(1d0);fingerprint=0_int64;ok=.false.
    message='streaming fragment Wannier closure requires MPI'
#endif
  end subroutine

  subroutine replicate_dg_fragment_wannier_representative(comm,fragment_id,values,gradients,&
      residual,correction,ok,message)
    integer,intent(in)::comm,fragment_id
    complex(real64),intent(inout)::values(:,:),gradients(:,:,:)
    real(real64),intent(out)::residual,correction
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::reference_values(:,:),reference_gradients(:,:,:)
    integer::rank,ierr,nwann,nbox,local_bad,global_bad,pair(2),reference_pair(2),&
      local_shape(2),minimum_shape(2),maximum_shape(2)
    real(real64)::local_correction,local_residual

    ok=.false.;message='';residual=huge(1d0);correction=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr)
    nwann=size(values,1);nbox=size(values,2);local_bad=0
    if(fragment_id<1.or.nwann<1.or.nbox<1.or.any(shape(gradients)/=[3,nwann,nbox]))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='invalid representative fragment Wannier replication contract';return
    endif
    local_shape=[nwann,nbox]
    call MPI_Allreduce(local_shape,minimum_shape,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_shape,maximum_shape,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(any(minimum_shape/=maximum_shape))then
      message='representative fragment Wannier shape differs across ranks';return
    endif
    pair=[fragment_id,rank]
    call MPI_Allreduce(pair,reference_pair,1,MPI_2INTEGER,MPI_MINLOC,comm,ierr)
    allocate(reference_values(nwann,nbox),reference_gradients(3,nwann,nbox))
    if(rank==reference_pair(2))then
      reference_values=values;reference_gradients=gradients
    endif
    call MPI_Bcast(reference_values,nwann*nbox,MPI_DOUBLE_COMPLEX,reference_pair(2),comm,ierr)
    call MPI_Bcast(reference_gradients,3*nwann*nbox,MPI_DOUBLE_COMPLEX,reference_pair(2),comm,ierr)
    local_correction=max(maxval(abs(values-reference_values)),&
      maxval(abs(gradients-reference_gradients)))
    call MPI_Allreduce(local_correction,correction,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    values=reference_values;gradients=reference_gradients
    local_residual=max(maxval(abs(values-reference_values)),&
      maxval(abs(gradients-reference_gradients)))
    call MPI_Allreduce(local_residual,residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    ok=ieee_is_finite(residual).and.ieee_is_finite(correction).and.residual==0d0
    if(ok)then;message='';else;message='representative fragment Wannier replication failed';endif
#else
    residual=huge(1d0);correction=huge(1d0);ok=.false.
    message='representative fragment Wannier replication requires MPI'
#endif
  end subroutine

  subroutine align_dg_fragment_wannier_gauge(comm,weights,values,gradients,tolerance,&
      residual,correction,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(inout)::values(:,:),gradients(:,:,:)
    real(real64),intent(out)::residual,correction
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::reference_values(:,:),weighted_local(:,:),overlap(:,:),gram(:,:),&
      vectors(:,:),inverse_root(:,:),unitary(:,:),rotated_values(:,:),rotated_gradient(:,:)
    real(real64),allocatable::spectrum(:)
    integer::rank,nproc,ierr,nwann,nbox,p,j,axis,local_bad,global_bad
    real(real64)::local_residual

    ok=.false.;message='';residual=huge(1d0);correction=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nwann=size(values,1);nbox=size(values,2);local_bad=0
    if(nproc<1.or.nwann<1.or.nbox<1.or.size(weights)/=nbox.or.&
        any(shape(gradients)/=[3,nwann,nbox]).or.tolerance<=0d0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid fragment Wannier gauge-alignment contract';return;endif
    allocate(reference_values(nwann,nbox),weighted_local(nbox,nwann),overlap(nwann,nwann),&
      gram(nwann,nwann),inverse_root(nwann,nwann),unitary(nwann,nwann),&
      rotated_values(nwann,nbox),rotated_gradient(nwann,nbox))
    if(rank==0)reference_values=values
    call MPI_Bcast(reference_values,nwann*nbox,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    do p=1,nbox
      weighted_local(p,:)=weights(p)*values(:,p)
    enddo
    call zgemm('C','T',nwann,nwann,nbox,(1d0,0d0),weighted_local,nbox,reference_values,nwann,&
      (0d0,0d0),overlap,nwann)
    gram=matmul(conjg(transpose(overlap)),overlap)
    gram=0.5d0*(gram+conjg(transpose(gram)))
    call hermitian_eigensystem(gram,spectrum,vectors,ok,message)
    if(.not.ok)return
    if(minval(spectrum)<=0d0)then;ok=.false.;message='singular fragment Wannier gauge overlap';return;endif
    inverse_root=vectors
    do j=1,nwann
      inverse_root(:,j)=inverse_root(:,j)/sqrt(spectrum(j))
    enddo
    inverse_root=matmul(inverse_root,conjg(transpose(vectors)))
    unitary=matmul(overlap,inverse_root)
    correction=maxval(abs(unitary-overlap))
    rotated_values=matmul(transpose(unitary),values)
    values=rotated_values
    do axis=1,3
      rotated_gradient=matmul(transpose(unitary),gradients(axis,:,:))
      gradients(axis,:,:)=rotated_gradient
    enddo
    local_residual=maxval(abs(values-reference_values))
    call MPI_Allreduce(local_residual,residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(.not.ieee_is_finite(residual).or.residual>tolerance)then
      ok=.false.;message='fragment Wannier gauges do not close under periodic translation';return
    endif
    ok=.true.
#else
    residual=huge(1d0);correction=huge(1d0);ok=.false.
    message='fragment Wannier gauge alignment requires MPI'
#endif
  end subroutine

  subroutine assemble_dg_distributed_candidate_symmetry(comm,local_candidate,weights,&
      symmetry_target_box_ids,symmetry_overlap,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::local_candidate(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    complex(real64),allocatable,intent(out)::symmetry_overlap(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::broadcast_candidate(:,:),mapped_candidate(:,:),&
      local_overlap(:,:),global_overlap(:,:),overlap_block(:,:)
    integer::rank,nproc,ierr,nlocal,ncandidate,nglobal,nsym,isym,owner,p,target_rank,target_point,&
      source_offset,target_offset,local_bad,global_bad,allocation_status

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    ncandidate=size(local_candidate,1);nlocal=size(local_candidate,2)
    nsym=size(symmetry_target_box_ids,2)
    local_bad=merge(0,1,ncandidate>0.and.nlocal>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.nsym>0)
    if(ncandidate>0.and.nproc>huge(nglobal)/ncandidate)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed candidate symmetry contract';return;endif
    nglobal=ncandidate*nproc
    allocate(symmetry_overlap(nglobal,nglobal,nsym),broadcast_candidate(ncandidate,nlocal),&
      mapped_candidate(nlocal,ncandidate),local_overlap(nglobal,nglobal),&
      global_overlap(nglobal,nglobal),overlap_block(ncandidate,ncandidate),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='cannot allocate distributed candidate symmetry workspace';return;endif
    symmetry_overlap=(0d0,0d0)
    source_offset=rank*ncandidate
    do isym=1,nsym
      target_rank=int((symmetry_target_box_ids(1,isym)-1_int64)/int(nlocal,int64))
      local_bad=merge(0,1,target_rank>=0.and.target_rank<nproc)
      do p=1,nlocal
        if(int((symmetry_target_box_ids(p,isym)-1_int64)/int(nlocal,int64))/=target_rank)local_bad=1
        target_point=int(modulo(symmetry_target_box_ids(p,isym)-1_int64,int(nlocal,int64)))+1
        if(target_point<1.or.target_point>nlocal)local_bad=1
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0)then;message='symmetry does not map a fragment box to one fragment box';return;endif
      mapped_candidate=(0d0,0d0)
      do owner=0,nproc-1
        if(rank==owner)broadcast_candidate=local_candidate
        call MPI_Bcast(broadcast_candidate,ncandidate*nlocal,MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        if(owner/=target_rank)cycle
        do p=1,nlocal
          target_point=int(modulo(symmetry_target_box_ids(p,isym)-1_int64,int(nlocal,int64)))+1
          mapped_candidate(p,:)=weights(p)*broadcast_candidate(:,target_point)
        enddo
      enddo
      local_overlap=(0d0,0d0);target_offset=target_rank*ncandidate
      call zgemm('C','T',ncandidate,ncandidate,nlocal,(1d0,0d0),mapped_candidate,nlocal,&
        local_candidate,ncandidate,(0d0,0d0),overlap_block,ncandidate)
      local_overlap(target_offset+1:target_offset+ncandidate,&
        source_offset+1:source_offset+ncandidate)=overlap_block
      call MPI_Allreduce(local_overlap,global_overlap,nglobal*nglobal,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      symmetry_overlap(:,:,isym)=global_overlap
    enddo
    ok=.true.
#else
    ok=.false.;message='distributed candidate symmetry requires MPI'
#endif
  end subroutine

  subroutine assemble_dg_distributed_basis_symmetry_overlap(comm,local_basis,weights,&
      symmetry_target_box_ids,symmetry_overlap,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::local_basis(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    complex(real64),allocatable,intent(out)::symmetry_overlap(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::target_basis(:,:),mapped_basis(:,:),local_overlap(:,:)
    integer::rank,nproc,ierr,nbasis,nlocal,nsym,isym,p,owner,target_rank,target_point,&
      local_bad,global_bad,allocation_status

    ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nbasis=size(local_basis,1);nlocal=size(local_basis,2);nsym=size(symmetry_target_box_ids,2)
    local_bad=merge(0,1,nbasis>0.and.nlocal>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.nsym>0.and.&
      all(ieee_is_finite(weights)).and.all(weights>=0d0).and.&
      all(ieee_is_finite(real(local_basis))).and.all(ieee_is_finite(aimag(local_basis))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed full-basis symmetry-overlap contract';return;endif
    allocate(symmetry_overlap(nbasis,nbasis,nsym),target_basis(nbasis,nlocal),&
      mapped_basis(nlocal,nbasis),local_overlap(nbasis,nbasis),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='cannot allocate distributed full-basis symmetry-overlap workspace';return;endif
    do isym=1,nsym
      local_bad=0
      do p=1,nlocal
        target_rank=int((symmetry_target_box_ids(p,isym)-1_int64)/int(nlocal,int64))
        target_point=int(modulo(symmetry_target_box_ids(p,isym)-1_int64,int(nlocal,int64)))+1
        if(target_rank<0.or.target_rank>=nproc.or.target_point<1.or.target_point>nlocal)local_bad=1
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0)then;message='symmetry core-point owner map is invalid';return;endif
      mapped_basis=(0d0,0d0)
      do owner=0,nproc-1
        if(rank==owner)target_basis=local_basis
        call MPI_Bcast(target_basis,size(target_basis),MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        do p=1,nlocal
          target_rank=int((symmetry_target_box_ids(p,isym)-1_int64)/int(nlocal,int64))
          if(target_rank/=owner)cycle
          target_point=int(modulo(symmetry_target_box_ids(p,isym)-1_int64,int(nlocal,int64)))+1
          mapped_basis(p,:)=weights(p)*target_basis(:,target_point)
        enddo
      enddo
      call zgemm('C','T',nbasis,nbasis,nlocal,(1d0,0d0),mapped_basis,nlocal,local_basis,nbasis,&
        (0d0,0d0),local_overlap,nbasis)
      call MPI_Allreduce(local_overlap,symmetry_overlap(:,:,isym),nbasis*nbasis,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    enddo
    ok=.true.
#else
    ok=.false.;message='distributed full-basis symmetry overlap requires MPI'
#endif
  end subroutine

  subroutine assemble_dg_distributed_basis_symmetry_overlap_rows(comm,local_basis,weights,&
      symmetry_target_box_ids,row_ids,symmetry_overlap_rows,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::local_basis(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    integer(int64),allocatable,intent(out)::row_ids(:)
    complex(real64),allocatable,intent(out)::symmetry_overlap_rows(:,:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::orbital_tile_size=64
    integer::rank,nproc,ierr,nbasis,nlocal,nsym,isym,tile_first,tile_count,&
      owner_first,owner_count,base,remainder,i,local_bad,global_bad,status
    integer(int64)::persistent_bytes,tile_bytes,complex_bytes
    integer,allocatable::receive_counts(:)
    complex(real64),allocatable::image_tile(:,:),packed_tile(:,:),reduced_tile(:,:)

    ok=.false.;message='';workspace_peak_bytes=0_int64
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nbasis=size(local_basis,1);nlocal=size(local_basis,2);nsym=size(symmetry_target_box_ids,2)
    local_bad=merge(0,1,nbasis>0.and.nlocal>0.and.nsym>0.and.size(weights)==nlocal.and.&
      size(symmetry_target_box_ids,1)==nlocal.and.all(weights>=0d0).and.&
      all(ieee_is_finite(weights)).and.all(ieee_is_finite(real(local_basis))).and.&
      all(ieee_is_finite(aimag(local_basis))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid row-owned symmetry-overlap contract';return
    endif
    base=nbasis/nproc;remainder=mod(nbasis,nproc)
    owner_count=base+merge(1,0,rank<remainder)
    owner_first=rank*base+min(rank,remainder)+1
    allocate(row_ids(owner_count),symmetry_overlap_rows(owner_count,nbasis,nsym),&
      receive_counts(nproc),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(row_ids))deallocate(row_ids)
      if(allocated(symmetry_overlap_rows))deallocate(symmetry_overlap_rows)
      if(allocated(receive_counts))deallocate(receive_counts)
      message='row-owned symmetry-overlap persistent allocation failed';return
    endif
    symmetry_overlap_rows=(0d0,0d0)
    do i=1,owner_count;row_ids(i)=int(owner_first+i-1,int64);enddo
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    persistent_bytes=complex_bytes*int(size(symmetry_overlap_rows),int64)
    tile_count=min(orbital_tile_size,nbasis)
    do i=1,nproc
      receive_counts(i)=(base+merge(1,0,i-1<remainder))*tile_count
    enddo
    allocate(image_tile(tile_count,nlocal),packed_tile(tile_count,nbasis),&
      reduced_tile(tile_count,owner_count),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='row-owned symmetry-overlap tile allocation failed';return
    endif
    tile_bytes=complex_bytes*int(size(image_tile)+size(packed_tile)+size(reduced_tile),int64)
    workspace_peak_bytes=persistent_bytes+tile_bytes
    do isym=1,nsym
      do tile_first=1,nbasis,orbital_tile_size
        tile_count=min(orbital_tile_size,nbasis-tile_first+1)
        call exchange_dg_point_permuted_orbital_rows(comm,&
          local_basis(tile_first:tile_first+tile_count-1,:),&
          symmetry_target_box_ids(:,isym),image_tile(1:tile_count,:),ok,message)
        if(.not.ok)return
        do i=1,nlocal
          image_tile(1:tile_count,i)=weights(i)*image_tile(1:tile_count,i)
        enddo
        call zgemm('N','C',tile_count,nbasis,nlocal,(1d0,0d0),image_tile,size(image_tile,1),&
          local_basis,nbasis,(0d0,0d0),packed_tile,size(packed_tile,1))
        if(tile_count<size(packed_tile,1))packed_tile(tile_count+1:,:)=(0d0,0d0)
        call MPI_Reduce_scatter(packed_tile,reduced_tile,receive_counts,&
          MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then
          ok=.false.;message='row-owned symmetry-overlap reduce-scatter failed';return
        endif
        symmetry_overlap_rows(:,tile_first:tile_first+tile_count-1,isym)=&
          transpose(reduced_tile(1:tile_count,:))
      enddo
    enddo
    ok=all(ieee_is_finite(real(symmetry_overlap_rows))).and.&
      all(ieee_is_finite(aimag(symmetry_overlap_rows)))
    if(ok)then;message='';else;message='row-owned symmetry overlap is nonfinite';endif
#else
    ok=.false.;message='row-owned symmetry overlap requires MPI';workspace_peak_bytes=0_int64
#endif
  end subroutine assemble_dg_distributed_basis_symmetry_overlap_rows

  subroutine gather_dg_single_symmetry_representation(comm,row_ids,row_values,operation,writer_rank,&
      representation,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,operation,writer_rank
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::row_values(:,:,:)
    complex(real64),allocatable,intent(out)::representation(:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::send_buffer(:),receive_buffer(:)
    integer,allocatable::receive_counts(:),receive_displacements(:)
    integer::rank,nproc,ierr,nstate,nstate_min,nstate_max,noperation,noperation_min,noperation_max,&
      operation_min,operation_max,writer_min,writer_max,local_rows,expected_rows,first_row
    integer::i,j,position,local_valid,global_valid,gather_error
    integer(int64)::element_count,local_workspace,global_workspace,complex_bytes,integer_bytes

    workspace_peak_bytes=0_int64;ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nstate=size(row_values,2);noperation=size(row_values,3);local_rows=size(row_values,1)
    local_valid=1
    if(writer_rank<0.or.writer_rank>=nproc.or.size(row_ids)/=local_rows.or.nstate<=0.or.&
        operation<1.or.operation>noperation) local_valid=0
    if(local_valid==1)then
      expected_rows=nstate/nproc+merge(1,0,rank<mod(nstate,nproc))
      first_row=rank*(nstate/nproc)+min(rank,mod(nstate,nproc))+1
      if(local_rows/=expected_rows) local_valid=0
      if(local_rows>0)then
        if(any(row_ids/=[(int(first_row+i-1,int64),i=1,local_rows)])) local_valid=0
      end if
      if(operation>=1.and.operation<=noperation)then
        if(.not.all(ieee_is_finite(real(row_values(:,:,operation)))).or.&
            .not.all(ieee_is_finite(aimag(row_values(:,:,operation))))) local_valid=0
      end if
    end if
    call MPI_Allreduce(nstate,nstate_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nstate,nstate_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(noperation,noperation_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(noperation,noperation_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(operation,operation_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(operation,operation_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(writer_rank,writer_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(writer_rank,writer_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nstate_min/=nstate_max.or.noperation_min/=noperation_max.or.&
        operation_min/=operation_max.or.writer_min/=writer_max) local_valid=0
    element_count=int(nstate,int64)*int(nstate,int64)
    if(element_count>int(huge(0),int64)) local_valid=0
    call MPI_Allreduce(local_valid,global_valid,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_valid/=1)then
      allocate(representation(0,0));message='invalid one-operation symmetry gather contract';return
    end if

    allocate(receive_counts(nproc),receive_displacements(nproc))
    do i=0,nproc-1
      expected_rows=nstate/nproc+merge(1,0,i<mod(nstate,nproc))
      first_row=i*(nstate/nproc)+min(i,mod(nstate,nproc))
      receive_counts(i+1)=expected_rows*nstate
      receive_displacements(i+1)=first_row*nstate
    end do
    allocate(send_buffer(local_rows*nstate));position=0
    do i=1,local_rows;do j=1,nstate
      position=position+1;send_buffer(position)=row_values(i,j,operation)
    end do;end do
    if(rank==writer_rank)then
      allocate(receive_buffer(nstate*nstate),representation(nstate,nstate))
    else
      allocate(receive_buffer(0),representation(0,0))
    end if
    call MPI_Gatherv(send_buffer,size(send_buffer),MPI_DOUBLE_COMPLEX,receive_buffer,receive_counts,&
      receive_displacements,MPI_DOUBLE_COMPLEX,writer_rank,comm,ierr)
    gather_error=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(MPI_IN_PLACE,gather_error,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gather_error/=0)then
      ok=.false.;message='one-operation symmetry gather failed';return
    end if
    if(rank==writer_rank)then
      position=0
      do i=1,nstate;do j=1,nstate
        position=position+1;representation(i,j)=receive_buffer(position)
      end do;end do
    end if
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    integer_bytes=int(storage_size(0)/8,int64)
    local_workspace=complex_bytes*int(size(send_buffer)+size(receive_buffer)+size(representation),int64)+&
      integer_bytes*int(size(receive_counts)+size(receive_displacements),int64)
    call MPI_Allreduce(local_workspace,global_workspace,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS
    if(ok)then
      workspace_peak_bytes=global_workspace;message=''
    else
      workspace_peak_bytes=0_int64;message='one-operation symmetry workspace reduction failed'
    end if
#else
    allocate(representation(0,0));workspace_peak_bytes=0_int64
    ok=.false.;message='one-operation symmetry gather requires MPI'
#endif
  end subroutine gather_dg_single_symmetry_representation

  subroutine validate_dg_row_owned_group_representation(comm,row_ids,representation_rows,&
      product_table,identity_operation,tolerance,identity_defect,unitarity_defect,closure_defect,&
      workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,product_table(:,:),identity_operation
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::representation_rows(:,:,:)
    real(real64),intent(in)::tolerance
    real(real64),intent(out)::identity_defect,unitarity_defect,closure_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,nstate,nlocal,nsym,base,remainder,owner,owner_first,owner_count,&
      operation,left,right,product,i,j,local_bad,global_bad
    complex(real64),allocatable::remote_rows(:,:),unitarity_tile(:,:),product_rows(:,:)
    real(real64)::local_identity,local_unitarity,local_closure,expected
    integer(int64)::complex_bytes
    ok=.false.;message='';identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);workspace_peak_bytes=0_int64
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nlocal=size(row_ids);nstate=size(representation_rows,2);nsym=size(representation_rows,3)
    local_bad=0
    if(ierr/=MPI_SUCCESS.or.nstate<=0.or.nlocal/=size(representation_rows,1).or.nsym<=0.or.&
        any(shape(product_table)/=[nsym,nsym]).or.identity_operation<1.or.identity_operation>nsym.or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.any(product_table<1).or.&
        any(product_table>nsym).or..not.all(ieee_is_finite(real(representation_rows))).or.&
        .not.all(ieee_is_finite(aimag(representation_rows))))local_bad=1
    base=nstate/nproc;remainder=mod(nstate,nproc)
    owner_count=base+merge(1,0,rank<remainder);owner_first=rank*base+min(rank,remainder)+1
    if(nlocal/=owner_count)local_bad=1
    do i=1,nlocal
      if(row_ids(i)/=int(owner_first+i-1,int64))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='invalid row-owned group contract';return;endif
    local_identity=0d0
    do i=1,nlocal;do j=1,nstate
      expected=merge(1d0,0d0,int(row_ids(i))==j)
      local_identity=max(local_identity,abs(representation_rows(i,j,identity_operation)-expected))
    enddo;enddo
    local_unitarity=0d0;local_closure=0d0;complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    allocate(product_rows(nlocal,nstate));workspace_peak_bytes=complex_bytes*int(size(product_rows),int64)
    do operation=1,nsym
      do owner=0,nproc-1
        owner_count=base+merge(1,0,owner<remainder);owner_first=owner*base+min(owner,remainder)+1
        allocate(remote_rows(owner_count,nstate),unitarity_tile(nlocal,owner_count))
        if(rank==owner)remote_rows=representation_rows(:,:,operation)
        call MPI_Bcast(remote_rows,size(remote_rows),MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        unitarity_tile=matmul(representation_rows(:,:,operation),conjg(transpose(remote_rows)))
        do j=1,owner_count;do i=1,nlocal
          expected=merge(1d0,0d0,int(row_ids(i))==owner_first+j-1)
          local_unitarity=max(local_unitarity,abs(unitarity_tile(i,j)-expected))
        enddo;enddo
        workspace_peak_bytes=max(workspace_peak_bytes,complex_bytes*&
          int(size(product_rows)+size(remote_rows)+size(unitarity_tile),int64))
        deallocate(remote_rows,unitarity_tile)
      enddo
    enddo
    do left=1,nsym;do right=1,nsym
      product=product_table(left,right);product_rows=(0d0,0d0)
      do owner=0,nproc-1
        owner_count=base+merge(1,0,owner<remainder);owner_first=owner*base+min(owner,remainder)+1
        allocate(remote_rows(owner_count,nstate))
        if(rank==owner)remote_rows=representation_rows(:,:,right)
        call MPI_Bcast(remote_rows,size(remote_rows),MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        product_rows=product_rows+matmul(representation_rows(:,owner_first:owner_first+owner_count-1,left),&
          remote_rows)
        workspace_peak_bytes=max(workspace_peak_bytes,complex_bytes*&
          int(size(product_rows)+size(remote_rows),int64));deallocate(remote_rows)
      enddo
      local_closure=max(local_closure,maxval(abs(product_rows-representation_rows(:,:,product))))
    enddo;enddo
    call MPI_Allreduce(local_identity,identity_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_unitarity,unitarity_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_closure,closure_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.max(identity_defect,max(unitarity_defect,closure_defect))<=tolerance
    if(ok)then;message='';else;message='row-owned group representation violates identity, unitarity, or closure';endif
#else
    ok=.false.;message='row-owned group validation requires MPI';identity_defect=huge(1d0)
    unitarity_defect=huge(1d0);closure_defect=huge(1d0);workspace_peak_bytes=0_int64
#endif
  end subroutine validate_dg_row_owned_group_representation

  subroutine validate_dg_streamed_affine_representation(comm,local_basis,weights,&
      symmetry_target_box_ids,identity_operation,subspace_defect,tolerance,identity_defect,&
      unitarity_defect,closure_defect,workspace_peak_bytes,ok,message,prepared_row_ids,prepared_rows)
    integer,intent(in)::comm,identity_operation
    complex(real64),intent(in)::local_basis(:,:)
    real(real64),intent(in)::weights(:),subspace_defect,tolerance
    integer(int64),intent(in)::symmetry_target_box_ids(:,:)
    real(real64),intent(out)::identity_defect,unitarity_defect,closure_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),allocatable,intent(out),optional::prepared_row_ids(:)
    complex(real64),allocatable,intent(out),optional::prepared_rows(:,:,:)
#ifdef USE_MPI
    integer::rank,nproc,ierr,nstate,nlocal,nsym,operation,owner,base,remainder,&
      owner_first,owner_count,i,j,local_bad,global_bad,status
    integer(int64),allocatable::row_ids(:)
    complex(real64),allocatable::rows(:,:,:),remote_rows(:,:),unitarity_tile(:,:)
    real(real64)::local_identity,local_unitarity,expected
    integer(int64)::operation_peak,complex_bytes,persistent_bytes
    ok=.false.;message='';identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);workspace_peak_bytes=0_int64
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nstate=size(local_basis,1);nlocal=size(local_basis,2);nsym=size(symmetry_target_box_ids,2)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.nstate>0.and.nlocal>0.and.nsym>0.and.&
      identity_operation>=1.and.identity_operation<=nsym.and.tolerance>0d0.and.&
      subspace_defect>=0d0.and.ieee_is_finite(tolerance).and.ieee_is_finite(subspace_defect))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='invalid streamed affine proof contract';return;endif
    local_identity=0d0;local_unitarity=0d0
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    base=nstate/nproc;remainder=mod(nstate,nproc)
    owner_count=base+merge(1,0,rank<remainder)
    if(present(prepared_row_ids).neqv.present(prepared_rows))then
      message='prepared affine representation outputs must be paired';return
    endif
    if(present(prepared_rows))then
      allocate(prepared_row_ids(owner_count),prepared_rows(owner_count,nstate,nsym),stat=status)
      call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        if(allocated(prepared_row_ids))deallocate(prepared_row_ids)
        if(allocated(prepared_rows))deallocate(prepared_rows)
        message='prepared affine representation allocation failed';return
      endif
    endif
    do operation=1,nsym
      call assemble_dg_distributed_basis_symmetry_overlap_rows(comm,local_basis,weights,&
        symmetry_target_box_ids(:,operation:operation),row_ids,rows,operation_peak,ok,message)
      if(.not.ok)return
      persistent_bytes=complex_bytes*int(size(rows),int64)
      workspace_peak_bytes=max(workspace_peak_bytes,operation_peak)
      if(present(prepared_rows))then
        if(operation==1)prepared_row_ids=row_ids
        if(any(prepared_row_ids/=row_ids))then
          message='prepared affine representation ownership changed';return
        endif
        prepared_rows(:,:,operation)=rows(:,:,1)
        workspace_peak_bytes=max(workspace_peak_bytes,operation_peak+complex_bytes*int(size(prepared_rows),int64))
      endif
      if(operation==identity_operation)then
        do i=1,size(row_ids);do j=1,nstate
          expected=merge(1d0,0d0,int(row_ids(i))==j)
          local_identity=max(local_identity,abs(rows(i,j,1)-expected))
        enddo;enddo
      endif
      do owner=0,nproc-1
        owner_count=base+merge(1,0,owner<remainder)
        owner_first=owner*base+min(owner,remainder)+1
        allocate(remote_rows(owner_count,nstate),unitarity_tile(size(row_ids),owner_count))
        if(rank==owner)remote_rows=rows(:,:,1)
        call MPI_Bcast(remote_rows,size(remote_rows),MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='streamed affine row broadcast failed';return;endif
        unitarity_tile=matmul(rows(:,:,1),conjg(transpose(remote_rows)))
        do j=1,owner_count;do i=1,size(row_ids)
          expected=merge(1d0,0d0,int(row_ids(i))==owner_first+j-1)
          local_unitarity=max(local_unitarity,abs(unitarity_tile(i,j)-expected))
        enddo;enddo
        workspace_peak_bytes=max(workspace_peak_bytes,persistent_bytes+complex_bytes*&
          int(size(remote_rows)+size(unitarity_tile),int64))
        deallocate(remote_rows,unitarity_tile)
      enddo
      deallocate(row_ids,rows)
    enddo
    call MPI_Allreduce(local_identity,identity_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_unitarity,unitarity_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    closure_defect=max(subspace_defect,2d0*subspace_defect+unitarity_defect)
    ok=ierr==MPI_SUCCESS.and.max(identity_defect,max(unitarity_defect,closure_defect))<=tolerance
    if(ok)then;message='';else;message='streamed affine representation violates proof tolerance';endif
#else
    ok=.false.;message='streamed affine proof requires MPI';identity_defect=huge(1d0)
    unitarity_defect=huge(1d0);closure_defect=huge(1d0);workspace_peak_bytes=0_int64
#endif
  end subroutine validate_dg_streamed_affine_representation

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
      periodic_localization_phase,candidate_axis_offset,center_representative_box_ids,&
      projection_seed_values)
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
    integer,intent(in),optional::candidate_axis_offset
    integer(int64),intent(in),optional::center_representative_box_ids(:)
    real(real64),intent(in),optional::projection_seed_values(:,:)
#ifdef USE_MPI
    complex(real64),allocatable::s_local(:,:),s(:,:),l_local(:,:),localizer(:,:),x(:,:),&
      occ_gram(:,:),occ_vectors(:,:),a_occ(:,:),overlap_occ_x(:,:),residual(:,:),&
      residual_gram(:,:),residual_vectors(:,:),q_comp(:,:),comp_localizer(:,:),&
      comp_vectors(:,:),transform(:,:),final_metric(:,:),metric_vectors(:,:),&
      occupied_reference(:,:),symmetry_product(:,:),orthogonal_transform(:,:),mapped_transform(:,:),&
      candidate_symmetry(:,:,:),all_candidate_flat(:),all_candidate(:,:),spatial_overlap(:,:),&
      link_local(:,:,:),link_global(:,:,:),symmetrized_localizer(:,:),all_phase_flat(:),&
      canonical_phase(:,:),all_wannier(:,:),mapped_candidate(:,:),distributed_spatial_overlap(:,:,:),&
      seed_overlap_local(:,:),seed_overlap(:,:),projected_seed(:,:)
    complex(real64),allocatable::seed_gram(:,:)
    complex(real64),allocatable::polar_vectors(:,:),polar_inverse(:,:)
    complex(real64),allocatable::symmetry_block(:,:)
    complex(real64),allocatable::block_vectors(:,:)
    complex(real64),allocatable::retained_projector(:,:),retained_projector_vectors(:,:)
    real(real64),allocatable::spectrum(:),occ_spectrum(:),residual_spectrum(:),comp_spectrum(:),&
      metric_spectrum(:),center_max_local(:),center_max_global(:),polar_spectrum(:)
    real(real64),allocatable::block_spectrum(:)
    real(real64),allocatable::retained_projector_spectrum(:)
    logical,allocatable::integration_core(:),all_core_unsorted(:),canonical_core(:)
    real(real64),allocatable::all_box_weights_unsorted(:),all_box_weights(:)
    integer(int64),allocatable::all_box_ids(:),all_target_ids(:,:),canonical_target_ids(:,:),&
      sorted_target_ids(:)
    integer(int64),allocatable::center_id_local(:),center_id_global(:)
    integer,allocatable::box_counts(:),box_displs(:),value_counts(:),value_displs(:),all_fragments(:),&
      canonical_fragments(:),canonical_ranks(:),&
      phase_counts(:),phase_displs(:)
    integer::nlocal,i,j,p,ierr,nproc,rank,local_bad,global_bad,&
      ncomp,nneed,nseed,selected_target,start_index,matrix_count,&
      scalar_min(4),scalar_max(4),scalar_local(4),&
      symmetry_present,symmetry_present_min,symmetry_present_max,nsym,isym,jsym,ksym,&
      nsym_min,nsym_max,core_mask_present,core_mask_present_min,core_mask_present_max,&
      phase_present,phase_present_min,phase_present_max,fragment_min,fragment_max
    integer::total_box,source,target,axis,source_axis,conjugate_flag,center_index,&
      local_candidate_count,candidate_begin,candidate_end
    integer::source_rank,middle_rank,final_rank,expected_rank,source_begin,middle_begin,final_begin
    real(real64)::largest,local_boundary_value,local_boundary_gradient
    real(real64)::tolerance_local(3),tolerance_min(3),tolerance_max(3)
    real(real64)::active_symmetry_tolerance,symmetry_tolerance_min,symmetry_tolerance_max,&
      symmetry_defect,product_residual,retained_projector_gap
    real(real64)::raw_symmetry_defect,polar_correction,best_product_residual,group_closure_defect
    complex(real64)::phase_ratio,trial_ratio
    real(real64),parameter::periodic_real_weight(3)=[sqrt(2d0),sqrt(3d0),sqrt(5d0)]
    real(real64),parameter::periodic_imag_weight(3)=[sqrt(7d0),sqrt(11d0),sqrt(13d0)]
    integer(int64)::local_fingerprint,global_fingerprint,quantized_density,matrix_count64,&
      expected_min,expected_max,total_box64,value_total64,running64
    integer(int64)::expected_box_min,expected_box_max
    integer(int64)::local_core_count,global_core_count
    logical::matched_product,phase_covariant,source_axis_used(3),distributed_candidates

    ok=.false.;message='';nlocal=size(physical_ids);nseed=0
    if(present(projection_seed_values))nseed=size(projection_seed_values,1)
    call MPI_Comm_size(comm,nproc,ierr)
    local_candidate_count=size(candidate_value,1)
    distributed_candidates=present(candidate_axis_offset)
    candidate_begin=1
    if(distributed_candidates)candidate_begin=candidate_axis_offset+1
    candidate_end=candidate_begin+local_candidate_count-1
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
    if(present(center_representative_box_ids))then
      if(size(center_representative_box_ids)/=nlocal)local_bad=1
      if(any(center_representative_box_ids<1_int64).or.&
          any(center_representative_box_ids>int(nlocal,int64)))local_bad=1
    endif
    if((.not.distributed_candidates.and.local_candidate_count/=ncandidate).or.&
        size(candidate_value,2)/=nlocal)local_bad=1
    if(distributed_candidates.and.(candidate_begin<1.or.candidate_end>ncandidate))local_bad=1
    if(size(candidate_gradient,1)/=3.or.size(candidate_gradient,2)/=local_candidate_count.or.&
        size(candidate_gradient,3)/=nlocal)local_bad=1
    if(size(occupied_coefficients,1)/=ncandidate.or.size(occupied_coefficients,2)/=noccupied)local_bad=1
    if(present(projection_seed_values))then
      if(nseed<1.or.nseed>ncandidate.or.size(projection_seed_values,2)/=nlocal)local_bad=1
      if(.not.all(ieee_is_finite(projection_seed_values)))local_bad=1
    endif
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
      allocate(all_target_ids(total_box,nsym),all_phase_flat(3*total_box))
      if(.not.distributed_candidates)allocate(all_candidate_flat(int(value_total64)))
      call MPI_Allgatherv(box_point_ids,nlocal,MPI_INTEGER8,all_box_ids,box_counts,box_displs,MPI_INTEGER8,comm,ierr)
      call MPI_Allgatherv(weights,nlocal,MPI_DOUBLE_PRECISION,all_box_weights_unsorted,box_counts,&
        box_displs,MPI_DOUBLE_PRECISION,comm,ierr)
      call MPI_Allgatherv(core_fragment,nlocal,MPI_INTEGER,all_fragments,box_counts,box_displs,MPI_INTEGER,comm,ierr)
      call MPI_Allgatherv(integration_core,nlocal,MPI_LOGICAL,all_core_unsorted,box_counts,box_displs,&
        MPI_LOGICAL,comm,ierr)
      if(.not.distributed_candidates)then
        call MPI_Allgatherv(candidate_value,ncandidate*nlocal,MPI_DOUBLE_COMPLEX,all_candidate_flat,&
          value_counts,value_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
      endif
      call MPI_Allgatherv(periodic_localization_phase,3*nlocal,MPI_DOUBLE_COMPLEX,all_phase_flat,&
        phase_counts,phase_displs,MPI_DOUBLE_COMPLEX,comm,ierr)
      do isym=1,nsym
        call MPI_Allgatherv(symmetry_target_box_ids(:,isym),nlocal,MPI_INTEGER8,all_target_ids(:,isym),&
          box_counts,box_displs,MPI_INTEGER8,comm,ierr)
      enddo
      if(int(total_box,int64)/=expected_box_count)then
        message='periodic-box construction call must contain complete fragment boxes';return
      endif
      allocate(all_box_weights(total_box),canonical_core(total_box),&
        canonical_phase(3,total_box),canonical_fragments(total_box),canonical_ranks(total_box))
      if(.not.distributed_candidates)allocate(all_candidate(ncandidate,total_box))
      allocate(canonical_target_ids(total_box,nsym))
      do source=1,total_box
        if(all_box_ids(source)<1_int64.or.all_box_ids(source)>int(total_box,int64))then
          message='invalid periodic-box point id';return
        endif
        target=int(all_box_ids(source))
        if(.not.distributed_candidates)&
          all_candidate(:,target)=all_candidate_flat((source-1)*ncandidate+1:source*ncandidate)
        all_box_weights(target)=all_box_weights_unsorted(source)
        canonical_core(target)=all_core_unsorted(source)
        canonical_fragments(target)=all_fragments(source)
        canonical_ranks(target)=count(box_displs<=source-1)-1
        canonical_phase(:,target)=all_phase_flat((source-1)*3+1:source*3)
        canonical_target_ids(target,:)=all_target_ids(source,:)
      enddo
      if(allocated(all_candidate_flat))deallocate(all_candidate_flat)
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
      do j=1,local_candidate_count;do i=1,local_candidate_count
        s_local(candidate_begin+i-1,candidate_begin+j-1)=&
          s_local(candidate_begin+i-1,candidate_begin+j-1)+&
          weights(p)*conjg(candidate_value(i,p))*candidate_value(j,p)
        if(nsym>0)then
          do axis=1,3
            link_local(axis,candidate_begin+i-1,candidate_begin+j-1)=&
              link_local(axis,candidate_begin+i-1,candidate_begin+j-1)+&
              weights(p)*conjg(candidate_value(i,p))*periodic_localization_phase(axis,p)*&
              candidate_value(j,p)
          enddo
        else
          l_local(candidate_begin+i-1,candidate_begin+j-1)=&
            l_local(candidate_begin+i-1,candidate_begin+j-1)+weights(p)*localization_coordinate(p)*&
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
    if(distributed_candidates)then
      allocate(spectrum(ncandidate),x(ncandidate,ncandidate));x=(0d0,0d0)
      do source_rank=0,nproc-1
        source_begin=source_rank*local_candidate_count+1
        call hermitian_eigensystem(s(source_begin:source_begin+local_candidate_count-1,&
          source_begin:source_begin+local_candidate_count-1),block_spectrum,block_vectors,ok,message)
        if(.not.ok)return
        spectrum(source_begin:source_begin+local_candidate_count-1)=block_spectrum
        x(source_begin:source_begin+local_candidate_count-1,&
          source_begin:source_begin+local_candidate_count-1)=block_vectors
        deallocate(block_spectrum,block_vectors)
      enddo
    else
      call hermitian_eigensystem(s,spectrum,x,ok,message)
      if(.not.ok)return
    endif
    largest=maxval(spectrum)
    if(largest<=0d0.or.count(spectrum>rank_tolerance*largest)<ncandidate)then
      ok=.false.;message='candidate rank loss in overlapping-Wannier construction';return
    endif
    do j=1,ncandidate
      x(:,j)=x(:,j)/sqrt(spectrum(j))
    enddo
    if(nsym>0)then
      allocate(candidate_symmetry(ncandidate,ncandidate,nsym),spatial_overlap(ncandidate,ncandidate))
      if(distributed_candidates)then
        call assemble_dg_distributed_candidate_symmetry(comm,candidate_value,weights,&
          symmetry_target_box_ids,distributed_spatial_overlap,ok,message)
        if(.not.ok)return
      else
        allocate(mapped_candidate(total_box,ncandidate))
      endif
      do isym=1,nsym
        if(distributed_candidates)then
          spatial_overlap=distributed_spatial_overlap(:,:,isym)
        else
        do source=1,total_box
          target=int(canonical_target_ids(source,isym))
          if(abs(all_box_weights(target)-all_box_weights(source))>&
              active_symmetry_tolerance*max(1d0,all_box_weights(source)))then
            ok=.false.;message='periodic-box symmetry does not preserve quadrature weights';return
          endif
          mapped_candidate(source,:)=all_box_weights(target)*all_candidate(:,target)
        enddo
        call zgemm('C','T',ncandidate,ncandidate,total_box,(1d0,0d0),mapped_candidate,total_box,&
          all_candidate,ncandidate,(0d0,0d0),spatial_overlap,ncandidate)
        endif
        candidate_symmetry(:,:,isym)=matmul(conjg(transpose(x)),matmul(spatial_overlap,x))
      enddo
      if(allocated(mapped_candidate))deallocate(mapped_candidate)
      if(allocated(distributed_spatial_overlap))deallocate(distributed_spatial_overlap)
      if(.not.all(ieee_is_finite(real(candidate_symmetry))).or.&
          .not.all(ieee_is_finite(aimag(candidate_symmetry))))then
        ok=.false.;message='nonfinite periodic-box candidate symmetry representation';return
      endif
      allocate(symmetry_product(ncandidate,ncandidate))
      if(distributed_candidates)then
        allocate(polar_vectors(local_candidate_count,local_candidate_count),&
          polar_inverse(local_candidate_count,local_candidate_count),&
          symmetry_block(local_candidate_count,local_candidate_count))
      else
        allocate(polar_vectors(ncandidate,ncandidate),polar_inverse(ncandidate,ncandidate))
      endif
      raw_symmetry_defect=0d0
      do isym=1,nsym
        raw_symmetry_defect=max(raw_symmetry_defect,&
          maxval(abs(matmul(conjg(transpose(candidate_symmetry(:,:,isym))),&
          candidate_symmetry(:,:,isym))-identity_complex(ncandidate))))
      enddo
      call MPI_Comm_rank(comm,rank,ierr)
      if(rank==0)write(*,'(a,es24.16)')&
        '[OW-GS-DIAGNOSTIC] candidate_symmetry_raw_unitarity_defect=',raw_symmetry_defect
      if(raw_symmetry_defect>max(100d0*active_symmetry_tolerance,&
          100d0*epsilon(1d0)*real(ncandidate,real64)))then
        ok=.false.;message='periodic-box candidate symmetry representation is not unitary';return
      endif
      polar_correction=0d0
      do isym=1,nsym
        if(distributed_candidates)then
          do source_rank=0,nproc-1
            source=source_rank*nlocal+1
            final_rank=int((canonical_target_ids(source,isym)-1_int64)/int(nlocal,int64))
            source_begin=source_rank*local_candidate_count+1
            final_begin=final_rank*local_candidate_count+1
            symmetry_block=candidate_symmetry(final_begin:final_begin+local_candidate_count-1,&
              source_begin:source_begin+local_candidate_count-1,isym)
            polar_inverse=matmul(conjg(transpose(symmetry_block)),symmetry_block)
            polar_inverse=0.5d0*(polar_inverse+conjg(transpose(polar_inverse)))
            call hermitian_eigensystem(polar_inverse,polar_spectrum,polar_vectors,ok,message)
            if(.not.ok)return
            if(minval(polar_spectrum)<=0d0)then
              ok=.false.;message='singular periodic-box candidate symmetry polar factor';return
            endif
            polar_inverse=polar_vectors
            do j=1,local_candidate_count
              polar_inverse(:,j)=polar_inverse(:,j)/sqrt(polar_spectrum(j))
            enddo
            polar_inverse=matmul(polar_inverse,conjg(transpose(polar_vectors)))
            candidate_symmetry(final_begin:final_begin+local_candidate_count-1,&
              source_begin:source_begin+local_candidate_count-1,isym)=matmul(symmetry_block,polar_inverse)
            polar_correction=max(polar_correction,maxval(abs(&
              candidate_symmetry(final_begin:final_begin+local_candidate_count-1,&
                source_begin:source_begin+local_candidate_count-1,isym)-symmetry_block)))
          enddo
        else
          symmetry_product=matmul(conjg(transpose(candidate_symmetry(:,:,isym))),&
            candidate_symmetry(:,:,isym))
          symmetry_product=0.5d0*(symmetry_product+conjg(transpose(symmetry_product)))
          call hermitian_eigensystem(symmetry_product,polar_spectrum,polar_vectors,ok,message)
          if(.not.ok)return
          if(minval(polar_spectrum)<=0d0)then
            ok=.false.;message='singular periodic-box candidate symmetry polar factor';return
          endif
          polar_inverse=polar_vectors
          do j=1,ncandidate
            polar_inverse(:,j)=polar_inverse(:,j)/sqrt(polar_spectrum(j))
          enddo
          polar_inverse=matmul(polar_inverse,conjg(transpose(polar_vectors)))
          spatial_overlap=candidate_symmetry(:,:,isym)
          candidate_symmetry(:,:,isym)=matmul(candidate_symmetry(:,:,isym),polar_inverse)
          polar_correction=max(polar_correction,maxval(abs(candidate_symmetry(:,:,isym)-spatial_overlap)))
        endif
      enddo
      symmetry_defect=0d0
      do isym=1,nsym
        symmetry_defect=max(symmetry_defect,maxval(abs(matmul(conjg(transpose(candidate_symmetry(:,:,isym))),&
          candidate_symmetry(:,:,isym))-identity_complex(ncandidate))))
      enddo
      if(rank==0)write(*,'(a,es24.16)')&
        '[OW-GS-DIAGNOSTIC] candidate_symmetry_polar_correction=',polar_correction
      if(symmetry_defect>active_symmetry_tolerance)then
        ok=.false.;message='polar-corrected periodic-box candidate symmetry is not unitary';return
      endif
      group_closure_defect=0d0
      do isym=1,nsym;do jsym=1,nsym
        best_product_residual=huge(1d0)
        do ksym=1,nsym
          product_residual=0d0
          if(distributed_candidates)then
            do source_rank=0,nproc-1
              source=source_rank*nlocal+1
              middle_rank=int((canonical_target_ids(source,jsym)-1_int64)/int(nlocal,int64))
              final_rank=int((canonical_target_ids(middle_rank*nlocal+1,isym)-1_int64)/int(nlocal,int64))
              expected_rank=int((canonical_target_ids(source,ksym)-1_int64)/int(nlocal,int64))
              if(final_rank/=expected_rank)then
                product_residual=huge(1d0);exit
              endif
              source_begin=source_rank*local_candidate_count+1
              middle_begin=middle_rank*local_candidate_count+1
              final_begin=final_rank*local_candidate_count+1
              symmetry_block=matmul(&
                candidate_symmetry(final_begin:final_begin+local_candidate_count-1,&
                  middle_begin:middle_begin+local_candidate_count-1,isym),&
                candidate_symmetry(middle_begin:middle_begin+local_candidate_count-1,&
                  source_begin:source_begin+local_candidate_count-1,jsym))
              product_residual=max(product_residual,maxval(abs(symmetry_block-&
                candidate_symmetry(final_begin:final_begin+local_candidate_count-1,&
                  source_begin:source_begin+local_candidate_count-1,ksym))))
            enddo
          else
            symmetry_product=matmul(candidate_symmetry(:,:,isym),candidate_symmetry(:,:,jsym))
            product_residual=maxval(abs(symmetry_product-candidate_symmetry(:,:,ksym)))
          endif
          do source=1,total_box
            target=int(canonical_target_ids(int(canonical_target_ids(source,jsym)),isym))
            if(target/=int(canonical_target_ids(source,ksym)))product_residual=huge(1d0)
          enddo
          best_product_residual=min(best_product_residual,product_residual)
        enddo
        group_closure_defect=max(group_closure_defect,best_product_residual)
      enddo;enddo
      if(rank==0)write(*,'(a,es24.16)')&
        '[OW-GS-DIAGNOSTIC] candidate_symmetry_group_closure_defect=',group_closure_defect
      if(group_closure_defect>active_symmetry_tolerance)then
        ok=.false.;message='periodic-box spatial and candidate symmetry representations are not homomorphic';return
      endif
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

    if(present(projection_seed_values))then
      allocate(seed_overlap_local(ncandidate,nseed),seed_overlap(ncandidate,nseed))
      seed_overlap_local=(0d0,0d0)
      do p=1,nlocal
        do j=1,nseed;do i=1,local_candidate_count
          seed_overlap_local(candidate_begin+i-1,j)=seed_overlap_local(candidate_begin+i-1,j)+&
            weights(p)*conjg(candidate_value(i,p))*projection_seed_values(j,p)
        enddo;enddo
      enddo
      call MPI_Allreduce(seed_overlap_local,seed_overlap,ncandidate*nseed,MPI_DOUBLE_COMPLEX,&
        MPI_SUM,comm,ierr)
      projected_seed=matmul(x,matmul(conjg(transpose(x)),seed_overlap))
      seed_gram=matmul(conjg(transpose(projected_seed)),matmul(s,projected_seed))
      seed_gram=0.5d0*(seed_gram+conjg(transpose(seed_gram)))
      largest=maxval(abs(seed_gram))
      if(largest<=tiny(1d0).or.&
          .not.positive_definite_above(seed_gram/largest,rank_tolerance))then
        ok=.false.;message='projected complete shell is rank deficient';return
      endif
      overlap_occ_x=matmul(conjg(transpose(a_occ)),matmul(s,projected_seed))
      residual=projected_seed-matmul(a_occ,overlap_occ_x)
    else
      overlap_occ_x=matmul(conjg(transpose(a_occ)),matmul(s,x))
      residual=x-matmul(a_occ,overlap_occ_x)
    endif
    residual_gram=matmul(conjg(transpose(residual)),matmul(s,residual))
    residual_gram=0.5d0*(residual_gram+conjg(transpose(residual_gram)))
    if(.not.finite_complex_matrix(residual_gram))then
      ok=.false.;message='nonfinite periodic localization complement Gram matrix';return
    endif
    if(present(projection_seed_values))then
      if(maxval(abs(residual_gram))<=&
          100d0*epsilon(1d0)*max(tiny(1d0),maxval(abs(seed_overlap))**2))then
        allocate(residual_spectrum(size(residual_gram,1)),&
          residual_vectors(size(residual_gram,1),size(residual_gram,1)))
        residual_spectrum=0d0;residual_vectors=identity_complex(size(residual_gram,1))
      else
        call hermitian_eigensystem(residual_gram,residual_spectrum,residual_vectors,ok,message)
        if(.not.ok)return
      endif
    else
      call hermitian_eigensystem(residual_gram,residual_spectrum,residual_vectors,ok,message)
      if(.not.ok)return
    endif
    if(present(projection_seed_values))then
      largest=maxval(residual_spectrum)
      if(largest<=0d0)then
        ncomp=0
      else
        ncomp=count(residual_spectrum>rank_tolerance*largest)
      endif
      selected_target=noccupied+ncomp
      nneed=ncomp
    else
      ncomp=count(residual_spectrum>rank_tolerance*max(1d0,maxval(residual_spectrum)))
      selected_target=ntarget
      nneed=selected_target-noccupied
    endif
    if(selected_target>ncandidate)then
      ok=.false.;message='occupied plus complete-shell direct sum exceeds candidate rank';return
    endif
    if(.not.present(projection_seed_values).and.ncomp<nneed)then
      ok=.false.;message='target rank loss in localization complement';return
    endif
    if(nneed>0)then
      if(present(projection_seed_values))then
        start_index=size(residual_spectrum)-nneed+1
      else
        start_index=size(residual_spectrum)-ncomp+1
      endif
      q_comp=matmul(residual,residual_vectors(:,start_index:size(residual_spectrum)))
      do j=1,size(q_comp,2)
        q_comp(:,j)=q_comp(:,j)/sqrt(residual_spectrum(start_index+j-1))
      enddo
      comp_localizer=matmul(conjg(transpose(q_comp)),matmul(localizer,q_comp))
      comp_localizer=0.5d0*(comp_localizer+conjg(transpose(comp_localizer)))
      if(.not.finite_complex_matrix(comp_localizer))then
        ok=.false.;message='nonfinite periodic localization complement matrix';return
      endif
      call hermitian_eigensystem(comp_localizer,comp_spectrum,comp_vectors,ok,message)
      if(.not.ok)return
      allocate(transform(ncandidate,selected_target))
      transform(:,1:noccupied)=a_occ
      transform(:,noccupied+1:selected_target)=matmul(q_comp,comp_vectors(:,1:nneed))
    else
      allocate(transform,source=a_occ)
    endif
    result%symmetry_closure_residual=0d0
    if(nsym>0)then
      orthogonal_transform=matmul(conjg(transpose(x)),matmul(s,transform))
      if(distributed_candidates)then
        retained_projector_gap=1d0
        if(selected_target<ncandidate)then
          allocate(retained_projector(ncandidate,ncandidate));retained_projector=(0d0,0d0)
          do isym=1,nsym
            mapped_transform=matmul(candidate_symmetry(:,:,isym),orthogonal_transform)
            retained_projector=retained_projector+&
              matmul(mapped_transform,conjg(transpose(mapped_transform)))
          enddo
          retained_projector=retained_projector/real(nsym,real64)
          retained_projector=0.5d0*(retained_projector+conjg(transpose(retained_projector)))
          call hermitian_eigensystem(retained_projector,retained_projector_spectrum,&
            retained_projector_vectors,ok,message)
          if(.not.ok)return
          retained_projector_gap=retained_projector_spectrum(ncandidate-selected_target+1)-&
            retained_projector_spectrum(ncandidate-selected_target)
        endif
        call MPI_Comm_rank(comm,rank,ierr)
        if(rank==0)write(*,'(a,es24.16)')&
          '[OW-GS-DIAGNOSTIC] retained_symmetry_projector_gap=',retained_projector_gap
        if(selected_target<ncandidate.and.&
            retained_projector_gap<=active_symmetry_tolerance)then
          ok=.false.;message='retained symmetry projector has no invariant-subspace gap';return
        endif
        if(selected_target<ncandidate)then
          transform=matmul(x,&
            retained_projector_vectors(:,ncandidate-selected_target+1:ncandidate))
          orthogonal_transform=matmul(conjg(transpose(x)),matmul(s,transform))
        endif
      endif
      allocate(result%symmetry_representation(selected_target,selected_target,nsym))
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

    result%candidate_rank=ncandidate;result%target_rank=selected_target
    result%retained_rank=selected_target;result%generation=generation
    allocate(result%physical_grid_ids,source=physical_ids)
    allocate(result%transform,source=transform)
    allocate(result%value(selected_target,nlocal),result%gradient(3,selected_target,nlocal))
    result%value=matmul(transpose(transform(candidate_begin:candidate_end,:)),candidate_value)
    do p=1,nlocal;do i=1,3
      result%gradient(i,:,p)=matmul(transpose(transform(candidate_begin:candidate_end,:)),&
        candidate_gradient(i,:,p))
    enddo;enddo
    if(.not.finite_complex_matrix(result%value).or.&
        .not.all(ieee_is_finite(real(result%gradient))).or.&
        .not.all(ieee_is_finite(aimag(result%gradient))))then
      ok=.false.;message='nonfinite periodic-box Wannier value or gradient tails';return
    endif
    block
      if(distributed_candidates)then
        allocate(result%center_box_point_ids(selected_target),center_max_local(selected_target),&
          center_max_global(selected_target),center_id_local(selected_target),&
          center_id_global(selected_target))
        center_max_local=0d0;center_id_local=huge(1_int64)
        do j=1,selected_target
          center_max_local(j)=maxval(abs(result%value(j,:))**2)
        enddo
        call MPI_Allreduce(center_max_local,center_max_global,selected_target,&
          MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
        do j=1,selected_target
          do p=1,nlocal
            if(center_max_global(j)-abs(result%value(j,p))**2<=&
                active_symmetry_tolerance*center_max_global(j))then
              center_id_local(j)=min(center_id_local(j),box_point_ids(p))
            endif
          enddo
        enddo
        call MPI_Allreduce(center_id_local,center_id_global,selected_target,&
          MPI_INTEGER8,MPI_MIN,comm,ierr)
        result%center_box_point_ids=center_id_global
        do j=1,selected_target
          if(center_max_global(j)<=0d0.or.center_id_global(j)==huge(1_int64).or.&
              .not.canonical_core(int(center_id_global(j))))then
            ok=.false.;message='distributed periodic-box Wannier center is not core owned';return
          endif
        enddo
      else
      if(nsym>0)then
        all_wannier=matmul(transpose(transform),all_candidate)
      else
        allocate(all_wannier(selected_target,nlocal));all_wannier=result%value
      endif
      if(.not.finite_complex_matrix(all_wannier))then
        ok=.false.;message='nonfinite periodic-box Wannier tails';return
      endif
      allocate(result%center_box_point_ids(selected_target))
      do j=1,selected_target
        if(nsym==0.and..not.present(center_representative_box_ids))then
          largest=0d0
          do p=1,nlocal
            if(integration_core(p))largest=max(largest,abs(all_wannier(j,p))**2)
          end do
        else
          largest=maxval(abs(all_wannier(j,:))**2)
        endif
        center_index=0
        if(largest<=0d0.or..not.ieee_is_finite(largest))then
          ok=.false.;message='periodic-box Wannier has no finite nonzero center density';return
        endif
        do source=1,size(all_wannier,2)
          if(nsym==0.and..not.present(center_representative_box_ids))then
            if(.not.integration_core(source))cycle
          endif
          if(largest-abs(all_wannier(j,source))**2<=&
              active_symmetry_tolerance*largest)then
            center_index=source;exit
          endif
        enddo
        if(present(center_representative_box_ids))then
          center_index=int(center_representative_box_ids(center_index))
        endif
        if(center_index==0)then
          ok=.false.;message='periodic-box Wannier center is not owned by the fragment core';return
        endif
        if(nsym>0)then
          if(.not.canonical_core(center_index))then
            ok=.false.;message='periodic-box Wannier center representative is not core owned';return
          endif
        else
          if(.not.integration_core(center_index))then
            ok=.false.;message='local periodic-box Wannier center representative is not core owned';return
          endif
        endif
        result%center_box_point_ids(j)=int(center_index,int64)
      enddo
      if(nsym>0)then
        call verify_dg_fragment_subspace_density_covariance(all_wannier,canonical_target_ids,&
          active_symmetry_tolerance,ok,message)
        if(.not.ok)return
      endif
      endif
    end block

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
    result%projection_inclusion_residual=0d0
    if(present(projection_seed_values))then
      residual_gram=metric_vectors
      do j=1,selected_target
        residual_gram(:,j)=residual_gram(:,j)/metric_spectrum(j)
      enddo
      residual_gram=matmul(residual_gram,conjg(transpose(metric_vectors)))
      seed_overlap=matmul(conjg(transpose(projected_seed)),matmul(s,projected_seed))
      overlap_occ_x=matmul(conjg(transpose(transform)),matmul(s,projected_seed))
      occ_gram=seed_overlap-matmul(conjg(transpose(overlap_occ_x)),&
        matmul(residual_gram,overlap_occ_x))
      result%projection_inclusion_residual=maxval(abs(occ_gram))/&
        max(tiny(1d0),maxval(abs(seed_overlap)))
      if(.not.ieee_is_finite(result%projection_inclusion_residual).or.&
          result%projection_inclusion_residual>rank_tolerance)then
        ok=.false.;message='complete projector shell inclusion tolerance exceeded';return
      endif
    endif
    local_boundary_value=0d0;local_boundary_gradient=0d0
    do p=1,nlocal
      if(.not.boundary_mask(p))cycle
      local_boundary_value=max(local_boundary_value,maxval(abs(result%value(:,p))))
      local_boundary_gradient=max(local_boundary_gradient,maxval(abs(result%gradient(:,:,p))))
    enddo
    call MPI_Allreduce(local_boundary_value,result%boundary_value_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_boundary_gradient,result%boundary_gradient_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Comm_rank(comm,rank,ierr)
    if(rank==0)then
      if(present(center_representative_box_ids))then
        write(*,'(a,es24.16,a,es24.16)')&
          '[OW-GS-DIAGNOSTIC] periodic_buffer_boundary_value_norm=',result%boundary_value_max,&
          ' periodic_buffer_boundary_gradient_norm=',result%boundary_gradient_max
      else
        write(*,'(a,es24.16,a,es24.16)')&
          '[OW-GS-DIAGNOSTIC] boundary_value_max=',result%boundary_value_max,&
          ' boundary_gradient_max=',result%boundary_gradient_max
      endif
    endif
    if(.not.present(center_representative_box_ids))then
      if(result%boundary_value_max>boundary_value_tolerance.or.&
          result%boundary_gradient_max>boundary_gradient_tolerance)then
        ok=.false.;message='buffer-boundary value or gradient tolerance exceeded';return
      endif
    endif

    call MPI_Comm_rank(comm,rank,ierr)
    allocate(result%center_owner_rank(selected_target),&
      result%center_owner_fragment(selected_target))
    if(nsym>0)then
      do j=1,selected_target
        result%center_owner_rank(j)=canonical_ranks(int(result%center_box_point_ids(j)))
        result%center_owner_fragment(j)=canonical_fragments(int(result%center_box_point_ids(j)))
      enddo
    else
      result%center_owner_rank=rank
      result%center_owner_fragment=fragment_min
    endif

    local_fingerprint=0_int64
    if(rank==0)local_fingerprint=ieor(int(generation,int64),&
      ishftc(int(selected_target,int64),11))
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
