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
  public::construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction
  public::verify_dg_overlapping_wannier_periodic_closure
  public::assemble_dg_distributed_candidate_symmetry
  public::assemble_dg_distributed_basis_symmetry_overlap
  public::assemble_dg_distributed_basis_symmetry_overlap_rows
  public::build_dg_pointwise_affine_owner_map
  public::select_dg_fixed_rank_symmetry_closed_subspace
  public::build_dg_distributed_symmetry_closed_basis
  public::align_dg_fragment_wannier_gauge
  public::replicate_dg_fragment_wannier_representative
  public::verify_dg_fragment_center_orbit
  public::verify_dg_fragment_subspace_density_covariance
  public::verify_dg_fragment_wannier_streaming_closure
  public::build_dg_core_owned_occupied_subspace
  public::verify_dg_uniform_fragment_target_rank
  public::assign_dg_overlapping_wannier_occupations
  public::find_dg_group_identity
  public::accumulate_dg_lcfo_buffer_contributions_to_core
  public::measure_dg_rank_fixed_symmetry_residuals
#ifdef USE_EIGENEXA
  public::measure_dg_rank_fixed_symmetry_residuals_eigenexa
#endif
  public::exchange_dg_point_permuted_orbital_rows
  public::accept_dg_boundary_calibrated_symmetry
  public::solve_dg_affine_common_fixed_point
  public::compute_dg_periodic_wannier_centers
  public::verify_dg_wannier_center_affine_orbits
  public::build_dg_balanced_orbital_ownership
  public::transpose_dg_spatial_cores_to_orbital_owners
  public::redistribute_dg_owned_orbitals_to_center_fragments
  public::assign_dg_periodic_centers_to_fragments
contains

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
      basis,weights,local_cyclic_metric,metric_peak,ok,detail)
    if(.not.ok)then;message='distributed rank-fixed metric: '//trim(detail);return;endif
    allocate(local_cyclic_vectors(info%nrow_local,info%ncol_local),metric_spectrum(size(basis,1)))
    call eigen_pdsyevd_ex_distributed_blocks(info,size(basis,1),local_cyclic_metric,metric_spectrum,&
      local_cyclic_vectors,eigen_ok,detail)
    if(.not.eigen_ok)then;message='distributed rank-fixed eigensystem: '//trim(detail);return;endif
    if(minval(metric_spectrum)<=epsilon(1d0)*max(1d0,maxval(metric_spectrum)))then
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
        if(mod(global_row-1,info%nprow)/=info%myrow)cycle
        local_row=(global_row-1)/info%nprow+1
        do global_col=1,nstate
          if(mod(global_col-1,info%npcol)/=info%mycol)cycle
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
      fractional_translations,tolerance,ok,message)
    real(real64),intent(in)::centers(:,:),fractional_translations(:,:),tolerance
    integer,intent(in)::integer_rotations(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::mapped_centers(:,:)
    integer,allocatable::matched_target(:)
    logical,allocatable::seen(:)
    integer::nwann,noperation,operation,source

    nwann=size(centers,2);noperation=size(integer_rotations,3)
    ok=nwann>0.and.size(centers,1)==3.and.noperation>0.and.size(integer_rotations,1)==3.and.&
      size(integer_rotations,2)==3.and.all(shape(fractional_translations)==[3,noperation]).and.&
      tolerance>0d0.and.all(ieee_is_finite(centers)).and.&
      all(ieee_is_finite(fractional_translations))
    if(.not.ok)then;message='invalid Wannier center affine-orbit contract';return;end if
    allocate(mapped_centers(3,nwann),matched_target(nwann),seen(nwann))
    do operation=1,noperation
      mapped_centers=modulo(matmul(real(integer_rotations(:,:,operation),real64),centers)+&
        spread(fractional_translations(:,operation),2,nwann),1d0)
      matched_target=0
      do source=1,nwann
        seen=.false.
        if(.not.augment_center_match(source))then
          ok=.false.;message='localized Wannier centers are not closed under full affine symmetry';return
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
      target_owner,target_point,pass,iw,rank_before,local_bad,global_bad,left,right,product,&
      source_owner,source_point,middle_owner,middle_point,effective_minimum_rank
    integer(int64)::middle_target,final_target
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
    do left=1,noperation;do right=1,noperation
      product=product_table(left,right)
      if(product<1.or.product>noperation)cycle
      do source_owner=1,nproc;do source_point=1,nlocal
        middle_target=all_maps(source_point,right,source_owner)
        middle_owner=int((middle_target-1_int64)/int(nlocal,int64))+1
        middle_point=int(modulo(middle_target-1_int64,int(nlocal,int64)))+1
        final_target=all_maps(middle_point,left,middle_owner)
        if(final_target/=all_maps(source_point,product,source_owner))local_bad=1
      end do;end do
    end do;end do
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
    integer,parameter::orbital_tile_size=32
    integer::rank,nproc,ierr,nbasis,nlocal,nsym,isym,tile_first,tile_count,owner,&
      owner_first,owner_count,base,remainder,i,j,local_bad,global_bad
    integer(int64)::persistent_bytes,tile_bytes,complex_bytes
    complex(real64),allocatable::image_tile(:,:),partial_rows(:,:),reduced_rows(:,:)

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
    allocate(row_ids(owner_count),symmetry_overlap_rows(owner_count,nbasis,nsym))
    symmetry_overlap_rows=(0d0,0d0)
    do i=1,owner_count;row_ids(i)=int(owner_first+i-1,int64);enddo
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    persistent_bytes=complex_bytes*int(size(symmetry_overlap_rows),int64)
    workspace_peak_bytes=persistent_bytes
    do isym=1,nsym
      do tile_first=1,nbasis,orbital_tile_size
        tile_count=min(orbital_tile_size,nbasis-tile_first+1)
        allocate(image_tile(tile_count,nlocal))
        call exchange_dg_point_permuted_orbital_rows(comm,&
          local_basis(tile_first:tile_first+tile_count-1,:),&
          symmetry_target_box_ids(:,isym),image_tile,ok,message)
        if(.not.ok)return
        do owner=0,nproc-1
          owner_count=base+merge(1,0,owner<remainder)
          owner_first=owner*base+min(owner,remainder)+1
          allocate(partial_rows(owner_count,tile_count),reduced_rows(owner_count,tile_count))
          partial_rows=(0d0,0d0);reduced_rows=(0d0,0d0)
          do j=1,tile_count;do i=1,owner_count
            partial_rows(i,j)=sum(weights*conjg(local_basis(owner_first+i-1,:))*image_tile(j,:))
          enddo;enddo
          call MPI_Reduce(partial_rows,reduced_rows,owner_count*tile_count,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,owner,comm,ierr)
          if(ierr/=MPI_SUCCESS)then
            ok=.false.;message='row-owned symmetry-overlap reduction failed';return
          endif
          if(rank==owner)symmetry_overlap_rows(:,tile_first:tile_first+tile_count-1,isym)=reduced_rows
          tile_bytes=complex_bytes*int(size(image_tile)+size(partial_rows)+size(reduced_rows),int64)
          workspace_peak_bytes=max(workspace_peak_bytes,persistent_bytes+tile_bytes)
          deallocate(partial_rows,reduced_rows)
        enddo
        deallocate(image_tile)
      enddo
    enddo
    ok=all(ieee_is_finite(real(symmetry_overlap_rows))).and.&
      all(ieee_is_finite(aimag(symmetry_overlap_rows)))
    if(ok)then;message='';else;message='row-owned symmetry overlap is nonfinite';endif
#else
    ok=.false.;message='row-owned symmetry overlap requires MPI';workspace_peak_bytes=0_int64
#endif
  end subroutine assemble_dg_distributed_basis_symmetry_overlap_rows

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
