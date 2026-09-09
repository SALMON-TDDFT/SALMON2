#include "config.h"
module rt_dg_hybrid_initialization
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_sparse_projection,only:validate_rt_dg_hybrid_sparse_hermiticity
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange,&
    exchange_rt_dg_sparse_matrix
  use rt_dg_hybrid_point_density,only:reconstruct_rt_dg_point_csr_density
  use rt_dg_hybrid_checkpoint_v4,only:s_rt_dg_hybrid_v4_shard,read_rt_dg_hybrid_checkpoint_v4
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter::salmon_xctype_none=0,salmon_xctype_pz=1,salmon_xctype_pzm=2,&
    salmon_xctype_pbe=3,salmon_xctype_pw=7
  integer(int64),parameter::cell_wrapped_position_convention_fingerprint=&
    int(z'43454C4C57524150',int64)
  type,public::s_rt_dg_hybrid_state
    logical::valid=.false.,initial_invariants_valid=.false.,density_freshly_reconstructed=.false.,&
      fixed_density_reference_valid=.false.
    real(real64)::startup_operator_covariance=huge(1d0),startup_projector_covariance=huge(1d0),&
      startup_orbital_residual=huge(1d0),startup_metric_defect=huge(1d0),reference_refresh_defect=0d0,&
      reference_refresh_scale=1d0
    integer::certified_rank=0,global_count=0,noccupied=0,operation_count=0,nonidentity_operation_count=0
    integer(int64)::payload_fingerprint=0_int64,operator_structure_fingerprint=0_int64,&
      operator_value_fingerprint=0_int64,scope_fingerprint=0_int64,density_workspace_peak_bytes=0_int64
    integer::density_payload_collective_count=0
    type(s_dg_hybrid_sparse_metric)::metric
    type(s_dg_hybrid_sparse_operators)::operators
    type(s_rt_dg_sparse_exchange)::basis_exchange
    integer(int64),allocatable::owned_row_ids(:),grid_ids(:)
    integer,allocatable::basis_point_offsets(:),basis_support_ids(:),basis_support_slots(:),basis_halo_ids(:)
    complex(real64),allocatable::coefficients(:,:),local_rows(:),kinetic_values(:),nonlocal_values(:),&
      sipg_values(:),basis_support_values(:),local_reference_correction(:),hamiltonian_reference_correction(:)
    real(real64),allocatable::density(:),grid_weights(:),occupations(:),eigenvalues(:),energy_receipt(:)
  end type s_rt_dg_hybrid_state
  public::initialize_rt_dg_hybrid_from_checkpoint,fingerprint_rt_dg_hybrid_scope,&
    fingerprint_rt_dg_hybrid_sparse_structure
contains
  subroutine initialize_rt_dg_hybrid_from_checkpoint(comm,path,theory,periodic,nspin,spinorbit,plus_u,hse,&
      fix_func,jm,xctype,tolerances,state,ok,message)
    integer,intent(in)::comm,nspin,xctype(:)
    character(*),intent(in)::path,theory
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    real(real64),intent(in)::tolerances(4)
    type(s_rt_dg_hybrid_state),intent(out)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_rt_dg_hybrid_v4_shard)::payload
    integer::i,j,edge,ierr,local_bad,global_bad,active_xc,nhalo,slot,payload_count
    integer(int64)::workspace_peak
    integer,allocatable::temporary_halo(:)
    complex(real64),allocatable::metric_edge_coefficients(:,:),operator_edge_coefficients(:,:),&
      s_coefficients(:,:),h_coefficients(:,:)
    real(real64),allocatable::reconstructed_density(:)
    real(real64)::local_value,global_value,local_scale,global_scale,local_residual,global_residual
    complex(real64)::local_inner,global_inner
    logical::manifest_exists,dense_v3_exists,exchange_ok
    character(256)::exchange_message
    ok=.false.;message='';state=s_rt_dg_hybrid_state();local_bad=0
    if((trim(theory)/='tddft_response'.and.trim(theory)/='tddft_pulse').or..not.periodic.or.nspin/=1.or.&
      spinorbit.or.plus_u.or.hse.or.fix_func.or.jm.or.size(xctype)<1)local_bad=1
    if(any(.not.ieee_is_finite(tolerances)).or.any(tolerances<=0d0))local_bad=1
    active_xc=0
    do i=1,size(xctype)
      if(xctype(i)/=salmon_xctype_none)active_xc=active_xc+1
      if(.not.(xctype(i)==salmon_xctype_none.or.xctype(i)==salmon_xctype_pz.or.&
        xctype(i)==salmon_xctype_pzm.or.xctype(i)==salmon_xctype_pbe.or.xctype(i)==salmon_xctype_pw))local_bad=1
    enddo
    if(active_xc==0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='unsupported local hybrid RT scope';return;endif
    inquire(file=trim(path)//'.manifest',exist=manifest_exists);inquire(file=trim(path),exist=dense_v3_exists)
    local_bad=merge(0,1,manifest_exists.or..not.dense_v3_exists)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='dense Hybrid v3 checkpoint is unsupported; regenerate distributed-native v4';return
    endif
    call read_rt_dg_hybrid_checkpoint_v4(comm,path,payload,ok,message);if(.not.ok)return
    local_bad=0
    if(size(payload%scope_selectors)/=8.or.any(payload%scope_selectors/=[1,1,1,0,0,0,0,0]).or.&
      size(payload%xc_types)/=size(xctype).or.any(payload%xc_types/=xctype).or.&
      payload%scope_fingerprint/=fingerprint_rt_dg_hybrid_scope(payload%scope_selectors,payload%xc_types))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='checkpoint/local hybrid RT scope mismatch';return;endif

    state%global_count=payload%global_count;state%certified_rank=payload%global_count
    state%noccupied=payload%nocc;state%operation_count=1;state%nonidentity_operation_count=0
    state%payload_fingerprint=payload%payload_fingerprint
    state%operator_structure_fingerprint=payload%operator_structure_fingerprint
    state%operator_value_fingerprint=payload%operator_fingerprint;state%scope_fingerprint=payload%scope_fingerprint
    allocate(state%owned_row_ids,source=payload%row_ids);allocate(state%coefficients,source=payload%initial_occupied_amplitudes)
    allocate(state%occupations,source=payload%occupations);allocate(state%eigenvalues,source=payload%eigenvalues)
    allocate(state%grid_ids,source=payload%grid_ids);allocate(state%grid_weights,source=payload%grid_weights)
    allocate(state%density,source=payload%density);allocate(state%energy_receipt,source=payload%energy_receipt)
    allocate(state%metric%owned_row_ids,source=payload%row_ids);allocate(state%metric%row_offsets,source=payload%metric_offsets)
    allocate(state%metric%column_ids,source=payload%metric_columns);allocate(state%metric%values,source=payload%metric_values)
    allocate(state%metric%active_rows(payload%global_count),state%metric%packet_ids(payload%global_count))
    state%metric%active_rows=.true.;state%metric%packet_ids=1;state%metric%global_count=payload%global_count
    state%metric%numerical_rank=payload%global_count;state%metric%fingerprint=payload%basis_fingerprint
    state%metric%condition_estimate=1d0;state%metric%maximum_value=0d0;state%metric%max_row_nnz=0
    if(size(payload%metric_values)>0)state%metric%maximum_value=maxval(abs(payload%metric_values))
    do i=1,size(payload%row_ids)
      state%metric%max_row_nnz=max(state%metric%max_row_nnz,payload%metric_offsets(i+1)-payload%metric_offsets(i))
    enddo
    state%metric%valid=.true.
    allocate(state%operators%owned_row_ids,source=payload%row_ids)
    allocate(state%operators%row_offsets,source=payload%operator_offsets)
    allocate(state%operators%column_ids,source=payload%operator_columns)
    allocate(state%operators%metric_values(size(payload%operator_columns)),source=(0d0,0d0))
    do i=1,size(payload%row_ids);do edge=payload%operator_offsets(i),payload%operator_offsets(i+1)-1
      do j=payload%metric_offsets(i),payload%metric_offsets(i+1)-1
        if(payload%metric_columns(j)==payload%operator_columns(edge))then
          state%operators%metric_values(edge)=payload%metric_values(j);exit
        endif
      enddo
    enddo;enddo
    allocate(state%operators%hamiltonian_values,source=payload%operator_values)
    allocate(state%operators%position_values,source=payload%position_values)
    state%operators%global_count=payload%global_count;state%operators%metric_fingerprint=payload%basis_fingerprint
    state%operators%position_convention_fingerprint=cell_wrapped_position_convention_fingerprint
    state%operators%fingerprint=payload%operator_fingerprint;state%operators%valid=.true.
    allocate(state%kinetic_values,source=payload%kinetic_values);allocate(state%nonlocal_values,source=payload%nonlocal_values)
    allocate(state%local_rows,source=payload%local_values);allocate(state%sipg_values,source=payload%sipg_values)
    allocate(state%basis_point_offsets,source=payload%basis_point_offsets)
    allocate(state%basis_support_ids,source=payload%basis_support_ids)
    allocate(state%basis_support_values,source=payload%basis_support_values)

    allocate(temporary_halo(size(payload%basis_support_ids)),state%basis_support_slots(size(payload%basis_support_ids)))
    nhalo=0
    do i=1,size(payload%basis_support_ids)
      slot=0
      do j=1,nhalo;if(temporary_halo(j)==payload%basis_support_ids(i))then;slot=j;exit;endif;enddo
      if(slot==0)then;nhalo=nhalo+1;temporary_halo(nhalo)=payload%basis_support_ids(i);slot=nhalo;endif
      state%basis_support_slots(i)=slot
    enddo
    allocate(state%basis_halo_ids(nhalo));state%basis_halo_ids=temporary_halo(:nhalo)
    call build_rt_dg_sparse_exchange(comm,payload%global_count,payload%basis_fingerprint,payload%row_ids,&
      state%basis_halo_ids,state%basis_exchange,ok,message)
    if(.not.ok)then;message='distributed-v4 point-support halo failed: '//trim(message);return;endif
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,payload%global_count,payload%row_ids,&
      payload%metric_offsets,payload%metric_columns,payload%metric_values,tolerances(1),ok,message)
    if(.not.ok)then;message='distributed-v4 metric Hermiticity failed: '//trim(message);return;endif
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,payload%global_count,payload%row_ids,&
      payload%operator_offsets,payload%operator_columns,payload%operator_values,tolerances(1),ok,message)
    if(.not.ok)then;message='distributed-v4 Hamiltonian Hermiticity failed: '//trim(message);return;endif

    call build_rt_dg_sparse_exchange(comm,payload%global_count,payload%basis_fingerprint,payload%row_ids,&
      payload%metric_columns,state%basis_exchange,ok,message);if(.not.ok)return
    allocate(metric_edge_coefficients(size(payload%metric_columns),payload%nocc))
    call exchange_rt_dg_sparse_matrix(comm,state%basis_exchange,state%coefficients,metric_edge_coefficients,&
      workspace_peak,payload_count,exchange_ok,exchange_message)
    if(.not.exchange_ok)then;ok=.false.;message='distributed-v4 metric coefficient halo failed';return;endif
    allocate(s_coefficients(size(payload%row_ids),payload%nocc));s_coefficients=(0d0,0d0)
    do i=1,size(payload%row_ids);do edge=payload%metric_offsets(i),payload%metric_offsets(i+1)-1
      s_coefficients(i,:)=s_coefficients(i,:)+payload%metric_values(edge)*metric_edge_coefficients(edge,:)
    enddo;enddo
    local_value=0d0;local_scale=1d0
    do i=1,payload%nocc;do j=1,payload%nocc
      local_inner=sum(conjg(state%coefficients(:,i))*s_coefficients(:,j))
      call MPI_Allreduce(local_inner,global_inner,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      local_value=max(local_value,abs(global_inner-merge((1d0,0d0),(0d0,0d0),i==j)))
      local_scale=max(local_scale,abs(global_inner))
    enddo;enddo
    state%startup_metric_defect=local_value/local_scale

    call build_rt_dg_sparse_exchange(comm,payload%global_count,payload%operator_structure_fingerprint,payload%row_ids,&
      payload%operator_columns,state%basis_exchange,ok,message);if(.not.ok)return
    allocate(operator_edge_coefficients(size(payload%operator_columns),payload%nocc))
    call exchange_rt_dg_sparse_matrix(comm,state%basis_exchange,state%coefficients,operator_edge_coefficients,&
      workspace_peak,payload_count,exchange_ok,exchange_message)
    if(.not.exchange_ok)then;ok=.false.;message='distributed-v4 operator coefficient halo failed';return;endif
    allocate(h_coefficients(size(payload%row_ids),payload%nocc));h_coefficients=(0d0,0d0)
    do i=1,size(payload%row_ids);do edge=payload%operator_offsets(i),payload%operator_offsets(i+1)-1
      h_coefficients(i,:)=h_coefficients(i,:)+payload%operator_values(edge)*operator_edge_coefficients(edge,:)
    enddo;enddo
    local_residual=0d0;local_scale=0d0
    do i=1,payload%nocc
      local_residual=local_residual+sum(abs(h_coefficients(:,i)-payload%eigenvalues(i)*s_coefficients(:,i))**2)
      local_scale=local_scale+sum(abs(h_coefficients(:,i))**2)
    enddo
    call MPI_Allreduce(local_residual,global_residual,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    state%startup_orbital_residual=sqrt(global_residual)/max(1d0,sqrt(global_scale))
    if(state%startup_metric_defect>tolerances(2).or.state%startup_orbital_residual>tolerances(3))then
      ok=.false.;message='distributed-v4 metric/stationarity certification failed';return
    endif
    call build_rt_dg_sparse_exchange(comm,payload%global_count,payload%basis_fingerprint,payload%row_ids,&
      state%basis_halo_ids,state%basis_exchange,ok,message);if(.not.ok)return
    allocate(reconstructed_density(size(payload%density)))
    call reconstruct_rt_dg_point_csr_density(comm,state%basis_exchange,state%basis_point_offsets,&
      state%basis_support_slots,state%basis_support_values,state%coefficients,state%occupations,&
      reconstructed_density,workspace_peak,payload_count,ok,message)
    if(.not.ok)return
    local_value=0d0;if(size(reconstructed_density)>0)local_value=maxval(abs(reconstructed_density-state%density))
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_value>tolerances(4))then
      ok=.false.;message='distributed-v4 checkpoint density reconstruction failed';return
    endif
    state%density_workspace_peak_bytes=workspace_peak;state%density_payload_collective_count=payload_count
    state%startup_operator_covariance=0d0;state%startup_projector_covariance=payload%acceptance_receipts(3)
    state%initial_invariants_valid=.true.;state%valid=.true.;ok=.true.;message=''
#else
    ok=.false.;message='hybrid RT initialization requires MPI'
#endif
  end subroutine initialize_rt_dg_hybrid_from_checkpoint

  integer(int64) function fingerprint_rt_dg_hybrid_scope(selectors,xctype) result(hash)
    integer,intent(in)::selectors(:),xctype(:)
    integer::i
    hash=int(z'A4093822299F31D0',int64)
    if(size(selectors)/=8)then;hash=0_int64;return;endif
    hash=ieor(ishftc(hash,7),int(selectors(1),int64));hash=ieor(ishftc(hash,7),int(selectors(3),int64))
    hash=ieor(ishftc(hash,7),int(selectors(2),int64))
    do i=4,8;hash=ieor(ishftc(hash,7),int(selectors(i),int64));enddo
    hash=ieor(ishftc(hash,7),int(size(xctype),int64))
    do i=1,size(xctype);hash=ieor(ishftc(hash,7),int(xctype(i),int64));enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_rt_dg_hybrid_scope

  subroutine fingerprint_rt_dg_hybrid_sparse_structure(comm,global_count,row_ids,row_offsets,column_ids,&
      ownership_fingerprint,position_convention_fingerprint,fingerprint,ok,message)
    integer,intent(in)::comm,global_count,row_offsets(:),column_ids(:)
    integer(int64),intent(in)::row_ids(:),ownership_fingerprint,position_convention_fingerprint
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,edge,ierr
    integer(int64)::local_hash,term
    local_hash=0_int64
    do i=1,size(row_ids)
      term=ieor(row_ids(i),ishftc(int(row_offsets(i),int64),17))
      term=ieor(term,ishftc(int(row_offsets(i+1),int64),29));local_hash=ieor(local_hash,ishftc(term,mod(7*i,63)))
      do edge=row_offsets(i),row_offsets(i+1)-1
        term=ieor(ishftc(row_ids(i),13),int(column_ids(edge),int64))
        local_hash=ieor(local_hash,ishftc(term,mod(11*edge+3*i,63)))
      enddo
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      fingerprint=0_int64;ok=.false.;message='sparse operator structure fingerprint reduction failed';return
    endif
    fingerprint=ieor(ishftc(fingerprint,11),ishftc(ownership_fingerprint,7))
    fingerprint=ieor(fingerprint,position_convention_fingerprint)
    fingerprint=ieor(fingerprint,ishftc(int(global_count,int64),37))
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.;message=''
#else
    fingerprint=0_int64;ok=.false.;message='sparse structure fingerprint requires MPI'
#endif
  end subroutine fingerprint_rt_dg_hybrid_sparse_structure
end module rt_dg_hybrid_initialization
