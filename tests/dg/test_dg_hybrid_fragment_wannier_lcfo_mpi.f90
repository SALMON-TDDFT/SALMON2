#include "config.h"
#ifdef DG_GENERALIZED_CONTRACT_SYNTAX
module dg_hybrid_generalized_red_contracts
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  type,public::s_dg_hybrid_dual_basis_catalog
    logical::valid=.false.
    type(s_dg_hybrid_fragment_basis),allocatable::fragment_bases(:)
    complex(real64),allocatable::union_to_complete(:,:)
    integer::uncompressed_rank=0
    integer::complete_rank=0
    integer(int64)::fragment_catalog_fingerprint=0_int64
    integer(int64)::complete_map_fingerprint=0_int64
    integer(int64)::complete_transform_binding_fingerprint=0_int64
    integer(int64),allocatable::uncompressed_global_basis_ids(:)
    integer,allocatable::uncompressed_owner_ranks(:),uncompressed_fragment_ids(:),&
      uncompressed_local_slots(:),uncompressed_sectors(:),uncompressed_generations(:)
  end type s_dg_hybrid_dual_basis_catalog
  interface
    subroutine compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_row_count,row_ids,weights,&
        wannier_values,pw_tile,wannier_fingerprint,packet_fingerprint,first_column,metric_tolerance,&
        coefficients,projected_pw,metric_rank,metric_condition,orthogonality_defect,workspace_peak_bytes,&
        fingerprint,ok,message)
      import::int64,real64
      integer,intent(in)::comm,global_row_count,first_column
      integer(int64),intent(in)::row_ids(:),wannier_fingerprint,packet_fingerprint
      real(real64),intent(in)::weights(:),metric_tolerance
      complex(real64),intent(in)::wannier_values(:,:),pw_tile(:,:)
      complex(real64),allocatable,intent(out)::coefficients(:,:),projected_pw(:,:)
      integer,intent(out)::metric_rank
      real(real64),intent(out)::metric_condition,orthogonality_defect
      integer(int64),intent(out)::workspace_peak_bytes,fingerprint
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine compute_dg_hybrid_generalized_wannier_projection_tile
    subroutine build_dg_hybrid_complete_union_map(comm,global_row_count,row_ids,weights,&
        uncompressed_basis_values,metric_tolerance,complete_basis_transform,complete_basis_values,&
        metric_rank,metric_condition,projector_fingerprint,ok,message)
      import::int64,real64
      integer,intent(in)::comm,global_row_count
      integer(int64),intent(in)::row_ids(:)
      real(real64),intent(in)::weights(:),metric_tolerance
      complex(real64),intent(in)::uncompressed_basis_values(:,:)
      complex(real64),allocatable,intent(out)::complete_basis_transform(:,:),complete_basis_values(:,:)
      integer,intent(out)::metric_rank
      real(real64),intent(out)::metric_condition
      integer(int64),intent(out)::projector_fingerprint
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine build_dg_hybrid_complete_union_map
    subroutine finalize_dg_hybrid_dual_basis_catalog(comm,global_row_count,row_ids,weights,fragment_bases,&
        uncompressed_global_basis_ids,uncompressed_basis_values,union_to_complete,expected_fragment_ranks,&
        maximum_terminal_rank_loss,seed_fragment_owner,seed_coefficients_in_uncompressed,metric_tolerance,&
        required_interface_fragment_ids,required_interface_row_ids,required_periodic_wrap_fragment_ids,&
        required_periodic_wrap_row_ids,required_projector_fragment_ids,required_projector_row_ids,packet_ids,&
        packet_neighbor_offsets,packet_neighbor_basis_ids,required_neighbor_packet_ids,&
        required_neighbor_basis_ids,tail_tolerance,catalog,preserved_seed_fragment_owner,&
        preserved_seed_coefficients_in_uncompressed,complete_seed_coefficients,seed_reconstruction_defect,ok,message,&
        expected_fragment_wannier_ranks)
      import::int64,real64,s_dg_hybrid_fragment_basis,s_dg_hybrid_dual_basis_catalog
      integer,intent(in)::comm,global_row_count,expected_fragment_ranks(:),maximum_terminal_rank_loss,&
        seed_fragment_owner(:),required_interface_fragment_ids(:),required_periodic_wrap_fragment_ids(:),&
        required_projector_fragment_ids(:),packet_ids(:),packet_neighbor_offsets(:),required_neighbor_packet_ids(:)
      integer(int64),intent(in)::row_ids(:),uncompressed_global_basis_ids(:),required_interface_row_ids(:),&
        required_periodic_wrap_row_ids(:),required_projector_row_ids(:),packet_neighbor_basis_ids(:),&
        required_neighbor_basis_ids(:)
      real(real64),intent(in)::weights(:),metric_tolerance,tail_tolerance
      type(s_dg_hybrid_fragment_basis),intent(in)::fragment_bases(:)
      complex(real64),intent(in)::uncompressed_basis_values(:,:),union_to_complete(:,:),&
        seed_coefficients_in_uncompressed(:,:)
      type(s_dg_hybrid_dual_basis_catalog),intent(out)::catalog
      integer,allocatable,intent(out)::preserved_seed_fragment_owner(:)
      complex(real64),allocatable,intent(out)::preserved_seed_coefficients_in_uncompressed(:,:),&
        complete_seed_coefficients(:,:)
      real(real64),intent(out)::seed_reconstruction_defect
      logical,intent(out)::ok
      character(*),intent(out)::message
      integer,intent(in),optional::expected_fragment_wannier_ranks(:)
    end subroutine finalize_dg_hybrid_dual_basis_catalog
  end interface
end module dg_hybrid_generalized_red_contracts
#endif

program test_dg_hybrid_fragment_wannier_lcfo_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
#ifdef DG_GENERALIZED_CONTRACT_SYNTAX
  use dg_hybrid_generalized_red_contracts,only:s_dg_hybrid_dual_basis_catalog,&
    compute_dg_hybrid_generalized_wannier_projection_tile,build_dg_hybrid_complete_union_map,&
    finalize_dg_hybrid_dual_basis_catalog
#else
  use dg_hybrid_wannier_complement,only:compute_dg_hybrid_generalized_wannier_projection_tile,&
    build_dg_hybrid_complete_union_map
  use dg_hybrid_projected_fragment_pipeline,only:s_dg_hybrid_dual_basis_catalog,&
    finalize_dg_hybrid_dual_basis_catalog
#endif
  implicit none
  integer,parameter::global_ngrid=8,nfragment=2,nwf=4,npw=2,nunion=nwf+npw,nseed=2,nocc=2
  integer,parameter::fragment_ranks(nfragment)=[4,2]
  real(real64),parameter::metric_tolerance=1d-8,comparison_tolerance=3d-7,&
    near_null_amplitude=1d-6,projected_tail_amplitude=1d-3,tail_tolerance=1d-5
  real(real64),parameter::global_weights(global_ngrid)=[0.50_real64,0.75_real64,1.00_real64,1.25_real64,&
    1.50_real64,1.75_real64,2.00_real64,2.25_real64]
  integer(int64),parameter::canonical_union_basis_ids(nunion)=[101_int64,205_int64,309_int64,450_int64,&
    701_int64,990_int64]
  integer,parameter::required_interface_fragment_ids(1)=[2],required_wrap_fragment_ids(1)=[2],&
    required_projector_fragment_ids(1)=[1]
  integer(int64),parameter::required_interface_rows(1)=[3_int64],required_wrap_rows(1)=[1_int64],&
    required_projector_rows(1)=[5_int64]
  integer,parameter::coverage_packet_ids(npw)=[41,73],required_neighbor_packet_ids(2)=[41,73]
  integer(int64),parameter::required_neighbor_basis_ids(2)=[309_int64,450_int64]

  type::s_case_result
    integer::projection_rank=0,metric_rank=0
    real(real64)::projection_condition=0d0,metric_condition=0d0
    real(real64)::projection_defect=huge(1d0),generalized_residual=huge(1d0),&
      generalized_orthogonality=huge(1d0)
    integer(int64)::projection_fingerprint=0_int64,map_fingerprint=0_int64,catalog_fingerprint=0_int64
    complex(real64),allocatable::span_projector(:,:),occupied_projector(:,:),seed_projector(:,:),projected_pw(:,:)
    real(real64),allocatable::eigenvalues(:),density(:),seed_density(:)
    integer,allocatable::basis_owner(:)
  end type s_case_result

  integer::comm,rank,nproc,ierr,nlocal,position,global_row,gauge_case
  integer(int64),allocatable::row_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::base_wf(:,:),base_pw(:,:),local_wf(:,:),local_pw(:,:),physical_seeds(:,:)
  type(s_case_result)::reference,current

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call require(nproc>=1.and.nproc<=global_ngrid,'fixture rank count exceeds dynamic row extent')
  nlocal=count([(mod(global_row-1,nproc)==rank,global_row=1,global_ngrid)])
  allocate(row_ids(nlocal),weights(nlocal),base_wf(nwf,global_ngrid),base_pw(npw,global_ngrid),&
    local_wf(nwf,nlocal),local_pw(npw,nlocal),physical_seeds(nseed,global_ngrid))
  call make_reference_functions(base_wf,base_pw,physical_seeds)
  position=0
  do global_row=global_ngrid,1,-1
    if(mod(global_row-1,nproc)/=rank)cycle
    position=position+1
    row_ids(position)=int(global_row,int64);weights(position)=global_weights(global_row)
    local_wf(:,position)=base_wf(:,global_row);local_pw(:,position)=base_pw(:,global_row)
  enddo

  call run_invariant_case(0,reference)
  call check_row_reorder_fingerprint_invariance(reference)
  do gauge_case=1,3
    call run_invariant_case(gauge_case,current)
    call compare_invariant_results(reference,current,gauge_case)
  enddo
  call run_full_rank_identity_case(reference%map_fingerprint)
  call run_negative_contracts
  if(rank==0)then
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)')'HYBRID_FRAGMENT_WANNIER_LCFO ranks=',nproc,&
      ' metric_rank=',reference%metric_rank,' generalized_fingerprint=',reference%projection_fingerprint,&
      ' map_fingerprint=',reference%map_fingerprint,' catalog_fingerprint=',reference%catalog_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid fragment Wannier LCFO on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains

  subroutine make_reference_functions(wf_values,pw_values,seed_values)
    complex(real64),intent(out)::wf_values(nwf,global_ngrid),pw_values(npw,global_ngrid),&
      seed_values(nseed,global_ngrid)
    real(real64)::normalization
    wf_values=(0d0,0d0)
    wf_values(1,1)=1d0/sqrt(global_weights(1));wf_values(2,2)=1d0/sqrt(global_weights(2))
    wf_values(3,1)=1d0/sqrt(2d0*global_weights(1))
    wf_values(3,3)=cmplx(0d0,1d0,real64)/sqrt(2d0*global_weights(3))
    normalization=sqrt(1d0+near_null_amplitude**2)
    wf_values(4,2)=1d0/(normalization*sqrt(global_weights(2)))
    wf_values(4,4)=near_null_amplitude/(normalization*sqrt(global_weights(4)))
    pw_values=(0d0,0d0)
    pw_values(1,:)=0.30d0*wf_values(1,:)+(0.20d0,-0.15d0)*wf_values(3,:)
    pw_values(1,5)=pw_values(1,5)+1d0/sqrt(global_weights(5))
    pw_values(1,7)=pw_values(1,7)+projected_tail_amplitude/sqrt(global_weights(7))
    pw_values(2,:)=(-0.25d0,0.10d0)*wf_values(2,:)+0.35d0*wf_values(4,:)
    pw_values(2,6)=pw_values(2,6)+(0.6d0,0.8d0)/sqrt(global_weights(6))
    pw_values(2,8)=pw_values(2,8)+cmplx(0d0,projected_tail_amplitude,real64)/sqrt(global_weights(8))
    seed_values(1,:)=wf_values(1,:);seed_values(2,:)=wf_values(3,:)
  end subroutine make_reference_functions

  subroutine make_block_gauge(gauge_kind,gauge)
    integer,intent(in)::gauge_kind
    complex(real64),intent(out)::gauge(nwf,nwf)
    real(real64)::c,s
    integer::i
    gauge=(0d0,0d0);do i=1,nwf;gauge(i,i)=1d0;enddo
    select case(gauge_kind)
    case(0)
    case(1)
      gauge(1,1)=exp(cmplx(0d0,0.31d0,real64));gauge(2,2)=exp(cmplx(0d0,-0.77d0,real64))
      gauge(3,3)=exp(cmplx(0d0,1.13d0,real64));gauge(4,4)=exp(cmplx(0d0,-0.48d0,real64))
    case(2)
      gauge=(0d0,0d0);gauge(2,1)=1d0;gauge(1,2)=1d0
      gauge(4,3)=(0d0,1d0);gauge(3,4)=(0d0,-1d0)
    case(3)
      gauge=(0d0,0d0);c=cos(0.37d0);s=sin(0.37d0)
      gauge(1,1)=c;gauge(2,1)=s;gauge(1,2)=-s;gauge(2,2)=c
      c=cos(0.63d0);s=sin(0.63d0)
      gauge(3,3)=c;gauge(4,3)=cmplx(0d0,s,real64)
      gauge(3,4)=cmplx(0d0,s,real64);gauge(4,4)=c
    case default
      error stop 'unknown block gauge'
    end select
  end subroutine make_block_gauge

  subroutine run_invariant_case(gauge_kind,result)
    integer,intent(in)::gauge_kind
    type(s_case_result),intent(out)::result
    complex(real64)::gauge(nwf,nwf),seed_coefficients(nunion,nseed),saved_seed_coefficients(nunion,nseed)
    complex(real64),allocatable::rotated_wf(:,:),coefficients(:,:),projected(:,:),hybrid_union(:,:),&
      complete_transform(:,:),complete_values(:,:),expected_complete_values(:,:),complete_seed_coefficients(:,:),&
      expected_complete_seed_coefficients(:,:),preserved_seed_coefficients(:,:),full_wf(:,:),full_projected(:,:),&
      full_hybrid(:,:),full_complete(:,:),expected_projected(:,:),wf_projector(:,:),raw_union_projector(:,:),&
      expected_coefficients(:,:),complete_projector(:,:),uncompressed_seed_local(:,:),complete_seed_local(:,:),full_uncompressed_seed(:,:),&
      full_complete_seed(:,:),overlap(:,:),hamiltonian(:,:),terminal_overlap(:,:),&
      terminal_hamiltonian(:,:),eigenvectors(:,:),psi(:,:),weighted_psi(:,:),weighted_seed(:,:)
    real(real64),allocatable::wf_metric_eigenvalues(:),union_metric_eigenvalues(:),complete_metric_eigenvalues(:)
    real(real64)::projection_condition,map_condition,projection_defect,seed_defect,&
      hdiag(global_ngrid),occupations(nocc)
    integer::projection_rank,metric_rank,oracle_rank,complete_oracle_rank,i,j
    integer(int64)::workspace,projection_fingerprint,map_fingerprint
    integer::expected_fragment_ranks(nfragment),seed_fragment_owner(nseed)
    integer,allocatable::preserved_seed_owner(:)
    type(s_dg_hybrid_fragment_basis),allocatable::fragment_bases(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    logical::ok
    character(256)::message

    call make_block_gauge(gauge_kind,gauge)
    allocate(rotated_wf(nwf,nlocal));rotated_wf=matmul(transpose(gauge),local_wf)
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,row_ids,weights,rotated_wf,&
      local_pw,1001_int64,2001_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok,'generalized WF projection failed: '//trim(message))
    call require(projection_rank==nwf-1,'WF-union near-null direction was not removed by G+')
    call require(all(shape(coefficients)==[nwf,npw]).and.all(shape(projected)==[npw,nlocal]),&
      'generalized projection returned invalid distributed shapes')
    call require(projection_condition>1d0.and.ieee_is_finite(projection_condition).and.&
      projection_defect<comparison_tolerance.and.workspace>0_int64.and.projection_fingerprint/=0_int64,&
      'generalized projection receipts are invalid')
    call check_generalized_projection_residual(rotated_wf,projected,projection_defect)

    call gather_state_rows(rotated_wf,full_wf);call gather_state_rows(projected,full_projected)
    call raw_projection_oracle(full_wf,base_pw,global_weights,metric_tolerance,expected_coefficients,&
      expected_projected,wf_projector,oracle_rank,wf_metric_eigenvalues)
    call require(oracle_rank==projection_rank,'independent weighted WF rank disagrees with receipt')
    call require(ieee_is_finite(projection_condition).and.&
      relative_scalar_defect(projection_condition,retained_condition(wf_metric_eigenvalues,metric_tolerance))<&
      comparison_tolerance,'generalized projection condition disagrees with the independent retained Gram')
    call require(maxval(abs(full_projected-expected_projected))<comparison_tolerance,&
      'generalized projection differs from independent weighted Moore-Penrose oracle')
    call require(maxval(abs(coefficients-expected_coefficients))<comparison_tolerance,&
      'generalized projection coefficients are not the minimum-norm G+K solution')

    allocate(hybrid_union(nunion,nlocal));hybrid_union(:nwf,:)=rotated_wf
    hybrid_union(nwf+1:nunion,:)=projected
    call gather_state_rows(hybrid_union,full_hybrid)
    call build_weighted_mp_projector(full_hybrid,global_weights,metric_tolerance,raw_union_projector,&
      oracle_rank,union_metric_eigenvalues)
    call require(oracle_rank==nunion-1,'full WF+projected-PW Hybrid union has wrong metric rank')
    call require(has_degenerate_retained_cluster(union_metric_eigenvalues,metric_tolerance),&
      'full Hybrid union lacks an exactly degenerate retained Gram cluster')

    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,weights,hybrid_union,metric_tolerance,&
      complete_transform,complete_values,metric_rank,map_condition,map_fingerprint,ok,message)
    call require(ok,'full Hybrid terminal union map failed: '//trim(message))
    call require(metric_rank==oracle_rank.and.all(shape(complete_transform)==[nunion,metric_rank]).and.&
      all(shape(complete_values)==[metric_rank,nlocal]),'terminal Hybrid map returned incorrect ranks or shapes')
    call require(ieee_is_finite(map_condition).and.&
      relative_scalar_defect(map_condition,retained_condition(union_metric_eigenvalues,metric_tolerance))<&
      comparison_tolerance,'terminal map condition disagrees with the independent retained Gram')
    call check_terminal_transform(full_hybrid,complete_transform)
    allocate(expected_complete_values(metric_rank,nlocal))
    expected_complete_values=matmul(transpose(complete_transform),hybrid_union)
    call require(maxval(abs(complete_values-expected_complete_values))<5d-14,&
      'complete_values does not equal transpose(T) times uncompressed Hybrid union')
    call gather_state_rows(complete_values,full_complete)
    call build_weighted_mp_projector(full_complete,global_weights,metric_tolerance,complete_projector,&
      complete_oracle_rank,complete_metric_eigenvalues)
    call require(complete_oracle_rank==metric_rank.and.&
      maxval(abs(complete_projector-raw_union_projector))<comparison_tolerance,&
      'terminal map changed independent raw weighted-union physical projector')

    seed_coefficients=(0d0,0d0);seed_coefficients(1,1)=1d0;seed_coefficients(3,2)=1d0
    seed_coefficients(:nwf,:)=matmul(conjg(transpose(gauge)),seed_coefficients(:nwf,:))
    saved_seed_coefficients=seed_coefficients;seed_fragment_owner=[1,2]
    expected_fragment_ranks=fragment_ranks
    call make_fragment_catalogs(hybrid_union,gauge_kind,fragment_bases)
    call assert_catalogs_are_sparse(fragment_bases)
    call finalize_with_evidence(fragment_bases,hybrid_union,complete_transform,expected_fragment_ranks,1,&
      seed_fragment_owner,seed_coefficients,0,catalog,preserved_seed_owner,preserved_seed_coefficients,&
      complete_seed_coefficients,seed_defect,ok,message)
    call require(ok,'dual Hybrid catalog finalization failed: '//trim(message))
    call require(catalog%valid.and.catalog%uncompressed_rank==nunion.and.catalog%complete_rank==metric_rank,&
      'dual catalog did not preserve full uncompressed Hybrid rank')
    call require(catalog%fragment_catalog_fingerprint/=0_int64.and.catalog%complete_map_fingerprint/=0_int64.and.&
      catalog%complete_transform_binding_fingerprint/=0_int64,&
      'dual catalog did not derive physical fingerprints from content')
    call require(allocated(catalog%union_to_complete).and.&
      maxval(abs(catalog%union_to_complete-complete_transform))<5d-14,&
      'dual catalog changed immutable terminal map')
    call require(catalog_metadata_equal(fragment_bases,catalog%fragment_bases),&
      'dual catalog rewrote fragment ID/generation/order/sector/points/provenance/values')
    call require(all(preserved_seed_owner==seed_fragment_owner).and.&
      bitwise_complex_equal(preserved_seed_coefficients,saved_seed_coefficients),&
      'dual catalog rewrote uncompressed seed-owner coordinates')
    call verify_catalog_tuples(catalog)
    result%catalog_fingerprint=catalog%fragment_catalog_fingerprint
    if(gauge_kind==0)call check_catalog_fingerprint_sensitivity(fragment_bases,hybrid_union,complete_transform,&
      expected_fragment_ranks,seed_fragment_owner,seed_coefficients,catalog)

    call verify_sparse_catalog_values(catalog%fragment_bases,full_hybrid)
    call metric_composed_seed_coefficients(full_hybrid,complete_transform,seed_coefficients,global_weights,&
      expected_complete_seed_coefficients)
    call require(maxval(abs(complete_seed_coefficients-expected_complete_seed_coefficients))<comparison_tolerance,&
      'complete seed coefficients do not use (T^H G T)^+ T^H G A')
    allocate(uncompressed_seed_local(nseed,nlocal),complete_seed_local(nseed,nlocal))
    uncompressed_seed_local=matmul(transpose(seed_coefficients),hybrid_union)
    complete_seed_local=matmul(transpose(complete_seed_coefficients),complete_values)
    call gather_state_rows(uncompressed_seed_local,full_uncompressed_seed)
    call gather_state_rows(complete_seed_local,full_complete_seed)
    call require(weighted_state_defect(full_uncompressed_seed,physical_seeds,global_weights)<comparison_tolerance,&
      'fragment embedding does not reconstruct weighted physical seeds')
    call require(weighted_state_defect(full_complete_seed,physical_seeds,global_weights)<comparison_tolerance.and.&
      weighted_state_defect(full_complete_seed,full_uncompressed_seed,global_weights)<comparison_tolerance.and.&
      seed_defect<comparison_tolerance,'terminal composition does not reconstruct weighted physical seeds')

    allocate(overlap(nunion,nunion),hamiltonian(nunion,nunion))
    call weighted_gram(full_hybrid,global_weights,overlap)
    hdiag=[0.35d0,0.82d0,1.47d0,2.31d0,3.40d0,4.80d0,6.10d0,7.70d0]
    hamiltonian=(0d0,0d0)
    do j=1,nunion;do i=1,nunion
      hamiltonian(i,j)=sum(global_weights*hdiag*conjg(full_hybrid(i,:))*full_hybrid(j,:))
    enddo;enddo
    call require(abs(overlap(1,3))>0.1d0.and.abs(hamiltonian(1,3))>0.01d0,&
      'fixture lacks nonzero cross-fragment Hybrid H/S blocks')
    terminal_overlap=matmul(conjg(transpose(complete_transform)),matmul(overlap,complete_transform))
    terminal_hamiltonian=matmul(conjg(transpose(complete_transform)),matmul(hamiltonian,complete_transform))
    call solve_generalized(terminal_hamiltonian,terminal_overlap,result%eigenvalues,eigenvectors,&
      result%generalized_residual,result%generalized_orthogonality)
    call require(result%generalized_residual<comparison_tolerance.and.&
      result%generalized_orthogonality<comparison_tolerance,&
      'terminal full-Hybrid generalized residual/orthogonality receipts are invalid')

    allocate(result%span_projector,source=raw_union_projector)
    allocate(result%projected_pw,source=full_projected)
    allocate(result%occupied_projector(global_ngrid,global_ngrid),result%density(global_ngrid),&
      result%seed_projector(global_ngrid,global_ngrid),result%seed_density(global_ngrid))
    occupations=[2d0,1d0];allocate(psi(nocc,global_ngrid),weighted_psi(nocc,global_ngrid))
    psi=matmul(transpose(eigenvectors(:,1:nocc)),full_complete)
    do i=1,global_ngrid;weighted_psi(:,i)=sqrt(global_weights(i))*psi(:,i);enddo
    result%occupied_projector=matmul(transpose(weighted_psi),conjg(weighted_psi))
    result%density=occupations(1)*abs(psi(1,:))**2+occupations(2)*abs(psi(2,:))**2
    allocate(weighted_seed(nseed,global_ngrid))
    do i=1,global_ngrid;weighted_seed(:,i)=sqrt(global_weights(i))*full_complete_seed(:,i);enddo
    result%seed_projector=matmul(transpose(weighted_seed),conjg(weighted_seed))
    result%seed_density=2d0*abs(full_complete_seed(1,:))**2+abs(full_complete_seed(2,:))**2
    result%projection_rank=projection_rank;result%metric_rank=metric_rank
    result%projection_condition=projection_condition;result%metric_condition=map_condition
    result%projection_defect=projection_defect;result%projection_fingerprint=projection_fingerprint
    result%map_fingerprint=map_fingerprint
    call derive_catalog_ownership(catalog%fragment_bases,result%basis_owner)
  end subroutine run_invariant_case

  subroutine compare_invariant_results(reference_result,rotated_result,gauge_kind)
    type(s_case_result),intent(in)::reference_result,rotated_result
    integer,intent(in)::gauge_kind
    character(32)::label
    write(label,'(a,i0)')'gauge case ',gauge_kind
    call require(rotated_result%projection_rank==reference_result%projection_rank.and.&
      rotated_result%metric_rank==reference_result%metric_rank,trim(label)//' changed retained ranks')
    call require(maxval(abs(rotated_result%span_projector-reference_result%span_projector))<comparison_tolerance,&
      trim(label)//' changed raw weighted-union physical projector')
    call require(maxval(abs(rotated_result%projected_pw-reference_result%projected_pw))<comparison_tolerance,&
      trim(label)//' changed generalized projected PWs')
    call require(maxval(abs(rotated_result%eigenvalues-reference_result%eigenvalues))<comparison_tolerance,&
      trim(label)//' changed terminal full-Hybrid generalized eigenvalues')
    call require(maxval(abs(rotated_result%occupied_projector-reference_result%occupied_projector))<&
      comparison_tolerance,trim(label)//' changed occupied physical projector')
    call require(maxval(abs(rotated_result%density-reference_result%density))<comparison_tolerance,&
      trim(label)//' changed reconstructed density')
    call require(maxval(abs(rotated_result%seed_projector-reference_result%seed_projector))<comparison_tolerance.and.&
      maxval(abs(rotated_result%seed_density-reference_result%seed_density))<comparison_tolerance,&
      trim(label)//' changed composed seed projector or density')
    call require(abs(rotated_result%generalized_residual-reference_result%generalized_residual)<comparison_tolerance.and.&
      abs(rotated_result%generalized_orthogonality-reference_result%generalized_orthogonality)<comparison_tolerance,&
      trim(label)//' changed generalized residual/orthogonality receipts')
    call require(all(rotated_result%basis_owner==reference_result%basis_owner),&
      trim(label)//' changed ownership derived from returned sparse catalogs')
    ! Degenerate retained clusters permit different transform columns and raw
    ! map fingerprints; only physical projectors and observables are compared.
    call require(rotated_result%projection_fingerprint/=0_int64.and.rotated_result%map_fingerprint/=0_int64,&
      trim(label)//' lost physical fingerprints')
  end subroutine compare_invariant_results

  subroutine check_row_reorder_fingerprint_invariance(reference_result)
    type(s_case_result),intent(in)::reference_result
    integer(int64),allocatable::permuted_ids(:)
    real(real64),allocatable::permuted_weights(:)
    complex(real64),allocatable::permuted_wf(:,:),permuted_pw(:,:),coefficients(:,:),projected(:,:),&
      union_values(:,:),transform(:,:),complete(:,:)
    integer::p,source,projection_rank,metric_rank
    integer(int64)::workspace,projection_fingerprint,map_fingerprint
    real(real64)::projection_condition,projection_defect,map_condition
    logical::ok
    character(256)::message
    allocate(permuted_ids(nlocal),permuted_weights(nlocal),permuted_wf(nwf,nlocal),permuted_pw(npw,nlocal))
    do p=1,nlocal
      source=nlocal-p+1;permuted_ids(p)=row_ids(source);permuted_weights(p)=weights(source)
      permuted_wf(:,p)=local_wf(:,source);permuted_pw(:,p)=local_pw(:,source)
    enddo
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,permuted_ids,permuted_weights,&
      permuted_wf,permuted_pw,1001_int64,2001_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok.and.projection_rank==reference_result%projection_rank.and.&
      projection_fingerprint==reference_result%projection_fingerprint,&
      'local spatial-row reordering changed generalized physical fingerprint or rank')
    do p=1,nlocal
      call require(maxval(abs(projected(:,p)-reference_result%projected_pw(:,int(permuted_ids(p)))))<&
        comparison_tolerance,'local row reorder changed projected PW row correspondence')
    enddo
    permuted_pw(1,:)=1.125d0*permuted_pw(1,:)
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,permuted_ids,permuted_weights,&
      permuted_wf,permuted_pw,1001_int64,2001_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok.and.projection_fingerprint/=reference_result%projection_fingerprint,&
      'generalized projection fingerprint is insensitive to changed PW payload content')
    do p=1,nlocal
      source=nlocal-p+1;permuted_pw(:,p)=local_pw(:,source)
    enddo
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,permuted_ids,permuted_weights,&
      permuted_wf,permuted_pw,1001_int64,2001_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok,'row-reorder projection rebuild failed: '//trim(message))
    allocate(union_values(nunion,nlocal));union_values(:nwf,:)=permuted_wf;union_values(nwf+1:,:)=projected
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,permuted_ids,permuted_weights,union_values,&
      metric_tolerance,transform,complete,metric_rank,map_condition,map_fingerprint,ok,message)
    call require(ok.and.metric_rank==reference_result%metric_rank.and.map_fingerprint==reference_result%map_fingerprint,&
      'local spatial-row reordering changed terminal-map physical fingerprint or rank')
  end subroutine check_row_reorder_fingerprint_invariance

  subroutine check_catalog_fingerprint_sensitivity(bases,union_values,transform,expected_ranks,seed_owner,&
      seed_coefficients,baseline)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::union_values(:,:),transform(:,:),seed_coefficients(:,:)
    integer,intent(in)::expected_ranks(:),seed_owner(:)
    type(s_dg_hybrid_dual_basis_catalog),intent(in)::baseline
    type(s_dg_hybrid_fragment_basis),allocatable::changed_bases(:)
    type(s_dg_hybrid_dual_basis_catalog)::changed_catalog
    complex(real64),allocatable::changed_transform(:,:),changed_union(:,:),changed_seed(:,:),&
      preserved_seed(:,:),complete_seed(:,:)
    integer,allocatable::preserved_owner(:)
    complex(real64)::phase
    real(real64)::seed_defect
    integer::b
    logical::ok
    character(256)::message

    changed_bases=bases
    do b=1,size(changed_bases)
      if(changed_bases(b)%fragment_id>0)changed_bases(b)%generation=changed_bases(b)%generation+1
    enddo
    call finalize_with_evidence(changed_bases,union_values,transform,expected_ranks,1,seed_owner,seed_coefficients,0,&
      changed_catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    call require(ok.and.changed_catalog%fragment_catalog_fingerprint/=baseline%fragment_catalog_fingerprint.and.&
      changed_catalog%complete_map_fingerprint==baseline%complete_map_fingerprint.and.&
      changed_catalog%complete_transform_binding_fingerprint==baseline%complete_transform_binding_fingerprint,&
      'dual fingerprints are insensitive to fragment generation metadata or couple it into the terminal map')

    allocate(changed_transform,source=transform);phase=cmplx(0d0,1d0,real64)
    changed_transform(:,1)=phase*changed_transform(:,1)
    call finalize_with_evidence(bases,union_values,changed_transform,expected_ranks,1,seed_owner,seed_coefficients,0,&
      changed_catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    call require(ok.and.changed_catalog%fragment_catalog_fingerprint==baseline%fragment_catalog_fingerprint.and.&
      changed_catalog%complete_map_fingerprint/=baseline%complete_map_fingerprint.and.&
      changed_catalog%complete_transform_binding_fingerprint/=baseline%complete_transform_binding_fingerprint,&
      'dual fingerprints are insensitive to terminal-map payload or couple it into fragment metadata')

    changed_bases=bases;allocate(changed_union,source=union_values);allocate(changed_seed,source=seed_coefficients)
    phase=exp(cmplx(0d0,0.19d0,real64));changed_union(1,:)=phase*changed_union(1,:)
    changed_seed(1,:)=conjg(phase)*changed_seed(1,:)
    call phase_catalog_column(changed_bases,canonical_union_basis_ids(1),phase)
    deallocate(changed_transform)
    allocate(changed_transform,source=transform);changed_transform(1,:)=conjg(phase)*changed_transform(1,:)
    call finalize_with_evidence(changed_bases,changed_union,changed_transform,expected_ranks,1,seed_owner,changed_seed,0,&
      changed_catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    call require(ok.and.changed_catalog%fragment_catalog_fingerprint/=baseline%fragment_catalog_fingerprint,&
      'fragment catalog fingerprint is insensitive to basis-value content')
  end subroutine check_catalog_fingerprint_sensitivity

  subroutine run_full_rank_identity_case(rank_deficient_map_fingerprint)
    integer(int64),intent(in)::rank_deficient_map_fingerprint
    complex(real64),allocatable::full_rank_wf(:,:),coefficients(:,:),projected(:,:),union_values(:,:),&
      complete_transform(:,:),complete_values(:,:),full_union(:,:),gram(:,:),preserved_seed_coefficients(:,:),&
      complete_seed_coefficients(:,:),changed_union(:,:)
    complex(real64)::seed_coefficients(nunion,nseed)
    integer::p,projection_rank,metric_rank,expected_ranks(nfragment),seed_owner(nseed)
    integer,allocatable::preserved_owner(:)
    integer(int64)::workspace,projection_fingerprint,map_fingerprint
    real(real64)::projection_condition,projection_defect,map_condition,seed_defect,local_magnitude(2),&
      global_magnitude(2)
    type(s_dg_hybrid_fragment_basis),allocatable::bases(:),changed_bases(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog,changed_catalog
    logical::ok,fingerprints_are_sensitive
    character(256)::message

    allocate(full_rank_wf(nwf,nlocal));full_rank_wf=local_wf
    full_rank_wf(4,:)=(0d0,0d0)
    do p=1,nlocal
      if(row_ids(p)==4_int64)full_rank_wf(4,p)=1d0/sqrt(global_weights(4))
    enddo
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,row_ids,weights,full_rank_wf,&
      local_pw,1101_int64,2101_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok,'full-rank generalized projection failed: '//trim(message))
    allocate(union_values(nunion,nlocal));union_values(:nwf,:)=full_rank_wf;union_values(nwf+1:,:)=projected
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,weights,union_values,metric_tolerance,&
      complete_transform,complete_values,metric_rank,map_condition,map_fingerprint,ok,message)
    call require(ok.and.metric_rank==nunion,'full-rank Hybrid union map was rejected')
    call require(map_fingerprint/=rank_deficient_map_fingerprint,&
      'terminal-map fingerprint is insensitive to a changed physical union payload')
    call gather_state_rows(union_values,full_union);allocate(gram(nunion,nunion))
    call weighted_gram(full_union,global_weights,gram)
    call require(maxval(abs(gram-diagonal_matrix(real(diagonal(gram),real64))))>0.1d0,&
      'full-rank identity-map fixture accidentally has an identity Gram')
    call require(identity_transform_bitwise(complete_transform),&
      'full-rank terminal map was whitened instead of returning bitwise identity')
    call require(bitwise_complex_equal(complete_values,union_values),&
      'full-rank terminal values are not a bitwise copy of the uncompressed Hybrid union')

    call make_fragment_catalogs(union_values,0,bases)
    call assert_catalogs_are_sparse(bases)
    seed_coefficients=(0d0,0d0);seed_coefficients(1,1)=1d0;seed_coefficients(3,2)=1d0
    seed_owner=[1,2];expected_ranks=fragment_ranks
    call finalize_with_evidence(bases,union_values,complete_transform,expected_ranks,0,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call require(ok.and.catalog%complete_rank==nunion.and.identity_transform_bitwise(catalog%union_to_complete),&
      'dual catalog did not retain the immutable bitwise identity map for a full-rank nonidentity Gram')

    ! Exercise values above unity on both fingerprint paths.  The publisher
    ! payload and canonical row-distributed union describe the same scaled
    ! basis, so this remains a valid full-rank identity-map catalog.
    allocate(changed_union,source=union_values);changed_union(1,:)=2d0*changed_union(1,:)
    local_magnitude=[maxval(abs(union_values(1,:))),maxval(abs(changed_union(1,:)))]
    call MPI_Allreduce(local_magnitude,global_magnitude,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call require(ierr==MPI_SUCCESS.and.all(global_magnitude>1d0).and.&
      global_magnitude(2)>global_magnitude(1),'magnitude fingerprint fixture did not exercise distinct values above one')
    call make_fragment_catalogs(changed_union,0,changed_bases)
    call finalize_with_evidence(changed_bases,changed_union,complete_transform,expected_ranks,0,seed_owner,&
      seed_coefficients,0,changed_catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call require(ok,'consistent magnitude-changed dual catalog was rejected: '//trim(message))
    fingerprints_are_sensitive=&
      changed_catalog%fragment_catalog_fingerprint/=catalog%fragment_catalog_fingerprint.and.&
      changed_catalog%complete_map_fingerprint/=catalog%complete_map_fingerprint
    if(rank==0.and..not.fingerprints_are_sensitive)write(*,'(a)')&
      'EXPECTED RED: magnitude>1 publisher/map payload changes collided in catalog fingerprints'
    call require(fingerprints_are_sensitive,&
      'dual catalog fingerprints are insensitive to distinct publisher/map payload magnitudes above one')
  end subroutine run_full_rank_identity_case

  subroutine run_negative_contracts
    complex(real64),allocatable::coefficients(:,:),projected(:,:),hybrid_union(:,:),complete_transform(:,:),&
      complete_values(:,:),preserved_seed_coefficients(:,:),complete_seed_coefficients(:,:),full_union(:,:),&
      lost_seed_values(:,:),lost_seed_full(:,:),seed_gram(:,:),seed_gram_inverse(:,:)
    complex(real64)::seed_coefficients(nunion,nseed),lost_seed_coefficients(nunion,nseed)
    real(real64),allocatable::seed_eigenvalues(:)
    real(real64)::projection_condition,projection_defect,map_condition,seed_defect,global_norm,local_norm
    integer::projection_rank,metric_rank,expected_ranks(nfragment),seed_owner(nseed),seed_rank
    integer,allocatable::preserved_owner(:)
    integer(int64)::workspace,projection_fingerprint,map_fingerprint
    integer(int64)::duplicate_basis_ids(nunion),unsorted_basis_ids(nunion)
    type(s_dg_hybrid_fragment_basis),allocatable::bases(:),damaged(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    logical::ok
    character(256)::message

    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,row_ids,weights,local_wf,&
      local_pw,1201_int64,2201_int64,1,metric_tolerance,coefficients,projected,projection_rank,&
      projection_condition,projection_defect,workspace,projection_fingerprint,ok,message)
    call require(ok,'negative-case setup projection failed: '//trim(message))
    allocate(hybrid_union(nunion,nlocal));hybrid_union(:nwf,:)=local_wf;hybrid_union(nwf+1:,:)=projected
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,weights,hybrid_union,metric_tolerance,&
      complete_transform,complete_values,metric_rank,map_condition,map_fingerprint,ok,message)
    call require(ok.and.metric_rank==nunion-1,'negative-case terminal-map setup failed')
    call make_fragment_catalogs(hybrid_union,0,bases)
    seed_coefficients=(0d0,0d0);seed_coefficients(1,1)=1d0;seed_coefficients(3,2)=1d0
    seed_owner=[1,2];expected_ranks=fragment_ranks

    call run_fragment_metadata_quality_negatives(bases,hybrid_union,complete_transform,expected_ranks,&
      seed_owner,seed_coefficients)
    call run_zero_valued_interface_coverage_positive(hybrid_union,seed_owner,seed_coefficients)
    call run_terminal_span_negative(bases,hybrid_union,complete_transform,expected_ranks,seed_owner,&
      seed_coefficients)
    call run_numerical_row_negatives(hybrid_union)
    call run_missing_wf_negative

    duplicate_basis_ids=canonical_union_basis_ids;duplicate_basis_ids(2)=duplicate_basis_ids(1)
    call finalize_with_evidence(bases,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message,duplicate_basis_ids)
    call expect_failure(ok,message,'duplicate','duplicate global Hybrid basis ID was accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'duplicate uncompressed ID failure partially published a catalog')

    unsorted_basis_ids=canonical_union_basis_ids
    unsorted_basis_ids(1:2)=[canonical_union_basis_ids(2),canonical_union_basis_ids(1)]
    call finalize_with_evidence(bases,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message,unsorted_basis_ids)
    call expect_failure_any(ok,message,'canonical','order','unsorted uncompressed Hybrid basis IDs were accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'unsorted uncompressed ID failure partially published a catalog')

    damaged=bases;call remove_catalog_column(damaged,1,205_int64)
    call finalize_with_evidence(damaged,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call expect_failure_any(ok,message,'rank','missing','missing local fragment WF column was silently accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'missing fragment column failure partially published a catalog')

    damaged=bases;call replace_catalog_id(damaged,2,450_int64,8081_int64)
    call finalize_with_evidence(damaged,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call expect_failure_any(ok,message,'basis','publisher','publisher ID inconsistent with canonical union was accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'publisher-ID mismatch partially published a catalog')

    damaged=bases;call replace_catalog_id(damaged,2,450_int64,101_int64)
    call finalize_with_evidence(damaged,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call expect_failure(ok,message,'duplicate','duplicate canonical publisher ID was accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'duplicate publisher failure partially published a catalog')

    call expect_coverage_failure(bases,hybrid_union,complete_transform,seed_owner,seed_coefficients,&
      2,3_int64,0,'interface')
    call expect_coverage_failure(bases,hybrid_union,complete_transform,seed_owner,seed_coefficients,&
      2,1_int64,0,'wrap')
    call expect_coverage_failure(bases,hybrid_union,complete_transform,seed_owner,seed_coefficients,&
      1,5_int64,0,'projector')
    call expect_coverage_failure(bases,hybrid_union,complete_transform,seed_owner,seed_coefficients,&
      0,0_int64,1,'neighbor')
    call prove_projected_tail_is_material(hybrid_union,1,7_int64)
    call expect_coverage_failure(bases,hybrid_union,complete_transform,seed_owner,seed_coefficients,&
      1,7_int64,0,'tail')

    ! The second seed is the normalized cross-fragment metric-null direction.
    ! Its physical Gram is full rank before terminal compression, so accepting
    ! it would prove that only an absolute (unnormalized) residual was checked.
    lost_seed_coefficients=(0d0,0d0);lost_seed_coefficients(1,1)=1d0
    lost_seed_coefficients(2,2)=1d0;lost_seed_coefficients(4,2)=-1d0
    allocate(lost_seed_values(nseed,nlocal))
    lost_seed_values=matmul(transpose(lost_seed_coefficients),hybrid_union)
    local_norm=sum(weights*abs(lost_seed_values(2,:))**2)
    call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call require(global_norm>0d0.and.ieee_is_finite(global_norm),&
      'cross-fragment lost seed has zero or nonfinite pre-normalization norm')
    lost_seed_coefficients(:,2)=lost_seed_coefficients(:,2)/sqrt(global_norm)
    lost_seed_values=matmul(transpose(lost_seed_coefficients),hybrid_union)
    call gather_state_rows(lost_seed_values,lost_seed_full)
    allocate(seed_gram(nseed,nseed));call weighted_gram(lost_seed_full,global_weights,seed_gram)
    call hermitian_pseudoinverse(seed_gram,metric_tolerance,seed_gram_inverse,seed_rank,seed_eigenvalues)
    call require(seed_rank==nseed.and.minval(seed_eigenvalues)>0.9d0,&
      'normalized seed-loss fixture is not full physical rank before compression')
    call finalize_with_evidence(bases,hybrid_union,complete_transform,expected_ranks,1,seed_owner,&
      lost_seed_coefficients,0,catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,ok,message)
    call expect_failure(ok,message,'seed','terminal compression accepted loss of a normalized physical seed')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'seed-rank failure partially published a catalog')

    call run_excessive_rank_loss_negative(seed_owner,seed_coefficients)
    call run_indefinite_metric_negative(hybrid_union)
    call run_huge_tail_tolerance_positive(bases,hybrid_union,complete_transform,expected_ranks,seed_owner,&
      seed_coefficients)
  end subroutine run_negative_contracts

  subroutine run_fragment_metadata_quality_negatives(bases,union_values,union_to_complete,expected_ranks,&
      seed_owner,seed_coefficients)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::union_values(:,:),union_to_complete(:,:),seed_coefficients(:,:)
    integer,intent(in)::expected_ranks(:),seed_owner(:)
    type(s_dg_hybrid_fragment_basis),allocatable::damaged(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer,allocatable::preserved_owner(:)
    complex(real64),allocatable::preserved_seed(:,:),complete_seed(:,:)
    real(real64)::seed_defect
    integer::b
    logical::ok,generation_rejected,sector_rejected,global_generation_rejected,global_sector_rejected
    character(256)::message

    damaged=bases
    do b=1,size(damaged)
      if(damaged(b)%fragment_id==1)damaged(b)%generation=damaged(b)%generation+1
    enddo
    call finalize_with_evidence(damaged,union_values,union_to_complete,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    generation_rejected=.not.ok.and.index(lowercase(message),'generation')>0.and.&
      catalog_outputs_unpublished(catalog,preserved_owner,preserved_seed,complete_seed,seed_defect)
    call MPI_Allreduce(generation_rejected,global_generation_rejected,1,MPI_LOGICAL,MPI_LAND,comm,ierr)

    damaged=bases
    call replace_catalog_sector(damaged,1,101_int64,2)
    call finalize_with_evidence(damaged,union_values,union_to_complete,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    sector_rejected=.not.ok.and.(index(lowercase(message),'wannier')>0.or.&
      index(lowercase(message),'sector')>0).and.&
      catalog_outputs_unpublished(catalog,preserved_owner,preserved_seed,complete_seed,seed_defect)
    call MPI_Allreduce(sector_rejected,global_sector_rejected,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(rank==0.and..not.global_generation_rejected)write(0,'(a)')&
      'EXPECTED RED: one-fragment generation mismatch was accepted or partially published'
    if(rank==0.and..not.global_sector_rejected)write(0,'(a)')&
      'EXPECTED RED: fragment Wannier-sector rank loss was accepted or partially published'
    call require(global_generation_rejected.and.global_sector_rejected,&
      'dual catalog fragment generation/Wannier-rank quality contracts are not implemented')
  end subroutine run_fragment_metadata_quality_negatives

  subroutine run_zero_valued_interface_coverage_positive(union_values,seed_owner,seed_coefficients)
    complex(real64),intent(in)::union_values(:,:),seed_coefficients(:,:)
    integer,intent(in)::seed_owner(:)
    complex(real64),allocatable::zero_union(:,:),zero_transform(:,:),zero_complete(:,:),full_zero(:,:),&
      preserved_seed(:,:),complete_seed(:,:)
    type(s_dg_hybrid_fragment_basis),allocatable::zero_bases(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer::p,b,point_position,metric_rank,expected_ranks(nfragment),local_publishers,global_publishers
    integer,allocatable::preserved_owner(:)
    integer(int64)::map_fingerprint
    real(real64)::map_condition,seed_defect
    logical::ok,local_evidence_ok,global_evidence_ok
    character(256)::message

    allocate(zero_union,source=union_values)
    do p=1,nlocal
      if(row_ids(p)==4_int64)zero_union(:,p)=(0d0,0d0)
    enddo
    call gather_state_rows(zero_union,full_zero)
    call require(maxval(abs(full_zero(:,4)))==0d0,&
      'zero-valued interface fixture retained a nonzero canonical basis value')
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,weights,zero_union,metric_tolerance,&
      zero_transform,zero_complete,metric_rank,map_condition,map_fingerprint,ok,message)
    call require(ok.and.metric_rank==nunion-1,&
      'zero-valued interface fixture changed the retained physical rank: '//trim(message))
    call make_fragment_catalogs(zero_union,0,zero_bases)
    local_publishers=0;local_evidence_ok=.true.
    do b=1,size(zero_bases)
      if(zero_bases(b)%fragment_id/=2)cycle
      local_publishers=local_publishers+1
      point_position=findloc(zero_bases(b)%buffer_point_ids,4_int64,dim=1)
      if(point_position<1)then
        local_evidence_ok=.false.
      elseif(any(zero_bases(b)%buffer_values(point_position,:)/=(0d0,0d0)))then
        local_evidence_ok=.false.
      endif
    enddo
    call MPI_Allreduce(local_publishers,global_publishers,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_evidence_ok,global_evidence_ok,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    call require(ierr==MPI_SUCCESS.and.global_publishers==1.and.global_evidence_ok,&
      'zero-valued interface fixture lacks exact publisher row-ID evidence')
    expected_ranks=fragment_ranks
    call finalize_with_evidence(zero_bases,zero_union,zero_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message,&
      interface_rows_override=[4_int64])
    if(rank==0.and..not.ok)write(0,'(a,a)')&
      'EXPECTED RED: an existing required interface row with an all-zero payload was rejected; diagnostic=',&
      trim(message)
    call require(ok.and.catalog%valid,&
      'required interface coverage incorrectly depends on a nonzero basis value: '//trim(message))
  end subroutine run_zero_valued_interface_coverage_positive

  subroutine run_terminal_span_negative(bases,union_values,union_to_complete,expected_ranks,seed_owner,&
      seed_coefficients)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::union_values(:,:),union_to_complete(:,:),seed_coefficients(:,:)
    integer,intent(in)::expected_ranks(:),seed_owner(:)
    complex(real64),allocatable::full_union(:,:),gram(:,:),gram_inverse(:,:),expected_projector(:,:),&
      original_projector(:,:),seed_coordinates(:,:),seed_coordinate_gram(:,:),seed_coordinate_inverse(:,:),&
      seed_orthogonal_projector(:,:),direction(:),rotation(:,:),rotated_transform(:,:),null_projector(:,:),&
      null_direction(:),bad_transform(:,:),column_gram(:,:),identity(:,:),retained_gram(:,:),&
      retained_inverse(:,:),bad_projector(:,:),bad_complete_values(:,:),bad_seed_coefficients(:,:),&
      original_seed_values(:,:),reconstructed_seed_values(:,:),preserved_seed(:,:),complete_seed(:,:)
    real(real64),allocatable::gram_eigenvalues(:),seed_coordinate_eigenvalues(:),retained_eigenvalues(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer::gram_rank,seed_coordinate_rank,retained_rank,ncomplete,nuncompressed,j,best_column
    integer,allocatable::preserved_owner(:)
    real(real64)::column_norm,best_norm,span_defect,reconstruction_defect,seed_defect,angle
    logical::ok,rejected
    character(256)::message

    call gather_state_rows(union_values,full_union)
    nuncompressed=size(union_to_complete,1);ncomplete=size(union_to_complete,2)
    allocate(gram(nuncompressed,nuncompressed));call weighted_gram(full_union,global_weights,gram)
    call hermitian_pseudoinverse(gram,metric_tolerance,gram_inverse,gram_rank,gram_eigenvalues)
    call require(gram_rank==ncomplete,'bad terminal-span fixture received the wrong union metric rank')
    expected_projector=matmul(gram,gram_inverse)
    original_projector=matmul(union_to_complete,conjg(transpose(union_to_complete)))
    call require(maxval(abs(original_projector-expected_projector))<comparison_tolerance,&
      'bad terminal-span fixture did not start from the physical retained projector')

    seed_coordinates=matmul(conjg(transpose(union_to_complete)),seed_coefficients)
    seed_coordinate_gram=matmul(conjg(transpose(seed_coordinates)),seed_coordinates)
    call hermitian_pseudoinverse(seed_coordinate_gram,metric_tolerance,seed_coordinate_inverse,&
      seed_coordinate_rank,seed_coordinate_eigenvalues)
    call require(seed_coordinate_rank==size(seed_coefficients,2),&
      'bad terminal-span fixture seeds are not independent retained coordinates')
    allocate(seed_orthogonal_projector(ncomplete,ncomplete));seed_orthogonal_projector=(0d0,0d0)
    do j=1,ncomplete;seed_orthogonal_projector(j,j)=1d0;enddo
    seed_orthogonal_projector=seed_orthogonal_projector-&
      matmul(seed_coordinates,matmul(seed_coordinate_inverse,conjg(transpose(seed_coordinates))))
    best_column=1;best_norm=-1d0
    do j=1,ncomplete
      column_norm=sqrt(max(0d0,real(dot_product(seed_orthogonal_projector(:,j),&
        seed_orthogonal_projector(:,j)),real64)))
      if(column_norm>best_norm)then;best_norm=column_norm;best_column=j;endif
    enddo
    allocate(direction(ncomplete));direction=seed_orthogonal_projector(:,best_column)/best_norm
    call require(best_norm>0.5d0.and.&
      maxval(abs(matmul(conjg(transpose(seed_coordinates)),direction)))<comparison_tolerance,&
      'cannot isolate a retained direction orthogonal to every physical seed')
    call complete_unitary_from_first(direction,rotation)
    rotated_transform=matmul(union_to_complete,rotation)

    allocate(null_projector(nuncompressed,nuncompressed));null_projector=(0d0,0d0)
    do j=1,nuncompressed;null_projector(j,j)=1d0;enddo
    null_projector=null_projector-original_projector
    best_column=1;best_norm=-1d0
    do j=1,nuncompressed
      column_norm=sqrt(max(0d0,real(dot_product(null_projector(:,j),null_projector(:,j)),real64)))
      if(column_norm>best_norm)then;best_norm=column_norm;best_column=j;endif
    enddo
    allocate(null_direction(nuncompressed));null_direction=null_projector(:,best_column)/best_norm
    call require(best_norm>0.5d0.and.&
      maxval(abs(matmul(conjg(transpose(union_to_complete)),null_direction)))<comparison_tolerance,&
      'cannot isolate the numerical-null coefficient direction')

    allocate(bad_transform,source=rotated_transform);angle=0.35d0
    bad_transform(:,1)=cos(angle)*rotated_transform(:,1)+sin(angle)*null_direction
    column_gram=matmul(conjg(transpose(bad_transform)),bad_transform)
    allocate(identity(ncomplete,ncomplete));identity=(0d0,0d0)
    do j=1,ncomplete;identity(j,j)=1d0;enddo
    retained_gram=matmul(conjg(transpose(bad_transform)),matmul(gram,bad_transform))
    call hermitian_pseudoinverse(retained_gram,metric_tolerance,retained_inverse,retained_rank,&
      retained_eigenvalues)
    bad_projector=matmul(bad_transform,conjg(transpose(bad_transform)))
    span_defect=maxval(abs(bad_projector-expected_projector))
    call require(maxval(abs(column_gram-identity))<comparison_tolerance.and.retained_rank==ncomplete.and.&
      span_defect>1d-2,'bad terminal-span fixture did not isolate TT^H completeness from rank/orthogonality')

    call metric_composed_seed_coefficients(full_union,bad_transform,seed_coefficients,global_weights,&
      bad_seed_coefficients)
    bad_complete_values=matmul(transpose(bad_transform),full_union)
    original_seed_values=matmul(transpose(seed_coefficients),full_union)
    reconstructed_seed_values=matmul(transpose(bad_seed_coefficients),bad_complete_values)
    reconstruction_defect=weighted_state_defect(original_seed_values,reconstructed_seed_values,global_weights)
    call require(reconstruction_defect<comparison_tolerance,&
      'bad terminal-span fixture accidentally loses a physical seed before the completeness gate')

    call finalize_with_evidence(bases,union_values,bad_transform,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message)
    rejected=.not.ok.and.len_trim(message)>0.and.&
      catalog_outputs_unpublished(catalog,preserved_owner,preserved_seed,complete_seed,seed_defect)
    if(rank==0.and..not.rejected)write(0,'(a,es12.4,a,a)')&
      'EXPECTED RED: terminal map with the wrong TT^H projector was accepted; defect=',span_defect,&
      ' diagnostic=',trim(message)
    call require(rejected,&
      'terminal map completeness failure was not collective and publication-safe')
  end subroutine run_terminal_span_negative

  subroutine complete_unitary_from_first(first_column,unitary)
    complex(real64),intent(in)::first_column(:)
    complex(real64),allocatable,intent(out)::unitary(:,:)
    complex(real64),allocatable::candidate(:)
    integer::n,column,basis_column,pass,j
    real(real64)::norm
    n=size(first_column);allocate(unitary(n,n),candidate(n));unitary=(0d0,0d0)
    norm=sqrt(max(0d0,real(dot_product(first_column,first_column),real64)))
    call require(norm>0d0,'cannot complete a zero direction to a unitary matrix')
    unitary(:,1)=first_column/norm;column=1
    do basis_column=1,n
      candidate=(0d0,0d0);candidate(basis_column)=1d0
      do pass=1,2
        do j=1,column
          candidate=candidate-unitary(:,j)*dot_product(unitary(:,j),candidate)
        enddo
      enddo
      norm=sqrt(max(0d0,real(dot_product(candidate,candidate),real64)))
      if(norm<=1d-10)cycle
      column=column+1;unitary(:,column)=candidate/norm
      if(column==n)exit
    enddo
    call require(column==n,'failed to complete the seed-orthogonal direction to a unitary matrix')
  end subroutine complete_unitary_from_first

  subroutine run_huge_tail_tolerance_positive(bases,union_values,union_to_complete,expected_ranks,&
      seed_owner,seed_coefficients)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::union_values(:,:),union_to_complete(:,:),seed_coefficients(:,:)
    integer,intent(in)::expected_ranks(:),seed_owner(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer,allocatable::preserved_owner(:)
    complex(real64),allocatable::preserved_seed(:,:),complete_seed(:,:)
    real(real64)::seed_defect,large_tail_tolerance
    logical::ok
    character(256)::message
    large_tail_tolerance=2d0*sqrt(huge(1d0))
    call require(ieee_is_finite(large_tail_tolerance).and.large_tail_tolerance>sqrt(huge(1d0)),&
      'huge-tail fixture did not cross the unsafe squaring threshold')
    call finalize_with_evidence(bases,union_values,union_to_complete,expected_ranks,1,seed_owner,&
      seed_coefficients,0,catalog,preserved_owner,preserved_seed,complete_seed,seed_defect,ok,message,&
      tail_tolerance_override=large_tail_tolerance)
    if(rank==0.and..not.ok)write(0,'(a,a)')&
      'EXPECTED RED: a huge finite tail tolerance was not processed safely; diagnostic=',trim(message)
    call require(ok.and.catalog%valid.and.ieee_is_finite(seed_defect),&
      'huge finite tail tolerance was rejected instead of being handled without squaring overflow')
  end subroutine run_huge_tail_tolerance_positive

  subroutine run_numerical_row_negatives(union_values)
    complex(real64),intent(in)::union_values(:,:)
    integer::kept,p
    integer(int64),allocatable::bad_ids(:)
    real(real64),allocatable::bad_weights(:)
    complex(real64),allocatable::bad_wf(:,:),bad_pw(:,:),bad_union(:,:),coefficients(:,:),projected(:,:),&
      transform(:,:),complete(:,:)
    integer::metric_rank
    integer(int64)::workspace,fingerprint
    real(real64)::condition,defect
    logical::ok
    character(256)::message

    kept=count(row_ids/=int(global_ngrid,int64))
    allocate(bad_ids(kept),bad_weights(kept),bad_wf(nwf,kept),bad_pw(npw,kept),bad_union(nunion,kept))
    kept=0
    do p=1,nlocal
      if(row_ids(p)==int(global_ngrid,int64))cycle
      kept=kept+1;bad_ids(kept)=row_ids(p);bad_weights(kept)=weights(p)
      bad_wf(:,kept)=local_wf(:,p);bad_pw(:,kept)=local_pw(:,p);bad_union(:,kept)=union_values(:,p)
    enddo
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,bad_ids,bad_weights,bad_wf,&
      bad_pw,1301_int64,2301_int64,1,metric_tolerance,coefficients,projected,metric_rank,condition,defect,&
      workspace,fingerprint,ok,message)
    call expect_failure(ok,message,'row','generalized projection accepted a missing spatial row')
    call expect_projection_unpublished(coefficients,projected,metric_rank,condition,defect,workspace,fingerprint,&
      'missing-row generalized projection exposed partial outputs')
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,bad_ids,bad_weights,bad_union,metric_tolerance,&
      transform,complete,metric_rank,condition,fingerprint,ok,message)
    call expect_failure(ok,message,'row','terminal map accepted a missing spatial row')
    call expect_map_unpublished(transform,complete,metric_rank,condition,fingerprint,&
      'missing-row terminal map exposed partial outputs')

    deallocate(bad_ids,bad_weights,bad_wf,bad_pw,bad_union)
    allocate(bad_ids(nlocal),bad_weights(nlocal),bad_wf(nwf,nlocal),bad_pw(npw,nlocal),bad_union(nunion,nlocal))
    bad_ids=row_ids;bad_weights=weights;bad_wf=local_wf;bad_pw=local_pw;bad_union=union_values
    do p=1,nlocal
      if(bad_ids(p)==int(global_ngrid,int64))bad_ids(p)=1_int64
    enddo
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,bad_ids,bad_weights,bad_wf,&
      bad_pw,1302_int64,2302_int64,1,metric_tolerance,coefficients,projected,metric_rank,condition,defect,&
      workspace,fingerprint,ok,message)
    call expect_failure(ok,message,'row','generalized projection accepted duplicate spatial rows')
    call expect_projection_unpublished(coefficients,projected,metric_rank,condition,defect,workspace,fingerprint,&
      'duplicate-row generalized projection exposed partial outputs')
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,bad_ids,bad_weights,bad_union,metric_tolerance,&
      transform,complete,metric_rank,condition,fingerprint,ok,message)
    call expect_failure(ok,message,'row','terminal map accepted duplicate spatial rows')
    call expect_map_unpublished(transform,complete,metric_rank,condition,fingerprint,&
      'duplicate-row terminal map exposed partial outputs')
  end subroutine run_numerical_row_negatives

  subroutine run_missing_wf_negative
    complex(real64),allocatable::bad_wf(:,:),coefficients(:,:),projected(:,:)
    integer::bad_count,metric_rank
    integer(int64)::workspace,fingerprint
    real(real64)::condition,defect
    logical::ok
    character(256)::message
    ! With one rank, nwf=3 is a valid self-describing API input because this
    ! numerical routine intentionally has no expected-WF catalog argument.
    ! Missing identity is tested below by the dual finalizer; here we test only
    ! the detectable collective rank disagreement.
    if(nproc==1)return
    bad_count=nwf-merge(1,0,rank==0);allocate(bad_wf(bad_count,nlocal))
    bad_wf=local_wf(:bad_count,:)
    call compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_ngrid,row_ids,weights,bad_wf,&
      local_pw,1303_int64,2303_int64,1,metric_tolerance,coefficients,projected,metric_rank,condition,defect,&
      workspace,fingerprint,ok,message)
    call expect_failure_any(ok,message,'wannier','rank','rank-disagreeing missing local WF column was accepted')
    call expect_projection_unpublished(coefficients,projected,metric_rank,condition,defect,workspace,fingerprint,&
      'missing-WF generalized projection exposed partial outputs')
  end subroutine run_missing_wf_negative

  subroutine run_indefinite_metric_negative(union_values)
    complex(real64),intent(in)::union_values(:,:)
    real(real64),allocatable::bad_weights(:)
    complex(real64),allocatable::transform(:,:),complete(:,:)
    integer::p,metric_rank
    integer(int64)::fingerprint
    real(real64)::condition
    logical::ok
    character(256)::message
    allocate(bad_weights,source=weights)
    do p=1,nlocal
      if(row_ids(p)==1_int64)bad_weights(p)=-bad_weights(p)
    enddo
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,bad_weights,union_values,metric_tolerance,&
      transform,complete,metric_rank,condition,fingerprint,ok,message)
    call expect_failure_any(ok,message,'indefinite','weight','an indefinite weighted union Gram was accepted')
    call expect_map_unpublished(transform,complete,metric_rank,condition,fingerprint,&
      'indefinite terminal map exposed partial outputs')
  end subroutine run_indefinite_metric_negative

  subroutine run_excessive_rank_loss_negative(seed_owner,seed_coefficients)
    integer,intent(in)::seed_owner(:)
    complex(real64),intent(in)::seed_coefficients(:,:)
    complex(real64),allocatable::collapsed(:,:),transform(:,:),complete(:,:),preserved_seed_coefficients(:,:),&
      complete_seed_coefficients(:,:)
    type(s_dg_hybrid_fragment_basis),allocatable::bases(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer::p,metric_rank,expected_ranks(nfragment)
    integer,allocatable::preserved_owner(:)
    integer(int64)::fingerprint
    real(real64)::condition,seed_defect
    logical::ok
    character(256)::message
    allocate(collapsed(nunion,nlocal))
    do p=1,nunion;collapsed(p,:)=local_wf(1,:);enddo
    call build_dg_hybrid_complete_union_map(comm,global_ngrid,row_ids,weights,collapsed,metric_tolerance,&
      transform,complete,metric_rank,condition,fingerprint,ok,message)
    call require(ok.and.metric_rank==1,'excessive-rank-loss setup did not produce a rank-one terminal union')
    call make_fragment_catalogs(collapsed,0,bases);expected_ranks=fragment_ranks
    call finalize_with_evidence(bases,collapsed,transform,expected_ranks,1,seed_owner,seed_coefficients,0,catalog,&
      preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,seed_defect,ok,message)
    call expect_failure(ok,message,'rank','terminal map accepted rank loss beyond the declared tolerance')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'excessive rank loss partially published a catalog')
  end subroutine run_excessive_rank_loss_negative

  subroutine make_fragment_catalogs(local_union,gauge_kind,bases)
    complex(real64),intent(in)::local_union(:,:)
    integer,intent(in)::gauge_kind
    type(s_dg_hybrid_fragment_basis),allocatable,intent(out)::bases(:)
    complex(real64),allocatable::full_union(:,:)
    integer::fragment,owner,j,p,slot,ncolumn
    integer,allocatable::slots(:),point_order(:)
    integer(int64)::bits
    call gather_state_rows(local_union,full_union);allocate(bases(nfragment))
    do fragment=1,nfragment
      owner=mod(fragment-1,nproc)
      if(rank/=owner)then
        allocate(bases(fragment)%global_ids(0),bases(fragment)%sector(0),&
          bases(fragment)%buffer_point_ids(0),bases(fragment)%buffer_values(0,0))
        cycle
      endif
      if(fragment==1)then
        allocate(slots,source=[6,5,2,1])
        allocate(point_order,source=[8,1,6,2,5,7])
      else
        allocate(slots,source=[4,3])
        allocate(point_order,source=[3,1,2,4])
      endif
      ncolumn=size(slots);bases(fragment)%fragment_id=fragment;bases(fragment)%generation=17
      allocate(bases(fragment)%global_ids(ncolumn),bases(fragment)%sector(ncolumn),&
        bases(fragment)%buffer_point_ids(size(point_order)),bases(fragment)%buffer_values(size(point_order),ncolumn))
      bases(fragment)%buffer_point_ids=int(point_order,int64)
      do j=1,ncolumn
        slot=slots(j);bases(fragment)%global_ids(j)=canonical_union_basis_ids(slot)
        bases(fragment)%sector(j)=merge(1,2,slot<=nwf)
        do p=1,size(point_order)
          bases(fragment)%buffer_values(p,j)=full_union(slot,point_order(p))
        enddo
      enddo
      bits=int(104729*fragment+1009*gauge_kind,int64)
      do j=1,ncolumn
        bits=ieor(bits,bases(fragment)%global_ids(j))
        bits=ieor(bits,ishft(int(bases(fragment)%sector(j),int64),modulo(7*j,31)))
      enddo
      if(bits==0_int64)bits=int(7919+fragment,int64)
      bases(fragment)%provenance_fingerprint=bits
      deallocate(slots,point_order)
    enddo
  end subroutine make_fragment_catalogs

  subroutine finalize_with_evidence(bases,union_values,union_to_complete,expected_ranks,maximum_rank_loss,&
      seed_owner,seed_coefficients,packet_mode,catalog,preserved_owner,preserved_seed_coefficients,&
      complete_seed_coefficients,seed_defect,ok,message,basis_ids_override,interface_rows_override,&
      tail_tolerance_override)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::union_values(:,:),union_to_complete(:,:),seed_coefficients(:,:)
    integer,intent(in)::expected_ranks(:),maximum_rank_loss,seed_owner(:),packet_mode
    type(s_dg_hybrid_dual_basis_catalog),intent(out)::catalog
    integer,allocatable,intent(out)::preserved_owner(:)
    complex(real64),allocatable,intent(out)::preserved_seed_coefficients(:,:),complete_seed_coefficients(:,:)
    real(real64),intent(out)::seed_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(in),optional::basis_ids_override(:)
    integer(int64),intent(in),optional::interface_rows_override(:)
    real(real64),intent(in),optional::tail_tolerance_override
    integer::packet_offsets(3)
    integer(int64),allocatable::neighbor_ids(:),interface_rows(:)
    integer(int64)::basis_ids(nunion)
    real(real64)::effective_tail_tolerance

    basis_ids=canonical_union_basis_ids
    if(present(basis_ids_override))basis_ids=basis_ids_override
    allocate(interface_rows,source=required_interface_rows)
    if(present(interface_rows_override))interface_rows=interface_rows_override
    effective_tail_tolerance=tail_tolerance
    if(present(tail_tolerance_override))effective_tail_tolerance=tail_tolerance_override
    if(packet_mode==0)then
      packet_offsets=[1,4,7]
      allocate(neighbor_ids,source=[101_int64,309_int64,701_int64,205_int64,450_int64,990_int64])
    else
      ! Still a structurally valid nonempty CSR graph.  Only the required
      ! packet-41 -> basis-309 edge is absent.
      packet_offsets=[1,3,6]
      allocate(neighbor_ids,source=[101_int64,701_int64,205_int64,450_int64,990_int64])
    endif
    call finalize_dg_hybrid_dual_basis_catalog(comm,global_ngrid,row_ids,weights,bases,basis_ids,union_values,&
      union_to_complete,expected_ranks,maximum_rank_loss,seed_owner,seed_coefficients,metric_tolerance,&
      required_interface_fragment_ids,interface_rows,required_wrap_fragment_ids,required_wrap_rows,&
      required_projector_fragment_ids,required_projector_rows,coverage_packet_ids,packet_offsets,neighbor_ids,&
      required_neighbor_packet_ids,required_neighbor_basis_ids,effective_tail_tolerance,catalog,preserved_owner,&
      preserved_seed_coefficients,complete_seed_coefficients,seed_defect,ok,message,&
      expected_fragment_wannier_ranks=[2,2])
  end subroutine finalize_with_evidence

  subroutine expect_coverage_failure(original_bases,union_values,union_to_complete,seed_owner,seed_coefficients,&
      fragment,row_id,packet_mode,cause)
    type(s_dg_hybrid_fragment_basis),intent(in)::original_bases(:)
    complex(real64),intent(in)::union_values(:,:),union_to_complete(:,:),seed_coefficients(:,:)
    integer,intent(in)::seed_owner(:),fragment,packet_mode
    integer(int64),intent(in)::row_id
    character(*),intent(in)::cause
    type(s_dg_hybrid_fragment_basis),allocatable::damaged(:)
    type(s_dg_hybrid_dual_basis_catalog)::catalog
    integer::expected_ranks(nfragment)
    integer,allocatable::preserved_owner(:)
    complex(real64),allocatable::preserved_seed_coefficients(:,:),complete_seed_coefficients(:,:)
    real(real64)::seed_defect
    logical::ok
    character(256)::message
    damaged=original_bases
    if(fragment>0)call remove_catalog_row(damaged,fragment,row_id)
    expected_ranks=fragment_ranks
    call finalize_with_evidence(damaged,union_values,union_to_complete,expected_ranks,1,seed_owner,&
      seed_coefficients,packet_mode,catalog,preserved_owner,preserved_seed_coefficients,&
      complete_seed_coefficients,seed_defect,ok,message)
    call expect_failure(ok,message,cause,'insufficient '//trim(cause)//' coverage was accepted')
    call expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,complete_seed_coefficients,&
      seed_defect,'coverage failure partially published a dual catalog')
  end subroutine expect_coverage_failure

  subroutine remove_catalog_row(bases,fragment,row_id)
    type(s_dg_hybrid_fragment_basis),intent(inout)::bases(:)
    integer,intent(in)::fragment
    integer(int64),intent(in)::row_id
    integer::b,p,kept
    integer(int64),allocatable::new_ids(:)
    complex(real64),allocatable::new_values(:,:)
    do b=1,size(bases)
      if(bases(b)%fragment_id/=fragment)cycle
      kept=count(bases(b)%buffer_point_ids/=row_id)
      allocate(new_ids(kept),new_values(kept,size(bases(b)%global_ids)));kept=0
      do p=1,size(bases(b)%buffer_point_ids)
        if(bases(b)%buffer_point_ids(p)==row_id)cycle
        kept=kept+1;new_ids(kept)=bases(b)%buffer_point_ids(p);new_values(kept,:)=bases(b)%buffer_values(p,:)
      enddo
      call move_alloc(new_ids,bases(b)%buffer_point_ids);call move_alloc(new_values,bases(b)%buffer_values)
      return
    enddo
  end subroutine remove_catalog_row

  subroutine remove_catalog_column(bases,fragment,global_id)
    type(s_dg_hybrid_fragment_basis),intent(inout)::bases(:)
    integer,intent(in)::fragment
    integer(int64),intent(in)::global_id
    integer::b,j,kept
    integer(int64),allocatable::new_ids(:)
    integer,allocatable::new_sector(:)
    complex(real64),allocatable::new_values(:,:)
    do b=1,size(bases)
      if(bases(b)%fragment_id/=fragment)cycle
      kept=count(bases(b)%global_ids/=global_id)
      allocate(new_ids(kept),new_sector(kept),new_values(size(bases(b)%buffer_point_ids),kept));kept=0
      do j=1,size(bases(b)%global_ids)
        if(bases(b)%global_ids(j)==global_id)cycle
        kept=kept+1;new_ids(kept)=bases(b)%global_ids(j);new_sector(kept)=bases(b)%sector(j)
        new_values(:,kept)=bases(b)%buffer_values(:,j)
      enddo
      call move_alloc(new_ids,bases(b)%global_ids);call move_alloc(new_sector,bases(b)%sector)
      call move_alloc(new_values,bases(b)%buffer_values);return
    enddo
  end subroutine remove_catalog_column

  subroutine replace_catalog_id(bases,fragment,old_id,new_id)
    type(s_dg_hybrid_fragment_basis),intent(inout)::bases(:)
    integer,intent(in)::fragment
    integer(int64),intent(in)::old_id,new_id
    integer::b,j
    do b=1,size(bases)
      if(bases(b)%fragment_id/=fragment)cycle
      do j=1,size(bases(b)%global_ids)
        if(bases(b)%global_ids(j)==old_id)bases(b)%global_ids(j)=new_id
      enddo
    enddo
  end subroutine replace_catalog_id

  subroutine replace_catalog_sector(bases,fragment,global_id,new_sector)
    type(s_dg_hybrid_fragment_basis),intent(inout)::bases(:)
    integer,intent(in)::fragment,new_sector
    integer(int64),intent(in)::global_id
    integer::b,j
    do b=1,size(bases)
      if(bases(b)%fragment_id/=fragment)cycle
      do j=1,size(bases(b)%global_ids)
        if(bases(b)%global_ids(j)==global_id)bases(b)%sector(j)=new_sector
      enddo
    enddo
  end subroutine replace_catalog_sector

  subroutine phase_catalog_column(bases,global_id,phase)
    type(s_dg_hybrid_fragment_basis),intent(inout)::bases(:)
    integer(int64),intent(in)::global_id
    complex(real64),intent(in)::phase
    integer::b,j
    do b=1,size(bases)
      do j=1,size(bases(b)%global_ids)
        if(bases(b)%global_ids(j)==global_id)bases(b)%buffer_values(:,j)=phase*bases(b)%buffer_values(:,j)
      enddo
    enddo
  end subroutine phase_catalog_column

  subroutine prove_projected_tail_is_material(union_values,fragment,row_id)
    complex(real64),intent(in)::union_values(:,:)
    integer,intent(in)::fragment
    integer(int64),intent(in)::row_id
    complex(real64),allocatable::full_union(:,:)
    integer::slot,row
    real(real64)::omitted_correction_norm
    slot=merge(nwf+1,nwf+2,fragment==1);row=int(row_id)
    call gather_state_rows(union_values,full_union)
    omitted_correction_norm=sqrt(global_weights(row)*abs(full_union(slot,row))**2)
    call require(omitted_correction_norm>tail_tolerance,&
      'tail negative did not omit a material projected-PW weighted correction norm')
  end subroutine prove_projected_tail_is_material

  subroutine gather_state_rows(local_values,global_values)
    complex(real64),intent(in)::local_values(:,:)
    complex(real64),allocatable,intent(out)::global_values(:,:)
    complex(real64),allocatable::contribution(:,:)
    integer::p,row
    allocate(contribution(size(local_values,1),global_ngrid),global_values(size(local_values,1),global_ngrid))
    contribution=(0d0,0d0)
    do p=1,size(row_ids)
      row=int(row_ids(p));contribution(:,row)=local_values(:,p)
    enddo
    call MPI_Allreduce(contribution,global_values,size(contribution),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(ierr==MPI_SUCCESS,'distributed row gather failed')
  end subroutine gather_state_rows

  subroutine verify_sparse_catalog_values(bases,full_union)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    complex(real64),intent(in)::full_union(:,:)
    integer::b,j,row,position,slot,local_seen(nunion),global_seen(nunion),local_bad,global_bad
    integer(int64)::expected_bits,actual_bits
    real(real64)::basis_omitted_squared,local_total_squared,global_total_squared,&
      local_max_omitted,global_max_omitted
    local_seen=0;local_bad=0;local_total_squared=0d0;local_max_omitted=0d0
    do b=1,size(bases)
      do j=1,size(bases(b)%global_ids)
        slot=canonical_slot(bases(b)%global_ids(j))
        if(slot==0)then;local_bad=1;cycle;endif
        local_seen(slot)=local_seen(slot)+1;basis_omitted_squared=0d0
        do row=1,global_ngrid
          position=findloc(bases(b)%buffer_point_ids,int(row,int64),dim=1)
          if(position==0)then
            basis_omitted_squared=basis_omitted_squared+global_weights(row)*abs(full_union(slot,row))**2
          else
            expected_bits=transfer(real(full_union(slot,row),real64),expected_bits)
            actual_bits=transfer(real(bases(b)%buffer_values(position,j),real64),actual_bits)
            if(expected_bits/=actual_bits)local_bad=1
            expected_bits=transfer(aimag(full_union(slot,row)),expected_bits)
            actual_bits=transfer(aimag(bases(b)%buffer_values(position,j)),actual_bits)
            if(expected_bits/=actual_bits)local_bad=1
          endif
        enddo
        local_total_squared=local_total_squared+basis_omitted_squared
        local_max_omitted=max(local_max_omitted,sqrt(basis_omitted_squared))
      enddo
    enddo
    call MPI_Allreduce(local_seen,global_seen,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_total_squared,global_total_squared,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_max_omitted,global_max_omitted,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call require(ierr==MPI_SUCCESS.and.global_bad==0.and.all(global_seen==1),&
      'present sparse catalog rows do not bitwise match their canonical Hybrid union values')
    call require(ieee_is_finite(global_max_omitted).and.ieee_is_finite(global_total_squared).and.&
      global_max_omitted<=tail_tolerance.and.sqrt(global_total_squared)<=tail_tolerance,&
      'sparse catalog omitted weighted norm exceeds the basis or aggregate tail tolerance')
  end subroutine verify_sparse_catalog_values

  subroutine assert_catalogs_are_sparse(bases)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    integer::b
    logical::local_sparse
    local_sparse=.true.
    do b=1,size(bases)
      if(bases(b)%fragment_id>0)local_sparse=local_sparse.and.&
        size(bases(b)%buffer_point_ids)>0.and.size(bases(b)%buffer_point_ids)<global_ngrid
    enddo
    call require(local_sparse,'fixture fragment publisher replicated every global spatial row')
  end subroutine assert_catalogs_are_sparse

  integer function canonical_slot(global_id)result(slot)
    integer(int64),intent(in)::global_id
    integer::i
    slot=0
    do i=1,nunion
      if(canonical_union_basis_ids(i)==global_id)then;slot=i;return;endif
    enddo
  end function canonical_slot

  subroutine verify_catalog_tuples(catalog)
    type(s_dg_hybrid_dual_basis_catalog),intent(in)::catalog
    integer::b,j,slot,local_count(nunion),global_count(nunion),local_owner(nunion),global_owner(nunion),&
      local_fragment(nunion),global_fragment(nunion),local_slot(nunion),global_slot(nunion),&
      local_sector(nunion),global_sector(nunion),local_generation(nunion),global_generation(nunion)
    logical::local_ok
    local_count=0;local_owner=0;local_fragment=0;local_slot=0;local_sector=0;local_generation=0;local_ok=.true.
    if(.not.allocated(catalog%fragment_bases))local_ok=.false.
    if(local_ok)then
      do b=1,size(catalog%fragment_bases)
        do j=1,size(catalog%fragment_bases(b)%global_ids)
          slot=canonical_slot(catalog%fragment_bases(b)%global_ids(j))
          if(slot==0)then;local_ok=.false.;cycle;endif
          local_count(slot)=local_count(slot)+1;local_owner(slot)=rank+1
          local_fragment(slot)=catalog%fragment_bases(b)%fragment_id
          local_slot(slot)=j;local_sector(slot)=catalog%fragment_bases(b)%sector(j)
          local_generation(slot)=catalog%fragment_bases(b)%generation
        enddo
      enddo
    endif
    call require(local_ok,'returned catalog tuple contains an unknown basis ID')
    call MPI_Allreduce(local_count,global_count,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_owner,global_owner,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_fragment,global_fragment,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_slot,global_slot,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_sector,global_sector,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_generation,global_generation,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call require(all(global_count==1),'returned sparse catalog does not give unique canonical basis ownership')
    call require(allocated(catalog%uncompressed_global_basis_ids),&
      'dual catalog omitted canonical uncompressed global basis IDs')
    call require(allocated(catalog%uncompressed_owner_ranks).and.allocated(catalog%uncompressed_fragment_ids).and.&
      allocated(catalog%uncompressed_local_slots).and.allocated(catalog%uncompressed_sectors),&
      'dual catalog omitted canonical ownership tuple fields')
    call require(allocated(catalog%uncompressed_generations),&
      'dual catalog omitted fragment generation from canonical ownership tuples')
    call require(size(catalog%uncompressed_global_basis_ids)==nunion.and.&
      all(catalog%uncompressed_global_basis_ids==canonical_union_basis_ids),&
      'dual catalog basis IDs are not sorted canonical noncontiguous IDs')
    call require(all(catalog%uncompressed_owner_ranks==global_owner-1).and.&
      all(catalog%uncompressed_fragment_ids==global_fragment).and.&
      all(catalog%uncompressed_local_slots==global_slot).and.&
      all(catalog%uncompressed_sectors==global_sector).and.&
      all(catalog%uncompressed_generations==global_generation),&
      'dual catalog owner/fragment/generation/local-slot/sector tuples were not derived from publishers')
  end subroutine verify_catalog_tuples

  subroutine derive_catalog_ownership(bases,owners)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    integer,allocatable,intent(out)::owners(:)
    integer::local_owner(nunion),global_owner(nunion),local_count(nunion),global_count(nunion),b,j,slot
    local_owner=0;local_count=0
    do b=1,size(bases)
      do j=1,size(bases(b)%global_ids)
        slot=canonical_slot(bases(b)%global_ids(j));if(slot==0)cycle
        local_owner(slot)=rank+1;local_count(slot)=local_count(slot)+1
      enddo
    enddo
    call MPI_Allreduce(local_owner,global_owner,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_count,global_count,nunion,MPI_INTEGER,MPI_SUM,comm,ierr)
    call require(all(global_count==1),'catalog ownership cannot be uniquely derived')
    allocate(owners(nunion));owners=global_owner-1
  end subroutine derive_catalog_ownership

  logical function catalog_metadata_equal(left,right)result(equal)
    type(s_dg_hybrid_fragment_basis),intent(in)::left(:),right(:)
    integer::b
    equal=size(left)==size(right);if(.not.equal)return
    do b=1,size(left)
      equal=left(b)%fragment_id==right(b)%fragment_id.and.left(b)%generation==right(b)%generation.and.&
        left(b)%provenance_fingerprint==right(b)%provenance_fingerprint
      if(.not.equal)return
      equal=allocated(left(b)%global_ids).eqv.allocated(right(b)%global_ids);if(.not.equal)return
      equal=allocated(left(b)%sector).eqv.allocated(right(b)%sector);if(.not.equal)return
      equal=allocated(left(b)%buffer_point_ids).eqv.allocated(right(b)%buffer_point_ids);if(.not.equal)return
      equal=allocated(left(b)%buffer_values).eqv.allocated(right(b)%buffer_values);if(.not.equal)return
      if(allocated(left(b)%global_ids))then
        if(any(shape(left(b)%global_ids)/=shape(right(b)%global_ids)))then;equal=.false.;return;endif
        if(any(left(b)%global_ids/=right(b)%global_ids))then;equal=.false.;return;endif
      endif
      if(allocated(left(b)%sector))then
        if(any(shape(left(b)%sector)/=shape(right(b)%sector)))then;equal=.false.;return;endif
        if(any(left(b)%sector/=right(b)%sector))then;equal=.false.;return;endif
      endif
      if(allocated(left(b)%buffer_point_ids))then
        if(any(shape(left(b)%buffer_point_ids)/=shape(right(b)%buffer_point_ids)))then;equal=.false.;return;endif
        if(any(left(b)%buffer_point_ids/=right(b)%buffer_point_ids))then;equal=.false.;return;endif
      endif
      if(allocated(left(b)%buffer_values))then
        if(any(shape(left(b)%buffer_values)/=shape(right(b)%buffer_values)))then;equal=.false.;return;endif
        if(.not.bitwise_complex_equal(left(b)%buffer_values,right(b)%buffer_values))then;equal=.false.;return;endif
      endif
    enddo
  end function catalog_metadata_equal

  logical function bitwise_complex_equal(left,right)result(equal)
    complex(real64),intent(in)::left(:,:),right(:,:)
    integer::i,j
    integer(int64)::left_bits,right_bits
    equal=all(shape(left)==shape(right));if(.not.equal)return
    do j=1,size(left,2);do i=1,size(left,1)
      left_bits=transfer(real(left(i,j),real64),left_bits);right_bits=transfer(real(right(i,j),real64),right_bits)
      if(left_bits/=right_bits)then;equal=.false.;return;endif
      left_bits=transfer(aimag(left(i,j)),left_bits);right_bits=transfer(aimag(right(i,j)),right_bits)
      if(left_bits/=right_bits)then;equal=.false.;return;endif
    enddo;enddo
  end function bitwise_complex_equal

  logical function identity_transform_bitwise(transform)result(equal)
    complex(real64),intent(in)::transform(:,:)
    complex(real64),allocatable::identity(:,:)
    integer::i
    if(size(transform,1)/=size(transform,2))then;equal=.false.;return;endif
    allocate(identity(size(transform,1),size(transform,2)));identity=(0d0,0d0)
    do i=1,size(identity,1);identity(i,i)=(1d0,0d0);enddo
    equal=bitwise_complex_equal(transform,identity)
  end function identity_transform_bitwise

  subroutine weighted_gram(values,row_weights,gram)
    complex(real64),intent(in)::values(:,:)
    real(real64),intent(in)::row_weights(:)
    complex(real64),intent(out)::gram(:,:)
    integer::i,j
    gram=(0d0,0d0)
    do j=1,size(values,1);do i=1,size(values,1)
      gram(i,j)=sum(row_weights*conjg(values(i,:))*values(j,:))
    enddo;enddo
  end subroutine weighted_gram

  subroutine hermitian_pseudoinverse(matrix,tolerance,inverse,retained_rank,eigenvalues)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::inverse(:,:)
    integer,intent(out)::retained_rank
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable::vectors(:,:),work(:)
    real(real64),allocatable::rwork(:)
    real(real64)::threshold
    integer::n,info,i,j,k
    external zheev
    n=size(matrix,1);allocate(vectors(n,n),eigenvalues(n),work(max(1,2*n-1)),rwork(max(1,3*n-2)),inverse(n,n))
    vectors=matrix
    call zheev('V','U',n,vectors,n,eigenvalues,work,size(work),rwork,info)
    call require(info==0,'fixture Hermitian eigensolver failed')
    threshold=tolerance*max(1d0,maxval(abs(eigenvalues)));retained_rank=count(eigenvalues>threshold)
    inverse=(0d0,0d0)
    do k=1,n
      if(eigenvalues(k)<=threshold)cycle
      do j=1,n;do i=1,n
        inverse(i,j)=inverse(i,j)+vectors(i,k)*conjg(vectors(j,k))/eigenvalues(k)
      enddo;enddo
    enddo
  end subroutine hermitian_pseudoinverse

  subroutine raw_projection_oracle(wannier,raw_pw,row_weights,tolerance,coefficients,projected,wannier_projector,&
      retained_rank,eigenvalues)
    complex(real64),intent(in)::wannier(:,:),raw_pw(:,:)
    real(real64),intent(in)::row_weights(:),tolerance
    complex(real64),allocatable,intent(out)::coefficients(:,:),projected(:,:),wannier_projector(:,:)
    integer,intent(out)::retained_rank
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable::gram(:,:),gram_inverse(:,:),cross(:,:),weighted_wannier(:,:)
    integer::i,j
    allocate(gram(size(wannier,1),size(wannier,1)));call weighted_gram(wannier,row_weights,gram)
    call hermitian_pseudoinverse(gram,tolerance,gram_inverse,retained_rank,eigenvalues)
    allocate(cross(size(wannier,1),size(raw_pw,1)))
    do j=1,size(raw_pw,1);do i=1,size(wannier,1)
      cross(i,j)=sum(row_weights*conjg(wannier(i,:))*raw_pw(j,:))
    enddo;enddo
    coefficients=matmul(gram_inverse,cross)
    allocate(projected(size(raw_pw,1),size(raw_pw,2)))
    projected=raw_pw-matmul(transpose(coefficients),wannier)
    allocate(weighted_wannier(size(wannier,1),size(wannier,2)))
    do i=1,size(wannier,2);weighted_wannier(:,i)=sqrt(row_weights(i))*wannier(:,i);enddo
    allocate(wannier_projector(size(wannier,2),size(wannier,2)))
    wannier_projector=matmul(transpose(weighted_wannier),matmul(gram_inverse,conjg(weighted_wannier)))
  end subroutine raw_projection_oracle

  subroutine build_weighted_mp_projector(values,row_weights,tolerance,projector,retained_rank,eigenvalues)
    complex(real64),intent(in)::values(:,:)
    real(real64),intent(in)::row_weights(:),tolerance
    complex(real64),allocatable,intent(out)::projector(:,:)
    integer,intent(out)::retained_rank
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable::gram(:,:),gram_inverse(:,:),weighted_values(:,:)
    integer::p
    allocate(gram(size(values,1),size(values,1)));call weighted_gram(values,row_weights,gram)
    call hermitian_pseudoinverse(gram,tolerance,gram_inverse,retained_rank,eigenvalues)
    allocate(weighted_values(size(values,1),size(values,2)))
    do p=1,size(values,2);weighted_values(:,p)=sqrt(row_weights(p))*values(:,p);enddo
    allocate(projector(size(values,2),size(values,2)))
    projector=matmul(transpose(weighted_values),matmul(gram_inverse,conjg(weighted_values)))
  end subroutine build_weighted_mp_projector

  subroutine check_generalized_projection_residual(wannier,projected,reported_defect)
    complex(real64),intent(in)::wannier(:,:),projected(:,:)
    real(real64),intent(in)::reported_defect
    complex(real64),allocatable::local_cross(:,:),global_cross(:,:)
    integer::i,j
    real(real64)::defect
    allocate(local_cross(size(wannier,1),size(projected,1)),global_cross(size(wannier,1),size(projected,1)))
    do j=1,size(projected,1);do i=1,size(wannier,1)
      local_cross(i,j)=sum(weights*conjg(wannier(i,:))*projected(j,:))
    enddo;enddo
    call MPI_Allreduce(local_cross,global_cross,size(local_cross),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    defect=maxval(abs(global_cross))
    call require(defect<comparison_tolerance.and.abs(defect-reported_defect)<comparison_tolerance,&
      'reported generalized orthogonality defect disagrees with direct weighted residual')
  end subroutine check_generalized_projection_residual

  logical function has_degenerate_retained_cluster(eigenvalues,tolerance)result(found)
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer::i,j
    real(real64)::threshold
    found=.false.;threshold=tolerance*max(1d0,maxval(abs(eigenvalues)))
    do j=2,size(eigenvalues);do i=1,j-1
      if(eigenvalues(i)>threshold.and.eigenvalues(j)>threshold.and.&
          abs(eigenvalues(i)-eigenvalues(j))<1d-12)found=.true.
    enddo;enddo
  end function has_degenerate_retained_cluster

  subroutine metric_composed_seed_coefficients(uncompressed,transform,seed_coefficients,row_weights,&
      complete_seed_coefficients)
    complex(real64),intent(in)::uncompressed(:,:),transform(:,:),seed_coefficients(:,:)
    real(real64),intent(in)::row_weights(:)
    complex(real64),allocatable,intent(out)::complete_seed_coefficients(:,:)
    complex(real64),allocatable::gram(:,:),complete_gram(:,:),complete_inverse(:,:),right_hand_side(:,:)
    real(real64),allocatable::eigenvalues(:)
    integer::retained_rank
    allocate(gram(size(uncompressed,1),size(uncompressed,1)));call weighted_gram(uncompressed,row_weights,gram)
    complete_gram=matmul(conjg(transpose(transform)),matmul(gram,transform))
    call hermitian_pseudoinverse(complete_gram,metric_tolerance,complete_inverse,retained_rank,eigenvalues)
    right_hand_side=matmul(conjg(transpose(transform)),matmul(gram,seed_coefficients))
    complete_seed_coefficients=matmul(complete_inverse,right_hand_side)
  end subroutine metric_composed_seed_coefficients

  real(real64) function weighted_state_defect(left,right,row_weights)result(defect)
    complex(real64),intent(in)::left(:,:),right(:,:)
    real(real64),intent(in)::row_weights(:)
    integer::state
    defect=0d0
    do state=1,size(left,1)
      defect=max(defect,sqrt(sum(row_weights*abs(left(state,:)-right(state,:))**2)))
    enddo
  end function weighted_state_defect

  subroutine solve_generalized(hamiltonian,overlap,eigenvalues,eigenvectors,residual,orthogonality)
    complex(real64),intent(in)::hamiltonian(:,:),overlap(:,:)
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable,intent(out)::eigenvectors(:,:)
    real(real64),intent(out)::residual,orthogonality
    complex(real64),allocatable::metric_vectors(:,:),metric_work(:),whitener(:,:),reduced_hamiltonian(:,:),&
      work(:),identity(:,:),residual_matrix(:,:)
    real(real64),allocatable::metric_eigenvalues(:),rwork(:)
    integer::n,info,i,j
    external zheev
    n=size(hamiltonian,1)
    allocate(metric_vectors(n,n),metric_eigenvalues(n),metric_work(max(1,2*n-1)),rwork(max(1,3*n-2)))
    metric_vectors=overlap
    call zheev('V','U',n,metric_vectors,n,metric_eigenvalues,metric_work,size(metric_work),rwork,info)
    call require(info==0.and.minval(metric_eigenvalues)>metric_tolerance,&
      'fixture generalized solver received a singular retained metric')
    allocate(whitener(n,n));whitener=metric_vectors
    do j=1,n;whitener(:,j)=whitener(:,j)/sqrt(metric_eigenvalues(j));enddo
    reduced_hamiltonian=matmul(conjg(transpose(whitener)),matmul(hamiltonian,whitener))
    allocate(eigenvalues(n),work(max(1,2*n-1)))
    call zheev('V','U',n,reduced_hamiltonian,n,eigenvalues,work,size(work),rwork,info)
    call require(info==0,'fixture reduced generalized eigensolver failed')
    allocate(eigenvectors(n,n));eigenvectors=matmul(whitener,reduced_hamiltonian)
    allocate(residual_matrix(n,n));residual_matrix=matmul(hamiltonian,eigenvectors)
    do j=1,n;residual_matrix(:,j)=residual_matrix(:,j)-eigenvalues(j)*matmul(overlap,eigenvectors(:,j));enddo
    residual=maxval(abs(residual_matrix))
    allocate(identity(n,n));identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;enddo
    orthogonality=maxval(abs(matmul(conjg(transpose(eigenvectors)),matmul(overlap,eigenvectors))-identity))
  end subroutine solve_generalized

  function diagonal(matrix)result(values)
    complex(real64),intent(in)::matrix(:,:)
    complex(real64)::values(min(size(matrix,1),size(matrix,2)))
    integer::i
    do i=1,size(values);values(i)=matrix(i,i);enddo
  end function diagonal

  function diagonal_matrix(values)result(matrix)
    real(real64),intent(in)::values(:)
    complex(real64)::matrix(size(values),size(values))
    integer::i
    matrix=(0d0,0d0);do i=1,size(values);matrix(i,i)=values(i);enddo
  end function diagonal_matrix

  real(real64) function retained_condition(eigenvalues,tolerance)result(condition)
    real(real64),intent(in)::eigenvalues(:),tolerance
    real(real64)::threshold,minimum_retained
    integer::i
    threshold=tolerance*max(1d0,maxval(abs(eigenvalues)));minimum_retained=huge(1d0)
    do i=1,size(eigenvalues)
      if(eigenvalues(i)>threshold)minimum_retained=min(minimum_retained,eigenvalues(i))
    enddo
    condition=maxval(eigenvalues)/minimum_retained
  end function retained_condition

  real(real64) function relative_scalar_defect(left,right)result(defect)
    real(real64),intent(in)::left,right
    defect=abs(left-right)/max(1d0,abs(left),abs(right))
  end function relative_scalar_defect

  subroutine check_terminal_transform(full_union,transform)
    complex(real64),intent(in)::full_union(:,:),transform(:,:)
    complex(real64),allocatable::gram(:,:),gram_inverse(:,:),expected_projector(:,:),actual_projector(:,:),&
      column_gram(:,:),identity(:,:)
    real(real64),allocatable::eigenvalues(:)
    integer::retained_rank,i
    allocate(gram(size(full_union,1),size(full_union,1)))
    call weighted_gram(full_union,global_weights,gram)
    call hermitian_pseudoinverse(gram,metric_tolerance,gram_inverse,retained_rank,eigenvalues)
    expected_projector=matmul(gram,gram_inverse)
    actual_projector=matmul(transform,conjg(transpose(transform)))
    column_gram=matmul(conjg(transpose(transform)),transform)
    allocate(identity(size(transform,2),size(transform,2)));identity=(0d0,0d0)
    do i=1,size(identity,1);identity(i,i)=1d0;enddo
    call require(size(transform,2)==retained_rank.and.maxval(abs(column_gram-identity))<comparison_tolerance,&
      'rank-deficient terminal transform columns are not Euclidean orthonormal')
    call require(maxval(abs(actual_projector-expected_projector))<comparison_tolerance,&
      'terminal transform does not span the independent retained coefficient-space projector')
  end subroutine check_terminal_transform

  logical function empty_complex(values)result(empty)
    complex(real64),allocatable,intent(in)::values(:,:)
    empty=.not.allocated(values)
    if(allocated(values))empty=size(values)==0
  end function empty_complex

  logical function empty_integer(values)result(empty)
    integer,allocatable,intent(in)::values(:)
    empty=.not.allocated(values)
    if(allocated(values))empty=size(values)==0
  end function empty_integer

  subroutine expect_projection_unpublished(coefficients,projected,metric_rank,condition,defect,workspace,&
      fingerprint,label)
    complex(real64),allocatable,intent(in)::coefficients(:,:),projected(:,:)
    integer,intent(in)::metric_rank
    real(real64),intent(in)::condition,defect
    integer(int64),intent(in)::workspace,fingerprint
    character(*),intent(in)::label
    call require(empty_complex(coefficients).and.empty_complex(projected).and.metric_rank==0.and.condition==0d0.and.&
      defect==0d0.and.workspace==0_int64.and.fingerprint==0_int64,label)
  end subroutine expect_projection_unpublished

  subroutine expect_map_unpublished(transform,complete,metric_rank,condition,fingerprint,label)
    complex(real64),allocatable,intent(in)::transform(:,:),complete(:,:)
    integer,intent(in)::metric_rank
    real(real64),intent(in)::condition
    integer(int64),intent(in)::fingerprint
    character(*),intent(in)::label
    call require(empty_complex(transform).and.empty_complex(complete).and.metric_rank==0.and.condition==0d0.and.&
      fingerprint==0_int64,label)
  end subroutine expect_map_unpublished

  subroutine expect_catalog_failure(catalog,preserved_owner,preserved_seed_coefficients,&
      complete_seed_coefficients,seed_defect,label)
    type(s_dg_hybrid_dual_basis_catalog),intent(in)::catalog
    integer,allocatable,intent(in)::preserved_owner(:)
    complex(real64),allocatable,intent(in)::preserved_seed_coefficients(:,:),complete_seed_coefficients(:,:)
    real(real64),intent(in)::seed_defect
    character(*),intent(in)::label
    logical::empty_fragments,empty_map,empty_ids,empty_tuples
    empty_fragments=.not.allocated(catalog%fragment_bases)
    if(allocated(catalog%fragment_bases))empty_fragments=size(catalog%fragment_bases)==0
    empty_map=empty_complex(catalog%union_to_complete)
    empty_ids=.not.allocated(catalog%uncompressed_global_basis_ids)
    if(allocated(catalog%uncompressed_global_basis_ids))empty_ids=size(catalog%uncompressed_global_basis_ids)==0
    empty_tuples=empty_integer(catalog%uncompressed_owner_ranks).and.&
      empty_integer(catalog%uncompressed_fragment_ids).and.empty_integer(catalog%uncompressed_local_slots).and.&
      empty_integer(catalog%uncompressed_sectors).and.empty_integer(catalog%uncompressed_generations)
    call require(.not.catalog%valid.and.catalog%uncompressed_rank==0.and.catalog%complete_rank==0.and.&
      catalog%fragment_catalog_fingerprint==0_int64.and.catalog%complete_map_fingerprint==0_int64.and.&
      catalog%complete_transform_binding_fingerprint==0_int64.and.&
      empty_fragments.and.empty_map.and.empty_ids.and.empty_tuples.and.&
      empty_integer(preserved_owner).and.&
      empty_complex(preserved_seed_coefficients).and.empty_complex(complete_seed_coefficients).and.seed_defect==0d0,label)
  end subroutine expect_catalog_failure

  logical function catalog_outputs_unpublished(catalog,preserved_owner,preserved_seed_coefficients,&
      complete_seed_coefficients,seed_defect)result(unpublished)
    type(s_dg_hybrid_dual_basis_catalog),intent(in)::catalog
    integer,allocatable,intent(in)::preserved_owner(:)
    complex(real64),allocatable,intent(in)::preserved_seed_coefficients(:,:),complete_seed_coefficients(:,:)
    real(real64),intent(in)::seed_defect
    logical::empty_fragments,empty_map,empty_ids,empty_tuples
    empty_fragments=.not.allocated(catalog%fragment_bases)
    if(allocated(catalog%fragment_bases))empty_fragments=size(catalog%fragment_bases)==0
    empty_map=empty_complex(catalog%union_to_complete)
    empty_ids=.not.allocated(catalog%uncompressed_global_basis_ids)
    if(allocated(catalog%uncompressed_global_basis_ids))empty_ids=size(catalog%uncompressed_global_basis_ids)==0
    empty_tuples=empty_integer(catalog%uncompressed_owner_ranks).and.&
      empty_integer(catalog%uncompressed_fragment_ids).and.empty_integer(catalog%uncompressed_local_slots).and.&
      empty_integer(catalog%uncompressed_sectors).and.empty_integer(catalog%uncompressed_generations)
    unpublished=.not.catalog%valid.and.catalog%uncompressed_rank==0.and.catalog%complete_rank==0.and.&
      catalog%fragment_catalog_fingerprint==0_int64.and.catalog%complete_map_fingerprint==0_int64.and.&
      catalog%complete_transform_binding_fingerprint==0_int64.and.&
      empty_fragments.and.empty_map.and.empty_ids.and.empty_tuples.and.empty_integer(preserved_owner).and.&
      empty_complex(preserved_seed_coefficients).and.empty_complex(complete_seed_coefficients).and.seed_defect==0d0
  end function catalog_outputs_unpublished

  subroutine expect_failure(ok,message,cause,label)
    logical,intent(in)::ok
    character(*),intent(in)::message,cause,label
    call require(.not.ok.and.index(lowercase(message),trim(lowercase(cause)))>0,&
      trim(label)//'; diagnostic='//trim(message))
  end subroutine expect_failure

  subroutine expect_failure_any(ok,message,cause1,cause2,label)
    logical,intent(in)::ok
    character(*),intent(in)::message,cause1,cause2,label
    call require(.not.ok.and.(index(lowercase(message),trim(lowercase(cause1)))>0.or.&
      index(lowercase(message),trim(lowercase(cause2)))>0),trim(label)//'; diagnostic='//trim(message))
  end subroutine expect_failure_any

  pure function lowercase(text)result(lowered)
    character(*),intent(in)::text
    character(len(text))::lowered
    integer::i,code
    lowered=text
    do i=1,len(text)
      code=iachar(text(i:i));if(code>=iachar('A').and.code<=iachar('Z'))lowered(i:i)=achar(code+32)
    enddo
  end function lowercase

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require

end program test_dg_hybrid_fragment_wannier_lcfo_mpi
