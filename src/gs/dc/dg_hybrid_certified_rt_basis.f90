#include "config.h"
module dg_hybrid_certified_rt_basis
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public::s_dg_hybrid_certified_rt_basis
    logical::valid=.false.
    logical::localization_converged=.false.
    logical::localization_symmetry_constrained=.false.
    integer::global_count=0,certified_rank=0,noccupied=0
    integer::localization_iterations=0
    integer(int64),allocatable::owned_row_ids(:)
    complex(real64),allocatable::c_cert(:,:),u_rt(:,:),b_rt(:,:)
    real(real64),allocatable::certified_eigenvalues(:)
    complex(real64),allocatable::initial_occupied_amplitudes(:,:)
    complex(real64),allocatable::metric_rt(:,:),hamiltonian_rt(:,:)
    complex(real64),allocatable::representation_rt(:,:,:)
    real(real64),allocatable::cartesian_rotations(:,:,:)
    complex(real64),allocatable::scalar_operators_rt(:,:,:)
    complex(real64),allocatable::vector_operators_rt(:,:,:,:)
    complex(real64),allocatable::tensor_operators_rt(:,:,:,:,:)
    real(real64),allocatable::centers(:,:),spreads_before(:),spreads_after(:)
    real(real64)::spread_before_total=0d0,spread_after_total=0d0,spread_improvement=0d0
    real(real64)::transform_unitarity_defect=huge(0d0)
    real(real64)::certified_metric_defect=huge(0d0),rt_metric_defect=huge(0d0)
    real(real64)::embedding_defect=huge(0d0),projector_invariance_defect=huge(0d0)
    real(real64)::target_symmetry_defect_before=huge(0d0)
    real(real64)::target_symmetry_defect_after=huge(0d0)
    real(real64)::energy_symmetry_defect_before=huge(0d0)
    real(real64)::energy_symmetry_defect_after=huge(0d0)
    real(real64)::symmetry_defect_invariance=huge(0d0)
    real(real64)::scalar_covariance_defect=huge(0d0)
    real(real64)::vector_covariance_defect=huge(0d0)
    real(real64)::tensor_covariance_defect=huge(0d0)
    integer(int64)::c_cert_fingerprint=0_int64
    integer(int64)::localization_fingerprint=0_int64
    integer(int64)::b_rt_fingerprint=0_int64
    integer(int64)::operator_fingerprint=0_int64
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_certified_rt_basis

  abstract interface
    subroutine dg_hybrid_certified_localizer(comm,global_count,certified_rank,row_ids,c_cert,&
        require_unconstrained,transform,centers,spreads_before,spreads_after,iterations,&
        symmetry_constrained,converged,ok,message)
      import int64,real64
      integer,intent(in)::comm,global_count,certified_rank
      integer(int64),intent(in)::row_ids(:)
      complex(real64),intent(in)::c_cert(:,:)
      logical,intent(in)::require_unconstrained
      complex(real64),allocatable,intent(out)::transform(:,:)
      real(real64),allocatable,intent(out)::centers(:,:),spreads_before(:),spreads_after(:)
      integer,intent(out)::iterations
      logical,intent(out)::symmetry_constrained,converged,ok
      character(*),intent(out)::message
    end subroutine dg_hybrid_certified_localizer
  end interface

  public::dg_hybrid_certified_localizer
  public::build_dg_hybrid_certified_rt_basis,validate_dg_hybrid_certified_rt_basis
contains
  subroutine build_dg_hybrid_certified_rt_basis(comm,global_count,row_ids,metric_rows,c_cert_rows,&
      certified_eigenvalues,noccupied,certified_representation,cartesian_rotations,scalar_operators,&
      vector_operators,tensor_operators,tolerance,localizer,result,ok,message)
    integer,intent(in)::comm,global_count,noccupied
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::metric_rows(:,:),c_cert_rows(:,:)
    real(real64),intent(in)::certified_eigenvalues(:),cartesian_rotations(:,:,:),tolerance
    complex(real64),intent(in)::certified_representation(:,:,:),scalar_operators(:,:,:)
    complex(real64),intent(in)::vector_operators(:,:,:,:),tensor_operators(:,:,:,:,:)
    procedure(dg_hybrid_certified_localizer)::localizer
    type(s_dg_hybrid_certified_rt_basis),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::transform(:,:)
    real(real64),allocatable::centers(:,:),spreads_before(:),spreads_after(:)
    integer::iterations
    logical::symmetry_constrained,converged,localizer_ok
    character(len(message))::localizer_message
#ifdef USE_MPI
    complex(real64),allocatable::gram(:,:),identity(:,:),energy(:,:),q_c_rows(:,:),q_b_rows(:,:)
    real(real64)::rotation_identity(3,3),rotation_gram(3,3),rotation_defect
    integer,allocatable::ownership(:),owner_codes(:),owner_ranks(:),local_positions(:)
    integer::certified_rank,noperation,nscalar,nvector,ntensor,i,operation,item,a,b,ierr,rank
    integer::ierr_ownership,ierr_owner,ierr_rank,ierr_failure,ierr_failure_agree,&
      ierr_bcast,ierr_bcast_agree
    integer::local_failure_rank,failure_rank,failure_status_bad,bcast_status_bad
    integer::local_bad,global_bad,minimum_integer,maximum_integer
    integer(int64)::minimum_hash,maximum_hash,tolerance_bits,input_hash
    real(real64)::metric_hermiticity_defect,metric_scale,q_c_scale,q_b_scale,scale,&
      spread_before_total,spread_after_total,safe_similarity_limit,safe_embedding_limit
    logical::fingerprint_ok,validation_ok,metric_positive_definite,before_sum_ok,after_sum_ok,&
      distributed_ok
    character(len(message))::validation_message

    result=s_dg_hybrid_certified_rt_basis();ok=.false.;message=''
    certified_rank=size(certified_eigenvalues)
    noperation=size(certified_representation,3)
    nscalar=size(scalar_operators,3);nvector=size(vector_operators,4);ntensor=size(tensor_operators,5)

    call agree_integer(comm,global_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT global size disagrees across ranks';return
    endif
    call agree_integer(comm,certified_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT rank disagrees across ranks';return
    endif
    call agree_integer(comm,noccupied,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT occupied rank disagrees across ranks';return
    endif
    call agree_integer(comm,noperation,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT operation count disagrees across ranks';return
    endif
    call agree_integer(comm,nscalar,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT scalar-operator count disagrees across ranks';return
    endif
    call agree_integer(comm,nvector,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT vector-operator count disagrees across ranks';return
    endif
    call agree_integer(comm,ntensor,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT tensor-operator count disagrees across ranks';return
    endif
    tolerance_bits=transfer(tolerance,0_int64)
    call agree_int64(comm,tolerance_bits,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT tolerance disagrees across ranks';return
    endif

    local_bad=0
    if(global_count<1.or.certified_rank<1.or.certified_rank>global_count.or.&
        certified_rank>ishft(huge(0),-1).or.&
        noccupied<1.or.noccupied>certified_rank.or.noperation<1.or.&
        nscalar<0.or.nvector<0.or.ntensor<0)then
      local_bad=1
    else if(.not.ieee_is_finite(tolerance))then
      local_bad=1
    else if(tolerance<=0d0.or.tolerance>=1d0)then
      local_bad=1
    else if(size(metric_rows,1)/=size(row_ids).or.size(metric_rows,2)/=global_count.or.&
        size(c_cert_rows,1)/=size(row_ids).or.size(c_cert_rows,2)/=certified_rank)then
      local_bad=1
    else if(any(shape(certified_representation)/=[certified_rank,certified_rank,noperation]).or.&
        any(shape(cartesian_rotations)/=[3,3,noperation]).or.&
        any(shape(scalar_operators)/=[certified_rank,certified_rank,nscalar]).or.&
        any(shape(vector_operators)/=[certified_rank,certified_rank,3,nvector]).or.&
        any(shape(tensor_operators)/=[certified_rank,certified_rank,3,3,ntensor]))then
      local_bad=1
    else if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)))then
      local_bad=1
    else if(.not.finite_complex_2(metric_rows).or..not.finite_complex_2(c_cert_rows).or.&
        .not.all(ieee_is_finite(certified_eigenvalues)).or.&
        .not.finite_complex_3(certified_representation).or.&
        .not.all(ieee_is_finite(cartesian_rotations)).or.&
        .not.finite_complex_3(scalar_operators).or..not.finite_complex_4(vector_operators).or.&
        .not.finite_complex_5(tensor_operators))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid certified RT construction contract';return
    endif

    local_bad=0
    safe_similarity_limit=huge(1d0)/(64d0*real(certified_rank,real64)**2)
    safe_embedding_limit=huge(1d0)/(8d0*real(certified_rank,real64)**2)
    if(certified_rank>1)then
      if(any(certified_eigenvalues(2:certified_rank)<certified_eigenvalues(1:certified_rank-1)))&
        local_bad=1
    endif
    if(maxval(abs(certified_eigenvalues))>safe_similarity_limit)local_bad=1
    if(maximum_complex_component_3(certified_representation)>1d0+tolerance)local_bad=1
    if(maximum_complex_component_3(scalar_operators)>safe_similarity_limit)local_bad=1
    if(maximum_complex_component_4(vector_operators)>safe_similarity_limit)local_bad=1
    if(maximum_complex_component_5(tensor_operators)>safe_similarity_limit)local_bad=1
    if(maximum_complex_component_2(c_cert_rows)>safe_embedding_limit)local_bad=1
    if(maxval(abs(cartesian_rotations))>1d0+tolerance)local_bad=1
    rotation_identity=0d0
    do i=1,3;rotation_identity(i,i)=1d0;enddo
    rotation_defect=0d0
    do operation=1,noperation
      rotation_gram=matmul(cartesian_rotations(:,:,operation),&
        transpose(cartesian_rotations(:,:,operation)))
      rotation_defect=max(rotation_defect,&
        sqrt(sum((rotation_gram-rotation_identity)**2))/sqrt(3d0))
    enddo
    if(.not.ieee_is_finite(rotation_defect).or.rotation_defect>tolerance)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='unordered certified spectrum or nonorthogonal Cartesian rotation';return
    endif

    input_hash=fingerprint_input_payload(certified_eigenvalues,certified_representation,&
      cartesian_rotations,scalar_operators,vector_operators,tensor_operators)
    call agree_int64(comm,input_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT replicated input payload disagrees across ranks';return
    endif

    rank=0;call MPI_Comm_rank(comm,rank,ierr_rank)
    local_bad=merge(0,1,ierr_rank==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT row-owner rank query failed';return
    endif
    allocate(ownership(global_count),owner_codes(global_count),owner_ranks(global_count),&
      local_positions(global_count),gram(certified_rank,certified_rank),&
      identity(certified_rank,certified_rank))
    ownership=0;owner_codes=0;local_positions=0
    do i=1,size(row_ids)
      ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1
      owner_codes(int(row_ids(i)))=owner_codes(int(row_ids(i)))+rank+1
      local_positions(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_count,MPI_INTEGER,MPI_SUM,comm,ierr_ownership)
    call MPI_Allreduce(MPI_IN_PLACE,owner_codes,global_count,MPI_INTEGER,MPI_SUM,comm,ierr_owner)
    if(ierr_ownership/=MPI_SUCCESS.or.ierr_owner/=MPI_SUCCESS.or.any(ownership/=1))then
      message='certified RT metric or coefficient rows are incomplete';return
    endif
    owner_ranks=owner_codes-1

    identity=(0d0,0d0)
    do i=1,certified_rank;identity(i,i)=1d0;enddo
    call inspect_row_owned_metric(comm,global_count,row_ids,metric_rows,owner_ranks,&
      local_positions,metric_scale,metric_hermiticity_defect,&
      metric_positive_definite,distributed_ok)
    if(.not.distributed_ok)then
      message='certified RT distributed metric inspection failed';return
    endif
    local_bad=0
    if(.not.ieee_is_finite(metric_hermiticity_defect))then
      local_bad=3
    else if(metric_hermiticity_defect>tolerance)then
      local_bad=3
    else if(.not.metric_positive_definite)then
      local_bad=2
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      result=s_dg_hybrid_certified_rt_basis()
      message='certified RT metric validation reduction failed';return
    else if(global_bad==3)then
      result=s_dg_hybrid_certified_rt_basis()
      message='certified RT construction metric is not Hermitian';return
    else if(global_bad==2)then
      result=s_dg_hybrid_certified_rt_basis()
      message='certified RT construction metric is singular or indefinite';return
    endif
    call build_row_owned_metric_dual_and_gram(comm,global_count,row_ids,metric_rows,&
      c_cert_rows,owner_ranks,local_positions,metric_scale,q_c_rows,q_c_scale,gram,distributed_ok)
    if(.not.distributed_ok)then
      message='certified RT coefficient metric projection failed';return
    endif
    result%certified_metric_defect=frobenius(gram-identity)/sqrt(real(certified_rank,real64))
    local_bad=0
    if(.not.ieee_is_finite(result%certified_metric_defect))then
      local_bad=1
    else if(result%certified_metric_defect>tolerance)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      result=s_dg_hybrid_certified_rt_basis()
      message='certified RT metric validation reduction failed';return
    else if(global_bad/=0)then
      result=s_dg_hybrid_certified_rt_basis()
      message='certified RT coefficients are not metric orthonormal';return
    endif

    localizer_message='';iterations=0;symmetry_constrained=.false.;converged=.false.;localizer_ok=.false.
    call localizer(comm,global_count,certified_rank,row_ids,c_cert_rows,.true.,transform,&
      centers,spreads_before,spreads_after,iterations,symmetry_constrained,converged,&
      localizer_ok,localizer_message)
    local_bad=merge(0,1,localizer_ok)
    if(local_bad==0)then
      if(.not.allocated(transform).or..not.allocated(centers).or.&
          .not.allocated(spreads_before).or..not.allocated(spreads_after))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='certified RT unconstrained localizer failure reduction failed';return
    endif
    if(global_bad/=0)then
      failure_rank=huge(0);ierr_bcast=MPI_SUCCESS
      local_failure_rank=huge(0)
      if(.not.localizer_ok.and.len_trim(localizer_message)>0)local_failure_rank=rank
      call MPI_Allreduce(local_failure_rank,failure_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr_failure)
      local_bad=merge(0,1,ierr_failure==MPI_SUCCESS);failure_status_bad=1
      call MPI_Allreduce(local_bad,failure_status_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_failure_agree)
      if(ierr_failure_agree==MPI_SUCCESS.and.failure_status_bad==0.and.failure_rank<huge(0))then
        if(rank/=failure_rank)localizer_message=''
        call MPI_Bcast(localizer_message,len(localizer_message),MPI_CHARACTER,&
          failure_rank,comm,ierr_bcast)
      endif
      local_bad=merge(0,1,ierr_bcast==MPI_SUCCESS);bcast_status_bad=1
      call MPI_Allreduce(local_bad,bcast_status_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bcast_agree)
      if(ierr_failure_agree==MPI_SUCCESS.and.failure_status_bad==0.and.&
          ierr_bcast_agree==MPI_SUCCESS.and.bcast_status_bad==0)then
        if(failure_rank<huge(0))then
          message='certified RT unconstrained localizer failed: '//trim(localizer_message)
        else
          message='certified RT unconstrained localizer failed'
        endif
      else
        message='certified RT unconstrained localizer failed'
      endif
      return
    endif

    local_bad=0
    if(.not.converged.or.symmetry_constrained.or.iterations<0)then
      local_bad=1
    else if(any(shape(transform)/=[certified_rank,certified_rank]).or.&
        any(shape(centers)/=[3,certified_rank]).or.size(spreads_before)/=certified_rank.or.&
        size(spreads_after)/=certified_rank)then
      local_bad=1
    else if(.not.finite_complex_2(transform).or..not.all(ieee_is_finite(centers)).or.&
        .not.all(ieee_is_finite(spreads_before)).or..not.all(ieee_is_finite(spreads_after)))then
      local_bad=1
    else if(maximum_complex_component_2(transform)>1d0+tolerance)then
      local_bad=1
    else if(any(spreads_before<0d0).or.any(spreads_after<0d0))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fixed-rank unconstrained localization result';return
    endif
    call sum_nonnegative_finite(spreads_before,spread_before_total,before_sum_ok)
    call sum_nonnegative_finite(spreads_after,spread_after_total,after_sum_ok)
    local_bad=merge(0,1,before_sum_ok.and.after_sum_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='localization spread aggregate is not finite';return
    endif
    call agree_integer(comm,iterations,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT localization iterations disagree across ranks';return
    endif
    result%localization_fingerprint=fingerprint_localization(transform,centers,&
      spreads_before,spreads_after,iterations,symmetry_constrained,converged)
    call agree_int64(comm,result%localization_fingerprint,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      result%localization_fingerprint=0_int64
      message='certified RT localization payload disagrees across ranks';return
    endif

    gram=matmul(conjg(transpose(transform)),transform)
    result%transform_unitarity_defect=frobenius(gram-identity)/sqrt(real(certified_rank,real64))
    local_bad=merge(0,1,ieee_is_finite(result%transform_unitarity_defect).and.&
      result%transform_unitarity_defect<=tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT localization transform is not unitary';return
    endif

    result%global_count=global_count;result%certified_rank=certified_rank;result%noccupied=noccupied
    result%localization_iterations=iterations;result%localization_converged=converged
    result%localization_symmetry_constrained=symmetry_constrained
    result%owned_row_ids=row_ids;result%c_cert=c_cert_rows;result%u_rt=transform
    result%b_rt=right_transform_rows(c_cert_rows,transform)
    result%certified_eigenvalues=certified_eigenvalues
    result%centers=centers;result%spreads_before=spreads_before;result%spreads_after=spreads_after
    result%spread_before_total=spread_before_total;result%spread_after_total=spread_after_total
    result%spread_improvement=result%spread_before_total-result%spread_after_total
    allocate(result%initial_occupied_amplitudes(certified_rank,noccupied),&
      result%metric_rt(certified_rank,certified_rank),result%hamiltonian_rt(certified_rank,certified_rank),&
      result%representation_rt(certified_rank,certified_rank,noperation),&
      result%scalar_operators_rt(certified_rank,certified_rank,nscalar),&
      result%vector_operators_rt(certified_rank,certified_rank,3,nvector),&
      result%tensor_operators_rt(certified_rank,certified_rank,3,3,ntensor))
    result%initial_occupied_amplitudes=conjg(transpose(transform(1:noccupied,:)))
    result%metric_rt=identity
    allocate(energy(certified_rank,certified_rank));energy=(0d0,0d0)
    do i=1,certified_rank;energy(i,i)=certified_eigenvalues(i);enddo
    result%hamiltonian_rt=similarity_transform(transform,energy)
    do operation=1,noperation
      result%representation_rt(:,:,operation)=similarity_transform(transform,&
        certified_representation(:,:,operation))
    enddo
    result%cartesian_rotations=cartesian_rotations
    do item=1,nscalar
      result%scalar_operators_rt(:,:,item)=similarity_transform(transform,scalar_operators(:,:,item))
    enddo
    do item=1,nvector;do a=1,3
      result%vector_operators_rt(:,:,a,item)=similarity_transform(transform,vector_operators(:,:,a,item))
    enddo;enddo
    do item=1,ntensor;do b=1,3;do a=1,3
      result%tensor_operators_rt(:,:,a,b,item)=similarity_transform(transform,tensor_operators(:,:,a,b,item))
    enddo;enddo;enddo

    local_bad=0
    if(maximum_complex_component_2(result%hamiltonian_rt)>safe_similarity_limit.or.&
        maximum_complex_component_3(result%representation_rt)>1d0+tolerance.or.&
        maximum_complex_component_3(result%scalar_operators_rt)>safe_similarity_limit.or.&
        maximum_complex_component_4(result%vector_operators_rt)>safe_similarity_limit.or.&
        maximum_complex_component_5(result%tensor_operators_rt)>safe_similarity_limit)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='unsafe transformed certified RT representation or operator payload';return
    endif

    call row_owned_embedding_defect(comm,result%c_cert,result%u_rt,result%b_rt,&
      result%embedding_defect,distributed_ok)
    if(.not.distributed_ok)then;message='certified RT embedding reduction failed';return;endif
    call build_row_owned_metric_dual_and_gram(comm,global_count,row_ids,metric_rows,&
      result%b_rt,owner_ranks,local_positions,metric_scale,q_b_rows,q_b_scale,gram,distributed_ok)
    if(.not.distributed_ok)then;message='certified RT basis metric projection failed';return;endif
    result%rt_metric_defect=frobenius(gram-identity)/sqrt(real(certified_rank,real64))
    call row_owned_projector_defect(comm,global_count,row_ids,result%c_cert,result%b_rt,&
      q_c_rows,q_b_rows,q_c_scale,q_b_scale,owner_ranks,local_positions,&
      result%projector_invariance_defect,distributed_ok)
    if(.not.distributed_ok)then;message='certified RT projector reduction failed';return;endif
    call compute_symmetry_defects(certified_representation,energy,&
      result%target_symmetry_defect_before,result%energy_symmetry_defect_before)
    call compute_symmetry_defects(result%representation_rt,result%hamiltonian_rt,&
      result%target_symmetry_defect_after,result%energy_symmetry_defect_after)
    result%symmetry_defect_invariance=max(&
      abs(result%target_symmetry_defect_after-result%target_symmetry_defect_before),&
      abs(result%energy_symmetry_defect_after-result%energy_symmetry_defect_before))
    local_bad=merge(0,1,all(ieee_is_finite([result%target_symmetry_defect_before,&
      result%target_symmetry_defect_after,result%energy_symmetry_defect_before,&
      result%energy_symmetry_defect_after,result%symmetry_defect_invariance])).and.&
      max(result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
      result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
      result%symmetry_defect_invariance)<=tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT symmetry validation failed before operator covariance';return
    endif
    call compute_covariance_defects(result%representation_rt,result%cartesian_rotations,&
      result%scalar_operators_rt,result%vector_operators_rt,result%tensor_operators_rt,&
      result%scalar_covariance_defect,result%vector_covariance_defect,result%tensor_covariance_defect)

    scale=tolerance
    local_bad=merge(0,1,all(ieee_is_finite([result%spread_before_total,result%spread_after_total,&
      result%spread_improvement,result%rt_metric_defect,result%embedding_defect,&
      result%projector_invariance_defect,result%target_symmetry_defect_before,&
      result%target_symmetry_defect_after,result%energy_symmetry_defect_before,&
      result%energy_symmetry_defect_after,result%symmetry_defect_invariance,&
      result%scalar_covariance_defect,result%vector_covariance_defect,result%tensor_covariance_defect])).and.&
      max(result%rt_metric_defect,result%embedding_defect,result%projector_invariance_defect,&
      result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
      result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
      result%symmetry_defect_invariance,result%scalar_covariance_defect,&
      result%vector_covariance_defect,result%tensor_covariance_defect)<=scale)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT projector, symmetry, or operator covariance validation failed';return
    endif

    call fingerprint_distributed_rows(comm,row_ids,result%c_cert,global_count,&
      result%c_cert_fingerprint,fingerprint_ok)
    if(.not.fingerprint_ok)then;message='certified RT coefficient fingerprint failed';return;endif
    call fingerprint_distributed_rows(comm,row_ids,result%b_rt,global_count,&
      result%b_rt_fingerprint,fingerprint_ok)
    if(.not.fingerprint_ok)then;message='certified RT basis fingerprint failed';return;endif
    result%operator_fingerprint=fingerprint_operator_payload(result%certified_eigenvalues,&
      result%initial_occupied_amplitudes,result%metric_rt,result%hamiltonian_rt,&
      result%representation_rt,result%cartesian_rotations,result%scalar_operators_rt,&
      result%vector_operators_rt,result%tensor_operators_rt)
    call agree_int64(comm,result%operator_fingerprint,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      result%operator_fingerprint=0_int64
      message='certified RT operator payload disagrees across ranks';return
    endif
    result%fingerprint=fingerprint_result_receipt(result)
    if(result%c_cert_fingerprint==0_int64.or.result%localization_fingerprint==0_int64.or.&
        result%b_rt_fingerprint==0_int64.or.result%operator_fingerprint==0_int64.or.&
        result%fingerprint==0_int64)then
      message='certified RT fingerprint construction failed';return
    endif

    result%valid=.true.
    call validate_dg_hybrid_certified_rt_basis(comm,result,tolerance,validation_ok,validation_message)
    if(.not.validation_ok)then
      result%valid=.false.;message='certified RT result self-validation failed: '//trim(validation_message);return
    endif
    ok=.true.;message=''
#else
    result=s_dg_hybrid_certified_rt_basis();ok=.false.
    message='certified RT basis construction requires MPI'
#endif
  end subroutine build_dg_hybrid_certified_rt_basis

  subroutine validate_dg_hybrid_certified_rt_basis(comm,result,tolerance,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_certified_rt_basis),intent(in)::result
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::identity(:,:),energy(:,:),expected_a(:,:),expected_h(:,:),&
      gram(:,:),representation_before(:,:,:)
    integer,allocatable::ownership(:)
    integer::r,noperation,nscalar,nvector,ntensor,i,operation,ierr,local_bad,global_bad
    integer::minimum_integer,maximum_integer
    integer(int64)::minimum_hash,maximum_hash,tolerance_bits,c_hash,u_hash,b_hash,operator_hash,final_hash
    real(real64)::unitarity_defect,embedding_defect,target_before,target_after,energy_before,&
      energy_after,invariance,&
      scalar_defect,vector_defect,tensor_defect,scale,spread_before_total,spread_after_total,&
      expected_improvement,spread_relation_scale,safe_similarity_limit,safe_embedding_limit
    logical::fingerprint_ok,before_sum_ok,after_sum_ok,distributed_ok

    ok=.false.;message='';local_bad=0
    if(.not.ieee_is_finite(tolerance))then
      local_bad=1
    else if(tolerance<=0d0.or.tolerance>=1d0)then
      local_bad=1
    endif
    if(.not.result%valid.or.result%global_count<1.or.result%certified_rank<1.or.&
        result%noccupied<1.or.result%noccupied>result%certified_rank.or.&
        .not.result%localization_converged.or.result%localization_symmetry_constrained)then
      local_bad=1
    endif
    if(.not.allocated(result%owned_row_ids).or..not.allocated(result%c_cert).or.&
        .not.allocated(result%u_rt).or..not.allocated(result%b_rt).or.&
        .not.allocated(result%certified_eigenvalues).or.&
        .not.allocated(result%initial_occupied_amplitudes).or..not.allocated(result%metric_rt).or.&
        .not.allocated(result%hamiltonian_rt).or..not.allocated(result%representation_rt).or.&
        .not.allocated(result%cartesian_rotations).or..not.allocated(result%scalar_operators_rt).or.&
        .not.allocated(result%vector_operators_rt).or..not.allocated(result%tensor_operators_rt).or.&
        .not.allocated(result%centers).or..not.allocated(result%spreads_before).or.&
        .not.allocated(result%spreads_after))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid certified RT result allocation or state';return
    endif

    r=result%certified_rank;noperation=size(result%representation_rt,3)
    nscalar=size(result%scalar_operators_rt,3);nvector=size(result%vector_operators_rt,4)
    ntensor=size(result%tensor_operators_rt,5)
    call agree_integer(comm,result%global_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result global size disagrees across ranks';return
    endif
    call agree_integer(comm,r,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result rank disagrees across ranks';return
    endif
    call agree_integer(comm,result%noccupied,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result occupied rank disagrees across ranks';return
    endif
    call agree_integer(comm,noperation,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result operation count disagrees across ranks';return
    endif
    call agree_integer(comm,nscalar,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result scalar count disagrees across ranks';return
    endif
    call agree_integer(comm,nvector,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result vector count disagrees across ranks';return
    endif
    call agree_integer(comm,ntensor,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='certified RT result tensor count disagrees across ranks';return
    endif
    tolerance_bits=transfer(tolerance,0_int64)
    call agree_int64(comm,tolerance_bits,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT result tolerance disagrees across ranks';return
    endif

    local_bad=0
    safe_similarity_limit=huge(1d0)/(64d0*real(r,real64)**2)
    safe_embedding_limit=huge(1d0)/(8d0*real(r,real64)**2)
    if(r>result%global_count.or.noperation<1.or.nscalar<0.or.nvector<0.or.ntensor<0.or.&
        size(result%c_cert,1)/=size(result%owned_row_ids).or.size(result%c_cert,2)/=r.or.&
        any(shape(result%u_rt)/=[r,r]).or.size(result%b_rt,1)/=size(result%owned_row_ids).or.&
        size(result%b_rt,2)/=r.or.size(result%certified_eigenvalues)/=r.or.&
        any(shape(result%initial_occupied_amplitudes)/=[r,result%noccupied]).or.&
        any(shape(result%metric_rt)/=[r,r]).or.any(shape(result%hamiltonian_rt)/=[r,r]).or.&
        any(shape(result%representation_rt)/=[r,r,noperation]).or.&
        any(shape(result%cartesian_rotations)/=[3,3,noperation]).or.&
        any(shape(result%scalar_operators_rt)/=[r,r,nscalar]).or.&
        any(shape(result%vector_operators_rt)/=[r,r,3,nvector]).or.&
        any(shape(result%tensor_operators_rt)/=[r,r,3,3,ntensor]).or.&
        any(shape(result%centers)/=[3,r]).or.size(result%spreads_before)/=r.or.&
        size(result%spreads_after)/=r)then
      local_bad=1
    else if(any(result%owned_row_ids<1_int64).or.&
        any(result%owned_row_ids>int(result%global_count,int64)).or.&
        .not.finite_complex_2(result%c_cert).or..not.finite_complex_2(result%u_rt).or.&
        .not.finite_complex_2(result%b_rt).or.&
        .not.all(ieee_is_finite(result%certified_eigenvalues)).or.&
        .not.finite_complex_2(result%initial_occupied_amplitudes).or.&
        .not.finite_complex_2(result%metric_rt).or..not.finite_complex_2(result%hamiltonian_rt).or.&
        .not.finite_complex_3(result%representation_rt).or.&
        .not.all(ieee_is_finite(result%cartesian_rotations)).or.&
        .not.finite_complex_3(result%scalar_operators_rt).or.&
        .not.finite_complex_4(result%vector_operators_rt).or.&
        .not.finite_complex_5(result%tensor_operators_rt).or.&
        .not.all(ieee_is_finite(result%centers)).or.&
        .not.all(ieee_is_finite(result%spreads_before)).or.&
        .not.all(ieee_is_finite(result%spreads_after)))then
      local_bad=1
    else if(any(result%spreads_before<0d0).or.any(result%spreads_after<0d0))then
      local_bad=1
    else if(maximum_complex_component_2(result%u_rt)>1d0+tolerance.or.&
        maximum_complex_component_3(result%representation_rt)>1d0+tolerance.or.&
        maxval(abs(result%cartesian_rotations))>1d0+tolerance)then
      local_bad=1
    else if(maxval(abs(result%certified_eigenvalues))>safe_similarity_limit.or.&
        maximum_complex_component_2(result%hamiltonian_rt)>safe_similarity_limit.or.&
        maximum_complex_component_3(result%scalar_operators_rt)>safe_similarity_limit.or.&
        maximum_complex_component_4(result%vector_operators_rt)>safe_similarity_limit.or.&
        maximum_complex_component_5(result%tensor_operators_rt)>safe_similarity_limit.or.&
        maximum_complex_component_2(result%c_cert)>safe_embedding_limit.or.&
        maximum_complex_component_2(result%b_rt)>safe_embedding_limit)then
      local_bad=1
    else if(.not.all(ieee_is_finite([result%spread_before_total,result%spread_after_total,&
        result%spread_improvement,result%transform_unitarity_defect,result%certified_metric_defect,&
        result%rt_metric_defect,result%embedding_defect,result%projector_invariance_defect,&
        result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
        result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
        result%symmetry_defect_invariance,result%scalar_covariance_defect,&
        result%vector_covariance_defect,result%tensor_covariance_defect])))then
      local_bad=1
    else if(result%spread_before_total<0d0.or.result%spread_after_total<0d0)then
      local_bad=1
    else if(any([result%transform_unitarity_defect,result%certified_metric_defect,&
        result%rt_metric_defect,result%embedding_defect,result%projector_invariance_defect,&
        result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
        result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
        result%symmetry_defect_invariance,result%scalar_covariance_defect,&
        result%vector_covariance_defect,result%tensor_covariance_defect]<0d0).or.&
        any([result%transform_unitarity_defect,result%certified_metric_defect,&
        result%rt_metric_defect,result%embedding_defect,result%projector_invariance_defect,&
        result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
        result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
        result%symmetry_defect_invariance,result%scalar_covariance_defect,&
        result%vector_covariance_defect,result%tensor_covariance_defect]>tolerance))then
      local_bad=1
    else if(result%localization_iterations<0.or.result%c_cert_fingerprint==0_int64.or.&
        result%localization_fingerprint==0_int64.or.result%b_rt_fingerprint==0_int64.or.&
        result%operator_fingerprint==0_int64.or.result%fingerprint==0_int64)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid certified RT result payload';return
    endif
    call sum_nonnegative_finite(result%spreads_before,spread_before_total,before_sum_ok)
    call sum_nonnegative_finite(result%spreads_after,spread_after_total,after_sum_ok)
    local_bad=merge(0,1,before_sum_ok.and.after_sum_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT result spread aggregate is not finite';return
    endif

    allocate(ownership(result%global_count));ownership=0
    do i=1,size(result%owned_row_ids)
      ownership(int(result%owned_row_ids(i)))=ownership(int(result%owned_row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,result%global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      message='certified RT result row ownership is incomplete';return
    endif

    allocate(identity(r,r),energy(r,r),expected_a(r,result%noccupied),expected_h(r,r),gram(r,r),&
      representation_before(r,r,noperation))
    identity=(0d0,0d0);energy=(0d0,0d0)
    do i=1,r;identity(i,i)=1d0;energy(i,i)=result%certified_eigenvalues(i);enddo
    gram=matmul(conjg(transpose(result%u_rt)),result%u_rt)
    unitarity_defect=frobenius(gram-identity)/sqrt(real(r,real64))
    expected_a=conjg(transpose(result%u_rt(1:result%noccupied,:)))
    expected_h=similarity_transform(result%u_rt,energy)
    call row_owned_embedding_defect(comm,result%c_cert,result%u_rt,result%b_rt,&
      embedding_defect,distributed_ok)
    if(.not.distributed_ok)then
      message='certified RT embedding reduction failed';return
    endif
    do operation=1,noperation
      representation_before(:,:,operation)=matmul(result%u_rt,&
        matmul(result%representation_rt(:,:,operation),conjg(transpose(result%u_rt))))
    enddo
    call compute_symmetry_defects(representation_before,energy,target_before,energy_before)
    call compute_symmetry_defects(result%representation_rt,result%hamiltonian_rt,target_after,energy_after)
    invariance=max(abs(target_after-target_before),abs(energy_after-energy_before))
    local_bad=merge(0,1,all(ieee_is_finite([target_before,target_after,energy_before,&
      energy_after,invariance])).and.max(target_before,target_after,energy_before,energy_after,&
      invariance)<=tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT result symmetry validation failed before operator covariance';return
    endif
    call compute_covariance_defects(result%representation_rt,result%cartesian_rotations,&
      result%scalar_operators_rt,result%vector_operators_rt,result%tensor_operators_rt,&
      scalar_defect,vector_defect,tensor_defect)
    expected_improvement=spread_before_total-spread_after_total
    spread_relation_scale=max(1d0,abs(spread_before_total),abs(spread_after_total),&
      abs(result%spread_improvement),abs(expected_improvement))
    scale=tolerance
    local_bad=merge(0,1,unitarity_defect<=scale.and.embedding_defect<=scale.and.&
      frobenius(result%metric_rt-identity)/sqrt(real(r,real64))<=scale.and.&
      frobenius(result%hamiltonian_rt-expected_h)/max(1d0,frobenius(expected_h))<=scale.and.&
      frobenius(result%initial_occupied_amplitudes-expected_a)/&
        max(1d0,frobenius(expected_a))<=scale.and.&
      max(target_before,target_after,energy_before,energy_after,invariance,scalar_defect,&
        vector_defect,tensor_defect,result%certified_metric_defect,result%rt_metric_defect,&
        result%projector_invariance_defect)<=scale.and.&
      abs(result%spread_before_total-spread_before_total)<=scale*max(1d0,abs(spread_before_total)).and.&
      abs(result%spread_after_total-spread_after_total)<=scale*max(1d0,abs(spread_after_total)).and.&
      abs(result%spread_improvement/spread_relation_scale-&
        expected_improvement/spread_relation_scale)<=scale.and.&
      abs(result%transform_unitarity_defect-unitarity_defect)<=scale.and.&
      abs(result%embedding_defect-embedding_defect)<=scale.and.&
      abs(result%target_symmetry_defect_before-target_before)<=scale.and.&
      abs(result%target_symmetry_defect_after-target_after)<=scale.and.&
      abs(result%energy_symmetry_defect_before-energy_before)<=scale.and.&
      abs(result%energy_symmetry_defect_after-energy_after)<=scale.and.&
      abs(result%symmetry_defect_invariance-invariance)<=scale.and.&
      abs(result%scalar_covariance_defect-scalar_defect)<=scale.and.&
      abs(result%vector_covariance_defect-vector_defect)<=scale.and.&
      abs(result%tensor_covariance_defect-tensor_defect)<=scale)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT result violates gauge, projector, symmetry, or covariance invariants';return
    endif

    call fingerprint_distributed_rows(comm,result%owned_row_ids,result%c_cert,result%global_count,&
      c_hash,fingerprint_ok)
    if(.not.fingerprint_ok)then;message='certified RT coefficient fingerprint validation failed';return;endif
    call fingerprint_distributed_rows(comm,result%owned_row_ids,result%b_rt,result%global_count,&
      b_hash,fingerprint_ok)
    if(.not.fingerprint_ok)then;message='certified RT basis fingerprint validation failed';return;endif
    u_hash=fingerprint_localization(result%u_rt,result%centers,result%spreads_before,&
      result%spreads_after,result%localization_iterations,result%localization_symmetry_constrained,&
      result%localization_converged)
    operator_hash=fingerprint_operator_payload(result%certified_eigenvalues,&
      result%initial_occupied_amplitudes,result%metric_rt,result%hamiltonian_rt,&
      result%representation_rt,result%cartesian_rotations,result%scalar_operators_rt,&
      result%vector_operators_rt,result%tensor_operators_rt)
    final_hash=fingerprint_result_receipt(result,c_hash,u_hash,b_hash,operator_hash)
    local_bad=merge(0,1,c_hash==result%c_cert_fingerprint.and.&
      u_hash==result%localization_fingerprint.and.b_hash==result%b_rt_fingerprint.and.&
      operator_hash==result%operator_fingerprint.and.final_hash==result%fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified RT result fingerprint mismatch';return
    endif
    call agree_int64(comm,u_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT localization differs across ranks';return
    endif
    call agree_int64(comm,operator_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT operators differ across ranks';return
    endif
    call agree_int64(comm,final_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='certified RT final receipt differs across ranks';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='certified RT basis validation requires MPI'
#endif
  end subroutine validate_dg_hybrid_certified_rt_basis

  pure function similarity_transform(transform,operator)result(transformed)
    complex(real64),intent(in)::transform(:,:),operator(:,:)
    complex(real64)::transformed(size(transform,2),size(transform,2))
    transformed=matmul(conjg(transpose(transform)),matmul(operator,transform))
  end function similarity_transform

  pure function right_transform_rows(rows,transform)result(transformed)
    complex(real64),intent(in)::rows(:,:),transform(:,:)
    complex(real64)::transformed(size(rows,1),size(transform,2))
    integer::local_row,column,k
    transformed=(0d0,0d0)
    do local_row=1,size(rows,1);do column=1,size(transform,2);do k=1,size(rows,2)
      transformed(local_row,column)=transformed(local_row,column)+rows(local_row,k)*transform(k,column)
    enddo;enddo;enddo
  end function right_transform_rows

  pure real(real64) function frobenius(values)result(norm)
    complex(real64),intent(in)::values(:,:)
    real(real64)::scale,sum_squares,component,root
    integer::i,j,part
    scale=0d0;sum_squares=1d0
    do j=1,size(values,2);do i=1,size(values,1);do part=1,2
      if(part==1)then;component=abs(real(values(i,j),real64));else;component=abs(aimag(values(i,j)));endif
      if(.not.ieee_is_finite(component))then;norm=huge(1d0);return;endif
      if(component>0d0)then
        if(scale<component)then
          sum_squares=1d0+sum_squares*(scale/component)**2;scale=component
        else
          sum_squares=sum_squares+(component/scale)**2
        endif
      endif
    enddo;enddo;enddo
    if(scale<=0d0)then;norm=0d0;return;endif
    root=sqrt(sum_squares)
    if(scale>huge(1d0)/root)then;norm=huge(1d0);else;norm=scale*root;endif
  end function frobenius

  subroutine sum_nonnegative_finite(values,total,ok)
    real(real64),intent(in)::values(:)
    real(real64),intent(out)::total
    logical,intent(out)::ok
    real(real64)::scale,scaled_total
    total=0d0;ok=.false.
    if(size(values)<1)return
    if(.not.all(ieee_is_finite(values)))return
    if(any(values<0d0))return
    scale=maxval(values)
    if(scale<=0d0)then;ok=.true.;return;endif
    scaled_total=sum(values/scale)
    if(.not.ieee_is_finite(scaled_total))return
    if(scaled_total<=0d0)return
    if(scale>huge(1d0)/scaled_total)return
    total=scale*scaled_total;ok=ieee_is_finite(total)
  end subroutine sum_nonnegative_finite

#ifdef USE_MPI
  subroutine inspect_row_owned_metric(comm,global_count,row_ids,metric_rows,owner_ranks,&
      local_positions,metric_scale,hermiticity_defect,positive_definite,ok)
    integer,intent(in)::comm,global_count,owner_ranks(:),local_positions(:)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::metric_rows(:,:)
    real(real64),intent(out)::metric_scale,hermiticity_defect
    logical,intent(out)::positive_definite,ok
    complex(real64),allocatable::row_buffer(:),lower_rows(:,:)
    real(real64),allocatable::row_norm_squares(:)
    complex(real64)::difference,remainder
    real(real64)::local_scale,local_squares(2),global_squares(2),pivot,threshold,&
      diagonal_budget,remaining_budget,remainder_bound,entry_square
    integer::rank,i,j,k,local_row,status,failed,local_bad,global_bad
    integer::ierr_rank,ierr_scale,ierr_bcast_status,ierr_bcast_row,ierr_update,ierr_squares,ierr_bad

    metric_scale=0d0;hermiticity_defect=huge(1d0);positive_definite=.false.;ok=.false.
    rank=0;call MPI_Comm_rank(comm,rank,ierr_rank)
    local_bad=merge(0,1,ierr_rank==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return

    local_scale=maximum_complex_component_2(metric_rows)
    call MPI_Allreduce(local_scale,metric_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr_scale)
    if(ierr_scale/=MPI_SUCCESS)return
    if(.not.ieee_is_finite(metric_scale).or.metric_scale<=0d0)then
      hermiticity_defect=0d0;positive_definite=.false.;ok=.true.;return
    endif

    allocate(row_buffer(global_count),lower_rows(size(row_ids),global_count),&
      row_norm_squares(size(row_ids)))
    local_squares=0d0;local_bad=0
    do j=1,global_count
      row_buffer=(0d0,0d0)
      if(rank==owner_ranks(j))then
        local_row=local_positions(j)
        row_buffer=metric_rows(local_row,:)/metric_scale
      endif
      call MPI_Bcast(row_buffer,global_count,MPI_DOUBLE_COMPLEX,owner_ranks(j),comm,ierr_bcast_row)
      if(ierr_bcast_row/=MPI_SUCCESS)local_bad=1
      do local_row=1,size(row_ids)
        i=int(row_ids(local_row))
        remainder=metric_rows(local_row,j)/metric_scale
        difference=remainder-conjg(row_buffer(i))
        local_squares(1)=local_squares(1)+real(difference,real64)**2+aimag(difference)**2
        local_squares(2)=local_squares(2)+real(remainder,real64)**2+aimag(remainder)**2
      enddo
    enddo
    global_squares=0d0
    call MPI_Allreduce(local_squares,global_squares,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr_squares)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_squares/=MPI_SUCCESS.or.ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return
    if(global_squares(2)<=0d0.or..not.all(ieee_is_finite(global_squares)))then
      hermiticity_defect=huge(1d0);positive_definite=.false.;ok=.true.;return
    endif
    hermiticity_defect=saturated_nonnegative_ratio(sqrt(global_squares(1)),sqrt(global_squares(2)))

    lower_rows=metric_rows/metric_scale;row_norm_squares=0d0
    threshold=64d0*epsilon(1d0);failed=0;local_bad=0
    do j=1,global_count
      status=failed;row_buffer=(0d0,0d0)
      if(rank==owner_ranks(j).and.status==0)then
        local_row=local_positions(j)
        pivot=real(lower_rows(local_row,j),real64)-row_norm_squares(local_row)
        if(.not.ieee_is_finite(pivot))then
          status=1
        else if(abs(aimag(lower_rows(local_row,j)))>threshold.or.pivot<=threshold)then
          status=1
        else
          if(j>1)row_buffer(1:j-1)=lower_rows(local_row,1:j-1)
          row_buffer(j)=cmplx(sqrt(pivot),0d0,real64)
          lower_rows(local_row,j)=row_buffer(j)
        endif
      endif
      call MPI_Bcast(status,1,MPI_INTEGER,owner_ranks(j),comm,ierr_bcast_status)
      call MPI_Bcast(row_buffer,global_count,MPI_DOUBLE_COMPLEX,owner_ranks(j),comm,ierr_bcast_row)
      if(ierr_bcast_status/=MPI_SUCCESS.or.ierr_bcast_row/=MPI_SUCCESS)local_bad=1
      global_bad=0
      if(status==0)then
        do local_row=1,size(row_ids)
          if(row_ids(local_row)<=int(j,int64))cycle
          remainder=lower_rows(local_row,j)
          do k=1,j-1
            remainder=remainder-lower_rows(local_row,k)*conjg(row_buffer(k))
          enddo
          i=int(row_ids(local_row));diagonal_budget=real(lower_rows(local_row,i),real64)
          remaining_budget=diagonal_budget-row_norm_squares(local_row)
          if(remaining_budget<=0d0.or..not.ieee_is_finite(remaining_budget))then
            global_bad=1
          else
            remainder_bound=sqrt(remaining_budget)*real(row_buffer(j),real64)
            if(abs(remainder)>remainder_bound)then
              global_bad=1
            else
              remainder=remainder/real(row_buffer(j),real64)
              entry_square=abs(remainder)**2
              lower_rows(local_row,j)=remainder
              row_norm_squares(local_row)=row_norm_squares(local_row)+entry_square
              if(.not.ieee_is_finite(real(remainder,real64)).or.&
                  .not.ieee_is_finite(aimag(remainder)).or.&
                  .not.ieee_is_finite(row_norm_squares(local_row)))global_bad=1
            endif
          endif
        enddo
      endif
      call MPI_Allreduce(MPI_IN_PLACE,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_update)
      if(ierr_update/=MPI_SUCCESS)local_bad=1
      failed=max(failed,status,global_bad)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return
    positive_definite=failed==0;ok=.true.
  end subroutine inspect_row_owned_metric

  subroutine build_row_owned_metric_dual_and_gram(comm,global_count,row_ids,metric_rows,&
      basis_rows,owner_ranks,local_positions,metric_scale,q_rows,q_scale,gram,ok)
    integer,intent(in)::comm,global_count,owner_ranks(:),local_positions(:)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::metric_rows(:,:),basis_rows(:,:)
    real(real64),intent(in)::metric_scale
    complex(real64),allocatable,intent(out)::q_rows(:,:)
    real(real64),intent(out)::q_scale
    complex(real64),intent(out)::gram(:,:)
    logical,intent(out)::ok
    complex(real64),allocatable::source_metric(:),source_basis(:),gram_pack(:)
    real(real64)::local_basis_scale,basis_scale,local_q_scale,normalized_q_scale,&
      metric_basis_scale,gram_scale,gram_component_bound
    integer::rank,r,source,j,a,b,local_row,local_bad,global_bad,ierr_rank,ierr_scale,&
      ierr_metric_bcast,ierr_basis_bcast,ierr_bad,ierr_q_scale,ierr_gram_bcast
    logical::product_ok

    ok=.false.;q_scale=0d0;gram=(0d0,0d0);r=size(basis_rows,2)
    allocate(q_rows(size(row_ids),r),source_metric(global_count),source_basis(r),gram_pack(2*r))
    q_rows=(0d0,0d0);rank=0;call MPI_Comm_rank(comm,rank,ierr_rank)
    local_bad=merge(0,1,ierr_rank==MPI_SUCCESS)
    local_basis_scale=maximum_complex_component_2(basis_rows)
    call MPI_Allreduce(local_basis_scale,basis_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr_scale)
    if(ierr_scale/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.local_bad/=0)return
    if(.not.ieee_is_finite(basis_scale).or.basis_scale<=0d0.or.&
        .not.ieee_is_finite(metric_scale).or.metric_scale<=0d0)return

    ! The dense metric is already streamed in global-row order for the distributed
    ! Cholesky check.  Reuse that deterministic order here so Gram/projector receipts
    ! remain bitwise invariant when the same rows are repartitioned across MPI ranks.
    local_bad=0
    do source=1,global_count
      source_metric=(0d0,0d0);source_basis=(0d0,0d0)
      if(rank==owner_ranks(source))then
        local_row=local_positions(source)
        source_metric=metric_rows(local_row,:)/metric_scale
        source_basis=basis_rows(local_row,:)/basis_scale
      endif
      call MPI_Bcast(source_metric,global_count,MPI_DOUBLE_COMPLEX,owner_ranks(source),comm,&
        ierr_metric_bcast)
      call MPI_Bcast(source_basis,r,MPI_DOUBLE_COMPLEX,owner_ranks(source),comm,ierr_basis_bcast)
      if(ierr_metric_bcast/=MPI_SUCCESS.or.ierr_basis_bcast/=MPI_SUCCESS)local_bad=1
      do local_row=1,size(row_ids)
        j=int(row_ids(local_row))
        do a=1,r
          q_rows(local_row,a)=q_rows(local_row,a)+conjg(source_basis(a))*source_metric(j)
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return

    local_q_scale=maximum_complex_component_2(q_rows);normalized_q_scale=0d0
    call MPI_Allreduce(local_q_scale,normalized_q_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr_q_scale)
    if(ierr_q_scale/=MPI_SUCCESS)return
    if(.not.ieee_is_finite(normalized_q_scale).or.normalized_q_scale<=0d0)return
    q_rows=q_rows/normalized_q_scale
    call safe_nonnegative_product(metric_scale,basis_scale,metric_basis_scale,product_ok)
    if(.not.product_ok)return
    call safe_nonnegative_product(metric_basis_scale,normalized_q_scale,q_scale,product_ok)
    if(.not.product_ok.or.q_scale<=0d0)return

    gram=(0d0,0d0);local_bad=0
    do j=1,global_count
      gram_pack=(0d0,0d0)
      if(rank==owner_ranks(j))then
        local_row=local_positions(j)
        gram_pack(1:r)=q_rows(local_row,:)
        gram_pack(r+1:2*r)=basis_rows(local_row,:)/basis_scale
      endif
      call MPI_Bcast(gram_pack,2*r,MPI_DOUBLE_COMPLEX,owner_ranks(j),comm,ierr_gram_bcast)
      if(ierr_gram_bcast/=MPI_SUCCESS)local_bad=1
      do b=1,r;do a=1,r
        gram(a,b)=gram(a,b)+gram_pack(a)*gram_pack(r+b)
      enddo;enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return
    call safe_nonnegative_product(q_scale,basis_scale,gram_scale,product_ok)
    if(.not.product_ok)return
    gram_component_bound=maximum_complex_component_2(gram)
    call safe_nonnegative_product(gram_scale,gram_component_bound,local_basis_scale,product_ok)
    if(.not.product_ok)return
    gram=gram*gram_scale
    if(.not.finite_complex_2(gram))return
    ok=.true.
  end subroutine build_row_owned_metric_dual_and_gram

  subroutine row_owned_embedding_defect(comm,c_rows,transform,b_rows,defect,ok)
    integer,intent(in)::comm
    complex(real64),intent(in)::c_rows(:,:),transform(:,:),b_rows(:,:)
    real(real64),intent(out)::defect
    logical,intent(out)::ok
    complex(real64)::expected,difference
    real(real64)::local_scales(2),global_scales(2),common_scale,local_squares(2),global_squares(2)
    integer::local_row,a,b,r,ierr_scales,ierr_squares

    defect=huge(1d0);ok=.false.;r=size(transform,1)
    local_scales=[maximum_complex_component_2(c_rows),maximum_complex_component_2(b_rows)]
    global_scales=0d0
    call MPI_Allreduce(local_scales,global_scales,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr_scales)
    if(ierr_scales/=MPI_SUCCESS.or..not.all(ieee_is_finite(global_scales)))return
    common_scale=maxval(global_scales)
    if(common_scale<=0d0)then;defect=0d0;ok=.true.;return;endif
    local_squares=0d0
    do local_row=1,size(c_rows,1);do b=1,r
      expected=(0d0,0d0)
      do a=1,r
        expected=expected+(c_rows(local_row,a)/common_scale)*transform(a,b)
      enddo
      difference=b_rows(local_row,b)/common_scale-expected
      local_squares(1)=local_squares(1)+real(difference,real64)**2+aimag(difference)**2
      difference=b_rows(local_row,b)/common_scale
      local_squares(2)=local_squares(2)+real(difference,real64)**2+aimag(difference)**2
    enddo;enddo
    global_squares=0d0
    call MPI_Allreduce(local_squares,global_squares,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr_squares)
    if(ierr_squares/=MPI_SUCCESS.or..not.all(ieee_is_finite(global_squares)))return
    defect=scaled_norm_ratio(common_scale,global_squares(1),global_squares(2));ok=.true.
  end subroutine row_owned_embedding_defect

  subroutine row_owned_projector_defect(comm,global_count,row_ids,c_rows,b_rows,q_c_rows,q_b_rows,&
      q_c_scale,q_b_scale,owner_ranks,local_positions,defect,ok)
    integer,intent(in)::comm,global_count,owner_ranks(:),local_positions(:)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::c_rows(:,:),b_rows(:,:),q_c_rows(:,:),q_b_rows(:,:)
    real(real64),intent(in)::q_c_scale,q_b_scale
    real(real64),intent(out)::defect
    logical,intent(out)::ok
    complex(real64),allocatable::q_pack(:)
    complex(real64)::projector_c,projector_b,difference
    real(real64)::local_scales(2),global_scales(2),projector_c_scale,projector_b_scale,&
      projector_scale,global_squares(2),row_pair(2)
    real(real64),allocatable::row_squares(:,:)
    integer::rank,r,i,j,a,local_row,local_bad,global_bad,ierr_rank,ierr_scales,ierr_bcast,&
      ierr_bad,ierr_row_bcast
    logical::product_ok

    defect=huge(1d0);ok=.false.;rank=0;r=size(c_rows,2)
    allocate(q_pack(2*r),row_squares(size(row_ids),2));row_squares=0d0
    call MPI_Comm_rank(comm,rank,ierr_rank)
    local_bad=merge(0,1,ierr_rank==MPI_SUCCESS)
    local_scales=[maximum_complex_component_2(c_rows),maximum_complex_component_2(b_rows)]
    global_scales=0d0
    call MPI_Allreduce(local_scales,global_scales,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr_scales)
    if(ierr_scales/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.local_bad/=0)return
    if(any(global_scales<=0d0).or..not.all(ieee_is_finite(global_scales)).or.&
        .not.ieee_is_finite(q_c_scale).or..not.ieee_is_finite(q_b_scale).or.&
        q_c_scale<=0d0.or.q_b_scale<=0d0)return
    call safe_nonnegative_product(global_scales(1),q_c_scale,projector_c_scale,product_ok)
    if(.not.product_ok)return
    call safe_nonnegative_product(global_scales(2),q_b_scale,projector_b_scale,product_ok)
    if(.not.product_ok)return
    projector_scale=max(projector_c_scale,projector_b_scale)
    if(projector_scale<=0d0)return

    local_bad=0
    do j=1,global_count
      q_pack=(0d0,0d0)
      if(rank==owner_ranks(j))then
        q_pack(1:r)=q_c_rows(local_positions(j),:)
        q_pack(r+1:2*r)=q_b_rows(local_positions(j),:)
      endif
      call MPI_Bcast(q_pack,2*r,MPI_DOUBLE_COMPLEX,owner_ranks(j),comm,ierr_bcast)
      if(ierr_bcast/=MPI_SUCCESS)local_bad=1
      do local_row=1,size(row_ids)
        projector_c=(0d0,0d0);projector_b=(0d0,0d0)
        do a=1,r
          projector_c=projector_c+(c_rows(local_row,a)/global_scales(1))*q_pack(a)
          projector_b=projector_b+(b_rows(local_row,a)/global_scales(2))*q_pack(r+a)
        enddo
        projector_c=projector_c*(projector_c_scale/projector_scale)
        projector_b=projector_b*(projector_b_scale/projector_scale)
        difference=projector_b-projector_c
        row_squares(local_row,1)=row_squares(local_row,1)+&
          real(difference,real64)**2+aimag(difference)**2
        row_squares(local_row,2)=row_squares(local_row,2)+&
          real(projector_c,real64)**2+aimag(projector_c)**2
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0)return
    global_squares=0d0;local_bad=0
    do i=1,global_count
      row_pair=0d0
      if(rank==owner_ranks(i))row_pair=row_squares(local_positions(i),:)
      call MPI_Bcast(row_pair,2,MPI_DOUBLE_PRECISION,owner_ranks(i),comm,ierr_row_bcast)
      if(ierr_row_bcast/=MPI_SUCCESS)local_bad=1
      global_squares=global_squares+row_pair
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr_bad)
    if(ierr_bad/=MPI_SUCCESS.or.global_bad/=0.or..not.all(ieee_is_finite(global_squares)))return
    defect=scaled_norm_ratio(projector_scale,global_squares(1),global_squares(2));ok=.true.
  end subroutine row_owned_projector_defect
#endif

  pure subroutine safe_nonnegative_product(left,right,product,ok)
    real(real64),intent(in)::left,right
    real(real64),intent(out)::product
    logical,intent(out)::ok
    product=0d0;ok=.false.
    if(.not.ieee_is_finite(left).or..not.ieee_is_finite(right))return
    if(left<0d0.or.right<0d0)return
    if(left<=0d0.or.right<=0d0)then;ok=.true.;return;endif
    if(right>=1d0)then
      if(left>huge(1d0)/right)return
    endif
    product=left*right;ok=ieee_is_finite(product)
  end subroutine safe_nonnegative_product

  pure real(real64) function saturated_nonnegative_product(left,right)result(product)
    real(real64),intent(in)::left,right
    if(left<=0d0.or.right<=0d0)then
      product=0d0
    else
      if(right>=1d0)then
        if(left>huge(1d0)/right)then
          product=huge(1d0);return
        endif
      endif
      product=left*right
    endif
  end function saturated_nonnegative_product

  pure real(real64) function saturated_nonnegative_ratio(numerator,denominator)result(ratio)
    real(real64),intent(in)::numerator,denominator
    if(numerator<=0d0)then
      ratio=0d0
    else if(denominator<=0d0)then
      ratio=huge(1d0)
    else
      if(denominator<1d0)then
        if(numerator>huge(1d0)*denominator)then
          ratio=huge(1d0);return
        endif
      endif
      ratio=numerator/denominator
    endif
  end function saturated_nonnegative_ratio

  pure real(real64) function scaled_norm_ratio(common_scale,numerator_squares,&
      denominator_squares)result(ratio)
    real(real64),intent(in)::common_scale,numerator_squares,denominator_squares
    real(real64)::numerator_root,denominator_root,denominator_norm
    if(common_scale<=0d0.or.numerator_squares<=0d0)then;ratio=0d0;return;endif
    numerator_root=sqrt(numerator_squares)
    denominator_root=sqrt(max(0d0,denominator_squares))
    denominator_norm=saturated_nonnegative_product(common_scale,denominator_root)
    if(denominator_norm>=1d0.and.denominator_root>0d0)then
      ratio=saturated_nonnegative_ratio(numerator_root,denominator_root)
    else
      ratio=saturated_nonnegative_product(common_scale,numerator_root)
    endif
  end function scaled_norm_ratio

  pure logical function finite_complex_2(values)result(finite)
    complex(real64),intent(in)::values(:,:)
    finite=all(ieee_is_finite(real(values,real64))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_2

  pure logical function finite_complex_3(values)result(finite)
    complex(real64),intent(in)::values(:,:,:)
    finite=all(ieee_is_finite(real(values,real64))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_3

  pure logical function finite_complex_4(values)result(finite)
    complex(real64),intent(in)::values(:,:,:,:)
    finite=all(ieee_is_finite(real(values,real64))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_4

  pure logical function finite_complex_5(values)result(finite)
    complex(real64),intent(in)::values(:,:,:,:,:)
    finite=all(ieee_is_finite(real(values,real64))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_5

  pure real(real64) function maximum_complex_component_2(values)result(maximum_value)
    complex(real64),intent(in)::values(:,:)
    maximum_value=0d0
    if(size(values)>0)maximum_value=max(maxval(abs(real(values,real64))),maxval(abs(aimag(values))))
  end function maximum_complex_component_2

  pure real(real64) function maximum_complex_component_3(values)result(maximum_value)
    complex(real64),intent(in)::values(:,:,:)
    maximum_value=0d0
    if(size(values)>0)maximum_value=max(maxval(abs(real(values,real64))),maxval(abs(aimag(values))))
  end function maximum_complex_component_3

  pure real(real64) function maximum_complex_component_4(values)result(maximum_value)
    complex(real64),intent(in)::values(:,:,:,:)
    maximum_value=0d0
    if(size(values)>0)maximum_value=max(maxval(abs(real(values,real64))),maxval(abs(aimag(values))))
  end function maximum_complex_component_4

  pure real(real64) function maximum_complex_component_5(values)result(maximum_value)
    complex(real64),intent(in)::values(:,:,:,:,:)
    maximum_value=0d0
    if(size(values)>0)maximum_value=max(maxval(abs(real(values,real64))),maxval(abs(aimag(values))))
  end function maximum_complex_component_5

  subroutine compute_symmetry_defects(representation,hamiltonian,target_defect,energy_defect)
    complex(real64),intent(in)::representation(:,:,:),hamiltonian(:,:)
    real(real64),intent(out)::target_defect,energy_defect
    complex(real64),allocatable::identity(:,:),work(:,:)
    integer::r,operation,i
    real(real64)::target_scale,energy_scale
    r=size(hamiltonian,1);allocate(identity(r,r),work(r,r));identity=(0d0,0d0)
    do i=1,r;identity(i,i)=1d0;enddo
    target_scale=sqrt(real(r,real64));energy_scale=max(1d0,frobenius(hamiltonian))
    target_defect=0d0;energy_defect=0d0
    do operation=1,size(representation,3)
      work=matmul(conjg(transpose(representation(:,:,operation))),representation(:,:,operation))-identity
      target_defect=max(target_defect,frobenius(work)/target_scale)
      work=matmul(conjg(transpose(representation(:,:,operation))),&
        matmul(hamiltonian,representation(:,:,operation)))-hamiltonian
      energy_defect=max(energy_defect,frobenius(work)/energy_scale)
    enddo
  end subroutine compute_symmetry_defects

  subroutine compute_covariance_defects(representation,rotations,scalar_operators,vector_operators,&
      tensor_operators,scalar_defect,vector_defect,tensor_defect)
    complex(real64),intent(in)::representation(:,:,:),scalar_operators(:,:,:)
    complex(real64),intent(in)::vector_operators(:,:,:,:),tensor_operators(:,:,:,:,:)
    real(real64),intent(in)::rotations(:,:,:)
    real(real64),intent(out)::scalar_defect,vector_defect,tensor_defect
    complex(real64),allocatable::transformed(:,:),expected(:,:)
    integer::r,operation,item,a,b,c,d
    real(real64)::scale
    r=size(representation,1);allocate(transformed(r,r),expected(r,r))
    scalar_defect=0d0;vector_defect=0d0;tensor_defect=0d0
    do operation=1,size(representation,3)
      do item=1,size(scalar_operators,3)
        transformed=matmul(conjg(transpose(representation(:,:,operation))),&
          matmul(scalar_operators(:,:,item),representation(:,:,operation)))
        scale=max(1d0,maxval(abs(scalar_operators(:,:,item))))
        scalar_defect=max(scalar_defect,maxval(abs(transformed-scalar_operators(:,:,item)))/scale)
      enddo
      do item=1,size(vector_operators,4)
        scale=max(1d0,maxval(abs(vector_operators(:,:,:,item))))
        do a=1,3
          transformed=matmul(conjg(transpose(representation(:,:,operation))),&
            matmul(vector_operators(:,:,a,item),representation(:,:,operation)))
          expected=(0d0,0d0)
          do b=1,3;expected=expected+rotations(a,b,operation)*vector_operators(:,:,b,item);enddo
          vector_defect=max(vector_defect,maxval(abs(transformed-expected))/scale)
        enddo
      enddo
      do item=1,size(tensor_operators,5)
        scale=max(1d0,maxval(abs(tensor_operators(:,:,:,:,item))))
        do b=1,3;do a=1,3
          transformed=matmul(conjg(transpose(representation(:,:,operation))),&
            matmul(tensor_operators(:,:,a,b,item),representation(:,:,operation)))
          expected=(0d0,0d0)
          do d=1,3;do c=1,3
            expected=expected+rotations(a,c,operation)*rotations(b,d,operation)*&
              tensor_operators(:,:,c,d,item)
          enddo;enddo
          tensor_defect=max(tensor_defect,maxval(abs(transformed-expected))/scale)
        enddo;enddo
      enddo
    enddo
  end subroutine compute_covariance_defects

#ifdef USE_MPI
  subroutine fingerprint_distributed_rows(comm,row_ids,rows,global_count,fingerprint,ok)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::rows(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr_hash,ierr_count,local_count,global_rows
    integer(int64)::entry,local_hash,global_hash
    local_hash=0_int64
    do i=1,size(row_ids)
      entry=mix_hash(int(z'243F6A8885A308D3',int64),row_ids(i))
      entry=mix_hash(entry,int(size(rows,2),int64))
      do j=1,size(rows,2);entry=mix_complex(entry,rows(i,j));enddo
      local_hash=ieor(local_hash,entry)
    enddo
    global_hash=0_int64;global_rows=0;local_count=size(row_ids)
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr_hash)
    call MPI_Allreduce(local_count,global_rows,1,MPI_INTEGER,MPI_SUM,comm,ierr_count)
    if(ierr_hash/=MPI_SUCCESS.or.ierr_count/=MPI_SUCCESS.or.global_rows/=global_count)then
      fingerprint=0_int64;ok=.false.;return
    endif
    fingerprint=mix_hash(int(z'13198A2E03707344',int64),int(global_count,int64))
    fingerprint=mix_hash(fingerprint,int(size(rows,2),int64));fingerprint=mix_hash(fingerprint,global_hash)
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
  end subroutine fingerprint_distributed_rows

  subroutine agree_integer(comm,value,minimum,maximum,ierr)
    integer,intent(in)::comm,value
    integer,intent(out)::minimum,maximum,ierr
    integer::ierr_minimum,ierr_maximum
    minimum=value;maximum=value
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr_minimum)
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr_maximum)
    if(ierr_minimum/=MPI_SUCCESS)then;ierr=ierr_minimum;else;ierr=ierr_maximum;endif
  end subroutine agree_integer

  subroutine agree_int64(comm,value,minimum,maximum,ierr)
    integer,intent(in)::comm
    integer(int64),intent(in)::value
    integer(int64),intent(out)::minimum,maximum
    integer,intent(out)::ierr
    integer::ierr_minimum,ierr_maximum
    minimum=value;maximum=value
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr_minimum)
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr_maximum)
    if(ierr_minimum/=MPI_SUCCESS)then;ierr=ierr_minimum;else;ierr=ierr_maximum;endif
  end subroutine agree_int64
#endif

  pure integer(int64) function mix_hash(hash,value)result(mixed)
    integer(int64),intent(in)::hash,value
    mixed=add_modulo_64(ishftc(hash,11),value)
    mixed=ieor(mixed,ishftc(mixed,25))
    mixed=add_modulo_64(mixed,int(z'9E3779B97F4A7C15',int64))
    mixed=ieor(mixed,ishft(mixed,-27))
    mixed=add_modulo_64(mixed,ishftc(value,17))
    mixed=ieor(mixed,ishftc(mixed,42))
  end function mix_hash

  pure integer(int64) function add_modulo_64(left,right)result(sum_value)
    integer(int64),intent(in)::left,right
    integer(int64),parameter::low_mask=int(z'00000000FFFFFFFF',int64)
    integer(int64)::low_value,high_value,carry
    low_value=iand(left,low_mask)+iand(right,low_mask)
    carry=ishft(low_value,-32)
    high_value=iand(ishft(left,-32),low_mask)+iand(ishft(right,-32),low_mask)+carry
    sum_value=ior(iand(low_value,low_mask),ishft(iand(high_value,low_mask),32))
  end function add_modulo_64

  pure integer(int64) function mix_real(hash,value)result(mixed)
    integer(int64),intent(in)::hash
    real(real64),intent(in)::value
    mixed=mix_hash(hash,transfer(value,0_int64))
  end function mix_real

  pure integer(int64) function mix_complex(hash,value)result(mixed)
    integer(int64),intent(in)::hash
    complex(real64),intent(in)::value
    mixed=mix_real(hash,real(value,real64));mixed=mix_real(mixed,aimag(value))
  end function mix_complex

  pure integer(int64) function hash_real_1(seed,values)result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:)
    integer::i
    hash=mix_hash(seed,int(size(values),int64))
    do i=1,size(values);hash=mix_real(hash,values(i));enddo
  end function hash_real_1

  pure integer(int64) function hash_real_2(seed,values)result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:,:)
    integer::i,j
    hash=mix_hash(mix_hash(seed,int(size(values,1),int64)),int(size(values,2),int64))
    do j=1,size(values,2);do i=1,size(values,1);hash=mix_real(hash,values(i,j));enddo;enddo
  end function hash_real_2

  pure integer(int64) function hash_real_3(seed,values)result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    hash=mix_hash(mix_hash(mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64))
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=mix_real(hash,values(i,j,k))
    enddo;enddo;enddo
  end function hash_real_3

  pure integer(int64) function hash_complex_2(seed,values)result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:)
    integer::i,j
    hash=mix_hash(mix_hash(seed,int(size(values,1),int64)),int(size(values,2),int64))
    do j=1,size(values,2);do i=1,size(values,1);hash=mix_complex(hash,values(i,j));enddo;enddo
  end function hash_complex_2

  pure integer(int64) function hash_complex_3(seed,values)result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    hash=mix_hash(mix_hash(mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64))
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=mix_complex(hash,values(i,j,k))
    enddo;enddo;enddo
  end function hash_complex_3

  pure integer(int64) function hash_complex_4(seed,values)result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:,:)
    integer::i,j,k,l
    hash=mix_hash(mix_hash(mix_hash(mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64)),int(size(values,4),int64))
    do l=1,size(values,4);do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=mix_complex(hash,values(i,j,k,l))
    enddo;enddo;enddo;enddo
  end function hash_complex_4

  pure integer(int64) function hash_complex_5(seed,values)result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:,:,:)
    integer::i,j,k,l,m
    hash=mix_hash(mix_hash(mix_hash(mix_hash(mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64)),int(size(values,4),int64)),&
      int(size(values,5),int64))
    do m=1,size(values,5);do l=1,size(values,4);do k=1,size(values,3)
      do j=1,size(values,2);do i=1,size(values,1)
        hash=mix_complex(hash,values(i,j,k,l,m))
      enddo;enddo
    enddo;enddo;enddo
  end function hash_complex_5

  pure integer(int64) function fingerprint_input_payload(eigenvalues,representation,rotations,&
      scalar_operators,vector_operators,tensor_operators)result(hash)
    real(real64),intent(in)::eigenvalues(:),rotations(:,:,:)
    complex(real64),intent(in)::representation(:,:,:),scalar_operators(:,:,:)
    complex(real64),intent(in)::vector_operators(:,:,:,:),tensor_operators(:,:,:,:,:)
    hash=hash_real_1(int(z'6A09E667F3BCC909',int64),eigenvalues)
    hash=hash_complex_3(hash,representation);hash=hash_real_3(hash,rotations)
    hash=hash_complex_3(hash,scalar_operators);hash=hash_complex_4(hash,vector_operators)
    hash=hash_complex_5(hash,tensor_operators)
    if(hash==0_int64)hash=1_int64
  end function fingerprint_input_payload

  pure integer(int64) function fingerprint_localization(transform,centers,spreads_before,&
      spreads_after,iterations,symmetry_constrained,converged)result(hash)
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::centers(:,:),spreads_before(:),spreads_after(:)
    integer,intent(in)::iterations
    logical,intent(in)::symmetry_constrained,converged
    hash=hash_complex_2(int(z'BB67AE8584CAA73B',int64),transform)
    hash=hash_real_2(hash,centers);hash=hash_real_1(hash,spreads_before)
    hash=hash_real_1(hash,spreads_after);hash=mix_hash(hash,int(iterations,int64))
    hash=mix_hash(hash,merge(1_int64,0_int64,symmetry_constrained))
    hash=mix_hash(hash,merge(1_int64,0_int64,converged))
    if(hash==0_int64)hash=1_int64
  end function fingerprint_localization

  pure integer(int64) function fingerprint_operator_payload(eigenvalues,occupied_amplitudes,&
      metric,hamiltonian,representation,rotations,scalar_operators,vector_operators,&
      tensor_operators)result(hash)
    real(real64),intent(in)::eigenvalues(:),rotations(:,:,:)
    complex(real64),intent(in)::occupied_amplitudes(:,:),metric(:,:),hamiltonian(:,:),&
      representation(:,:,:),scalar_operators(:,:,:),vector_operators(:,:,:,:),&
      tensor_operators(:,:,:,:,:)
    hash=hash_real_1(int(z'3C6EF372FE94F82B',int64),eigenvalues)
    hash=hash_complex_2(hash,occupied_amplitudes);hash=hash_complex_2(hash,metric)
    hash=hash_complex_2(hash,hamiltonian);hash=hash_complex_3(hash,representation)
    hash=hash_real_3(hash,rotations);hash=hash_complex_3(hash,scalar_operators)
    hash=hash_complex_4(hash,vector_operators);hash=hash_complex_5(hash,tensor_operators)
    if(hash==0_int64)hash=1_int64
  end function fingerprint_operator_payload

  pure integer(int64) function fingerprint_result_receipt(result,c_hash,u_hash,b_hash,&
      operator_hash)result(hash)
    type(s_dg_hybrid_certified_rt_basis),intent(in)::result
    integer(int64),intent(in),optional::c_hash,u_hash,b_hash,operator_hash
    integer(int64)::c_value,u_value,b_value,operator_value
    c_value=result%c_cert_fingerprint;u_value=result%localization_fingerprint
    b_value=result%b_rt_fingerprint;operator_value=result%operator_fingerprint
    if(present(c_hash))c_value=c_hash;if(present(u_hash))u_value=u_hash
    if(present(b_hash))b_value=b_hash;if(present(operator_hash))operator_value=operator_hash
    hash=mix_hash(int(z'A54FF53A5F1D36F1',int64),int(result%global_count,int64))
    hash=mix_hash(hash,int(result%certified_rank,int64));hash=mix_hash(hash,int(result%noccupied,int64))
    hash=mix_hash(hash,c_value);hash=mix_hash(hash,u_value);hash=mix_hash(hash,b_value)
    hash=mix_hash(hash,operator_value);hash=mix_hash(hash,int(result%localization_iterations,int64))
    hash=mix_real(hash,result%spread_before_total);hash=mix_real(hash,result%spread_after_total)
    hash=mix_real(hash,result%spread_improvement);hash=mix_real(hash,result%transform_unitarity_defect)
    hash=mix_real(hash,result%certified_metric_defect);hash=mix_real(hash,result%rt_metric_defect)
    hash=mix_real(hash,result%embedding_defect);hash=mix_real(hash,result%projector_invariance_defect)
    hash=mix_real(hash,result%target_symmetry_defect_before)
    hash=mix_real(hash,result%target_symmetry_defect_after)
    hash=mix_real(hash,result%energy_symmetry_defect_before)
    hash=mix_real(hash,result%energy_symmetry_defect_after)
    hash=mix_real(hash,result%symmetry_defect_invariance)
    hash=mix_real(hash,result%scalar_covariance_defect)
    hash=mix_real(hash,result%vector_covariance_defect)
    hash=mix_real(hash,result%tensor_covariance_defect)
    if(hash==0_int64)hash=1_int64
  end function fingerprint_result_receipt
end module dg_hybrid_certified_rt_basis
