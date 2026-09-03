#include "config.h"
module dg_hybrid_wannier_complement
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_generalized_metric_factor
    logical::valid=.false.
    integer::global_row_count=0
    integer::wannier_rank=0
    integer::metric_rank=0
    real(real64)::metric_tolerance=0d0
    real(real64)::metric_cutoff=0d0
    real(real64)::metric_condition=0d0
    integer(int64)::wannier_fingerprint=0_int64
    integer(int64)::semantic_fingerprint=0_int64
    integer(int64)::local_binding_fingerprint=0_int64
    complex(real64),allocatable::gram(:,:)
    complex(real64),allocatable::inverse(:,:)
    complex(real64),allocatable::eigenvectors(:,:)
    complex(real64),allocatable::retained_projector(:,:)
    real(real64),allocatable::eigenvalues(:)
  end type s_dg_hybrid_generalized_metric_factor
  public::project_dg_hybrid_wannier_complement,compute_dg_hybrid_wannier_projection_tile,&
    materialize_dg_hybrid_projected_pw_tile,prepare_dg_hybrid_generalized_wannier_metric,&
    apply_dg_hybrid_generalized_wannier_projection_tile,&
    compute_dg_hybrid_generalized_wannier_projection_tile,build_dg_hybrid_complete_union_map
  interface
    subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
      import::real64
      character(1),intent(in)::jobz,uplo
      integer,intent(in)::n,lda,lwork
      complex(real64),intent(inout)::a(lda,*)
      real(real64),intent(out)::w(*)
      complex(real64),intent(inout)::work(*)
      real(real64),intent(inout)::rwork(*)
      integer,intent(out)::info
    end subroutine zheev
  end interface
contains
  subroutine materialize_dg_hybrid_projected_pw_tile(wannier_buffer,raw_pw_buffer,projection_coefficients,&
      projected_pw_buffer,ok,message)
    complex(real64),intent(in)::wannier_buffer(:,:),raw_pw_buffer(:,:),projection_coefficients(:,:)
    complex(real64),allocatable,intent(out)::projected_pw_buffer(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nw,width,npoint,allocation_status
    ok=.false.;message='';nw=size(wannier_buffer,1);npoint=size(wannier_buffer,2);width=size(raw_pw_buffer,1)
    if(nw<1.or.width<1.or.npoint<1.or.size(raw_pw_buffer,2)/=npoint.or.&
        any(shape(projection_coefficients)/=[nw,width]))then
      message='invalid projected PW buffer shape';return
    endif
    if(.not.finite_complex(wannier_buffer).or..not.finite_complex(raw_pw_buffer).or.&
        .not.finite_complex(projection_coefficients))then
      message='nonfinite projected PW buffer input';return
    endif
    allocate(projected_pw_buffer(width,npoint),stat=allocation_status)
    if(allocation_status/=0)then
      message='cannot allocate projected PW buffer';return
    endif
    projected_pw_buffer=raw_pw_buffer-matmul(transpose(projection_coefficients),wannier_buffer)
    if(.not.finite_complex(projected_pw_buffer))then
      deallocate(projected_pw_buffer);message='nonfinite projected PW buffer output';return
    endif
    ok=.true.;message=''
  end subroutine materialize_dg_hybrid_projected_pw_tile

  subroutine compute_dg_hybrid_wannier_projection_tile(comm,global_row_count,row_ids,weights,wannier_values,&
      pw_tile,wannier_fingerprint,packet_fingerprint,first_column,tolerance,coefficients,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count,first_column
    integer(int64),intent(in)::row_ids(:),wannier_fingerprint,packet_fingerprint
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(in)::wannier_values(:,:),pw_tile(:,:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nw,width,nlocal,i,j,k,ierr,local_bad,global_bad,minimum_integer,maximum_integer
    integer,allocatable::ownership(:)
    integer(int64)::minimum_bits,maximum_bits,bits
    complex(real64),allocatable::local_coefficients(:,:),gram_local(:,:),gram_global(:,:)
    real(real64)::gram_defect
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;local_bad=0
    nlocal=size(row_ids);nw=size(wannier_values,1);width=size(pw_tile,1)
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent projection row extent';return;endif
    call agree_integer(first_column,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent projection tile origin';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent projection tolerance';return;endif
    call agree_int64(wannier_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.wannier_fingerprint==0_int64)then
      message='invalid projection Wannier provenance';return
    endif
    call agree_int64(packet_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.packet_fingerprint==0_int64)then
      message='invalid projection packet provenance';return
    endif
    if(global_row_count<1.or.first_column<1.or.nw<1.or.width<1.or.&
        size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal.or.size(pw_tile,2)/=nlocal)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_row_count),int64)))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or..not.ieee_is_finite(tolerance).or.&
        tolerance<1d-15.or.tolerance>1d-2.or..not.finite_complex(wannier_values).or.&
        .not.finite_complex(pw_tile))local_bad=1
    allocate(ownership(max(1,global_row_count)));ownership=0
    do i=1,nlocal
      if(row_ids(i)>=1_int64.and.row_ids(i)<=int(max(0,global_row_count),int64))&
        ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,max(0,global_row_count),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
    elseif(global_row_count>0)then
      if(any(ownership(:global_row_count)/=1))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid projection tile contract';return;endif
    allocate(local_coefficients(nw,width),coefficients(nw,width),gram_local(nw,nw),gram_global(nw,nw))
    do j=1,width;do i=1,nw
      local_coefficients(i,j)=sum(weights*conjg(wannier_values(i,:))*pw_tile(j,:))
    enddo;enddo
    do j=1,nw;do i=1,nw
      gram_local(i,j)=sum(weights*conjg(wannier_values(i,:))*wannier_values(j,:))
    enddo;enddo
    call MPI_Allreduce(local_coefficients,coefficients,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projection coefficient reduction failed';return;endif
    call MPI_Allreduce(gram_local,gram_global,nw*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projection Gram reduction failed';return;endif
    gram_defect=0d0
    do j=1,nw;do i=1,nw
      if(i==j)then;gram_defect=max(gram_defect,abs(gram_global(i,j)-1d0))
      else;gram_defect=max(gram_defect,abs(gram_global(i,j)));endif
    enddo;enddo
    if(gram_defect>100d0*tolerance)then;message='projection Wannier frame is not orthonormal';return;endif
    workspace_peak_bytes=16_int64*int(2*nw*width+2*nw*nw,int64)+4_int64*int(global_row_count,int64)
    fingerprint=ieor(wannier_fingerprint,ishftc(packet_fingerprint,13))
    fingerprint=ieor(fingerprint,int(first_column,int64))
    do j=1,width;do i=1,nw
      bits=transfer(real(coefficients(i,j)),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=transfer(aimag(coefficients(i,j)),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
    enddo;enddo
    if(fingerprint==0_int64)fingerprint=1877_int64
    ok=.true.;message=''
#else
    ok=.false.;message='hybrid projection tile requires MPI';workspace_peak_bytes=0_int64;fingerprint=0_int64
    allocate(coefficients(0,0))
#endif
  end subroutine compute_dg_hybrid_wannier_projection_tile

  subroutine prepare_dg_hybrid_generalized_wannier_metric(comm,global_row_count,row_ids,weights,&
      wannier_values,wannier_fingerprint,metric_tolerance,factor,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:),wannier_fingerprint
    real(real64),intent(in)::weights(:),metric_tolerance
    complex(real64),intent(in)::wannier_values(:,:)
    type(s_dg_hybrid_generalized_metric_factor),intent(out)::factor
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nw,i,j,k,ierr,rank,minimum_integer,maximum_integer
    integer::local_code,global_code,allocation_status,info,lwork,lrwork
    integer,allocatable::ownership(:)
    integer(int64)::bits,minimum_bits,maximum_bits,local_row_xor,global_row_xor,n2,candidate_fingerprint
    integer(int64)::complex_elements,real_elements,integer_elements
    real(real64)::hermitian_defect,hermitian_scale,spectrum_scale,roundoff_floor,negative_limit
    real(real64)::smallest_retained,metric_cutoff,candidate_condition
    integer::candidate_rank
    integer(int64)::candidate_workspace
    complex(real64),allocatable::gram_local(:,:),gram_global(:,:),eigenvectors(:,:),inverse(:,:),projector(:,:),work(:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    logical::workspace_ok

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    nlocal=size(row_ids);nw=size(wannier_values,1)
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent generalized metric global row count';return
    endif
    call agree_integer(nw,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent generalized Wannier rank';return
    endif
    bits=transfer(metric_tolerance,bits)
    call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent generalized metric tolerance';return
    endif
    call agree_int64(wannier_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.wannier_fingerprint==0_int64)then
      message='invalid or inconsistent generalized Wannier provenance';return
    endif

    local_code=0
    if(global_row_count<1.or.nw<1.or..not.product_fits_default_integer(nw,nw))local_code=1
    if(size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal)local_code=max(local_code,1)
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_row_count),int64)))local_code=max(local_code,2)
    if(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0))local_code=max(local_code,3)
    if(.not.finite_complex(wannier_values))local_code=max(local_code,4)
    if(.not.ieee_is_finite(metric_tolerance).or.metric_tolerance<1d-15.or.metric_tolerance>1d-2)&
      local_code=max(local_code,5)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      select case(global_code)
      case(2);message='invalid generalized metric spatial row ID'
      case(3);message='invalid generalized metric row weight'
      case(4);message='nonfinite generalized Wannier values'
      case(5);message='invalid generalized metric tolerance'
      case default;message='invalid generalized Wannier rank or shape'
      end select
      return
    endif

    allocate(ownership(global_row_count),stat=allocation_status)
    local_code=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      if(allocated(ownership))deallocate(ownership)
      message='cannot allocate generalized metric row ownership';return
    endif
    ownership=0
    do i=1,nlocal
      ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    local_code=0
    if(ierr/=MPI_SUCCESS)then
      local_code=1
    elseif(any(ownership/=1))then
      local_code=1
    endif
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      deallocate(ownership);message='duplicate or missing generalized metric spatial row';return
    endif

    local_row_xor=local_basis_row_xor(row_ids,weights,wannier_values)
    call MPI_Allreduce(local_row_xor,global_row_xor,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      deallocate(ownership);message='generalized metric row fingerprint reduction failed';return
    endif
    candidate_fingerprint=int(z'9E3779B97F4A7C15',int64)
    call mix_hash_word(candidate_fingerprint,int(global_row_count,int64))
    call mix_hash_word(candidate_fingerprint,int(nw,int64))
    call mix_hash_word(candidate_fingerprint,wannier_fingerprint)
    call mix_hash_word(candidate_fingerprint,bits)
    call mix_hash_word(candidate_fingerprint,global_row_xor)
    if(candidate_fingerprint==0_int64)candidate_fingerprint=1_int64

    n2=int(nw,int64)*int(nw,int64);lwork=max(1,2*nw-1);lrwork=max(1,3*nw-2)
    allocate(gram_local(nw,nw),gram_global(nw,nw),eigenvectors(nw,nw),inverse(nw,nw),projector(nw,nw),&
      eigenvalues(nw),work(lwork),rwork(lrwork),stat=allocation_status)
    local_code=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_prepare_workspace();message='cannot allocate generalized metric factor workspace';return
    endif
    gram_local=(0d0,0d0)
    do j=1,nw;do i=1,nw
      gram_local(i,j)=sum(weights*conjg(wannier_values(i,:))*wannier_values(j,:))
    enddo;enddo
    call MPI_Allreduce(gram_local,gram_global,nw*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_prepare_workspace();message='generalized Wannier Gram reduction failed';return
    endif
    hermitian_scale=max(1d0,maxval(abs(gram_global)))
    hermitian_defect=maxval(abs(gram_global-conjg(transpose(gram_global))))
    if(.not.ieee_is_finite(hermitian_defect).or.hermitian_defect>10d0*metric_tolerance*hermitian_scale)then
      call cleanup_prepare_workspace();message='generalized Wannier Gram is not Hermitian';return
    endif
    gram_global=0.5d0*(gram_global+conjg(transpose(gram_global)))
    do i=1,nw;gram_global(i,i)=cmplx(real(gram_global(i,i),real64),0d0,real64);enddo

    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_prepare_workspace();message='generalized metric communicator query failed';return
    endif
    eigenvectors=gram_global;info=0
    if(rank==0)call zheev('V','U',nw,eigenvectors,nw,eigenvalues,work,lwork,rwork,info)
    call MPI_Bcast(info,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.info/=0)then
      call cleanup_prepare_workspace();message='generalized Wannier Gram diagonalization failed';return
    endif
    call MPI_Bcast(eigenvalues,nw,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Bcast(eigenvectors,nw*nw,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_prepare_workspace();message='generalized metric eigensystem broadcast failed';return
    endif
    if(.not.all(ieee_is_finite(eigenvalues)).or..not.finite_complex(eigenvectors))then
      call cleanup_prepare_workspace();message='nonfinite generalized metric eigensystem';return
    endif

    spectrum_scale=max(1d0,maxval(abs(eigenvalues)))
    roundoff_floor=64d0*epsilon(1d0)*spectrum_scale*real(max(1,nw),real64)
    negative_limit=roundoff_floor
    metric_cutoff=max(metric_tolerance*spectrum_scale,roundoff_floor)
    if(any(eigenvalues< -negative_limit))then
      call cleanup_prepare_workspace();message='indefinite generalized Wannier Gram matrix';return
    endif
    if(any(abs(eigenvalues-metric_cutoff)<=16d0*roundoff_floor))then
      call cleanup_prepare_workspace();message='ambiguous generalized Wannier metric rank at cutoff';return
    endif
    candidate_rank=count(eigenvalues>metric_cutoff)
    if(candidate_rank<1)then
      call cleanup_prepare_workspace();message='generalized Wannier metric has zero retained rank';return
    endif
    smallest_retained=huge(1d0)
    do k=1,nw
      if(eigenvalues(k)>metric_cutoff)smallest_retained=min(smallest_retained,eigenvalues(k))
    enddo
    candidate_condition=maxval(eigenvalues)/smallest_retained
    if(.not.ieee_is_finite(candidate_condition).or.candidate_condition<1d0)then
      call cleanup_prepare_workspace();message='invalid generalized Wannier metric condition';return
    endif

    inverse=(0d0,0d0);projector=(0d0,0d0)
    do k=1,nw
      if(eigenvalues(k)<=metric_cutoff)cycle
      do j=1,nw;do i=1,nw
        projector(i,j)=projector(i,j)+eigenvectors(i,k)*conjg(eigenvectors(j,k))
        inverse(i,j)=inverse(i,j)+eigenvectors(i,k)*conjg(eigenvectors(j,k))/eigenvalues(k)
      enddo;enddo
    enddo
    if(.not.finite_complex(inverse).or..not.finite_complex(projector))then
      call cleanup_prepare_workspace();message='nonfinite generalized metric inverse';return
    endif

    complex_elements=5_int64*n2+int(lwork,int64)
    real_elements=int(nw+lrwork,int64);integer_elements=int(global_row_count,int64)
    call calculate_workspace_bytes(complex_elements,real_elements,integer_elements,candidate_workspace,workspace_ok)
    local_code=merge(0,1,workspace_ok)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_prepare_workspace();message='generalized metric workspace receipt overflow';return
    endif

    factor%global_row_count=global_row_count;factor%wannier_rank=nw;factor%metric_rank=candidate_rank
    factor%metric_cutoff=metric_cutoff;factor%metric_condition=candidate_condition
    factor%metric_tolerance=metric_tolerance;factor%wannier_fingerprint=wannier_fingerprint
    factor%semantic_fingerprint=candidate_fingerprint
    factor%local_binding_fingerprint=local_binding_fingerprint(row_ids,weights,wannier_values)
    call move_alloc(gram_global,factor%gram);call move_alloc(inverse,factor%inverse)
    call move_alloc(eigenvectors,factor%eigenvectors);call move_alloc(projector,factor%retained_projector)
    call move_alloc(eigenvalues,factor%eigenvalues)
    factor%valid=.true.;workspace_peak_bytes=candidate_workspace;fingerprint=candidate_fingerprint
    ok=.true.;message=''
    call cleanup_prepare_workspace()
#else
    ok=.false.;message='generalized Wannier metric requires MPI'
    workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup_prepare_workspace()
      if(allocated(ownership))deallocate(ownership)
      if(allocated(gram_local))deallocate(gram_local)
      if(allocated(gram_global))deallocate(gram_global)
      if(allocated(eigenvectors))deallocate(eigenvectors)
      if(allocated(inverse))deallocate(inverse)
      if(allocated(projector))deallocate(projector)
      if(allocated(eigenvalues))deallocate(eigenvalues)
      if(allocated(work))deallocate(work)
      if(allocated(rwork))deallocate(rwork)
    end subroutine cleanup_prepare_workspace
#endif
  end subroutine prepare_dg_hybrid_generalized_wannier_metric

  subroutine apply_dg_hybrid_generalized_wannier_projection_tile(comm,row_ids,weights,wannier_values,&
      pw_tile,packet_fingerprint,first_column,factor,coefficients,projected_pw,orthogonality_defect,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,first_column
    integer(int64),intent(in)::row_ids(:),packet_fingerprint
    real(real64),intent(in)::weights(:)
    complex(real64),intent(in)::wannier_values(:,:),pw_tile(:,:)
    type(s_dg_hybrid_generalized_metric_factor),intent(in)::factor
    complex(real64),allocatable,intent(out)::coefficients(:,:),projected_pw(:,:)
    real(real64),intent(out)::orthogonality_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nw,width,i,j,ierr,minimum_integer,maximum_integer
    integer::local_code,global_code,allocation_status
    integer(int64)::minimum_bits,maximum_bits,local_tile_hash,global_tile_hash,n2,candidate_fingerprint
    integer(int64)::complex_elements,real_elements
    integer(int64)::candidate_workspace
    real(real64)::cross_scale,acceptance_limit,candidate_defect
    complex(real64),allocatable::cross_local(:,:),cross_global(:,:),candidate_coefficients(:,:),candidate_projected(:,:)
    logical::workspace_ok

    ok=.false.;message='';orthogonality_defect=0d0;workspace_peak_bytes=0_int64;fingerprint=0_int64
    nlocal=size(row_ids);nw=size(wannier_values,1);width=size(pw_tile,1)
    call agree_integer(nw,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent generalized Wannier rank in projection tile';return
    endif
    call agree_integer(width,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent generalized projection tile width';return
    endif
    call agree_integer(first_column,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent generalized projection tile origin';return
    endif
    call agree_int64(packet_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.packet_fingerprint==0_int64)then
      message='invalid or inconsistent generalized projection packet provenance';return
    endif

    local_code=0
    if(.not.factor%valid.or.first_column<1.or.nw<1.or.width<1)local_code=1
    if(size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal.or.size(pw_tile,2)/=nlocal)local_code=max(local_code,1)
    if(.not.finite_complex(wannier_values).or..not.finite_complex(pw_tile))local_code=max(local_code,2)
    if(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0))local_code=max(local_code,3)
    if(.not.allocated(factor%gram).or..not.allocated(factor%inverse).or..not.allocated(factor%eigenvectors).or.&
        .not.allocated(factor%retained_projector).or..not.allocated(factor%eigenvalues))local_code=max(local_code,4)
    if(local_code<4)then
      if(factor%global_row_count<1.or.factor%wannier_rank/=nw.or.factor%metric_rank<1.or.&
          factor%metric_rank>nw.or.factor%wannier_fingerprint==0_int64.or.factor%semantic_fingerprint==0_int64)&
        local_code=max(local_code,4)
      if(any(row_ids<1_int64).or.any(row_ids>int(factor%global_row_count,int64)))local_code=max(local_code,5)
      if(any(shape(factor%gram)/=[nw,nw]).or.any(shape(factor%inverse)/=[nw,nw]).or.&
          any(shape(factor%eigenvectors)/=[nw,nw]).or.any(shape(factor%retained_projector)/=[nw,nw]).or.&
          size(factor%eigenvalues)/=nw)local_code=max(local_code,4)
    endif
    if(local_code==0)then
      if(local_binding_fingerprint(row_ids,weights,wannier_values)/=factor%local_binding_fingerprint)local_code=6
      if(.not.finite_complex(factor%gram).or..not.finite_complex(factor%inverse).or.&
          .not.finite_complex(factor%eigenvectors).or..not.finite_complex(factor%retained_projector).or.&
          .not.all(ieee_is_finite(factor%eigenvalues)))local_code=4
    endif
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      select case(global_code)
      case(2);message='nonfinite generalized projection tile input'
      case(3);message='invalid generalized projection row weight'
      case(4);message='invalid generalized metric factor'
      case(5);message='invalid generalized projection spatial row'
      case(6);message='generalized metric factor does not match Wannier rows or weights'
      case default;message='invalid generalized projection tile shape'
      end select
      return
    endif
    call agree_int64(factor%semantic_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent generalized metric factor provenance';return
    endif
    if(.not.product_fits_default_integer(nw,width))then
      message='generalized projection tile MPI count overflow';return
    endif

    allocate(cross_local(nw,width),cross_global(nw,width),candidate_coefficients(nw,width),&
      candidate_projected(width,nlocal),stat=allocation_status)
    local_code=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_apply_workspace();message='cannot allocate generalized projection tile workspace';return
    endif
    do j=1,width;do i=1,nw
      cross_local(i,j)=sum(weights*conjg(wannier_values(i,:))*pw_tile(j,:))
    enddo;enddo
    call MPI_Allreduce(cross_local,cross_global,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_apply_workspace();message='generalized Wannier-PW overlap reduction failed';return
    endif
    cross_scale=max(1d0,maxval(abs(cross_global)))
    candidate_coefficients=matmul(factor%inverse,cross_global)
    candidate_projected=pw_tile-matmul(transpose(candidate_coefficients),wannier_values)
    local_code=merge(0,1,finite_complex(candidate_coefficients).and.finite_complex(candidate_projected))
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_apply_workspace();message='nonfinite generalized projected PW values';return
    endif
    do j=1,width;do i=1,nw
      cross_local(i,j)=sum(weights*conjg(wannier_values(i,:))*candidate_projected(j,:))
    enddo;enddo
    call MPI_Allreduce(cross_local,cross_global,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_apply_workspace();message='generalized projection residual reduction failed';return
    endif
    cross_local=matmul(factor%retained_projector,cross_global)
    candidate_defect=maxval(abs(cross_local))
    acceptance_limit=max(100d0*factor%metric_tolerance*cross_scale,&
      512d0*epsilon(1d0)*real(max(1,nw),real64)*cross_scale)
    if(.not.ieee_is_finite(candidate_defect).or.candidate_defect>acceptance_limit)then
      call cleanup_apply_workspace();message='generalized projected PW is not metric orthogonal';return
    endif

    local_tile_hash=local_basis_row_xor(row_ids,weights,pw_tile)
    call MPI_Allreduce(local_tile_hash,global_tile_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup_apply_workspace();message='generalized projection fingerprint reduction failed';return
    endif
    candidate_fingerprint=int(z'D1B54A32D192ED03',int64)
    call mix_hash_word(candidate_fingerprint,factor%semantic_fingerprint)
    call mix_hash_word(candidate_fingerprint,packet_fingerprint)
    call mix_hash_word(candidate_fingerprint,int(first_column,int64))
    call mix_hash_word(candidate_fingerprint,int(width,int64))
    call mix_hash_word(candidate_fingerprint,int(factor%metric_rank,int64))
    call mix_hash_word(candidate_fingerprint,global_tile_hash)
    if(candidate_fingerprint==0_int64)candidate_fingerprint=2_int64

    n2=int(nw,int64)*int(nw,int64)
    complex_elements=4_int64*n2+3_int64*int(nw,int64)*int(width,int64)+&
      int(width,int64)*int(nlocal,int64)
    real_elements=int(nw,int64)
    call calculate_workspace_bytes(complex_elements,real_elements,0_int64,candidate_workspace,workspace_ok)
    local_code=merge(0,1,workspace_ok)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_apply_workspace();message='generalized projection workspace receipt overflow';return
    endif
    call move_alloc(candidate_coefficients,coefficients);call move_alloc(candidate_projected,projected_pw)
    orthogonality_defect=candidate_defect;workspace_peak_bytes=candidate_workspace
    fingerprint=candidate_fingerprint
    call cleanup_apply_workspace();ok=.true.;message=''
#else
    ok=.false.;message='generalized Wannier projection requires MPI'
    orthogonality_defect=0d0;workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup_apply_workspace()
      if(allocated(cross_local))deallocate(cross_local)
      if(allocated(cross_global))deallocate(cross_global)
      if(allocated(candidate_coefficients))deallocate(candidate_coefficients)
      if(allocated(candidate_projected))deallocate(candidate_projected)
    end subroutine cleanup_apply_workspace
#endif
  end subroutine apply_dg_hybrid_generalized_wannier_projection_tile

  subroutine compute_dg_hybrid_generalized_wannier_projection_tile(comm,global_row_count,row_ids,weights,&
      wannier_values,pw_tile,wannier_fingerprint,packet_fingerprint,first_column,metric_tolerance,&
      coefficients,projected_pw,metric_rank,metric_condition,orthogonality_defect,workspace_peak_bytes,&
      fingerprint,ok,message)
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
    type(s_dg_hybrid_generalized_metric_factor)::factor
    integer::candidate_rank
    integer(int64)::prepare_workspace,prepare_fingerprint,candidate_workspace,candidate_fingerprint
    real(real64)::candidate_condition,candidate_defect
    logical::stage_ok
    character(256)::stage_message

    ok=.false.;message='';metric_rank=0;metric_condition=0d0
    orthogonality_defect=0d0;workspace_peak_bytes=0_int64;fingerprint=0_int64
    call prepare_dg_hybrid_generalized_wannier_metric(comm,global_row_count,row_ids,weights,wannier_values,&
      wannier_fingerprint,metric_tolerance,factor,prepare_workspace,prepare_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    candidate_rank=factor%metric_rank;candidate_condition=factor%metric_condition
    call apply_dg_hybrid_generalized_wannier_projection_tile(comm,row_ids,weights,wannier_values,pw_tile,&
      packet_fingerprint,first_column,factor,coefficients,projected_pw,candidate_defect,candidate_workspace,&
      candidate_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    metric_rank=candidate_rank;metric_condition=candidate_condition;orthogonality_defect=candidate_defect
    workspace_peak_bytes=max(prepare_workspace,candidate_workspace);fingerprint=candidate_fingerprint
    ok=.true.;message=''
  end subroutine compute_dg_hybrid_generalized_wannier_projection_tile

  subroutine build_dg_hybrid_complete_union_map(comm,global_row_count,row_ids,weights,&
      uncompressed_basis_values,metric_tolerance,complete_basis_transform,complete_basis_values,&
      metric_rank,metric_condition,projector_fingerprint,ok,message)
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
#ifdef USE_MPI
    type(s_dg_hybrid_generalized_metric_factor)::factor
    complex(real64),allocatable::candidate_transform(:,:),candidate_values(:,:),discarded_metric(:,:)
    integer::n,nlocal,i,j,k,retained_column,ierr,local_code,global_code,allocation_status,candidate_metric_rank
    integer(int64)::workspace,source_fingerprint,quantized,bits,prepare_fingerprint,candidate_fingerprint
    real(real64)::orthogonality_defect,retained_defect,discarded_defect,scale,quantum,candidate_condition
    logical::stage_ok
    character(256)::stage_message

    ok=.false.;message='';metric_rank=0;metric_condition=0d0;projector_fingerprint=0_int64
    n=size(uncompressed_basis_values,1);nlocal=size(row_ids)
    source_fingerprint=int(z'94D049BB133111EB',int64)
    call mix_hash_word(source_fingerprint,int(n,int64));call mix_hash_word(source_fingerprint,int(global_row_count,int64))
    bits=transfer(metric_tolerance,bits);call mix_hash_word(source_fingerprint,bits)
    if(source_fingerprint==0_int64)source_fingerprint=3_int64
    call prepare_dg_hybrid_generalized_wannier_metric(comm,global_row_count,row_ids,weights,&
      uncompressed_basis_values,source_fingerprint,metric_tolerance,factor,workspace,prepare_fingerprint,&
      stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    candidate_metric_rank=factor%metric_rank;candidate_condition=factor%metric_condition
    allocate(candidate_transform(n,candidate_metric_rank),candidate_values(candidate_metric_rank,nlocal),&
      discarded_metric(n,n),&
      stat=allocation_status)
    local_code=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_union_workspace();message='cannot allocate complete Hybrid union map';return
    endif
    candidate_transform=(0d0,0d0)
    if(candidate_metric_rank==n)then
      do i=1,n;candidate_transform(i,i)=(1d0,0d0);enddo
      candidate_values=uncompressed_basis_values
    else
      retained_column=0
      do k=1,n
        if(factor%eigenvalues(k)<=factor%metric_cutoff)cycle
        retained_column=retained_column+1
        candidate_transform(:,retained_column)=factor%eigenvectors(:,k)
      enddo
      if(retained_column/=candidate_metric_rank)then
        call cleanup_union_workspace();message='complete Hybrid union retained-rank construction failed';return
      endif
      candidate_values=matmul(transpose(candidate_transform),uncompressed_basis_values)
    endif
    local_code=merge(0,1,finite_complex(candidate_transform).and.finite_complex(candidate_values))
    call MPI_Allreduce(local_code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      call cleanup_union_workspace();message='nonfinite complete Hybrid union map';return
    endif

    orthogonality_defect=0d0;retained_defect=0d0
    do j=1,candidate_metric_rank;do i=1,candidate_metric_rank
      if(i==j)then
        orthogonality_defect=max(orthogonality_defect,&
          abs(sum(conjg(candidate_transform(:,i))*candidate_transform(:,j))-1d0))
      else
        orthogonality_defect=max(orthogonality_defect,&
          abs(sum(conjg(candidate_transform(:,i))*candidate_transform(:,j))))
      endif
      retained_defect=max(retained_defect,abs(sum(conjg(candidate_transform(:,i))*&
        matmul(factor%gram,candidate_transform(:,j)))))
    enddo;enddo
    if(orthogonality_defect>100d0*metric_tolerance)then
      call cleanup_union_workspace();message='complete Hybrid union transform is not orthonormal';return
    endif
    discarded_metric=factor%gram-matmul(factor%retained_projector,factor%gram)
    discarded_defect=maxval(abs(discarded_metric));scale=max(1d0,maxval(abs(factor%gram)))
    if(discarded_defect>max(10d0*factor%metric_cutoff*real(max(1,n),real64),&
        512d0*epsilon(1d0)*scale*real(max(1,n),real64)))then
      call cleanup_union_workspace();message='complete Hybrid union discarded metric is not null';return
    endif
    if(.not.ieee_is_finite(retained_defect).or.retained_defect<=factor%metric_cutoff)then
      call cleanup_union_workspace();message='complete Hybrid union retained metric is singular';return
    endif

    candidate_fingerprint=int(z'BF58476D1CE4E5B9',int64)
    call mix_hash_word(candidate_fingerprint,factor%semantic_fingerprint)
    call mix_hash_word(candidate_fingerprint,int(candidate_metric_rank,int64))
    quantum=max(100d0*metric_tolerance,1024d0*epsilon(1d0))
    do j=1,n;do i=1,n
      if(candidate_metric_rank==n)then
        quantized=merge(nint(1d0/quantum,int64),0_int64,i==j)
      else
        quantized=nint(real(factor%retained_projector(i,j),real64)/quantum,int64)
      endif
      call mix_hash_word(candidate_fingerprint,quantized)
      if(candidate_metric_rank==n)then
        quantized=0_int64
      else
        quantized=nint(aimag(factor%retained_projector(i,j))/quantum,int64)
      endif
      call mix_hash_word(candidate_fingerprint,quantized)
    enddo;enddo
    if(candidate_fingerprint==0_int64)candidate_fingerprint=4_int64
    call move_alloc(candidate_transform,complete_basis_transform)
    call move_alloc(candidate_values,complete_basis_values)
    metric_rank=candidate_metric_rank;metric_condition=candidate_condition
    projector_fingerprint=candidate_fingerprint
    call cleanup_union_workspace();ok=.true.;message=''
#else
    ok=.false.;message='complete Hybrid union map requires MPI'
    metric_rank=0;metric_condition=0d0;projector_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup_union_workspace()
      if(allocated(candidate_transform))deallocate(candidate_transform)
      if(allocated(candidate_values))deallocate(candidate_values)
      if(allocated(discarded_metric))deallocate(discarded_metric)
    end subroutine cleanup_union_workspace
#endif
  end subroutine build_dg_hybrid_complete_union_map

  subroutine project_dg_hybrid_wannier_complement(comm,global_row_count,row_ids,weights,wannier_values,&
      pw_values,wannier_fingerprint,packet_fingerprint,packet_ids,near_offsets,near_wannier_ids,&
      diagnose_full_tail,tolerance,projected_values,&
      omitted_tail,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(in)::wannier_values(:,:),pw_values(:,:)
    integer(int64),intent(in)::wannier_fingerprint,packet_fingerprint
    integer,intent(in)::packet_ids(:),near_offsets(:),near_wannier_ids(:)
    logical,intent(in)::diagnose_full_tail
    complex(real64),allocatable,intent(out)::projected_values(:,:)
    real(real64),intent(out)::omitted_tail
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::block_width=8
    integer::rank,ierr,nlocal,nw,np,i,j,j0,j1,width,p,k,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer,root,position
    integer,allocatable::ownership_count(:),owner(:),owner_position(:)
    integer(int64)::bits,minimum_bits,maximum_bits,complex_elements,integer_elements,quantized
    complex(real64),allocatable::local_block(:,:),global_block(:,:),remote_row(:)
    complex(real64)::local_scalar,global_scalar
    real(real64)::local_max,global_value_scale,global_weight_scale,safe_weight,gram_defect,tail_square,&
      packet_tail_square,cross_defect,quantization_limit,value_bound_scale,global_quantization_scale
    logical::diagnose_min,diagnose_max
    ok=.false.;message='';omitted_tail=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64
    local_bad=0;nlocal=size(row_ids);nw=size(wannier_values,1);np=size(pw_values,1)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement spatial extent';return
    endif
    call agree_integer(nw,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement Wannier extent';return
    endif
    call agree_integer(np,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement PW extent';return
    endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent hybrid complement tolerance';return
    endif
    call agree_int64(wannier_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.wannier_fingerprint==0_int64)then
      message='invalid or inconsistent retained-Wannier provenance';return
    endif
    call agree_int64(packet_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.packet_fingerprint==0_int64)then
      message='invalid or inconsistent windowed-PW packet provenance';return
    endif
    call agree_logical(diagnose_full_tail,diagnose_min,diagnose_max,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.diagnose_min.neqv.diagnose_max)then
      message='inconsistent hybrid omitted-tail diagnostic mode';return
    endif
    if(global_row_count<1.or.nw<1.or.np<1)local_bad=1
    if(size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal.or.size(pw_values,2)/=nlocal)local_bad=1
    if(size(packet_ids)/=np.or.size(near_offsets)/=np+1)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.finite_complex(wannier_values).or.&
      .not.finite_complex(pw_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid complement shape or finite contract';return
    endif
    if(tolerance<1d-15.or.tolerance>1d-2.or.any(weights<=0d0).or.any(packet_ids<1))local_bad=1
    do i=1,np
      call agree_integer(packet_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid PW packet IDs';return
      endif
    enddo
    do i=1,np+1
      call agree_integer(near_offsets(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid complement neighbor offsets';return
      endif
    enddo
    do i=1,size(near_wannier_ids)
      call agree_integer(near_wannier_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid complement neighbor IDs';return
      endif
    enddo
    if(near_offsets(1)/=1.or.near_offsets(np+1)/=size(near_wannier_ids)+1)local_bad=1
    if(any(near_offsets(2:np+1)<near_offsets(1:np)))local_bad=1
    if(any(near_wannier_ids<1).or.any(near_wannier_ids>nw))local_bad=1
    do p=1,np
      do k=near_offsets(p)+1,near_offsets(p+1)-1
        if(near_wannier_ids(k)<=near_wannier_ids(k-1))local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid complement sparse-neighbor catalog';return
    endif
    allocate(ownership_count(global_row_count),owner(global_row_count),owner_position(global_row_count),&
      local_block(nw,min(block_width,nw)),global_block(nw,min(block_width,nw)),&
      remote_row(np),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid complement workspace';return
    endif
    ownership_count=0;owner=-1;owner_position=0
    do i=1,nlocal
      ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1
      owner(int(row_ids(i)))=rank;owner_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement ownership count failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement position reduction failed';return;endif
    if(any(ownership_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='duplicate or missing hybrid complement spatial row';return
    endif
    local_max=0d0
    if(nlocal>0)local_max=max(maxval(abs(wannier_values)),maxval(abs(pw_values)))
    call MPI_Allreduce(local_max,global_value_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement value scale reduction failed';return;endif
    local_max=0d0;if(nlocal>0)local_max=maxval(weights)
    call MPI_Allreduce(local_max,global_weight_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement weight scale reduction failed';return;endif
    if(global_value_scale>0d0)then
      value_bound_scale=max(1d0,global_value_scale)
      safe_weight=huge(1d0)/64d0
      safe_weight=safe_weight/real(global_row_count,real64)/real(max(1,nw),real64)
      safe_weight=safe_weight/value_bound_scale/value_bound_scale/value_bound_scale
      if(global_weight_scale>safe_weight)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid complement multiplication magnitude is unsafe';return
    endif
    gram_defect=0d0
    do j0=1,nw,block_width
      j1=min(nw,j0+block_width-1);width=j1-j0+1;local_block(:,1:width)=(0d0,0d0)
      do j=1,width;do i=1,nw
        local_block(i,j)=sum(weights*conjg(wannier_values(i,:))*wannier_values(j0+j-1,:))
      enddo;enddo
      call MPI_Allreduce(local_block,global_block,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid Wannier Gram reduction failed';return;endif
      do j=1,width;do i=1,nw
        if(i==j0+j-1)then
          gram_defect=max(gram_defect,abs(global_block(i,j)-1d0))
        else
          gram_defect=max(gram_defect,abs(global_block(i,j)))
        endif
      enddo;enddo
    enddo
    if(gram_defect>100d0*tolerance)then
      call cleanup();message='retained Wannier frame is not orthonormal';return
    endif
    complex_elements=int(np,int64)*int(nlocal,int64)+&
      2_int64*int(nw,int64)*int(min(block_width,nw),int64)+&
      int(np,int64)
    integer_elements=3_int64*int(global_row_count,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.&
      integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>&
      huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid complement workspace receipt overflow';return
    endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(projected_values(np,nlocal),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate projected hybrid PW values';return
    endif
    projected_values=pw_values;tail_square=0d0
    do p=1,np
      packet_tail_square=0d0
      if(diagnose_full_tail)then
        do j0=1,nw,block_width
          j1=min(nw,j0+block_width-1);width=j1-j0+1
          do j=1,width
            local_block(j,1)=sum(weights*conjg(wannier_values(j0+j-1,:))*pw_values(p,:))
          enddo
          call MPI_Allreduce(local_block(1:width,1),global_block(1:width,1),width,&
            MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement overlap reduction failed';return;endif
          do j=1,width
            i=j0+j-1
            if(is_near(i,p,near_offsets,near_wannier_ids))then
              projected_values(p,:)=projected_values(p,:)-wannier_values(i,:)*global_block(j,1)
            else
              packet_tail_square=packet_tail_square+abs(global_block(j,1))**2
            endif
          enddo
        enddo
        tail_square=max(tail_square,packet_tail_square)
      else
        width=near_offsets(p+1)-near_offsets(p)
        do j=1,width
          i=near_wannier_ids(near_offsets(p)+j-1)
          local_scalar=sum(weights*conjg(wannier_values(i,:))*pw_values(p,:))
          call MPI_Allreduce(local_scalar,global_scalar,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid local overlap reduction failed';return;endif
          projected_values(p,:)=projected_values(p,:)-wannier_values(i,:)*global_scalar
        enddo
      endif
    enddo
    omitted_tail=sqrt(tail_square)
    if(diagnose_full_tail.and.omitted_tail>tolerance)then
      call cleanup();message='omitted Wannier projection tail exceeds tolerance';return
    endif
    cross_defect=0d0
    do p=1,np;do i=1,nw
      local_scalar=sum(weights*conjg(wannier_values(i,:))*projected_values(p,:))
      call MPI_Allreduce(local_scalar,global_scalar,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement residual reduction failed';return;endif
      cross_defect=max(cross_defect,abs(global_scalar))
    enddo;enddo
    if(cross_defect>2d0*tolerance)then
      call cleanup();message='projected PW is not orthogonal to retained Wanniers';return
    endif
    quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*tolerance
    local_max=0d0
    if(nlocal>0)local_max=max(maxval(abs(real(projected_values))),maxval(abs(aimag(projected_values))))
    call MPI_Allreduce(local_max,global_quantization_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup();message='hybrid complement quantization scale reduction failed';return
    endif
    if(global_quantization_scale>quantization_limit)then
      call cleanup();message='hybrid complement fingerprint range is unsafe';return
    endif
    fingerprint=ieor(int(z'A54FF53A5F1D36F1',int64),wannier_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),packet_fingerprint)
    do i=1,np
      fingerprint=ieor(ishftc(fingerprint,9),int(packet_ids(i),int64))
      fingerprint=ieor(ishftc(fingerprint,9),int(near_offsets(i+1)-near_offsets(i),int64))
      do k=near_offsets(i),near_offsets(i+1)-1
        fingerprint=ieor(ishftc(fingerprint,9),int(near_wannier_ids(k),int64))
      enddo
    enddo
    do i=1,global_row_count
      root=owner(i);remote_row=(0d0,0d0)
      if(rank==root)remote_row=projected_values(:,owner_position(i))
      call MPI_Bcast(remote_row,np,MPI_DOUBLE_COMPLEX,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement fingerprint broadcast failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),int(i,int64))
      do p=1,np
        quantized=nint(real(remote_row(p))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
        quantized=nint(aimag(remote_row(p))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='hybrid Wannier complement requires MPI';omitted_tail=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;allocate(projected_values(0,0))
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup()
      if(allocated(projected_values))deallocate(projected_values)
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position)
      if(allocated(local_block))deallocate(local_block)
      if(allocated(global_block))deallocate(global_block)
      if(allocated(remote_row))deallocate(remote_row)
    end subroutine cleanup
#endif
  end subroutine project_dg_hybrid_wannier_complement

  logical function is_near(wannier_id,pw_id,offsets,ids)
    integer,intent(in)::wannier_id,pw_id,offsets(:),ids(:)
    integer::k
    is_near=.false.
    do k=offsets(pw_id),offsets(pw_id+1)-1
      if(ids(k)==wannier_id)then;is_near=.true.;return;endif
    enddo
  end function is_near

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex

  logical function product_fits_default_integer(left,right)
    integer,intent(in)::left,right
    product_fits_default_integer=left>=0.and.right>=0
    if(.not.product_fits_default_integer)return
    if(right==0)return
    product_fits_default_integer=left<=huge(0)/right
  end function product_fits_default_integer

  integer(int64) function local_basis_row_xor(row_ids,weights,values)result(fingerprint)
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::weights(:)
    complex(real64),intent(in)::values(:,:)
    integer::row
    fingerprint=0_int64
    do row=1,size(row_ids)
      fingerprint=ieor(fingerprint,basis_row_fingerprint(row_ids(row),weights(row),values(:,row)))
    enddo
  end function local_basis_row_xor

  integer(int64) function local_binding_fingerprint(row_ids,weights,values)result(fingerprint)
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::weights(:)
    complex(real64),intent(in)::values(:,:)
    fingerprint=int(z'243F6A8885A308D3',int64)
    call mix_hash_word(fingerprint,int(size(row_ids),int64))
    call mix_hash_word(fingerprint,local_basis_row_xor(row_ids,weights,values))
    if(fingerprint==0_int64)fingerprint=5_int64
  end function local_binding_fingerprint

  integer(int64) function basis_row_fingerprint(row_id,weight,values)result(fingerprint)
    integer(int64),intent(in)::row_id
    real(real64),intent(in)::weight
    complex(real64),intent(in)::values(:)
    integer::i
    integer(int64)::bits
    fingerprint=int(z'13198A2E03707344',int64)
    call mix_hash_word(fingerprint,row_id)
    bits=transfer(weight,bits);call mix_hash_word(fingerprint,bits)
    call mix_hash_word(fingerprint,int(size(values),int64))
    do i=1,size(values)
      bits=transfer(real(values(i),real64),bits);call mix_hash_word(fingerprint,bits)
      bits=transfer(aimag(values(i)),bits);call mix_hash_word(fingerprint,bits)
    enddo
    if(fingerprint==0_int64)fingerprint=ieor(row_id,6_int64)
  end function basis_row_fingerprint

  subroutine mix_hash_word(fingerprint,value)
    integer(int64),intent(inout)::fingerprint
    integer(int64),intent(in)::value
    integer(int64)::mixed,cross_bits
    mixed=ieor(value,ishftc(value,23))
    cross_bits=iand(ishftc(fingerprint,17),ishftc(mixed,41))
    fingerprint=ieor(ieor(ishftc(fingerprint,29),mixed),cross_bits)
    cross_bits=iand(not(fingerprint),ishftc(mixed,7))
    fingerprint=ieor(ieor(fingerprint,cross_bits),shiftr(mixed,19))
  end subroutine mix_hash_word

  subroutine calculate_workspace_bytes(complex_elements,real_elements,integer_elements,bytes,ok)
    integer(int64),intent(in)::complex_elements,real_elements,integer_elements
    integer(int64),intent(out)::bytes
    logical,intent(out)::ok
    bytes=0_int64;ok=.true.
    call add_workspace_bytes(bytes,complex_elements,16_int64,ok)
    call add_workspace_bytes(bytes,real_elements,8_int64,ok)
    call add_workspace_bytes(bytes,integer_elements,4_int64,ok)
  end subroutine calculate_workspace_bytes

  subroutine add_workspace_bytes(total,elements,element_bytes,ok)
    integer(int64),intent(inout)::total
    integer(int64),intent(in)::elements,element_bytes
    logical,intent(inout)::ok
    if(.not.ok)return
    if(elements<0_int64.or.element_bytes<0_int64)then;ok=.false.;return;endif
    if(element_bytes>0_int64)then
      if(elements>(huge(total)-total)/element_bytes)then;ok=.false.;return;endif
    endif
    total=total+elements*element_bytes
  end subroutine add_workspace_bytes

#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value
    integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
  subroutine agree_logical(value,minimum_value,maximum_value,comm,ierr)
    logical,intent(in)::value
    logical,intent(out)::minimum_value,maximum_value
    integer,intent(in)::comm
    integer,intent(out)::ierr
    integer::input,minimum_integer,maximum_integer
    input=merge(1,0,value)
    call agree_integer(input,minimum_integer,maximum_integer,comm,ierr)
    minimum_value=minimum_integer==1;maximum_value=maximum_integer==1
  end subroutine agree_logical
#endif
end module dg_hybrid_wannier_complement
