#include "config.h"
module dg_overlapping_wannier_construction
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_overlapping_wannier_metric,only:assemble_dg_overlapping_wannier_metric
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_overlapping_wannier_construction
    integer::candidate_rank=0,target_rank=0,retained_rank=0,generation=0
    integer(int64)::transform_fingerprint=0_int64
    integer,allocatable::center_owner_rank(:),center_owner_fragment(:)
    integer(int64),allocatable::physical_grid_ids(:)
    complex(real64),allocatable::value(:,:),gradient(:,:,:),transform(:,:)
    real(real64)::occupied_inclusion_residual=huge(1d0)
    real(real64)::boundary_value_max=huge(1d0),boundary_gradient_max=huge(1d0)
    real(real64)::metric_minimum_eigenvalue=0d0,metric_condition_number=huge(1d0)
  end type
  public::construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction
contains
  subroutine construct_dg_overlapping_wannier_basis(comm,ncandidate,ntarget,noccupied,physical_ids,&
      core_fragment,weights,localization_coordinate,boundary_mask,candidate_value,candidate_gradient,&
      occupied_coefficients,expected_core_count,generation,boundary_value_tolerance,&
      boundary_gradient_tolerance,rank_tolerance,result,ok,message)
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
#ifdef USE_MPI
    complex(real64),allocatable::s_local(:,:),s(:,:),l_local(:,:),localizer(:,:),x(:,:),&
      occ_gram(:,:),occ_vectors(:,:),a_occ(:,:),overlap_occ_x(:,:),residual(:,:),&
      residual_gram(:,:),residual_vectors(:,:),q_comp(:,:),comp_localizer(:,:),&
      comp_vectors(:,:),transform(:,:),final_metric(:,:),metric_vectors(:,:)
    real(real64),allocatable::spectrum(:),occ_spectrum(:),residual_spectrum(:),comp_spectrum(:),&
      metric_spectrum(:),local_score(:),all_score(:,:)
    integer(int64),allocatable::local_center_id(:),all_center_id(:,:)
    integer,allocatable::local_center_fragment(:),all_center_fragment(:,:)
    logical,allocatable::complete_pairs(:,:)
    integer::nlocal,i,j,p,ierr,nproc,rank,local_bad,global_bad,metric_rejected,ownership_count,&
      ncomp,nneed,start_index,owner,r,matrix_count
    real(real64)::largest,local_boundary_value,local_boundary_gradient,occ_scale
    integer(int64)::local_fingerprint,global_fingerprint,quantized_density,matrix_count64

    ok=.false.;message='';nlocal=size(physical_ids)
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
    if(any(physical_ids<=0_int64).or.any(core_fragment<=0).or.any(weights<=0d0))local_bad=1
    if(ncandidate>0)then
      if(int(ncandidate,int64)>huge(1_int64)/int(ncandidate,int64))then
        local_bad=1;matrix_count64=0_int64
      else
        matrix_count64=int(ncandidate,int64)*int(ncandidate,int64)
        if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
      endif
    else
      matrix_count64=0_int64
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid overlapping-Wannier construction metadata';return;endif

    allocate(s_local(ncandidate,ncandidate),l_local(ncandidate,ncandidate))
    s_local=(0d0,0d0);l_local=(0d0,0d0)
    do p=1,nlocal
      do j=1,ncandidate;do i=1,ncandidate
        s_local(i,j)=s_local(i,j)+weights(p)*conjg(candidate_value(i,p))*candidate_value(j,p)
        l_local(i,j)=l_local(i,j)+weights(p)*localization_coordinate(p)*&
          conjg(candidate_value(i,p))*candidate_value(j,p)
      enddo;enddo
    enddo
    allocate(s(ncandidate,ncandidate),localizer(ncandidate,ncandidate))
    matrix_count=int(matrix_count64)
    call MPI_Allreduce(s_local,s,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(l_local,localizer,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    s=0.5d0*(s+conjg(transpose(s)));localizer=0.5d0*(localizer+conjg(transpose(localizer)))
    if(.not.positive_definite_above(s,rank_tolerance))then
      message='candidate rank loss in overlapping-Wannier construction';return
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

    occ_gram=matmul(conjg(transpose(occupied_coefficients)),matmul(s,occupied_coefficients))
    occ_gram=0.5d0*(occ_gram+conjg(transpose(occ_gram)))
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

    overlap_occ_x=matmul(conjg(transpose(a_occ)),matmul(s,x))
    residual=x-matmul(a_occ,overlap_occ_x)
    residual_gram=matmul(conjg(transpose(residual)),matmul(s,residual))
    residual_gram=0.5d0*(residual_gram+conjg(transpose(residual_gram)))
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
      call hermitian_eigensystem(comp_localizer,comp_spectrum,comp_vectors,ok,message)
      if(.not.ok)return
      allocate(transform(ncandidate,ntarget))
      transform(:,1:noccupied)=a_occ
      transform(:,noccupied+1:ntarget)=matmul(q_comp,comp_vectors(:,1:nneed))
    else
      allocate(transform,source=a_occ)
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

    allocate(complete_pairs(ntarget,nlocal));complete_pairs=.true.
    call assemble_dg_overlapping_wannier_metric(comm,ntarget,physical_ids,weights,result%value,&
      complete_pairs,expected_core_count,rank_tolerance,final_metric,metric_vectors,metric_spectrum,&
      result%metric_minimum_eigenvalue,result%metric_condition_number,metric_rejected,&
      ownership_count,ok,message)
    if(.not.ok)return
    if(metric_rejected/=0.or.ownership_count/=expected_core_count)then
      ok=.false.;message='retained overlapping-Wannier metric rank or ownership loss';return
    endif

    occ_gram=matmul(conjg(transpose(occupied_coefficients)),matmul(s,occupied_coefficients))
    overlap_occ_x=matmul(conjg(transpose(transform)),matmul(s,occupied_coefficients))
    result%occupied_inclusion_residual=maxval(abs(occ_gram-&
      matmul(conjg(transpose(overlap_occ_x)),overlap_occ_x)))/max(1d0,maxval(abs(occ_gram)))
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

    call MPI_Comm_size(comm,nproc,ierr);call MPI_Comm_rank(comm,rank,ierr)
    allocate(local_score(ntarget),local_center_id(ntarget),local_center_fragment(ntarget))
    local_score=-1d0;local_center_id=huge(1_int64);local_center_fragment=0
    do j=1,ntarget;do p=1,nlocal
      occ_scale=abs(result%value(j,p))**2
      if(occ_scale>local_score(j).or.(occ_scale==local_score(j).and.physical_ids(p)<local_center_id(j)))then
        local_score(j)=occ_scale;local_center_id(j)=physical_ids(p);local_center_fragment(j)=core_fragment(p)
      endif
    enddo;enddo
    allocate(all_score(ntarget,nproc),all_center_id(ntarget,nproc),all_center_fragment(ntarget,nproc))
    call MPI_Allgather(local_score,ntarget,MPI_DOUBLE_PRECISION,all_score,ntarget,MPI_DOUBLE_PRECISION,comm,ierr)
    call MPI_Allgather(local_center_id,ntarget,MPI_INTEGER8,all_center_id,ntarget,MPI_INTEGER8,comm,ierr)
    call MPI_Allgather(local_center_fragment,ntarget,MPI_INTEGER,all_center_fragment,ntarget,MPI_INTEGER,comm,ierr)
    allocate(result%center_owner_rank(ntarget),result%center_owner_fragment(ntarget))
    do j=1,ntarget
      owner=1
      do r=2,nproc
        if(all_score(j,r)>all_score(j,owner).or.(all_score(j,r)==all_score(j,owner).and.&
            all_center_id(j,r)<all_center_id(j,owner)))owner=r
      enddo
      result%center_owner_rank(j)=owner-1
      result%center_owner_fragment(j)=all_center_fragment(j,owner)
    enddo

    local_fingerprint=0_int64
    if(rank==0)local_fingerprint=ieor(int(generation,int64),ishftc(int(ntarget,int64),11))
    do p=1,nlocal
      quantized_density=nint(sum(abs(result%value(:,p))**2)*1d10,int64)
      local_fingerprint=ieor(local_fingerprint,ieor(ishftc(physical_ids(p),17),quantized_density))
    enddo
    call MPI_Allreduce(local_fingerprint,global_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    result%transform_fingerprint=global_fingerprint
    if(result%transform_fingerprint==0_int64)result%transform_fingerprint=1_int64
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier construction requires MPI'
#endif
  end subroutine

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
    if(allocated(result%value))deallocate(result%value)
    if(allocated(result%gradient))deallocate(result%gradient)
    if(allocated(result%transform))deallocate(result%transform)
  end subroutine
end module dg_overlapping_wannier_construction
