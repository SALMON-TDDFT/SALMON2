#include "config.h"
module dg_overlapping_wannier_w90
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::estimate_dg_w90_coordinator_bytes,validate_dg_w90_result
  public::setup_dg_w90_gamma_library,run_dg_w90_gamma_library
  public::assemble_dg_w90_gamma_matrices
  public::apply_dg_w90_gamma_transform
  public::inherit_dg_w90_affine_receipts
  public::validate_dg_w90_convergence_log
  public::align_dg_w90_character_sector_gauge
  public::validate_dg_w90_localization_cluster
contains
  subroutine validate_dg_w90_localization_cluster(eigenvalues,selected_count,tolerance,ok,message)
    ! The localization eigensolver must supply this spectrum in ascending order.
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer,intent(in)::selected_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i
    real(real64)::scale
    ok=.false.;message=''
    if(size(eigenvalues)<1.or.selected_count<1.or.selected_count>size(eigenvalues).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or..not.all(ieee_is_finite(eigenvalues)))then
      message='invalid Wannier90 localization-cluster contract';return
    endif
    do i=2,size(eigenvalues)
      if(eigenvalues(i)<eigenvalues(i-1))then;message='Wannier90 localization spectrum is not ordered';return;endif
    enddo
    if(selected_count<size(eigenvalues))then
      scale=max(1d0,maxval(abs(eigenvalues)))
      if(abs(eigenvalues(selected_count+1)-eigenvalues(selected_count))<=tolerance*scale)then
        message='Wannier90 selection splits a degenerate localization cluster';return
      endif
    endif
    ok=.true.
  end subroutine validate_dg_w90_localization_cluster

#ifdef USE_MPI
  subroutine align_dg_w90_character_sector_gauge(comm,row_ids,sector_rows,reference_rows,&
      gamma_rows,conjugate_rows,gamma_sewing_defect,tolerance,&
      aligned_rows,aligned_conjugate_rows,singular_values,&
      canonical_channel_keys,polar_defect,gamma_defect,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::sector_rows(:,:),reference_rows(:,:),gamma_rows(:,:),conjugate_rows(:,:)
    real(real64),intent(in)::gamma_sewing_defect,tolerance
    complex(real64),allocatable,intent(out)::aligned_rows(:,:),aligned_conjugate_rows(:,:)
    real(real64),allocatable,intent(out)::singular_values(:)
    integer(int64),allocatable,intent(out)::canonical_channel_keys(:)
    real(real64),intent(out)::polar_defect,gamma_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::link(:,:),polar(:,:),svd_left(:,:),svd_right(:,:),gram(:,:),&
      remote_row(:),generated(:,:),svd_work(:)
    complex(real64),allocatable::projector_row(:)
    complex(real64),allocatable::ordered_reference(:,:)
    real(real64),allocatable::svd_rwork(:)
    integer,allocatable::owner(:),position(:),ownership_count(:)
    integer,allocatable::channel_order(:)
    integer::nlocal,n,m,i,j,k,rank,ierr,svd_info,svd_lwork,local_bad,global_bad,allocation_status
    integer::minimum_n,maximum_n,minimum_m,maximum_m,retained_singular_count
    integer(int64)::complex_elements,real_elements,integer_elements,byte_term
    real(real64)::minimum_tolerance,maximum_tolerance,minimum_gamma_receipt,maximum_gamma_receipt,&
      singular_scale,pivot_magnitude,candidate_magnitude
    complex(real64)::pivot_value,stream_value,phase_factor
    logical::receipt_valid
    complex(real64)::projector_value
    interface
      subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
        character,intent(in)::jobu,jobvt
        integer,intent(in)::m,n,lda,ldu,ldvt,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::s(*),rwork(*)
        complex(8),intent(out)::u(ldu,*),vt(ldvt,*)
        integer,intent(out)::info
      end subroutine
    end interface
    ok=.false.;message='';polar_defect=huge(1d0);gamma_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64
    nlocal=size(row_ids);m=size(sector_rows,2);n=size(gamma_rows,2)
    local_bad=merge(0,1,n>=1.and.m>=1.and.nlocal>=0.and.size(sector_rows,1)==nlocal.and.&
        all(shape(reference_rows)==[nlocal,m]).and.all(shape(conjugate_rows)==[nlocal,m]).and.&
        gamma_sewing_defect>=0d0.and.gamma_sewing_defect<=tolerance.and.&
        ieee_is_finite(gamma_sewing_defect).and.&
        size(gamma_rows,1)==nlocal.and.tolerance>=1d-15.and.tolerance<=1d-2.and.&
        ieee_is_finite(tolerance).and.all(row_ids>=1_int64).and.all(row_ids<=int(n,int64)).and.&
        all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
        all(ieee_is_finite(real(reference_rows))).and.all(ieee_is_finite(aimag(reference_rows))).and.&
        all(ieee_is_finite(real(gamma_rows))).and.all(ieee_is_finite(aimag(gamma_rows))).and.&
        all(ieee_is_finite(real(conjugate_rows))).and.all(ieee_is_finite(aimag(conjugate_rows))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid Wannier90 character-sector alignment contract';return
    endif
    call MPI_Allreduce(n,minimum_n,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(n,maximum_n,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(m,minimum_m,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(m,maximum_m,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_n/=maximum_n.or.minimum_m/=maximum_m)then
      message='Wannier90 alignment metadata disagree across ranks';return
    endif
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(gamma_sewing_defect,minimum_gamma_receipt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(gamma_sewing_defect,maximum_gamma_receipt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_tolerance/=maximum_tolerance.or.&
        minimum_gamma_receipt/=maximum_gamma_receipt)then
      message='Wannier90 alignment metadata disagree across ranks';return
    endif
    if(m>0.and.m>huge(0)/m)then;message='Wannier90 sector-link MPI count overflows';return;endif
    if(m>huge(0)/5)then;message='Wannier90 SVD workspace extent overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 communicator rank query failed';return;endif
    allocate(owner(n),position(n),ownership_count(n),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 ownership allocation failed';return;endif
    owner=0;position=0;ownership_count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;ownership_count(int(row_ids(i)))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 ownership rank reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 ownership position reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      message='Wannier90 sector rows are not uniquely owned';return
    endif
    allocate(link(m,m),polar(m,m),svd_left(m,m),svd_right(m,m),gram(m,m),singular_values(m),&
      ordered_reference(nlocal,m),channel_order(m),canonical_channel_keys(m),&
      aligned_rows(nlocal,m),aligned_conjugate_rows(nlocal,m),generated(nlocal,m),remote_row(m),&
      projector_row(n),&
      stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 sector workspace allocation failed';return;endif
    ordered_reference=reference_rows
    do j=1,m
      candidate_magnitude=maxval(abs(ordered_reference(:,j)))
      call MPI_Allreduce(candidate_magnitude,pivot_magnitude,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference pivot reduction failed';return;endif
      pivot_value=(0d0,0d0)
      do k=1,n
        stream_value=(0d0,0d0)
        if(rank==owner(k)-1)stream_value=ordered_reference(position(k),j)
        call MPI_Bcast(stream_value,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference phase stream failed';return;endif
        if(abs(stream_value)>=pivot_magnitude-10d0*tolerance*max(1d0,pivot_magnitude))then
          pivot_value=stream_value;exit
        endif
      enddo
      pivot_magnitude=abs(pivot_value)
      if(pivot_magnitude<=tolerance)then;message='singular Wannier90 reference channel';return;endif
      phase_factor=conjg(pivot_value)/pivot_magnitude
      ordered_reference(:,j)=ordered_reference(:,j)*phase_factor
      canonical_channel_keys(j)=int(z'BB67AE8584CAA73B',int64)
      do k=1,n
        stream_value=(0d0,0d0)
        if(rank==owner(k)-1)stream_value=ordered_reference(position(k),j)
        call MPI_Bcast(stream_value,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference key stream failed';return;endif
        canonical_channel_keys(j)=ieor(ishftc(canonical_channel_keys(j),11),&
          nint(real(stream_value,real64)/(1000d0*tolerance),int64))
        canonical_channel_keys(j)=ieor(ishftc(canonical_channel_keys(j),11),&
          nint(aimag(stream_value)/(1000d0*tolerance),int64))
      enddo
    enddo
    channel_order=[(i,i=1,m)]
    do i=2,m
      k=channel_order(i);j=i-1
      do while(j>=1)
        if(canonical_channel_keys(channel_order(j))<=canonical_channel_keys(k))exit
        channel_order(j+1)=channel_order(j);j=j-1
      enddo
      channel_order(j+1)=k
    enddo
    do i=2,m
      if(canonical_channel_keys(channel_order(i))==canonical_channel_keys(channel_order(i-1)))then
        message='Wannier90 canonical reference-channel fingerprints collide';return
      endif
    enddo
    ordered_reference=ordered_reference(:,channel_order)
    canonical_channel_keys=canonical_channel_keys(channel_order)
    gram=matmul(conjg(transpose(sector_rows)),sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 input sector frame is not orthonormal';return
    endif
    gram=matmul(conjg(transpose(ordered_reference)),ordered_reference)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 reference frame is not orthonormal';return
    endif
    gram=matmul(conjg(transpose(conjugate_rows)),conjugate_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 conjugate sector frame is not orthonormal';return
    endif
    link=matmul(conjg(transpose(sector_rows)),ordered_reference)
    call MPI_Allreduce(MPI_IN_PLACE,link,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 localization-link reduction failed';return;endif
    allocate(svd_rwork(max(1,5*m)),svd_work(1),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 SVD query allocation failed';return;endif
    svd_left=link
    svd_lwork=-1
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,&
      svd_lwork,svd_rwork,svd_info)
    if(svd_info/=0.or..not.ieee_is_finite(real(svd_work(1))))then
      message='Wannier90 localization-link SVD workspace query failed';return
    endif
    if(real(svd_work(1),real64)>real(huge(0),real64))then
      message='Wannier90 SVD workspace extent overflows';return
    endif
    svd_lwork=max(1,ceiling(real(svd_work(1),real64)));deallocate(svd_work)
    allocate(svd_work(svd_lwork),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 SVD workspace allocation failed';return;endif
    svd_left=link
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,&
      svd_lwork,svd_rwork,svd_info)
    if(svd_info/=0.or..not.all(ieee_is_finite(singular_values)))then
      message='Wannier90 localization-link SVD failed';return
    endif
    singular_scale=max(1d0,maxval(singular_values))
    retained_singular_count=count(singular_values>tolerance*singular_scale)
    if(retained_singular_count>0.and.retained_singular_count<m)then
      if(abs(singular_values(retained_singular_count)-singular_values(retained_singular_count+1))<=&
          10d0*tolerance*singular_scale)then
        message='Wannier90 rank threshold splits a degenerate singular-value cluster';return
      endif
    endif
    if(minval(singular_values)<=tolerance*singular_scale)then
      message='singular Wannier90 localization link';return
    endif
    polar=matmul(polar,svd_right)
    gram=matmul(conjg(transpose(polar)),polar)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    polar_defect=maxval(abs(gram));if(polar_defect>10d0*tolerance)then
      message='Wannier90 localization-link polar factor is not unitary';return
    endif
    aligned_rows=matmul(sector_rows,polar)
    generated=(0d0,0d0)
    do k=1,n
      remote_row=(0d0,0d0)
      if(rank==owner(k)-1)remote_row=conjg(aligned_rows(position(k),:))
      call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma row stream failed';return;endif
      do i=1,nlocal;generated(i,:)=generated(i,:)+gamma_rows(i,k)*remote_row;enddo
    enddo
    aligned_conjugate_rows=generated
    gram=matmul(conjg(transpose(generated)),generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 Gamma image leaks outside a unitary conjugate-sector frame';return
    endif
    gram=matmul(conjg(transpose(conjugate_rows)),generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma overlap reduction failed';return;endif
    gram=matmul(conjg(transpose(gram)),gram)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    gamma_defect=maxval(abs(gram));if(gamma_defect>10d0*tolerance)then
      message='Wannier90 aligned sectors violate Gamma conjugate pairing';return
    endif
    fingerprint=int(z'6A09E667F3BCC909',int64)
    do k=1,n
      remote_row=(0d0,0d0)
      if(rank==owner(k)-1)remote_row=aligned_rows(position(k),:)
      call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 aligned-projector row stream failed';return;endif
      projector_row=(0d0,0d0)
      do i=1,nlocal
        projector_row(int(row_ids(i)))=sum(remote_row*conjg(aligned_rows(i,:)))
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,projector_row,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 aligned-projector fingerprint reduction failed';return;endif
      if(rank==0)then
        do j=1,n
          projector_value=projector_row(j)
          fingerprint=ieor(ishftc(fingerprint,13),&
            nint(real(projector_value,real64)/(100d0*tolerance),int64))
          fingerprint=ieor(ishftc(fingerprint,13),nint(aimag(projector_value)/(100d0*tolerance),int64))
        enddo
      endif
    enddo
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 alignment fingerprint broadcast failed';return;endif
    complex_elements=0_int64;real_elements=0_int64;integer_elements=0_int64;receipt_valid=.true.
    call checked_add(complex_elements,size(sector_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(reference_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(gamma_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(conjugate_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(link,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(polar,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_left,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_right,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_work,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(gram,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(ordered_reference,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(aligned_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(aligned_conjugate_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(generated,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(remote_row,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(projector_row,kind=int64),receipt_valid)
    ! Conservative allowance for the largest MATMUL/LAPACK temporary owned by this routine.
    call checked_add(complex_elements,max(int(m,int64)*int(m,int64),int(nlocal,int64)*int(m,int64)),receipt_valid)
    call checked_add(real_elements,size(singular_values,kind=int64),receipt_valid)
    call checked_add(real_elements,size(svd_rwork,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(row_ids,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(canonical_channel_keys,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(owner,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(position,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(ownership_count,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(channel_order,kind=int64),receipt_valid)
    if(receipt_valid)call checked_product([complex_elements,16_int64],workspace_peak_bytes,receipt_valid)
    if(receipt_valid)call checked_product([real_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(receipt_valid)call checked_product([integer_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(.not.receipt_valid.or.workspace_peak_bytes<=0_int64)then
      workspace_peak_bytes=0_int64;message='Wannier90 alignment workspace receipt overflows';return
    endif
    ok=.true.;message=''
  end subroutine align_dg_w90_character_sector_gauge
#endif

  subroutine validate_dg_w90_convergence_log(path,maximum_iterations,iterations,ok,message)
    character(*),intent(in)::path
    integer,intent(in)::maximum_iterations
    integer,intent(out)::iterations
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::unit,io,parsed_iteration
    character(1024)::line
    logical::exists,have_final,have_convergence
    iterations=-1;ok=.false.;message='';have_final=.false.;have_convergence=.false.
    if(len_trim(path)==0.or.maximum_iterations<1)then
      message='invalid Wannier90 convergence-log contract';return
    endif
    inquire(file=trim(path),exist=exists)
    if(.not.exists)then;message='Wannier90 convergence log is missing';return;endif
    open(newunit=unit,file=trim(path),status='old',action='read',iostat=io)
    if(io/=0)then;message='Wannier90 convergence log cannot be opened';return;endif
    do
      read(unit,'(a)',iostat=io)line
      if(io/=0)exit
      if(index(line,'<-- CONV')>0)then
        read(line,*,iostat=io)parsed_iteration
        if(io==0)iterations=max(iterations,parsed_iteration)
        io=0
      endif
      if(index(line,'Wannierisation convergence criteria satisfied')>0)have_convergence=.true.
      if(index(adjustl(line),'Final State')==1)have_final=.true.
    enddo
    close(unit)
    if(.not.have_final)then;message='Wannier90 convergence log has no final state';return;endif
    if(.not.have_convergence)then;message='Wannier90 did not report convergence';return;endif
    if(iterations<0.or.iterations>=maximum_iterations)then
      message='Wannier90 exhausted its iteration limit';return
    endif
    ok=.true.;message=''
  end subroutine validate_dg_w90_convergence_log

  subroutine inherit_dg_w90_affine_receipts(transform,affine_subspace_defect,tolerance,&
      identity_defect,unitarity_defect,closure_defect,workspace_peak_bytes,ok,message)
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::affine_subspace_defect,tolerance
    real(real64),intent(out)::identity_defect,unitarity_defect,closure_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:)
    integer::nstate,i
    ok=.false.;message='';identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);workspace_peak_bytes=0_int64;nstate=size(transform,1)
    if(nstate<1.or.size(transform,2)/=nstate.or.tolerance<=0d0.or.&
       affine_subspace_defect<0d0.or..not.ieee_is_finite(tolerance).or.&
       .not.ieee_is_finite(affine_subspace_defect).or.&
       .not.all(ieee_is_finite(real(transform))).or.&
       .not.all(ieee_is_finite(aimag(transform))))then
      message='invalid MLWF affine-receipt inheritance contract';return
    endif
    allocate(gram(nstate,nstate));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nstate;gram(i,i)=gram(i,i)-1d0;enddo
    unitarity_defect=maxval(abs(gram))
    identity_defect=affine_subspace_defect;closure_defect=affine_subspace_defect
    workspace_peak_bytes=int(storage_size((0d0,0d0))/8,int64)*int(size(gram),int64)
    ok=max(identity_defect,max(unitarity_defect,closure_defect))<=tolerance
    if(.not.ok)message='MLWF gauge cannot inherit the accepted affine proof'
  end subroutine inherit_dg_w90_affine_receipts

  subroutine apply_dg_w90_gamma_transform(comm,physical_ids,values,gradients,transform,centers,&
      tolerance,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::physical_ids(:)
    complex(real64),intent(inout)::values(:,:),gradients(:,:,:),transform(:,:)
    real(real64),intent(inout)::centers(:,:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,npoint,i,j,k,axis,ierr,status
    integer,allocatable::order(:)
    logical,allocatable::used(:)
    complex(real64),allocatable::ordered_transform(:,:),new_values(:,:),new_gradients(:,:,:),gram(:,:)
    real(real64),allocatable::ordered_centers(:,:),local_maximum(:),global_maximum(:)
    integer(int64),allocatable::local_id(:),global_id(:)
    complex(real64),allocatable::local_pivot(:),global_pivot(:)
    real(real64)::scale
    logical::precedes
    ok=.false.;message='';status=0;nstate=size(values,1);npoint=size(values,2)
    if(nstate<=0.or.size(values,2)/=size(physical_ids).or.&
        any(shape(gradients)/=[3,nstate,npoint]).or.any(shape(transform)/=[nstate,nstate]).or.&
        any(shape(centers)/=[3,nstate]).or.tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(real(gradients))).or..not.all(ieee_is_finite(aimag(gradients))).or.&
        .not.all(ieee_is_finite(real(transform))).or..not.all(ieee_is_finite(aimag(transform))).or.&
        .not.all(ieee_is_finite(centers)).or.any(physical_ids<=0_int64))status=1
    if(maxval(abs(aimag(transform)))>tolerance*max(1d0,maxval(abs(transform))))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma MLWF transform contract';return;endif
    allocate(gram(nstate,nstate));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nstate;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>tolerance*max(1d0,real(nstate,real64)))status=2
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='Gamma MLWF transform is not unitary';return;endif
    allocate(order(nstate),used(nstate));used=.false.
    do i=1,nstate
      order(i)=0
      do j=1,nstate
        if(used(j))cycle
        if(order(i)==0)then
          order(i)=j
        else
          precedes=.false.
          do axis=1,3
            if(modulo(centers(axis,j),1d0)<modulo(centers(axis,order(i)),1d0)-tolerance)then
              precedes=.true.;exit
            else if(modulo(centers(axis,j),1d0)>modulo(centers(axis,order(i)),1d0)+tolerance)then
              exit
            endif
          enddo
          if(.not.precedes.and.all(abs(modulo(centers(:,j),1d0)-&
              modulo(centers(:,order(i)),1d0))<=tolerance))then
            do k=1,nstate
              if(abs(transform(k,j))>abs(transform(k,order(i)))+tolerance)then
                precedes=.true.;exit
              else if(abs(transform(k,j))<abs(transform(k,order(i)))-tolerance)then
                exit
              endif
            enddo
          endif
          if(precedes)order(i)=j
        endif
      enddo
      used(order(i))=.true.
    enddo
    allocate(ordered_transform(nstate,nstate),ordered_centers(3,nstate))
    ordered_transform=transform(:,order);ordered_centers=centers(:,order)
    allocate(new_values(nstate,npoint),new_gradients(3,nstate,npoint))
    new_values=matmul(transpose(ordered_transform),values)
    do axis=1,3;new_gradients(axis,:,:)=matmul(transpose(ordered_transform),gradients(axis,:,:));enddo
    allocate(local_maximum(nstate),global_maximum(nstate),local_id(nstate),global_id(nstate),&
      local_pivot(nstate),global_pivot(nstate))
    do i=1,nstate
      if(npoint>0)then
        j=maxloc(abs(new_values(i,:)),dim=1);local_maximum(i)=abs(new_values(i,j))
      else
        local_maximum(i)=-1d0
      endif
    enddo
    call MPI_Allreduce(local_maximum,global_maximum,nstate,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    local_id=huge(0_int64)
    do i=1,nstate;do j=1,npoint
      scale=max(1d0,global_maximum(i))
      if(abs(abs(new_values(i,j))-global_maximum(i))<=tolerance*scale)&
        local_id(i)=min(local_id(i),physical_ids(j))
    enddo;enddo
    call MPI_Allreduce(local_id,global_id,nstate,MPI_INTEGER8,MPI_MIN,comm,ierr)
    local_pivot=(0d0,0d0)
    do i=1,nstate;do j=1,npoint
      if(physical_ids(j)==global_id(i))local_pivot(i)=new_values(i,j)
    enddo;enddo
    call MPI_Allreduce(local_pivot,global_pivot,nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_maximum<=tolerance).or.&
        maxval(abs(aimag(global_pivot)))>tolerance*max(1d0,maxval(abs(global_pivot))))then
      message='cannot determine canonical Gamma MLWF signs';return
    endif
    do i=1,nstate
      if(real(global_pivot(i),real64)<0d0)then
        ordered_transform(:,i)=-ordered_transform(:,i);new_values(i,:)=-new_values(i,:)
        new_gradients(:,i,:)=-new_gradients(:,i,:)
      endif
    enddo
    transform=ordered_transform;centers=ordered_centers;values=new_values;gradients=new_gradients
    ok=.true.
#else
    ok=.false.;message='Gamma MLWF transform application requires MPI'
#endif
  end subroutine apply_dg_w90_gamma_transform

  subroutine assemble_dg_w90_gamma_matrices(comm,values,anchors,weights,fractional,nncell,&
      coordinator_byte_limit,m_matrix,a_matrix,coordinator_bytes,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,nncell(:,:)
    complex(real64),intent(in)::values(:,:),anchors(:,:)
    real(real64),intent(in)::weights(:),fractional(:,:)
    integer(int64),intent(in)::coordinator_byte_limit
    complex(real64),allocatable,intent(out)::m_matrix(:,:,:),a_matrix(:,:)
    integer(int64),intent(out)::coordinator_bytes,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::tile_size=32
    integer::rank,ierr,status,nband,nwann,npoint,nntot,m0,m1,n0,n1,b,p,m,n,count,allocation_status
    integer::local_dimensions(3),minimum_dimensions(3),maximum_dimensions(3)
    integer(int64)::output_elements,output_bytes,tile_bytes,complex_bytes,peak
    integer(int64)::minimum_limit,maximum_limit
    real(real64)::angle
    complex(real64)::phase
    complex(real64),allocatable::local_tile(:,:),reduced_tile(:,:)
    logical::arithmetic_ok
    ok=.false.;message='';coordinator_bytes=0_int64;workspace_peak_bytes=0_int64;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    nband=size(values,1);npoint=size(values,2);nwann=size(anchors,1);nntot=size(nncell,2)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.npoint<0.or.&
        size(anchors,2)/=npoint.or.size(weights)/=npoint.or.&
        any(shape(fractional)/=[3,npoint]).or.size(nncell,1)/=3.or.nntot<=0.or.&
        coordinator_byte_limit<0_int64.or..not.all(ieee_is_finite(real(values))).or.&
        .not.all(ieee_is_finite(aimag(values))).or..not.all(ieee_is_finite(real(anchors))).or.&
        .not.all(ieee_is_finite(aimag(anchors))).or..not.all(ieee_is_finite(weights)).or.&
        .not.all(ieee_is_finite(fractional)).or.any(weights<0d0))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));message='invalid distributed Wannier90 matrix contract';return
    endif
    local_dimensions=[nband,nwann,nntot]
    call MPI_Allreduce(local_dimensions,minimum_dimensions,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_dimensions,maximum_dimensions,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(coordinator_byte_limit,minimum_limit,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(coordinator_byte_limit,maximum_limit,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions).or.minimum_limit/=maximum_limit)status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));message='rank-inconsistent Wannier90 matrix contract';return
    endif
    call estimate_dg_w90_coordinator_bytes(nband,nwann,nntot,1,coordinator_bytes,&
      arithmetic_ok,message)
    if(.not.arithmetic_ok)status=2
    if(arithmetic_ok.and.coordinator_bytes>coordinator_byte_limit)status=3
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    call checked_product([int(nband,int64),int(nband,int64),int(nntot,int64)],&
      output_elements,arithmetic_ok)
    if(arithmetic_ok)call checked_add(output_elements,int(nband,int64)*int(nwann,int64),arithmetic_ok)
    if(arithmetic_ok)call checked_product([output_elements,complex_bytes],output_bytes,arithmetic_ok)
    if(.not.arithmetic_ok)status=2
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));workspace_peak_bytes=0_int64
      if(status==3)message='Wannier90 coordinator byte limit exceeded'
      if(status==2)message='Wannier90 matrix byte estimate overflow'
      return
    endif
    allocation_status=0
    if(rank==0)then
      allocate(m_matrix(nband,nband,nntot),a_matrix(nband,nwann),stat=allocation_status)
      if(allocation_status==0)then;m_matrix=(0d0,0d0);a_matrix=(0d0,0d0);endif
    else
      allocate(m_matrix(0,0,0),a_matrix(0,0),stat=allocation_status)
    endif
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(allocation_status/=0.or.ierr/=MPI_SUCCESS)then
      if(.not.allocated(m_matrix))allocate(m_matrix(0,0,0))
      if(.not.allocated(a_matrix))allocate(a_matrix(0,0))
      message='cannot allocate Wannier90 coordinator matrices';return
    endif
    peak=merge(output_bytes,0_int64,rank==0)
    do b=1,nntot
      do n0=1,nband,tile_size
        n1=min(n0+tile_size-1,nband)
        do m0=1,nband,tile_size
          m1=min(m0+tile_size-1,nband);allocate(local_tile(m1-m0+1,n1-n0+1));local_tile=(0d0,0d0)
          do p=1,npoint
            angle=-2d0*acos(-1d0)*dot_product(real(nncell(:,b),real64),fractional(:,p))
            phase=cmplx(cos(angle),sin(angle),real64)
            do n=n0,n1;do m=m0,m1
              local_tile(m-m0+1,n-n0+1)=local_tile(m-m0+1,n-n0+1)+&
                weights(p)*conjg(values(m,p))*phase*values(n,p)
            enddo;enddo
          enddo
          count=size(local_tile);allocate(reduced_tile(size(local_tile,1),size(local_tile,2)))
          call MPI_Reduce(local_tile,reduced_tile,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
          call checked_product([2_int64,int(count,int64),complex_bytes],tile_bytes,arithmetic_ok)
          if(arithmetic_ok)peak=max(peak,merge(output_bytes,0_int64,rank==0)+tile_bytes)
          if(rank==0.and.ierr==MPI_SUCCESS)m_matrix(m0:m1,n0:n1,b)=reduced_tile
          deallocate(local_tile,reduced_tile);if(ierr/=MPI_SUCCESS)status=4
        enddo
      enddo
    enddo
    do n0=1,nwann,tile_size
      n1=min(n0+tile_size-1,nwann)
      do m0=1,nband,tile_size
        m1=min(m0+tile_size-1,nband);allocate(local_tile(m1-m0+1,n1-n0+1));local_tile=(0d0,0d0)
        do p=1,npoint;do n=n0,n1;do m=m0,m1
          local_tile(m-m0+1,n-n0+1)=local_tile(m-m0+1,n-n0+1)+&
            weights(p)*conjg(values(m,p))*anchors(n,p)
        enddo;enddo;enddo
        count=size(local_tile);allocate(reduced_tile(size(local_tile,1),size(local_tile,2)))
        call MPI_Reduce(local_tile,reduced_tile,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
        call checked_product([2_int64,int(count,int64),complex_bytes],tile_bytes,arithmetic_ok)
        if(arithmetic_ok)peak=max(peak,merge(output_bytes,0_int64,rank==0)+tile_bytes)
        if(rank==0.and.ierr==MPI_SUCCESS)a_matrix(m0:m1,n0:n1)=reduced_tile
        deallocate(local_tile,reduced_tile);if(ierr/=MPI_SUCCESS)status=4
      enddo
    enddo
    workspace_peak_bytes=peak
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=status==0.and.ierr==MPI_SUCCESS
    if(.not.ok)message='Wannier90 distributed matrix reduction failed'
#else
    ok=.false.;message='Wannier90 matrix assembly requires MPI';coordinator_bytes=0_int64
    workspace_peak_bytes=0_int64;allocate(m_matrix(0,0,0),a_matrix(0,0))
#endif
  end subroutine assemble_dg_w90_gamma_matrices

  subroutine setup_dg_w90_gamma_library(comm,seed,real_lattice,reciprocal_lattice,atom_symbols,&
      atoms_cart,nband,nwann,nntot,nncell,ok,message)
    integer,intent(in)::comm,nband,nwann
    character(*),intent(in)::seed
    real(real64),intent(in)::real_lattice(3,3),reciprocal_lattice(3,3),atoms_cart(:,:)
    character(*),intent(in)::atom_symbols(:)
    integer,intent(out)::nntot
    integer,allocatable,intent(out)::nncell(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#if defined(USE_MPI) && defined(USE_WANNIER90)
    integer,parameter::num_nnmax=12
    integer::rank,ierr,unit,io,axis,atom,num_bands_out,num_wann_out,status
    logical::dmn_exists
    integer::mp_grid(3),nnlist(1,num_nnmax),nncell_max(3,1,num_nnmax),exclude_bands(max(1,nband))
    integer::proj_l(max(1,nband)),proj_m(max(1,nband)),proj_radial(max(1,nband))
    integer::proj_s(max(1,nband))
    real(real64)::kpoint(3,1),proj_site(3,max(1,nband)),proj_z(3,max(1,nband)),&
      proj_x(3,max(1,nband)),proj_zona(max(1,nband)),proj_s_qaxis(3,max(1,nband))
    interface
      subroutine wannier_setup(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
          recip_lattice_loc,kpt_latt_loc,num_bands_tot,num_atoms_loc,atom_symbols_loc,&
          atoms_cart_loc,gamma_only_loc,spinors_loc,nntot_loc,nnlist_loc,nncell_loc,&
          num_bands_loc,num_wann_loc,proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,&
          proj_z_loc,proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)
        import real64,num_nnmax
        character(*),intent(in)::seed_name
        integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_tot,num_atoms_loc
        real(real64),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),kpt_latt_loc(3,num_kpts_loc)
        character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
        real(real64),intent(in)::atoms_cart_loc(3,num_atoms_loc)
        logical,intent(in)::gamma_only_loc,spinors_loc
        integer,intent(out)::nntot_loc,nnlist_loc(num_kpts_loc,num_nnmax),&
          nncell_loc(3,num_kpts_loc,num_nnmax),num_bands_loc,num_wann_loc
        real(real64),intent(out)::proj_site_loc(3,num_bands_tot),proj_z_loc(3,num_bands_tot),&
          proj_x_loc(3,num_bands_tot),proj_zona_loc(num_bands_tot)
        integer,intent(out)::proj_l_loc(num_bands_tot),proj_m_loc(num_bands_tot),&
          proj_radial_loc(num_bands_tot),exclude_bands_loc(num_bands_tot),proj_s_loc(num_bands_tot)
        real(real64),intent(out)::proj_s_qaxis_loc(3,num_bands_tot)
      end subroutine wannier_setup
    end interface
    ok=.false.;message='';nntot=0;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.size(atom_symbols)<=0.or.&
        any(shape(atoms_cart)/=[3,size(atom_symbols)]).or..not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or..not.all(ieee_is_finite(atoms_cart)))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma Wannier90 setup contract';return;endif
    mp_grid=[1,1,1];kpoint=0d0
    if(rank==0)then
      inquire(file=trim(seed)//'.dmn',exist=dmn_exists)
      if(.not.dmn_exists)status=4
    endif
    if(rank==0.and.status==0)then
      open(newunit=unit,file=trim(seed)//'.win',status='replace',action='write',iostat=io)
      if(io/=0)then
        status=2
      else
        write(unit,'(a,i0)')'num_bands = ',nband
        write(unit,'(a,i0)')'num_wann = ',nwann
        write(unit,'(a)')'num_iter = 200'
        write(unit,'(a)')'conv_tol = 1.d-12'
        write(unit,'(a)')'conv_window = 5'
        write(unit,'(a)')'gamma_only = true'
        write(unit,'(a)')'site_symmetry = .true.'
        write(unit,'(a)')'symmetrize_eps = 1.d-10'
        write(unit,'(a)')'begin unit_cell_cart';write(unit,'(a)')'bohr'
        do axis=1,3;write(unit,'(3(es24.16,1x))')real_lattice(:,axis);enddo
        write(unit,'(a)')'end unit_cell_cart'
        write(unit,'(a)')'begin atoms_cart';write(unit,'(a)')'bohr'
        do atom=1,size(atom_symbols)
          write(unit,'(a,1x,3(es24.16,1x))')trim(atom_symbols(atom)),atoms_cart(:,atom)
        enddo
        write(unit,'(a)')'end atoms_cart'
        write(unit,'(a)')'begin projections';write(unit,'(a)')'random'
        write(unit,'(a)')'end projections'
        write(unit,'(a)')'mp_grid = 1 1 1'
        write(unit,'(a)')'begin kpoints';write(unit,'(a)')'0.0 0.0 0.0'
        write(unit,'(a)')'end kpoints';close(unit)
        call wannier_setup(trim(seed),mp_grid,1,real_lattice,reciprocal_lattice,kpoint,nband,&
          size(atom_symbols),atom_symbols,atoms_cart,.true.,.false.,nntot,nnlist,nncell_max,&
          num_bands_out,num_wann_out,proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,&
          proj_zona,exclude_bands,proj_s,proj_s_qaxis)
        if(nntot<1.or.nntot>num_nnmax.or.num_bands_out/=nband.or.num_wann_out/=nwann)status=3
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(nntot,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(nncell_max,3*num_nnmax,MPI_INTEGER,0,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma library setup failed';return;endif
    allocate(nncell(3,nntot));nncell=nncell_max(:,1,1:nntot);ok=.true.
#else
    ok=.false.;message='Wannier90 Gamma setup requires MPI and USE_WANNIER90';nntot=0
#endif
  end subroutine setup_dg_w90_gamma_library

  subroutine run_dg_w90_gamma_library(comm,seed,real_lattice,reciprocal_lattice,atom_symbols,&
      atoms_cart,m_matrix,a_matrix,eigenvalues,initial_gauge_spread,tolerance,transform,centers,&
      spreads,spread,ok,message,convergence_iterations_out)
    integer,intent(in)::comm
    character(*),intent(in)::seed
    real(real64),intent(in)::real_lattice(3,3),reciprocal_lattice(3,3),atoms_cart(:,:),&
      eigenvalues(:),initial_gauge_spread,tolerance
    character(*),intent(in)::atom_symbols(:)
    complex(real64),intent(in)::m_matrix(:,:,:),a_matrix(:,:)
    complex(real64),allocatable,intent(out)::transform(:,:)
    real(real64),allocatable,intent(out)::centers(:,:),spreads(:)
    real(real64),intent(out)::spread(3)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(out),optional::convergence_iterations_out
#if defined(USE_MPI) && defined(USE_WANNIER90)
    integer::rank,ierr,nband,nwann,nntot,status,mp_grid(3),matrix_dimensions(3),convergence_iterations
    real(real64)::kpoint(3,1)
    complex(real64),allocatable::u(:,:,:),uopt(:,:,:),m4(:,:,:,:),a3(:,:,:)
    real(real64),allocatable::e2(:,:)
    logical,allocatable::lwindow(:,:)
    interface
      subroutine wannier_run(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
          recip_lattice_loc,kpt_latt_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc,&
          atom_symbols_loc,atoms_cart_loc,gamma_only_loc,m_matrix_loc,a_matrix_loc,&
          eigenvalues_loc,u_matrix_loc,u_matrix_opt_loc,lwindow_loc,wann_centres_loc,&
          wann_spreads_loc,spread_loc)
        import real64
        character(*),intent(in)::seed_name
        integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc
        real(real64),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),kpt_latt_loc(3,num_kpts_loc)
        character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
        real(real64),intent(in)::atoms_cart_loc(3,num_atoms_loc)
        logical,intent(in)::gamma_only_loc
        complex(real64),intent(in)::m_matrix_loc(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc),&
          a_matrix_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
        real(real64),intent(in)::eigenvalues_loc(num_bands_loc,num_kpts_loc)
        complex(real64),intent(out)::u_matrix_loc(num_wann_loc,num_wann_loc,num_kpts_loc),&
          u_matrix_opt_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
        logical,intent(out)::lwindow_loc(num_bands_loc,num_kpts_loc)
        real(real64),intent(out)::wann_centres_loc(3,num_wann_loc),wann_spreads_loc(num_wann_loc),spread_loc(3)
      end subroutine wannier_run
    end interface
    ok=.false.;message='';spread=0d0;status=0;convergence_iterations=-1
    if(present(convergence_iterations_out))convergence_iterations_out=-1
    call MPI_Comm_rank(comm,rank,ierr)
    matrix_dimensions=0
    if(rank==0)matrix_dimensions=[size(m_matrix,1),size(a_matrix,2),size(m_matrix,3)]
    call MPI_Bcast(matrix_dimensions,3,MPI_INTEGER,0,comm,ierr)
    nband=matrix_dimensions(1);nwann=matrix_dimensions(2);nntot=matrix_dimensions(3)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.size(eigenvalues)/=nband.or.nntot<=0.or.&
        any(shape(atoms_cart)/=[3,size(atom_symbols)]).or.&
        .not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or..not.all(ieee_is_finite(atoms_cart)))status=1
    if(rank==0)then
      if(size(m_matrix,2)/=nband.or.size(a_matrix,1)/=nband.or.&
          .not.all(ieee_is_finite(real(m_matrix))).or..not.all(ieee_is_finite(aimag(m_matrix))).or.&
          .not.all(ieee_is_finite(real(a_matrix))).or..not.all(ieee_is_finite(aimag(a_matrix))))status=1
    else if(size(m_matrix)/=0.or.size(a_matrix)/=0)then
      status=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma Wannier90 run contract';return;endif
    allocate(transform(nwann,nwann),centers(3,nwann),spreads(nwann));transform=(0d0,0d0)
    centers=0d0;spreads=0d0;mp_grid=[1,1,1];kpoint=0d0
    if(rank==0)then
      allocate(u(nwann,nwann,1),uopt(nband,nwann,1),lwindow(nband,1),&
        m4(nband,nband,nntot,1),a3(nband,nwann,1),e2(nband,1))
      m4(:,:,:,1)=m_matrix;a3(:,:,1)=a_matrix;e2(:,1)=eigenvalues
      call wannier_run(trim(seed),mp_grid,1,real_lattice,reciprocal_lattice,kpoint,nband,nwann,&
        nntot,size(atom_symbols),atom_symbols,atoms_cart,.true.,m4,a3,e2,u,uopt,lwindow,&
        centers,spreads,spread)
      transform=matmul(uopt(:,:,1),u(:,:,1))
      call validate_dg_w90_convergence_log(trim(seed)//'.wout',200,convergence_iterations,ok,message)
      if(ok)call validate_dg_w90_result(transform,centers,spreads,spread,initial_gauge_spread,&
        tolerance,ok,message)
      status=merge(0,2,ok)
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(convergence_iterations,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(transform,size(transform),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(centers,size(centers),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spreads,size(spreads),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spread,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ok=status==0.and.ierr==MPI_SUCCESS
    if(present(convergence_iterations_out))convergence_iterations_out=convergence_iterations
    if(ok)then;message='';else;message='Wannier90 Gamma library run failed validation';endif
#else
    ok=.false.;message='Wannier90 Gamma run requires MPI and USE_WANNIER90';spread=0d0
    if(present(convergence_iterations_out))convergence_iterations_out=-1
#endif
  end subroutine run_dg_w90_gamma_library

  subroutine estimate_dg_w90_coordinator_bytes(nband,nwann,nntot,nkpoint,nbytes,ok,message)
    integer,intent(in)::nband,nwann,nntot,nkpoint
    integer(int64),intent(out)::nbytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::nb,nw,nn,nk,complex_elements,real_elements,term
    ok=.false.;message='';nbytes=0_int64
    if(nband<=0.or.nwann<=0.or.nwann>nband.or.nntot<=0.or.nkpoint<=0)then
      message='invalid Wannier90 coordinator dimensions';return
    endif
    nb=int(nband,int64);nw=int(nwann,int64);nn=int(nntot,int64);nk=int(nkpoint,int64)
    complex_elements=0_int64;real_elements=0_int64
    ! Input plus the library-owned original/projected overlap copies.
    call checked_product([2_int64,nb,nb,nn,nk],term,ok)
    if(.not.ok)then;message='Wannier90 M-matrix byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    ! A input, optimized subspace, returned U/Uopt, and SALMON scatter copy.
    call checked_product([3_int64,nb,nw,nk],term,ok)
    if(.not.ok)then;message='Wannier90 A/Uopt byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    call checked_product([2_int64,nw,nw,nk],term,ok)
    if(.not.ok)then;message='Wannier90 U-matrix byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    call checked_product([nb,nk],term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,term,ok);if(.not.ok)goto 900
    call checked_product([4_int64,nw],term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,3_int64,ok);if(.not.ok)goto 900
    call checked_product([complex_elements,int(storage_size((0d0,0d0))/8,int64)],term,ok)
    if(.not.ok)goto 900
    nbytes=term
    call checked_product([real_elements,int(storage_size(0d0)/8,int64)],term,ok)
    if(.not.ok)goto 900
    call checked_add(nbytes,term,ok);if(.not.ok)goto 900
    ok=.true.;return
900 message='Wannier90 coordinator byte estimate overflow';nbytes=0_int64
  end subroutine estimate_dg_w90_coordinator_bytes

  subroutine validate_dg_w90_result(transform,centers,spreads,spread,initial_gauge_spread,&
      tolerance,ok,message)
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::centers(:,:),spreads(:),spread(:),initial_gauge_spread,tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:)
    real(real64)::scale,defect,imaginary_defect
    integer::i,nwann
    ok=.false.;message='';nwann=size(transform,1)
    if(nwann<=0.or.size(transform,2)/=nwann.or.any(shape(centers)/=[3,nwann]).or.&
        size(spreads)/=nwann.or.size(spread)/=3.or.tolerance<=0d0.or.&
        .not.ieee_is_finite(tolerance).or..not.ieee_is_finite(initial_gauge_spread).or.&
        initial_gauge_spread<0d0)then
      message='invalid Wannier90 result dimensions or tolerance';return
    endif
    if(.not.all(ieee_is_finite(real(transform))).or..not.all(ieee_is_finite(aimag(transform))).or.&
        .not.all(ieee_is_finite(centers)).or..not.all(ieee_is_finite(spreads)).or.&
        .not.all(ieee_is_finite(spread)))then
      message='Wannier90 result is nonfinite';return
    endif
    scale=max(1d0,max(maxval(abs(spreads)),maxval(abs(spread))))
    if(any(spreads < -tolerance*scale).or.any(spread < -tolerance*scale))then
      message='Wannier90 result has a physically negative spread';return
    endif
    scale=max(1d0,maxval(abs(transform)))
    imaginary_defect=maxval(abs(aimag(transform)))
    if(imaginary_defect>tolerance*scale)then
      message='Wannier90 transform violates the Gamma-real gauge';return
    endif
    allocate(gram(nwann,nwann));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nwann;gram(i,i)=gram(i,i)-1d0;enddo
    defect=maxval(abs(gram))
    if(defect>tolerance*max(1d0,real(nwann,real64)))then
      message='Wannier90 transform is not unitary';return
    endif
    if(spread(3)>initial_gauge_spread+tolerance*max(1d0,initial_gauge_spread))then
      message='Wannier90 increased the gauge-dependent spread';return
    endif
    ok=.true.
  end subroutine validate_dg_w90_result

  subroutine checked_product(factors,value,ok)
    integer(int64),intent(in)::factors(:)
    integer(int64),intent(out)::value
    logical,intent(out)::ok
    integer::i
    value=1_int64;ok=all(factors>=0_int64)
    if(.not.ok)return
    do i=1,size(factors)
      if(factors(i)==0_int64)then;value=0_int64;return;endif
      if(value>huge(value)/factors(i))then;value=0_int64;ok=.false.;return;endif
      value=value*factors(i)
    enddo
  end subroutine checked_product

  subroutine checked_add(value,increment,ok)
    integer(int64),intent(inout)::value
    integer(int64),intent(in)::increment
    logical,intent(inout)::ok
    if(.not.ok)return
    ok=value>=0_int64.and.increment>=0_int64.and.value<=huge(value)-increment
    if(ok)value=value+increment
  end subroutine checked_add
end module dg_overlapping_wannier_w90
