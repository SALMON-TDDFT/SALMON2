module dg_wpw_s_orthogonal_complement
  use mpi,only:MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_Comm_rank,MPI_Comm_size,MPI_INTEGER,MPI_INTEGER8,&
    MPI_MAX,MPI_MIN,&
    MPI_SUM,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUCCESS,MPI_COMM_NULL
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_wpw_owner_exchange,only:s_dg_wpw_owner_schedule,initialize_dg_wpw_owner_schedule,&
    fetch_rows_from_owners,reduce_w_partial_to_owners,release_dg_wpw_owner_schedule
  use dg_wpw_bounded_operator,only:s_dg_wpw_bounded_operator,apply_h_dg_wpw_bounded,&
    apply_s_dg_wpw_bounded,global_gram_dg_wpw_bounded
  implicit none
  private

  type,public::s_dg_wpw_s_orthogonal_complement
    logical::valid=.false.
    integer(8)::metric_fingerprint=0_8
    real(8)::relative_cutoff=0d0,solve_residual=huge(1d0),cross_metric_defect=huge(1d0)
    real(8)::hermitian_defect=huge(1d0),minimum_rayleigh=0d0,condition_estimate=huge(1d0),&
      cutoff_low_weight=0d0
    integer::numerical_p_rank=0,iteration_count=0
    integer,allocatable::global_p_ids(:)
    complex(8),allocatable::a_owned_w_global_p(:,:)
    type(s_dg_wpw_owner_schedule)::global_p_schedule
  end type

  public::initialize_dg_wpw_s_orthogonal_complement
  public::apply_h_dg_wpw_s_orthogonal_complement,apply_s_dg_wpw_s_orthogonal_complement
  public::map_dg_wpw_complement_to_original,map_dg_wpw_original_to_complement
  public::validate_dg_wpw_s_orthogonal_complement,release_dg_wpw_s_orthogonal_complement
contains

  subroutine initialize_dg_wpw_s_orthogonal_complement(op,relative_cutoff,transform,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    real(8),intent(in)::relative_cutoff
    type(s_dg_wpw_s_orthogonal_complement),intent(inout)::transform
    integer,intent(out)::info
    type(s_dg_wpw_s_orthogonal_complement)::candidate
    integer,allocatable::counts(:),displs(:),all_p(:),all_w(:),all_peers(:)
    complex(8),allocatable::rhs(:,:),solution(:,:),partial(:,:),zero_p(:,:)
    real(8),allocatable::rhs_norm(:),residual(:)
    real(8)::max_residual,max_cross,hermitian_defect,global_min,global_max
    integer::comm,nrank,rank,i,j,k,kend,nblock,block_width,nall,n_global_w,astat,local_bad,global_bad,ierr,&
      solve_info,max_iter
    integer::powner,iterations,total_iterations

    info=1;local_bad=0;comm=op%w_schedule%comm
    if(comm==MPI_COMM_NULL)return
    if(.not.op%valid.or..not.ieee_is_finite(relative_cutoff).or.relative_cutoff<=0d0.or.&
      relative_cutoff>=1d0.or..not.allocated(op%owned_w_ids).or..not.allocated(op%owned_p_ids).or.&
      .not.allocated(op%required_w_ids).or..not.allocated(op%required_p_ids).or.&
      .not.allocated(op%ww_s_dense).or..not.allocated(op%wp_s_dense).or.&
      .not.allocated(op%pp_s_dense).or..not.finite2(op%ww_s_dense).or.&
      .not.finite2(op%wp_s_dense).or..not.finite2(op%pp_s_dense))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(size(op%owned_w_ids),n_global_w,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(size(op%owned_p_ids),nall,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.n_global_w<=0.or.nall<=0)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(counts(nrank),displs(nrank),all_p(nall),all_peers(max(0,nrank-1)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allgather(size(op%owned_p_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    displs(1)=0
    do i=2,nrank;displs(i)=displs(i-1)+counts(i-1);enddo
    call MPI_Allgatherv(op%owned_p_ids,size(op%owned_p_ids),MPI_INTEGER,all_p,counts,displs,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call sort_int(all_p)
    local_bad=merge(0,1,strictly_increasing(all_p))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call gather_global_ids(comm,op%owned_w_ids,all_w,solve_info)
    if(solve_info/=0.or.size(all_w)/=n_global_w)return
    call diagnose_sww(op,all_w,relative_cutoff,hermitian_defect,global_min,global_max,solve_info)
    local_bad=merge(0,1,solve_info==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    candidate%global_p_ids=all_p
    j=0
    do i=0,nrank-1
      if(i/=rank)then;j=j+1;all_peers(j)=i;endif
    enddo
    call initialize_dg_wpw_owner_schedule(candidate%global_p_schedule,comm,all_peers,&
      op%owned_p_ids,candidate%global_p_ids,solve_info)
    local_bad=merge(0,1,solve_info==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%global_p_schedule);return
    endif
    block_width=min(8,nall)
    allocate(candidate%a_owned_w_global_p(size(op%owned_w_ids),nall),rhs(size(op%owned_w_ids),block_width),&
      solution(size(op%owned_w_ids),block_width),partial(size(op%required_w_ids),block_width),&
      zero_p(size(op%owned_p_ids),block_width),rhs_norm(block_width),residual(block_width),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%global_p_schedule);return
    endif
    candidate%a_owned_w_global_p=0;zero_p=0
    max_residual=0d0;max_cross=0d0;total_iterations=0;max_iter=max(200,20*n_global_w)
    do k=1,nall,block_width
      kend=min(nall,k+block_width-1);nblock=kend-k+1;partial=0;rhs=0;solution=0;zero_p=0
      do j=1,nblock
        powner=find_id(op%owned_p_ids,all_p(k+j-1))
        if(powner>0)partial(:,j)=op%wp_s_dense(:,powner)
      enddo
      call reduce_w_partial_to_owners(op%w_schedule,partial(:,1:nblock),rhs(:,1:nblock),solve_info)
      if(solve_info/=0)then;local_bad=1;exit;endif
      call solve_sww_deflated_block(op,rhs(:,1:nblock),zero_p(:,1:nblock),relative_cutoff,max_iter,&
        solution(:,1:nblock),rhs_norm(1:nblock),residual(1:nblock),iterations,solve_info)
      if(solve_info/=0)then
        if(rank==0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-PW-COMPLEMENT-FAIL] block_pcg first_p=',k,&
          ' info=',solve_info
        local_bad=1;exit
      endif
      candidate%a_owned_w_global_p(:,k:kend)=solution(:,1:nblock)
      max_residual=max(max_residual,maxval(residual(1:nblock)/max(1d-300,rhs_norm(1:nblock))))
      max_cross=max(max_cross,maxval(residual(1:nblock)/max(1d-300,rhs_norm(1:nblock))))
      total_iterations=max(total_iterations,iterations)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%global_p_schedule);return
    endif
    call diagnose_sperp(op,candidate%global_p_schedule,candidate%global_p_ids,&
      candidate%a_owned_w_global_p,relative_cutoff,candidate%numerical_p_rank,&
      candidate%cutoff_low_weight,solve_info)
    local_bad=merge(0,1,solve_info==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%global_p_schedule);return
    endif
    candidate%metric_fingerprint=compute_metric_fingerprint(op)
    candidate%relative_cutoff=relative_cutoff;candidate%solve_residual=max_residual
    candidate%cross_metric_defect=max_cross;candidate%minimum_rayleigh=global_min
    candidate%hermitian_defect=hermitian_defect
    candidate%condition_estimate=global_max/global_min;candidate%iteration_count=total_iterations
    local_bad=merge(0,1,candidate%metric_fingerprint/=0_8.and.max_residual<=1d-11.and.max_cross<=1d-11)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call release_dg_wpw_owner_schedule(candidate%global_p_schedule);return
    endif
    candidate%valid=.true.
    call move_transform(candidate,transform);info=0
  end subroutine initialize_dg_wpw_s_orthogonal_complement

  subroutine solve_sww_deflated_block(op,b,zero_p,tolerance,max_iter,x,bnorm,residual,iterations,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    complex(8),intent(in)::b(:,:),zero_p(:,:);real(8),intent(in)::tolerance;integer,intent(in)::max_iter
    complex(8),intent(out)::x(:,:);real(8),intent(out)::bnorm(:),residual(:)
    integer,intent(out)::iterations,info
    complex(8),allocatable::gram(:,:),eigenvectors(:,:),work(:),active_b(:,:),active_zero_p(:,:),&
      active_x(:,:),rw(:,:),partial(:,:),owned(:,:),ax(:,:),r(:,:),active_coeff(:,:)
    complex(8)::work_query(1);real(8),allocatable::eigenvalues(:),rwork(:),active_bnorm(:),active_residual(:)
    integer::nrhs,nactive,first_active,i,astat,local_bad,global_bad,ierr,gram_info,lwork,lapack_info
    external::zheev
    info=1;x=0;bnorm=0;residual=huge(1d0);iterations=0;nrhs=size(b,2);local_bad=0
    allocate(gram(nrhs,nrhs),eigenvectors(nrhs,nrhs),eigenvalues(nrhs),rwork(max(1,3*nrhs-2)),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call global_gram_dg_wpw_bounded(op,b,zero_p,b,zero_p,gram,gram_info);if(gram_info/=0)return
    do i=1,nrhs;bnorm(i)=sqrt(max(0d0,real(gram(i,i),8)));enddo
    eigenvectors=gram;lwork=-1
    call zheev('V','U',nrhs,eigenvectors,nrhs,eigenvalues,work_query,lwork,rwork,lapack_info)
    if(lapack_info/=0)return;lwork=max(1,int(real(work_query(1),8)));allocate(work(lwork),stat=astat)
    if(astat/=0)return
    call zheev('V','U',nrhs,eigenvectors,nrhs,eigenvalues,work,lwork,rwork,lapack_info)
    if(lapack_info/=0.or..not.all(ieee_is_finite(eigenvalues)))return
    first_active=1
    do while(first_active<=nrhs.and.eigenvalues(first_active)<=&
      100d0*epsilon(1d0)*max(1d-300,eigenvalues(nrhs)))
      first_active=first_active+1
    enddo
    nactive=nrhs-first_active+1
    if(nactive==0)then;residual=0;info=0;return;endif
    allocate(active_b(size(b,1),nactive),active_zero_p(size(zero_p,1),nactive),&
      active_x(size(x,1),nactive),active_bnorm(nactive),active_residual(nactive),&
      active_coeff(nrhs,nactive),rw(size(op%required_w_ids),nrhs),&
      partial(size(op%required_w_ids),nrhs),owned(size(op%owned_w_ids),nrhs),&
      ax(size(op%owned_w_ids),nrhs),r(size(op%owned_w_ids),nrhs),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    active_coeff=eigenvectors(:,first_active:nrhs)
    active_b=matmul(b,active_coeff);active_zero_p=0
    call solve_sww_block_pcg(op,active_b,active_zero_p,tolerance,max_iter,active_x,&
      active_bnorm,active_residual,iterations,info)
    if(info/=0)return
    x=matmul(active_x,conjg(transpose(active_coeff)))
    call apply_sww(op,x,rw,partial,owned,ax,gram_info);if(gram_info/=0)then;info=1;return;endif
    r=ax-b;call global_gram_dg_wpw_bounded(op,r,zero_p,r,zero_p,gram,gram_info)
    if(gram_info/=0)then;info=1;return;endif
    do i=1,nrhs;residual(i)=sqrt(max(0d0,real(gram(i,i),8)));enddo
    info=0
  end subroutine solve_sww_deflated_block

  subroutine solve_sww_block_pcg(op,b,zero_p,tolerance,max_iter,x,bnorm,residual,iterations,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    complex(8),intent(in)::b(:,:),zero_p(:,:)
    real(8),intent(in)::tolerance;integer,intent(in)::max_iter
    complex(8),intent(out)::x(:,:);real(8),intent(out)::bnorm(:),residual(:)
    integer,intent(out)::iterations,info
    complex(8),allocatable::r(:,:),z(:,:),p(:,:),ap(:,:),rw(:,:),partial(:,:),owned(:,:)
    complex(8),allocatable::rz(:,:),rz_new(:,:),pap(:,:),rr(:,:),alpha(:,:),beta(:,:)
    real(8),allocatable::diag(:)
    integer::i,nrhs,pos,astat,local_bad,global_bad,ierr,apply_info,solve_info
    info=1;x=0;residual=huge(1d0);bnorm=0d0;iterations=0;local_bad=0;nrhs=size(b,2)
    allocate(r(size(b,1),nrhs),z(size(b,1),nrhs),p(size(b,1),nrhs),ap(size(b,1),nrhs),&
      rw(size(op%required_w_ids),nrhs),partial(size(op%required_w_ids),nrhs),owned(size(b,1),nrhs),&
      diag(size(b,1)),rz(nrhs,nrhs),rz_new(nrhs,nrhs),pap(nrhs,nrhs),rr(nrhs,nrhs),&
      alpha(nrhs,nrhs),beta(nrhs,nrhs),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    do i=1,size(diag)
      pos=find_id(op%required_w_ids,op%owned_w_ids(i));diag(i)=0d0
      if(pos>0)diag(i)=real(op%ww_s_dense(pos,pos),8)
      if(.not.ieee_is_finite(diag(i)).or.diag(i)<=0d0)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    r=b;z=r/spread(diag,2,nrhs);p=z
    call global_gram_dg_wpw_bounded(op,r,zero_p,z,zero_p,rz,apply_info);if(apply_info/=0)return
    call global_gram_dg_wpw_bounded(op,b,zero_p,b,zero_p,rr,apply_info);if(apply_info/=0)return
    do i=1,nrhs;bnorm(i)=sqrt(max(1d-300,real(rr(i,i),8)));enddo
    do iterations=1,max_iter
      call apply_sww(op,p,rw,partial,owned,ap,apply_info);if(apply_info/=0)return
      call global_gram_dg_wpw_bounded(op,p,zero_p,ap,zero_p,pap,apply_info);if(apply_info/=0)return
      alpha=rz;call solve_hermitian_rank_revealing(pap,alpha,solve_info);if(solve_info/=0)return
      x=x+matmul(p,alpha);r=r-matmul(ap,alpha)
      call global_gram_dg_wpw_bounded(op,r,zero_p,r,zero_p,rr,apply_info);if(apply_info/=0)return
      do i=1,nrhs;residual(i)=sqrt(max(0d0,real(rr(i,i),8)));enddo
      if(all(residual<=min(1d-12,tolerance)*bnorm))then;info=0;return;endif
      z=r/spread(diag,2,nrhs)
      call global_gram_dg_wpw_bounded(op,r,zero_p,z,zero_p,rz_new,apply_info);if(apply_info/=0)return
      beta=rz_new;call solve_hermitian_rank_revealing(rz,beta,solve_info);if(solve_info/=0)return
      p=z+matmul(p,beta);rz=rz_new
    enddo
  end subroutine solve_sww_block_pcg

  subroutine solve_hermitian_rank_revealing(matrix,rhs,info)
    complex(8),intent(in)::matrix(:,:);complex(8),intent(inout)::rhs(:,:);integer,intent(out)::info
    complex(8),allocatable::eigenvectors(:,:),work(:),spectral(:,:);real(8),allocatable::eigenvalues(:),rwork(:)
    complex(8)::work_query(1);real(8)::threshold;integer::n,i,astat,lwork,lapack_info
    external::zheev
    info=1;n=size(matrix,1)
    if(size(matrix,2)/=n.or.size(rhs,1)/=n)return
    allocate(eigenvectors(n,n),eigenvalues(n),rwork(max(1,3*n-2)),spectral(n,size(rhs,2)),stat=astat)
    if(astat/=0)return;eigenvectors=0.5d0*(matrix+conjg(transpose(matrix)))
    lwork=-1;call zheev('V','U',n,eigenvectors,n,eigenvalues,work_query,lwork,rwork,lapack_info)
    if(lapack_info/=0)return;lwork=max(1,int(real(work_query(1),8)));allocate(work(lwork),stat=astat)
    if(astat/=0)return
    call zheev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,rwork,lapack_info)
    if(lapack_info/=0.or..not.all(ieee_is_finite(eigenvalues)))return
    threshold=100d0*epsilon(1d0)*max(1d-300,maxval(abs(eigenvalues)))
    if(any(eigenvalues < -threshold))return
    spectral=matmul(conjg(transpose(eigenvectors)),rhs)
    do i=1,n
      if(eigenvalues(i)>threshold)then;spectral(i,:)=spectral(i,:)/eigenvalues(i)
      else;spectral(i,:)=0;endif
    enddo
    rhs=matmul(eigenvectors,spectral);info=0
  end subroutine solve_hermitian_rank_revealing

  subroutine apply_sww(op,x,rw,partial,owned,y,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    complex(8),intent(in)::x(:,:);complex(8),intent(out)::rw(:,:),partial(:,:),owned(:,:),y(:,:)
    integer,intent(out)::info
    call fetch_rows_from_owners(op%w_schedule,x,rw,info);if(info/=0)return
    partial=matmul(op%ww_s_dense,rw)
    call reduce_w_partial_to_owners(op%w_schedule,partial,owned,info);if(info/=0)return
    y=owned
  end subroutine apply_sww

  subroutine apply_sperp(op,p_schedule,a,xp,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_owner_schedule),intent(in)::p_schedule
    complex(8),intent(in)::a(:,:),xp(:,:);complex(8),intent(out)::yp(:,:);integer,intent(out)::info
    complex(8),allocatable::global_xp(:,:),ow(:,:),bw(:,:),bp(:,:),partial(:,:),correction(:,:)
    integer::nrhs,astat,local_bad,global_bad,ierr
    info=1;yp=0;nrhs=size(xp,2);local_bad=0
    if(size(xp,1)/=size(op%owned_p_ids).or.size(a,1)/=size(op%owned_w_ids).or.&
      size(a,2)/=size(p_schedule%required_ids).or.any(shape(yp)/=shape(xp)).or..not.finite2(xp))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(global_xp(size(p_schedule%required_ids),nrhs),ow(size(op%owned_w_ids),nrhs),&
      bw(size(op%owned_w_ids),nrhs),bp(size(op%owned_p_ids),nrhs),&
      partial(size(p_schedule%required_ids),nrhs),correction(size(op%owned_p_ids),nrhs),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call fetch_rows_from_owners(p_schedule,xp,global_xp,info);if(info/=0)return
    ow=-matmul(a,global_xp)
    call apply_s_dg_wpw_bounded(op,op%operator_epoch,op%layout_fingerprint,ow,xp,bw,bp,info)
    if(info/=0)return
    partial=matmul(conjg(transpose(a)),bw)
    call reduce_w_partial_to_owners(p_schedule,partial,correction,info);if(info/=0)return
    yp=bp-correction;local_bad=merge(0,1,finite2(yp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine apply_sperp

  subroutine diagnose_sperp(op,p_schedule,global_p_ids,a,relative_cutoff,numerical_rank,low_weight,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_owner_schedule),intent(in)::p_schedule
    integer,intent(in)::global_p_ids(:);complex(8),intent(in)::a(:,:);real(8),intent(in)::relative_cutoff
    integer,intent(out)::numerical_rank,info;real(8),intent(out)::low_weight
    complex(8),allocatable::xp(:,:),yp(:,:),factor(:,:),pivot_row(:),pivot_row_local(:)
    real(8),allocatable::diagonal(:),initial_diagonal(:)
    real(8)::local_max,global_max,initial_max,local_trace,global_trace,local_low,global_low,pivot_value
    complex(8)::correction
    integer::i,k,pos,pivot_id,local_pivot,astat,local_bad,global_bad,ierr,apply_info
    info=1;numerical_rank=0;low_weight=0d0;local_bad=0
    allocate(xp(size(op%owned_p_ids),1),yp(size(op%owned_p_ids),1),&
      factor(size(op%owned_p_ids),size(global_p_ids)),diagonal(size(op%owned_p_ids)),&
      initial_diagonal(size(op%owned_p_ids)),pivot_row(size(global_p_ids)),&
      pivot_row_local(size(global_p_ids)),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    factor=0;diagonal=0
    do k=1,size(global_p_ids)
      xp=0;pos=find_id(op%owned_p_ids,global_p_ids(k));if(pos>0)xp(pos,1)=1
      call apply_sperp(op,p_schedule,a,xp,yp,apply_info);if(apply_info/=0)then;local_bad=1;exit;endif
      if(pos>0)then
        if(abs(aimag(yp(pos,1)))>1d-11*max(1d0,abs(yp(pos,1))))local_bad=1
        diagonal(pos)=real(yp(pos,1),8)
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    initial_diagonal=diagonal;local_max=maxval(diagonal);local_trace=sum(max(0d0,diagonal))
    call MPI_Allreduce(local_max,initial_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,op%w_schedule%comm,ierr)
    call MPI_Allreduce(local_trace,global_trace,1,MPI_DOUBLE_PRECISION,MPI_SUM,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.initial_max<=0d0.or.global_trace<=0d0)return
    do k=1,size(global_p_ids)
      local_max=maxval(diagonal)
      call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,op%w_schedule%comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(global_max<=relative_cutoff*initial_max)exit
      local_pivot=huge(0)
      do i=1,size(op%owned_p_ids)
        if(abs(diagonal(i)-global_max)<=1d-14*max(1d0,global_max))local_pivot=min(local_pivot,op%owned_p_ids(i))
      enddo
      call MPI_Allreduce(local_pivot,pivot_id,1,MPI_INTEGER,MPI_MIN,op%w_schedule%comm,ierr)
      if(ierr/=MPI_SUCCESS.or.pivot_id==huge(0))return
      pivot_row=0;pivot_row_local=0;pos=find_id(op%owned_p_ids,pivot_id)
      if(pos>0.and.k>1)pivot_row_local(1:k-1)=factor(pos,1:k-1)
      if(k>1)then
        call MPI_Allreduce(pivot_row_local,pivot_row,k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,op%w_schedule%comm,ierr)
        if(ierr/=MPI_SUCCESS)return
      endif
      pivot_value=sqrt(global_max);xp=0;if(pos>0)xp(pos,1)=1
      call apply_sperp(op,p_schedule,a,xp,yp,apply_info);if(apply_info/=0)return
      do i=1,size(op%owned_p_ids)
        correction=0
        if(k>1)correction=sum(factor(i,1:k-1)*conjg(pivot_row(1:k-1)))
        factor(i,k)=(yp(i,1)-correction)/pivot_value
        diagonal(i)=diagonal(i)-abs(factor(i,k))**2
        if(diagonal(i)<-1d-10*initial_max)local_bad=1
        diagonal(i)=max(0d0,diagonal(i))
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      numerical_rank=k
    enddo
    local_low=sum(max(0d0,diagonal));call MPI_Allreduce(local_low,global_low,1,MPI_DOUBLE_PRECISION,&
      MPI_SUM,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    low_weight=global_low/global_trace
    if(.not.ieee_is_finite(low_weight))return
    info=0
  end subroutine diagnose_sperp

  subroutine map_dg_wpw_complement_to_original(op,transform,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::info
    complex(8),allocatable::global_xp(:,:)
    integer::nvec,astat,local_bad,global_bad,ierr
    yw=0;yp=0;info=1;nvec=size(xw,2);local_bad=0
    call validate_shapes(op,transform,xw,xp,yw,yp,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)return
    call collective_nvec(op%w_schedule%comm,nvec,info);if(info/=0)return
    allocate(global_xp(size(transform%global_p_ids),nvec),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call fetch_rows_from_owners(transform%global_p_schedule,xp,global_xp,info);if(info/=0)return
    yw=xw-matmul(transform%a_owned_w_global_p,global_xp);yp=xp
    local_bad=merge(0,1,finite2(yw).and.finite2(yp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine map_dg_wpw_complement_to_original

  subroutine map_dg_wpw_original_to_complement(op,transform,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::info
    complex(8),allocatable::global_xp(:,:)
    integer::nvec,astat,local_bad,global_bad,ierr
    yw=0;yp=0;info=1;nvec=size(xw,2);local_bad=0
    call validate_shapes(op,transform,xw,xp,yw,yp,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)return
    call collective_nvec(op%w_schedule%comm,nvec,info);if(info/=0)return
    allocate(global_xp(size(transform%global_p_ids),nvec),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call fetch_rows_from_owners(transform%global_p_schedule,xp,global_xp,info);if(info/=0)return
    yw=xw+matmul(transform%a_owned_w_global_p,global_xp);yp=xp
    local_bad=merge(0,1,finite2(yw).and.finite2(yp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine map_dg_wpw_original_to_complement

  subroutine apply_h_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:);integer,intent(out)::info
    call apply_transformed(op,transform,.true.,xw,xp,yw,yp,info)
  end subroutine apply_h_dg_wpw_s_orthogonal_complement

  subroutine apply_s_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:);integer,intent(out)::info
    call apply_transformed(op,transform,.false.,xw,xp,yw,yp,info)
  end subroutine apply_s_dg_wpw_s_orthogonal_complement

  subroutine apply_transformed(op,transform,use_h,xw,xp,yw,yp,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform;logical,intent(in)::use_h
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:);integer,intent(out)::info
    complex(8),allocatable::ow(:,:),opv(:,:),bw(:,:),bp(:,:),partial(:,:),correction(:,:)
    integer::nvec,astat,local_bad,global_bad,ierr
    yw=0;yp=0;info=1;nvec=size(xw,2);local_bad=0
    if(op%w_schedule%comm==MPI_COMM_NULL)return
    call validate_shapes(op,transform,xw,xp,yw,yp,local_bad)
    if(.not.allocated(transform%global_p_ids).or.&
      .not.allocated(transform%a_owned_w_global_p))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call validate_dg_wpw_s_orthogonal_complement(op,transform,info);if(info/=0)return
    call collective_nvec(op%w_schedule%comm,nvec,info);if(info/=0)return
    allocate(ow(size(xw,1),nvec),opv(size(xp,1),nvec),bw(size(xw,1),nvec),bp(size(xp,1),nvec),&
      partial(size(transform%global_p_ids),nvec),correction(size(op%owned_p_ids),nvec),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call map_dg_wpw_complement_to_original(op,transform,xw,xp,ow,opv,info);if(info/=0)return
    if(use_h)then
      call apply_h_dg_wpw_bounded(op,op%operator_epoch,op%layout_fingerprint,ow,opv,bw,bp,info)
    else
      call apply_s_dg_wpw_bounded(op,op%operator_epoch,op%layout_fingerprint,ow,opv,bw,bp,info)
    endif
    if(info/=0)return
    partial=matmul(conjg(transpose(transform%a_owned_w_global_p)),bw)
    call reduce_w_partial_to_owners(transform%global_p_schedule,partial,correction,info);if(info/=0)return
    yw=bw;yp=bp-correction
    local_bad=merge(0,1,finite2(yw).and.finite2(yp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine apply_transformed

  subroutine validate_dg_wpw_s_orthogonal_complement(op,transform,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform;integer,intent(out)::info
    integer::local_bad,global_bad,ierr;integer(8)::current_fingerprint
    local_bad=0;info=1
    current_fingerprint=compute_metric_fingerprint(op)
    if(.not.op%valid.or..not.transform%valid.or..not.allocated(transform%global_p_ids).or.&
      .not.allocated(transform%a_owned_w_global_p))local_bad=1
    if(local_bad==0)then
      if(any(shape(transform%a_owned_w_global_p)/=[size(op%owned_w_ids),size(transform%global_p_ids)]).or.&
        .not.finite2(transform%a_owned_w_global_p).or.&
        transform%metric_fingerprint/=current_fingerprint)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,op%w_schedule%comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine validate_dg_wpw_s_orthogonal_complement

  subroutine release_dg_wpw_s_orthogonal_complement(transform)
    type(s_dg_wpw_s_orthogonal_complement),intent(inout)::transform
    call release_dg_wpw_owner_schedule(transform%global_p_schedule)
    if(allocated(transform%global_p_ids))deallocate(transform%global_p_ids)
    if(allocated(transform%a_owned_w_global_p))deallocate(transform%a_owned_w_global_p)
    transform=s_dg_wpw_s_orthogonal_complement()
  end subroutine release_dg_wpw_s_orthogonal_complement

  subroutine move_transform(source,destination)
    type(s_dg_wpw_s_orthogonal_complement),intent(inout)::source,destination
    call release_dg_wpw_s_orthogonal_complement(destination)
    destination%valid=source%valid;destination%metric_fingerprint=source%metric_fingerprint
    destination%relative_cutoff=source%relative_cutoff;destination%solve_residual=source%solve_residual
    destination%cross_metric_defect=source%cross_metric_defect
    destination%hermitian_defect=source%hermitian_defect
    destination%minimum_rayleigh=source%minimum_rayleigh
    destination%condition_estimate=source%condition_estimate
    destination%cutoff_low_weight=source%cutoff_low_weight
    destination%numerical_p_rank=source%numerical_p_rank
    destination%iteration_count=source%iteration_count
    destination%global_p_schedule=source%global_p_schedule
    source%global_p_schedule=s_dg_wpw_owner_schedule()
    call move_alloc(source%global_p_ids,destination%global_p_ids)
    call move_alloc(source%a_owned_w_global_p,destination%a_owned_w_global_p)
    source%valid=.false.
  end subroutine move_transform

  subroutine validate_shapes(op,transform,xw,xp,yw,yp,bad)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    type(s_dg_wpw_s_orthogonal_complement),intent(in)::transform
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(in)::yw(:,:),yp(:,:);integer,intent(out)::bad
    bad=0
    if(.not.transform%valid.or.size(xw,1)/=size(op%owned_w_ids).or.&
      size(xp,1)/=size(op%owned_p_ids).or.size(xp,2)/=size(xw,2).or.&
      any(shape(yw)/=shape(xw)).or.any(shape(yp)/=shape(xp)).or.&
      .not.finite2(xw).or..not.finite2(xp))bad=1
  end subroutine validate_shapes

  subroutine collective_nvec(comm,nvec,info)
    integer,intent(in)::comm,nvec;integer,intent(out)::info
    integer::lo,hi,ierr
    info=1;call MPI_Allreduce(nvec,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nvec,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.lo==hi)info=0
  end subroutine collective_nvec

  subroutine gather_global_ids(comm,owned_ids,global_ids,info)
    integer,intent(in)::comm,owned_ids(:);integer,allocatable,intent(out)::global_ids(:);integer,intent(out)::info
    integer,allocatable::counts(:),displs(:);integer::nrank,nall,i,astat,ierr,local_bad,global_bad
    info=1;call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(counts(nrank),displs(nrank),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allgather(size(owned_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    displs(1)=0
    do i=2,nrank;displs(i)=displs(i-1)+counts(i-1);enddo
    nall=sum(counts);allocate(global_ids(nall),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allgatherv(owned_ids,size(owned_ids),MPI_INTEGER,global_ids,counts,displs,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call sort_int(global_ids);local_bad=merge(0,1,strictly_increasing(global_ids))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine gather_global_ids

  subroutine diagnose_sww(op,global_w_ids,relative_cutoff,hermitian_defect,min_rayleigh,max_rayleigh,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer,intent(in)::global_w_ids(:);real(8),intent(in)::relative_cutoff
    real(8),intent(out)::hermitian_defect,min_rayleigh,max_rayleigh;integer,intent(out)::info
    complex(8),allocatable::basis(:,:),y(:,:),z(:,:),rw(:,:),partial(:,:),owned(:,:),zero_p(:,:)
    complex(8)::norm_gram(1,1),trace_gram(1,1),trace_s2
    real(8)::norm_s2,defect2,probe_min,probe_max
    integer::k,pos,seed,astat,local_bad,global_bad,ierr,apply_info
    info=1;hermitian_defect=huge(1d0);min_rayleigh=huge(1d0);max_rayleigh=0d0
    allocate(basis(size(op%owned_w_ids),1),y(size(op%owned_w_ids),1),z(size(op%owned_w_ids),1),&
      rw(size(op%required_w_ids),1),partial(size(op%required_w_ids),1),owned(size(op%owned_w_ids),1),&
      zero_p(size(op%owned_p_ids),1),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    zero_p=0;norm_s2=0d0;trace_s2=0
    do k=1,size(global_w_ids)
      basis=0;pos=find_id(op%owned_w_ids,global_w_ids(k));if(pos>0)basis(pos,1)=1
      call apply_sww(op,basis,rw,partial,owned,y,apply_info);if(apply_info/=0)return
      call global_gram_dg_wpw_bounded(op,y,zero_p,y,zero_p,norm_gram,apply_info);if(apply_info/=0)return
      call apply_sww(op,y,rw,partial,owned,z,apply_info);if(apply_info/=0)return
      call global_gram_dg_wpw_bounded(op,basis,zero_p,z,zero_p,trace_gram,apply_info);if(apply_info/=0)return
      norm_s2=norm_s2+real(norm_gram(1,1),8);trace_s2=trace_s2+trace_gram(1,1)
    enddo
    defect2=max(0d0,2d0*norm_s2-2d0*real(trace_s2,8))
    hermitian_defect=sqrt(defect2)/max(1d-300,sqrt(norm_s2))
    if(.not.ieee_is_finite(hermitian_defect).or.hermitian_defect>relative_cutoff)return
    do seed=1,2
      call lanczos_sww(op,global_w_ids,seed,probe_min,probe_max,apply_info)
      if(apply_info/=0)return
      min_rayleigh=min(min_rayleigh,probe_min);max_rayleigh=max(max_rayleigh,probe_max)
    enddo
    if(.not.ieee_is_finite(min_rayleigh).or..not.ieee_is_finite(max_rayleigh).or.&
      max_rayleigh<=0d0.or.min_rayleigh<=relative_cutoff*max_rayleigh)return
    info=0
  end subroutine diagnose_sww

  subroutine lanczos_sww(op,global_w_ids,seed,min_eval,max_eval,info)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer,intent(in)::global_w_ids(:),seed
    real(8),intent(out)::min_eval,max_eval;integer,intent(out)::info
    complex(8),allocatable::q(:,:),z(:,:),rw(:,:),partial(:,:),owned(:,:),zero_p(:,:),tri(:,:),work(:)
    complex(8)::gram(1,1),work_query(1);real(8),allocatable::eval(:),rwork(:)
    real(8)::alpha,beta,normq;integer::i,j,k,pos,m,mactual,astat,local_bad,global_bad,ierr,apply_info,lwork,lapack_info
    external::zheev
    info=1;min_eval=0d0;max_eval=0d0;m=size(global_w_ids)
    allocate(q(size(op%owned_w_ids),m),z(size(op%owned_w_ids),1),rw(size(op%required_w_ids),1),&
      partial(size(op%required_w_ids),1),owned(size(op%owned_w_ids),1),zero_p(size(op%owned_p_ids),1),&
      tri(m,m),eval(m),rwork(max(1,3*m-2)),stat=astat)
    local_bad=merge(0,1,astat==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,&
      op%w_schedule%comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    zero_p=0;q=0;tri=0
    do i=1,size(op%owned_w_ids)
      q(i,1)=cmplx(sin(real((seed+2)*op%owned_w_ids(i),8)),cos(real((2*seed+1)*op%owned_w_ids(i),8)),8)
    enddo
    call global_gram_dg_wpw_bounded(op,q(:,1:1),zero_p,q(:,1:1),zero_p,gram,apply_info)
    if(apply_info/=0)return
    normq=sqrt(max(0d0,real(gram(1,1),8)));if(normq<=1d-300)return;q(:,1)=q(:,1)/normq
    beta=0d0;mactual=m
    do j=1,m
      call apply_sww(op,q(:,j:j),rw,partial,owned,z,apply_info);if(apply_info/=0)return
      call global_gram_dg_wpw_bounded(op,q(:,j:j),zero_p,z,zero_p,gram,apply_info);if(apply_info/=0)return
      alpha=real(gram(1,1),8);tri(j,j)=alpha
      z=z-alpha*q(:,j:j)
      if(j>1)z=z-beta*q(:,j-1:j-1)
      do i=1,j
        call global_gram_dg_wpw_bounded(op,q(:,i:i),zero_p,z,zero_p,gram,apply_info)
        if(apply_info/=0)return
        z=z-q(:,i:i)*gram(1,1)
      enddo
      if(j==m)exit
      call global_gram_dg_wpw_bounded(op,z,zero_p,z,zero_p,gram,apply_info);if(apply_info/=0)return
      beta=sqrt(max(0d0,real(gram(1,1),8)))
      if(beta<=100d0*epsilon(1d0)*max(1d0,abs(alpha)))then
        beta=0d0
        restart_basis: do k=1,size(global_w_ids)
          z=0;pos=find_id(op%owned_w_ids,global_w_ids(k));if(pos>0)z(pos,1)=1
          do i=1,j
            call global_gram_dg_wpw_bounded(op,q(:,i:i),zero_p,z,zero_p,gram,apply_info)
            if(apply_info/=0)return
            z=z-q(:,i:i)*gram(1,1)
          enddo
          call global_gram_dg_wpw_bounded(op,z,zero_p,z,zero_p,gram,apply_info);if(apply_info/=0)return
          normq=sqrt(max(0d0,real(gram(1,1),8)))
          if(normq>100d0*epsilon(1d0))then;q(:,j+1)=z(:,1)/normq;exit restart_basis;endif
        enddo restart_basis
        if(normq<=100d0*epsilon(1d0))then;mactual=j;exit;endif
      else
        tri(j,j+1)=beta;tri(j+1,j)=beta;q(:,j+1)=z(:,1)/beta
      endif
    enddo
    lwork=-1;call zheev('N','U',mactual,tri,m,eval,work_query,lwork,rwork,lapack_info)
    if(lapack_info/=0)return;lwork=max(1,int(real(work_query(1),8)));allocate(work(lwork),stat=astat)
    if(astat/=0)return
    call zheev('N','U',mactual,tri,m,eval,work,lwork,rwork,lapack_info);if(lapack_info/=0)return
    min_eval=eval(1);max_eval=eval(mactual);info=0
  end subroutine lanczos_sww

  integer(8) function compute_metric_fingerprint(op)result(fingerprint)
    type(s_dg_wpw_bounded_operator),intent(in)::op
    integer(8)::local_hash;integer::i,ierr,rank
    local_hash=88172645463325252_8
    call MPI_Comm_rank(op%w_schedule%comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_8;return;endif
    call hash_word(local_hash,int(rank+1,8))
    do i=1,len_trim(op%metric_convention);call hash_word(local_hash,int(iachar(op%metric_convention(i:i)),8));enddo
    do i=1,size(op%owned_w_ids);call hash_word(local_hash,int(op%owned_w_ids(i),8));enddo
    do i=1,size(op%owned_p_ids);call hash_word(local_hash,int(op%owned_p_ids(i),8));enddo
    do i=1,size(op%required_w_ids);call hash_word(local_hash,int(op%required_w_ids(i),8));enddo
    do i=1,size(op%required_p_ids);call hash_word(local_hash,int(op%required_p_ids(i),8));enddo
    call hash_complex2(local_hash,op%ww_s_dense);call hash_complex2(local_hash,op%wp_s_dense)
    call hash_complex2(local_hash,op%pp_s_dense)
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_SUM,op%w_schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS)fingerprint=0_8
  end function compute_metric_fingerprint

  subroutine hash_complex2(hash,x)
    integer(8),intent(inout)::hash;complex(8),intent(in)::x(:,:);integer::i,j
    do j=1,size(x,2);do i=1,size(x,1)
      call hash_word(hash,transfer(real(x(i,j),8),hash));call hash_word(hash,transfer(aimag(x(i,j)),hash))
    enddo;enddo
  end subroutine hash_complex2
  subroutine hash_word(hash,word)
    integer(8),intent(inout)::hash;integer(8),intent(in)::word
    hash=ieor(ishftc(hash,13),word);hash=ieor(hash,ishft(hash,-7))
  end subroutine hash_word
  integer function find_id(ids,target)result(pos)
    integer,intent(in)::ids(:),target;integer::i
    pos=0;do i=1,size(ids);if(ids(i)==target)then;pos=i;return;endif;enddo
  end function find_id
  subroutine sort_int(x)
    integer,intent(inout)::x(:);integer::i,j,t
    do i=2,size(x);t=x(i);j=i-1;do while(j>=1);if(x(j)<=t)exit;x(j+1)=x(j);j=j-1;enddo;x(j+1)=t;enddo
  end subroutine sort_int
  logical function strictly_increasing(x)result(ok)
    integer,intent(in)::x(:);integer::i
    ok=size(x)>0.and.all(x>0);do i=2,size(x);if(x(i)<=x(i-1))ok=.false.;enddo
  end function strictly_increasing
  logical function finite2(x)result(ok)
    complex(8),intent(in)::x(:,:);ok=all(ieee_is_finite(real(x,8))).and.all(ieee_is_finite(aimag(x)))
  end function finite2
end module dg_wpw_s_orthogonal_complement
