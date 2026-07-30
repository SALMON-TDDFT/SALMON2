#include "config.h"
module dg_overlapping_wannier_solver
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::solve_dg_overlapping_wannier_coefficients
contains
  subroutine solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,nstate,max_iterations,&
      tolerance,metric_tolerance,coefficients,eigenvalues,maximum_residual,orthogonality_defect,&
      metric_condition,ok,message,initial_coefficients)
    integer,intent(in)::comm,nstate,max_iterations
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    real(8),intent(in)::tolerance,metric_tolerance
    complex(8),intent(out)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:),maximum_residual,orthogonality_defect,metric_condition
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),intent(in),optional::initial_coefficients(:,:)
#ifdef USE_MPI
    complex(8),allocatable::qlocal(:,:),vector(:),hvector(:),svector(:),residual(:),direction(:),&
      hdirection(:),sdirection(:),gram(:,:),block_basis(:,:),block_coefficients(:,:),&
      block_h(:,:),block_s(:,:),block_residual(:,:),candidate_block(:,:),previous_direction(:,:),old_q(:,:)
    complex(8),allocatable::reference_initial(:,:)
    complex(8)::value,reduced_h(2,2),reduced_s(2,2),ritz_vector(2)
    real(8),allocatable::block_eigenvalues(:)
    real(8)::ritz_value,pivot,residual_norm,block_maximum_residual,warm_defect,global_warm_defect
    integer,allocatable::ownership(:)
    integer::n,nlocal,i,j,band,iteration,ierr,local_bad,global_bad,block_size,maximum_block_size,&
      warm_local,warm_min,warm_max
    integer::scalar_local(3),scalar_min(3),scalar_max(3)
    real(8)::real_local(2),real_min(2),real_max(2)

    ok=.false.;message='';coefficients=(0d0,0d0);eigenvalues=0d0
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);metric_condition=huge(1d0)
    n=size(hrows,2);nlocal=size(row_ids)
    scalar_local=[n,nstate,max_iterations]
    call MPI_Allreduce(scalar_local,scalar_min,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(scalar_local,scalar_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    real_local=[tolerance,metric_tolerance]
    call MPI_Allreduce(real_local,real_min,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(real_local,real_max,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    warm_local=merge(1,0,present(initial_coefficients))
    call MPI_Allreduce(warm_local,warm_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(warm_local,warm_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=0
    if(any(scalar_min/=scalar_max).or.any(real_min/=real_max).or.warm_min/=warm_max)local_bad=1
    if(n<1.or.nstate<1.or.nstate>n.or.max_iterations<1)local_bad=1
    if(tolerance<=0d0.or.metric_tolerance<=0d0.or.metric_tolerance>=1d0)local_bad=1
    if(.not.all(ieee_is_finite(real_local)))local_bad=1
    if(size(hrows,1)/=nlocal.or.any(shape(srows)/=shape(hrows)))local_bad=1
    if(any(shape(coefficients)/=[n,nstate]).or.size(eigenvalues)/=nstate)local_bad=1
    if(any(row_ids<1_8).or.any(row_ids>int(n,8)))local_bad=1
    if(.not.finite_matrix(hrows).or..not.finite_matrix(srows))local_bad=1
    if(present(initial_coefficients))then
      if(any(shape(initial_coefficients)/=[n,nstate]).or..not.finite_matrix(initial_coefficients))local_bad=1
    endif
    if(int(n,int64)>0_int64.and.int(n,int64)>int(huge(1),int64)/int(n,int64))local_bad=1
    if(int(nstate,int64)>0_int64.and.int(nstate,int64)>int(huge(1),int64)/int(nstate,int64))local_bad=1
    if(int(nstate,int64)>int(huge(1),int64)/3_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed overlapping-Wannier coefficient solve contract';return;endif
    if(warm_max==1)then
      allocate(reference_initial,source=initial_coefficients)
      call MPI_Bcast(reference_initial,n*nstate,MPI_DOUBLE_COMPLEX,0,comm,ierr)
      warm_defect=maxval(abs(initial_coefficients-reference_initial))
      call MPI_Allreduce(warm_defect,global_warm_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(global_warm_defect>metric_tolerance)then
        message='rank-inconsistent warm-start coefficient block';return
      endif
    endif

    allocate(ownership(n));ownership=0
    do i=1,nlocal;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(ownership/=1))then;message='coefficient rows do not form a unique global partition';return;endif
    call distributed_hermiticity(comm,row_ids,hrows,tolerance,local_bad)
    call distributed_hermiticity(comm,row_ids,srows,metric_tolerance,global_bad)
    if(local_bad/=0.or.global_bad/=0)then;message='non-Hermitian distributed H or S';return;endif

    allocate(qlocal(nlocal,n),vector(n),hvector(nlocal),svector(nlocal),residual(nlocal),&
      direction(nlocal),hdirection(nlocal),sdirection(nlocal),gram(nstate,nstate))
    qlocal=(0d0,0d0)
    do j=1,n
      do i=1,nlocal
        if(row_ids(i)==int(j,8))qlocal(i,j)=1d0
      enddo
      call s_orthogonalize(comm,row_ids,srows,qlocal(:,j),qlocal(:,1:j-1),pivot,ok)
      if(.not.ok)then;message='overlap metric is indefinite or numerically singular';return;endif
    enddo
    call estimate_metric_condition(comm,row_ids,srows,qlocal,metric_condition,ok)
    if(.not.ok.or.metric_condition*metric_tolerance>=1d0)then
      ok=.false.;message='overlap metric condition gate failed';return
    endif
    if(present(initial_coefficients))then
      do j=1,nstate
        do i=1,nlocal
          qlocal(i,j)=initial_coefficients(int(row_ids(i)),j)
        enddo
        call s_orthogonalize(comm,row_ids,srows,qlocal(:,j),qlocal(:,1:j-1),pivot,ok)
        if(.not.ok)then;message='warm-start coefficient block lost metric rank';return;endif
      enddo
    endif

    maximum_block_size=min(n,3*nstate)
    allocate(block_basis(nlocal,maximum_block_size),block_coefficients(n,maximum_block_size),&
      block_eigenvalues(maximum_block_size),block_h(nlocal,nstate),block_s(nlocal,nstate),&
      block_residual(nlocal,nstate),candidate_block(nlocal,2*nstate),&
      previous_direction(nlocal,nstate),old_q(nlocal,nstate))
    previous_direction=(0d0,0d0)
    do iteration=1,max_iterations
      do band=1,nstate
        call gather_vector(comm,row_ids,qlocal(:,band),coefficients(:,band))
      enddo
      call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal(:,1:nstate),coefficients,eigenvalues,ok)
      if(.not.ok)then;message='coefficient block Rayleigh-Ritz failed';return;endif
      block_h=matmul(hrows,coefficients);block_s=matmul(srows,coefficients)
      block_maximum_residual=0d0
      do band=1,nstate
        block_residual(:,band)=block_h(:,band)-eigenvalues(band)*block_s(:,band)
        call relative_residual(comm,block_residual(:,band),block_h(:,band),block_s(:,band),&
          eigenvalues(band),residual_norm)
        if(.not.ieee_is_finite(residual_norm))then
          ok=.false.;message='nonfinite coefficient block residual';return
        endif
        block_maximum_residual=max(block_maximum_residual,residual_norm)
      enddo
      if(block_maximum_residual<=tolerance)exit
      block_basis=(0d0,0d0);block_basis(:,1:nstate)=qlocal(:,1:nstate);block_size=nstate
      old_q=qlocal(:,1:nstate)
      candidate_block(:,1:nstate)=-block_residual
      candidate_block(:,nstate+1:2*nstate)=previous_direction
      call append_s_orthonormal_block(comm,row_ids,srows,&
        candidate_block(:,1:merge(2*nstate,nstate,iteration>1)),block_basis,block_size,ok)
      if(.not.ok)then;message='coefficient block residual orthogonalization failed';return;endif
      if(block_size==nstate)then
        ok=.false.;message='coefficient block residual space lost metric rank';return
      endif
      do j=1,block_size
        call gather_vector(comm,row_ids,block_basis(:,j),block_coefficients(:,j))
      enddo
      call block_rayleigh_ritz(comm,row_ids,hrows,srows,block_basis(:,1:block_size),&
        block_coefficients(:,1:block_size),block_eigenvalues(1:block_size),ok)
      if(.not.ok)then;message='coefficient expanded block Rayleigh-Ritz failed';return;endif
      previous_direction=old_q
      qlocal(:,1:nstate)=block_basis(:,1:nstate)
      if(block_size==n)then
        eigenvalues=block_eigenvalues(1:nstate)
        coefficients=block_coefficients(:,1:nstate)
        call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
          maximum_residual,orthogonality_defect)
        if(maximum_residual<=10d0*tolerance.and.orthogonality_defect<=10d0*tolerance)exit
        ok=.false.;message='full-space Ritz residual exceeds numerical quality gate';return
      endif
    enddo
    if(iteration>max_iterations)then
      ok=.false.
      write(message,'(a,i0,a,es12.4,a,es12.4)')&
        'coefficient block iteration did not converge: iterations=',max_iterations,&
        ' residual=',block_maximum_residual,' metric_condition=',metric_condition
      return
    endif
    do band=1,nstate
      call gather_vector(comm,row_ids,qlocal(:,band),coefficients(:,band))
    enddo
    call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal(:,1:nstate),coefficients,eigenvalues,ok)
    if(.not.ok)then;message='final block Rayleigh-Ritz failed';return;endif
    call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
      maximum_residual,orthogonality_defect)
    if(maximum_residual>10d0*tolerance.or.orthogonality_defect>10d0*tolerance)then
      ok=.false.;message='coefficient residual or S-orthonormality gate failed';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier coefficient solve requires MPI'
#endif
  end subroutine

#ifdef USE_MPI
  subroutine gather_vector(comm,row_ids,local,global)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::local(:)
    complex(8),intent(out)::global(:)
    complex(8),allocatable::staged(:)
    integer::i,ierr
    allocate(staged(size(global)));staged=(0d0,0d0)
    do i=1,size(local);staged(int(row_ids(i)))=local(i);enddo
    call MPI_Allreduce(staged,global,size(global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine inner(comm,left,right,value)
    integer,intent(in)::comm
    complex(8),intent(in)::left(:),right(:)
    complex(8),intent(out)::value
    complex(8)::local
    integer::ierr
    local=sum(conjg(left)*right)
    call MPI_Allreduce(local,value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine s_orthogonalize(comm,row_ids,srows,vector,previous,pivot,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    real(8),intent(out)::pivot
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
    call inner(comm,vector,svector,value);pivot=real(value,8)
    ok=ieee_is_finite(pivot).and.pivot>0d0.and.abs(aimag(value))<=1d-10*max(1d0,abs(value))
    if(ok)vector=vector/sqrt(pivot)
  end subroutine

  subroutine s_project_out(comm,row_ids,srows,vector,previous,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    ok=finite_matrix(reshape(vector,[size(vector),1]))
  end subroutine

  subroutine distributed_hermiticity(comm,row_ids,rows,tolerance,bad)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::rows(:,:)
    real(8),intent(in)::tolerance
    integer,intent(out)::bad
    integer,parameter::row_batch_size=32
    complex(8),allocatable::local_batch(:,:),global_batch(:,:)
    real(8)::local_scale,scale,local_defect,defect
    integer::i,j,ierr,n,first,last,nbatch
    n=size(rows,2)
    allocate(local_batch(min(row_batch_size,n),n),global_batch(min(row_batch_size,n),n))
    local_scale=1d0
    if(size(rows)>0)local_scale=max(local_scale,maxval(abs(rows)))
    local_defect=0d0
    do first=1,n,row_batch_size
      last=min(n,first+row_batch_size-1);nbatch=last-first+1
      local_batch=(0d0,0d0)
      do i=1,size(row_ids)
        if(row_ids(i)>=int(first,8).and.row_ids(i)<=int(last,8))&
          local_batch(int(row_ids(i))-first+1,:)=rows(i,:)
      enddo
      call MPI_Allreduce(local_batch,global_batch,size(local_batch),&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      do i=1,size(row_ids)
        do j=first,last
          local_defect=max(local_defect,&
            abs(rows(i,j)-conjg(global_batch(j-first+1,int(row_ids(i))))))
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_scale,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_defect,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    bad=merge(1,0,defect>tolerance*scale)
  end subroutine

  subroutine relative_residual(comm,residual,hvector,svector,eigenvalue,value)
    integer,intent(in)::comm
    complex(8),intent(in)::residual(:),hvector(:),svector(:)
    real(8),intent(in)::eigenvalue
    real(8),intent(out)::value
    complex(8)::z
    real(8)::nr,nh,ns
    call inner(comm,residual,residual,z);nr=sqrt(max(0d0,real(z,8)))
    call inner(comm,hvector,hvector,z);nh=sqrt(max(0d0,real(z,8)))
    call inner(comm,svector,svector,z);ns=sqrt(max(0d0,real(z,8)))
    value=nr/max(tiny(1d0),nh+abs(eigenvalue)*ns)
  end subroutine

  subroutine estimate_metric_condition(comm,row_ids,srows,q,condition,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),q(:,:)
    real(8),intent(out)::condition
    logical,intent(out)::ok
    complex(8),allocatable::gathered_column(:),root_q(:,:)
    integer(8),allocatable::gathered_ids(:)
    complex(8)::gram_value
    real(8)::squares_s,squares_inverse,local_squares,local_scale,scale_s,&
      norm_s,norm_inverse
    real(8)::factor_s
    integer,allocatable::counts(:),displacements(:)
    integer::ierr,n,i,j,k,rank,nproc,nlocal
    n=size(q,2)
    nlocal=size(row_ids)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    local_scale=0d0
    if(size(srows)>0)local_scale=maxval(abs(srows))
    call MPI_Allreduce(local_scale,scale_s,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(scale_s<=0d0.or..not.ieee_is_finite(scale_s))then;ok=.false.;return;endif
    local_squares=sum((abs(srows)/scale_s)**2)
    call MPI_Allreduce(local_squares,squares_s,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    allocate(counts(nproc),displacements(nproc))
    call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm,ierr)
    if(rank==0)then
      displacements(1)=0
      do i=2,nproc;displacements(i)=displacements(i-1)+counts(i-1);enddo
      allocate(gathered_ids(n),gathered_column(n),root_q(n,n))
    else
      allocate(gathered_ids(1),gathered_column(1),root_q(1,1))
    endif
    call MPI_Gatherv(row_ids,nlocal,MPI_INTEGER8,gathered_ids,counts,displacements,&
      MPI_INTEGER8,0,comm,ierr)
    do j=1,n
      call MPI_Gatherv(q(:,j),nlocal,MPI_DOUBLE_COMPLEX,gathered_column,counts,displacements,&
        MPI_DOUBLE_COMPLEX,0,comm,ierr)
      if(rank==0)then
        do k=1,n;root_q(int(gathered_ids(k)),j)=gathered_column(k);enddo
      endif
    enddo
    squares_inverse=0d0
    if(rank==0)then
      do j=1,n;do i=1,n
        gram_value=dot_product(root_q(:,i),root_q(:,j))
        squares_inverse=squares_inverse+abs(gram_value)**2
      enddo;enddo
    endif
    call MPI_Bcast(squares_inverse,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    factor_s=sqrt(squares_s)
    if(factor_s<=0d0)then;ok=.false.;return;endif
    if(squares_inverse<=0d0.or..not.ieee_is_finite(squares_inverse))then;ok=.false.;return;endif
    if(factor_s>1d0)then
      if(scale_s>huge(1d0)/factor_s)then
        condition=huge(1d0);ok=.false.;return
      endif
    endif
    norm_s=scale_s*factor_s;norm_inverse=sqrt(squares_inverse)
    if(norm_inverse>1d0)then
      if(norm_s>huge(1d0)/norm_inverse)then
        condition=huge(1d0);ok=.false.;return
      else
        condition=norm_s*norm_inverse
      endif
    else
      condition=norm_s*norm_inverse
    endif
    ok=ieee_is_finite(condition).and.condition>=1d0
  end subroutine

  subroutine block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal,coefficients,eigenvalues,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    complex(8),intent(inout)::qlocal(:,:)
    complex(8),intent(inout)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:)
    logical,intent(out)::ok
    complex(8),allocatable::hq(:,:),projected(:,:),local_projected(:,:),rotated_q(:,:),work(:)
    real(8),allocatable::rwork(:)
    integer::p,ierr,n,lapack_info,lwork
    external::zheev
    n=size(eigenvalues)
    hq=matmul(hrows,coefficients)
    local_projected=matmul(conjg(transpose(qlocal)),hq)
    allocate(projected(size(eigenvalues),size(eigenvalues)))
    call MPI_Allreduce(local_projected,projected,size(eigenvalues)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    projected=0.5d0*(projected+conjg(transpose(projected)))
    if(.not.finite_matrix(projected))then;ok=.false.;return;endif
    allocate(rwork(max(1,3*n-2)))
    lwork=max(1,2*n-1)
    allocate(work(lwork))
    call zheev('V','U',n,projected,n,eigenvalues,work,lwork,rwork,lapack_info)
    if(lapack_info/=0.or..not.finite_matrix(projected))then;ok=.false.;return;endif
    allocate(rotated_q(size(qlocal,1),n))
    rotated_q=matmul(qlocal,projected)
    qlocal=rotated_q
    do p=1,size(eigenvalues);call gather_vector(comm,row_ids,qlocal(:,p),coefficients(:,p));enddo
    ok=all(ieee_is_finite(eigenvalues)).and.finite_matrix(coefficients)
  end subroutine

  subroutine append_s_orthonormal_block(comm,row_ids,srows,candidates,basis,basis_size,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),candidates(:,:)
    complex(8),intent(inout)::basis(:,:)
    integer,intent(inout)::basis_size
    logical,intent(out)::ok
    complex(8),allocatable::work_block(:,:),global_block(:,:),sblock(:,:),cross(:,:),local_cross(:,:),&
      gram(:,:),local_gram(:,:),work(:),rotated(:,:)
    real(8),allocatable::spectrum(:),rwork(:)
    integer::ncandidate,n,keep,i,j,ierr,info,lwork,projection_pass
    external::zheev
    ncandidate=size(candidates,2);n=size(srows,2)
    allocate(work_block(size(candidates,1),ncandidate),global_block(n,ncandidate),&
      sblock(size(candidates,1),ncandidate))
    work_block=candidates
    if(basis_size>0)then
      allocate(local_cross(basis_size,ncandidate),cross(basis_size,ncandidate))
      do projection_pass=1,2
        call gather_block(comm,row_ids,work_block,global_block)
        sblock=matmul(srows,global_block)
        local_cross=matmul(conjg(transpose(basis(:,1:basis_size))),sblock)
        call MPI_Allreduce(local_cross,cross,basis_size*ncandidate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        work_block=work_block-matmul(basis(:,1:basis_size),cross)
      enddo
    endif
    call gather_block(comm,row_ids,work_block,global_block)
    sblock=matmul(srows,global_block)
    allocate(local_gram(ncandidate,ncandidate),gram(ncandidate,ncandidate),spectrum(ncandidate))
    local_gram=matmul(conjg(transpose(work_block)),sblock)
    call MPI_Allreduce(local_gram,gram,ncandidate*ncandidate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    gram=0.5d0*(gram+conjg(transpose(gram)))
    allocate(rwork(max(1,3*ncandidate-2)))
    lwork=max(1,2*ncandidate-1);allocate(work(lwork))
    call zheev('V','U',ncandidate,gram,ncandidate,spectrum,work,lwork,rwork,info)
    if(info/=0.or..not.all(ieee_is_finite(spectrum)))then;ok=.false.;return;endif
    keep=count(spectrum>max(tiny(1d0),1d-12*max(0d0,maxval(spectrum))))
    keep=min(keep,size(basis,2)-basis_size)
    if(keep<1)then;ok=.true.;return;endif
    allocate(rotated(size(work_block,1),keep))
    do j=1,keep
      i=ncandidate-keep+j
      rotated(:,j)=matmul(work_block,gram(:,i))/sqrt(spectrum(i))
    enddo
    basis(:,basis_size+1:basis_size+keep)=rotated
    basis_size=basis_size+keep
    ok=finite_matrix(rotated)
  end subroutine

  subroutine gather_block(comm,row_ids,local,global)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::local(:,:)
    complex(8),intent(out)::global(:,:)
    complex(8),allocatable::staged(:,:)
    integer::i,ierr
    allocate(staged(size(global,1),size(global,2)));staged=(0d0,0d0)
    do i=1,size(local,1);staged(int(row_ids(i)),:)=local(i,:);enddo
    call MPI_Allreduce(staged,global,size(global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine coefficient_diagnostics(comm,row_ids,hrows,srows,c,e,residual,orthogonality)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:),c(:,:)
    real(8),intent(in)::e(:)
    real(8),intent(out)::residual,orthogonality
    complex(8),allocatable::rlocal(:,:),gram(:,:),local_gram(:,:),hv(:),sv(:)
    real(8)::band_residual
    integer::j,ierr
    allocate(rlocal(size(row_ids),size(e)),gram(size(e),size(e)),local_gram(size(e),size(e)))
    allocate(hv(size(row_ids)),sv(size(row_ids)));residual=0d0
    do j=1,size(e)
      hv=matmul(hrows,c(:,j));sv=matmul(srows,c(:,j));rlocal(:,j)=hv-e(j)*sv
      call relative_residual(comm,rlocal(:,j),hv,sv,e(j),band_residual)
      residual=max(residual,band_residual)
    enddo
    local_gram=matmul(conjg(transpose(c(int(row_ids),:))),matmul(srows,c))
    call MPI_Allreduce(local_gram,gram,size(e)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do j=1,size(e);gram(j,j)=gram(j,j)-1d0;enddo
    orthogonality=maxval(abs(gram))
    if(.not.all(ieee_is_finite([residual,orthogonality])))then
      residual=huge(1d0);orthogonality=huge(1d0)
    endif
  end subroutine

  subroutine lowest_generalized_2x2(h,s,eigenvalue,eigenvector,ok)
    complex(8),intent(in)::h(2,2),s(2,2)
    real(8),intent(out)::eigenvalue
    complex(8),intent(out)::eigenvector(2)
    logical,intent(out)::ok
    complex(8)::l(2,2),linv(2,2),a(2,2),u(2),phase,hs(2,2),temporary(2,2)
    real(8)::l11,l22,s22_root,disc_scaled,norm,scale,delta,t,c,absb_scaled,&
      difference_scaled,eigenvalue_scaled,xdisc,ydisc,metric_ratio,phase_threshold
    real(8)::hamiltonian_scale,total_scale
    complex(8)::quotient
    integer::i,j,k
    ok=.false.;eigenvalue=0d0;eigenvector=(0d0,0d0)
    if(.not.finite_matrix(h).or..not.finite_matrix(s))return
    if(real(s(1,1),8)<=0d0.or.real(s(2,2),8)<=0d0)return
    l11=sqrt(real(s(1,1),8))
    l=(0d0,0d0);l(1,1)=l11
    call safe_complex_divide(s(2,1),l11,l(2,1),ok);if(.not.ok)return
    s22_root=sqrt(real(s(2,2),8))
    if(safe_complex_abs(l(2,1))>=s22_root)return
    metric_ratio=safe_complex_abs(l(2,1))/s22_root
    l22=s22_root*sqrt((1d0-metric_ratio)*(1d0+metric_ratio));l(2,2)=l22
    if(.not.all(ieee_is_finite([l11,l22])).or.min(l11,l22)<=0d0)return
    if(l11<1d0/huge(1d0).or.l22<1d0/huge(1d0))return
    linv=(0d0,0d0);linv(1,1)=1d0/l11;linv(2,2)=1d0/l22
    call safe_complex_divide(l(2,1),l11,quotient,ok);if(.not.ok)return
    call safe_complex_divide(-quotient,l22,linv(2,1),ok);if(.not.ok)return
    hamiltonian_scale=max(tiny(1d0),max(maxval(abs(real(h,8))),maxval(abs(aimag(h)))))
    hs=h/hamiltonian_scale
    temporary=(0d0,0d0);a=(0d0,0d0)
    do j=1,2;do i=1,2;do k=1,2
      temporary(i,j)=temporary(i,j)+linv(i,k)*hs(k,j)
    enddo;enddo;enddo
    do j=1,2;do i=1,2;do k=1,2
      a(i,j)=a(i,j)+temporary(i,k)*conjg(linv(j,k))
    enddo;enddo;enddo
    a=0.5d0*(a+conjg(transpose(a)))
    if(.not.finite_matrix(a))return
    scale=max(tiny(1d0),max(maxval(abs(real(a,8))),maxval(abs(aimag(a)))))
    delta=real(a(1,1),8)/scale-real(a(2,2),8)/scale
    xdisc=abs(delta)
    absb_scaled=safe_real_hypot(real(a(1,2),8)/scale,aimag(a(1,2))/scale)
    ydisc=2d0*absb_scaled
    disc_scaled=sqrt(xdisc*xdisc+ydisc*ydisc)
    eigenvalue_scaled=0.5d0*(real(a(1,1),8)/scale+real(a(2,2),8)/scale-disc_scaled)
    if(scale>1d0)then
      if(hamiltonian_scale>huge(1d0)/scale)return
    endif
    total_scale=scale*hamiltonian_scale
    if(abs(eigenvalue_scaled)>1d0)then
      if(total_scale>huge(1d0)/abs(eigenvalue_scaled))return
    endif
    eigenvalue=total_scale*eigenvalue_scaled
    phase_threshold=epsilon(1d0)*max(1d0,1d0/scale)
    if(absb_scaled>phase_threshold)then
      phase=(a(1,2)/scale)/absb_scaled
      if(real(a(1,1),8)<=real(a(2,2),8))then
        difference_scaled=real(a(2,2),8)/scale-real(a(1,1),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[cmplx(c,0d0,8),-t*c*conjg(phase)]
      else
        difference_scaled=real(a(1,1),8)/scale-real(a(2,2),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[-t*c*phase,cmplx(c,0d0,8)]
      endif
    else
      u=(0d0,0d0);u(merge(1,2,real(a(1,1),8)<=real(a(2,2),8)))=1d0
    endif
    scale=max(safe_complex_abs(u(1)),safe_complex_abs(u(2)))
    if(scale<=0d0.or..not.ieee_is_finite(scale))return
    norm=scale*safe_real_hypot(safe_complex_abs(u(1))/scale,safe_complex_abs(u(2))/scale)
    if(norm<=0d0.or..not.ieee_is_finite(norm))return
    u=u/norm;eigenvector=matmul(conjg(transpose(linv)),u)
    ok=all(ieee_is_finite(real(eigenvector,8))).and.all(ieee_is_finite(aimag(eigenvector)))
  end subroutine
#endif

  real(8) function safe_complex_abs(value)
    complex(8),intent(in)::value
    safe_complex_abs=safe_real_hypot(real(value,8),aimag(value))
  end function

  real(8) function safe_real_hypot(x,y)
    real(8),intent(in)::x,y
    real(8)::larger,ratio,factor
    larger=max(abs(x),abs(y))
    if(larger==0d0)then
      safe_real_hypot=0d0
      return
    endif
    ratio=min(abs(x),abs(y))/larger
    factor=sqrt(1d0+ratio*ratio)
    if(larger>huge(1d0)/factor)then
      safe_real_hypot=huge(1d0)
    else
      safe_real_hypot=larger*factor
    endif
  end function

  subroutine safe_complex_divide(value,divisor,result,ok)
    complex(8),intent(in)::value
    real(8),intent(in)::divisor
    complex(8),intent(out)::result
    logical,intent(out)::ok
    real(8)::magnitude
    result=(0d0,0d0);ok=.false.
    if(divisor<=0d0.or..not.ieee_is_finite(divisor))return
    magnitude=safe_complex_abs(value)
    if(divisor<1d0)then
      if(magnitude>huge(1d0)*divisor)return
    endif
    result=value/divisor
    ok=ieee_is_finite(real(result,8))
    if(ok)ok=ieee_is_finite(aimag(result))
  end subroutine

  logical function finite_matrix(matrix)
    complex(8),intent(in)::matrix(:,:)
    integer::i,j
    real(8)::real_part,imaginary_part
    finite_matrix=.true.
    do j=1,size(matrix,2)
      do i=1,size(matrix,1)
        real_part=real(matrix(i,j),8)
        imaginary_part=aimag(matrix(i,j))
        if(real_part/=real_part)then
          finite_matrix=.false.
          return
        endif
        if(imaginary_part/=imaginary_part)then
          finite_matrix=.false.
          return
        endif
        if(abs(real_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
        if(abs(imaginary_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
      enddo
    enddo
  end function
end module
