module dg_spatial_mlwf_seed
  use iso_fortran_env, only: real64,int64
  use mpi
  implicit none
  private
  public::spatial_seed_options,spatial_seed_diagnostics,spatial_seed_sparse_matrix,&
    make_spatial_seed_sparse_matrix,construct_spatial_mlwf_seed,apply_s_projector
  public::certified_sparse_metric_solve

  type spatial_seed_sparse_matrix
    integer::extent=0
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::values(:)
  end type

  type spatial_seed_options
    real(real64)::lower_bound=0d0,upper_bound=0d0
    real(real64)::occupied_upper=0d0,empty_lower=0d0,empty_upper=0d0
    integer::filter_degree=8,max_metric_iterations=80,fill_limit=huge(1)
    real(real64)::tolerance=1d-10,rank_tolerance=1d-12,condition_limit=1d8
  end type
  type spatial_seed_diagnostics
    real(real64)::metric_defect=huge(1d0),metric_condition=huge(1d0)
    real(real64)::occupied_projector_defect=huge(1d0),empty_projector_defect=huge(1d0)
    real(real64)::mutual_projector_defect=huge(1d0),ritz_residual=huge(1d0)
    real(real64)::intertwining_defect=huge(1d0),map_metric_defect=huge(1d0)
    real(real64)::map_linearity_defect=huge(1d0),map_equivariance_defect=huge(1d0)
    real(real64)::back_transform_defect=huge(1d0)
    real(real64)::metric_solve_residual=huge(1d0)
    integer::metric_iterations=0,filter_fill=0
    integer(int64)::dense_seed_workspace_bytes=0_int64
    integer::sparse_nnz_peak=0
    integer(int64)::sparse_workspace_peak_bytes=0_int64
  end type
  interface construct_spatial_mlwf_seed
    module procedure construct_spatial_mlwf_seed_dense
    module procedure construct_spatial_mlwf_seed_sparse
  end interface
contains
  subroutine make_spatial_seed_sparse_matrix(matrix,drop_tolerance,sparse,ok,message)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),intent(in)::drop_tolerance
    type(spatial_seed_sparse_matrix),intent(out)::sparse
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,k,n,nnz
    ok=.false.;message='';n=size(matrix,1)
    if(size(matrix,2)/=n.or.drop_tolerance<0d0)then;message='invalid sparse matrix fixture';return;end if
    nnz=count(abs(matrix)>drop_tolerance);sparse%extent=n
    allocate(sparse%row_offsets(n+1),sparse%column_ids(nnz),sparse%values(nnz));k=1;sparse%row_offsets(1)=1
    do i=1,n
      do j=1,n
        if(abs(matrix(i,j))>drop_tolerance)then
          sparse%column_ids(k)=j;sparse%values(k)=matrix(i,j);k=k+1
        end if
      end do
      sparse%row_offsets(i+1)=k
    end do
    ok=.true.;message='ok'
  end subroutine

  subroutine construct_spatial_mlwf_seed_sparse(comm,h_sparse,s_sparse,yocc,yempty,bg,docc,dempty,a,tg,opt,&
      cocc,cempty,diag,fingerprint,ok,message)
    integer,intent(in)::comm
    type(spatial_seed_sparse_matrix),intent(in)::h_sparse,s_sparse
    complex(real64),intent(in)::yocc(:,:),yempty(:,:),bg(:,:,:),docc(:,:,:),dempty(:,:,:),a(:,:),tg(:,:,:)
    type(spatial_seed_options),intent(in)::opt
    complex(real64),allocatable,intent(out)::cocc(:,:),cempty(:,:)
    type(spatial_seed_diagnostics),intent(out)::diag
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(spatial_seed_sparse_matrix)::x,xx,sxx,update,xnew,xadj,metric_product
    complex(real64),allocatable::q(:,:),trial(:,:),sc(:,:),hc(:,:),gram(:,:),cross(:,:),lambda(:,:),unitary(:,:)
    real(real64)::s_lower,s_upper,h_lower,h_upper,radius,scale,defect
    integer::i,k,n,iteration,ierr,local_bad,global_bad
    logical::local_ok
    character(256)::local_message
    call validate_sparse_matrix(h_sparse,ok,message);if(.not.ok)return
    call validate_sparse_matrix(s_sparse,ok,message);if(.not.ok)return
    if(h_sparse%extent/=s_sparse%extent)then;ok=.false.;message='sparse H/S extent mismatch';return;end if
    ok=.false.;message='';diag=spatial_seed_diagnostics();n=s_sparse%extent;diag%dense_seed_workspace_bytes=0_int64
    if(size(yocc,1)/=n.or.size(yempty,1)/=n)then;message='sparse seed trial shape mismatch';return;end if
    call sparse_hermitian_bounds(s_sparse,s_lower,s_upper,local_ok)
    if(.not.local_ok.or.s_lower<=0d0)then;message='S is not certified Hermitian positive definite';return;end if
    diag%metric_condition=s_upper/s_lower
    if(diag%metric_condition>opt%condition_limit)then;message='S condition limit exceeded';return;end if
    scale=s_upper;call sparse_identity(n,1d0/sqrt(scale),x)
    do iteration=1,opt%max_metric_iterations
      call sparse_multiply(x,x,opt%rank_tolerance*0.01d0,opt%fill_limit,xx,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;end if
      call sparse_multiply(s_sparse,xx,opt%rank_tolerance*0.01d0,opt%fill_limit,sxx,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;end if
      call sparse_three_identity_minus(sxx,update)
      call sparse_multiply(x,update,opt%rank_tolerance*0.01d0,opt%fill_limit,xnew,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;end if
      call sparse_scale(xnew,0.5d0);x=xnew
      call sparse_adjoint(x,xadj)
      call sparse_multiply(xadj,s_sparse,opt%rank_tolerance*0.01d0,opt%fill_limit,sxx,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;end if
      call sparse_multiply(sxx,x,opt%rank_tolerance*0.01d0,opt%fill_limit,metric_product,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;end if
      diag%metric_defect=sparse_identity_defect(metric_product);diag%metric_iterations=iteration
      diag%sparse_nnz_peak=max(diag%sparse_nnz_peak,size(x%values),size(xx%values),size(sxx%values),&
        size(update%values),size(metric_product%values))
      if(diag%metric_defect<=opt%tolerance)exit
    end do
    if(diag%metric_defect>opt%tolerance)then;message='metric purification did not converge';return;end if
    call sparse_hermitian_bounds(h_sparse,h_lower,h_upper,local_ok)
    if(.not.local_ok)then;message='H is not Hermitian';return;end if
    h_lower=h_lower/s_upper;h_upper=h_upper/s_lower
    if(opt%lower_bound>h_lower+opt%tolerance.or.opt%upper_bound<h_upper-opt%tolerance)then
      message='spectral bounds do not enclose verified sparse generalized bounds';return
    end if
    if(opt%empty_lower<=opt%occupied_upper)then;message='occupied-empty gap is closed';return;end if
    call sparse_filtered_block(h_sparse,x,xadj,yocc,opt%occupied_upper,opt%lower_bound,opt%filter_degree,q)
    call sparse_apply(x,q,trial);diag%back_transform_defect=0d0
    call reynolds_sparse(s_sparse,trial,bg,docc,opt%rank_tolerance,cocc,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    call sparse_filtered_block(h_sparse,x,xadj,yempty,opt%upper_bound,opt%empty_lower,opt%filter_degree,q)
    call sparse_apply(x,q,trial);call sparse_apply(s_sparse,trial,sc)
    cross=matmul(conjg(transpose(cocc)),sc);trial=trial-matmul(cocc,cross)
    call reynolds_sparse(s_sparse,trial,bg,dempty,opt%rank_tolerance,cempty,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    diag%filter_fill=size(x%values)+count(abs(cocc)>opt%rank_tolerance)+count(abs(cempty)>opt%rank_tolerance)
    diag%sparse_workspace_peak_bytes=16_int64*int(5*n+2*diag%sparse_nnz_peak,int64)+&
      4_int64*int(5*n+2*diag%sparse_nnz_peak,int64)
    if(diag%filter_fill>opt%fill_limit)then;message='filter fill limit exceeded';return;end if
    call sparse_gram(s_sparse,cocc,cocc,gram);call subtract_identity(gram)
    diag%occupied_projector_defect=maxval(abs(gram))
    call sparse_gram(s_sparse,cempty,cempty,gram);call subtract_identity(gram)
    diag%empty_projector_defect=maxval(abs(gram))
    call sparse_gram(s_sparse,cocc,cempty,cross);diag%mutual_projector_defect=maxval(abs(cross))
    diag%ritz_residual=max(sparse_ritz_defect(h_sparse,s_sparse,cocc),sparse_ritz_defect(h_sparse,s_sparse,cempty))
    diag%intertwining_defect=max(intertwining_defect_sparse(cocc,bg,docc),intertwining_defect_sparse(cempty,bg,dempty))
    diag%map_metric_defect=max(gram_defect_sparse(s_sparse,cocc,a),gram_defect_sparse(s_sparse,cempty,a))
    allocate(unitary(max(size(cocc,2),size(cempty,2)),max(size(cocc,2),size(cempty,2))))
    unitary=(0d0,0d0);do i=1,size(unitary,1);unitary(i,i)=1d0;end do
    diag%map_linearity_defect=max(linearity_defect(a,cocc,unitary(:size(cocc,2),:size(cocc,2))),&
      linearity_defect(a,cempty,unitary(:size(cempty,2),:size(cempty,2))))
    diag%map_equivariance_defect=max(equivariance_defect(a,cocc,bg,tg),equivariance_defect(a,cempty,bg,tg))
    defect=max(diag%occupied_projector_defect,diag%empty_projector_defect,diag%mutual_projector_defect,&
      diag%ritz_residual,diag%intertwining_defect,diag%map_metric_defect,diag%map_linearity_defect,&
      diag%map_equivariance_defect)
    if(diag%map_metric_defect>opt%tolerance)then;message='map metric pullback failed';return;end if
    if(diag%map_equivariance_defect>opt%tolerance)then;message='map equivariance failed';return;end if
    if(defect>10d0*opt%tolerance)then;message='sparse seed projector, Ritz, or symmetry gate failed';return;end if
    fingerprint=1469598103934665603_int64
    do i=1,size(cocc,1);fingerprint=ieor(ishftc(fingerprint,7),transfer(real(cocc(i,1)),fingerprint));end do
    local_bad=0;call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='MPI sparse seed agreement failed';return;end if
    ok=.true.;message='ok'
  end subroutine

  subroutine validate_sparse_matrix(matrix,ok,message)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    logical,intent(out)::ok;character(*),intent(out)::message
    integer::i,k
    ok=.false.;message='invalid sparse CSR contract'
    if(matrix%extent<1.or..not.allocated(matrix%row_offsets).or..not.allocated(matrix%column_ids).or.&
      .not.allocated(matrix%values))return
    if(size(matrix%row_offsets)/=matrix%extent+1.or.size(matrix%column_ids)/=size(matrix%values))return
    if(matrix%row_offsets(1)/=1.or.matrix%row_offsets(matrix%extent+1)/=size(matrix%values)+1)return
    if(any(matrix%row_offsets(2:)<matrix%row_offsets(:matrix%extent)))return
    do i=1,matrix%extent
      do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
        if(matrix%column_ids(k)<1.or.matrix%column_ids(k)>matrix%extent)return
        if(k>matrix%row_offsets(i))then;if(matrix%column_ids(k)<=matrix%column_ids(k-1))return;end if
      end do
    end do
    ok=.true.;message='ok'
  end subroutine

  subroutine sparse_to_dense(sparse,matrix)
    type(spatial_seed_sparse_matrix),intent(in)::sparse
    complex(real64),allocatable,intent(out)::matrix(:,:)
    integer::i,k
    allocate(matrix(sparse%extent,sparse%extent));matrix=(0d0,0d0)
    do i=1,sparse%extent;do k=sparse%row_offsets(i),sparse%row_offsets(i+1)-1
      matrix(i,sparse%column_ids(k))=sparse%values(k)
    end do;end do
  end subroutine

  subroutine sparse_identity(n,scale,matrix)
    integer,intent(in)::n;real(real64),intent(in)::scale
    type(spatial_seed_sparse_matrix),intent(out)::matrix
    integer::i
    matrix%extent=n;allocate(matrix%row_offsets(n+1),matrix%column_ids(n),matrix%values(n))
    do i=1,n;matrix%row_offsets(i)=i;matrix%column_ids(i)=i;matrix%values(i)=scale;end do
    matrix%row_offsets(n+1)=n+1
  end subroutine

  subroutine sparse_multiply(left,right,drop_tolerance,fill_limit,result,ok,message)
    type(spatial_seed_sparse_matrix),intent(in)::left,right
    real(real64),intent(in)::drop_tolerance
    integer,intent(in)::fill_limit
    type(spatial_seed_sparse_matrix),intent(out)::result
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::accumulator(:),values(:)
    logical,allocatable::marked(:)
    integer,allocatable::touched(:),columns(:),offsets(:)
    integer::n,i,ka,kb,j,ntouched,nnz,capacity,p,q,tmp
    n=left%extent;ok=.false.;message=''
    if(right%extent/=n)then;message='sparse multiply extent mismatch';return;end if
    capacity=max(1,min(fill_limit,max(n,size(left%values)+size(right%values))))
    allocate(accumulator(n),marked(n),touched(n),columns(capacity),values(capacity),offsets(n+1))
    accumulator=(0d0,0d0);marked=.false.;nnz=0;offsets(1)=1
    do i=1,n
      ntouched=0
      do ka=left%row_offsets(i),left%row_offsets(i+1)-1
        p=left%column_ids(ka)
        do kb=right%row_offsets(p),right%row_offsets(p+1)-1
          j=right%column_ids(kb)
          if(.not.marked(j))then;ntouched=ntouched+1;touched(ntouched)=j;marked(j)=.true.;end if
          accumulator(j)=accumulator(j)+left%values(ka)*right%values(kb)
        end do
      end do
      do p=2,ntouched
        tmp=touched(p);q=p-1
        do while(q>=1.and.touched(q)>tmp);touched(q+1)=touched(q);q=q-1;end do
        touched(q+1)=tmp
      end do
      do p=1,ntouched
        j=touched(p)
        if(abs(accumulator(j))>drop_tolerance)then
          if(nnz>=fill_limit)then;message='sparse multiplication fill limit exceeded';return;end if
          if(nnz==capacity)then;message='sparse multiplication workspace bound exceeded';return;end if
          nnz=nnz+1;columns(nnz)=j;values(nnz)=accumulator(j)
        end if
        accumulator(j)=(0d0,0d0);marked(j)=.false.
      end do
      offsets(i+1)=nnz+1
    end do
    result%extent=n;allocate(result%row_offsets(n+1),result%column_ids(nnz),result%values(nnz))
    result%row_offsets=offsets;result%column_ids=columns(:nnz);result%values=values(:nnz)
    ok=.true.;message='ok'
  end subroutine

  subroutine sparse_three_identity_minus(matrix,result)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    type(spatial_seed_sparse_matrix),intent(out)::result
    complex(real64),allocatable::dense_row(:)
    integer::i,k,j,nnz
    allocate(dense_row(matrix%extent));nnz=size(matrix%values)+matrix%extent
    result%extent=matrix%extent;allocate(result%row_offsets(matrix%extent+1),&
      result%column_ids(nnz),result%values(nnz));nnz=0;result%row_offsets(1)=1
    do i=1,matrix%extent
      dense_row=(0d0,0d0);dense_row(i)=3d0
      do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
        dense_row(matrix%column_ids(k))=dense_row(matrix%column_ids(k))-matrix%values(k)
      end do
      do j=1,matrix%extent
        if(abs(dense_row(j))>0d0)then;nnz=nnz+1;result%column_ids(nnz)=j;result%values(nnz)=dense_row(j);end if
      end do
      result%row_offsets(i+1)=nnz+1
    end do
    if(nnz<size(result%values))then
      result%column_ids=result%column_ids(:nnz);result%values=result%values(:nnz)
    end if
  end subroutine

  subroutine sparse_scale(matrix,scale)
    type(spatial_seed_sparse_matrix),intent(inout)::matrix;real(real64),intent(in)::scale
    matrix%values=scale*matrix%values
  end subroutine

  subroutine sparse_apply(matrix,x,y)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    complex(real64),intent(in)::x(:,:)
    complex(real64),allocatable,intent(out)::y(:,:)
    integer::i,k
    allocate(y(matrix%extent,size(x,2)));y=(0d0,0d0)
    do i=1,matrix%extent;do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
      y(i,:)=y(i,:)+matrix%values(k)*x(matrix%column_ids(k),:)
    end do;end do
  end subroutine

  subroutine sparse_apply_vector(matrix,x,y)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    complex(real64),intent(in)::x(:)
    complex(real64),allocatable,intent(out)::y(:)
    integer::i,k
    allocate(y(matrix%extent));y=(0d0,0d0)
    do i=1,matrix%extent;do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
      y(i)=y(i)+matrix%values(k)*x(matrix%column_ids(k))
    end do;end do
  end subroutine

  subroutine certified_sparse_metric_solve(s,rhs,tolerance,max_iterations,solution,iterations,residual,ok,message)
    type(spatial_seed_sparse_matrix),intent(in)::s
    complex(real64),intent(in)::rhs(:)
    real(real64),intent(in)::tolerance
    integer,intent(in)::max_iterations
    complex(real64),allocatable,intent(out)::solution(:)
    integer,intent(out)::iterations
    real(real64),intent(out)::residual
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::r(:),p(:),ap(:)
    real(real64)::rr,rr_new,alpha,beta,rhs_norm
    ok=.false.;message='';iterations=0;residual=huge(1d0)
    if(size(rhs)/=s%extent.or.max_iterations<1.or.tolerance<=0d0)then;message='invalid sparse metric solve contract';return;end if
    allocate(solution(s%extent));solution=(0d0,0d0);r=rhs;p=r
    rr=real(dot_product(r,r),real64);rhs_norm=max(1d0,sqrt(rr))
    if(rr==0d0)then;residual=0d0;ok=.true.;message='ok';return;end if
    do iterations=1,max_iterations
      call sparse_apply_vector(s,p,ap)
      alpha=rr/real(dot_product(p,ap),real64)
      if(alpha<=0d0)then;message='sparse metric solve lost positive definiteness';return;end if
      solution=solution+alpha*p;r=r-alpha*ap;rr_new=real(dot_product(r,r),real64)
      residual=sqrt(rr_new)/rhs_norm
      if(residual<=tolerance)then;ok=.true.;message='ok';return;end if
      beta=rr_new/rr;p=r+beta*p;rr=rr_new
    end do
    message='certified sparse metric solve iteration bound exceeded'
  end subroutine

  subroutine sparse_adjoint(matrix,adjoint)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    type(spatial_seed_sparse_matrix),intent(out)::adjoint
    integer,allocatable::counts(:),next(:);integer::i,k,j
    allocate(counts(matrix%extent),next(matrix%extent));counts=0
    do k=1,size(matrix%values);counts(matrix%column_ids(k))=counts(matrix%column_ids(k))+1;end do
    adjoint%extent=matrix%extent;allocate(adjoint%row_offsets(matrix%extent+1),&
      adjoint%column_ids(size(matrix%values)),adjoint%values(size(matrix%values)))
    adjoint%row_offsets(1)=1
    do i=1,matrix%extent;adjoint%row_offsets(i+1)=adjoint%row_offsets(i)+counts(i);end do
    next=adjoint%row_offsets(:matrix%extent)
    do i=1,matrix%extent;do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
      j=matrix%column_ids(k);adjoint%column_ids(next(j))=i;adjoint%values(next(j))=conjg(matrix%values(k));next(j)=next(j)+1
    end do;end do
  end subroutine

  real(real64) function sparse_identity_defect(matrix)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    integer::i,k,j;real(real64)::value
    sparse_identity_defect=0d0
    do i=1,matrix%extent
      value=0d0
      do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
        j=matrix%column_ids(k)
        if(j==i)then;value=abs(matrix%values(k)-1d0);else;value=abs(matrix%values(k));end if
        sparse_identity_defect=max(sparse_identity_defect,value)
      end do
      if(.not.any(matrix%column_ids(matrix%row_offsets(i):matrix%row_offsets(i+1)-1)==i))&
        sparse_identity_defect=max(sparse_identity_defect,1d0)
    end do
  end function

  subroutine sparse_hermitian_bounds(matrix,lower,upper,ok)
    type(spatial_seed_sparse_matrix),intent(in)::matrix
    real(real64),intent(out)::lower,upper
    logical,intent(out)::ok
    integer::i,k,j,reverse
    real(real64)::diagonal,radius,scale
    lower=huge(1d0);upper=-huge(1d0);scale=max(1d0,maxval(abs(matrix%values)));ok=.true.
    do i=1,matrix%extent
      diagonal=0d0;radius=0d0
      do k=matrix%row_offsets(i),matrix%row_offsets(i+1)-1
        j=matrix%column_ids(k)
        if(j==i)then
          if(abs(aimag(matrix%values(k)))>1d-12*scale)ok=.false.
          diagonal=real(matrix%values(k),real64)
        else
          radius=radius+abs(matrix%values(k));reverse=find_sparse_entry(matrix,j,i)
          if(reverse==0.or.abs(matrix%values(k)-conjg(matrix%values(reverse)))>1d-12*scale)ok=.false.
        end if
      end do
      lower=min(lower,diagonal-radius);upper=max(upper,diagonal+radius)
    end do
  end subroutine

  integer function find_sparse_entry(matrix,row,column)
    type(spatial_seed_sparse_matrix),intent(in)::matrix;integer,intent(in)::row,column
    integer::k
    find_sparse_entry=0
    do k=matrix%row_offsets(row),matrix%row_offsets(row+1)-1
      if(matrix%column_ids(k)==column)then;find_sparse_entry=k;return;end if
    end do
  end function

  subroutine sparse_filtered_block(h,x,xadj,y,edge,opposite,degree,result)
    type(spatial_seed_sparse_matrix),intent(in)::h,x,xadj
    complex(real64),intent(in)::y(:,:)
    real(real64),intent(in)::edge,opposite
    integer,intent(in)::degree
    complex(real64),allocatable,intent(out)::result(:,:)
    complex(real64),allocatable::first(:,:),second(:,:),third(:,:)
    integer::iteration
    result=y
    do iteration=1,degree
      call sparse_apply(x,result,first);call sparse_apply(h,first,second);call sparse_apply(xadj,second,third)
      result=(edge*result-third)/max(1d0,abs(edge-opposite))
    end do
    if(maxval(abs(result))>0d0)result=result/maxval(abs(result))
  end subroutine

  subroutine sparse_gram(s,left,right,gram)
    type(spatial_seed_sparse_matrix),intent(in)::s
    complex(real64),intent(in)::left(:,:),right(:,:)
    complex(real64),allocatable,intent(out)::gram(:,:)
    complex(real64),allocatable::sright(:,:)
    call sparse_apply(s,right,sright);gram=matmul(conjg(transpose(left)),sright)
  end subroutine

  subroutine subtract_identity(matrix)
    complex(real64),intent(inout)::matrix(:,:);integer::i
    do i=1,min(size(matrix,1),size(matrix,2));matrix(i,i)=matrix(i,i)-1d0;end do
  end subroutine

  subroutine reynolds_sparse(s,y,bg,d,rank_tolerance,c,ok,message)
    type(spatial_seed_sparse_matrix),intent(in)::s
    complex(real64),intent(in)::y(:,:),bg(:,:,:),d(:,:,:)
    real(real64),intent(in)::rank_tolerance
    complex(real64),allocatable,intent(out)::c(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::z(:,:),gram(:,:),vectors(:,:),work(:)
    real(real64),allocatable::eigenvalue(:),rwork(:)
    integer::g,m,info,lwork,i
    m=size(y,2);allocate(z(size(y,1),m));z=(0d0,0d0)
    do g=1,size(bg,3);z=z+matmul(bg(:,:,g),matmul(y,conjg(transpose(d(:,:,g)))));end do
    z=z/real(size(bg,3),real64);call sparse_gram(s,z,z,gram)
    allocate(eigenvalue(m),rwork(max(1,3*m-2)));lwork=max(1,2*m);allocate(work(lwork))
    call zheev('V','U',m,gram,m,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or.minval(eigenvalue)<=rank_tolerance)then;ok=.false.;message='complete irrep rank revelation failed';return;end if
    vectors=gram;do i=1,m;vectors(:,i)=vectors(:,i)/sqrt(eigenvalue(i));end do
    gram=matmul(vectors,conjg(transpose(gram)));c=matmul(z,gram);ok=.true.;message='ok'
  end subroutine

  real(real64) function sparse_ritz_defect(h,s,c)
    type(spatial_seed_sparse_matrix),intent(in)::h,s
    complex(real64),intent(in)::c(:,:)
    complex(real64),allocatable::hc(:,:),sc(:,:),lambda(:,:)
    call sparse_apply(h,c,hc);call sparse_apply(s,c,sc);lambda=matmul(conjg(transpose(c)),hc)
    sparse_ritz_defect=maxval(abs(hc-matmul(sc,lambda)))/max(1d0,maxval(abs(hc)))
  end function

  real(real64) function intertwining_defect_sparse(c,bg,d)
    complex(real64),intent(in)::c(:,:),bg(:,:,:),d(:,:,:);integer::g
    intertwining_defect_sparse=0d0
    do g=1,size(bg,3);intertwining_defect_sparse=max(intertwining_defect_sparse,&
      maxval(abs(matmul(bg(:,:,g),c)-matmul(c,d(:,:,g)))));end do
  end function

  real(real64) function gram_defect_sparse(s,c,a)
    type(spatial_seed_sparse_matrix),intent(in)::s
    complex(real64),intent(in)::c(:,:),a(:,:)
    complex(real64),allocatable::gram(:,:),ac(:,:)
    call sparse_gram(s,c,c,gram);ac=matmul(a,c)
    gram_defect_sparse=maxval(abs(gram-matmul(conjg(transpose(ac)),ac)))
  end function

  subroutine construct_spatial_mlwf_seed_dense(comm,h,s,yocc,yempty,bg,docc,dempty,a,tg,opt,&
      cocc,cempty,diag,fingerprint,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::h(:,:),s(:,:),yocc(:,:),yempty(:,:),bg(:,:,:),docc(:,:,:),dempty(:,:,:),a(:,:),tg(:,:,:)
    type(spatial_seed_options),intent(in)::opt
    complex(real64),allocatable,intent(out)::cocc(:,:),cempty(:,:)
    type(spatial_seed_diagnostics),intent(out)::diag
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::n,info,i,j,g,ierr,local_bad,global_bad
    complex(real64),allocatable::x(:,:),hbar(:,:),work(:,:),qfiltered(:,:),identity(:,:),projector(:,:),unitary(:,:)
    real(real64)::scale,defect,minimum_pivot,maximum_pivot,verified_lower,verified_upper,radius
    logical::local_ok
    character(256)::local_message
    ok=.false.;message='';fingerprint=0_int64;diag=spatial_seed_diagnostics();n=size(s,1)
    if(size(s,2)/=n.or.size(h,1)/=n.or.size(h,2)/=n.or.size(a,1)/=n.or.size(a,2)/=n)then
      message='seed matrix shape mismatch';return
    end if
    if(maxval(abs(s-conjg(transpose(s))))>opt%tolerance)then;message='S is not Hermitian positive definite';return;end if
    if(maxval(abs(h-conjg(transpose(h))))>opt%tolerance)then;message='H is not Hermitian';return;end if
    if(opt%empty_lower<=opt%occupied_upper)then;message='occupied-empty gap is closed';return;end if
    if(opt%filter_degree<1)then;message='filter degree must be positive';return;end if
    allocate(work(n,n));work=s
    call zpotrf('U',n,work,n,info)
    if(info/=0)then;message='S is not Hermitian positive definite';return;end if
    minimum_pivot=minval([(abs(work(i,i)),i=1,n)]);maximum_pivot=maxval([(abs(work(i,i)),i=1,n)])
    diag%metric_condition=(maximum_pivot/minimum_pivot)**2
    if(diag%metric_condition>opt%condition_limit)then;message='S condition limit exceeded';return;end if
    allocate(x(n,n),identity(n,n));identity=(0d0,0d0)
    do i=1,n;identity(i,i)=1d0;end do
    scale=maxval(sum(abs(s),dim=2));x=identity/sqrt(scale)
    do i=1,opt%max_metric_iterations
      x=0.5d0*matmul(x,3d0*identity-matmul(s,matmul(x,x)))
      diag%metric_defect=maxval(abs(matmul(conjg(transpose(x)),matmul(s,x))-identity))
      diag%metric_iterations=i
      if(diag%metric_defect<=opt%tolerance)exit
    end do
    if(diag%metric_defect>opt%tolerance)then;message='metric purification did not converge';return;end if
    if(count(abs(x)>opt%rank_tolerance)>opt%fill_limit)then;message='metric purification fill limit exceeded';return;end if
    allocate(hbar(n,n));hbar=matmul(conjg(transpose(x)),matmul(h,x))
    verified_lower=huge(1d0);verified_upper=-huge(1d0)
    do i=1,n
      radius=sum(abs(hbar(i,:)))-abs(hbar(i,i))
      verified_lower=min(verified_lower,real(hbar(i,i),real64)-radius)
      verified_upper=max(verified_upper,real(hbar(i,i),real64)+radius)
    end do
    if(opt%lower_bound>verified_lower+opt%tolerance.or.opt%upper_bound<verified_upper-opt%tolerance)then
      message='spectral bounds do not enclose verified H_bar bounds';return
    end if
    call filtered_block(hbar,yocc,opt%occupied_upper,opt%lower_bound,opt%filter_degree,qfiltered)
    work=matmul(x,qfiltered)
    diag%back_transform_defect=maxval(abs(work-matmul(x,qfiltered)))
    call reynolds_and_orthonormalize(s,work,bg,docc,opt%rank_tolerance,cocc,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    call filtered_block(hbar,yempty,opt%upper_bound,opt%empty_lower,opt%filter_degree,qfiltered)
    work=matmul(x,qfiltered)
    if(size(cocc,2)>0)work=work-matmul(cocc,matmul(conjg(transpose(cocc)),matmul(s,work)))
    call reynolds_and_orthonormalize(s,work,bg,dempty,opt%rank_tolerance,cempty,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    diag%filter_fill=count(abs(cocc)>opt%rank_tolerance)+count(abs(cempty)>opt%rank_tolerance)
    if(diag%filter_fill>opt%fill_limit)then;message='filter fill limit exceeded';return;end if
    call apply_s_projector(cocc,s,projector,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    diag%occupied_projector_defect=maxval(abs(matmul(projector,projector)-projector))
    call apply_s_projector(cempty,s,work,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;end if
    diag%empty_projector_defect=maxval(abs(matmul(work,work)-work))
    diag%mutual_projector_defect=maxval(abs(matmul(projector,work)))
    diag%ritz_residual=max(ritz_defect(h,s,cocc),ritz_defect(h,s,cempty))
    diag%intertwining_defect=max(intertwining_defect(s,cocc,bg,docc),intertwining_defect(s,cempty,bg,dempty))
    diag%map_metric_defect=max(gram_defect(s,cocc,a),gram_defect(s,cempty,a))
    allocate(unitary(max(size(cocc,2),size(cempty,2)),max(size(cocc,2),size(cempty,2))))
    unitary=(0d0,0d0);do i=1,size(unitary,1);unitary(i,i)=1d0;end do
    diag%map_linearity_defect=max(linearity_defect(a,cocc,unitary(:size(cocc,2),:size(cocc,2))),&
      linearity_defect(a,cempty,unitary(:size(cempty,2),:size(cempty,2))))
    diag%map_equivariance_defect=max(equivariance_defect(a,cocc,bg,tg),equivariance_defect(a,cempty,bg,tg))
    defect=max(diag%occupied_projector_defect,diag%empty_projector_defect,diag%mutual_projector_defect,&
      diag%ritz_residual,diag%intertwining_defect,diag%map_metric_defect,diag%map_linearity_defect,&
      diag%map_equivariance_defect)
    if(diag%map_metric_defect>opt%tolerance)then;message='map metric pullback failed';return;end if
    if(diag%map_linearity_defect>opt%tolerance)then;message='map linearity failed';return;end if
    if(diag%map_equivariance_defect>opt%tolerance)then;message='map equivariance failed';return;end if
    if(defect>10d0*opt%tolerance)then
      write(message,'(a,8es11.3)')'seed projector, Ritz, or symmetry gate failed: ',&
        diag%occupied_projector_defect,diag%empty_projector_defect,diag%mutual_projector_defect,&
        diag%ritz_residual,diag%intertwining_defect,diag%map_metric_defect,&
        diag%map_linearity_defect,diag%map_equivariance_defect
      return
    end if
    fingerprint=1469598103934665603_int64
    do j=1,size(cocc,2);do i=1,size(cocc,1)
      fingerprint=ieor(fingerprint,transfer(real(cocc(i,j),real64),fingerprint))
      fingerprint=fingerprint*1099511628211_int64
    end do;end do
    local_bad=0;call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='MPI seed agreement failed';return;end if
    ok=.true.;message='ok'
  end subroutine

  subroutine filtered_block(hbar,y,edge,opposite,degree,result)
    complex(real64),intent(in)::hbar(:,:),y(:,:)
    real(real64),intent(in)::edge,opposite
    integer,intent(in)::degree
    complex(real64),allocatable,intent(out)::result(:,:)
    integer::i,n
    complex(real64),allocatable::operator(:,:)
    n=size(hbar,1);allocate(operator(n,n));operator=-hbar
    do i=1,n;operator(i,i)=operator(i,i)+edge;end do
    operator=operator/max(abs(edge-opposite),1d0);result=y
    do i=1,degree;result=matmul(operator,result);end do
    if(maxval(abs(result))>0d0)result=result/maxval(abs(result))
  end subroutine

  subroutine reynolds_and_orthonormalize(s,y,bg,d,rank_tolerance,c,ok,message)
    complex(real64),intent(in)::s(:,:),y(:,:),bg(:,:,:),d(:,:,:)
    real(real64),intent(in)::rank_tolerance
    complex(real64),allocatable,intent(out)::c(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::z(:,:),gram(:,:),vectors(:,:),work(:)
    real(real64),allocatable::eigenvalue(:),rwork(:)
    integer::g,m,n,info,lwork,i
    n=size(y,1);m=size(y,2);allocate(z(n,m));z=(0d0,0d0)
    do g=1,size(bg,3);z=z+matmul(bg(:,:,g),matmul(y,conjg(transpose(d(:,:,g)))));end do
    z=z/real(size(bg,3),real64);allocate(gram(m,m));gram=matmul(conjg(transpose(z)),matmul(s,z))
    allocate(eigenvalue(m),rwork(max(1,3*m-2)));lwork=max(1,2*m);allocate(work(lwork))
    call zheev('V','U',m,gram,m,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or.minval(eigenvalue)<=rank_tolerance)then;ok=.false.;message='complete irrep rank revelation failed';return;end if
    vectors=gram
    do i=1,m;vectors(:,i)=vectors(:,i)/sqrt(eigenvalue(i));end do
    gram=matmul(vectors,conjg(transpose(gram)))
    c=matmul(z,gram)
    ok=.true.;message='ok'
  end subroutine

  subroutine apply_s_projector(c,s,p,ok,message)
    complex(real64),intent(in)::c(:,:),s(:,:)
    complex(real64),allocatable,intent(out)::p(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:)
    integer::m,info
    m=size(c,2);allocate(gram(m,m));gram=matmul(conjg(transpose(c)),matmul(s,c))
    call zpotrf('U',m,gram,m,info)
    if(info/=0)then;ok=.false.;message='projector Gram is singular';return;end if
    call zpotri('U',m,gram,m,info)
    if(info/=0)then;ok=.false.;message='projector Gram inverse failed';return;end if
    call fill_lower(gram);p=matmul(c,matmul(gram,matmul(conjg(transpose(c)),s)))
    ok=.true.;message='ok'
  end subroutine

  subroutine fill_lower(matrix)
    complex(real64),intent(inout)::matrix(:,:);integer::i,j
    do j=1,size(matrix,2);do i=j+1,size(matrix,1);matrix(i,j)=conjg(matrix(j,i));end do;end do
  end subroutine

  real(real64) function ritz_defect(h,s,c)
    complex(real64),intent(in)::h(:,:),s(:,:),c(:,:)
    complex(real64),allocatable::lambda(:,:),residual(:,:)
    lambda=matmul(conjg(transpose(c)),matmul(h,c));residual=matmul(h,c)-matmul(s,matmul(c,lambda))
    ritz_defect=maxval(abs(residual))/max(1d0,maxval(abs(matmul(h,c))))
  end function
  real(real64) function intertwining_defect(s,c,bg,d)
    complex(real64),intent(in)::s(:,:),c(:,:),bg(:,:,:),d(:,:,:);integer::g
    intertwining_defect=0d0
    do g=1,size(bg,3);intertwining_defect=max(intertwining_defect,maxval(abs(matmul(bg(:,:,g),c)-matmul(c,d(:,:,g)))));end do
  end function
  real(real64) function gram_defect(s,c,a)
    complex(real64),intent(in)::s(:,:),c(:,:),a(:,:)
    gram_defect=maxval(abs(matmul(conjg(transpose(c)),matmul(s,c))-&
      matmul(conjg(transpose(matmul(a,c))),matmul(a,c))))
  end function
  real(real64) function linearity_defect(a,c,u)
    complex(real64),intent(in)::a(:,:),c(:,:),u(:,:)
    linearity_defect=maxval(abs(matmul(a,matmul(c,u))-matmul(matmul(a,c),u)))
  end function
  real(real64) function equivariance_defect(a,c,bg,tg)
    complex(real64),intent(in)::a(:,:),c(:,:),bg(:,:,:),tg(:,:,:);integer::g
    equivariance_defect=0d0
    do g=1,size(bg,3);equivariance_defect=max(equivariance_defect,&
      maxval(abs(matmul(a,matmul(bg(:,:,g),c))-matmul(tg(:,:,g),matmul(a,c)))));end do
  end function
end module
