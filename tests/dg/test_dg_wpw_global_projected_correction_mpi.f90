program test_dg_wpw_global_projected_correction_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_global_projected_correction
  implicit none
  type s_fixture_context
    integer::rank=0,epoch=7,apply_calls=0,gram_calls=0,fail_apply_call=0,stale_apply_call=0,&
      stale_fingerprint_apply_call=0,disagree_gram_call=0
    integer(8)::fingerprint=41_8
  end type
  type(s_fixture_context)::ctx
  type(s_dg_wpw_global_correction_controls)::controls
  type(s_dg_wpw_global_correction_diagnostics)::diagnostics
  complex(8)::qw(1,2),qp(1,2),rw(1,3),rp(1,3),zw(1,3),zp(1,3)
  real(8)::eigenvalues(3)
  integer::ierr,info

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,ctx%rank,ierr)
  call initialize_problem(qw,qp,rw,rp,eigenvalues)
  controls%restart=4;controls%max_iterations=5;controls%state_batch=2
  controls%relative_tolerance=1d-11
  call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
    global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
  if(info/=0.or.any(.not.diagnostics%converged).or.maxval(diagnostics%s_orthogonality)>1d-11.or.&
    maxval(diagnostics%equation_defect)>1d-9)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  call check_dense_oracle
  if(maxval(diagnostics%iterations)>5)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  if(.not.any(diagnostics%breakdown_status==DG_WPW_GLOBAL_BREAKDOWN_HAPPY))&
    call MPI_Abort(MPI_COMM_WORLD,25,ierr)
  call check_matching_diagnostics
  controls%restart=1;controls%max_iterations=5;controls%relative_tolerance=1d-2
  call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
    global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
  if(info/=0.or..not.any(diagnostics%restart_count>0))call MPI_Abort(MPI_COMM_WORLD,27,ierr)
  controls%restart=4;controls%max_iterations=5;controls%relative_tolerance=1d-11

  ! Zero projected RHS and happy breakdown are accepted using explicit residuals.
  call apply_s(ctx,qw(:,1:1),qp(:,1:1),rw(:,1:1),rp(:,1:1),info)
  call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
    global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
  if(info/=0.or..not.diagnostics%zero_rhs(1).or.maxval(abs(zw(:,1)))+maxval(abs(zp(:,1)))>0d0)&
    call MPI_Abort(MPI_COMM_WORLD,3,ierr)

  call check_collective_failures
  call release_dg_wpw_global_correction_diagnostics(diagnostics)
  call release_dg_wpw_global_correction_diagnostics(diagnostics)
  if(allocated(diagnostics%converged))call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  if(ctx%rank==0)write(*,'(a)')"PASS global projected correction MPI fixture"
  call MPI_Finalize(ierr)
contains
  subroutine check_dense_oracle
    complex(8)::hmat(4,4),smat(4,4),qfull(4,2),qlocal(4,2),rfull(4,3),rlocal(4,3),&
      zfull(4,3),zlocal(4,3),identity(4,4),pr(4,4),pl(4,4),sq(4,2),op(4,4),&
      kkt(6,6),rhs_kkt(6,1),b(4),x(4),right_formula(4),left_formula(4),operator_formula(4)
    integer::s,j,lapack_info,mpi_err,ipiv(6)
    hmat=reshape([&
      (1.5d0,0d0),(0.2d0,-0.1d0),(0.3d0,0.2d0),(0.0d0,0.1d0),&
      (0.2d0,0.1d0),(2.2d0,0d0),(-0.15d0,0.05d0),(0.25d0,0d0),&
      (0.3d0,-0.2d0),(-0.15d0,-0.05d0),(3.1d0,0d0),(0.4d0,-0.1d0),&
      (0.0d0,-0.1d0),(0.25d0,0d0),(0.4d0,0.1d0),(4.0d0,0d0)],[4,4])
    smat=reshape([&
      (2d0,0d0),(0.1d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0.1d0,0d0),(1.7d0,0d0),(0.05d0,0d0),(0d0,0d0),&
      (0d0,0d0),(0.05d0,0d0),(2.3d0,0d0),(0.1d0,0d0),&
      (0d0,0d0),(0d0,0d0),(0.1d0,0d0),(1.9d0,0d0)],[4,4])
    qlocal=0;rlocal=0;zlocal=0
    qlocal(1+2*ctx%rank,:)=qw(1,:);qlocal(2+2*ctx%rank,:)=qp(1,:)
    rlocal(1+2*ctx%rank,:)=rw(1,:);rlocal(2+2*ctx%rank,:)=rp(1,:)
    zlocal(1+2*ctx%rank,:)=zw(1,:);zlocal(2+2*ctx%rank,:)=zp(1,:)
    call MPI_Allreduce(qlocal,qfull,8,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(rlocal,rfull,12,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    call MPI_Allreduce(zlocal,zfull,12,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    identity=0;do j=1,4;identity(j,j)=1;enddo
    sq=matmul(smat,qfull)
    pr=identity-matmul(qfull,matmul(conjg(transpose(qfull)),smat))
    pl=identity-matmul(sq,conjg(transpose(qfull)))
    x=[(0.3d0,0.2d0),(-0.4d0,0.1d0),(0.8d0,-0.2d0),(0.5d0,0.6d0)]
    right_formula=x-matmul(qfull,matmul(conjg(transpose(qfull)),matmul(smat,x)))
    left_formula=x-matmul(sq,matmul(conjg(transpose(qfull)),x))
    if(maxval(abs(right_formula-matmul(pr,x)))>1d-13.or.&
      maxval(abs(left_formula-matmul(pl,x)))>1d-13.or.&
      maxval(abs(matmul(conjg(transpose(qfull)),matmul(smat,pr))))>1d-13)&
      call MPI_Abort(MPI_COMM_WORLD,21,ierr)
    do s=1,3
      op=matmul(pl,matmul(hmat-eigenvalues(s)*smat,pr))
      operator_formula=matmul(pl,matmul(hmat-eigenvalues(s)*smat,matmul(pr,x)))
      if(maxval(abs(operator_formula-matmul(op,x)))>1d-13)call MPI_Abort(MPI_COMM_WORLD,22,ierr)
      b=-matmul(pl,rfull(:,s))
      if(sqrt(sum(abs(b-matmul(op,zfull(:,s)))**2))>1d-9)&
        call MPI_Abort(MPI_COMM_WORLD,23,ierr)
      kkt=0;kkt(1:4,1:4)=hmat-eigenvalues(s)*smat;kkt(1:4,5:6)=sq
      kkt(5:6,1:4)=matmul(conjg(transpose(qfull)),smat)
      rhs_kkt=0;rhs_kkt(1:4,1)=b
      call zgesv(6,1,kkt,6,ipiv,rhs_kkt,6,lapack_info)
      if(lapack_info/=0.or.maxval(abs(rhs_kkt(1:4,1)-zfull(:,s)))>1d-9)&
        call MPI_Abort(MPI_COMM_WORLD,24,ierr)
    enddo
  end subroutine

  subroutine initialize_problem(ow,op,brw,brp,eval)
    complex(8),intent(out)::ow(1,2),op(1,2),brw(1,3),brp(1,3)
    real(8),intent(out)::eval(3)
    ! Distributed dense oracle: global row order is W0,P0,W1,P1.
    ow=0;op=0
    if(ctx%rank==0)then
      ow(1,1)=1d0/sqrt(2d0)
    else
      op(1,2)=1d0/sqrt(1.9d0)
    endif
    brw(1,:)=[cmplx(0.3d0+ctx%rank,0.2d0,8),cmplx(-0.7d0,0.1d0*ctx%rank,8),&
      cmplx(0.9d0-0.2d0*ctx%rank,-0.3d0,8)]
    brp(1,:)=[cmplx(-0.2d0,0.4d0,8),cmplx(0.8d0-ctx%rank,0.2d0,8),&
      cmplx(-0.4d0+0.3d0*ctx%rank,0.5d0,8)]
    eval=[0.2d0,0.65d0,1.1d0]
  end subroutine

  subroutine check_collective_failures
    complex(8),allocatable::badw(:,:),badp(:,:),badqw(:,:),badqp(:,:),&
      count_rw(:,:),count_rp(:,:),count_zw(:,:),count_zp(:,:)
    real(8),allocatable::count_eval(:)
    complex(8)::saved
    integer::selected,ncount
    controls%restart=0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,10);controls%restart=4
    controls%restart=17
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,101);controls%restart=4
    controls%max_iterations=0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,102);controls%max_iterations=5
    controls%state_batch=17
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,103);controls%state_batch=2
    controls%relative_tolerance=1d0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,104);controls%relative_tolerance=1d-11
    allocate(badw(merge(2,1,ctx%rank==0),3),badp(1,3));badw=0;badp=0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,badw,badp,controls,badw,badp,diagnostics,info)
    call require_failure(info,11);deallocate(badw,badp)
    allocate(badqw(1,merge(1,2,ctx%rank==0)),badqp(1,merge(1,2,ctx%rank==0)))
    badqw=0;badqp=0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,badqw,badqp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,111);deallocate(badqw,badqp)
    eigenvalues(2)=ieee_value(0d0,ieee_quiet_nan)
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,12);eigenvalues(2)=0.65d0
    saved=qw(1,1);if(ctx%rank==1)qw(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,121);qw(1,1)=saved
    saved=rp(1,2);if(ctx%rank==0)rp(1,2)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_failure(info,122);rp(1,2)=saved
    ncount=merge(2,3,ctx%rank==0)
    allocate(count_eval(ncount),count_rw(1,ncount),count_rp(1,ncount),&
      count_zw(1,ncount),count_zp(1,ncount))
    count_eval=0;count_rw=0;count_rp=0
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,count_eval,count_rw,count_rp,controls,&
      count_zw,count_zp,diagnostics,info)
    call require_failure(info,123)
    deallocate(count_eval,count_rw,count_rp,count_zw,count_zp)
    ctx%apply_calls=0;ctx%fail_apply_call=1
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_initialized_failure(info,13);ctx%fail_apply_call=0
    ctx%apply_calls=0;ctx%stale_apply_call=2
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_initialized_failure(info,14);ctx%stale_apply_call=0;ctx%epoch=7
    ctx%apply_calls=0;ctx%stale_fingerprint_apply_call=2
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_initialized_failure(info,141)
    ctx%stale_fingerprint_apply_call=0;ctx%fingerprint=41_8
    controls%restart=2;controls%max_iterations=3;controls%relative_tolerance=1d-15
    call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
      global_gram,validate_state,7,41_8,qw(:,1:1),qp(:,1:1),eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
    call require_initialized_failure(info,15)
    if(.not.any(diagnostics%iterations==3).or..not.any(diagnostics%restart_count==1))&
      call MPI_Abort(MPI_COMM_WORLD,151,ierr)
    controls%restart=4;controls%max_iterations=5;controls%relative_tolerance=1d-11
    do selected=1,80
      call initialize_problem(qw,qp,rw,rp,eigenvalues)
      ctx%gram_calls=0;ctx%disagree_gram_call=selected
      call solve_dg_wpw_global_projected_correction(ctx,MPI_COMM_WORLD,apply_h,apply_s,&
        global_gram,validate_state,7,41_8,qw,qp,eigenvalues,rw,rp,controls,zw,zp,diagnostics,info)
      if(info/=0.and.allocated(diagnostics%state_status))then
        if(any(diagnostics%state_status==DG_WPW_GLOBAL_STATE_INVALID_RESULT))exit
      endif
    enddo
    if(selected>80)call MPI_Abort(MPI_COMM_WORLD,16,ierr)
    call require_initialized_failure(info,16);ctx%disagree_gram_call=0
  end subroutine check_collective_failures
  subroutine check_matching_diagnostics
    call check_real_match(diagnostics%initial_residual,261)
    call check_real_match(diagnostics%final_residual,262)
    call check_real_match(diagnostics%relative_residual,263)
    call check_real_match(diagnostics%s_orthogonality,264)
    call check_real_match(diagnostics%equation_defect,265)
    call check_real_match(diagnostics%projected_fraction,266)
    call check_real_match(diagnostics%correction_norm,267)
    call check_real_match(diagnostics%amplification,268)
    call check_int_match(diagnostics%iterations,269)
    call check_int_match(diagnostics%restart_count,270)
    call check_int_match(diagnostics%breakdown_status,271)
    call check_int_match(diagnostics%state_status,272)
  end subroutine
  subroutine check_real_match(values,code)
    real(8),intent(in)::values(:);integer,intent(in)::code
    real(8)::min_value,max_value;integer::s,mpi_err
    do s=1,size(values)
      call MPI_Allreduce(values(s),min_value,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,mpi_err)
      call MPI_Allreduce(values(s),max_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_err)
      if(abs(max_value-min_value)>1d-13)call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    enddo
  end subroutine
  subroutine check_int_match(values,code)
    integer,intent(in)::values(:),code
    integer::min_value,max_value,s,mpi_err
    do s=1,size(values)
      call MPI_Allreduce(values(s),min_value,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_err)
      call MPI_Allreduce(values(s),max_value,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_err)
      if(min_value/=max_value)call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    enddo
  end subroutine
  subroutine require_failure(value,code)
    integer,intent(in)::value,code
    integer::same_info,min_info,max_info
    call MPI_Allreduce(value,min_info,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,same_info)
    call MPI_Allreduce(value,max_info,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,same_info)
    if(value==0.or.min_info/=max_info)call MPI_Abort(MPI_COMM_WORLD,code,ierr)
  end subroutine
  subroutine require_initialized_failure(value,code)
    integer,intent(in)::value,code
      call require_failure(value,code)
      if(.not.allocated(diagnostics%state_status).or.any(zw/=(0d0,0d0)).or.any(zp/=(0d0,0d0)))&
        call MPI_Abort(MPI_COMM_WORLD,code+100,ierr)
      call check_matching_diagnostics
  end subroutine

  subroutine apply_h(context,xw,xp,yw,yp,apply_info)
    class(*),intent(inout)::context
    complex(8),intent(in)::xw(:,:),xp(:,:)
    complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::apply_info
    call dense_apply(context,xw,xp,yw,yp,.true.,apply_info)
  end subroutine
  subroutine apply_s(context,xw,xp,yw,yp,apply_info)
    class(*),intent(inout)::context
    complex(8),intent(in)::xw(:,:),xp(:,:)
    complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::apply_info
    call dense_apply(context,xw,xp,yw,yp,.false.,apply_info)
  end subroutine
  subroutine dense_apply(context,xw,xp,yw,yp,use_h,apply_info)
    class(*),intent(inout)::context
    complex(8),intent(in)::xw(:,:),xp(:,:)
    complex(8),intent(out)::yw(:,:),yp(:,:)
    logical,intent(in)::use_h
    integer,intent(out)::apply_info
    complex(8)::local(4,size(xw,2)),global(4,size(xw,2)),matrix(4,4)
    integer::mpi_info
    select type(c=>context);type is(s_fixture_context)
      c%apply_calls=c%apply_calls+1
      if(c%stale_apply_call>0.and.c%apply_calls==c%stale_apply_call)c%epoch=c%epoch+1
      if(c%stale_fingerprint_apply_call>0.and.c%apply_calls==c%stale_fingerprint_apply_call)&
        c%fingerprint=c%fingerprint+1_8
      apply_info=merge(7,0,c%fail_apply_call>0.and.c%apply_calls==c%fail_apply_call)
    class default;apply_info=9;return;end select
    local=0;local(1+2*ctx%rank,:)=xw(1,:);local(2+2*ctx%rank,:)=xp(1,:)
    call MPI_Allreduce(local,global,4*size(xw,2),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_info)
    if(use_h)then
      matrix=reshape([&
        (1.5d0,0d0),(0.2d0,-0.1d0),(0.3d0,0.2d0),(0.0d0,0.1d0),&
        (0.2d0,0.1d0),(2.2d0,0d0),(-0.15d0,0.05d0),(0.25d0,0d0),&
        (0.3d0,-0.2d0),(-0.15d0,-0.05d0),(3.1d0,0d0),(0.4d0,-0.1d0),&
        (0.0d0,-0.1d0),(0.25d0,0d0),(0.4d0,0.1d0),(4.0d0,0d0)],[4,4])
    else
      matrix=reshape([&
        (2d0,0d0),(0.1d0,0d0),(0d0,0d0),(0d0,0d0),&
        (0.1d0,0d0),(1.7d0,0d0),(0.05d0,0d0),(0d0,0d0),&
        (0d0,0d0),(0.05d0,0d0),(2.3d0,0d0),(0.1d0,0d0),&
        (0d0,0d0),(0d0,0d0),(0.1d0,0d0),(1.9d0,0d0)],[4,4])
    endif
    local=matmul(matrix,global);yw(1,:)=local(1+2*ctx%rank,:);yp(1,:)=local(2+2*ctx%rank,:)
    if(mpi_info/=MPI_SUCCESS)apply_info=8
  end subroutine
  subroutine global_gram(left,right,nrow,nleft,nright,gram,gram_info)
    integer,intent(in)::nrow,nleft,nright
    complex(8),intent(in)::left(nrow,nleft),right(nrow,nright)
    complex(8),intent(out)::gram(nleft,nright)
    integer,intent(out)::gram_info
    complex(8)::local(nleft,nright)
    local=matmul(conjg(transpose(left)),right);ctx%gram_calls=ctx%gram_calls+1
    if(ctx%disagree_gram_call/=0.and.(ctx%disagree_gram_call<0.or.&
      ctx%gram_calls==ctx%disagree_gram_call))then
      gram=local;gram_info=0
    else
      call MPI_Allreduce(local,gram,nleft*nright,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,gram_info)
    endif
  end subroutine
  subroutine validate_state(context,expected_epoch,expected_fingerprint,validate_info)
    class(*),intent(inout)::context
    integer,intent(in)::expected_epoch
    integer(8),intent(in)::expected_fingerprint
    integer,intent(out)::validate_info
    select type(c=>context);type is(s_fixture_context)
      validate_info=merge(0,1,c%epoch==expected_epoch.and.c%fingerprint==expected_fingerprint)
    class default;validate_info=1;end select
  end subroutine
end program
