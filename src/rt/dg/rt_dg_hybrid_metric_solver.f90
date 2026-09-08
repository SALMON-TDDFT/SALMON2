#include "config.h"
module rt_dg_hybrid_metric_solver
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange,&
    exchange_rt_dg_sparse_values,clear_rt_dg_sparse_exchange
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::solve_rt_dg_hybrid_metric
contains
    subroutine solve_rt_dg_hybrid_metric(comm,metric,rhs_owned,relative_tolerance,max_iterations,solution_owned,&
      iteration_count,relative_residual,workspace_peak_bytes,fingerprint,ok,message,cached_exchange_plan)
    integer,intent(in)::comm,max_iterations
    type(s_dg_hybrid_sparse_metric),intent(in)::metric
    complex(real64),intent(in)::rhs_owned(:,:)
    real(real64),intent(in)::relative_tolerance
    complex(real64),allocatable,intent(out)::solution_owned(:,:)
    integer,intent(out)::iteration_count
    real(real64),intent(out)::relative_residual
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_rt_dg_sparse_exchange),optional,intent(inout)::cached_exchange_plan
#ifdef USE_MPI
    integer::i,j,k,row,nowned,nrhs,n,iter,ierr,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer
    type(s_rt_dg_sparse_exchange)::exchange_plan
    integer(int64)::bits,minimum_bits,maximum_bits,complex_count,real_count,integer_count,quantized,local_hash,global_hash,row_hash
    complex(real64),allocatable::r(:),z(:),p(:),ap(:),edge_values(:)
    complex(real64)::local_dot,global_dot
    real(real64),allocatable::diagonal(:)
    real(real64)::rho,rho_new,denominator,alpha,beta,bnorm,resnorm,target,local_real,solution_scale,&
      quantization_limit,rhs_scale,safe_rhs_scale,maximum_residual
    logical::converged,exchange_ok,owns_exchange
    character(256)::exchange_message
    ok=.false.;message='';iteration_count=0;relative_residual=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64
    owns_exchange=.true.
    nowned=size(metric%owned_row_ids);nrhs=size(rhs_owned,2);n=metric%global_count;local_bad=0
    call agree_integer(n,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid metric solver extent';return;endif
    call agree_integer(nrhs,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid metric RHS count';return;endif
    call agree_integer(max_iterations,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid metric iteration cap';return;endif
    bits=transfer(relative_tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent hybrid metric solver tolerance';return;endif
    call agree_int64(metric%fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.metric%fingerprint==0_int64)then
      message='inconsistent hybrid metric solver provenance';return
    endif
    if(.not.metric%valid.or.n<1.or.nrhs<1.or.max_iterations<1)local_bad=1
    if(nrhs>huge(0)-4)local_bad=1
    if(size(rhs_owned,1)/=nowned.or.size(metric%row_offsets)/=nowned+1)local_bad=1
    if(size(metric%active_rows)/=n.or.size(metric%packet_ids)/=n)local_bad=1
    if(size(metric%column_ids)/=size(metric%values))local_bad=1
    if(any(metric%owned_row_ids<1_int64).or.any(metric%owned_row_ids>int(n,int64)))local_bad=1
    if(.not.ieee_is_finite(relative_tolerance).or.relative_tolerance<1d-15.or.relative_tolerance>1d-2)local_bad=1
    if(.not.ieee_is_finite(metric%condition_estimate).or.metric%condition_estimate<1d0)local_bad=1
    if(.not.finite_matrix(rhs_owned).or..not.finite_vector(metric%values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid metric solver contract';return;endif
    complex_count=0_int64;real_count=0_int64;integer_count=0_int64
    call add_product(complex_count,int(nowned,int64),int(nrhs+4,int64),local_bad)
    call add_count(complex_count,int(size(metric%column_ids),int64),local_bad)
    call add_count(real_count,int(nowned,int64),local_bad)
    if(local_bad==0)then
      if(complex_count>huge(workspace_peak_bytes)/16_int64.or.real_count>huge(workspace_peak_bytes)/8_int64.or.&
        integer_count>huge(workspace_peak_bytes)/4_int64)local_bad=1
    endif
    if(local_bad==0)then
      if(16_int64*complex_count>huge(workspace_peak_bytes)-8_int64*real_count)local_bad=1
    endif
    if(local_bad==0)then
      if(16_int64*complex_count+8_int64*real_count>huge(workspace_peak_bytes)-4_int64*integer_count)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid metric solver workspace overflow';return;endif
    workspace_peak_bytes=16_int64*complex_count+8_int64*real_count+4_int64*integer_count
    allocate(solution_owned(nowned,nrhs),r(nowned),z(nowned),p(nowned),ap(nowned),edge_values(size(metric%column_ids)),&
      diagonal(nowned),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid metric solver workspace';return;endif
    owns_exchange=.not.present(cached_exchange_plan)
    if(present(cached_exchange_plan))then
      local_bad=merge(0,1,cached_exchange_plan%valid.and.&
        cached_exchange_plan%catalog_fingerprint==metric%fingerprint)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='cached metric exchange agreement failed';return;endif
      if(global_bad/=0)then;call cleanup();message='stale cached metric exchange plan';return;endif
    else
      call build_rt_dg_sparse_exchange(comm,n,metric%fingerprint,metric%owned_row_ids,metric%column_ids,&
        exchange_plan,exchange_ok,exchange_message)
      if(.not.exchange_ok)then;call cleanup();message='hybrid metric exchange failed: '//trim(exchange_message);return;endif
      local_bad=merge(0,1,exchange_plan%workspace_peak_bytes<=huge(workspace_peak_bytes)-workspace_peak_bytes)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid metric exchange receipt overflow';return;endif
      workspace_peak_bytes=workspace_peak_bytes+exchange_plan%workspace_peak_bytes
    endif
    local_real=0d0;if(nowned>0)local_real=maxval(abs(rhs_owned))
    call MPI_Allreduce(local_real,rhs_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_rhs_scale=sqrt(huge(1d0))/(16d0*sqrt(real(n,real64)))
    if(ierr/=MPI_SUCCESS.or.rhs_scale>safe_rhs_scale)local_bad=1
    diagonal=0d0
    do i=1,nowned
      row=int(metric%owned_row_ids(i))
      if(.not.metric%active_rows(row))cycle
      do k=metric%row_offsets(i),metric%row_offsets(i+1)-1
        if(metric%column_ids(k)==row)diagonal(i)=real(metric%values(k))
      enddo
      if(.not.ieee_is_finite(diagonal(i)).or.diagonal(i)<=relative_tolerance*metric%maximum_value)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid metric solver has invalid ownership or diagonal';return;endif
    target=relative_tolerance/max(1d0,metric%condition_estimate)
    solution_owned=(0d0,0d0);fingerprint=metric%fingerprint;maximum_residual=0d0
    do j=1,nrhs
      r=rhs_owned(:,j);z=(0d0,0d0)
      do i=1,nowned
        row=int(metric%owned_row_ids(i))
        if(metric%active_rows(row))then
          z(i)=r(i)/diagonal(i)
        else if(abs(r(i))>relative_tolerance)then
          local_bad=1
        else
          r(i)=(0d0,0d0)
        endif
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='RHS has support outside active metric rank';return;endif
      p=z
      call global_real_dot(r,z,rho,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric residual dot failed';return;endif
      local_real=sum(abs(r)**2);call MPI_Allreduce(local_real,bnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric RHS norm reduction failed';return;endif
      if(bnorm==0d0)then
        relative_residual=0d0;cycle
      endif
      converged=.false.
      do iter=1,max_iterations
        if(present(cached_exchange_plan))then
          call exchange_rt_dg_sparse_values(comm,cached_exchange_plan,p,edge_values,ierr)
        else
          call exchange_rt_dg_sparse_values(comm,exchange_plan,p,edge_values,ierr)
        endif
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric search-vector exchange failed';return;endif
        call sparse_apply(edge_values,ap)
        local_dot=sum(conjg(p)*ap);call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric curvature reduction failed';return;endif
        denominator=real(global_dot)
        if(.not.ieee_is_finite(denominator).or.denominator<=tiny(1d0).or.&
          abs(aimag(global_dot))>100d0*relative_tolerance*max(1d0,abs(denominator)))local_bad=1
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='singular or non-Hermitian metric curvature';return;endif
        alpha=rho/denominator
        if(.not.ieee_is_finite(alpha))local_bad=1
        solution_owned(:,j)=solution_owned(:,j)+alpha*p;r=r-alpha*ap
        if(.not.finite_vector(solution_owned(:,j)).or..not.finite_vector(r))local_bad=1
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='nonfinite hybrid metric iteration';return;endif
        local_real=sum(abs(r)**2);call MPI_Allreduce(local_real,resnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric residual norm reduction failed';return;endif
        relative_residual=sqrt(resnorm/bnorm);iteration_count=max(iteration_count,iter)
        if(relative_residual<=target)then;converged=.true.;exit;endif
        z=(0d0,0d0)
        do i=1,nowned
          row=int(metric%owned_row_ids(i))
          if(metric%active_rows(row))then
            z(i)=r(i)/diagonal(i)
          else if(abs(r(i))>relative_tolerance)then
            local_bad=1
          else
            r(i)=(0d0,0d0)
          endif
        enddo
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
          call cleanup();message='metric iteration leaked into inactive rank';return
        endif
        call global_real_dot(r,z,rho_new,ierr)
        if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(rho_new).or.rho<=tiny(1d0))then
          call cleanup();message='metric preconditioned residual failed';return
        endif
        beta=rho_new/rho
        if(.not.ieee_is_finite(beta))then;call cleanup();message='nonfinite metric recurrence';return;endif
        p=z+beta*p;rho=rho_new
      enddo
      if(.not.converged)then;call cleanup();message='hybrid metric solver iteration cap reached';return;endif
      maximum_residual=max(maximum_residual,relative_residual)
      local_real=0d0;if(nowned>0)local_real=maxval(abs(solution_owned(:,j)))
      call MPI_Allreduce(local_real,solution_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric solution scale reduction failed';return;endif
      quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*relative_tolerance
      if(solution_scale>quantization_limit)then;call cleanup();message='metric solution fingerprint range is unsafe';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),int(j,int64))
      local_hash=0_int64
      do i=1,nowned
        row_hash=ieor(ishftc(int(metric%owned_row_ids(i),int64),11),int(j,int64))
        quantized=nint(real(solution_owned(i,j))/(100d0*relative_tolerance),int64)
        row_hash=ieor(ishftc(row_hash,9),quantized)
        quantized=nint(aimag(solution_owned(i,j))/(100d0*relative_tolerance),int64)
        row_hash=ieor(ishftc(row_hash,9),quantized);local_hash=ieor(local_hash,row_hash)
      enddo
      call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric fingerprint reduction failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),global_hash)
    enddo
    relative_residual=maximum_residual
    if(fingerprint==0_int64)fingerprint=1_int64;ok=.true.
#else
    ok=.false.;message='hybrid metric solver requires MPI';iteration_count=0
    relative_residual=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine sparse_apply(values_by_edge,local_values)
      complex(real64),intent(in)::values_by_edge(:);complex(real64),intent(out)::local_values(:)
      local_values=(0d0,0d0)
      do i=1,nowned
        row=int(metric%owned_row_ids(i));if(.not.metric%active_rows(row))cycle
        do k=metric%row_offsets(i),metric%row_offsets(i+1)-1
          if(metric%active_rows(metric%column_ids(k)))&
            local_values(i)=local_values(i)+metric%values(k)*values_by_edge(k)
        enddo
      enddo
    end subroutine sparse_apply
    subroutine global_real_dot(left,right,value,status)
      complex(real64),intent(in)::left(:),right(:);real(real64),intent(out)::value;integer,intent(out)::status
      local_dot=sum(conjg(left)*right)
      call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,status)
      value=real(global_dot)
      if(status==MPI_SUCCESS)then
        if(abs(aimag(global_dot))>100d0*relative_tolerance*max(1d0,abs(value)))status=MPI_ERR_OTHER
      endif
    end subroutine global_real_dot
    subroutine cleanup()
      if(allocated(solution_owned))deallocate(solution_owned)
      if(allocated(r))deallocate(r)
      if(allocated(z))deallocate(z)
      if(allocated(p))deallocate(p)
      if(allocated(ap))deallocate(ap)
      if(allocated(edge_values))deallocate(edge_values)
      if(allocated(diagonal))deallocate(diagonal)
      if(owns_exchange)call clear_rt_dg_sparse_exchange(exchange_plan)
    end subroutine cleanup
#endif
  end subroutine solve_rt_dg_hybrid_metric

  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  subroutine add_count(total,value,bad)
    integer(int64),intent(inout)::total;integer(int64),intent(in)::value;integer,intent(inout)::bad
    if(bad/=0)return
    if(value<0_int64.or.total>huge(total)-value)then;bad=1;else;total=total+value;endif
  end subroutine add_count
  subroutine add_product(total,left,right,bad)
    integer(int64),intent(inout)::total;integer(int64),intent(in)::left,right;integer,intent(inout)::bad
    if(bad/=0)return
    if(left<0_int64.or.right<0_int64)then;bad=1;return;endif
    if(left/=0_int64)then;if(right>huge(total)/left)then;bad=1;return;endif;endif
    call add_count(total,left*right,bad)
  end subroutine add_product
#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm;integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value;integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module rt_dg_hybrid_metric_solver
